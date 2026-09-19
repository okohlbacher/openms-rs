#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Check that C++ source citations point at the lines whose code they quote.

Every port document, manifest and Rust module cites the source it reproduces as
``<File>.<ext>:<a>`` or ``:<a>-<b>``, and every one of those line numbers was
read off by a human. Package A6 alone shipped four citation defects; one of
them - a range quoted for code ten lines below it - survived two review rounds
in seven places at once, in a round whose own finding was about a wrong line
range. Nothing in the gate battery reads a line number back, so nothing catches
them.

This does. It resolves each citation against the pins the repository already
declares - the core SDK checkouts under ``.reference/``, and the TOPP and CLI
packages read out of their git objects at the pinned revisions, never out of a
working tree - and then, at four levels of evidence:

* a cited range must exist: it may not run backwards or end past the end of
  its file;
* a citation that names its file and one line may not name a blank one;
* a code fragment quoted beside the citation, if the cited file contains it at
  all, must be inside the cited lines - a fragment the file does not contain is
  a paraphrase, a Rust name or a proposed fix, and is ignored;
* a transcribed code block's ``// :NNN`` annotations must match those lines,
  read against the file the block's own entry names.

What makes the third check quiet enough to be worth running is that it asks
only about quotations, and only about the one a citation is written beside. A
backticked *name* - ``writeHeader_``, ``MzTabFile::load`` - is a reference, and
what it names is usually the function the cited lines sit inside rather than
anything on them. A backticked *fragment* with syntax in it reproduces a line.
And a quotation on the far side of a full stop or a semicolon belongs to
another clause, so it is not attached.

*Beside* is meant strictly. Each quotation answers for one citation and each
citation for one quotation, paired over the whole paragraph at once, so that a
paragraph citing one file eight times has each of its citations read against
what stands next to it. A correct citation cannot vouch for a wrong one beside
it - which is the defect this was written for - though a span the same
paragraph cites *around* a line does answer for a quotation of that span.

  python3 tools/check_source_citations.py            # gate
  python3 tools/check_source_citations.py --report   # show the pins and the count
  python3 tools/check_source_citations.py --verbose  # and what could not be checked

Citations it cannot check are counted, never guessed at: a file no reachable
pin contains, a bare ``:a-b`` that fits no file the same paragraph cites, a
file name that more than one pin carries and nothing else narrows, and a
manifest too malformed to parse. ``--report`` also says which pin answered how
many citations and confirmed how many, because a citation confirmed against the
wrong file is worse than an unchecked one.

What this does not catch, stated plainly so that a green run is not read for
more than it says. Of the 2,979 citations it resolves, 97 are confirmed against
code quoted beside them; the rest are checked only for existing, because most
citations in this repository paraphrase the source instead of reproducing it,
and a paraphrase cannot be read back. The A6 defect that this tool was written
for - eight places citing ``:290-293`` for code that sits at ``:280-283`` - is
still not caught, for exactly that reason: those places write
"iter = getFirstChild()" where the source has
``xercesc::DOMNode* firstChild = currentNode->getFirstChild();``. The replay is
on the record and was repeated after every change made here. What the tool does
cover is the class of defect, not that instance: the block in the issue log that
transcribes the same code is read line by line, and shifting it fires.

So the lever that would raise the confirmed fraction is a convention rather than
a cleverer checker - quote the source verbatim in the code span beside the
citation, and this reads it back against the pinned file. Whoever picks this up
next should spend the effort there, and on the population named in the lane's
``not_done``: a bare range written under a file named without a line number,
which is where the issue log puts most of its line numbers, and which needs its
own measurement pass over all the ranges that are currently left unresolved.
"""

import argparse
import collections
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROVENANCE = ROOT / "SOURCE_PROVENANCE.json"
PACKAGE_PINS = ROOT / "tests/data/topp_cli_provenance.json"
ISSUE_LOG = "OpenMS_CPP_ISSUES.md"

SOURCE_EXTENSIONS = ("h", "hpp", "cpp", "cxx")
EXTENSIONS = "|".join(SOURCE_EXTENSIONS)

# A citation names a C++ file - optionally by its full path - and a line or range.
CITATION = re.compile(
    r"(?P<path>(?:[A-Za-z0-9_.-]+/)*)(?P<file>[A-Za-z_][A-Za-z0-9_]*\.(?:" + EXTENSIONS + r"))"
    r":(?P<first>\d+)(?:[-\u2013\u2014](?P<last>\d+))?\b"
)
# A C++ file named on its own, with no line number: how an entry introduces the
# file whose code it goes on to transcribe.
FILE_NAMED = re.compile(
    r"(?P<path>(?:[A-Za-z0-9_.-]+/)*)(?P<file>[A-Za-z_][A-Za-z0-9_]*\.(?:" + EXTENSIONS + r"))\b"
)
# A cheap test for whether a span could hold either of those at all. Both open
# with a long character class, which costs time quadratic in a run of word
# characters that turns out not to end in a source extension - and a fixture
# table can hold a single 84 kB cell of captured stdout. Asking for the
# extension first is linear in the span, and a span without one holds neither.
HAS_SOURCE = re.compile(r"\.(?:" + EXTENSIONS + r")\b")
# A bare range continues the file named before it: "Decoder.cpp:165-168, :179-181".
CONTINUATION = re.compile(r"(?<![\w.:/-]):(?P<first>\d+)(?:[-\u2013\u2014](?P<last>\d+))?\b")
# A transcribed code block annotates its lines: "DOMNode* iter = firstChild;  // :282".
ANNOTATED = re.compile(r"^(?P<code>.*?)\s*//\s*:(?P<line>\d+)\s*$")
# A code span may be broken over two lines by the document's own wrapping.
CODE_SPAN = re.compile(r"`([^`]+?)`")
# A manifest has no code spans, so a quotation in one has to be recognised by
# its shape: a name, or a chain of them, related by an operator to another name
# or call. Requiring a call or a member access on one side keeps "t_wait = 0.2 s"
# - prose about a value - from being read as a line of the source.
NAME = r"[A-Za-z_][\w:.]*(?:\s*(?:->|::|\.)\s*[A-Za-z_][\w:.]*)*(?:\([^)]{0,40}\))?"
BARE_QUOTATION = re.compile(NAME + r"\s*(?:==|!=|<=|>=|<<|=)\s*" + NAME)
IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]{2,}")
# What separates a quotation from a reference: a name, a call or a signature
# says which symbol is meant, while an assignment, a comparison, a statement, a
# stream write or a directive reproduces something that was written on a line.
SYNTAX = re.compile(r"[=;!#]|<<|>>|->|\*|\+|&&|\|\||\b(?:if|for|while|return|switch|case|throw)\b")
# A qualified name, possibly of an operator: "MzTabFile::load",
# "MSChromatogram::operator==". Punctuation makes it look like code without
# making it a quotation of anything.
REFERENCE = re.compile(r"[\w:.~]+(?:operator\s*(?:\[\]|\(\)|[^\w\s]{1,3}|\s+[\w:]+))?")
# A full stop, a semicolon, a colon or an em dash between a quotation and a
# citation ends the clause the citation belongs to, and with it the claim that
# they go together: what follows a colon is a remark about the citation -
# "...:1779-1795: still inside the source's `if (...)`" - and not a quotation of
# those lines. A comma is not one of them: this repository writes a quotation
# and its own citation as "`++window_count`, `Estimator.h:365`". Markup that
# closes after the stop - the bold of a heading, a closing backtick or bracket -
# does not keep it from being one.
BOUNDARY = re.compile(r"[.;:\u2014][*_`)\]]*\s")
# How this repository writes a path inside one of the pinned packages: from the
# package's own root, "OpenMS4-topp/src/FileInfo.cpp", and from the checkout
# that holds them all, "OpenMS4-tests/packages/cli/source/APPLICATIONS/". The
# pin's own paths begin below that prefix, so it has to come off before the path
# can narrow anything - and while it is off, it says which pin is meant.
PACKAGE_PATH = re.compile(r"^(?:OpenMS4-tests/packages/|OpenMS4-)(?P<package>[A-Za-z][\w-]*)/")
REVISION = re.compile(r"\b[0-9a-f]{40}\b")
ISSUE_HEADING = re.compile(r"^##\s+CPP-\d+\b")

# Shorter than this, a fragment such as "f.end)" is not specific enough to
# locate anything, unless it is long enough in words to be a statement.
MIN_QUOTATION = 14

# Everything in the repository that cites the C++: the documents, the manifests
# and the Rust sources, whose module and item documentation cites it too, plus
# the fixture tables and probe scripts under tests/data/, which carry the source
# a row was read off in a column of their own. Not tools/: this checker and its
# tests are full of invented citations of an invented Decoder.cpp, and reading
# them would be reading a fixture as a claim.
DOCUMENT_GLOBS = (
    "*.md", "docs/**/*.md", "docs/**/*.json",
    "tests/data/**/*.md", "tests/data/**/*.json", "tests/data/**/*.tsv", "tests/data/**/*.py",
    "src/**/*.rs", "tests/**/*.rs", "examples/*.rs", "build.rs",
)


class Directory:
    """A pinned source revision unpacked on disk under .reference/."""

    def __init__(self, path):
        self.path = path

    def paths(self):
        return [
            str(item.relative_to(self.path))
            for item in self.path.rglob("*")
            if item.suffix[1:] in SOURCE_EXTENSIONS and item.is_file()
        ]

    def text(self, path):
        return (self.path / path).read_text(errors="replace")

    def __str__(self):
        return str(self.path.relative_to(ROOT) if self.path.is_relative_to(ROOT) else self.path)


class Objects:
    """A pinned source revision read out of a repository's git objects.

    The working tree is never touched: a package checkout may sit at any
    revision, and the pin is the only thing the citations were written against.
    """

    def __init__(self, repository, revision):
        self.repository = repository
        self.revision = revision

    def paths(self):
        listing = git(self.repository, "ls-tree", "-r", "--name-only", self.revision)
        return [p for p in listing.splitlines() if p.rsplit(".", 1)[-1] in SOURCE_EXTENSIONS]

    def text(self, path):
        return git(self.repository, "show", f"{self.revision}:{path}") or ""

    def __str__(self):
        return f"{self.repository.name} {self.revision[:7]}"


def git(repository, *arguments):
    """Run git in a repository, or return None when it is not there or fails."""
    if not repository.is_dir():
        return None
    done = subprocess.run(
        ["git", "-C", str(repository), *arguments],
        capture_output=True, text=True, errors="replace", check=False,
    )
    return done.stdout if done.returncode == 0 else None


def checkout_roots():
    """Where an unpacked pin may sit: beside this checkout, and beside the main one.

    A port lane works in a ``git worktree``, which carries the tracked files but
    not the gitignored ``.reference/`` checkouts, so the main worktree is tried
    too. That keeps the tool usable in the place where citations get written.
    """
    roots = [ROOT]
    common = git(ROOT, "rev-parse", "--git-common-dir")
    if common:
        main = (ROOT / common.strip()).resolve().parent
        if main != ROOT:
            roots.append(main)
    return roots


class Pins:
    """Every pinned revision the repository declares, and the files in each."""

    def __init__(self, reference=None, packages=None):
        self.declared = {"core": json.loads(PROVENANCE.read_text())["target_sdk"]["commit"]}
        self.declared.update(json.loads(PACKAGE_PINS.read_text())["package_revisions"])
        roots = checkout_roots()
        self.reference = pathlib.Path(reference) if reference else first_directory(
            [root / ".reference" for root in roots]
        )
        self.packages = pathlib.Path(packages) if packages else first_directory(
            [root.parent / "OpenMS4-tests" / "packages" for root in roots]
        )
        self.sources = {}
        self.index = collections.defaultdict(list)
        self._lines = {}
        self._collect()

    def _collect(self):
        for path in sorted(self.reference.glob("*")) if self.reference.is_dir() else []:
            top = git(path, "rev-parse", "--show-toplevel")
            if top is None or pathlib.Path(top.strip()).resolve() != path.resolve():
                continue  # Not its own checkout; do not let git walk up to a parent.
            head = git(path, "rev-parse", "HEAD")
            if head:
                self.sources.setdefault(head.strip(), Directory(path))
        for package in ("topp", "cli", "test_data"):
            revision = self.declared.get(package)
            repository = self.packages / package.replace("_", "-")
            if revision and git(repository, "cat-file", "-e", f"{revision}^{{commit}}") is not None:
                self.sources.setdefault(revision, Objects(repository, revision))
        for revision, source in self.sources.items():
            for path in source.paths():
                self.index[(revision, path.rsplit("/", 1)[-1])].append(path)

    def missing(self):
        """The declared pins that are not reachable here."""
        return {name: sha for name, sha in self.declared.items() if sha not in self.sources}

    def paths_named(self, revision, name):
        """Every path in one revision whose file name is this one."""
        return self.index.get((revision, name), [])

    def packaged(self, inside, name, besides):
        """A package pin whose layout this path fits, when none of ``besides`` does.

        The issue log checks each entry against the revisions the entry itself
        names, because its entries predate the current pin - and those are core
        revisions, so a TOPP or CLI path written in an old entry found nothing
        to narrow in and fell back to the core file of the same name. A package
        has exactly one pin in this repository and no older revision of it to
        prefer, so a path that fits its layout names it whatever the entry says:
        ``src/PeakPickerHiRes.cpp:170-186`` is ``doLowMemAlgorithm`` in the TOPP
        tool, and nothing at all in the core algorithm of the same name.
        """
        return [
            (revision, path)
            for revision in dict.fromkeys(self.declared.values())
            if revision not in besides and isinstance(self.sources.get(revision), Objects)
            for path in self.paths_named(revision, name)
            if path.endswith(inside + name)
        ]

    def lines(self, revision, path):
        key = (revision, path)
        if key not in self._lines:
            self._lines[key] = self.sources[revision].text(path).splitlines()
        return self._lines[key]

    def where(self, revision, path):
        return f"{self.sources[revision]}:{path}"

    def label(self, revision):
        """A pin's short name for a tally: what the repository calls it, if anything."""
        for name, sha in sorted(self.declared.items()):
            if sha == revision:
                return f"{name} {revision[:7]}"
        source = self.sources[revision]
        return source.path.name if isinstance(source, Directory) else str(source)


class Unreadable(Exception):
    """A document that cannot be parsed at all, and so cannot be checked."""


def split_package(directory):
    """Split a cited path into the package its prefix names, if any, and the rest.

    ``OpenMS4-topp/src/`` is the TOPP package's own ``src/``; ``FORMAT/`` is
    nobody's package and comes back unchanged.
    """
    match = PACKAGE_PATH.match(directory)
    if not match:
        return None, directory
    return match.group("package").replace("-", "_"), directory[match.end():]


def first_directory(paths):
    for path in paths:
        if path.is_dir():
            return path
    return paths[0]


def units(path, text, default_revisions):
    """Cut a document into spans, each paired with the revisions it may cite.

    Prose - Markdown, and the documentation comments of a Rust module, which
    wrap the same way - is cut at blank lines, at list bullets and at table
    rows, because a citation is quoted by its own bullet or row and not by its
    neighbours. A manifest is cut into its string values, which is exactly how
    its prose is written. The issue log additionally declares a source revision
    per entry, and its entries predate the current pin, so each is checked
    against the revisions that entry itself names.
    """
    if path.suffix == ".json":
        found = []
        try:
            tree = json.loads(text)
        except ValueError as problem:
            raise Unreadable(f"not valid JSON ({problem})") from problem

        def walk(node):
            if isinstance(node, str):
                found.append((node, default_revisions))
            elif isinstance(node, dict):
                for value in node.values():
                    walk(value)
            elif isinstance(node, list):
                for value in node:
                    walk(value)

        walk(tree)
        return found

    if path.suffix == ".tsv":
        # A fixture table: one row is one record, and the column that names the
        # upstream test a row was read off belongs to that row alone. Cutting at
        # rows also keeps CITATION off a whole table at once - thousands of
        # characters of tab-separated values it cannot match, which its optional
        # leading path directory backtracks through at a cost out of all
        # proportion to the handful of citations these files hold.
        return [(line, default_revisions) for line in text.splitlines() if line.strip()]

    issue_log = path.name == ISSUE_LOG
    if issue_log:
        # Entries predate the current pin and name the revision they were read
        # at, in a "Source revision" line, in their evidence, or in both; an
        # entry is checked against every revision it names, and against the
        # log's own header revision when it names none.
        header = tuple(REVISION.findall(text[: text.find("\n## ")]))
        sections = {}
        current_id = None
        for line in text.splitlines():
            if ISSUE_HEADING.match(line):
                current_id = line
                sections[current_id] = []
            if current_id is not None:
                sections[current_id].extend(REVISION.findall(line))
        named = {key: tuple(dict.fromkeys(value)) or header for key, value in sections.items()}
        revisions = header or default_revisions
    else:
        header, named, revisions = (), {}, default_revisions

    spans, current, fenced = [], [], False
    for line in text.splitlines():
        if line.lstrip().startswith("```"):
            fenced = not fenced
            spans.append(("\n".join(current), revisions))
            current = []
            continue
        if fenced:
            continue  # A fenced block quotes without citing.
        if issue_log and ISSUE_HEADING.match(line):
            revisions = named.get(line) or header or default_revisions
        stripped = line.strip()
        alone = stripped.startswith(("|", "#"))
        if alone or not stripped or re.match(r"[-*+]\s|\d+[.)]\s", stripped):
            spans.append(("\n".join(current), revisions))
            current = [] if alone or not stripped else [line]
            if alone:
                spans.append((line, revisions))
        else:
            current.append(line)
    spans.append(("\n".join(current), revisions))
    return [(span, revisions) for span, revisions in spans if span.strip()]


def code_spans(unit):
    """The character ranges of the unit's inline code spans."""
    return [(m.start(1), m.end(1)) for m in CODE_SPAN.finditer(unit)]


def citations_in(unit, quoted_only):
    """The unit's citations: those naming a file, and the bare ranges continuing one.

    ``Decoder.cpp:165-168, :179-181`` writes the second range without repeating
    the file. In prose a bare range is only read as a continuation when it is
    set in a code span, because ``startProgress (:288)`` in a table cell refers
    to a class the cell names without an extension, not to the file cited beside
    it; a manifest has no code spans, so there the shape has to stand on its own.
    """
    named, bare = [], []
    if not HAS_SOURCE.search(unit):
        return [], []
    for match in CITATION.finditer(unit):
        named.append(
            (match.start(), match.group("path"), match.group("file"),
             int(match.group("first")), int(match.group("last") or match.group("first")),
             unit[match.start():match.end()])
        )
    if not named:
        return [], []
    covered = [(start, start + len(text)) for start, _, _, _, _, text in named]
    spans = code_spans(unit) if quoted_only else [(0, len(unit))]
    for match in CONTINUATION.finditer(unit):
        if any(start <= match.start() < end for start, end in covered):
            continue
        if not any(start <= match.start() and match.end() <= end for start, end in spans):
            continue
        bare.append(
            (match.start(), int(match.group("first")),
             int(match.group("last") or match.group("first")), match.group(0))
        )
    return sorted(named), sorted(bare)


def flatten(text):
    """Collapse whitespace so a quotation can be matched across a line break."""
    return re.sub(r"\s+", " ", text).strip()


def quotations(unit, spans=CODE_SPAN):
    """The code fragments the unit reproduces verbatim, with where they sit.

    A name - ``writeHeader_``, ``MzTabFile::load``, ``updateRanges()``,
    ``MSChromatogram::operator==`` - is a reference to a symbol, and what it
    names is almost always the enclosing function or the class the cited lines
    belong to rather than anything written on them, so checking one against a
    line range produces nothing but noise. A fragment that reproduces a line
    instead of naming a symbol - an assignment, a comparison, a statement, a
    stream write, a directive - is a quotation, and a quotation is the thing a
    line number is supposed to point at.

    ``spans`` says how a quotation is delimited: backticks in prose, and in a
    manifest, which has none, the shape of the code itself. Either pattern
    yields the fragment in its last group, or in the whole match when it has no
    groups.
    """
    found, group = [], spans.groups
    for match in spans.finditer(unit):
        if match.group(0).count("\n") > 1:
            continue  # An unpaired backtick, not a span.
        span = flatten(CITATION.sub(" ", match.group(group)))
        if len(span) < MIN_QUOTATION and span.count(" ") < 2:
            continue
        if not SYNTAX.search(span) or not IDENTIFIER.search(span):
            continue
        if REFERENCE.fullmatch(span):
            continue
        if spans is BARE_QUOTATION and not re.search(r"\(|->|::", span):
            continue  # Prose about a value, not a line of the source.
        found.append((match.start(group), match.end(group), span))
    return found


def attach(unit, citations, quoted):
    """Pair each citation with the one quotation it is written beside.

    A paragraph quotes more than it cites: a residual is cited at its own line
    while the constant it reuses, quoted in the same sentence, lives ten lines
    higher and is named without a citation of its own; a sentence later another
    paragraph quotes something else entirely. Attaching every quotation in the
    paragraph to every citation in it would report all of those. A quotation may
    answer for a citation only when it is not separated from it by a full stop,
    a semicolon, a colon or an em dash, which is as far as one clause reaches.

    Among the pairings that leaves, the one taken is the pairing of the whole
    unit: every citation answers for at most one quotation and every quotation
    for at most one citation, as many are paired as can be, and among those the
    closest overall. Both halves matter. A paragraph that cites one file eight
    times and quotes the line of one of them - "fills them under ``#pragma omp
    parallel for`` (``:409``, loop body to ``:487``)" - must give that quotation
    to one citation and not to all eight, or seven correct citations are
    reported. And "``a = f();`` at ``:214-218``, ``b = g();`` at ``:228-232``"
    puts the second quotation four characters from the first citation and its
    own six away, so pairing each citation with whatever is nearest crosses
    them over; pairing the unit as a whole costs 12 that way against 49, and
    reads it as it is written.
    """
    allowed = {}
    for order, (position, length, _) in enumerate(citations):
        finish = position + length
        for index, (start, end, fragment) in enumerate(quoted):
            if end > position and start < finish:
                allowed[(order, index)] = 0
                continue
            between = unit[end:position] if end <= position else unit[finish:start]
            if not BOUNDARY.search(between):
                allowed[(order, index)] = len(between)
    if not allowed:
        return {}
    chosen = pair_up(allowed, len(citations), len(quoted))
    return {citations[order][2]: quoted[index][2] for order, index in chosen}


# Beyond this many quotations in one unit the pairing is taken greedily rather
# than as a whole: the search is exponential in them, and the most any document
# here puts in one paragraph is twelve.
MOST_QUOTATIONS_PAIRED = 14


def pair_up(allowed, citations, quotations):
    """The pairing of a unit: most pairs first, then least distance overall.

    ``allowed`` gives the distance of every pair that may be made at all. The
    search walks the citations in order, carrying for each set of quotations
    already spoken for the best way to have reached it, which is exponential in
    the quotations and linear in the citations. A unit with more quotations
    than :data:`MOST_QUOTATIONS_PAIRED` is paired greedily instead, closest
    pair first, which can pair both fewer of them and *differently*: greedy
    takes a close pair that the whole-unit pairing would have split in order to
    do better overall, so on such a unit a quotation can be attached to a
    citation it was not written beside. Nothing in this repository reaches that
    many - twelve is the most any unit holds - so the greedy branch is a guard
    against a document that grows rather than a path anything here takes, and
    ``test_the_greedy_fallback_may_pair_differently`` pins what it does.
    """
    if quotations > MOST_QUOTATIONS_PAIRED:
        taken, spoken, greedy = set(), set(), []
        for distance, order, index in sorted((d, o, i) for (o, i), d in allowed.items()):
            if order not in taken and index not in spoken:
                taken.add(order)
                spoken.add(index)
                greedy.append((order, index))
        return sorted(greedy)
    states = {0: (0, 0, ())}
    for order in range(citations):
        moves = {}
        for mask, value in states.items():
            for candidate, reached in [(value, mask)] + [
                ((value[0] - 1, value[1] + allowed[(order, index)], value[2] + ((order, index),)),
                 mask | 1 << index)
                for index in range(quotations)
                if not mask >> index & 1 and (order, index) in allowed
            ]:
                if reached not in moves or candidate < moves[reached]:
                    moves[reached] = candidate
        states = moves
    return min(states.values())[2]


def annotations_in(unit, certain=True):
    """The ``code  // :NNN`` line annotations of a transcribed code block.

    ``certain`` says whether the file the block belongs to was cited in the
    block's own span or only inferred from the entry it sits under; an inferred
    file is checked more cautiously, see :func:`problems_with`.
    """
    found = []
    for line in unit.splitlines():
        match = ANNOTATED.match(line)
        if match and match.group("code").strip():
            found.append((int(match.group("line")), match.group("code").strip(), certain))
    return found


def problems_with(pins, revision, path, ranges, quoted, annotated):
    """Everything wrong with one file's citations in one unit, against one revision.

    ``ranges`` are that file's cited ranges, ``quoted`` the quotations attached
    to them as ``(fragment, first, last, text)`` - each read against the one
    range it was written beside - and ``annotated`` the ``// :NNN`` lines of a
    transcribed block. Returns ``(problems, confirmed)``: the caller adds the
    confirmations only once it has settled on a revision, so a candidate that
    is rejected does not leave its count behind.
    """
    lines = pins.lines(revision, path)
    found, confirmed = [], 0
    for first, last, text, _ in ranges:
        if last < first:
            found.append(f"{text}: the range runs backwards")
        elif last > len(lines):
            found.append(f"{text}: the file has {len(lines)} lines")
    if found:
        return found, 0
    for first, last, text, named in ranges:
        # Only for a citation that names its file: a bare range's file is
        # inferred, and a blank line is too weak a signal to report on a guess.
        if named and first == last and not lines[first - 1].strip():
            found.append(f"{text}: that line is blank")
    if found:
        return found, 0
    whole = flatten(" ".join(lines))
    for fragment, first, last, text in sorted(quoted):
        # The lines the quotation may sit on: the ones it is written beside,
        # and any span the same unit cites around them. A document that cites a
        # loop whole - "the banded dynamic program, `SpectrumAlignment.h:79-176`"
        # - and then single lines inside it quotes the loop, not the line, so
        # the enclosing span answers for the quotation as much as the line does.
        # Two ranges that merely sit beside each other do not, which is what
        # keeps a right citation from vouching for a wrong one next to it.
        allowed = [(first, last)] + [
            (a, b) for a, b, _, _ in ranges if a <= first and last <= b and (a, b) != (first, last)
        ]
        if any(fragment in flatten(" ".join(lines[a - 1:b])) for a, b in allowed):
            confirmed += 1
            continue
        if fragment not in whole:
            continue
        found.append(f"{text}: `{fragment}` is in the file but not on the cited lines")
    for number, code, certain in annotated:
        body = flatten(code).rstrip(". ")
        if not certain and body not in whole:
            continue  # The entry's file is not the one this block was taken from.
        if number > len(lines):
            found.append(f"// :{number}: the file has {len(lines)} lines")
        elif body not in flatten(lines[number - 1]):
            found.append(f"// :{number} is `{lines[number - 1].strip()}`, not `{code}`")
        else:
            confirmed += 1
    return found, confirmed


def resolvable(pins, revisions, directory, name):
    """Every (revision, path) a cited file name resolves to, across the revisions.

    Two pins carry the same file name more often than one would think -
    ``FileInfo.cpp`` is a core SDK source *and* a TOPP tool, and fourteen names
    collide that way here - so answering such a citation with whichever pin
    happens to be declared first reads the document against a file it never
    meant. The path the citation writes decides instead, in two steps.

    A path that opens with one of this repository's package prefixes names its
    pin outright, and the citation is resolved in that pin and nowhere else; a
    prefix naming a package with no pin resolves nowhere, and is counted as
    unreachable rather than answered by something else. Otherwise every revision
    whose own layout the rest of the path fits is kept and the revisions it does
    not fit are dropped, instead of each of them falling back to its own file of
    that name. If the path fits none of the revisions offered but does fit a
    package pin outside them, that pin answers - see :meth:`Pins.packaged`. Only
    a name with no usable path, or one whose path fits nothing anywhere,
    resolves everywhere it exists, and :func:`check_file` then has to tell those
    apart or count the citation ambiguous.
    """
    package, inside = split_package(directory)
    if package is not None:
        pinned = pins.declared.get(package)
        revisions = (pinned,) if pinned else ()
    narrowed, anywhere = [], []
    for revision in revisions:
        if revision not in pins.sources:
            continue
        for path in pins.paths_named(revision, name):
            anywhere.append((revision, path))
            if inside and path.endswith(inside + name):
                narrowed.append((revision, path))
    if inside and not narrowed and package is None:
        narrowed = pins.packaged(inside, name, revisions)
    return narrowed or anywhere


def holding(pins, attempts, quoted):
    """The candidates whose own text holds the most of the quotations beside them.

    What tells a tool's ``FileInfo.cpp`` from the SDK's, when the citation wrote
    no path, is the code the document quotes next to it: the pin whose file
    contains that line is the one the citation was written against. Nothing
    quoted, or nothing found in any candidate, settles nothing - the caller then
    keeps every candidate and counts the citation as ambiguous.
    """
    fragments = {fragment for fragment, _, _, _ in quoted}
    if not fragments:
        return []
    scored = collections.defaultdict(list)
    for revision, path in attempts:
        whole = flatten(" ".join(pins.lines(revision, path)))
        scored[sum(fragment in whole for fragment in fragments)].append((revision, path))
    best = max(scored)
    return scored[best] if best else []


def owner_of(pins, revisions, bare, named):
    """Which cited file a bare range continues.

    The nearest citation before it is the usual answer, but documents also write
    a bare range for a file named earlier in the sentence, or - the reason this
    has to be careful - for one they never name with its extension at all. So
    the file is only accepted when the range could exist in it, and a bare range
    that fits nothing the unit cites is left unresolved rather than guessed at.
    """
    position, _, last, _ = bare
    before = [item for item in named if item[0] < position]
    after = [item for item in named if item[0] >= position]
    for _, directory, name, _, _, _ in reversed(before) if before else []:
        for revision, path in resolvable(pins, revisions, directory, name):
            if last <= len(pins.lines(revision, path)):
                return (directory, name)
    for _, directory, name, _, _, _ in after:
        for revision, path in resolvable(pins, revisions, directory, name):
            if last <= len(pins.lines(revision, path)):
                return (directory, name)
    return None


def check_file(pins, revisions, key, cited, quoted, annotated, report):
    """Check one cited file against every revision and path it resolves to.

    A file that no reachable pin contains is counted, not guessed at. Where the
    name still reaches into more than one pin - the same header under two core
    checkouts, or a TOPP tool and an SDK source that share a name - the pin
    whose file holds the code quoted beside the citation answers it; where
    nothing quoted settles it the citation is counted as ambiguous, so that a
    number stands against the chance it was read in the wrong file, instead of
    the first pin declared quietly taking it. Among the candidates that remain,
    one that has nothing wrong with it settles the matter; when they all
    disagree the one that disagrees least is reported, since a citation is
    written against one file, not all of them. Either way the pin that answered
    is tallied, because which pin confirmed how much is the thing a reader has
    to be able to see.
    """
    directory, name = key
    attempts = resolvable(pins, revisions, directory, name)
    if not attempts:
        report["skipped"] += len(cited)
        report["skipped_files"][name] += len(cited)
        return []
    if len({revision for revision, _ in attempts}) > 1:
        attempts = holding(pins, attempts, quoted) or attempts
        if len({revision for revision, _ in attempts}) > 1:
            report["ambiguous"] += len(cited)
            report["ambiguous_files"][name] += len(cited)
    disagreements = []
    for revision, path in attempts:
        wrong, confirmed = problems_with(pins, revision, path, cited, quoted, annotated)
        if not wrong:
            report["checked"] += len(cited)
            report["quoted"] += confirmed
            report["answered"][pins.label(revision)] += len(cited)
            report["confirmed_by"][pins.label(revision)] += confirmed
            return []
        disagreements.append((revision, path, wrong))
    # Every candidate disagrees; report the one that disagrees least.
    revision, path, wrong = min(disagreements, key=lambda item: len(item[2]))
    report["answered"][pins.label(revision)] += len(cited)
    return [(pins.where(revision, path), wrong)]


def check_unit(pins, unit, revisions, quoted_only, report, context=None):
    """Check one unit against the revisions it may cite; return its findings.

    ``context`` is the source file the document named most recently before this
    unit. A transcribed code block is written under the entry that names its
    file and is separated from it by a blank line, so it is a unit of its own
    with no citation in it; without the context its ``// :NNN`` annotations
    would have nothing to be read against, which is how they went unchecked.
    """
    named, bare = citations_in(unit, quoted_only)
    annotated = annotations_in(unit, certain=bool(named))
    if not named:
        if not annotated or context is None:
            return []
        return check_file(pins, revisions, context, [], (), annotated, report)
    ranges = collections.defaultdict(list)
    placed = []
    for index, (position, directory, name, first, last, text) in enumerate(named):
        ranges[(directory, name)].append((first, last, text, True, index))
        placed.append((position, len(text), index))
    for offset, (position, first, last, text) in enumerate(bare):
        owner = owner_of(pins, revisions, (position, first, last, text), named)
        if owner is None:
            report["unresolved"] += 1
            continue
        index = len(named) + offset
        ranges[owner].append((first, last, text, False, index))
        placed.append((position, len(text), index))
    attached = attach(unit, placed, quotations(unit, CODE_SPAN if quoted_only else BARE_QUOTATION))
    findings = []
    for key, cited in ranges.items():
        quoted = [
            (attached[index], first, last, text)
            for first, last, text, _, index in cited if index in attached
        ]
        ranged = [item[:4] for item in cited]
        findings.extend(check_file(pins, revisions, key, ranged, quoted, annotated, report))
        annotated = ()  # One block belongs to one file: the one cited first.
    return findings


def named_file(unit):
    """The last source file the unit names, cited or not, or None.

    An entry names its file once - in the issue log, on an ``Affected
    file/function`` line - and transcribes it further down, so the name has to
    be carried forward from the span that gives it to the span that uses it.
    """
    last = None
    if not HAS_SOURCE.search(unit):
        return None
    for match in FILE_NAMED.finditer(unit):
        last = (match.group("path"), match.group("file"))
    return last


def tally():
    """The counters one run fills in, in one place so that a caller cannot miss one."""
    return {
        "checked": 0, "quoted": 0, "skipped": 0, "unresolved": 0, "ambiguous": 0,
        "skipped_files": collections.Counter(), "ambiguous_files": collections.Counter(),
        "answered": collections.Counter(), "confirmed_by": collections.Counter(),
    }


def documents():
    seen = []
    for pattern in DOCUMENT_GLOBS:
        for path in sorted(ROOT.glob(pattern)):
            if path.is_file() and path not in seen:
                seen.append(path)
    return seen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", action="store_true", help="Print the pins and what was checked")
    parser.add_argument("--verbose", action="store_true", help="Also print what could not be checked")
    parser.add_argument("--reference", help="The directory of pinned core checkouts (default .reference/)")
    parser.add_argument("--packages", help="Where the pinned TOPP and CLI package repositories are")
    parser.add_argument(
        "--require-pins", action="store_true",
        help="Fail if a declared pin is unreachable instead of skipping its citations",
    )
    arguments = parser.parse_args()

    pins = Pins(arguments.reference, arguments.packages)
    default = tuple(pins.declared.values())
    report = tally()
    failures, unreadable = [], []
    for document in documents():
        name = str(document.relative_to(ROOT))
        # Markdown, and a Rust doc comment, which is Markdown; a manifest has no
        # code spans, so there a bare range has to stand on its own shape.
        quoted_only = document.suffix in (".md", ".rs")
        try:
            spans = units(document, document.read_text(errors="replace"), default)
        except Unreadable as problem:
            unreadable.append(f"{name}: {problem}; its citations are skipped")
            continue
        context = None
        for unit, revisions in spans:
            context = named_file(unit) or context
            for where, problems in check_unit(pins, unit, revisions, quoted_only, report, context):
                failures.extend((name, where, problem) for problem in problems)

    missing = pins.missing()
    if missing:
        listed = ", ".join(f"{k} {v[:7]}" for k, v in sorted(missing.items()))
        print(f"Pinned sources not reachable, their citations are skipped: {listed}")
    for line in unreadable:
        print(line)
    for name, where, problem in failures:
        print(f"{name}: {where}\n    {problem}")
    if arguments.verbose and report["skipped_files"]:
        print("Cited files that no reachable pin contains:")
        for name, count in report["skipped_files"].most_common():
            print(f"  {name} ({count})")
    if arguments.verbose and report["ambiguous_files"]:
        print("Cited file names that more than one pin carries, with nothing to tell them apart:")
        for name, count in report["ambiguous_files"].most_common():
            print(f"  {name} ({count})")
    if arguments.report or arguments.verbose:
        print("Pins: " + ", ".join(
            f"{name} {sha[:7]}" + ("" if sha in pins.sources else " (unreachable)")
            for name, sha in sorted(pins.declared.items())
        ))
        # Which pin answered how much: a citation confirmed against the wrong
        # file of the right name is worse than one nothing was checked against,
        # so the split has to be visible and not only the total.
        print("Answered by: " + (", ".join(
            f"{where} {count} ({report['confirmed_by'][where]} confirmed)"
            for where, count in report["answered"].most_common()
        ) or "no pin"))
    summary = (
        f"{report['checked']} citations checked against the pins, "
        f"{report['quoted']} of them confirmed against code quoted beside them; "
        f"{report['skipped']} skipped for an unreachable file, "
        f"{report['ambiguous']} that more than one pin could answer and "
        f"{report['unresolved']} bare ranges left unresolved"
    )
    if unreadable:
        summary += f", and {len(unreadable)} document(s) that could not be read"
    if missing and arguments.require_pins:
        print(f"{summary}; --require-pins was given.")
        return 1
    if failures:
        print(f"\n{len(failures)} citation problem(s). {summary}.")
        return 1
    print(f"{summary}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

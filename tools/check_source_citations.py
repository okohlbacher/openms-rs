#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Check that C++ source citations point at the lines whose code they quote.

Every port document cites the source it reproduces as ``<File>.<ext>:<a>`` or
``:<a>-<b>``, and every one of those line numbers was read by a human. Package
A6 alone shipped four citation defects; one of them - a range quoted ten lines
above the code it described - survived two review rounds in seven places at
once, in a round whose own finding was about a wrong line range. Nothing in the
gate battery reads a line number, so nothing catches them.

This does. For each citation it resolves the file against the pins the
repository already declares - the core SDK checkouts under ``.reference/``, and
the TOPP and CLI packages read out of their git objects at the pinned
revisions, never out of a working tree - and checks that the source the
surrounding text quotes really is inside the cited lines.

What counts as a quotation is decided by the cited file, not by a word list. A
token spelled out in a code span beside the citation is treated as a quotation
only when the cited file contains it at all, and contains it on few enough
lines to locate something; a token the file never contains is prose, a Rust
name or a test name, and is ignored. So the check fires on exactly the defect
it is for: the quoted identifier is in the file the citation names, but not in
the lines the citation gives. It also rejects a range that runs backwards or
ends past the end of its file, and verifies the ``// :NNN`` line annotations
that the issue log's transcribed code blocks carry.

  python3 tools/check_source_citations.py            # gate
  python3 tools/check_source_citations.py --report   # show what was checked
  python3 tools/check_source_citations.py --verbose  # and what could not be
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
    r":(?P<first>\d+)(?:-(?P<last>\d+))?\b"
)
# A bare range continues the file named before it: "Decoder.cpp:165-168, :179-181".
CONTINUATION = re.compile(r"(?<![\w.:/-]):(?P<first>\d+)(?:-(?P<last>\d+))?\b")
# A transcribed code block annotates its lines: "DOMNode* iter = firstChild;  // :282".
ANNOTATED = re.compile(r"^(?P<code>.*?)\s*//\s*:(?P<line>\d+)\s*$")
CODE_SPAN = re.compile(r"`([^`\n]+)`")
IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]{2,}")
# What separates a quotation from a reference: a name, a call or a signature
# says which symbol is meant, while an assignment, a comparison, a statement, a
# stream write or a directive reproduces something that was written on a line.
SYNTAX = re.compile(r"[=;!#]|<<|>>|->|\*|\+|&&|\|\||\b(?:if|for|while|return|switch|case|throw)\b")
# A full stop or a semicolon between a quotation and a citation ends the
# clause the citation belongs to, and with it the claim that they go together.
BOUNDARY = re.compile(r"[.;]\s")
REVISION = re.compile(r"\b[0-9a-f]{40}\b")
ISSUE_HEADING = re.compile(r"^##\s+CPP-\d+\b")

# Shorter than this, a fragment such as "f.end)" is not specific enough to
# locate anything, unless it is long enough in words to be a statement.
MIN_QUOTATION = 14

# Documents and manifests are scanned. The Rust sources are not: their citations
# sit in doc comments whose surrounding code is Rust, so there is nothing there
# that the C++ could be quoted by.
DOCUMENT_GLOBS = ("*.md", "docs/*.md", "docs/*.json", "tests/data/*.json", "tests/data/*/*.json")


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

    def candidates(self, revision, name, directory):
        """Paths in one revision whose file name matches, narrowed by any given path."""
        paths = self.index.get((revision, name), [])
        if directory:
            narrowed = [p for p in paths if p.endswith(directory + name)]
            if narrowed:
                return narrowed
        return paths

    def lines(self, revision, path):
        key = (revision, path)
        if key not in self._lines:
            self._lines[key] = self.sources[revision].text(path).splitlines()
        return self._lines[key]

    def where(self, revision, path):
        return f"{self.sources[revision]}:{path}"


def first_directory(paths):
    for path in paths:
        if path.is_dir():
            return path
    return paths[0]


def units(path, text, default_revisions):
    """Cut a document into spans, each paired with the revisions it may cite.

    Markdown is cut at blank lines, at list bullets and at table rows, because a
    citation is quoted by its own bullet or row and not by its neighbours. A
    manifest is cut into its string values, which is exactly how its prose is
    written. The issue log additionally declares a source revision per entry,
    and entries predate the current pin, so each of its entries is checked
    against the revisions that entry names.
    """
    if path.suffix == ".json":
        found = []

        def walk(node):
            if isinstance(node, str):
                found.append((node, default_revisions))
            elif isinstance(node, dict):
                for value in node.values():
                    walk(value)
            elif isinstance(node, list):
                for value in node:
                    walk(value)

        walk(json.loads(text))
        return found

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
        if issue_log and ISSUE_HEADING.match(line):
            revisions = named.get(line) or header or default_revisions
        stripped = line.strip()
        breaks = (
            not stripped
            or stripped.startswith(("|", "#"))
            or bool(re.match(r"[-*+]\s|\d+[.)]\s", stripped))
        )
        if breaks and not fenced:
            spans.append(("\n".join(current), revisions))
            current = [line] if stripped else []
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


def quotations(unit):
    """The code fragments the unit reproduces verbatim, with where they sit.

    A backticked name - ``writeHeader_``, ``MzTabFile::load``, ``updateRanges()``
    - is a reference to a symbol, and what it names is almost always the
    enclosing function or the class the cited lines belong to rather than
    anything written on them, so checking one against a line range produces
    nothing but noise. A fragment that reproduces a line instead of naming a
    symbol - an assignment, a comparison, a statement, a stream write, a
    directive - is a quotation, and a quotation is the thing a line number is
    supposed to point at.
    """
    found = []
    for match in CODE_SPAN.finditer(unit):
        span = flatten(CITATION.sub(" ", match.group(1)))
        if len(span) < MIN_QUOTATION and span.count(" ") < 2:
            continue
        if not SYNTAX.search(span) or not IDENTIFIER.search(span):
            continue
        found.append((match.start(1), match.end(1), span))
    return found


def attach(unit, citations, quoted):
    """Give each citation the one quotation it is written beside.

    A paragraph quotes more than it cites: a residual is cited at its own line
    while the constant it reuses, quoted in the same sentence, lives ten lines
    higher and is named without a citation of its own; a sentence later another
    paragraph quotes something else entirely. Attaching every quotation in the
    paragraph to every citation in it would report all of those. The quotation a
    citation answers for is the nearest one it is not separated from by a full
    stop or a semicolon, which is as far as one clause reaches.
    """
    attached = collections.defaultdict(set)
    for position, length, key in citations:
        finish = position + length
        near = []
        for start, end, fragment in quoted:
            between = unit[end:position] if end <= position else unit[finish:start]
            if end > position and start < finish:
                near.append((0, fragment))
            elif not BOUNDARY.search(between):
                near.append((len(between), fragment))
        if near:
            attached[key].add(min(near)[1])
    return attached


def annotations_in(unit):
    """The ``code  // :NNN`` line annotations of a transcribed code block."""
    found = []
    for line in unit.splitlines():
        match = ANNOTATED.match(line)
        if match and match.group("code").strip():
            found.append((int(match.group("line")), match.group("code").strip()))
    return found


def compress(numbers):
    """Render line numbers as ranges: [1, 2, 3, 9] -> '1-3, 9'."""
    spans, start, previous = [], numbers[0], numbers[0]
    for number in numbers[1:]:
        if number == previous + 1:
            previous = number
            continue
        spans.append((start, previous))
        start = previous = number
    spans.append((start, previous))
    return ", ".join(str(a) if a == b else f"{a}-{b}" for a, b in spans)


def problems_with(pins, revision, path, ranges, quoted, annotated):
    """Everything wrong with one file's citations in one unit, against one revision."""
    lines = pins.lines(revision, path)
    found = []
    for first, last, text in ranges:
        if last < first:
            found.append(f"{text}: the range runs backwards")
        elif last > len(lines):
            found.append(f"{text}: the file has {len(lines)} lines")
    if found:
        return found
    named = " / ".join(text for _, _, text in ranges)
    whole = flatten(" ".join(lines))
    cited = flatten(" ".join(line for first, last, _ in ranges for line in lines[first - 1:last]))
    for fragment in sorted(quoted):
        if fragment in cited or fragment not in whole:
            continue
        found.append(f"{named}: `{fragment}` is in the file but not on the cited lines")
    for number, code in annotated:
        if number > len(lines):
            found.append(f"// :{number}: the file has {len(lines)} lines")
        elif flatten(code).rstrip(". ") not in flatten(lines[number - 1]):
            found.append(f"// :{number} is `{lines[number - 1].strip()}`, not `{code}`")
    return found


def resolvable(pins, revisions, directory, name):
    """Every (revision, path) a cited file name resolves to, across the revisions."""
    return [
        (revision, path)
        for revision in revisions
        if revision in pins.sources
        for path in pins.candidates(revision, name, directory)
    ]


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


def check_unit(pins, unit, revisions, quoted_only, report):
    """Check one unit against the revisions it may cite; return its findings."""
    named, bare = citations_in(unit, quoted_only)
    if not named:
        return []
    ranges = collections.defaultdict(list)
    placed = []
    for position, directory, name, first, last, text in named:
        ranges[(directory, name)].append((first, last, text))
        placed.append((position, len(text), (directory, name)))
    for position, first, last, text in bare:
        owner = owner_of(pins, revisions, (position, first, last, text), named)
        if owner is None:
            report["unresolved"] += 1
            continue
        ranges[owner].append((first, last, text))
        placed.append((position, len(text), owner))
    quoted = attach(unit, placed, quotations(unit)) if quoted_only else {}
    annotated = annotations_in(unit)
    findings = []
    for key, cited in ranges.items():
        directory, name = key
        attempts = resolvable(pins, revisions, directory, name)
        if not attempts:
            report["skipped"] += len(cited)
            report["skipped_files"][name] += len(cited)
            continue
        disagreements = []
        for revision, path in attempts:
            wrong = problems_with(pins, revision, path, cited, quoted.get(key, ()), annotated)
            if not wrong:
                report["checked"] += len(cited)
                break
            disagreements.append((revision, path, wrong))
        else:
            # Every candidate disagrees; report the one that disagrees least.
            revision, path, wrong = min(disagreements, key=lambda item: len(item[2]))
            findings.append((pins.where(revision, path), wrong))
    return findings


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
    report = {"checked": 0, "skipped": 0, "unresolved": 0, "skipped_files": collections.Counter()}
    failures = []
    for document in documents():
        name = str(document.relative_to(ROOT))
        quoted_only = document.suffix == ".md"
        for unit, revisions in units(document, document.read_text(errors="replace"), default):
            for where, problems in check_unit(pins, unit, revisions, quoted_only, report):
                failures.extend((name, where, problem) for problem in problems)

    missing = pins.missing()
    if missing:
        listed = ", ".join(f"{k} {v[:7]}" for k, v in sorted(missing.items()))
        print(f"Pinned sources not reachable, their citations are skipped: {listed}")
    for name, where, problem in failures:
        print(f"{name}: {where}\n    {problem}")
    if arguments.verbose and report["skipped_files"]:
        print("Cited files that no reachable pin contains:")
        for name, count in report["skipped_files"].most_common():
            print(f"  {name} ({count})")
    if arguments.report or arguments.verbose:
        print("Pins: " + ", ".join(
            f"{name} {sha[:7]}" + ("" if sha in pins.sources else " (unreachable)")
            for name, sha in sorted(pins.declared.items())
        ))
    summary = (
        f"{report['checked']} citations checked against the pins, "
        f"{report['skipped']} skipped for an unreachable file and "
        f"{report['unresolved']} bare ranges left unresolved"
    )
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

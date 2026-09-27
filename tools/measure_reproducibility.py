#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Measure which recorded sha256 digests a third party could recompute, and which CI does.

The repository records sha256 digests all over: manifests pin the C++ source a
port was read from, the fixtures a test compares against and the oracle outputs
it was derived from, and prose repeats them. A digest is only evidence to a
reader who can recompute it. This counts every digest the tracked files record,
works out what each one names, and tries to recompute it the way a third party
who clones this repository could. ``docs/REPRODUCIBILITY.md`` reports the result
and says what it means.

What the text names decides the *family* of an occurrence:

* ``in_repo`` - a tracked file here. Recomputed from the file as it is now
  (``match``), else from every version of that path in the history of ``HEAD``
  (``earlier revision``: the digest was right when it was written and the file
  has moved on since), else ``mismatch``.
* ``upstream`` - a file of the C++ packages (core, cli, topp, test-data), which
  are public repositories at the pinned revisions. Recomputed from their git
  objects, at the revisions the manifest itself names first and the declared
  pins after them.
* ``crate`` - a ``Cargo.lock`` checksum. Cargo recomputes it on every fetch.
* ``oracle`` - an artifact under ``../oracle/``, which is outside this repository
  and published nowhere. Recomputed here when the file is on this machine, which
  says what the lab can recompute, not what a third party can.
* ``build`` - a binary or an install (``libOpenMS.so``, a driver binary, a path
  under ``/ceph``): only an identical build would reproduce it.
* ``log`` - the digest of a gate's or a probe's console log (``log_sha256``,
  ``stderr_sha256``), which was never kept anywhere.
* ``unbound`` - a digest with no artifact the text names that this could find:
  prose, a table cell, a key that names a quantity rather than a file. A digest
  literal in Rust source with no path beside it lands here too; no Rust code in
  this crate computes a sha256, so no Rust test recomputes one.

Independently of what the text names, every digest is also looked up *by
content*: is it the sha256 of a file a third party can get - a file tracked
here now (``tree``), any blob in the history of ``HEAD`` (``history``), or any
file of the four packages at their pins (``pinned``)? A digest the text does not
bind (or binds to the wrong path) but that some public file hashes to is not
lost; it is recomputable by whoever searches for it, and is reported apart from
the ones the text binds.

With ``--ci`` it also runs, in place, every Python gate
``.github/workflows/rust.yml`` runs that can compute a sha256 at all (the gate or
a ``tools/`` module it imports names ``hashlib``; a gate that does not, such as
the one that shells out to cargo, cannot recompute a digest and is not run),
with every sha256 they compute logged. A recorded digest one of them computes,
or a ``Cargo.lock`` checksum, is one CI recomputes on every push. The gates read
nothing outside the checkout, so this is what CI sees.

An occurrence is one digest written in one place; the same pair of artifact and
digest written in four manifests is four occurrences and one distinct pair.
Both are reported. Hex strings that are not digests - the hex-encoded terms of
the controlled-vocabulary projection tables, for one - are not counted: in a
manifest a digest has to sit under a key that says so (``sha256``, ``digest``,
``checksum``, ``hash``, or a key that is itself a file name), and in text the
line has to say so.

    python3 tools/measure_reproducibility.py              # the summary
    python3 tools/measure_reproducibility.py --ci --json  # the numbers docs/REPRODUCIBILITY.md cites
    python3 tools/measure_reproducibility.py --ci --source ../OpenMS4-tests/packages/core
        # and what CI would recompute with check_core_sdk.py --source added (core clean at the pin)
    python3 tools/measure_reproducibility.py --list mismatch   # what did not recompute

This measures; it never fails a build. It needs the package repositories under
``../OpenMS4-tests/packages/`` (beside this checkout or the main one) for the
``upstream`` family and the ``pinned`` index, and says so when they are missing:
``git clone`` the four public repositories there (their URLs are in
``SOURCE_PROVENANCE.json`` and ``tests/data/topp_cli_provenance.json``).
"""

import argparse
import collections
import hashlib
import json
import pathlib
import re
import subprocess
import sys
import tempfile
import threading

ROOT = pathlib.Path(__file__).resolve().parents[1]
HEX = re.compile(r"(?<![0-9A-Za-z])[0-9a-f]{64}(?![0-9A-Za-z])")
WHOLE = re.compile(r"[0-9a-f]{64}")
DIGEST_KEY = re.compile(r"sha256|sha-256|digest|checksum|hash", re.I)
DIGEST_TEXT = re.compile(r"sha-?256|digest|checksum|hash", re.I)
PATH_KEYS = ("path", "file", "header", "fixture", "artifact", "input", "output", "source", "location", "name")
PACKAGES = {"core": "core", "cli": "cli", "topp": "topp", "test_data": "test-data"}
# A manifest names the revisions it read under keys that end in one of these:
# ``source_revision``, ``old_package_revision``, ``commit``.
REVISION_KEY = re.compile(r"(?:^|_)(?:revision|commit)$")
EMPTY = hashlib.sha256(b"").hexdigest()
LOG_KEY = re.compile(r"(?:^|_)(?:log|stdout|stderr)_sha256$")
BUILD_WORDS = re.compile(r"libOpenMS|dylib|\.so\b|binary|executable|/ceph/|product-sdk|install", re.I)
UPSTREAM_PREFIXES = (
    "src/openms/", "src/openswathalgo/", "src/tests/", "src/topp/", "src/utils/", "src/openms_gui/",
    "share/OpenMS/", "cmake/", "source/", "include/", "tests/topp/", "topp/", "contrib/",
    "OpenMS4-core/", "OpenMS4-cli/", "OpenMS4-topp/", "OpenMS4-test-data/", "OpenMS4-tests/", "packages/",
)
PACKAGE_PREFIX = re.compile(r"^(?:OpenMS4-tests/packages/|OpenMS4-|packages/)(core|cli|topp|test-data)/")
WORKFLOW = ROOT / ".github/workflows/rust.yml"
BOUND = {("in_repo", "match"), ("in_repo", "earlier revision"), ("upstream", "match"), ("crate", "cargo")}


def git(repository, *arguments, stdin=None):
    done = subprocess.run(
        ["git", "-C", str(repository), *arguments], input=stdin,
        capture_output=True, check=False,
    )
    return done.stdout if done.returncode == 0 else None


def tracked():
    out = git(ROOT, "ls-files", "-z")
    if out is None:
        sys.exit("cannot list the tracked files")
    return [p for p in out.decode().split("\0") if p]


def package_root():
    roots = [ROOT]
    common = git(ROOT, "rev-parse", "--git-common-dir")
    if common:
        main = (ROOT / common.decode().strip()).resolve().parent
        if main != ROOT:
            roots.append(main)
    for root in roots:
        candidate = root.parent / "OpenMS4-tests" / "packages"
        if candidate.is_dir():
            return candidate, roots
    return None, roots


def blob_digests(repository, blobs):
    """sha256 of each blob id, read in one ``git cat-file --batch`` stream."""
    blobs = list(dict.fromkeys(blobs))
    if not blobs:
        return {}
    process = subprocess.Popen(["git", "-C", str(repository), "cat-file", "--batch"],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    def feed():  # From a thread: git blocks on a full stdout pipe while we would block on stdin.
        process.stdin.write("".join(f"{b}\n" for b in blobs).encode())
        process.stdin.close()

    writer = threading.Thread(target=feed)
    writer.start()
    result = {}
    out = process.stdout
    for blob in blobs:
        header = out.readline().split()
        if len(header) < 3 or header[1] != b"blob":
            continue
        size = int(header[2])
        digest = hashlib.sha256()
        remaining = size
        while remaining:
            chunk = out.read(min(remaining, 1 << 20))
            digest.update(chunk)
            remaining -= len(chunk)
        out.read(1)
        result[blob] = digest.hexdigest()
    writer.join()
    process.wait()
    return result


class Content:
    """Which public bytes hash to a digest: the tree, the history of HEAD, the pinned packages."""

    def __init__(self, files, packages, pins):
        self.tree = collections.defaultdict(set)
        for name in files:
            path = ROOT / name
            if path.is_file():
                self.tree[hashlib.sha256(path.read_bytes()).hexdigest()].add(name)
        self.history = collections.defaultdict(set)
        listing = git(ROOT, "rev-list", "--objects", "HEAD") or b""
        named = {}
        for line in listing.decode().splitlines():
            sha, _, path = line.partition(" ")
            if path:
                named.setdefault(sha, path)
        kinds = git(ROOT, "cat-file", "--batch-check=%(objectname) %(objecttype)",
                    stdin="".join(f"{s}\n" for s in named).encode()) or b""
        blobs = [s for s, kind in (l.split() for l in kinds.decode().splitlines()) if kind == "blob"]
        for blob, digest in blob_digests(ROOT, blobs).items():
            self.history[digest].add(named[blob])
        self.pinned = collections.defaultdict(set)
        if packages is not None:
            for package, revision in pins.items():
                repository = packages / PACKAGES[package]
                listing = git(repository, "ls-tree", "-r", "-z", revision) or b""
                entries = {}
                for entry in listing.decode().split("\0"):
                    meta, _, path = entry.partition("\t")
                    parts = meta.split()
                    if len(parts) == 3 and parts[1] == "blob":
                        entries.setdefault(parts[2], path)
                for blob, digest in blob_digests(repository, entries).items():
                    self.pinned[digest].add(f"{PACKAGES[package]}@{revision[:7]}:{entries[blob]}")

    def where(self, digest):
        if digest == EMPTY:
            return None  # Every empty file hashes to it: finding one says nothing about which was meant.
        if digest in self.tree:
            return "tree"
        if digest in self.history:
            return "history"
        if digest in self.pinned:
            return "pinned"
        return None


class Upstream:
    """Blob digests of the package repositories, read from their git objects."""

    def __init__(self, packages):
        self.packages = packages
        provenance = json.loads((ROOT / "SOURCE_PROVENANCE.json").read_text())
        self.pins = {"core": provenance["target_sdk"]["commit"]}
        self.pins.update(json.loads((ROOT / "tests/data/topp_cli_provenance.json").read_text())["package_revisions"])
        self.cache = {}

    def available(self):
        return self.packages is not None

    def digest(self, package, revision, path):
        key = (package, revision, path)
        if key not in self.cache:
            repository = self.packages / PACKAGES[package]
            blob = git(repository, "cat-file", "blob", f"{revision}:{path}") if self.packages else None
            self.cache[key] = hashlib.sha256(blob).hexdigest() if blob is not None else None
        return self.cache[key]

    def lookup(self, path, revisions, expected):
        """(found, matched): whether any package holds the path at a tried revision, and whether one matches."""
        prefix = PACKAGE_PREFIX.match(path)
        order = ("core", "cli", "topp", "test_data")
        if prefix:
            named = prefix.group(1).replace("-", "_")
            order = (named,) + tuple(p for p in order if p != named)
            path = path[prefix.end():]
        found = False
        for package in order:
            tried = [r for r in revisions if package == "core"] + [self.pins.get(package)]
            for revision in dict.fromkeys(r for r in tried if r):
                digest = self.digest(package, revision, path)
                if digest is None:
                    continue
                found = True
                if digest == expected:
                    return True, True
        return found, False


def path_history(path, cache={}):
    """sha256 of every version of one tracked path in the history of HEAD."""
    if path not in cache:
        raw = git(ROOT, "log", "--full-history", "--format=", "--raw", "--no-abbrev", "--no-renames", "HEAD", "--", path) or b""
        blobs = set()
        for line in raw.decode().splitlines():
            parts = line.split()
            if line.startswith(":") and len(parts) >= 5:
                blobs.update(p for p in parts[2:4] if set(p) != {"0"})
        cache[path] = set(blob_digests(ROOT, blobs).values())
    return cache[path]


def manifest_revisions(tree):
    found = []

    def walk(node):
        if isinstance(node, dict):
            for key, value in node.items():
                if REVISION_KEY.search(key) and isinstance(value, str) and re.fullmatch(r"[0-9a-f]{40}", value):
                    found.append(value)
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)

    walk(tree)
    return list(dict.fromkeys(found))


def looks_like_path(text):
    return bool(re.search(r"[/\\]|\.[A-Za-z0-9]{1,6}$", text)) and " " not in text and len(text) < 300


def json_occurrences(document, tree):
    """(digest, candidate path or None, key, structured) for every digest in a manifest."""
    found = []
    sums_file = "sha256" in document.name.lower()

    def walk(node, key, parent, trail):
        if isinstance(node, dict):
            for k, v in node.items():
                walk(v, k, node, trail + [k])
        elif isinstance(node, list):
            for v in node:
                walk(v, key, parent, trail)
        elif isinstance(node, str):
            if WHOLE.fullmatch(node):
                if key is None:
                    return
                path = None
                if isinstance(parent, dict):
                    stem = re.sub(r"_?sha256$", "", key)
                    # A digest under a prefixed key - ``source_sha256``,
                    # ``derived_from_sha256`` - names the artifact its prefix
                    # names, not the entry's own ``path``: a fixture's
                    # ``source_sha256`` is the upstream file it was derived from.
                    names = (stem, stem + "_path", stem + "_file") if stem and stem != key else PATH_KEYS
                    for name in names:
                        if isinstance(parent.get(name), str) and looks_like_path(parent[name]):
                            path = parent[name]
                            break
                if looks_like_path(key):
                    path = "/".join(trail) if sums_file else key
                if path is None and not (DIGEST_KEY.search(key) or sums_file):
                    return  # A 64-hex string under a key that does not say it is a digest.
                found.append((node, path, key, True))
            else:
                for match in HEX.finditer(node):
                    window = node[max(0, match.start() - 160):match.start()]
                    if DIGEST_TEXT.search(window) or DIGEST_TEXT.search(node[match.end():match.end() + 20]):
                        found.append((match.group(0), None, key, False))

    walk(tree, None, None, [])
    return found


def text_occurrences(document, text):
    found = []
    if document.name == "Cargo.lock":
        for match in re.finditer(r'checksum = "([0-9a-f]{64})"', text):
            found.append((match.group(1), None, "crate", True))
        return found
    sums = document.name in ("SHA256SUMS",) or document.suffix == ".sha256"
    for line in text.splitlines():
        if sums:
            match = re.match(r"([0-9a-f]{64})\s+\*?(.+)$", line.strip())
            if match:
                found.append((match.group(1), match.group(2).strip(), "sums", True))
            continue
        for match in HEX.finditer(line):
            if not DIGEST_TEXT.search(line):
                continue
            path = None
            for token in re.findall(r"[`\"']?([\w./-]+/[\w.-]+\.[A-Za-z0-9]{1,8})[`\"']?", line):
                path = token
                break
            found.append((match.group(0), path, "text", False))
    return found


def resolve(document, path, key, digest, revisions, tracked_set, upstream, roots):
    """The family and the recompute outcome of one occurrence."""
    if key == "crate":
        return "crate", "cargo"
    if path is None:
        if key and BUILD_WORDS.search(key):
            return "build", "not recomputable"
        if key and LOG_KEY.search(key):
            return "log", "not retained"
        return "unbound", "not recomputable"
    raw = path.strip()
    if raw.startswith("/") or BUILD_WORDS.search(raw) and not raw.startswith(("../oracle", "oracle/")):
        if raw.startswith("/") or re.search(r"\.(so|dylib|a)$|/ceph/|product-sdk", raw):
            return "build", "not recomputable"
    if raw.startswith(("../oracle/", "oracle/")):
        local = None
        for root in roots:
            candidate = root.parent / "oracle" / raw.split("oracle/", 1)[1]
            if candidate.is_file():
                local = hashlib.sha256(candidate.read_bytes()).hexdigest() == digest
                break
        return "oracle", {None: "absent here", True: "matches here", False: "differs here"}[local]
    here = document.parent.relative_to(ROOT)
    for candidate in (raw, str(here / raw), raw.lstrip("./")):
        norm = str(pathlib.PurePosixPath(candidate))
        if norm in tracked_set:
            if hashlib.sha256((ROOT / norm).read_bytes()).hexdigest() == digest:
                return "in_repo", "match"
            if digest in path_history(norm):
                return "in_repo", "earlier revision"
            return "in_repo", "mismatch"
    if raw.startswith(UPSTREAM_PREFIXES) or not raw.startswith(("tests/", "docs/", "resources/")):
        if upstream.available():
            found, matched = upstream.lookup(raw, revisions, digest)
            if matched:
                return "upstream", "match"
            if found:
                return "upstream", "mismatch"
        if raw.startswith(UPSTREAM_PREFIXES):
            return "upstream", "not found" if upstream.available() else "packages absent"
    return "unbound", "path not found"


def hashing_gates():
    """The Python gates rust.yml runs that can compute a sha256: the gate or a tools module it imports names hashlib."""
    gates = []
    for line in WORKFLOW.read_text().splitlines():
        match = re.match(r"\s*-\s*run:\s*python3\s+(tools/\S+\.py)(.*)$", line)
        if not match:
            continue
        script, arguments = match.group(1), match.group(2).split()
        seen, pending, hashing = set(), [ROOT / script], False
        while pending:
            source = pending.pop()
            if source in seen or not source.is_file():
                continue
            seen.add(source)
            text = source.read_text()
            hashing = hashing or "hashlib" in text
            for name in re.findall(r"^\s*(?:from\s+(\w+)\s+import|import\s+(\w+))", text, re.M):
                pending.append(source.parent / f"{name[0] or name[1]}.py")
        if hashing:
            gates.append((script, arguments))
    return gates


TRACER = r'''
import hashlib, os, runpy, sys
log_path, script, *arguments = sys.argv[1:]
log = open(log_path, "a")
real_sha256, real_new = hashlib.sha256, hashlib.new
class Traced:
    def __init__(self, inner): self._inner = inner
    def update(self, data): self._inner.update(data)
    def hexdigest(self):
        value = self._inner.hexdigest(); log.write(value + "\n"); return value
    def digest(self):
        log.write(self._inner.hexdigest() + "\n"); return self._inner.digest()
    def copy(self): return Traced(self._inner.copy())
    def __getattr__(self, name): return getattr(self._inner, name)
hashlib.sha256 = lambda *a, **k: Traced(real_sha256(*a, **k))
hashlib.new = lambda name, *a, **k: (Traced(real_new(name, *a, **k))
                                     if name.lower().replace("-", "") == "sha256" else real_new(name, *a, **k))
sys.argv = [script, *arguments]
sys.path.insert(0, os.path.dirname(os.path.abspath(script)))
try:
    runpy.run_path(script, run_name="__main__")
finally:
    log.close()
'''


def ci_digests(gates):
    """Every sha256 these gates compute when run in place, and each gate's exit status."""
    with tempfile.TemporaryDirectory() as scratch:
        tracer = pathlib.Path(scratch) / "tracer.py"
        tracer.write_text(TRACER)
        log = pathlib.Path(scratch) / "digests.log"
        log.touch()
        statuses = {}
        for script, arguments in gates:
            done = subprocess.run([sys.executable, str(tracer), str(log), script, *arguments],
                                  cwd=ROOT, capture_output=True, check=False)
            statuses[" ".join([script, *arguments])] = done.returncode
        return set(log.read_text().split()), statuses


def measure(with_ci=False, source=None):
    packages, roots = package_root()
    upstream = Upstream(packages)
    files = tracked()
    tracked_set = set(files)
    content = Content(files, packages, upstream.pins)
    computed, statuses = ci_digests(hashing_gates()) if with_ci else (None, {})
    projected = None
    if with_ci and source is not None:
        # The one step CI does not take: the SDK check against a clean checkout of the pin.
        extra, more = ci_digests([("tools/check_core_sdk.py", ["--source", str(source)])])
        projected = computed | extra
        statuses.update(more)
    found = []
    for name in files:
        document = ROOT / name
        if document.suffix in (".png", ".gz", ".bz2", ".zip", ".mzML", ".idXML", ".featureXML", ".consensusXML"):
            continue
        try:
            text = document.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        if not HEX.search(text):
            continue
        revisions = []
        if document.suffix == ".json":
            try:
                tree = json.loads(text)
            except ValueError:
                continue
            revisions = manifest_revisions(tree)
            occurrences = json_occurrences(document, tree)
        else:
            occurrences = text_occurrences(document, text)
        found.extend((document, revisions, *occurrence) for occurrence in occurrences)
    # Prose that repeats a digest a manifest records names what the manifest
    # binds it to - the line's first path is often the manifest itself, or a
    # neighbour. The first manifest to bind the digest to a path decides.
    manifest = {}
    for document, revisions, digest, path, key, structured in found:
        if structured and path is not None:
            manifest.setdefault(digest, (document, revisions, path))
    rows = []
    for document, revisions, digest, path, key, structured in found:
        name = str(document.relative_to(ROOT))
        where_bound = document
        if not structured and digest in manifest:
            where_bound, revisions, path = manifest[digest]
        family, outcome = resolve(where_bound, path, key, digest, revisions, tracked_set, upstream, roots)
        where = content.where(digest)
        if (family, outcome) in BOUND:
            tier = "bound"
        elif where:
            tier = "content only"
        elif family in ("oracle", "build"):
            tier = "lab only"
        else:
            tier = "none"
        rows.append({
            "document": name, "digest": digest, "path": path, "key": key,
            "structured": structured, "family": family, "outcome": outcome, "content": where,
            "tier": tier,
            # A gate that hashes zero bytes anywhere computes EMPTY; that recomputes no particular file.
            "ci": None if computed is None else (family == "crate" or (digest in computed and digest != EMPTY)),
            "ci_with_source": None if projected is None else (
                family == "crate" or (digest in projected and digest != EMPTY)),
        })
    return rows, packages, statuses


def count(rows, field):
    occurrences = collections.Counter(r[field] for r in rows)
    distinct = collections.Counter()
    seen = set()
    for r in rows:
        pair = (r[field], r["path"] or r["document"], r["digest"])
        if pair not in seen:
            seen.add(pair)
            distinct[r[field]] += 1
    return occurrences, distinct


def distinct_pairs(rows):
    return len({(r["path"] or r["document"], r["digest"]) for r in rows})


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", action="store_true", help="Print the counts as JSON")
    parser.add_argument("--ci", action="store_true",
                        help="Also run the hashing Python gates of rust.yml in place and count what they recompute")
    parser.add_argument("--source", type=pathlib.Path, metavar="CORE_CHECKOUT",
                        help="With --ci, also count what CI would recompute if it ran check_core_sdk.py --source "
                             "against this clean checkout of the pinned core (a clone of the public repository)")
    parser.add_argument("--list", metavar="VALUE",
                        help="List the occurrences whose family, outcome, content or tier is VALUE")
    arguments = parser.parse_args()
    if arguments.source and not arguments.ci:
        parser.error("--source needs --ci")
    rows, packages, statuses = measure(arguments.ci, arguments.source)
    for row in rows:
        row["family/outcome"] = f"{row['family']}/{row['outcome']}"
    if arguments.list:
        for r in rows:
            if arguments.list in (r["outcome"], r["family"], r["content"], r["tier"]):
                print(f"{r['document']}\t{r['family']}\t{r['outcome']}\t{r['content']}\t{r['key']}\t{r['path']}\t{r['digest'][:12]}")
        return 0
    total, total_distinct = len(rows), distinct_pairs(rows)
    families, families_distinct = count(rows, "family/outcome")
    tiers, tiers_distinct = count(rows, "tier")
    contents, contents_distinct = count(rows, "content")
    ci_rows = [r for r in rows if r["ci"]]
    summary = {
        "occurrences": total,
        "distinct": total_distinct,
        "structured_occurrences": sum(1 for r in rows if r["structured"]),
        "empty_digest": dict(collections.Counter(r["tier"] for r in rows if r["digest"] == EMPTY)),
        "by_family": {k: {"occurrences": families[k], "distinct": families_distinct[k]} for k in sorted(families)},
        "by_tier": {k: {"occurrences": tiers[k], "distinct": tiers_distinct[k]} for k in sorted(tiers)},
        "by_content": {str(k): {"occurrences": contents[k], "distinct": contents_distinct[k]} for k in contents},
        "packages": str(packages) if packages else None,
    }
    if arguments.ci:
        summary["ci"] = {
            "gates": statuses,
            "recomputed": {"occurrences": len(ci_rows), "distinct": distinct_pairs(ci_rows)},
            "recomputed_by_family": dict(collections.Counter(r["family"] for r in ci_rows)),
        }
    projected_rows = [r for r in rows if r["ci_with_source"]]
    if arguments.source:
        summary["ci_with_source"] = {
            "recomputed": {"occurrences": len(projected_rows), "distinct": distinct_pairs(projected_rows)},
            "recomputed_by_family": dict(collections.Counter(r["family"] for r in projected_rows)),
        }
    if arguments.json:
        print(json.dumps(summary, indent=1))
        return 0
    print(f"{total} digest occurrences ({total_distinct} distinct artifact/digest pairs) in the tracked files.")
    for key, value in sorted(summary["by_family"].items(), key=lambda item: -item[1]["occurrences"]):
        print(f"  {key:30} {value['occurrences']:6} occurrences {value['distinct']:6} distinct")
    print("By tier (bound: the text names bytes a third party can get and they recompute; content only: it does "
          "not, but some public file hashes to the digest; lab only: an oracle or a build, no public bytes):")
    for key, value in sorted(summary["by_tier"].items(), key=lambda item: -item[1]["occurrences"]):
        print(f"  {key:30} {value['occurrences']:6} ({100 * value['occurrences'] / total:5.1f}%)"
              f" {value['distinct']:6} distinct ({100 * value['distinct'] / total_distinct:5.1f}%)")
    if arguments.ci:
        failed = [gate for gate, status in statuses.items() if status]
        print(f"Recomputed by CI ({len(statuses)} hashing gates run{', FAILED: ' + ', '.join(failed) if failed else ''}"
              f", plus Cargo.lock): {len(ci_rows)} of {total} occurrences ({100 * len(ci_rows) / total:.1f}%), "
              f"{distinct_pairs(ci_rows)} of {total_distinct} distinct ({100 * distinct_pairs(ci_rows) / total_distinct:.1f}%).")
    if arguments.source:
        print(f"With check_core_sdk.py --source as well: {len(projected_rows)} occurrences "
              f"({100 * len(projected_rows) / total:.1f}%), {distinct_pairs(projected_rows)} distinct "
              f"({100 * distinct_pairs(projected_rows) / total_distinct:.1f}%).")
    if packages is None:
        print("The package repositories were not found: the upstream family and the pinned index are not recomputed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

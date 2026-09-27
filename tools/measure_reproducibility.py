#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Measure which recorded sha256 digests a third party could recompute.

The repository records sha256 digests all over: manifests pin the C++ source a
port was read from, the fixtures a test compares against and the oracle outputs
it was derived from, and prose repeats them. A digest is only evidence to a
reader who can recompute it. This counts every digest the tracked files record,
works out what each one names, and tries to recompute it the way a third party
who clones this repository could:

* ``in_repo`` - the digest names a tracked file here. Recomputed now; the
  result is ``match`` or ``mismatch``.
* ``upstream`` - it names a file of the C++ packages (core, cli, topp,
  test-data), whose repositories are public at the pinned revisions. Recomputed
  here from the package repositories' git objects, at the revisions the
  manifest itself names first and the declared pins after them. A third party
  can do the same after cloning those repositories; nothing in CI does.
* ``crate`` - a ``Cargo.lock`` checksum. Cargo recomputes it on every fetch.
* ``rust_test`` - a digest literal in Rust source. Whether a test recomputes it
  is not decided here; it is counted apart so that it is not claimed.
* ``oracle`` - it names an artifact under ``../oracle/``, which is outside this
  repository and published nowhere. Recomputed here when the file is on this
  machine, which says what the lab can recompute, not what a third party can.
* ``build`` - a binary or an install (``libOpenMS.so``, a driver binary, a
  path under ``/ceph``): only an identical build would reproduce it.
* ``unbound`` - a digest with no artifact the text names that this could find:
  prose, a table cell, a key that names a quantity rather than a file. It may
  well be right; nobody can recompute it from what is written.

An occurrence is one digest written in one place; the same pair of artifact and
digest written in four manifests is four occurrences and one distinct pair.
Both are reported. Hex strings that are not digests - the hex-encoded terms of
the controlled-vocabulary projection tables, for one - are not counted: in a
manifest a digest has to sit under a key that says so (``sha256``, ``digest``,
``checksum``, ``hash``, or a key that is itself a file name), and in text the
line has to say so.

    python3 tools/measure_reproducibility.py            # the summary
    python3 tools/measure_reproducibility.py --json     # the numbers, for docs/REPRODUCIBILITY.md
    python3 tools/measure_reproducibility.py --list mismatch   # what did not recompute

This measures; it never fails a build. It needs the package repositories under
``../OpenMS4-tests/packages/`` (beside this checkout or the main one) for the
``upstream`` family and says so when they are missing.
"""

import argparse
import collections
import hashlib
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
HEX = re.compile(r"(?<![0-9A-Za-z])[0-9a-f]{64}(?![0-9A-Za-z])")
WHOLE = re.compile(r"[0-9a-f]{64}")
DIGEST_KEY = re.compile(r"sha256|sha-256|digest|checksum|hash", re.I)
DIGEST_TEXT = re.compile(r"sha-?256|digest|checksum|hash", re.I)
PATH_KEYS = ("path", "file", "fixture", "artifact", "input", "output", "source", "location", "name")
PACKAGES = {"core": "core", "cli": "cli", "topp": "topp", "test_data": "test-data"}
REVISION_KEYS = ("source_revision", "revision", "reference_revision", "commit", "upstream_revision",
                 "sdk_revision", "core_revision")
BUILD_WORDS = re.compile(r"libOpenMS|dylib|\.so\b|binary|executable|/ceph/|product-sdk|install", re.I)
UPSTREAM_PREFIXES = (
    "src/openms/", "src/openswathalgo/", "src/tests/", "src/topp/", "src/utils/", "src/openms_gui/",
    "share/OpenMS/", "cmake/", "source/", "include/", "tests/topp/", "topp/", "contrib/",
    "OpenMS4-core/", "OpenMS4-cli/", "OpenMS4-topp/", "OpenMS4-test-data/", "OpenMS4-tests/",
)


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
        path = re.sub(r"^(?:OpenMS4-tests/packages/|OpenMS4-)(core|cli|topp|test-data)/", "", path)
        found = False
        for package in ("core", "cli", "topp", "test_data"):
            tried = [r for r in revisions if package == "core"] + [self.pins.get(package)]
            for revision in dict.fromkeys(r for r in tried if r):
                digest = self.digest(package, revision, path)
                if digest is None:
                    continue
                found = True
                if digest == expected:
                    return True, True
        return found, False


def manifest_revisions(tree):
    found = []

    def walk(node):
        if isinstance(node, dict):
            for key, value in node.items():
                if key in REVISION_KEYS and isinstance(value, str) and re.fullmatch(r"[0-9a-f]{40}", value):
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
                    for name in PATH_KEYS:
                        if isinstance(parent.get(name), str) and looks_like_path(parent[name]):
                            path = parent[name]
                            break
                    stem = re.sub(r"_?sha256$", "", key)
                    if path is None and stem != key:
                        for name in (stem, stem + "_path", stem + "_file"):
                            if isinstance(parent.get(name), str) and looks_like_path(parent[name]):
                                path = parent[name]
                                break
                if looks_like_path(key):
                    path = "/".join(t for t in trail if looks_like_path(t) or t == key) if sums_file else key
                    if sums_file:
                        path = "/".join(trail)
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
    if document.suffix == ".rs" and path is None:
        return "rust_test", "not decided"
    if path is None:
        if key and BUILD_WORDS.search(key):
            return "build", "not recomputable"
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
            ok = hashlib.sha256((ROOT / norm).read_bytes()).hexdigest() == digest
            return "in_repo", "match" if ok else "mismatch"
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


def measure():
    packages, roots = package_root()
    upstream = Upstream(packages)
    files = tracked()
    tracked_set = set(files)
    rows = []
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
        for digest, path, key, structured in occurrences:
            family, outcome = resolve(document, path, key, digest, revisions, tracked_set, upstream, roots)
            rows.append({
                "document": name, "digest": digest, "path": path, "key": key,
                "structured": structured, "family": family, "outcome": outcome,
            })
    return rows, packages


def summarise(rows):
    families = collections.Counter((r["family"], r["outcome"]) for r in rows)
    distinct = collections.Counter()
    seen = set()
    for r in rows:
        pair = (r["family"], r["outcome"], r["path"] or r["document"], r["digest"])
        if pair not in seen:
            seen.add(pair)
            distinct[(r["family"], r["outcome"])] += 1
    return families, distinct


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", action="store_true", help="Print the counts as JSON")
    parser.add_argument("--list", metavar="OUTCOME_OR_FAMILY", help="List the occurrences with this outcome or family")
    arguments = parser.parse_args()
    rows, packages = measure()
    families, distinct = summarise(rows)
    third_party = {("in_repo", "match"), ("upstream", "match"), ("crate", "cargo")}
    total = len(rows)
    recomputable = sum(v for k, v in families.items() if k in third_party)
    distinct_total = sum(distinct.values())
    distinct_recomputable = sum(v for k, v in distinct.items() if k in third_party)
    if arguments.list:
        for r in rows:
            if arguments.list in (r["outcome"], r["family"]):
                print(f"{r['document']}\t{r['family']}\t{r['outcome']}\t{r['key']}\t{r['path']}\t{r['digest'][:12]}")
        return 0
    if arguments.json:
        print(json.dumps({
            "occurrences": total,
            "distinct": distinct_total,
            "recomputable_by_a_third_party": {"occurrences": recomputable, "distinct": distinct_recomputable},
            "by_family": {f"{a}/{b}": {"occurrences": families[(a, b)], "distinct": distinct[(a, b)]}
                          for a, b in sorted(families)},
            "packages": str(packages) if packages else None,
        }, indent=1))
        return 0
    print(f"{total} digest occurrences ({distinct_total} distinct artifact/digest pairs) in the tracked files.")
    for (family, outcome), count in sorted(families.items(), key=lambda item: -item[1]):
        print(f"  {family:10} {outcome:18} {count:6} occurrences {distinct[(family, outcome)]:6} distinct")
    print(f"Recomputable by a third party (in_repo match, upstream match, crate): {recomputable} of {total} "
          f"occurrences ({100 * recomputable / total:.1f}%), {distinct_recomputable} of {distinct_total} distinct "
          f"({100 * distinct_recomputable / distinct_total:.1f}%).")
    if packages is None:
        print("The package repositories were not found: the upstream family is not recomputed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

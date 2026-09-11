#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Advance docs/core-sdk-update.json to a new Core SDK target.

The record is updated from the *measured delta* against the current target
rather than rebuilt from scratch. Build registration comes from the CMake
`sources.cmake` union, whose classification is intricate and already recorded
per file; re-deriving it would risk silently reclassifying files that did not
change. So every unchanged path keeps its recorded classification verbatim, and
the tool refuses to run when a registration input changed, because then the
classification genuinely has to be redone by hand.

    python3 tools/core_sdk_retarget.py --target .reference/openms4-core-<short> \
        --review docs/CORE_SDK_<SHORT>_REVIEW.md --write

Historical fixture pins are never rewritten. Carried-forward reference paths
take target hashes; a path whose bytes changed keeps its original hash and needs
an explicit behaviour review, which the tool reports and refuses to invent.
"""

import argparse
import hashlib
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "docs/core-sdk-update.json"
ARCHIVE = ROOT / ".reference/openms4-core"

REGISTRATION_INPUTS = ("sources.cmake", "CMakeLists.txt", "includes.cmake", "OpenSwathAlgoFiles.cmake")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def physical_lines(path):
    data = path.read_bytes()
    return 0 if not data else data.count(b"\n") + (0 if data.endswith(b"\n") else 1)


def scientific_path(path):
    return (
        path.startswith("src/openms/include/OpenMS/") and path.endswith(".h")
        or path.startswith("src/openms/source/") and path.endswith((".cpp", ".h"))
        or path.startswith("src/openswathalgo/include/") and path.endswith(".h")
        or path.startswith("src/openswathalgo/source/") and path.endswith(".cpp")
    )


def scientific_set(root):
    tracked = subprocess.check_output(["git", "-C", str(root), "ls-files"], text=True).splitlines()
    return {p for p in tracked if scientific_path(p)}


def measure(root, path, template):
    """Refresh the measurable fields of a file record, keeping its classification."""
    target = root / path
    record = dict(template)
    record["bytes"] = target.stat().st_size
    record["physical_lines"] = physical_lines(target)
    record["sha256"] = digest(target)
    return record


def summarise(files, root):
    by_kind = {}
    for item in files:
        entry = by_kind.setdefault(item["kind"], {"files": 0, "physical_lines": 0, "bytes": 0})
        entry["files"] += 1
        entry["physical_lines"] += item["physical_lines"]
        entry["bytes"] += item["bytes"]
    tests = sorted((root / "src/tests/class_tests/openms/source").rglob("*.cpp"))
    return {
        "by_kind": {k: by_kind[k] for k in sorted(by_kind)},
        "registration": dict(Counter(item["registration"] for item in files).most_common()),
        "class_test_cpp_files": len(tests),
        "class_test_cpp_lines": sum(physical_lines(p) for p in tests),
    }


def archive_entry(path, template):
    """Build the original-archive side of a delta entry for a newly changed path."""
    source = ARCHIVE / path
    if not source.is_file():
        return None
    record = {k: template[k] for k in ("path", "kind", "domain")}
    record["bytes"] = source.stat().st_size
    record["physical_lines"] = physical_lines(source)
    record["sha256"] = digest(source)
    for key in ("registration", "registration_file", "configuration_note"):
        if key in template:
            record[key] = template[key]
    return record


def rebuild_delta(previous_delta, files, changed_paths, removed_paths, added_paths):
    """Fold newly changed paths into the delta against the original archive."""
    entries = {item["path"]: item for item in previous_delta["changed_or_removed"]}
    current = {item["path"]: item for item in files}
    for path in sorted(changed_paths | added_paths):
        item = current[path]
        if path in entries:
            entries[path]["current"] = item
            entries[path]["status"] = "changed" if entries[path]["old"] else "added"
            continue
        old = archive_entry(path, item)
        entries[path] = {
            "path": path,
            "status": "changed" if old else "added",
            "kind": item["kind"],
            "domain": item["domain"],
            "old": old,
            "current": item,
        }
    for path in sorted(removed_paths):
        entry = entries.get(path)
        if entry and entry["old"]:
            entry["status"] = "removed"
            entry["current"] = None
        elif entry:
            del entries[path]
    changed = [entries[p] for p in sorted(entries)]
    counts = Counter(item["status"] for item in changed)
    unchanged = previous_delta["summary"]["unchanged"] - len(
        [p for p in changed_paths | removed_paths if p not in
         {i["path"] for i in previous_delta["changed_or_removed"]}]
    )
    summary = {"unchanged": unchanged, **{k: counts[k] for k in ("changed", "removed", "added") if counts[k]}}
    domains = defaultdict(lambda: {"added": 0, "removed": 0, "changed": 0, "unchanged": 0})
    for domain, value in previous_delta["by_domain"].items():
        domains[domain] = dict(value)
    for path in changed_paths:
        if path not in {i["path"] for i in previous_delta["changed_or_removed"]}:
            domain = current[path]["domain"]
            domains[domain]["changed"] += 1
            domains[domain]["unchanged"] -= 1
    return {
        "old_summary": previous_delta["old_summary"],
        "summary": summary,
        "by_domain": {d: domains[d] for d in sorted(domains)},
        "changed_or_removed": changed,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--review", default="")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()

    target = args.target.resolve()
    git = lambda *a: subprocess.check_output(["git", "-C", str(target), *a], text=True).strip()
    assert Path(git("rev-parse", "--show-toplevel")).resolve() == target, "Target is not its own Git root"
    assert not git("status", "--porcelain"), "Target checkout must be clean"
    revision = git("rev-parse", "HEAD")

    previous = json.loads(RECORD.read_text())
    old_revision = previous["identity"]["current_package_revision"]
    assert revision != old_revision, "Target already recorded"

    changed_upstream = git("diff", "--name-only", f"{old_revision}..{revision}").splitlines()
    touched_registration = [p for p in changed_upstream if p.endswith(REGISTRATION_INPUTS)]
    if touched_registration:
        raise SystemExit(
            "Build registration inputs changed; the per-file classification must be redone by hand:\n  "
            + "\n  ".join(touched_registration))

    recorded = {item["path"]: item for item in previous["files"]}
    present = scientific_set(target)
    removed = set(recorded) - present
    added = present - set(recorded)
    if added:
        raise SystemExit(
            "New scientific files need a registration classification, which this tool does not derive:\n  "
            + "\n  ".join(sorted(added)))

    files, changed = [], set()
    for path in sorted(present):
        item = measure(target, path, recorded[path])
        if item["sha256"] != recorded[path]["sha256"]:
            changed.add(path)
        files.append(item)

    refs, needs_review = [], []
    for item in previous["carried_forward_references"]:
        path, source_hash = item["path"], item.get("source_sha256", item["sha256"])
        reference = target / path
        if not reference.is_file():
            refs.append(dict(item))
            continue
        target_hash = digest(reference)
        record = {"path": path, "sha256": target_hash, "original_manifests": item["original_manifests"]}
        if target_hash != source_hash:
            record["source_sha256"] = source_hash
            record["review"] = item.get("review", "")
            if item["sha256"] != target_hash:
                needs_review.append(path)
        refs.append(record)

    summary = summarise(files, target)
    identity = dict(previous["identity"])
    identity["current_package_revision"] = revision
    identity["current_checkout_kind"] = "clean detached Git worktree of the packaged Core SDK repository"
    identity["current_worktree_status"] = "clean"
    identity["retrieved_on"] = subprocess.check_output(["date", "+%Y-%m-%d"], text=True).strip()

    record = {
        "schema_version": previous["schema_version"],
        "identity": identity,
        "methodology": previous["methodology"],
        "summary": summary,
        "files": files,
        "registration_evidence": [
            {"path": item["path"], "sha256": digest(target / item["path"])}
            for item in previous["registration_evidence"]
        ],
        "physical_files_by_domain": previous["physical_files_by_domain"],
        "scope_delta": rebuild_delta(previous["scope_delta"], files, changed, removed, added),
        "carried_forward_references": refs,
        "carried_forward_reference_count": len(refs),
        "reference_policy": previous["reference_policy"],
        "target_update_history": previous["target_update_history"] + [previous["previous_target_update"]],
        "previous_target_update": {
            "revision": old_revision,
            "target_revision": revision,
            "comparison": f"https://github.com/okohlbacher/OpenMS4-core/compare/{old_revision}...{revision}",
            "review": args.review,
            "registered_public_header_scope_unchanged":
                summary["registration"].get("registered_public_header")
                == previous["summary"]["registration"].get("registered_public_header"),
            "changed_historical_reference_paths": sum(1 for r in refs if "source_sha256" in r),
        },
    }

    for item in record["registration_evidence"]:
        was = next(e for e in previous["registration_evidence"] if e["path"] == item["path"])
        assert item["sha256"] == was["sha256"], f"Registration evidence changed: {item['path']}"

    print(f"{old_revision[:7]} -> {revision[:7]}")
    print(f"  scientific files: {len(files)}  changed: {sorted(changed)}  removed: {sorted(removed)}")
    print(f"  registered public headers: {summary['registration'].get('registered_public_header')}"
          f" (unchanged: {record['previous_target_update']['registered_public_header_scope_unchanged']})")
    print(f"  class tests: {summary['class_test_cpp_files']} files, {summary['class_test_cpp_lines']} lines")
    if needs_review:
        print("  carried-forward reference bytes changed and need an explicit review:")
        for path in needs_review:
            print("    ", path)
    if args.write:
        RECORD.write_text(json.dumps(record, indent=1) + "\n")
        print(f"  wrote {RECORD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

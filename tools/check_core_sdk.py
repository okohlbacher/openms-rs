#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Check the native port's SDK target; optionally verify a clean source checkout."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def scientific_path(path):
    return (
        path.startswith("src/openms/include/OpenMS/") and path.endswith(".h")
        or path.startswith("src/openms/source/") and path.endswith((".cpp", ".h"))
        or path.startswith("src/openswathalgo/include/") and path.endswith(".h")
        or path.startswith("src/openswathalgo/source/") and path.endswith(".cpp")
    )


def verify(source=None):
    record = json.loads((ROOT / "docs/core-sdk-update.json").read_text())
    identity = record["identity"]
    provenance = json.loads((ROOT / "SOURCE_PROVENANCE.json").read_text())
    target = provenance["target_sdk"]
    revision = identity["current_package_revision"]
    assert re.fullmatch(r"[0-9a-f]{40}", revision), "Invalid SDK revision"
    assert target["commit"] == revision
    assert target["version"] == identity["sdk_version"]
    graph = provenance["identification_graph"]
    assert graph["target_revision"] == revision
    graph_refs = json.loads((ROOT / graph["reference_provenance"]).read_text())
    assert graph_refs["revision"] == revision
    graph_sources = graph_refs["sources"]
    current_sources = []
    for manifest in provenance.get("current_sdk_reference_manifests", []):
        data = json.loads((ROOT / manifest).read_text())
        pins = [data[key] for key in ["revision", "source_revision", "reference_revision", "commit"] if key in data]
        assert pins and all(pin == revision for pin in pins), manifest
        current_sources.extend(data.get("sources", data.get("source_files", [])))
        for item in [*data.get("fixtures", []), *data.get("files", [])]:
            path = Path(item["path"])
            assert not path.is_absolute() and ".." not in path.parts, item["path"]
            assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == item["sha256"], item["path"]
    for item in current_sources:
        path = Path(item["path"])
        assert not path.is_absolute() and ".." not in path.parts, item["path"]
        assert re.fullmatch(r"[0-9a-f]{64}", item["sha256"])
    assert len(graph_sources) == len({item["path"] for item in graph_sources})
    for item in graph_sources:
        assert scientific_path(item["path"]) or item["path"].startswith("src/tests/")
        assert re.fullmatch(r"[0-9a-f]{64}", item["sha256"])
    rust = (ROOT / "src/lib.rs").read_text()
    for name, value in [("CORE_SDK_REVISION", revision), ("CORE_SDK_VERSION", target["version"])]:
        assert re.search(rf'pub const {name}: &str = "([^"]+)";', rust)[1] == value, name

    files = record["files"]
    by_path = {item["path"]: item for item in files}
    assert len(by_path) == len(files), "Duplicate inventory path"
    for item in files:
        assert scientific_path(item["path"]), item["path"]
        assert re.fullmatch(r"[0-9a-f]{64}", item["sha256"])
    for kind, summary in record["summary"]["by_kind"].items():
        selected = [item for item in files if item["kind"] == kind]
        assert len(selected) == summary["files"], kind
        for field in ["bytes", "physical_lines"]:
            assert sum(item[field] for item in selected) == summary[field], (kind, field)
    assert Counter(item["registration"] for item in files) == record["summary"]["registration"]
    refs = record["carried_forward_references"]
    assert len(refs) == len({item["path"] for item in refs}) == record["carried_forward_reference_count"]
    assert len(refs) == target["carried_forward_reference_paths"]
    for item in refs:
        assert not Path(item["path"]).is_absolute()
        for manifest in item["original_manifests"]:
            assert (ROOT / manifest).is_file(), manifest
    changes = record["scope_delta"]["changed_or_removed"]
    for item in changes:
        if item["status"] == "removed":
            assert item["path"] not in by_path
        else:
            assert item["status"] == "changed"
            assert by_path[item["path"]]["sha256"] == item["current"]["sha256"]
    assert Counter(item["status"] for item in changes) == {
        key: value for key, value in record["scope_delta"]["summary"].items() if key != "unchanged"
    }

    if source is not None:
        source = source.resolve()

        def git(*args):
            return subprocess.check_output(["git", "-C", str(source), *args], text=True).strip()

        # Prevent an archive path from accidentally selecting a parent Git repository.
        assert Path(git("rev-parse", "--show-toplevel")).resolve() == source
        assert git("rev-parse", "HEAD") == revision, "Source revision differs from target"
        assert not git("status", "--porcelain"), "Source checkout must be clean"
        tracked = git("ls-files").splitlines()
        assert {path for path in tracked if scientific_path(path)} == set(by_path)
        verified = set()
        for item in [*files, *refs, *record["registration_evidence"], *graph_sources, *current_sources]:
            path = source / item["path"]
            assert hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"], item["path"]
            verified.add(item["path"])
        print(f"Verified {len(verified)} distinct current source/registration/reference files at {revision}.")
    print(f"SDK {target['version']} target, {len(files)} scientific files, {len(refs)} historical reference paths, {len(graph_sources)} graph and {len(current_sources)} added source references agree.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, help="Clean independent checkout of the pinned Core SDK")
    verify(parser.parse_args().source)

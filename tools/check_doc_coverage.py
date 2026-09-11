#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Report and ratchet rustdoc coverage of the crate's public items.

The C++ SDK documents nearly every public member; this port does not yet, and
`missing_docs` cannot simply be switched on without failing the build in
thousands of places. This gate measures coverage per module against a recorded
floor and fails when a module regresses, so documentation can only improve.

    python3 tools/check_doc_coverage.py              # check against the floors
    python3 tools/check_doc_coverage.py --report     # print every module
    python3 tools/check_doc_coverage.py --write      # record current as the floor

A module at 100% is pinned there: it may never lose documentation again.
"""

import argparse
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FLOORS = ROOT / "docs/doc-coverage.json"

ITEM = re.compile(r"\s*pub (?:fn|struct|enum|trait|const|type|mod) ")
ATTRIBUTE = re.compile(r"\s*#\[")


def measure(path):
    """Public items and how many carry a doc comment, ignoring attributes."""
    lines = path.read_text(encoding="utf8", errors="replace").splitlines()
    items = documented = 0
    for index, line in enumerate(lines):
        if not ITEM.match(line):
            continue
        items += 1
        previous = index - 1
        while previous >= 0 and ATTRIBUTE.match(lines[previous]):
            previous -= 1
        if previous >= 0 and lines[previous].strip().startswith("///"):
            documented += 1
    return items, documented


def survey():
    result = {}
    for path in sorted((ROOT / "src").rglob("*.rs")):
        items, documented = measure(path)
        if items:
            result[path.relative_to(ROOT).as_posix()] = {"items": items, "documented": documented}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--write", action="store_true", help="Record current coverage as the new floor")
    parser.add_argument("--report", action="store_true", help="Print per-module coverage")
    args = parser.parse_args()

    current = survey()
    items = sum(v["items"] for v in current.values())
    documented = sum(v["documented"] for v in current.values())
    percent = 100 * documented / items if items else 100.0

    if args.report:
        for name, value in sorted(current.items(), key=lambda kv: (kv[1]["documented"] / kv[1]["items"], -kv[1]["items"])):
            share = 100 * value["documented"] / value["items"]
            print(f"  {share:5.1f}%  {value['documented']:4d}/{value['items']:4d}  {name}")

    if args.write:
        FLOORS.write_text(json.dumps(
            {"note": "Per-module rustdoc coverage floors. Modules may only improve; "
                     "regenerate with tools/check_doc_coverage.py --write after documenting a module.",
             "total": {"items": items, "documented": documented},
             "modules": current}, indent=1) + "\n")
        print(f"Recorded {documented}/{items} = {percent:.1f}% as the floor.")
        return

    if not FLOORS.is_file():
        raise SystemExit("No recorded floors; run with --write once to establish them.")
    floors = json.loads(FLOORS.read_text())["modules"]

    regressed = []
    for name, value in current.items():
        floor = floors.get(name)
        if floor is None:
            # A new module must be fully documented; there is no precedent to regress from.
            if value["documented"] < value["items"]:
                regressed.append(f"{name}: new module documents {value['documented']}/{value['items']}, needs all")
            continue
        was = floor["documented"] / floor["items"]
        now = value["documented"] / value["items"]
        if now < was - 1e-9:
            regressed.append(
                f"{name}: {100 * now:.1f}% ({value['documented']}/{value['items']}) "
                f"below floor {100 * was:.1f}% ({floor['documented']}/{floor['items']})")

    print(f"rustdoc coverage: {documented}/{items} public items = {percent:.1f}%")
    if regressed:
        print("Documentation regressed:")
        for line in regressed:
            print("  ", line)
        raise SystemExit(1)
    recorded = json.loads(FLOORS.read_text())["total"]
    if documented > recorded["documented"] or items != recorded["items"]:
        print("Coverage improved or the surface changed; run --write to record the new floor.")


if __name__ == "__main__":
    main()

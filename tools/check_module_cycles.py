#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Ratchet the crate's cross-module dependency graph.

Rust permits module cycles inside one crate and forbids them between crates,
which is exactly why they accumulate here unnoticed: nothing fails. They matter
because a workspace split - separating the scientific core from the CLI, the
TOPP tools and the heavy optional adapters - cannot happen while they exist.
The count has already grown from 12 two-cycles to 14 while porting continued.

This does not break the existing cycles; unpicking them is its own project. It
freezes the module-pair edges that exist today and fails on a NEW one, so the
graph can only improve. An edge is a pair of top-level `src/` modules, so adding
a file to an existing module costs nothing as long as it reaches for modules
that module already reaches for.

  python3 tools/check_module_cycles.py            # gate
  python3 tools/check_module_cycles.py --report   # show the graph and its cycles
  python3 tools/check_module_cycles.py --write    # re-record, only ever narrowing
"""

import argparse
import collections
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
BASELINE = ROOT / "docs/module-cycles.json"


def modules():
    return {p.name for p in (ROOT / "src").iterdir() if p.is_dir()}


def edges():
    """Map each top-level module to the other top-level modules it names."""
    tops = modules()
    found = collections.defaultdict(set)
    for path in sorted((ROOT / "src").rglob("*.rs")):
        parts = path.relative_to(ROOT / "src").parts
        owner = parts[0] if parts[0] in tops else path.stem
        if owner not in tops:
            continue
        for target in re.findall(r"\bcrate::(\w+)", path.read_text(errors="ignore")):
            if target in tops and target != owner:
                found[owner].add(target)
    return {k: sorted(v) for k, v in sorted(found.items())}


def two_cycles(graph):
    return sorted(
        {tuple(sorted((a, b))) for a, targets in graph.items() for b in targets if a in graph.get(b, ())}
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Re-record the baseline")
    parser.add_argument("--report", action="store_true", help="Print the graph")
    args = parser.parse_args()

    graph = edges()
    pairs = two_cycles(graph)
    flat = {(a, b) for a, targets in graph.items() for b in targets}

    if args.report:
        for owner, targets in graph.items():
            print(f"{owner:<18} -> {', '.join(targets)}")
        print(f"\n{len(flat)} edges, {len(pairs)} mutually-dependent pairs:")
        for a, b in pairs:
            print(f"  {a} <-> {b}")
        return

    if args.write:
        recorded = {"edges": graph, "mutually_dependent_pairs": [list(p) for p in pairs]}
        if BASELINE.exists():
            before = json.loads(BASELINE.read_text())
            was = {(a, b) for a, t in before["edges"].items() for b in t}
            if flat - was:
                print(f"Refusing to widen the baseline with {sorted(flat - was)}")
                return 1
        BASELINE.write_text(json.dumps(recorded, indent=1) + "\n")
        print(f"Recorded {len(flat)} edges and {len(pairs)} mutually-dependent pairs.")
        return 0

    before = json.loads(BASELINE.read_text())
    was = {(a, b) for a, t in before["edges"].items() for b in t}
    added = sorted(flat - was)
    if added:
        print("New cross-module dependencies introduce or deepen cycles:")
        for a, b in added:
            note = " (creates a mutual dependency)" if (a in graph.get(b, ())) else ""
            print(f"  {a} -> {b}{note}")
        print("\nEither route the call through an existing edge, or - if the new edge is")
        print("deliberate - record it with: python3 tools/check_module_cycles.py --write")
        return 1
    removed = len(was) - len(flat)
    print(f"{len(flat)} cross-module edges, {len(pairs)} mutually-dependent pairs"
          + (f"; {removed} edge(s) removed since the baseline - re-record with --write." if removed > 0 else "."))
    return 0


if __name__ == "__main__":
    sys.exit(main())

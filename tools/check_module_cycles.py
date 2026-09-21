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
records the module-pair edges that exist today and fails on any new edge that
closes a cycle, that is, whenever the target module already reaches the source.
A new edge that closes no cycle cannot block the split and is allowed; record it
with --write, which refuses to record a cycle-closing edge. An edge is a pair of
top-level `src/` modules, so adding a file to an existing module costs nothing
as long as it reaches for modules that module already reaches for.

  python3 tools/check_module_cycles.py            # gate
  python3 tools/check_module_cycles.py --report   # show the graph and its cycles
  python3 tools/check_module_cycles.py --write    # re-record; refuses cycle-closing edges
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


def brace_group_end(text, opening):
    """Index of the `}` matching the `{` at `opening`, or the end of `text`."""
    depth = 0
    for index in range(opening, len(text)):
        if text[index] == "{":
            depth += 1
        elif text[index] == "}":
            depth -= 1
            if depth == 0:
                return index
    return len(text)


def group_items(body):
    """Split a brace group's body on the commas that no nested group encloses."""
    items, depth, start = [], 0, 0
    for index, char in enumerate(body):
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
        elif char == "," and depth == 0:
            items.append(body[start:index])
            start = index + 1
    items.append(body[start:])
    return items


def named_after_crate(text):
    """Yield the first path segment of every `crate::` reference in `text`.

    Only that first segment can name a top-level module, so this is the whole
    edge and nothing deeper: `crate::metadata::MetaValue` names `metadata`.

    A braced group used to hide every module it names, because nothing follows
    `crate::` but the `{` and a bare `crate::(\\w+)` match therefore found
    nothing. Two real edges were invisible that way - `param -> metadata`,
    which closed a cycle, and `interfaces -> metadata`, which did not - and an
    edge the gate cannot see is a cycle the gate cannot refuse. Each item of a
    group starts its own path, so `use crate::{Result, metadata::MetaValue}`
    names `Result` and `metadata`, and a group may nest to any depth.
    """
    for match in re.finditer(r"\bcrate::", text):
        rest = text[match.end() :]
        if rest.startswith("{"):
            end = brace_group_end(text, match.end())
            items = group_items(text[match.end() + 1 : end])
        else:
            items = [rest]
        for item in items:
            leading = re.match(r"\s*(\w+)", item)
            if leading:
                yield leading.group(1)


def edges():
    """Map each top-level module to the other top-level modules it names."""
    tops = modules()
    found = collections.defaultdict(set)
    for path in sorted((ROOT / "src").rglob("*.rs")):
        parts = path.relative_to(ROOT / "src").parts
        owner = parts[0] if parts[0] in tops else path.stem
        if owner not in tops:
            continue
        for target in named_after_crate(path.read_text(errors="ignore")):
            if target in tops and target != owner:
                found[owner].add(target)
    return {k: sorted(v) for k, v in sorted(found.items())}


def two_cycles(graph):
    return sorted(
        {tuple(sorted((a, b))) for a, targets in graph.items() for b in targets if a in graph.get(b, ())}
    )


def reaches(graph, start, goal):
    """Whether `goal` can be reached from `start` along recorded edges."""
    seen, stack = set(), [start]
    while stack:
        node = stack.pop()
        if node == goal:
            return True
        if node not in seen:
            seen.add(node)
            stack.extend(graph.get(node, ()))
    return False


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
            closing = [(a, b) for a, b in sorted(flat - was) if reaches(graph, b, a)]
            if closing:
                print(f"Refusing to record edges that close a cycle: {closing}")
                return 1
        BASELINE.write_text(json.dumps(recorded, indent=1) + "\n")
        print(f"Recorded {len(flat)} edges and {len(pairs)} mutually-dependent pairs.")
        return 0

    before = json.loads(BASELINE.read_text())
    was = {(a, b) for a, t in before["edges"].items() for b in t}
    added = sorted(flat - was)
    closing = [(a, b) for a, b in added if reaches(graph, b, a)]
    if closing:
        print("New cross-module dependencies close a cycle (the target already reaches the source):")
        for a, b in closing:
            note = " (mutual dependency)" if a in graph.get(b, ()) else ""
            print(f"  {a} -> {b}{note}")
        print("\nRoute the call through an existing edge instead.")
        return 1
    notes = []
    if added:
        notes.append(f"new acyclic edge(s) {added} not yet recorded - record with --write")
    removed = len(was - flat)
    if removed:
        notes.append(f"{removed} edge(s) removed since the baseline - re-record with --write")
    print(f"{len(flat)} cross-module edges, {len(pairs)} mutually-dependent pairs"
          + ("; " + "; ".join(notes) if notes else "."))
    return 0


if __name__ == "__main__":
    sys.exit(main())

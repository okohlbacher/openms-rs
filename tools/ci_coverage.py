#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Say what the Rust workflow runs, and refuse to let it run less.

A CI change that silently runs less looks exactly like a green one. Nothing
compared one version of `.github/workflows/rust.yml` with the next, so a
dropped `--test`, a narrowed feature slice or a job that stopped running on a
push would have gone unnoticed; and the pre-push sweep
(`~/.local/bin/openms-ci-sweep.sh`) only runs what it can read out of the file,
so a cargo invocation it cannot see is one nobody runs before a push.

This reads the workflow the way the sweep does - every `run:` line of every
step, leading environment assignments stripped - and expands each cargo line
into the units it covers: one per job, runner, toolchain, environment, command,
feature set and target, where a line that names no target covers every target
cargo would select for it (all 332 integration tests for a bare `cargo test`).
Then:

* `--check` refuses a workflow that covers anything but the recorded baseline
  (`tools/ci_coverage_baseline.json`) - less is a drop, and more is a unit a
  later change could drop unseen until it is recorded - that hides a cargo
  invocation from the sweep, or that runs cargo on every push in a job the
  sweep does not sweep;
* `--write` re-records the baseline, and refuses to record a drop unless told
  `--allow-drop`, so a smaller baseline is a line in a reviewed diff
  (`tools/ci_matrix.py --write` records additions itself);
* `--superset OLD NEW` proves one workflow file covers every unit another does;
* `--report` prints the reduced-feature audit of docs/CI_COVERAGE.md.

It is deliberately strict. A construct it does not model - a step `if:` other
than a matrix leg, `continue-on-error`, a matrix `include`, a `cargo test`
filter - is an error, not something to guess about: a guess is how a gate
comes to skip a line.

  python3 tools/ci_coverage.py --check
  python3 tools/ci_coverage.py --report
  python3 tools/ci_coverage.py --superset old-rust.yml .github/workflows/rust.yml
  python3 tools/ci_coverage.py --write [--allow-drop]
"""

import argparse
import collections
import itertools
import json
import pathlib
import re
import shlex
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github/workflows/rust.yml"
BASELINE = ROOT / "tools/ci_coverage_baseline.json"

# The jobs ~/.local/bin/openms-ci-sweep.sh sweeps when it is not told
# otherwise. A job that runs cargo on every push and is not among them is a
# job whose lines nobody runs before the push.
SWEPT_JOBS = ("quality", "portable-feature-graph", "test", "minimum-rust")

# The sweep's own rule for a leading environment assignment, verbatim.
ENV_ASSIGNMENT = re.compile(r"^[A-Za-z_][A-Za-z_0-9]*=(\"[^\"]*\"|'[^']*'|\S*)\s+")
LEG_CONDITION = re.compile(r"^\s*(?:\$\{\{\s*)?matrix\.([A-Za-z_][\w-]*)\s*==\s*'([^']*)'\s*(?:\}\})?\s*$")


class Unmodelled(Exception):
    """The workflow uses something this tool does not model; it will not guess."""


def load_yaml(path):
    try:
        import yaml
    except ImportError:  # pragma: no cover - the environment decides
        raise SystemExit("ci_coverage needs PyYAML (apt: python3-yaml, pip: pyyaml)")
    with open(path) as handle:
        return yaml.safe_load(handle)


# --------------------------------------------------------------------------
# The crate: its integration tests and its features
# --------------------------------------------------------------------------

def test_targets(root=ROOT):
    """Every integration-test target cargo discovers: tests/*.rs and tests/*/main.rs."""
    manifest = (root / "Cargo.toml").read_text()
    if re.search(r"^\s*\[\[test\]\]", manifest, re.M) or re.search(r"^\s*autotests\s*=", manifest, re.M):
        raise Unmodelled("Cargo.toml declares [[test]] targets or autotests; teach test_targets() about them")
    names = {p.stem for p in (root / "tests").glob("*.rs")}
    names |= {p.parent.name for p in (root / "tests").glob("*/main.rs")}
    return sorted(names)


def cargo_features(root=ROOT):
    """The [features] table: feature -> the features and dependencies it enables."""
    text = (root / "Cargo.toml").read_text()
    section = re.search(r"^\[features\]\s*\n(.*?)(?=^\[|\Z)", text, re.S | re.M)
    table = {}
    for match in re.finditer(r"^([A-Za-z0-9_-]+)\s*=\s*\[(.*?)\]", section.group(1), re.M | re.S):
        table[match.group(1)] = re.findall(r'"([^"]+)"', match.group(2))
    return table


def enabled(features, table):
    """The crate features a feature list turns on, with everything they imply."""
    out, todo = set(), list(features)
    while todo:
        feature = todo.pop()
        if feature in out or "/" in feature or feature.startswith("dep:"):
            continue
        out.add(feature)
        todo.extend(table.get(feature, ()))
    return out


# --------------------------------------------------------------------------
# cfg predicates in the tests
# --------------------------------------------------------------------------

def closing_paren(text, opening):
    """Index of the `)` matching the `(` at `opening`."""
    depth = 0
    for index in range(opening, len(text)):
        if text[index] == "(":
            depth += 1
        elif text[index] == ")":
            depth -= 1
            if depth == 0:
                return index
    raise Unmodelled("unbalanced parentheses in a cfg predicate")


def split_top_level(body):
    parts, depth, start = [], 0, 0
    for index, char in enumerate(body):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == "," and depth == 0:
            parts.append(body[start:index])
            start = index + 1
    parts.append(body[start:])
    return [p for p in parts if p.strip()]


def parse_predicate(text):
    """A cfg predicate as ('feature', name) | ('all'|'any', [..]) | ('not', p) | ('other', text)."""
    text = text.strip()
    match = re.fullmatch(r'feature\s*=\s*"([^"]+)"', text)
    if match:
        return ("feature", match.group(1))
    match = re.match(r"(all|any|not)\s*\(", text)
    if match and closing_paren(text, match.end() - 1) == len(text) - 1:
        args = [parse_predicate(p) for p in split_top_level(text[match.end() : -1])]
        if match.group(1) == "not":
            if len(args) != 1:
                raise Unmodelled(f"not() takes one predicate: {text!r}")
            return ("not", args[0])
        return (match.group(1), args)
    return ("other", re.sub(r"\s+", " ", text))


# What the non-feature predicates the tests use are on the Linux runners every
# push job runs on. A predicate that is not here is an error, not a guess.
LINUX_X86_64 = {
    "unix": True,
    "windows": False,
    "test": True,
    "debug_assertions": True,
    'target_os = "linux"': True,
    'target_os = "macos"': False,
    'target_os = "windows"': False,
    'target_arch = "x86_64"': True,
    'target_arch = "aarch64"': False,
}


def holds(predicate, features):
    kind = predicate[0]
    if kind == "feature":
        return predicate[1] in features
    if kind == "all":
        return all(holds(p, features) for p in predicate[1])
    if kind == "any":
        return any(holds(p, features) for p in predicate[1])
    if kind == "not":
        return not holds(predicate[1], features)
    if predicate[1] in LINUX_X86_64:
        return LINUX_X86_64[predicate[1]]
    raise Unmodelled(f"cfg predicate {predicate[1]!r} is not modelled")


def features_named(predicate):
    kind = predicate[0]
    if kind == "feature":
        return {predicate[1]}
    if kind in ("all", "any"):
        return set().union(set(), *(features_named(p) for p in predicate[1]))
    if kind == "not":
        return features_named(predicate[1])
    return set()


def without_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def module_sources(path, seen=None):
    """A test file and every file it pulls in with `#[path = ...] mod x;`."""
    seen = set() if seen is None else seen
    path = path.resolve()
    if path in seen:
        return []
    seen.add(path)
    text = path.read_text()
    found = [text]
    for match in re.finditer(r'#\[path\s*=\s*"([^"]+)"\]\s*(?:pub\s+)?mod\s+\w+\s*;', text):
        found += module_sources(path.parent / match.group(1), seen)
    return found


class TestGates:
    """What a test target's cfgs say about the features it exercises.

    `gate` is its crate-level `#![cfg(...)]`, the predicate without which the
    binary compiles to nothing; `item_features` are the features named by every
    other cfg in it - items, `cfg!` branches, `cfg_attr` - and in the modules it
    pulls in.
    """

    def __init__(self, name, root=ROOT):
        path = root / "tests" / f"{name}.rs"
        if not path.exists():
            path = root / "tests" / name / "main.rs"
        self.name = name
        self.gate = None
        self.item_features = set()
        for index, text in enumerate(module_sources(path)):
            text = without_comments(text)
            if index == 0:
                match = re.search(r"#!\[cfg\(", text)
                if match:
                    end = closing_paren(text, match.end() - 1)
                    self.gate = parse_predicate(text[match.end() : end])
                    text = text[: match.start()] + text[end + 1 :]
            for match in re.finditer(r"(#!?\[cfg\(|#\[cfg_attr\(|\bcfg!\()", text):
                body = text[match.end() : closing_paren(text, match.end() - 1)]
                if match.group(1).startswith("#[cfg_attr"):
                    body = split_top_level(body)[0]
                self.item_features |= features_named(parse_predicate(body))

    def runs_under(self, features):
        """Whether the binary holds any test at all under these enabled features."""
        return self.gate is None or holds(self.gate, features)

    def minimal_feature_sets(self, table):
        """The smallest feature lists that satisfy the gate; [()] when there is none."""
        if self.gate is None:
            return [()]
        named = sorted(features_named(self.gate))
        found = []
        for size in range(len(named) + 1):
            for combination in itertools.combinations(named, size):
                if any(set(smaller) <= set(combination) for smaller in found):
                    continue
                if holds(self.gate, enabled(combination, table)):
                    found.append(combination)
        return found


# --------------------------------------------------------------------------
# The workflow: cargo lines and the units they cover
# --------------------------------------------------------------------------

def strip_assignments(line):
    """(assignments, rest): the sweep's rule for a line's leading VAR=value words."""
    assignments = []
    while ENV_ASSIGNMENT.match(line):
        assignments.append(ENV_ASSIGNMENT.match(line).group(0).strip())
        line = ENV_ASSIGNMENT.sub("", line, count=1)
    return assignments, line


def matrix_combinations(job_name, job):
    strategy = job.get("strategy") or {}
    matrix = strategy.get("matrix") or {}
    if not isinstance(matrix, dict):
        raise Unmodelled(f"{job_name}: a matrix that is not a mapping")
    for key in ("include", "exclude"):
        if key in matrix:
            raise Unmodelled(f"{job_name}: matrix {key} is not modelled")
    keys = sorted(matrix)
    for key in keys:
        if not isinstance(matrix[key], list) or not matrix[key]:
            raise Unmodelled(f"{job_name}: matrix.{key} must be a non-empty list")
    return [dict(zip(keys, values)) for values in itertools.product(*(matrix[k] for k in keys))]


def substitute(text, combination):
    def value(match):
        key = match.group(1)
        if key not in combination:
            raise Unmodelled(f"${{{{ matrix.{key} }}}} names no matrix key")
        return str(combination[key])

    return re.sub(r"\$\{\{\s*matrix\.([A-Za-z_][\w-]*)\s*\}\}", value, text)


def step_runs_in(job_name, step, combination):
    condition = step.get("if")
    if condition is None:
        return True
    match = LEG_CONDITION.match(str(condition))
    if not match:
        raise Unmodelled(f"{job_name}: step if {condition!r} is not a matrix leg condition")
    key, wanted = match.groups()
    if key not in combination:
        raise Unmodelled(f"{job_name}: step if names matrix.{key}, which the job does not define")
    return str(combination[key]) == wanted


class CargoLine:
    """One cargo invocation of the workflow, as the sweep would run it."""

    def __init__(self, job, runner, toolchain, job_condition, env, text):
        self.job = job
        self.runner = runner
        self.toolchain = toolchain
        self.job_condition = job_condition
        self.env = env
        self.text = text


def cargo_lines(spec):
    """Every cargo line of every job, once per matrix combination it runs in."""
    lines = []
    if not isinstance(spec.get("jobs"), dict):
        raise Unmodelled("the workflow has no jobs mapping")
    for job_name, job in spec["jobs"].items():
        for key in ("continue-on-error", "defaults", "container", "services"):
            if key in job:
                raise Unmodelled(f"{job_name}: job-level {key} is not modelled")
        job_condition = str(job.get("if", "")).strip()
        steps = job.get("steps") or []
        toolchains = [s["uses"].split("@", 1)[1] for s in steps if str(s.get("uses", "")).startswith("dtolnay/rust-toolchain@")]
        if len(toolchains) > 1:
            raise Unmodelled(f"{job_name}: more than one rust-toolchain step")
        toolchain = toolchains[0] if toolchains else ""
        combinations = matrix_combinations(job_name, job)
        for index, step in enumerate(steps):
            for key in ("continue-on-error", "working-directory", "shell"):
                if key in step:
                    raise Unmodelled(f"{job_name}: step {key} is not modelled")
            run = step.get("run")
            if run is None:
                continue
            ran_somewhere = False
            for combination in combinations:
                if not step_runs_in(job_name, step, combination):
                    continue
                ran_somewhere = True
                runner = substitute(str(job.get("runs-on", "")), combination)
                step_env = [f"{k}={v}" for k, v in sorted((step.get("env") or {}).items())]
                for raw in str(run).splitlines():
                    raw = substitute(raw.strip(), combination)
                    assignments, rest = strip_assignments(raw)
                    if rest.startswith("cargo "):
                        env = " ".join(sorted(assignments + step_env))
                        lines.append(CargoLine(job_name, runner, toolchain, job_condition, env, rest))
                    elif re.search(r"\bcargo\b", raw) and not re.match(r"\s*#", raw):
                        raise Unmodelled(
                            f"{job_name}: {raw!r} runs cargo in a form the pre-push sweep cannot see; "
                            "write it as its own line that starts with `cargo `"
                        )
            if not ran_somewhere:
                raise Unmodelled(f"{job_name}: step {index} runs in no matrix combination")
    return lines


TARGET_FLAGS = {"--lib", "--bins", "--tests", "--examples", "--benches", "--all-targets", "--doc"}
PLAIN_FLAGS = {"--locked", "--all", "--workspace", "--no-deps", "--frozen", "--offline"}


def parse_cargo(text):
    """(subcommand, command key, feature key, selection) for one cargo line."""
    words = shlex.split(text)
    if words[0] != "cargo":
        raise Unmodelled(f"not a cargo line: {text!r}")
    index = 1
    toolchain = ""
    if words[index].startswith("+"):
        toolchain = words[index][1:]
        index += 1
    subcommand = words[index]
    index += 1
    if subcommand not in ("test", "check", "clippy", "build", "doc", "fmt"):
        raise Unmodelled(f"cargo {subcommand} is not modelled")
    all_features, no_default, features = False, False, []
    tests, kinds, plain, rest = [], set(), [], []
    while index < len(words):
        word = words[index]
        if word == "--":
            rest = words[index + 1 :]
            break
        if word == "--all-features":
            all_features = True
        elif word == "--no-default-features":
            no_default = True
        elif word in ("--features", "-F"):
            index += 1
            features += [f for f in re.split(r"[\s,]+", words[index]) if f]
        elif word.startswith("--features="):
            features += [f for f in re.split(r"[\s,]+", word.split("=", 1)[1]) if f]
        elif word == "--test":
            index += 1
            tests.append(words[index])
        elif word in TARGET_FLAGS:
            kinds.add(word)
        elif word in PLAIN_FLAGS:
            plain.append(word)
        else:
            raise Unmodelled(f"cargo flag {word!r} is not modelled: {text!r}")
        index += 1
    if subcommand == "test" and rest:
        raise Unmodelled(f"a cargo test filter or harness argument is not modelled: {text!r}")
    if all_features:
        feature_key = "ALL"
    else:
        feature_key = " ".join(sorted(set(features) | (set() if no_default else {"default"}))) or "{}"
    command = " ".join([subcommand] + sorted(plain) + (["--"] + rest if rest else []))
    return toolchain, subcommand, command, feature_key, tests, kinds


def selected_targets(subcommand, tests, kinds, all_tests):
    """The targets a line selects, with every integration test named."""
    everything = ["lib", "bins", "examples", "benches"] + [f"test:{t}" for t in all_tests]
    if subcommand == "fmt":
        return ["sources"]
    if subcommand == "doc":
        return ["api-docs"]
    chosen = [f"test:{t}" for t in tests]
    if "--all-targets" in kinds:
        chosen += everything
    for flag, target in (("--lib", "lib"), ("--bins", "bins"), ("--examples", "examples"), ("--benches", "benches")):
        if flag in kinds:
            chosen.append(target)
    if "--tests" in kinds:
        chosen += [f"test:{t}" for t in all_tests]
    if "--doc" in kinds:
        chosen.append("doctests")
    if chosen:
        return sorted(set(chosen))
    if subcommand == "test":
        return sorted(set(["lib", "bins", "examples", "doctests"] + [f"test:{t}" for t in all_tests]))
    return ["bins", "lib"]


def units(lines, all_tests):
    """Every unit the lines cover: (job, when, runner, toolchain, env, command, features, target)."""
    unknown = set()
    covered = set()
    for line in lines:
        toolchain, subcommand, command, feature_key, tests, kinds = parse_cargo(line.text)
        unknown |= set(tests) - set(all_tests)
        for target in selected_targets(subcommand, tests, kinds, all_tests):
            covered.add(
                (line.job, line.job_condition, line.runner, toolchain or line.toolchain,
                 line.env, command, feature_key, target)
            )
    if unknown:
        raise Unmodelled(f"--test names no integration test: {sorted(unknown)}")
    return covered


def check_sweep_visibility(spec, lines):
    """A push job that runs cargo must be one the sweep sweeps."""
    problems = []
    for job in sorted({l.job for l in lines if not l.job_condition}):
        if job not in SWEPT_JOBS:
            problems.append(f"job {job!r} runs cargo on every push, but the pre-push sweep does not sweep it")
    triggers = spec.get("on", spec.get(True))
    if not isinstance(triggers, dict) or "push" not in triggers:
        problems.append("the workflow no longer runs on push")
    return problems


def workflow_units(path, all_tests):
    spec = load_yaml(path)
    lines = cargo_lines(spec)
    return spec, lines, units(lines, all_tests)


# --------------------------------------------------------------------------
# Baseline
# --------------------------------------------------------------------------

FIELDS = ("job", "when", "runner", "toolchain", "env", "command", "features")


def to_record(covered):
    """The baseline's shape: one entry per distinct line context, its targets listed."""
    grouped = collections.defaultdict(list)
    for unit in sorted(covered):
        grouped[unit[:7]].append(unit[7])
    return [dict(zip(FIELDS, key), targets=targets) for key, targets in sorted(grouped.items())]


def from_record(record):
    return {tuple(entry[f] for f in FIELDS) + (target,) for entry in record for target in entry["targets"]}


def describe(missing, limit=40):
    out = []
    record = to_record(missing)
    for entry in record[:limit]:
        context = " | ".join(entry[f] or "-" for f in FIELDS)
        targets = entry["targets"]
        shown = ", ".join(targets[:8]) + (f", ... ({len(targets)} targets)" if len(targets) > 8 else "")
        out.append(f"  {context}: {shown}")
    if len(record) > limit:
        out.append(f"  ... and {len(record) - limit} more groups")
    return "\n".join(out)


def counts(lines, covered):
    tests = {u for u in covered if u[7].startswith("test:")}
    return {
        "cargo lines": len(lines),
        "units": len(covered),
        "(job, features, test) pairs": len({(u[0], u[6], u[7]) for u in tests}),
        "(features, test) pairs": len({(u[6], u[7]) for u in tests}),
    }


# --------------------------------------------------------------------------
# The reduced-feature audit
# --------------------------------------------------------------------------

def audit(lines, covered, all_tests, table, jobs=("test", "minimum-rust")):
    """Per job: which targets run non-empty under a reduced feature set."""
    gates = {name: TestGates(name) for name in all_tests}
    reduced = collections.defaultdict(set)
    for unit in covered:
        if unit[7].startswith("test:") and unit[6] != "ALL":
            reduced[(unit[0], unit[7][5:])].add(unit[6])
    table_rows = {}
    for name in all_tests:
        row = {}
        for job in jobs:
            keys = reduced.get((job, name), set())
            live = sorted(k for k in keys if gates[name].runs_under(enabled(() if k == "{}" else k.split(), table)))
            row[job] = "runs" if live else ("empty" if keys else "none")
        table_rows[name] = row
    return gates, table_rows


def report(path):
    all_tests = test_targets()
    table = cargo_features()
    spec, lines, covered = workflow_units(path, all_tests)
    gates, rows = audit(lines, covered, all_tests, table)
    print(json.dumps(counts(lines, covered)))
    for job in ("test", "minimum-rust"):
        tally = collections.Counter(row[job] for row in rows.values())
        print(f"\n{job}: runs non-empty under a reduced set {tally['runs']}, "
              f"only as an empty binary {tally['empty']}, no reduced line {tally['none']}")
        for state in ("empty", "none"):
            names = [n for n, r in rows.items() if r[job] == state]
            gated = [n for n in names if gates[n].gate is not None]
            items = [n for n in names if gates[n].gate is None and gates[n].item_features]
            plain = [n for n in names if gates[n].gate is None and not gates[n].item_features]
            if names:
                print(f"  {state}: {len(names)} = {len(gated)} crate-gated + {len(items)} with item-level "
                      f"feature cfgs + {len(plain)} feature-agnostic")
                for label, group in (("crate-gated", gated), ("item-level cfgs", items)):
                    for name in group:
                        g = gates[name]
                        what = sorted(features_named(g.gate)) if g.gate else sorted(g.item_features)
                        print(f"    {label:16} {name:45} {' '.join(what)}")
    never = [n for n, r in rows.items() if all(r[j] != "runs" for j in ("test", "minimum-rust"))]
    print(f"\nnever non-empty under a reduced set in either job: {len(never)}")
    for name in never:
        print(f"  {name:45} gate {' '.join(sorted(features_named(gates[name].gate or ('other', ''))))}")


# --------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true", help="gate against the recorded baseline")
    parser.add_argument("--write", action="store_true", help="re-record the baseline")
    parser.add_argument("--allow-drop", action="store_true", help="let --write record a smaller baseline")
    parser.add_argument("--report", action="store_true", help="print the reduced-feature audit")
    parser.add_argument("--superset", nargs=2, metavar=("OLD", "NEW"), help="prove NEW covers every unit of OLD")
    parser.add_argument("--workflow", default=str(WORKFLOW))
    args = parser.parse_args(argv)
    try:
        return run(args)
    except Unmodelled as problem:
        print(f"ci_coverage: {problem}")
        return 2


def record(workflow=WORKFLOW, allow_drop=False):
    """Re-record the baseline from `workflow`; (exit status, message). Refuses a drop."""
    spec, lines, covered = workflow_units(workflow, test_targets())
    problems = check_sweep_visibility(spec, lines)
    if problems:
        return 1, "\n".join(problems)
    if BASELINE.exists():
        missing = from_record(json.loads(BASELINE.read_text())) - covered
        if missing and not allow_drop:
            return 1, (f"Refusing to record a baseline that drops {len(missing)} unit(s):\n{describe(missing)}\n"
                       "Re-run tools/ci_coverage.py --write --allow-drop only if running less is the intent.")
    BASELINE.write_text(json.dumps(to_record(covered), indent=1) + "\n")
    return 0, f"Recorded {len(covered)} units from {len(lines)} cargo lines."


def run(args):
    all_tests = test_targets()
    if args.superset:
        old_spec, old_lines, old = workflow_units(args.superset[0], all_tests)
        new_spec, new_lines, new = workflow_units(args.superset[1], all_tests)
        print("before:", json.dumps(counts(old_lines, old)))
        print("after: ", json.dumps(counts(new_lines, new)))
        missing = old - new
        if missing:
            print(f"NOT a superset: {len(missing)} unit(s) of {args.superset[0]} are not covered:")
            print(describe(missing))
            return 1
        print(f"superset: every one of the {len(old)} units before is covered after "
              f"({len(new - old)} units added)")
        return 0
    if args.report:
        report(args.workflow)
        return 0
    if args.write:
        status, message = record(args.workflow, args.allow_drop)
        print(message)
        return status
    spec, lines, covered = workflow_units(args.workflow, all_tests)
    problems = check_sweep_visibility(spec, lines)
    recorded = from_record(json.loads(BASELINE.read_text()))
    missing, added = recorded - covered, covered - recorded
    if missing:
        problems.append(f"{len(missing)} recorded unit(s) no longer run:\n{describe(missing)}")
    if added:
        # Recorded exactly, so that a unit added today is one a later change
        # cannot drop unseen.
        problems.append(f"{len(added)} unit(s) run that the baseline does not record - record them with "
                        f"python3 tools/ci_matrix.py --write (or tools/ci_coverage.py --write):\n{describe(added, 8)}")
    if problems:
        print("\n".join(problems))
        return 1
    print(f"{len(lines)} cargo lines cover {len(covered)} units, exactly the recorded ones.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""What the coverage gate counts as running, and what it refuses to guess.

The gate exists because a CI change that runs less looks green. These pin the
cases that decide whether it can see a drop: a `--test` taken off a line, a
line no leg runs, a job the sweep does not sweep, a cargo call the sweep cannot
read - and the changes that must NOT count as a drop, such as moving a line to
another leg or folding two lines with the same features into one.
"""

import pathlib
import tempfile
import textwrap
import unittest

import yaml

from ci_coverage import (
    TestGates,
    Unmodelled,
    cargo_lines,
    check_sweep_visibility,
    enabled,
    holds,
    parse_cargo,
    parse_predicate,
    selected_targets,
    units,
)

TESTS = ["alpha", "beta", "gamma"]


def workflow(jobs):
    return yaml.safe_load("on:\n  push:\njobs:\n" + textwrap.indent(textwrap.dedent(jobs), "  "))


def covered(jobs):
    return units(cargo_lines(workflow(jobs)), TESTS)


def job(runs, name="minimum-rust", extra=""):
    steps = "".join(f"    - run: {r}\n" for r in runs)
    return f"{name}:\n  runs-on: ubuntu-latest\n{extra}  steps:\n    - uses: dtolnay/rust-toolchain@1.85.0\n{steps}"


class Selection(unittest.TestCase):
    def test_a_bare_cargo_test_selects_every_integration_test_and_the_doctests(self):
        chosen = selected_targets("test", [], set(), TESTS)
        self.assertEqual(chosen, ["bins", "doctests", "examples", "lib", "test:alpha", "test:beta", "test:gamma"])

    def test_all_targets_adds_benches_and_leaves_out_the_doctests(self):
        chosen = selected_targets("test", [], {"--all-targets"}, TESTS)
        self.assertIn("benches", chosen)
        self.assertNotIn("doctests", chosen)

    def test_a_named_test_selects_only_itself(self):
        self.assertEqual(selected_targets("test", ["beta"], set(), TESTS), ["test:beta"])

    def test_doc_selects_only_the_doctests(self):
        self.assertEqual(selected_targets("test", [], {"--doc"}, TESTS), ["doctests"])


class Features(unittest.TestCase):
    def test_no_default_features_with_a_quoted_list(self):
        _, _, command, key, tests, _ = parse_cargo('cargo test --locked --no-default-features --features "b a" --test x')
        self.assertEqual((command, key, tests), ("test --locked", "a b", ["x"]))

    def test_the_default_set_is_part_of_the_key_when_it_is_not_turned_off(self):
        self.assertEqual(parse_cargo("cargo test --locked --features x")[3], "default x")

    def test_all_features_is_its_own_key(self):
        self.assertEqual(parse_cargo("cargo test --locked --all-features --features x")[3], "ALL")

    def test_no_features_at_all(self):
        self.assertEqual(parse_cargo("cargo test --locked --no-default-features")[3], "{}")

    def test_an_unknown_flag_is_refused_rather_than_ignored(self):
        with self.assertRaises(Unmodelled):
            parse_cargo("cargo test --locked --release")

    def test_a_test_filter_is_refused_because_it_runs_less(self):
        with self.assertRaises(Unmodelled):
            parse_cargo("cargo test --locked -- only_this")


class Drops(unittest.TestCase):
    def test_taking_a_test_off_a_line_is_a_drop(self):
        before = covered(job(["cargo test --locked --no-default-features --test alpha --test beta"]))
        after = covered(job(["cargo test --locked --no-default-features --test alpha"]))
        self.assertEqual({u[-1] for u in before - after}, {"test:beta"})

    def test_dropping_locked_is_a_drop(self):
        before = covered(job(["cargo test --locked --no-default-features --test alpha"]))
        after = covered(job(["cargo test --no-default-features --test alpha"]))
        self.assertTrue(before - after)

    def test_a_different_toolchain_is_a_drop(self):
        before = covered(job(["cargo test --locked --test alpha"]))
        after = covered(job(["cargo +stable test --locked --test alpha"]))
        self.assertTrue(before - after)

    def test_folding_two_lines_with_the_same_features_is_not_a_drop(self):
        before = covered(job(["cargo test --locked --no-default-features --features x --test alpha",
                              "cargo test --locked --no-default-features --features x --test beta"]))
        after = covered(job(['cargo test --locked --no-default-features --features "x" --test beta --test alpha']))
        self.assertEqual(before, after)

    def test_a_bare_line_covers_every_named_test_under_the_same_features(self):
        before = covered(job(["cargo test --locked --no-default-features --test alpha --test gamma"]))
        after = covered(job(["cargo test --locked --no-default-features"]))
        self.assertFalse(before - after)

    def test_a_job_level_if_that_stops_a_push_job_is_a_drop(self):
        before = covered(job(["cargo test --locked --test alpha"]))
        after = covered(job(["cargo test --locked --test alpha"], extra="  if: github.event_name == 'workflow_dispatch'\n"))
        self.assertTrue(before - after)

    def test_a_runner_taken_out_of_the_matrix_is_a_drop(self):
        two = "  strategy:\n    matrix:\n      os: [macos-latest, windows-latest]\n"
        before = covered(job(["cargo test --locked --test alpha"], extra=two).replace("ubuntu-latest", "${{ matrix.os }}"))
        one = "  strategy:\n    matrix:\n      os: [macos-latest]\n"
        after = covered(job(["cargo test --locked --test alpha"], extra=one).replace("ubuntu-latest", "${{ matrix.os }}"))
        self.assertEqual({u[2] for u in before - after}, {"windows-latest"})

    def test_a_job_environment_is_part_of_what_a_line_runs(self):
        before = covered(job(["cargo test --locked --test alpha"]))
        after = covered(job(["cargo test --locked --test alpha"], extra="  env:\n    RUSTFLAGS: -C x\n"))
        self.assertTrue(before - after)

    def test_an_environment_assignment_is_kept_apart_as_the_sweep_keeps_it(self):
        lines = cargo_lines(workflow(job(['RUSTFLAGS="-C x" cargo check --locked --all-features'])))
        self.assertEqual((lines[0].env, lines[0].text), ('RUSTFLAGS="-C x"', "cargo check --locked --all-features"))


LEGS = "  strategy:\n    matrix:\n      leg: [a, b]\n"


def legged(steps):
    body = "".join(f"    - if: matrix.leg == '{leg}'\n      run: {run}\n" for leg, run in steps)
    return f"minimum-rust:\n  runs-on: ubuntu-latest\n{LEGS}  steps:\n{body}"


class Legs(unittest.TestCase):
    def test_moving_a_line_to_another_leg_is_not_a_change(self):
        before = covered(legged([("a", "cargo test --locked --test alpha"), ("b", "cargo test --locked --test beta")]))
        after = covered(legged([("b", "cargo test --locked --test alpha"), ("a", "cargo test --locked --test beta")]))
        self.assertEqual(before, after)

    def test_each_legged_line_runs_once_not_once_per_leg(self):
        lines = cargo_lines(workflow(legged([("a", "cargo test --locked --test alpha"), ("b", "cargo test --locked --test beta")])))
        self.assertEqual(len(lines), 2)

    def test_a_line_no_leg_runs_is_refused(self):
        with self.assertRaises(Unmodelled):
            covered(legged([("a", "cargo test --locked --test alpha"), ("c", "cargo test --locked --test beta")]))

    def test_any_other_step_condition_is_refused(self):
        steps = "minimum-rust:\n  runs-on: ubuntu-latest\n  steps:\n    - if: github.ref == 'refs/heads/main'\n      run: cargo test --locked\n"
        with self.assertRaises(Unmodelled):
            covered(steps)

    def test_continue_on_error_is_refused_because_a_failure_would_stop_counting(self):
        steps = "minimum-rust:\n  runs-on: ubuntu-latest\n  steps:\n    - continue-on-error: true\n      run: cargo test --locked\n"
        with self.assertRaises(Unmodelled):
            covered(steps)


class Sweep(unittest.TestCase):
    def test_a_cargo_call_the_sweep_cannot_read_is_refused(self):
        with self.assertRaises(Unmodelled):
            covered(job(["cd sub && cargo test --locked"]))

    def test_a_push_job_the_sweep_does_not_sweep_is_reported(self):
        spec = workflow(job(["cargo test --locked"], name="new-slices"))
        self.assertTrue(check_sweep_visibility(spec, cargo_lines(spec)))

    def test_a_tag_only_job_outside_the_sweep_is_fine(self):
        spec = workflow(job(["cargo test --locked"], name="cross-platform", extra="  if: startsWith(github.ref, 'refs/tags/')\n"))
        self.assertEqual(check_sweep_visibility(spec, cargo_lines(spec)), [])

    def test_a_workflow_that_stops_running_on_push_is_reported(self):
        spec = workflow(job(["cargo test --locked"]))
        spec[True] = {"workflow_dispatch": None}
        self.assertTrue(check_sweep_visibility(spec, cargo_lines(spec)))

    def test_the_list_form_of_on_is_read(self):
        spec = workflow(job(["cargo test --locked"]))
        spec[True] = ["push", "pull_request"]
        self.assertEqual(check_sweep_visibility(spec, cargo_lines(spec)), [])

    def test_the_swept_jobs_are_fine(self):
        spec = workflow(job(["cargo test --locked"], name="minimum-rust"))
        self.assertEqual(check_sweep_visibility(spec, cargo_lines(spec)), [])


TABLE = {"mzml": ["numpress", "dep:x"], "numpress": [], "paramxml": [], "featurexml": [], "sqmass": ["sqlite", "mzml"], "sqlite": []}


class Gates(unittest.TestCase):
    def test_implied_features_are_enabled(self):
        self.assertEqual(enabled(["sqmass"], TABLE), {"sqmass", "sqlite", "mzml", "numpress"})

    def test_all_needs_every_feature(self):
        gate = parse_predicate('all(feature = "mzml", feature = "paramxml")')
        self.assertFalse(holds(gate, {"mzml"}))
        self.assertTrue(holds(gate, {"mzml", "paramxml"}))

    def test_an_unmodelled_predicate_is_refused(self):
        with self.assertRaises(Unmodelled):
            holds(parse_predicate('target_os = "freebsd"'), set())

    def smallest(self, gate):
        probe = TestGates.__new__(TestGates)
        probe.gate = parse_predicate(gate) if gate else None
        return probe.minimal_feature_sets(TABLE)

    def test_no_gate_is_the_empty_set(self):
        self.assertEqual(self.smallest(None), [()])

    def test_any_yields_each_alternative(self):
        self.assertEqual(self.smallest('any(feature = "featurexml", feature = "paramxml")'), [("featurexml",), ("paramxml",)])

    def test_a_feature_implied_by_another_is_not_needed_twice(self):
        self.assertEqual(self.smallest('all(feature = "mzml", feature = "numpress")'), [("mzml",)])

    def test_the_three_feature_gate_of_the_picked_finder(self):
        gate = 'all(feature = "mzml", feature = "paramxml", feature = "featurexml")'
        self.assertEqual(self.smallest(gate), [("featurexml", "mzml", "paramxml")])


class CrateGate(unittest.TestCase):
    def gates_of(self, source):
        with tempfile.TemporaryDirectory() as root:
            (pathlib.Path(root) / "tests").mkdir()
            (pathlib.Path(root) / "tests" / "t.rs").write_text(source)
            return TestGates("t", root=pathlib.Path(root))

    def test_the_gate_at_the_top_is_the_crate_gate(self):
        gates = self.gates_of(
            '//! doc\n#![allow(x)]\n#![cfg(feature = "mzml")]\nuse a::b;\n#[cfg(feature = "idxml")]\nfn f() {}\n'
        )
        self.assertEqual((gates.gate, gates.item_features), (("feature", "mzml"), {"idxml"}))

    def test_an_inner_cfg_below_an_item_is_refused(self):
        with self.assertRaises(Unmodelled):
            self.gates_of('use a::b;\nmod m {\n#![cfg(feature = "mzml")]\n}\n')

    def test_two_crate_gates_are_refused(self):
        with self.assertRaises(Unmodelled):
            self.gates_of('#![cfg(feature = "a")]\n#![cfg(feature = "b")]\n')

    def test_a_module_pulled_in_by_path_counts(self):
        with tempfile.TemporaryDirectory() as root:
            tests = pathlib.Path(root) / "tests"
            (tests / "support").mkdir(parents=True)
            (tests / "support" / "s.rs").write_text('#[cfg(feature = "network")]\nfn g() {}\n')
            (tests / "t.rs").write_text('#[path = "support/s.rs"]\nmod s;\n')
            self.assertEqual(TestGates("t", root=pathlib.Path(root)).item_features, {"network"})


class Repository(unittest.TestCase):
    def test_the_real_gates_parse(self):
        from ci_coverage import test_targets
        for name in test_targets():
            TestGates(name)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""The generated jobs cover what the generator promises, read back as the sweep reads them.

`ci_matrix.py --check` only says the committed text is what the generator
prints. These read that text back through `ci_coverage.py`, which parses the
workflow the way the pre-push sweep does, and check the promise itself: every
integration test runs in the minimum-Rust job under each smallest feature set
its gate admits, every line runs in exactly one leg, and the spec lists no test
the derivation would have added anyway.
"""

import collections
import unittest

import ci_coverage
import ci_matrix
from ci_matrix import Bare, Job, Raw, Slice, SpecError, resolve


class Committed(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tests = ci_coverage.test_targets()
        cls.table = ci_coverage.cargo_features()
        cls.spec = ci_coverage.load_yaml(ci_coverage.WORKFLOW)
        cls.lines = ci_coverage.cargo_lines(cls.spec)
        cls.units = ci_coverage.units(cls.lines, cls.tests)

    def test_the_committed_jobs_are_what_the_generator_prints(self):
        text = ci_coverage.WORKFLOW.read_text()
        self.assertEqual(ci_matrix.split(text)[1], ci_matrix.render())

    def test_every_test_runs_in_minimum_rust_under_each_smallest_admitted_set(self):
        reduced = collections.defaultdict(set)
        for job, when, runner, toolchain, env, command, features, target in self.units:
            if job == "minimum-rust" and toolchain == "1.85.0" and command == "test --locked":
                reduced[target].add(features)
        missing = []
        for name in self.tests:
            for smallest in ci_coverage.TestGates(name).minimal_feature_sets(self.table):
                key = " ".join(sorted(smallest)) or "{}"
                if key not in reduced[f"test:{name}"]:
                    missing.append((name, key))
        self.assertEqual(missing, [])

    def test_every_test_runs_non_empty_under_some_reduced_set_of_minimum_rust(self):
        _, rows = ci_coverage.audit(self.lines, self.units, self.tests, self.table)
        self.assertEqual([n for n, r in rows.items() if r["minimum-rust"] != "runs"], [])

    def test_each_generated_line_runs_in_exactly_one_leg(self):
        for job in ("test", "minimum-rust"):
            spec_job = self.spec["jobs"][job]
            legs = spec_job["strategy"]["matrix"]["leg"]
            for step in spec_job["steps"]:
                if "run" in step:
                    runs_in = [leg for leg in legs if ci_coverage.step_runs_in(job, step, {"leg": leg})]
                    self.assertEqual(len(runs_in), 1, step)

    def test_each_leg_keeps_its_own_cache(self):
        # Legs of one job share GITHUB_JOB; without a per-leg key they would race
        # to save one cache and the loser's dependencies would never be kept.
        for job in ("test", "minimum-rust"):
            cache = [s for s in self.spec["jobs"][job]["steps"] if str(s.get("uses", "")).startswith("Swatinem/rust-cache")]
            self.assertEqual(cache[0]["with"]["key"], "${{ matrix.leg }}")


TESTS = ["a", "b"]


class Spec(unittest.TestCase):
    def resolve(self, lines, derive=False):
        return resolve(Job("j", "1.85.0", ("x", "y"), lines, derive=derive, derived_leg="x"), {}, TESTS)

    def test_an_unknown_leg_is_refused(self):
        with self.assertRaises(SpecError):
            self.resolve([Raw("z", "--all-features --all-targets")])

    def test_a_leg_with_nothing_to_run_is_refused(self):
        with self.assertRaises(SpecError):
            self.resolve([Raw("x", "--all-features --all-targets")])

    def test_two_slices_with_the_same_features_are_refused(self):
        with self.assertRaises(SpecError):
            self.resolve([Slice("x", "m p", "a"), Slice("y", "p m", "b")])

    def test_a_slice_a_bare_line_already_runs_is_refused(self):
        with self.assertRaises(SpecError):
            self.resolve([Bare("x", "m"), Slice("y", "m", "a")])

    def test_an_unknown_test_is_refused(self):
        with self.assertRaises(SpecError):
            self.resolve([Slice("x", "m", "nope"), Bare("y")])

    def test_a_bare_line_spells_its_features(self):
        out = self.resolve([Bare("x"), Bare("y", "m p")])
        self.assertEqual(out, [("x", "cargo test --locked --no-default-features"),
                               ("y", 'cargo test --locked --no-default-features --features "m p"')])


class Derivation(unittest.TestCase):
    def test_the_repository_spec_resolves(self):
        table, tests = ci_coverage.cargo_features(), ci_coverage.test_targets()
        for job in ci_matrix.JOBS:
            resolve(job, table, tests)

    def test_a_listed_test_the_derivation_adds_anyway_is_refused(self):
        table, tests = ci_coverage.cargo_features(), ci_coverage.test_targets()
        gated = next(n for n in tests if ci_coverage.TestGates(n).gate == ("feature", "mzml"))
        job = Job("j", "1.85.0", ("x",), [Bare("x"), Slice("x", "mzml", gated)], derive=True, derived_leg="x")
        with self.assertRaises(SpecError):
            resolve(job, table, tests)


if __name__ == "__main__":
    unittest.main()

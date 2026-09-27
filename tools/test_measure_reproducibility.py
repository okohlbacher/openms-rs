#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Regression checks for how measure_reproducibility.py reads what a digest names.

The numbers in docs/REPRODUCIBILITY.md are only as good as the binding of each
digest to the artifact it names. These pin the decisions that moved them: a
prefixed key names its prefix's artifact and not the entry's own path, a log's
digest is a log's, the digest of zero bytes proves nothing by content, and prose
that repeats a manifest's digest is bound as the manifest binds it.
"""

import pathlib
import unittest

import measure_reproducibility as measure

MANIFEST = pathlib.Path("tests/data/example_provenance.json")
A, B = "a" * 64, "b" * 64


def bound(tree):
    return {(digest, path, key) for digest, path, key, _ in measure.json_occurrences(MANIFEST, tree)}


class BindingTests(unittest.TestCase):
    def test_a_plain_digest_names_the_entry_path(self):
        self.assertEqual(bound({"files": [{"path": "tests/data/x.tsv", "sha256": A}]}),
                         {(A, "tests/data/x.tsv", "sha256")})

    def test_a_prefixed_digest_names_its_prefix_and_not_the_entry_path(self):
        # A fixture derived from an upstream file records that file's digest as
        # ``source_sha256`` beside its own ``path``; it is not the fixture's.
        entry = {"path": "tests/data/x.tsv", "source": "src/tests/data/x.mzML", "source_sha256": A}
        self.assertEqual(bound({"files": [entry]}), {(A, "src/tests/data/x.mzML", "source_sha256")})
        entry = {"path": "tests/data/x.tsv", "source_sha256": A}
        self.assertEqual(bound({"files": [entry]}), {(A, None, "source_sha256")})

    def test_a_hex_string_under_a_key_that_is_no_digest_is_not_counted(self):
        self.assertEqual(bound({"term": A}), set())

    def test_revisions_are_read_from_any_revision_or_commit_key(self):
        tree = {"old_package_revision": "1" * 40, "nested": {"commit": "2" * 40}, "revisionist": "3" * 40}
        self.assertEqual(measure.manifest_revisions(tree), ["1" * 40, "2" * 40])


class FamilyTests(unittest.TestCase):
    def resolve(self, path, key):
        return measure.resolve(measure.ROOT / MANIFEST, path, key, A, [], set(), None, [])

    def test_a_log_digest_with_no_path_is_a_log_never_kept(self):
        self.assertEqual(self.resolve(None, "log_sha256"), ("log", "not retained"))
        self.assertEqual(self.resolve(None, "stderr_sha256"), ("log", "not retained"))

    def test_a_crate_checksum_is_cargos(self):
        self.assertEqual(self.resolve(None, "crate"), ("crate", "cargo"))

    def test_an_install_path_is_a_build(self):
        self.assertEqual(self.resolve("/ceph/x/lib/libOpenMS.so", "sha256"), ("build", "not recomputable"))


class ContentTests(unittest.TestCase):
    def test_the_digest_of_zero_bytes_is_never_found_by_content(self):
        content = object.__new__(measure.Content)
        content.tree = {measure.EMPTY: {"tests/data/empty.txt"}, A: {"tests/data/a.txt"}}
        content.history, content.pinned = {}, {B: {"core@bc9cc12:x"}}
        self.assertIsNone(content.where(measure.EMPTY))
        self.assertEqual(content.where(A), "tree")
        self.assertEqual(content.where(B), "pinned")


class GateTests(unittest.TestCase):
    def test_only_gates_that_can_hash_are_traced(self):
        gates = [script for script, _ in measure.hashing_gates()]
        self.assertIn("tools/check_core_sdk.py", gates)
        # It shells out to cargo and hashes nothing; running it would build the crate.
        self.assertNotIn("tools/check_schema_feature_graph.py", gates)


if __name__ == "__main__":
    unittest.main()

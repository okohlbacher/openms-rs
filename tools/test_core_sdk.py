#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Regression checks for preserving historical source provenance on SDK updates."""

from copy import deepcopy
import unittest

from check_core_sdk import reference_sources


class SourceTargetReviewTests(unittest.TestCase):
    def setUp(self):
        self.revision = "2" * 40
        self.data = {
            "source_revision": "1" * 40,
            "sources": [{"path": "source.cpp", "sha256": "a" * 64}],
            "target_verification": {
                "revision": self.revision,
                "source_revisions": ["1" * 40],
                "source_changes": [{
                    "path": "source.cpp",
                    "source_sha256": "a" * 64,
                    "target_sha256": "b" * 64,
                    "review": "Reviewed change to unconsumed diagnostics.",
                }],
            },
        }

    def test_review_preserves_origin_and_selects_target_bytes(self):
        original = deepcopy(self.data)
        selected = reference_sources(self.data, self.revision, "fixture.json")
        self.assertEqual(selected[0]["sha256"], "b" * 64)
        self.assertEqual(self.data, original)
        self.data["target_verification"]["source_changes"] = []
        self.assertEqual(reference_sources(self.data, self.revision, "fixture.json"), self.data["sources"])

    def test_rejects_stale_target_or_unreviewed_origin(self):
        invalid = [deepcopy(self.data) for _ in range(4)]
        del invalid[0]["target_verification"]
        invalid[1]["target_verification"]["revision"] = "3" * 40
        invalid[2]["source_revision"] = "4" * 40
        invalid[3]["target_verification"]["source_revisions"] = []
        for data in invalid:
            with self.subTest(data=data), self.assertRaises(AssertionError):
                reference_sources(data, self.revision, "fixture.json")

    def test_rejects_wrong_source_hash_and_untraceable_changes(self):
        invalid = [deepcopy(self.data) for _ in range(5)]
        invalid[0]["target_verification"]["source_changes"][0]["source_sha256"] = "c" * 64
        invalid[1]["target_verification"]["source_changes"][0]["path"] = "unrelated.cpp"
        invalid[2]["target_verification"]["source_changes"][0]["review"] = " "
        invalid[3]["target_verification"]["source_changes"][0]["target_sha256"] = "invalid"
        invalid[4]["target_verification"]["source_changes"] *= 2
        for data in invalid:
            with self.subTest(data=data), self.assertRaises(AssertionError):
                reference_sources(data, self.revision, "fixture.json")


if __name__ == "__main__":
    unittest.main()

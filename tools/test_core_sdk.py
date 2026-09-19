#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Regression checks for preserving historical source provenance on SDK updates."""

from copy import deepcopy
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

import check_core_sdk
from check_core_sdk import check_external_artifacts, reference_sources


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


class ExternalArtifactTests(unittest.TestCase):
    """Where a C++ probe artifact moved out of this repository may not be found.

    These pin the *scope* of the search, because the scope was what was wrong:
    it walked the working tree, and the main worktree carries the gitignored
    ``.reference/`` core checkouts while a fresh one does not, so the same
    check failed from one and passed from the other. The repository is what git
    tracks, so each case here builds a small repository and asks it.
    """

    ARTIFACT = {
        "external_reference_note": "The probe sources live under ../oracle/.",
        "external_reference_artifacts": [{
            "path": "../oracle/datetime/include/OpenMS/CONCEPT/Exception.h",
            "sha256": "a" * 64,
            "origin_key": "datetime",
        }],
    }

    def setUp(self):
        self.tree = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tree, True)
        subprocess.run(["git", "-C", str(self.tree), "init", "-q"], check=True)
        (self.tree / ".gitignore").write_text(".reference/\n")
        self.add("src/lib.rs", "pub const CORE_SDK_REVISION: &str = \"\";\n")

    def add(self, path, text="x\n", tracked=True):
        item = self.tree / path
        item.parent.mkdir(parents=True, exist_ok=True)
        item.write_text(text)
        if tracked:
            subprocess.run(["git", "-C", str(self.tree), "add", "-f", path], check=True)
        return item

    def check(self, data=None):
        check_core_sdk.tracked_paths.cache_clear()
        self.addCleanup(check_core_sdk.tracked_paths.cache_clear)
        with mock.patch.object(check_core_sdk, "ROOT", self.tree):
            check_external_artifacts(deepcopy(data or self.ARTIFACT), "fixture.json")

    def test_a_pinned_checkout_beside_the_repository_is_not_a_copy_in_it(self):
        # The regression: on disk, ignored by git, and present in one worktree
        # and not another - so the verdict may not turn on it.
        self.add(".reference/openms4-core/src/openms/include/OpenMS/CONCEPT/Exception.h",
                 tracked=False)
        self.check()

    def test_a_committed_file_of_the_same_name_is_still_refused(self):
        self.add("src/openms/include/OpenMS/CONCEPT/Exception.h")
        with self.assertRaises(AssertionError):
            self.check()

    def test_a_copy_at_the_mirrored_path_is_refused_committed_or_not(self):
        # Not a C++ suffix, so only the mirrored-path check can catch it; that
        # one place is asked of the repository and of the disk both, because a
        # file sitting exactly there is unambiguous however it got there.
        data = deepcopy(self.ARTIFACT)
        data["external_reference_artifacts"][0]["path"] = "../oracle/datetime/report.json"
        self.add("docs/datetime/report.json")
        self.check(data)
        self.add("oracle/datetime/report.json", tracked=False)
        with self.assertRaises(AssertionError):
            self.check(data)

    def test_an_artifact_must_be_under_oracle_and_carry_its_origin_and_note(self):
        outside = deepcopy(self.ARTIFACT)
        outside["external_reference_artifacts"][0]["path"] = "../elsewhere/Exception.h"
        unsourced = deepcopy(self.ARTIFACT)
        del unsourced["external_reference_artifacts"][0]["origin_key"]
        unnoted = deepcopy(self.ARTIFACT)
        del unnoted["external_reference_note"]
        for data in (outside, unsourced, unnoted):
            with self.subTest(data=data), self.assertRaises(AssertionError):
                self.check(data)

    def test_tracked_paths_lists_the_repository_and_not_what_sits_beside_it(self):
        self.add(".reference/openms4-core/src/openms/CONCEPT/Exception.h", tracked=False)
        check_core_sdk.tracked_paths.cache_clear()
        self.addCleanup(check_core_sdk.tracked_paths.cache_clear)
        with mock.patch.object(check_core_sdk, "ROOT", self.tree):
            listed = check_core_sdk.tracked_paths()
        self.assertIn(Path("src/lib.rs"), listed)
        self.assertFalse([path for path in listed if path.parts[0] == ".reference"])


if __name__ == "__main__":
    unittest.main()

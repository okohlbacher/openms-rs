#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Regression checks for which module references the cycle gate can see.

The gate is a ratchet, so an edge it cannot see is a cycle it cannot refuse.
It once extracted edges with a bare ``crate::(\\w+)``, which matches nothing
after ``use crate::{`` - and most of this crate's ``crate::`` references are
braced groups. One real edge hid there long enough to close a cycle.

These pin what a reference is: the first path segment after ``crate::``, and
the first segment of each item of a braced group at any nesting depth. Nothing
deeper counts, because nothing deeper can name a top-level module.
"""

import unittest

from check_module_cycles import brace_group_end, group_items, named_after_crate


def named(text):
    return sorted(named_after_crate(text))


class FirstSegment(unittest.TestCase):
    def test_a_plain_path_names_its_first_segment_only(self):
        self.assertEqual(named("crate::metadata::MetaValue::new()"), ["metadata"])

    def test_several_plain_paths_are_all_named(self):
        self.assertEqual(
            named("use crate::kernel::X;\nuse crate::format::Y;"), ["format", "kernel"]
        )

    def test_a_longer_identifier_ending_in_crate_is_not_a_reference(self):
        self.assertEqual(named("mycrate::metadata::X"), [])

    def test_a_trailing_crate_with_nothing_after_it_is_ignored(self):
        self.assertEqual(named("something crate::"), [])


class BracedGroups(unittest.TestCase):
    def test_a_braced_group_names_the_module_a_bare_pattern_could_not_see(self):
        # The defect: nothing follows `crate::` but `{`, so `crate::(\w+)`
        # matched nothing and `param -> metadata` stayed invisible.
        self.assertEqual(
            named("use crate::{Result, metadata::{MetaInfo, MetaValue}};"),
            ["Result", "metadata"],
        )

    def test_a_module_imported_directly_in_a_group_is_named(self):
        # No `::` follows the module at all, so matching `(\w+)::` inside the
        # group would still miss both of these.
        self.assertEqual(
            named("use crate::{Error, metadata, param};"), ["Error", "metadata", "param"]
        )

    def test_each_item_of_a_group_starts_its_own_path(self):
        self.assertEqual(
            named("use crate::{format::{a, b}, kernel::c, math};"),
            ["format", "kernel", "math"],
        )

    def test_a_segment_below_the_first_is_not_named(self):
        # `metadata` here is a submodule of `identification`, not the
        # top-level module of the same name, and must not be recorded as one.
        self.assertEqual(named("use crate::{identification::metadata::X};"), ["identification"])

    def test_a_group_may_span_lines(self):
        self.assertEqual(
            named("use crate::{\n    Result,\n    metadata::MetaInfo,\n    system,\n};"),
            ["Result", "metadata", "system"],
        )

    def test_a_nested_group_does_not_end_the_outer_one(self):
        self.assertEqual(
            named("use crate::{a::{b, c}, kernel::D};"), ["a", "kernel"]
        )

    def test_an_unterminated_group_is_read_to_the_end_rather_than_hanging(self):
        self.assertEqual(named("use crate::{kernel::A, format"), ["format", "kernel"])

    def test_an_empty_group_names_nothing(self):
        self.assertEqual(named("use crate::{};"), [])


class Helpers(unittest.TestCase):
    def test_brace_group_end_finds_the_matching_close(self):
        text = "use crate::{a::{b}, c};"
        self.assertEqual(text[brace_group_end(text, text.index("{"))], "}")
        self.assertEqual(brace_group_end(text, text.index("{")), text.rindex("}"))

    def test_brace_group_end_returns_the_length_when_unbalanced(self):
        text = "use crate::{a"
        self.assertEqual(brace_group_end(text, text.index("{")), len(text))

    def test_group_items_splits_only_on_commas_the_group_itself_owns(self):
        self.assertEqual(group_items("a::{b, c}, d"), ["a::{b, c}", " d"])


if __name__ == "__main__":
    unittest.main()

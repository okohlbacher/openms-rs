#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Regression checks for the source-citation checker's reading of a document.

The checker is only worth running if it is quiet, and it is quiet because of a
handful of decisions about what a citation is, which quotation it answers for
and what may be reported on an inferred file. These pin those decisions on
synthetic text and a stand-in for the pinned source, so that they hold on a
machine with no pinned C++ checkout to hand - which is every CI runner.
"""

import pathlib
import unittest

from check_source_citations import annotations_in, attach, citations_in, problems_with, quotations, units

# A stand-in for one pinned file, with a blank line at 4 and a walk at 5-8.
LINES = [
    "void Decoder::walk()",                    # 1
    "{",                                       # 2
    "  DOMNode* iter = getFirstChild();",      # 3
    "",                                        # 4
    "  while (iter != lastChild)",             # 5
    "  {",                                     # 6
    "    iter = iter->getNextSibling();",      # 7
    "  }",                                     # 8
    "}",                                       # 9
]


class Pinned:
    """Stands in for a pinned checkout: one revision, one file."""

    @staticmethod
    def lines(_revision, _path):
        return LINES


def check(ranges, quoted=(), annotated=()):
    return problems_with(Pinned(), "r", "Decoder.cpp", ranges, quoted, annotated)


class CitationReadingTests(unittest.TestCase):
    def test_reads_a_path_a_single_line_and_an_en_dash_range(self):
        named, bare = citations_in("see HANDLERS/Decoder.cpp:3 and Other.h:10–12", True)
        self.assertEqual(
            [(item[1], item[2], item[3], item[4]) for item in named],
            [("HANDLERS/", "Decoder.cpp", 3, 3), ("", "Other.h", 10, 12)],
        )
        self.assertEqual(bare, [])

    def test_a_bare_range_continues_a_file_only_where_it_is_quoted(self):
        quoted = "`Decoder.cpp:3`, `:5-8`"
        self.assertEqual([item[1:] for item in citations_in(quoted, True)[1]], [(5, 8, ":5-8")])
        # A table cell writing "startProgress (:288)" of a class it names without
        # an extension must not be read as a range of the file cited beside it.
        plain = "`Decoder.cpp:3`, ProgressLogger (:288)"
        self.assertEqual(citations_in(plain, True)[1], [])
        self.assertEqual([item[1:] for item in citations_in(plain, False)[1]], [(288, 288, ":288")])

    def test_a_version_or_a_time_is_not_a_bare_range(self):
        self.assertEqual(citations_in("`Decoder.cpp:3`, version 1:2, gcc-14:4", True)[1], [])

    def test_a_bare_range_needs_a_file_before_it(self):
        self.assertEqual(citations_in("`:5-8` alone", True), ([], []))


class QuotationTests(unittest.TestCase):
    def fragments(self, unit):
        return {item[2] for item in quotations(unit)}

    def test_a_name_is_a_reference_and_a_statement_is_a_quotation(self):
        self.assertEqual(
            self.fragments(
                "`writeHeader_`, `MzTabFile::load` and `updateRanges()` reach "
                "`iter = iter->getNextSibling();` here"
            ),
            {"iter = iter->getNextSibling();"},
        )

    def test_a_citation_is_not_part_of_the_quotation_beside_it(self):
        self.assertEqual(self.fragments("`Decoder.cpp:3-7`"), set())

    def test_a_quotation_is_flattened_across_a_line_break(self):
        self.assertEqual(
            self.fragments("`while (iter !=\n   lastChild)`"), {"while (iter != lastChild)"}
        )


class AttachmentTests(unittest.TestCase):
    def attached(self, unit):
        named, _ = citations_in(unit, True)
        placed = [(item[0], len(item[5]), item[2]) for item in named]
        return {key: sorted(value) for key, value in attach(unit, placed, quotations(unit)).items()}

    def test_the_nearest_quotation_is_the_one_the_citation_answers_for(self):
        unit = "reuses `iter = getFirstChild();` and writes `iter = iter->getNextSibling();` (`Decoder.cpp:7`)"
        self.assertEqual(self.attached(unit), {"Decoder.cpp": ["iter = iter->getNextSibling();"]})

    def test_a_full_stop_or_a_semicolon_detaches_a_quotation(self):
        for glue in (". Elsewhere,", "; elsewhere,"):
            unit = f"`Decoder.cpp:3`{glue} the other branch writes `iter = iter->getNextSibling();`"
            with self.subTest(glue=glue):
                self.assertEqual(self.attached(unit), {})

    def test_a_quotation_after_the_citation_still_attaches(self):
        unit = "`Decoder.cpp:5-8` runs `iter = iter->getNextSibling();` once"
        self.assertEqual(self.attached(unit), {"Decoder.cpp": ["iter = iter->getNextSibling();"]})


class UnitTests(unittest.TestCase):
    def spans(self, name, text):
        return [span for span, _ in units(pathlib.Path(name), text, ("r",))]

    def test_a_bullet_a_row_and_a_paragraph_are_separate_units(self):
        text = "- first `A.cpp:1`\n- second `B.cpp:2`\n\n| a | `C.cpp:3` |\nprose `D.cpp:4`\n"
        self.assertEqual(
            self.spans("doc.md", text),
            ["- first `A.cpp:1`", "- second `B.cpp:2`", "| a | `C.cpp:3` |", "prose `D.cpp:4`"],
        )

    def test_a_fenced_block_is_not_a_unit(self):
        self.assertEqual(self.spans("doc.md", "before\n\n```\n`A.cpp:1`\n```\n\nafter"), ["before", "after"])

    def test_a_manifest_is_cut_into_its_string_values(self):
        self.assertEqual(
            self.spans("m.json", '{"why": "A.cpp:1", "how": ["B.cpp:2"]}'), ["A.cpp:1", "B.cpp:2"]
        )


class ProblemTests(unittest.TestCase):
    def test_a_range_may_not_run_backwards_or_end_past_the_file(self):
        self.assertIn("runs backwards", check([(8, 5, ":8-5", True)])[0])
        self.assertIn("has 9 lines", check([(5, 99, ":5-99", True)])[0])

    def test_a_named_single_line_may_not_be_blank_and_an_inferred_one_may(self):
        self.assertIn("blank", check([(4, 4, "Decoder.cpp:4", True)])[0])
        self.assertEqual(check([(4, 4, ":4", False)]), [])

    def test_a_quotation_must_be_inside_the_cited_lines(self):
        quoted = {"iter = iter->getNextSibling();"}
        self.assertEqual(check([(5, 8, ":5-8", True)], quoted), [])
        self.assertIn("not on the cited lines", check([(1, 4, ":1-4", True)], quoted)[0])

    def test_several_ranges_for_one_file_are_taken_together(self):
        quoted = {"iter = iter->getNextSibling();"}
        self.assertEqual(check([(1, 3, ":1-3", True), (7, 7, ":7", True)], quoted), [])

    def test_a_quotation_the_file_does_not_contain_is_not_a_finding(self):
        self.assertEqual(check([(1, 4, ":1-4", True)], {"iter = next(iter);"}), [])

    def test_an_annotation_must_match_the_line_it_names(self):
        self.assertEqual(check([(3, 3, ":3", True)], (), annotations_in("  DOMNode* iter = getFirstChild();   // :3")), [])
        wrong = annotations_in("  DOMNode* iter = getFirstChild();   // :5")
        self.assertIn("// :5 is", check([(3, 3, ":3", True)], (), wrong)[0])


if __name__ == "__main__":
    unittest.main()

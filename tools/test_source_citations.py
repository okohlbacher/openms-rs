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

import collections
import pathlib
import unittest

from check_source_citations import (
    Unreadable, annotations_in, attach, check_unit, citations_in, named_file, problems_with,
    quotations, units,
)

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

    sources = {"r": "a checkout"}

    @staticmethod
    def lines(_revision, _path):
        return LINES

    @staticmethod
    def candidates(_revision, name, _directory):
        return ["HANDLERS/" + name] if name == "Decoder.cpp" else []

    @staticmethod
    def where(_revision, path):
        return "pinned:" + path


def check(ranges, quoted=(), annotated=()):
    return problems_with(Pinned(), "r", "Decoder.cpp", ranges, quoted, annotated)[0]


def quoted_at(fragment, first, last, text=None):
    """A quotation as :func:`attach` hands it on: with the range it stands beside."""
    return (fragment, first, last, text or f":{first}-{last}")


def findings_for(unit, context=None):
    report = {"checked": 0, "quoted": 0, "skipped": 0, "unresolved": 0, "skipped_files": collections.Counter()}
    found = check_unit(Pinned(), unit, ("r",), True, report, context)
    return [problem for _, problems in found for problem in problems]


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

    def test_an_operator_name_is_a_reference_and_an_operator_call_is_not(self):
        self.assertEqual(self.fragments("`MSChromatogram::operator==` and `std::operator==`"), set())
        self.assertEqual(
            self.fragments("`f.BaseFeature::operator=(c)`"), {"f.BaseFeature::operator=(c)"}
        )

    def test_a_manifest_quotation_is_recognised_by_its_shape(self):
        from check_source_citations import BARE_QUOTATION
        found = {item[2] for item in quotations(
            "the walk is iter = iter->getNextSibling(); bounded against t_wait = 0.2 s",
            BARE_QUOTATION,
        )}
        self.assertEqual(found, {"iter = iter->getNextSibling()"})

    def test_a_citation_is_not_part_of_the_quotation_beside_it(self):
        self.assertEqual(self.fragments("`Decoder.cpp:3-7`"), set())

    def test_a_quotation_is_flattened_across_a_line_break(self):
        self.assertEqual(
            self.fragments("`while (iter !=\n   lastChild)`"), {"while (iter != lastChild)"}
        )


class AttachmentTests(unittest.TestCase):
    def attached(self, unit):
        named, _ = citations_in(unit, True)
        placed = [(item[0], len(item[5]), item[5]) for item in named]
        return attach(unit, placed, quotations(unit))

    def test_the_nearest_quotation_is_the_one_the_citation_answers_for(self):
        unit = "reuses `iter = getFirstChild();` and writes `iter = iter->getNextSibling();` (`Decoder.cpp:7`)"
        self.assertEqual(self.attached(unit), {"Decoder.cpp:7": "iter = iter->getNextSibling();"})

    def test_a_full_stop_a_semicolon_or_a_colon_detaches_a_quotation(self):
        for glue in (". Elsewhere,", "; elsewhere,", ": still inside", ".** Elsewhere,"):
            unit = f"`Decoder.cpp:3`{glue} the other branch writes `iter = iter->getNextSibling();`"
            with self.subTest(glue=glue):
                self.assertEqual(self.attached(unit), {})

    def test_a_quotation_after_the_citation_still_attaches(self):
        unit = "`Decoder.cpp:5-8` runs `iter = iter->getNextSibling();` once"
        self.assertEqual(self.attached(unit), {"Decoder.cpp:5-8": "iter = iter->getNextSibling();"})

    def test_one_quotation_among_several_citations_answers_for_only_one(self):
        # Seven of these eight would be reported if the quotation were shared out.
        unit = ("allocates at `Decoder.cpp:1-2`, walks under `iter = iter->getNextSibling();` "
                "(`Decoder.cpp:7`, body to `Decoder.cpp:8`) and returns at `Decoder.cpp:9`")
        self.assertEqual(self.attached(unit), {"Decoder.cpp:7": "iter = iter->getNextSibling();"})

    def test_two_citations_each_keep_their_own_quotation(self):
        # The second quotation stands four characters from the first citation and
        # six from its own; pairing the paragraph as a whole does not cross them.
        unit = ("`DOMNode* iter = getFirstChild();` at `Decoder.cpp:3`, "
                "`iter = iter->getNextSibling();` at `Decoder.cpp:7`")
        self.assertEqual(
            self.attached(unit),
            {
                "Decoder.cpp:3": "DOMNode* iter = getFirstChild();",
                "Decoder.cpp:7": "iter = iter->getNextSibling();",
            },
        )


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
        walk = "iter = iter->getNextSibling();"
        self.assertEqual(check([(5, 8, ":5-8", True)], [quoted_at(walk, 5, 8)]), [])
        self.assertIn(
            "not on the cited lines", check([(1, 4, ":1-4", True)], [quoted_at(walk, 1, 4)])[0]
        )

    def test_a_quotation_is_read_against_the_range_it_stands_beside(self):
        # The range it stands beside, not the union of everything the unit cites
        # for that file: a right citation may not vouch for a wrong one next to it.
        walk = "iter = iter->getNextSibling();"
        ranges = [(1, 3, ":1-3", True), (7, 7, ":7", True)]
        self.assertEqual(check(ranges, [quoted_at(walk, 7, 7, ":7")]), [])
        self.assertIn("not on the cited lines", check(ranges, [quoted_at(walk, 1, 3, ":1-3")])[0])

    def test_a_span_the_unit_cites_around_the_line_answers_for_it_too(self):
        # "the walk, `Decoder.cpp:1-9`" and then single lines inside it: a
        # quotation of the span is a quotation of the span, not of the line.
        walk = "iter = iter->getNextSibling();"
        ranges = [(1, 9, ":1-9", True), (5, 5, ":5", True)]
        self.assertEqual(check(ranges, [quoted_at(walk, 5, 5, ":5")]), [])

    def test_a_quotation_the_file_does_not_contain_is_not_a_finding(self):
        self.assertEqual(check([(1, 4, ":1-4", True)], [quoted_at("iter = next(iter);", 1, 4)]), [])

    def test_an_annotation_must_match_the_line_it_names(self):
        self.assertEqual(check([(3, 3, ":3", True)], (), annotations_in("  DOMNode* iter = getFirstChild();   // :3")), [])
        wrong = annotations_in("  DOMNode* iter = getFirstChild();   // :5")
        self.assertIn("// :5 is", check([(3, 3, ":3", True)], (), wrong)[0])

    def test_an_annotation_on_an_inferred_file_is_reported_only_if_it_fits_it(self):
        # The file came from the entry the block sits under, not from a citation
        # in the block, so code the file does not hold at all is another file's.
        here = annotations_in("  DOMNode* iter = getFirstChild();   // :5", certain=False)
        self.assertIn("// :5 is", check([], (), here)[0])
        elsewhere = annotations_in("  int unrelated_zz = 1;   // :5", certain=False)
        self.assertEqual(check([], (), elsewhere), [])


class EntryTests(unittest.TestCase):
    """The checks that need a whole entry, not one span of it."""

    ENTRY = (
        "## CPP-999 - Decoder skips the first child\n"
        "\n"
        "**Affected file/function:** `HANDLERS/Decoder.cpp`, `Decoder::walk` (`:1`)\n"
        "\n"
        "**Issue:** The child walk is\n"
        "\n"
        "    DOMNode* iter = getFirstChild();                      // :3\n"
        "    while (iter != lastChild)                             // :5\n"
    )

    def entry(self, text):
        found, context = [], None
        for unit, _ in units(pathlib.Path("log.md"), text, ("r",)):
            context = named_file(unit) or context
            found.extend(findings_for(unit, context))
        return found

    def test_the_file_named_by_the_entry_reaches_the_transcribed_block(self):
        self.assertEqual(self.entry(self.ENTRY), [])

    def test_a_wrong_annotation_in_a_transcribed_block_is_reported(self):
        # The block is its own span - a blank line separates it from the line
        # that names its file - so nothing in it cites anything.
        wrong = self.ENTRY.replace("// :5", "// :7")
        self.assertIn("// :7 is", "".join(self.entry(wrong)))

    def test_named_file_carries_the_last_file_the_span_names(self):
        self.assertEqual(
            named_file("**Affected:** `OTHER/First.cpp` and `HANDLERS/Decoder.cpp`"),
            ("HANDLERS/", "Decoder.cpp"),
        )
        self.assertIsNone(named_file("no source file here"))

    def test_a_manifest_that_does_not_parse_is_named_rather_than_crashing(self):
        with self.assertRaises(Unreadable):
            units(pathlib.Path("m.json"), '{"why": "A.cpp:1",}', ("r",))


if __name__ == "__main__":
    unittest.main()

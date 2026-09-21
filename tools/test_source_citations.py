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
import hashlib
import pathlib
import tempfile
import unittest
from unittest import mock

import check_source_citations as checker
from check_source_citations import (
    Directory, Objects, Pins, Retained, Unreadable, annotations_in, attach, check_file,
    check_unit, citations_in, holding, named_file, pair_up, problems_with, quotations,
    resolvable, split_package, tally, units,
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

    declared = {}

    @staticmethod
    def paths_named(_revision, name):
        return ["HANDLERS/" + name] if name == "Decoder.cpp" else []

    # `resolvable` asks these two after the offered revisions; this stand-in has
    # neither a package pin nor a retained directory, so both answer nothing.
    @staticmethod
    def packaged(_inside, _name, _besides):
        return []

    @staticmethod
    def retained(_inside, _name):
        return []

    @staticmethod
    def where(_revision, path):
        return "pinned:" + path

    @staticmethod
    def label(revision):
        return revision


def check(ranges, quoted=(), annotated=()):
    return problems_with(Pinned(), "r", "Decoder.cpp", ranges, quoted, annotated)[0]


def quoted_at(fragment, first, last, text=None):
    """A quotation as :func:`attach` hands it on: with the range it stands beside."""
    return (fragment, first, last, text or f":{first}-{last}")


def findings_for(unit, context=None):
    report = tally()
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

    def test_a_name_with_a_plus_in_it_is_read_whole(self):
        # `bits/c++locale.h` is a real libstdc++ header. Matching from the
        # second `+` reported it as `locale.h` - a file the citation does not
        # name and no pin has - under "cited files that no reachable pin
        # contains", which is the one place the report must not be wrong.
        named, _ = citations_in("stripped at `bits/c++locale.h:74`", True)
        self.assertEqual(
            [(item[1], item[2], item[3]) for item in named], [("bits/", "c++locale.h", 74)]
        )
        self.assertEqual(named_file("see `bits/c++locale.h`"), ("bits/", "c++locale.h"))

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
        for glue in (". Elsewhere,", "; elsewhere,", ": still inside", ".** Elsewhere,",
                     "\u2014 Elsewhere,"):
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

    def test_the_greedy_fallback_may_pair_differently(self):
        # Above MOST_QUOTATIONS_PAIRED the pairing is taken greedily, and greedy
        # is not merely a smaller version of the whole-unit answer: it takes the
        # closest pair, which here costs it both of the others. Nothing in this
        # repository has that many quotations in one unit; this pins what the
        # branch does if one ever grows to.
        allowed = {(0, 0): 1, (0, 1): 10, (1, 0): 2, (1, 1): 100}
        self.assertEqual(pair_up(allowed, 2, 2), ((0, 1), (1, 0)))
        with mock.patch.object(checker, "MOST_QUOTATIONS_PAIRED", 1):
            self.assertEqual(pair_up(allowed, 2, 2), [(0, 0), (1, 1)])

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


CORE, TOPP = "c" * 40, "t" * 40


class TwoPins:
    """Two pins that both carry a ``FileInfo.cpp``: the core SDK, and a TOPP tool.

    Fourteen file names collide across the pins this repository declares, and
    this is the shape of all of them - one file deep in the SDK's tree, one at
    the top of a tool's ``src/``, sharing nothing but a name.
    """

    declared = {"core": CORE, "topp": TOPP}
    sources = {CORE: "core bc9cc12", TOPP: "topp 174b576"}
    paths = {CORE: ["src/openms/source/FORMAT/FileInfo.cpp"], TOPP: ["src/FileInfo.cpp"]}
    text = {
        "src/openms/source/FORMAT/FileInfo.cpp": [
            "void FileInfo::run()",                      # 1
            "{",                                         # 2
            '  os << "Number of peaks: " << count;',     # 3
            "",                                          # 4
            "  for (const auto& s : exp)",               # 5
            "  {",                                       # 6
            "    total += s.size();",                    # 7
            "  }",                                       # 8
            "}",                                         # 9
        ],
        "src/FileInfo.cpp": ['  registerFlag_("c", "Check for corrupt data");'],
    }

    def paths_named(self, revision, name):
        return self.paths[revision] if name == "FileInfo.cpp" else []

    def packaged(self, inside, name, besides):
        # Only the tool is a package pin here; the SDK checkout is not.
        return [] if TOPP in besides else [
            (TOPP, path) for path in self.paths_named(TOPP, name)
            if path.endswith(inside + name)
        ]

    def retained(self, _inside, _name):
        return []  # Nothing retained here; `RetainedSourceTests` covers that limb.

    def lines(self, _revision, path):
        return self.text[path]

    def label(self, revision):
        return self.sources[revision]

    def where(self, revision, path):
        return f"{self.sources[revision]}:{path}"


class OnePinTwoPaths:
    """One pin that carries ``Macros.h`` at two paths, as the core checkout does.

    A file name does not need two *pins* to be ambiguous. At ``bc9cc12`` the
    core SDK carries ``Macros.h`` under both ``CONCEPT/`` and
    ``OPENSWATHALGO/``, so a bare ``Macros.h:1`` reaches two different files
    inside a single revision and nothing but the quotation beside it can say
    which one the document meant.
    """

    declared = {"core": CORE}
    sources = {CORE: "core bc9cc12"}
    paths = {
        CORE: [
            "src/openms/include/OpenMS/CONCEPT/Macros.h",
            "src/openswathalgo/include/OpenMS/OPENSWATHALGO/Macros.h",
        ]
    }
    text = {
        "src/openms/include/OpenMS/CONCEPT/Macros.h": ["#define OPENMS_PRECONDITION(x, m)"],
        "src/openswathalgo/include/OpenMS/OPENSWATHALGO/Macros.h": ["#define OPENSWATH_PRECONDITION(x, m)"],
    }

    def paths_named(self, revision, name):
        return self.paths[revision] if name == "Macros.h" else []

    def packaged(self, inside, name, besides):
        return []

    def retained(self, _inside, _name):
        return []

    def lines(self, _revision, path):
        return self.text[path]

    def label(self, _revision):
        return "core bc9cc12"

    def where(self, revision, path):
        return f"core bc9cc12:{path}"


FLAG = 'registerFlag_("c", "Check for corrupt data");'
PEAKS = 'os << "Number of peaks: " << count;'


class PinResolutionTests(unittest.TestCase):
    """Which pin answers a citation, when two of them carry the file name."""

    def resolved(self, directory):
        return [revision for revision, _ in resolvable(TwoPins(), (CORE, TOPP), directory, "FileInfo.cpp")]

    def test_a_bare_name_both_pins_carry_resolves_into_both(self):
        self.assertEqual(self.resolved(""), [CORE, TOPP])

    def test_a_package_prefix_resolves_in_that_pin_and_nowhere_else(self):
        self.assertEqual(self.resolved("OpenMS4-topp/src/"), [TOPP])
        self.assertEqual(self.resolved("OpenMS4-tests/packages/topp/src/"), [TOPP])

    def test_a_prefix_naming_a_package_with_no_pin_resolves_nowhere(self):
        self.assertEqual(self.resolved("OpenMS4-flash/src/"), [])

    def test_a_path_that_fits_one_pin_drops_the_pins_it_does_not_fit(self):
        # The defect this closes: "src/" is the tool's whole layout and fits
        # nothing in the SDK, so the SDK may not answer with its own file of
        # that name - it did, because it is the revision declared first.
        self.assertEqual(self.resolved("src/"), [TOPP])
        self.assertEqual(self.resolved("src/openms/source/FORMAT/"), [CORE])

    def test_a_package_pin_claims_a_path_no_offered_revision_fits(self):
        # An issue-log entry is checked against the revisions it names, and
        # those are core revisions; a TOPP path in such an entry used to find
        # nothing to narrow in and fall back to the SDK file of that name.
        self.assertEqual(
            [revision for revision, _ in resolvable(TwoPins(), (CORE,), "src/", "FileInfo.cpp")],
            [TOPP],
        )
        # But an offered revision that does fit keeps it.
        self.assertEqual(
            [revision for revision, _ in
             resolvable(TwoPins(), (CORE,), "src/openms/source/FORMAT/", "FileInfo.cpp")],
            [CORE],
        )

    def test_a_path_that_fits_no_pin_at_all_still_resolves_by_name(self):
        # Some citations write an absolute path into a .reference/ checkout.
        self.assertEqual(
            self.resolved(".reference/openms4-core-bc9cc12/src/openms/source/FORMAT/"),
            [CORE, TOPP],
        )

    def test_split_package_leaves_an_ordinary_path_alone(self):
        self.assertEqual(split_package("FORMAT/"), (None, "FORMAT/"))
        self.assertEqual(split_package("OpenMS4-topp/src/"), ("topp", "src/"))
        self.assertEqual(
            split_package("OpenMS4-tests/packages/test-data/topp/"), ("test_data", "topp/")
        )

    def test_a_tool_citation_is_read_against_the_tool_and_not_the_sdk(self):
        # ":1-9" exists in the SDK's file of the same name and not in the
        # tool's, so before this the SDK answered and the citation passed.
        report = tally()
        found = check_file(
            TwoPins(), (CORE, TOPP), ("OpenMS4-topp/src/", "FileInfo.cpp"),
            [(1, 9, "OpenMS4-topp/src/FileInfo.cpp:1-9", True)], (), (), report,
        )
        self.assertEqual([where for where, _ in found], ["topp 174b576:src/FileInfo.cpp"])
        self.assertIn("the file has 1 lines", found[0][1][0])
        self.assertEqual(report["checked"], 0)

    def test_the_pin_whose_file_holds_the_quotation_answers_a_bare_name(self):
        report = tally()
        found = check_file(
            TwoPins(), (CORE, TOPP), ("", "FileInfo.cpp"), [(1, 1, "FileInfo.cpp:1", True)],
            [(FLAG, 1, 1, "FileInfo.cpp:1")], (), report,
        )
        self.assertEqual(found, [])
        self.assertEqual(report["answered"], collections.Counter({"topp 174b576": 1}))
        self.assertEqual(report["confirmed_by"], collections.Counter({"topp 174b576": 1}))
        self.assertEqual(report["ambiguous"], 0)

    def test_choosing_the_pin_by_its_quotation_does_not_excuse_the_line(self):
        # The SDK holds this line, so the SDK answers - and then says it is at
        # :3 and not at the :5 the citation named.
        report = tally()
        found = check_file(
            TwoPins(), (CORE, TOPP), ("", "FileInfo.cpp"), [(5, 5, "FileInfo.cpp:5", True)],
            [(PEAKS, 5, 5, "FileInfo.cpp:5")], (), report,
        )
        self.assertEqual([where for where, _ in found], ["core bc9cc12:src/openms/source/FORMAT/FileInfo.cpp"])
        self.assertIn("not on the cited lines", found[0][1][0])

    def test_a_bare_name_with_nothing_quoted_beside_it_is_counted_ambiguous(self):
        report = tally()
        self.assertEqual(
            check_file(TwoPins(), (CORE, TOPP), ("", "FileInfo.cpp"),
                       [(1, 1, "FileInfo.cpp:1", True)], (), (), report),
            [],
        )
        self.assertEqual(report["ambiguous"], 1)
        self.assertEqual(report["ambiguous_files"]["FileInfo.cpp"], 1)
        self.assertEqual(report["answered"], collections.Counter({"core bc9cc12": 1}))
        # An ambiguous citation is still checked - one of the pins answered it -
        # so the count is a subset of "checked" and not a fourth kind of
        # unchecked. The summary line has to be read that way, and says so.
        self.assertEqual(report["checked"], 1)
        self.assertEqual(report["skipped"], 0)

    def test_a_name_only_one_pin_carries_is_not_ambiguous(self):
        report = tally()
        check_file(TwoPins(), (CORE,), ("", "FileInfo.cpp"),
                   [(1, 1, "FileInfo.cpp:1", True)], (), (), report)
        self.assertEqual(report["ambiguous"], 0)
        self.assertEqual(report["checked"], 1)

    def test_a_name_one_pin_carries_at_two_paths_is_counted_ambiguous_too(self):
        # The guard used to ask how many REVISIONS were still standing, so a
        # name one pin carried at two paths was settled silently by whichever
        # file agreed first. Path granularity is the same risk one level down.
        report = tally()
        self.assertEqual(
            check_file(OnePinTwoPaths(), (CORE,), ("", "Macros.h"),
                       [(1, 1, "Macros.h:1", True)], (), (), report),
            [],
        )
        self.assertEqual(report["ambiguous"], 1)
        self.assertEqual(report["ambiguous_files"]["Macros.h"], 1)
        self.assertEqual(report["checked"], 1)

    def test_a_quotation_tells_two_paths_of_one_pin_apart(self):
        # And when the document does quote the file, nothing is ambiguous:
        # `holding` narrows inside the revision exactly as it does across pins.
        report = tally()
        quoted = [("#define OPENSWATH_PRECONDITION(x, m)", 1, 1, "Macros.h:1")]
        self.assertEqual(
            check_file(OnePinTwoPaths(), (CORE,), ("", "Macros.h"),
                       [(1, 1, "Macros.h:1", True)], quoted, (), report),
            [],
        )
        self.assertEqual(report["ambiguous"], 0)
        self.assertEqual(report["checked"], 1)
        self.assertEqual(report["quoted"], 1)

    def test_a_citation_no_candidate_agrees_with_is_reported_not_counted_ambiguous(self):
        # The other direction of the subset relation. When every candidate
        # disagrees the citation is REPORTED by name, and `checked` does not
        # move - so `ambiguous` must not move either, or the summary would say
        # "N checked, M of them answered by a file that may be the wrong one"
        # with M not a subset of N.
        report = tally()
        found = check_file(TwoPins(), (CORE, TOPP), ("", "FileInfo.cpp"),
                           [(99, 100, "FileInfo.cpp:99-100", True)], (), (), report)
        self.assertEqual(len(found), 1)
        self.assertIn("the file has", found[0][1][0])
        self.assertEqual(report["checked"], 0)
        self.assertEqual(report["ambiguous"], 0)
        self.assertEqual(report["ambiguous_files"], collections.Counter())
        self.assertEqual(report["skipped"], 0)
        # It was still answered by one of them, which is what `--report` shows.
        self.assertEqual(sum(report["answered"].values()), 1)

    def test_holding_decides_nothing_when_no_candidate_has_the_quotation(self):
        attempts = [(CORE, "src/openms/source/FORMAT/FileInfo.cpp"), (TOPP, "src/FileInfo.cpp")]
        self.assertEqual(holding(TwoPins(), attempts, [("nothing = here();", 1, 1, ":1")]), [])
        self.assertEqual(holding(TwoPins(), attempts, ()), [])


class PackagedPinTests(unittest.TestCase):
    """`Pins.packaged` itself, on a real `Pins` rather than a stand-in.

    Every other test here supplies its own `packaged`, so the shipped method
    was never executed by the suite: mutating its comprehension so that it can
    never match left all of it green while moving two citations on the real
    tree. These call the real one.

    It is the third limb of "resolve a citation in the pin its own path names":
    the issue log checks each entry against the revisions that entry names, and
    those are core revisions, so a TOPP path written in an old entry found
    nothing to narrow in and fell back to the core file of the same name. A
    package has exactly one pin here and no older revision to prefer, so a path
    that fits its layout names it whatever the entry says.
    """

    CORE_PATH = "src/openms/source/PROCESSING/CENTROIDING/PeakPickerHiRes.cpp"
    TOPP_PATH = "src/PeakPickerHiRes.cpp"

    def pins(self):
        """A `Pins` with one `Directory` core pin and one `Objects` package pin."""
        pins = object.__new__(Pins)
        pins.declared = {"core": CORE, "topp": TOPP}
        pins.sources = {
            CORE: Directory(pathlib.Path("/nowhere/openms4-core-bc9cc12")),
            TOPP: Objects(pathlib.Path("/nowhere/packages/topp"), TOPP),
        }
        pins.index = {
            (CORE, "PeakPickerHiRes.cpp"): [self.CORE_PATH],
            (TOPP, "PeakPickerHiRes.cpp"): [self.TOPP_PATH],
        }
        return pins

    def test_a_package_pin_claims_a_path_only_its_layout_fits(self):
        self.assertEqual(
            self.pins().packaged("src/", "PeakPickerHiRes.cpp", (CORE,)),
            [(TOPP, self.TOPP_PATH)],
        )

    def test_a_core_shaped_path_is_not_claimed_by_the_package_pin(self):
        self.assertEqual(
            self.pins().packaged(
                "src/openms/source/PROCESSING/CENTROIDING/", "PeakPickerHiRes.cpp", (CORE,)
            ),
            [],
        )

    def test_a_revision_already_offered_is_not_claimed_again(self):
        # `besides` is what `resolvable` already tried; claiming one of those
        # would answer twice with the same file.
        self.assertEqual(
            self.pins().packaged("src/", "PeakPickerHiRes.cpp", (CORE, TOPP)), []
        )

    def test_only_a_package_pin_may_claim_it_never_an_unpacked_checkout(self):
        # The core pin is a `Directory`, not an `Objects`, so even a path that
        # fits its layout is not claimed this way - `resolvable` had its chance.
        pins = self.pins()
        self.assertEqual(
            pins.packaged("src/openms/source/PROCESSING/CENTROIDING/", "PeakPickerHiRes.cpp", ()),
            [],
        )
        # And with both offered as package pins, only the one whose layout fits.
        pins.sources[CORE] = Objects(pathlib.Path("/nowhere/core"), CORE)
        self.assertEqual(
            pins.packaged("src/", "PeakPickerHiRes.cpp", ()), [(TOPP, self.TOPP_PATH)]
        )

    def test_a_name_the_package_pin_does_not_carry_is_not_claimed(self):
        self.assertEqual(self.pins().packaged("src/", "FileInfo.cpp", (CORE,)), [])


class RetainedSourceTests(unittest.TestCase):
    """`Retained` and `Pins.retained`: the digest is what admits a file.

    A ``libstdc++`` header is not in any pin and cannot be, so a citation of
    ``bits/stl_algo.h:1806`` was counted unreachable and nothing read it back.
    Reading a retained copy instead is only worth anything if the copy is the
    one the port measured, which is the whole content of these: a file whose
    bytes do not hash to the declared digest, or that is not there, must go on
    counting as unreachable rather than quietly answering a citation. The
    failure this guards against is a later toolchain leaving a *different*
    ``stl_algo.h`` at that path, whose line 1806 is something else entirely -
    which would confirm a citation against a file the port never read.
    """

    BODY = b"enum { _S_threshold = 16 };\nwhile (__last - __first > 16)\n"
    DIGEST = hashlib.sha256(BODY).hexdigest()

    def retained(self, body=BODY, digest=None, under="bits/", write=True):
        # `addCleanup`, not `enterContext`: the latter is Python 3.11 and up,
        # and nothing else here asks for a floor that high.
        holder = tempfile.TemporaryDirectory()
        self.addCleanup(holder.cleanup)
        directory = pathlib.Path(holder.name)
        if write:
            (directory / "stl_algo.h").write_bytes(body)
        return Retained(directory, {"stl_algo.h": digest or self.DIGEST}, under)

    def test_a_file_at_the_declared_digest_is_admitted_under_its_cited_path(self):
        source = self.retained()
        self.assertEqual(source.paths(), ["bits/stl_algo.h"])
        self.assertEqual(source.refused, {})
        self.assertIn("_S_threshold", source.text("bits/stl_algo.h"))

    def test_a_file_whose_bytes_changed_is_refused_and_stays_unreachable(self):
        source = self.retained(body=b"something else entirely\n")
        self.assertEqual(source.paths(), [])
        self.assertIn("stl_algo.h", source.refused)
        self.assertIn("sha256", source.refused["stl_algo.h"])

    def test_a_file_that_is_not_retained_here_is_refused(self):
        source = self.retained(write=False)
        self.assertEqual(source.paths(), [])
        self.assertEqual(source.refused, {"stl_algo.h": "not retained here"})

    def test_a_declared_digest_no_file_matches_admits_nothing(self):
        source = self.retained(digest="0" * 64)
        self.assertEqual(source.paths(), [])

    def pins(self, under="bits/", body=BODY):
        """A `Pins` with one core pin and one retained directory."""
        pins = object.__new__(Pins)
        pins.declared = {"core": CORE}
        pins.sources = {CORE: Directory(pathlib.Path("/nowhere/openms4-core-bc9cc12"))}
        pins.index = collections.defaultdict(list)
        pins.index[(CORE, "Macros.h")] = ["src/openms/include/OpenMS/CONCEPT/Macros.h"]
        source = self.retained(body=body, under=under)
        key = "retained:../oracle/libstdcxx"
        pins.sources[key] = source
        pins.retained_keys = [key]
        for path in source.paths():
            pins.index[(key, path.rsplit("/", 1)[-1])].append(path)
        pins._lines = {}
        return pins, key

    def test_a_retained_file_answers_a_name_no_pin_carries(self):
        pins, key = self.pins()
        self.assertEqual(
            resolvable(pins, (CORE,), "bits/", "stl_algo.h"), [(key, "bits/stl_algo.h")]
        )
        # And with no directory written at all, which is how half of them are cited.
        self.assertEqual(
            resolvable(pins, (CORE,), "", "stl_algo.h"), [(key, "bits/stl_algo.h")]
        )

    def test_a_refused_file_answers_nothing(self):
        pins, _ = self.pins(body=b"a different header\n")
        self.assertEqual(resolvable(pins, (CORE,), "bits/", "stl_algo.h"), [])

    def test_a_citation_writing_another_directory_is_not_answered_by_it(self):
        pins, _ = self.pins()
        self.assertEqual(resolvable(pins, (CORE,), "ext/", "stl_algo.h"), [])

    def test_a_retained_copy_never_takes_a_name_a_pin_carries(self):
        # The limb is last for this reason: retaining a file may add an answer,
        # never move one. A pin's `Macros.h` still answers `Macros.h`, even
        # when a retained directory happens to carry that name too.
        pins, _ = self.pins()
        source = pins.sources["retained:../oracle/libstdcxx"]
        source.admitted["bits/Macros.h"] = source.admitted["bits/stl_algo.h"]
        pins.index[("retained:../oracle/libstdcxx", "Macros.h")] = ["bits/Macros.h"]
        self.assertEqual(
            resolvable(pins, (CORE,), "", "Macros.h"),
            [(CORE, "src/openms/include/OpenMS/CONCEPT/Macros.h")],
        )

    def test_two_retained_directories_of_one_name_are_told_apart(self):
        # Both bundles retain a directory called `libstdcxx`. Labelling each by
        # its last component alone put both under one name in the tally - a
        # checker reporting about itself exactly the confusion it exists to
        # report about the documents.
        first = self.retained()
        second = self.retained()
        first.name = "sne-completion/libstdcxx"
        second.name = "a2-textfmt-linux/libstdcxx"
        self.assertNotEqual(str(first), str(second))
        self.assertEqual(str(first), "retained sne-completion/libstdcxx")

    def test_a_retained_file_is_checked_like_any_other(self):
        # It resolves, so it is read back: a range past its end is a finding,
        # which is the point of resolving it rather than skipping it.
        pins, _ = self.pins()
        report = tally()
        found = check_file(
            pins, (CORE,), ("bits/", "stl_algo.h"),
            [(1, 99, "bits/stl_algo.h:1-99", True)], (), (), report,
        )
        self.assertIn("the file has 2 lines", "".join(p for _, ps in found for p in ps))
        self.assertEqual(report["skipped"], 0)


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

    def test_a_fixture_table_is_cut_into_its_rows(self):
        # One row is one record: a citation in it answers for that row, and the
        # pattern is never put to a whole table at once.
        text = "case\torigin\nfirst\tA.cpp:1\n\nsecond\tB.cpp:2\n"
        self.assertEqual(
            self.spans("t.tsv", text), ["case\torigin", "first\tA.cpp:1", "second\tB.cpp:2"]
        )

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

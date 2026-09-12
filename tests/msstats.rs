// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FORMAT/MSstatsFile.h`: the long-format MSstats and MSstatsTMT writers.
//!
//! The two class-test sections carry no assertion at all — both bodies are the
//! comment "tested via MSstatsConverter tool" — so the expectations here come
//! from the retained C++ outputs of the upstream `TOPP_MSstatsConverter_*`
//! tests instead, which is stronger evidence than the class test could give.
//! See `tests/data/msstats_provenance.json`.

#![cfg(feature = "consensusxml")]

use openms::format::msstats::{
    self, IsoOptions, LfqOptions, RetentionTimeSummarization, assemble_run_map,
    check_condition_iso, check_condition_lfq,
};
use openms::format::{consensusxml, experimental_design_file};
use openms::metadata::{ExperimentalDesign, MSFileSectionEntry, SampleSection};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// The shipped consensusXML read budget is not enough for these fixtures: the
/// 471-row TMT map exhausts `max_work` while parsing its modified peptide
/// sequences (each parse charges the modification registry's size), and the
/// twelve-run label-free map exhausts `max_payload_bytes`. Both are
/// caller-supplied ceilings, so the test raises them rather than shrinking the
/// retained C++ input.
fn load_map(name: &str) -> openms::kernel::ConsensusMap {
    consensusxml::load_with_options(
        data(name),
        &consensusxml::ReadOptions {
            max_payload_bytes: 4 * 1024 * 1024 * 1024,
            max_work: 20_000_000_000,
            ..Default::default()
        },
    )
    .unwrap()
}

fn design(name: &str) -> ExperimentalDesign {
    experimental_design_file::load(data(name), &Default::default()).unwrap()
}

/// The upstream suite compares MSstatsConverter output with `FuzzyDiff`
/// (`packages/test-data/topp/CMakeLists.txt:1643`, `:1647`, `:1651`), whose
/// configured tolerance is a 1% ratio or a 0.01 absolute difference. A byte
/// comparison is therefore *not* the upstream contract, and indeed is not
/// satisfiable: the retained `Intensity` fields are written in a fixed
/// six-decimal form that the current `StringUtils::toStr(float)` no longer
/// produces (it emits shortest-round-trip scientific notation above 1e4), so
/// the reference file and today's C++ already disagree textually while naming
/// the same numbers.
///
/// This comparison is much tighter than upstream's: every field is compared
/// exactly unless both sides parse as numbers, in which case the relative
/// difference must stay below 1e-6, with a 1e-9 absolute floor for values near
/// zero. Single-precision intensities are the only quantity that needs any
/// tolerance at all, and 1e-6 is the f32 decimal resolution.
fn assert_matches_reference(produced: &[String], reference: &str) {
    let expected: Vec<&str> = reference.lines().collect();
    assert_eq!(
        produced.len(),
        expected.len(),
        "row count differs: produced {} rows, reference has {}",
        produced.len(),
        expected.len()
    );
    for (row, (got, want)) in produced.iter().zip(&expected).enumerate() {
        let got_fields: Vec<&str> = got.split(',').collect();
        let want_fields: Vec<&str> = want.split(',').collect();
        assert_eq!(
            got_fields.len(),
            want_fields.len(),
            "row {row} has {} fields, reference has {}\n  got:  {got}\n  want: {want}",
            got_fields.len(),
            want_fields.len()
        );
        for (column, (a, b)) in got_fields.iter().zip(&want_fields).enumerate() {
            if a == b {
                continue;
            }
            match (a.parse::<f64>(), b.parse::<f64>()) {
                (Ok(x), Ok(y)) => {
                    let difference = (x - y).abs();
                    let scale = x.abs().max(y.abs());
                    assert!(
                        difference <= 1e-9 || difference / scale <= 1e-6,
                        "row {row} column {column}: {a} != {b}\n  got:  {got}\n  want: {want}"
                    );
                }
                _ => panic!(
                    "row {row} column {column}: {a:?} != {b:?}\n  got:  {got}\n  want: {want}"
                ),
            }
        }
    }
}

/// Tier 1: the retained output of `TOPP_MSstatsConverter_2`, an isobaric
/// (MSstatsTMT) conversion of a ten-channel TMT consensusXML with two
/// fractions in one fraction group.
#[test]
fn iso_reproduces_the_retained_cpp_output() {
    let map = load_map("msstats_iso_in.consensusXML.gz");
    let report = msstats::prepare_iso(
        &map,
        &design("msstats_iso_design.tsv"),
        &IsoOptions::default(),
    )
    .unwrap();
    let reference = std::fs::read_to_string(data("msstats_iso_expected.csv")).unwrap();
    assert_matches_reference(&report.lines, &reference);
    // The source declares its MSstats-run to fraction-group map in storeISO and
    // never fills it, so nothing is reported.
    assert!(report.run_to_fraction_group.is_empty());
}

/// Tier 1: the retained output of `TOPP_MSstatsConverter_3`, the same
/// consensusXML under a design that splits the two files into two fraction
/// groups of one fraction each, which changes `Run` and `TechRepMixture`.
#[test]
fn iso_reproduces_the_retained_cpp_output_for_a_second_design() {
    let map = load_map("msstats_iso_in.consensusXML.gz");
    let report = msstats::prepare_iso(
        &map,
        &design("msstats_iso_fractiongroup_design.tsv"),
        &IsoOptions::default(),
    )
    .unwrap();
    let reference =
        std::fs::read_to_string(data("msstats_iso_fractiongroup_expected.csv")).unwrap();
    assert_matches_reference(&report.lines, &reference);
}

/// Tier 1: the retained output of `TOPP_MSstatsConverter_1`, a label-free
/// conversion of a twelve-run consensusXML with `max` summarization.
#[test]
fn lfq_reproduces_the_retained_cpp_output() {
    let map = load_map("msstats_lfq_in.consensusXML.gz");
    let report = msstats::prepare_lfq(
        &map,
        &design("msstats_lfq_design.tsv"),
        &LfqOptions {
            retention_time_summarization: RetentionTimeSummarization::Max,
            ..Default::default()
        },
    )
    .unwrap();
    let reference = std::fs::read_to_string(data("msstats_lfq_expected.csv")).unwrap();
    assert_matches_reference(&report.lines, &reference);
    // Twelve single-fraction runs, each its own fraction group.
    assert_eq!(report.run_to_fraction_group.len(), 12);
    assert_eq!(report.run_to_fraction_group.get(&1), Some(&1));
}

/// The header the label-free layout writes, and the two optional columns.
#[test]
fn lfq_header_matches_the_retained_reference_header() {
    let reference = std::fs::read_to_string(data("msstats_lfq_expected.csv")).unwrap();
    let map = load_map("msstats_lfq_in.consensusXML.gz");
    let design = design("msstats_lfq_design.tsv");
    let report = msstats::prepare_lfq(&map, &design, &LfqOptions::default()).unwrap();
    assert_eq!(report.lines[0], reference.lines().next().unwrap());
    assert_eq!(
        report.lines[0],
        "ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,\
         IsotopeLabelType,Condition,BioReplicate,Run,Intensity,Reference"
    );
    // The design is not fractionated, so no Fraction column is written; manual
    // summarization adds the leading RetentionTime column.
    let manual = msstats::prepare_lfq(
        &map,
        &design,
        &LfqOptions {
            retention_time_summarization: RetentionTimeSummarization::Manual,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(manual.lines[0].starts_with("RetentionTime,ProteinName,"));
    assert!(
        manual
            .warnings
            .iter()
            .any(|w| w.contains("rt_summarization set to manual"))
    );
}

/// `is_isotope_label_type` selects `H` for the whole file, where the default
/// is the `L` MSstats documents for endogenous peptides.
#[test]
fn the_isotope_label_type_column_switches_between_l_and_h() {
    let map = load_map("msstats_lfq_in.consensusXML.gz");
    let design = design("msstats_lfq_design.tsv");
    let column = |report: &msstats::MSstatsReport| -> Vec<String> {
        report.lines[1..]
            .iter()
            .map(|line| line.split(',').nth(5).unwrap().to_owned())
            .collect()
    };
    let endogenous = msstats::prepare_lfq(&map, &design, &LfqOptions::default()).unwrap();
    assert!(column(&endogenous).iter().all(|value| value == "L"));
    let labelled = msstats::prepare_lfq(
        &map,
        &design,
        &LfqOptions {
            is_isotope_label_type: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(column(&labelled).iter().all(|value| value == "H"));
    assert_eq!(labelled.lines.len(), endogenous.lines.len());
}

/// `reannotate_filenames` replaces the consensus map's own run paths, and is
/// checked against the design like any other set of run names.
///
/// One name is consumed per *column header*, not per distinct file: the source
/// walks the column headers and pops one raw path for each, so a ten-channel
/// TMT map with two files needs twenty names. A shorter list leaves the
/// remaining columns with an empty name, which then fails the design check.
#[test]
fn reannotated_filenames_replace_the_maps_own_run_paths() {
    let map = load_map("msstats_iso_in.consensusXML.gz");
    let design = design("msstats_iso_design.tsv");
    // The same basenames reached by absolute paths: the source takes the
    // basename of every reannotated name, so the result is unchanged.
    let relocated: Vec<String> = map
        .primary_ms_run_path()
        .iter()
        .map(|name| format!("/elsewhere/{name}"))
        .collect();
    assert_eq!(relocated.len(), 20);
    let report = msstats::prepare_iso(
        &map,
        &design,
        &IsoOptions {
            reannotate_filenames: relocated,
            ..Default::default()
        },
    )
    .unwrap();
    let reference = std::fs::read_to_string(data("msstats_iso_expected.csv")).unwrap();
    assert_matches_reference(&report.lines, &reference);
    // A reannotated name the design does not declare is refused.
    let error = msstats::prepare_iso(
        &map,
        &design,
        &IsoOptions {
            reannotate_filenames: vec!["not_in_the_design.mzML".into()],
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(
        format!("{error}").contains("not the same as in the experimental design"),
        "unexpected error: {error}"
    );
}

/// The isobaric header, and the source's forced reversion to manual.
#[test]
fn iso_reverts_a_non_manual_summarization_with_a_warning() {
    let map = load_map("msstats_iso_in.consensusXML.gz");
    let report = msstats::prepare_iso(
        &map,
        &design("msstats_iso_design.tsv"),
        &IsoOptions {
            retention_time_summarization: RetentionTimeSummarization::Sum,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        report.lines[0],
        "RetentionTime,ProteinName,PeptideSequence,Charge,Channel,Condition,BioReplicate,Run,\
         Mixture,TechRepMixture,Fraction,Intensity,Reference"
    );
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("Reverting to 'manual'"))
    );
    // Reverting means the rows are identical to the manual conversion.
    let reference = std::fs::read_to_string(data("msstats_iso_expected.csv")).unwrap();
    assert_matches_reference(&report.lines, &reference);
}

/// A label-free conversion refuses a design with more than one label.
#[test]
fn lfq_refuses_a_multi_label_design() {
    let map = load_map("msstats_iso_in.consensusXML.gz");
    let error = msstats::prepare_lfq(
        &map,
        &design("msstats_iso_design.tsv"),
        &LfqOptions::default(),
    )
    .unwrap_err();
    assert!(
        format!("{error}").contains("Too many labels"),
        "unexpected error: {error}"
    );
}

/// Keeping shared peptides adds rows for peptides that map into more than one
/// indistinguishable protein group.
#[test]
fn keeping_shared_peptides_adds_rows() {
    let map = load_map("msstats_iso_in.consensusXML.gz");
    let design = design("msstats_iso_design.tsv");
    let dropped = msstats::prepare_iso(&map, &design, &IsoOptions::default()).unwrap();
    let kept = msstats::prepare_iso(
        &map,
        &design,
        &IsoOptions {
            remove_shared_peptides: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(dropped.shared_peptides_dropped > 0);
    assert_eq!(kept.shared_peptides_dropped, 0);
    assert!(kept.lines.len() > dropped.lines.len());
    assert!(
        dropped
            .warnings
            .iter()
            .any(|w| w.contains("shared peptides"))
    );
}

fn factor_design(factors: &[&str], values: &[&[&str]]) -> ExperimentalDesign {
    let mut columns = BTreeMap::new();
    for (index, factor) in factors.iter().enumerate() {
        columns.insert((*factor).to_owned(), index);
    }
    let mut rows = BTreeMap::new();
    let mut content = Vec::new();
    let mut files = Vec::new();
    for (index, row) in values.iter().enumerate() {
        rows.insert((index + 1).to_string(), index);
        content.push(row.iter().map(|value| (*value).to_owned()).collect());
        files.push(MSFileSectionEntry {
            fraction_group: u32::try_from(index + 1).unwrap(),
            fraction: 1,
            path: format!("run{}.mzML", index + 1),
            label: 1,
            sample: u32::try_from(index).unwrap(),
            sample_name: (index + 1).to_string(),
        });
    }
    let section = SampleSection::from_table(content, rows, columns).unwrap();
    ExperimentalDesign::from_sections(files, section).unwrap()
}

/// A missing MSstats column is refused, naming the column.
#[test]
fn a_missing_condition_or_bioreplicate_column_is_refused() {
    let design = factor_design(
        &["MSstats_Condition", "MSstats_BioReplicate"],
        &[&["1", "1"], &["2", "2"]],
    );
    assert!(
        check_condition_lfq(
            design.sample_section(),
            "MSstats_BioReplicate",
            "MSstats_Condition"
        )
        .unwrap()
        .is_empty()
    );
    let error =
        check_condition_lfq(design.sample_section(), "Absent", "MSstats_Condition").unwrap_err();
    assert!(format!("{error}").contains("Absent"), "{error}");
    // The isobaric layout additionally needs the mixture column.
    let error = check_condition_iso(
        design.sample_section(),
        "MSstats_BioReplicate",
        "MSstats_Condition",
        "MSstats_Mixture",
    )
    .unwrap_err();
    assert!(format!("{error}").contains("MSstats_Mixture"), "{error}");
}

/// A biological replicate that recurs under two conditions is warned about,
/// not refused.
#[test]
fn a_recurring_bioreplicate_is_warned_about() {
    let design = factor_design(
        &["MSstats_Condition", "MSstats_BioReplicate"],
        &[&["1", "7"], &["2", "7"]],
    );
    let warnings = check_condition_lfq(
        design.sample_section(),
        "MSstats_BioReplicate",
        "MSstats_Condition",
    )
    .unwrap();
    assert_eq!(warnings.len(), 1);
    assert!(
        warnings[0].contains("'7' occurs under 2 different"),
        "{:?}",
        warnings
    );
    assert!(warnings[0].contains("paired design"), "{:?}", warnings);
}

/// A non-ASCII factor value survives into the condition and biological
/// replicate columns, and into the recurring-replicate warning, without being
/// split mid-codepoint.
#[test]
fn non_ascii_factor_values_are_carried_through_whole() {
    let design = factor_design(
        &["MSstats_Condition", "MSstats_BioReplicate"],
        &[
            &["\u{6761}\u{4ef6}", "\u{8907}\u{88fd}"],
            &["control", "\u{8907}\u{88fd}"],
        ],
    );
    let warnings = check_condition_lfq(
        design.sample_section(),
        "MSstats_BioReplicate",
        "MSstats_Condition",
    )
    .unwrap();
    assert_eq!(warnings.len(), 1);
    assert!(warnings[0].contains('\u{8907}'), "{:?}", warnings);
    assert!(warnings[0].contains('\u{6761}'), "{:?}", warnings);
    assert_eq!(
        design
            .sample_section()
            .factor_value_by_row(0, "MSstats_Condition")
            .unwrap(),
        "\u{6761}\u{4ef6}"
    );
}

/// MSstats numbers runs by `(file basename, fraction)` in design order, from
/// one, while OpenMS splits one run into fractions.
#[test]
fn run_numbers_enumerate_file_and_fraction_pairs() {
    let iso = assemble_run_map(&design("msstats_iso_design.tsv"));
    assert_eq!(iso.len(), 2);
    assert_eq!(
        iso.get(&("QExactiveHF02_03904.mzML".to_owned(), 1)),
        Some(&1)
    );
    assert_eq!(
        iso.get(&("QExactiveHF02_03905.mzML".to_owned(), 2)),
        Some(&2)
    );
    let lfq = assemble_run_map(&design("msstats_lfq_design.tsv"));
    assert_eq!(lfq.len(), 12);
    assert_eq!(
        lfq.get(&("JD_06232014_sample1-A.mzML".to_owned(), 1)),
        Some(&1)
    );
}

/// Summarization names, and the refusal the source lacks.
#[test]
fn summarization_names_round_trip_and_unknown_names_are_refused() {
    for name in ["manual", "max", "min", "mean", "sum"] {
        assert_eq!(
            RetentionTimeSummarization::from_name(name).unwrap().name(),
            name
        );
    }
    assert!(RetentionTimeSummarization::from_name("median").is_err());
    // The source would have written the intensity 0 for this name.
    assert!(RetentionTimeSummarization::from_name("").is_err());
}

/// Writing to a path produces exactly the prepared rows, newline terminated.
#[test]
fn stored_output_is_the_prepared_rows() {
    let map = load_map("msstats_iso_in.consensusXML.gz");
    let design = design("msstats_iso_design.tsv");
    let directory = std::env::temp_dir().join(format!("openms-msstats-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("iso.csv");
    let report = msstats::store_iso(&path, &map, &design, &IsoOptions::default()).unwrap();
    let written = std::fs::read_to_string(&path).unwrap();
    assert_eq!(written.lines().count(), report.lines.len());
    assert!(written.ends_with('\n'));
    assert_eq!(written.lines().next().unwrap(), report.lines[0]);
    std::fs::remove_file(&path).unwrap();
    std::fs::remove_dir_all(&directory).unwrap();
}

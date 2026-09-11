// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source `ExperimentalDesign_test.cpp` literals over the retained fixtures,
//! plus the mapping and error paths that test only exercises indirectly.

use openms::format::experimental_design_file::{self as design_file, ReadOptions};
use openms::identification::ProteinIdentification;
use openms::kernel::{ColumnHeader, ConsensusMap, FeatureMap};
use openms::metadata::{ExperimentalDesign, MSFileSectionEntry, MetaValue, SampleSection};
use std::collections::{BTreeMap, BTreeSet};

fn fixture(name: &str) -> String {
    format!("tests/data/experimental_design_{name}.tsv")
}
fn load(name: &str) -> ExperimentalDesign {
    design_file::load(fixture(name), &ReadOptions::default()).unwrap()
}
fn labelfree() -> [ExperimentalDesign; 3] {
    [
        load("input_1"),
        load("input_1_single_table"),
        load("input_3_single_table"),
    ]
}
fn fourplex() -> [ExperimentalDesign; 2] {
    [load("input_2"), load("input_2_single_table")]
}
fn names(set: &BTreeSet<String>) -> Vec<&str> {
    set.iter().map(String::as_str).collect()
}

#[test]
fn source_counts_over_both_table_layouts() {
    for design in labelfree() {
        assert_eq!(design.number_of_samples(), 12);
        assert_eq!(design.number_of_fractions(), 1);
        assert_eq!(design.number_of_labels(), 1);
        assert_eq!(design.number_of_ms_files(), 12);
        assert_eq!(design.number_of_fraction_groups(), 12);
        assert!(!design.is_fractionated());
        assert!(design.same_nr_of_ms_files_per_fraction());
        assert_eq!(design.sample(1, 1).unwrap(), 0);
        assert_eq!(design.sample(12, 1).unwrap(), 11);
    }
    for design in fourplex() {
        assert_eq!(design.number_of_samples(), 8);
        assert_eq!(design.number_of_fractions(), 3);
        assert_eq!(design.number_of_labels(), 4);
        assert_eq!(design.number_of_ms_files(), 6);
        assert_eq!(design.number_of_fraction_groups(), 2);
        assert!(design.is_fractionated());
        assert!(design.same_nr_of_ms_files_per_fraction());
        assert_eq!(design.sample(1, 1).unwrap(), 0);
        assert_eq!(design.sample(2, 4).unwrap(), 7);
    }
}

#[test]
fn unknown_fraction_group_and_label_is_reported() {
    assert!(load("input_1").sample(99999, 99999).is_err());
}

#[test]
fn source_fraction_to_ms_files_mapping() {
    for design in labelfree() {
        let mapping = design.fraction_to_ms_files_mapping();
        assert_eq!(mapping.len(), 1);
        assert_eq!(mapping[&1].len(), 12);
    }
    for design in fourplex() {
        let mapping = design.fraction_to_ms_files_mapping();
        assert_eq!(mapping.len(), 3);
        for fraction in 1..=3 {
            assert_eq!(mapping[&fraction].len(), 8);
        }
    }
}

#[test]
fn source_path_label_mappings() {
    for design in labelfree() {
        assert_eq!(design.path_label_to_sample_mapping(true).unwrap().len(), 12);
        let fractions = design.path_label_to_fraction_mapping(true).unwrap();
        assert_eq!(fractions.len(), 12);
        assert!(fractions.values().all(|fraction| *fraction == 1));
    }
    for design in fourplex() {
        let samples = design.path_label_to_sample_mapping(true).unwrap();
        assert_eq!(samples.len(), 24);
        assert!(samples.values().all(|sample| *sample <= 7));
        let fractions = design.path_label_to_fraction_mapping(true).unwrap();
        assert_eq!(fractions.len(), 24);
        assert!(fractions.values().all(|f| (1..=3).contains(f)));
        // Fraction group follows the technical-replicate suffix of the file.
        for ((path, _), group) in design.path_label_to_fraction_group_mapping(true).unwrap() {
            assert_eq!(group, if path.contains("TR2") { 2 } else { 1 });
        }
    }
    // Twelve unfractionated label-free files are twelve consecutive groups, in
    // the lexical order of their basenames.
    for design in [load("input_1"), load("input_1_single_table")] {
        let groups = design.path_label_to_fraction_group_mapping(true).unwrap();
        assert_eq!(groups.len(), 12);
        assert!(groups.values().copied().eq(1..=12));
    }
}

#[test]
fn ambiguous_basename_and_label_is_rejected() {
    let rows = vec![
        MSFileSectionEntry {
            fraction_group: 1,
            path: "/tmp/run_a/shared_name.mzML".into(),
            sample: 0,
            sample_name: "S1".into(),
            ..Default::default()
        },
        MSFileSectionEntry {
            fraction_group: 2,
            path: "/tmp/run_b/shared_name.mzML".into(),
            sample: 1,
            sample_name: "S2".into(),
            ..Default::default()
        },
    ];
    let mut section = SampleSection::new();
    section.add_sample("S1", Vec::new());
    section.add_sample("S2", Vec::new());
    let design = ExperimentalDesign::from_sections(rows, section).unwrap();
    assert!(design.path_label_to_sample_mapping(true).is_err());
    // The full paths stay distinct, so the same design maps without basenames.
    assert_eq!(design.path_label_to_sample_mapping(false).unwrap().len(), 2);
}

#[test]
fn source_factors_and_factor_values_agree_across_layouts() {
    let designs = labelfree();
    for design in &designs {
        assert_eq!(
            design.sample_section().factors().collect::<Vec<_>>(),
            ["MSstats_BioReplicate", "MSstats_Condition", "Sample"]
        );
    }
    // A missing Sample column is filled from the fraction group, so the three
    // layouts describe the same twelve samples.
    for row in 0..12 {
        for factor in ["MSstats_BioReplicate", "MSstats_Condition", "Sample"] {
            let values: Vec<&str> = designs
                .iter()
                .map(|design| {
                    design
                        .sample_section()
                        .factor_value_by_row(row, factor)
                        .unwrap()
                })
                .collect();
            assert_eq!(values[0], values[1], "row {row} factor {factor}");
            assert_eq!(values[0], values[2], "row {row} factor {factor}");
        }
    }
    assert_eq!(
        load("input_1")
            .sample_section()
            .factor_value("4", "MSstats_Condition")
            .unwrap(),
        "2"
    );
    for design in fourplex() {
        assert_eq!(
            design.sample_section().factors().collect::<Vec<_>>(),
            ["Sample"]
        );
    }
}

#[test]
fn conditions_group_samples_and_ignore_replicate_columns() {
    let design = load("input_1");
    // Four conditions of three biological replicates each; the replicate column
    // is excluded by name, the condition column is not.
    let conditions = design.condition_to_sample_mapping().unwrap();
    assert_eq!(conditions.len(), 4);
    assert_eq!(
        conditions.keys().cloned().collect::<Vec<_>>(),
        [vec!["1"], vec!["2"], vec!["3"], vec!["4"]]
    );
    assert_eq!(
        conditions[&vec!["1".to_string()]],
        BTreeSet::from([0, 1, 2])
    );
    assert_eq!(
        conditions[&vec!["4".to_string()]],
        BTreeSet::from([9, 10, 11])
    );

    let by_name = design.sample_to_condition_mapping().unwrap();
    assert_eq!(by_name.len(), 12);
    assert_eq!(by_name["1"], 0);
    assert_eq!(by_name["12"], 3);

    // Every factor counts for prefractionation, so each sample stays its own.
    let prefractionation = design.sample_to_prefractionation_mapping().unwrap();
    assert_eq!(prefractionation.len(), 12);
    assert_eq!(prefractionation.values().collect::<BTreeSet<_>>().len(), 12);

    let unique = design.unique_sample_row_to_sample_mapping().unwrap();
    assert_eq!(unique.len(), 12);
    assert_eq!(
        names(&unique[&vec!["1".to_string(), "1".to_string()]]),
        ["1"]
    );

    let per_condition = design.condition_to_path_label_vector().unwrap();
    assert_eq!(per_condition.len(), 4);
    for entries in &per_condition {
        assert_eq!(entries.len(), 3);
    }
    assert!(per_condition[0][0].0.ends_with("JD_06232014_sample1-A.raw"));

    let paths = design.path_label_to_condition_mapping(true).unwrap();
    assert_eq!(paths.len(), 12);
    assert_eq!(paths[&("JD_06232014_sample1-A.raw".to_string(), 1)], 0);
    assert_eq!(paths[&("JD_06232014_sample4_C.raw".to_string(), 1)], 3);
    assert_eq!(
        design
            .path_label_to_prefractionation_mapping(true)
            .unwrap()
            .len(),
        12
    );
}

#[test]
fn a_factorless_section_gives_every_sample_its_own_group() {
    // The shape produced by every from_* constructor.
    let design = ExperimentalDesign::from_identifications(&[run(&["a.mzML", "b.mzML"])]).unwrap();
    assert_eq!(design.sample_section().factors().len(), 0);
    assert_eq!(
        design.sample_to_condition_mapping().unwrap(),
        BTreeMap::from([("0".to_string(), 0), ("1".to_string(), 1)])
    );
    assert_eq!(
        design.sample_to_prefractionation_mapping().unwrap(),
        BTreeMap::from([("0".to_string(), 0), ("1".to_string(), 1)])
    );
    // The other three collapse a factor-less section into one condition.
    assert_eq!(design.condition_to_sample_mapping().unwrap().len(), 1);
    assert_eq!(design.condition_to_path_label_vector().unwrap().len(), 1);
    assert_eq!(
        design.path_label_to_condition_mapping(true).unwrap().len(),
        2
    );
}

#[test]
fn source_validation_rules() {
    // Missing fractions and unsorted input rows are accepted.
    load("input_2_wrong");
    load("input_2_wrong_3");
    // Fraction groups must still be consecutive from 1.
    assert!(design_file::load(fixture("input_2_wrong_2"), &ReadOptions::default()).is_err());

    let duplicate = vec![MSFileSectionEntry::default(), MSFileSectionEntry::default()];
    assert!(ExperimentalDesign::from_sections(duplicate, SampleSection::new()).is_err());

    // One fraction group may not hold two samples when there is one label.
    let two_samples = vec![
        MSFileSectionEntry {
            fraction: 1,
            path: "a.mzML".into(),
            sample: 0,
            ..Default::default()
        },
        MSFileSectionEntry {
            fraction: 2,
            path: "b.mzML".into(),
            sample: 1,
            ..Default::default()
        },
    ];
    assert!(ExperimentalDesign::from_sections(two_samples.clone(), SampleSection::new()).is_err());
    // The same rows are valid once the two samples are separate channels.
    let mut labeled = two_samples;
    labeled[1].label = 2;
    labeled[1].fraction = 1;
    assert!(ExperimentalDesign::from_sections(labeled, SampleSection::new()).is_ok());

    // An empty file section is valid and unchecked.
    assert!(ExperimentalDesign::from_sections(Vec::new(), SampleSection::new()).is_ok());
}

#[test]
fn comment_lines_are_skipped_and_rows_are_sorted() {
    let design = load("BSA_design_onetable_nonconsec");
    assert_eq!(design.ms_file_section().len(), 6);
    assert_eq!(design.number_of_fraction_groups(), 4);
    assert_eq!(design.number_of_samples(), 3);
    let groups: Vec<u32> = design
        .ms_file_section()
        .iter()
        .map(|row| row.fraction_group)
        .collect();
    assert_eq!(groups, [1, 2, 2, 3, 3, 4]);
    // A sample reused across fraction groups keeps one sample row.
    let mapping = design.path_label_to_fraction_group_mapping(true).unwrap();
    assert_eq!(mapping[&("BSA1_F1.mzML".to_string(), 1)], 1);
    assert_eq!(mapping[&("BSA1_F2.mzML".to_string(), 1)], 4);
    assert_eq!(
        design.sample_section().samples().collect::<Vec<_>>(),
        ["BSA1", "BSA2", "BSA3"]
    );
    // "Rep1" is a replicate column by name, so it does not split conditions.
    assert_eq!(design.condition_to_sample_mapping().unwrap().len(), 3);
}

#[test]
fn filtering_by_basename_rebuilds_the_sample_section() {
    let mut design = load("input_1");
    let keep = BTreeSet::from([
        "JD_06232014_sample1-A.raw".to_string(),
        "JD_06232014_sample1_B.raw".to_string(),
    ]);
    assert_eq!(design.filter_by_basenames(&keep).unwrap(), 10);
    assert_eq!(design.ms_file_section().len(), 2);
    assert_eq!(design.number_of_samples(), 2);
    // Sample indices are renumbered to the surviving rows, values retained.
    assert_eq!(design.ms_file_section()[0].sample, 0);
    assert_eq!(design.ms_file_section()[1].sample, 1);
    assert_eq!(
        design
            .sample_section()
            .factor_value("2", "MSstats_Condition")
            .unwrap(),
        "1"
    );
    // A disjoint set empties the file section rather than failing.
    let mut empty = load("input_1");
    assert_eq!(
        empty
            .filter_by_basenames(&BTreeSet::from(["nothing.raw".to_string()]))
            .unwrap(),
        12
    );
    assert!(empty.ms_file_section().is_empty());
}

#[test]
fn one_table_needs_a_sample_column_for_multiplexed_data() {
    let text = "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\n1\t1\ta.mzML\t2\n";
    assert!(parse(text).is_err());
    // Label 1 keeps working: the fraction group becomes the sample name.
    let design = parse("Fraction_Group\tFraction\tSpectra_Filepath\tLabel\n1\t1\ta.mzML\t1\n")
        .unwrap()
        .0;
    assert_eq!(design.ms_file_section()[0].sample_name, "1");
}

#[test]
fn two_table_rejects_unknown_file_section_headers_and_unknown_samples() {
    let unknown_header = "Fraction_Group\tFraction\tSpectra_Filepath\tSample\tExtra\n\
         1\t1\ta.mzML\t1\tx\n\nSample\tCondition\n1\tc\n";
    assert!(parse(unknown_header).is_err());

    let unknown_sample = "Fraction_Group\tFraction\tSpectra_Filepath\tSample\n\
         1\t1\ta.mzML\t7\n\nSample\tCondition\n1\tc\n";
    assert!(parse(unknown_sample).is_err());

    let duplicate_sample = "Fraction_Group\tFraction\tSpectra_Filepath\tSample\n\
         1\t1\ta.mzML\t1\n\nSample\tCondition\n1\tc\n1\td\n";
    assert!(parse(duplicate_sample).is_err());
}

#[test]
fn ragged_rows_are_rejected_rather_than_read_out_of_bounds() {
    // Source reads past the end of a short row in both layouts; see CPP-059.
    let short_one_table =
        "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\n1\t1\ta.mzML\n";
    assert!(parse(short_one_table).is_err());
    let short_sample_row = "Fraction_Group\tFraction\tSpectra_Filepath\tSample\n\
         1\t1\ta.mzML\t1\n\nSample\tCondition\n1\n";
    assert!(parse(short_sample_row).is_err());
    // A negative index wraps to a huge unsigned value in the source; see CPP-060.
    let negative = "Fraction_Group\tFraction\tSpectra_Filepath\tSample\n1\t-1\ta.mzML\t1\n";
    assert!(parse(negative).is_err());
    let duplicate_header =
        "Fraction_Group\tFraction\tSpectra_Filepath\tFraction\n1\t1\ta.mzML\t1\n";
    assert!(parse(duplicate_header).is_err());
    let missing_header = "Fraction\tSpectra_Filepath\n1\ta.mzML\n";
    assert!(parse(missing_header).is_err());
}

#[test]
fn mismatched_factors_for_one_sample_are_reported_not_fatal() {
    let text = "Fraction_Group\tFraction\tSpectra_Filepath\tSample\tCondition\n\
         1\t1\ta.mzML\t1\tcontrol\n\
         2\t1\tb.mzML\t1\ttreated\n";
    let (design, warnings) = parse(text).unwrap();
    assert_eq!(warnings.len(), 1);
    assert!(warnings[0].contains("factors for sample '1'"));
    // The first row's values are the ones retained.
    assert_eq!(
        design
            .sample_section()
            .factor_value("1", "Condition")
            .unwrap(),
        "control"
    );
}

#[test]
fn a_missing_spectra_file_only_fails_when_required() {
    let text = "Fraction_Group\tFraction\tSpectra_Filepath\tSample\n1\t1\tabsent.mzML\t1\n";
    let design = parse(text).unwrap().0;
    // Nothing resolved, so the string is kept as written.
    assert_eq!(design.ms_file_section()[0].path, "absent.mzML");
    let text_file = openms::format::text::TextFile::from_reader(
        text.as_bytes(),
        &openms::format::text::ReadOptions {
            trim_lines: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(
        design_file::load_text(
            &text_file,
            "tests/data/design.tsv",
            &ReadOptions {
                require_spectra_files: true,
                ..Default::default()
            },
        )
        .is_err()
    );
}

#[test]
fn an_empty_file_loads_as_an_empty_design() {
    let design = parse("").unwrap().0;
    assert!(design.ms_file_section().is_empty());
    assert_eq!(design.number_of_samples(), 0);
    assert_eq!(design.number_of_labels(), 0);
}

#[test]
fn sample_section_answers_names_from_either_store() {
    let mut added = SampleSection::new();
    added.add_sample("BSA1", Vec::new());
    added.add_sample("BSA2", Vec::new());
    assert_eq!(added.sample_name(0).unwrap(), "BSA1");
    assert_eq!(added.sample_name(1).unwrap(), "BSA2");
    assert_eq!(added.sample_row("BSA2").unwrap(), 1);
    assert!(added.sample_name(7).is_err());
    assert!(added.sample_row("nope").is_err());

    let from_file = SampleSection::from_table(
        vec![
            vec!["S_a".into(), "control".into()],
            vec!["S_b".into(), "treated".into()],
        ],
        BTreeMap::from([("S_a".to_string(), 0), ("S_b".to_string(), 1)]),
        BTreeMap::from([
            ("Sample".to_string(), 0),
            ("MSstats_Condition".to_string(), 1),
        ]),
    )
    .unwrap();
    assert_eq!(from_file.sample_name(0).unwrap(), "S_a");
    assert_eq!(from_file.sample_name(1).unwrap(), "S_b");
    assert_eq!(
        from_file.factor_column_index("MSstats_Condition").unwrap(),
        1
    );
    assert!(from_file.has_factor("Sample") && from_file.has_sample("S_b"));
    assert!(from_file.factor_value("S_a", "absent").is_err());
    assert!(from_file.factor_value("absent", "Sample").is_err());
    // A row shorter than the column map is refused up front.
    assert!(
        SampleSection::from_table(
            vec![vec!["S_a".into()]],
            BTreeMap::from([("S_a".to_string(), 0)]),
            BTreeMap::from([("Sample".to_string(), 0), ("Condition".to_string(), 1)]),
        )
        .is_err()
    );
}

#[test]
fn sample_names_that_are_not_row_indices_keep_the_mappings_working() {
    for (a, b) in [("1", "2"), ("BSA1", "BSA2")] {
        let rows = (0..2u32)
            .map(|index| MSFileSectionEntry {
                fraction_group: index + 1,
                path: format!("/data/run_{index}.mzML"),
                sample: index,
                sample_name: if index == 0 { a.into() } else { b.into() },
                ..Default::default()
            })
            .collect();
        let mut section = SampleSection::new();
        section.add_sample(a, Vec::new());
        section.add_sample(b, Vec::new());
        let design = ExperimentalDesign::from_sections(rows, section).unwrap();
        for mapping in [
            design.sample_to_prefractionation_mapping().unwrap(),
            design.sample_to_condition_mapping().unwrap(),
        ] {
            assert_eq!(mapping.len(), design.number_of_samples() as usize);
            assert!(mapping.contains_key(a) && mapping.contains_key(b));
        }
        assert_eq!(
            design
                .path_label_to_prefractionation_mapping(false)
                .unwrap()
                .len(),
            2
        );
        assert_eq!(
            design.path_label_to_condition_mapping(false).unwrap().len(),
            2
        );
    }
}

fn run(paths: &[&str]) -> ProteinIdentification {
    ProteinIdentification {
        primary_ms_run_paths: paths.iter().map(|p| p.to_string()).collect(),
        ..Default::default()
    }
}

#[test]
fn from_identifications_numbers_groups_from_one_and_samples_from_zero() {
    let design =
        ExperimentalDesign::from_identifications(&[run(&["a.mzML", "b.mzML", "c.mzML"])]).unwrap();
    assert_eq!(design.ms_file_section().len(), 3);
    for (index, row) in design.ms_file_section().iter().enumerate() {
        assert_eq!(row.sample as usize, index);
        assert_eq!(row.fraction_group as usize, index + 1);
    }
    // An inferred design must satisfy the same validator as a loaded one.
    let revalidated = ExperimentalDesign::from_sections(
        design.ms_file_section().to_vec(),
        design.sample_section().clone(),
    )
    .unwrap();
    assert_eq!(revalidated.number_of_fraction_groups(), 3);
    assert!(
        ExperimentalDesign::from_identifications(&[])
            .unwrap()
            .ms_file_section()
            .is_empty()
    );
}

fn header(filename: &str, label: &str, channel: Option<i64>) -> ColumnHeader {
    let mut header = ColumnHeader {
        filename: filename.into(),
        label: label.into(),
        ..Default::default()
    };
    if let Some(id) = channel {
        header.metadata.insert("channel_id".into(), id.into());
    }
    header
}

#[test]
fn from_consensus_map_derives_labels_groups_and_samples() {
    // TMT10-plex: one file, ten channels, ten samples, one fraction group.
    let mut map = ConsensusMap::new();
    map.experiment_type = "labeled_MS2".into();
    for channel in 0..10u64 {
        map.column_headers.insert(
            channel,
            header("/data/TMTTenPlex.mzML", "tmt10plex", Some(channel as i64)),
        );
    }
    let design = ExperimentalDesign::from_consensus_map(&map).unwrap();
    assert_eq!(design.number_of_labels(), 10);
    assert_eq!(design.number_of_ms_files(), 1);
    let rows = design.ms_file_section();
    assert_eq!((rows[0].label, rows[9].label), (1, 10));
    assert_eq!((rows[0].fraction_group, rows[9].fraction_group), (1, 1));
    assert_eq!((rows[0].fraction, rows[9].fraction), (1, 1));
    assert_eq!((rows[0].sample, rows[9].sample), (0, 9));
    assert_eq!(rows[0].path, "/data/TMTTenPlex.mzML");

    // Label-free: no channel ids, so each file becomes its own fraction group.
    let mut map = ConsensusMap::new();
    map.column_headers
        .insert(0, header("raw_file1.mzML", "", None));
    map.column_headers
        .insert(1, header("raw_file2.mzML", "", None));
    let design = ExperimentalDesign::from_consensus_map(&map).unwrap();
    assert_eq!(design.number_of_labels(), 1);
    assert_eq!(design.number_of_ms_files(), 2);
    let rows = design.ms_file_section();
    assert_eq!((rows[0].label, rows[1].label), (1, 1));
    assert_eq!((rows[0].fraction_group, rows[1].fraction_group), (1, 2));
    assert_eq!((rows[0].sample, rows[1].sample), (0, 1));
    assert_eq!(rows[0].path, "raw_file1.mzML");
    assert_eq!(rows[1].path, "raw_file2.mzML");

    // Annotated sample names are interned and kept in the sample section.
    let mut map = ConsensusMap::new();
    for (index, name) in ["Sample_A", "Sample_B"].iter().enumerate() {
        let mut column = header(&format!("file_{index}.mzML"), "label-free", None);
        column
            .metadata
            .insert("sample_name".into(), (*name).to_string().into());
        map.column_headers.insert(index as u64, column);
    }
    let design = ExperimentalDesign::from_consensus_map(&map).unwrap();
    assert!(design.sample_section().has_sample("Sample_A"));
    assert!(design.sample_section().has_sample("Sample_B"));
}

#[test]
fn annotated_fractions_need_a_fraction_group() {
    let mut map = ConsensusMap::new();
    let mut column = header("a.mzML", "", None);
    column.metadata.insert("fraction".into(), 2_i64.into());
    map.column_headers.insert(0, column);
    assert!(ExperimentalDesign::from_consensus_map(&map).is_err());

    let mut column = header("a.mzML", "", None);
    column.metadata.insert("fraction".into(), 2_i64.into());
    column
        .metadata
        .insert("fraction_group".into(), 1_i64.into());
    map.column_headers.insert(0, column);
    let design = ExperimentalDesign::from_consensus_map(&map).unwrap();
    assert_eq!(design.ms_file_section()[0].fraction, 2);

    // A nonnumeric index is refused rather than cast from an inactive union.
    let mut column = header("a.mzML", "", None);
    column
        .metadata
        .insert("fraction".into(), MetaValue::from("two".to_string()));
    column
        .metadata
        .insert("fraction_group".into(), 1_i64.into());
    map.column_headers.insert(0, column);
    assert!(ExperimentalDesign::from_consensus_map(&map).is_err());
}

#[test]
fn from_feature_map_needs_exactly_one_primary_run() {
    let mut map = FeatureMap::new();
    assert!(ExperimentalDesign::from_feature_map(&map).is_err());
    map.metadata.insert(
        "spectra_data".into(),
        MetaValue::from(vec!["file://C:/raw_file1.mzML".to_string()]),
    );
    let design = ExperimentalDesign::from_feature_map(&map).unwrap();
    assert_eq!(design.number_of_labels(), 1);
    assert_eq!(design.number_of_ms_files(), 1);
    let row = &design.ms_file_section()[0];
    assert_eq!(
        (row.label, row.fraction_group, row.fraction, row.sample),
        (1, 1, 1, 0)
    );
    assert_eq!(row.path, "file://C:/raw_file1.mzML");

    map.metadata.insert(
        "spectra_data".into(),
        MetaValue::from(vec!["a.mzML".to_string(), "b.mzML".to_string()]),
    );
    assert!(ExperimentalDesign::from_feature_map(&map).is_err());
}

#[test]
fn annotate_column_headers_stamps_fraction_structure() {
    // Two fraction files of one quantification unit, two channels each.
    let mut map = ConsensusMap::new();
    map.experiment_type = "labeled_MS2".into();
    for file in 0..2u64 {
        for channel in 0..2u64 {
            map.column_headers.insert(
                file * 2 + channel,
                header(
                    &format!("/data/run{}.mzML", file + 1),
                    "tmt10plex",
                    Some(channel as i64),
                ),
            );
        }
    }
    let design = parse(
        "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\n\
         1\t1\trun1.mzML\t1\t1\n1\t1\trun1.mzML\t2\t2\n\
         1\t2\trun2.mzML\t1\t1\n1\t2\trun2.mzML\t2\t2\n",
    )
    .unwrap()
    .0;
    assert_eq!(design.annotate_column_headers(&mut map).unwrap(), 0);
    for header in map.column_headers.values() {
        assert_eq!(header.metadata["fraction_group"].as_i64().unwrap(), 1);
    }
    let fraction = |index: u64| {
        map.column_headers[&index].metadata["fraction"]
            .as_i64()
            .unwrap()
    };
    assert_eq!(
        [fraction(0), fraction(1), fraction(2), fraction(3)],
        [1, 1, 2, 2]
    );

    // A header the design does not describe is counted, not silently skipped.
    map.column_headers
        .insert(99, header("/data/not_in_design.mzML", "tmt10plex", Some(0)));
    assert_eq!(design.annotate_column_headers(&mut map).unwrap(), 1);
    assert!(
        !map.column_headers[&99]
            .metadata
            .contains_key("fraction_group")
    );
}

#[test]
fn annotate_column_headers_handles_inferred_and_ambiguous_headers() {
    // A design inferred from the map itself annotates, sample name included.
    let mut map = ConsensusMap::new();
    map.column_headers
        .insert(0, header("/data/run_0.mzML", "label-free", None));
    map.column_headers
        .insert(1, header("/data/run_1.mzML", "label-free", None));
    let inferred = ExperimentalDesign::from_consensus_map(&map).unwrap();
    assert_eq!(inferred.annotate_column_headers(&mut map).unwrap(), 0);
    for index in [0, 1] {
        let metadata = &map.column_headers[&index].metadata;
        assert!(metadata.contains_key("fraction_group"));
        assert!(metadata.contains_key("fraction"));
        assert!(metadata.contains_key("sample_name"));
    }

    // Two headers of one file without channel ids both resolve to label 1, so
    // neither can be attributed to the row they collapse onto.
    let mut map = ConsensusMap::new();
    map.experiment_type = "labeled_MS2".into();
    map.column_headers
        .insert(0, header("/data/run1.mzML", "light", None));
    map.column_headers
        .insert(1, header("/data/run1.mzML", "heavy", None));
    let design = parse(
        "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\n\
         1\t1\trun1.mzML\t1\t1\n1\t1\trun1.mzML\t2\t2\n",
    )
    .unwrap()
    .0;
    assert_eq!(design.annotate_column_headers(&mut map).unwrap(), 2);
    for index in [0, 1] {
        assert!(
            !map.column_headers[&index]
                .metadata
                .contains_key("sample_name")
        );
    }
}

#[test]
fn resource_limits_reject_oversized_sections() {
    let rows = vec![
        MSFileSectionEntry {
            path: "a".repeat(1024),
            ..Default::default()
        };
        ExperimentalDesign::MAX_ROWS + 1
    ];
    assert!(ExperimentalDesign::from_sections(rows, SampleSection::new()).is_err());
}

/// Parse an in-memory design, as the source `load(TextFile, ...)` overload does.
fn parse(text: &str) -> openms::Result<(ExperimentalDesign, Vec<String>)> {
    let lines = openms::format::text::TextFile::from_reader(
        text.as_bytes(),
        &openms::format::text::ReadOptions {
            trim_lines: true,
            ..Default::default()
        },
    )?;
    design_file::load_text(&lines, "tests/data/design.tsv", &ReadOptions::default())
}

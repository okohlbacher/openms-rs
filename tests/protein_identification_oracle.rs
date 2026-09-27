// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source: METADATA/ProteinIdentification.h/.cpp and ProteinIdentification_test.cpp
// at core bc9cc12. Expected values come from the upstream class test or from the
// Release oracle `protein-identification-closure` (Linux x86_64 Release build on
// ibminode06), whose records are tests/data/protein_identification_oracle.tsv; see
// tests/data/protein_identification_oracle_provenance.json.

use openms::chemistry::AASequence;
use openms::identification::*;
use openms::kernel::{ConsensusFeature, ConsensusMap, DataArray};
use openms::metadata::SourceFile;
use openms::{Error, MSExperiment};
use std::collections::{BTreeMap, BTreeSet};

const ORACLE: &str = include_str!("data/protein_identification_oracle.tsv");

/// The Release records of one oracle case, keyed by record key.
fn oracle(case: &str) -> BTreeMap<&'static str, &'static str> {
    let records: BTreeMap<_, _> = ORACLE
        .lines()
        .filter_map(|line| {
            let mut fields = line.splitn(3, '\t');
            let (name, key, value) = (fields.next()?, fields.next()?, fields.next()?);
            (name == case).then_some((key, value))
        })
        .collect();
    assert_eq!(records.get("exit"), Some(&"0"), "{case} ran to completion");
    records
}

fn flag(value: bool) -> &'static str {
    if value { "true" } else { "false" }
}

fn list(values: &[String]) -> String {
    format!("[{}]", values.join("|"))
}

#[test]
fn charge_ranges_match_every_release_record() {
    let mut compared = 0;
    for case in ["charge_range", "charge_range2"] {
        for (key, expected) in oracle(case) {
            let Some(text) = key
                .strip_prefix("charges<")
                .and_then(|rest| rest.strip_suffix('>'))
            else {
                continue;
            };
            let parameters = SearchParameters {
                charges: text.into(),
                ..Default::default()
            };
            let actual = match parameters.charge_range() {
                Ok((low, high)) => format!("{low},{high}"),
                Err(Error::InvalidValue(_)) => "throws ConversionError".into(),
                Err(Error::MissingInformation(_)) => "throws MissingInformation".into(),
                Err(other) => format!("unexpected {other}"),
            };
            assert_eq!(actual, expected, "charges {text:?}");
            compared += 1;
        }
    }
    assert_eq!(compared, 71);
}

#[test]
fn upstream_charge_range_goldens() {
    for (text, expected) in [
        ("1,2,3", (1, 3)),
        ("+2-+5", (2, 5)),
        ("-1,-2,-3", (-3, -1)),
        ("2", (2, 2)),
    ] {
        let parameters = SearchParameters {
            charges: text.into(),
            ..Default::default()
        };
        assert_eq!(parameters.charge_range().unwrap(), expected, "{text}");
    }
}

#[test]
fn unreadable_charges_neither_invalidate_nor_block_merging() {
    // test-data topp/SequestAdapter_4_output.idXML stores charges="+1+3"; the
    // source reads it as (0, 0) and never rejects the settings.
    let mut run = ProteinIdentification::new();
    run.search_parameters.charges = "+1+3".into();
    run.validate().unwrap();
    assert_eq!(run.search_parameters.charge_range().unwrap(), (0, 0));
    run.sort().unwrap();
    let records = oracle("mergeable_invalid_charges");
    let mut a = ProteinIdentification::new();
    a.search_engine = "Comet".into();
    let mut c = a.clone();
    c.search_parameters.charges = "invalid".into();
    assert_eq!(
        flag(
            a.search_parameters
                .mergeable(&c.search_parameters, "labeled_MS1")
                .unwrap()
        ),
        records["sp_mergeable"]
    );
    assert_eq!(
        flag(a.peptide_ids_mergeable(&c, "labeled_MS1").unwrap()),
        records["pi_mergeable"]
    );
    let d = c.clone();
    assert_eq!(
        flag(c.peptide_ids_mergeable(&d, "label-free").unwrap()),
        records["pi_mergeable_same_invalid"]
    );
}

#[cfg(feature = "idxml")]
#[test]
fn idxml_with_a_source_unreadable_charge_string_loads() {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/protein_identification_order_comet_first.idXML"
    ))
    .unwrap()
    .replace("charges=\"\"", "charges=\"+1+3\"");
    let document = openms::format::idxml::read(text.as_bytes()).unwrap();
    let run = &document.protein_identifications[0];
    assert_eq!(run.search_parameters.charges, "+1+3");
    assert_eq!(run.search_parameters.charge_range().unwrap(), (0, 0));
}

#[test]
fn settings_equality_leaves_out_settings_metadata_like_release() {
    let records = oracle("sp_equality");
    let a = SearchParameters::default();
    let mut c = SearchParameters::default();
    c.metadata.insert("SE:Comet".into(), "x".into());
    assert_eq!(flag(a.source_equal(&c)), records["sp_meta_only_eq"]);
    assert_eq!(flag(!a.source_equal(&c)), records["sp_meta_only_ne"]);
    assert_ne!(a, c, "the derived PartialEq stays structural");
    let d = SearchParameters {
        mass_type: PeakMassType::Average,
        ..Default::default()
    };
    assert_eq!(flag(a.source_equal(&d)), records["sp_mass_type_eq"]);
    let e = SearchParameters {
        missed_cleavages: 3,
        ..Default::default()
    };
    assert_eq!(flag(a.source_equal(&e)), records["sp_missed_cleavages_eq"]);

    let p1 = ProteinIdentification::new();
    let mut p2 = ProteinIdentification::new();
    p2.set_inference_engine("Epifany");
    assert_eq!(
        flag(p1.source_equal(&p2)),
        records["pi_search_meta_only_eq"]
    );
    assert_eq!(p1.inference_engine().unwrap(), records["p1_inference"]);
    assert_eq!(p2.inference_engine().unwrap(), records["p2_inference"]);
    let mut p3 = ProteinIdentification::new();
    p3.metadata.insert("x".into(), "1".into());
    assert_eq!(flag(p1.source_equal(&p3)), records["pi_run_meta_eq"]);
    let mut p4 = ProteinIdentification::new();
    p4.set_primary_ms_run_path(&["a.mzML".into()], false)
        .unwrap();
    assert_eq!(flag(p1.source_equal(&p4)), records["pi_run_path_eq"]);
    assert!(p1.source_equal(&p1.clone()));
}

#[test]
fn name_tables_match_release() {
    let records = oracle("names");
    assert_eq!(
        list(&PeakMassType::all_names()),
        records["peak_mass_type_names"]
    );
    assert_eq!(PeakMassType::Monoisotopic.name(), "Monoisotopic");
    assert_eq!(PeakMassType::Average.name(), "Average");
    // Only the four specificities the native enum represents; no-cterm (8) and
    // no-nterm (9) remain an open gap.
    for (specificity, index) in [
        (EnzymeTermSpecificity::None, 0),
        (EnzymeTermSpecificity::Semi, 1),
        (EnzymeTermSpecificity::Full, 2),
        (EnzymeTermSpecificity::Unknown, 3),
    ] {
        let mut run = ProteinIdentification::new();
        run.search_parameters.enzyme_specificity = specificity;
        let pairs = run.search_engine_settings_as_pairs("").unwrap();
        let exported = pairs
            .iter()
            .find(|(key, _)| key == "enzyme_term_specificity")
            .map(|(_, value)| value.as_str());
        assert_eq!(
            exported,
            Some(records[format!("settings_specificity{index}").as_str()])
        );
    }
}

#[test]
fn run_path_lists_match_release() {
    let records = oracle("run_paths");
    let mut run = ProteinIdentification::new();
    assert_eq!(
        run.nr_primary_ms_run_paths(false).to_string(),
        records["nr_initial"]
    );
    run.set_primary_ms_run_path(&[], false).unwrap();
    assert_eq!(
        list(run.primary_ms_run_path(false)),
        records["get_after_empty_set"]
    );
    run.set_primary_ms_run_path(&["a.mzML".into(), "b.raw".into()], false)
        .unwrap();
    assert_eq!(
        list(run.primary_ms_run_path(false)),
        records["get_after_set"]
    );
    run.add_primary_ms_run_path(&["c.mzML".into()], false)
        .unwrap();
    run.add_primary_ms_run_path(&["a.mzML".into()], false)
        .unwrap();
    assert_eq!(
        list(run.primary_ms_run_path(false)),
        records["get_after_add"]
    );
    assert_eq!(
        run.nr_primary_ms_run_paths(false).to_string(),
        records["nr_after_add"]
    );
    run.add_primary_ms_run_path(&["d".into(), "e".into()], true)
        .unwrap();
    assert_eq!(
        list(run.primary_ms_run_path(true)),
        records["get_raw_after_add"]
    );
    assert_eq!(
        run.nr_primary_ms_run_paths(true).to_string(),
        records["nr_raw"]
    );
    run.set_primary_ms_run_path(&["x.raw".into()], true)
        .unwrap();
    assert_eq!(
        list(run.primary_ms_run_path(true)),
        records["get_raw_after_set"]
    );
    assert_eq!(
        list(run.primary_ms_run_path(false)),
        records["get_after_raw_set"]
    );
}

/// An oracle label, the experiment's (path, name) source files and the fallback paths.
type Scenario<'a> = (&'a str, Vec<(String, String)>, Vec<String>);

#[test]
fn run_paths_from_an_experiment_match_release() {
    let records = oracle("run_paths_experiment");
    let dir = std::env::temp_dir().join(format!("openms-protein-run-paths-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for name in [
        "exists.mzML",
        "exists.mzml",
        "UPPER.MZML",
        "exists.mzML.gz",
        "exists.raw",
        "RUN.RAW",
    ] {
        std::fs::write(dir.join(name), "x").unwrap();
    }
    std::fs::create_dir_all(dir.join("folder.mzML")).unwrap();
    let base = dir.to_str().unwrap().to_string();
    let file_uri = format!("file://{base}");
    let one = |path: &str, name: &str| vec![(path.to_string(), name.to_string())];
    let fallback = vec!["fallback.mzML".to_string()];
    let scenarios: Vec<Scenario> = vec![
        ("existing_mzML", one(&base, "exists.mzML"), fallback.clone()),
        ("missing_mzML", one(&base, "missing.mzML"), fallback.clone()),
        (
            "existing_lower_mzml",
            one(&base, "exists.mzml"),
            fallback.clone(),
        ),
        (
            "existing_upper_MZML",
            one(&base, "UPPER.MZML"),
            fallback.clone(),
        ),
        (
            "existing_mzML_gz",
            one(&base, "exists.mzML.gz"),
            fallback.clone(),
        ),
        (
            "directory_named_mzML",
            one(&base, "folder.mzML"),
            fallback.clone(),
        ),
        ("existing_raw", one(&base, "exists.raw"), fallback.clone()),
        ("upper_RAW", one(&base, "RUN.RAW"), fallback.clone()),
        ("missing_raw", one(&base, "missing.raw"), fallback.clone()),
        (
            "two_files",
            vec![
                (base.clone(), "exists.mzML".into()),
                (base.clone(), "exists.mzML".into()),
            ],
            fallback.clone(),
        ),
        ("empty_name", one(&base, ""), fallback.clone()),
        ("file_uri", one(&file_uri, "exists.mzML"), fallback.clone()),
        ("no_files", Vec::new(), fallback.clone()),
        (
            "missing_mzML_empty_fallback",
            one(&base, "missing.mzML"),
            Vec::new(),
        ),
        ("raw_empty_fallback", one(&base, "exists.raw"), Vec::new()),
        ("other_type", one(&base, "exists.mzXML"), fallback.clone()),
    ];
    // The source keeps the paths in metadata; an absent key and an empty list
    // are one state in the dedicated native lists (a fresh run starts empty).
    let show = |paths: &[String], expected: &str| {
        if paths.is_empty() && expected == "absent" {
            "absent".to_string()
        } else {
            list(paths).replace(&base, "<dir>")
        }
    };
    for (label, files, paths) in &scenarios {
        let mut experiment = MSExperiment::default();
        for (path, name) in files {
            experiment.settings.source_files.push(SourceFile {
                path: path.clone(),
                name: name.clone(),
                ..Default::default()
            });
        }
        let mut run = ProteinIdentification::new();
        run.set_primary_ms_run_path_from_experiment(paths, &experiment)
            .unwrap();
        let mzml = records[format!("{label}.spectra_data").as_str()];
        let raw = records[format!("{label}.spectra_data_raw").as_str()];
        assert_eq!(show(run.primary_ms_run_path(false), mzml), mzml, "{label}");
        assert_eq!(show(run.primary_ms_run_path(true), raw), raw, "{label}");
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn run_path_name_rule_agrees_with_the_file_type_table() {
    // The run-path overload restricts FileNameUtils::getTypeByFileName to mzML
    // and raw; format::file_types::type_by_file_name is the full port of it.
    use openms::format::FileType;
    use openms::format::file_types::type_by_file_name;
    let dir = std::env::temp_dir().join(format!("openms-protein-run-kinds-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let base = dir.to_str().unwrap().to_string();
    let names = [
        "a.mzML",
        "a.MZML",
        "a.mzml.GZ",
        "a.mzML.bz2",
        "a.mzML.zip",
        "a.mzML.gz.gz",
        "a.raw",
        "a.RAW.gz",
        "a.pep.xml",
        "a.mzXML",
        "mzML",
        "a.",
        ".mzML",
        "a.mzML.tar",
        "a.gz",
        "fid",
        "a.idXML",
        "a.raw.mzML",
        "a.mzML.raw",
    ];
    for name in names {
        std::fs::write(dir.join(name), "x").unwrap();
        let mut experiment = MSExperiment::default();
        experiment.settings.source_files.push(SourceFile {
            path: base.clone(),
            name: name.into(),
            ..Default::default()
        });
        let mut run = ProteinIdentification::new();
        run.set_primary_ms_run_path_from_experiment(&["fb".into()], &experiment)
            .unwrap();
        let full = format!("{base}/{name}");
        let (mzml, raw): (Vec<String>, Vec<String>) = match type_by_file_name(&full) {
            FileType::MzMl => (vec![full.clone()], Vec::new()),
            FileType::Raw => (vec!["fb".into()], vec![full.clone()]),
            _ => (vec!["fb".into()], Vec::new()),
        };
        assert_eq!(run.primary_ms_run_path(false), mzml, "{name}");
        assert_eq!(run.primary_ms_run_path(true), raw, "{name}");
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn an_existing_mzml_leaves_the_raw_paths_alone() {
    // Source: the mzML branch returns before anything else is set.
    let dir = std::env::temp_dir().join(format!("openms-protein-raw-kept-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("run.mzML"), "x").unwrap();
    let mut experiment = MSExperiment::default();
    experiment.settings.source_files.push(SourceFile {
        path: dir.to_str().unwrap().into(),
        name: "run.mzML".into(),
        ..Default::default()
    });
    let mut run = ProteinIdentification::new();
    run.raw_ms_run_paths = vec!["kept.raw".into()];
    run.set_primary_ms_run_path_from_experiment(&["ignored.mzML".into()], &experiment)
        .unwrap();
    assert_eq!(run.raw_ms_run_paths, ["kept.raw"]);
    assert_eq!(run.primary_ms_run_path(false).len(), 1);
    std::fs::remove_dir_all(&dir).unwrap();
}

fn protein(accession: &str, sequence: &str) -> ProteinHit {
    ProteinHit::new(1.0, 0, accession, sequence).unwrap()
}

fn evidence(accession: &str, start: Option<usize>, end: Option<usize>) -> PeptideEvidence {
    PeptideEvidence {
        protein_accession: accession.into(),
        start,
        end,
        ..Default::default()
    }
}

fn identifications(sequence: &str, evidences: Vec<PeptideEvidence>) -> Vec<PeptideIdentification> {
    let mut hit = PeptideHit::new(0.0, 0, 1, AASequence::parse(sequence).unwrap()).unwrap();
    hit.evidences = evidences;
    vec![PeptideIdentification {
        hits: vec![hit],
        ..Default::default()
    }]
}

fn coverage_record(run: &ProteinIdentification) -> String {
    let values: Vec<String> = run
        .hits
        .iter()
        .map(|hit| format!("{}:{}", hit.accession, hit.coverage.unwrap()))
        .collect();
    list(&values)
}

#[test]
fn coverage_reads_only_what_it_uses_like_release() {
    let records = oracle("coverage");
    let p1 = || protein("P1", "ABCDEFGHIJ");
    let cases: Vec<(&str, Vec<ProteinHit>, Vec<PeptideEvidence>)> = vec![
        (
            "unrelated_unknown_positions",
            vec![p1()],
            vec![
                evidence("P1", Some(2), Some(4)),
                evidence("OTHER", None, None),
            ],
        ),
        (
            "unrelated_reversed_positions",
            vec![p1()],
            vec![
                evidence("P1", Some(2), Some(4)),
                evidence("OTHER", Some(5), Some(3)),
            ],
        ),
        (
            "own_unknown_positions",
            vec![p1()],
            vec![evidence("P1", None, None)],
        ),
        (
            "own_reversed_positions",
            vec![p1()],
            vec![evidence("P1", Some(5), Some(3))],
        ),
        (
            "empty_sequence_without_evidence",
            vec![p1(), protein("P2", "")],
            vec![evidence("P1", Some(2), Some(4))],
        ),
        (
            "last_residue",
            vec![p1()],
            vec![evidence("P1", Some(9), Some(9))],
        ),
        (
            "nan_protein_score",
            vec![p1()],
            vec![evidence("P1", Some(2), Some(4))],
        ),
    ];
    for (label, hits, evidences) in cases {
        let mut run = ProteinIdentification {
            hits,
            ..Default::default()
        };
        if label == "nan_protein_score" {
            run.hits[0].score = f64::NAN;
        }
        let before = run.clone();
        let actual = match run.compute_coverage(&identifications("ABC", evidences)) {
            Ok(()) => coverage_record(&run),
            Err(Error::MissingInformation(_)) => {
                assert_eq!(
                    format!("{run:?}"),
                    format!("{before:?}"),
                    "{label} is atomic"
                );
                "throws MissingInformation".into()
            }
            Err(other) => format!("unexpected {other}"),
        };
        assert_eq!(actual, records[label], "{label}");
    }
}

#[test]
fn modifications_and_coverage_ignore_scores_and_foreign_evidence_like_release() {
    let records = oracle("validation_free_helpers");
    let mut run = ProteinIdentification {
        hits: vec![protein("P1", "ABPEPMKXX")],
        ..Default::default()
    };
    run.hits[0].score = f64::NAN;
    let mut ids = identifications(
        "PEPM(Oxidation)K",
        vec![
            evidence("P1", Some(2), Some(6)),
            evidence("OTHER", Some(5), Some(3)),
        ],
    );
    ids[0].hits[0].score = f64::NAN;
    run.compute_modifications(&ids, &BTreeSet::new()).unwrap();
    let mods: Vec<String> = run.hits[0]
        .modifications
        .iter()
        .map(|m| format!("{}:{}", m.position, m.modification.full_id()))
        .collect();
    assert_eq!(list(&mods), records["modifications"]);
    run.compute_coverage(&ids).unwrap();
    let expected: f64 = records["coverage"].parse().unwrap();
    assert_eq!(run.hits[0].coverage, Some(expected));
}

#[test]
fn consensus_map_overloads_visit_features_then_unassigned() {
    // Source computeCoverage/computeModifications(const ConsensusMap&, ...):
    // every feature's identifications, then the unassigned ones on request.
    let run = ProteinIdentification {
        hits: vec![protein("P1", "ABPEPMKXXAB")],
        ..Default::default()
    };
    let assigned = identifications("PEPM(Oxidation)K", vec![evidence("P1", Some(2), Some(6))]);
    let unassigned = identifications("A", vec![evidence("P1", Some(9), Some(10))]);
    let mut feature = ConsensusFeature::default();
    feature.base.peptide_identifications = assigned.clone();
    let map = ConsensusMap {
        features: vec![feature],
        unassigned_peptide_identifications: unassigned.clone(),
        ..Default::default()
    };
    let all: Vec<PeptideIdentification> = assigned.iter().chain(&unassigned).cloned().collect();
    for (use_unassigned, ids) in [(false, &assigned), (true, &all)] {
        let mut from_map = run.clone();
        from_map
            .compute_coverage_from_consensus_map(&map, use_unassigned)
            .unwrap();
        let mut from_list = run.clone();
        from_list.compute_coverage(ids).unwrap();
        assert_eq!(from_map.hits, from_list.hits);
        let mut from_map = run.clone();
        from_map
            .compute_modifications_from_consensus_map(&map, &BTreeSet::new(), use_unassigned)
            .unwrap();
        let mut from_list = run.clone();
        from_list
            .compute_modifications(ids, &BTreeSet::new())
            .unwrap();
        assert_eq!(from_map.hits, from_list.hits);
    }
    let mut covered = run.clone();
    covered
        .compute_coverage_from_consensus_map(&map, true)
        .unwrap();
    assert_eq!(covered.hits[0].coverage, Some(100.0 * 7.0 / 11.0));
}

#[test]
fn original_engine_follows_release_where_key_orders_coincide() {
    // Release takes the first SE: key (without "percolator") in registry order,
    // which the recorded key lists show. Where that key is also the lexically
    // first one the native answer must match; the other cases are the
    // documented registry-order divergence and are only checked on the
    // Release side here.
    let mut agreeing = 0;
    for case in [
        "idxml_comet_first",
        "idxml_msgf_first",
        "engine_fresh",
        "engine_global",
        "engine_consensus",
    ] {
        let records = oracle(case);
        let keys: Vec<&str> = records["keys"]
            .trim_start_matches('[')
            .trim_end_matches(']')
            .split('|')
            .collect();
        let candidates: Vec<&str> = keys
            .iter()
            .copied()
            .filter(|key| key.starts_with("SE:") && !key.contains("percolator"))
            .collect();
        let registry_first = candidates[0];
        assert_eq!(
            &registry_first[3..],
            records["original"],
            "{case}: the source rule"
        );
        let lexical_first = *candidates.iter().min().unwrap();
        if registry_first != lexical_first {
            continue;
        }
        let mut run = ProteinIdentification::new();
        run.search_engine = if case == "engine_consensus" {
            "ConsensusID".into()
        } else {
            "Percolator".into()
        };
        for key in &keys {
            run.search_parameters
                .metadata
                .insert((*key).into(), "v".into());
        }
        assert_eq!(
            run.original_search_engine_name(),
            records["original"],
            "{case}"
        );
        agreeing += 1;
    }
    assert_eq!(agreeing, 1);
}

#[cfg(feature = "idxml")]
#[test]
fn idxml_run_names_the_release_engine_when_orders_coincide() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/protein_identification_order_comet_first.idXML"
    );
    let document = openms::format::idxml::load(path).unwrap();
    let run = &document.protein_identifications[0];
    assert_eq!(
        run.original_search_engine_name(),
        oracle("idxml_comet_first")["original"]
    );
}

#[test]
fn upstream_protein_group_order_and_arrays() {
    // ProteinIdentification_test.cpp: ProteinGroup operator< and get*DataArrayByName.
    let mut p1 = ProteinGroup::default();
    let mut p2 = ProteinGroup::default();
    assert!(!p1.source_less(&p2) && !p2.source_less(&p1));
    p1.probability = 0.1;
    p2.probability = 0.2;
    assert!(p2.source_less(&p1));
    assert!(!p1.source_less(&p2));
    p2.probability = 0.1;
    p1.accessions.push("bla".into());
    p2.accessions.push("bla".into());
    assert!(!p1.source_less(&p2) && !p2.source_less(&p1));
    p2.accessions.push("blubb".into());
    assert!(p1.source_less(&p2));
    assert!(!p2.source_less(&p1));
    p1.accessions.push("laber".into());
    assert!(!p1.source_less(&p2));
    assert!(p2.source_less(&p1));

    let mut group = ProteinGroup::default();
    group
        .integer_data_arrays
        .push(DataArray::new("int_da", Vec::new()));
    group
        .string_data_arrays
        .push(DataArray::new("str_da", Vec::new()));
    group
        .float_data_arrays
        .push(DataArray::new("float_da", Vec::new()));
    assert_eq!(
        group.integer_data_array_by_name("int_da").unwrap().name,
        "int_da"
    );
    assert_eq!(
        group.string_data_array_by_name("str_da").unwrap().name,
        "str_da"
    );
    assert_eq!(
        group.float_data_array_by_name("float_da").unwrap().name,
        "float_da"
    );
    assert!(group.integer_data_array_by_name("missing_int").is_err());
    assert!(group.string_data_array_by_name("missing_str").is_err());
    assert!(group.float_data_array_by_name("missing_float").is_err());
    assert!(group.integer_data_array_by_name_mut("missing_int").is_err());
    assert!(group.string_data_array_by_name_mut("missing_str").is_err());
    assert!(group.float_data_array_by_name_mut("missing_float").is_err());
    group
        .float_data_array_by_name_mut("float_da")
        .unwrap()
        .data
        .push(2.5);
    group
        .integer_data_array_by_name_mut("int_da")
        .unwrap()
        .data
        .push(3);
    group
        .string_data_array_by_name_mut("str_da")
        .unwrap()
        .data
        .push("s".into());
    assert_eq!(group.float_data_arrays[0].data, [2.5]);
    assert_eq!(group.integer_data_arrays[0].data, [3]);
    assert_eq!(group.string_data_arrays[0].data, ["s"]);
}

#[test]
fn upstream_find_hit_and_peak_mass_type_names() {
    // ProteinIdentification_test.cpp: findHit and getAllNamesOfPeakMassType.
    let mut run = ProteinIdentification::new();
    run.hits.push(protein("test1", ""));
    run.hits.push(protein("test2", ""));
    assert_eq!(run.find_hit_mut("test1").unwrap().accession, "test1");
    assert_eq!(run.find_hit_mut("test2").unwrap().accession, "test2");
    assert!(run.find_hit_mut("test3").is_none());
    run.find_hit_mut("test2").unwrap().score = 5.0;
    assert_eq!(run.find_hit("test2").unwrap().score, 5.0);
    assert_eq!(PeakMassType::all_names(), ["Monoisotopic", "Average"]);
    assert_eq!(PeakMassType::NAMES.len(), 2);
}

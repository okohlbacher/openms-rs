// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `TransitionParquetFile` (`.oswpq` assay libraries).
//!
//! Two sources of expected values:
//!
//! - `TransitionParquetFile_test.cpp` at core `bc9cc12`, both sections
//!   transcribed (the class test's `convertTargetedExp` of
//!   `MRMAssay_detectingTransistionCompound_input.TraML` is replaced by the
//!   light experiment the Release build produced from it, record
//!   `traml_input`, because `OpenSwathDataAccessHelper::convertTargetedExp`
//!   is not ported).
//! - The executed Release oracle `../oracle/O5c-parquet` (ibminode06,
//!   `openms4-release-bc9cc12-c19e494-174b576`, Arrow C++ 23.0.1, run twice
//!   with identical records and files): `tests/data/transition_parquet_file/`
//!   holds its record stream and every table and metadata file it wrote.
//!   The reader is replayed on the C++-written containers and compared field
//!   by field with the C++ reader's dump; the writer is run on the same
//!   inputs and its tables compared value by value with the C++ writer's.
//!   See `tests/data/transition_parquet_file_provenance.json`.
#![cfg(all(feature = "arrow", feature = "file-compression"))]

use openms::Error;
use openms::analysis::openswath::openswath_library_id_normalizer::{
    OpenSwathLibraryIDNormalizer, SourceIDMapping,
};
use openms::analysis::openswath::transition_experiment::{
    LightCompound, LightModification, LightProtein, LightTargetedExperiment, LightTransition,
};
use openms::analysis::openswath::transition_parquet_file::TransitionParquetFile;
use openms::chemistry::AASequence;
use openms::format::arrow_schema_registry::{OSWPrecursorSchema, OSWTransitionSchema};
use openms::format::parquet_file::ParquetFile;
use openms::format::zip_archive_file::ZipArchiveFile;
use openms::system::file::TempDir;

use arrow_array::{
    ArrayRef, BooleanArray, Float64Array, Int32Array, Int64Array, RecordBatch, StringArray,
};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::sync::Arc;

const DATA: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/transition_parquet_file"
);

// ------------------------------------------------------------------ records

fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c == '\\' {
            match chars.next() {
                Some('t') => out.push('\t'),
                Some('n') => out.push('\n'),
                Some('\\') => out.push('\\'),
                Some(other) => {
                    out.push('\\');
                    out.push(other);
                }
                None => out.push('\\'),
            }
        } else {
            out.push(c);
        }
    }
    out
}

fn records() -> Vec<Vec<String>> {
    let text = std::fs::read_to_string(Path::new(DATA).join("records.tsv")).unwrap();
    text.lines()
        .map(|line| line.split('\t').map(unescape).collect())
        .collect()
}

fn of_case<'a>(all: &'a [Vec<String>], kind: &str, case: &str) -> Vec<&'a Vec<String>> {
    all.iter()
        .filter(|r| r[0] == kind && r.get(1).map(String::as_str) == Some(case))
        .collect()
}

fn outcome<'a>(all: &'a [Vec<String>], case: &str) -> &'a Vec<String> {
    all.iter()
        .find(|r| (r[0] == "OK" || r[0] == "E") && r[1] == case)
        .unwrap_or_else(|| panic!("no outcome record for {case}"))
}

fn num(s: &str) -> f64 {
    match s {
        "nan" | "-nan" => f64::NAN,
        "inf" => f64::INFINITY,
        "-inf" => f64::NEG_INFINITY,
        _ => s.parse().unwrap_or_else(|_| panic!("bad number {s}")),
    }
}

fn same(a: f64, b: f64) -> bool {
    (a.is_nan() && b.is_nan()) || a.to_bits() == b.to_bits()
}

fn split_list(s: &str) -> Vec<String> {
    if s.is_empty() {
        Vec::new()
    } else {
        s.split('|').map(str::to_owned).collect()
    }
}

fn mods(s: &str) -> Vec<LightModification> {
    s.split(',')
        .filter(|m| !m.is_empty())
        .map(|m| {
            let (location, unimod) = m.split_once(':').unwrap();
            LightModification {
                location: location.parse().unwrap(),
                unimod_id: unimod.parse().unwrap(),
            }
        })
        .collect()
}

fn compound_of(r: &[String]) -> LightCompound {
    LightCompound {
        id: r[3].clone(),
        rt: num(&r[4]),
        drift_time: num(&r[5]),
        charge: r[6].parse().unwrap(),
        sequence: r[7].clone(),
        protein_refs: split_list(&r[8]),
        compound_name: r[9].clone(),
        modifications: mods(&r[10]),
        rt_start: num(&r[11]),
        rt_end: num(&r[12]),
        peptide_group_label: r[13].clone(),
        gene_name: r[14].clone(),
        sum_formula: r[15].clone(),
        label_type: r[16].clone(),
        smiles: r[17].clone(),
        adducts: r[18].clone(),
    }
}

fn transition_of(r: &[String]) -> LightTransition {
    let mut t = LightTransition {
        transition_name: r[3].clone(),
        peptide_ref: r[4].clone(),
        library_intensity: num(&r[5]),
        product_mz: num(&r[6]),
        precursor_mz: num(&r[7]),
        precursor_im: num(&r[8]),
        fragment_charge: r[9].parse().unwrap(),
        fragment_nr: r[14].parse().unwrap(),
        ..LightTransition::default()
    };
    t.set_decoy(r[10] == "1");
    t.set_detecting_transition(r[11] == "1");
    t.set_identifying_transition(r[12] == "1");
    t.set_quantifying_transition(r[13] == "1");
    t.set_fragment_type(&r[15]);
    t
}

fn protein_of(r: &[String]) -> LightProtein {
    LightProtein {
        id: r[3].clone(),
        sequence: r[4].clone(),
        uniprot_id: r[5].clone(),
    }
}

/// The light experiment a record dump describes, in dump order.
fn experiment(all: &[Vec<String>], case: &str) -> LightTargetedExperiment {
    LightTargetedExperiment {
        compounds: of_case(all, "C", case)
            .into_iter()
            .map(|r| compound_of(r))
            .collect(),
        transitions: of_case(all, "T", case)
            .into_iter()
            .map(|r| transition_of(r))
            .collect(),
        proteins: of_case(all, "P", case)
            .into_iter()
            .map(|r| protein_of(r))
            .collect(),
    }
}

fn id_maps(all: &[Vec<String>], case: &str) -> BTreeMap<String, BTreeMap<String, String>> {
    let mut maps: BTreeMap<String, BTreeMap<String, String>> = BTreeMap::new();
    for name in [
        "precursor_source_to_canonical",
        "precursor_canonical_to_source",
        "transition_canonical_to_source",
    ] {
        maps.insert(name.to_owned(), BTreeMap::new());
    }
    for r in of_case(all, "M", case) {
        maps.get_mut(&r[2])
            .unwrap()
            .insert(r[3].clone(), r[4].clone());
    }
    maps
}

fn maps_of(ids: &SourceIDMapping) -> BTreeMap<String, BTreeMap<String, String>> {
    let sorted = |m: &HashMap<String, String>| {
        m.iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect::<BTreeMap<_, _>>()
    };
    BTreeMap::from([
        (
            "precursor_source_to_canonical".to_owned(),
            sorted(&ids.precursor_source_to_canonical),
        ),
        (
            "precursor_canonical_to_source".to_owned(),
            sorted(&ids.precursor_canonical_to_source),
        ),
        (
            "transition_canonical_to_source".to_owned(),
            sorted(&ids.transition_canonical_to_source),
        ),
    ])
}

// ------------------------------------------------------------------ comparisons

fn assert_compound_eq(case: &str, got: &LightCompound, want: &LightCompound) {
    let ctx = format!("{case} compound {}", want.id);
    assert_eq!(got.id, want.id, "{ctx}");
    assert!(same(got.rt, want.rt), "{ctx} rt {} {}", got.rt, want.rt);
    assert!(same(got.drift_time, want.drift_time), "{ctx} drift");
    assert_eq!(got.charge, want.charge, "{ctx} charge");
    assert_eq!(got.sequence, want.sequence, "{ctx} sequence");
    assert_eq!(got.protein_refs, want.protein_refs, "{ctx} protein_refs");
    assert_eq!(got.compound_name, want.compound_name, "{ctx} compound_name");
    assert_eq!(got.modifications, want.modifications, "{ctx} modifications");
    assert!(same(got.rt_start, want.rt_start), "{ctx} rt_start");
    assert!(same(got.rt_end, want.rt_end), "{ctx} rt_end");
    assert_eq!(got.peptide_group_label, want.peptide_group_label, "{ctx}");
    assert_eq!(got.gene_name, want.gene_name, "{ctx}");
    assert_eq!(got.sum_formula, want.sum_formula, "{ctx}");
    assert_eq!(got.label_type, want.label_type, "{ctx}");
    assert_eq!(got.smiles, want.smiles, "{ctx}");
    assert_eq!(got.adducts, want.adducts, "{ctx}");
}

fn assert_transition_eq(case: &str, got: &LightTransition, want: &LightTransition) {
    let ctx = format!("{case} transition {}", want.transition_name);
    assert_eq!(got.transition_name, want.transition_name, "{ctx}");
    assert_eq!(got.peptide_ref, want.peptide_ref, "{ctx}");
    assert!(same(got.library_intensity, want.library_intensity), "{ctx}");
    assert!(same(got.product_mz, want.product_mz), "{ctx}");
    assert!(same(got.precursor_mz, want.precursor_mz), "{ctx}");
    assert!(same(got.precursor_im, want.precursor_im), "{ctx}");
    assert_eq!(got.fragment_charge, want.fragment_charge, "{ctx}");
    assert_eq!(got.flags, want.flags, "{ctx}");
    assert_eq!(got.fragment_nr, want.fragment_nr, "{ctx}");
    assert_eq!(got.fragment_type, want.fragment_type, "{ctx}");
    assert_eq!(got.annotation(), want.annotation(), "{ctx}");
    assert_eq!(got.peptidoforms, want.peptidoforms, "{ctx}");
}

/// Compare a Rust-read experiment with the C++ reader's dump. Transitions
/// must agree in order; compounds and proteins are compared by identity,
/// since the source emits them in `std::unordered_map` order (see the
/// module notes), and their Rust order is checked against the precursor
/// table's row order separately.
fn assert_experiment_matches(case: &str, got: &LightTargetedExperiment, all: &[Vec<String>]) {
    let want = experiment(all, case);
    assert_eq!(
        got.compounds.len(),
        want.compounds.len(),
        "{case} compounds"
    );
    let by_id: HashMap<&str, &LightCompound> =
        got.compounds.iter().map(|c| (c.id.as_str(), c)).collect();
    for w in &want.compounds {
        assert_compound_eq(case, by_id[w.id.as_str()], w);
    }
    assert_eq!(
        got.transitions.len(),
        want.transitions.len(),
        "{case} transitions"
    );
    for (g, w) in got.transitions.iter().zip(&want.transitions) {
        assert_transition_eq(case, g, w);
    }
    let got_p: BTreeSet<(String, String, String)> = got
        .proteins
        .iter()
        .map(|p| (p.id.clone(), p.sequence.clone(), p.uniprot_id.clone()))
        .collect();
    let want_p: BTreeSet<(String, String, String)> = want
        .proteins
        .iter()
        .map(|p| (p.id.clone(), p.sequence.clone(), p.uniprot_id.clone()))
        .collect();
    assert_eq!(got.proteins.len(), want.proteins.len(), "{case} proteins");
    assert_eq!(got_p, want_p, "{case} proteins");
}

/// The error a C++ exception record corresponds to.
fn assert_error_matches(case: &str, err: &Error, record: &[String]) {
    let kind = record[2].as_str();
    let what = record[3].as_str();
    // `Exception::InvalidValue::what()` is "the value '<v>' was used but is
    // not valid; <message>"; the port carries message and value in one string.
    let (value, message) = match what.strip_prefix("the value '") {
        Some(rest) => {
            let (value, message) = rest.split_once("' was used but is not valid; ").unwrap();
            (value, message)
        }
        None => ("", what),
    };
    match (kind, err) {
        ("OpenMS::Exception::InvalidValue", Error::InvalidValue(m)) => {
            assert!(m.contains(message), "{case}: '{m}' lacks '{message}'");
            assert!(m.contains(value), "{case}: '{m}' lacks '{value}'");
        }
        ("OpenMS::Exception::MissingInformation", Error::MissingInformation(m)) => {
            assert_eq!(m, message, "{case}");
        }
        ("OpenMS::Exception::FileNotFound", Error::Io(e)) => {
            assert_eq!(e.kind(), std::io::ErrorKind::NotFound, "{case}: {e}");
        }
        _ => panic!("{case}: C++ threw {kind} ({what}), Rust returned {err:?}"),
    }
}

fn assert_tables_equal(case: &str, got: &Path, want: &Path) {
    let g = ParquetFile::read_table(got).unwrap();
    let w = ParquetFile::read_table(want).unwrap();
    assert_eq!(g.num_rows(), w.num_rows(), "{case} {}", want.display());
    let gs = g.schema();
    let ws = w.schema();
    let names = |s: &arrow_schema::Schema| {
        s.fields()
            .iter()
            .map(|f| (f.name().clone(), f.data_type().clone(), f.is_nullable()))
            .collect::<Vec<_>>()
    };
    assert_eq!(names(&gs), names(&ws), "{case} schema {}", want.display());
    for (i, field) in ws.fields().iter().enumerate() {
        assert_eq!(
            g.column(i).to_data(),
            w.column(i).to_data(),
            "{case} column {} of {}",
            field.name(),
            want.display()
        );
    }
}

fn assert_metadata_equal(case: &str, got: &Path, want: &Path) {
    let mask = |text: String| {
        text.lines()
            .map(|line| {
                if line.trim_start().starts_with("\"build_time\"") {
                    "    \"build_time\": <masked>,".to_owned()
                } else {
                    line.to_owned()
                }
            })
            .collect::<Vec<_>>()
            .join("\n")
    };
    let g = std::fs::read_to_string(got).unwrap();
    let w = std::fs::read_to_string(want).unwrap();
    assert!(g.contains("    \"build_time\": \"\",\n"), "{case}: {g}");
    assert!(g.ends_with("}\n") && w.ends_with("}\n"));
    assert_eq!(mask(g), mask(w), "{case} metadata.json");
}

fn assert_library_matches(case: &str, dir: &Path) {
    let files = Path::new(DATA).join("files");
    for name in ["precursors.parquet", "transitions.parquet"] {
        assert_tables_equal(
            case,
            &dir.join("library").join(name),
            &files.join(format!("{case}__{name}")),
        );
    }
    assert_metadata_equal(
        case,
        &dir.join("library/metadata.json"),
        &files.join(format!("{case}__metadata.json")),
    );
}

/// A `.oswpq` directory assembled from the oracle's flattened files of `source`.
fn container(tmp: &Path, source: &str, name: &str) -> PathBuf {
    let dir = tmp.join(format!("{name}.oswpq"));
    std::fs::create_dir_all(dir.join("library")).unwrap();
    let files = Path::new(DATA).join("files");
    for entry in ["precursors.parquet", "transitions.parquet", "metadata.json"] {
        let from = files.join(format!("{source}__{entry}"));
        if from.exists() {
            std::fs::copy(&from, dir.join("library").join(entry)).unwrap();
        }
    }
    dir
}

fn tmp() -> TempDir {
    TempDir::new_in(std::env::temp_dir(), false).unwrap()
}

// ------------------------------------------------------------------ class test

/// `TransitionParquetFile_test.cpp`, section
/// `convertParquetToTargetedExperiment`.
#[test]
fn class_test_convert_parquet_to_targeted_experiment() {
    let all = records();
    let light = experiment(&all, "traml_input");
    let compound_count = light.compounds.len().min(2);
    assert_eq!(compound_count, 2);

    let mut compound_to_precursor: BTreeMap<String, i64> = BTreeMap::new();
    for (i, compound) in light.compounds.iter().take(compound_count).enumerate() {
        compound_to_precursor.insert(compound.id.clone(), i as i64 + 1);
    }
    let mut transitions: Vec<LightTransition> = Vec::new();
    let mut precursor_mz: HashMap<String, f64> = HashMap::new();
    for t in &light.transitions {
        if compound_to_precursor.contains_key(&t.peptide_ref) {
            transitions.push(t.clone());
            precursor_mz
                .entry(t.peptide_ref.clone())
                .or_insert(t.precursor_mz);
        }
    }
    assert!(!transitions.is_empty());

    let tmp = tmp();
    let base_dir = tmp.path().join("test.oswpq");
    std::fs::create_dir_all(base_dir.join("library")).unwrap();

    let compounds = &light.compounds[..compound_count];
    let unmodified = |s: &str| {
        if s.is_empty() {
            String::new()
        } else {
            AASequence::parse(s)
                .map(|a| a.as_str().to_owned())
                .unwrap_or_default()
        }
    };
    let precursor_columns = |ids: Vec<i64>| -> Vec<ArrayRef> {
        vec![
            Arc::new(Int64Array::from(ids)),
            Arc::new(Float64Array::from_iter_values(
                compounds.iter().map(|c| precursor_mz[&c.id]),
            )),
            Arc::new(Int32Array::from_iter_values(
                compounds.iter().map(|c| c.charge),
            )),
            Arc::new(Float64Array::from_iter_values(
                compounds.iter().map(|c| c.rt),
            )),
            Arc::new(Float64Array::from_iter_values(
                compounds.iter().map(|c| c.drift_time),
            )),
            Arc::new(BooleanArray::from(vec![false; compound_count])),
            Arc::new(StringArray::from_iter_values(
                compounds.iter().map(|c| c.id.clone()),
            )),
            Arc::new(StringArray::from_iter_values(
                compounds.iter().map(|c| c.sequence.clone()),
            )),
            Arc::new(StringArray::from_iter_values(
                compounds.iter().map(|c| unmodified(&c.sequence)),
            )),
            Arc::new(StringArray::from_iter_values(
                compounds.iter().map(|c| c.protein_refs.join(";")),
            )),
        ]
    };
    let transition_columns = |precursor_refs: Vec<i64>| -> Vec<ArrayRef> {
        vec![
            Arc::new(Int64Array::from_iter_values(1..=transitions.len() as i64)),
            Arc::new(Int64Array::from(precursor_refs)),
            Arc::new(StringArray::from_iter_values(
                transitions.iter().map(|t| t.transition_name.clone()),
            )),
            Arc::new(Float64Array::from_iter_values(
                transitions.iter().map(|t| t.product_mz),
            )),
            Arc::new(Int32Array::from_iter_values(
                transitions.iter().map(|t| i32::from(t.fragment_charge)),
            )),
            Arc::new(StringArray::from_iter_values(
                transitions.iter().map(|t| t.fragment_type_str()),
            )),
            Arc::new(StringArray::from_iter_values(
                transitions.iter().map(|t| t.annotation()),
            )),
            Arc::new(Int32Array::from_iter_values(
                transitions.iter().map(|t| i32::from(t.fragment_nr)),
            )),
            Arc::new(BooleanArray::from_iter(
                transitions
                    .iter()
                    .map(|t| Some(t.is_detecting_transition())),
            )),
            Arc::new(BooleanArray::from_iter(
                transitions
                    .iter()
                    .map(|t| Some(t.is_identifying_transition())),
            )),
            Arc::new(BooleanArray::from_iter(
                transitions
                    .iter()
                    .map(|t| Some(t.is_quantifying_transition())),
            )),
            Arc::new(Float64Array::from_iter_values(
                transitions.iter().map(|t| t.library_intensity),
            )),
            Arc::new(BooleanArray::from_iter(
                transitions.iter().map(|t| Some(t.decoy())),
            )),
        ]
    };
    let ids: Vec<i64> = compounds
        .iter()
        .map(|c| compound_to_precursor[&c.id])
        .collect();
    let refs: Vec<i64> = transitions
        .iter()
        .map(|t| compound_to_precursor[&t.peptide_ref])
        .collect();
    let precursors =
        RecordBatch::try_new(OSWPrecursorSchema::schema(), precursor_columns(ids)).unwrap();
    ParquetFile::write_table(&precursors, base_dir.join("library/precursors.parquet")).unwrap();
    let trans =
        RecordBatch::try_new(OSWTransitionSchema::schema(), transition_columns(refs)).unwrap();
    ParquetFile::write_table(&trans, base_dir.join("library/transitions.parquet")).unwrap();

    let reader = TransitionParquetFile::new();
    let mut out = LightTargetedExperiment::default();
    let mut source_ids = SourceIDMapping::default();
    reader
        .convert_parquet_to_targeted_experiment_with_source_ids(
            &base_dir,
            &mut out,
            Some(&mut source_ids),
        )
        .unwrap();
    assert_eq!(out.compounds.len(), compound_count);
    assert_eq!(out.transitions.len(), transitions.len());

    let expected_precursor_ids: BTreeSet<String> = compound_to_precursor
        .values()
        .map(|v| v.to_string())
        .collect();
    for (source, canonical) in &compound_to_precursor {
        let canonical = canonical.to_string();
        assert_eq!(source_ids.precursor_source_to_canonical[source], canonical);
        assert_eq!(
            source_ids.precursor_canonical_to_source[&canonical],
            *source
        );
    }
    let mut compound_refs = BTreeSet::new();
    for compound in &out.compounds {
        assert!(expected_precursor_ids.contains(&compound.id));
        compound_refs.insert(compound.id.clone());
    }
    assert_eq!(compound_refs.len(), expected_precursor_ids.len());
    for (i, t) in out.transitions.iter().enumerate() {
        assert_eq!(t.transition_name, (i + 1).to_string());
        assert_eq!(
            t.peptide_ref,
            compound_to_precursor[&transitions[i].peptide_ref].to_string()
        );
        assert!(compound_refs.contains(&t.peptide_ref));
    }

    // Duplicate persistent precursor ids are rejected at the reader boundary.
    let duplicate_dir = tmp.path().join("duplicate.oswpq");
    std::fs::create_dir_all(duplicate_dir.join("library")).unwrap();
    let dup_precursors =
        RecordBatch::try_new(OSWPrecursorSchema::schema(), precursor_columns(vec![1; 2])).unwrap();
    ParquetFile::write_table(
        &dup_precursors,
        duplicate_dir.join("library/precursors.parquet"),
    )
    .unwrap();
    let dup_trans = RecordBatch::try_new(
        OSWTransitionSchema::schema(),
        transition_columns(vec![1; transitions.len()]),
    )
    .unwrap();
    ParquetFile::write_table(
        &dup_trans,
        duplicate_dir.join("library/transitions.parquet"),
    )
    .unwrap();
    let mut duplicate_out = LightTargetedExperiment::default();
    let err = reader
        .convert_parquet_to_targeted_experiment(&duplicate_dir, &mut duplicate_out)
        .unwrap_err();
    assert!(matches!(err, Error::InvalidValue(_)), "{err:?}");
    assert!(duplicate_out.compounds.is_empty());
}

/// `TransitionParquetFile_test.cpp`, section
/// `convertLightTargetedExperimentToParquet`.
#[test]
fn class_test_convert_light_targeted_experiment_to_parquet() {
    let all = records();
    let mut light = experiment(&all, "traml_input");
    assert!(!light.compounds.is_empty());
    assert!(!light.transitions.is_empty());

    let tmp = tmp();
    let out_dir = tmp.path().join("roundtrip.oswpq");
    std::fs::create_dir_all(&out_dir).unwrap();
    let source_ids = OpenSwathLibraryIDNormalizer::normalize_source_ids(&mut light).unwrap();
    let writer = TransitionParquetFile::new();
    writer
        .convert_light_targeted_experiment_to_parquet_with_source_ids(
            &out_dir,
            &light,
            Some(&source_ids),
        )
        .unwrap();
    assert!(out_dir.join("library/precursors.parquet").exists());
    assert!(out_dir.join("library/transitions.parquet").exists());
    assert!(out_dir.join("library/metadata.json").exists());

    let reader = TransitionParquetFile::new();
    let mut roundtrip = LightTargetedExperiment::default();
    let mut roundtrip_ids = SourceIDMapping::default();
    reader
        .convert_parquet_to_targeted_experiment_with_source_ids(
            &out_dir,
            &mut roundtrip,
            Some(&mut roundtrip_ids),
        )
        .unwrap();
    assert_eq!(
        roundtrip_ids.precursor_source_to_canonical.len(),
        source_ids.precursor_source_to_canonical.len()
    );
    for (source, canonical) in &source_ids.precursor_source_to_canonical {
        assert_eq!(
            roundtrip_ids.precursor_source_to_canonical[source],
            *canonical
        );
    }
    assert_eq!(
        roundtrip_ids.transition_canonical_to_source.len(),
        source_ids.transition_canonical_to_source.len()
    );
    for (canonical, source) in &source_ids.transition_canonical_to_source {
        assert_eq!(
            roundtrip_ids.transition_canonical_to_source[canonical],
            *source
        );
    }
    OpenSwathLibraryIDNormalizer::validate_canonical_ids(&roundtrip).unwrap();
    assert_eq!(roundtrip.compounds.len(), light.compounds.len());
    assert_eq!(roundtrip.transitions.len(), light.transitions.len());
    assert_eq!(roundtrip.proteins.len(), light.proteins.len());
    let ids: BTreeSet<&str> = roundtrip.compounds.iter().map(|c| c.id.as_str()).collect();
    for t in &roundtrip.transitions {
        assert!(ids.contains(t.peptide_ref.as_str()));
    }
    for (i, t) in roundtrip.transitions.iter().enumerate() {
        assert_eq!(t.transition_name, i.to_string());
        assert!(
            (t.product_mz - light.transitions[i].product_mz).abs() <= 1e-5 * t.product_mz.abs()
        );
    }

    // A canonical-looking but malformed experiment is refused, not renumbered.
    assert!(light.transitions.len() > 1);
    let mut malformed = light.clone();
    malformed.transitions[1].transition_name = malformed.transitions[0].transition_name.clone();
    let malformed_out = tmp.path().join("malformed.oswpq");
    let err = writer
        .convert_light_targeted_experiment_to_parquet(&malformed_out, &malformed)
        .unwrap_err();
    assert!(matches!(err, Error::InvalidValue(_)), "{err:?}");
    assert!(!malformed_out.exists());
}

// ------------------------------------------------------------------ oracle: reader

/// Reader cases: (oracle case, container source, read with source ids).
const READ_CASES: &[(&str, &str, bool)] = &[
    ("traml_dir_read", "traml_dir_write", true),
    ("traml_dir_read2", "traml_dir_write", false),
    ("traml_norm_read", "traml_norm_write", true),
    ("syn_read", "syn_write", true),
    ("w_empty_read", "w_empty", true),
    ("r_basic", "r_basic", true),
    ("r_list_proteins", "r_list_proteins", true),
    ("r_extra_column", "r_extra_column", true),
    ("r_missing_traml_col", "r_missing_traml_col", true),
    ("r_missing_annotation_col", "r_missing_annotation_col", true),
    ("r_missing_type_col", "r_missing_type_col", true),
    ("r_charge_int64", "r_charge_int64", true),
    ("r_dup_precursor", "r_dup_precursor", true),
    ("r_unknown_precursor", "r_unknown_precursor", true),
    ("r_null_mz", "r_null_mz", true),
    ("r_ambiguous_traml", "r_ambiguous_traml", true),
    ("r_ambiguous_traml_noids", "r_ambiguous_traml", false),
    ("r_dup_transition", "r_dup_transition", true),
    ("r_negative_precursor", "r_negative_precursor", true),
    ("r_wrapping_casts", "r_wrapping_casts", true),
    (
        "r_missing_transitions_file",
        "r_missing_transitions_file",
        true,
    ),
];

fn read(dir: &Path, with_ids: bool) -> (Result<()>, LightTargetedExperiment, SourceIDMapping) {
    let mut exp = LightTargetedExperiment::default();
    exp.compounds.push(LightCompound {
        id: "stale".into(),
        ..LightCompound::default()
    });
    let mut ids = SourceIDMapping::default();
    ids.precursor_source_to_canonical
        .insert("stale".into(), "stale".into());
    let reader = TransitionParquetFile::new();
    let result = if with_ids {
        reader.convert_parquet_to_targeted_experiment_with_source_ids(dir, &mut exp, Some(&mut ids))
    } else {
        reader.convert_parquet_to_targeted_experiment(dir, &mut exp)
    };
    (result, exp, ids)
}

type Result<T> = std::result::Result<T, Error>;

/// Rust compound order is precursor row order.
fn assert_row_order(case: &str, dir: &Path, exp: &LightTargetedExperiment) {
    let table = ParquetFile::read_table(dir.join("library/precursors.parquet")).unwrap();
    let ids = ParquetFile::get_column(&table, "precursor_id").unwrap();
    let rows: Vec<String> = (0..table.num_rows())
        .map(|r| {
            ParquetFile::get_int64(Some(ids.as_ref()), r, 0, false)
                .unwrap()
                .to_string()
        })
        .collect();
    let got: Vec<String> = exp.compounds.iter().map(|c| c.id.clone()).collect();
    assert_eq!(got, rows, "{case} compound order");
}

fn check_read_case(all: &[Vec<String>], case: &str, dir: &Path, with_ids: bool) {
    let (result, exp, ids) = read(dir, with_ids);
    let record = outcome(all, case);
    if record[0] == "OK" {
        result.unwrap_or_else(|e| panic!("{case}: {e:?}"));
        assert_experiment_matches(case, &exp, all);
        assert_row_order(case, dir, &exp);
        if with_ids {
            assert_eq!(maps_of(&ids), id_maps(all, case), "{case} source ids");
        }
    } else {
        let err = result.unwrap_err();
        assert_error_matches(case, &err, record);
        assert!(exp.compounds.is_empty() && exp.transitions.is_empty());
        if with_ids {
            assert_eq!(ids, SourceIDMapping::default(), "{case}");
        }
    }
}

#[test]
fn oracle_reader_on_release_written_containers() {
    let all = records();
    let tmp = tmp();
    for &(case, source, with_ids) in READ_CASES {
        let dir = container(tmp.path(), source, case);
        check_read_case(&all, case, &dir, with_ids);
    }
    // Zipped containers: the same library read through a store-only archive.
    for (case, source) in [("traml_zip_read", "traml_zip_write"), ("r_zip", "r_basic")] {
        let dir = container(tmp.path(), source, &format!("{case}_src"));
        let zip = tmp.path().join(format!("{case}.oswpq"));
        ZipArchiveFile::zip_directory(&dir, &zip).unwrap();
        check_read_case(&all, case, &zip, true);
    }
    // A path that does not exist.
    let missing = tmp.path().join("does_not_exist.oswpq");
    let (result, _, _) = read(&missing, true);
    assert_error_matches(
        "r_missing_dir",
        &result.unwrap_err(),
        outcome(&all, "r_missing_dir"),
    );
}

/// Documented divergence: the source lets `AASequence::fromString`'s
/// `ParseError` (here a missing `)`) escape and abort the read; the port
/// treats every unparseable sequence like an unknown modification and keeps
/// the compound without modifications.
#[test]
fn oracle_reader_parse_error_divergence() {
    let all = records();
    assert_eq!(
        outcome(&all, "r_parse_error")[2],
        "OpenMS::Exception::ParseError"
    );
    let tmp = tmp();
    let dir = container(tmp.path(), "r_parse_error", "r_parse_error");
    let (result, exp, _) = read(&dir, true);
    result.unwrap();
    let compound = exp.compounds.iter().find(|c| c.id == "3").unwrap();
    assert_eq!(compound.sequence, "PEPT(UniMod:35");
    assert!(compound.modifications.is_empty());
    // Every other record equals r_basic, which the Release build read.
    let basic = experiment(&all, "r_basic");
    for want in basic.compounds.iter().filter(|c| c.id != "3") {
        let got = exp.compounds.iter().find(|c| c.id == want.id).unwrap();
        assert_compound_eq("r_parse_error", got, want);
    }
}

// ------------------------------------------------------------------ oracle: writer

fn syn_ids() -> SourceIDMapping {
    let mut ids = SourceIDMapping::default();
    ids.precursor_canonical_to_source
        .insert("5".into(), "PEP_five".into());
    ids.precursor_canonical_to_source
        .insert("17".into(), "DECOY_PEP".into());
    ids.transition_canonical_to_source
        .insert("4".into(), "tr_four".into());
    ids.transition_canonical_to_source
        .insert("0".into(), "tr_zero".into());
    ids
}

#[test]
fn oracle_writer_tables_match_release() {
    let all = records();
    let tmp = tmp();
    let writer = TransitionParquetFile::new();
    let traml = experiment(&all, "traml_input");

    // 2-arg compatibility path on source-style ids.
    assert_eq!(outcome(&all, "traml_dir_write")[0], "OK");
    let dir = tmp.path().join("traml_dir.oswpq");
    std::fs::create_dir_all(&dir).unwrap();
    writer
        .convert_light_targeted_experiment_to_parquet(&dir, &traml)
        .unwrap();
    assert_library_matches("traml_dir_write", &dir);
    assert!(!tmp.path().join("traml_dir.oswpq.idx.json").exists());

    // Explicit normalisation, then the 3-arg writer.
    let mut normalized = traml.clone();
    let ids = OpenSwathLibraryIDNormalizer::normalize_source_ids(&mut normalized).unwrap();
    let want = experiment(&all, "traml_norm_input");
    for (g, w) in normalized.compounds.iter().zip(&want.compounds) {
        assert_compound_eq("traml_norm_input", g, w);
    }
    for (g, w) in normalized.transitions.iter().zip(&want.transitions) {
        assert_transition_eq("traml_norm_input", g, w);
    }
    assert_eq!(maps_of(&ids), id_maps(&all, "traml_norm_input"));
    let dir = tmp.path().join("traml_norm.oswpq");
    std::fs::create_dir_all(&dir).unwrap();
    writer
        .convert_light_targeted_experiment_to_parquet_with_source_ids(&dir, &normalized, Some(&ids))
        .unwrap();
    assert_library_matches("traml_norm_write", &dir);

    // Synthetic canonical library, with and without provenance.
    let syn = experiment(&all, "syn_input");
    for (case, provenance) in [("syn_write", Some(syn_ids())), ("syn_write_noids", None)] {
        assert_eq!(outcome(&all, case)[0], "OK");
        let dir = tmp.path().join(format!("{case}.oswpq"));
        std::fs::create_dir_all(&dir).unwrap();
        match &provenance {
            Some(ids) => writer
                .convert_light_targeted_experiment_to_parquet_with_source_ids(&dir, &syn, Some(ids))
                .unwrap(),
            None => writer
                .convert_light_targeted_experiment_to_parquet(&dir, &syn)
                .unwrap(),
        }
        assert_library_matches(case, &dir);
    }

    // Empty library.
    let dir = tmp.path().join("w_empty.oswpq");
    std::fs::create_dir_all(&dir).unwrap();
    writer
        .convert_light_targeted_experiment_to_parquet(&dir, &LightTargetedExperiment::default())
        .unwrap();
    assert_library_matches("w_empty", &dir);
}

#[test]
fn oracle_writer_zip_layout_matches_release() {
    let all = records();
    let tmp = tmp();
    let writer = TransitionParquetFile::new();
    let traml = experiment(&all, "traml_input");

    let zip = tmp.path().join("traml.oswpq");
    writer
        .convert_light_targeted_experiment_to_parquet(&zip, &traml)
        .unwrap();
    let want: Vec<String> = of_case(&all, "Z", "traml_zip_write")
        .into_iter()
        .filter(|r| r[2].parse::<usize>().is_ok())
        .map(|r| r[3].clone())
        .collect();
    assert_eq!(ZipArchiveFile::list_entries(&zip), want);
    assert!(!tmp.path().join("traml.oswpq.tmp").exists());

    // The index names the three library entries (sizes are the encoder's).
    let mut extract = None;
    let idx =
        ZipArchiveFile::extract_entry_to_temp_file(&zip, "traml.oswpq.tmp.idx.json", &mut extract)
            .unwrap();
    let idx_text = std::fs::read_to_string(idx).unwrap();
    let keys = |text: &str| {
        text.lines()
            .filter_map(|l| l.trim().strip_prefix('"'))
            .map(|l| l.split('"').next().unwrap().to_owned())
            .collect::<Vec<_>>()
    };
    let want_idx = of_case(&all, "Z", "traml_zip_write")
        .into_iter()
        .find(|r| r[2] == "idx_json")
        .unwrap()[3]
        .clone();
    assert_eq!(keys(&idx_text), keys(&want_idx));

    // The archived tables are the ones the Release build archived.
    let unpacked = tmp.path().join("unpacked");
    for entry in ["precursors.parquet", "transitions.parquet", "metadata.json"] {
        let path = ZipArchiveFile::extract_entry_to_temp_file(
            &zip,
            &format!("library/{entry}"),
            &mut extract,
        )
        .unwrap();
        std::fs::create_dir_all(unpacked.join("library")).unwrap();
        std::fs::copy(path, unpacked.join("library").join(entry)).unwrap();
    }
    assert_library_matches("traml_zip_write", &unpacked);

    // An existing non-archive file at the destination is replaced.
    let over = tmp.path().join("over.oswpq");
    std::fs::write(&over, "not a zip").unwrap();
    writer
        .convert_light_targeted_experiment_to_parquet(&over, &traml)
        .unwrap();
    let want: Vec<String> = of_case(&all, "Z", "zip_over_file")
        .into_iter()
        .map(|r| r[3].clone())
        .collect();
    assert_eq!(ZipArchiveFile::list_entries(&over), want);
}

#[test]
fn oracle_writer_errors_match_release() {
    let all = records();
    let tmp = tmp();
    let writer = TransitionParquetFile::new();
    let syn = experiment(&all, "syn_input");
    let traml = experiment(&all, "traml_input");

    let mut malformed = syn.clone();
    malformed.transitions[1].transition_name = malformed.transitions[0].transition_name.clone();
    let dir = tmp.path().join("w_malformed.oswpq");
    std::fs::create_dir_all(&dir).unwrap();
    let err = writer
        .convert_light_targeted_experiment_to_parquet(&dir, &malformed)
        .unwrap_err();
    assert_error_matches(
        "w_malformed_canonical",
        &err,
        outcome(&all, "w_malformed_canonical"),
    );
    let zip = tmp.path().join("w_malformed.zip.oswpq");
    let err = writer
        .convert_light_targeted_experiment_to_parquet(&zip, &malformed)
        .unwrap_err();
    assert_error_matches(
        "w_malformed_canonical_zip",
        &err,
        outcome(&all, "w_malformed_canonical_zip"),
    );
    assert!(!zip.exists());

    let mut no_transition = syn.clone();
    no_transition.compounds.push(LightCompound {
        id: "99".into(),
        sequence: "AAAK".into(),
        charge: 2,
        ..LightCompound::default()
    });
    let dir = tmp.path().join("w_no_transition.oswpq");
    std::fs::create_dir_all(&dir).unwrap();
    let err = writer
        .convert_light_targeted_experiment_to_parquet_with_source_ids(&dir, &no_transition, None)
        .unwrap_err();
    assert_error_matches("w_no_transition", &err, outcome(&all, "w_no_transition"));

    let dir = tmp.path().join("w_src3.oswpq");
    std::fs::create_dir_all(&dir).unwrap();
    let err = writer
        .convert_light_targeted_experiment_to_parquet_with_source_ids(&dir, &traml, None)
        .unwrap_err();
    assert_error_matches(
        "w_source_style_noncanonical_3arg",
        &err,
        outcome(&all, "w_source_style_noncanonical_3arg"),
    );
}

/// Documented divergence (see [`oracle_reader_parse_error_divergence`]): the
/// source's writer aborts on a `ParseError` from the unmodified-sequence
/// conversion; the port writes the row with an empty `unmodified_sequence`.
#[test]
fn oracle_writer_parse_error_divergence() {
    let all = records();
    assert_eq!(
        outcome(&all, "w_parse_error")[2],
        "OpenMS::Exception::ParseError"
    );
    let mut syn = experiment(&all, "syn_input");
    syn.compounds[4].sequence = "PEPT(UniMod:35".into();
    let tmp = tmp();
    let dir = tmp.path().join("w_parse_error.oswpq");
    std::fs::create_dir_all(&dir).unwrap();
    TransitionParquetFile::new()
        .convert_light_targeted_experiment_to_parquet(&dir, &syn)
        .unwrap();
    let table = ParquetFile::read_table(dir.join("library/precursors.parquet")).unwrap();
    let modified = ParquetFile::get_column(&table, "modified_sequence").unwrap();
    let unmodified = ParquetFile::get_column(&table, "unmodified_sequence").unwrap();
    assert_eq!(
        ParquetFile::get_string(Some(modified.as_ref()), 4).unwrap(),
        "PEPT(UniMod:35"
    );
    assert_eq!(
        ParquetFile::get_string(Some(unmodified.as_ref()), 4).unwrap(),
        ""
    );
}

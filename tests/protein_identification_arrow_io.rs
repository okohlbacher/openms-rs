// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

#![cfg(feature = "arrow")]

//! `ProteinIdentificationArrowIO` against the executed C++ Release build.
//!
//! Expected values are the records of the executed Linux x86_64 Release
//! install `openms4-release-bc9cc12-c19e494-174b576` (Arrow C++ 23.0.1), run
//! twice on ibminode06 with identical records (evidence tier 1):
//! `tests/data/arrow_parquet/cpp_release_records_3.tsv` (sha256
//! `b5bc983e77cff80caeaeb3ff9f488ae2def3e7263e776874f8c94cd9eb8c9493`).
//! The driver imported the protein tables of six C++-written bundles from
//! test-data `0cb15f2` (copied under `tests/data/arrow_parquet/fixtures/`)
//! and exported a constructed set of runs to the `pia_*.parquet` files
//! beside the records. Every run, hit and group field is compared, with
//! metadata sorted by key (the source's order follows its process-wide
//! registry).

use arrow_array::cast::AsArray;
use arrow_array::{Array, RecordBatch};
use openms::chemistry::{ModificationsDB, SequenceModification};
use openms::comparison::Tolerance;
use openms::format::arrow_io_helpers::ParquetWriteConfig;
use openms::format::arrow_schema_registry::validation::{Mode, validate};
use openms::format::arrow_schema_registry::{
    ConsensusFeatureSchema, FeatureSchema, PSMSchema, ProteinGroupSchema, ProteinSchema,
    SearchParamsSchema,
};
use openms::format::parquet_file::ParquetFile;
use openms::format::protein_identification_arrow_io::ProteinIdentificationArrowIO as Pia;
use openms::identification::{
    EnzymeTermSpecificity, PeakMassType, PeptideIdentification, ProteinGroup, ProteinHit,
    ProteinIdentification, ProteinModification,
};
use openms::kernel::DataArray;
use openms::metadata::{MetaInfo, MetaValue, MetaValueData};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

const RECORDS: &str = include_str!("data/arrow_parquet/cpp_release_records_3.tsv");

const BUNDLES: [&str; 6] = [
    "topp/IDMerger_idparquet_input1.idparquet",
    "topp/IDMerger_idparquet_input2.idparquet",
    "topp/IDRipper_4_input.idparquet",
    "topp/THIRDPARTY/PercolatorAdapter_idparquet_in.idparquet",
    "topp/FileConverter_9_output.featureparquet",
    "topp/FileConverter_14_output.consensusparquet",
];

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/arrow_parquet")
        .join(name)
}

fn fixture(bundle: &str, file: &str) -> PathBuf {
    data("fixtures").join(bundle).join(file)
}

fn unescape(text: &str) -> String {
    let raw = text.as_bytes();
    let mut bytes = Vec::new();
    let mut i = 0;
    while i < raw.len() {
        if raw[i] == b'\\' && i + 1 < raw.len() {
            match raw[i + 1] {
                b't' => bytes.push(b'\t'),
                b'n' => bytes.push(b'\n'),
                b'\\' => bytes.push(b'\\'),
                b'x' => {
                    let hex = std::str::from_utf8(&raw[i + 2..i + 4]).unwrap();
                    bytes.push(u8::from_str_radix(hex, 16).unwrap());
                    i += 2;
                }
                other => panic!("unknown escape {other}"),
            }
            i += 2;
        } else {
            bytes.push(raw[i]);
            i += 1;
        }
    }
    String::from_utf8(bytes).unwrap()
}

/// Records under `prefix` (kind, key, value), in order.
fn records_with_prefix(prefix: &str) -> Vec<(String, String, String)> {
    RECORDS
        .lines()
        .map(|line| {
            let mut parts = line.splitn(3, '\t');
            (
                parts.next().unwrap().to_owned(),
                unescape(parts.next().unwrap()),
                unescape(parts.next().unwrap_or("")),
            )
        })
        .filter(|(_, key, _)| key.starts_with(prefix))
        .collect()
}

fn rec(kind: &str, key: &str) -> String {
    RECORDS
        .lines()
        .find_map(|line| {
            let mut parts = line.splitn(3, '\t');
            let k = parts.next()?;
            let kk = unescape(parts.next()?);
            (k == kind && kk == key).then(|| unescape(parts.next().unwrap_or("")))
        })
        .unwrap_or_else(|| panic!("no record {kind}\t{key}"))
}

/// C `printf("%a")` as glibc prints a double.
fn hexd(v: f64) -> String {
    if v.is_nan() {
        return if v.is_sign_negative() { "-nan" } else { "nan" }.into();
    }
    let sign = if v.is_sign_negative() { "-" } else { "" };
    if v.is_infinite() {
        return format!("{sign}inf");
    }
    let bits = v.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i64;
    let fraction = bits & ((1u64 << 52) - 1);
    if exponent == 0 && fraction == 0 {
        return format!("{sign}0x0p+0");
    }
    let (lead, e) = if exponent == 0 {
        (0, -1022)
    } else {
        (1, exponent - 1023)
    };
    let digits = format!("{fraction:013x}");
    let digits = digits.trim_end_matches('0');
    let mantissa = if digits.is_empty() {
        format!("{lead}")
    } else {
        format!("{lead}.{digits}")
    };
    let plus = if e >= 0 { "+" } else { "" };
    format!("{sign}0x{mantissa}p{plus}{e}")
}

/// The driver's `describe(DataValue)`.
fn describe(value: &MetaValue) -> String {
    let list = |items: Vec<String>| format!("[{}]", items.join(", "));
    match value.data() {
        MetaValueData::String(text) => format!("STRING_VALUE|{text}"),
        MetaValueData::Integer(n) => format!("INT_VALUE|{n}"),
        MetaValueData::Float(x) => format!("DOUBLE_VALUE|{}", hexd(*x)),
        MetaValueData::StringList(v) => format!("STRING_LIST|{}", list(v.clone())),
        MetaValueData::IntegerList(v) => {
            format!("INT_LIST|{}", list(v.iter().map(i64::to_string).collect()))
        }
        MetaValueData::FloatList(v) => {
            let mut text = "DOUBLE_LIST".to_owned();
            for x in v {
                text.push('|');
                text.push_str(&hexd(*x));
            }
            text
        }
        MetaValueData::Empty => "EMPTY_VALUE|".to_owned(),
    }
}

/// The driver's `metas()`: entries sorted by key.
fn metas(meta: &MetaInfo) -> String {
    meta.iter()
        .map(|(k, v)| format!("{k}={};", describe(v)))
        .collect()
}

/// A recorded `metas()` without the `spectra_data` entries, which the port
/// keeps as the run's path fields instead of metadata.
fn without_run_paths(recorded: &str) -> String {
    recorded
        .split_inclusive(';')
        .filter(|entry| {
            !entry.starts_with("spectra_data=") && !entry.starts_with("spectra_data_raw=")
        })
        .collect()
}

/// `ProteinIdentification::getInferenceEngine` and `...Version`, the rule
/// the source's getters follow.
fn inference(run: &ProteinIdentification) -> (String, String) {
    let as_engine = matches!(
        run.search_engine.as_str(),
        "Fido" | "BayesianProteinInference" | "Epifany" | "ProteinInference"
    ) || (run.search_engine == "Percolator"
        && !run.indistinguishable_groups.is_empty());
    let text = |key: &str| {
        run.search_parameters
            .metadata
            .get(key)
            .map(ToString::to_string)
    };
    let engine = text("InferenceEngine").unwrap_or_else(|| {
        if as_engine {
            run.search_engine.clone()
        } else {
            String::new()
        }
    });
    let version = text("InferenceEngineVersion").unwrap_or_else(|| {
        if engine.is_empty() {
            String::new()
        } else {
            run.search_engine_version.clone()
        }
    });
    (engine, version)
}

/// `EnzymaticDigestion::Specificity` as the driver prints it.
fn specificity(value: EnzymeTermSpecificity) -> &'static str {
    match value {
        EnzymeTermSpecificity::None => "0",
        EnzymeTermSpecificity::Semi => "1",
        EnzymeTermSpecificity::Full => "2",
        EnzymeTermSpecificity::Unknown => "3",
    }
}

fn tolerance(t: Tolerance) -> String {
    match t {
        Tolerance::Absolute(v) => format!("{}|0", hexd(v)),
        Tolerance::Ppm(v) => format!("{}|1", hexd(v)),
    }
}

fn data_arrays(group: &ProteinGroup) -> String {
    let mut text = String::new();
    for a in &group.float_data_arrays {
        text.push_str(&format!("F:{}=", a.name));
        for v in &a.data {
            text.push_str(&format!("{},", hexd(f64::from(*v))));
        }
        text.push(';');
    }
    for a in &group.string_data_arrays {
        text.push_str(&format!("S:{}=", a.name));
        for v in &a.data {
            text.push_str(&format!("{v},"));
        }
        text.push(';');
    }
    for a in &group.integer_data_arrays {
        text.push_str(&format!("I:{}=", a.name));
        for v in &a.data {
            text.push_str(&format!("{v},"));
        }
        text.push(';');
    }
    text
}

/// Compare the port's runs with the driver's `dumpRuns(tag, ...)` records.
fn assert_runs_match(tag: &str, runs: &[ProteinIdentification]) {
    assert_eq!(runs.len().to_string(), rec("runs", tag), "{tag}");
    for (r, run) in runs.iter().enumerate() {
        let k = format!("{tag}/{r}");
        let field = |name: &str| rec("run", &format!("{k}/{name}"));
        let sp = &run.search_parameters;
        assert_eq!(run.identifier, field("identifier"), "{k}");
        assert_eq!(run.search_engine, field("search_engine"), "{k}");
        assert_eq!(
            run.search_engine_version,
            field("search_engine_version"),
            "{k}"
        );
        assert_eq!(
            run.date_time.clone().unwrap_or_else(|| "<invalid>".into()),
            field("date"),
            "{k}"
        );
        assert_eq!(run.score_type, field("score_type"), "{k}");
        assert_eq!(
            u8::from(run.higher_score_better).to_string(),
            field("higher_score_better"),
            "{k}"
        );
        assert_eq!(
            hexd(run.significance_threshold),
            field("significance_threshold"),
            "{k}"
        );
        let (engine, version) = inference(run);
        assert_eq!(engine, field("inference_engine"), "{k}");
        assert_eq!(version, field("inference_engine_version"), "{k}");
        // `<exists>|<path>|<path>...`: the port has no "set but empty" state.
        let paths = field("primary_ms_run_paths");
        let recorded: Vec<&str> = paths.split('|').skip(1).filter(|p| !p.is_empty()).collect();
        assert_eq!(run.primary_ms_run_paths, recorded, "{k}");
        let raw = field("raw_ms_run_paths");
        let recorded: Vec<&str> = raw.split('|').skip(1).filter(|p| !p.is_empty()).collect();
        assert_eq!(run.raw_ms_run_paths, recorded, "{k}");
        assert_eq!(
            metas(&run.metadata),
            without_run_paths(&field("metavalues")),
            "{k}"
        );
        assert_eq!(sp.database, field("db"), "{k}");
        assert_eq!(sp.database_version, field("db_version"), "{k}");
        assert_eq!(sp.taxonomy, field("taxonomy"), "{k}");
        assert_eq!(sp.charges, field("charges"), "{k}");
        assert_eq!(
            match sp.mass_type {
                PeakMassType::Monoisotopic => "MONOISOTOPIC",
                PeakMassType::Average => "AVERAGE",
            },
            field("mass_type"),
            "{k}"
        );
        assert_eq!(tolerance(sp.precursor_tolerance), field("precursor"), "{k}");
        assert_eq!(tolerance(sp.fragment_tolerance), field("fragment"), "{k}");
        assert_eq!(sp.digestion_enzyme, field("enzyme"), "{k}");
        assert_eq!(
            specificity(sp.enzyme_specificity),
            field("specificity"),
            "{k}"
        );
        assert_eq!(
            sp.missed_cleavages.to_string(),
            field("missed_cleavages"),
            "{k}"
        );
        assert_eq!(sp.fixed_modifications.join("|"), field("fixed"), "{k}");
        assert_eq!(
            sp.variable_modifications.join("|"),
            field("variable"),
            "{k}"
        );
        assert_eq!(metas(&sp.metadata), field("sp_metavalues"), "{k}");
        assert_eq!(run.hits.len().to_string(), field("hits"), "{k}");
        for (h, hit) in run.hits.iter().enumerate() {
            let mods: String = hit
                .modifications
                .iter()
                .map(|m| {
                    format!(
                        "{}:{};",
                        m.position,
                        m.modification.known().map(|k| k.full_id()).unwrap_or("?")
                    )
                })
                .collect();
            let text = format!(
                "{}|{}|{}|{}|{}|{}|{}",
                hit.accession,
                hexd(hit.score),
                hit.rank,
                hexd(hit.coverage.unwrap_or(-1.0)),
                hit.sequence,
                mods,
                metas(&hit.metadata)
            );
            assert_eq!(text, rec("hit", &format!("{k}/hit/{h}")), "{k}/hit/{h}");
        }
        for (what, groups) in [
            ("protein_groups", &run.protein_groups),
            ("indistinguishable", &run.indistinguishable_groups),
        ] {
            assert_eq!(groups.len().to_string(), field(what), "{k}/{what}");
            for (g, group) in groups.iter().enumerate() {
                let text = format!(
                    "{}|{}|{}",
                    hexd(group.probability),
                    group.accessions.join(","),
                    data_arrays(group)
                );
                assert_eq!(
                    text,
                    rec("group", &format!("{k}/{what}/{g}")),
                    "{k}/{what}/{g}"
                );
            }
        }
    }
    // Every record of the tag was consulted by one of the lookups above.
    let consulted = records_with_prefix(&format!("{tag}/")).len();
    let expected: usize = runs
        .iter()
        .map(|run| {
            // 26 run fields, 2 group counts, the hits and groups.
            28 + run.hits.len() + run.protein_groups.len() + run.indistinguishable_groups.len()
        })
        .sum();
    assert_eq!(consulted, expected, "{tag}: unexpected extra records");
}

/// Records `import <bundle>` and `run`/`hit`/`group` `<bundle>/...`: every
/// retained C++-written bundle imports as the executed source imports it.
#[test]
fn retained_bundles_import_as_in_the_release_build() {
    for bundle in BUNDLES {
        let mut runs = Vec::new();
        let warnings = Pia::import_from_parquet(
            fixture(bundle, "proteins.parquet"),
            fixture(bundle, "protein_groups.parquet"),
            fixture(bundle, "search_params.parquet"),
            &mut runs,
        )
        .unwrap_or_else(|e| panic!("{bundle}: {e}"));
        assert_eq!(rec("import", bundle), "1");
        assert!(warnings.is_empty(), "{bundle}: {warnings:?}");
        assert_runs_match(bundle, &runs);
    }
}

/// Every table of every retained bundle validates strictly against the
/// ported registry (the promise of the crate-decision entry).
#[test]
fn retained_bundle_tables_validate_strictly_against_the_registry() {
    for bundle in BUNDLES {
        let mut files = vec![
            ("proteins.parquet", ProteinSchema::schema()),
            ("protein_groups.parquet", ProteinGroupSchema::schema()),
            ("search_params.parquet", SearchParamsSchema::schema()),
            ("psms.parquet", PSMSchema::schema()),
        ];
        if bundle.ends_with("featureparquet") {
            files.push(("features.parquet", FeatureSchema::schema()));
        }
        if bundle.ends_with("consensusparquet") {
            files.push((
                "consensus_features.parquet",
                ConsensusFeatureSchema::schema(),
            ));
        }
        for (file, schema) in files {
            let table = ParquetFile::read_table(fixture(bundle, file)).unwrap();
            let result = validate(&table, &schema, Mode::Strict);
            assert!(result.valid, "{bundle}/{file}: {result}");
        }
    }
}

fn oxidation() -> SequenceModification {
    SequenceModification::Known(
        ModificationsDB::global()
            .get_modification_handle("Oxidation (M)", None, None)
            .unwrap(),
    )
}

fn data_array<T>(name: &str, data: Vec<T>) -> DataArray<T> {
    DataArray {
        name: name.into(),
        data,
        metadata: MetaInfo::new(),
        data_processing: Vec::new(),
    }
}

/// The driver's `constructed()` runs.
fn constructed() -> Vec<ProteinIdentification> {
    let mut a = ProteinIdentification {
        identifier: "run_A".into(),
        search_engine: "Comet".into(),
        search_engine_version: "2024.01".into(),
        date_time: Some("2023-11-14T22:13:20".into()),
        score_type: "q-value".into(),
        higher_score_better: false,
        significance_threshold: 0.05,
        primary_ms_run_paths: vec!["/data/a.mzML".into(), "/data/b.mzML".into()],
        ..ProteinIdentification::default()
    };
    a.metadata.insert("run_note".into(), "first".into());
    a.metadata.insert("run_count".into(), 7i64.into());
    let sp = &mut a.search_parameters;
    sp.database = "/db/human.fasta".into();
    sp.database_version = "2024_01".into();
    sp.taxonomy = "9606".into();
    sp.charges = "2:4".into();
    sp.mass_type = PeakMassType::Average;
    sp.precursor_tolerance = Tolerance::Ppm(10.0);
    sp.fragment_tolerance = Tolerance::Absolute(0.02);
    sp.digestion_enzyme = "Trypsin".into();
    sp.enzyme_specificity = EnzymeTermSpecificity::Semi;
    sp.missed_cleavages = 2;
    sp.fixed_modifications = vec!["Carbamidomethyl (C)".into()];
    sp.variable_modifications = vec!["Oxidation (M)".into(), "Acetyl (N-term)".into()];
    sp.metadata.insert("sp_int".into(), 3i64.into());
    sp.metadata.insert(
        "sp_list".into(),
        vec!["x".to_owned(), "y".to_owned()].into(),
    );
    sp.metadata
        .insert("InferenceEngine".into(), "Epifany".into());
    sp.metadata
        .insert("InferenceEngineVersion".into(), "1.2".into());
    let mut h1 = ProteinHit {
        score: 12.5,
        rank: 1,
        accession: "P1".into(),
        sequence: "MKCM".into(),
        coverage: Some(42.5),
        modifications: vec![
            ProteinModification {
                position: 0,
                modification: oxidation(),
            },
            ProteinModification {
                position: 3,
                modification: oxidation(),
            },
        ],
        ..ProteinHit::default()
    };
    h1.set_description("first protein");
    h1.metadata.insert("target_decoy".into(), "target".into());
    h1.metadata
        .insert("hit_double".into(), MetaValue::try_from(0.25).unwrap());
    h1.metadata.insert("hit_ints".into(), vec![1i64, -2].into());
    h1.metadata.insert(
        "hit_doubles".into(),
        MetaValue::try_from(vec![1.5, -0.5]).unwrap(),
    );
    let mut h2 = ProteinHit {
        score: 3.0,
        rank: 2,
        accession: "DECOY_P2".into(),
        ..ProteinHit::default()
    };
    h2.metadata.insert("target_decoy".into(), "decoy".into());
    h2.set_description("");
    let h3 = ProteinHit {
        score: -1.0,
        rank: 3,
        accession: "P3".into(),
        sequence: "PEPTIDE".into(),
        ..ProteinHit::default()
    };
    a.hits = vec![h1, h2, h3];
    a.protein_groups = vec![
        ProteinGroup {
            probability: 0.9,
            accessions: vec!["P1".into(), "DECOY_P2".into()],
            float_data_arrays: vec![data_array("abundances", vec![1.5f32, 0.1f32])],
            string_data_arrays: vec![data_array("labels", vec!["a".to_owned(), String::new()])],
            integer_data_arrays: vec![data_array("counts", vec![4, -7])],
        },
        ProteinGroup {
            probability: 0.4,
            accessions: vec!["P3".into()],
            ..ProteinGroup::default()
        },
    ];
    a.indistinguishable_groups = vec![ProteinGroup {
        probability: 1.0,
        accessions: vec!["P1".into()],
        ..ProteinGroup::default()
    }];
    let b = ProteinIdentification {
        identifier: "run_B".into(),
        search_engine: "Fido".into(),
        score_type: "Posterior Probability".into(),
        higher_score_better: true,
        hits: vec![ProteinHit {
            score: 0.5,
            rank: 1,
            accession: "Q1".into(),
            ..ProteinHit::default()
        }],
        ..ProteinIdentification::default()
    };
    vec![a, b]
}

/// A table with its `metavalues`-style columns' entries sorted by name
/// within each row, and without schema metadata, so that a C++-written table
/// (whose entry order follows its registry) and the port's compare equal.
fn normalized(table: &RecordBatch) -> Vec<Vec<String>> {
    let mut rows = vec![Vec::new(); table.num_rows()];
    for (index, field) in table.schema().fields().iter().enumerate() {
        let column = table.column(index);
        for (row, cells) in rows.iter_mut().enumerate() {
            let mut cell = if column.is_null(row) {
                "<null>".to_owned()
            } else {
                arrow_cell(column.as_ref(), row)
            };
            if field.name().ends_with("metavalues") {
                let mut entries: Vec<&str> = cell.split(";;").collect();
                entries.sort();
                cell = entries.join(";;");
            }
            cells.push(format!("{}={cell}", field.name()));
        }
    }
    rows
}

/// One cell as text, recursing into lists and structs.
fn arrow_cell(array: &dyn Array, row: usize) -> String {
    use arrow_array::types::*;
    use arrow_schema::DataType;
    if array.is_null(row) {
        return "<null>".into();
    }
    match array.data_type() {
        DataType::Utf8 => array.as_string::<i32>().value(row).to_owned(),
        DataType::Float64 => hexd(array.as_primitive::<Float64Type>().value(row)),
        DataType::Int32 => array.as_primitive::<Int32Type>().value(row).to_string(),
        DataType::Int64 => array.as_primitive::<Int64Type>().value(row).to_string(),
        DataType::Boolean => array.as_boolean().value(row).to_string(),
        DataType::Timestamp(unit, _) => format!(
            "{:?}:{}",
            unit,
            array.as_primitive::<TimestampMillisecondType>().value(row)
        ),
        DataType::List(_) => {
            let list = array.as_list::<i32>();
            let values = list.value(row);
            (0..values.len())
                .map(|i| arrow_cell(values.as_ref(), i))
                .collect::<Vec<_>>()
                .join(";;")
        }
        DataType::Struct(_) => {
            let s = array.as_struct();
            format!(
                "{{{}}}",
                (0..s.num_columns())
                    .map(|c| arrow_cell(s.column(c).as_ref(), row))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        }
        other => panic!("unexpected type {other}"),
    }
}

/// Records `export ...`: the three tables the port exports equal the ones
/// the executed source wrote for the same runs, cell by cell, and carry the
/// same file-level metadata keys.
#[test]
fn exported_tables_match_the_release_build() {
    let runs = constructed();
    assert_runs_match("constructed", &runs);
    let dir = std::env::temp_dir().join(format!("openms-pia-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let config = ParquetWriteConfig::default();
    let none = BTreeMap::new();
    let defs = BTreeMap::from([("run_B".to_owned(), "custom definitions text".to_owned())]);
    for (file, record) in [
        ("pia_proteins.parquet", "proteins"),
        ("pia_protein_groups.parquet", "protein_groups"),
        ("pia_search_params.parquet", "search_params"),
        ("pia_search_params_defs.parquet", "search_params_defs"),
    ] {
        assert_eq!(rec("export", record), "1");
        let path = dir.join(file);
        match record {
            "proteins" => Pia::export_proteins_to_parquet(&runs, &path, &config),
            "protein_groups" => Pia::export_protein_groups_to_parquet(&runs, &path, &config),
            "search_params" => Pia::export_search_params_to_parquet(&runs, &path, &config, &none),
            _ => Pia::export_search_params_to_parquet(&runs, &path, &config, &defs),
        }
        .unwrap();
        let ours = ParquetFile::read_table(&path).unwrap();
        let theirs = ParquetFile::read_table(data(file)).unwrap();
        assert_eq!(ours.schema().fields(), theirs.schema().fields(), "{file}");
        assert_eq!(normalized(&ours), normalized(&theirs), "{file}");
        let mut ours_keys: Vec<&String> = ours.schema().metadata().keys().collect();
        let mut theirs_keys: Vec<&String> = theirs.schema().metadata().keys().collect();
        ours_keys.sort();
        theirs_keys.sort();
        assert_eq!(ours_keys, theirs_keys, "{file}");
        for key in ["file_type", "creator", "software_provider"] {
            assert_eq!(
                ours.schema().metadata().get(key),
                theirs.schema().metadata().get(key),
                "{file} {key}"
            );
        }
    }
    assert_eq!(rec("export", "proteins file_type"), "proteins");
}

/// Records `roundtrip ...`: the port's export read back by the port, and
/// the C++-written files read by the port, reconstruct the runs the executed
/// source reconstructs.
#[test]
fn exported_runs_round_trip_as_in_the_release_build() {
    let dir = std::env::temp_dir().join(format!("openms-pia-rt-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let runs = constructed();
    let config = ParquetWriteConfig::default();
    Pia::export_proteins_to_parquet(&runs, dir.join("p.parquet"), &config).unwrap();
    Pia::export_protein_groups_to_parquet(&runs, dir.join("g.parquet"), &config).unwrap();
    Pia::export_search_params_to_parquet(&runs, dir.join("s.parquet"), &config, &BTreeMap::new())
        .unwrap();
    assert_eq!(rec("roundtrip", "ok"), "1");
    for (proteins, groups, search) in [
        (
            dir.join("p.parquet"),
            dir.join("g.parquet"),
            dir.join("s.parquet"),
        ),
        (
            data("pia_proteins.parquet"),
            data("pia_protein_groups.parquet"),
            data("pia_search_params.parquet"),
        ),
    ] {
        let mut back = Vec::new();
        let warnings = Pia::import_from_parquet(proteins, groups, search, &mut back).unwrap();
        assert!(warnings.is_empty(), "{warnings:?}");
        assert_runs_match("roundtrip", &back);
    }
}

/// Records `checkUniqueIdentifiers`, `synthesize ...` and `rename ...`.
#[test]
fn run_identifier_handling_matches_the_release_build() {
    let duplicate = vec![
        ProteinIdentification {
            identifier: "same".into(),
            ..ProteinIdentification::default()
        },
        ProteinIdentification {
            identifier: "same".into(),
            ..ProteinIdentification::default()
        },
    ];
    let record = rec("checkUniqueIdentifiers", "duplicate");
    let detail = record.rsplit("; ").next().unwrap();
    match Pia::check_unique_identifiers(&duplicate) {
        Err(openms::Error::InvalidValue(message)) => {
            assert!(record.starts_with("InvalidValue|"));
            assert!(message.contains(detail), "{message}");
            assert!(message.ends_with("same"));
        }
        other => panic!("{other:?}"),
    }
    assert!(Pia::check_unique_identifiers(&constructed()).is_ok());
    let mut runs = constructed();
    runs[1].identifier = "run_A".into();
    runs[1].search_engine = String::new();
    let (rename, warnings) = Pia::synthesize_run_identifiers(&mut runs);
    assert_eq!(rename.len().to_string(), rec("synthesize", "map size"));
    assert_eq!(warnings.len(), 1);
    let prefix = |id: &str| id[..id.rfind('_').unwrap()].to_owned();
    assert_eq!(
        prefix(&runs[0].identifier),
        rec("synthesize", "run 0 prefix")
    );
    assert_eq!(
        prefix(&runs[1].identifier),
        rec("synthesize", "run 1 prefix")
    );
    assert_eq!(
        u8::from(rename["run_A"] == runs[1].identifier).to_string(),
        rec("synthesize", "map value is run 1")
    );
    assert_ne!(runs[0].identifier, runs[1].identifier);
    let mut peptides = vec![
        PeptideIdentification {
            identifier: "run_A".into(),
            ..PeptideIdentification::default()
        },
        PeptideIdentification {
            identifier: "orphan".into(),
            ..PeptideIdentification::default()
        },
    ];
    Pia::apply_run_identifier_rename(&rename, &mut peptides);
    assert_eq!(
        u8::from(peptides[0].identifier == runs[1].identifier).to_string(),
        rec("rename", "matched")
    );
    assert_eq!(peptides[1].identifier, rec("rename", "orphan"));
}

/// Native: an unknown enzyme and an unknown modification are warnings, as
/// the source logs them, and a table of the wrong shape is refused without
/// touching the runs.
#[test]
fn import_warnings_and_refusals() {
    let dir = std::env::temp_dir().join(format!("openms-pia-warn-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let mut runs = constructed();
    runs[0].search_parameters.digestion_enzyme = "NoSuchEnzyme".into();
    let table = Pia::export_search_params_to_arrow(&runs, &BTreeMap::new()).unwrap();
    let mut back = Vec::new();
    let warnings = Pia::import_search_params_from_arrow(&table, &mut back).unwrap();
    assert_eq!(
        warnings,
        vec!["ProteinIdentificationArrowIO: Unknown enzyme 'NoSuchEnzyme'"]
    );
    assert_eq!(back[0].search_parameters.digestion_enzyme, "unknown_enzyme");
    // A hit whose run is missing creates the run and warns.
    let proteins = Pia::export_proteins_to_arrow(&runs).unwrap();
    let mut only = Vec::new();
    let warnings = Pia::import_proteins_from_arrow(&proteins, &mut only).unwrap();
    assert_eq!(only.len(), 2);
    assert_eq!(warnings.len(), 2);
    // A wrong-typed table is refused and leaves the runs alone.
    let before = only.clone();
    assert!(Pia::import_proteins_from_arrow(&table, &mut only).is_err());
    assert_eq!(only, before);
}

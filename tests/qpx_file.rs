// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `QPXFile` (`FORMAT/QPXFile.h`).
//!
//! Two kinds of evidence:
//!
//! * **Executed differential.** `tests/data/qpx/cpp_release_records.tsv` and
//!   the `psm_*.parquet` files beside it were produced by the oracle driver
//!   `F1c-qpx` linking the Linux x86_64 Release install
//!   `openms4-release-bc9cc12-c19e494-174b576` (g++ 14.4 `-O2 -DNDEBUG`,
//!   Arrow C++ 23.0.1) on ibminode06, run twice with identical records. The
//!   driver builds the identifications [`runs`] and [`peptides`] build here
//!   and records every cell of both export layouts, the warnings, the import,
//!   the Parquet files' row groups and footers, and the exceptions.
//! * **The upstream class test** `QPXFile_test.cpp` (core `bc9cc12`),
//!   section by section, where it needs only `QPXFile`.
#![cfg(feature = "arrow")]

use arrow_array::cast::AsArray;
use arrow_array::types::{Float32Type, Float64Type, Int16Type, Int32Type, Int64Type};
use arrow_array::{Array, RecordBatch};
use arrow_schema::DataType;
use openms::Error;
use openms::chemistry::AASequence;
use openms::format::arrow_io_helpers::{ParquetCompression, ParquetWriteConfig};
use openms::format::arrow_schema_registry::validation::cpp_field_string;
use openms::format::parquet_file::ParquetFile;
use openms::format::qpx_file::{QpxFile, StreamingOptions};
use openms::format::qpx_identity::FeatureLinks;
use openms::identification::{
    FlankingResidue, PeakAnnotation, PeptideEvidence, PeptideHit, PeptideIdentification,
    ProteinIdentification,
};
use openms::metadata::{MetaValue, MetaValueData};
use parquet::file::reader::FileReader;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

// ---------------------------------------------------------------------------
// Records
// ---------------------------------------------------------------------------

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/qpx")
        .join(name)
}

/// The driver's `esc` undone.
fn unescape(text: &str) -> String {
    let bytes = text.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'\\' && i + 1 < bytes.len() {
            match bytes[i + 1] {
                b't' => out.push(b'\t'),
                b'n' => out.push(b'\n'),
                b'\\' => out.push(b'\\'),
                b'x' if i + 3 < bytes.len() => {
                    let hex = std::str::from_utf8(&bytes[i + 2..i + 4]).unwrap();
                    out.push(u8::from_str_radix(hex, 16).unwrap());
                    i += 4;
                    continue;
                }
                other => {
                    out.push(b'\\');
                    out.push(other);
                }
            }
            i += 2;
        } else {
            out.push(bytes[i]);
            i += 1;
        }
    }
    String::from_utf8(out).unwrap()
}

/// `(section, key) -> value`, in file order.
fn records() -> &'static Vec<(String, String, String)> {
    static RECORDS: OnceLock<Vec<(String, String, String)>> = OnceLock::new();
    RECORDS.get_or_init(|| {
        let text = std::fs::read_to_string(data("cpp_release_records.tsv")).unwrap()
            + &std::fs::read_to_string(data("cpp_release_import_records.tsv")).unwrap();
        text.lines()
            .map(|line| {
                let mut parts = line.splitn(3, '\t');
                let section = parts.next().unwrap().to_owned();
                let key = unescape(parts.next().unwrap());
                let value = unescape(parts.next().unwrap_or(""));
                (section, key, value)
            })
            .collect()
    })
}

fn rec(section: &str, key: &str) -> &'static str {
    records()
        .iter()
        .find(|(s, k, _)| s == section && k == key)
        .map(|(_, _, v)| v.as_str())
        .unwrap_or_else(|| panic!("no record {section} {key}"))
}

fn recs_with_prefix(section: &str, prefix: &str) -> Vec<(&'static str, &'static str)> {
    records()
        .iter()
        .filter(|(s, k, _)| s == section && k.starts_with(prefix))
        .map(|(_, k, v)| (k.as_str(), v.as_str()))
        .collect()
}

/// The warnings the driver captured for `tag`, without the log colour codes,
/// the log's repeat summaries (`<...> occurred N times`) and the sequence
/// parser's own warning, which is not QPXFile's.
fn cpp_warnings(tag: &str) -> Vec<String> {
    let count: usize = rec("warn", &format!("{tag}/count")).parse().unwrap();
    (0..count)
        .map(|i| {
            rec("warn", &format!("{tag}/{i}"))
                .replace("\x1b[93m", "")
                .replace("\x1b[39m", "")
        })
        .filter(|line| !line.starts_with('<') && !line.starts_with("Warning: unknown modification"))
        .collect()
}

// ---------------------------------------------------------------------------
// The data set the driver builds
// ---------------------------------------------------------------------------

fn run(
    identifier: &str,
    score_type: &str,
    higher_better: bool,
    paths: &[&str],
) -> ProteinIdentification {
    ProteinIdentification {
        identifier: identifier.into(),
        score_type: score_type.into(),
        higher_score_better: higher_better,
        primary_ms_run_paths: paths.iter().map(|p| (*p).to_owned()).collect(),
        ..ProteinIdentification::default()
    }
}

fn runs() -> Vec<ProteinIdentification> {
    vec![
        run("run1", "q-value", false, &["/data/run1.mzML"]),
        run(
            "merged",
            "q-value",
            false,
            &["/a/fracA.mzML", "/b/fracB.mzML.gz"],
        ),
        run("nopath", "hyperscore", true, &[]),
    ]
}

fn hit(sequence: &str, charge: i32, score: f64) -> PeptideHit {
    PeptideHit {
        sequence: AASequence::parse(sequence).unwrap(),
        charge,
        score,
        ..PeptideHit::default()
    }
}

fn float(value: f64) -> MetaValue {
    MetaValue::try_from(value).unwrap()
}

fn peptide(identifier: &str, score_type: &str, higher_better: bool) -> PeptideIdentification {
    PeptideIdentification {
        identifier: identifier.into(),
        score_type: score_type.into(),
        higher_score_better: higher_better,
        ..PeptideIdentification::default()
    }
}

fn peptides(with_unattributable: bool) -> Vec<PeptideIdentification> {
    let mut peptides = Vec::new();

    let mut p0 = peptide("run1", "q-value", false);
    p0.rt = Some(1234.5);
    p0.mz = Some(500.25);
    p0.set_spectrum_reference("controllerType=0 controllerNumber=1 scan=1000");
    p0.metadata.insert("ion_mobility".into(), float(0.75));
    p0.metadata.insert("pep_note".into(), "hello".into());
    let mut h0 = hit("PEM(Oxidation)TIDER", 2, 0.01);
    h0.rank = 1;
    h0.metadata.insert("target_decoy".into(), "target".into());
    h0.metadata
        .insert("Posterior Error Probability".into(), float(0.02));
    h0.metadata.insert("predicted_RT".into(), float(1200.0));
    h0.metadata.insert("missed_cleavages".into(), 1i64.into());
    h0.metadata.insert("COMET:deltaCn".into(), float(0.5));
    h0.metadata.insert("q-value".into(), float(0.01));
    h0.metadata.insert("expect".into(), float(0.001));
    h0.metadata.insert("hyperscore".into(), 35i64.into());
    h0.metadata.insert(
        "custom_list".into(),
        vec!["a".to_owned(), "b".to_owned()].into(),
    );
    h0.evidences = vec![
        PeptideEvidence {
            protein_accession: "PROT_A".into(),
            aa_before: FlankingResidue::Residue('K'),
            aa_after: FlankingResidue::Residue('A'),
            start: Some(10),
            end: Some(17),
        },
        PeptideEvidence {
            protein_accession: "PROT_B".into(),
            ..PeptideEvidence::default()
        },
    ];
    h0.peak_annotations = vec![
        PeakAnnotation {
            mz: 175.119,
            intensity: 1000.0,
            charge: 1,
            annotation: "y1".into(),
        },
        PeakAnnotation {
            mz: 276.155,
            intensity: 500.0,
            charge: 1,
            annotation: "y2".into(),
        },
    ];
    p0.hits = vec![h0, hit("SECONDK", 1, 0.5)];
    peptides.push(p0);

    let mut p1 = peptide("merged", "Posterior Error Probability", false);
    p1.mz = Some(600.3);
    p1.set_spectrum_reference("scan=42");
    p1.metadata.insert("id_merge_index".into(), 1i64.into());
    let mut h = hit("(Acetyl)PEPTIDEC(Carbamidomethyl)K", 3, 0.05);
    h.metadata.insert("target_decoy".into(), "decoy".into());
    h.metadata.insert("IM".into(), float(1.1));
    h.metadata.insert("predicted_rt".into(), float(88.0));
    h.evidences = vec![PeptideEvidence {
        protein_accession: "PROT_C".into(),
        aa_before: FlankingResidue::NTerminus,
        aa_after: FlankingResidue::CTerminus,
        start: Some(0),
        end: Some(9),
    }];
    p1.hits = vec![h];
    peptides.push(p1);

    let mut p2 = peptide("merged", "q-value", true);
    p2.set_spectrum_reference("scan=43");
    p2.metadata.insert("id_merge_index".into(), 0i64.into());
    peptides.push(p2);

    if with_unattributable {
        let mut p3 = peptide("nopath", "hyperscore", true);
        p3.rt = Some(50.0);
        p3.mz = Some(400.0);
        p3.set_spectrum_reference("index=5");
        p3.hits = vec![hit("PEPTIDEK", 0, 12.5)];
        peptides.push(p3);
    }

    let mut p4 = peptide("run1", "q-value", false);
    p4.rt = Some(1300.25);
    p4.mz = Some(450.5);
    p4.set_spectrum_reference("controllerType=0 controllerNumber=1 scan=1001");
    p4.metadata.insert(
        "reference_file_name".into(),
        "/elsewhere/other_run.mzML".into(),
    );
    p4.hits = vec![
        hit("PEPTIDEK.(Amidated)", 2, 0.2),
        hit("PEPTM[+123.456]IDEK", 2, 0.3),
    ];
    peptides.push(p4);
    peptides
}

// ---------------------------------------------------------------------------
// Cell rendering, as the driver's `cell`
// ---------------------------------------------------------------------------

fn cell(array: &dyn Array, i: usize) -> String {
    if array.is_null(i) {
        return "null".into();
    }
    match array.data_type() {
        DataType::Utf8 => format!("\"{}\"", array.as_string::<i32>().value(i)),
        DataType::Int16 => format!("i:{}", array.as_primitive::<Int16Type>().value(i)),
        DataType::Int32 => format!("i:{}", array.as_primitive::<Int32Type>().value(i)),
        DataType::Int64 => format!("i:{}", array.as_primitive::<Int64Type>().value(i)),
        DataType::Float32 => format!(
            "f32:{:08x}",
            array.as_primitive::<Float32Type>().value(i).to_bits()
        ),
        DataType::Float64 => format!(
            "f64:{:016x}",
            array.as_primitive::<Float64Type>().value(i).to_bits()
        ),
        DataType::Boolean => format!("b:{}", array.as_boolean().value(i)),
        DataType::List(_) => {
            let list = array.as_list::<i32>();
            let offsets = list.value_offsets();
            let values = list.values();
            let items: Vec<String> = (offsets[i] as usize..offsets[i + 1] as usize)
                .map(|k| cell(values.as_ref(), k))
                .collect();
            format!("[{}]", items.join("\x1e"))
        }
        DataType::Struct(fields) => {
            let structure = array.as_struct();
            let items: Vec<String> = fields
                .iter()
                .enumerate()
                .map(|(f, field)| {
                    format!("{}={}", field.name(), cell(structure.column(f).as_ref(), i))
                })
                .collect();
            format!("{{{}}}", items.join(","))
        }
        other => format!("<type {other}>"),
    }
}

/// List columns whose entries follow metadata order, which differs by design
/// (the source's registry order, key order here): compared as sorted lists.
const UNORDERED: [&str; 3] = ["additional_scores", "psm_metavalues", "spectrum_metavalues"];

fn sorted_list(text: &str) -> String {
    let Some(body) = text.strip_prefix('[').and_then(|t| t.strip_suffix(']')) else {
        return text.to_owned();
    };
    if body.is_empty() {
        return text.to_owned();
    }
    let mut items: Vec<&str> = body.split('\x1e').collect();
    items.sort_unstable();
    format!("[{}]", items.join("\x1e"))
}

/// Relative tolerance of the internal layout's `calculated_mz`, a double.
///
/// `AASequence`'s port sums a modified residue's mass as its internal
/// formula's mass plus the modification's, where the source subtracts water
/// from the full residue's formula mass (`Residue::getMonoWeight(Internal)`),
/// so the sum can differ in the last bits (measured: 2 ULP, 4.6e-16
/// relative, on `PEM(Oxidation)TIDER` at charge 2). This is accumulation
/// order in the sequence arithmetic QPXFile calls, not a QPXFile difference;
/// the QPX view narrows the value to `float` and is compared exactly.
const CALCULATED_MZ_RELATIVE: f64 = 1e-15;

/// Whether two rendered cells agree: exactly, as sorted lists for the
/// metadata-ordered columns, and within [`CALCULATED_MZ_RELATIVE`] for a
/// double `calculated_mz`.
fn cells_agree(column: &str, actual: &str, expected: &str) -> bool {
    if UNORDERED.contains(&column) {
        return sorted_list(actual) == sorted_list(expected);
    }
    if column == "calculated_mz" {
        if let (Some(a), Some(b)) = (actual.strip_prefix("f64:"), expected.strip_prefix("f64:")) {
            let a = f64::from_bits(u64::from_str_radix(a, 16).unwrap());
            let b = f64::from_bits(u64::from_str_radix(b, 16).unwrap());
            return (a - b).abs() <= CALCULATED_MZ_RELATIVE * b.abs();
        }
    }
    actual == expected
}

/// Compare `table` with every `cell <tag>/...` record, and its row and
/// column counts and schema with the `table`/`schema` records.
fn assert_table_matches(tag: &str, table: &RecordBatch) {
    assert_eq!(
        rec("table", tag),
        format!("rows={} cols={}", table.num_rows(), table.num_columns()),
        "{tag}"
    );
    let schema: Vec<String> = table
        .schema()
        .fields()
        .iter()
        .map(|field| cpp_field_string(field))
        .collect();
    assert_eq!(rec("schema", tag), schema.join("\n"), "{tag} schema");
    let mut compared = 0;
    for (key, expected) in recs_with_prefix("cell", &format!("{tag}/")) {
        let mut parts = key.rsplitn(3, '/');
        let row: usize = parts.next().unwrap().parse().unwrap();
        let column = parts.next().unwrap();
        let array = table.column_by_name(column).unwrap();
        let actual = cell(array.as_ref(), row);
        assert!(
            cells_agree(column, &actual, expected),
            "{key}\n  ours: {actual}\ntheirs: {expected}"
        );
        compared += 1;
    }
    assert_eq!(compared, table.num_rows() * table.num_columns(), "{tag}");
}

// ---------------------------------------------------------------------------
// Differential against the Release build
// ---------------------------------------------------------------------------

#[test]
fn both_export_layouts_match_the_release_build_cell_by_cell() {
    let proteins = runs();
    let peptides = peptides(true);

    let internal_best = QpxFile::export_to_arrow(&proteins, &peptides, false).unwrap();
    assert_table_matches("internal_best", &internal_best.value);
    assert!(internal_best.warnings.is_empty());
    assert_eq!(cpp_warnings("internal_best"), Vec::<String>::new());

    let internal_all = QpxFile::export_to_arrow(&proteins, &peptides, true).unwrap();
    assert_table_matches("internal_all", &internal_all.value);
    assert_eq!(internal_all.warnings, cpp_warnings("internal_all"));

    let qpx_best = QpxFile::export_psms_to_qpx_arrow(&proteins, &peptides, false, None).unwrap();
    assert_table_matches("qpx_best", &qpx_best.value);
    assert_eq!(qpx_best.warnings, cpp_warnings("qpx_best"));

    let qpx_all = QpxFile::export_psms_to_qpx_arrow(&proteins, &peptides, true, None).unwrap();
    assert_table_matches("qpx_all", &qpx_all.value);
    // The source logs the same warning once per call; its log suppresses the
    // repeat, so the record holds none for this second call.
    assert_eq!(qpx_all.warnings.len(), 1);
    assert_eq!(qpx_all.warnings, qpx_best.warnings);

    // feature_id is looked up by psm_id.
    let first = qpx_best.value.column_by_name("psm_id").unwrap();
    let first = first.as_primitive::<Int64Type>().value(0);
    let links: FeatureLinks = [(first, 777)].into_iter().collect();
    let linked = QpxFile::export_psms_to_qpx_arrow(&proteins, &peptides, false, Some(&links))
        .unwrap()
        .value;
    let feature_id = linked.column_by_name("feature_id").unwrap();
    for row in 0..linked.num_rows() {
        assert_eq!(
            cell(feature_id.as_ref(), row),
            rec("cell", &format!("qpx_links/feature_id/{row}"))
        );
    }
}

#[test]
fn warnings_match_the_release_build() {
    let mut clash = runs();
    clash[1].primary_ms_run_paths = vec!["/x/run1.mzML".into(), "/y/fracB.mzML".into()];
    let collisions =
        QpxFile::export_psms_to_qpx_arrow(&clash, &peptides(false), false, None).unwrap();
    assert_eq!(collisions.warnings, cpp_warnings("collisions"));

    let mut duplicate = runs();
    duplicate[2].primary_ms_run_paths = vec!["/data/run1.mzML".into()];
    let exported = QpxFile::export_to_arrow(&duplicate, &peptides(false), false).unwrap();
    assert_eq!(exported.warnings, cpp_warnings("duplicate_paths"));
    assert_eq!(rec("table", "duplicate_paths"), "ok");
}

// ---------------------------------------------------------------------------
// Import
// ---------------------------------------------------------------------------

fn bits(value: f64) -> String {
    format!("f64:{:016x}", value.to_bits())
}

fn describe(value: &MetaValue) -> String {
    let list = |items: Vec<String>| format!("[{}]", items.join(", "));
    match value.data() {
        MetaValueData::String(text) => format!("STRING_VALUE|{text}"),
        MetaValueData::Integer(number) => format!("INT_VALUE|{number}"),
        MetaValueData::Float(number) => format!("DOUBLE_VALUE|{}", bits(*number)),
        MetaValueData::StringList(items) => format!("STRING_LIST|{}", list(items.clone())),
        MetaValueData::IntegerList(items) => format!(
            "INT_LIST|{}",
            list(items.iter().map(i64::to_string).collect())
        ),
        other => panic!("not in the data set: {other:?}"),
    }
}

fn metas<'a>(entries: impl Iterator<Item = (&'a String, &'a MetaValue)>, rank: u32) -> String {
    let mut sorted: BTreeMap<String, String> = entries
        .map(|(key, value)| (key.clone(), describe(value)))
        .collect();
    if rank != 0 {
        // The source keeps the rank as this metadata value.
        sorted.insert("rank".into(), format!("INT_VALUE|{rank}"));
    }
    sorted
        .into_iter()
        .map(|(key, value)| format!("{key}={value};"))
        .collect()
}

fn flank(residue: FlankingResidue) -> char {
    residue.code()
}

fn position(value: Option<usize>) -> String {
    value.map_or("-1".into(), |at| at.to_string())
}

fn assert_import_matches(
    tag: &str,
    proteins: &[ProteinIdentification],
    peptides: &[PeptideIdentification],
) {
    assert_eq!(
        rec("import", &format!("{tag}/proteins")),
        proteins.len().to_string()
    );
    for (i, run) in proteins.iter().enumerate() {
        let paths: String = run
            .primary_ms_run_paths
            .iter()
            .map(|path| format!("{path}|"))
            .collect();
        assert_eq!(
            rec("import", &format!("{tag}/protein/{i}")),
            format!(
                "{}|{}|{}|{}",
                run.identifier,
                run.score_type,
                u8::from(run.higher_score_better),
                paths
            )
        );
    }
    assert_eq!(
        rec("import", &format!("{tag}/peptides")),
        peptides.len().to_string()
    );
    for (i, peptide) in peptides.iter().enumerate() {
        let key = format!("{tag}/peptide/{i}");
        assert_eq!(
            rec("import", &key),
            format!(
                "{}|{}|{}|{}|{}|{}",
                peptide.identifier,
                peptide.score_type,
                u8::from(peptide.higher_score_better),
                peptide.rt.map_or("nort".into(), bits),
                peptide.mz.map_or("nomz".into(), bits),
                metas(peptide.metadata.iter(), 0)
            ),
            "{key}"
        );
        for (h, hit) in peptide.hits.iter().enumerate() {
            let evidences: String = hit
                .evidences
                .iter()
                .map(|e| {
                    format!(
                        "{},{},{},{},{};",
                        e.protein_accession,
                        flank(e.aa_before),
                        flank(e.aa_after),
                        position(e.start),
                        position(e.end)
                    )
                })
                .collect();
            assert_eq!(
                rec("import", &format!("{key}/hit/{h}")),
                format!(
                    "{}|{}|{}|{}|{}|{}",
                    hit.sequence,
                    hit.charge,
                    bits(hit.score),
                    hit.rank,
                    evidences,
                    metas(hit.metadata.iter(), hit.rank)
                ),
                "{key}/hit/{h}"
            );
        }
    }
}

/// `import_all` is recorded in the process that also built the data set, so
/// the C++ import resolves the mass delta of `PEPTM[+123.4560]IDEK` against a
/// `ModificationsDB` that `AASequence::fromString("PEPTM[+123.456]IDEK")` has
/// already extended ("adding it to the database"). The port keeps no such
/// process-wide state; it matches `import_fresh`, the same table imported by
/// a fresh process (read from the C++-written `internal_all.parquet`). Every
/// record of the in-process import that does not involve that hit is
/// compared too.
#[test]
fn import_matches_the_release_build() {
    let internal = QpxFile::export_to_arrow(&runs(), &peptides(true), true)
        .unwrap()
        .value;
    let mut proteins = Vec::new();
    let mut imported = Vec::new();
    let warnings = QpxFile::import_from_arrow(&internal, &mut proteins, &mut imported).unwrap();
    assert_eq!(rec("import", "import_fresh/ok"), "true");
    assert_eq!(warnings, cpp_warnings("import_fresh"));
    assert_import_matches("import_fresh", &proteins, &imported);

    // The C++-written table imports to the same identifications.
    let theirs = ParquetFile::read_table(data("internal_all.parquet")).unwrap();
    let mut their_proteins = Vec::new();
    let mut their_imported = Vec::new();
    let their_warnings =
        QpxFile::import_from_arrow(&theirs, &mut their_proteins, &mut their_imported).unwrap();
    assert_eq!(their_warnings, warnings);
    assert_eq!(their_proteins, proteins);
    assert_eq!(their_imported, imported);

    // The in-process record differs only in the mass-tag hit and its warnings.
    assert_eq!(rec("import", "import_all/ok"), "true");
    assert_eq!(cpp_warnings("import_all"), Vec::<String>::new());
    let fresh: Vec<(&str, &str)> = recs_with_prefix("import", "import_fresh/")
        .into_iter()
        .filter(|(key, _)| !key.ends_with("/ok"))
        .collect();
    let in_process: Vec<(&str, &str)> = recs_with_prefix("import", "import_all/")
        .into_iter()
        .filter(|(key, _)| !key.ends_with("/ok"))
        .collect();
    assert_eq!(fresh.len(), in_process.len());
    let mass_tag_hit = "peptide/3/hit/1";
    for ((fresh_key, fresh_value), (key, value)) in fresh.iter().zip(&in_process) {
        assert_eq!(
            &fresh_key["import_fresh/".len()..],
            &key["import_all/".len()..]
        );
        if key.ends_with(mass_tag_hit) {
            assert!(value.starts_with("PEPTM[+123.456]IDEK|"), "{value}");
            assert!(fresh_value.starts_with("PEPTMIDEK|"), "{fresh_value}");
        } else {
            assert_eq!(fresh_value, value, "{key}");
        }
    }

    // Re-exported as QPX: the per-PSM origin survives the round trip. The
    // mass-tag row is the one whose sequence the fresh import could not
    // reconstruct, so its peptidoform, modifications and m/z differ from the
    // in-process record; every other cell is compared.
    let again = QpxFile::export_psms_to_qpx_arrow(&proteins, &imported, true, None).unwrap();
    let skipped = [
        "psm_id",
        "peptidoform",
        "modifications",
        "calculated_mz",
        "mass_error_ppm",
    ];
    for (key, expected) in recs_with_prefix("cell", "qpx_reimport/") {
        let mut parts = key.rsplitn(3, '/');
        let row: usize = parts.next().unwrap().parse().unwrap();
        let column = parts.next().unwrap();
        if row == 5 && skipped.contains(&column) {
            continue;
        }
        let array = again.value.column_by_name(column).unwrap();
        let actual = cell(array.as_ref(), row);
        assert!(
            cells_agree(column, &actual, expected),
            "{key}\n  ours: {actual}\ntheirs: {expected}"
        );
    }
}

/// The class test's lossy peptidoform, and one that cannot be parsed.
#[test]
fn import_warnings_match_the_release_build() {
    let proteins = vec![run("run_1", "score", true, &[])];
    let mut identification = peptide("run_1", "score", true);
    identification.rt = Some(17.8);
    identification.mz = Some(712.7828);
    identification.set_spectrum_reference("scan=9");
    identification.hits = vec![hit("AEADNLDDKK", 2, 0.1)];
    let table = QpxFile::export_to_arrow(&proteins, &[identification], true)
        .unwrap()
        .value;
    let replace = |table: &RecordBatch, name: &str, values: arrow_array::ArrayRef| {
        let index = table.schema().index_of(name).unwrap();
        let mut columns = table.columns().to_vec();
        columns[index] = values;
        RecordBatch::try_new(table.schema(), columns).unwrap()
    };
    let nuxl = replace(
        &table,
        "peptidoform",
        std::sync::Arc::new(arrow_array::StringArray::from(vec![
            "AEADNLDDK[NuXL:U-H2O]K",
        ])),
    );
    let nuxl = replace(
        &nuxl,
        "calculated_mz",
        std::sync::Arc::new(arrow_array::Float64Array::from(vec![712.7824985834709])),
    );
    for (tag, table) in [
        ("import_nuxl", nuxl.clone()),
        (
            "import_unparseable",
            replace(
                &nuxl,
                "peptidoform",
                std::sync::Arc::new(arrow_array::StringArray::from(vec!["AEAD[[K"])),
            ),
        ),
    ] {
        let mut out_proteins = proteins.clone();
        let mut out = Vec::new();
        let warnings = QpxFile::import_from_arrow(&table, &mut out_proteins, &mut out).unwrap();
        assert_eq!(rec("import", &format!("{tag}/ok")), "true");
        // The first warning quotes the ProForma parser's own message, whose
        // text is the parser port's (chemistry::proforma), not QPXFile's:
        // the source's reads "Expected ']' in: <input>" for the NuXL case.
        // QPXFile's framing around it is compared exactly.
        let expected = cpp_warnings(tag);
        assert_eq!(warnings.len(), expected.len(), "{tag}");
        let peptidoform = if tag == "import_nuxl" {
            "AEADNLDDK[NuXL:U-H2O]K"
        } else {
            "AEAD[[K"
        };
        let prefix = format!("QPXFile: failed to parse peptidoform '{peptidoform}' (");
        let suffix = "); falling back to the unmodified sequence column.";
        for text in [&warnings[0], &expected[0]] {
            assert!(
                text.starts_with(&prefix) && text.ends_with(suffix),
                "{text}"
            );
        }
        assert_eq!(warnings[1..], expected[1..], "{tag}");
        assert_import_matches(tag, &out_proteins, &out);
    }
}

// ---------------------------------------------------------------------------
// Parquet files
// ---------------------------------------------------------------------------

fn temp(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("qpx_file_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

/// `rows=<n>`, the row-group layout and the footer pairs (dynamic values
/// masked) of a Parquet file, as the driver's `dumpParquet` records them.
fn parquet_layout(path: &Path) -> (String, String, Vec<(String, String)>) {
    let reader = parquet::file::serialized_reader::SerializedFileReader::new(
        std::fs::File::open(path).unwrap(),
    )
    .unwrap();
    let metadata = reader.metadata();
    let groups: String = metadata
        .row_groups()
        .iter()
        .map(|group| format!("{},", group.num_rows()))
        .collect();
    let mut pairs: Vec<(String, String)> = metadata
        .file_metadata()
        .key_value_metadata()
        .map(|kv| {
            kv.iter()
                .filter(|kv| kv.key != "ARROW:schema")
                .map(|kv| {
                    let value = match kv.key.as_str() {
                        "uuid" | "creation_date" => "<dynamic>".to_owned(),
                        _ => kv.value.clone().unwrap_or_default(),
                    };
                    (kv.key.clone(), value)
                })
                .collect()
        })
        .unwrap_or_default();
    pairs.sort();
    (
        format!("{}:{groups}", metadata.num_row_groups()),
        metadata.file_metadata().num_rows().to_string(),
        pairs,
    )
}

/// The Rust file has the C++ file's rows, row groups and footer, the pairs
/// compared as a set (Parquet C++ keeps insertion order, the port key
/// order), and both files read back to the same table.
fn assert_file_matches(tag: &str, ours: &Path, theirs: &str) {
    let (groups, rows, pairs) = parquet_layout(ours);
    assert_eq!(
        groups,
        rec("parquet", &format!("{tag}/row_groups")),
        "{tag}"
    );
    assert_eq!(rows, rec("parquet", &format!("{tag}/rows")), "{tag}");
    let mut expected: Vec<(String, String)> = recs_with_prefix("parquet", &format!("{tag}/kv/"))
        .into_iter()
        .map(|(key, value)| (key[tag.len() + 4..].to_owned(), value.to_owned()))
        .collect();
    expected.sort();
    assert_eq!(pairs, expected, "{tag}");
    let (their_groups, their_rows, _) = parquet_layout(&data(theirs));
    assert_eq!((their_groups, their_rows), (groups, rows));

    let ours = ParquetFile::read_table(ours).unwrap();
    let theirs = ParquetFile::read_table(data(theirs)).unwrap();
    assert_eq!(ours.schema().fields(), theirs.schema().fields(), "{tag}");
    assert_eq!(ours.num_rows(), theirs.num_rows(), "{tag}");
    for (index, field) in ours.schema().fields().iter().enumerate() {
        for row in 0..ours.num_rows() {
            let actual = cell(ours.column(index).as_ref(), row);
            let expected = cell(theirs.column(index).as_ref(), row);
            assert!(
                cells_agree(field.name(), &actual, &expected),
                "{tag}/{}/{row}\n  ours: {actual}\ntheirs: {expected}",
                field.name()
            );
        }
    }
}

#[test]
fn parquet_files_match_the_release_build() {
    let proteins = runs();
    let peptides = peptides(false);
    let borrowed: Vec<&PeptideIdentification> = peptides.iter().collect();

    let oneshot = temp("psm_oneshot.parquet");
    let warnings = QpxFile::export_to_parquet(
        &proteins,
        &peptides,
        &oneshot,
        true,
        &ParquetWriteConfig::default(),
        None,
    )
    .unwrap();
    assert_eq!(rec("parquet", "oneshot/ok"), "true");
    assert_eq!(warnings, cpp_warnings("oneshot"));
    assert_file_matches("oneshot", &oneshot, "psm_oneshot.parquet");

    let stream = temp("psm_stream.parquet");
    let options = StreamingOptions {
        export_all_psms: true,
        batch_size: 2,
        config: ParquetWriteConfig {
            row_group_size: 2,
            ..ParquetWriteConfig::default()
        },
        n_threads: 1,
    };
    let warnings =
        QpxFile::export_to_parquet_streaming(&proteins, &borrowed, &stream, &options, None)
            .unwrap();
    assert_eq!(rec("parquet", "stream/ok"), "true");
    assert_eq!(warnings, cpp_warnings("stream"));
    assert_file_matches("stream", &stream, "psm_stream.parquet");

    let empty = temp("psm_empty.parquet");
    QpxFile::export_to_parquet_streaming(
        &proteins,
        &[],
        &empty,
        &StreamingOptions::default(),
        None,
    )
    .unwrap();
    assert_eq!(rec("parquet", "empty/ok"), "true");
    assert_file_matches("empty", &empty, "psm_empty.parquet");

    // A batch without rows (an identification without hits) is an empty
    // row group, as Parquet C++'s WriteTable writes one.
    let gap = temp("psm_gap.parquet");
    let options = StreamingOptions {
        batch_size: 1,
        ..StreamingOptions::default()
    };
    QpxFile::export_to_parquet_streaming(
        &proteins,
        &[borrowed[0], borrowed[2], borrowed[1]],
        &gap,
        &options,
        None,
    )
    .unwrap();
    assert_eq!(rec("parquet", "gap/ok"), "true");
    assert_file_matches("gap", &gap, "psm_gap.parquet");

    // An unattributable PSM (empty run_file_name) is refused by the value
    // validation before anything is written.
    let bad = temp("psm_bad.parquet");
    let with_bad = self::peptides(true);
    assert!(
        QpxFile::export_to_parquet(
            &proteins,
            &with_bad,
            &bad,
            false,
            &ParquetWriteConfig::default(),
            None
        )
        .is_err()
    );
    assert_eq!(rec("parquet", "unattributable/ok"), "false");
    assert!(!bad.exists());
    assert_eq!(rec("parquet", "unattributable"), "missing");

    // LZ4 has no QPX compression token; the internal table is not a QPX one.
    let lz4 = ParquetWriteConfig {
        compression: ParquetCompression::Lz4,
        ..ParquetWriteConfig::default()
    };
    assert!(
        QpxFile::export_to_parquet(&proteins, &peptides, temp("lz4.parquet"), false, &lz4, None)
            .is_err()
    );
    assert_eq!(rec("parquet", "lz4/ok"), "false");
    let internal = QpxFile::export_to_arrow(&proteins, &peptides, false)
        .unwrap()
        .value;
    assert!(
        QpxFile::export_table_to_parquet(
            &internal,
            temp("internal.parquet"),
            &ParquetWriteConfig::default(),
            ""
        )
        .is_err()
    );
    assert_eq!(rec("parquet", "internal_table/ok"), "false");
}

// ---------------------------------------------------------------------------
// Exceptions
// ---------------------------------------------------------------------------

fn message(error: Error) -> (String, String) {
    match error {
        Error::InvalidValue(message) => ("InvalidValue".into(), message),
        Error::MissingInformation(message) => ("MissingInformation".into(), message),
        other => ("other".into(), other.to_string()),
    }
}

#[test]
fn refusals_match_the_release_build() {
    let proteins = runs();
    let base = peptides(false);

    let mut duplicate = proteins.clone();
    duplicate.push(proteins[0].clone());
    let (kind, text) =
        message(QpxFile::require_resolvable_merge_indices(&duplicate, &base).unwrap_err());
    assert_eq!(
        format!("{kind}|{text}"),
        rec("exception", "duplicate_identifier")
    );

    let mut missing = base.clone();
    missing[1].metadata.remove("id_merge_index");
    let (kind, text) =
        message(QpxFile::require_resolvable_merge_indices(&proteins, &missing).unwrap_err());
    assert_eq!(format!("{kind}|{text}"), rec("exception", "missing_index"));
    let borrowed: Vec<&PeptideIdentification> = missing.iter().collect();
    let (kind, text) =
        message(QpxFile::require_resolvable_merge_indices_of(&proteins, &borrowed).unwrap_err());
    assert_eq!(format!("{kind}|{text}"), rec("exception", "missing_index"));
    let (kind, text) =
        message(QpxFile::export_psms_to_qpx_arrow(&proteins, &missing, false, None).unwrap_err());
    assert_eq!(
        format!("{kind}|{text}"),
        rec("exception", "export_missing_index")
    );

    let mut out_of_range = base.clone();
    out_of_range[1]
        .metadata
        .insert("id_merge_index".into(), 2i64.into());
    let (kind, text) =
        message(QpxFile::require_resolvable_merge_indices(&proteins, &out_of_range).unwrap_err());
    assert_eq!(format!("{kind}|{text}"), rec("exception", "out_of_range"));

    let mut not_integer = base.clone();
    not_integer[1]
        .metadata
        .insert("id_merge_index".into(), "1".into());
    let (kind, text) =
        message(QpxFile::require_resolvable_merge_indices(&proteins, &not_integer).unwrap_err());
    assert_eq!(format!("{kind}|{text}"), rec("exception", "not_integer"));

    // The source throws Exception::ConversionError; the port has no such
    // variant and returns the same message as Error::InvalidValue.
    let mut double = base.clone();
    double[0].hits[0]
        .metadata
        .insert("missed_cleavages".into(), float(1.5));
    let (_, text) =
        message(QpxFile::export_psms_to_qpx_arrow(&proteins, &double, false, None).unwrap_err());
    assert_eq!(
        format!("ConversionError|{text}"),
        rec("exception", "missed_cleavages_double")
    );
}

// ---------------------------------------------------------------------------
// The upstream class test, QPXFile_test.cpp
// ---------------------------------------------------------------------------

fn test_run(paths: &[&str]) -> ProteinIdentification {
    let mut run = run("test_search", "TestScore", true, paths);
    run.search_engine = "TestEngine".into();
    run
}

/// `static std::shared_ptr<arrow::Table> exportPSMsToQPXArrow(...)`.
#[test]
fn class_test_export_psms_to_qpx_arrow() {
    let proteins = vec![test_run(&[])];
    let mut peptides = Vec::new();
    for (i, sequence) in ["PEPTIDER", "PEM(Oxidation)TIDER", "DFPIANGER"]
        .iter()
        .enumerate()
    {
        let mut identification = peptide("test_search", "TestScore", true);
        identification.rt = Some(1234.5 + i as f64 * 100.0);
        identification.mz = Some(500.25 + i as f64 * 50.0);
        identification.set_spectrum_reference(format!(
            "controllerType=0 controllerNumber=1 scan={}",
            1000 + i
        ));
        let mut first = hit(sequence, 2 + i as i32, 0.95 - i as f64 * 0.1);
        first.metadata.insert(
            "target_decoy".into(),
            if i % 2 == 0 { "target" } else { "decoy" }.into(),
        );
        first.evidences = vec![PeptideEvidence {
            protein_accession: format!("TEST_PROTEIN_{i}"),
            ..PeptideEvidence::default()
        }];
        identification.hits = vec![first];
        if i == 0 {
            identification.hits.push(hit("SECOND", 1, 0.5));
        }
        peptides.push(identification);
    }
    let table = QpxFile::export_psms_to_qpx_arrow(&proteins, &peptides, false, None)
        .unwrap()
        .value;
    assert_eq!(table.num_rows(), 3);
    assert_eq!(table.num_columns(), 26);
    let names: Vec<String> = table
        .schema()
        .fields()
        .iter()
        .map(|f| f.name().clone())
        .collect();
    assert_eq!(
        names,
        [
            "psm_id",
            "sequence",
            "peptidoform",
            "modifications",
            "charge",
            "posterior_error_probability",
            "is_decoy",
            "calculated_mz",
            "observed_mz",
            "mass_error_ppm",
            "additional_scores",
            "predicted_rt",
            "run_file_name",
            "cv_params",
            "scan",
            "rt",
            "ion_mobility",
            "missed_cleavages",
            "protein_accessions",
            "cross_links",
            "mz_array",
            "intensity_array",
            "charge_array",
            "ion_type_array",
            "ion_mobility_array",
            "feature_id",
        ]
    );
    let types: Vec<&DataType> = [0, 25, 4, 5, 7, 8, 9]
        .iter()
        .map(|i| table.schema_ref().field(*i).data_type())
        .collect();
    assert_eq!(
        types,
        [
            &DataType::Int64,
            &DataType::Int64,
            &DataType::Int16,
            &DataType::Float64,
            &DataType::Float32,
            &DataType::Float32,
            &DataType::Float32
        ]
    );
    for i in [14, 18, 3] {
        assert!(matches!(
            table.schema_ref().field(i).data_type(),
            DataType::List(_)
        ));
    }
    let text = |name: &str, row: usize| {
        table
            .column_by_name(name)
            .unwrap()
            .as_string::<i32>()
            .value(row)
            .to_owned()
    };
    assert_eq!(text("sequence", 0), "PEPTIDER");
    assert_eq!(text("sequence", 1), "PEMTIDER");
    assert_eq!(text("sequence", 2), "DFPIANGER");
    assert_eq!(text("peptidoform", 0), "PEPTIDER");
    assert!(text("peptidoform", 1).contains('['));
    let decoy = table.column_by_name("is_decoy").unwrap();
    let decoy = decoy.as_boolean();
    assert_eq!(
        (decoy.value(0), decoy.value(1), decoy.value(2)),
        (false, true, false)
    );
    let scan = table.column_by_name("scan").unwrap();
    for row in 0..3 {
        assert_eq!(cell(scan.as_ref(), row), format!("[i:{}]", 1000 + row));
    }
}

/// `[EXTRA] QPX value validation rejects empty and duplicate PSM keys before writing`.
#[test]
fn class_test_value_validation_before_writing() {
    use openms::format::qpx_value_validation::{QpxValueValidation, View};
    let proteins = vec![run("merged", "q-value", false, &["/data/runA.mzML"])];
    let peptides = vec![merged_psm("PEPTIDEK", Some(0))];
    let table = QpxFile::export_psms_to_qpx_arrow(&proteins, &peptides, false, None)
        .unwrap()
        .value;
    let mut validator = QpxValueValidation::new(View::Psm);
    assert!(validator.validate(&table).valid);
    let again = validator.validate(&table);
    assert!(!again.valid);
    assert!(again.to_string().contains("primary key"));
    validator.reset();
    assert!(validator.validate(&table).valid);
    let doubled = arrow_select::concat::concat_batches(&table.schema(), [&table, &table]).unwrap();
    assert!(!QpxValueValidation::new(View::Psm).validate(&doubled).valid);

    let index = table.schema().index_of("run_file_name").unwrap();
    let mut columns = table.columns().to_vec();
    columns[index] = std::sync::Arc::new(arrow_array::StringArray::from(vec![""]));
    let empty_run = RecordBatch::try_new(table.schema(), columns).unwrap();
    assert!(
        !QpxValueValidation::new(View::Psm)
            .validate(&empty_run)
            .valid
    );
    // (A physical null under the non-nullable run_file_name cannot be put
    // into an arrow-rs RecordBatch: try_new refuses it.)

    let output = temp("empty_run.parquet");
    assert!(
        QpxFile::export_table_to_parquet(&empty_run, &output, &ParquetWriteConfig::default(), "")
            .is_err()
    );
    assert!(!output.exists());
}

/// `exportPSMsToQPXArrow(...) with export_all_psms`.
#[test]
fn class_test_export_all_psms() {
    let proteins = vec![test_run(&[])];
    let mut identification = peptide("test_search", "TestScore", true);
    identification.rt = Some(1234.5);
    identification.mz = Some(500.25);
    for (i, sequence) in ["PEPTIDER", "ALTERNATIVE", "THIRDPSM"].iter().enumerate() {
        let mut h = hit(sequence, 2 + i as i32, 0.95 - i as f64 * 0.1);
        h.metadata.insert("target_decoy".into(), "target".into());
        h.evidences = vec![PeptideEvidence {
            protein_accession: format!("TEST_PROTEIN_{i}"),
            ..PeptideEvidence::default()
        }];
        identification.hits.push(h);
    }
    let table = QpxFile::export_psms_to_qpx_arrow(&proteins, &[identification], true, None)
        .unwrap()
        .value;
    assert_eq!(table.num_rows(), 3);
    let charge = table.column_by_name("charge").unwrap();
    let charge = charge.as_primitive::<Int16Type>();
    assert_eq!(
        (charge.value(0), charge.value(1), charge.value(2)),
        (2, 3, 4)
    );
    let sequence = table.column_by_name("sequence").unwrap();
    let sequence = sequence.as_string::<i32>();
    assert_eq!(
        (sequence.value(0), sequence.value(1), sequence.value(2)),
        ("PEPTIDER", "ALTERNATIVE", "THIRDPSM")
    );
}

/// `static bool exportToParquet(...)`.
#[test]
fn class_test_export_to_parquet() {
    let proteins = vec![test_run(&["test_run.mzML"])];
    let mut identification = peptide("test_search", "TestScore", true);
    identification.rt = Some(1234.5);
    identification.mz = Some(500.25);
    identification.set_spectrum_reference("controllerType=0 controllerNumber=1 scan=1234");
    let mut h = hit("PEM(Oxidation)TIDER", 2, 0.95);
    h.metadata.insert("target_decoy".into(), "target".into());
    h.evidences = vec![
        PeptideEvidence {
            protein_accession: "TEST_PROTEIN".into(),
            ..PeptideEvidence::default()
        },
        PeptideEvidence {
            protein_accession: "TEST_PROTEIN_2".into(),
            ..PeptideEvidence::default()
        },
    ];
    identification.hits = vec![h];
    let output = temp("class_export.parquet");
    QpxFile::export_to_parquet(
        &proteins,
        &[identification],
        &output,
        false,
        &ParquetWriteConfig::default(),
        None,
    )
    .unwrap();
    let table = ParquetFile::read_table(&output).unwrap();
    assert_eq!(table.num_rows(), 1);
    assert_eq!(table.num_columns(), 26);
    assert!(matches!(
        table.column_by_name("modifications").unwrap().data_type(),
        DataType::List(_)
    ));
    assert!(matches!(
        table
            .column_by_name("protein_accessions")
            .unwrap()
            .data_type(),
        DataType::List(_)
    ));
    assert_eq!(
        table.column_by_name("mass_error_ppm").unwrap().data_type(),
        &DataType::Float32
    );
    let schema = table.schema();
    let metadata = schema.metadata();
    assert_eq!(metadata.get("qpx_version").unwrap(), "1.1");
    assert_eq!(metadata.get("creator").unwrap(), "OpenMS");
    assert_eq!(metadata.get("file_type").unwrap(), "psm_file");
    assert_eq!(metadata.get("compression_format").unwrap(), "zstd");
    assert_eq!(metadata.get("scan_format").unwrap(), "scan");
    assert!(
        metadata
            .get("software_provider")
            .unwrap()
            .starts_with("OpenMS ")
    );
    assert!(metadata.contains_key("creation_date"));
    assert!(metadata.contains_key("uuid"));
}

/// `[EXTRA] QPX scan_format is derived from the spectrum native IDs`.
#[test]
fn class_test_scan_format() {
    let scan_format_for = |reference: &str| -> String {
        let proteins = vec![run("sf_search", "", true, &["scan_format_run.mzML"])];
        let mut identification = peptide("sf_search", "TestScore", true);
        if !reference.is_empty() {
            identification.set_spectrum_reference(reference);
        }
        identification.hits = vec![hit("PEPTIDER", 2, 1.0)];
        let output = temp(&format!("scan_format_{}.parquet", reference.len()));
        if QpxFile::export_to_parquet(
            &proteins,
            &[identification],
            &output,
            false,
            &ParquetWriteConfig::default(),
            None,
        )
        .is_err()
        {
            return "<write failed>".into();
        }
        let table = ParquetFile::read_table(&output).unwrap();
        let schema = table.schema();
        schema
            .metadata()
            .get("scan_format")
            .cloned()
            .unwrap_or_default()
    };
    assert_eq!(
        scan_format_for("controllerType=0 controllerNumber=1 scan=1234"),
        "scan"
    );
    assert_eq!(scan_format_for("scan=42"), "scan");
    assert_eq!(scan_format_for("index=7"), "index");
    assert_eq!(scan_format_for("some_opaque_identifier"), "<write failed>");
    assert_eq!(scan_format_for(""), "<write failed>");
}

/// `Test modifications structured output`.
#[test]
fn class_test_modifications() {
    let proteins = vec![run("test_search", "", true, &[])];
    let export = |sequence: &str| {
        let mut identification = peptide("test_search", "TestScore", true);
        identification.rt = Some(100.0);
        identification.mz = Some(500.0);
        identification.hits = vec![hit(sequence, 2, 0.9)];
        let table = QpxFile::export_psms_to_qpx_arrow(&proteins, &[identification], false, None)
            .unwrap()
            .value;
        let modifications = table.column_by_name("modifications").unwrap();
        let list = modifications.as_list::<i32>();
        list.value_length(0)
    };
    assert!(export("PEM(Oxidation)TIDER") > 0);
    assert_eq!(export("PEPTIDER"), 0);
}

/// `[EXTRA] the three QPX views agree on their run names`, psm view only:
/// the feature and pg views are `ConsensusMapArrowExport` and
/// `ProteinGroupArrowExport`, not ported.
#[test]
fn class_test_run_names_of_the_psm_view() {
    let paths = ["/data/frac_a/RUN_ONE.mzML", "/data/frac_b/RUN_TWO.mzML"];
    let proteins = vec![run("merged", "q-value", false, &paths)];
    let mut peptides = Vec::new();
    for m in 0..2i64 {
        let mut identification = peptide("merged", "q-value", false);
        identification.set_spectrum_reference(format!(
            "controllerType=0 controllerNumber=1 scan={}",
            10 + m
        ));
        identification
            .metadata
            .insert("id_merge_index".into(), m.into());
        identification.hits = vec![hit("PEPTIDEK", 2, 0.01)];
        peptides.push(identification);
    }
    let table = QpxFile::export_psms_to_qpx_arrow(&proteins, &peptides, false, None)
        .unwrap()
        .value;
    let names = table.column_by_name("run_file_name").unwrap();
    let names = names.as_string::<i32>();
    let set: std::collections::BTreeSet<&str> = (0..names.len()).map(|r| names.value(r)).collect();
    assert_eq!(set, ["RUN_ONE", "RUN_TWO"].into_iter().collect());
}

/// `[EXTRA] importFromArrow_round_trip`.
#[test]
fn class_test_import_round_trip() {
    let proteins = vec![run("run_1", "score", true, &[])];
    let mut first = peptide("run_1", "score", true);
    first.rt = Some(123.4);
    first.mz = Some(567.89);
    first.set_spectrum_reference("scan=42");
    let mut h = hit("PEPTIDE", 2, 0.95);
    h.metadata.insert("target_decoy".into(), "target".into());
    h.metadata.insert("COMET:deltaCn".into(), float(0.5));
    h.evidences = vec![PeptideEvidence {
        protein_accession: "sp|P12345|EXAMPLE".into(),
        aa_before: FlankingResidue::Residue('K'),
        aa_after: FlankingResidue::Residue('R'),
        start: Some(42),
        end: Some(48),
    }];
    first.hits = vec![h];
    let mut second = peptide("run_1", "score", true);
    second.rt = Some(234.5);
    second.mz = Some(678.90);
    second.set_spectrum_reference("scan=99");
    let mut h2 = hit("ACDEFGHIK", 3, 0.80);
    h2.metadata.insert("target_decoy".into(), "target".into());
    second.hits = vec![h2];

    let table = QpxFile::export_to_arrow(&proteins, &[first, second], true)
        .unwrap()
        .value;
    assert_eq!(table.num_rows(), 2);
    let mut out_proteins = proteins.clone();
    let mut out = Vec::new();
    QpxFile::import_from_arrow(&table, &mut out_proteins, &mut out).unwrap();
    assert_eq!(out.len(), 2);
    let hit0 = &out[0].hits[0];
    assert_eq!(hit0.sequence.to_string(), "PEPTIDE");
    assert_eq!(hit0.charge, 2);
    assert_eq!(hit0.score, 0.95);
    assert_eq!(
        hit0.metadata.get("target_decoy").unwrap().to_string(),
        "target"
    );
    assert_eq!(
        hit0.metadata
            .get("COMET:deltaCn")
            .unwrap()
            .as_f64()
            .unwrap(),
        0.5
    );
    assert_eq!(out[0].rt, Some(123.4));
    assert_eq!(out[0].mz, Some(567.89));
    assert_eq!(out[0].identifier, "run_1");
    assert_eq!(out[0].spectrum_reference(), "scan=42");
    assert_eq!(
        hit0.evidences,
        [PeptideEvidence {
            protein_accession: "sp|P12345|EXAMPLE".into(),
            aa_before: FlankingResidue::Residue('K'),
            aa_after: FlankingResidue::Residue('R'),
            start: Some(42),
            end: Some(48),
        }]
    );
    assert_eq!(out[1].hits[0].sequence.to_string(), "ACDEFGHIK");
    assert_eq!(out[1].hits[0].charge, 3);
    assert_eq!(out[1].hits[0].score, 0.80);
    assert_eq!(out[1].rt, Some(234.5));
    assert_eq!(out[1].mz, Some(678.90));
    assert_eq!(out[1].spectrum_reference(), "scan=99");
}

/// `[EXTRA] importFromArrow_round_trip_peptide_evidence_unknown_sentinels`.
#[test]
fn class_test_evidence_sentinels_round_trip() {
    let proteins = vec![run("run_1", "score", true, &[])];
    let mut identification = peptide("run_1", "score", true);
    identification.rt = Some(11.0);
    identification.mz = Some(222.0);
    identification.set_spectrum_reference("scan=1");
    let mut h = hit("PEPTIDE", 2, 0.5);
    h.metadata.insert("target_decoy".into(), "target".into());
    let unknown = PeptideEvidence {
        protein_accession: "PROT_X".into(),
        ..PeptideEvidence::default()
    };
    let terminal = PeptideEvidence {
        protein_accession: "PROT_TERM".into(),
        aa_before: FlankingResidue::NTerminus,
        aa_after: FlankingResidue::CTerminus,
        start: Some(0),
        end: Some(6),
    };
    h.evidences = vec![unknown.clone(), terminal.clone()];
    identification.hits = vec![h];
    let table = QpxFile::export_to_arrow(&proteins, &[identification], true)
        .unwrap()
        .value;
    let mut out_proteins = proteins.clone();
    let mut out = Vec::new();
    QpxFile::import_from_arrow(&table, &mut out_proteins, &mut out).unwrap();
    assert_eq!(out[0].hits[0].evidences, [unknown, terminal]);
}

/// `[EXTRA] importFromArrow_appends_shell_for_unknown_run_identifier`.
#[test]
fn class_test_import_appends_shell_run() {
    let mut identification = peptide("run_unknown", "hyperscore", true);
    identification.rt = Some(99.0);
    identification.mz = Some(400.0);
    identification.set_spectrum_reference("scan=7");
    let mut h = hit("TESTVAK", 2, 0.70);
    h.metadata.insert("target_decoy".into(), "target".into());
    identification.hits = vec![h];
    let table = QpxFile::export_to_arrow(&[], &[identification], true)
        .unwrap()
        .value;
    let mut proteins = Vec::new();
    let mut out = Vec::new();
    QpxFile::import_from_arrow(&table, &mut proteins, &mut out).unwrap();
    assert_eq!(proteins.len(), 1);
    assert_eq!(proteins[0].identifier, "run_unknown");
    assert_eq!(proteins[0].score_type, "hyperscore");
    assert!(proteins[0].higher_score_better);
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].identifier, "run_unknown");
    assert_eq!(out[0].hits[0].sequence.to_string(), "TESTVAK");
}

fn streaming_data(
    m: usize,
    modified: &str,
) -> (Vec<ProteinIdentification>, Vec<PeptideIdentification>) {
    let proteins = vec![run("test_search", "TestScore", true, &["/data/run1.mzML"])];
    let peptides = (0..m)
        .map(|i| {
            let mut identification = peptide("test_search", "TestScore", true);
            identification.rt = Some(100.0 + i as f64);
            identification.mz = Some(400.0 + (i % 500) as f64);
            identification.set_spectrum_reference(format!(
                "controllerType=0 controllerNumber=1 scan={}",
                1000 + i
            ));
            let mut h = hit(
                if i % 3 == 0 { modified } else { "PEPTIDEK" },
                2 + (i % 3) as i32,
                0.99 - (i % 100) as f64 * 0.001,
            );
            h.metadata.insert(
                "target_decoy".into(),
                if i % 2 == 0 { "target" } else { "decoy" }.into(),
            );
            h.evidences = vec![PeptideEvidence {
                protein_accession: format!("PROT_{}", i % 50),
                ..PeptideEvidence::default()
            }];
            identification.hits = vec![h];
            identification
        })
        .collect();
    (proteins, peptides)
}

/// The columns the class test compares, in row order.
fn stable_columns(table: &RecordBatch) -> Vec<String> {
    (0..table.num_rows())
        .map(|row| {
            [
                "sequence",
                "peptidoform",
                "charge",
                "rt",
                "posterior_error_probability",
            ]
            .iter()
            .map(|name| cell(table.column_by_name(name).unwrap().as_ref(), row))
            .collect::<Vec<_>>()
            .join("|")
        })
        .collect()
}

/// `static bool exportToParquetStreaming(...)`.
#[test]
fn class_test_streaming() {
    let (proteins, peptides) = streaming_data(2500, "PEM(Oxidation)TIDER");
    let borrowed: Vec<&PeptideIdentification> = peptides.iter().collect();
    let stream = temp("stream_class.parquet");
    let options = StreamingOptions {
        batch_size: 1000,
        ..StreamingOptions::default()
    };
    QpxFile::export_to_parquet_streaming(&proteins, &borrowed, &stream, &options, None).unwrap();
    let streamed = ParquetFile::read_table(&stream).unwrap();
    assert_eq!(streamed.num_rows(), 2500);
    assert_eq!(streamed.num_columns(), 26);
    let schema = streamed.schema();
    assert_eq!(schema.metadata().get("file_type").unwrap(), "psm_file");
    assert_eq!(schema.metadata().get("qpx_version").unwrap(), "1.1");

    let reference = temp("stream_ref.parquet");
    QpxFile::export_to_parquet(
        &proteins,
        &peptides,
        &reference,
        false,
        &ParquetWriteConfig::default(),
        None,
    )
    .unwrap();
    let reference = ParquetFile::read_table(&reference).unwrap();
    assert_eq!(stable_columns(&streamed), stable_columns(&reference));

    let empty = temp("stream_empty.parquet");
    QpxFile::export_to_parquet_streaming(
        &proteins,
        &[],
        &empty,
        &StreamingOptions::default(),
        None,
    )
    .unwrap();
    let empty = ParquetFile::read_table(&empty).unwrap();
    assert_eq!((empty.num_rows(), empty.num_columns()), (0, 26));

    let one = temp("stream_one.parquet");
    let options = StreamingOptions {
        batch_size: 1,
        ..StreamingOptions::default()
    };
    QpxFile::export_to_parquet_streaming(&proteins, &borrowed[..1], &one, &options, None).unwrap();
    assert_eq!(ParquetFile::read_table(&one).unwrap().num_rows(), 1);

    // A key repeated after a batch boundary is refused, and nothing is left.
    let duplicate = temp("stream_duplicate.parquet");
    assert!(
        QpxFile::export_to_parquet_streaming(
            &proteins,
            &[borrowed[0], borrowed[0]],
            &duplicate,
            &options,
            None
        )
        .is_err()
    );
    assert!(!duplicate.exists());

    let zero = temp("stream_zero.parquet");
    let options = StreamingOptions {
        batch_size: 0,
        ..StreamingOptions::default()
    };
    QpxFile::export_to_parquet_streaming(&proteins, &borrowed, &zero, &options, None).unwrap();
    assert_eq!(ParquetFile::read_table(&zero).unwrap().num_rows(), 2500);
}

/// `[EXTRA] exportToParquetStreaming parallel build (n_threads)`, and the
/// port's stronger contract: the file is byte-identical for every thread
/// count once its dynamic metadata is equal.
#[test]
fn class_test_streaming_threads() {
    let (proteins, peptides) = streaming_data(2500, "PEM(Oxidation)TIDEK");
    let borrowed: Vec<&PeptideIdentification> = peptides.iter().collect();
    let reference = temp("threads_ref.parquet");
    QpxFile::export_to_parquet(
        &proteins,
        &peptides,
        &reference,
        false,
        &ParquetWriteConfig::default(),
        None,
    )
    .unwrap();
    let reference = stable_columns(&ParquetFile::read_table(&reference).unwrap());
    let mut layouts = Vec::new();
    for threads in [1, 2, 8] {
        let file = temp(&format!("threads_{threads}.parquet"));
        let options = StreamingOptions {
            batch_size: 1000,
            n_threads: threads,
            ..StreamingOptions::default()
        };
        QpxFile::export_to_parquet_streaming(&proteins, &borrowed, &file, &options, None).unwrap();
        let table = ParquetFile::read_table(&file).unwrap();
        assert_eq!(table.num_rows(), 2500);
        assert_eq!(stable_columns(&table), reference);
        layouts.push(parquet_layout(&file).0);
    }
    assert!(layouts.windows(2).all(|pair| pair[0] == pair[1]));

    // export_all_psms with 1..3 hits per identification
    let mut multi = Vec::new();
    let mut expected_rows = 0;
    for i in 0..777 {
        let mut identification = peptide("test_search", "TestScore", true);
        identification.rt = Some(10.0 + i as f64);
        identification.mz = Some(300.0 + (i % 100) as f64);
        identification
            .set_spectrum_reference(format!("controllerType=0 controllerNumber=1 scan={i}"));
        let hits = 1 + i % 3;
        for h in 0..hits {
            identification.hits.push(hit(
                if h == 0 { "PEPTIDEK" } else { "DFPIANGER" },
                2 + h as i32,
                0.9 - 0.01 * h as f64,
            ));
        }
        expected_rows += hits;
        multi.push(identification);
    }
    let all_reference = temp("threads_all_ref.parquet");
    QpxFile::export_to_parquet(
        &proteins,
        &multi,
        &all_reference,
        true,
        &ParquetWriteConfig::default(),
        None,
    )
    .unwrap();
    let all_reference = ParquetFile::read_table(&all_reference).unwrap();
    assert_eq!(all_reference.num_rows(), expected_rows);
    let all = temp("threads_all.parquet");
    let borrowed_multi: Vec<&PeptideIdentification> = multi.iter().collect();
    let options = StreamingOptions {
        export_all_psms: true,
        batch_size: 50,
        n_threads: 8,
        ..StreamingOptions::default()
    };
    QpxFile::export_to_parquet_streaming(&proteins, &borrowed_multi, &all, &options, None).unwrap();
    let all = ParquetFile::read_table(&all).unwrap();
    assert_eq!(stable_columns(&all), stable_columns(&all_reference));

    // No hits at all, with 8 threads
    let empty_hits: Vec<PeptideIdentification> =
        (0..50).map(|_| PeptideIdentification::default()).collect();
    let borrowed_empty: Vec<&PeptideIdentification> = empty_hits.iter().collect();
    let file = temp("threads_no_hits.parquet");
    let options = StreamingOptions {
        batch_size: 10,
        n_threads: 8,
        ..StreamingOptions::default()
    };
    QpxFile::export_to_parquet_streaming(&proteins, &borrowed_empty, &file, &options, None)
        .unwrap();
    let table = ParquetFile::read_table(&file).unwrap();
    assert_eq!((table.num_rows(), table.num_columns()), (0, 26));

    for (name, slice) in [
        ("t0", &borrowed[..0]),
        ("t1", &borrowed[..1]),
        ("t3", &borrowed[..3]),
    ] {
        let file = temp(&format!("threads_{name}.parquet"));
        let options = StreamingOptions {
            n_threads: 8,
            ..StreamingOptions::default()
        };
        QpxFile::export_to_parquet_streaming(&proteins, slice, &file, &options, None).unwrap();
        assert_eq!(
            ParquetFile::read_table(&file).unwrap().num_rows(),
            slice.len()
        );
    }
}

/// `[EXTRA] exportToParquetStreaming dedicated-column values (index path + precedence)`.
#[test]
fn class_test_dedicated_columns() {
    let proteins = vec![run("run1", "q-value", false, &["/data/fallback.mzML"])];
    let mut template = peptide("run1", "q-value", false);
    template.rt = Some(42.0);
    template.mz = Some(555.5);
    let mut h = hit("PEPTIDEK", 3, 0.001);
    h.metadata.insert("target_decoy".into(), "decoy".into());
    h.metadata
        .insert("Posterior Error Probability".into(), float(0.25));
    h.metadata.insert("predicted_RT".into(), float(40.0));
    h.metadata.insert("predicted_rt".into(), float(999.0));
    h.metadata.insert("ion_mobility".into(), float(1.23));
    h.metadata.insert("IM".into(), float(9.99));
    h.metadata.insert("missed_cleavages".into(), 2i64.into());
    h.metadata
        .insert("reference_file_name".into(), "/data/real.mzML".into());
    template.hits = vec![h];
    let peptides: Vec<PeptideIdentification> = (0..10)
        .map(|i| {
            let mut replicate = template.clone();
            replicate.set_spectrum_reference(format!(
                "controllerType=0 controllerNumber=1 scan={}",
                4242 + i
            ));
            replicate
        })
        .collect();
    let borrowed: Vec<&PeptideIdentification> = peptides.iter().collect();
    let streamed = temp("dedicated_stream.parquet");
    let options = StreamingOptions {
        batch_size: 3,
        n_threads: 4,
        ..StreamingOptions::default()
    };
    QpxFile::export_to_parquet_streaming(&proteins, &borrowed, &streamed, &options, None).unwrap();
    let oneshot = temp("dedicated_oneshot.parquet");
    QpxFile::export_to_parquet(
        &proteins,
        &peptides,
        &oneshot,
        false,
        &ParquetWriteConfig::default(),
        None,
    )
    .unwrap();
    for file in [streamed, oneshot] {
        let table = ParquetFile::read_table(&file).unwrap();
        assert_eq!(table.num_rows(), 10);
        for row in [0, 9] {
            let get = |name: &str| cell(table.column_by_name(name).unwrap().as_ref(), row);
            assert_eq!(get("is_decoy"), "b:true");
            assert_eq!(
                get("predicted_rt"),
                format!("f32:{:08x}", 40.0f32.to_bits())
            );
            assert_eq!(
                get("ion_mobility"),
                format!("f32:{:08x}", 1.23f32.to_bits())
            );
            assert_eq!(get("missed_cleavages"), "i:2");
            assert_eq!(get("run_file_name"), "\"real\"");
            assert_eq!(get("posterior_error_probability"), bits(0.25));
        }
    }
}

fn merged_psm(sequence: &str, merge_index: Option<i64>) -> PeptideIdentification {
    let mut identification = peptide("merged", "q-value", false);
    identification.rt = Some(42.0);
    identification.mz = Some(555.5);
    identification.set_spectrum_reference("controllerType=0 controllerNumber=1 scan=17");
    if let Some(index) = merge_index {
        identification
            .metadata
            .insert("id_merge_index".into(), index.into());
    }
    identification.hits = vec![hit(sequence, 2, 0.001)];
    identification
}

fn run_names(path: &Path) -> Vec<String> {
    let table = ParquetFile::read_table(path).unwrap();
    let names = table.column_by_name("run_file_name").unwrap();
    let names = names.as_string::<i32>();
    (0..names.len())
        .map(|r| names.value(r).to_owned())
        .collect()
}

/// `[EXTRA] run_file_name resolves per PSM on a merged run`.
#[test]
fn class_test_merged_run_names() {
    let proteins = vec![run(
        "merged",
        "q-value",
        false,
        &["/data/runA.mzML", "/data/runB.mzML.gz"],
    )];
    let peptides = vec![
        merged_psm("PEPTIDEK", Some(0)),
        merged_psm("TESTPEPTIDER", Some(1)),
        merged_psm("ANOTHERPEPK", Some(1)),
    ];
    let borrowed: Vec<&PeptideIdentification> = peptides.iter().collect();
    let oneshot = temp("merged_oneshot.parquet");
    QpxFile::export_to_parquet(
        &proteins,
        &peptides,
        &oneshot,
        false,
        &ParquetWriteConfig::default(),
        None,
    )
    .unwrap();
    let streamed = temp("merged_stream.parquet");
    let options = StreamingOptions {
        batch_size: 2,
        n_threads: 4,
        ..StreamingOptions::default()
    };
    QpxFile::export_to_parquet_streaming(&proteins, &borrowed, &streamed, &options, None).unwrap();
    for file in [oneshot, streamed] {
        assert_eq!(run_names(&file), ["runA", "runB", "runB"]);
    }
}

/// `[EXTRA] duplicate identification-run identifiers are refused by both merge-index preflights`.
#[test]
fn class_test_duplicate_run_identifiers() {
    let mut proteins = vec![run("merged", "q-value", false, &["/data/runA.mzML"])];
    let mut duplicate = proteins[0].clone();
    duplicate.primary_ms_run_paths = vec!["/data/runB.mzML".into()];
    proteins.push(duplicate);
    let peptides = vec![merged_psm("PEPTIDEK", Some(0))];
    assert!(matches!(
        QpxFile::require_resolvable_merge_indices(&proteins, &peptides),
        Err(Error::InvalidValue(_))
    ));
    let borrowed: Vec<&PeptideIdentification> = peptides.iter().collect();
    assert!(matches!(
        QpxFile::require_resolvable_merge_indices_of(&proteins, &borrowed),
        Err(Error::InvalidValue(_))
    ));
}

/// `[EXTRA] merged run without a usable id_merge_index is refused`.
#[test]
fn class_test_unusable_merge_index() {
    let paths = ["/data/runA.mzML", "/data/runB.mzML"];
    let proteins = vec![run("merged", "q-value", false, &paths)];
    let config = ParquetWriteConfig::default();
    let options = StreamingOptions {
        batch_size: 2,
        n_threads: 4,
        ..StreamingOptions::default()
    };
    let missing =
        |result: openms::Result<Vec<String>>| matches!(result, Err(Error::MissingInformation(_)));

    let peptides = vec![
        merged_psm("PEPTIDEK", Some(0)),
        merged_psm("TESTPEPTIDER", None),
    ];
    let borrowed: Vec<&PeptideIdentification> = peptides.iter().collect();
    assert!(matches!(
        QpxFile::require_resolvable_merge_indices_of(&proteins, &borrowed),
        Err(Error::MissingInformation(_))
    ));
    assert!(missing(QpxFile::export_to_parquet(
        &proteins,
        &peptides,
        temp("mi_a1.parquet"),
        false,
        &config,
        None
    )));
    let streamed = temp("mi_a2.parquet");
    assert!(missing(QpxFile::export_to_parquet_streaming(
        &proteins, &borrowed, &streamed, &options, None
    )));
    assert!(!streamed.exists());

    for (name, index) in [
        ("mi_b", MetaValue::from(2i64)),
        ("mi_c", MetaValue::from(-1i64)),
        ("mi_d", MetaValue::from("not_an_index")),
    ] {
        let mut identification = merged_psm("PEPTIDEK", Some(0));
        identification
            .metadata
            .insert("id_merge_index".into(), index);
        let peptides = vec![identification];
        assert!(missing(QpxFile::export_to_parquet(
            &proteins,
            &peptides,
            temp(&format!("{name}.parquet")),
            false,
            &config,
            None
        )));
        let borrowed: Vec<&PeptideIdentification> = peptides.iter().collect();
        assert!(missing(QpxFile::export_to_parquet_streaming(
            &proteins,
            &borrowed,
            temp(&format!("{name}_s.parquet")),
            &options,
            None
        )));
    }

    let single = vec![run("merged", "q-value", false, &["/data/only.mzML"])];
    let file = temp("mi_e.parquet");
    QpxFile::export_to_parquet(
        &single,
        &[merged_psm("PEPTIDEK", None)],
        &file,
        false,
        &config,
        None,
    )
    .unwrap();
    assert_eq!(run_names(&file), ["only"]);
}

/// `[EXTRA] run_file_name stem collisions warn but do not fail the export`.
#[test]
fn class_test_stem_collisions() {
    let proteins = vec![run(
        "merged",
        "q-value",
        false,
        &["/a/run1.mzML", "/b/run1.mzML"],
    )];
    let peptides = vec![
        merged_psm("PEPTIDEK", Some(0)),
        merged_psm("TESTPEPTIDER", Some(1)),
    ];
    let file = temp("collision.parquet");
    let warnings = QpxFile::export_to_parquet(
        &proteins,
        &peptides,
        &file,
        false,
        &ParquetWriteConfig::default(),
        None,
    )
    .unwrap();
    assert_eq!(run_names(&file), ["run1", "run1"]);
    assert_eq!(
        warnings,
        [
            "QPXFile: several MS runs share the run_file_name 'run1' after stripping path and \
          extension: /a/run1.mzML, /b/run1.mzML. They cannot be told apart in the exported QPX \
          tables."
        ]
    );
}

/// `[EXTRA] per-PSM origin survives an internal-format round trip`.
#[test]
fn class_test_origin_survives_round_trip() {
    let proteins = vec![run(
        "merged",
        "q-value",
        false,
        &["/data/runA.mzML", "/data/runB.mzML"],
    )];
    let peptides = vec![
        merged_psm("PEPTIDEK", Some(0)),
        merged_psm("TESTPEPTIDER", Some(1)),
    ];
    let internal = QpxFile::export_to_arrow(&proteins, &peptides, false)
        .unwrap()
        .value;
    let reference = internal.column_by_name("reference_file_name").unwrap();
    let reference = reference.as_string::<i32>();
    assert_eq!(
        (reference.value(0), reference.value(1)),
        ("/data/runA.mzML", "/data/runB.mzML")
    );
    let mut proteins_in = Vec::new();
    let mut peptides_in = Vec::new();
    QpxFile::import_from_arrow(&internal, &mut proteins_in, &mut peptides_in).unwrap();
    let qpx = QpxFile::export_psms_to_qpx_arrow(&proteins_in, &peptides_in, false, None)
        .unwrap()
        .value;
    assert_eq!(qpx.num_rows(), 2);
    let names = qpx.column_by_name("run_file_name").unwrap();
    let names = names.as_string::<i32>();
    assert_eq!((names.value(0), names.value(1)), ("runA", "runB"));
}

// ---------------------------------------------------------------------------
// Native checks
// ---------------------------------------------------------------------------

/// The import refuses a table of the wrong layout, and leaves its outputs
/// untouched on error.
#[test]
fn import_refuses_incompatible_tables_atomically() {
    let qpx = QpxFile::export_psms_to_qpx_arrow(&runs(), &peptides(false), false, None)
        .unwrap()
        .value;
    let mut proteins = runs();
    let mut out = vec![PeptideIdentification::default()];
    match QpxFile::import_from_arrow(&qpx, &mut proteins, &mut out) {
        Err(Error::InvalidValue(message)) => {
            assert!(message.starts_with("QPXFile::importFromArrow: Incompatible PSM schema: "))
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(out.len(), 1);
    assert_eq!(proteins.len(), 3);

    // a sequence that cannot be parsed (and no peptidoform to prefer)
    let internal = QpxFile::export_to_arrow(&runs(), &peptides(false), false)
        .unwrap()
        .value;
    let index = internal.schema().index_of("sequence").unwrap();
    let mut columns = internal.columns().to_vec();
    columns[index] = std::sync::Arc::new(arrow_array::StringArray::from(vec![
        Some("PEP[[TIDE"),
        None,
        None,
    ]));
    let bad = RecordBatch::try_new(internal.schema(), columns).unwrap();
    let index = bad.schema().index_of("peptidoform").unwrap();
    let mut columns = bad.columns().to_vec();
    columns[index] = arrow_array::new_null_array(&DataType::Utf8, 3);
    let bad = RecordBatch::try_new(bad.schema(), columns).unwrap();
    let before = out.clone();
    assert!(QpxFile::import_from_arrow(&bad, &mut proteins, &mut out).is_err());
    assert_eq!(out, before);
}

/// Every metadata value the numeric columns read must be a number: the
/// source reinterprets a string's storage as a double.
#[test]
fn non_numeric_column_metadata_is_refused() {
    let mut peptides = peptides(false);
    peptides[0].hits[0]
        .metadata
        .insert("predicted_RT".into(), "soon".into());
    match QpxFile::export_to_arrow(&runs(), &peptides, false) {
        Err(Error::InvalidValue(message)) => assert_eq!(
            message,
            "QPXFile: metadata value 'predicted_RT' of type 'String' and value 'soon' is not numeric"
        ),
        other => panic!("{other:?}"),
    }
}

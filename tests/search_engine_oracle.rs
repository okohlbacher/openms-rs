// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Differential tests of the search-engine readers against the executed C++
//! Release build (tier 1).
//!
//! `tests/data/search_engine_oracle/*.out` were written by
//! `tests/data/search_engine_oracle/driver.cpp.txt`, built against the Linux
//! x86_64 Release install `openms4-release-bc9cc12-c19e494-174b576` and run on
//! ibminode06 (`host.txt`, `compiler.txt`, `hashes.txt`), twice, with
//! byte-identical output once the run date is masked. Each dump lists every
//! identification, hit, meta value and evidence the C++ reader produced from
//! the upstream class-test fixture; the Rust reader's result is dumped in the
//! same format and compared line by line.
//!
//! Normalisations, each a representational difference, not a tolerance:
//! the run date (`DateTime::now()` in the source, pinned here) is masked on
//! both sides; the C++ `significance_threshold` and `rank` meta values are
//! dropped, because OpenMS4 stores those two fields as meta values while the
//! Rust types hold them as fields, which the `id`/`hit` lines already print;
//! a C++ sequence is re-rendered through the Rust `AASequence` so notation
//! differences in `toString` do not hide or fake a difference; Sequest's `-`
//! terminus flanking residue (45) is the Rust terminus marker (`[` 91 before,
//! `]` 93 after); and PepNovo's `Parsed N ids` log line is dropped.
#![cfg(feature = "idxml")]

use openms::chemistry::{AASequence, ModificationDefinitionsSet};
use openms::comparison::Tolerance;
use openms::data_structures::DateTime;
use openms::format::omssa_xml_file::{OMSSALoadOptions, OMSSAXMLFile};
use openms::format::pepnovo_infile::PepNovoInfile;
use openms::format::pepnovo_outfile::{IndexPosMapping, PepNovoOutfile};
use openms::format::sequest_outfile::SequestOutfile;
use openms::format::xtandem_xml_file::XTandemXMLFile;
use openms::identification::{PeptideIdentification, ProteinIdentification};
use openms::metadata::{MetaInfo, MetaValueData};
use std::collections::BTreeMap;
use std::fmt::Write;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn oracle(name: &str) -> String {
    std::fs::read_to_string(data(&format!("search_engine_oracle/{name}"))).unwrap()
}

fn run_date() -> DateTime {
    "2026-09-27 12:00:00".parse().unwrap()
}

/// `printf("%.17g")`, the driver's number format.
fn g17(value: f64) -> String {
    if value.is_nan() {
        return "nan".into();
    }
    if value.is_infinite() {
        return if value > 0.0 { "inf" } else { "-inf" }.into();
    }
    if value == 0.0 {
        return if value.is_sign_negative() { "-0" } else { "0" }.into();
    }
    let scientific = format!("{value:.16e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap();
    let exponent: i32 = exponent.parse().unwrap();
    let strip = |text: &str| -> String {
        if text.contains('.') {
            text.trim_end_matches('0').trim_end_matches('.').to_owned()
        } else {
            text.to_owned()
        }
    };
    if !(-4..17).contains(&exponent) {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{}e{sign}{:02}", strip(mantissa), exponent.abs())
    } else {
        let decimals = (16 - exponent) as usize;
        strip(&format!("{value:.decimals$}"))
    }
}

fn tolerance(value: Tolerance) -> f64 {
    match value {
        Tolerance::Absolute(value) | Tolerance::Ppm(value) => value,
    }
}

fn dump_meta(out: &mut String, meta: &MetaInfo, indent: &str) {
    for (key, value) in meta {
        let (kind, text) = match value.data() {
            MetaValueData::String(text) => (0, text.clone()),
            MetaValueData::Integer(value) => (1, value.to_string()),
            MetaValueData::Float(value) => (2, g17(*value)),
            other => panic!("unexpected meta value {other:?}"),
        };
        writeln!(out, "{indent}meta {key} type={kind} value={text}").unwrap();
    }
}

fn dump_protein(out: &mut String, protein: &ProteinIdentification) {
    writeln!(
        out,
        "protein identifier={} engine={} version={} score_type={} higher={}",
        protein.identifier,
        protein.search_engine,
        protein.search_engine_version,
        protein.score_type,
        u8::from(protein.higher_score_better)
    )
    .unwrap();
    let parameters = &protein.search_parameters;
    write!(
        out,
        "  params fragment_tol={} precursor_tol={} variable=",
        g17(tolerance(parameters.fragment_tolerance)),
        g17(tolerance(parameters.precursor_tolerance))
    )
    .unwrap();
    for modification in &parameters.variable_modifications {
        write!(out, "[{modification}]").unwrap();
    }
    writeln!(out).unwrap();
    writeln!(out, "  protein_hits {}", protein.hits.len()).unwrap();
    for hit in &protein.hits {
        writeln!(
            out,
            "  protein_hit acc={} score={}",
            hit.accession,
            g17(hit.score)
        )
        .unwrap();
    }
}

fn dump_peptides(out: &mut String, ids: &[PeptideIdentification]) {
    writeln!(out, "peptide_ids {}", ids.len()).unwrap();
    for id in ids {
        writeln!(
            out,
            "id identifier={} score_type={} higher={} threshold={} rt={} mz={}",
            id.identifier,
            id.score_type,
            u8::from(id.higher_score_better),
            g17(id.significance_threshold),
            g17(id.rt.unwrap_or(f64::NAN)),
            g17(id.mz.unwrap_or(f64::NAN))
        )
        .unwrap();
        dump_meta(out, &id.metadata, "  ");
        writeln!(out, "  hits {}", id.hits.len()).unwrap();
        for hit in &id.hits {
            writeln!(
                out,
                "  hit seq={} score={} rank={} charge={}",
                hit.sequence,
                g17(hit.score),
                hit.rank,
                hit.charge
            )
            .unwrap();
            dump_meta(out, &hit.metadata, "    ");
            for evidence in &hit.evidences {
                let position = |value: Option<usize>| value.map_or(-1, |v| v as i64);
                writeln!(
                    out,
                    "    evidence acc={} start={} end={} before={} after={}",
                    evidence.protein_accession,
                    position(evidence.start),
                    position(evidence.end),
                    u32::from(evidence.aa_before.code()),
                    u32::from(evidence.aa_after.code())
                )
                .unwrap();
            }
        }
    }
}

/// Replace every `YYYY-MM-DD` date, with an optional ` hh:mm:ss` time, by
/// `<DATE>`.
fn mask_dates(text: &str) -> String {
    let bytes = text.as_bytes();
    let digits = |from: usize, count: usize| {
        bytes
            .get(from..from + count)
            .is_some_and(|slice| slice.iter().all(u8::is_ascii_digit))
    };
    let mut out = String::with_capacity(text.len());
    let mut i = 0;
    while i < bytes.len() {
        let date = digits(i, 4)
            && bytes.get(i + 4) == Some(&b'-')
            && digits(i + 5, 2)
            && bytes.get(i + 7) == Some(&b'-')
            && digits(i + 8, 2);
        if date {
            let mut end = i + 10;
            let time = bytes.get(end) == Some(&b' ')
                && digits(end + 1, 2)
                && bytes.get(end + 3) == Some(&b':')
                && digits(end + 4, 2)
                && bytes.get(end + 6) == Some(&b':')
                && digits(end + 7, 2);
            if time {
                end += 9;
            }
            out.push_str("<DATE>");
            i = end;
        } else {
            let c = text[i..].chars().next().unwrap();
            out.push(c);
            i += c.len_utf8();
        }
    }
    out
}

/// The normalisations listed in the module documentation.
fn normalise_cpp(text: &str) -> Vec<String> {
    mask_dates(text)
        .lines()
        .filter(|line| {
            !line.starts_with("  meta significance_threshold ")
                && !line.starts_with("    meta rank ")
                && !line.starts_with("Parsed ")
        })
        .map(|line| {
            if let Some(rest) = line.strip_prefix("  hit seq=") {
                let (sequence, tail) = rest.split_once(" score=").unwrap();
                let rendered = AASequence::parse(sequence)
                    .map(|parsed| parsed.to_string())
                    .unwrap_or_else(|_| sequence.to_owned());
                return format!("  hit seq={rendered} score={tail}");
            }
            if line.starts_with("    evidence ") {
                return line
                    .replace("before=45", "before=91")
                    .replace("after=45", "after=93");
            }
            line.to_owned()
        })
        .collect()
}

fn compare(rust: &str, cpp: &str) {
    let rust: Vec<String> = mask_dates(rust).lines().map(str::to_owned).collect();
    let cpp = normalise_cpp(cpp);
    for (index, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        assert_eq!(r, c, "first difference at line {}", index + 1);
    }
    assert_eq!(rust.len(), cpp.len(), "line counts differ");
}

#[test]
fn xtandem_xml_matches_the_release_build() {
    let mods = ModificationDefinitionsSet::from_comma_separated(
        "",
        "Carbamidomethyl (C),Oxidation (M),Carboxymethyl (C)",
    )
    .unwrap();
    for (fixture, expected) in [
        ("XTandemXMLFile_test.xml", "xtandem_xml_1.out"),
        ("XTandemXMLFile_test_2.xml", "xtandem_xml_2.out"),
    ] {
        let reader = XTandemXMLFile {
            run_date: Some(run_date()),
            ..XTandemXMLFile::new()
        };
        let result = reader.load(data(fixture), &mods).unwrap();
        let mut out = String::new();
        dump_protein(&mut out, &result.protein_identification);
        dump_peptides(&mut out, &result.peptide_identifications);
        write!(out, "fixed").unwrap();
        for name in result.modification_definitions.fixed_names() {
            write!(out, " [{name}]").unwrap();
        }
        write!(out, "\nvariable").unwrap();
        for name in result.modification_definitions.variable_names() {
            write!(out, " [{name}]").unwrap();
        }
        writeln!(out).unwrap();
        compare(&out, &oracle(expected));
    }
}

#[test]
fn omssa_xml_matches_the_release_build() {
    let mut reader = OMSSAXMLFile::new().unwrap();
    reader.run_date = Some(run_date());
    let result = reader
        .load(data("OMSSAXMLFile_test_1.xml"), OMSSALoadOptions::default())
        .unwrap();
    let mut out = String::new();
    dump_protein(&mut out, &result.protein_identification);
    dump_peptides(&mut out, &result.peptide_identifications);
    compare(&out, &oracle("omssa_xml.out"));
}

#[test]
fn sequest_outfile_matches_the_release_build() {
    let mut protein = ProteinIdentification::default();
    let id = SequestOutfile::new()
        .load(data("SequestOutfile.out"), &mut protein, 1.0, &[], false)
        .unwrap();
    let mut out = String::new();
    dump_protein(&mut out, &protein);
    dump_peptides(&mut out, id.as_slice());
    compare(&out, &oracle("sequest_out.out"));
}

#[test]
fn pepnovo_outfile_matches_the_release_build() {
    let keys: BTreeMap<String, String> = [("K+42", "Acetyl (K)"), ("Y+42", "Acetyl (Y)")]
        .into_iter()
        .map(|(key, id)| (key.to_owned(), id.to_owned()))
        .collect();
    let mut index = IndexPosMapping::new();
    index.insert(0, (1510.5732421875, 747.761901855469));
    index.insert(1, (1530.11535644531, 549.856262207031));
    index.insert(2, (1533.16589355469, 358.174530029297));
    index.insert(3, (1111.0, 2222.0));
    let reader = PepNovoOutfile {
        run_date: Some(run_date()),
    };
    for (mapping, expected) in [
        (index, "pepnovo_out_index.out"),
        (IndexPosMapping::new(), "pepnovo_out_title.out"),
    ] {
        let result = reader
            .load(data("PepNovoOutfile.out"), -4.0, &mapping, &keys)
            .unwrap();
        let mut out = String::new();
        dump_protein(&mut out, &result.protein_identification);
        dump_peptides(&mut out, &result.peptide_identifications);
        compare(&out, &oracle(expected));
    }
}

#[test]
fn pepnovo_infile_matches_the_release_build() {
    let mut infile = PepNovoInfile::new();
    infile
        .set_modifications(&["Phospho (C)"], &["Phospho (D)", "Ethanolamine (C-term)"])
        .unwrap();
    let mut written = Vec::new();
    infile.write(&mut written).unwrap();
    let written = String::from_utf8(written).unwrap();
    let mut out = format!("file_bytes {}\n{written}\n<<EOF\n", written.len());
    for (key, id) in infile.modifications() {
        writeln!(out, "key {key} -> {id}").unwrap();
    }
    assert_eq!(out, oracle("pepnovo_infile.out"));
}

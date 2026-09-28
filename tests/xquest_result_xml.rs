// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/xquest_result_xml.rs` (`FORMAT/XQuestResultXMLFile.h`,
//! `FORMAT/HANDLERS/XQuestResultXMLHandler.h`).
//!
//! Evidence, strongest first (see `docs/XQUEST_RESULT_XML_SUPPORT.md` and
//! `tests/data/xquest_result_xml_provenance.json`):
//!
//! - tier 1, executed differential: an oracle driver linked against the Linux
//!   x86_64 Release libOpenMS of openms4-release-bc9cc12-c19e494-174b576 ran on
//!   ibminode06, twice, identically (`../oracle/F3a-xquest`). For every input it
//!   recorded `load` (every field of both out-parameters and the three
//!   getters), `store` of that load, `IdXMLFile::load` + `store` for the
//!   upstream idXML files, and `splitByNth` / `splitByMiddle`;
//!   `tests/data/xquest_result_xml/expected/` holds the records;
//! - tier 3: the literals of `XQuestResultXMLFile_test.cpp`.
//!
//! The reader needs the `idxml` feature.
#![cfg(feature = "idxml")]

use openms::Error;
use openms::chemistry::{AASequence, SequenceModification};
use openms::format::file_info::text_format::{ostream_g, to_str};
use openms::format::idxml;
use openms::format::xquest_result_xml::{
    ReadLimits, XQuestDocument, XQuestResultXmlFile, XQuestResultXmlHandler,
};
use openms::identification::{EnzymeTermSpecificity, PeakMassType};
use openms::metadata::{MetaInfo, MetaValueData};
use std::fs;
use std::path::PathBuf;

const DATE: &str = "2000-01-01 00:00:00";

fn data() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/xquest_result_xml")
}

fn input(name: &str) -> PathBuf {
    data().join("inputs").join(name)
}

fn expected(name: &str) -> String {
    fs::read_to_string(data().join("expected").join(name)).unwrap()
}

fn g17(value: f64) -> String {
    ostream_g(value, 17)
}

/// The source's `AASequence::toString`: each annotation is `(` + the
/// modification's id + `)` (`ResidueModification::toString`).
fn source_to_string(sequence: &AASequence) -> String {
    let annotation = |modification: &SequenceModification| match modification.known() {
        Some(known) => format!("({})", known.name()),
        None => panic!("the xQuest reader sets named modifications only"),
    };
    let mut out = String::new();
    if let Some(modification) = sequence.n_terminal_modification() {
        out.push('.');
        out.push_str(&annotation(modification));
    }
    for (index, residue) in sequence.as_str().chars().enumerate() {
        out.push(residue);
        if let Some(modification) = sequence.residue_modification(index).unwrap() {
            out.push_str(&annotation(modification));
        }
    }
    if let Some(modification) = sequence.c_terminal_modification() {
        out.push('.');
        out.push_str(&annotation(modification));
    }
    out
}

/// The driver's `value`: doubles with `%.17g`, lists as it joins them.
fn render_meta(out: &mut Vec<String>, meta: &MetaInfo, indent: &str) {
    for (key, value) in meta {
        let (kind, text) = match value.data() {
            MetaValueData::Empty => ("empty", String::new()),
            MetaValueData::String(text) => ("string", text.clone()),
            MetaValueData::Integer(number) => ("int", number.to_string()),
            MetaValueData::Float(number) => ("double", g17(*number)),
            MetaValueData::StringList(list) => ("stringList", format!("[{}]", list.join("|"))),
            MetaValueData::IntegerList(list) => (
                "intList",
                format!(
                    "[{}]",
                    list.iter()
                        .map(i64::to_string)
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
            ),
            MetaValueData::FloatList(list) => (
                "doubleList",
                format!(
                    "[{}]",
                    list.iter().map(|v| g17(*v)).collect::<Vec<_>>().join(",")
                ),
            ),
        };
        out.push(format!("{indent}meta {key} {kind} {text}"));
    }
}

fn tolerance(tolerance: openms::comparison::Tolerance) -> (f64, u8) {
    match tolerance {
        openms::comparison::Tolerance::Ppm(value) => (value, 1),
        openms::comparison::Tolerance::Absolute(value) => (value, 0),
    }
}

/// The oracle driver's rendering of a successful load.
fn render(file: &XQuestResultXmlFile, document: &XQuestDocument) -> Vec<String> {
    let mut out: Vec<String> = document.log.clone();
    out.push("status ok".into());
    out.push(format!(
        "n_hits {} min_score {} max_score {}",
        file.number_of_hits(),
        g17(file.min_score()),
        g17(file.max_score())
    ));
    out.push(format!(
        "protein_ids {}",
        document.protein_identifications.len()
    ));
    for run in &document.protein_identifications {
        let date = run
            .date_time
            .as_deref()
            .map(|date| date.replace('T', " "))
            .unwrap_or_else(|| "0000-00-00 00:00:00".into());
        out.push(format!(
            "prot identifier=[{}] engine=[{}] version=[{}] date=[{date}] score_type=[{}] higher={}",
            run.identifier,
            run.search_engine,
            run.search_engine_version,
            run.score_type,
            u8::from(run.higher_score_better)
        ));
        out.push(format!(" runs {}", run.primary_ms_run_paths.len()));
        for path in &run.primary_ms_run_paths {
            out.push(format!(" run [{path}]"));
        }
        let sp = &run.search_parameters;
        let (precursor, precursor_ppm) = tolerance(sp.precursor_tolerance);
        let (fragment, fragment_ppm) = tolerance(sp.fragment_tolerance);
        let mass_type = match sp.mass_type {
            PeakMassType::Monoisotopic => 0,
            PeakMassType::Average => 1,
        };
        // `EnzymaticDigestion::Specificity`: NONE 0, SEMI 1, FULL 2, UNKNOWN 3.
        let specificity = match sp.enzyme_specificity {
            EnzymeTermSpecificity::None => 0,
            EnzymeTermSpecificity::Semi => 1,
            EnzymeTermSpecificity::Full => 2,
            EnzymeTermSpecificity::Unknown => 3,
        };
        out.push(format!(
            " sp db=[{}] charges=[{}] enzyme=[{}] missed={} prec_tol={} prec_ppm={precursor_ppm} frag_tol={} \
             frag_ppm={fragment_ppm} mass_type={mass_type} specificity={specificity}",
            sp.database,
            sp.charges,
            sp.digestion_enzyme,
            sp.missed_cleavages,
            g17(precursor),
            g17(fragment),
        ));
        let mut fixed = " fixed".to_owned();
        for modification in &sp.fixed_modifications {
            fixed.push_str(&format!(" [{modification}]"));
        }
        out.push(fixed);
        let mut variable = " variable".to_owned();
        for modification in &sp.variable_modifications {
            variable.push_str(&format!(" [{modification}]"));
        }
        out.push(variable);
        render_meta(&mut out, &sp.metadata, "  sp ");
        render_meta(&mut out, &run.metadata, " ");
        out.push(format!(" hits {}", run.hits.len()));
        for hit in &run.hits {
            out.push(format!(
                " phit [{}] score={} rank={}",
                hit.accession,
                g17(hit.score),
                hit.rank
            ));
            render_meta(&mut out, &hit.metadata, "  ");
        }
    }
    out.push(format!(
        "peptide_ids {}",
        document.peptide_identifications.len()
    ));
    for id in &document.peptide_identifications {
        let optional = |value: Option<f64>| value.map_or_else(|| "none".to_owned(), g17);
        out.push(format!(
            "id identifier=[{}] score_type=[{}] higher={} rt={} mz={} hits={}",
            id.identifier,
            id.score_type,
            u8::from(id.higher_score_better),
            optional(id.rt),
            optional(id.mz),
            id.hits.len()
        ));
        render_meta(&mut out, &id.metadata, " ");
        for hit in &id.hits {
            out.push(format!(
                " hit seq=[{}] score={} rank={} charge={} evidences={}",
                source_to_string(&hit.sequence),
                g17(hit.score),
                hit.rank,
                hit.charge,
                hit.evidences.len()
            ));
            for evidence in &hit.evidences {
                let position =
                    |p: Option<usize>| p.map_or_else(|| "-1".to_owned(), |p| p.to_string());
                out.push(format!(
                    "  ev [{}] start={} end={} before={} after={}",
                    evidence.protein_accession,
                    position(evidence.start),
                    position(evidence.end),
                    evidence.aa_before.code(),
                    evidence.aa_after.code()
                ));
            }
            render_meta(&mut out, &hit.metadata, "  ");
        }
    }
    out
}

/// The crate's error kind for each C++ exception the oracle recorded.
fn same_kind(cpp: &str, error: &Error) -> bool {
    let name = cpp.split(':').next().unwrap_or_default();
    match name {
        // `fatalError` and `AASequence::fromString`.
        "status exception Parse Error" => matches!(error, Error::Parse { .. }),
        // `ConversionError` (`toDouble`, `toInt32`, `DataValue` casts),
        // `ElementNotFound` (the enzyme), `IllegalArgument` (`splitByMiddle`)
        // and the bare `std::exception` of an unsupported link type.
        "status exception ConversionError"
        | "status exception ElementNotFound"
        | "status exception IllegalArgument"
        | "status std_exception std::exception" => matches!(error, Error::InvalidValue(_)),
        // `std::string::at` and `IndexOverflow`.
        "status std_exception basic_string::at" | "status exception IndexOverflow" => {
            matches!(error, Error::InvalidRange(_))
        }
        _ => false,
    }
}

fn status_line(text: &str) -> &str {
    text.lines()
        .find(|line| line.starts_with("status "))
        .expect("the driver prints a status")
}

/// Every case the oracle recorded a `load` for.
fn cases(suffix: &str) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(data().join("expected"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter_map(|name| name.strip_suffix(suffix).map(str::to_owned))
        .collect();
    names.sort();
    names
}

// ---------------------------------------------------------------------------
// Executed differential: load
// ---------------------------------------------------------------------------

#[test]
fn release_load_differential() {
    let cases = cases(".load.out");
    assert_eq!(cases.len(), 20, "{cases:?}");
    for case in &cases {
        let recorded = expected(&format!("{case}.load.out"));
        assert_eq!(expected(&format!("{case}.load.exit")).trim(), "0", "{case}");
        let mut file = XQuestResultXmlFile::new();
        let result = file.load(input(&format!("{case}.xquest.xml")));
        let status = status_line(&recorded);
        match result {
            Ok(document) => {
                assert_eq!(status, "status ok", "{case}: the source threw");
                let lines: Vec<&str> = recorded.lines().collect();
                let rendered = render(&file, &document);
                let rendered: Vec<&str> = rendered.iter().map(String::as_str).collect();
                if rendered != lines {
                    let first = rendered
                        .iter()
                        .zip(&lines)
                        .position(|(a, b)| a != b)
                        .unwrap_or(rendered.len().min(lines.len()));
                    panic!(
                        "{case}: first difference at line {first}:\n port:   {:?}\n source: {:?}\n ({} vs {} lines)",
                        rendered.get(first),
                        lines.get(first),
                        rendered.len(),
                        lines.len()
                    );
                }
            }
            Err(error) => {
                assert!(
                    same_kind(status, &error),
                    "{case}: source {status:?}, port {error:?}"
                );
                // A failed load leaves the adapter as it was.
                assert_eq!(file.number_of_hits(), -1, "{case}");
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Executed differential: store
// ---------------------------------------------------------------------------

/// The attributes `writeTo` writes by name, in its order; every other
/// attribute of a `<search_hit>` is a meta value written in the order of the
/// source process's `MetaInfoRegistry` (byte order in the port), so the
/// comparison sorts those.
const FIXED: &[&str] = &[
    "search_hit_rank",
    "id",
    "type",
    "structure",
    "seq1",
    "seq2",
    "prot1",
    "prot2",
    "topology",
    "xlinkposition",
    "Mr",
    "mz",
    "charge",
    "xlinkermass",
    "measured_mass",
    "error",
    "error_rel",
    "xlinkions_matched",
    "backboneions_matched",
    "xcorrx",
    "xcorrb",
    "match_odds",
    "prescore",
    "num_of_matched_ions_alpha",
    "num_of_matched_ions_beta",
    "num_of_matched_common_ions_alpha",
    "num_of_matched_common_ions_beta",
    "num_of_matched_xlink_ions_alpha",
    "num_of_matched_xlink_ions_beta",
    "TIC",
    "wTIC",
    "intsum",
    "xprophet_f",
    "fdr_type",
    "fdr",
    "annotated_spec",
    "score",
];

/// Attributes as written: ` key="value"` pairs, where a key may hold spaces.
fn attributes(body: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut rest = body;
    while let Some(eq) = rest.find("=\"") {
        let key = rest[..eq].trim_start().to_owned();
        let after = &rest[eq + 2..];
        let end = after.find('"').expect("closed attribute value");
        out.push((key, after[..end].to_owned()));
        rest = &after[end + 1..];
    }
    out
}

fn normalize(text: &str) -> Vec<String> {
    text.lines()
        .map(|line| {
            if let Some(body) = line.strip_prefix("<search_hit ") {
                let mut fixed = Vec::new();
                let mut extra = Vec::new();
                for (key, value) in attributes(body.trim_end_matches(" >")) {
                    if FIXED.contains(&key.as_str()) {
                        fixed.push(format!("{key}={value}"));
                    } else {
                        extra.push(format!("{key}={value}"));
                    }
                }
                extra.sort();
                format!("<search_hit {} | {}", fixed.join(" "), extra.join(" "))
            } else if line.starts_with("<xquest_results ") {
                let start = line.find(" date=\"").unwrap() + 7;
                let end = start + line[start..].find('"').unwrap();
                format!("{}<DATE>{}", &line[..start], &line[end..])
            } else {
                line.to_owned()
            }
        })
        .collect()
}

fn compare_written(case: &str, written: &str, recorded: &str) {
    let written = normalize(written);
    let recorded = normalize(recorded);
    if written != recorded {
        let first = written
            .iter()
            .zip(&recorded)
            .position(|(a, b)| a != b)
            .unwrap_or(written.len().min(recorded.len()));
        panic!(
            "{case}: first difference at line {first}:\n port:   {:?}\n source: {:?}\n ({} vs {} lines)",
            written.get(first),
            recorded.get(first),
            written.len(),
            recorded.len()
        );
    }
}

#[test]
fn release_store_differential() {
    let cases = cases(".store.out");
    assert_eq!(cases.len(), 20, "{cases:?}");
    for case in &cases {
        let recorded = expected(&format!("{case}.store.out"));
        let status = status_line(&recorded);
        let loaded = XQuestResultXmlFile::new().load(input(&format!("{case}.xquest.xml")));
        let written = loaded.and_then(|document| {
            let mut text = Vec::new();
            XQuestResultXmlHandler::write_to(
                &mut text,
                &document.protein_identifications,
                &document.peptide_identifications,
                DATE,
            )?;
            Ok(String::from_utf8(text).unwrap())
        });
        match written {
            Ok(text) => {
                assert_eq!(status, "status ok", "{case}: the source threw");
                let stored = expected(&format!("{case}.stored.xquest.xml"));
                compare_written(case, &text, &stored);
            }
            Err(error) => assert!(
                same_kind(status, &error),
                "{case}: source {status:?}, port {error:?}"
            ),
        }
    }
}

#[test]
fn release_convert_differential() {
    let cases = cases(".convert.out");
    assert_eq!(cases.len(), 2, "{cases:?}");
    for case in &cases {
        let recorded = expected(&format!("{case}.convert.out"));
        assert_eq!(status_line(&recorded), "status ok", "{case}");
        let document = idxml::load(input(&format!("{case}.idXML"))).unwrap();
        let mut text = Vec::new();
        XQuestResultXmlHandler::write_to(
            &mut text,
            &document.protein_identifications,
            &document.peptide_identifications,
            DATE,
        )
        .unwrap();
        let stored = expected(&format!("{case}.converted.xquest.xml"));
        compare_written(case, &String::from_utf8(text).unwrap(), &stored);
    }
}

// ---------------------------------------------------------------------------
// Executed differential: splitByNth / splitByMiddle
// ---------------------------------------------------------------------------

#[test]
fn release_split_differential() {
    let table = fs::read_to_string(data().join("splits.tsv")).unwrap();
    let mut count = 0;
    for (index, row) in table.lines().enumerate() {
        let fields: Vec<&str> = row.split('\t').collect();
        let (how, text, separator, n) = (fields[0], fields[1], fields[2], fields[3]);
        let separator = separator.chars().next().unwrap();
        let result = if how == "middle" {
            XQuestResultXmlHandler::split_by_middle(text, separator)
        } else {
            XQuestResultXmlHandler::split_by_nth(text, separator, n.parse().unwrap())
        };
        let recorded = expected(&format!("split{}.out", index + 1));
        let status = status_line(&recorded);
        match result {
            Ok(parts) => {
                assert_eq!(status, "status ok", "{row}");
                let rendered: Vec<String> = parts.iter().map(|p| format!("part [{p}]")).collect();
                let lines: Vec<&str> = recorded.lines().skip(1).collect();
                assert_eq!(rendered, lines, "{row}");
            }
            Err(error) => {
                assert!(same_kind(status, &error), "{row}: {status:?} vs {error:?}");
                if let Error::InvalidRange(message) = &error {
                    // The message of `std::string::at` itself.
                    assert_eq!(format!("status std_exception {message}"), status, "{row}");
                }
            }
        }
        count += 1;
    }
    assert_eq!(count, 14);
}

// ---------------------------------------------------------------------------
// START_SECTION(void store(...)) of XQuestResultXMLFile_test.cpp, which
// loads the fixture, stores it, loads it back and compares.
// ---------------------------------------------------------------------------

#[test]
fn class_test_store_and_reload() {
    let mut file = XQuestResultXmlFile::new();
    let original = file
        .load(input("XQuestResultXMLFile_test_data.xquest.xml"))
        .unwrap();
    let directory = std::env::temp_dir().join(format!("xquest_rt_{}", std::process::id()));
    fs::create_dir_all(&directory).unwrap();
    let path = directory.join("roundtrip.xquest.xml");
    file.store(
        &path,
        &original.protein_identifications,
        &original.peptide_identifications,
    )
    .unwrap();
    let reloaded = XQuestResultXmlFile::new().load(&path).unwrap();
    fs::remove_dir_all(&directory).unwrap();

    let before = &original.peptide_identifications;
    let after = &reloaded.peptide_identifications;
    let text = |meta: &MetaInfo, key: &str| meta.get(key).map(|value| value.to_string());
    let close = |a: f64, b: f64| (a - b).abs() <= 1e-5 * a.abs().max(b.abs()).max(1e-5);
    for i in (0..after.len()).step_by(20) {
        for (k, hit) in after[i].hits.iter().enumerate() {
            let source = &before[i].hits[k];
            assert!(close(hit.score, source.score));
            assert_eq!(hit.charge, source.charge);
            for key in ["xl_rank", "xl_type", "xl_pos1", "xl_pos2", "sequence_beta"] {
                assert_eq!(
                    text(&hit.metadata, key),
                    text(&source.metadata, key),
                    "{i} {key}"
                );
            }
            assert_eq!(hit.sequence.to_string(), source.sequence.to_string());
            for key in ["OpenPepXL:match-odds", "OpenPepXL:intsum"] {
                assert!(close(
                    hit.metadata[key].as_f64().unwrap(),
                    source.metadata[key].as_f64().unwrap()
                ));
            }
        }
    }

    assert_eq!(after.len(), 296);
    let check = |index: usize,
                 hits: usize,
                 pos1: i64,
                 pos2: i64,
                 alpha: &str,
                 beta_spec: &str,
                 seq: &str,
                 beta: &str,
                 spec: &str| {
        let id = &after[index];
        assert_eq!(id.hits.len(), hits, "{index}");
        let hit = &id.hits[0];
        assert_eq!(hit.metadata["xl_type"].as_str().unwrap(), "cross-link");
        assert_eq!(hit.metadata["xl_pos1"].as_i64().unwrap(), pos1);
        assert_eq!(hit.metadata["xl_pos2"].as_i64().unwrap(), pos2);
        assert_eq!(hit.metadata["xl_term_spec_alpha"].as_str().unwrap(), alpha);
        assert_eq!(
            hit.metadata["xl_term_spec_beta"].as_str().unwrap(),
            beta_spec
        );
        assert_eq!(hit.sequence.to_string(), seq);
        assert_eq!(hit.metadata["sequence_beta"].as_str().unwrap(), beta);
        let _ = spec;
    };
    check(
        0,
        1,
        14,
        5,
        "ANYWHERE",
        "ANYWHERE",
        "LTEIISHDPNIELHKK",
        "VEGCPKHPK",
        "",
    );
    check(
        17,
        1,
        15,
        11,
        "C_TERM",
        "ANYWHERE",
        "VILHLKEDQTEYLEER",
        "EYGCAPWPMVEKLIK",
        "",
    );
    check(
        289,
        1,
        15,
        0,
        "ANYWHERE",
        "ANYWHERE",
        "DYHFVNATEESDALAKLR",
        "KETFDDLPK",
        "",
    );
    check(
        279,
        2,
        0,
        6,
        "N_TERM",
        "ANYWHERE",
        "MASGSCQGCEEDEETLKK",
        "NTEGTQKQK",
        "",
    );

    // XFDR attributes survive the load and the store/load round trip, on
    // exactly one hit, found by its id.
    for ids in [before, after] {
        let mut matches = 0;
        for id in ids {
            for hit in &id.hits {
                if text(&hit.metadata, "OpenPepXL:id").as_deref()
                    != Some("LTEIISHDPNIELHKK-VEGCPKHPK-a15-b6")
                {
                    continue;
                }
                matches += 1;
                assert_eq!(hit.metadata["XFDR:used_for_FDR"].as_i64().unwrap(), 1);
                assert_eq!(hit.metadata["XFDR:fdr_type"].as_str().unwrap(), "q-value");
            }
        }
        assert_eq!(matches, 1);
    }
}

// ---------------------------------------------------------------------------
// Native behaviour
// ---------------------------------------------------------------------------

#[test]
fn getters_before_a_load_and_after_a_failed_one() {
    let mut file = XQuestResultXmlFile::new();
    assert_eq!(file.number_of_hits(), -1);
    assert_eq!(file.min_score(), 0.0);
    assert_eq!(file.max_score(), 0.0);
    file.load(input("x01_opxl_types.xquest.xml")).unwrap();
    assert_eq!(file.number_of_hits(), 4);
    assert!(file.load(input("x06_unsupported_type.xquest.xml")).is_err());
    assert_eq!(file.number_of_hits(), 4);
    assert_eq!(file.min_score(), -2.5);
    assert_eq!(file.max_score(), 10.5);
}

#[test]
fn store_refuses_another_known_extension() {
    let error = XQuestResultXmlFile::new()
        .store("out.idXML", &[], &[])
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
}

#[test]
fn store_without_a_run_is_refused() {
    let mut text = Vec::new();
    let error = XQuestResultXmlHandler::write_to(&mut text, &[], &[], DATE).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert!(text.is_empty());
}

#[test]
fn ceilings_refuse_instead_of_truncating() {
    let path = input("x12_many_ties.xquest.xml");
    let read = |limits: ReadLimits| {
        XQuestResultXmlFile::new().read_with_limits(
            std::io::BufReader::new(fs::File::open(&path).unwrap()),
            &limits,
        )
    };
    assert!(read(ReadLimits::default()).is_ok());
    let hits = read(ReadLimits {
        max_hits: 20,
        ..ReadLimits::default()
    });
    assert!(matches!(hits, Err(Error::InvalidValue(_))), "{hits:?}");
    let accessions = read(ReadLimits {
        max_accessions_per_hit: 1,
        ..ReadLimits::default()
    });
    assert!(
        matches!(accessions, Err(Error::InvalidValue(_))),
        "{accessions:?}"
    );
    let sequence = read(ReadLimits {
        max_sequence_bytes: 4,
        ..ReadLimits::default()
    });
    assert!(
        matches!(sequence, Err(Error::InvalidValue(_))),
        "{sequence:?}"
    );
    let bytes = read(ReadLimits {
        max_bytes: 100,
        ..ReadLimits::default()
    });
    assert!(bytes.is_err());
}

#[test]
fn attribute_text_the_source_indexes_past_is_refused() {
    let document = |spectrum: &str| {
        format!(
            "<?xml version=\"1.0\"?>\n<xquest_results xquest_version=\"OpenPepXL 1.0\" date=\"2018-03-06 11:52:13\" \
             tolerancemeasure_ms1=\"ppm\" tolerancemeasure_ms2=\"Da\" ms1tolerance=\"10\" ms2tolerance=\"0.2\" \
             xlink_ms2tolerance=\"0.3\" crosslinkername=\"DSS\" xlinkermw=\"138.0680796\" database=\"db\" \
             AArequired1=\"K\" AArequired2=\"K\" cp_isotopediff=\"\" enzyme_name=\"Trypsin\" missed_cleavages=\"2\" \
             ntermxlinkable=\"0\" charges=\"2\" run_path=\"\">\n{spectrum}\n</xquest_results>\n"
        )
    };
    // One retention time where the source reads two.
    let text = document(
        "<spectrum_search spectrum=\"a.1\" rtsecscans=\"100\" mzscans=\"\" mz_precursor=\"500\" \
         charge_precursor=\"2\" scan_index_light=\"1\" scan_index_heavy=\"1\"></spectrum_search>",
    );
    let result = XQuestResultXmlFile::new().read(text.as_bytes());
    assert!(matches!(result, Err(Error::InvalidRange(_))), "{result:?}");
    // A well-formed minimal file without hits loads.
    let result = XQuestResultXmlFile::new().read(document("").as_bytes());
    let document = result.unwrap();
    assert!(document.peptide_identifications.is_empty());
    assert_eq!(
        document.protein_identifications[0]
            .search_parameters
            .charges,
        "2"
    );
}

#[test]
fn to_str_matches_the_recorded_rt_text() {
    // `rtsecscans="200.0:200.0"` in the recorded store of x01.
    assert_eq!(to_str(200.0), "200.0");
}

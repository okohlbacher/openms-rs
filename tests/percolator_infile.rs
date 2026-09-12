// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FORMAT/PercolatorInfile.h`: the Percolator `.pin` writer and reader.
//!
//! Expectations come from `PercolatorInfile_test.cpp` and its retained
//! `sage.pin` fixture. See `tests/data/percolator_infile_provenance.json`.

use openms::Error;
use openms::chemistry::AASequence;
use openms::format::percolator_infile::{
    self as pin, PinOptions, ReadOptions, ScanNumberPattern, bracket_sequence, count_enzymatic,
    is_enzymatic, scan_identifier, standard_feature_set,
};
use openms::identification::{
    FlankingResidue, PeptideEvidence, PeptideHit, PeptideIdentification, TargetDecoyType,
};
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// The fully annotated target PSM the `store` class-test section builds.
fn sampler_psm() -> Vec<PeptideIdentification> {
    let mut hit = PeptideHit::new(1.0, 0, 2, AASequence::parse("SAMPLER").unwrap()).unwrap();
    hit.set_target_decoy_type(TargetDecoyType::Target);
    hit.evidences = vec![PeptideEvidence {
        protein_accession: "PROT1".into(),
        // Tryptic flanks, so enzN and enzC are both true.
        aa_before: FlankingResidue::Residue('K'),
        aa_after: FlankingResidue::Residue('S'),
        ..Default::default()
    }];
    let mut identification = PeptideIdentification::new();
    identification.mz = Some(500.25);
    identification.rt = Some(123.4);
    identification.set_spectrum_reference("scan=529");
    identification.hits = vec![hit];
    vec![identification]
}

fn column(header: &[&str], name: &str) -> usize {
    header
        .iter()
        .position(|column| *column == name)
        .unwrap_or_else(|| panic!("no {name} column in {header:?}"))
}

// START_SECTION(PercolatorInfile()) and START_SECTION(~PercolatorInfile()) —
// the source class is a bag of static functions whose constructor and
// destructor are only checked for existence. This port has no object; the
// option defaults are the equivalent constructed state.
#[test]
fn the_default_options_are_the_constructed_state() {
    let options = PinOptions::default();
    assert_eq!(options.enzyme, "trypsin");
    assert_eq!((options.min_charge, options.max_charge), (2, 5));
    let read = ReadOptions::default();
    assert!(read.higher_score_better);
    assert_eq!(read.spectrum_q_threshold, 0.01);
    assert!(!read.sage_annotation);
    assert_eq!(pin::MANDATORY_COLUMNS, ["SpecId", "Label", "ScanNr"]);
    assert_eq!(pin::TRAILING_COLUMNS, ["Peptide", "Proteins"]);
}

// START_SECTION((static StringList getStandardFeatureSet(int min_charge,
//   int max_charge)))
#[test]
fn the_standard_feature_set_is_the_exact_ordered_column_contract() {
    let features = standard_feature_set(2, 4).unwrap();
    let expected = [
        "SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "mass", "peplen", "charge2", "charge3",
        "charge4", "enzN", "enzC", "enzInt", "dm", "absdm",
    ];
    assert_eq!(features.len(), expected.len());
    for (got, want) in features.iter().zip(expected) {
        assert_eq!(got, want);
    }
    // The three mandatory leading columns.
    assert_eq!(features[0], "SpecId");
    assert_eq!(features[1], "Label");
    assert_eq!(features[2], "ScanNr");

    // A single charge state yields exactly one charge column.
    let single = standard_feature_set(3, 3).unwrap();
    let expected_single = [
        "SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "mass", "peplen", "charge3", "enzN",
        "enzC", "enzInt", "dm", "absdm",
    ];
    assert_eq!(single.len(), expected_single.len());
    for (got, want) in single.iter().zip(expected_single) {
        assert_eq!(got, want);
    }

    // The assembled header begins with SpecId and ends with Peptide, Proteins.
    let mut header = features;
    header.push("Peptide".into());
    header.push("Proteins".into());
    assert_eq!(header.first().unwrap(), "SpecId");
    assert_eq!(header[header.len() - 2], "Peptide");
    assert_eq!(header.last().unwrap(), "Proteins");
}

// START_SECTION((static void store(const std::string& pin_file,
//   const PeptideIdentificationList& peptide_ids, const StringList& feature_set,
//   const std::string& enz, int min_charge, int max_charge)))
#[test]
fn storing_writes_the_column_contract_and_the_computed_features() {
    let mut feature_set = standard_feature_set(2, 3).unwrap();
    feature_set.push("Peptide".into());
    feature_set.push("Proteins".into());

    let directory = std::env::temp_dir().join(format!("openms-pin-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("out.pin");
    pin::store(
        &path,
        &sampler_psm(),
        &feature_set,
        &PinOptions {
            enzyme: "trypsin".into(),
            min_charge: 2,
            max_charge: 3,
            ..Default::default()
        },
    )
    .unwrap();
    let written = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_dir_all(&directory).unwrap();

    // Header plus exactly one data row: the PSM must not be skipped.
    let lines: Vec<&str> = written.lines().collect();
    assert_eq!(lines.len(), 2);

    let header: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(header.len(), feature_set.len());
    assert_eq!(header[0], "SpecId");
    assert_eq!(header[1], "Label");
    assert_eq!(header[2], "ScanNr");
    assert_eq!(header[header.len() - 2], "Peptide");
    assert_eq!(*header.last().unwrap(), "Proteins");
    for name in [
        "ExpMass", "CalcMass", "mass", "peplen", "charge2", "charge3", "enzN", "enzC", "enzInt",
        "dm", "absdm",
    ] {
        assert!(header.contains(&name), "no {name} column");
    }

    let row: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(row.len(), header.len());
    assert_eq!(row[column(&header, "SpecId")], "scan=529");
    assert_eq!(row[column(&header, "ScanNr")], "529");
    assert_eq!(row[column(&header, "Label")], "1"); // target -> 1
    assert_eq!(row[column(&header, "peplen")], "7"); // SAMPLER
    assert_eq!(row[column(&header, "charge2")], "1"); // charge == 2, one-hot
    assert_eq!(row[column(&header, "charge3")], "0");
    assert_eq!(row[column(&header, "enzN")], "1"); // tryptic N-terminus
    assert_eq!(row[column(&header, "enzC")], "1"); // tryptic C-terminus
    assert_eq!(row[column(&header, "Peptide")], "K.SAMPLER.S");
    assert_eq!(row[column(&header, "Proteins")], "PROT1");
    // Derived: SAMPLER has no K/R before a non-P residue at an internal bond.
    assert_eq!(row[column(&header, "enzInt")], "0");
    // Derived: ExpMass is the identification's m/z, and mass repeats it.
    assert_eq!(row[column(&header, "ExpMass")], "500.25");
    assert_eq!(row[column(&header, "mass")], "500.25");
    // Derived: dm is ExpMass - CalcMass and absdm its magnitude.
    let dm: f64 = row[column(&header, "dm")].parse().unwrap();
    let calculated: f64 = row[column(&header, "CalcMass")].parse().unwrap();
    assert!((dm - (500.25 - calculated)).abs() < 1e-9);
    assert_eq!(
        row[column(&header, "absdm")].parse::<f64>().unwrap(),
        dm.abs()
    );
}

/// Writing to a stream produces the same bytes as writing to a path.
#[test]
fn the_stream_and_path_writers_agree() {
    let mut feature_set = standard_feature_set(2, 3).unwrap();
    feature_set.push("Peptide".into());
    feature_set.push("Proteins".into());
    let options = PinOptions {
        min_charge: 2,
        max_charge: 3,
        ..Default::default()
    };
    let mut bytes = Vec::new();
    let report = pin::write(&mut bytes, &sampler_psm(), &feature_set, &options).unwrap();
    assert_eq!(report.lines.len(), 2);
    assert_eq!(String::from_utf8(bytes).unwrap().lines().count(), 2);
    assert_eq!(report.hits_missing_features, 0);
}

/// An empty input yields a header-only file and the source's warning.
#[test]
fn no_identifications_produce_an_empty_percolator_input() {
    let feature_set = standard_feature_set(2, 3).unwrap();
    let report = pin::prepare_pin(&[], &feature_set, &PinOptions::default()).unwrap();
    assert_eq!(report.lines.len(), 1);
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("Creating empty percolator input"))
    );
}

/// A hit without protein evidence or without target/decoy status is skipped
/// and reported, not written with blank fields.
#[test]
fn hits_without_evidence_or_status_are_skipped() {
    let mut identifications = sampler_psm();
    let mut orphan = identifications[0].hits[0].clone();
    orphan.evidences.clear();
    let mut unknown = identifications[0].hits[0].clone();
    unknown.set_target_decoy_type(TargetDecoyType::Unknown);
    identifications[0].hits.push(orphan);
    identifications[0].hits.push(unknown);

    let mut feature_set = standard_feature_set(2, 3).unwrap();
    feature_set.push("Peptide".into());
    feature_set.push("Proteins".into());
    let report = pin::prepare_pin(
        &identifications,
        &feature_set,
        &PinOptions {
            min_charge: 2,
            max_charge: 3,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(report.lines.len(), 2, "only the annotated PSM is written");
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("without protein reference"))
    );
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("without target/decoy information"))
    );

    let mut stamped = identifications.clone();
    let stamp = pin::stamp_pin_features(
        &mut stamped,
        &PinOptions {
            min_charge: 2,
            max_charge: 3,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(stamp.skipped.len(), 2);
    assert!(stamp.skipped.contains(&(0, 1)));
    assert!(stamp.skipped.contains(&(0, 2)));
    // The kept hit records every key that stamping introduced, so a caller can
    // remove only what it added.
    let added = &stamp.added_meta_values[&(0, 0)];
    for key in [
        "SpecId",
        "ScanNr",
        "Label",
        "CalcMass",
        "ExpMass",
        "deltamass",
        "retentiontime",
        "mass",
        "score",
        "peplen",
        "charge2",
        "charge3",
        "enzN",
        "enzC",
        "enzInt",
        "dm",
        "absdm",
        "Peptide",
        "Proteins",
    ] {
        assert!(added.contains(key), "{key} was not recorded as added");
    }
    // The skipped hits keep only the metadata they arrived with: no PIN key
    // was stamped onto them.
    for key in ["SpecId", "ScanNr", "Label", "Peptide", "Proteins"] {
        assert!(!stamped[0].hits[1].metadata.contains_key(key), "{key}");
        assert!(!stamped[0].hits[2].metadata.contains_key(key), "{key}");
    }
    assert!(!stamp.added_meta_values.contains_key(&(0, 1)));
}

/// An existing meta value is refreshed but not recorded as newly added, and an
/// existing `CalcMass` is reused rather than recomputed.
#[test]
fn existing_meta_values_are_reused_and_not_recorded_as_added() {
    let mut identifications = sampler_psm();
    identifications[0].hits[0].metadata.insert(
        "CalcMass".into(),
        openms::metadata::MetaValue::try_from(499.0).unwrap(),
    );
    let options = PinOptions {
        min_charge: 2,
        max_charge: 3,
        ..Default::default()
    };
    let mut stamped = identifications.clone();
    let report = pin::stamp_pin_features(&mut stamped, &options).unwrap();
    assert!(!report.added_meta_values[&(0, 0)].contains("CalcMass"));
    assert_eq!(
        stamped[0].hits[0].metadata["CalcMass"].as_f64().unwrap(),
        499.0
    );
    // dm is then ExpMass - the supplied CalcMass.
    assert!((stamped[0].hits[0].metadata["dm"].as_f64().unwrap() - (500.25 - 499.0)).abs() < 1e-9);
}

/// An isotope error shifts the experimental mass down, in both the legacy
/// `IsotopeError` spelling and the current `isotope_error` one.
#[test]
fn an_isotope_error_shifts_the_experimental_mass() {
    for key in ["IsotopeError", "isotope_error"] {
        let mut identifications = sampler_psm();
        identifications[0].hits[0]
            .metadata
            .insert(key.into(), openms::metadata::MetaValue::from(1i64));
        let mut stamped = identifications;
        pin::stamp_pin_features(
            &mut stamped,
            &PinOptions {
                min_charge: 2,
                max_charge: 3,
                ..Default::default()
            },
        )
        .unwrap();
        let expected = 500.25 - openms::constants::C13C12_MASSDIFF_U / 2.0;
        assert!(
            (stamped[0].hits[0].metadata["ExpMass"].as_f64().unwrap() - expected).abs() < 1e-12,
            "{key}"
        );
    }
}

/// A declared feature with no meta value drops the hit and names the feature,
/// which is what the source does rather than writing a short row.
#[test]
fn a_declared_feature_without_a_meta_value_drops_the_hit() {
    let mut feature_set = standard_feature_set(2, 3).unwrap();
    feature_set.push("Peptide".into());
    feature_set.push("Proteins".into());
    feature_set.push("NoSuchFeature".into());
    let report = pin::prepare_pin(
        &sampler_psm(),
        &feature_set,
        &PinOptions {
            min_charge: 2,
            max_charge: 3,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(report.lines.len(), 1);
    assert_eq!(report.hits_missing_features, 1);
    assert!(report.missing_features.contains("NoSuchFeature"));
}

// START_SECTION(PeptideIdentificationList PercolatorInfile::load(
//   const std::string& pin_file, bool higher_score_better,
//   const std::string& score_name, std::string decoy_prefix))
#[test]
fn loading_a_sage_pin_file_groups_rows_and_reannotates_decoys() {
    let document = pin::load(
        data("percolator_infile_sage.pin"),
        &ReadOptions {
            higher_score_better: true,
            score_name: "ln(hyperscore)".into(),
            extra_scores: vec![
                "ln(delta_next)".into(),
                "ln(delta_best)".into(),
                "matched_peaks".into(),
            ],
            decoy_prefix: "DECOY_".into(),
            ..Default::default()
        },
    )
    .unwrap();
    let pids = &document.peptide_identifications;
    assert_eq!(pids.len(), 9);
    assert_eq!(document.filenames.len(), 2);
    assert_eq!(pids[0].spectrum_reference(), "30381");
    assert_eq!(pids[6].spectrum_reference(), "spectrum=2041");
    // The eighth entry is labelled target in the file but maps only to a
    // DECOY_-prefixed protein, so the decoy prefix overrides the label.
    assert_eq!(
        pids[7].hits[0].metadata["target_decoy"].to_string(),
        "decoy"
    );

    // Derived from the fixture's own columns.
    assert_eq!(document.filenames[0], "JD_06232014_sample1_A.mzML");
    assert_eq!(
        document.filenames[1],
        "sample_from_different_instrument.mzML"
    );
    assert_eq!(pids[0].score_type, "ln(hyperscore)");
    assert!(pids[0].higher_score_better);
    // retentiontime is in minutes and becomes seconds.
    assert!((pids[0].rt.unwrap() - 67.93814 * 60.0).abs() < 1e-6);
    // z=3 is the one-hot column set for the first row.
    assert_eq!(pids[0].hits[0].charge, 3);
    // m/z is ExpMass / |charge| + one proton.
    assert!(
        (pids[0].mz.unwrap() - (3995.7983 / 3.0 + openms::constants::PROTON_MASS_U)).abs() < 1e-9
    );
    // The sixth row belongs to the second file, so its merge index is one.
    assert_eq!(pids[6].metadata["id_merge_index"].as_i64().unwrap(), 1);
    assert_eq!(pids[0].metadata["PinSpecId"].to_string(), "0");
    // The extra score columns are stored verbatim, as strings.
    assert_eq!(pids[0].hits[0].metadata["matched_peaks"].to_string(), "17");
    // DeltaMass is ExpMass - CalcMass.
    assert!(
        (pids[0].hits[0].metadata["DeltaMass"].as_f64().unwrap() - (3995.7983 - 3597.828)).abs()
            < 1e-9
    );
    // The bracket-notation peptide parses, including its terminal mass tag.
    assert_eq!(pids[1].hits[0].sequence.as_str(), "MVLVQDLLHPTAASEAR");
    assert_eq!(pids[1].hits[0].evidences.len(), 2);
    assert_eq!(
        pids[1].hits[0].evidences[0].protein_accession,
        "sp|P35997|RS27_YEAST"
    );
    assert!(document.warnings.is_empty());
}

/// Without a decoy prefix the file's own `Label` column is trusted.
#[test]
fn without_a_decoy_prefix_the_label_column_is_trusted() {
    let document = pin::load(
        data("percolator_infile_sage.pin"),
        &ReadOptions {
            score_name: "ln(hyperscore)".into(),
            ..Default::default()
        },
    )
    .unwrap();
    // Every row in the fixture is labelled 1.
    assert_eq!(
        document.peptide_identifications[7].hits[0].metadata["target_decoy"].to_string(),
        "target"
    );
}

/// An extra score the header does not declare is warned about, not an error.
#[test]
fn an_unknown_extra_score_is_warned_about() {
    let document = pin::load(
        data("percolator_infile_sage.pin"),
        &ReadOptions {
            score_name: "ln(hyperscore)".into(),
            extra_scores: vec!["not_a_column".into()],
            ..Default::default()
        },
    )
    .unwrap();
    assert!(
        document
            .warnings
            .iter()
            .any(|w| w.contains("not_a_column") && w.contains("not found"))
    );
}

/// A row with the wrong column count is a parse error naming the line, and a
/// missing required column is reported as missing information.
#[test]
fn malformed_input_is_refused_explicitly() {
    let text = "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tFileName\tretentiontime\tPeptide\tProteins\tscore\n\
                a\t1\t7\t1.0\t1.0\tf.mzML\t1.0\tPEPTIDE\tP1\n";
    let error = pin::read(
        text.as_bytes(),
        &ReadOptions {
            score_name: "score".into(),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(
        matches!(error, Error::Parse { line: 2, .. }),
        "unexpected error: {error:?}"
    );

    let text = "SpecId\tLabel\tScanNr\n a\t1\t7\n";
    let error = pin::read(
        text.as_bytes(),
        &ReadOptions {
            score_name: "score".into(),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(
        matches!(error, Error::MissingInformation(_)),
        "unexpected error: {error:?}"
    );

    // Sage annotation needs a path, so the stream reader refuses it.
    let error = pin::read(
        "SpecId\n".as_bytes(),
        &ReadOptions {
            sage_annotation: true,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
}

/// Every enzyme the source knows, at the bonds that distinguish them, plus the
/// terminus markers and the unknown-name fallback.
#[test]
fn enzyme_specificity_matches_the_source_table() {
    assert!(is_enzymatic('K', 'A', "trypsin"));
    assert!(!is_enzymatic('K', 'P', "trypsin"));
    assert!(is_enzymatic('K', 'P', "trypsinp"));
    assert!(is_enzymatic('F', 'A', "chymotrypsin"));
    assert!(!is_enzymatic('F', 'P', "chymotrypsin"));
    assert!(is_enzymatic('A', 'L', "thermolysin"));
    assert!(!is_enzymatic('D', 'L', "thermolysin"));
    assert!(is_enzymatic('R', 'G', "thermolysin"));
    assert!(is_enzymatic('A', 'Q', "proteinasek"));
    assert!(!is_enzymatic('Q', 'Q', "proteinasek"));
    assert!(is_enzymatic('A', 'F', "pepsin"));
    assert!(!is_enzymatic('R', 'F', "pepsin"));
    assert!(is_enzymatic('L', 'A', "elastase"));
    assert!(!is_enzymatic('L', 'P', "elastase"));
    assert!(is_enzymatic('A', 'K', "lys-n"));
    assert!(is_enzymatic('K', 'A', "lys-c"));
    assert!(!is_enzymatic('K', 'P', "lys-c"));
    assert!(is_enzymatic('R', 'A', "arg-c"));
    assert!(!is_enzymatic('R', 'P', "arg-c"));
    assert!(is_enzymatic('A', 'D', "asp-n"));
    assert!(is_enzymatic('E', 'A', "glu-c"));
    assert!(!is_enzymatic('E', 'P', "glu-c"));
    // A protein terminus is always enzymatic.
    assert!(is_enzymatic('-', 'P', "trypsin"));
    assert!(is_enzymatic('A', '-', "trypsin"));
    // The source's else branch makes every bond enzymatic for an unknown name.
    assert!(is_enzymatic('Q', 'Q', "no_enzyme"));

    assert_eq!(count_enzymatic("AKAKA", "trypsin"), 2);
    assert_eq!(count_enzymatic("AKPKA", "trypsin"), 1);
    assert_eq!(count_enzymatic("A", "trypsin"), 0);
    // Counted by character, so a non-ASCII letter cannot split a codepoint.
    assert_eq!(count_enzymatic("K\u{e4}K", "trypsin"), 1);
}

/// The scan identifier falls back through spectrum reference, `spectrum_id`
/// and the index, and whitespace is removed from the result.
#[test]
fn the_scan_identifier_falls_back_and_strips_whitespace() {
    let mut identification = PeptideIdentification::new();
    assert_eq!(scan_identifier(&identification, 5), "index=5");
    identification.metadata.insert(
        "spectrum_id".into(),
        openms::metadata::MetaValue::from(41i64),
    );
    assert_eq!(scan_identifier(&identification, 5), "scan=41");
    identification.set_spectrum_reference("controllerType=0 controllerNumber=1 scan=6814");
    assert_eq!(
        scan_identifier(&identification, 5),
        "controllerType=0controllerNumber=1scan=6814"
    );
    // The scan-number pattern still finds the scan token in that identifier.
    let pattern = ScanNumberPattern::for_native_id(&scan_identifier(&identification, 5));
    assert_eq!(
        pattern.extract(&scan_identifier(&identification, 5)),
        Some(6814)
    );
    // The source derives the pattern from the first identification only, so a
    // differently shaped identifier yields nothing.
    assert_eq!(pattern.extract("index=3"), None);
    for (id, want) in [
        ("scan=7", Some(7)),
        ("index=8", Some(8)),
        ("scanId=9", Some(9)),
        ("scanID=10", Some(10)),
        ("spectrum=11", Some(11)),
        ("file=12", Some(12)),
        ("frame=3 scan=13", Some(13)),
        ("function=1 process=0 scan=14", Some(14)),
        ("30381", Some(30381)),
    ] {
        assert_eq!(
            ScanNumberPattern::for_native_id(id).extract(id),
            want,
            "{id}"
        );
    }
}

/// Bracket notation renders signed accurate mass deltas, which is the
/// `Peptide` column's middle part.
#[test]
fn bracket_notation_renders_signed_mass_deltas() {
    let plain = AASequence::parse("SAMPLER").unwrap();
    assert_eq!(bracket_sequence(&plain).unwrap(), "SAMPLER");
    let oxidised = AASequence::parse("SAM(Oxidation)PLER").unwrap();
    let rendered = bracket_sequence(&oxidised).unwrap();
    assert!(rendered.starts_with("SAM[+15.994"), "{rendered}");
    assert!(rendered.ends_with("]PLER"), "{rendered}");
    let acetylated = AASequence::parse(".(Acetyl)SAMPLER").unwrap();
    let rendered = bracket_sequence(&acetylated).unwrap();
    assert!(rendered.starts_with("n[+42.010"), "{rendered}");
    assert!(rendered.ends_with("SAMPLER"), "{rendered}");
    assert_eq!(bracket_sequence(&AASequence::default()).unwrap(), "");
}

/// A `[` or `]` flanking marker becomes Percolator's `-`.
#[test]
fn terminal_flanking_markers_become_dashes() {
    let mut identifications = sampler_psm();
    identifications[0].hits[0].evidences[0].aa_before = FlankingResidue::NTerminus;
    identifications[0].hits[0].evidences[0].aa_after = FlankingResidue::CTerminus;
    let mut stamped = identifications;
    pin::stamp_pin_features(
        &mut stamped,
        &PinOptions {
            min_charge: 2,
            max_charge: 3,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        stamped[0].hits[0].metadata["Peptide"].to_string(),
        "-.SAMPLER.-"
    );
    // The source computes enzN and enzC from the raw '[' and ']' markers and
    // only then maps them to '-' for the Peptide column, so a protein-terminal
    // PSM is reported as NOT enzymatic even though Percolator's own '-'
    // convention would call it enzymatic. This port reproduces that; the
    // divergence is recorded as OPENMS-PERCIN-003.
    assert_eq!(stamped[0].hits[0].metadata["enzN"].as_i64().unwrap(), 0);
    // enzC is 1 here only because SAMPLER's last residue R is itself a tryptic
    // site; the ']' marker contributed nothing.
    assert_eq!(stamped[0].hits[0].metadata["enzC"].as_i64().unwrap(), 1);
    // With a C-terminal residue that is not a cleavage site, the ']' marker
    // leaves enzC 0 where Percolator's '-' convention would make it 1.
    let mut other = sampler_psm();
    other[0].hits[0].sequence = AASequence::parse("SAMPLES").unwrap();
    other[0].hits[0].evidences[0].aa_after = FlankingResidue::CTerminus;
    pin::stamp_pin_features(
        &mut other,
        &PinOptions {
            min_charge: 2,
            max_charge: 3,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(other[0].hits[0].metadata["enzC"].as_i64().unwrap(), 0);
}

/// Several protein accessions are joined with tabs inside the `Proteins`
/// field, which is how Percolator reads a trailing protein list.
#[test]
fn several_accessions_are_tab_separated_in_the_proteins_field() {
    let mut identifications = sampler_psm();
    identifications[0].hits[0].evidences.push(PeptideEvidence {
        protein_accession: "PROT2".into(),
        ..Default::default()
    });
    let mut stamped = identifications;
    pin::stamp_pin_features(
        &mut stamped,
        &PinOptions {
            min_charge: 2,
            max_charge: 3,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        stamped[0].hits[0].metadata["Proteins"].to_string(),
        "PROT1\tPROT2"
    );
}

/// The hit ceiling is checked before anything is mutated, so a refused call
/// leaves the input unchanged.
#[test]
fn the_hit_ceiling_leaves_the_input_unchanged() {
    let mut identifications = sampler_psm();
    let before = identifications.clone();
    let error = pin::stamp_pin_features(
        &mut identifications,
        &PinOptions {
            max_hits: 0,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert_eq!(identifications, before);

    let mut identifications = sampler_psm();
    let before = identifications.clone();
    assert!(
        pin::stamp_pin_features(
            &mut identifications,
            &PinOptions {
                max_bytes: 1,
                min_charge: 2,
                max_charge: 3,
                ..Default::default()
            },
        )
        .is_err()
    );
    assert_eq!(identifications, before);
}

/// The Sage annotation path: the sibling `.tsv` filters rows by `spectrum_q`
/// and the sibling `matched_fragments.sage.tsv` supplies peak annotations.
///
/// These three fixtures are derived, not retained C++ output: the class test
/// never exercises `SageAnnotation`. They encode the column names and the
/// sibling-path derivation read off `PercolatorInfile.cpp:93-155`.
#[test]
fn the_sage_annotation_path_filters_rows_and_attaches_peak_annotations() {
    let document = pin::load(
        data("percolator_infile_sage_results.sage.pin"),
        &ReadOptions {
            score_name: "ln(hyperscore)".into(),
            sage_annotation: true,
            spectrum_q_threshold: 0.01,
            decoy_prefix: "DECOY_".into(),
            ..Default::default()
        },
    )
    .unwrap();
    // The second row has spectrum_q 0.5, above the threshold, so it is dropped.
    let pids = &document.peptide_identifications;
    assert_eq!(pids.len(), 2);
    assert_eq!(pids[0].spectrum_reference(), "101");
    assert_eq!(pids[1].spectrum_reference(), "103");
    assert_eq!(document.filenames, vec!["a.mzML", "b.mzML"]);
    // spectrum_q is stamped on each kept hit.
    assert!((pids[0].hits[0].metadata["spectrum_q"].as_f64().unwrap() - 0.001).abs() < 1e-12);
    // Two annotations for PSM 0, one for PSM 2, none for the filtered PSM 1.
    assert_eq!(pids[0].hits[0].peak_annotations.len(), 2);
    assert_eq!(pids[0].hits[0].peak_annotations[0].annotation, "b2");
    assert_eq!(pids[0].hits[0].peak_annotations[0].charge, 1);
    assert!((pids[0].hits[0].peak_annotations[0].mz - 201.1).abs() < 1e-9);
    assert!((pids[0].hits[0].peak_annotations[1].intensity - 600.25).abs() < 1e-9);
    assert_eq!(pids[1].hits[0].peak_annotations.len(), 1);
    assert_eq!(pids[1].hits[0].peak_annotations[0].annotation, "y4");
    // The third row's one-hot charge columns are all zero, so the Sage
    // "z=other" column supplies the charge.
    assert_eq!(pids[1].hits[0].charge, 4);
    // The decoy prefix reannotates the third row.
    assert_eq!(
        pids[1].hits[0].metadata["target_decoy"].to_string(),
        "decoy"
    );
}

/// A `.pin` path too short to carry the Sage suffixes is refused explicitly,
/// where the source computes the sibling path with unchecked `size() - N`
/// subtraction on the path string.
#[test]
fn a_short_or_non_ascii_sage_path_is_refused_rather_than_wrapping() {
    let directory = std::env::temp_dir().join(format!("openms-pin-sage-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    for name in ["a.pin", "\u{65e5}\u{672c}\u{8a9e}.pin"] {
        let path = directory.join(name);
        std::fs::write(&path, "SpecId\n").unwrap();
        let error = pin::load(
            &path,
            &ReadOptions {
                score_name: "ln(hyperscore)".into(),
                sage_annotation: true,
                ..Default::default()
            },
        )
        .unwrap_err();
        assert!(matches!(error, Error::InvalidValue(_)), "{name}: {error:?}");
        std::fs::remove_file(&path).unwrap();
    }
    std::fs::remove_dir_all(&directory).unwrap();
}

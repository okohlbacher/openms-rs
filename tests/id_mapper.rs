// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `ANALYSIS/ID/IDMapper` against the upstream class test
//! (`IDMapper_test.cpp` at core `bc9cc12`) and the retained outputs of the six
//! `TOPP_IDMapper_*` tests (test-data `0cb15f2`). See
//! `tests/data/id_mapper_provenance.json`.
#![cfg(all(
    feature = "idxml",
    feature = "featurexml",
    feature = "consensusxml",
    feature = "mzml"
))]

use openms::analysis::id::id_mapper::{
    IDMapper, IDMapperLimits, Measure, MzReference, SpectraAnnotation,
    UNKNOWN_SEARCH_RUN_IDENTIFIER,
};
use openms::chemistry::AASequence;
use openms::format::{consensusxml, featurexml, idxml, mzml};
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::kernel::{ConsensusFeature, ConsensusMap, MSExperiment, MSSpectrum, Precursor};
use openms::param::ParamValue;
use openms::{Error, Result};

fn data(name: &str) -> String {
    format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))
}

fn seq(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}

fn mapper(rt: f64, mz: f64, measure: Measure, ignore_charge: bool) -> IDMapper {
    IDMapper {
        rt_tolerance: rt,
        mz_tolerance: mz,
        measure,
        ignore_charge,
        ..IDMapper::default()
    }
}

fn spectrum(rt: f64, precursor_mz: f64, native_id: &str) -> MSSpectrum {
    MSSpectrum {
        rt,
        native_id: native_id.into(),
        precursors: vec![Precursor {
            mz: precursor_mz,
            ..Precursor::default()
        }],
        ..MSSpectrum::default()
    }
}

// ---------------------------------------------------------------- parameters

#[test]
fn parameters_follow_the_source_defaults() {
    let defaults = IDMapper::defaults().unwrap();
    assert_eq!(
        defaults.value("rt_tolerance").unwrap(),
        &ParamValue::Float(5.0)
    );
    assert_eq!(
        defaults.value("mz_tolerance").unwrap(),
        &ParamValue::Float(20.0)
    );
    assert_eq!(
        defaults.value("mz_measure").unwrap().as_str().unwrap(),
        "ppm"
    );
    assert_eq!(
        defaults.value("mz_reference").unwrap().as_str().unwrap(),
        "precursor"
    );
    assert_eq!(
        defaults.value("ignore_charge").unwrap().as_str().unwrap(),
        "false"
    );
    assert_eq!(defaults.valid_strings("mz_measure").unwrap(), ["ppm", "Da"]);
    assert_eq!(
        defaults.valid_strings("mz_reference").unwrap(),
        ["precursor", "peptide"]
    );
    assert_eq!(defaults.entry("rt_tolerance").unwrap().min_float, 0.0);
    assert_eq!(defaults.entry("mz_tolerance").unwrap().min_float, 0.0);
    assert_eq!(
        IDMapper::from_param(&defaults).unwrap(),
        IDMapper::default()
    );
}

/// Class test sections `IDMapper(const IDMapper&)` and `operator=`: the
/// parameters survive a copy.
#[test]
fn parameters_round_trip() {
    let mut p = IDMapper::defaults().unwrap();
    p.set_value("rt_tolerance", ParamValue::Float(0.5), "", &[])
        .unwrap();
    p.set_value("mz_tolerance", ParamValue::Float(0.05), "", &[])
        .unwrap();
    p.set_value("mz_measure", ParamValue::String("ppm".into()), "", &[])
        .unwrap();
    let mapper = IDMapper::from_param(&p).unwrap();
    let copy = mapper.clone();
    assert_eq!(
        copy.to_param().unwrap().value("rt_tolerance").unwrap(),
        &ParamValue::Float(0.5)
    );
    assert_eq!(
        copy.to_param().unwrap().value("mz_tolerance").unwrap(),
        &ParamValue::Float(0.05)
    );
    assert_eq!(copy.rt_tolerance, 0.5);
    assert_eq!(copy.measure, Measure::Ppm);
}

#[test]
fn parameters_are_checked_and_unknown_ones_warned() {
    let mut p = IDMapper::defaults().unwrap();
    p.set_value(
        "consensus:use_subelements",
        ParamValue::String("true".into()),
        "",
        &[],
    )
    .unwrap();
    p.set_value(
        "mz_reference",
        ParamValue::String("peptide".into()),
        "",
        &[],
    )
    .unwrap();
    p.set_value("ignore_charge", ParamValue::String("true".into()), "", &[])
        .unwrap();
    let (mapper, warnings) = IDMapper::from_param_with_warnings(&p).unwrap();
    assert_eq!(mapper.mz_reference, MzReference::Peptide);
    assert!(mapper.ignore_charge);
    assert!(!warnings.is_empty());

    let mut bad = IDMapper::defaults().unwrap();
    bad.set_value("rt_tolerance", ParamValue::Float(-1.0), "", &[])
        .unwrap();
    assert!(IDMapper::from_param(&bad).is_err());
    let mut bad = IDMapper::defaults().unwrap();
    bad.set_value("mz_measure", ParamValue::String("Th".into()), "", &[])
        .unwrap();
    assert!(IDMapper::from_param(&bad).is_err());
}

/// Class test `[EXTRA] getAbsoluteMZTolerance_`.
#[test]
fn absolute_mz_tolerance_matches_class_test() {
    let mut m = IDMapper {
        mz_tolerance: 1.0,
        ..IDMapper::default()
    };
    assert!((m.absolute_mz_tolerance(1000.0) - 0.001).abs() < 1e-12);
    m.mz_tolerance = 3.0;
    assert!((m.absolute_mz_tolerance(1000.0) - 0.003).abs() < 1e-12);
    m.measure = Measure::Da;
    assert_eq!(m.absolute_mz_tolerance(1000.0), 3.0);
}

/// Class test `[EXTRA] isMatch_`.
#[test]
fn is_match_matches_class_test() {
    let mut m = IDMapper::default();
    assert!(m.is_match(1.0, 1000.0, 1000.001));
    m.mz_tolerance = 3.0;
    assert!(m.is_match(4.0, 1000.0, 1000.0028));
    assert!(!m.is_match(4.0, 1000.0, 1000.004));
    assert!(m.is_match(4.0, 1000.0, 999.9972));
    assert!(!m.is_match(4.0, 1000.0, 999.996));
    m.measure = Measure::Da;
    assert!(m.is_match(5.0, 999.0, 1002.0));
    assert!(!m.is_match(5.0, 999.0, 1002.1));
    // Native: a zero reference and an RT beyond tolerance never match.
    assert!(!IDMapper::default().is_match(0.0, 0.0, 0.0));
    assert!(!m.is_match(5.5, 999.0, 999.0));
}

// ---------------------------------------------------------------- spectra

fn three_spectra(native: bool) -> MSExperiment {
    let ids = if native {
        ["spectrum=1234", "spectrum=6666", "spectrum=4321"]
    } else {
        ["", "", ""]
    };
    MSExperiment {
        spectra: vec![
            spectrum(60.0, 0.0, ids[0]),
            spectrum(181.0, 20.0, ids[1]),
            spectrum(120.0001, 11.0, ids[2]),
        ],
        ..MSExperiment::default()
    }
}

/// Class test section `annotate(AnnotatedMSRun&, ...)`: RT mapping and
/// native-ID mapping.
#[test]
fn annotate_spectra_matches_class_test() -> Result<()> {
    let doc1 = idxml::load(data("IDMapper_1.idXML"))?;
    let doc7 = idxml::load(data("IDMapper_7.idXML"))?;
    for doc in [&doc1, &doc7] {
        assert_eq!(doc.peptide_identifications.len(), 3);
        let sizes: Vec<_> = doc
            .peptide_identifications
            .iter()
            .map(|p| p.hits.len())
            .collect();
        assert_eq!(sizes, [2, 1, 2]);
        assert_eq!(doc.protein_identifications.len(), 1);
        assert_eq!(doc.protein_identifications[0].hits.len(), 2);
    }
    let m = mapper(0.5, 0.05, Measure::Da, true);
    for (doc, native) in [(&doc1, false), (&doc7, true)] {
        let experiment = three_spectra(native);
        let mut run = SpectraAnnotation::default();
        m.annotate_spectra(
            &experiment,
            &mut run,
            &doc.peptide_identifications,
            &doc.protein_identifications,
            false,
            false,
        )?;
        assert_eq!(run.protein_identifications.len(), 1);
        let accessions: Vec<_> = run.protein_identifications[0]
            .hits
            .iter()
            .map(|h| h.accession.as_str())
            .collect();
        assert_eq!(accessions, ["ABCDE", "FGHIJ"]);
        let ids = &run.peptide_identifications;
        assert_eq!(ids.len(), 3);
        assert_eq!(ids[0].hits.len(), 2);
        assert_eq!(ids[0].hits[0].sequence, seq("LHASGITVTEIPVTATNFK"));
        assert_eq!(ids[0].hits[1].sequence, seq("MRSLGYVAVISAVATDTDK"));
        assert_eq!(ids[1].hits.len(), 0);
        assert_eq!(ids[2].hits.len(), 1);
        assert_eq!(ids[2].hits[0].sequence, seq("HSKLSAK"));
    }
    Ok(())
}

/// Semantics of the commented-out class test for the FeatureMap overload
/// (`IDMapper_test.cpp:85-162`, disabled upstream): hits of every feature
/// identification, with missing coordinates taken from the feature. The
/// disabled test's own expectations use an older one-ID-per-spectrum layout
/// and are not transcribed; these are derived from the source body.
#[test]
fn annotate_spectra_from_features_fills_coordinates() -> Result<()> {
    use openms::kernel::{Feature, FeatureMap};
    let mut feature = Feature::new(9.0, 900.0, 1.0);
    let mut pid = PeptideIdentification {
        identifier: "myID".into(),
        hits: vec![PeptideHit::default(); 4],
        ..PeptideIdentification::default()
    };
    feature.peptide_identifications.push(pid.clone());
    pid.mz = Some(800.0);
    pid.rt = Some(9.05);
    feature.peptide_identifications.push(pid);
    let mut fm = FeatureMap::from_features(vec![feature]);
    fm.protein_identifications = vec![Default::default(); 2];
    let experiment = MSExperiment {
        spectra: vec![
            spectrum(8.9, 0.0, ""),
            spectrum(9.1, 20.0, ""),
            spectrum(12.0, 11.0, ""),
        ],
        ..MSExperiment::default()
    };
    let m = mapper(0.3, 0.05, Measure::Da, true);
    let mut run = SpectraAnnotation::default();
    let stats = m.annotate_spectra_from_features(&experiment, &mut run, &fm, true, true)?;
    assert_eq!(run.protein_identifications.len(), 2);
    // map_ms1: RT alone decides; both IDs (RT 9.0 and 9.05) lie strictly within
    // 0.3 s of 8.9 and 9.1, so each of these spectra receives 4 + 4 hits.
    assert_eq!(run.peptide_identifications[0].hits.len(), 8);
    assert_eq!(run.peptide_identifications[1].hits.len(), 8);
    assert_eq!(run.peptide_identifications[2].hits.len(), 0);
    assert_eq!(stats.assigned, 2);
    // Without MS1 mapping, no precursor m/z (0, 20, 11) is near 800 or 900.
    let stats = m.annotate_spectra_from_features(&experiment, &mut run, &fm, true, false)?;
    assert_eq!(run.protein_identifications.len(), 2);
    assert!(
        run.peptide_identifications
            .iter()
            .all(|p| p.hits.is_empty())
    );
    assert_eq!(stats.assigned, 0);
    assert_eq!(stats.unassigned, 2);
    Ok(())
}

#[test]
fn annotate_spectra_edge_cases() -> Result<()> {
    let m = mapper(0.5, 0.05, Measure::Da, true);
    let experiment = three_spectra(false);
    let mut run = SpectraAnnotation {
        peptide_identifications: vec![PeptideIdentification::default()],
        protein_identifications: vec![Default::default()],
    };
    // No identifications: nothing but the optional clearing happens.
    m.annotate_spectra(&experiment, &mut run, &[], &[], false, false)?;
    assert_eq!(run.peptide_identifications.len(), 1);
    m.annotate_spectra(&experiment, &mut run, &[], &[], true, false)?;
    assert!(run.peptide_identifications.is_empty());
    assert!(run.protein_identifications.is_empty());
    // A missing RT is the source MissingInformation; nothing is changed.
    let before = run.clone();
    let missing = PeptideIdentification {
        mz: Some(1.0),
        ..PeptideIdentification::default()
    };
    assert!(matches!(
        m.annotate_spectra(&experiment, &mut run, &[missing], &[], true, false),
        Err(Error::MissingInformation(_))
    ));
    assert_eq!(run, before);
    // Zero RT tolerance: the strict source window never matches by RT.
    let doc = idxml::load(data("IDMapper_1.idXML"))?;
    let strict = mapper(0.0, 1e9, Measure::Da, true);
    let stats = strict.annotate_spectra(
        &experiment,
        &mut run,
        &doc.peptide_identifications,
        &[],
        true,
        true,
    )?;
    assert_eq!(stats.assigned, 0);
    Ok(())
}

// ---------------------------------------------------------------- features

fn sequences(ids: &[PeptideIdentification]) -> Vec<String> {
    ids.iter().map(|p| p.hits[0].sequence.to_string()).collect()
}

/// Class test section `annotate(FeatureMap&, ...)`.
#[test]
fn annotate_features_matches_class_test() -> Result<()> {
    let doc = idxml::load(data("IDMapper_2.idXML"))?;
    let (ids, proteins) = (&doc.peptide_identifications, &doc.protein_identifications);
    let none = MSExperiment::default();

    // Convex hulls.
    let mut fm = featurexml::load(data("IDMapper_2.featureXML"))?;
    mapper(0.0, 0.0, Measure::Da, true)
        .annotate_features(&mut fm, ids, proteins, false, false, &none)?;
    assert_eq!(fm.protein_identifications.len(), 1);
    let acc: Vec<_> = fm.protein_identifications[0]
        .hits
        .iter()
        .map(|h| h.accession.as_str())
        .collect();
    assert_eq!(acc, ["ABCDE", "FGHIJ"]);
    let pids = &fm.features[0].peptide_identifications;
    assert_eq!(pids.len(), 7);
    assert!(pids.iter().all(|p| p.hits.len() == 1));
    assert_eq!(sequences(pids), ["A", "K", "C", "D", "E", "F", "I"]);
    assert_eq!(
        sequences(&fm.unassigned_peptide_identifications),
        ["G", "H", "L"]
    );

    // Centroids.
    let mut fm2 = featurexml::load(data("IDMapper_2.featureXML"))?;
    mapper(4.0, 1.5, Measure::Da, true)
        .annotate_features(&mut fm2, ids, proteins, true, true, &none)?;
    assert_eq!(fm2.protein_identifications.len(), 1);
    assert_eq!(
        sequences(&fm2.features[0].peptide_identifications),
        ["A", "K"]
    );
    assert_eq!(
        sequences(&fm2.unassigned_peptide_identifications),
        ["C", "D", "E", "F", "G", "H", "I", "L"]
    );

    // Charge-specific matching.
    let mut fm = featurexml::load(data("IDMapper_2.featureXML"))?;
    mapper(0.0, 0.0, Measure::Da, false)
        .annotate_features(&mut fm, ids, proteins, false, false, &none)?;
    assert_eq!(fm.protein_identifications.len(), 1);
    assert_eq!(
        sequences(&fm.features[0].peptide_identifications),
        ["A", "K", "C"]
    );
    assert_eq!(fm.unassigned_peptide_identifications.len(), 7);

    // ppm.
    let doc4 = idxml::load(data("IDMapper_4.idXML"))?;
    let mut fm_ppm = featurexml::load(data("IDMapper_4.featureXML"))?;
    mapper(4.0, 3.0, Measure::Ppm, true).annotate_features(
        &mut fm_ppm,
        &doc4.peptide_identifications,
        &doc4.protein_identifications,
        false,
        false,
        &none,
    )?;
    let f = &fm_ppm.features;
    assert_eq!(f[0].peptide_identifications.len(), 1);
    assert_eq!(f[0].peptide_identifications[0].hits.len(), 2);
    assert_eq!(
        f[0].peptide_identifications[0].hits[0].sequence,
        seq("LHASGITVTEIPVTATNFK")
    );
    assert_eq!(f[1].peptide_identifications.len(), 0);
    assert_eq!(f[2].peptide_identifications.len(), 1);
    assert_eq!(f[2].peptide_identifications[0].hits.len(), 1);
    assert_eq!(
        f[2].peptide_identifications[0].hits[0].sequence,
        seq("HSKLSAK")
    );
    assert_eq!(f[3].peptide_identifications.len(), 0);
    assert_eq!(f[4].peptide_identifications.len(), 1);
    assert_eq!(f[4].peptide_identifications[0].hits.len(), 2);
    assert_eq!(
        f[4].peptide_identifications[0].hits[0].sequence,
        seq("RASNSPQDPQSATAHSFR")
    );
    assert_eq!(f[5].peptide_identifications.len(), 0);
    let un = &fm_ppm.unassigned_peptide_identifications;
    assert_eq!(un.len(), 2);
    assert_eq!(un[0].hits[0].sequence, seq("DEAD"));
    assert_eq!(un[0].hits[1].sequence, seq("DEADA"));
    assert_eq!(un[1].hits[0].sequence, seq("DEADAA"));
    assert_eq!(un[1].hits[1].sequence, seq("DEADAAA"));
    Ok(())
}

#[test]
fn annotate_features_is_atomic_and_bounded() -> Result<()> {
    let doc = idxml::load(data("IDMapper_2.idXML"))?;
    let mut fm = featurexml::load(data("IDMapper_2.featureXML"))?;
    let before = fm.clone();
    let mut ids = doc.peptide_identifications.clone();
    ids.last_mut().unwrap().rt = None;
    let m = mapper(0.0, 0.0, Measure::Da, true);
    assert!(matches!(
        m.annotate_features(
            &mut fm,
            &ids,
            &doc.protein_identifications,
            false,
            false,
            &MSExperiment::default()
        ),
        Err(Error::MissingInformation(_))
    ));
    assert_eq!(fm, before);
    let tight = IDMapper {
        limits: IDMapperLimits {
            max_rt_bins: 1,
            ..IDMapperLimits::default()
        },
        ..m.clone()
    };
    assert!(matches!(
        tight.annotate_features(
            &mut fm,
            &doc.peptide_identifications,
            &[],
            false,
            false,
            &MSExperiment::default()
        ),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(fm, before);
    // An empty map leaves every identification unassigned, with a warning.
    let mut empty = openms::kernel::FeatureMap::new();
    let stats = m.annotate_features(
        &mut empty,
        &doc.peptide_identifications,
        &[],
        false,
        false,
        &MSExperiment::default(),
    )?;
    assert_eq!(
        empty.unassigned_peptide_identifications.len(),
        doc.peptide_identifications.len()
    );
    assert_eq!(stats.warnings.len(), 1);
    Ok(())
}

// ---------------------------------------------------------------- consensus

fn close(a: Option<f64>, b: Option<f64>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => (a - b).abs() <= 1e-9 * a.abs().max(1.0),
        (None, None) => true,
        _ => false,
    }
}

/// Compare the identifications the mapper produced with those a C++ build
/// stored: coordinates, hits (sequence and charge) and the meta values the
/// mapper writes. Run identifiers are not compared: the readers name runs
/// differently from the C++ writer.
fn same_ids(context: &str, got: &[PeptideIdentification], want: &[PeptideIdentification]) {
    assert_eq!(got.len(), want.len(), "{context}: identification count");
    for (k, (g, w)) in got.iter().zip(want).enumerate() {
        assert!(
            close(g.rt, w.rt),
            "{context} #{k}: RT {:?} vs {:?}",
            g.rt,
            w.rt
        );
        assert!(
            close(g.mz, w.mz),
            "{context} #{k}: m/z {:?} vs {:?}",
            g.mz,
            w.mz
        );
        assert_eq!(g.hits.len(), w.hits.len(), "{context} #{k}: hit count");
        for (gh, wh) in g.hits.iter().zip(&w.hits) {
            assert_eq!(gh.sequence, wh.sequence, "{context} #{k}: sequence");
            assert_eq!(gh.charge, wh.charge, "{context} #{k}: charge");
        }
        for key in ["map_index", "spectrum_index", "spectrum_reference"] {
            assert_eq!(
                g.metadata.get(key).map(ToString::to_string),
                w.metadata.get(key).map(ToString::to_string),
                "{context} #{k}: meta value {key}"
            );
        }
    }
}

fn same_consensus(context: &str, got: &ConsensusMap, want: &ConsensusMap) {
    assert_eq!(
        got.features.len(),
        want.features.len(),
        "{context}: feature count"
    );
    for (i, (g, w)) in got.features.iter().zip(&want.features).enumerate() {
        same_ids(
            &format!("{context} feature {i}"),
            &g.peptide_identifications,
            &w.peptide_identifications,
        );
    }
    same_ids(
        &format!("{context} unassigned"),
        &got.unassigned_peptide_identifications,
        &want.unassigned_peptide_identifications,
    );
}

/// Class test section `annotate(ConsensusMap&, ...)`, against the consensusXML
/// files the C++ class test compares with (`TEST_FILE_SIMILAR`).
#[test]
fn annotate_consensus_matches_class_test_outputs() -> Result<()> {
    let doc3 = idxml::load(data("IDMapper_3.idXML"))?;
    let doc5 = idxml::load(data("IDMapper_5.idXML"))?;
    let m = IDMapper {
        mz_tolerance: 0.01,
        measure: Measure::Da,
        ignore_charge: true,
        ..IDMapper::default()
    };
    let none = MSExperiment::default();
    for (subelements, expected) in [
        (false, "IDMapper_3_out1.consensusXML"),
        (true, "IDMapper_3_out2.consensusXML"),
    ] {
        let mut cm = consensusxml::load(data("IDMapper_3.consensusXML"))?;
        m.annotate_consensus(
            &mut cm,
            &doc3.peptide_identifications,
            &doc3.protein_identifications,
            subelements,
            false,
            &none,
        )?;
        same_consensus(expected, &cm, &consensusxml::load(data(expected))?);
        assert_eq!(
            cm.protein_identifications.len(),
            doc3.protein_identifications.len()
        );
    }
    let m5 = mapper(20.0, 20.0, Measure::Ppm, true);
    let mut cm = consensusxml::load(data("IDMapper_5.consensusXML"))?;
    m5.annotate_consensus(
        &mut cm,
        &doc5.peptide_identifications,
        &doc5.protein_identifications,
        true,
        true,
        &none,
    )?;
    same_consensus(
        "IDMapper_5_out1",
        &cm,
        &consensusxml::load(data("IDMapper_5_out1.consensusXML"))?,
    );
    Ok(())
}

/// Class test: charge-specific matching on a one-feature map.
#[test]
fn annotate_consensus_checks_charges() -> Result<()> {
    let doc3 = idxml::load(data("IDMapper_3.idXML"))?;
    let ids = &doc3.peptide_identifications;
    let mut m = IDMapper {
        mz_tolerance: 0.01,
        measure: Measure::Da,
        ignore_charge: true,
        ..IDMapper::default()
    };
    let mut feature = ConsensusFeature::new();
    feature.rt = 4101.48;
    feature.mz = 117.1;
    feature.charge = 2;
    let mut cm = ConsensusMap::from_features(vec![feature]);
    m.annotate_consensus(
        &mut cm,
        ids,
        &doc3.protein_identifications,
        false,
        false,
        &MSExperiment::default(),
    )?;
    assert_eq!(cm.features[0].peptide_identifications.len(), 1);
    assert_eq!(
        cm.features[0].peptide_identifications[0].hits[0].sequence,
        seq("ACSF")
    );
    assert_eq!(cm.unassigned_peptide_identifications.len(), ids.len() - 1);

    cm.features[0].peptide_identifications.clear();
    cm.unassigned_peptide_identifications.clear();
    m.ignore_charge = false;
    m.annotate_consensus(
        &mut cm,
        ids,
        &doc3.protein_identifications,
        false,
        false,
        &MSExperiment::default(),
    )?;
    assert_eq!(cm.features[0].peptide_identifications.len(), 0);
    assert_eq!(cm.unassigned_peptide_identifications.len(), ids.len());
    Ok(())
}

/// Class test: annotation of precursors without identification.
#[test]
fn annotate_consensus_unidentified_precursors_match_class_test_outputs() -> Result<()> {
    let m = IDMapper {
        mz_tolerance: 0.01,
        measure: Measure::Da,
        ignore_charge: true,
        ..IDMapper::default()
    };
    let rts = [
        306.58, 306.58, 312.738, 312.738, 3112.53, 3840.95, 3849.22, 3870.67, 3880.9, 3892.26,
    ];
    let exact = [
        426.849, 405.85, 506.815, 484.83, 496.244, 430.212, 446.081, 453.233, 400.172, 437.227,
    ];
    let mismatch = [
        426.85899, 405.0, 506.815, 484.85, 496.244, 430.0, 446.081, 453.0, 400.172, 437.239,
    ];
    let experiment = |mzs: &[f64], rts: &[f64]| MSExperiment {
        spectra: mzs
            .iter()
            .zip(rts)
            .map(|(&mz, &rt)| spectrum(rt, mz, ""))
            .collect(),
        ..MSExperiment::default()
    };
    let multiple_rts = [306.58, 305.58, 307.58, 304.58, 308.58];
    let cases = [
        (experiment(&exact, &rts), "IDMapper_6_out1.consensusXML"),
        (experiment(&mismatch, &rts), "IDMapper_6_out2.consensusXML"),
        (
            experiment(&[426.849; 5], &multiple_rts),
            "IDMapper_6_out3.consensusXML",
        ),
    ];
    for (spectra, expected) in cases {
        let mut cm = consensusxml::load(data("IDMapper_3.consensusXML"))?;
        m.annotate_consensus(&mut cm, &[], &[], false, false, &spectra)?;
        same_consensus(expected, &cm, &consensusxml::load(data(expected))?);
        // No run existed, so the unknown run was added and every empty ID uses it.
        assert_eq!(cm.protein_identifications.len(), 1);
        assert_eq!(
            cm.protein_identifications[0].identifier,
            UNKNOWN_SEARCH_RUN_IDENTIFIER
        );
        assert!(
            cm.features
                .iter()
                .flat_map(|f| &f.peptide_identifications)
                .all(|p| p.identifier == UNKNOWN_SEARCH_RUN_IDENTIFIER)
        );
    }
    Ok(())
}

// ---------------------------------------------------------------- TOPP retained outputs

/// The TOPP tool's option defaults: `mz_reference` defaults to its second
/// valid string, `peptide`, and `feature:use_centroid_mz` to `true`
/// (`IDMapper.cpp` at topp `174b576`).
fn topp_mapper(
    rt: f64,
    mz: f64,
    measure: Measure,
    reference: MzReference,
    ignore_charge: bool,
) -> IDMapper {
    IDMapper {
        rt_tolerance: rt,
        mz_tolerance: mz,
        measure,
        mz_reference: reference,
        ignore_charge,
        ..IDMapper::default()
    }
}

fn topp_features(
    number: u32,
    mapper: &IDMapper,
    use_centroid_mz: bool,
    with_spectra: bool,
) -> Result<()> {
    let doc = idxml::load(data(&format!("IDMapper_{number}_input.idXML")))?;
    let mut map = featurexml::load(data(&format!("IDMapper_{number}_input.featureXML")))?;
    let spectra = if with_spectra {
        mzml::load(data(&format!("IDMapper_{number}_input.mzML")))?
    } else {
        MSExperiment::default()
    };
    mapper.annotate_features(
        &mut map,
        &doc.peptide_identifications,
        &doc.protein_identifications,
        false,
        use_centroid_mz,
        &spectra,
    )?;
    let want = featurexml::load(data(&format!("IDMapper_{number}_output.featureXML")))?;
    assert_eq!(map.features.len(), want.features.len());
    for (i, (g, w)) in map.features.iter().zip(&want.features).enumerate() {
        same_ids(
            &format!("TOPP_IDMapper_{number} feature {i}"),
            &g.peptide_identifications,
            &w.peptide_identifications,
        );
    }
    same_ids(
        &format!("TOPP_IDMapper_{number} unassigned"),
        &map.unassigned_peptide_identifications,
        &want.unassigned_peptide_identifications,
    );
    Ok(())
}

fn topp_consensus(number: u32, mapper: &IDMapper) -> Result<()> {
    let doc = idxml::load(data(&format!("IDMapper_{number}_input.idXML")))?;
    let mut map = consensusxml::load(data(&format!("IDMapper_{number}_input.consensusXML")))?;
    mapper.annotate_consensus(
        &mut map,
        &doc.peptide_identifications,
        &doc.protein_identifications,
        false,
        true,
        &MSExperiment::default(),
    )?;
    map.sort_peptide_identifications_by_map_index()?;
    let want = consensusxml::load(data(&format!("IDMapper_{number}_output.consensusXML")))?;
    same_consensus(&format!("TOPP_IDMapper_{number}"), &map, &want);
    Ok(())
}

/// `TOPP_IDMapper_1`: `-mz_measure Da -mz_tolerance 1 -ignore_charge
/// -mz_reference precursor -feature:use_centroid_mz false`.
#[test]
fn topp_idmapper_1() -> Result<()> {
    topp_features(
        1,
        &topp_mapper(5.0, 1.0, Measure::Da, MzReference::Precursor, true),
        false,
        false,
    )
}

/// `TOPP_IDMapper_2`: `-mz_measure Da -mz_tolerance 1 -mz_reference precursor`.
#[test]
fn topp_idmapper_2() -> Result<()> {
    topp_consensus(
        2,
        &topp_mapper(5.0, 1.0, Measure::Da, MzReference::Precursor, false),
    )
}

/// `TOPP_IDMapper_3`: `-mz_measure ppm -mz_tolerance 3 -rt_tolerance 4
/// -ignore_charge -mz_reference precursor`.
#[test]
fn topp_idmapper_3() -> Result<()> {
    topp_features(
        3,
        &topp_mapper(4.0, 3.0, Measure::Ppm, MzReference::Precursor, true),
        true,
        false,
    )
}

/// `TOPP_IDMapper_4`: `-spectra:in … -mz_measure ppm -mz_tolerance 20
/// -rt_tolerance 10`.
#[test]
fn topp_idmapper_4() -> Result<()> {
    topp_features(
        4,
        &topp_mapper(10.0, 20.0, Measure::Ppm, MzReference::Peptide, false),
        true,
        true,
    )
}

/// `TOPP_IDMapper_5`: `-spectra:in … -mz_measure ppm -mz_tolerance 10
/// -rt_tolerance 20 -feature:use_centroid_mz false`.
#[test]
fn topp_idmapper_5() -> Result<()> {
    topp_features(
        5,
        &topp_mapper(20.0, 10.0, Measure::Ppm, MzReference::Peptide, false),
        false,
        true,
    )
}

/// `TOPP_IDMapper_6`: isobaric (TMT) consensus map, mapped by scan ID.
#[test]
fn topp_idmapper_6() -> Result<()> {
    topp_consensus(
        6,
        &topp_mapper(5.0, 20.0, Measure::Ppm, MzReference::Peptide, false),
    )
}

#[test]
fn isobaric_mapping_requires_spectrum_references() -> Result<()> {
    let mut map = consensusxml::load(data("IDMapper_6_input.consensusXML"))?;
    let before = map.clone();
    let mut doc = idxml::load(data("IDMapper_6_input.idXML"))?;
    for id in &mut doc.peptide_identifications {
        id.metadata.remove("spectrum_reference");
    }
    let m = IDMapper::default();
    assert!(matches!(
        m.annotate_consensus(
            &mut map,
            &doc.peptide_identifications,
            &doc.protein_identifications,
            false,
            false,
            &MSExperiment::default()
        ),
        Err(Error::MissingInformation(_))
    ));
    assert_eq!(map, before);
    Ok(())
}

#[test]
fn map_precursors_partitions_spectra() -> Result<()> {
    let doc = idxml::load(data("IDMapper_1.idXML"))?;
    let first = &doc.peptide_identifications[0];
    let mut experiment = MSExperiment::default();
    experiment.spectra.push(MSSpectrum::default()); // no precursor
    experiment
        .spectra
        .push(spectrum(first.rt.unwrap(), first.mz.unwrap(), "")); // identified
    experiment
        .spectra
        .push(spectrum(first.rt.unwrap() + 1.0, first.mz.unwrap(), "")); // off in RT
    let state = IDMapper::map_precursors_to_identifications_default(
        &experiment,
        &doc.peptide_identifications,
    )?;
    assert_eq!(state.no_precursors, [0]);
    assert_eq!(state.identified, [1]);
    assert_eq!(state.unidentified, [2]);
    Ok(())
}

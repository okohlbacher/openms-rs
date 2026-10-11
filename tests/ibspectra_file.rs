// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/ibspectra.rs`, the port of `FORMAT/IBSpectraFile.h`.
//!
//! Every `START_SECTION` of the upstream `IBSpectraFile_test.cpp` (core
//! `bc9cc12`) is reproduced. Its real example compares with the retained
//! `IBSpectraFile.ibspectra.csv` through `TEST_FILE_SIMILAR`, that is the
//! fuzzy comparator with the class-test defaults (relative 1 + 1e-5, absolute
//! 1e-5, `ClassTest.cpp:35-38`); it needs the `consensusxml` feature. Fixtures
//! and hashes are in `tests/data/ibspectra_file_provenance.json`.

use openms::Error;
use openms::chemistry::AASequence;
use openms::format::ibspectra::IBSpectraFile;
use openms::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use openms::kernel::{BaseFeature, ColumnHeader, ConsensusFeature, ConsensusMap, FeatureHandle};
use std::path::{Path, PathBuf};

// Only the consensusxml-gated tests read a fixture; without the feature this would be dead code.
#[cfg_attr(not(feature = "consensusxml"), allow(dead_code))]
fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/ibspectra_file")
        .join(name)
}

// START_SECTION((IBSpectraFile()))
#[test]
fn default_construction() {
    let file = IBSpectraFile::new();
    assert_eq!(file, IBSpectraFile);
}

// START_SECTION((IBSpectraFile(const IBSpectraFile& other)))
// START_SECTION((IBSpectraFile& operator=(const IBSpectraFile& rhs)))
#[test]
fn copy_and_assignment() {
    let original = IBSpectraFile::new();
    let copy = original;
    let mut assigned = IBSpectraFile::new();
    assert_eq!(assigned, copy);
    assigned = original;
    assert_eq!(assigned, original);
}

// START_SECTION((void store(const std::string& filename, const ConsensusMap& cm)))
#[test]
fn store_refuses_maps_without_isobaric_data() {
    let mut cm_no_ms2quant = ConsensusMap::new();
    cm_no_ms2quant.experiment_type = "labeled_MS1".into();
    match IBSpectraFile::new().store("not-a-file-name", &cm_no_ms2quant) {
        Err(Error::InvalidValue(message)) => assert_eq!(
            message,
            "Given ConsensusMap does not hold any isobaric quantification data."
        ),
        other => panic!("{other:?}"),
    }

    let mut cm_wrong_channel_count = ConsensusMap::new();
    cm_wrong_channel_count.experiment_type = "labeled_MS2".into();
    for index in 0..3 {
        cm_wrong_channel_count
            .column_headers
            .insert(index, ColumnHeader::default());
    }
    match IBSpectraFile::new().store("not-a-file-name", &cm_wrong_channel_count) {
        Err(Error::InvalidValue(message)) => assert_eq!(
            message,
            "Could not guess isobaric quantification data from ConsensusMap due to non-matching number of input maps."
        ),
        other => panic!("{other:?}"),
    }
    assert!(!Path::new("not-a-file-name").exists());
}

#[cfg(feature = "consensusxml")]
#[test]
fn store_matches_the_retained_example() {
    use openms::concept::fuzzy_string_comparator::{FuzzyStringComparator, LogDestination};
    let cm = openms::format::consensusxml::load(data("IBSpectraFile.consensusXML")).unwrap();
    let dir = openms::system::file::TempDir::new(false).unwrap();
    let tmp_filename = dir.path().join("out.ibspectra.csv");
    IBSpectraFile::new().store(&tmp_filename, &cm).unwrap();

    let mut comparator = FuzzyStringComparator::new();
    comparator.set_log_destination(LogDestination::Buffer);
    comparator.set_acceptable_relative(1.0 + 1e-5);
    comparator.set_acceptable_absolute(1e-5);
    let similar = comparator.compare_files(&tmp_filename, &data("IBSpectraFile.ibspectra.csv"));
    assert!(
        similar,
        "{}\n{}",
        String::from_utf8_lossy(comparator.log()),
        std::fs::read_to_string(&tmp_filename).unwrap()
    );
}

// ----------------------------------------------------------------------------
// Native regressions derived from IBSpectraFile.cpp
// ----------------------------------------------------------------------------

fn handle(map_index: u64, mz: f64, intensity: f32) -> FeatureHandle {
    let mut feature = BaseFeature::new(10.0, mz, intensity);
    feature.unique_id = map_index + 1;
    FeatureHandle::new(map_index, &feature)
}

fn four_plex() -> ConsensusMap {
    let mut cm = ConsensusMap::new();
    cm.experiment_type = "itraq".into();
    for index in 0..4 {
        cm.column_headers.insert(index, ColumnHeader::default());
    }
    cm
}

fn table(cm: &ConsensusMap) -> String {
    let mut out = Vec::new();
    IBSpectraFile::new().write(&mut out, cm).unwrap();
    String::from_utf8(out).unwrap()
}

#[test]
fn unidentified_features_and_channel_intensities() {
    let mut cm = four_plex();
    let mut cf = ConsensusFeature::new();
    cf.base = BaseFeature::new(600.5, 500.25, 1000.0);
    cf.base.unique_id = 99;
    cf.insert(handle(0, 114.11, 5.0)).unwrap();
    cf.insert(handle(1, 115.9, 0.000_001)).unwrap(); // below the threshold
    cf.insert(handle(2, 117.2, 7.0)).unwrap();
    cm.features.push(cf);
    let mut zero = ConsensusFeature::new();
    zero.base = BaseFeature::new(1.0, 2.0, 0.0);
    cm.features.push(zero);

    let text = table(&cm);
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(
        lines[0],
        "accession\tpeptide\tmodif\tcharge\ttheo.mass\texp.mass\tparent.intens\tretention.time\tspectrum\tsearch.engine\tX114_mass\tX115_mass\tX116_mass\tX117_mass\tX114_ions\tX115_ions\tX116_ions\tX117_ions"
    );
    // The zero-intensity feature writes nothing.
    assert_eq!(lines.len(), 2);
    let fields: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(
        &fields[..10],
        &[
            "UNIDENTIFIED_PROTEIN",
            "UNIDENTIFIED_PEPTIDE",
            "",
            "0",
            "-1.0",
            "500.25",
            "1000.0",
            "600.5",
            "99",
            "open-ms-generic"
        ]
    );
    // Channel 115 is below the threshold, channel 116 has no handle, and the
    // handle at 117.2 truncates to channel 117.
    assert_eq!(&fields[14..], &["5.0", "0.0", "0.0", "7.0"]);
}

#[test]
fn identified_features_write_one_line_per_known_accession() {
    let mut cm = four_plex();
    let mut proteins = ProteinIdentification::default();
    for accession in ["P2", "P1"] {
        proteins.hits.push(ProteinHit {
            accession: accession.into(),
            ..ProteinHit::default()
        });
    }
    cm.protein_identifications.push(proteins);

    let sequence = AASequence::parse("PEPTIDEK").unwrap();
    let mut hit = PeptideHit::new(1.0, 1, 2, sequence.clone()).unwrap();
    for accession in ["P2", "P1", "P9"] {
        hit.evidences.push(PeptideEvidence {
            protein_accession: accession.into(),
            ..PeptideEvidence::default()
        });
    }
    let mut identification = PeptideIdentification::default();
    identification.hits.push(hit);
    let mut cf = ConsensusFeature::new();
    cf.base = BaseFeature::new(1.0, 2.0, 3.0);
    cf.base.peptide_identifications.push(identification);
    cm.features.push(cf);

    let text = table(&cm);
    let lines: Vec<&str> = text.lines().collect();
    // P9 is not in the protein identification and is skipped; the others come
    // in sorted order.
    assert_eq!(lines.len(), 3);
    let theo = openms::format::file_info::text_format::to_str(
        sequence
            .mono_mass_for(openms::chemistry::PeptideFragmentType::Full, 2)
            .unwrap(),
    );
    for (line, accession) in lines[1..].iter().zip(["P1", "P2"]) {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[0], accession);
        assert_eq!(fields[1], "PEPTIDEK");
        // No N-terminal modification, then ":" for each of the eight residues.
        assert_eq!(fields[2], "::::::::");
        assert_eq!(fields[3], "2");
        assert_eq!(fields[4], theo);
    }
}

#[test]
fn a_peptide_identification_without_hits_is_refused() {
    let mut cm = four_plex();
    cm.protein_identifications
        .push(ProteinIdentification::default());
    let mut cf = ConsensusFeature::new();
    cf.base = BaseFeature::new(1.0, 2.0, 3.0);
    cf.base
        .peptide_identifications
        .push(PeptideIdentification::default());
    cm.features.push(cf);
    assert!(matches!(
        IBSpectraFile::new().write(Vec::new(), &cm),
        Err(Error::InvalidValue(_))
    ));
}

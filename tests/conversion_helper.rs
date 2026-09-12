// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `ConversionHelper_test.cpp` (3 sections) at Core SDK
//! `bc9cc12514c768385ce121d6ca4bb710fe1983c4`, plus native checks for the
//! conversions in `src/kernel/conversion_helper.rs`.
//!
//! `docs/CONVERSION_HELPER_SUPPORT.md` holds the section-to-test table. All
//! transcribed expectations are tier 3 (source review): the literals come from
//! the pinned class test, no C++ was built or executed.

use openms::concept::{HasUniqueId, UniqueIdGenerator};
use openms::identification::{PeptideIdentification, ProteinIdentification};
use openms::kernel::conversion_helper::MapConversion;
use openms::kernel::features::{ConsensusFeature, Feature, FeatureMap};
use openms::kernel::{MSChromatogram, MSExperiment, MSSpectrum, Peak1D};

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}
fn generator() -> UniqueIdGenerator {
    UniqueIdGenerator::from_seed(20_260_912)
}

/// `fm` of the first section (lines 28-36): three features whose retention
/// time, m/z and unique ID all derive from the loop index.
fn source_feature_map() -> FeatureMap {
    let mut map = FeatureMap::new();
    for i in 0..3u64 {
        let mut feature = Feature::new(i as f64 * 77.7, i as f64 + 100.35, 0.0);
        feature.base.unique_id = i * 33 + 17;
        map.features.push(feature);
    }
    map
}

/// `mse` of the second section (lines 62-78): three spectra of four peaks each,
/// with `mz = 10 * m + i + 100.35`, `intensity = 900 + 7 * m + 5 * i` and
/// `rt = m * 5`.
fn source_peak_map() -> MSExperiment {
    let mut experiment = MSExperiment::new();
    for m in 0..3u32 {
        let mut spectrum = MSSpectrum::default();
        for i in 0..4u32 {
            spectrum.peaks.push(Peak1D::new(
                f64::from(10 * m + i) + 100.35,
                (900 + 7 * m + 5 * i) as f32,
            ));
        }
        spectrum.rt = f64::from(m * 5);
        experiment.spectra.push(spectrum);
    }
    experiment
}

/// The `cm` fixture built at file scope before the third section (line 102):
/// `MapConversion::convert(33, mse, cm, 8)`.
///
/// The generator is threaded through rather than created here: the source draws
/// from one process-wide singleton, so two consecutive conversions cannot
/// produce the same container ID, and a caller-owned generator only matches
/// that when it is reused rather than re-seeded.
fn source_consensus_map(
    generator: &mut UniqueIdGenerator,
) -> openms::kernel::features::ConsensusMap {
    MapConversion::peak_map_to_consensus(33, &source_peak_map(), Some(8), generator).unwrap()
}

/// `static void convert(UInt64 const input_map_index, FeatureMap const&
/// input_map, ConsensusMap& output_map, Size n = -1)` (line 25).
#[test]
fn feature_map_to_consensus() {
    let fm = source_feature_map();
    let cm = MapConversion::feature_map_to_consensus(33, &fm, None).unwrap();

    assert_eq!(cm.features.len(), 3);
    assert_eq!(cm.column_headers[&33].size, 3);
    for i in 0..3usize {
        assert_eq!(cm.features[i].len(), 1);
        let handle = cm.features[i].handles()[0];
        assert_eq!(handle.map_index, 33);
        assert_eq!(handle.unique_id, i as u64 * 33 + 17);
        close(handle.rt, i as f64 * 77.7);
        close(handle.mz, i as f64 + 100.35);
    }

    // Truncated copy: the header size still reports the full input size.
    let cm = MapConversion::feature_map_to_consensus(33, &fm, Some(2)).unwrap();
    assert_eq!(cm.features.len(), 2);
    assert_eq!(cm.column_headers[&33].size, 3);
}

/// Native: the feature-map conversion carries the identification records,
/// stamps `map_index` on the copied peptide identifications and takes the
/// container unique ID from the input.
#[test]
fn feature_map_to_consensus_carries_records_and_stamps_map_index() {
    let mut fm = source_feature_map();
    fm.unique_id = 4242;
    fm.protein_identifications
        .resize(2, ProteinIdentification::default());
    fm.unassigned_peptide_identifications
        .resize(3, PeptideIdentification::default());
    fm.features[0]
        .base
        .peptide_identifications
        .push(PeptideIdentification::default());

    let cm = MapConversion::feature_map_to_consensus(7, &fm, None).unwrap();
    assert_eq!(cm.unique_id, 4242);
    assert_eq!(cm.protein_identifications.len(), 2);
    assert_eq!(cm.unassigned_peptide_identifications.len(), 3);
    assert_eq!(
        cm.features[0].base.peptide_identifications[0].metadata["map_index"]
            .as_i64()
            .unwrap(),
        7
    );
    // The input is untouched.
    assert!(
        !fm.features[0].base.peptide_identifications[0]
            .metadata
            .contains_key("map_index")
    );
}

/// `static void convert(UInt64 const input_map_index, PeakMap& input_map,
/// ConsensusMap& output_map, Size n = -1)` (line 80).
#[test]
fn peak_map_to_consensus() {
    let mse = source_peak_map();
    let cm = MapConversion::peak_map_to_consensus(33, &mse, Some(8), &mut generator()).unwrap();

    assert_eq!(cm.features.len(), 8);
    // Intensities run 900, 905, 910, 915, 907, 912, 917, 922, 914, 919, 924, 929;
    // the eight most intense end at 912.
    assert_eq!(cm.features.last().unwrap().intensity, 912.0);
    // Descending intensity order, one handle each, all tagged with the index.
    assert_eq!(cm.features[0].intensity, 929.0);
    for (position, feature) in cm.features.iter().enumerate() {
        assert_eq!(feature.len(), 1);
        assert_eq!(feature.handles()[0].map_index, 33);
        assert_eq!(feature.handles()[0].unique_id, position as u64);
    }
    // The header size reports what was written, unlike the feature-map overload.
    assert_eq!(cm.column_headers[&33].size, 8);
    // A peak map has no container unique ID, so a fresh one is drawn.
    assert!(cm.unique_id.has_valid_unique_id());
}

/// Native: `None` keeps every MS1 peak, and the clamp uses the number of points
/// actually collected rather than the source's `MSExperiment::getSize()`, which
/// counts MS2 spectra and chromatograms the conversion never reads.
#[test]
fn peak_map_to_consensus_clamps_against_collected_points() {
    let mut mse = source_peak_map();
    let all = MapConversion::peak_map_to_consensus(0, &mse, None, &mut generator()).unwrap();
    assert_eq!(all.features.len(), 12);

    let mut ms2 = MSSpectrum {
        ms_level: 2,
        ..MSSpectrum::default()
    };
    ms2.peaks.push(Peak1D::new(50.0, 1.0));
    mse.spectra.push(ms2);
    let mut chromatogram = MSChromatogram::default();
    chromatogram
        .peaks
        .push(openms::kernel::ChromatogramPeak::new(1.0, 1.0));
    mse.chromatograms.push(chromatogram);

    // Source: n = Size(-1) is clamped to getSize() == 14 and the partial sort
    // then runs past the end of the 12-element MS1 vector.
    let clamped = MapConversion::peak_map_to_consensus(0, &mse, None, &mut generator()).unwrap();
    assert_eq!(clamped.features.len(), 12);
    assert_eq!(clamped.column_headers[&0].size, 12);
    // An empty run yields an empty map rather than an out-of-range access.
    let empty =
        MapConversion::peak_map_to_consensus(0, &MSExperiment::new(), None, &mut generator())
            .unwrap();
    assert!(empty.features.is_empty());
    assert_eq!(empty.column_headers[&0].size, 0);
}

/// `static void convert(ConsensusMap const& input_map, const bool keep_uids,
/// FeatureMap& output_map)` (line 104).
#[test]
fn consensus_to_feature_map() {
    let mut generator = generator();
    let cm = source_consensus_map(&mut generator);
    let out_fm = MapConversion::consensus_to_feature_map(&cm, true, &mut generator).unwrap();

    assert_eq!(cm.unique_id, out_fm.unique_id);
    assert_eq!(
        cm.protein_identifications.len(),
        out_fm.protein_identifications.len()
    );
    assert_eq!(
        cm.unassigned_peptide_identifications.len(),
        out_fm.unassigned_peptide_identifications.len()
    );
    assert_eq!(cm.features.len(), out_fm.features.len());
    for i in 0..cm.features.len() {
        // The source compares the two elements; only the BaseFeature part of a
        // consensus feature survives the conversion, and that part is equal.
        assert_eq!(cm.features[i].base, out_fm.features[i].base);
    }

    let out_fm = MapConversion::consensus_to_feature_map(&cm, false, &mut generator).unwrap();
    assert_ne!(cm.unique_id, out_fm.unique_id);
    for i in 0..cm.features.len() {
        close(cm.features[i].rt, out_fm.features[i].rt);
        close(cm.features[i].mz, out_fm.features[i].mz);
        close(
            f64::from(cm.features[i].intensity),
            f64::from(out_fm.features[i].intensity),
        );
        assert_ne!(cm.features[i].unique_id, out_fm.features[i].unique_id);
    }
    // Fresh IDs are unique, so the resulting map validates.
    out_fm.validate().unwrap();
}

/// Native: the consensus conversion keeps the document identity and the
/// identification records, and drops the handles, ratios and hulls.
#[test]
fn consensus_to_feature_map_keeps_records_and_drops_handles() {
    let mut generator = generator();
    let mut cm = source_consensus_map(&mut generator);
    cm.identifier = "lsid".into();
    cm.loaded_file_path = "input.consensusXML".into();
    cm.loaded_file_type = openms::format::FileType::ConsensusXml;
    cm.protein_identifications
        .resize(2, ProteinIdentification::default());
    cm.unassigned_peptide_identifications
        .resize(1, PeptideIdentification::default());
    cm.features[0].base.quality = 0.5;
    cm.features[0].base.charge = 3;
    cm.features[0]
        .base
        .metadata
        .insert("note".into(), "kept".into());

    let out = MapConversion::consensus_to_feature_map(&cm, true, &mut generator).unwrap();
    assert_eq!(out.identifier, "lsid");
    assert_eq!(out.loaded_file_path, "input.consensusXML");
    assert_eq!(out.loaded_file_type, openms::format::FileType::ConsensusXml);
    assert_eq!(out.protein_identifications.len(), 2);
    assert_eq!(out.unassigned_peptide_identifications.len(), 1);
    assert_eq!(out.features[0].quality, 0.5);
    assert_eq!(out.features[0].charge, 3);
    assert_eq!(out.features[0].metadata["note"].as_str().unwrap(), "kept");
    // The consensus grouping and every Feature-only member are empty.
    assert!(out.features[0].convex_hulls.is_empty());
    assert!(out.features[0].subordinates.is_empty());
    assert_eq!(out.features[0].quality_rt, 0.0);
    assert_eq!(out.features[0].quality_mz, 0.0);
}

/// Native: the declared ceilings are checked and a rejected conversion
/// allocates nothing.
#[test]
fn conversion_ceilings() {
    assert_eq!(MapConversion::MAX_ITEMS, 10_000_000);
    let empty = MapConversion::feature_map_to_consensus(0, &FeatureMap::new(), Some(5)).unwrap();
    assert!(empty.features.is_empty());
    assert_eq!(empty.column_headers[&0].size, 0);
    // A non-finite coordinate is rejected before any element is produced; the
    // source copies it unchecked.
    let mut broken = FeatureMap::new();
    broken.features.push(Feature::new(f64::NAN, 1.0, 1.0));
    assert!(MapConversion::feature_map_to_consensus(0, &broken, None).is_err());

    // The consensus conversion accepts an empty map.
    let out = MapConversion::consensus_to_feature_map(
        &openms::kernel::features::ConsensusMap::new(),
        false,
        &mut generator(),
    )
    .unwrap();
    assert!(out.features.is_empty());
    assert!(ConsensusFeature::new().is_empty());
}

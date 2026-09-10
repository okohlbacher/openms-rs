// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::format::{FileHandler, FileType, file_handler::type_by_content};
use openms::{Error, MSExperiment, MSSpectrum, Peak1D};

struct Directory(std::path::PathBuf);
impl Directory {
    fn new() -> Self {
        static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let id = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("openms-rs-handler-{}-{id}", std::process::id()));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn experiment() -> MSExperiment {
    MSExperiment {
        spectra: vec![MSSpectrum {
            ms_level: 2,
            peaks: vec![Peak1D::new(100.5, 7.0)],
            ..Default::default()
        }],
        ..Default::default()
    }
}

#[test]
fn capabilities_distinguish_recognition_from_real_build_support() {
    assert_eq!(FileType::from_name("featureXML"), FileType::FeatureXml);
    assert!(!FileHandler::can_read_experiment(FileType::FeatureXml));
    assert!(!FileHandler::can_write_experiment(FileType::Raw));
    assert_eq!(
        FileHandler::can_read_experiment(FileType::MzMl),
        cfg!(feature = "mzml")
    );
    assert!(matches!(
        FileHandler::read_experiment(b"".as_slice(), FileType::Raw),
        Err(Error::Unsupported(_))
    ));
    assert!(FileHandler::can_read_experiment(FileType::Ms2));
    assert!(FileHandler::can_write_experiment(FileType::Ms2));
    assert!(
        !FileType::Ms2
            .source_properties()
            .contains(&openms::format::FileProperty::Writeable)
    );
}

#[test]
fn content_sniffing_finds_source_markers_without_claiming_validity() {
    for (text, expected) in [
        ("BEGIN IONS\nPEPMASS=500\nEND IONS\n", FileType::Mgf),
        ("# comment\n>p\nPEPTIDE\n", FileType::Fasta),
        ("<?xml version=\"1.0\"?>\n<mzML>", FileType::MzMl),
        (
            "<mzML>\n<cv id=\"IMS\" name=\"Imaging MS Ontology\"/>",
            FileType::ImzMl,
        ),
        ("<featureMap>", FileType::FeatureXml),
        ("500 2\n100 1\n200 2\n300 3\n400 4\n", FileType::Dta),
        (
            "#SEC\n1 100 1\n1 200 2\n2 100 3\n2 200 4\n",
            FileType::Dta2d,
        ),
        ("H\tCreationDate\tdate", FileType::Ms2),
        ("", FileType::Unknown),
    ] {
        assert_eq!(type_by_content(text.as_bytes()), expected, "{text}");
    }
    let mut text = b"<mzML>\n".to_vec();
    text.extend_from_slice(&vec![b' '; 65_536]);
    text.extend_from_slice(b"IMS:1000050");
    assert_eq!(type_by_content(&text), FileType::MzMl);
}

#[test]
fn dispatcher_roundtrip_and_unknown_suffix_detection() {
    let dir = Directory::new();
    let path = dir.0.join("peaks.data");
    FileHandler::store_experiment(&path, &experiment(), Some(FileType::Mgf)).unwrap();
    let copy = FileHandler::load_experiment(&path, &[FileType::Mgf]).unwrap();
    assert_eq!(copy.spectra[0].peaks, experiment().spectra[0].peaks);
    assert!(FileHandler::load_experiment(&path, &[FileType::MzMl]).is_err());
    assert!(
        FileHandler::store_experiment(dir.0.join("wrong.mzML"), &experiment(), Some(FileType::Mgf))
            .is_err()
    );
}

#[test]
fn invalid_serialization_preserves_existing_output_and_removes_temporary() {
    let dir = Directory::new();
    let path = dir.0.join("existing.mgf");
    std::fs::write(&path, "original").unwrap();
    let mut bad = experiment();
    bad.spectra[0].peaks[0].intensity = f32::NAN;
    assert!(FileHandler::store_experiment(&path, &bad, None).is_err());
    assert_eq!(std::fs::read_to_string(&path).unwrap(), "original");
    assert_eq!(std::fs::read_dir(&dir.0).unwrap().count(), 1);
    let bad_path = dir.0.join("output.mgf.zip");
    assert!(FileHandler::store_experiment(&bad_path, &experiment(), None).is_err());
    assert!(!bad_path.exists());
}

#[test]
fn dta_dispatch_rejects_extra_spectra_and_metadata_without_emitting_output() {
    let mut value = experiment();
    value.spectra.push(value.spectra[0].clone());
    let mut bytes = Vec::new();
    assert!(FileHandler::write_experiment(&mut bytes, &value, FileType::Dta).is_err());
    assert!(bytes.is_empty());
}

#[cfg(feature = "mzml")]
#[test]
fn gzip_roundtrip_and_content_detection_replay_decompressed_bytes() {
    let dir = Directory::new();
    let path = dir.0.join("peaks.mgf.gz");
    FileHandler::store_experiment(&path, &experiment(), None).unwrap();
    assert!(std::fs::read(&path).unwrap().starts_with(&[0x1f, 0x8b]));
    let unknown = dir.0.join("peaks.opaque");
    std::fs::rename(path, &unknown).unwrap();
    let copy = FileHandler::load_experiment(&unknown, &[]).unwrap();
    assert_eq!(copy.spectra[0].peaks, experiment().spectra[0].peaks);
}

#[cfg(feature = "mzml")]
#[test]
fn baseline_tool_core_chain_preserves_peaks_arrays_and_mzml_identity() {
    use openms::kernel::{DataArray, SpectrumType};
    use openms::processing::{
        SpectrumFilter,
        baseline::{MorphologicalFilter, MorphologicalMethod, StructuringElement},
    };
    let dir = Directory::new();
    let input = dir.0.join("profile.mzML");
    let output = dir.0.join("corrected.mzML.gz");
    let mut profile = experiment();
    profile.spectra[0].ms_level = 1;
    profile.spectra[0].rt = 15.;
    profile.spectra[0].native_id = "scan=1".into();
    profile.spectra[0].spectrum_type = SpectrumType::Profile;
    profile.spectra[0].peaks = [2., 2., 6., 2., 2.]
        .into_iter()
        .enumerate()
        .map(|(i, y)| Peak1D::new(100. + i as f64, y))
        .collect();
    profile.spectra[0]
        .integer_data_arrays
        .push(DataArray::new("index", vec![0, 1, 2, 3, 4]));
    FileHandler::store_experiment(&input, &profile, None).unwrap();
    let mut processed = FileHandler::load_experiment(&input, &[FileType::MzMl]).unwrap();
    MorphologicalFilter::new(
        MorphologicalMethod::TopHat,
        StructuringElement::DataPoints(3),
    )
    .unwrap()
    .filter_experiment(&mut processed)
    .unwrap();
    FileHandler::store_experiment(&output, &processed, None).unwrap();
    let result = FileHandler::load_experiment(&output, &[FileType::MzMl]).unwrap();
    let s = &result.spectra[0];
    assert_eq!(
        s.peaks.iter().map(|p| p.intensity).collect::<Vec<_>>(),
        [0., 0., 4., 0., 0.]
    );
    assert_eq!(s.native_id, "scan=1");
    assert_eq!(s.rt, 15.);
    assert_eq!(
        s.integer_data_arrays,
        profile.spectra[0].integer_data_arrays
    );
    assert_eq!(s.spectrum_type, SpectrumType::Profile);
}

#[cfg(feature = "file-compression")]
#[test]
fn bzip2_roundtrip_sniffing_and_truncated_stream_detection() {
    let dir = Directory::new();
    let path = dir.0.join("peaks.mgf.bz2");
    FileHandler::store_experiment(&path, &experiment(), None).unwrap();
    let compressed = std::fs::read(&path).unwrap();
    assert!(compressed.starts_with(b"BZh"));
    let unknown = dir.0.join("compressed.opaque");
    std::fs::rename(&path, &unknown).unwrap();
    assert_eq!(
        FileHandler::load_experiment(&unknown, &[]).unwrap().spectra[0].peaks,
        experiment().spectra[0].peaks
    );
    for cut in [3, compressed.len() - 1] {
        std::fs::write(&path, &compressed[..cut]).unwrap();
        assert!(FileHandler::load_experiment(&path, &[]).is_err());
    }
}

#[cfg(not(feature = "file-compression"))]
#[test]
fn unavailable_compression_is_reported_before_touching_destination() {
    let dir = Directory::new();
    for suffix in ["gz", "bz2"] {
        let path = dir.0.join(format!("peaks.mgf.{suffix}"));
        std::fs::write(&path, "original").unwrap();
        assert!(FileHandler::store_experiment(&path, &experiment(), None).is_err());
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "original");
    }
}

#[cfg(all(feature = "featurexml", feature = "consensusxml"))]
#[test]
fn feature_and_consensus_dispatch_preserve_typed_maps_through_unknown_suffixes() {
    let dir = Directory::new();
    let source = openms::format::featurexml::read(
        include_bytes!("data/featurexml_source_1.featureXML").as_slice(),
    )
    .unwrap();
    let path = dir.0.join("features.opaque.gz");
    FileHandler::store_feature_map(&path, &source, Some(FileType::FeatureXml)).unwrap();
    let mut copy = FileHandler::load_feature_map(&path, &[FileType::FeatureXml]).unwrap();
    assert_eq!(copy.loaded_file_type, FileType::FeatureXml);
    copy.loaded_file_path.clear();
    copy.loaded_file_type = FileType::Unknown;
    assert_eq!(copy, source);
    assert!(FileHandler::load_feature_map(&path, &[FileType::ConsensusXml]).is_err());
    let source = openms::format::consensusxml::read(
        include_bytes!("data/consensusxml/ConsensusXMLFile_1.consensusXML").as_slice(),
    )
    .unwrap();
    let path = dir.0.join("consensus.opaque.bz2");
    FileHandler::store_consensus_map(&path, &source, Some(FileType::ConsensusXml)).unwrap();
    let mut copy = FileHandler::load_consensus_map(&path, &[FileType::ConsensusXml]).unwrap();
    assert_eq!(copy.loaded_file_type, FileType::ConsensusXml);
    copy.loaded_file_path.clear();
    copy.loaded_file_type = FileType::Unknown;
    assert_eq!(copy, source);
    assert!(FileHandler::load_feature_map(&path, &[]).is_err());
}

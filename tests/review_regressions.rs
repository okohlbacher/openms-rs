// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::format::mgf;
use openms::processing::{LinearResamplerAlign, SpectrumFilter, ThresholdMower};
use openms::{MSExperiment, MSSpectrum, Peak1D, Precursor};

#[test]
fn threshold_rejects_nonfinite_peaks_without_mutation() {
    for intensity in [f32::INFINITY, f32::NEG_INFINITY, f32::NAN] {
        let mut spectrum =
            MSSpectrum::from_peaks(vec![Peak1D::new(100.0, intensity), Peak1D::new(200.0, 1.0)]);
        let before = spectrum.peaks[0].intensity.to_bits();
        assert!(
            ThresholdMower::default()
                .filter_spectrum(&mut spectrum)
                .is_err()
        );
        assert_eq!(spectrum.len(), 2);
        assert_eq!(spectrum.peaks[0].intensity.to_bits(), before);
    }
    let mut spectrum = MSSpectrum::from_peaks(vec![Peak1D::new(f64::INFINITY, 1.0)]);
    assert!(
        ThresholdMower::default()
            .filter_spectrum(&mut spectrum)
            .is_err()
    );
    assert_eq!(spectrum.peaks[0].mz, f64::INFINITY);
}

#[test]
fn aligned_resampling_does_not_invent_peaks_for_an_empty_spectrum() {
    // LinearResamplerAlign.h::raster_align explicitly returns for empty input.
    let mut spectrum = MSSpectrum {
        name: "empty scan".into(),
        ..MSSpectrum::default()
    };
    let before = spectrum.clone();
    LinearResamplerAlign::new(0.5)
        .unwrap()
        .raster_align(&mut spectrum, 100.0, 101.0)
        .unwrap();
    assert_eq!(spectrum, before);
}

#[test]
fn mgf_charge_roundtrip_covers_entire_i32_domain() {
    for charge in [i32::MIN, -3, 0, 2, i32::MAX] {
        let experiment = MSExperiment {
            spectra: vec![MSSpectrum {
                ms_level: 2,
                precursors: vec![Precursor::new(500.0, charge)],
                ..MSSpectrum::default()
            }],
            ..MSExperiment::default()
        };
        let mut bytes = Vec::new();
        mgf::write(&mut bytes, &experiment).unwrap();
        let copy = mgf::read(bytes.as_slice()).unwrap();
        assert_eq!(copy.spectra[0].precursors, experiment.spectra[0].precursors);
    }
    for charge in ["2147483648+", "2147483649-", "-2-", "2+ and 3+"] {
        let text = format!("BEGIN IONS\nCHARGE={charge}\nEND IONS\n");
        assert!(mgf::read(text.as_bytes()).is_err());
    }
}

#[test]
fn mgf_writer_rejects_metadata_that_would_be_dropped_as_a_comment() {
    for key in ["#COMMENT", ";COMMENT", "!COMMENT", "/COMMENT"] {
        let mut spectrum = MSSpectrum::default();
        spectrum
            .metadata
            .insert(key.into(), "retained value".into());
        let experiment = MSExperiment {
            spectra: vec![MSSpectrum::default(), spectrum],
            ..MSExperiment::default()
        };
        let mut bytes = Vec::new();
        assert!(
            mgf::write(&mut bytes, &experiment).is_err(),
            "accepted {key}"
        );
        assert!(bytes.is_empty(), "semantic validation must precede output");
    }
}

#[test]
fn mgf_writer_rejects_case_insensitive_duplicate_keys_before_output() {
    let mut spectrum = MSSpectrum::default();
    spectrum.metadata.insert("SCANS".into(), "1".into());
    spectrum.metadata.insert("scans".into(), "2".into());
    let experiment = MSExperiment {
        spectra: vec![spectrum],
        ..MSExperiment::default()
    };
    let mut bytes = Vec::new();
    assert!(mgf::write(&mut bytes, &experiment).is_err());
    assert!(bytes.is_empty());
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::kernel::{BaseFeature, ConsensusFeature, Feature, PeakIndex};
use openms::{MSExperiment, MSSpectrum, Peak1D};

#[test]
fn constructor_validity_and_clear_keep_source_component_rules() {
    let invalid = PeakIndex::default();
    assert_eq!((invalid.spectrum, invalid.peak), (usize::MAX, usize::MAX));
    assert!(!invalid.is_valid());
    let feature = PeakIndex::for_feature(17);
    assert_eq!((feature.spectrum, feature.peak), (usize::MAX, 17));
    assert!(feature.is_valid());
    assert!(!PeakIndex::for_feature(usize::MAX).is_valid());
    let mut peak = PeakIndex::new(5, 17);
    assert_eq!((peak.spectrum, peak.peak), (5, 17));
    assert!(peak.is_valid());
    let copied = peak;
    assert_eq!(peak, copied);
    assert_ne!(peak, feature);
    assert_ne!(peak, PeakIndex::new(2, 5));
    assert!(!PeakIndex::new(2, usize::MAX).is_valid());
    peak.clear();
    assert_eq!(peak, invalid);
}

#[test]
fn source_feature_and_peak_access_literals() {
    let features: Vec<_> = (1..=5)
        .map(|mz| Feature::from(BaseFeature::new(0., f64::from(mz), 0.)))
        .collect();
    let consensus: Vec<_> = (1..=5)
        .map(|mz| ConsensusFeature::from(BaseFeature::new(0., f64::from(mz) + 0.1, 0.)))
        .collect();
    for (index, mz) in [(4, 5.), (0, 1.)] {
        assert_eq!(
            PeakIndex::for_feature(index)
                .get_feature(&features)
                .unwrap()
                .mz,
            mz
        );
        assert_eq!(
            PeakIndex::for_feature(index)
                .get_feature(&consensus)
                .unwrap()
                .mz,
            mz + 0.1
        );
    }
    assert!(PeakIndex::for_feature(5).get_feature(&features).is_err());
    let mut exp = MSExperiment {
        spectra: vec![MSSpectrum::new(), MSSpectrum::new(), MSSpectrum::new()],
        ..Default::default()
    };
    for (i, s) in exp.spectra.iter_mut().enumerate() {
        s.rt = (i + 1) as f64;
    }
    exp.spectra[0].peaks.resize(15, Peak1D::default());
    exp.spectra[2].peaks = (1..=3).map(|mz| Peak1D::new(f64::from(mz), 0.)).collect();
    assert_eq!(
        PeakIndex::new(0, usize::MAX).get_spectrum(&exp).unwrap().rt,
        1.
    );
    assert_eq!(
        PeakIndex::new(2, usize::MAX).get_spectrum(&exp).unwrap().rt,
        3.
    );
    assert_eq!(PeakIndex::new(0, 0).get_peak(&exp).unwrap().mz, 0.);
    assert_eq!(PeakIndex::new(2, 0).get_peak(&exp).unwrap().mz, 1.);
    assert_eq!(PeakIndex::new(2, 2).get_peak(&exp).unwrap().mz, 3.);
    assert!(PeakIndex::new(2, 16).get_peak(&exp).is_err());
    assert!(PeakIndex::new(3, 0).get_peak(&exp).is_err());
}

#[test]
fn borrowed_access_checks_only_the_relevant_dimensions() {
    let exp = MSExperiment {
        spectra: vec![
            MSSpectrum::new(),
            MSSpectrum::from_peaks(vec![Peak1D::new(10., 20.), Peak1D::new(30., 40.)]),
        ],
        ..Default::default()
    };
    let index = PeakIndex::new(1, 1);
    assert!(std::ptr::eq(
        index.get_peak(&exp).unwrap(),
        &exp.spectra[1].peaks[1]
    ));
    assert!(std::ptr::eq(
        index.get_spectrum(&exp).unwrap(),
        &exp.spectra[1]
    ));
    // Source getSpectrum does not inspect isValid or the peak dimension.
    assert!(PeakIndex::new(1, usize::MAX).get_spectrum(&exp).is_ok());
    assert!(PeakIndex::new(0, 0).get_peak(&exp).is_err());
    assert!(PeakIndex::new(1, 2).get_peak(&exp).is_err());
    assert!(PeakIndex::new(2, 0).get_spectrum(&exp).is_err());
    assert!(PeakIndex::for_feature(0).get_peak(&exp).is_err());
    assert!(PeakIndex::default().get_feature::<u8>(&[]).is_err());
    let features = vec![Feature::default(), Feature::default()];
    assert!(std::ptr::eq(
        index.get_feature(&features).unwrap(),
        &features[1]
    ));
    let consensus = vec![ConsensusFeature::default()];
    assert!(std::ptr::eq(
        PeakIndex::new(999, 0).get_feature(&consensus).unwrap(),
        &consensus[0]
    ));
}

#[test]
fn slice_oracle_covers_every_dimension_and_boundary_without_index_arithmetic() {
    let exp = MSExperiment {
        spectra: (0..4)
            .map(|s| {
                MSSpectrum::from_peaks((0..s).map(|p| Peak1D::new(p as f64, s as f32)).collect())
            })
            .collect(),
        ..Default::default()
    };
    for s in [0, 1, 2, 3, 4, usize::MAX] {
        for p in [0, 1, 2, 3, 4, usize::MAX] {
            let index = PeakIndex::new(s, p);
            assert_eq!(index.is_valid(), p != usize::MAX);
            assert_eq!(index.get_spectrum(&exp).ok(), exp.spectra.get(s));
            assert_eq!(
                index.get_peak(&exp).ok(),
                exp.spectra.get(s).and_then(|s| s.peaks.get(p))
            );
        }
    }
}

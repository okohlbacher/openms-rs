// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{AASequence, IonSeries, ProteaseDigestion};
use openms::format::{dta, fasta, mgf};
use openms::processing::{NLargest, Normalizer, SpectrumFilter, ThresholdMower};
use openms::{MSExperiment, MSSpectrum, Peak1D, Precursor};

#[test]
fn fasta_to_fragments_to_mgf_to_processing() {
    let proteins = fasta::read(b">P1 test protein\nDFPIANGER\n".as_slice()).unwrap();
    let protein = AASequence::parse(&proteins[0].sequence).unwrap();
    let peptides = ProteaseDigestion::default().digest(&protein).unwrap();
    assert_eq!(peptides.len(), 1);
    let peptide = &peptides[0].sequence;
    assert!((peptide.mono_mass().unwrap() - 1017.487964).abs() < 0.00001);
    let fragments = peptide.fragment_ions(1).unwrap();
    assert_eq!(
        fragments
            .iter()
            .filter(|f| f.series == IonSeries::B)
            .count(),
        8
    );
    assert_eq!(
        fragments
            .iter()
            .filter(|f| f.series == IonSeries::Y)
            .count(),
        8
    );
    let mut spec = MSSpectrum {
        ms_level: 2,
        precursors: vec![Precursor::new(peptide.mz(2).unwrap(), 2)],
        peaks: fragments.iter().map(|f| Peak1D::new(f.mz, 100.0)).collect(),
        ..Default::default()
    };
    spec.sort_by_position().unwrap();
    let exp = MSExperiment {
        spectra: vec![spec],
        ..Default::default()
    };
    let mut output = Vec::new();
    mgf::write(&mut output, &exp).unwrap();
    let mut parsed = mgf::read(output.as_slice()).unwrap();
    assert_eq!(parsed.spectra[0].peaks, exp.spectra[0].peaks);
    Normalizer::default()
        .filter_experiment(&mut parsed)
        .unwrap();
    assert_eq!(parsed.spectra[0].calculate_tic(), 16.0);
}

#[test]
fn upstream_dta_filter_and_mgf_roundtrip() {
    let mut spec = dta::read(&include_bytes!("data/Transformers_tests.dta")[..]).unwrap();
    ThresholdMower { threshold: 10.0 }
        .filter_spectrum(&mut spec)
        .unwrap();
    NLargest { n: 10 }.filter_spectrum(&mut spec).unwrap();
    Normalizer::default().filter_spectrum(&mut spec).unwrap();
    spec.sort_by_position().unwrap();
    let experiment = MSExperiment {
        spectra: vec![spec],
        ..Default::default()
    };
    let mut bytes = Vec::new();
    mgf::write(&mut bytes, &experiment).unwrap();
    let parsed = mgf::read(bytes.as_slice()).unwrap();
    assert_eq!(parsed.spectra[0].peaks, experiment.spectra[0].peaks);
    assert_eq!(
        parsed.spectra[0].precursors,
        experiment.spectra[0].precursors
    );
    assert_eq!(parsed.spectra[0].len(), 10);
}

#[test]
fn modified_peptide_envelopes_collapse_to_annotated_fragment_spectrum() {
    use openms::chemistry::{TheoreticalIsotopeModel, TheoreticalSpectrumGenerator};
    use openms::comparison::{SpectrumAlignment, Tolerance};
    use openms::processing::deisotoping::Deisotoper;
    let peptide = AASequence::parse("(Acetyl)AC(Carbamidomethyl)M(Oxidation)K").unwrap();
    let generator = TheoreticalSpectrumGenerator {
        add_metainfo: true,
        isotope_model: TheoreticalIsotopeModel::Coarse { max_peaks: 4 },
        ..Default::default()
    };
    let envelopes = generator.generate(&peptide, 1, 1, Some(3)).unwrap();
    assert_eq!(envelopes.len(), 20);
    let result = Deisotoper {
        max_charge: 1,
        keep_only_deisotoped: true,
        use_decreasing_model: false,
        add_up_intensity: true,
        annotate_charge: true,
        annotate_isotope_peak_count: true,
        ..Default::default()
    }
    .deisotope(&envelopes)
    .unwrap();
    assert_eq!(result.clusters.len(), 5);
    assert_eq!(result.spectrum.len(), 5);
    assert!(
        result
            .clusters
            .iter()
            .all(|c| c.charge == 1 && c.peak_indices.len() == 4)
    );
    result.spectrum.validate().unwrap();
    assert!(
        result
            .spectrum
            .peaks
            .iter()
            .all(|p| (p.intensity - 1.0).abs() < 2e-7)
    );
    let monoisotopic = TheoreticalSpectrumGenerator {
        isotope_model: TheoreticalIsotopeModel::None,
        ..generator
    }
    .generate(&peptide, 1, 1, Some(3))
    .unwrap();
    // The source's coarse generator adds neutral H rather than proton mass,
    // plus its formula anchor differs slightly from declared terminal deltas.
    let alignment = SpectrumAlignment {
        tolerance: Tolerance::Absolute(0.002),
        ..Default::default()
    };
    assert_eq!(
        alignment.align(&result.spectrum, &monoisotopic).unwrap(),
        [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)]
    );
    for (coarse, mono) in result.spectrum.peaks.iter().zip(&monoisotopic.peaks) {
        assert!((0.00054..0.00056).contains(&(coarse.mz - mono.mz)));
    }
}

#[test]
fn profile_centroids_become_consensus_measurements() {
    use openms::kernel::{BaseFeature, ConsensusFeature, FeatureHandle};
    use openms::processing::peak_picking::PeakPickerHiRes;
    let profile = MSSpectrum::from_peaks(
        [200., 250., 450., 250., 200.]
            .into_iter()
            .enumerate()
            .map(|(i, intensity)| Peak1D::new(100.0 + i as f64 * 0.01, intensity))
            .collect(),
    );
    let picked = PeakPickerHiRes::default().pick_spectrum(&profile).unwrap();
    assert_eq!(picked.spectrum.len(), 1);
    let peak = picked.spectrum.peaks[0];
    let a = BaseFeature {
        rt: 100.0,
        mz: peak.mz,
        intensity: peak.intensity,
        unique_id: 1,
        ..Default::default()
    };
    let b = BaseFeature {
        rt: 110.0,
        intensity: 2.0 * peak.intensity,
        unique_id: 2,
        ..a.clone()
    };
    let mut consensus = ConsensusFeature::from_feature(0, &a).unwrap();
    consensus.insert(FeatureHandle::new(1, &b)).unwrap();
    consensus.compute_consensus().unwrap();
    assert_eq!(consensus.rt, 105.0);
    assert!((consensus.mz - 100.02).abs() < 1e-6);
    assert!((consensus.intensity - 675.0).abs() < 1e-3);
    assert_eq!(consensus.len(), 2);
}

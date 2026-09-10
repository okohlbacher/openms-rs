// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{
    AASequence, EmpiricalFormula, FineIsotopePatternGenerator, FineIsotopeStop,
    TheoreticalIonSeries, TheoreticalIsotopeModel, TheoreticalSpectrumGenerator,
};
use openms::{MSSpectrum, Peak1D};

fn settings() -> TheoreticalSpectrumGenerator {
    TheoreticalSpectrumGenerator {
        isotope_model: TheoreticalIsotopeModel::Fine {
            unexplained_probability: 0.01,
        },
        add_metainfo: true,
        add_losses: true,
        add_precursor_peaks: true,
        add_all_precursor_charges: true,
        ..Default::default()
    }
}

#[test]
fn modified_peptide_fine_formulas_and_annotation_alignment_survive_selection() {
    let peptide = AASequence::parse("AC(Carbamidomethyl)M(Oxidation)K").unwrap();
    let generator = FineIsotopePatternGenerator::default();
    let charged = peptide.formula().unwrap().with_charge(2);
    let explicit = peptide
        .formula()
        .unwrap()
        .checked_add(&EmpiricalFormula::parse("H2").unwrap())
        .unwrap();
    assert_eq!(
        generator.run(&charged).unwrap(),
        generator.run(&explicit).unwrap()
    );
    assert_eq!(charged.charge(), 2);

    let spectrum = settings().generate(&peptide, 1, 2, Some(3)).unwrap();
    assert!(spectrum.is_sorted());
    assert!(spectrum.peaks.len() > 30);
    assert_eq!(
        spectrum.string_data_arrays[0].data.len(),
        spectrum.peaks.len()
    );
    assert_eq!(
        spectrum.integer_data_arrays[0].data.len(),
        spectrum.peaks.len()
    );
    let chosen: Vec<_> = spectrum.integer_data_arrays[0]
        .data
        .iter()
        .enumerate()
        .filter_map(|(index, &charge)| (charge == 2).then_some(index))
        .collect();
    let mut selected = spectrum.clone();
    selected.select(&chosen).unwrap();
    selected.validate().unwrap();
    assert!(!selected.is_empty());
    for (index, &original) in chosen.iter().enumerate() {
        assert_eq!(selected.peaks[index], spectrum.peaks[original]);
        assert_eq!(
            selected.string_data_arrays[0].data[index],
            spectrum.string_data_arrays[0].data[original]
        );
        assert_eq!(selected.integer_data_arrays[0].data[index], 2);
    }
}

#[test]
fn isotope_requests_require_known_chemistry_and_append_errors_preserve_the_destination() {
    let peptide = AASequence::parse("AX[150.0]K").unwrap();
    let mut target = MSSpectrum::from_peaks(vec![Peak1D::new(50.0, 7.0)]);
    target.name = "preserved observed data".into();
    let before = target.clone();
    assert!(
        settings()
            .append_to(&mut target, &peptide, 1, 2, None)
            .is_err()
    );
    assert_eq!(target, before);

    // The same observed mass remains usable without a formula-based request.
    let mono = TheoreticalSpectrumGenerator {
        ion_series: vec![TheoreticalIonSeries::B],
        add_first_prefix_ion: true,
        ..Default::default()
    };
    assert!(mono.generate(&peptide, 1, 1, None).is_ok());
    let invalid = FineIsotopePatternGenerator {
        stop: FineIsotopeStop::UnexplainedProbability(f64::NAN),
    };
    assert!(invalid.run(&EmpiricalFormula::default()).is_err());
}

#[cfg(feature = "mzml")]
#[test]
fn fine_theoretical_spectra_roundtrip_through_mzml_with_annotations_and_acquisition() {
    use openms::MSExperiment;
    use openms::format::mzml;
    use openms::metadata::{ActivationMethod, DriftTimeUnit};
    let peptide = AASequence::parse("AC(Carbamidomethyl)M(Oxidation)K").unwrap();
    let mut spectrum = settings().generate(&peptide, 1, 2, Some(3)).unwrap();
    spectrum.native_id = "scan=1".into();
    spectrum.name = "modified peptide fine isotope spectrum".into();
    let p = &mut spectrum.precursors[0];
    p.activation_methods.insert(ActivationMethod::Hcd);
    p.isolation_window_lower_offset = 0.4;
    p.isolation_window_upper_offset = 0.5;
    p.drift_time = Some(0.9);
    p.drift_time_unit = DriftTimeUnit::InverseReducedMobility;
    let e = MSExperiment {
        spectra: vec![spectrum],
        ..Default::default()
    };
    for compressed in [false, true] {
        let mut output = Vec::new();
        mzml::write_with_options(
            &mut output,
            &e,
            &mzml::WriteOptions {
                zlib_compression: compressed,
            },
        )
        .unwrap();
        assert_eq!(mzml::read(output.as_slice()).unwrap(), e);
    }
}

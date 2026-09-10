// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::MSSpectrum;
use openms::chemistry::{FineIsotopeIterator, FineIsotopePatternGenerator, FineIsotopeStop};
use openms::kernel::DataArray;

fn populations() -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    // Two independent single-atom populations can have the same total mass
    // while retaining distinct configuration probabilities.
    (
        vec![vec![12.0, 13.0]; 2],
        vec![vec![0.25, 0.75], vec![0.5, 0.5]],
    )
}

fn annotated_spectrum() -> MSSpectrum {
    let (masses, probabilities) = populations();
    let configurations: Vec<_> =
        FineIsotopeIterator::from_isotopes(&[1, 1], &masses, &probabilities)
            .unwrap()
            .collect::<openms::Result<_>>()
            .unwrap();
    let mut spectrum = MSSpectrum::from_peaks(
        configurations
            .iter()
            .map(|c| c.to_peak().unwrap())
            .collect(),
    );
    spectrum.integer_data_arrays.push(DataArray::new(
        "ConfigurationIndex",
        (0..configurations.len()).map(|i| i as i32).collect(),
    ));
    spectrum.sort_by_position().unwrap();
    for (peak, &index) in spectrum
        .peaks
        .iter()
        .zip(&spectrum.integer_data_arrays[0].data)
    {
        assert_eq!(*peak, configurations[index as usize].to_peak().unwrap());
    }
    spectrum
}

#[test]
fn enriched_streams_materialize_without_merging_and_keep_annotations_through_selection() {
    let (masses, probabilities) = populations();
    let distribution = FineIsotopePatternGenerator {
        stop: FineIsotopeStop::AbsoluteThreshold(0.0),
    }
    .run_with_isotopes(&[1, 1], &masses, &probabilities)
    .unwrap();
    let mut spectrum = annotated_spectrum();
    assert_eq!(distribution.len(), 4);
    for (isotope, peak) in distribution.peaks().iter().zip(&spectrum.peaks) {
        assert_eq!(isotope.mass, peak.mz);
        assert_eq!(isotope.probability, f64::from(peak.intensity));
    }
    let mut middle_weights: Vec<_> = spectrum
        .peaks
        .iter()
        .filter(|p| p.mz == 25.0)
        .map(|p| p.intensity)
        .collect();
    middle_weights.sort_by(f32::total_cmp);
    assert_eq!(middle_weights, [0.125, 0.375]);

    let chosen: Vec<_> = spectrum
        .peaks
        .iter()
        .enumerate()
        .filter_map(|(i, peak)| (peak.intensity > 0.2).then_some(i))
        .collect();
    let before = spectrum.clone();
    spectrum.select(&chosen).unwrap();
    spectrum.validate().unwrap();
    assert_eq!(spectrum.peaks.len(), 2);
    for (current, &original) in chosen.iter().enumerate() {
        assert_eq!(spectrum.peaks[current], before.peaks[original]);
        assert_eq!(
            spectrum.integer_data_arrays[0].data[current],
            before.integer_data_arrays[0].data[original]
        );
    }
}

#[test]
fn a_caller_can_stop_at_raw_coverage_and_resume_the_same_owned_stream() {
    let (masses, probabilities) = populations();
    let mut stream = FineIsotopeIterator::from_isotopes(&[1, 1], &masses, &probabilities).unwrap();
    drop(masses);
    drop(probabilities);
    let mut coverage = 0.0;
    let mut prefix = Vec::new();
    while coverage < 0.7 {
        let config = stream.next().unwrap().unwrap();
        coverage += config.probability;
        prefix.push(config);
    }
    assert_eq!(prefix.len(), 2);
    assert!((coverage - 0.75).abs() < 1e-15);
    let remaining = stream.collect::<openms::Result<Vec<_>>>().unwrap();
    assert_eq!(remaining.len(), 2);
    assert!((remaining.iter().map(|c| c.probability).sum::<f64>() - 0.25).abs() < 1e-15);
    // Yielded values own their data and remain valid after exhausting the stream.
    assert_eq!(prefix[0].to_peak().unwrap().intensity, 0.375);
}

#[cfg(feature = "mzml")]
#[test]
fn duplicate_mass_enrichment_configurations_survive_mzml_as_distinct_annotated_peaks() {
    use openms::MSExperiment;
    use openms::format::mzml;
    let mut spectrum = annotated_spectrum();
    spectrum.native_id = "scan=1".into();
    spectrum.name = "independent enriched populations".into();
    let experiment = MSExperiment {
        spectra: vec![spectrum],
        ..Default::default()
    };
    for zlib_compression in [false, true] {
        let mut output = Vec::new();
        mzml::write_with_options(
            &mut output,
            &experiment,
            &mzml::WriteOptions { zlib_compression },
        )
        .unwrap();
        let restored = mzml::read(output.as_slice()).unwrap();
        assert_eq!(restored, experiment);
        assert_eq!(
            restored.spectra[0]
                .peaks
                .iter()
                .filter(|p| p.mz == 25.0)
                .count(),
            2
        );
    }
}

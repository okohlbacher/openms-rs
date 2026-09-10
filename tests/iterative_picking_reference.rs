// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Independent source-derived iterative-picker checks; source tests have no numerical goldens.

use openms::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, SpectrumType,
};
use openms::processing::SpectrumFilter;
use openms::processing::iterative::PeakPickerIterative;

fn symmetric(offset: f64) -> MSSpectrum {
    MSSpectrum::from_peaks(
        [0., 1., 4., 10., 4., 1., 0.]
            .into_iter()
            .enumerate()
            .map(|(i, y)| Peak1D::new(offset + i as f64, y))
            .collect(),
    )
}
fn picker() -> PeakPickerIterative {
    PeakPickerIterative {
        signal_to_noise: 0.,
        ..Default::default()
    }
}

#[test]
fn symmetric_public_picker_reports_sample_sum_and_full_precision_borders() {
    for offset in [100., 100_000_000.25] {
        let mut input = symmetric(offset);
        input.name = "manual symmetric trace".into();
        input.native_id = "scan=7".into();
        input.metadata.insert("sample".into(), "reference".into());
        input
            .float_data_arrays
            .push(DataArray::new("profile labels", vec![7.; 7]));
        input
            .integer_data_arrays
            .push(DataArray::new("indices", vec![1; 7]));
        input
            .string_data_arrays
            .push(DataArray::new("notes", vec!["point".to_owned(); 7]));
        let original = input.clone();
        let result = picker().pick_spectrum(&input).unwrap();
        assert_eq!(input, original);
        assert_eq!(
            result.picked.spectrum.peaks,
            [Peak1D::new(f64::from((offset + 3.) as f32), 20.)]
        );
        assert_eq!(result.regions.len(), 1);
        assert_eq!(
            (
                result.regions[0].center_index,
                result.regions[0].left_index,
                result.regions[0].right_index
            ),
            (3, 0, 6)
        );
        assert_eq!(result.picked.boundaries[0].min, offset);
        assert_eq!(result.picked.boundaries[0].max, offset + 6.);
        let out = &result.picked.spectrum;
        assert_eq!(
            out.float_data_arrays[0],
            DataArray::new("IntegratedIntensity", vec![20.])
        );
        assert_eq!(
            out.float_data_arrays[1],
            DataArray::new("leftWidth", vec![offset as f32])
        );
        assert_eq!(
            out.float_data_arrays[2],
            DataArray::new("rightWidth", vec![(offset + 6.) as f32])
        );
        assert_eq!(out.spectrum_type, SpectrumType::Centroid);
        assert_eq!(out.metadata, input.metadata);
        assert_eq!(out.name, input.name);
        assert_eq!(out.native_id, input.native_id);
        assert_eq!(
            result.picked.omitted_arrays,
            ["profile labels", "indices", "notes"]
        );
        assert!(out.integer_data_arrays.is_empty() && out.string_data_arrays.is_empty());
        out.validate().unwrap();
    }
}

#[test]
fn experiment_ms1_selection_clears_only_generated_arrays_and_keeps_chromatograms() {
    let mut ms1 = symmetric(100.);
    ms1.float_data_arrays
        .push(DataArray::new("input annotations", vec![1.; 7]));
    let mut ms2 = ms1.clone();
    ms2.ms_level = 2;
    let chrom = MSChromatogram {
        peaks: vec![ChromatogramPeak::new(1., 2.)],
        ..Default::default()
    };
    let input = MSExperiment {
        spectra: vec![ms1, ms2.clone()],
        chromatograms: vec![chrom.clone()],
        ..Default::default()
    };
    let original = input.clone();
    let p = PeakPickerIterative {
        ms1_only: true,
        clear_meta_data: true,
        ..picker()
    };
    let result = p.pick_experiment(&input).unwrap();
    assert_eq!(input, original);
    assert_eq!(result.experiment.spectra[0].peaks, [Peak1D::new(103., 20.)]);
    assert!(result.experiment.spectra[0].float_data_arrays.is_empty());
    assert_eq!(result.omitted_spectrum_arrays[0], ["input annotations"]);
    assert_eq!(result.experiment.spectra[1], ms2);
    assert_eq!(result.spectrum_regions[0].as_ref().unwrap().len(), 1);
    assert_eq!(result.spectrum_regions[1], None);
    assert_eq!(result.experiment.chromatograms, [chrom]);
    // clear_meta_data is a source experiment option; direct picking keeps arrays.
    assert_eq!(
        p.pick_spectrum(&input.spectra[0])
            .unwrap()
            .picked
            .spectrum
            .float_data_arrays
            .len(),
        3
    );
}

#[test]
fn invalid_late_records_and_resource_failures_are_atomic() {
    let input = symmetric(100.);
    let mut wrong = input.clone();
    wrong.peaks.swap(2, 3);
    let experiment = MSExperiment {
        spectra: vec![input.clone(), wrong],
        ..Default::default()
    };
    let original = experiment.clone();
    assert!(picker().pick_experiment(&experiment).is_err());
    assert_eq!(experiment, original);
    let mut mutable = experiment.clone();
    assert!(picker().filter_experiment(&mut mutable).is_err());
    assert_eq!(mutable, original);
    for p in [
        PeakPickerIterative {
            max_points: 6,
            ..picker()
        },
        PeakPickerIterative {
            max_work: 1,
            ..picker()
        },
        PeakPickerIterative {
            iterations: usize::MAX,
            ..picker()
        },
        PeakPickerIterative {
            peak_width: f64::NAN,
            ..picker()
        },
    ] {
        let mut value = input.clone();
        assert!(p.filter_spectrum(&mut value).is_err());
        assert_eq!(value, input);
    }
    let default_result = PeakPickerIterative::default()
        .pick_spectrum(&input)
        .unwrap();
    // Seven samples leave the source's ten-element noise windows sparse.
    assert!(default_result.picked.spectrum.is_empty());
}

#[test]
fn realistic_profile_results_integrate_reported_raw_regions() {
    // These are original HiRes INPUT traces. Their existing HiRes output files
    // are deliberately not used as iterative centroid reference values.
    for text in [
        include_str!("data/peak_picking_orbitrap.tsv"),
        include_str!("data/peak_picking_ftms.tsv"),
    ] {
        let input = MSSpectrum::from_peaks(
            text.lines()
                .filter(|l| !l.starts_with('#'))
                .map(|line| {
                    let fields: Vec<_> = line.split_whitespace().collect();
                    Peak1D::new(fields[0].parse().unwrap(), fields[1].parse().unwrap())
                })
                .collect(),
        );
        let result = PeakPickerIterative::default()
            .pick_spectrum(&input)
            .unwrap();
        assert!(!result.picked.spectrum.is_empty());
        assert_eq!(result.regions.len(), result.picked.spectrum.len());
        assert!(result.picked.spectrum.is_sorted());
        for (region, peak) in result.regions.iter().zip(&result.picked.spectrum.peaks) {
            let raw = &input.peaks[region.left_index..=region.right_index];
            let sum = raw.iter().fold(0., |v, p| v + f64::from(p.intensity));
            let weighted = raw
                .iter()
                .fold(0., |v, p| v + p.mz * f64::from(p.intensity));
            assert_eq!(peak.intensity, sum as f32);
            assert_eq!(peak.mz, f64::from((weighted / sum) as f32));
            assert!(
                region.left_index < region.center_index && region.center_index < region.right_index
            );
        }
        result.picked.spectrum.validate().unwrap();
    }
}

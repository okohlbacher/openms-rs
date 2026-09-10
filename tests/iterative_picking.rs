// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, SpectrumType,
};
use openms::processing::SpectrumFilter;
use openms::processing::iterative::PeakPickerIterative;

fn profile() -> MSSpectrum {
    MSSpectrum {
        peaks: [0.0, 1.0, 4.0, 10.0, 4.0, 1.0, 0.0]
            .into_iter()
            .enumerate()
            .map(|(i, intensity)| Peak1D::new(100.0 + i as f64, intensity))
            .collect(),
        spectrum_type: SpectrumType::Profile,
        ..Default::default()
    }
}
fn picker() -> PeakPickerIterative {
    PeakPickerIterative {
        signal_to_noise: 0.0,
        ..Default::default()
    }
}

#[test]
fn analytic_symmetric_peak_has_raw_area_and_original_boundaries() {
    let input = profile();
    let output = picker().pick_spectrum(&input).unwrap();
    assert_eq!(output.picked.spectrum.peaks, [Peak1D::new(103.0, 20.0)]);
    assert_eq!(output.picked.spectrum.spectrum_type, SpectrumType::Centroid);
    assert_eq!(output.regions.len(), 1);
    let region = output.regions[0];
    assert_eq!(
        (region.center_index, region.left_index, region.right_index),
        (3, 0, 6)
    );
    assert_eq!(output.picked.boundaries[0].min, 100.0);
    assert_eq!(output.picked.boundaries[0].max, 106.0);
    let arrays = &output.picked.spectrum.float_data_arrays;
    assert_eq!(
        arrays.iter().map(|a| a.name.as_str()).collect::<Vec<_>>(),
        ["IntegratedIntensity", "leftWidth", "rightWidth"]
    );
    assert_eq!(arrays[0].data, [20.0]);
    assert_eq!(arrays[1].data, [100.0]);
    assert_eq!(arrays[2].data, [106.0]);
    assert_eq!(input, profile());
}

#[test]
fn default_seed_noise_rejects_sparse_peak_and_width_check_can_remove_candidate() {
    let input = profile();
    assert!(
        PeakPickerIterative::default()
            .pick_spectrum(&input)
            .unwrap()
            .picked
            .spectrum
            .is_empty()
    );
    let checked = PeakPickerIterative {
        check_width_internally: true,
        ..picker()
    };
    assert!(
        checked
            .pick_spectrum(&input)
            .unwrap()
            .picked
            .spectrum
            .is_empty()
    );
    // Equality is permitted by the source's internal width sanity check.
    let checked = PeakPickerIterative {
        peak_width: 1.0,
        ..checked
    };
    assert_eq!(
        checked.pick_spectrum(&input).unwrap().picked.spectrum.len(),
        1
    );
}

#[test]
fn profile_arrays_are_reported_and_all_record_metadata_is_preserved() {
    let mut input = profile();
    input.native_id = "scan=13".into();
    input.name = "TOF profile".into();
    input.rt = 123.5;
    input.metadata.insert("source".into(), "synthetic".into());
    input
        .float_data_arrays
        .push(DataArray::new("Ion Mobility", vec![1.0; 7]));
    input
        .integer_data_arrays
        .push(DataArray::new("index", (0..7).collect()));
    input
        .string_data_arrays
        .push(DataArray::new("labels", vec!["sample".into(); 7]));
    let result = picker().pick_spectrum(&input).unwrap();
    assert_eq!(
        result.picked.omitted_arrays,
        ["Ion Mobility", "index", "labels"]
    );
    assert_eq!(result.picked.spectrum.native_id, input.native_id);
    assert_eq!(result.picked.spectrum.name, input.name);
    assert_eq!(result.picked.spectrum.rt, input.rt);
    assert_eq!(result.picked.spectrum.metadata, input.metadata);
    assert!(result.picked.spectrum.integer_data_arrays.is_empty());
    assert!(result.picked.spectrum.string_data_arrays.is_empty());
    result.picked.spectrum.validate().unwrap();
}

#[test]
fn experiment_flags_copy_excluded_spectra_and_preserve_chromatograms() {
    let mut ms2 = profile();
    ms2.ms_level = 2;
    ms2.float_data_arrays
        .push(DataArray::new("keep", vec![3.0; 7]));
    let input = MSExperiment {
        spectra: vec![profile(), ms2.clone()],
        chromatograms: vec![MSChromatogram {
            peaks: vec![ChromatogramPeak::new(1.0, 12.0)],
            ..Default::default()
        }],
        ..Default::default()
    };
    let picker = PeakPickerIterative {
        ms1_only: true,
        clear_meta_data: true,
        ..picker()
    };
    // Experiment flags have no effect on the single-spectrum overload.
    assert_eq!(
        picker
            .pick_spectrum(&ms2)
            .unwrap()
            .picked
            .spectrum
            .float_data_arrays
            .len(),
        3
    );
    let output = picker.pick_experiment(&input).unwrap();
    assert_eq!(
        output.experiment.spectra[0].peaks,
        [Peak1D::new(103.0, 20.0)]
    );
    assert!(output.experiment.spectra[0].float_data_arrays.is_empty());
    assert_eq!(output.experiment.spectra[1], ms2);
    assert_eq!(output.experiment.chromatograms, input.chromatograms);
    assert!(output.spectrum_regions[0].is_some());
    assert!(output.spectrum_regions[1].is_none());
    assert!(output.omitted_spectrum_arrays[1].is_empty());
    let mut mutable = input.clone();
    picker.filter_experiment(&mut mutable).unwrap();
    assert_eq!(mutable, output.experiment);
}

#[test]
fn errors_are_atomic_for_spectrum_and_experiment() {
    let mut invalid = profile();
    invalid.peaks[2].intensity = -1.0;
    let before = invalid.clone();
    assert!(picker().filter_spectrum(&mut invalid).is_err());
    assert_eq!(invalid, before);
    let mut experiment = MSExperiment {
        spectra: vec![profile(), invalid],
        ..Default::default()
    };
    let before = experiment.clone();
    assert!(picker().filter_experiment(&mut experiment).is_err());
    assert_eq!(experiment, before);
    let mut overflowing = profile();
    for p in &mut overflowing.peaks {
        p.intensity *= f32::MAX / 10.0;
    }
    let before = overflowing.clone();
    assert!(picker().filter_spectrum(&mut overflowing).is_err());
    assert_eq!(overflowing, before);
}

#[test]
fn malformed_profiles_and_options_are_rejected() {
    let base = profile();
    for change in 0..6 {
        let mut invalid = base.clone();
        match change {
            0 => invalid.peaks[2].mz = invalid.peaks[1].mz,
            1 => invalid.peaks.swap(1, 2),
            2 => invalid.peaks[0].mz = -1.0,
            3 => invalid.peaks[2].mz = f64::INFINITY,
            4 => invalid.peaks[2].intensity = f32::NAN,
            _ => invalid
                .float_data_arrays
                .push(DataArray::new("bad", vec![1.0])),
        }
        assert!(picker().pick_spectrum(&invalid).is_err());
    }
    for invalid in [
        PeakPickerIterative {
            iterations: 0,
            ..picker()
        },
        PeakPickerIterative {
            signal_to_noise: -1.0,
            ..picker()
        },
        PeakPickerIterative {
            peak_width: f64::NAN,
            ..picker()
        },
        PeakPickerIterative {
            spacing_difference: -1.0,
            ..picker()
        },
        PeakPickerIterative {
            max_points: 0,
            ..picker()
        },
        PeakPickerIterative {
            max_work: 0,
            ..picker()
        },
    ] {
        assert!(invalid.pick_spectrum(&base).is_err());
    }
}

#[test]
fn point_iteration_noise_setup_and_refinement_work_are_bounded() {
    let input = profile();
    for limited in [
        PeakPickerIterative {
            max_points: 6,
            ..picker()
        },
        PeakPickerIterative {
            iterations: usize::MAX,
            ..picker()
        },
        PeakPickerIterative {
            max_work: 2,
            ..picker()
        },
        PeakPickerIterative {
            max_work: 30,
            ..Default::default()
        },
    ] {
        assert!(limited.pick_spectrum(&input).is_err());
    }
    let mut limited = PeakPickerIterative {
        max_work: 100,
        ..Default::default()
    };
    limited.noise_estimator.bin_count = 101;
    // Boundary histogram setup is checked even when no seeds survive.
    assert!(limited.pick_spectrum(&input).is_err());
}

#[test]
fn empty_short_and_flat_profiles_return_empty_annotated_centroids() {
    for n in [0, 1, 2, 3, 4, 20] {
        let input = MSSpectrum {
            peaks: (0..n).map(|i| Peak1D::new(i as f64, 1.0)).collect(),
            name: "short or flat".into(),
            ..Default::default()
        };
        let output = picker().pick_spectrum(&input).unwrap();
        assert!(output.picked.spectrum.is_empty());
        assert_eq!(output.picked.spectrum.name, input.name);
        assert_eq!(output.picked.spectrum.spectrum_type, SpectrumType::Centroid);
        assert_eq!(output.picked.spectrum.float_data_arrays.len(), 3);
        assert!(output.regions.is_empty());
    }
}

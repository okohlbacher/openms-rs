// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `IMDataArrayUtils_test.cpp` and of the `setIMUnit`,
//! `reshapeIMFrameToMany` and `reshapeIMFrameToSingle` sections of
//! `IMDataConverter_test.cpp` (core SDK bc9cc12), cited by source line.
//! Literals are tier 3, transcribed from the class tests; values derived from
//! the source instead say so.

use openms::concept::constants::user_param;
use openms::concept::math_functions::ppm_to_mass;
use openms::kernel::DataArray;
use openms::kernel::im_data_array_utils::{
    ImDataArrayUtils, MEAN_ION_MOBILITY_ARRAY, RAW_INVERSE_REDUCED_ION_MOBILITY_ARRAY,
};
use openms::kernel::im_data_converter::ImDataConverter;
use openms::metadata::DriftTimeUnit;
use openms::{Error, MSExperiment, MSSpectrum, Peak1D};

// IMDataArrayUtils_test.cpp L19 "unit round trips and compatibility entry points"
#[test]
fn section_unit_round_trips_and_compatibility_entry_points() {
    for expected in [
        DriftTimeUnit::Millisecond,
        DriftTimeUnit::InverseReducedMobility,
    ] {
        let mut direct = DataArray::new("", vec![1.25f32]);
        let mut legacy = DataArray::new("", Vec::<f32>::new());
        ImDataArrayUtils::set_im_unit(&mut direct, expected).unwrap();
        ImDataConverter::set_im_unit(&mut legacy, expected).unwrap();
        assert_eq!(direct.name, legacy.name);
        assert_eq!(direct.data[0], 1.25f32);
        assert_eq!(ImDataArrayUtils::im_unit(&direct), Some(expected));
        assert_eq!(ImDataConverter::im_unit(&direct), Some(expected));
    }
    let mut array = DataArray::new("", Vec::<f32>::new());
    for unit in [
        DriftTimeUnit::None,
        DriftTimeUnit::FaimsCompensationVoltage,
        DriftTimeUnit::CollisionCrossSection,
    ] {
        assert!(matches!(
            ImDataArrayUtils::set_im_unit(&mut array, unit),
            Err(Error::InvalidValue(_))
        ));
    }
    // Derived: a refused unit leaves the array's name alone.
    assert_eq!(array.name, "");
}

// IMDataArrayUtils_test.cpp L41 "ontology and vendor names retain their units"
#[test]
fn section_ontology_and_vendor_names_retain_their_units() {
    let mut array = DataArray::new("", Vec::<f32>::new());
    array.name = "mean inverse reduced ion mobility array".into();
    assert_eq!(
        ImDataArrayUtils::im_unit(&array),
        Some(DriftTimeUnit::InverseReducedMobility)
    );
    array.name = user_param::ION_MOBILITY.into();
    assert_eq!(
        ImDataArrayUtils::im_unit(&array),
        Some(DriftTimeUnit::Millisecond)
    );
    array.name = user_param::INVERSE_REDUCED_ION_MOBILITY.into();
    assert_eq!(
        ImDataArrayUtils::im_unit(&array),
        Some(DriftTimeUnit::InverseReducedMobility)
    );
    array.name = format!("{} MS:1002954", user_param::ION_MOBILITY);
    assert_eq!(
        ImDataArrayUtils::im_unit(&array),
        Some(DriftTimeUnit::CollisionCrossSection)
    );
    array.name = "unrelated auxiliary data".into();
    assert_eq!(ImDataArrayUtils::im_unit(&array), None);
}

/// Derived from `IMDataArrayUtils::setIMUnit`: the two CV names it writes.
#[test]
fn set_im_unit_writes_the_two_cv_names() {
    let mut array = DataArray::new("old", vec![0.5f32]);
    ImDataArrayUtils::set_im_unit(&mut array, DriftTimeUnit::Millisecond).unwrap();
    assert_eq!(array.name, MEAN_ION_MOBILITY_ARRAY);
    assert_eq!(array.name, "mean ion mobility array");
    ImDataArrayUtils::set_im_unit(&mut array, DriftTimeUnit::InverseReducedMobility).unwrap();
    assert_eq!(array.name, RAW_INVERSE_REDUCED_ION_MOBILITY_ARRAY);
    assert_eq!(array.name, "raw inverse reduced ion mobility array");
    assert_eq!(array.data, [0.5f32]);
}

// IMDataConverter_test.cpp L79 static void setIMUnit(...), L100 getIMUnit (NOT_TESTABLE there)
#[test]
fn converter_section_set_im_unit() {
    let mut fda = DataArray::new("", Vec::<f32>::new());
    assert!(matches!(
        ImDataConverter::set_im_unit(&mut fda, DriftTimeUnit::FaimsCompensationVoltage),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        ImDataConverter::set_im_unit(&mut fda, DriftTimeUnit::None),
        Err(Error::InvalidValue(_))
    ));
    ImDataConverter::set_im_unit(&mut fda, DriftTimeUnit::Millisecond).unwrap();
    assert_eq!(
        ImDataConverter::im_unit(&fda),
        Some(DriftTimeUnit::Millisecond)
    );
    ImDataConverter::set_im_unit(&mut fda, DriftTimeUnit::InverseReducedMobility).unwrap();
    assert_eq!(
        ImDataConverter::im_unit(&fda),
        Some(DriftTimeUnit::InverseReducedMobility)
    );
    fda.name = "mean inverse reduced ion mobility array".into();
    assert_eq!(
        ImDataConverter::im_unit(&fda),
        Some(DriftTimeUnit::InverseReducedMobility)
    );
}

/// The class test's shared `frame` (IMDataConverter_test.cpp:105-122).
fn frame() -> MSSpectrum {
    let peak = |mz: f64, intensity: f32| Peak1D { mz, intensity };
    let mut frame = MSSpectrum {
        rt: 1.0,
        ..MSSpectrum::default()
    };
    frame.peaks = vec![
        peak(1.0, 11.0),
        peak(1.0 + ppm_to_mass(4.0, 1.0).unwrap(), 12.0),
        peak(1.2, 13.0),
        peak(2.0, 20.0),
        peak(3.0 - ppm_to_mass(3.0, 3.0).unwrap(), 32.0),
        peak(3.0, 31.0),
        peak(4.0, 40.0),
        peak(5.0, 50.0),
        peak(6.0, 60.0),
        peak(7.0, 70.0),
    ];
    // `afa.assign({1.1, 1.11, ...})` narrows each double literal to float.
    let values: Vec<f32> = [1.1, 1.11, 1.11, 2.2, 3.2, 3.22, 4.4, 5.6, 6.6, 7.7]
        .iter()
        .map(|value: &f64| *value as f32)
        .collect();
    let mut afa = DataArray::new("", values);
    ImDataConverter::set_im_unit(&mut afa, DriftTimeUnit::Millisecond).unwrap();
    frame.float_data_arrays.push(afa);
    frame
}

/// The class test's shared `spec` (IMDataConverter_test.cpp:124-128).
fn spec() -> MSSpectrum {
    MSSpectrum {
        rt: 2.0,
        peaks: vec![
            Peak1D {
                mz: 111.0,
                intensity: -1.0,
            },
            Peak1D {
                mz: 222.0,
                intensity: -2.0,
            },
            Peak1D {
                mz: 333.0,
                intensity: -3.0,
            },
        ],
        ..MSSpectrum::default()
    }
}

// IMDataConverter_test.cpp L131 static MSExperiment reshapeIMFrameToMany(MSSpectrum im_frame),
// with L197 reshapeIMFrameToSingle ("tested above")
#[test]
fn converter_section_reshape_im_frame_to_many_and_back() {
    // Not an IM frame: Exception::MissingInformation.
    assert!(matches!(
        ImDataConverter::reshape_im_frame_to_many(&spec()),
        Err(Error::MissingInformation(_))
    ));

    let frame = frame();
    let exp = ImDataConverter::reshape_im_frame_to_many(&frame).unwrap();
    assert_eq!(exp.len(), 9); // nine different IM values
    assert_eq!(exp.spectra[0].peaks.len(), 1);
    assert_eq!(exp.spectra[1].peaks.len(), 2);
    assert_eq!(exp.spectra[1].peaks[0].intensity, 12.0f32);
    assert_eq!(exp.spectra[1].peaks[1].intensity, 13.0f32);
    // TEST_EQUAL(exp[0].getDriftTime(), 1.1f): the float widened to double.
    assert_eq!(exp.spectra[0].drift_time, f64::from(1.1f32));
    assert_eq!(exp.spectra[0].drift_time_unit, DriftTimeUnit::Millisecond);
    assert_eq!(exp.spectra[0].rt, 1.0);
    assert_eq!(exp.spectra[1].drift_time, f64::from(1.11f32));
    assert_eq!(exp.spectra[8].drift_time, f64::from(7.7f32));
    assert!(exp.is_im_frame().unwrap());

    let frame_reconstruct = ImDataConverter::reshape_im_frame_to_single(&exp).unwrap();
    assert_eq!(frame_reconstruct.len(), 1);
    assert_eq!(frame_reconstruct.spectra[0], frame);
}

/// Derived from `reshapeIMFrameToMany`: an unsorted frame is sorted by ion
/// mobility first (stably), the output is m/z-sorted within each value, and the
/// input is untouched; an empty frame needs no ion-mobility array.
#[test]
fn reshape_to_many_sorts_an_unsorted_frame_and_leaves_the_input_alone() {
    let mut frame = MSSpectrum {
        rt: 5.0,
        ..MSSpectrum::default()
    };
    frame.peaks = vec![
        Peak1D {
            mz: 300.0,
            intensity: 1.0,
        },
        Peak1D {
            mz: 100.0,
            intensity: 2.0,
        },
        Peak1D {
            mz: 200.0,
            intensity: 3.0,
        },
    ];
    let mut mobility = DataArray::new("", vec![0.9f32, 0.7, 0.9]);
    ImDataConverter::set_im_unit(&mut mobility, DriftTimeUnit::InverseReducedMobility).unwrap();
    frame.float_data_arrays.push(mobility);
    let before = frame.clone();
    let many = ImDataConverter::reshape_im_frame_to_many(&frame).unwrap();
    assert_eq!(frame, before);
    assert_eq!(many.len(), 2);
    assert_eq!(
        many.spectra[0].peaks,
        [Peak1D {
            mz: 100.0,
            intensity: 2.0
        }]
    );
    assert_eq!(
        many.spectra[1].peaks,
        [
            Peak1D {
                mz: 200.0,
                intensity: 3.0
            },
            Peak1D {
                mz: 300.0,
                intensity: 1.0
            }
        ]
    );
    assert_eq!(
        many.spectra[1].drift_time_unit,
        DriftTimeUnit::InverseReducedMobility
    );
    assert!(many.spectra.iter().all(|s| s.float_data_arrays.is_empty()));

    assert!(
        ImDataConverter::reshape_im_frame_to_many(&MSSpectrum::default())
            .unwrap()
            .is_empty()
    );

    // An ion-mobility array shorter than the peak list is refused, where the
    // source would silently drop the surplus peaks.
    let mut short = frame.clone();
    short.float_data_arrays[0].data.pop();
    assert!(matches!(
        ImDataConverter::reshape_im_frame_to_many(&short),
        Err(Error::InvalidValue(_))
    ));
}

/// Derived from `reshapeIMFrameToSingle`: spectra without a drift time are
/// copied and end the frame, frames split on retention time, and a FAIMS unit
/// cannot name an ion-mobility array.
#[test]
fn reshape_to_single_frames_by_retention_time() {
    let with_drift = |rt: f64, drift: f64, mz: f64| MSSpectrum {
        rt,
        drift_time: drift,
        drift_time_unit: DriftTimeUnit::Millisecond,
        peaks: vec![Peak1D { mz, intensity: 1.0 }],
        ..MSSpectrum::default()
    };
    let mut experiment = MSExperiment::new();
    experiment.spectra = vec![
        with_drift(1.0, 0.5, 100.0),
        with_drift(1.0, 0.6, 200.0),
        spec(),
        with_drift(3.0, 0.5, 300.0),
    ];
    let before = experiment.clone();
    let single = ImDataConverter::reshape_im_frame_to_single(&experiment).unwrap();
    assert_eq!(experiment, before);
    assert_eq!(single.len(), 3);
    assert_eq!(single.spectra[0].peaks.len(), 2);
    assert_eq!(single.spectra[0].float_data_arrays[0].data, [0.5f32, 0.6]);
    assert_eq!(single.spectra[0].drift_time, -1.0);
    assert_eq!(single.spectra[0].drift_time_unit, DriftTimeUnit::None);
    assert_eq!(single.spectra[1], spec());
    assert_eq!(single.spectra[2].float_data_arrays[0].data, [0.5f32]);

    let mut faims = MSExperiment::new();
    let mut spectrum = with_drift(1.0, -45.0, 100.0);
    spectrum.drift_time_unit = DriftTimeUnit::FaimsCompensationVoltage;
    faims.spectra.push(spectrum);
    assert!(matches!(
        ImDataConverter::reshape_im_frame_to_single(&faims),
        Err(Error::InvalidValue(_))
    ));
    assert!(
        ImDataConverter::reshape_im_frame_to_single(&MSExperiment::new())
            .unwrap()
            .is_empty()
    );
}

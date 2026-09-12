// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `SpectrumHelper_test.cpp` (nine sections) plus native checks.
//! Literals are transcribed from the pinned class test (evidence tier 3).

use openms::kernel::spectrum_helper::{
    IntensityAveragingMethod as Method, PeakContainer, SpectrumHelperLimits, UniquePositionOptions,
    copy_spectrum_meta, data_array_by_name, data_array_by_name_mut, data_array_index_by_name,
    make_peak_position_unique, make_peak_position_unique_with, subtract_minimum_intensity,
    subtract_minimum_intensity_with_limits,
};
use openms::kernel::{ChromatogramPeak, DataArray, MSChromatogram, MSSpectrum, Peak1D};

fn close(actual: f32, expected: f64) {
    let actual = f64::from(actual);
    assert!(
        (actual - expected).abs() <= 1e-5 * expected.abs().max(1.0),
        "{actual:.17} != {expected:.17}"
    );
}

/// Bitwise peak comparison plus ordinary equality of every other field, for
/// fixtures that deliberately contain NaN (which is never `==` to itself).
fn assert_unchanged<C: PeakContainer + Clone + PartialEq + std::fmt::Debug>(
    actual: &C,
    before: &C,
) {
    let bits = |c: &C| -> Vec<(u64, u32)> {
        c.peaks()
            .iter()
            .map(|p| (C::position(p).to_bits(), C::intensity(p).to_bits()))
            .collect()
    };
    assert_eq!(bits(actual), bits(before));
    let mut actual = actual.clone();
    actual.peaks_mut().clear();
    let mut before = before.clone();
    before.peaks_mut().clear();
    assert_eq!(actual, before);
}

/// Source: three copies of the same array named f1/f2/f3.
fn three<T: Clone>(values: Vec<T>) -> Vec<DataArray<T>> {
    ["f1", "f2", "f3"]
        .iter()
        .map(|name| DataArray::new(*name, values.clone()))
        .collect()
}

fn check_lookup<T: Clone + PartialEq + std::fmt::Debug>(arrays: &mut [DataArray<T>]) {
    assert!(data_array_index_by_name(arrays, "f2").is_some());
    assert_eq!(data_array_index_by_name(arrays, "NOT_THERE"), None);
    assert_eq!(data_array_index_by_name(arrays, "f1"), Some(0));
    assert_eq!(data_array_index_by_name(arrays, "f2"), Some(1));
    assert_eq!(data_array_index_by_name(arrays, "f3"), Some(2));
    // Const overload: borrowed lookup.
    assert_eq!(
        data_array_by_name(arrays, "f2").map(|a| a.name.as_str()),
        Some("f2")
    );
    assert_eq!(
        data_array_by_name(arrays, "f1").map(|a| a.name.as_str()),
        Some("f1")
    );
    assert!(data_array_by_name(arrays, "NOT_THERE").is_none());
    assert_eq!(data_array_by_name(arrays, "f3").unwrap().data.len(), 5);
    // Mutable overload: the borrow is the same element as the index.
    let expected = arrays[2].data.clone();
    assert_eq!(data_array_by_name_mut(arrays, "f3").unwrap().data, expected);
    assert!(data_array_by_name_mut(arrays, "NOT_THERE").is_none());
    // Name comparison is exact; the first match wins.
    arrays[1].name = "f1".into();
    assert_eq!(data_array_index_by_name(arrays, "f1"), Some(0));
    assert_eq!(data_array_index_by_name(arrays, "F1"), None);
}

// Section 1: MSSpectrum::FloatDataArrays getDataArrayByName
#[test]
fn spectrum_float_data_array_by_name() {
    let mut ds = MSSpectrum::new();
    ds.float_data_arrays = three(vec![56.0_f32, 201.0, 31.0, 201.0, 201.0]);
    check_lookup(&mut ds.float_data_arrays);
}

// Section 2: MSSpectrum::StringDataArrays getDataArrayByName
#[test]
fn spectrum_string_data_array_by_name() {
    let mut ds = MSSpectrum::new();
    ds.string_data_arrays = three(
        ["56", "201", "31", "201", "201"]
            .iter()
            .map(|s| s.to_string())
            .collect(),
    );
    check_lookup(&mut ds.string_data_arrays);
}

// Section 3: MSSpectrum::IntegerDataArrays getDataArrayByName
#[test]
fn spectrum_integer_data_array_by_name() {
    let mut ds = MSSpectrum::new();
    ds.integer_data_arrays = three(vec![56_i32, 201, 31, 201, 201]);
    check_lookup(&mut ds.integer_data_arrays);
}

// Section 4: MSChromatogram::FloatDataArrays getDataArrayByName
#[test]
fn chromatogram_float_data_array_by_name() {
    let mut ds = MSChromatogram::new();
    ds.float_data_arrays = three(vec![56.0_f32, 201.0, 31.0, 201.0, 201.0]);
    check_lookup(&mut ds.float_data_arrays);
}

// Section 5: MSChromatogram::StringDataArrays getDataArrayByName
#[test]
fn chromatogram_string_data_array_by_name() {
    let mut ds = MSChromatogram::new();
    ds.string_data_arrays = three(
        ["56", "201", "31", "201", "201"]
            .iter()
            .map(|s| s.to_string())
            .collect(),
    );
    check_lookup(&mut ds.string_data_arrays);
}

// Section 6: MSChromatogram::IntegerDataArrays getDataArrayByName
#[test]
fn chromatogram_integer_data_array_by_name() {
    let mut ds = MSChromatogram::new();
    ds.integer_data_arrays = three(vec![56_i32, 201, 31, 201, 201]);
    check_lookup(&mut ds.integer_data_arrays);
}

// Section 7: removePeaks(p, pos_start, pos_end). The source function maps onto
// the existing `retain_peaks` with the inclusive PosBegin/PosEnd range; the
// section's literals are asserted through that mapping.
#[test]
fn remove_peaks_maps_to_retain_peaks_with_inclusive_range() {
    let mut s = MSSpectrum::new();
    let mut c = MSChromatogram::new();
    let mut ida = DataArray::new("ida", Vec::new());
    for i in 5..15_i32 {
        // RTs: [5 14]
        s.peaks.push(Peak1D::new(f64::from(i), 0.0));
        c.peaks.push(ChromatogramPeak::new(f64::from(i), 0.0));
        ida.data.push(i);
    }
    s.integer_data_arrays.push(ida.clone());
    c.integer_data_arrays.push(ida);

    // start rt (3) is lower than the minimum (5) within the spectrum
    let mut s1 = s.clone();
    s1.retain_peaks(|p| p.mz >= 3.0 && p.mz <= 6.0).unwrap();
    assert_eq!(s1.len(), 2);
    assert_eq!(s1.integer_data_arrays[0].data.len(), 2);
    assert_eq!(s1.peaks[0].mz, 5.0);
    assert_eq!(s1.peaks[1].mz, 6.0);
    assert_eq!(s1.integer_data_arrays[0].data, vec![5, 6]);

    // no peak within the requested range
    let mut s2 = s.clone();
    s2.retain_peaks(|p| p.mz >= 0.0 && p.mz <= 4.0).unwrap();
    assert_eq!(s2.len(), 0);
    assert_eq!(s2.integer_data_arrays[0].data.len(), 0);

    // end rt (16) is higher than the maximum (14) within the chromatogram
    let mut c1 = c.clone();
    c1.retain_peaks(|p| p.rt >= 12.0 && p.rt <= 16.0).unwrap();
    assert_eq!(c1.len(), 3);
    assert_eq!(c1.integer_data_arrays[0].data.len(), 3);
    assert_eq!(c1.peaks[0].rt, 12.0);
    assert_eq!(c1.peaks[1].rt, 13.0);
    assert_eq!(c1.peaks[2].rt, 14.0);
    assert_eq!(c1.integer_data_arrays[0].data, vec![12, 13, 14]);

    // all within the range
    let mut c2 = c.clone();
    c2.retain_peaks(|p| p.rt >= 9.0 && p.rt <= 12.0).unwrap();
    assert_eq!(c2.len(), 4);
    assert_eq!(c2.integer_data_arrays[0].data.len(), 4);
    assert_eq!(c2.peaks[0].rt, 9.0);
    assert_eq!(c2.peaks[1].rt, 10.0);
    assert_eq!(c2.peaks[2].rt, 11.0);
    assert_eq!(c2.peaks[3].rt, 12.0);

    let mut s_empty = MSSpectrum::new();
    s_empty
        .retain_peaks(|p| p.mz >= 9.0 && p.mz <= 12.0)
        .unwrap();
    assert_eq!(s_empty.len(), 0);

    // Native difference: the source silently skips arrays whose length differs
    // from the peak count; `retain_peaks` rejects them and leaves everything
    // unchanged.
    let mut s3 = s.clone();
    s3.integer_data_arrays[0].data.push(99);
    let before = s3.clone();
    assert!(s3.retain_peaks(|p| p.mz >= 3.0 && p.mz <= 6.0).is_err());
    assert_eq!(s3, before);
}

// Section 8: subtractMinimumIntensity(p)
#[test]
fn subtract_minimum_intensity_source_literals() {
    let mut s = MSSpectrum::new();
    let mut c = MSChromatogram::new();
    for i in -5..5_i32 {
        // Intensities: [-5 4]
        s.peaks.push(Peak1D::new(0.0, i as f32));
        c.peaks.push(ChromatogramPeak::new(0.0, i as f32));
    }

    subtract_minimum_intensity(&mut s).unwrap();
    close(s.peaks[0].intensity, 0.0);
    close(s.peaks[1].intensity, 1.0);
    close(s.peaks[9].intensity, 9.0);

    subtract_minimum_intensity(&mut c).unwrap();
    close(c.peaks[0].intensity, 0.0);
    close(c.peaks[1].intensity, 1.0);
    close(c.peaks[9].intensity, 9.0);

    let mut c_empty = MSChromatogram::new();
    subtract_minimum_intensity(&mut c_empty).unwrap();
    assert_eq!(c_empty.len(), 0);

    s.clear(true);
    c.clear(true);
    for i in 5..15_i32 {
        // Intensities: [5 14]
        s.peaks.push(Peak1D::new(0.0, i as f32));
        c.peaks.push(ChromatogramPeak::new(0.0, i as f32));
    }

    subtract_minimum_intensity(&mut s).unwrap();
    close(s.peaks[0].intensity, 0.0);
    close(s.peaks[1].intensity, 1.0);
    close(s.peaks[9].intensity, 9.0);

    subtract_minimum_intensity(&mut c).unwrap();
    close(c.peaks[0].intensity, 0.0);
    close(c.peaks[1].intensity, 1.0);
    close(c.peaks[9].intensity, 9.0);
}

#[test]
fn subtract_minimum_intensity_native_boundaries() {
    // Data arrays and metadata are untouched, as the source note states.
    let mut s = MSSpectrum::from_peaks(vec![Peak1D::new(1.0, 3.0), Peak1D::new(2.0, 7.0)]);
    s.float_data_arrays
        .push(DataArray::new("fda", vec![3.0, 7.0]));
    s.rt = 42.0;
    subtract_minimum_intensity(&mut s).unwrap();
    assert_eq!(s.peaks, vec![Peak1D::new(1.0, 0.0), Peak1D::new(2.0, 4.0)]);
    assert_eq!(s.float_data_arrays[0].data, vec![3.0, 7.0]);
    assert_eq!(s.rt, 42.0);

    // Non-finite intensities and exceeded limits are rejected atomically.
    let mut bad = MSSpectrum::from_peaks(vec![Peak1D::new(1.0, f32::NAN), Peak1D::new(2.0, 1.0)]);
    let before = bad.clone();
    assert!(subtract_minimum_intensity(&mut bad).is_err());
    assert_unchanged(&bad, &before);

    let mut overflow = MSSpectrum::from_peaks(vec![
        Peak1D::new(1.0, -f32::MAX),
        Peak1D::new(2.0, f32::MAX),
    ]);
    let before = overflow.clone();
    assert!(subtract_minimum_intensity(&mut overflow).is_err());
    assert_eq!(overflow, before);

    let mut limited = MSSpectrum::from_peaks(vec![Peak1D::new(1.0, 3.0), Peak1D::new(2.0, 7.0)]);
    let before = limited.clone();
    let limits = SpectrumHelperLimits {
        max_peaks: 1,
        ..Default::default()
    };
    assert!(subtract_minimum_intensity_with_limits(&mut limited, limits).is_err());
    assert_eq!(limited, before);
    let limits = SpectrumHelperLimits {
        max_work: 5,
        ..Default::default()
    };
    assert!(subtract_minimum_intensity_with_limits(&mut limited, limits).is_err());
    assert_eq!(limited, before);
    let limits = SpectrumHelperLimits {
        max_peaks: 2,
        max_work: 6,
    };
    subtract_minimum_intensity_with_limits(&mut limited, limits).unwrap();
    assert_eq!(limited.peaks[1].intensity, 4.0);
}

fn s_template() -> MSSpectrum {
    MSSpectrum::from_peaks(vec![
        Peak1D::new(1.0, 1.0),
        Peak1D::new(2.0, 4.0),
        Peak1D::new(2.0, 8.0),
        Peak1D::new(3.0, 9.0),
        Peak1D::new(4.0, 7.0),
        Peak1D::new(2.0, 10.0),
    ])
}

fn c_template() -> MSChromatogram {
    MSChromatogram::from_peaks(vec![
        ChromatogramPeak::new(1.0, 1.0),
        ChromatogramPeak::new(2.0, 4.0),
        ChromatogramPeak::new(2.0, 8.0),
        ChromatogramPeak::new(3.0, 9.0),
        ChromatogramPeak::new(4.0, 7.0),
        ChromatogramPeak::new(2.0, 10.0),
    ])
}

fn intensities<C: PeakContainer>(container: &C) -> Vec<f32> {
    container.peaks().iter().map(C::intensity).collect()
}

fn positions<C: PeakContainer>(container: &C) -> Vec<f64> {
    container.peaks().iter().map(C::position).collect()
}

// Section 9: makePeakPositionUnique(p, m)
#[test]
fn make_peak_position_unique_source_literals() {
    let cases: [(Method, [f64; 4]); 5] = [
        (Method::Median, [1.0, 8.0, 9.0, 7.0]),
        (Method::Mean, [1.0, (4.0 + 8.0 + 10.0) / 3.0, 9.0, 7.0]),
        (Method::Sum, [1.0, 4.0 + 8.0 + 10.0, 9.0, 7.0]),
        (Method::Min, [1.0, 4.0, 9.0, 7.0]),
        (Method::Max, [1.0, 10.0, 9.0, 7.0]),
    ];
    for (method, expected) in cases {
        let mut s = s_template();
        make_peak_position_unique(&mut s, method).unwrap();
        assert_eq!(s.len(), 4, "{method:?}");
        for (actual, expected) in intensities(&s).into_iter().zip(expected) {
            close(actual, expected);
        }
        assert_eq!(positions(&s), vec![1.0, 2.0, 3.0, 4.0]);

        let mut c = c_template();
        make_peak_position_unique(&mut c, method).unwrap();
        assert_eq!(c.len(), 4, "{method:?}");
        for (actual, expected) in intensities(&c).into_iter().zip(expected) {
            close(actual, expected);
        }
        assert_eq!(positions(&c), vec![1.0, 2.0, 3.0, 4.0]);
    }

    // Source default argument is MEDIAN.
    assert_eq!(Method::default(), Method::Median);

    let mut s_empty = MSSpectrum::new();
    make_peak_position_unique(&mut s_empty, Method::default()).unwrap();
    assert_eq!(s_empty.len(), 0);

    let mut c_empty = MSChromatogram::new();
    make_peak_position_unique(&mut c_empty, Method::default()).unwrap();
    assert_eq!(c_empty.len(), 0);
}

#[test]
fn make_peak_position_unique_median_and_sum_conventions() {
    // Even group size averages the two middle values (Math::median).
    let mut s = MSSpectrum::from_peaks(vec![
        Peak1D::new(5.0, 9.0),
        Peak1D::new(5.0, 1.0),
        Peak1D::new(5.0, 4.0),
        Peak1D::new(5.0, 2.0),
    ]);
    make_peak_position_unique(&mut s, Method::Median).unwrap();
    assert_eq!(s.peaks, vec![Peak1D::new(5.0, 3.0)]);

    // Groups are separated by a strict `>` after a stable sort, so -0.0 and
    // 0.0 share a group and the group keeps the first position seen.
    let mut z = MSSpectrum::from_peaks(vec![
        Peak1D::new(0.0, 2.0),
        Peak1D::new(-0.0, 4.0),
        Peak1D::new(1.0, 1.0),
    ]);
    make_peak_position_unique(&mut z, Method::Sum).unwrap();
    assert_eq!(z.len(), 2);
    assert_eq!(z.peaks[0].intensity, 6.0);
    assert!(z.peaks[0].mz == 0.0 && z.peaks[0].mz.is_sign_positive());

    // The sum accumulates in f64 in storage order and narrows once: three
    // f32 values whose pairwise f32 sums would round differently.
    let mut f = MSSpectrum::from_peaks(vec![
        Peak1D::new(1.0, 16_777_216.0),
        Peak1D::new(1.0, 1.0),
        Peak1D::new(1.0, 1.0),
    ]);
    make_peak_position_unique(&mut f, Method::Sum).unwrap();
    assert_eq!(f.peaks[0].intensity, 16_777_218.0);
    // Mean of the same group: 16777218 / 3 in f64, then narrowed.
    let mut m = MSSpectrum::from_peaks(vec![
        Peak1D::new(1.0, 16_777_216.0),
        Peak1D::new(1.0, 1.0),
        Peak1D::new(1.0, 1.0),
    ]);
    make_peak_position_unique(&mut m, Method::Mean).unwrap();
    assert_eq!(m.peaks[0].intensity, (16_777_218.0_f64 / 3.0) as f32);
}

#[test]
fn make_peak_position_unique_keeps_metadata_and_refuses_arrays_by_default() {
    let mut s = s_template();
    s.rt = 12.5;
    s.ms_level = 2;
    s.name = "merged".into();
    s.native_id = "scan=7".into();
    make_peak_position_unique(&mut s, Method::Max).unwrap();
    assert_eq!(s.rt, 12.5);
    assert_eq!(s.ms_level, 2);
    assert_eq!(s.name, "merged");
    assert_eq!(s.native_id, "scan=7");
    assert_eq!(intensities(&s), vec![1.0, 10.0, 9.0, 7.0]);

    // Attached arrays (even empty placeholders) are refused; nothing changes.
    let mut with_arrays = s_template();
    with_arrays
        .string_data_arrays
        .push(DataArray::new("labels", Vec::new()));
    let before = with_arrays.clone();
    assert!(make_peak_position_unique(&mut with_arrays, Method::Median).is_err());
    assert_eq!(with_arrays, before);

    // An empty container returns before the policy check, as in the source.
    let mut empty = MSSpectrum::new();
    empty
        .float_data_arrays
        .push(DataArray::new("fda", Vec::new()));
    empty.rt = 3.0;
    let before = empty.clone();
    make_peak_position_unique(&mut empty, Method::Median).unwrap();
    assert_eq!(empty, before);

    // Explicit discard keeps metadata and drops the arrays.
    let mut discard = s_template();
    discard.rt = 9.0;
    discard
        .integer_data_arrays
        .push(DataArray::new("ida", vec![1, 2, 3, 4, 5, 6]));
    let options = UniquePositionOptions {
        discard_data_arrays: true,
        ..Default::default()
    };
    make_peak_position_unique_with(&mut discard, Method::Sum, options).unwrap();
    assert_eq!(discard.rt, 9.0);
    assert!(discard.integer_data_arrays.is_empty());
    assert_eq!(intensities(&discard), vec![1.0, 22.0, 9.0, 7.0]);

    // Source option set: default-constructed container with merged peaks only.
    let mut source = s_template();
    source.rt = 9.0;
    source.ms_level = 3;
    source.name = "gone".into();
    source
        .float_data_arrays
        .push(DataArray::new("fda", vec![0.0; 6]));
    make_peak_position_unique_with(&mut source, Method::Min, UniquePositionOptions::source())
        .unwrap();
    let mut expected = MSSpectrum::from_peaks(vec![
        Peak1D::new(1.0, 1.0),
        Peak1D::new(2.0, 4.0),
        Peak1D::new(3.0, 9.0),
        Peak1D::new(4.0, 7.0),
    ]);
    expected.rt = MSSpectrum::default().rt;
    assert_eq!(source, expected);
    assert_eq!(source.rt, -1.0);
    assert_eq!(source.ms_level, 1);
    assert_eq!(source.name, "");

    let mut chrom = c_template();
    chrom.name = "chrom".into();
    chrom
        .string_data_arrays
        .push(DataArray::new("s", vec![String::new(); 6]));
    make_peak_position_unique_with(&mut chrom, Method::Median, UniquePositionOptions::source())
        .unwrap();
    assert_eq!(chrom.name, "");
    assert!(chrom.string_data_arrays.is_empty());
    assert_eq!(intensities(&chrom), vec![1.0, 8.0, 9.0, 7.0]);
}

#[test]
fn make_peak_position_unique_rejects_invalid_input_atomically() {
    let mut nan = s_template();
    nan.peaks[2].mz = f64::NAN;
    let before = nan.clone();
    assert!(make_peak_position_unique(&mut nan, Method::Median).is_err());
    assert_unchanged(&nan, &before);

    let mut inf = c_template();
    inf.peaks[1].intensity = f32::INFINITY;
    let before = inf.clone();
    assert!(make_peak_position_unique(&mut inf, Method::Sum).is_err());
    assert_eq!(inf, before);

    // A sum overflowing f32 is a non-finite merged intensity.
    let mut overflow =
        MSSpectrum::from_peaks(vec![Peak1D::new(1.0, f32::MAX), Peak1D::new(1.0, f32::MAX)]);
    let before = overflow.clone();
    assert!(make_peak_position_unique(&mut overflow, Method::Sum).is_err());
    assert_eq!(overflow, before);
    make_peak_position_unique(&mut overflow, Method::Max).unwrap();
    assert_eq!(overflow.peaks, vec![Peak1D::new(1.0, f32::MAX)]);

    // Limits are checked before mutation.
    let mut limited = s_template();
    let before = limited.clone();
    let mut options = UniquePositionOptions::default();
    options.limits.max_peaks = 5;
    assert!(make_peak_position_unique_with(&mut limited, Method::Median, options).is_err());
    assert_eq!(limited, before);
    options.limits = SpectrumHelperLimits {
        max_peaks: 6,
        max_work: 100,
    };
    assert!(make_peak_position_unique_with(&mut limited, Method::Median, options).is_err());
    assert_eq!(limited, before);
    // 2 * 6 + 8 * 6 * bit_length(6) = 12 + 144 = 156 units suffice.
    options.limits.max_work = 156;
    make_peak_position_unique_with(&mut limited, Method::Median, options).unwrap();
    assert_eq!(limited.len(), 4);
}

// The source header declares copySpectrumMeta but the class test has no section
// for it; this is a native test of the .cpp behaviour.
#[test]
fn copy_spectrum_meta_copies_metadata_only() {
    let mut input = s_template();
    input.rt = 100.5;
    input.ms_level = 2;
    input.name = "in".into();
    input.native_id = "scan=1".into();
    input.metadata.insert("key".into(), "value".into());
    // Source copySpectrumMeta assigns drift time and its unit explicitly.
    input.drift_time = 12.5;
    input.drift_time_unit = openms::metadata::DriftTimeUnit::Millisecond;
    input
        .float_data_arrays
        .push(DataArray::new("fda", vec![1.0; 6]));

    // clear_spectrum = true: output holds only the input's metadata.
    let mut output = spectrum_with_data();
    copy_spectrum_meta(&input, &mut output, true);
    assert!(output.peaks.is_empty());
    assert!(output.float_data_arrays.is_empty());
    assert!(output.integer_data_arrays.is_empty());
    assert!(output.string_data_arrays.is_empty());
    assert_eq!(output.rt, 100.5);
    assert_eq!(output.ms_level, 2);
    assert_eq!(output.name, "in");
    assert_eq!(output.native_id, "scan=1");
    assert_eq!(output.metadata["key"].as_str().unwrap(), "value");
    assert_eq!(output.drift_time, 12.5);
    assert_eq!(
        output.drift_time_unit,
        openms::metadata::DriftTimeUnit::Millisecond
    );
    let mut expected = input.clone();
    expected.peaks.clear();
    expected.float_data_arrays.clear();
    assert_eq!(output, expected);

    // clear_spectrum = false: output keeps its peaks and arrays.
    let mut output = spectrum_with_data();
    let kept_peaks = output.peaks.clone();
    let kept_arrays = output.integer_data_arrays.clone();
    copy_spectrum_meta(&input, &mut output, false);
    assert_eq!(output.peaks, kept_peaks);
    assert_eq!(output.integer_data_arrays, kept_arrays);
    assert!(output.float_data_arrays.is_empty());
    assert_eq!(output.rt, 100.5);
    assert_eq!(output.ms_level, 2);
    assert_eq!(output.name, "in");
    assert_eq!(output.metadata["key"].as_str().unwrap(), "value");

    // The input is not modified and the source default (clear) is explicit.
    assert_eq!(input.peaks.len(), 6);
    assert_eq!(input.float_data_arrays.len(), 1);
}

/// Output fixture for the copy test: one peak, one integer array, own metadata.
fn spectrum_with_data() -> MSSpectrum {
    let mut spectrum = MSSpectrum::from_peaks(vec![Peak1D::new(7.0, 7.0)]);
    spectrum.rt = 1.0;
    spectrum.ms_level = 1;
    spectrum.name = "out".into();
    spectrum
        .integer_data_arrays
        .push(DataArray::new("ida", vec![7]));
    spectrum
}

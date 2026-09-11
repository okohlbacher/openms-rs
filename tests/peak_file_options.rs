// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::format::peak_options::{
    DEFAULT_NUMPRESS_ERROR_TOLERANCE, EMPTY_PEAK_FILE_RANGE, MAX_PEAK_FILE_MS_LEVELS,
    NUMPRESS_MASS_TIME_WARNING, NumpressCompression, NumpressConfig, PeakFileOptions,
};
use openms::kernel::NumericRange;

fn range(min: f64, max: f64) -> NumericRange {
    NumericRange { min, max }
}

#[test]
fn complete_header_default_state() {
    let options = PeakFileOptions::new();
    assert!(!options.metadata_only);
    assert!(!options.force_mq_compatibility);
    assert!(!options.force_tpp_compatibility);
    assert!(options.write_supplemental_data);
    assert!(!options.mz_32_bit);
    assert!(options.intensity_32_bit);
    assert!(!options.zlib_compression);
    assert!(!options.always_append_data);
    assert!(!options.skip_xml_checks);
    assert!(options.sort_spectra_by_mz);
    assert!(options.sort_chromatograms_by_rt);
    assert!(options.fill_data);
    assert!(options.write_index);
    assert_eq!(options.max_data_pool_size, 100);
    assert!(options.precursor_mz_selected_ion);
    assert!(!options.skip_chromatograms);
    assert!(!options.has_filters());
    assert!(!options.has_rt_range());
    assert!(!options.has_mz_range());
    assert!(!options.has_intensity_range());
    assert!(!options.has_precursor_mz_range());
    for actual in [
        options.rt_range(),
        options.mz_range(),
        options.intensity_range(),
        options.precursor_mz_range(),
    ] {
        assert_eq!(actual.min.to_bits(), f64::MAX.to_bits());
        assert_eq!(actual.max.to_bits(), (-f64::MAX).to_bits());
    }
    assert!(!options.has_ms_levels());
    assert!(!options.contains_ms_level(1));
    assert!(options.ms_levels().is_empty());
    assert_eq!(
        options.numpress_configuration_mass_time(),
        NumpressConfig::default()
    );
    assert_eq!(
        options.numpress_configuration_intensity(),
        NumpressConfig::default()
    );
    assert_eq!(
        options.numpress_configuration_float_data_array(),
        NumpressConfig::default()
    );
    assert_eq!(options, PeakFileOptions::default());
}

#[test]
fn source_scalar_setters_are_public_fields_and_copy_independently() {
    let mut changed = PeakFileOptions::default();
    changed.metadata_only = true;
    changed.force_mq_compatibility = true;
    changed.force_tpp_compatibility = true;
    changed.write_supplemental_data = false;
    changed.mz_32_bit = true;
    changed.intensity_32_bit = false;
    changed.zlib_compression = true;
    changed.always_append_data = true;
    changed.skip_xml_checks = true;
    changed.sort_spectra_by_mz = false;
    changed.sort_chromatograms_by_rt = false;
    changed.fill_data = false;
    changed.write_index = false;
    changed.max_data_pool_size = 250; // literal source class-test value
    changed.precursor_mz_selected_ion = false;
    changed.skip_chromatograms = true;
    let copied = changed.clone();
    assert_eq!(copied, changed);
    assert!(
        copied.metadata_only && copied.force_mq_compatibility && copied.force_tpp_compatibility
    );
    assert!(!copied.write_supplemental_data && copied.mz_32_bit && !copied.intensity_32_bit);
    assert!(copied.zlib_compression && copied.always_append_data && copied.skip_xml_checks);
    assert!(!copied.sort_spectra_by_mz && !copied.sort_chromatograms_by_rt && !copied.fill_data);
    assert!(!copied.write_index && !copied.precursor_mz_selected_ion && copied.skip_chromatograms);
    assert_eq!(copied.max_data_pool_size, 250);
    // Every switch above is excluded from the source's narrower hasFilters.
    assert!(!copied.has_filters());
    changed.max_data_pool_size = 0;
    assert_eq!(changed.max_data_pool_size, 0);
    changed.max_data_pool_size = usize::MAX;
    assert_eq!(changed.max_data_pool_size, usize::MAX);
    assert_eq!(copied.max_data_pool_size, 250);
}

#[test]
fn literal_source_range_and_has_filters_cases() {
    let mut options = PeakFileOptions::default();
    options.set_rt_range(range(2.0, 4.0));
    assert!(options.has_rt_range());
    assert_eq!(options.rt_range(), range(2.0, 4.0));
    options.set_mz_range(range(3.0, 5.0));
    assert!(options.has_mz_range());
    assert_eq!(options.mz_range(), range(3.0, 5.0));
    options.set_intensity_range(range(3.0, 5.0));
    assert!(options.has_intensity_range());
    assert_eq!(options.intensity_range(), range(3.0, 5.0));
    options.set_precursor_mz_range(range(400.0, 1200.0));
    assert!(options.has_precursor_mz_range());
    assert_eq!(options.precursor_mz_range(), range(400.0, 1200.0));
    assert_eq!(options.clone(), options);
    for kind in 0..3 {
        let mut options = PeakFileOptions::default();
        assert!(!options.has_filters());
        match kind {
            0 => options.set_rt_range(range(10.0, 100.0)),
            1 => options.add_ms_level(2).unwrap(),
            _ => options.set_precursor_mz_range(range(400.0, 1200.0)),
        }
        assert!(options.has_filters());
    }
}

#[test]
fn exact_empty_sentinel_activation_is_not_geometrical_emptiness() {
    let mut options = PeakFileOptions::default();
    options.set_rt_range(range(5.0, 5.0));
    assert!(options.has_rt_range());
    assert!(options.has_filters());
    options.set_rt_range(EMPTY_PEAK_FILE_RANGE);
    assert!(!options.has_rt_range());
    assert!(!options.has_filters());
    options.set_mz_range(EMPTY_PEAK_FILE_RANGE);
    options.set_intensity_range(EMPTY_PEAK_FILE_RANGE);
    assert!(options.has_mz_range() && options.has_intensity_range());
    assert_eq!(options.mz_range(), EMPTY_PEAK_FILE_RANGE);
    assert_eq!(options.intensity_range(), EMPTY_PEAK_FILE_RANGE);
    assert!(!options.has_filters());
    options.set_precursor_mz_range(EMPTY_PEAK_FILE_RANGE);
    assert!(options.has_precursor_mz_range());
    assert!(options.has_filters());
    // Range getter equality does not erase the source's independent enabled bits.
    let default = PeakFileOptions::default();
    assert_eq!(default.precursor_mz_range(), options.precursor_mz_range());
    assert_ne!(default, options);
}

#[test]
fn range_value_storage_preserves_signed_zero_nonfinite_and_raw_endpoints() {
    let nan = f64::from_bits(0x7ff8_0000_0000_0042);
    for bounds in [
        range(-8.0, -1.0),
        range(-0.0, 0.0),
        range(f64::NEG_INFINITY, f64::INFINITY),
        range(nan, 5.0),
        range(2.0, nan),
        range(10.0, -5.0),
    ] {
        let mut options = PeakFileOptions::default();
        options.set_rt_range(bounds);
        options.set_mz_range(bounds);
        options.set_intensity_range(bounds);
        options.set_precursor_mz_range(bounds);
        assert!(options.has_rt_range());
        assert!(options.has_mz_range());
        assert!(options.has_intensity_range());
        assert!(options.has_precursor_mz_range());
        for actual in [
            options.rt_range(),
            options.mz_range(),
            options.intensity_range(),
            options.precursor_mz_range(),
        ] {
            assert_eq!(actual.min.to_bits(), bounds.min.to_bits());
            assert_eq!(actual.max.to_bits(), bounds.max.to_bits());
        }
    }
}

#[test]
fn literal_ms_levels_and_signed_ordered_duplicate_membership() {
    let mut options = PeakFileOptions::default();
    options.set_ms_levels(&[1, 3, 5]).unwrap();
    assert!(options.has_ms_levels());
    assert!(options.contains_ms_level(3));
    assert!(!options.contains_ms_level(2));
    assert_eq!(options.ms_levels(), [1, 3, 5]);
    options.clear_ms_levels();
    assert!(!options.has_ms_levels());
    assert!(!options.has_filters());
    for level in [1, 3, 5] {
        options.add_ms_level(level).unwrap();
    }
    assert_eq!(options.ms_levels(), [1, 3, 5]);
    let mut input = vec![3, -1, 0, 3, i32::MIN, i32::MAX];
    options.set_ms_levels(&input).unwrap();
    input[0] = 42;
    assert_eq!(options.ms_levels(), [3, -1, 0, 3, i32::MIN, i32::MAX]);
    for value in [-1, 0, i32::MIN, i32::MAX] {
        assert!(options.contains_ms_level(value));
    }
    let mut copied = options.clone();
    copied.clear_ms_levels();
    assert!(!copied.contains_ms_level(3));
    assert!(options.contains_ms_level(3));
    options.set_ms_levels(&[]).unwrap();
    assert!(!options.has_ms_levels());
    assert!(!options.has_filters());
}

#[test]
fn native_ms_level_allocation_bound_is_checked_atomically() {
    let mut options = PeakFileOptions::default();
    options.add_ms_level(2).unwrap();
    let oversized = vec![1; MAX_PEAK_FILE_MS_LEVELS + 1];
    assert!(options.set_ms_levels(&oversized).is_err());
    assert_eq!(options.ms_levels(), [2]);
    options
        .set_ms_levels(&oversized[..MAX_PEAK_FILE_MS_LEVELS])
        .unwrap();
    assert!(options.add_ms_level(3).is_err());
    assert_eq!(options.ms_levels().len(), MAX_PEAK_FILE_MS_LEVELS);
    assert!(options.ms_levels().iter().all(|v| *v == 1));
    options.clear_ms_levels();
    options.add_ms_level(-5).unwrap();
    assert_eq!(options.ms_levels(), [-5]);
}

#[test]
fn source_numpress_configuration_defaults_and_exact_mode_names() {
    let default = NumpressConfig::default();
    assert_eq!(default.fixed_point, 0.0);
    assert_eq!(default.error_tolerance.to_bits(), 0.0001_f64.to_bits());
    assert_eq!(DEFAULT_NUMPRESS_ERROR_TOLERANCE, 0.0001);
    assert_eq!(default.compression, NumpressCompression::None);
    assert!(default.estimate_fixed_point);
    assert_eq!(default.linear_fp_mass_acc, -1.0);
    for (index, (mode, name)) in NumpressCompression::ALL
        .into_iter()
        .zip(["none", "linear", "pic", "slof"])
        .enumerate()
    {
        assert_eq!(usize::from(mode as u8), index);
        assert_eq!(mode.name(), name);
        assert_eq!(name.parse::<NumpressCompression>().unwrap(), mode);
        let mut config = default;
        config.set_compression(name).unwrap();
        assert_eq!(config.compression, mode);
        for invalid in [
            "", "Linear", "PIC", " slof", "none ", "linear\n", "unknown", "π",
        ] {
            assert!(config.set_compression(invalid).is_err());
            assert_eq!(config.compression, mode);
        }
    }
}

#[test]
fn numpress_dimensions_copy_all_values_and_return_literal_source_warning() {
    let mut options = PeakFileOptions::default();
    let source_warning = "Warning, compression of m/z or time dimension with pic or slof algorithms can lead to data loss";
    assert_eq!(NUMPRESS_MASS_TIME_WARNING, source_warning);
    for mode in NumpressCompression::ALL {
        let config = NumpressConfig {
            fixed_point: -123.0,
            error_tolerance: 0.0,
            compression: mode,
            estimate_fixed_point: false,
            linear_fp_mass_acc: 0.0001,
        };
        let warning = options.set_numpress_configuration_mass_time(config);
        assert_eq!(
            warning,
            match mode {
                NumpressCompression::Pic | NumpressCompression::Slof => Some(source_warning),
                _ => None,
            }
        );
        assert_eq!(options.numpress_configuration_mass_time(), config);
        assert_eq!(
            options.numpress_configuration_intensity(),
            NumpressConfig::default()
        );
        assert_eq!(
            options.numpress_configuration_float_data_array(),
            NumpressConfig::default()
        );
    }
    let unusual = NumpressConfig {
        fixed_point: f64::NEG_INFINITY,
        error_tolerance: f64::INFINITY,
        compression: NumpressCompression::Slof,
        estimate_fixed_point: false,
        linear_fp_mass_acc: f64::from_bits(0x7ff8_0000_0000_0012),
    };
    options.set_numpress_configuration_intensity(unusual);
    options.set_numpress_configuration_float_data_array(unusual);
    for copy in [
        options.numpress_configuration_intensity(),
        options.numpress_configuration_float_data_array(),
    ] {
        assert_eq!(copy.fixed_point.to_bits(), unusual.fixed_point.to_bits());
        assert_eq!(
            copy.error_tolerance.to_bits(),
            unusual.error_tolerance.to_bits()
        );
        assert_eq!(
            copy.linear_fp_mass_acc.to_bits(),
            unusual.linear_fp_mass_acc.to_bits()
        );
        assert_eq!(copy.compression, unusual.compression);
        assert!(!copy.estimate_fixed_point);
    }
    let mut retrieved = options.numpress_configuration_mass_time();
    retrieved.fixed_point = 999.0;
    assert_ne!(
        retrieved.fixed_point,
        options.numpress_configuration_mass_time().fixed_point
    );
}

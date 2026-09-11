// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Scientific peak-file option values. These options do not execute filters or
//! enable codecs by themselves; adapter resource limits are separate.

use crate::kernel::NumericRange;
use crate::{Error, Result};
use std::str::FromStr;

/// Exact default DRange<1> sentinel, using finite extrema rather than infinities.
pub const EMPTY_PEAK_FILE_RANGE: NumericRange = NumericRange {
    min: f64::MAX,
    max: f64::MIN,
};
/// Native bound on the owned MS-level vector; duplicates count separately.
pub const MAX_PEAK_FILE_MS_LEVELS: usize = 1_000_000;
/// Source default maximum relative error for Numpress configuration.
pub const DEFAULT_NUMPRESS_ERROR_TOLERANCE: f64 = 0.0001;
/// Returned by the mass/time setter for PIC and SLOF, without global output.
pub const NUMPRESS_MASS_TIME_WARNING: &str = "Warning, compression of m/z or time dimension with pic or slof algorithms can lead to data loss";

/// Numpress mode identity only. Codec execution is a separate API.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum NumpressCompression {
    #[default]
    None = 0,
    Linear = 1,
    Pic = 2,
    Slof = 3,
}
impl NumpressCompression {
    pub const ALL: [Self; 4] = [Self::None, Self::Linear, Self::Pic, Self::Slof];
    pub const fn name(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Linear => "linear",
            Self::Pic => "pic",
            Self::Slof => "slof",
        }
    }
}
impl FromStr for NumpressCompression {
    type Err = Error;
    fn from_str(text: &str) -> Result<Self> {
        match text {
            "none" => Ok(Self::None),
            "linear" => Ok(Self::Linear),
            "pic" => Ok(Self::Pic),
            "slof" => Ok(Self::Slof),
            _ => Err(Error::InvalidValue(
                "invalid Numpress compression scheme".into(),
            )),
        }
    }
}

/// Complete source Numpress configuration. Scalar values are stored verbatim,
/// including nonfinite values; a future executing codec must validate its inputs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NumpressConfig {
    pub fixed_point: f64,
    pub error_tolerance: f64,
    pub compression: NumpressCompression,
    pub estimate_fixed_point: bool,
    pub linear_fp_mass_acc: f64,
}
impl Default for NumpressConfig {
    fn default() -> Self {
        Self {
            fixed_point: 0.0,
            error_tolerance: DEFAULT_NUMPRESS_ERROR_TOLERANCE,
            compression: NumpressCompression::None,
            estimate_fixed_point: true,
            linear_fp_mass_acc: -1.0,
        }
    }
}
impl NumpressConfig {
    /// Exact case-sensitive source string mapping; an error leaves state unchanged.
    pub fn set_compression(&mut self, text: &str) -> Result<()> {
        self.compression = text.parse()?;
        Ok(())
    }
}

/// Complete PeakFileOptions state. Public scalar fields replace mechanical C++
/// getter/setter pairs. No options are automatically applied to any format API.
#[derive(Clone, Debug, PartialEq)]
pub struct PeakFileOptions {
    pub metadata_only: bool,
    /// mzXML-only source instrument/index compatibility switch.
    pub force_mq_compatibility: bool,
    /// mzML-only source isolation-window writing compatibility switch.
    pub force_tpp_compatibility: bool,
    /// Supplemental peak data in mzData.
    pub write_supplemental_data: bool,
    /// Source flag controls both spectrum m/z and chromatogram RT precision.
    pub mz_32_bit: bool,
    pub intensity_32_bit: bool,
    pub zlib_compression: bool,
    pub always_append_data: bool,
    pub skip_xml_checks: bool,
    pub sort_spectra_by_mz: bool,
    pub sort_chromatograms_by_rt: bool,
    pub fill_data: bool,
    pub write_index: bool,
    /// Processing batch size, not a resource cap. Source accepts zero.
    pub max_data_pool_size: usize,
    pub precursor_mz_selected_ion: bool,
    pub skip_chromatograms: bool,
    rt_range: NumericRange,
    mz_range: Option<NumericRange>,
    intensity_range: Option<NumericRange>,
    precursor_mz_range: Option<NumericRange>,
    ms_levels: Vec<i32>,
    np_mass_time: NumpressConfig,
    np_intensity: NumpressConfig,
    np_float_data_array: NumpressConfig,
}
impl Default for PeakFileOptions {
    fn default() -> Self {
        Self {
            metadata_only: false,
            force_mq_compatibility: false,
            force_tpp_compatibility: false,
            write_supplemental_data: true,
            mz_32_bit: false,
            intensity_32_bit: true,
            zlib_compression: false,
            always_append_data: false,
            skip_xml_checks: false,
            sort_spectra_by_mz: true,
            sort_chromatograms_by_rt: true,
            fill_data: true,
            write_index: true,
            max_data_pool_size: 100,
            precursor_mz_selected_ion: true,
            skip_chromatograms: false,
            rt_range: EMPTY_PEAK_FILE_RANGE,
            mz_range: None,
            intensity_range: None,
            precursor_mz_range: None,
            ms_levels: Vec::new(),
            np_mass_time: NumpressConfig::default(),
            np_intensity: NumpressConfig::default(),
            np_float_data_array: NumpressConfig::default(),
        }
    }
}
impl PeakFileOptions {
    pub fn new() -> Self {
        Self::default()
    }
    /// Copies the supplied source-style bounds exactly. Only the default sentinel
    /// clears RT filtering; equal endpoints, infinities and NaNs do not clear it.
    pub fn set_rt_range(&mut self, range: NumericRange) {
        self.rt_range = range;
    }
    pub fn rt_range(&self) -> NumericRange {
        self.rt_range
    }
    pub fn has_rt_range(&self) -> bool {
        self.rt_range != EMPTY_PEAK_FILE_RANGE
    }
    /// Unlike RT, setting even the empty sentinel enables this option.
    pub fn set_mz_range(&mut self, range: NumericRange) {
        self.mz_range = Some(range);
    }
    pub fn mz_range(&self) -> NumericRange {
        self.mz_range.unwrap_or(EMPTY_PEAK_FILE_RANGE)
    }
    pub fn has_mz_range(&self) -> bool {
        self.mz_range.is_some()
    }
    pub fn set_intensity_range(&mut self, range: NumericRange) {
        self.intensity_range = Some(range);
    }
    pub fn intensity_range(&self) -> NumericRange {
        self.intensity_range.unwrap_or(EMPTY_PEAK_FILE_RANGE)
    }
    pub fn has_intensity_range(&self) -> bool {
        self.intensity_range.is_some()
    }
    pub fn set_precursor_mz_range(&mut self, range: NumericRange) {
        self.precursor_mz_range = Some(range);
    }
    pub fn precursor_mz_range(&self) -> NumericRange {
        self.precursor_mz_range.unwrap_or(EMPTY_PEAK_FILE_RANGE)
    }
    pub fn has_precursor_mz_range(&self) -> bool {
        self.precursor_mz_range.is_some()
    }

    /// Copies in caller order, without sorting, deduplication or sign restrictions.
    pub fn set_ms_levels(&mut self, levels: &[i32]) -> Result<()> {
        if levels.len() > MAX_PEAK_FILE_MS_LEVELS {
            return Err(Error::InvalidValue(
                "MS-level option count limit exceeded".into(),
            ));
        }
        self.ms_levels = levels.to_vec();
        Ok(())
    }
    pub fn add_ms_level(&mut self, level: i32) -> Result<()> {
        if self.ms_levels.len() >= MAX_PEAK_FILE_MS_LEVELS {
            return Err(Error::InvalidValue(
                "MS-level option count limit exceeded".into(),
            ));
        }
        self.ms_levels.push(level);
        Ok(())
    }
    pub fn clear_ms_levels(&mut self) {
        self.ms_levels.clear();
    }
    pub fn has_ms_levels(&self) -> bool {
        !self.ms_levels.is_empty()
    }
    /// Exact membership. Empty means false, not an implicit wildcard.
    pub fn contains_ms_level(&self, level: i32) -> bool {
        self.ms_levels.contains(&level)
    }
    pub fn ms_levels(&self) -> &[i32] {
        &self.ms_levels
    }

    /// Source predicate for whole-spectrum filters. Peak m/z/intensity filters,
    /// metadata_only and skip_chromatograms are deliberately not included.
    pub fn has_filters(&self) -> bool {
        self.has_rt_range() || self.has_ms_levels() || self.has_precursor_mz_range()
    }

    /// Retains the configuration even when returning the source loss warning.
    /// Native callers choose where to route the warning instead of implicit stderr.
    pub fn set_numpress_configuration_mass_time(
        &mut self,
        config: NumpressConfig,
    ) -> Option<&'static str> {
        self.np_mass_time = config;
        matches!(
            config.compression,
            NumpressCompression::Pic | NumpressCompression::Slof
        )
        .then_some(NUMPRESS_MASS_TIME_WARNING)
    }
    pub fn numpress_configuration_mass_time(&self) -> NumpressConfig {
        self.np_mass_time
    }
    pub fn set_numpress_configuration_intensity(&mut self, config: NumpressConfig) {
        self.np_intensity = config;
    }
    pub fn numpress_configuration_intensity(&self) -> NumpressConfig {
        self.np_intensity
    }
    pub fn set_numpress_configuration_float_data_array(&mut self, config: NumpressConfig) {
        self.np_float_data_array = config;
    }
    pub fn numpress_configuration_float_data_array(&self) -> NumpressConfig {
        self.np_float_data_array
    }
}

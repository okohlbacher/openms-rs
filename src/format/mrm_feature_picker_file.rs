// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Loading component and component group picker parameters from a CSV file.
//!
//! Ports `FORMAT/MRMFeaturePickerFile.h` and its `.cpp`. The file has two
//! required columns, `component_name` and `component_group_name`. Any other
//! column is read by its header:
//!
//! - `TransitionGroupPicker:PeakPickerChromatogram:<param>` feeds the
//!   component's parameter set;
//! - `TransitionGroupPicker:<param>` (without the `PeakPickerChromatogram:`
//!   prefix) feeds the component group's parameter set;
//! - every other column is ignored.
//!
//! Cells that are empty are not written; rows with an empty `component_name`
//! or `component_group_name` are skipped. A reduced example:
//!
//! ```text
//! component_name,component_group_name,TransitionGroupPicker:stop_after_feature,TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length
//! arg-L.arg-L_1.Heavy,arg-L,2,15
//! arg-L.arg-L_1.Light,arg-L,2,17
//! orn.orn_1.Heavy,orn,3,21
//! ```
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/MRM_FEATURE_QC_SUPPORT.md).

use crate::analysis::openswath::mrm_feature_picker::{ComponentGroupParams, ComponentParams};
use crate::concept::boost_regex::BoostRegex;
use crate::format::csv::{CsvFile, ReadOptions};
use crate::format::mrm_feature_qc_file::{string_utils_to_double, string_utils_to_int32};
use crate::param::{Param, ParamValue};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::path::Path;

const GROUP_PATTERN: &str = "TransitionGroupPicker:(?!PeakPickerChromatogram:)(.+)";
const COMPONENT_PATTERN: &str = "TransitionGroupPicker:PeakPickerChromatogram:(.+)";

const DOUBLE_HEADERS: [&str; 9] = [
    "gauss_width",
    "peak_width",
    "signal_to_noise",
    "sn_win_len",
    "stop_after_intensity_ratio",
    "min_peak_width",
    "recalculate_peaks_max_z",
    "minimal_quality",
    "resample_boundary",
];
const BOOL_HEADERS: [&str; 7] = [
    "use_gauss",
    "write_sn_log_messages",
    "remove_overlapping_peaks",
    "recalculate_peaks",
    "use_precursors",
    "compute_peak_quality",
    "compute_peak_shape_metrics",
];
const UINT_HEADERS: [&str; 3] = [
    "sgolay_frame_length",
    "sgolay_polynomial_order",
    "sn_bin_count",
];
const INT_HEADERS: [&str; 1] = ["stop_after_feature"];

/// Reads [`ComponentParams`] and [`ComponentGroupParams`] from a CSV file.
///
/// Source `MRMFeaturePickerFile`, which derives from `CsvFile`; the CSV buffer
/// is local to [`load`](Self::load) here, so the type holds no state.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MRMFeaturePickerFile;

impl MRMFeaturePickerFile {
    /// Read `filename`, returning one [`ComponentParams`] per row with
    /// non-empty names and one [`ComponentGroupParams`] the first time each
    /// `component_group_name` is seen.
    ///
    /// The file is parsed as comma separated and its first line is the header.
    /// Later rows of an already seen group extend the component list but add
    /// no group entry, so a group's parameters are those of its first row.
    /// With zero or one line no header is validated and both lists are empty.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when `filename` cannot be opened;
    /// [`Error::MissingInformation`] when the `component_name` or
    /// `component_group_name` column is absent;
    /// [`Error::InvalidValue`] when a numeric cell of a typed parameter cannot
    /// be converted, or a row is shorter than a header's column index (the
    /// source reads past the end of the row).
    pub fn load(
        &self,
        filename: impl AsRef<Path>,
    ) -> Result<(Vec<ComponentParams>, Vec<ComponentGroupParams>)> {
        let mut csv = CsvFile::new();
        csv.load(filename, &ReadOptions::default())?;
        let mut headers: BTreeMap<String, usize> = BTreeMap::new();
        if csv.row_count() >= 2 {
            for (i, name) in csv.row(0)?.1.into_iter().enumerate() {
                headers.insert(name, i);
            }
            if !(headers.contains_key("component_name")
                && headers.contains_key("component_group_name"))
            {
                return Err(Error::MissingInformation(
                    "Columns component_name and/or component_group_name not found.".into(),
                ));
            }
        }
        let mut cp_list = Vec::new();
        let mut cgp_list: Vec<ComponentGroupParams> = Vec::new();
        for i in 1..csv.row_count() {
            let line = csv.row(i)?.1;
            if let Some((cp, cgp)) = self.extract_params_from_line(&line, &headers)? {
                cp_list.push(cp);
                if !cgp_list
                    .iter()
                    .any(|c| c.component_group_name == cgp.component_group_name)
                {
                    cgp_list.push(cgp);
                }
            }
        }
        Ok((cp_list, cgp_list))
    }

    /// Extract the parameters of one parsed line.
    ///
    /// Returns `None` when `component_name` or `component_group_name` is empty.
    /// Every column whose header contains `TransitionGroupPicker:` feeds the
    /// group parameters, or the component parameters when the prefix is
    /// `TransitionGroupPicker:PeakPickerChromatogram:`; see
    /// [`set_cast_value`](Self::set_cast_value) for the value types.
    ///
    /// # Errors
    ///
    /// As [`load`](Self::load); also [`Error::MissingInformation`] when
    /// `headers` has no name column (the source dereferences `end()`).
    pub fn extract_params_from_line(
        &self,
        line: &[String],
        headers: &BTreeMap<String, usize>,
    ) -> Result<Option<(ComponentParams, ComponentGroupParams)>> {
        let cell = |header: &str, index: usize| -> Result<&str> {
            line.get(index).map(String::as_str).ok_or_else(|| {
                Error::InvalidValue(format!(
                    "row has {} cell(s) but column '{header}' is number {}; the source reads past the end of the row",
                    line.len(),
                    index + 1
                ))
            })
        };
        let column = |name: &str| -> Result<usize> {
            headers
                .get(name)
                .copied()
                .ok_or_else(|| Error::MissingInformation(format!("Column {name} not found.")))
        };
        let component_name = cell("component_name", column("component_name")?)?.to_owned();
        let component_group_name =
            cell("component_group_name", column("component_group_name")?)?.to_owned();
        if component_name.is_empty() || component_group_name.is_empty() {
            return Ok(None);
        }
        let mut cp = ComponentParams {
            component_name,
            component_group_name: component_group_name.clone(),
            params: Param::new(),
        };
        let mut cgp = ComponentGroupParams {
            component_group_name,
            params: Param::new(),
        };
        let group = BoostRegex::new(GROUP_PATTERN)?;
        let component = BoostRegex::new(COMPONENT_PATTERN)?;
        for (header, &index) in headers {
            if let Some(found) = group.search(header.as_bytes())? {
                if let Some(key) = found.get(1) {
                    self.set_cast_value(&header[key], cell(header, index)?, &mut cgp.params)?;
                }
            } else if let Some(found) = component.search(header.as_bytes())? {
                if let Some(key) = found.get(1) {
                    self.set_cast_value(&header[key], cell(header, index)?, &mut cp.params)?;
                }
            }
        }
        Ok(Some((cp, cgp)))
    }

    /// Insert `value` into `params` under `key`, converted to the type the
    /// parameter name has in the source's fixed table.
    ///
    /// An empty `value` is dropped. The types are:
    ///
    /// - `gauss_width`, `peak_width`, `signal_to_noise`, `sn_win_len`,
    ///   `stop_after_intensity_ratio`, `min_peak_width`,
    ///   `recalculate_peaks_max_z`, `minimal_quality`, `resample_boundary`:
    ///   floating point;
    /// - `use_gauss`, `write_sn_log_messages`, `remove_overlapping_peaks`,
    ///   `recalculate_peaks`, `use_precursors`, `compute_peak_quality`,
    ///   `compute_peak_shape_metrics`: the string `true` when the value is
    ///   exactly `true` or `TRUE` and `false` for anything else (the source
    ///   header calls this case-insensitive; the code is not);
    /// - `sgolay_frame_length`, `sgolay_polynomial_order`, `sn_bin_count`: an
    ///   unsigned integer from the value read as a double and truncated;
    /// - `stop_after_feature`: a 32-bit integer;
    /// - anything else: the string unchanged.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a numeric value cannot be converted
    /// (`StringUtils::toDouble`/`toInt32` throw `ConversionError`).
    ///
    /// The unsigned truncation of a negative or out-of-range double is
    /// undefined in the source; this port saturates to `0` and `u32::MAX`.
    pub fn set_cast_value(&self, key: &str, value: &str, params: &mut Param) -> Result<()> {
        if value.is_empty() {
            return Ok(());
        }
        let converted = if DOUBLE_HEADERS.contains(&key) {
            ParamValue::Float(string_utils_to_double(value)?)
        } else if BOOL_HEADERS.contains(&key) {
            ParamValue::String(if value == "true" || value == "TRUE" {
                "true".into()
            } else {
                "false".into()
            })
        } else if UINT_HEADERS.contains(&key) {
            ParamValue::Integer(i64::from(string_utils_to_double(value)? as u32))
        } else if INT_HEADERS.contains(&key) {
            ParamValue::Integer(i64::from(string_utils_to_int32(value)?))
        } else {
            ParamValue::String(value.to_owned())
        };
        params.set_value(key, converted, "", &[])
    }
}

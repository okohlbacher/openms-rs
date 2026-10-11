// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Loading and storing [`MRMFeatureQC`] criteria as a comma-separated file.
//!
//! Ports `FORMAT/MRMFeatureQCFile.h` and its `.cpp`. A component file has the
//! header `component_name, retention_time_l, retention_time_u, intensity_l,
//! intensity_u, overall_quality_l, overall_quality_u`; a component group file
//! starts with `component_group_name` and adds the `n_*_l`/`n_*_u` counts, the
//! ion ratio columns and the same bounds. Any further column named
//! `metaValue_<name>_l` or `metaValue_<name>_u` sets a bound of the metavalue
//! `<name>`. Cells left empty take the defaults of [`MRMFeatureQC`]'s members;
//! rows with an empty component or component group name are skipped.
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/MRM_FEATURE_QC_SUPPORT.md).
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::mrm_feature_qc::MRMFeatureQC;
//! use openms::format::mrm_feature_qc_file::MRMFeatureQCFile;
//!
//! let dir = std::env::temp_dir().join("openms_mrm_feature_qc_file_doctest");
//! std::fs::create_dir_all(&dir).unwrap();
//! let path = dir.join("components.csv");
//! std::fs::write(
//!     &path,
//!     "component_name,retention_time_l,retention_time_u,metaValue_sn_score_l\nglucose,1.5,2.5,10\n",
//! )
//! .unwrap();
//! let mut qc = MRMFeatureQC::default();
//! MRMFeatureQCFile::default().load(&path, &mut qc, false).unwrap();
//! assert_eq!(qc.component_qcs[0].component_name, "glucose");
//! assert_eq!(qc.component_qcs[0].retention_time_l, 1.5);
//! assert_eq!(qc.component_qcs[0].meta_value_qc["sn_score"], (10.0, 1e12));
//! ```
//!
//! [`MRMFeatureQC`]: crate::analysis::openswath::mrm_feature_qc::MRMFeatureQC

use crate::analysis::openswath::mrm_feature_qc::{
    ComponentGroupQCs, ComponentQCs, DEFAULT_COUNT_UPPER_BOUND, DEFAULT_UPPER_BOUND, MRMFeatureQC,
};
use crate::concept::boost_regex::BoostRegex;
use crate::format::csv::{CsvFile, ReadOptions};
use crate::format::file_info::text_format::to_str;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::path::Path;

/// The source's `std::map<std::string, Size>` of header name to column index.
pub type HeaderMap = BTreeMap<String, usize>;

const META_VALUE_PATTERN: &str = "metaValue_(.+)_(l|u)";

/// C `isspace` in the "C" locale, as `strtod` and `strtol` skip it.
fn is_c_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

fn conversion(function: &str, text: &str) -> Error {
    Error::InvalidValue(format!("{function}: no conversion of '{text}'"))
}

/// `std::stod`: the longest `strtod` prefix, which may be followed by anything.
///
/// Leading C whitespace is skipped. `inf`, `infinity` and `nan` are accepted
/// in any case. A decimal literal that overflows to infinity or underflows to
/// zero is the source's `std::out_of_range`, and a string with no number is its
/// `std::invalid_argument`; both are [`Error::InvalidValue`] here. Hexadecimal
/// floating literals, which `strtod` also reads, are refused as
/// [`Error::Unsupported`].
pub(crate) fn std_stod(text: &str) -> Result<f64> {
    let bytes = text.as_bytes();
    let mut i = 0;
    while bytes.get(i).is_some_and(|&b| is_c_space(b)) {
        i += 1;
    }
    let start = i;
    if matches!(bytes.get(i), Some(b'+' | b'-')) {
        i += 1;
    }
    let ci = |at: usize, word: &str| {
        bytes
            .get(at..at + word.len())
            .is_some_and(|w| w.eq_ignore_ascii_case(word.as_bytes()))
    };
    if ci(i, "infinity") || ci(i, "inf") {
        let negative = bytes[start] == b'-';
        return Ok(if negative {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        });
    }
    if ci(i, "nan") {
        return Ok(f64::NAN);
    }
    if bytes.get(i) == Some(&b'0') && matches!(bytes.get(i + 1), Some(b'x' | b'X')) {
        let digit_at = |at: usize| bytes.get(at).is_some_and(u8::is_ascii_hexdigit);
        if digit_at(i + 2) || (bytes.get(i + 2) == Some(&b'.') && digit_at(i + 3)) {
            return Err(Error::Unsupported(
                "std::stod: hexadecimal floating literals are not read".into(),
            ));
        }
    }
    let int_start = i;
    while bytes.get(i).is_some_and(u8::is_ascii_digit) {
        i += 1;
    }
    let mut digits = i - int_start;
    if bytes.get(i) == Some(&b'.') {
        let frac_start = i + 1;
        let mut j = frac_start;
        while bytes.get(j).is_some_and(u8::is_ascii_digit) {
            j += 1;
        }
        if digits > 0 || j > frac_start {
            digits += j - frac_start;
            i = j;
        }
    }
    if digits == 0 {
        return Err(conversion("stod", text));
    }
    if matches!(bytes.get(i), Some(b'e' | b'E')) {
        let mut j = i + 1;
        if matches!(bytes.get(j), Some(b'+' | b'-')) {
            j += 1;
        }
        let exponent_start = j;
        while bytes.get(j).is_some_and(u8::is_ascii_digit) {
            j += 1;
        }
        if j > exponent_start {
            i = j;
        }
    }
    let literal = &text[start..i];
    let value: f64 = literal.parse().map_err(|_| conversion("stod", text))?;
    let nonzero_mantissa = literal
        .bytes()
        .take_while(|b| !matches!(b, b'e' | b'E'))
        .any(|b| (b'1'..=b'9').contains(&b));
    if value.is_infinite() || (value == 0.0 && nonzero_mantissa) {
        return Err(Error::InvalidValue(format!(
            "stod: '{text}' is out of the range of double"
        )));
    }
    Ok(value)
}

/// `std::stoi`: the longest `strtol` prefix, which may be followed by anything.
///
/// Leading C whitespace and one sign are skipped. A string with no digits is
/// the source's `std::invalid_argument` and a value outside `i32` its
/// `std::out_of_range`; both are [`Error::InvalidValue`].
pub(crate) fn std_stoi(text: &str) -> Result<i32> {
    let bytes = text.as_bytes();
    let mut i = 0;
    while bytes.get(i).is_some_and(|&b| is_c_space(b)) {
        i += 1;
    }
    let negative = bytes.get(i) == Some(&b'-');
    if matches!(bytes.get(i), Some(b'+' | b'-')) {
        i += 1;
    }
    let digits_start = i;
    let mut magnitude: i128 = 0;
    while let Some(&b) = bytes.get(i).filter(|b| b.is_ascii_digit()) {
        magnitude = (magnitude * 10 + i128::from(b - b'0')).min(i128::from(i64::MAX) + 1);
        i += 1;
    }
    if i == digits_start {
        return Err(conversion("stoi", text));
    }
    let value = if negative { -magnitude } else { magnitude };
    i32::try_from(value)
        .map_err(|_| Error::InvalidValue(format!("stoi: '{text}' is out of the range of int")))
}

/// `StringUtils::toDouble`, the strict conversion `MRMFeaturePickerFile` uses:
/// surrounding whitespace is allowed and the rest must be one number.
///
/// One leading `+` is skipped, `nan` (with an optional parenthesised payload)
/// and `inf`/`infinity` are read in any case, a finite literal that overflows
/// is an error, and hexadecimal floats are not read.
pub(crate) fn string_utils_to_double(text: &str) -> Result<f64> {
    let failure = || {
        Error::InvalidValue(format!(
            "Could not convert string '{text}' to a double value"
        ))
    };
    let trimmed = text.trim_matches(|c: char| c.is_ascii() && is_c_space(c as u8));
    if trimmed.is_empty() {
        return Err(failure());
    }
    let lower = trimmed.to_ascii_lowercase();
    if lower == "nan" || (lower.starts_with("nan(") && lower.ends_with(')')) {
        return Ok(f64::NAN);
    }
    let body = trimmed.strip_prefix('+').unwrap_or(trimmed);
    let is_infinity_word = {
        let word = body.strip_prefix('-').unwrap_or(body).to_ascii_lowercase();
        word == "inf" || word == "infinity"
    };
    if body.starts_with(|c: char| c.is_ascii_alphabetic()) && !is_infinity_word {
        return Err(failure());
    }
    let value: f64 = body.parse().map_err(|_| failure())?;
    if value.is_infinite() && !is_infinity_word {
        return Err(failure());
    }
    Ok(value)
}

/// `StringUtils::toInt32`, the strict conversion `MRMFeaturePickerFile` uses.
pub(crate) fn string_utils_to_int32(text: &str) -> Result<i32> {
    let trimmed = text.trim_matches(|c: char| c.is_ascii() && is_c_space(c as u8));
    let body = trimmed.strip_prefix('+').unwrap_or(trimmed);
    body.parse::<i32>().map_err(|_| {
        Error::InvalidValue(format!(
            "Could not convert string '{text}' to an integer value"
        ))
    })
}

/// Loads and stores [`MRMFeatureQC`] files.
///
/// Source `MRMFeatureQCFile` privately inherits `CsvFile` and publicly
/// `ProgressLogger`; neither is used by `load` or `store`, so this type holds
/// no state. The methods the source declares `protected` are public here
/// because the class test reaches them through a facade.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MRMFeatureQCFile;

fn cell<'a>(line: &'a [String], index: usize, header: &str) -> Result<&'a str> {
    line.get(index).map(String::as_str).ok_or_else(|| {
        Error::InvalidValue(format!(
            "row has {} cell(s) but column '{header}' is number {}; the source reads past the end of the row",
            line.len(),
            index + 1
        ))
    })
}

impl MRMFeatureQCFile {
    /// Load a QC file, replacing the component QCs (`is_component_group`
    /// false) or the component group QCs (true) of `qc`.
    ///
    /// The first row names the columns; every later row is one QC. The other
    /// list of `qc` is left untouched.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when the file cannot be opened and
    /// [`Error::InvalidValue`] when a number cannot be converted
    /// (`std::stod`/`std::stoi` throw in the source) or a row has fewer cells
    /// than the header has columns (the source reads out of bounds). Unlike the
    /// source, which clears the list before parsing, `qc` is left unchanged on
    /// error.
    pub fn load(
        &self,
        filename: impl AsRef<Path>,
        qc: &mut MRMFeatureQC,
        is_component_group: bool,
    ) -> Result<()> {
        let csv = CsvFile::from_path(filename, &ReadOptions::default())?;
        let mut headers = HeaderMap::new();
        let mut first_row: Vec<String> = Vec::new();
        if csv.row_count() > 0 {
            first_row = csv.row(0)?.1;
        }
        for (i, name) in first_row.iter().enumerate() {
            headers.insert(name.clone(), i);
        }
        if is_component_group {
            let mut loaded = Vec::new();
            for i in 1..csv.row_count() {
                let line = csv.row(i)?.1;
                self.push_values_from_line_component_group(&line, &headers, &mut loaded)?;
            }
            qc.component_group_qcs = loaded;
        } else {
            let mut loaded = Vec::new();
            for i in 1..csv.row_count() {
                let line = csv.row(i)?.1;
                self.push_values_from_line_component(&line, &headers, &mut loaded)?;
            }
            qc.component_qcs = loaded;
        }
        Ok(())
    }

    /// Store the component QCs (`is_component_group` false) or the component
    /// group QCs (true) of `qc` as a comma-separated file.
    ///
    /// The metavalue columns are those of the first QC, in key order. When the
    /// list is empty no file is written, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when the file cannot be created and
    /// [`Error::InvalidValue`] when a QC's metavalue names differ from the first
    /// QC's. The source writes the values positionally, which misaligns the
    /// columns when the names differ and writes past the end of the row when a
    /// QC has more metavalues.
    pub fn store(
        &self,
        filename: impl AsRef<Path>,
        qc: &MRMFeatureQC,
        is_component_group: bool,
    ) -> Result<()> {
        let mut csv = CsvFile::new();
        if is_component_group {
            let Some(first) = qc.component_group_qcs.first() else {
                return Ok(());
            };
            let mut headers: Vec<String> = [
                "component_group_name",
                "retention_time_l",
                "retention_time_u",
                "intensity_l",
                "intensity_u",
                "overall_quality_l",
                "overall_quality_u",
                "n_heavy_l",
                "n_heavy_u",
                "n_light_l",
                "n_light_u",
                "n_detecting_l",
                "n_detecting_u",
                "n_quantifying_l",
                "n_quantifying_u",
                "n_identifying_l",
                "n_identifying_u",
                "n_transitions_l",
                "n_transitions_u",
                "ion_ratio_pair_name_1",
                "ion_ratio_pair_name_2",
                "ion_ratio_l",
                "ion_ratio_u",
                "ion_ratio_feature_name",
            ]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
            for name in first.meta_value_qc.keys() {
                headers.push(format!("metaValue_{name}_l"));
                headers.push(format!("metaValue_{name}_u"));
            }
            csv.add_row(&headers)?;
            for q in &qc.component_group_qcs {
                let mut row = vec![
                    q.component_group_name.clone(),
                    to_str(q.retention_time_l),
                    to_str(q.retention_time_u),
                    to_str(q.intensity_l),
                    to_str(q.intensity_u),
                    to_str(q.overall_quality_l),
                    to_str(q.overall_quality_u),
                    q.n_heavy_l.to_string(),
                    q.n_heavy_u.to_string(),
                    q.n_light_l.to_string(),
                    q.n_light_u.to_string(),
                    q.n_detecting_l.to_string(),
                    q.n_detecting_u.to_string(),
                    q.n_quantifying_l.to_string(),
                    q.n_quantifying_u.to_string(),
                    q.n_identifying_l.to_string(),
                    q.n_identifying_u.to_string(),
                    q.n_transitions_l.to_string(),
                    q.n_transitions_u.to_string(),
                    q.ion_ratio_pair_name_1.clone(),
                    q.ion_ratio_pair_name_2.clone(),
                    to_str(q.ion_ratio_l),
                    to_str(q.ion_ratio_u),
                    q.ion_ratio_feature_name.clone(),
                ];
                append_meta_cells(&mut row, &q.meta_value_qc, &first.meta_value_qc)?;
                csv.add_row(&row)?;
            }
        } else {
            let Some(first) = qc.component_qcs.first() else {
                return Ok(());
            };
            let mut headers: Vec<String> = [
                "component_name",
                "retention_time_l",
                "retention_time_u",
                "intensity_l",
                "intensity_u",
                "overall_quality_l",
                "overall_quality_u",
            ]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
            for name in first.meta_value_qc.keys() {
                headers.push(format!("metaValue_{name}_l"));
                headers.push(format!("metaValue_{name}_u"));
            }
            csv.add_row(&headers)?;
            for q in &qc.component_qcs {
                let mut row = vec![
                    q.component_name.clone(),
                    to_str(q.retention_time_l),
                    to_str(q.retention_time_u),
                    to_str(q.intensity_l),
                    to_str(q.intensity_u),
                    to_str(q.overall_quality_l),
                    to_str(q.overall_quality_u),
                ];
                append_meta_cells(&mut row, &q.meta_value_qc, &first.meta_value_qc)?;
                csv.add_row(&row)?;
            }
        }
        csv.store(filename)
    }

    /// Save the values of one row to a new [`ComponentQCs`] appended to
    /// `c_qcs`.
    ///
    /// A row whose `component_name` is missing or empty is skipped. Every
    /// header of the form `metaValue_<name>_l` or `_u` (matched, as in the
    /// source, anywhere in the header with `metaValue_(.+)_(l|u)`) sets a
    /// metavalue bound through [`set_pair_value`](Self::set_pair_value); an
    /// empty cell there takes the default of its bound.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] on an unconvertible number or when `line` is
    /// shorter than a header's column index.
    pub fn push_values_from_line_component(
        &self,
        line: &[String],
        headers: &HeaderMap,
        c_qcs: &mut Vec<ComponentQCs>,
    ) -> Result<()> {
        let name = self.get_cast_value_string(headers, line, "component_name", "")?;
        if name.is_empty() {
            return Ok(());
        }
        let mut c = ComponentQCs {
            component_name: name,
            ..ComponentQCs::default()
        };
        c.retention_time_l = self.get_cast_value_f64(headers, line, "retention_time_l", 0.0)?;
        c.retention_time_u =
            self.get_cast_value_f64(headers, line, "retention_time_u", DEFAULT_UPPER_BOUND)?;
        c.intensity_l = self.get_cast_value_f64(headers, line, "intensity_l", 0.0)?;
        c.intensity_u =
            self.get_cast_value_f64(headers, line, "intensity_u", DEFAULT_UPPER_BOUND)?;
        c.overall_quality_l = self.get_cast_value_f64(headers, line, "overall_quality_l", 0.0)?;
        c.overall_quality_u =
            self.get_cast_value_f64(headers, line, "overall_quality_u", DEFAULT_UPPER_BOUND)?;
        self.read_meta_values(line, headers, &mut c.meta_value_qc)?;
        c_qcs.push(c);
        Ok(())
    }

    /// Save the values of one row to a new [`ComponentGroupQCs`] appended to
    /// `cg_qcs`; as [`push_values_from_line_component`](Self::push_values_from_line_component)
    /// with `component_group_name` as the required column.
    ///
    /// # Errors
    ///
    /// As [`push_values_from_line_component`](Self::push_values_from_line_component).
    pub fn push_values_from_line_component_group(
        &self,
        line: &[String],
        headers: &HeaderMap,
        cg_qcs: &mut Vec<ComponentGroupQCs>,
    ) -> Result<()> {
        let name = self.get_cast_value_string(headers, line, "component_group_name", "")?;
        if name.is_empty() {
            return Ok(());
        }
        let mut cg = ComponentGroupQCs {
            component_group_name: name,
            ..ComponentGroupQCs::default()
        };
        let upper = DEFAULT_UPPER_BOUND;
        let n_upper = DEFAULT_COUNT_UPPER_BOUND;
        cg.retention_time_l = self.get_cast_value_f64(headers, line, "retention_time_l", 0.0)?;
        cg.retention_time_u = self.get_cast_value_f64(headers, line, "retention_time_u", upper)?;
        cg.intensity_l = self.get_cast_value_f64(headers, line, "intensity_l", 0.0)?;
        cg.intensity_u = self.get_cast_value_f64(headers, line, "intensity_u", upper)?;
        cg.overall_quality_l = self.get_cast_value_f64(headers, line, "overall_quality_l", 0.0)?;
        cg.overall_quality_u =
            self.get_cast_value_f64(headers, line, "overall_quality_u", upper)?;
        cg.n_heavy_l = self.get_cast_value_i32(headers, line, "n_heavy_l", 0)?;
        cg.n_heavy_u = self.get_cast_value_i32(headers, line, "n_heavy_u", n_upper)?;
        cg.n_light_l = self.get_cast_value_i32(headers, line, "n_light_l", 0)?;
        cg.n_light_u = self.get_cast_value_i32(headers, line, "n_light_u", n_upper)?;
        cg.n_detecting_l = self.get_cast_value_i32(headers, line, "n_detecting_l", 0)?;
        cg.n_detecting_u = self.get_cast_value_i32(headers, line, "n_detecting_u", n_upper)?;
        cg.n_quantifying_l = self.get_cast_value_i32(headers, line, "n_quantifying_l", 0)?;
        cg.n_quantifying_u = self.get_cast_value_i32(headers, line, "n_quantifying_u", n_upper)?;
        cg.n_identifying_l = self.get_cast_value_i32(headers, line, "n_identifying_l", 0)?;
        cg.n_identifying_u = self.get_cast_value_i32(headers, line, "n_identifying_u", n_upper)?;
        cg.n_transitions_l = self.get_cast_value_i32(headers, line, "n_transitions_l", 0)?;
        cg.n_transitions_u = self.get_cast_value_i32(headers, line, "n_transitions_u", n_upper)?;
        cg.ion_ratio_pair_name_1 =
            self.get_cast_value_string(headers, line, "ion_ratio_pair_name_1", "")?;
        cg.ion_ratio_pair_name_2 =
            self.get_cast_value_string(headers, line, "ion_ratio_pair_name_2", "")?;
        cg.ion_ratio_l = self.get_cast_value_f64(headers, line, "ion_ratio_l", 0.0)?;
        cg.ion_ratio_u = self.get_cast_value_f64(headers, line, "ion_ratio_u", upper)?;
        cg.ion_ratio_feature_name =
            self.get_cast_value_string(headers, line, "ion_ratio_feature_name", "")?;
        self.read_meta_values(line, headers, &mut cg.meta_value_qc)?;
        cg_qcs.push(cg);
        Ok(())
    }

    fn read_meta_values(
        &self,
        line: &[String],
        headers: &HeaderMap,
        meta_values_qc: &mut BTreeMap<String, (f64, f64)>,
    ) -> Result<()> {
        let regex = BoostRegex::new(META_VALUE_PATTERN)?;
        for (header, &index) in headers {
            let Some(found) = regex.search(header.as_bytes())? else {
                continue;
            };
            let (Some(name), Some(boundary)) = (found.get(1), found.get(2)) else {
                continue;
            };
            let value = cell(line, index, header)?;
            self.set_pair_value(&header[name], value, &header[boundary], meta_values_qc)?;
        }
        Ok(())
    }

    /// Set one of the two bounds of a metavalue.
    ///
    /// When `key` is present its pair is updated at the `boundary` position
    /// (`"l"` lower, anything else upper). When it is absent a new pair is
    /// created whose other element is the default: `1e12` for a missing upper
    /// bound, `0.0` for a missing lower bound. An empty `value` is the default
    /// of its own bound (`0.0` lower, `1e12` upper); otherwise it is read with
    /// `std::stod` semantics.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `value` has no leading number.
    pub fn set_pair_value(
        &self,
        key: &str,
        value: &str,
        boundary: &str,
        meta_values_qc: &mut BTreeMap<String, (f64, f64)>,
    ) -> Result<()> {
        let lower = boundary == "l";
        let cast = if value.is_empty() {
            if lower { 0.0 } else { DEFAULT_UPPER_BOUND }
        } else {
            std_stod(value)?
        };
        match meta_values_qc.get_mut(key) {
            Some(pair) => {
                if lower {
                    pair.0 = cast;
                } else {
                    pair.1 = cast;
                }
            }
            None => {
                meta_values_qc.insert(
                    key.to_owned(),
                    if lower {
                        (cast, DEFAULT_UPPER_BOUND)
                    } else {
                        (0.0, cast)
                    },
                );
            }
        }
        Ok(())
    }

    /// Read a column of a row as an `Int` (source `getCastValue_`, `Int`
    /// overload): the cell converted with `std::stoi` when the column exists
    /// and the cell is not empty, otherwise `default_value`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] on an unconvertible number or when the column's
    /// index lies past the end of `line`.
    pub fn get_cast_value_i32(
        &self,
        headers: &HeaderMap,
        line: &[String],
        header: &str,
        default_value: i32,
    ) -> Result<i32> {
        match headers.get(header) {
            Some(&index) => {
                let text = cell(line, index, header)?;
                if text.is_empty() {
                    Ok(default_value)
                } else {
                    std_stoi(text)
                }
            }
            None => Ok(default_value),
        }
    }

    /// Read a column of a row as a `double` (source `getCastValue_`, `double`
    /// overload), converted with `std::stod` semantics.
    ///
    /// # Errors
    ///
    /// As [`get_cast_value_i32`](Self::get_cast_value_i32).
    pub fn get_cast_value_f64(
        &self,
        headers: &HeaderMap,
        line: &[String],
        header: &str,
        default_value: f64,
    ) -> Result<f64> {
        match headers.get(header) {
            Some(&index) => {
                let text = cell(line, index, header)?;
                if text.is_empty() {
                    Ok(default_value)
                } else {
                    std_stod(text)
                }
            }
            None => Ok(default_value),
        }
    }

    /// Read a column of a row as a string (source `getCastValue_`, `String`
    /// overload).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the column's index lies past the end of
    /// `line`.
    pub fn get_cast_value_string(
        &self,
        headers: &HeaderMap,
        line: &[String],
        header: &str,
        default_value: &str,
    ) -> Result<String> {
        match headers.get(header) {
            Some(&index) => {
                let text = cell(line, index, header)?;
                Ok(if text.is_empty() {
                    default_value.to_owned()
                } else {
                    text.to_owned()
                })
            }
            None => Ok(default_value.to_owned()),
        }
    }
}

fn append_meta_cells(
    row: &mut Vec<String>,
    meta: &BTreeMap<String, (f64, f64)>,
    template: &BTreeMap<String, (f64, f64)>,
) -> Result<()> {
    if !meta.keys().eq(template.keys()) {
        return Err(Error::InvalidValue(
            "every QC of a file must carry the same metavalue names as the first".into(),
        ));
    }
    for (lower, upper) in meta.values() {
        row.push(to_str(*lower));
        row.push(to_str(*upper));
    }
    Ok(())
}

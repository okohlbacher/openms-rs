// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The `DefaultParamHandler` contract of `SpectrumAnnotator`: the eight
//! parameters its constructor declares (`SpectrumAnnotator.cpp:24-78`) and the
//! members its `updateMembers_` reads from them (`:475-485`).

use super::SpectrumAnnotator;
use crate::param::handler::{set_flag, unsigned_parameter};
use crate::param::{DefaultParamHandling, Param, ParamValue};
use crate::{Error, Result};

/// The seven boolean parameters, in declaration order around
/// `topNmatch_fragmenterrors` (declared sixth): key and description.
const BEFORE_TOP_N: [(&str, &str); 5] = [
    (
        "basic_statistics",
        "If set, meta values for peak_number, sum_intensity, matched_ion_number, matched_intensity are added",
    ),
    (
        "list_of_ions_matched",
        "If set, meta values for matched_ions are added",
    ),
    (
        "max_series",
        "If set, meta values for max_series_type, max_series_size are added",
    ),
    // The source's description, copied from an older isotope option, says
    // nothing about the signal-to-noise meta values the parameter controls.
    (
        "S/N_statistics",
        "If set to 1 isotope peaks of the product ion peaks are added",
    ),
    (
        "precursor_statistics",
        "If set, meta values for precursor_in_ms2 are added",
    ),
];

const TOP_N: &str = "topNmatch_fragmenterrors";

const AFTER_TOP_N: [(&str, &str); 2] = [
    (
        "fragmenterror_statistics",
        "If set, meta values for median_fragment_error, IQR_fragment_error are added",
    ),
    (
        "terminal_series_match_ratio",
        "If set, meta values for NTermIonCurrentRatio, CTermIonCurrentRatio are added",
    ),
];

impl SpectrumAnnotator {
    /// The constructor's declarations carrying this annotator's values.
    fn parameter_tree(&self) -> Result<Param> {
        let before = [
            self.basic_statistics,
            self.list_of_ions_matched,
            self.max_series,
            self.sn_statistics,
            self.precursor_statistics,
        ];
        let after = [
            self.fragment_error_statistics,
            self.terminal_series_match_ratio,
        ];
        let top_n = u32::try_from(self.top_n_fragment_errors).map_err(|_| {
            Error::InvalidValue(format!(
                "SpectrumAnnotator: '{TOP_N}' = {} does not fit the source's unsigned int",
                self.top_n_fragment_errors
            ))
        })?;
        let mut param = Param::new();
        for ((key, description), value) in BEFORE_TOP_N.into_iter().zip(before) {
            set_flag(&mut param, key, value, description, &[])?;
        }
        // `unsigned(7)`: an unsigned value is an integer parameter with no
        // restriction; `updateMembers_` refuses a negative one.
        param.set_value(
            TOP_N,
            ParamValue::Integer(i64::from(top_n)),
            "If set n > 0,  meta values for topN_meanfragmenterror, topN_MSEfragmenterror, topN_stddevfragmenterror are added",
            &[],
        )?;
        for ((key, description), value) in AFTER_TOP_N.into_iter().zip(after) {
            set_flag(&mut param, key, value, description, &[])?;
        }
        Ok(param)
    }
}

/// Source `SpectrumAnnotator`'s `DefaultParamHandler` base.
///
/// | Parameter | Default | Field |
/// |---|---|---|
/// | `basic_statistics` | `true` | [`basic_statistics`](SpectrumAnnotator::basic_statistics) |
/// | `list_of_ions_matched` | `true` | [`list_of_ions_matched`](SpectrumAnnotator::list_of_ions_matched) |
/// | `max_series` | `true` | [`max_series`](SpectrumAnnotator::max_series) |
/// | `S/N_statistics` | `true` | [`sn_statistics`](SpectrumAnnotator::sn_statistics) |
/// | `precursor_statistics` | `true` | [`precursor_statistics`](SpectrumAnnotator::precursor_statistics) |
/// | `topNmatch_fragmenterrors` (int) | 7 | [`top_n_fragment_errors`](SpectrumAnnotator::top_n_fragment_errors) |
/// | `fragmenterror_statistics` | `true` | [`fragment_error_statistics`](SpectrumAnnotator::fragment_error_statistics) |
/// | `terminal_series_match_ratio` | `true` | [`terminal_series_match_ratio`](SpectrumAnnotator::terminal_series_match_ratio) |
impl DefaultParamHandling for SpectrumAnnotator {
    const NAME: &'static str = "SpectrumAnnotator";

    fn defaults() -> Result<Param> {
        Self::default().parameter_tree()
    }

    /// # Errors
    ///
    /// [`Error::InvalidValue`] when
    /// [`top_n_fragment_errors`](SpectrumAnnotator::top_n_fragment_errors)
    /// exceeds the source's `unsigned`.
    fn to_param(&self) -> Result<Param> {
        self.parameter_tree()
    }

    /// Source `updateMembers_`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a negative `topNmatch_fragmenterrors`, which
    /// the source's `(unsigned)` conversion refuses with
    /// `Exception::ConversionError` (the parameter has no restriction, so
    /// `setParameters` itself accepts it). A value above
    /// [`MAX_ANNOTATION_TOP_N`](super::MAX_ANNOTATION_TOP_N) is accepted here and
    /// refused by [`add_ion_match_statistics`](SpectrumAnnotator::add_ion_match_statistics).
    fn update_members(param: &Param) -> Result<Self> {
        let flag = |key: &str| -> Result<bool> { param.value(key)?.to_bool() };
        let top_n = unsigned_parameter(param, TOP_N)?;
        Ok(Self {
            basic_statistics: flag("basic_statistics")?,
            list_of_ions_matched: flag("list_of_ions_matched")?,
            max_series: flag("max_series")?,
            sn_statistics: flag("S/N_statistics")?,
            precursor_statistics: flag("precursor_statistics")?,
            top_n_fragment_errors: usize::try_from(top_n).map_err(|_| {
                Error::InvalidValue(format!(
                    "SpectrumAnnotator: '{TOP_N}' = {top_n} does not fit this platform"
                ))
            })?,
            fragment_error_statistics: flag("fragmenterror_statistics")?,
            terminal_series_match_ratio: flag("terminal_series_match_ratio")?,
        })
    }
}

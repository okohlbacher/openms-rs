// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Flagging and filtering of MRM/SRM features against QC criteria.
//!
//! Ports `ANALYSIS/OPENSWATH/MRMFeatureFilter.h` and its `.cpp`.
//! [`MRMFeatureFilter`] checks the features of a
//! [`FeatureMap`](crate::kernel::FeatureMap) against the bounds of an
//! [`MRMFeatureQC`] and either records the verdict in metavalues (`flag`) or
//! removes the features and subordinates that fail (`filter`). It also
//! estimates QC bounds from standards, `%RSD` limits from pooled QC samples
//! and background interference from blanks.
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/MRM_FEATURE_QC_SUPPORT.md).
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::mrm_feature_filter::MRMFeatureFilter;
//! use openms::kernel::Feature;
//!
//! let filter = MRMFeatureFilter::default();
//! let mut a = Feature::default();
//! let mut b = Feature::default();
//! a.metadata.insert("native_id".into(), "a".into());
//! b.metadata.insert("native_id".into(), "b".into());
//! a.intensity = 3.0;
//! b.intensity = 4.0;
//! assert_eq!(filter.calculate_ion_ratio(&a, &b, "intensity").unwrap(), 0.75);
//! ```
//!
//! # Notes
//!
//! The source writes `OPENMS_LOG_DEBUG` messages when a metavalue or the
//! internal standard is missing; debug logging has no equivalent here and the
//! messages are not produced.
//!
//! [`MRMFeatureFilter`]: crate::analysis::openswath::mrm_feature_filter::MRMFeatureFilter
//! [`MRMFeatureQC`]: crate::analysis::openswath::mrm_feature_qc::MRMFeatureQC

use crate::analysis::openswath::mrm_feature_qc::{
    ComponentGroupQCs, ComponentQCs, MRMFeatureQC, MetaValueBounds,
};
use crate::analysis::openswath::mrm_feature_support::{meta_double, meta_string, native_id};
use crate::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition;
use crate::analysis::openswath::targeted_experiment::TargetedExperiment;
use crate::kernel::{Feature, FeatureMap};
use crate::metadata::{MetaValue, MetaValueData};
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use std::collections::BTreeMap;

/// The source `DefaultParamHandler` name.
pub const MRM_FEATURE_FILTER_NAME: &str = "MRMFeatureFilter";

/// What [`MRMFeatureFilter`] does with a feature or subordinate that fails.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum FlagOrFilter {
    /// Keep everything and record `..._pass` and `..._message` metavalues.
    #[default]
    Flag,
    /// Remove the failing components and component groups.
    Filter,
}

impl FlagOrFilter {
    /// The parameter string, `flag` or `filter`.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Flag => "flag",
            Self::Filter => "filter",
        }
    }
}

/// A number the range helpers of the filter work on: `double` and `Int`.
pub trait RangeValue: Copy + PartialOrd {
    /// The additive zero, `T(0)` in the source.
    const ZERO: Self;
}
impl RangeValue for f64 {
    const ZERO: Self = 0.0;
}
impl RangeValue for i32 {
    const ZERO: Self = 0;
}

/// An absolute quantitation method as far as
/// [`MRMFeatureFilter::transfer_lloq_and_uloq_to_calculated_concentration_bounds`]
/// reads it.
///
/// Source `AbsoluteQuantitationMethod` (`ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h`)
/// is not ported yet; its port implements this trait with its `getComponentName`,
/// `getLLOQ` and `getULOQ`.
pub trait QuantitationBounds {
    /// The component the method quantifies (`getComponentName`).
    fn component_name(&self) -> &str;
    /// The lower limit of quantitation (`getLLOQ`).
    fn lloq(&self) -> f64;
    /// The upper limit of quantitation (`getULOQ`).
    fn uloq(&self) -> f64;
}

/// Counts of heavy and light labels and of transition types in a component
/// group, as returned by
/// [`MRMFeatureFilter::count_labels_and_transition_types`] (the source's
/// `std::map<std::string, int>` with keys `n_heavy`, `n_light`,
/// `n_quantifying`, `n_identifying`, `n_detecting` and `n_transitions`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct LabelCounts {
    /// Number of subordinates whose `LabelType` is `Heavy`.
    pub n_heavy: i32,
    /// Number of subordinates whose `LabelType` is `Light`.
    pub n_light: i32,
    /// Number of quantifying transitions.
    pub n_quantifying: i32,
    /// Number of identifying transitions.
    pub n_identifying: i32,
    /// Number of detecting transitions.
    pub n_detecting: i32,
    /// Number of subordinates.
    pub n_transitions: i32,
}

/// Flags or filters features that do not pass QC criteria.
///
/// Source `MRMFeatureFilter`, a `DefaultParamHandler` with the single
/// parameter `flag_or_filter` (`flag` by default, tagged `advanced`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MRMFeatureFilter {
    /// Flag or filter (remove) components and component groups that fail.
    pub flag_or_filter: FlagOrFilter,
}

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

fn at<'a, T>(items: &'a [T], index: usize, what: &str) -> Result<&'a T> {
    items.get(index).ok_or_else(|| {
        bad(format!(
            "{what}: the filter values hold {} entries, but index {index} was requested (the source's std::vector::at throws)",
            items.len()
        ))
    })
}

fn bounds_at<'a>(
    map: &'a BTreeMap<String, MetaValueBounds>,
    key: &str,
) -> Result<&'a MetaValueBounds> {
    map.get(key).ok_or_else(|| {
        bad(format!(
            "metavalue QC '{key}' is missing from the filter values (the source's std::map::at throws)"
        ))
    })
}

/// Source `checkRange`: `value >= value_l && value <= value_u`.
fn in_range<T: RangeValue>(value: T, value_l: T, value_u: T) -> bool {
    value >= value_l && value <= value_u
}

fn component_group_scalars(q: &ComponentGroupQCs) -> ([f64; 8], [i32; 12]) {
    (
        [
            q.retention_time_l,
            q.retention_time_u,
            q.intensity_l,
            q.intensity_u,
            q.overall_quality_l,
            q.overall_quality_u,
            q.ion_ratio_l,
            q.ion_ratio_u,
        ],
        [
            q.n_heavy_l,
            q.n_heavy_u,
            q.n_light_l,
            q.n_light_u,
            q.n_detecting_l,
            q.n_detecting_u,
            q.n_quantifying_l,
            q.n_quantifying_u,
            q.n_identifying_l,
            q.n_identifying_u,
            q.n_transitions_l,
            q.n_transitions_u,
        ],
    )
}

fn set_component_group_scalars(q: &mut ComponentGroupQCs, doubles: [f64; 8], ints: [i32; 12]) {
    q.retention_time_l = doubles[0];
    q.retention_time_u = doubles[1];
    q.intensity_l = doubles[2];
    q.intensity_u = doubles[3];
    q.overall_quality_l = doubles[4];
    q.overall_quality_u = doubles[5];
    q.ion_ratio_l = doubles[6];
    q.ion_ratio_u = doubles[7];
    q.n_heavy_l = ints[0];
    q.n_heavy_u = ints[1];
    q.n_light_l = ints[2];
    q.n_light_u = ints[3];
    q.n_detecting_l = ints[4];
    q.n_detecting_u = ints[5];
    q.n_quantifying_l = ints[6];
    q.n_quantifying_u = ints[7];
    q.n_identifying_l = ints[8];
    q.n_identifying_u = ints[9];
    q.n_transitions_l = ints[10];
    q.n_transitions_u = ints[11];
}

fn component_scalars(q: &ComponentQCs) -> [f64; 6] {
    [
        q.retention_time_l,
        q.retention_time_u,
        q.intensity_l,
        q.intensity_u,
        q.overall_quality_l,
        q.overall_quality_u,
    ]
}

fn set_component_scalars(q: &mut ComponentQCs, doubles: [f64; 6]) {
    q.retention_time_l = doubles[0];
    q.retention_time_u = doubles[1];
    q.intensity_l = doubles[2];
    q.intensity_u = doubles[3];
    q.overall_quality_l = doubles[4];
    q.overall_quality_u = doubles[5];
}

/// `Int / size_t` as C++ evaluates it: the `Int` converts to `unsigned long`
/// (so a negative value wraps), the quotient converts back to `Int`.
fn int_div(value: i32, divisor: u64, what: &str) -> Result<i32> {
    if divisor == 0 {
        return Err(bad(format!(
            "{what}: the source divides an integer count by zero here, which terminates the process"
        )));
    }
    Ok(((i64::from(value) as u64) / divisor) as i32)
}

/// `mean != 0 ? std::sqrt(var) / mean * 100 : 0` over doubles.
fn perc_rsd(mean: f64, var: f64) -> f64 {
    if mean != 0.0 {
        var.sqrt() / mean * 100.0
    } else {
        0.0
    }
}

/// The same expression over `Int` operands, truncated back to `Int`.
fn perc_rsd_int(mean: i32, var: i32) -> i32 {
    if mean != 0 {
        (f64::from(var).sqrt() / f64::from(mean) * 100.0) as i32
    } else {
        0
    }
}

impl MRMFeatureFilter {
    /// A filter with the source defaults (`flag`).
    pub fn new() -> Self {
        Self::default()
    }

    /// The default parameters (source `getDefaultParameters`): `flag_or_filter`
    /// with the valid strings `flag` and `filter`.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which this constant tree cannot trigger.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The current parameters in the layout of [`defaults`](Self::defaults).
    ///
    /// # Errors
    ///
    /// As [`defaults`](Self::defaults).
    pub fn to_param(&self) -> Result<Param> {
        let entry = ParamEntry {
            name: "flag_or_filter".into(),
            description:
                "Flag or Filter (i.e., remove) Components or transitions that do not pass the QC."
                    .into(),
            value: ParamValue::String(self.flag_or_filter.name().into()),
            tags: ["advanced".to_owned()].into_iter().collect(),
            valid_strings: vec!["flag".into(), "filter".into()],
            ..ParamEntry::default()
        };
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: vec![entry],
            nodes: Vec::new(),
        })
    }

    /// Build a filter from parameters (source `setParameters` followed by
    /// `updateMembers_`); unknown parameters are ignored.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `flag_or_filter` is neither `flag` nor
    /// `filter` or is not a string.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`from_param`](Self::from_param), also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`from_param`](Self::from_param).
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(MRM_FEATURE_FILTER_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |complete| {
            let flag_or_filter = match complete.value("flag_or_filter")?.as_str()? {
                "flag" => FlagOrFilter::Flag,
                "filter" => FlagOrFilter::Filter,
                other => {
                    return Err(bad(format!(
                        "flag_or_filter must be 'flag' or 'filter', not '{other}'"
                    )));
                }
            };
            Ok(Self { flag_or_filter })
        })
    }

    /// Source `checkRange`: whether `value` lies in `[value_l, value_u]`.
    pub fn check_range<T: RangeValue>(&self, value: T, value_l: T, value_u: T) -> bool {
        in_range(value, value_l, value_u)
    }

    /// Source `updateRange`: widen `[value_l, value_u]` to include `value`.
    pub fn update_range<T: RangeValue>(&self, value: T, value_l: &mut T, value_u: &mut T) {
        if value < *value_l {
            *value_l = value;
        }
        if value > *value_u {
            *value_u = value;
        }
    }

    /// Source `setRange`: set the range to `[0, value]`, or `[value, 0]` for a
    /// negative `value` (a NaN `value` takes the second form).
    pub fn set_range<T: RangeValue>(&self, value: T, value_l: &mut T, value_u: &mut T) {
        if value >= T::ZERO {
            *value_l = T::ZERO;
            *value_u = value;
        } else {
            *value_l = value;
            *value_u = T::ZERO;
        }
    }

    /// Source `initRange`: set both bounds to `value`.
    pub fn init_range<T: RangeValue>(&self, value: T, value_l: &mut T, value_u: &mut T) {
        *value_l = value;
        *value_u = value;
    }

    /// Source `getUniqueSorted`: the messages sorted, without duplicates.
    pub fn get_unique_sorted(&self, messages: &[String]) -> Vec<String> {
        let mut unique = messages.to_vec();
        unique.sort();
        unique.dedup();
        unique
    }

    /// Count the heavy and light labels and the quantifying, identifying and
    /// detecting transitions of a component group.
    ///
    /// Each subordinate is matched to the first transition of `transitions`
    /// whose native ID equals the subordinate's `native_id` metavalue; a
    /// subordinate without a match is counted with the flags of a
    /// default-constructed transition (detecting and quantifying, not
    /// identifying). `LabelType` is read leniently as a string.
    pub fn count_labels_and_transition_types(
        &self,
        component_group: &Feature,
        transitions: &TargetedExperiment,
    ) -> LabelCounts {
        let mut counts = LabelCounts::default();
        let default_transition = ReactionMonitoringTransition::new();
        for subordinate in &component_group.subordinates {
            let wanted = subordinate
                .metadata
                .get("native_id")
                .and_then(|v| match v.data() {
                    MetaValueData::String(text) => Some(text.as_str()),
                    _ => None,
                });
            let transition = wanted
                .and_then(|id| {
                    transitions
                        .transitions()
                        .iter()
                        .find(|t| t.native_id() == id)
                })
                .unwrap_or(&default_transition);
            match meta_string(&subordinate.metadata, "LabelType").as_str() {
                "Heavy" => counts.n_heavy += 1,
                "Light" => counts.n_light += 1,
                _ => {}
            }
            if transition.is_quantifying_transition() {
                counts.n_quantifying += 1;
            }
            if transition.is_identifying_transition() {
                counts.n_identifying += 1;
            }
            if transition.is_detecting_transition() {
                counts.n_detecting += 1;
            }
            counts.n_transitions += 1;
        }
        counts
    }

    /// The ratio between two transitions on a feature: the `intensity` member
    /// or a metavalue.
    ///
    /// With `feature_name == "intensity"` the ratio is
    /// `component_1.intensity / component_2.intensity` when both carry a
    /// `native_id`, the intensity of `component_1` alone when only it does,
    /// and `0` otherwise. For any other name the metavalue of that name takes
    /// the place of the intensity, and `0` is returned when `component_1` does
    /// not have it. A zero denominator yields the IEEE infinity, as in the
    /// source.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a metavalue read is neither an integer nor
    /// a float (the source's `(double)` cast throws on an empty value and
    /// reads garbage for other types).
    pub fn calculate_ion_ratio(
        &self,
        component_1: &Feature,
        component_2: &Feature,
        feature_name: &str,
    ) -> Result<f64> {
        let has = |f: &Feature, key: &str| f.metadata.contains_key(key);
        if feature_name == "intensity" {
            if has(component_1, "native_id") && has(component_2, "native_id") {
                Ok(f64::from(component_1.intensity) / f64::from(component_2.intensity))
            } else if has(component_1, "native_id") {
                Ok(f64::from(component_1.intensity))
            } else {
                Ok(0.0)
            }
        } else if has(component_1, feature_name) && has(component_2, feature_name) {
            let feature_1 = meta_double(&component_1.metadata, feature_name)?;
            let feature_2 = meta_double(&component_2.metadata, feature_name)?;
            Ok(feature_1 / feature_2)
        } else if has(component_1, feature_name) {
            meta_double(&component_1.metadata, feature_name)
        } else {
            Ok(0.0)
        }
    }

    /// The absolute retention time difference between two features.
    pub fn calculate_rt_difference(&self, component_1: &Feature, component_2: &Feature) -> f64 {
        (component_1.rt - component_2.rt).abs()
    }

    /// The resolution between two features, `2 * |RT2 - RT1| / (W1 + W2)`.
    ///
    /// A width is `width_at_50` (the FWHM) times `1.7` when the feature has
    /// that metavalue, otherwise the feature's `width` member. A
    /// non-positive width sum yields `0`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `width_at_50` is neither an integer nor a float.
    pub fn calculate_resolution(
        &self,
        component_1: &Feature,
        component_2: &Feature,
    ) -> Result<f64> {
        let rt_diff = (component_1.rt - component_2.rt).abs();
        let width = |f: &Feature| -> Result<f64> {
            if f.metadata.contains_key("width_at_50") {
                Ok(meta_double(&f.metadata, "width_at_50")? * 1.7)
            } else {
                Ok(f64::from(f.width))
            }
        };
        let width_sum = width(component_1)? + width(component_2)?;
        if width_sum <= 0.0 {
            return Ok(0.0);
        }
        Ok(2.0 * rt_diff / width_sum)
    }

    /// Whether the metavalue `meta_value_key` of `component` lies in
    /// `[meta_value_l, meta_value_u]`.
    ///
    /// Returns `(passes, key_exists)`; a missing key passes.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the metavalue is neither an integer nor a float.
    pub fn check_meta_value(
        &self,
        component: &Feature,
        meta_value_key: &str,
        meta_value_l: f64,
        meta_value_u: f64,
    ) -> Result<(bool, bool)> {
        if component.metadata.contains_key(meta_value_key) {
            let value = meta_double(&component.metadata, meta_value_key)?;
            Ok((in_range(value, meta_value_l, meta_value_u), true))
        } else {
            Ok((true, false))
        }
    }

    /// Widen the range to include the metavalue (source `updateMetaValue`).
    /// Returns whether the key exists.
    ///
    /// # Errors
    ///
    /// As [`check_meta_value`](Self::check_meta_value).
    pub fn update_meta_value(
        &self,
        component: &Feature,
        meta_value_key: &str,
        meta_value_l: &mut f64,
        meta_value_u: &mut f64,
    ) -> Result<bool> {
        if component.metadata.contains_key(meta_value_key) {
            let value = meta_double(&component.metadata, meta_value_key)?;
            self.update_range(value, meta_value_l, meta_value_u);
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Set the range to `[0, value]` or `[value, 0]` from the metavalue
    /// (source `setMetaValue`). Returns whether the key exists.
    ///
    /// # Errors
    ///
    /// As [`check_meta_value`](Self::check_meta_value).
    pub fn set_meta_value(
        &self,
        component: &Feature,
        meta_value_key: &str,
        meta_value_l: &mut f64,
        meta_value_u: &mut f64,
    ) -> Result<bool> {
        if component.metadata.contains_key(meta_value_key) {
            let value = meta_double(&component.metadata, meta_value_key)?;
            self.set_range(value, meta_value_l, meta_value_u);
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Set both bounds to the metavalue (source `initMetaValue`). Returns
    /// whether the key exists.
    ///
    /// # Errors
    ///
    /// As [`check_meta_value`](Self::check_meta_value).
    pub fn init_meta_value(
        &self,
        component: &Feature,
        meta_value_key: &str,
        meta_value_l: &mut f64,
        meta_value_u: &mut f64,
    ) -> Result<bool> {
        if component.metadata.contains_key(meta_value_key) {
            let value = meta_double(&component.metadata, meta_value_key)?;
            self.init_range(value, meta_value_l, meta_value_u);
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Flag or filter the features and subordinates of `features` against
    /// `filter_criteria`.
    ///
    /// Features are matched to component group QCs by their `PeptideRef`
    /// metavalue and subordinates to component QCs by their `native_id`.
    /// Each subordinate is checked against the component QC with its name; the
    /// component group QC of the feature is evaluated once for **every**
    /// subordinate, as in the source, so a group with `k` subordinates
    /// contributes `k` times as many tests and failure messages to the group
    /// score. The group checks are retention time, intensity, overall quality,
    /// the six label and transition counts, each ion ratio pair present among
    /// the subordinates, and the metavalue bounds.
    ///
    /// In flag mode a subordinate gets `QC_transition_score`,
    /// `QC_transition_pass` (integer `1` or `0`) and `QC_transition_message`
    /// (the sorted unique failure names), and the feature gets
    /// `QC_transition_group_score`, `QC_transition_group_pass` and
    /// `QC_transition_group_message`. In filter mode the failing subordinates
    /// are removed, features that fail or keep no subordinate are removed, and
    /// only the scores are written.
    ///
    /// The source replaces the whole map with a fresh one in filter mode,
    /// which also discards its identifier, protein identifications, data
    /// processing and the other map-level members. This port replaces only the
    /// features and keeps the map-level members.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a metavalue a bound applies to is neither
    /// an integer nor a float.
    pub fn filter_feature_map(
        &self,
        features: &mut FeatureMap,
        filter_criteria: &MRMFeatureQC,
        transitions: &TargetedExperiment,
    ) -> Result<()> {
        let mut staged = features.features.clone();
        let mut filtered: Vec<Feature> = Vec::new();
        for feature in &mut staged {
            let component_group_name = meta_string(&feature.metadata, "PeptideRef");
            let counts = self.count_labels_and_transition_types(feature, transitions);
            let mut subordinates_filtered: Vec<Feature> = Vec::new();
            let mut cg_qc_pass = true;
            let mut cg_fail: Vec<String> = Vec::new();
            let mut cg_tests_count: u32 = 0;

            for sub_it in 0..feature.subordinates.len() {
                let component_name = native_id(&feature.subordinates[sub_it]);
                let mut c_qc_pass = true;
                let mut c_fail: Vec<String> = Vec::new();

                for cg in &filter_criteria.component_group_qcs {
                    if cg.component_group_name != component_group_name {
                        continue;
                    }
                    let checks: [(bool, &str); 9] = [
                        (
                            in_range(feature.rt, cg.retention_time_l, cg.retention_time_u),
                            "retention_time",
                        ),
                        (
                            in_range(f64::from(feature.intensity), cg.intensity_l, cg.intensity_u),
                            "intensity",
                        ),
                        (
                            in_range(
                                f64::from(feature.quality),
                                cg.overall_quality_l,
                                cg.overall_quality_u,
                            ),
                            "overall_quality",
                        ),
                        (
                            in_range(counts.n_heavy, cg.n_heavy_l, cg.n_heavy_u),
                            "n_heavy",
                        ),
                        (
                            in_range(counts.n_light, cg.n_light_l, cg.n_light_u),
                            "n_light",
                        ),
                        (
                            in_range(counts.n_detecting, cg.n_detecting_l, cg.n_detecting_u),
                            "n_detecting",
                        ),
                        (
                            in_range(counts.n_quantifying, cg.n_quantifying_l, cg.n_quantifying_u),
                            "n_quantifying",
                        ),
                        (
                            in_range(counts.n_identifying, cg.n_identifying_l, cg.n_identifying_u),
                            "n_identifying",
                        ),
                        (
                            in_range(counts.n_transitions, cg.n_transitions_l, cg.n_transitions_u),
                            "n_transitions",
                        ),
                    ];
                    for (pass, name) in checks {
                        if !pass {
                            cg_qc_pass = false;
                            cg_fail.push(name.to_owned());
                        }
                    }
                    cg_tests_count += 9;

                    // ion ratio QC
                    for sub_it2 in 0..feature.subordinates.len() {
                        let component_name2 = native_id(&feature.subordinates[sub_it2]);
                        if !cg.ion_ratio_pair_name_1.is_empty()
                            && !cg.ion_ratio_pair_name_2.is_empty()
                            && cg.ion_ratio_pair_name_1 == component_name
                            && cg.ion_ratio_pair_name_2 == component_name2
                        {
                            let ion_ratio = self.calculate_ion_ratio(
                                &feature.subordinates[sub_it],
                                &feature.subordinates[sub_it2],
                                &cg.ion_ratio_feature_name,
                            )?;
                            if !in_range(ion_ratio, cg.ion_ratio_l, cg.ion_ratio_u) {
                                cg_qc_pass = false;
                                cg_fail.push(format!(
                                    "ion_ratio_pair[{component_name}/{component_name2}]"
                                ));
                            }
                            cg_tests_count += 1;
                        }
                    }

                    for (key, (lower, upper)) in &cg.meta_value_qc {
                        let (pass, exists) = self.check_meta_value(feature, key, *lower, *upper)?;
                        if !pass {
                            cg_qc_pass = false;
                            cg_fail.push(key.clone());
                        }
                        if exists {
                            cg_tests_count += 1;
                        }
                    }
                }

                let mut c_tests_count: u32 = 0;
                for c in &filter_criteria.component_qcs {
                    if c.component_name != component_name {
                        continue;
                    }
                    let sub = &feature.subordinates[sub_it];
                    let checks: [(bool, &str); 3] = [
                        (
                            in_range(sub.rt, c.retention_time_l, c.retention_time_u),
                            "retention_time",
                        ),
                        (
                            in_range(f64::from(sub.intensity), c.intensity_l, c.intensity_u),
                            "intensity",
                        ),
                        (
                            in_range(
                                f64::from(sub.quality),
                                c.overall_quality_l,
                                c.overall_quality_u,
                            ),
                            "overall_quality",
                        ),
                    ];
                    for (pass, name) in checks {
                        if !pass {
                            c_qc_pass = false;
                            c_fail.push(name.to_owned());
                        }
                    }
                    c_tests_count += 3;
                    for (key, (lower, upper)) in &c.meta_value_qc {
                        let (pass, exists) = self.check_meta_value(sub, key, *lower, *upper)?;
                        if !pass {
                            c_qc_pass = false;
                            c_fail.push(key.clone());
                        }
                        if exists {
                            c_tests_count += 1;
                        }
                    }
                }

                let c_score = score(c_fail.len(), c_tests_count);
                let sub = &mut feature.subordinates[sub_it];
                put_score(sub, "QC_transition_score", c_score)?;
                match (c_qc_pass, self.flag_or_filter) {
                    (true, FlagOrFilter::Filter) => subordinates_filtered.push(sub.clone()),
                    (true, FlagOrFilter::Flag) => {
                        put_verdict(sub, "QC_transition", true, Vec::new());
                    }
                    (false, FlagOrFilter::Filter) => {}
                    (false, FlagOrFilter::Flag) => {
                        let unique = self.get_unique_sorted(&c_fail);
                        put_verdict(sub, "QC_transition", false, unique);
                    }
                }
            }

            let cg_score = score(cg_fail.len(), cg_tests_count);
            put_score(feature, "QC_transition_group_score", cg_score)?;
            match (cg_qc_pass, self.flag_or_filter) {
                (true, FlagOrFilter::Filter) => {
                    if !subordinates_filtered.is_empty() {
                        let mut kept = feature.clone();
                        kept.subordinates = subordinates_filtered;
                        filtered.push(kept);
                    }
                }
                (true, FlagOrFilter::Flag) => {
                    put_verdict(feature, "QC_transition_group", true, Vec::new());
                }
                (false, FlagOrFilter::Filter) => {}
                (false, FlagOrFilter::Flag) => {
                    let unique = self.get_unique_sorted(&cg_fail);
                    put_verdict(feature, "QC_transition_group", false, unique);
                }
            }
        }
        features.features = match self.flag_or_filter {
            FlagOrFilter::Filter => filtered,
            FlagOrFilter::Flag => staged,
        };
        Ok(())
    }

    /// Flag or filter against the `%RSD` limits in `filter_criteria`, using
    /// the `%RSD` values of [`estimate_perc_rsd`](Self::estimate_perc_rsd) in
    /// `filter_values`.
    ///
    /// Each criterion is compared with the **upper** bound of the matching
    /// entry of `filter_values` (the position in the list, not the name,
    /// pairs them). Labels and transition counts are not checked. Metavalues
    /// written are `QC_transition_%RSD_score`, `QC_transition_%RSD_pass`,
    /// `QC_transition_%RSD_message` on subordinates and
    /// `QC_transition_group_%RSD_score`, `..._pass`, `..._message` on
    /// features. The map-level members are kept in filter mode, as in
    /// [`filter_feature_map`](Self::filter_feature_map).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `filter_values` holds fewer entries than
    /// `filter_criteria` or lacks a metavalue of a criterion (the source's
    /// `.at` throws `std::out_of_range`).
    pub fn filter_feature_map_perc_rsd(
        &self,
        features: &mut FeatureMap,
        filter_criteria: &MRMFeatureQC,
        filter_values: &MRMFeatureQC,
    ) -> Result<()> {
        let mut staged = features.features.clone();
        let mut filtered: Vec<Feature> = Vec::new();
        for feature in &mut staged {
            let component_group_name = meta_string(&feature.metadata, "PeptideRef");
            let mut subordinates_filtered: Vec<Feature> = Vec::new();
            let mut cg_qc_pass = true;
            let mut cg_fail: Vec<String> = Vec::new();
            let mut cg_tests_count: u32 = 0;

            for sub_it in 0..feature.subordinates.len() {
                let component_name = native_id(&feature.subordinates[sub_it]);
                let mut c_qc_pass = true;
                let mut c_fail: Vec<String> = Vec::new();

                for (cg_qc_it, cg) in filter_criteria.component_group_qcs.iter().enumerate() {
                    if cg.component_group_name != component_group_name {
                        continue;
                    }
                    let fv = at(
                        &filter_values.component_group_qcs,
                        cg_qc_it,
                        "component group",
                    )?;
                    let checks: [(bool, &str); 3] = [
                        (
                            in_range(
                                fv.retention_time_u,
                                cg.retention_time_l,
                                cg.retention_time_u,
                            ),
                            "retention_time",
                        ),
                        (
                            in_range(fv.intensity_u, cg.intensity_l, cg.intensity_u),
                            "intensity",
                        ),
                        (
                            in_range(
                                fv.overall_quality_u,
                                cg.overall_quality_l,
                                cg.overall_quality_u,
                            ),
                            "overall_quality",
                        ),
                    ];
                    for (pass, name) in checks {
                        if !pass {
                            cg_qc_pass = false;
                            cg_fail.push(name.to_owned());
                        }
                    }
                    cg_tests_count += 3;
                    if !cg.ion_ratio_pair_name_1.is_empty() && !cg.ion_ratio_pair_name_2.is_empty()
                    {
                        if !in_range(fv.ion_ratio_u, cg.ion_ratio_l, cg.ion_ratio_u) {
                            cg_qc_pass = false;
                            cg_fail.push(format!(
                                "ion_ratio_pair[{}/{}]",
                                cg.ion_ratio_pair_name_1, cg.ion_ratio_pair_name_2
                            ));
                        }
                        cg_tests_count += 1;
                    }
                    for (key, (lower, upper)) in &cg.meta_value_qc {
                        let observed = bounds_at(&fv.meta_value_qc, key)?.1;
                        if !in_range(observed, *lower, *upper) {
                            cg_qc_pass = false;
                            cg_fail.push(key.clone());
                        }
                        cg_tests_count += 1;
                    }
                }

                let mut c_tests_count: u32 = 0;
                for (c_qc_it, c) in filter_criteria.component_qcs.iter().enumerate() {
                    if c.component_name != component_name {
                        continue;
                    }
                    let fv = at(&filter_values.component_qcs, c_qc_it, "component")?;
                    let checks: [(bool, &str); 3] = [
                        (
                            in_range(fv.retention_time_u, c.retention_time_l, c.retention_time_u),
                            "retention_time",
                        ),
                        (
                            in_range(fv.intensity_u, c.intensity_l, c.intensity_u),
                            "intensity",
                        ),
                        (
                            in_range(
                                fv.overall_quality_u,
                                c.overall_quality_l,
                                c.overall_quality_u,
                            ),
                            "overall_quality",
                        ),
                    ];
                    for (pass, name) in checks {
                        if !pass {
                            c_qc_pass = false;
                            c_fail.push(name.to_owned());
                        }
                    }
                    c_tests_count += 3;
                    for (key, (lower, upper)) in &c.meta_value_qc {
                        let observed = bounds_at(&fv.meta_value_qc, key)?.1;
                        if !in_range(observed, *lower, *upper) {
                            c_qc_pass = false;
                            c_fail.push(key.clone());
                        }
                        c_tests_count += 1;
                    }
                }

                let c_score = score(c_fail.len(), c_tests_count);
                let sub = &mut feature.subordinates[sub_it];
                put_score(sub, "QC_transition_%RSD_score", c_score)?;
                match (c_qc_pass, self.flag_or_filter) {
                    (true, FlagOrFilter::Filter) => subordinates_filtered.push(sub.clone()),
                    (true, FlagOrFilter::Flag) => {
                        put_verdict(sub, "QC_transition_%RSD", true, Vec::new());
                    }
                    (false, FlagOrFilter::Filter) => {}
                    (false, FlagOrFilter::Flag) => {
                        let unique = self.get_unique_sorted(&c_fail);
                        put_verdict(sub, "QC_transition_%RSD", false, unique);
                    }
                }
            }

            let cg_score = score(cg_fail.len(), cg_tests_count);
            put_score(feature, "QC_transition_group_%RSD_score", cg_score)?;
            match (cg_qc_pass, self.flag_or_filter) {
                (true, FlagOrFilter::Filter) => {
                    if !subordinates_filtered.is_empty() {
                        let mut kept = feature.clone();
                        kept.subordinates = subordinates_filtered;
                        filtered.push(kept);
                    }
                }
                (true, FlagOrFilter::Flag) => {
                    put_verdict(feature, "QC_transition_group_%RSD", true, Vec::new());
                }
                (false, FlagOrFilter::Filter) => {}
                (false, FlagOrFilter::Flag) => {
                    let unique = self.get_unique_sorted(&cg_fail);
                    put_verdict(feature, "QC_transition_group_%RSD", false, unique);
                }
            }
        }
        features.features = match self.flag_or_filter {
            FlagOrFilter::Filter => filtered,
            FlagOrFilter::Flag => staged,
        };
        Ok(())
    }

    /// Flag or filter on the percent background interference, only on
    /// intensity.
    ///
    /// The interference of a component or group is the upper intensity bound
    /// of the matching entry of `filter_values` (from
    /// [`estimate_background_interferences`](Self::estimate_background_interferences))
    /// divided by the feature intensity times 100, which is compared with the
    /// intensity bounds of `filter_criteria`. Metavalues written use the
    /// prefix `QC_transition_%BackgroundInterference` and
    /// `QC_transition_group_%BackgroundInterference`. A feature intensity of
    /// zero yields an infinite or NaN interference, as in the source.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `filter_values` holds fewer entries than
    /// `filter_criteria`.
    pub fn filter_feature_map_background_interference(
        &self,
        features: &mut FeatureMap,
        filter_criteria: &MRMFeatureQC,
        filter_values: &MRMFeatureQC,
    ) -> Result<()> {
        let mut staged = features.features.clone();
        let mut filtered: Vec<Feature> = Vec::new();
        for feature in &mut staged {
            let component_group_name = meta_string(&feature.metadata, "PeptideRef");
            let mut subordinates_filtered: Vec<Feature> = Vec::new();
            let mut cg_qc_pass = true;
            let mut cg_fail: Vec<String> = Vec::new();
            let mut cg_tests_count: u32 = 0;

            for sub_it in 0..feature.subordinates.len() {
                let component_name = native_id(&feature.subordinates[sub_it]);
                let mut c_qc_pass = true;
                let mut c_fail: Vec<String> = Vec::new();

                for (cg_qc_it, cg) in filter_criteria.component_group_qcs.iter().enumerate() {
                    if cg.component_group_name != component_group_name {
                        continue;
                    }
                    let fv = at(
                        &filter_values.component_group_qcs,
                        cg_qc_it,
                        "component group",
                    )?;
                    let perc = fv.intensity_u / f64::from(feature.intensity) * 100.0;
                    if !in_range(perc, cg.intensity_l, cg.intensity_u) {
                        cg_qc_pass = false;
                        cg_fail.push("intensity".to_owned());
                    }
                    cg_tests_count += 1;
                }

                let mut c_tests_count: u32 = 0;
                for (c_qc_it, c) in filter_criteria.component_qcs.iter().enumerate() {
                    if c.component_name != component_name {
                        continue;
                    }
                    let fv = at(&filter_values.component_qcs, c_qc_it, "component")?;
                    let perc =
                        fv.intensity_u / f64::from(feature.subordinates[sub_it].intensity) * 100.0;
                    if !in_range(perc, c.intensity_l, c.intensity_u) {
                        c_qc_pass = false;
                        c_fail.push("intensity".to_owned());
                    }
                    c_tests_count += 1;
                }

                let c_score = score(c_fail.len(), c_tests_count);
                let sub = &mut feature.subordinates[sub_it];
                put_score(sub, "QC_transition_%BackgroundInterference_score", c_score)?;
                match (c_qc_pass, self.flag_or_filter) {
                    (true, FlagOrFilter::Filter) => subordinates_filtered.push(sub.clone()),
                    (true, FlagOrFilter::Flag) => {
                        put_verdict(
                            sub,
                            "QC_transition_%BackgroundInterference",
                            true,
                            Vec::new(),
                        );
                    }
                    (false, FlagOrFilter::Filter) => {}
                    (false, FlagOrFilter::Flag) => {
                        let unique = self.get_unique_sorted(&c_fail);
                        put_verdict(sub, "QC_transition_%BackgroundInterference", false, unique);
                    }
                }
            }

            let cg_score = score(cg_fail.len(), cg_tests_count);
            put_score(
                feature,
                "QC_transition_group_%BackgroundInterference_score",
                cg_score,
            )?;
            match (cg_qc_pass, self.flag_or_filter) {
                (true, FlagOrFilter::Filter) => {
                    if !subordinates_filtered.is_empty() {
                        let mut kept = feature.clone();
                        kept.subordinates = subordinates_filtered;
                        filtered.push(kept);
                    }
                }
                (true, FlagOrFilter::Flag) => {
                    put_verdict(
                        feature,
                        "QC_transition_group_%BackgroundInterference",
                        true,
                        Vec::new(),
                    );
                }
                (false, FlagOrFilter::Filter) => {}
                (false, FlagOrFilter::Flag) => {
                    let unique = self.get_unique_sorted(&cg_fail);
                    put_verdict(
                        feature,
                        "QC_transition_group_%BackgroundInterference",
                        false,
                        unique,
                    );
                }
            }
        }
        features.features = match self.flag_or_filter {
            FlagOrFilter::Filter => filtered,
            FlagOrFilter::Flag => staged,
        };
        Ok(())
    }

    /// Estimate the lower and upper bounds of the QC template from samples
    /// (typically standards).
    ///
    /// For every sample, feature and subordinate whose names appear in
    /// `filter_template`, the bounds are widened by the feature's retention
    /// time, intensity, overall quality, label and transition counts and
    /// ion ratios, and by the metavalues the template names. With
    /// `init_template_values` the first sample overwrites the template bounds
    /// instead of widening them.
    ///
    /// The component group metavalues use `update_meta_value` for the first
    /// initialising sample and `init_meta_value` for every other one, the
    /// reverse of the component metavalues; the source does the same and this
    /// port reproduces it (see the support document).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a metavalue a bound applies to is neither
    /// an integer nor a float.
    pub fn estimate_default_mrm_feature_qc_values(
        &self,
        samples: &[FeatureMap],
        filter_template: &mut MRMFeatureQC,
        transitions: &TargetedExperiment,
        init_template_values: bool,
    ) -> Result<()> {
        for (sample_it, sample) in samples.iter().enumerate() {
            let first = sample_it == 0 && init_template_values;
            for feature in &sample.features {
                let component_group_name = meta_string(&feature.metadata, "PeptideRef");
                let counts = self.count_labels_and_transition_types(feature, transitions);
                for (sub_it, sub) in feature.subordinates.iter().enumerate() {
                    let component_name = native_id(sub);
                    for cg in &mut filter_template.component_group_qcs {
                        if cg.component_group_name != component_group_name {
                            continue;
                        }
                        let apply_f64 = |value: f64, l: &mut f64, u: &mut f64| {
                            if first {
                                self.init_range(value, l, u);
                            } else {
                                self.update_range(value, l, u);
                            }
                        };
                        let apply_i32 = |value: i32, l: &mut i32, u: &mut i32| {
                            if first {
                                self.init_range(value, l, u);
                            } else {
                                self.update_range(value, l, u);
                            }
                        };
                        apply_f64(
                            feature.rt,
                            &mut cg.retention_time_l,
                            &mut cg.retention_time_u,
                        );
                        apply_f64(
                            f64::from(feature.intensity),
                            &mut cg.intensity_l,
                            &mut cg.intensity_u,
                        );
                        apply_f64(
                            f64::from(feature.quality),
                            &mut cg.overall_quality_l,
                            &mut cg.overall_quality_u,
                        );
                        apply_i32(counts.n_heavy, &mut cg.n_heavy_l, &mut cg.n_heavy_u);
                        apply_i32(counts.n_light, &mut cg.n_light_l, &mut cg.n_light_u);
                        apply_i32(
                            counts.n_detecting,
                            &mut cg.n_detecting_l,
                            &mut cg.n_detecting_u,
                        );
                        apply_i32(
                            counts.n_quantifying,
                            &mut cg.n_quantifying_l,
                            &mut cg.n_quantifying_u,
                        );
                        apply_i32(
                            counts.n_identifying,
                            &mut cg.n_identifying_l,
                            &mut cg.n_identifying_u,
                        );
                        apply_i32(
                            counts.n_transitions,
                            &mut cg.n_transitions_l,
                            &mut cg.n_transitions_u,
                        );

                        for sub2 in &feature.subordinates {
                            let component_name2 = native_id(sub2);
                            if !cg.ion_ratio_pair_name_1.is_empty()
                                && !cg.ion_ratio_pair_name_2.is_empty()
                                && cg.ion_ratio_pair_name_1 == component_name
                                && cg.ion_ratio_pair_name_2 == component_name2
                            {
                                let ion_ratio = self.calculate_ion_ratio(
                                    sub,
                                    sub2,
                                    &cg.ion_ratio_feature_name,
                                )?;
                                apply_f64(ion_ratio, &mut cg.ion_ratio_l, &mut cg.ion_ratio_u);
                            }
                        }

                        for (key, bounds) in &mut cg.meta_value_qc {
                            if first {
                                self.update_meta_value(feature, key, &mut bounds.0, &mut bounds.1)?;
                            } else {
                                self.init_meta_value(feature, key, &mut bounds.0, &mut bounds.1)?;
                            }
                        }
                    }

                    for c in &mut filter_template.component_qcs {
                        if c.component_name != component_name {
                            continue;
                        }
                        let apply = |value: f64, l: &mut f64, u: &mut f64| {
                            if first {
                                self.init_range(value, l, u);
                            } else {
                                self.update_range(value, l, u);
                            }
                        };
                        apply(sub.rt, &mut c.retention_time_l, &mut c.retention_time_u);
                        apply(
                            f64::from(sub.intensity),
                            &mut c.intensity_l,
                            &mut c.intensity_u,
                        );
                        apply(
                            f64::from(sub.quality),
                            &mut c.overall_quality_l,
                            &mut c.overall_quality_u,
                        );
                        for (key, bounds) in &mut c.meta_value_qc {
                            if first {
                                self.init_meta_value(sub, key, &mut bounds.0, &mut bounds.1)?;
                            } else {
                                self.update_meta_value(sub, key, &mut bounds.0, &mut bounds.1)?;
                            }
                        }
                    }
                    let _ = sub_it;
                }
            }
        }
        Ok(())
    }

    /// Copy the LLOQ and ULOQ of the quantitation methods into the
    /// `calculated_concentration` metavalue bounds of the component QCs of
    /// the same name. Methods with `LLOQ == 0` and `ULOQ == 0` are skipped.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a matching component QC has no
    /// `calculated_concentration` bound (the source's `.at` throws).
    pub fn transfer_lloq_and_uloq_to_calculated_concentration_bounds<M: QuantitationBounds>(
        &self,
        quantitation_methods: &[M],
        filter_template: &mut MRMFeatureQC,
    ) -> Result<()> {
        for method in quantitation_methods {
            if method.lloq() == 0.0 && method.uloq() == 0.0 {
                continue;
            }
            for c in &mut filter_template.component_qcs {
                if c.component_name == method.component_name() {
                    let bounds = c
                        .meta_value_qc
                        .get_mut("calculated_concentration")
                        .ok_or_else(|| {
                            bad("component QC has no 'calculated_concentration' metavalue bound")
                        })?;
                    bounds.0 = method.lloq();
                    bounds.1 = method.uloq();
                }
            }
        }
        Ok(())
    }

    /// Estimate the variability (`%RSD`) of the template's values over several
    /// samples (pooled QC or replicate unknown samples).
    ///
    /// The values of each sample are accumulated, their mean and variance
    /// (divided by `n - 1`) computed and `sqrt(var) / mean * 100` stored in
    /// `filter_template`; use it with
    /// [`filter_feature_map_perc_rsd`](Self::filter_feature_map_perc_rsd).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a metavalue is neither an integer nor a
    /// float, and when fewer than two samples are given while the template has
    /// component group QCs (the source divides an `Int` by zero).
    pub fn estimate_perc_rsd(
        &self,
        samples: &[FeatureMap],
        filter_template: &mut MRMFeatureQC,
        transitions: &TargetedExperiment,
    ) -> Result<()> {
        let mut filter_values = Vec::new();
        self.accumulate_filter_values(&mut filter_values, samples, filter_template, transitions)?;
        let mut filter_mean = MRMFeatureQC::default();
        self.calculate_filter_values_mean(&mut filter_mean, &filter_values, filter_template)?;
        let mut filter_var = MRMFeatureQC::default();
        self.calculate_filter_values_var(
            &mut filter_var,
            &filter_values,
            &filter_mean,
            filter_template,
        )?;
        self.calculate_filter_values_perc_rsd(filter_template, &filter_mean, &filter_var)
    }

    /// Estimate the background interference from blank samples: the average
    /// of the template's values, stored in `filter_template`, for use with
    /// [`filter_feature_map_background_interference`](Self::filter_feature_map_background_interference).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a metavalue is neither an integer nor a
    /// float, and when no sample is given while the template has component
    /// group QCs.
    pub fn estimate_background_interferences(
        &self,
        samples: &[FeatureMap],
        filter_template: &mut MRMFeatureQC,
        transitions: &TargetedExperiment,
    ) -> Result<()> {
        let mut filter_values = Vec::new();
        self.accumulate_filter_values(&mut filter_values, samples, filter_template, transitions)?;
        // The source passes the template as both output and template; the
        // clone has the same names and is zeroed the same way.
        let template = filter_template.clone();
        self.calculate_filter_values_mean(filter_template, &filter_values, &template)
    }

    /// Accumulate the values of each sample, one [`MRMFeatureQC`] per sample.
    ///
    /// Each value starts as a copy of `filter_template`; every matching
    /// feature then overwrites the retention time, intensity, quality, count,
    /// ion ratio and metavalue ranges with `[0, value]` or `[value, 0]`
    /// (`setRange`), so the result holds the **last** matching feature of
    /// the sample.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a metavalue is neither an integer nor a float.
    pub fn accumulate_filter_values(
        &self,
        filter_values: &mut Vec<MRMFeatureQC>,
        samples: &[FeatureMap],
        filter_template: &MRMFeatureQC,
        transitions: &TargetedExperiment,
    ) -> Result<()> {
        let mut produced = Vec::with_capacity(samples.len());
        for sample in samples {
            let mut filter_value = filter_template.clone();
            for feature in &sample.features {
                let component_group_name = meta_string(&feature.metadata, "PeptideRef");
                let counts = self.count_labels_and_transition_types(feature, transitions);
                for sub in &feature.subordinates {
                    let component_name = native_id(sub);
                    for cg in &mut filter_value.component_group_qcs {
                        if cg.component_group_name != component_group_name {
                            continue;
                        }
                        self.set_range(
                            feature.rt,
                            &mut cg.retention_time_l,
                            &mut cg.retention_time_u,
                        );
                        self.set_range(
                            f64::from(feature.intensity),
                            &mut cg.intensity_l,
                            &mut cg.intensity_u,
                        );
                        self.set_range(
                            f64::from(feature.quality),
                            &mut cg.overall_quality_l,
                            &mut cg.overall_quality_u,
                        );
                        self.set_range(counts.n_heavy, &mut cg.n_heavy_l, &mut cg.n_heavy_u);
                        self.set_range(counts.n_light, &mut cg.n_light_l, &mut cg.n_light_u);
                        self.set_range(
                            counts.n_detecting,
                            &mut cg.n_detecting_l,
                            &mut cg.n_detecting_u,
                        );
                        self.set_range(
                            counts.n_quantifying,
                            &mut cg.n_quantifying_l,
                            &mut cg.n_quantifying_u,
                        );
                        self.set_range(
                            counts.n_identifying,
                            &mut cg.n_identifying_l,
                            &mut cg.n_identifying_u,
                        );
                        self.set_range(
                            counts.n_transitions,
                            &mut cg.n_transitions_l,
                            &mut cg.n_transitions_u,
                        );
                        for sub2 in &feature.subordinates {
                            let component_name2 = native_id(sub2);
                            if !cg.ion_ratio_pair_name_1.is_empty()
                                && !cg.ion_ratio_pair_name_2.is_empty()
                                && cg.ion_ratio_pair_name_1 == component_name
                                && cg.ion_ratio_pair_name_2 == component_name2
                            {
                                let ion_ratio = self.calculate_ion_ratio(
                                    sub,
                                    sub2,
                                    &cg.ion_ratio_feature_name,
                                )?;
                                self.set_range(ion_ratio, &mut cg.ion_ratio_l, &mut cg.ion_ratio_u);
                            }
                        }
                        for (key, bounds) in &mut cg.meta_value_qc {
                            self.set_meta_value(feature, key, &mut bounds.0, &mut bounds.1)?;
                        }
                    }
                    for c in &mut filter_value.component_qcs {
                        if c.component_name != component_name {
                            continue;
                        }
                        self.set_range(sub.rt, &mut c.retention_time_l, &mut c.retention_time_u);
                        self.set_range(
                            f64::from(sub.intensity),
                            &mut c.intensity_l,
                            &mut c.intensity_u,
                        );
                        self.set_range(
                            f64::from(sub.quality),
                            &mut c.overall_quality_l,
                            &mut c.overall_quality_u,
                        );
                        for (key, bounds) in &mut c.meta_value_qc {
                            self.set_meta_value(sub, key, &mut bounds.0, &mut bounds.1)?;
                        }
                    }
                }
            }
            produced.push(filter_value);
        }
        filter_values.extend(produced);
        Ok(())
    }

    /// Copy `filter_template` into `filter_zeros` and set every bound, count
    /// and metavalue bound to zero.
    pub fn zero_filter_values(
        &self,
        filter_zeros: &mut MRMFeatureQC,
        filter_template: &MRMFeatureQC,
    ) {
        *filter_zeros = filter_template.clone();
        for cg in &mut filter_zeros.component_group_qcs {
            set_component_group_scalars(cg, [0.0; 8], [0; 12]);
            for bounds in cg.meta_value_qc.values_mut() {
                *bounds = (0.0, 0.0);
            }
        }
        for c in &mut filter_zeros.component_qcs {
            set_component_scalars(c, [0.0; 6]);
            for bounds in c.meta_value_qc.values_mut() {
                *bounds = (0.0, 0.0);
            }
        }
    }

    /// The mean of every QC value over `filter_values`, in `filter_mean`.
    ///
    /// Sums are taken sample by sample and divided by the number of samples.
    /// Counts are `Int` members: they are added with wrapping arithmetic and
    /// divided as the source does, an `Int` converted to `unsigned long`
    /// divided by the sample count.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a sample lacks a QC or metavalue entry of
    /// `filter_template`, and when `filter_values` is empty while the template
    /// has component group QCs (an integer division by zero in the source).
    pub fn calculate_filter_values_mean(
        &self,
        filter_mean: &mut MRMFeatureQC,
        filter_values: &[MRMFeatureQC],
        filter_template: &MRMFeatureQC,
    ) -> Result<()> {
        let mut mean = MRMFeatureQC::default();
        self.zero_filter_values(&mut mean, filter_template);
        for filter in filter_values {
            for (cg_it, m) in mean.component_group_qcs.iter_mut().enumerate() {
                let f = at(&filter.component_group_qcs, cg_it, "component group")?;
                let (mut doubles, mut ints) = component_group_scalars(m);
                let (fd, fi) = component_group_scalars(f);
                for k in 0..doubles.len() {
                    doubles[k] += fd[k];
                }
                for k in 0..ints.len() {
                    ints[k] = ints[k].wrapping_add(fi[k]);
                }
                set_component_group_scalars(m, doubles, ints);
                for (key, bounds) in &mut m.meta_value_qc {
                    let fb = bounds_at(&f.meta_value_qc, key)?;
                    bounds.0 += fb.0;
                    bounds.1 += fb.1;
                }
            }
            for (c_it, m) in mean.component_qcs.iter_mut().enumerate() {
                let f = at(&filter.component_qcs, c_it, "component")?;
                let mut doubles = component_scalars(m);
                let fd = component_scalars(f);
                for k in 0..doubles.len() {
                    doubles[k] += fd[k];
                }
                set_component_scalars(m, doubles);
                for (key, bounds) in &mut m.meta_value_qc {
                    let fb = bounds_at(&f.meta_value_qc, key)?;
                    bounds.0 += fb.0;
                    bounds.1 += fb.1;
                }
            }
        }
        let n = filter_values.len() as u64;
        for m in &mut mean.component_group_qcs {
            let (mut doubles, mut ints) = component_group_scalars(m);
            for value in &mut doubles {
                *value /= n as f64;
            }
            for value in &mut ints {
                *value = int_div(*value, n, "mean of the labels and transition counts")?;
            }
            set_component_group_scalars(m, doubles, ints);
            for bounds in m.meta_value_qc.values_mut() {
                bounds.0 /= n as f64;
                bounds.1 /= n as f64;
            }
        }
        for m in &mut mean.component_qcs {
            let mut doubles = component_scalars(m);
            for value in &mut doubles {
                *value /= n as f64;
            }
            set_component_scalars(m, doubles);
            for bounds in m.meta_value_qc.values_mut() {
                bounds.0 /= n as f64;
                bounds.1 /= n as f64;
            }
        }
        *filter_mean = mean;
        Ok(())
    }

    /// The sample variance (divided by `n - 1`) of every QC value over
    /// `filter_values`, about `filter_mean`, in `filter_var`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a sample or the mean lacks an entry of
    /// `filter_template`, and when exactly one sample is given while the
    /// template has component group QCs (an integer division by zero in the
    /// source). With no samples the divisor `n - 1` wraps to `2^64 - 1`, as in
    /// the source.
    pub fn calculate_filter_values_var(
        &self,
        filter_var: &mut MRMFeatureQC,
        filter_values: &[MRMFeatureQC],
        filter_mean: &MRMFeatureQC,
        filter_template: &MRMFeatureQC,
    ) -> Result<()> {
        let mut var = MRMFeatureQC::default();
        self.zero_filter_values(&mut var, filter_template);
        for filter in filter_values {
            for (cg_it, v) in var.component_group_qcs.iter_mut().enumerate() {
                let f = at(&filter.component_group_qcs, cg_it, "component group")?;
                let m = at(
                    &filter_mean.component_group_qcs,
                    cg_it,
                    "component group mean",
                )?;
                let (mut doubles, mut ints) = component_group_scalars(v);
                let (fd, fi) = component_group_scalars(f);
                let (md, mi) = component_group_scalars(m);
                for k in 0..doubles.len() {
                    let d = fd[k] - md[k];
                    doubles[k] += d * d;
                }
                for k in 0..ints.len() {
                    let d = f64::from(fi[k].wrapping_sub(mi[k]));
                    ints[k] = (f64::from(ints[k]) + d * d) as i32;
                }
                set_component_group_scalars(v, doubles, ints);
                for (key, bounds) in &mut v.meta_value_qc {
                    let fb = bounds_at(&f.meta_value_qc, key)?;
                    let mb = bounds_at(&m.meta_value_qc, key)?;
                    let d0 = fb.0 - mb.0;
                    let d1 = fb.1 - mb.1;
                    bounds.0 += d0 * d0;
                    bounds.1 += d1 * d1;
                }
            }
            for (c_it, v) in var.component_qcs.iter_mut().enumerate() {
                let f = at(&filter.component_qcs, c_it, "component")?;
                let m = at(&filter_mean.component_qcs, c_it, "component mean")?;
                let mut doubles = component_scalars(v);
                let fd = component_scalars(f);
                let md = component_scalars(m);
                for k in 0..doubles.len() {
                    let d = fd[k] - md[k];
                    doubles[k] += d * d;
                }
                set_component_scalars(v, doubles);
                for (key, bounds) in &mut v.meta_value_qc {
                    let fb = bounds_at(&f.meta_value_qc, key)?;
                    let mb = bounds_at(&m.meta_value_qc, key)?;
                    let d0 = fb.0 - mb.0;
                    let d1 = fb.1 - mb.1;
                    bounds.0 += d0 * d0;
                    bounds.1 += d1 * d1;
                }
            }
        }
        let n_minus_1 = (filter_values.len() as u64).wrapping_sub(1);
        for v in &mut var.component_group_qcs {
            let (mut doubles, mut ints) = component_group_scalars(v);
            for value in &mut doubles {
                *value /= n_minus_1 as f64;
            }
            for value in &mut ints {
                *value = int_div(
                    *value,
                    n_minus_1,
                    "variance of the labels and transition counts",
                )?;
            }
            set_component_group_scalars(v, doubles, ints);
            for bounds in v.meta_value_qc.values_mut() {
                bounds.0 /= n_minus_1 as f64;
                bounds.1 /= n_minus_1 as f64;
            }
        }
        for v in &mut var.component_qcs {
            let mut doubles = component_scalars(v);
            for value in &mut doubles {
                *value /= n_minus_1 as f64;
            }
            set_component_scalars(v, doubles);
            for bounds in v.meta_value_qc.values_mut() {
                bounds.0 /= n_minus_1 as f64;
                bounds.1 /= n_minus_1 as f64;
            }
        }
        *filter_var = var;
        Ok(())
    }

    /// The relative standard deviation in percent, `sqrt(var) / mean * 100`,
    /// of every QC value, written over the values of `filter_rsd`; a zero
    /// mean gives `0`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `filter_mean` or `filter_var` lacks an
    /// entry of `filter_rsd`.
    pub fn calculate_filter_values_perc_rsd(
        &self,
        filter_rsd: &mut MRMFeatureQC,
        filter_mean: &MRMFeatureQC,
        filter_var: &MRMFeatureQC,
    ) -> Result<()> {
        for (cg_it, r) in filter_rsd.component_group_qcs.iter_mut().enumerate() {
            let m = at(
                &filter_mean.component_group_qcs,
                cg_it,
                "component group mean",
            )?;
            let v = at(
                &filter_var.component_group_qcs,
                cg_it,
                "component group variance",
            )?;
            let (mut doubles, mut ints) = component_group_scalars(r);
            let (md, mi) = component_group_scalars(m);
            let (vd, vi) = component_group_scalars(v);
            for k in 0..doubles.len() {
                doubles[k] = perc_rsd(md[k], vd[k]);
            }
            for k in 0..ints.len() {
                ints[k] = perc_rsd_int(mi[k], vi[k]);
            }
            set_component_group_scalars(r, doubles, ints);
            for (key, bounds) in &mut r.meta_value_qc {
                let mb = bounds_at(&m.meta_value_qc, key)?;
                let vb = bounds_at(&v.meta_value_qc, key)?;
                bounds.0 = perc_rsd(mb.0, vb.0);
                bounds.1 = perc_rsd(mb.1, vb.1);
            }
        }
        for (c_it, r) in filter_rsd.component_qcs.iter_mut().enumerate() {
            let m = at(&filter_mean.component_qcs, c_it, "component mean")?;
            let v = at(&filter_var.component_qcs, c_it, "component variance")?;
            let mut doubles = component_scalars(r);
            let md = component_scalars(m);
            let vd = component_scalars(v);
            for k in 0..doubles.len() {
                doubles[k] = perc_rsd(md[k], vd[k]);
            }
            set_component_scalars(r, doubles);
            for (key, bounds) in &mut r.meta_value_qc {
                let mb = bounds_at(&m.meta_value_qc, key)?;
                let vb = bounds_at(&v.meta_value_qc, key)?;
                bounds.0 = perc_rsd(mb.0, vb.0);
                bounds.1 = perc_rsd(mb.1, vb.1);
            }
        }
        Ok(())
    }
}

/// `tests ? 1.0 - failures / (double)tests : 1.0`.
fn score(failures: usize, tests: u32) -> f64 {
    if tests != 0 {
        1.0 - failures as f64 / f64::from(tests)
    } else {
        1.0
    }
}

fn put_score(feature: &mut Feature, key: &str, value: f64) -> Result<()> {
    feature
        .metadata
        .insert(key.to_owned(), MetaValue::try_from(value)?);
    Ok(())
}

/// `setMetaValue(<prefix>_pass, bool)` and `(<prefix>_message, StringList)`.
/// A `bool` converts to `DataValue(int)` in the source, so the flag is the
/// integer `1` or `0`.
fn put_verdict(feature: &mut Feature, prefix: &str, pass: bool, messages: Vec<String>) {
    feature
        .metadata
        .insert(format!("{prefix}_pass"), MetaValue::from(i32::from(pass)));
    feature
        .metadata
        .insert(format!("{prefix}_message"), MetaValue::from(messages));
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Quality-control bounds for MRM/SRM components and component groups.
//!
//! Ports `ANALYSIS/OPENSWATH/MRMFeatureQC.h` (the `.cpp` is empty). An
//! [`MRMFeatureQC`] holds the lower and upper bounds that
//! [`MRMFeatureFilter`](crate::analysis::openswath::mrm_feature_filter::MRMFeatureFilter)
//! checks features against; it is read and written by
//! [`MRMFeatureQCFile`](crate::format::mrm_feature_qc_file::MRMFeatureQCFile).
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/MRM_FEATURE_QC_SUPPORT.md).
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::mrm_feature_qc::{ComponentGroupQCs, ComponentQCs, MRMFeatureQC};
//!
//! let mut qc = MRMFeatureQC::default();
//! qc.component_qcs.push(ComponentQCs {
//!     component_name: "heavy_13C6_glucose".into(),
//!     retention_time_l: 5.0,
//!     retention_time_u: 7.0,
//!     intensity_l: 1000.0,
//!     ..ComponentQCs::default()
//! });
//! qc.component_group_qcs.push(ComponentGroupQCs {
//!     component_group_name: "glucose".into(),
//!     ion_ratio_pair_name_1: "quantifier".into(),
//!     ion_ratio_pair_name_2: "qualifier".into(),
//!     ion_ratio_l: 0.8,
//!     ion_ratio_u: 1.2,
//!     ..ComponentGroupQCs::default()
//! });
//! assert_eq!(qc.component_qcs[0].retention_time_u, 7.0);
//! assert_eq!(qc.component_group_qcs[0].n_heavy_u, 100);
//! ```

use std::collections::BTreeMap;

/// The "unbounded" upper bound the source gives every `double` upper bound, `1e12`.
pub const DEFAULT_UPPER_BOUND: f64 = 1e12;

/// The default upper bound of every integer count (`Int ... { 100 }` in source).
pub const DEFAULT_COUNT_UPPER_BOUND: i32 = 100;

/// Lower and upper bound of a metavalue QC, `(lower, upper)`; both inclusive.
pub type MetaValueBounds = (f64, f64);

/// Quality controls for one component: a transition or sub-feature that
/// corresponds to a single peptide or metabolite.
///
/// Source `MRMFeatureQC::ComponentQCs`. The defaults are the source member
/// initialisers: every lower bound `0.0` and every upper bound `1e12`.
///
/// # Equality
///
/// `PartialEq` compares every member, the metavalue map included. The source
/// `operator==` compares the maps with `std::equal(first1, last1, first2, ...)`,
/// which reads `other.meta_value_qc` past its end when `other` holds fewer
/// entries (undefined behaviour) and reports equality when `other` merely
/// *starts* with this map's entries; the exact map comparison here is what the
/// source intends, and agrees with it whenever both maps have the same size.
#[derive(Clone, Debug, PartialEq)]
pub struct ComponentQCs {
    /// Name of the component.
    pub component_name: String,
    /// Retention time lower bound.
    pub retention_time_l: f64,
    /// Retention time upper bound.
    pub retention_time_u: f64,
    /// Intensity lower bound.
    pub intensity_l: f64,
    /// Intensity upper bound.
    pub intensity_u: f64,
    /// Overall quality lower bound.
    pub overall_quality_l: f64,
    /// Overall quality upper bound.
    pub overall_quality_u: f64,
    /// Bounds on feature metavalues, by metavalue name.
    pub meta_value_qc: BTreeMap<String, MetaValueBounds>,
}

impl Default for ComponentQCs {
    fn default() -> Self {
        Self {
            component_name: String::new(),
            retention_time_l: 0.0,
            retention_time_u: DEFAULT_UPPER_BOUND,
            intensity_l: 0.0,
            intensity_u: DEFAULT_UPPER_BOUND,
            overall_quality_l: 0.0,
            overall_quality_u: DEFAULT_UPPER_BOUND,
            meta_value_qc: BTreeMap::new(),
        }
    }
}

/// Quality controls within a component group: a transition group or feature
/// holding every transition of one component (peptide or metabolite).
///
/// Source `MRMFeatureQC::ComponentGroupQCs`; the defaults are the source member
/// initialisers (bounds `0.0` and `1e12`, counts `0` and `100`, empty names).
/// Equality is as for [`ComponentQCs`].
#[derive(Clone, Debug, PartialEq)]
pub struct ComponentGroupQCs {
    /// Name of the component group.
    pub component_group_name: String,
    /// Retention time lower bound.
    pub retention_time_l: f64,
    /// Retention time upper bound.
    pub retention_time_u: f64,
    /// Intensity lower bound.
    pub intensity_l: f64,
    /// Intensity upper bound.
    pub intensity_u: f64,
    /// Overall quality lower bound.
    pub overall_quality_l: f64,
    /// Overall quality upper bound.
    pub overall_quality_u: f64,
    /// Lower bound of the number of heavy-labelled transitions.
    pub n_heavy_l: i32,
    /// Upper bound of the number of heavy-labelled transitions.
    pub n_heavy_u: i32,
    /// Lower bound of the number of light-labelled transitions.
    pub n_light_l: i32,
    /// Upper bound of the number of light-labelled transitions.
    pub n_light_u: i32,
    /// Lower bound of the number of detecting transitions.
    pub n_detecting_l: i32,
    /// Upper bound of the number of detecting transitions.
    pub n_detecting_u: i32,
    /// Lower bound of the number of quantifying transitions.
    pub n_quantifying_l: i32,
    /// Upper bound of the number of quantifying transitions.
    pub n_quantifying_u: i32,
    /// Lower bound of the number of identifying transitions.
    pub n_identifying_l: i32,
    /// Upper bound of the number of identifying transitions.
    pub n_identifying_u: i32,
    /// Lower bound of the number of transitions.
    pub n_transitions_l: i32,
    /// Upper bound of the number of transitions.
    pub n_transitions_u: i32,
    /// Name (native ID) of the numerator component of the ion ratio. Empty
    /// disables the ion ratio check.
    pub ion_ratio_pair_name_1: String,
    /// Name (native ID) of the denominator component of the ion ratio.
    pub ion_ratio_pair_name_2: String,
    /// Ion ratio lower bound.
    pub ion_ratio_l: f64,
    /// Ion ratio upper bound.
    pub ion_ratio_u: f64,
    /// Name of the feature (`intensity` or a metavalue such as `peak_apex_int`)
    /// the ion ratio is computed on.
    pub ion_ratio_feature_name: String,
    /// Bounds on feature metavalues, by metavalue name.
    pub meta_value_qc: BTreeMap<String, MetaValueBounds>,
}

impl Default for ComponentGroupQCs {
    fn default() -> Self {
        Self {
            component_group_name: String::new(),
            retention_time_l: 0.0,
            retention_time_u: DEFAULT_UPPER_BOUND,
            intensity_l: 0.0,
            intensity_u: DEFAULT_UPPER_BOUND,
            overall_quality_l: 0.0,
            overall_quality_u: DEFAULT_UPPER_BOUND,
            n_heavy_l: 0,
            n_heavy_u: DEFAULT_COUNT_UPPER_BOUND,
            n_light_l: 0,
            n_light_u: DEFAULT_COUNT_UPPER_BOUND,
            n_detecting_l: 0,
            n_detecting_u: DEFAULT_COUNT_UPPER_BOUND,
            n_quantifying_l: 0,
            n_quantifying_u: DEFAULT_COUNT_UPPER_BOUND,
            n_identifying_l: 0,
            n_identifying_u: DEFAULT_COUNT_UPPER_BOUND,
            n_transitions_l: 0,
            n_transitions_u: DEFAULT_COUNT_UPPER_BOUND,
            ion_ratio_pair_name_1: String::new(),
            ion_ratio_pair_name_2: String::new(),
            ion_ratio_l: 0.0,
            ion_ratio_u: DEFAULT_UPPER_BOUND,
            ion_ratio_feature_name: String::new(),
            meta_value_qc: BTreeMap::new(),
        }
    }
}

/// Quality controls between two component groups, for parameters that involve
/// more than one of them, such as a minimum retention time difference.
///
/// Source `MRMFeatureQC::ComponentGroupPairQCs`. Its numeric members have no
/// initialiser in the source, so a default-constructed value holds
/// indeterminate numbers; [`Default`] here zeroes them. The source declares no
/// `operator==`; the derived one compares every member.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ComponentGroupPairQCs {
    /// Name of the component group.
    pub component_group_name: String,
    /// Name of the component group to calculate the resolution or retention
    /// time difference against.
    pub resolution_pair_name: String,
    /// Resolution lower bound.
    pub resolution_l: f64,
    /// Resolution upper bound.
    pub resolution_u: f64,
    /// Retention time difference lower bound.
    pub rt_diff_l: f64,
    /// Retention time difference upper bound.
    pub rt_diff_u: f64,
}

/// The QC criteria of a targeted experiment: bounds per component, per
/// component group and per pair of component groups.
///
/// Source `MRMFeatureQC`: a plain aggregate of three lists with no behaviour of
/// its own (the class declares only its default constructor and destructor).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MRMFeatureQC {
    /// The component QCs.
    pub component_qcs: Vec<ComponentQCs>,
    /// The component group QCs.
    pub component_group_qcs: Vec<ComponentGroupQCs>,
    /// The component group pair QCs.
    pub component_group_pair_qcs: Vec<ComponentGroupPairQCs>,
}

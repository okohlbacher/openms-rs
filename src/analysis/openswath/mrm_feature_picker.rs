// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Per-component and per-component-group parameters for the transition group picker.
//!
//! Ports `ANALYSIS/OPENSWATH/MRMFeaturePicker.h` (the `.cpp` is empty): the two
//! parameter records that
//! [`MRMTransitionGroupPicker`](crate::analysis::openswath::mrm_transition_group_picker::MRMTransitionGroupPicker)
//! consumers load from a file with
//! [`MRMFeaturePickerFile`](crate::format::mrm_feature_picker_file::MRMFeaturePickerFile).
//! Examples of the parameter names are `TransitionGroupPicker:compute_peak_quality`,
//! `TransitionGroupPicker:stop_after_feature`,
//! `TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise` and
//! `TransitionGroupPicker:PeakPickerChromatogram:sn_bin_count`.
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/MRM_FEATURE_QC_SUPPORT.md).

use crate::param::Param;

/// A single component with its parameters (source
/// `MRMFeaturePicker::ComponentParams`).
///
/// Neither name is empty when produced by
/// [`MRMFeaturePickerFile`](crate::format::mrm_feature_picker_file::MRMFeaturePickerFile).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ComponentParams {
    /// The name of the component; not an empty string.
    pub component_name: String,
    /// The name of the component group the component belongs to; not an empty string.
    pub component_group_name: String,
    /// The parameters pertaining to a single component.
    pub params: Param,
}

/// A component group with its parameters (source
/// `MRMFeaturePicker::ComponentGroupParams`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ComponentGroupParams {
    /// The name of the component group; not an empty string.
    pub component_group_name: String,
    /// The parameters pertaining to a component group.
    pub params: Param,
}

/// The source class `MRMFeaturePicker`: an empty class that only scopes the
/// two parameter records, which are free types in this module.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MRMFeaturePicker;

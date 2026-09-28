// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Test doubles for the [`itransition`](super::itransition) access traits.
//!
//! Ports `OPENSWATHALGO/DATAACCESS/MockObjects.h` with `MockObjects.cpp`. Each
//! mock returns exactly what its public fields hold. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md).
//!
//! The source members carry an `m_` prefix (`m_rt_vec`, `m_intensity`, ...);
//! the port drops it.

use super::itransition::{
    FeatureAccess, FeatureAccessPtr, MRMFeatureAccess, SignalToNoiseAccess, TransitionGroupAccess,
};
use crate::Result;
use std::collections::BTreeMap;
use std::sync::Arc;

/// Mock implementing [`FeatureAccess`]. Source `MockFeature`.
///
/// The source constructor leaves the scalar members uninitialised; the port
/// defaults them to 0.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MockFeature {
    /// Returned by [`FeatureAccess::rt_values`].
    pub rt_vec: Vec<f64>,
    /// Returned by [`FeatureAccess::intensity_values`].
    pub intensity_vec: Vec<f64>,
    /// Returned by [`FeatureAccess::intensity`].
    pub intensity: f32,
    /// Returned by [`FeatureAccess::rt`].
    pub rt: f64,
}

impl FeatureAccess for MockFeature {
    fn rt_values(&self) -> Vec<f64> {
        self.rt_vec.clone()
    }
    fn intensity_values(&self) -> Vec<f64> {
        self.intensity_vec.clone()
    }
    fn intensity(&self) -> f32 {
        self.intensity
    }
    fn rt(&self) -> f64 {
        self.rt
    }
}

/// Mock implementing [`MRMFeatureAccess`]. Source `MockMRMFeature`.
///
/// Native IDs come back in the sorted order of the maps, as in the source
/// (`std::map`).
#[derive(Clone, Debug, Default)]
pub struct MockMRMFeature {
    /// Transition features by native ID.
    pub features: BTreeMap<String, Arc<MockFeature>>,
    /// Precursor features by native ID.
    pub precursor_features: BTreeMap<String, Arc<MockFeature>>,
    /// Returned by [`MRMFeatureAccess::intensity`].
    pub intensity: f32,
    /// Returned by [`MRMFeatureAccess::rt`].
    pub rt: f64,
    /// Returned by [`MRMFeatureAccess::meta_value`] for every name.
    pub metavalue: f64,
}

impl MRMFeatureAccess for MockMRMFeature {
    /// The feature, or `None`. The source looks the ID up with
    /// `std::map::operator[]`, which inserts and returns a null pointer for an
    /// unknown ID; the port neither inserts nor needs a null.
    fn feature(&self, native_id: &str) -> Option<FeatureAccessPtr> {
        self.features
            .get(native_id)
            .map(|f| Arc::clone(f) as FeatureAccessPtr)
    }
    fn precursor_feature(&self, native_id: &str) -> Option<FeatureAccessPtr> {
        self.precursor_features
            .get(native_id)
            .map(|f| Arc::clone(f) as FeatureAccessPtr)
    }
    fn native_ids(&self) -> Vec<String> {
        self.features.keys().cloned().collect()
    }
    fn precursor_ids(&self) -> Vec<String> {
        self.precursor_features.keys().cloned().collect()
    }
    fn intensity(&self) -> f32 {
        self.intensity
    }
    fn rt(&self) -> f64 {
        self.rt
    }
    /// [`Self::metavalue`] whatever `name` is; never fails.
    fn meta_value(&self, _name: &str) -> Result<f64> {
        Ok(self.metavalue)
    }
    fn len(&self) -> usize {
        self.features.len()
    }
}

/// Mock implementing [`TransitionGroupAccess`]. Source `MockTransitionGroup`.
///
/// [`TransitionGroupAccess::len`] returns [`Self::size`], which is not tied to
/// the length of [`Self::native_ids`], as in the source.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MockTransitionGroup {
    /// Returned by [`TransitionGroupAccess::len`].
    pub size: usize,
    /// Returned by [`TransitionGroupAccess::native_ids`].
    pub native_ids: Vec<String>,
    /// Returned by [`TransitionGroupAccess::library_intensities`].
    pub library_intensities: Vec<f64>,
}

impl TransitionGroupAccess for MockTransitionGroup {
    fn len(&self) -> usize {
        self.size
    }
    fn native_ids(&self) -> Vec<String> {
        self.native_ids.clone()
    }
    fn library_intensities(&self) -> Vec<f64> {
        self.library_intensities.clone()
    }
}

/// Mock implementing [`SignalToNoiseAccess`]. Source `MockSignalToNoise`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MockSignalToNoise {
    /// Returned by [`SignalToNoiseAccess::value_at_rt`] for every RT.
    pub sn_value: f64,
}

impl SignalToNoiseAccess for MockSignalToNoise {
    fn value_at_rt(&self, _rt: f64) -> f64 {
        self.sn_value
    }
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The access interfaces OpenSWATH scoring reads features, peak groups,
//! transition groups and signal-to-noise estimates through.
//!
//! Ports `OPENSWATHALGO/DATAACCESS/ITransition.h` (header-only). The source's
//! four abstract classes become four traits; their `I` prefix becomes an
//! `Access` suffix so the traits do not shadow the kernel types
//! ([`Feature`](crate::kernel::Feature), [`MRMFeature`](crate::kernel::mrm::MRMFeature))
//! they are implemented over:
//!
//! | Source | Port |
//! |---|---|
//! | `IFeature` (and `std::shared_ptr<IFeature>`, now a borrow) | [`FeatureAccess`] |
//! | `IMRMFeature` | [`MRMFeatureAccess`] |
//! | `ITransitionGroup` | [`TransitionGroupAccess`] |
//! | `ISignalToNoise`, `ISignalToNoisePtr` | [`SignalToNoiseAccess`], [`SignalToNoisePtr`] |
//!
//! Implementations over OpenMS types live in
//! [`mrm_feature_access`](super::mrm_feature_access); test doubles in
//! [`mock_objects`](super::mock_objects). See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md).
//!
//! Out-parameters (`getRT(std::vector<double>&)`) become returned vectors, and
//! the non-`const` getters of the source take `&self`: no implementation in the
//! source or the port mutates state to answer them (the source's own comment on
//! `ISignalToNoise::getValueAtRT` says it is non-`const` only "due to OpenMS
//! implementation", which the port's implementation does not need).

use crate::Result;
use std::sync::Arc;

/// Shared handle to a [`SignalToNoiseAccess`], as the source `ISignalToNoisePtr`.
pub type SignalToNoisePtr = Arc<dyn SignalToNoiseAccess>;

/// A single-transition feature: an elution profile with its apex RT and
/// intensity. Source `IFeature`.
pub trait FeatureAccess: Send + Sync {
    /// Retention times of the feature's points (source `getRT(std::vector<double>&)`).
    ///
    /// # Errors
    ///
    /// Implementation-specific; the OpenMS implementation fails when the
    /// feature does not have exactly one convex hull.
    fn rt_values(&self) -> Result<Vec<f64>>;
    /// Intensities of the feature's points, parallel to [`Self::rt_values`]
    /// (source `getIntensity(std::vector<double>&)`).
    ///
    /// # Errors
    ///
    /// As [`Self::rt_values`].
    fn intensity_values(&self) -> Result<Vec<f64>>;
    /// The feature intensity.
    fn intensity(&self) -> f32;
    /// The feature (apex) retention time.
    fn rt(&self) -> f64;
}

/// A peak group across several transitions (and optional precursors). Source
/// `IMRMFeature`.
pub trait MRMFeatureAccess: Send + Sync {
    /// The feature of the transition with native ID `native_id`, or `None`.
    ///
    /// The source returns a `std::shared_ptr<IFeature>`; the port borrows the
    /// feature from the peak group, which every source caller only uses
    /// locally. For an unknown ID the OpenMS implementation throws
    /// `Exception::ElementNotFound` and the mock returns a null pointer (after
    /// inserting it into its map); both become `None`.
    fn feature(&self, native_id: &str) -> Option<&dyn FeatureAccess>;
    /// The precursor feature with native ID `native_id`, or `None`; see
    /// [`Self::feature`].
    fn precursor_feature(&self, native_id: &str) -> Option<&dyn FeatureAccess>;
    /// Native IDs of the transition features.
    fn native_ids(&self) -> Vec<String>;
    /// Native IDs of the precursor features.
    fn precursor_ids(&self) -> Vec<String>;
    /// The peak-group intensity.
    fn intensity(&self) -> f32;
    /// The peak-group retention time.
    fn rt(&self) -> f64;
    /// The numeric meta value `name` of the peak group.
    ///
    /// # Errors
    ///
    /// Implementation-specific; the OpenMS implementation fails when the value
    /// is missing or not numeric, where the source throws on the conversion.
    fn meta_value(&self, name: &str) -> Result<f64>;
    /// Number of transition features.
    fn len(&self) -> usize;
    /// Whether there are no transition features.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// A transition group: the transitions of one assay. Source `ITransitionGroup`.
pub trait TransitionGroupAccess: Send + Sync {
    /// Number of transitions.
    fn len(&self) -> usize;
    /// Whether there are no transitions.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Native IDs of the transitions.
    fn native_ids(&self) -> Vec<String>;
    /// Library intensities of the transitions (source
    /// `getLibraryIntensities(std::vector<double>&)`).
    ///
    /// # Errors
    ///
    /// Implementation-specific; the OpenMS implementation reports the resource
    /// ceiling of [`MRMTransitionGroup::library_intensity`](crate::kernel::mrm::MRMTransitionGroup::library_intensity).
    fn library_intensities(&self) -> Result<Vec<f64>>;
}

/// A signal-to-noise estimate along retention time. Source `ISignalToNoise`.
pub trait SignalToNoiseAccess: Send + Sync {
    /// The signal-to-noise value at retention time `rt`.
    fn value_at_rt(&self, rt: f64) -> f64;
}

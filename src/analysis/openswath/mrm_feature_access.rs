// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Implementations of the OpenSWATH access traits over OpenMS kernel types.
//!
//! Ports `ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h` with
//! `MRMFeatureAccessOpenMS.cpp`: minimal implementations of the
//! [`itransition`](crate::analysis::openswath::itransition) interfaces.
//!
//! | Source | Port | Implements |
//! |---|---|---|
//! | `FeatureOpenMS` | [`FeatureOpenMS`](crate::analysis::openswath::mrm_feature_access::FeatureOpenMS) | [`FeatureAccess`](crate::analysis::openswath::itransition::FeatureAccess) |
//! | `MRMFeatureOpenMS` | [`MRMFeatureOpenMS`](crate::analysis::openswath::mrm_feature_access::MRMFeatureOpenMS) | [`MRMFeatureAccess`](crate::analysis::openswath::itransition::MRMFeatureAccess) |
//! | `TransitionGroupOpenMS<SpectrumT, TransitionT>` | [`TransitionGroupOpenMS`](crate::analysis::openswath::mrm_feature_access::TransitionGroupOpenMS) | [`TransitionGroupAccess`](crate::analysis::openswath::itransition::TransitionGroupAccess) |
//! | `SignalToNoiseOpenMS<ContainerT>` | [`SignalToNoiseOpenMS`](crate::analysis::openswath::mrm_feature_access::SignalToNoiseOpenMS) | [`SignalToNoiseAccess`](crate::analysis::openswath::itransition::SignalToNoiseAccess) |
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md).
//!
//! The wrappers borrow the kernel records they wrap, as the source's hold
//! references and raw pointers; the borrow checker now enforces the lifetime
//! the source leaves to the caller. The source's file-scope "default instances"
//! (template instantiations kept for the linker) have no counterpart.

use super::itransition::{
    FeatureAccess, MRMFeatureAccess, SignalToNoiseAccess, TransitionGroupAccess,
};
use crate::kernel::features::Feature;
use crate::kernel::mrm::{MRMFeature, MRMTransitionGroup, NativeIdentified, Transition};
use crate::kernel::{MSChromatogram, MSSpectrum};
use crate::processing::peak_picking::SignalToNoiseEstimatorMedian;
use crate::{Error, Result};

fn single_hull_points(feature: &Feature) -> Result<Vec<crate::kernel::Point2D>> {
    if feature.convex_hulls.len() != 1 {
        return Err(Error::InvalidValue(format!(
            "there needs to be exactly one convex hull per feature, found {}",
            feature.convex_hulls.len()
        )));
    }
    Ok(feature.convex_hulls[0].hull_points())
}

/// An implementation of the OpenSWATH feature access interface over an OpenMS
/// [`Feature`]: the elution profile is the feature's single convex hull, read
/// as `(RT, intensity)` points.
#[derive(Clone, Copy, Debug)]
pub struct FeatureOpenMS<'a> {
    feature: &'a Feature,
}

impl<'a> FeatureOpenMS<'a> {
    /// Wrap `feature`.
    pub fn new(feature: &'a Feature) -> Self {
        Self { feature }
    }

    /// The wrapped feature.
    pub fn feature(&self) -> &'a Feature {
        self.feature
    }
}

impl FeatureAccess for FeatureOpenMS<'_> {
    /// The x coordinates (RT) of the hull points.
    ///
    /// The source appends to the caller's vector; the port returns a fresh one.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] unless the feature has exactly one convex hull
    /// (source `Exception::IllegalArgument`).
    fn rt_values(&self) -> Result<Vec<f64>> {
        Ok(single_hull_points(self.feature)?
            .iter()
            .map(|p| p.rt)
            .collect())
    }

    /// The y coordinates (intensity) of the hull points; errors as
    /// [`FeatureAccess::rt_values`].
    fn intensity_values(&self) -> Result<Vec<f64>> {
        Ok(single_hull_points(self.feature)?
            .iter()
            .map(|p| p.mz)
            .collect())
    }

    fn intensity(&self) -> f32 {
        self.feature.intensity
    }

    fn rt(&self) -> f64 {
        self.feature.rt
    }
}

/// One set of sub-features of a peak group with its native IDs (the source's
/// private `FeatureLane_`). `features[i]` belongs to `ids[i]`.
#[derive(Clone, Debug)]
struct FeatureLane<'a> {
    features: Vec<FeatureOpenMS<'a>>,
    ids: Vec<String>,
    /// Whether `ids` was supplied by the caller as aligned with the storage
    /// order, which enables the index-hint fast paths.
    aligned: bool,
}

impl<'a> FeatureLane<'a> {
    fn from_lookup(ids: Vec<String>, lookup: impl Fn(&str) -> Result<&'a Feature>) -> Result<Self> {
        let features = ids
            .iter()
            .map(|id| lookup(id).map(FeatureOpenMS::new))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            features,
            ids,
            aligned: false,
        })
    }

    fn find_index(&self, native_id: &str) -> Option<usize> {
        self.ids.iter().position(|id| id == native_id)
    }

    fn not_found(native_id: &str) -> Error {
        Error::MissingInformation(format!("no feature with native ID '{native_id}'"))
    }

    fn intensities_by_ids(&self, native_ids: &[String]) -> Result<Vec<Vec<f64>>> {
        if self.aligned && native_ids == self.ids.as_slice() {
            return self.features.iter().map(|f| f.intensity_values()).collect();
        }
        native_ids
            .iter()
            .map(|id| {
                let index = self.find_index(id).ok_or_else(|| Self::not_found(id))?;
                self.features[index].intensity_values()
            })
            .collect()
    }

    fn intensity(&self, native_id: &str, expected_index: usize) -> Result<f32> {
        if self.aligned && self.ids.get(expected_index).map(String::as_str) == Some(native_id) {
            return Ok(self.features[expected_index].intensity());
        }
        let index = self
            .find_index(native_id)
            .ok_or_else(|| Self::not_found(native_id))?;
        Ok(self.features[index].intensity())
    }
}

/// An implementation of the OpenSWATH MRM feature access interface over an
/// OpenMS [`MRMFeature`] (peak group).
///
/// Fragment (transition) features and precursor features are held in two
/// independent lanes, each a list of wrapped features with their native IDs.
/// Native-ID lookups are linear and return the first match.
///
/// The source builds the feature wrappers lazily on the first `getFeature`;
/// the port builds them in the constructor (they are two words each).
#[derive(Clone, Debug)]
pub struct MRMFeatureOpenMS<'a> {
    mrmfeature: &'a MRMFeature,
    feature_lane: FeatureLane<'a>,
    precursor_feature_lane: FeatureLane<'a>,
}

impl<'a> MRMFeatureOpenMS<'a> {
    /// Build both lanes from the native-ID maps stored in `mrmfeature`, in the
    /// maps' (lexicographic) key order.
    pub fn new(mrmfeature: &'a MRMFeature) -> Self {
        let lane = |ids: Vec<String>, features: Vec<FeatureOpenMS<'a>>| FeatureLane {
            features,
            ids,
            aligned: false,
        };
        let ids: Vec<String> = mrmfeature.feature_ids().map(str::to_string).collect();
        let features = ids
            .iter()
            .filter_map(|id| mrmfeature.feature(id).ok().map(FeatureOpenMS::new))
            .collect();
        let precursor_ids: Vec<String> = mrmfeature
            .precursor_feature_ids()
            .map(str::to_string)
            .collect();
        let precursor_features = precursor_ids
            .iter()
            .filter_map(|id| {
                mrmfeature
                    .precursor_feature(id)
                    .ok()
                    .map(FeatureOpenMS::new)
            })
            .collect();
        Self {
            mrmfeature,
            feature_lane: lane(ids, features),
            precursor_feature_lane: lane(precursor_ids, precursor_features),
        }
    }

    /// Build both lanes from explicit fragment and precursor native-ID lists,
    /// looking each ID up in `mrmfeature`. The lists define the order of
    /// [`MRMFeatureAccess::native_ids`] and [`MRMFeatureAccess::precursor_ids`].
    ///
    /// # Errors
    ///
    /// The error of [`MRMFeature::feature`] or [`MRMFeature::precursor_feature`]
    /// for an ID the peak group does not hold (the source's lookup throws too).
    pub fn with_ids(
        mrmfeature: &'a MRMFeature,
        feature_ids: &[String],
        precursor_feature_ids: &[String],
    ) -> Result<Self> {
        Ok(Self {
            mrmfeature,
            feature_lane: FeatureLane::from_lookup(feature_ids.to_vec(), |id| {
                mrmfeature.feature(id)
            })?,
            precursor_feature_lane: FeatureLane::from_lookup(
                precursor_feature_ids.to_vec(),
                |id| mrmfeature.precursor_feature(id),
            )?,
        })
    }

    /// Build both lanes with an externally aligned native-ID order.
    ///
    /// `feature_ids` and `precursor_feature_ids` define the aligned order used
    /// by the index-hint fast paths ([`Self::feature_intensity`],
    /// [`Self::feature_intensities`]). Fragment features are read directly
    /// from the storage of `mrmfeature` ([`MRMFeature::features`]), so the
    /// caller must ensure that `feature_ids` matches that storage order; the
    /// lookup-ID lists are only used to validate sizes and, for precursors, to
    /// find the features.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] when `feature_ids` and `feature_lookup_ids`, or
    ///   `precursor_feature_ids` and `precursor_feature_lookup_ids`, differ in
    ///   length, or when the fragment storage is not as long as
    ///   `feature_lookup_ids` (the source's three `IllegalArgument`s);
    /// - the error of [`MRMFeature::precursor_feature`] for an unknown
    ///   precursor lookup ID.
    pub fn with_aligned_ids(
        mrmfeature: &'a MRMFeature,
        feature_ids: &[String],
        precursor_feature_ids: &[String],
        feature_lookup_ids: &[String],
        precursor_feature_lookup_ids: &[String],
    ) -> Result<Self> {
        if feature_ids.len() != feature_lookup_ids.len() {
            return Err(Error::InvalidValue(
                "feature id cache needs to match feature lookup id cache".into(),
            ));
        }
        if precursor_feature_ids.len() != precursor_feature_lookup_ids.len() {
            return Err(Error::InvalidValue(
                "precursor feature id cache needs to match precursor feature lookup id cache"
                    .into(),
            ));
        }
        let storage = mrmfeature.features();
        if storage.len() != feature_lookup_ids.len() {
            return Err(Error::InvalidValue(
                "feature storage needs to match feature lookup id cache".into(),
            ));
        }
        let feature_lane = FeatureLane {
            features: storage.iter().map(FeatureOpenMS::new).collect(),
            ids: feature_ids.to_vec(),
            aligned: true,
        };
        let mut precursor_feature_lane =
            FeatureLane::from_lookup(precursor_feature_lookup_ids.to_vec(), |id| {
                mrmfeature.precursor_feature(id)
            })?;
        precursor_feature_lane.ids = precursor_feature_ids.to_vec();
        precursor_feature_lane.aligned = true;
        Ok(Self {
            mrmfeature,
            feature_lane,
            precursor_feature_lane,
        })
    }

    /// The wrapped peak group.
    pub fn mrm_feature(&self) -> &'a MRMFeature {
        self.mrmfeature
    }

    /// The hull intensities of the fragment features named by `native_ids`, one
    /// vector per ID. Source `getFeatureIntensities`.
    ///
    /// When the lanes were built by [`Self::with_aligned_ids`] and `native_ids`
    /// equals the aligned ID list, features are taken by position; otherwise
    /// each ID is looked up.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] for an unknown ID (source
    /// `Exception::ElementNotFound`), and the errors of
    /// [`FeatureAccess::intensity_values`].
    pub fn feature_intensities(&self, native_ids: &[String]) -> Result<Vec<Vec<f64>>> {
        self.feature_lane.intensities_by_ids(native_ids)
    }

    /// The hull intensities of the precursor features named by `native_ids`;
    /// rules and errors as [`Self::feature_intensities`]. Source
    /// `getPrecursorFeatureIntensities`.
    ///
    /// # Errors
    ///
    /// As [`Self::feature_intensities`].
    pub fn precursor_feature_intensities(&self, native_ids: &[String]) -> Result<Vec<Vec<f64>>> {
        self.precursor_feature_lane.intensities_by_ids(native_ids)
    }

    /// The intensity of the fragment feature `native_id`.
    ///
    /// `expected_index` is an optional fast-path hint, trusted only when the
    /// lanes were built by [`Self::with_aligned_ids`] and the aligned ID at
    /// that index is `native_id`; the default native-ID order is map-based and
    /// is never assumed to match the storage order. Otherwise the ID is looked
    /// up. Source `getFeatureIntensity`.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] for an unknown ID.
    pub fn feature_intensity(&self, native_id: &str, expected_index: usize) -> Result<f32> {
        self.feature_lane.intensity(native_id, expected_index)
    }
}

impl MRMFeatureAccess for MRMFeatureOpenMS<'_> {
    fn feature(&self, native_id: &str) -> Option<&dyn FeatureAccess> {
        let index = self.feature_lane.find_index(native_id)?;
        Some(&self.feature_lane.features[index] as &dyn FeatureAccess)
    }

    fn precursor_feature(&self, native_id: &str) -> Option<&dyn FeatureAccess> {
        let index = self.precursor_feature_lane.find_index(native_id)?;
        Some(&self.precursor_feature_lane.features[index] as &dyn FeatureAccess)
    }

    fn native_ids(&self) -> Vec<String> {
        self.feature_lane.ids.clone()
    }

    fn precursor_ids(&self) -> Vec<String> {
        self.precursor_feature_lane.ids.clone()
    }

    fn intensity(&self) -> f32 {
        self.mrmfeature.intensity
    }

    fn rt(&self) -> f64 {
        self.mrmfeature.rt
    }

    /// The peak group's meta value `name` as a number.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when the value is absent and
    /// [`Error::InvalidValue`] when it is not numeric; the source converts the
    /// `DataValue` to `double`, which throws in both cases.
    fn meta_value(&self, name: &str) -> Result<f64> {
        self.mrmfeature
            .metadata
            .get(name)
            .ok_or_else(|| Error::MissingInformation(format!("no meta value '{name}'")))?
            .as_f64()
    }

    /// Number of fragment features in the lane (the aligned ID count for
    /// [`MRMFeatureOpenMS::with_aligned_ids`]).
    fn len(&self) -> usize {
        self.feature_lane.ids.len()
    }
}

/// An implementation of the OpenSWATH transition group access interface over
/// an OpenMS [`MRMTransitionGroup`].
#[derive(Clone, Copy, Debug)]
pub struct TransitionGroupOpenMS<'a, C, T> {
    trgroup: &'a MRMTransitionGroup<C, T>,
}

impl<'a, C, T> TransitionGroupOpenMS<'a, C, T> {
    /// Wrap `trgroup`.
    pub fn new(trgroup: &'a MRMTransitionGroup<C, T>) -> Self {
        Self { trgroup }
    }
}

impl<C, T> TransitionGroupAccess for TransitionGroupOpenMS<'_, C, T>
where
    C: NativeIdentified + Sync,
    T: Transition + Sync,
{
    /// The number of fragment-ion chromatograms (source `MRMTransitionGroup::size`).
    fn len(&self) -> usize {
        self.trgroup.size()
    }

    /// The native IDs of the fragment-ion chromatograms, in storage order.
    fn native_ids(&self) -> Vec<String> {
        self.trgroup
            .chromatograms()
            .iter()
            .map(|c| c.native_id().to_string())
            .collect()
    }

    /// The transitions' library intensities, negatives clamped to zero; see
    /// [`MRMTransitionGroup::library_intensity`].
    fn library_intensities(&self) -> Result<Vec<f64>> {
        self.trgroup.library_intensity()
    }
}

/// An implementation of the OpenSWATH signal-to-noise access interface: a
/// median signal-to-noise estimate over a spectrum or chromatogram, read at the
/// point nearest to a retention time.
///
/// The source class is a template over the container and keeps a reference to
/// it; the port copies the positions and computes the estimates once, in the
/// constructor, as the source's constructor does (`sn_.init(chromatogram_)`).
#[derive(Clone, Debug, PartialEq)]
pub struct SignalToNoiseOpenMS {
    positions: Vec<f64>,
    signal_to_noise: Vec<f64>,
}

impl SignalToNoiseOpenMS {
    fn estimator(
        sn_win_len: f64,
        sn_bin_count: usize,
        write_log_messages: bool,
    ) -> SignalToNoiseEstimatorMedian {
        SignalToNoiseEstimatorMedian {
            window_length: sn_win_len,
            bin_count: sn_bin_count,
            write_log_messages,
            ..SignalToNoiseEstimatorMedian::default()
        }
    }

    /// Estimate over the points of `positions`/`intensities` (sorted by
    /// position), with the estimator defaults except the window length
    /// `sn_win_len`, the bin count `sn_bin_count` and `write_log_messages`.
    ///
    /// # Errors
    ///
    /// The errors of [`SignalToNoiseEstimatorMedian::estimate`] (unequal
    /// lengths, unsorted positions, invalid parameters, resource ceilings).
    pub fn from_points(
        positions: Vec<f64>,
        intensities: &[f64],
        sn_win_len: f64,
        sn_bin_count: usize,
        write_log_messages: bool,
    ) -> Result<Self> {
        let estimates = Self::estimator(sn_win_len, sn_bin_count, write_log_messages)
            .estimate(&positions, intensities)?;
        Ok(Self {
            positions,
            signal_to_noise: estimates.signal_to_noise,
        })
    }

    /// Estimate over a chromatogram; see [`Self::from_points`].
    ///
    /// # Errors
    ///
    /// As [`Self::from_points`].
    pub fn from_chromatogram(
        chromatogram: &MSChromatogram,
        sn_win_len: f64,
        sn_bin_count: usize,
        write_log_messages: bool,
    ) -> Result<Self> {
        let positions = chromatogram.peaks.iter().map(|p| p.rt).collect();
        let intensities: Vec<f64> = chromatogram
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect();
        Self::from_points(
            positions,
            &intensities,
            sn_win_len,
            sn_bin_count,
            write_log_messages,
        )
    }

    /// Estimate over a spectrum (positions are m/z); see [`Self::from_points`].
    ///
    /// # Errors
    ///
    /// As [`Self::from_points`].
    pub fn from_spectrum(
        spectrum: &MSSpectrum,
        sn_win_len: f64,
        sn_bin_count: usize,
        write_log_messages: bool,
    ) -> Result<Self> {
        let positions = spectrum.peaks.iter().map(|p| p.mz).collect();
        let intensities: Vec<f64> = spectrum
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect();
        Self::from_points(
            positions,
            &intensities,
            sn_win_len,
            sn_bin_count,
            write_log_messages,
        )
    }
}

impl SignalToNoiseAccess for SignalToNoiseOpenMS {
    /// The estimate at the point nearest to `rt`; -1 for an empty container.
    ///
    /// The candidate is the first point with a position greater than `rt` (the
    /// last point when there is none) and the point before it; the one before
    /// wins only when it is strictly closer, as in the source.
    fn value_at_rt(&self, rt: f64) -> f64 {
        let n = self.positions.len();
        if n == 0 {
            return -1.0;
        }
        let mut iter = self.positions.partition_point(|&p| p <= rt);
        if iter == n {
            iter -= 1;
        }
        let prev = iter.saturating_sub(1);
        if (self.positions[prev] - rt).abs() < (self.positions[iter] - rt).abs() {
            self.signal_to_noise[prev]
        } else {
            self.signal_to_noise[iter]
        }
    }
}

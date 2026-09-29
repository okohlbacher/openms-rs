// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! All peaks of an experiment that passed the filters for one peak pattern.
//!
//! Native equivalent of `FEATUREFINDER/MultiplexFilteredMSExperiment.h`. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`.

use super::filtered_peak::MultiplexFilteredPeak;
use crate::concept::unique_id::UniqueIdGenerator;
use crate::kernel::{
    BaseFeature, ColumnHeader, ConsensusFeature, ConsensusMap, FeatureHandle, MSExperiment,
};
use crate::{Error, Result};

/// The peaks (and optionally their raw data points) of an experiment that
/// correspond to one specific peak pattern (source
/// `MultiplexFilteredMSExperiment`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MultiplexFilteredMSExperiment {
    /// Peaks which passed the peak pattern filter (source `result_`).
    result: Vec<MultiplexFilteredPeak>,
}

impl MultiplexFilteredMSExperiment {
    /// An empty result (source constructor).
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a peak to the results (source `addPeak`).
    pub fn add_peak(&mut self, peak: MultiplexFilteredPeak) {
        self.result.push(peak);
    }

    /// Peak `i`, or `None` past the end (source `getPeak`, which returns a
    /// copy and is undefined past the end).
    pub fn peak(&self, i: usize) -> Option<&MultiplexFilteredPeak> {
        self.result.get(i)
    }

    /// All peaks in insertion order.
    pub fn peaks(&self) -> &[MultiplexFilteredPeak] {
        &self.result
    }

    /// m/z of peak `i`, or `None` past the end (source `getMZ(size_t)`).
    pub fn mz_at(&self, i: usize) -> Option<f64> {
        self.result.get(i).map(MultiplexFilteredPeak::mz)
    }

    /// m/z of all peaks (source `getMZ()`).
    pub fn mz(&self) -> Vec<f64> {
        self.result.iter().map(MultiplexFilteredPeak::mz).collect()
    }

    /// RT of peak `i` widened to `f64`, or `None` past the end (source
    /// `getRT(size_t)`).
    pub fn rt_at(&self, i: usize) -> Option<f64> {
        self.result.get(i).map(|p| f64::from(p.rt()))
    }

    /// RT of all peaks widened to `f64` (source `getRT()`).
    pub fn rt(&self) -> Vec<f64> {
        self.result.iter().map(|p| f64::from(p.rt())).collect()
    }

    /// Number of peaks (source `size`).
    pub fn len(&self) -> usize {
        self.result.len()
    }

    /// Whether no peak passed.
    pub fn is_empty(&self) -> bool {
        self.result.is_empty()
    }

    /// The consensus map source `writeDebugOutput` stores: one consensus
    /// feature per peak at the peak's RT and m/z (intensity 1, charge 1,
    /// quality 1) with one handle per satellite, in satellite order.
    ///
    /// Handle `k` of a feature lies in map `k`, at the satellite peak's RT and
    /// m/z in `exp_picked`, intensity 1 and "charge" the satellite's pattern
    /// index, as the source abuses the charge field. Column header `k` counts
    /// the handles in map `k` and is labelled `satellite_k` with file name
    /// `satellites`. The features are sorted by position; then the map and each
    /// feature in turn receive a unique ID from `ids`, where the source draws
    /// from its process-wide generator; the experiment type is `label-free`.
    ///
    /// `exp_picked` is the original (not white) centroided experiment the
    /// satellite indices refer to.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a satellite index lies outside
    /// `exp_picked` (undefined behaviour in the source) or a pattern index does
    /// not fit a charge, and the errors of the consensus map operations.
    pub fn debug_consensus_map(
        &self,
        exp_picked: &MSExperiment,
        ids: &mut UniqueIdGenerator,
    ) -> Result<ConsensusMap> {
        let mut map = ConsensusMap::new();
        for peak in &self.result {
            let mut base = BaseFeature::new(f64::from(peak.rt()), peak.mz(), 1.0);
            base.charge = 1;
            base.quality = 1.0;
            let mut consensus = ConsensusFeature::from(base);
            for (count, (pattern_idx, satellite)) in peak.satellite_entries().enumerate() {
                let spectrum = exp_picked.spectra.get(satellite.rt_idx).ok_or_else(|| {
                    Error::InvalidValue("satellite spectrum index outside the experiment".into())
                })?;
                let point = spectrum.peaks.get(satellite.mz_idx).ok_or_else(|| {
                    Error::InvalidValue("satellite peak index outside the spectrum".into())
                })?;
                let mut handle_base = BaseFeature::new(spectrum.rt, point.mz, 1.0);
                handle_base.charge = i32::try_from(pattern_idx).map_err(|_| {
                    Error::InvalidValue("pattern index does not fit a charge".into())
                })?;
                let map_index = count as u64;
                consensus.insert(FeatureHandle::new(map_index, &handle_base))?;
                let header = map
                    .column_headers
                    .entry(map_index)
                    .or_insert_with(ColumnHeader::default);
                header.size += 1;
                header.label = format!("satellite_{count}");
                header.filename = "satellites".into();
            }
            map.features.push(consensus);
        }
        map.sort_by_position()?;
        // `applyMemberFunction(&UniqueIdInterface::setUniqueId)` assigns the
        // map's own ID first, then each feature's in order.
        map.unique_id = ids.get_unique_id();
        for feature in &mut map.features {
            feature.unique_id = ids.get_unique_id();
        }
        map.experiment_type = "label-free".into();
        Ok(map)
    }

    /// Write all peaks to a consensusXML file (source `writeDebugOutput`), the
    /// map of [`Self::debug_consensus_map`].
    ///
    /// # Errors
    ///
    /// As [`Self::debug_consensus_map`], and the errors of
    /// [`crate::format::consensusxml::store`].
    #[cfg(feature = "consensusxml")]
    pub fn write_debug_output(
        &self,
        exp_picked: &MSExperiment,
        debug_out: impl AsRef<std::path::Path>,
        ids: &mut UniqueIdGenerator,
    ) -> Result<()> {
        let map = self.debug_consensus_map(exp_picked, ids)?;
        crate::format::consensusxml::store(debug_out, &map)
    }
}

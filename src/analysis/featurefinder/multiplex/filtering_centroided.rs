// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Filtering centroided data for multiplex peak patterns.
//!
//! Native equivalent of `FEATUREFINDER/MultiplexFilteringCentroided.h`. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`.

use std::ops::{Deref, DerefMut};

use super::filtered_ms_experiment::MultiplexFilteredMSExperiment;
use super::filtered_peak::MultiplexFilteredPeak;
use super::filtering::{MultiplexFiltering, rt_begin, rt_end};
use super::isotopic_peak_pattern::MultiplexIsotopicPeakPattern;
use crate::Result;
use crate::concept::progress_logger::ProgressReporter;
use crate::kernel::MSExperiment;

/// Filters centroided data for peak patterns (source
/// `MultiplexFilteringCentroided`).
///
/// Searches the centroided data for the patterns of multiple peptides and
/// generates one filter result per pattern. The shared state and filters are
/// those of [`MultiplexFiltering`], reachable through [`Deref`] as the
/// source's base class members are.
#[derive(Clone, Debug)]
pub struct MultiplexFilteringCentroided {
    base: MultiplexFiltering,
}

impl Deref for MultiplexFilteringCentroided {
    type Target = MultiplexFiltering;
    fn deref(&self) -> &MultiplexFiltering {
        &self.base
    }
}

impl DerefMut for MultiplexFilteringCentroided {
    fn deref_mut(&mut self) -> &mut MultiplexFiltering {
        &mut self.base
    }
}

impl MultiplexFilteringCentroided {
    /// Prepare filtering (source constructor); the arguments are those of
    /// [`MultiplexFiltering::new`], with `averagine_type` one of `peptide`,
    /// `RNA` and `DNA` (the source's default is `peptide`).
    ///
    /// # Errors
    ///
    /// As [`MultiplexFiltering::new`].
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        exp_centroided: &MSExperiment,
        patterns: Vec<MultiplexIsotopicPeakPattern>,
        isotopes_per_peptide_min: usize,
        isotopes_per_peptide_max: usize,
        intensity_cutoff: f64,
        rt_band: f64,
        mz_tolerance: f64,
        mz_tolerance_unit_in_ppm: bool,
        peptide_similarity: f64,
        averagine_similarity: f64,
        averagine_similarity_scaling: f64,
        averagine_type: &str,
    ) -> Result<Self> {
        Ok(Self {
            base: MultiplexFiltering::new(
                exp_centroided,
                patterns,
                isotopes_per_peptide_min,
                isotopes_per_peptide_max,
                intensity_cutoff,
                rt_band,
                mz_tolerance,
                mz_tolerance_unit_in_ppm,
                peptide_similarity,
                averagine_similarity,
                averagine_similarity_scaling,
                averagine_type,
            )?,
        })
    }

    /// Filter for the patterns: one result per pattern, in pattern order
    /// (source `filter`).
    ///
    /// # Errors
    ///
    /// As [`Self::filter_with_progress`].
    pub fn filter(&mut self) -> Result<Vec<MultiplexFilteredMSExperiment>> {
        self.filter_with_progress(&mut ProgressReporter::silent())
    }

    /// [`Self::filter`], reporting progress with the label "filtering LC-MS
    /// data" over `patterns * spectra` steps, one step per non-empty white
    /// spectrum, as the source does.
    ///
    /// For each pattern the white experiment is rebuilt, then every peak of
    /// every non-empty white spectrum is tried as a primary peak, with the RT
    /// band `[rt - rt_band / 2, rt + rt_band / 2]`, through
    /// [`MultiplexFiltering::filter_peak_positions`],
    /// [`MultiplexFiltering::filter_averagine_model`] and
    /// [`MultiplexFiltering::filter_peptide_correlation`]. A peak passing all
    /// three is added to the result and its mass traces are blacklisted at
    /// once, so later primary peaks of the same pattern already see the new
    /// marks (through the blacklist, not the white experiment, which is only
    /// rebuilt for the next pattern).
    ///
    /// The source runs the peak loop under `#pragma omp parallel for`, with the
    /// blacklist read by every thread while a `critical` section writes it; its
    /// result then depends on thread scheduling. This port runs the loop
    /// serially, which is the source's result with one thread (`-threads 1`,
    /// as the upstream tests run).
    ///
    /// # Errors
    ///
    /// The errors of the three filters and of the blacklisting, and of
    /// `progress`. An error leaves the blacklist part-way updated.
    pub fn filter_with_progress(
        &mut self,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<Vec<MultiplexFilteredMSExperiment>> {
        let steps = self
            .base
            .patterns
            .len()
            .saturating_mul(self.base.exp_centroided.spectra.len());
        progress.start_count(steps, "filtering LC-MS data")?;
        let mut done = 0usize;
        let mut filter_results = Vec::with_capacity(self.base.patterns.len());
        for pattern_idx in 0..self.base.patterns.len() {
            let pattern = self.base.patterns[pattern_idx].clone();
            let mut result = MultiplexFilteredMSExperiment::new();
            self.base.update_white_ms_experiment();
            for idx_rt in 0..self.base.exp_centroided_white.spectra.len() {
                let white = &self.base.exp_centroided_white.spectra;
                if white[idx_rt].peaks.is_empty() {
                    continue;
                }
                done += 1;
                progress.set_count(done)?;
                let rt = white[idx_rt].rt;
                let band = rt_begin(white, rt - self.base.rt_band / 2.0)
                    ..rt_end(white, rt + self.base.rt_band / 2.0);
                for s in 0..white[idx_rt].peaks.len() {
                    let mz = self.base.exp_centroided_white.spectra[idx_rt].peaks[s].mz;
                    let mut peak = MultiplexFilteredPeak::new(
                        mz,
                        rt as f32,
                        self.base.exp_centroided_mapping[idx_rt][s],
                        idx_rt,
                    );
                    if !self
                        .base
                        .filter_peak_positions(mz, band.clone(), &pattern, &mut peak)?
                    {
                        continue;
                    }
                    if !self.base.filter_averagine_model(&pattern, &peak)? {
                        continue;
                    }
                    if !self.base.filter_peptide_correlation(&pattern, &peak)? {
                        continue;
                    }
                    self.base.blacklist_peak(&peak, pattern_idx)?;
                    result.add_peak(peak);
                }
            }
            filter_results.push(result);
        }
        progress.end()?;
        Ok(filter_results)
    }
}

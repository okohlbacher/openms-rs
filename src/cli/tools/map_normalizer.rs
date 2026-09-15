// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Normalizes peak intensities to a percentage of the run maximum.
//!
//! Ports `OpenMS4-topp/src/MapNormalizer.cpp`. The output carries the source's
//! `intensity normalization` processing record.
//!
//! The run maximum is the *combined* one: `main_` calls `exp.updateRanges()`
//! and then `exp.getMaxIntensity()`, and `MSExperiment::updateRanges()` folds
//! the spectrum ranges over every MS level **and** the chromatogram ranges into
//! `combined_ranges_` (`MSExperiment.cpp:665-719`), which is what
//! `getMaxIntensity()` reads (`MSExperiment.h:1059`). A TIC chromatogram sums a
//! whole spectrum, so it routinely outranks the most intense single peak and
//! then sets the scale on its own. See [`MapNormalizer::run_maximum_intensity`].

use crate::cli::{ExitCode, Tool, ToolContext, ToolSpec};
use crate::format::file_handler::FileHandler;
use crate::format::file_types::FileType;
use crate::kernel::{MSExperiment, SummaryLimits};
use crate::metadata::ProcessingAction;
use crate::{Error, Result};

/// The `MapNormalizer` TOPP tool.
pub struct MapNormalizer;

impl Tool for MapNormalizer {
    const NAME: &'static str = "MapNormalizer";
    const DESCRIPTION: &'static str = "Normalizes peak intensities in an MS run.";

    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_input_file("in", "<file>", "", "input file ", true, false, &[])?;
        spec.set_valid_formats("in", &["mzML"])?;
        spec.register_output_file("out", "<file>", "", "output file ", true, false)?;
        spec.set_valid_formats("out", &["mzML"])?;
        Ok(())
    }

    /// Source `main_`, run on the worker pool that `-threads` sizes, as
    /// `TOPPBase::main` applies the setting before `main_`
    /// (`TOPPBase.cpp:408-415`). See [`ToolContext::in_thread_pool`].
    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        ctx.in_thread_pool(|| Self::run_in_pool(ctx))?
    }
}

impl MapNormalizer {
    /// The work [`MSExperiment::combined_ranges_with_limits`] charges for this
    /// map, as its own ceiling for the call.
    ///
    /// That query charges one unit per spectrum, per peak, per chromatogram and
    /// per chromatogram point against [`SummaryLimits::max_work`], whose default
    /// of 50,000,000 no real LC-MS run fits: the benchmark's 1.2 GB Velos run
    /// charges 88,521,983 and would be refused. The map is fully materialized by
    /// the time the tool asks for its ranges — what bounds the input is the
    /// reader's own size-derived ceilings in `src/format/mzml_scaling.rs` — so
    /// the map's own size is the honest ceiling here: still bounded work, and it
    /// can never refuse a map the reader already admitted. The other two limits
    /// are left at their defaults; this query allocates nothing and produces no
    /// output points, so neither is charged.
    fn range_limits(experiment: &MSExperiment) -> SummaryLimits {
        let mut work = experiment
            .spectra
            .len()
            .saturating_add(experiment.chromatograms.len());
        for spectrum in &experiment.spectra {
            work = work.saturating_add(spectrum.peaks.len());
        }
        for chromatogram in &experiment.chromatograms {
            work = work.saturating_add(chromatogram.peaks.len());
        }
        SummaryLimits {
            max_work: work,
            ..SummaryLimits::default()
        }
    }

    /// The run maximum the source normalizes against: `exp.getMaxIntensity()`
    /// after `exp.updateRanges()`.
    ///
    /// `MSExperiment::updateRanges()` extends `combined_ranges_` first with the
    /// spectrum ranges over every MS level and then with the chromatogram
    /// ranges (`MSExperiment.cpp:665-719`), and `getMaxIntensity()` returns
    /// that combined maximum (`MSExperiment.h:1059`).
    /// [`MSExperiment::combined_ranges_with_limits`] is that fold, in that
    /// order. Taking [`MSExperiment::ranges`] instead, which covers spectra
    /// alone, understates the maximum whenever a chromatogram carries the most
    /// intense value and scales every MS1 peak by the wrong factor.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the combined intensity range is
    /// empty — no spectrum peak and no chromatogram point anywhere — because
    /// the source refuses such a run rather than writing it through. `main_`
    /// asks for `getMaxIntensity()` unconditionally, one line before its peak
    /// loop, and `RangeBase::getMax()` (`RangeManager.h:139-146`) throws
    /// `Exception::InvalidRange` on an empty range with no assertion guard, so
    /// a Release build throws too. Executed at the pinned Release build on both
    /// shapes that reach it — scans that carry a retention time but no point,
    /// and an empty `spectrumList` — the C++ exits 8 with *Empty or
    /// uninitialized range object. Did you forget to call updateRanges()?* and
    /// writes no output file.
    ///
    /// The port refuses the same runs and likewise writes nothing, but exits 6
    /// rather than 8: `Error::InvalidValue` and [`Error::InvalidRange`] map to
    /// `ILLEGAL_PARAMETERS` crate-wide, where the source's exceptions of those
    /// names derive from `BaseException` and reach its `UNKNOWN_ERROR` arm. See
    /// the mapping documented at `run_failure` in `src/cli.rs`, and
    /// `tests/data/topp_map_normalizer_provenance.json` for the executed
    /// commands and digests.
    ///
    /// Also returns [`Error::InvalidValue`] on a nonfinite retention time, m/z
    /// or intensity, as the query does.
    fn run_maximum_intensity(experiment: &MSExperiment) -> Result<f64> {
        experiment
            .combined_ranges_with_limits(Self::range_limits(experiment))?
            .intensity
            .map(|range| range.max)
            .ok_or_else(|| Error::InvalidValue("run has no intensities to normalize".into()))
    }

    /// The tool body, as the source `main_`.
    fn run_in_pool(ctx: &ToolContext) -> Result<ExitCode> {
        let mut experiment = FileHandler::load_experiment(ctx.string("in")?, &[FileType::MzMl])?;

        // Source takes the combined maximum over every MS level and every
        // chromatogram, then divides by a hundredth of it, so the most intense
        // value in the run becomes 100. An empty combined intensity range is
        // refused here, as the source refuses it.
        let maximum = Self::run_maximum_intensity(&experiment)?;

        // The source has no guard here and would emit NaN (an all-zero run) or
        // a sign flip (an all-negative run) instead; refusing is the one
        // deliberate deviation recorded in docs/TOPP_CLI_SUPPORT.md.
        let scale = maximum / 100.0;
        if !(scale.is_finite() && scale > 0.0) {
            return Err(Error::InvalidValue(
                "run maximum intensity must be finite and positive".into(),
            ));
        }

        // Only MS1 is scaled; the source leaves higher levels untouched and
        // its commented-out chromatogram branch is not ported.
        for spectrum in &mut experiment.spectra {
            if spectrum.ms_level < 2 {
                for peak in &mut spectrum.peaks {
                    peak.intensity = (f64::from(peak.intensity) / scale) as f32;
                }
            }
        }

        // Source addDataProcessing_(exp, getProcessingInfo_(NORMALIZATION)).
        let processing = ctx.processing_info(&[ProcessingAction::IntensityNormalization])?;
        ctx.add_data_processing(&mut experiment, &processing);
        FileHandler::store_experiment(ctx.string("out")?, &experiment, Some(FileType::MzMl))?;
        Ok(ExitCode::ExecutionOk)
    }
}

#[cfg(test)]
mod tests {
    use super::{MapNormalizer, SummaryLimits};
    use crate::kernel::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};

    /// Two spectra of two and three peaks, one chromatogram of four points:
    /// 2 spectra + 5 peaks + 1 chromatogram + 4 points = 12 units of work.
    fn experiment() -> MSExperiment {
        let mut experiment = MSExperiment::default();
        experiment.spectra.push(MSSpectrum::from_peaks(vec![
            Peak1D::new(100.0, 1.0),
            Peak1D::new(200.0, 2.0),
        ]));
        experiment.spectra.push(MSSpectrum::from_peaks(vec![
            Peak1D::new(300.0, 3.0),
            Peak1D::new(400.0, 4.0),
            Peak1D::new(500.0, 5.0),
        ]));
        experiment
            .chromatograms
            .push(MSChromatogram::from_peaks(vec![
                ChromatogramPeak::new(1.0, 10.0),
                ChromatogramPeak::new(2.0, 20.0),
                ChromatogramPeak::new(3.0, 30.0),
                ChromatogramPeak::new(4.0, 40.0),
            ]));
        experiment
    }

    /// The derived ceiling is exactly what the query charges: one unit less is
    /// refused, so it is not slack, and the exact figure passes, so it can never
    /// refuse a map that is already in memory.
    #[test]
    fn the_derived_ceiling_is_exactly_what_the_query_charges() {
        let experiment = experiment();
        let limits = MapNormalizer::range_limits(&experiment);
        assert_eq!(limits.max_work, 12);
        assert!(experiment.combined_ranges_with_limits(limits).is_ok());
        assert!(
            experiment
                .combined_ranges_with_limits(SummaryLimits {
                    max_work: limits.max_work - 1,
                    ..limits
                })
                .is_err(),
            "the ceiling must be tight, or it is not measuring the real charge"
        );
    }

    /// The maximum is the combined one, and a chromatogram can set it alone.
    #[test]
    fn a_chromatogram_can_carry_the_run_maximum() {
        let experiment = experiment();
        assert_eq!(
            MapNormalizer::run_maximum_intensity(&experiment).unwrap(),
            40.0,
            "the chromatogram's 40 outranks the most intense peak, 5"
        );
        let spectra_only = MSExperiment {
            chromatograms: Vec::new(),
            ..experiment
        };
        assert_eq!(
            MapNormalizer::run_maximum_intensity(&spectra_only).unwrap(),
            5.0
        );
    }

    /// An empty combined intensity range is refused, as the source refuses it:
    /// `getMaxIntensity()` is asked for before the peak loop and throws on an
    /// empty range. The pinned C++ Release build exits 8 and writes nothing on
    /// both shapes; this port exits 6 and writes nothing. See
    /// `an_empty_combined_intensity_range_is_refused` in
    /// `tests/topp_map_normalizer.rs` for the end-to-end case.
    #[test]
    fn an_empty_run_is_refused() {
        let empty = MSExperiment::default();
        assert!(MapNormalizer::run_maximum_intensity(&empty).is_err());
        let mut without_peaks = MSExperiment::default();
        without_peaks.spectra.push(MSSpectrum::new());
        assert!(
            MapNormalizer::run_maximum_intensity(&without_peaks).is_err(),
            "a spectrum with no peaks contributes an RT but no intensity, and \
             the source throws on that intensity range just the same"
        );
    }
}

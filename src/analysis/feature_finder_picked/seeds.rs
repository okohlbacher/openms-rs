// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Isotope-pattern precalculation and seed selection of the picked feature
//! finder (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`).
//!
//! [`SeedStage`](crate::analysis::feature_finder_picked::seeds::SeedStage) runs
//! the source `run_` up to and including seed selection:
//!
//! - step 0: the parameters, the abundance overrides and the score arrays;
//! - step 1: the intensity quantiles and intensity scores
//!   ([`IntensityThresholds`](crate::analysis::feature_finder_picked::scoring::IntensityThresholds));
//! - step 2: the trace scores and local-maximum flags;
//! - step 2.5: one averagine isotope pattern per mass window
//!   ([`IsotopeWindows`](crate::analysis::feature_finder_picked::seeds::IsotopeWindows));
//! - step 3.1: the pattern score of every peak, per charge;
//! - step 3.2: the overall score and the seeds of each charge
//!   ([`ChargeSeeds`](crate::analysis::feature_finder_picked::seeds::ChargeSeeds)).
//!
//! The source runs steps 3.1 to 3.3 charge by charge. Step 3.3, the extension of
//! the seeds, only reads the arrays of its own charge, so computing every
//! charge's seeds first
//! ([`SeedStage::compute`](crate::analysis::feature_finder_picked::seeds::SeedStage::compute))
//! gives the same arrays and seeds. The extension is
//! [`extension`](crate::analysis::feature_finder_picked::extension), and
//! [`feature_stage`](crate::analysis::feature_finder_picked::algorithm::feature_stage)
//! puts the log lines back into the source's order, each charge's seed count
//! followed by its candidate count. The algorithm instance
//! ([`FeatureFinderAlgorithmPicked`](crate::analysis::feature_finder_picked::instance::FeatureFinderAlgorithmPicked))
//! keeps the source's order instead, selecting one charge's seeds right before
//! extending them, so its log lines, progress calls and debug files come in
//! the source's order.
//!
//! Isotope patterns are computed in the source's binary32 arithmetic
//! ([`ProbabilityPrecision::SourceSingle`](crate::chemistry::isotopes::ProbabilityPrecision::SourceSingle))
//! with the source `trimLeft`. See
//! `docs/FEATURE_FINDER_PICKED_SUPPORT.md` and
//! `docs/ISOTOPE_SOURCE_PRECISION_SUPPORT.md`.

use crate::analysis::feature_finder_picked::algorithm::{
    AbundanceOverride, Limits, Options, Settings, validate_input,
};
use crate::analysis::feature_finder_picked::debug::{LOG_PRECALCULATING, LogSink, NoLog};
use crate::analysis::feature_finder_picked::helper_structs::{
    IsotopePattern, PatternPeak, Seed, TheoreticalIsotopePattern,
};
use crate::analysis::feature_finder_picked::instance::Progress;
use crate::analysis::feature_finder_picked::scoring::{
    IntensityThresholds, ScoreArrays, Work, fill_intensity_scores, fill_trace_scores,
    find_isotope_logged, isotope_score_logged, libstdcxx, ms1_ranges, nearest, reset_pattern,
    x86_64,
};
use crate::analysis::feature_finder_picked::source_sort::source_sort_reversed_by;
use crate::chemistry::isotopes::{
    CoarseIsotopePatternGenerator, CoarseMassMode, IsotopeDistribution, IsotopePeak,
    ProbabilityPrecision,
};
use crate::kernel::{FeatureMap, MSExperiment};
use crate::param::Param;
use crate::{Error, Result};

/// The precalculated theoretical isotope patterns, one per mass window (source
/// member `isotope_distributions_`, filled in step 2.5 of `run_`).
#[derive(Clone, Debug, PartialEq)]
pub struct IsotopeWindows {
    mass_window_width: f64,
    patterns: Vec<TheoreticalIsotopePattern>,
}

impl IsotopeWindows {
    /// Precalculate the patterns for masses up to `max_mz * charge_high`: step
    /// 2.5 of source `run_`.
    ///
    /// There are `ceil(max_mz * charge_high / width) + 1` windows, converted to
    /// `Size` as the Linux x86_64 Release build converts it (crate-private
    /// `x86_64::truncate_to_u64`, the `cvttsd2si`/`btc` sequence at
    /// `libOpenMS.so` `0x18e46f4`): an infinite maximum m/z, or one whose count
    /// reaches `2^64`, gives **no** window, and the first pattern lookup of
    /// step 3.1 then fails ([`Self::get`]); the executed build throws the same
    /// `Exception::InvalidValue` there. A count in `[2^63, 2^64)` makes the
    /// source's `resize` throw `std::length_error` and is refused here by
    /// [`Limits::max_isotope_windows`] first. Window `i`
    /// holds the averagine estimate for the peptide mass `0.5 * width + i *
    /// width` from a coarse generator limited to `settings.max_isotopes()` peaks,
    /// in source binary32 precision and with the abundance overrides of
    /// `options`. The estimate is trimmed on the left with the source `trimLeft`
    /// and then on the right, both at `intensity_percentage_optional`; the number
    /// of isotopes removed on the left is kept as `trimmed_left`. The leading
    /// isotopes below `intensity_percentage` are optional at the beginning, and
    /// those below it after the first required isotope are optional at the end.
    /// The intensities are finally divided by their maximum, which is kept as
    /// `max`.
    ///
    /// A pattern whose every isotope lies below the trimming cutoff is kept whole
    /// by the left trim and then emptied by the right trim, as in the source; its
    /// `max` is zero.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the converted window count exceeds
    /// [`Limits::max_isotope_windows`], when the windows could hold more
    /// than [`Limits::max_pattern_values`] values, or when the isotope generator
    /// fails (for example because every retained bin underflows in binary32,
    /// where the source produces NaN weights). Returns [`Error::Unsupported`] for
    /// a changed abundance under [`AbundanceOverride::Refuse`].
    pub fn precalculate(max_mz: f64, settings: &Settings, options: &Options) -> Result<Self> {
        Self::precalculate_onto(None, max_mz, settings, options)
    }

    /// [`Self::precalculate`] on the windows an earlier run of the same
    /// algorithm object left behind: source step 2.5 on a member that `run_`
    /// never clears.
    ///
    /// The source resizes `isotope_distributions_` to the new window count,
    /// which keeps the earlier windows below that count, and then *appends*
    /// each window's new intensities to the kept ones
    /// (`FeatureFinderAlgorithmPicked.cpp:361`, `:378`). The optional isotopes,
    /// the maximum and the normalisation are then computed over the whole,
    /// longer vector, and `trimmed_left` is overwritten. A reused object
    /// therefore works with longer patterns whose earlier part is already
    /// normalised, and usually a maximum of 1, from its second run on; the
    /// executed Release build finds different features in the second run on
    /// the same input for that reason. `None` gives the first run.
    ///
    /// # Errors
    ///
    /// As [`Self::precalculate`]; the pattern-value ceiling counts the kept
    /// values too.
    pub fn precalculate_onto(
        previous: Option<&IsotopeWindows>,
        max_mz: f64,
        settings: &Settings,
        options: &Options,
    ) -> Result<Self> {
        let width = settings.mass_window_width;
        let max_mass = max_mz * f64::from(settings.charge_high);
        // Source `Size num_isotopes = std::ceil(max_mass / mass_window_width_) + 1`.
        let count = x86_64::truncate_to_u64((max_mass / width).ceil() + 1.0);
        let limits = options.limits;
        let count = usize::try_from(count)
            .ok()
            .filter(|&count| count <= limits.max_isotope_windows)
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "{count} isotope windows for a maximum mass of {max_mass} Da and a window \
                     width of {width} Da exceed the limit of {}",
                    limits.max_isotope_windows
                ))
            })?;
        let max_isotopes = settings.max_isotopes();
        let kept: &[TheoreticalIsotopePattern] = previous.map_or(&[], |windows| {
            &windows.patterns[..count.min(windows.patterns.len())]
        });
        let kept_values = kept
            .iter()
            .fold(0usize, |total, pattern| total.saturating_add(pattern.len()));
        if count
            .saturating_mul(max_isotopes)
            .saturating_add(kept_values)
            > limits.max_pattern_values
        {
            return Err(Error::InvalidValue(format!(
                "{count} isotope windows of up to {max_isotopes} isotopes, and {kept_values} \
                 values an earlier run left, exceed the limit of {} pattern values",
                limits.max_pattern_values
            )));
        }
        let generator = pattern_generator(settings, options.abundance_override)?;
        let mut patterns = Vec::new();
        patterns
            .try_reserve_exact(count)
            .map_err(|_| Error::InvalidValue("cannot allocate the isotope windows".into()))?;
        patterns.extend_from_slice(kept);
        patterns.resize_with(count, TheoreticalIsotopePattern::default);
        for (index, pattern) in patterns.iter_mut().enumerate() {
            let mut distribution =
                generator.estimate_from_peptide_weight(0.5 * width + index as f64 * width)?;
            theoretical_pattern_onto(pattern, &mut distribution, settings)?;
        }
        Ok(Self {
            mass_window_width: width,
            patterns,
        })
    }

    /// The pattern for `mass`: source `getIsotopeDistribution_`, window
    /// `(Size) floor(mass / width)`.
    ///
    /// The conversion to `Size` is undefined in C++ for a NaN, infinite,
    /// negative or too large quotient; this follows the instructions the
    /// Linux x86_64 Release build emits for it (crate-private
    /// `x86_64::truncate_to_u64`): a NaN mass gives index `2^63`, an infinite
    /// one index 0, a negative one a wrapped index. The index is then compared
    /// with the window count as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the window was not precalculated,
    /// with the `what()` text of the source's `Exception::InvalidValue`:
    /// `the value '<index>' was used but is not valid; IsotopeDistribution not
    /// precalculated. Maximum allowed index is <count>`. The executed Release
    /// build gives exactly this text for a NaN m/z (index
    /// `9223372036854775808`), an infinite m/z (no window) and an m/z of
    /// `1e300` (no window).
    pub fn get(&self, mass: f64) -> Result<&TheoreticalIsotopePattern> {
        let index = x86_64::truncate_to_u64((mass / self.mass_window_width).floor());
        usize::try_from(index)
            .ok()
            .and_then(|index| self.patterns.get(index))
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "the value '{index}' was used but is not valid; IsotopeDistribution not \
                     precalculated. Maximum allowed index is {}",
                    self.patterns.len()
                ))
            })
    }

    /// The window width in Da.
    pub fn mass_window_width(&self) -> f64 {
        self.mass_window_width
    }

    /// The patterns, by window index.
    pub fn patterns(&self) -> &[TheoreticalIsotopePattern] {
        &self.patterns
    }
}

/// The coarse generator of step 2.5, with the overrides of step 0.
fn pattern_generator(
    settings: &Settings,
    policy: AbundanceOverride,
) -> Result<CoarseIsotopePatternGenerator> {
    let mut generator = CoarseIsotopePatternGenerator::new(
        Some(settings.max_isotopes()),
        CoarseMassMode::Approximate,
    )?
    .with_precision(ProbabilityPrecision::SourceSingle);
    let overrides = [
        (
            settings.abundance_12c_changed,
            "isotopic_pattern:abundance_12C",
            "C",
            12.0,
            settings.abundance_12c,
        ),
        (
            settings.abundance_14n_changed,
            "isotopic_pattern:abundance_14N",
            "N",
            14.0,
            settings.abundance_14n,
        ),
    ];
    for (changed, key, symbol, light, abundance) in overrides {
        if !changed {
            continue;
        }
        if policy == AbundanceOverride::Refuse {
            return Err(Error::Unsupported(format!(
                "{key} = {abundance}: the source builds this isotope override with a stray \
                 (0, 1) peak and the port does not reproduce it; select \
                 AbundanceOverride::Intended for the intended two-isotope override"
            )));
        }
        let distribution = IsotopeDistribution::from_peaks(vec![
            IsotopePeak {
                mass: light,
                probability: abundance / 100.0,
            },
            IsotopePeak {
                mass: light + 1.0,
                probability: 1.0 - (abundance / 100.0),
            },
        ])?;
        generator.set_isotope_override(symbol, distribution)?;
    }
    Ok(generator)
}

/// Trim, classify and normalise one estimate into `target`: the body of the
/// step 2.5 loop, appending to what `target` already holds as the source
/// appends to its member.
fn theoretical_pattern_onto(
    target: &mut TheoreticalIsotopePattern,
    distribution: &mut IsotopeDistribution,
    settings: &Settings,
) -> Result<()> {
    let size_before = distribution.len();
    distribution.trim_left_source(settings.intensity_percentage_optional)?;
    let trimmed_left = size_before - distribution.len();
    distribution.trim_right(settings.intensity_percentage_optional)?;
    let mut intensity = std::mem::take(&mut target.intensity);
    intensity
        .try_reserve(distribution.len())
        .map_err(|_| Error::InvalidValue("cannot allocate an isotope window".into()))?;
    intensity.extend(distribution.peaks().iter().map(|peak| peak.probability));
    let mut begin = 0;
    let mut end = 0;
    let mut is_begin = true;
    let mut is_end = false;
    for &value in &intensity {
        if value < settings.intensity_percentage {
            if !is_end && !is_begin {
                is_end = true;
            }
            if is_begin {
                begin += 1;
            } else if is_end {
                end += 1;
            }
        } else if is_begin {
            is_begin = false;
        }
    }
    let mut max = 0.0;
    for &value in &intensity {
        if value > max {
            max = value;
        }
    }
    for value in &mut intensity {
        *value /= max;
    }
    *target = TheoreticalIsotopePattern {
        intensity,
        optional_begin: begin,
        optional_end: end,
        max,
        trimmed_left,
    };
    Ok(())
}

/// A user-specified seed position, the part of a seed feature the source reads.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UserSeed {
    /// m/z of the seed feature.
    pub mz: f64,
    /// Retention time of the seed feature in seconds.
    pub rt: f64,
}

/// The seeds of one charge, most intense first.
#[derive(Clone, Debug, PartialEq)]
pub struct ChargeSeeds {
    /// The charge.
    pub charge: i32,
    /// The seeds, in the order step 3.3 extends them.
    pub seeds: Vec<Seed>,
}

/// The state of a run after seed selection: steps 0 to 3.2 of source `run_`.
///
/// It owns the (possibly sorted) experiment, the typed settings, the user seeds
/// sorted by m/z, the intensity quantiles, the isotope windows, the per-peak
/// score arrays, the seeds of every charge and the log so far. Seed extension
/// builds on it.
#[derive(Clone, Debug, PartialEq)]
pub struct SeedStage {
    experiment: MSExperiment,
    settings: Settings,
    user_seeds: Vec<UserSeed>,
    thresholds: IntensityThresholds,
    windows: IsotopeWindows,
    scores: ScoreArrays,
    charges: Vec<ChargeSeeds>,
    log: Vec<String>,
    /// Number of charges the stage selects seeds for.
    charge_count: usize,
    /// The scoring work budget left for the charges not yet selected.
    work: Work,
}

impl SeedStage {
    /// Validate the input, apply the parameters and select seeds, with default
    /// [`Options`].
    ///
    /// Returns `Ok(None)` for an experiment without spectra, where source `run`
    /// returns before touching the parameters.
    ///
    /// # Errors
    ///
    /// As [`Self::run_with_options`].
    pub fn run(
        experiment: MSExperiment,
        user_seeds: &FeatureMap,
        parameters: &Param,
    ) -> Result<Option<Self>> {
        Self::run_with_options(experiment, user_seeds, parameters, &Options::default())
    }

    /// Validate the input ([`validate_input`]), apply the parameters
    /// ([`Settings::from_parameters`]) and select seeds ([`Self::compute`]), in
    /// the source's order.
    ///
    /// The log starts with the unsorted-input warning, if any, followed by the
    /// unknown-parameter warnings.
    ///
    /// # Errors
    ///
    /// Every error of the three steps.
    pub fn run_with_options(
        mut experiment: MSExperiment,
        user_seeds: &FeatureMap,
        parameters: &Param,
        options: &Options,
    ) -> Result<Option<Self>> {
        let mut log = Vec::new();
        if !validate_input(&mut experiment, &mut log)? {
            return Ok(None);
        }
        let (settings, warnings) = Settings::from_parameters(parameters)?;
        log.extend(warnings);
        Self::compute(experiment, user_seeds, settings, options, log).map(Some)
    }

    /// Steps 0 to 3.2 of source `run_` on validated input.
    ///
    /// `experiment` must hold MS1 spectra sorted as [`validate_input`] leaves
    /// them; `log` is extended with one line
    /// `Found <n> seeds for charge <c>.` per charge, which the source prints to
    /// `std::cout`.
    ///
    /// `mass_trace:min_spectra = 1` makes the number of scans inspected per side
    /// zero, which the source does not check: it inspects no neighbouring scan,
    /// divides every trace score by zero, scores every peak NaN and therefore
    /// finds no seed, returning an empty map. That is defined behaviour and the
    /// port reproduces it (`CPP-271`); nothing later in the algorithm is
    /// reached, so the source's `size_t(-1)` delta buffer in
    /// `extendMassTrace_` stays unreachable.
    ///
    /// Seeds are the local trace maxima whose overall score, the `f32` cube root
    /// `powf(trace * intensity * pattern, 1/3)`, reaches `seed:min_score`. With
    /// user seeds, the threshold is `user-seed:min_score` instead, and a peak
    /// must also lie strictly within `user-seed:mz_tolerance` and
    /// `user-seed:rt_tolerance` of some user seed. The cube root is
    /// [`overall_score`]. Seeds are sorted by descending `f32` intensity,
    /// the order [`Seed::is_less_intense_than`] defines. The source sorts with
    /// `std::sort(seeds.rbegin(), seeds.rend())`, which leaves seeds of equal
    /// intensity in an order the standard does not specify; the port puts them
    /// where the Linux x86_64 Release build's libstdc++ introsort does
    /// ([`crate::analysis::feature_finder_picked::source_sort`]).
    ///
    /// `write_debug` does not change what this stage computes. Its log lines
    /// and seed maps are produced by the algorithm instance, which runs the same
    /// steps with the debug output attached.
    ///
    /// # Errors
    ///
    /// - [`Error::Unsupported`] for a changed isotope abundance under
    ///   [`AbundanceOverride::Refuse`], which is not the default; see that type.
    /// - [`Error::InvalidValue`] when `charge_low` exceeds `charge_high` by more
    ///   than one ([`Settings::charge_count`]), for a zero or infinite
    ///   intensity bin step read by the seed loop under
    ///   [`DegenerateBinStep::Refuse`](crate::analysis::feature_finder_picked::algorithm::DegenerateBinStep::Refuse),
    ///   which is not the default, for a NaN user-seed m/z among different
    ///   seed m/z values (the source's `std::sort` is then undefined), for the
    ///   undefined cases of [`IntensityThresholds::compute`], and when a [`Limits`]
    ///   ceiling is exceeded, checked before the allocation or computation it
    ///   bounds. A zero or infinite step under the default
    ///   [`DegenerateBinStep::Source`](crate::analysis::feature_finder_picked::algorithm::DegenerateBinStep::Source)
    ///   is computed as the Linux x86_64 Release build computes it.
    /// - [`Error::InvalidRange`] when every retention time or every m/z is
    ///   NaN, and [`Error::InvalidValue`] with the source's text when a pattern
    ///   lookup of step 3.1 finds no precalculated window
    ///   ([`IsotopeWindows::get`]): for a NaN m/z, and for every input when an
    ///   infinite or huge maximum m/z left no window.
    pub fn compute(
        experiment: MSExperiment,
        user_seeds: &FeatureMap,
        settings: Settings,
        options: &Options,
        log: Vec<String>,
    ) -> Result<Self> {
        let mut progress = Progress::silent();
        let mut stage = Self::prepare(
            experiment,
            user_seeds,
            settings,
            options,
            log,
            None,
            &mut NoLog,
            &mut progress,
        )?;
        while stage.select_next_charge(&mut NoLog, &mut progress)? {
            progress.end()?;
            stage.log_seed_count();
        }
        Ok(stage)
    }

    /// Steps 0 to 2.5 of source `run_`, leaving the charges for
    /// [`Self::select_next_charge`].
    ///
    /// `log` receives the source's first `log_` line and `progress` the
    /// source's `startProgress`/`setProgress`/`endProgress` calls of steps 1, 2
    /// and 2.5, with their labels and values, in the source's order. The
    /// `setProgress` calls of steps 1 and 2 are made after each step's loop
    /// rather than inside it; the progress logger shows them only when a
    /// wall-clock second has passed, so which of them it prints depends on
    /// timing in the source as here.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn prepare<L: LogSink>(
        experiment: MSExperiment,
        user_seeds: &FeatureMap,
        settings: Settings,
        options: &Options,
        log: Vec<String>,
        previous_windows: Option<&IsotopeWindows>,
        debug_log: &mut L,
        progress: &mut Progress<'_>,
    ) -> Result<Self> {
        let limits = options.limits;
        // The source sorts the user seeds first (`FeatureFinderAlgorithmPicked.cpp:190`).
        let user_seeds = sorted_user_seeds(user_seeds)?;
        preflight(&experiment, &settings, options)?;
        let charge_count = settings.charge_count()?;
        if settings.abundance_12c_changed || settings.abundance_14n_changed {
            // Fails early under the refusing policy, before any work.
            pattern_generator(&settings, options.abundance_override)?;
        }
        let mut scores = ScoreArrays::new(
            &experiment,
            settings.charge_low,
            charge_count,
            limits.max_score_bytes,
        )?;
        let mut work = Work::new(limits.max_work);

        // Step 1: intensity quantiles and scores.
        if debug_log.enabled() {
            debug_log.put(LOG_PRECALCULATING);
        }
        let bins = settings.intensity_bins;
        let cells = i64::try_from(bins.saturating_mul(bins))
            .map_err(|_| Error::InvalidValue("intensity bin count overflow".into()))?;
        progress.start(0, cells, "Precalculating intensity scores")?;
        let thresholds = IntensityThresholds::compute_with_work(
            &experiment,
            settings.intensity_bins,
            &mut work,
        )?;
        progress.set_each(0..cells)?;
        fill_intensity_scores(&experiment, &thresholds, &mut scores, &mut work)?;
        progress.end()?;

        // Step 2: trace scores and local maxima.
        let spectra = experiment.spectra.len();
        let end_iteration = spectra - settings.min_spectra.min(spectra);
        progress.start(
            progress_value(settings.min_spectra)?,
            progress_value(end_iteration)?,
            "Precalculating mass trace scores",
        )?;
        fill_trace_scores(
            &experiment,
            settings.min_spectra,
            settings.trace_tolerance,
            &mut scores,
            &mut work,
        )?;
        progress.set_each(progress_value(settings.min_spectra)?..progress_value(end_iteration)?)?;
        progress.end()?;

        // Step 2.5: isotope patterns per mass window.
        let (_, mz_range) = ms1_ranges(&experiment)?;
        let windows =
            IsotopeWindows::precalculate_onto(previous_windows, mz_range.max, &settings, options)?;
        progress.start(
            0,
            progress_value(windows.patterns().len())?,
            "Precalculating isotope distributions",
        )?;
        progress.end()?;

        Ok(Self {
            experiment,
            settings,
            user_seeds,
            thresholds,
            windows,
            scores,
            charges: Vec::with_capacity(charge_count),
            log,
            charge_count,
            work,
        })
    }

    /// Steps 3.1 and 3.2 of source `run_` for the next charge: its pattern
    /// scores and its seeds, with the `Found <n> seeds for charge <c>.` line.
    ///
    /// Returns `false` when every charge is selected. `log` receives the
    /// source's per-peak `log_` lines of step 3.1, and `progress` the progress
    /// calls of both steps (the `setProgress` calls of step 3.1 after its
    /// loop). Step 3.2's progress is left open: the source stores the debug
    /// seed map before it calls `endProgress`, and prints the seed count after
    /// that ([`Self::log_seed_count`]), so the caller does both.
    pub(crate) fn select_next_charge<L: LogSink>(
        &mut self,
        debug_log: &mut L,
        progress: &mut Progress<'_>,
    ) -> Result<bool> {
        let charge_index = self.charges.len();
        if charge_index >= self.charge_count {
            return Ok(false);
        }
        let charge = self.settings.charge_low + charge_index as i32;
        let spectra = progress_value(self.experiment.spectra.len())?;
        progress.start(
            0,
            spectra,
            &format!("Calculating isotope pattern scores for charge {charge}"),
        )?;
        fill_pattern_scores(
            &self.experiment,
            &self.windows,
            &self.settings,
            charge,
            charge_index,
            &mut self.scores,
            &mut self.work,
            debug_log,
        )?;
        progress.set_each(0..spectra)?;
        progress.end()?;
        let end_iteration = self.experiment.spectra.len()
            - self.settings.min_spectra.min(self.experiment.spectra.len());
        let begin = progress_value(self.settings.min_spectra)?;
        let end = progress_value(end_iteration)?;
        progress.start(begin, end, &format!("Finding seeds for charge {charge}"))?;
        let seeds = select_seeds(
            &self.experiment,
            &self.settings,
            &self.user_seeds,
            charge_index,
            &mut self.scores,
            &mut |s| progress.set(progress_value(s)?),
        )?;
        self.charges.push(ChargeSeeds { charge, seeds });
        Ok(true)
    }

    /// Record the source's `std::cout` line of step 3.2 for the charge selected
    /// last, which the source prints after the debug seed file.
    pub(crate) fn log_seed_count(&mut self) {
        if let Some(last) = self.charges.last() {
            self.log.push(format!(
                "Found {} seeds for charge {}.",
                last.seeds.len(),
                last.charge
            ));
        }
    }

    /// Number of charges the stage selects seeds for.
    pub fn charge_count(&self) -> usize {
        self.charge_count
    }

    /// The experiment, sorted when the input was not.
    pub fn experiment(&self) -> &MSExperiment {
        &self.experiment
    }

    /// The typed settings.
    pub fn settings(&self) -> &Settings {
        &self.settings
    }

    /// The user seeds, sorted by m/z (source `seeds_.sortByMZ()`).
    pub fn user_seeds(&self) -> &[UserSeed] {
        &self.user_seeds
    }

    /// The intensity quantiles of step 1.
    pub fn thresholds(&self) -> &IntensityThresholds {
        &self.thresholds
    }

    /// The isotope windows of step 2.5.
    pub fn windows(&self) -> &IsotopeWindows {
        &self.windows
    }

    /// The per-peak scores.
    pub fn scores(&self) -> &ScoreArrays {
        &self.scores
    }

    /// The seeds of every charge, from `charge_low` upwards.
    pub fn charges(&self) -> &[ChargeSeeds] {
        &self.charges
    }

    /// The log lines so far.
    pub fn log(&self) -> &[String] {
        &self.log
    }

    /// The experiment, consuming the stage.
    pub fn into_experiment(self) -> MSExperiment {
        self.experiment
    }
}

/// Size ceilings that bound every allocation of the stage, and the
/// [`DegenerateBinStep::Refuse`](crate::analysis::feature_finder_picked::algorithm::DegenerateBinStep::Refuse)
/// check.
fn preflight(experiment: &MSExperiment, settings: &Settings, options: &Options) -> Result<()> {
    use crate::analysis::feature_finder_picked::algorithm::DegenerateBinStep;
    let limits: &Limits = &options.limits;
    if experiment.spectra.len() > limits.max_spectra {
        return Err(Error::InvalidValue(format!(
            "{} spectra exceed the limit of {}",
            experiment.spectra.len(),
            limits.max_spectra
        )));
    }
    let peaks = experiment
        .spectra
        .iter()
        .try_fold(0usize, |total, spectrum| {
            total.checked_add(spectrum.peaks.len())
        })
        .ok_or_else(|| Error::InvalidValue("peak count overflow".into()))?;
    if peaks > limits.max_peaks {
        return Err(Error::InvalidValue(format!(
            "{peaks} peaks exceed the limit of {}",
            limits.max_peaks
        )));
    }
    let charges = settings.charge_count()?;
    if charges > limits.max_charges {
        return Err(Error::InvalidValue(format!(
            "{charges} charges exceed the limit of {}",
            limits.max_charges
        )));
    }
    if settings.intensity_bins > limits.max_intensity_bins {
        return Err(Error::InvalidValue(format!(
            "intensity:bins {} exceeds the limit of {}",
            settings.intensity_bins, limits.max_intensity_bins
        )));
    }
    if settings.max_isotopes() > IsotopePattern::MAX_SIZE {
        return Err(Error::InvalidValue(
            "isotope pattern size exceeds IsotopePattern::MAX_SIZE".into(),
        ));
    }
    if options.degenerate_bin_step == DegenerateBinStep::Refuse {
        refuse_degenerate_bin_step(experiment, settings)?;
    }
    Ok(())
}

/// [`DegenerateBinStep::Refuse`](crate::analysis::feature_finder_picked::algorithm::DegenerateBinStep::Refuse):
/// an [`Error::InvalidValue`] when a bin step is zero or infinite and the seed
/// loop reads the resulting intensity scores.
///
/// The seed loop (`FeatureFinderAlgorithmPicked.cpp:493-498`) visits the scans
/// `min_spectra .. n - min(min_spectra, n)`, which is empty for `n <= 2 *
/// min_spectra`. There the source computes every intensity score but never
/// reads one, and its result, an empty feature map, does not depend on them;
/// such an input is not refused.
fn refuse_degenerate_bin_step(experiment: &MSExperiment, settings: &Settings) -> Result<()> {
    use crate::analysis::feature_finder_picked::scoring::{bin_steps, degenerate_step};
    let spectra = experiment.spectra.len();
    let loop_end = spectra - settings.min_spectra.min(spectra);
    if settings.min_spectra >= loop_end {
        return Ok(());
    }
    let (rt, mz) = ms1_ranges(experiment)?;
    let (rt_step, mz_step) = bin_steps(&rt, &mz, settings.intensity_bins);
    if degenerate_step(rt_step) || degenerate_step(mz_step) {
        return Err(Error::InvalidValue(format!(
            "FeatureFinderAlgorithmPicked: the intensity bin steps are {rt_step} (RT {} to {}) \
             and {mz_step} (m/z {} to {}) for {} bins; the source converts floor(NaN) or \
             floor(inf) to UInt for every peak here, which is undefined behaviour, and \
             DegenerateBinStep::Refuse is selected",
            rt.min, rt.max, mz.min, mz.max, settings.intensity_bins
        )));
    }
    Ok(())
}

/// The user seed positions sorted by m/z: source `seeds_.sortByMZ()`
/// (`FeatureFinderAlgorithmPicked.cpp:190`, `std::sort` with `Feature::MZLess`).
///
/// Infinite positions and a NaN retention time are ordinary values here: the
/// retention time is no sort key, and `<` orders infinities. A NaN m/z is
/// equivalent to every other m/z under `<`, so it is a strict weak ordering
/// only while every other m/z is equivalent to every other; with two different
/// non-NaN m/z values the source's `std::sort` is undefined, and that is
/// [`Error::InvalidValue`] here. Otherwise the order of the equivalent seeds is
/// unspecified in the source; this sorts stably by `f64::total_cmp`, and
/// [`near_user_seed`] gives the same answer for every order of such seeds.
fn sorted_user_seeds(seeds: &FeatureMap) -> Result<Vec<UserSeed>> {
    let mut sorted = Vec::with_capacity(seeds.features.len());
    let mut has_nan = false;
    let mut low = f64::INFINITY;
    let mut high = f64::NEG_INFINITY;
    for feature in &seeds.features {
        if feature.mz.is_nan() {
            has_nan = true;
        } else {
            low = low.min(feature.mz);
            high = high.max(feature.mz);
        }
        sorted.push(UserSeed {
            mz: feature.mz,
            rt: feature.rt,
        });
    }
    if has_nan && low < high {
        return Err(Error::InvalidValue(format!(
            "FeatureFinderAlgorithmPicked: a user seed has a NaN m/z among seeds with the \
             different m/z values {low} and {high}; the source sorts them with std::sort, whose \
             comparator is then not a strict weak ordering (undefined behaviour)"
        )));
    }
    sorted.sort_by(|a, b| a.mz.total_cmp(&b.mz));
    Ok(sorted)
}

/// Pattern scores of one charge: step 3.1 of source `run_`.
///
/// For every peak at m/z `x`, the pattern of window `x * charge` is placed so
/// that its most intense isotope (the first of equal maxima) lies on the peak:
/// isotope `i` is searched at `x + (i - max_isotope) / charge` with
/// [`find_isotope`], starting from the peak nearest to `x - (len + 1) /
/// charge`. The [`isotope_score`](crate::analysis::feature_finder_picked::scoring::isotope_score)
/// with m/z distances raises the stored pattern score of every matched,
/// non-removed peak to that score when it is higher, compared in `f64` and
/// stored as `f32`.
#[allow(clippy::too_many_arguments)]
fn fill_pattern_scores<L: LogSink>(
    experiment: &MSExperiment,
    windows: &IsotopeWindows,
    settings: &Settings,
    charge: i32,
    charge_index: usize,
    scores: &mut ScoreArrays,
    work: &mut Work,
    log: &mut L,
) -> Result<()> {
    let spectra = &experiment.spectra;
    let offsets: Vec<usize> = (0..spectra.len()).map(|s| scores.offset(s)).collect();
    let target = scores.pattern_mut(charge_index);
    let c = f64::from(charge);
    let mut pattern = IsotopePattern::default();
    for (s, spectrum) in spectra.iter().enumerate() {
        for peak in &spectrum.peaks {
            let mz = peak.mz;
            let isotopes = windows.get(mz * c)?;
            let size = isotopes.len();
            let mut max_isotope = 0;
            for (i, &value) in isotopes.intensity.iter().enumerate() {
                if isotopes.intensity[max_isotope] < value {
                    max_isotope = i;
                }
            }
            let mut peak_index = nearest(&spectrum.peaks, mz - (size + 1) as f64 / c, |p| p.mz)
                .ok_or_else(|| Error::InvalidValue("empty spectrum in step 3.1".into()))?;
            reset_pattern(&mut pattern, size);
            let mut units = 1u64;
            for i in 0..size {
                let isotope_pos = mz + (i as f64 - max_isotope as f64) / c;
                units += find_isotope_logged(
                    spectra,
                    isotope_pos,
                    s,
                    &mut pattern,
                    i,
                    &mut peak_index,
                    settings.pattern_tolerance,
                    log,
                )?;
            }
            let (pattern_score, score_units) = isotope_score_logged(
                isotopes,
                &mut pattern,
                true,
                settings.min_isotope_fit,
                settings.optional_fit_improvement,
                log,
            )?;
            work.consume(units + score_units)?;
            if pattern_score > 0.0 {
                for (found, &found_spectrum) in pattern.peak.iter().zip(&pattern.spectrum) {
                    let PatternPeak::Found(found_peak) = *found else {
                        continue;
                    };
                    let slot = &mut target[offsets[found_spectrum] + found_peak];
                    if pattern_score > f64::from(*slot) {
                        *slot = pattern_score as f32;
                    }
                }
            }
        }
    }
    Ok(())
}

/// The overall score of step 3.2: `pow(trace * intensity * pattern, 1.0f / 3.0f)`,
/// the `f32` product formed left to right.
///
/// The source calls `std::pow(float, float)`, the `float` overload, which is
/// the platform's `powf` and is not correctly rounded everywhere: the executed
/// Apple `powf` of the C++ oracle misrounds 99 of the 30,840 overall scores of
/// the retained executions, by one binary32 step each (12 of 3,084 on
/// FeatureFinderCentroided_1). This evaluates the power in `f64` with
/// `libm::pow`, a pure-Rust port of musl's, and rounds once to `f32`. That is
/// the correctly rounded binary32 power for every one of those scores (checked
/// against 60-digit decimal arithmetic) and the same on every machine, so it
/// agrees with the oracle everywhere except on the scores the oracle misrounds.
/// `libm::powf`, the direct binary32 port, disagreed with the oracle on 226 of
/// the 3,084 FeatureFinderCentroided_1 scores. Against the Linux x86_64 Release
/// build, whose glibc 2.39 `powf` is not correctly rounded either, 8 of the
/// 30,840 retained scores are one binary32 step apart (`CPP-272`).
///
/// A NaN product, which a zero or infinite intensity bin step produces
/// ([`DegenerateBinStep`](crate::analysis::feature_finder_picked::algorithm::DegenerateBinStep)),
/// gives a NaN score whose bits are the product's, as
/// the executed `powf` returns them.
pub fn overall_score(trace: f32, intensity: f32, pattern: f32) -> f32 {
    let product = trace * intensity * pattern;
    crate::analysis::feature_finder_picked::scoring::x86_64::narrow(libm::pow(
        f64::from(product),
        f64::from(1.0f32 / 3.0f32),
    ))
}

/// Overall scores and seeds of one charge: step 3.2 of source `run_`.
///
/// `progress` is called with every scan index the loop visits, as the source
/// calls `setProgress`.
fn select_seeds(
    experiment: &MSExperiment,
    settings: &Settings,
    user_seeds: &[UserSeed],
    charge_index: usize,
    scores: &mut ScoreArrays,
    progress: &mut dyn FnMut(usize) -> Result<()>,
) -> Result<Vec<Seed>> {
    let spectra = &experiment.spectra;
    let end = spectra.len() - settings.min_spectra.min(spectra.len());
    let starts: Vec<usize> = (0..spectra.len()).map(|s| scores.offset(s)).collect();
    let (trace, intensity, local_max, pattern, overall) = scores.seed_inputs(charge_index);
    let use_user_seeds = !user_seeds.is_empty();
    let mut seeds = Vec::new();
    for s in settings.min_spectra..end {
        progress(s)?;
        let spectrum = &spectra[s];
        for (p, peak) in spectrum.peaks.iter().enumerate() {
            let slot = starts[s] + p;
            let score = overall_score(trace[slot], intensity[slot], pattern[slot]);
            overall[slot] = score;
            if local_max[slot] == 0.0 {
                continue;
            }
            let score = f64::from(score);
            let seed = Seed::new(s, p, peak.intensity);
            if !use_user_seeds {
                if score >= settings.seed_min_score {
                    seeds.push(seed);
                }
            } else if score >= settings.user_seed_min_score
                && near_user_seed(user_seeds, peak.mz, spectrum.rt, settings)
            {
                seeds.push(seed);
            }
        }
    }
    // Source: `std::sort(seeds.rbegin(), seeds.rend())` with `Seed::operator<`.
    // The intensities are finite (validated input), so the order is strict weak
    // and the sort cannot fail.
    source_sort_reversed_by(&mut seeds, Seed::is_less_intense_than)?;
    Ok(seeds)
}

/// Whether a user seed lies strictly within the user-seed tolerances of a peak.
///
/// Source: `std::lower_bound` by m/z at `mz - mz_tolerance`, then a forward walk
/// that stops above `mz + mz_tolerance` or at the first seed with `|Δm/z| <
/// mz_tolerance` and `|ΔRT| < rt_tolerance`.
///
/// Infinite and NaN differences fail both tolerance tests. With NaN seed m/z
/// values, which [`sorted_user_seeds`] admits only when every other seed m/z
/// is one value `x`, no NaN seed ever matches, and `x` matches only when `x >=
/// mz - mz_tolerance`, where the search starts at the first seed whatever the
/// order; the answer is therefore the same for every order of the seeds. The
/// search follows libstdc++'s probes ([`libstdcxx::lower_bound`]).
fn near_user_seed(user_seeds: &[UserSeed], mz: f64, rt: f64, settings: &Settings) -> bool {
    let tolerance = settings.user_seed_mz_tolerance;
    let start = libstdcxx::lower_bound(user_seeds, |seed| seed.mz < mz - tolerance);
    for seed in &user_seeds[start..] {
        if seed.mz > mz + tolerance {
            break;
        }
        if (seed.mz - mz).abs() < tolerance
            && (seed.rt - rt).abs() < settings.user_seed_rt_tolerance
        {
            return true;
        }
    }
    false
}

/// A scan or window count as a progress value (`SignedSize` in the source).
fn progress_value(value: usize) -> Result<i64> {
    i64::try_from(value).map_err(|_| Error::InvalidValue("progress value overflow".into()))
}

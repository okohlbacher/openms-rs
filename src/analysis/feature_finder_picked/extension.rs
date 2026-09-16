// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Isotope fit and mass-trace extension of the picked feature finder
//! (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`, step 3.3.1 of `run_`).
//!
//! A seed becomes a feature candidate in two steps, both here:
//!
//! 1. [`find_best_isotope_fit`](crate::analysis::feature_finder_picked::extension::find_best_isotope_fit) slides the seed's averagine pattern over the
//!    peaks of the seed's own scan and keeps the placement with the best
//!    isotope score that still contains the seed (source `findBestIsotopeFit_`);
//! 2. [`extend_mass_traces`](crate::analysis::feature_finder_picked::extension::extend_mass_traces) turns each matched isotope into a *mass trace* by
//!    walking the neighbouring scans in both retention-time directions (source
//!    `extendMassTraces_` with `extendMassTrace_`).
//!
//! Both read the per-peak scores of
//! [`ScoreArrays`](crate::analysis::feature_finder_picked::scoring::ScoreArrays)
//! that the seed stage filled; [`OverallScores`](crate::analysis::feature_finder_picked::extension::OverallScores) is the one array the extension
//! needs, the overall score of the charge being extended.
//!
//! Arithmetic follows the source operation by operation, including where it
//! reads a `float` into a `double` expression and where an unsigned difference
//! wraps. Two source defects are reproduced deliberately and are listed in
//! `docs/FEATURE_FINDER_PICKED_SUPPORT.md` (B7 candidates 1 and 2): the
//! trace-index comparison of [`extend_mass_traces`](crate::analysis::feature_finder_picked::extension::extend_mass_traces) and the better-seed search
//! that re-reads the moving seed's m/z while comparing against the original
//! one.
//!
//! The module is serial: the source parallelises one level above it, over
//! seeds.

use crate::analysis::feature_finder_picked::algorithm::Settings;
use crate::analysis::feature_finder_picked::debug::{LogSink, NoLog, g, g32, put_all};
use crate::analysis::feature_finder_picked::helper_structs::{
    IsotopePattern, MassTrace, MassTraces, PatternPeak, Seed, TracePeak,
};
use crate::analysis::feature_finder_picked::scoring::{
    ScoreArrays, find_isotope_logged, isotope_score_logged, position_score, reset_pattern,
};
use crate::analysis::feature_finder_picked::seeds::IsotopeWindows;
use crate::format::file_info::text_format::to_str;
use crate::kernel::{MSSpectrum, nearest};
use crate::{Error, Result};

/// The overall-score array of one charge, the only score array the extension
/// reads: source `map_[s].getFloatDataArrays()[meta_index_overall]`.
///
/// The source keeps the array on each spectrum; this port keeps it in
/// [`ScoreArrays`], so the view pairs the arrays with the charge index.
#[derive(Clone, Copy, Debug)]
pub struct OverallScores<'a> {
    scores: &'a ScoreArrays,
    charge_index: usize,
}

impl<'a> OverallScores<'a> {
    /// The overall scores of charge index `charge_index` (`charge -
    /// charge_low`).
    pub fn new(scores: &'a ScoreArrays, charge_index: usize) -> Self {
        Self {
            scores,
            charge_index,
        }
    }

    /// The overall score of one peak, or `None` outside the arrays.
    ///
    /// The source indexes both arrays without a check; every call below stays
    /// in range, and an out-of-range index becomes [`Error::InvalidValue`]
    /// rather than undefined behaviour.
    pub fn get(&self, spectrum: usize, peak: usize) -> Option<f32> {
        self.scores
            .overall(self.charge_index, spectrum)?
            .get(peak)
            .copied()
    }

    /// The charge index this view reads.
    pub fn charge_index(&self) -> usize {
        self.charge_index
    }
}

fn out_of_range(what: &str) -> Error {
    Error::InvalidValue(format!(
        "FeatureFinderAlgorithmPicked seed extension: {what} is out of range"
    ))
}

/// Whether `pattern` matched the seed's own peak: the source's `seed_contained`
/// test, `pattern.peak[iso] == center.peak && pattern.spectrum[iso] ==
/// center.spectrum`.
fn contains_seed(pattern: &IsotopePattern, seed: Seed) -> bool {
    pattern
        .peak
        .iter()
        .zip(&pattern.spectrum)
        .any(|(peak, &spectrum)| {
            *peak == PatternPeak::Found(seed.peak) && spectrum == seed.spectrum
        })
}

/// The best averagine placement containing `seed`: source
/// `findBestIsotopeFit_`.
///
/// The theoretical pattern is the precalculated window of `seed_mz * charge`.
/// Its m/z reach, `(isotopes + 1) / charge`, bounds a linear search around the
/// seed peak in the seed's own scan: every peak of that window is tried as the
/// pattern's first isotope, its isotopes are matched with
/// [`find_isotope`] (which also inspects the neighbouring scans), and the
/// placement is scored with [`isotope_score`] *without* the m/z-distance factor.
/// A placement counts only while it matches the seed peak itself, checked both
/// before scoring and again afterwards, because the scoring may drop optional
/// isotopes. The highest score wins, compared with a strict `>`, so the
/// earliest of equal scores is kept.
///
/// Returns the best score, `0.0` when no placement qualified, and the matching
/// pattern, whose `theoretical_pattern` is the window's pattern in either case,
/// as the source assigns it unconditionally at the end.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `seed` does not address a peak of
/// `spectra`, when the window for `seed_mz * charge` was not precalculated
/// ([`IsotopeWindows::get`]), and from the helpers. The source indexes without
/// those checks.
pub fn find_best_isotope_fit(
    spectra: &[MSSpectrum],
    windows: &IsotopeWindows,
    settings: &Settings,
    seed: Seed,
    charge: i32,
) -> Result<(f64, IsotopePattern)> {
    find_best_isotope_fit_logged(spectra, windows, settings, seed, charge, &mut NoLog)
}

/// [`find_best_isotope_fit`] writing the source's debug lines to `log`.
pub(crate) fn find_best_isotope_fit_logged<L: LogSink>(
    spectra: &[MSSpectrum],
    windows: &IsotopeWindows,
    settings: &Settings,
    seed: Seed,
    charge: i32,
    log: &mut L,
) -> Result<(f64, IsotopePattern)> {
    if log.enabled() {
        put_all(
            log,
            &[
                "Testing isotope patterns for charge ",
                &charge.to_string(),
                ": \n",
            ],
        );
    }
    let spectrum = spectra
        .get(seed.spectrum)
        .ok_or_else(|| out_of_range("the seed's spectrum"))?;
    let center_mz = spectrum
        .peaks
        .get(seed.peak)
        .ok_or_else(|| out_of_range("the seed's peak"))?
        .mz;
    let c = f64::from(charge);
    let isotopes = windows.get(center_mz * c)?;
    if log.enabled() {
        put_all(
            log,
            &[
                " - Seed: ",
                &seed.peak.to_string(),
                " (mz:",
                &g(center_mz),
                ")\n",
            ],
        );
    }
    let size = isotopes.len();
    let mass_window = (size + 1) as f64 / c;
    if log.enabled() {
        put_all(log, &[" - Mass window: ", &g(mass_window), "\n"]);
    }

    // Source: `end` walks up while the peak stays below the window, then steps
    // back by one. The first test always succeeds (`center_mz < center_mz +
    // mass_window` for a positive window), so `end` never underflows.
    let mut end = seed.peak;
    while end < spectrum.peaks.len() && spectrum.peaks[end].mz < center_mz + mass_window {
        end += 1;
    }
    let end = end
        .checked_sub(1)
        .ok_or_else(|| out_of_range("the isotope search window"))?;
    // Source: `begin` is a `SignedSize` that may reach -1 before stepping back up.
    let mut begin = seed.peak as i64;
    while begin >= 0 && spectrum.peaks[begin as usize].mz > center_mz - mass_window {
        begin -= 1;
    }
    let begin = (begin + 1) as usize;
    if log.enabled() {
        for (name, index) in [(" - Begin: ", begin), (" - End: ", end)] {
            put_all(
                log,
                &[
                    name,
                    &index.to_string(),
                    " (mz:",
                    &g(spectrum.peaks[index].mz),
                    ")\n",
                ],
            );
        }
    }

    let mut max_score = 0.0;
    let mut best = IsotopePattern::new(0)?;
    let mut pattern = IsotopePattern::new(size)?;
    for start in begin..=end {
        let start_mz = spectrum.peaks[start].mz;
        let mut peak_index = start;
        reset_pattern(&mut pattern, size);
        if log.enabled() {
            put_all(
                log,
                &[
                    " - Fitting at ",
                    &start.to_string(),
                    " (mz:",
                    &g(start_mz),
                    ")\n",
                ],
            );
        }
        for iso in 0..size {
            let pos = start_mz + iso as f64 / c;
            find_isotope_logged(
                spectra,
                pos,
                seed.spectrum,
                &mut pattern,
                iso,
                &mut peak_index,
                settings.pattern_tolerance,
                log,
            )?;
        }
        if !contains_seed(&pattern, seed) {
            put_all(log, &["   - aborting: seed is not contained!\n"]);
            continue;
        }
        let (score, _) = isotope_score_logged(
            isotopes,
            &mut pattern,
            false,
            settings.min_isotope_fit,
            settings.optional_fit_improvement,
            log,
        )?;
        if !contains_seed(&pattern, seed) {
            put_all(
                log,
                &["   - aborting: seed was removed during isotope fit!\n"],
            );
            continue;
        }
        if log.enabled() {
            put_all(log, &["   - final score: ", &g(score), "\n"]);
        }
        if score > max_score {
            max_score = score;
            best.peak.clone_from(&pattern.peak);
            best.spectrum.clone_from(&pattern.spectrum);
            best.intensity.clone_from(&pattern.intensity);
            best.mz_score.clone_from(&pattern.mz_score);
            best.theoretical_mz.clone_from(&pattern.theoretical_mz);
        }
    }
    if log.enabled() {
        put_all(log, &[" - best score              : ", &g(max_score), "\n"]);
    }
    best.theoretical_pattern = isotopes.clone();
    Ok((max_score, best))
}

/// The mass traces of one feature candidate: source `extendMassTraces_`.
///
/// The most intense matched isotope starts the *maximum trace*, which is
/// extended without retention-time bounds; its first and last retention time
/// then bound every other trace. A candidate is dropped, by returning no trace
/// at all, when that maximum trace has fewer than three peaks or fewer than
/// `2 * min_spectra - max_missing` peaks.
///
/// Every other matched isotope first looks for a stronger start peak within
/// `min_spectra` scans (see below), is then extended inside the maximum trace's
/// bounds, and is appended. [`MassTraces::max_trace`] is the position of the
/// maximum trace in the result.
///
/// # Preserved source defect
///
/// An isotope whose trace stays below three peaks is handled by comparing the
/// *pattern* index `p` with [`MassTraces::max_trace`], which is an index into
/// the traces collected so far and is still `0` before the maximum trace is
/// reached. `p < max_trace` clears the traces collected so far and skips the
/// isotope, `p > max_trace` stops the extension, and `p == max_trace` falls
/// through and appends the invalid trace. Before the maximum trace is reached
/// that means an invalid trace at `p == 0` is kept and any later invalid trace
/// stops the extension, which is not what the comment ("Missing traces in the
/// middle of a pattern are not acceptable") describes. The behaviour is
/// reproduced; see B7 candidate 1 in `docs/FEATURE_FINDER_PICKED_SUPPORT.md`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a matched isotope does not address a
/// peak of `spectra`, when `pattern` has no matched isotope at all (the source
/// then dereferences the pattern's first entry, which is undefined behaviour;
/// it is unreachable from `run_`, where the pattern always contains the seed),
/// when the theoretical pattern is shorter than the matched pattern, and from
/// [`extend_mass_trace`].
pub fn extend_mass_traces(
    spectra: &[MSSpectrum],
    overall: OverallScores<'_>,
    settings: &Settings,
    pattern: &IsotopePattern,
) -> Result<MassTraces> {
    extend_mass_traces_logged(spectra, overall, settings, pattern, &mut NoLog)
}

/// The source's `" - extending from: RT / mz (int: I)"`-style line of a trace
/// start.
fn log_trace_point<L: LogSink>(log: &mut L, label: &str, peak: &TracePeak) {
    if log.enabled() {
        put_all(
            log,
            &[
                label,
                &g(peak.rt),
                " / ",
                &g(peak.mz),
                " (int: ",
                &g32(peak.intensity),
                ")\n",
            ],
        );
    }
}

/// [`extend_mass_traces`] writing the source's debug lines to `log`.
pub(crate) fn extend_mass_traces_logged<L: LogSink>(
    spectra: &[MSSpectrum],
    overall: OverallScores<'_>,
    settings: &Settings,
    pattern: &IsotopePattern,
    log: &mut L,
) -> Result<MassTraces> {
    let peak_of = |p: usize| -> Result<Option<TracePeak>> {
        let PatternPeak::Found(peak) = pattern.peak[p] else {
            return Ok(None);
        };
        let spectrum_index = pattern.spectrum[p];
        let spectrum = spectra
            .get(spectrum_index)
            .ok_or_else(|| out_of_range("a matched isotope's spectrum"))?;
        let found = spectrum
            .peaks
            .get(peak)
            .ok_or_else(|| out_of_range("a matched isotope's peak"))?;
        Ok(Some(TracePeak::new(
            spectrum_index,
            peak,
            spectrum.rt,
            found.mz,
            found.intensity,
        )))
    };

    // Source: the strongest matched isotope, compared as `float > double` with
    // the running maximum starting at 0.0, so a zero-intensity peak never wins.
    let mut max_int = 0.0;
    let mut max_trace_index = 0usize;
    for p in 0..pattern.peak.len() {
        let Some(found) = peak_of(p)? else { continue };
        if f64::from(found.intensity) > max_int {
            max_int = f64::from(found.intensity);
            max_trace_index = p;
        }
    }
    let matched = if pattern.peak.is_empty() {
        None
    } else {
        peak_of(max_trace_index)?
    };
    let Some(start) = matched else {
        return Err(Error::InvalidValue(
            "FeatureFinderAlgorithmPicked seed extension: the isotope pattern matched no peak; \
             the source reads its first entry here"
                .into(),
        ));
    };
    if pattern.theoretical_pattern.len() < pattern.peak.len() {
        return Err(Error::InvalidValue(
            "FeatureFinderAlgorithmPicked seed extension: the theoretical pattern is shorter \
             than the matched pattern"
                .into(),
        ));
    }

    if log.enabled() {
        put_all(
            log,
            &[
                " - Trace ",
                &max_trace_index.to_string(),
                " (maximum intensity)\n",
            ],
        );
    }
    log_trace_point(log, "   - extending from: ", &start);
    let mut max_trace = MassTrace::default();
    max_trace.peaks.push(start);
    extend_mass_trace_logged(
        &mut max_trace,
        spectra,
        overall,
        start.spectrum,
        start.mz,
        false,
        settings,
        0.0,
        0.0,
        log,
    )?;
    extend_mass_trace_logged(
        &mut max_trace,
        spectra,
        overall,
        start.spectrum,
        start.mz,
        true,
        settings,
        0.0,
        0.0,
        log,
    )?;
    // Both are set, because the trace holds the start peak.
    let rt_min = max_trace.peaks[0].rt;
    let rt_max = max_trace.peaks[max_trace.peaks.len() - 1].rt;
    if log.enabled() {
        put_all(
            log,
            &["   - rt bounds: ", &g(rt_min), "-", &g(rt_max), "\n"],
        );
    }

    let mut traces = MassTraces::new();
    traces.reserve(pattern.peak.len())?;
    // Source: `2 * min_spectra_ - max_missing_trace_peaks_` in `UInt`, which
    // wraps when more missing peaks than twice the half-window are allowed; the
    // comparison against a `Size` then always holds and the candidate is
    // dropped.
    let min_spectra = u32::try_from(settings.min_spectra).unwrap_or(u32::MAX);
    let required = u64::from(
        min_spectra
            .wrapping_mul(2)
            .wrapping_sub(settings.max_missing_trace_peaks),
    );
    if !max_trace.is_valid() || (max_trace.peaks.len() as u64) < required {
        put_all(
            log,
            &["   - could not extend trace with maximum intensity => abort\n"],
        );
        return Ok(traces);
    }

    let mut max_trace = Some(max_trace);
    for p in 0..pattern.peak.len() {
        if log.enabled() {
            put_all(log, &[" - Trace ", &p.to_string(), "\n"]);
        }
        if p == max_trace_index {
            put_all(log, &["   - previously extended maximum trace\n"]);
            let mut trace = max_trace.take().ok_or_else(|| {
                Error::InvalidValue("the maximum trace was already consumed".into())
            })?;
            trace.theoretical_int = pattern.theoretical_pattern.intensity[p];
            traces.push(trace);
            traces.max_trace = traces.len() - 1;
            continue;
        }
        // Source: `-2` (removed during the isotope fit) and `-1` (missing) are
        // both skipped.
        match pattern.peak[p] {
            PatternPeak::Removed => {
                put_all(log, &["   - removed during isotope fit\n"]);
                continue;
            }
            PatternPeak::NotFound => {
                put_all(log, &["   - missing\n"]);
                continue;
            }
            PatternPeak::Found(_) => {}
        }
        let Some(found) = peak_of(p)? else { continue };
        log_trace_point(log, "   - trace seed: ", &found);

        // Source: look for a stronger start peak in the surrounding
        // `min_spectra` scans. `begin` is an unsigned difference that wraps
        // when the isotope sits in one of the first `min_spectra` scans, which
        // makes the range empty and skips the search; `end` is exclusive and
        // clamped to the number of scans.
        let mut seed_spectrum = found.spectrum;
        let mut seed_peak = found.peak;
        let begin = seed_spectrum.wrapping_sub(settings.min_spectra);
        let end = seed_spectrum
            .saturating_add(settings.min_spectra)
            .min(spectra.len());
        let mz = found.mz;
        let mut intensity = f64::from(found.intensity);
        for spectrum_index in begin..end {
            let spectrum = &spectra[spectrum_index];
            if spectrum.peaks.is_empty() {
                continue;
            }
            // Source: the search m/z is re-read from the *current* start peak,
            // which moves, while the tolerance below compares against the
            // original `mz`.
            let current = spectra[seed_spectrum].peaks[seed_peak].mz;
            let Some(peak_index) = nearest(&spectrum.peaks, current, |peak| peak.mz) else {
                continue;
            };
            let candidate = spectrum.peaks[peak_index];
            if f64::from(candidate.intensity) <= intensity
                || (mz - candidate.mz).abs() >= settings.pattern_tolerance
            {
                continue;
            }
            seed_spectrum = spectrum_index;
            seed_peak = peak_index;
            intensity = f64::from(candidate.intensity);
        }

        let spectrum = &spectra[seed_spectrum];
        let peak = spectrum.peaks[seed_peak];
        let start = TracePeak::new(
            seed_spectrum,
            seed_peak,
            spectrum.rt,
            peak.mz,
            peak.intensity,
        );
        log_trace_point(log, "   - extending from: ", &start);
        let mut trace = MassTrace::default();
        trace.peaks.push(start);
        extend_mass_trace_logged(
            &mut trace,
            spectra,
            overall,
            seed_spectrum,
            peak.mz,
            false,
            settings,
            rt_min,
            rt_max,
            log,
        )?;
        extend_mass_trace_logged(
            &mut trace,
            spectra,
            overall,
            seed_spectrum,
            peak.mz,
            true,
            settings,
            rt_min,
            rt_max,
            log,
        )?;

        if !trace.is_valid() {
            put_all(log, &["   - could not extend trace \n"]);
            if p < traces.max_trace {
                // Source: clears the traces collected so far; `max_trace` keeps
                // its value and is then stale.
                traces.clear();
                continue;
            } else if p > traces.max_trace {
                break;
            }
            // p == traces.max_trace: the source falls through and appends the
            // invalid trace (see the preserved defect above).
        }
        trace.theoretical_int = pattern.theoretical_pattern.intensity[p];
        traces.push(trace);
    }
    Ok(traces)
}

/// Extend one mass trace into the neighbouring scans: source
/// `extendMassTrace_`.
///
/// `spectrum_index` is the scan the trace starts in; the walk begins at the
/// next scan below (`increase_rt` false) or above it. The upward call first
/// reverses the peaks collected so far, so the trace ends in chronological
/// order. In each scan the peak nearest to `mz` is taken unless it is missing,
/// its overall score is below `0.01`, or its m/z score against
/// `mass_trace:mz_tolerance` is zero; more than `mass_trace:max_missing`
/// consecutive such scans stop the walk.
///
/// After every accepted peak the mean of the last `min_spectra` relative
/// intensity changes is compared with the slope bound, which is doubled when
/// retention-time bounds are given (`max_rt != min_rt`). Exceeding it removes
/// the last `min_spectra - 1` peaks added by *this* call and stops the walk.
///
/// # Arguments
///
/// * `overall` - overall scores of the charge being extended.
/// * `mz` - m/z the trace follows, the start peak's.
/// * `min_rt`, `max_rt` - retention-time bounds in seconds; equal values (the
///   source's `0, 0` default) mean no bounds, as the source decides.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the trace is empty, when
/// `mass_trace:min_spectra` is below 2 so that the source's delta buffer has
/// `size_t(-1)` entries (unreachable, because such a run finds no seed), when a
/// peak's overall score is missing, and when the delta buffer cannot be
/// allocated.
#[allow(clippy::too_many_arguments)]
pub fn extend_mass_trace(
    trace: &mut MassTrace,
    spectra: &[MSSpectrum],
    overall: OverallScores<'_>,
    spectrum_index: usize,
    mz: f64,
    increase_rt: bool,
    settings: &Settings,
    min_rt: f64,
    max_rt: f64,
) -> Result<()> {
    extend_mass_trace_logged(
        trace,
        spectra,
        overall,
        spectrum_index,
        mz,
        increase_rt,
        settings,
        min_rt,
        max_rt,
        &mut NoLog,
    )
}

/// [`extend_mass_trace`] writing the source's `   - Added <n> peaks (abort:
/// <reason>)` line to `log`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn extend_mass_trace_logged<L: LogSink>(
    trace: &mut MassTrace,
    spectra: &[MSSpectrum],
    overall: OverallScores<'_>,
    spectrum_index: usize,
    mz: f64,
    increase_rt: bool,
    settings: &Settings,
    min_rt: f64,
    max_rt: f64,
    log: &mut L,
) -> Result<()> {
    let mut index = spectrum_index as i64;
    if increase_rt {
        index += 1;
        trace.peaks.reverse();
    } else {
        index -= 1;
    }
    let boundaries = max_rt != min_rt;
    let current_slope_bound = (1.0 + f64::from(u8::from(boundaries))) * settings.slope_bound;
    let delta_count = settings.min_spectra;
    let pad = delta_count.checked_sub(1).ok_or_else(|| {
        Error::InvalidValue(
            "FeatureFinderAlgorithmPicked seed extension: mass_trace:min_spectra below 2 leaves \
             no scan on either side; the source allocates size_t(-1) intensity deltas here"
                .into(),
        )
    })?;
    let mut deltas: Vec<f64> = Vec::new();
    deltas
        .try_reserve(pad + 1)
        .map_err(|_| Error::InvalidValue("cannot allocate the intensity deltas".into()))?;
    deltas.resize(pad, 0.0);
    let mut last_observed_intensity = f64::from(
        trace
            .peaks
            .last()
            .ok_or_else(|| {
                Error::InvalidValue(
                    "FeatureFinderAlgorithmPicked seed extension: an empty trace has no last \
                     intensity"
                        .into(),
                )
            })?
            .intensity,
    );
    let mut missing_peaks: u32 = 0;
    let peaks_before_extension = trace.peaks.len();
    // The source's `abort_reason`, empty when the walk runs out of scans.
    let mut abort_reason = String::new();

    while (!increase_rt && index >= 0) || (increase_rt && index < spectra.len() as i64) {
        let s = index as usize;
        let spectrum = &spectra[s];
        if boundaries
            && ((!increase_rt && spectrum.rt < min_rt) || (increase_rt && spectrum.rt > max_rt))
        {
            if log.enabled() {
                abort_reason = "Hit upper/lower boundary".to_owned();
            }
            break;
        }
        let found = if spectrum.peaks.is_empty() {
            None
        } else {
            nearest(&spectrum.peaks, mz, |peak| peak.mz)
        };
        // Source: a peak is "missing" when none was found, its overall score is
        // below 0.01, or its m/z score is exactly zero. All three comparisons
        // are false for NaN, which therefore counts as present, as in the
        // source.
        let accepted = match found {
            None => None,
            Some(peak_index) => {
                let score = overall
                    .get(s, peak_index)
                    .ok_or_else(|| out_of_range("an overall score"))?;
                let peak = spectrum.peaks[peak_index];
                if f64::from(score) < 0.01
                    || position_score(mz, peak.mz, settings.trace_tolerance) == 0.0
                {
                    None
                } else {
                    Some((peak_index, peak))
                }
            }
        };
        match accepted {
            None => {
                missing_peaks += 1;
                if missing_peaks > settings.max_missing_trace_peaks {
                    if log.enabled() {
                        abort_reason = "too many peaks missing".to_owned();
                    }
                    break;
                }
            }
            Some((peak_index, peak)) => {
                missing_peaks = 0;
                trace.peaks.push(TracePeak::new(
                    s,
                    peak_index,
                    spectrum.rt,
                    peak.mz,
                    peak.intensity,
                ));
                let intensity = f64::from(peak.intensity);
                deltas.push((intensity - last_observed_intensity) / last_observed_intensity);
                last_observed_intensity = intensity;
                // Source: the mean of the last `delta_count` deltas, summed
                // left to right from 0.0.
                let start = deltas.len() - delta_count;
                let average_delta =
                    deltas[start..].iter().fold(0.0, |sum, d| sum + d) / delta_count as f64;
                if average_delta > current_slope_bound {
                    if log.enabled() {
                        // Source: `std::string + double`, `StringUtils::toStr`.
                        abort_reason = format!(
                            "Average delta above threshold: {}/{}",
                            to_str(average_delta),
                            to_str(current_slope_bound)
                        );
                    }
                    let remove = (trace.peaks.len() - peaks_before_extension).min(delta_count - 1);
                    trace.peaks.truncate(trace.peaks.len() - remove);
                    break;
                }
            }
        }
        index += if increase_rt { 1 } else { -1 };
    }
    if log.enabled() {
        put_all(
            log,
            &[
                "   - Added ",
                &(trace.peaks.len() - peaks_before_extension).to_string(),
                " peaks (abort: ",
                &abort_reason,
                ")\n",
            ],
        );
    }
    Ok(())
}

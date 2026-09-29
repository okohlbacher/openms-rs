// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Correlation of individual mass traces into pseudo-spectra.
//!
//! Ported from `ANALYSIS/OPENSWATH/MasstraceCorrelator.h` and `.cpp` at
//! OpenMS4-core `bc9cc12`, the algorithm behind the `ClusterMassTraces` and
//! `ClusterMassTracesByPrecursor` tools. Mass traces are the (RT, intensity)
//! points of the consensus features of one map; each unused trace seeds a
//! pseudo-spectrum, and every later (lower-intensity) trace that co-elutes
//! with it (Pearson correlation and cross-correlation lag) joins it. See
//! `docs/MRM_TRANSITION_GROUP_PICKER_SUPPORT.md` (section "Mass trace
//! correlation").

use crate::analysis::openswath::scoring::{normalized_cross_correlation, xcorr_array_get_max_peak};
use crate::concept::progress_logger::{ProgressLogger, ProgressReporter, progress_value};
use crate::kernel::{ConsensusMap, MSSpectrum, Peak1D};
use crate::math::statistic_functions::pearson_correlation_coefficient;
use crate::param::{DefaultParamHandler, Param, ParamNode};
use crate::processing::peak_picking::{int_entry, integer};
use crate::{Error, Result};

/// The handler name the source constructor passes,
/// `DefaultParamHandler("MRMFeatureFinderScoring")` — not the class name, a
/// source slip reproduced here because unknown-parameter warnings carry it.
pub const MASSTRACE_CORRELATOR_NAME: &str = "MRMFeatureFinderScoring";

/// Label of the cache-building progress section.
pub const CACHE_PROGRESS_LABEL: &str = "create consensus map cache";

/// Label of the correlation progress section (with the source's trailing
/// space).
pub const CORRELATION_PROGRESS_LABEL: &str = "correlating masstraces ";

/// A mass trace: `(RT, intensity)` points (source typedef
/// `MasstracePointsType`).
pub type MasstracePoints = Vec<(f64, f64)>;

/// The scores of one pair of mass traces (source `scoreHullpoints`'
/// out-parameters).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HullpointScores {
    /// Pearson correlation of the matched intensities; NaN when a trace is
    /// constant.
    pub pearson_score: f64,
    /// The cross-correlation lag at its maximum (the co-elution score) and the
    /// correlation there (the shape score), computed only when the Pearson
    /// score is above the minimum; the source then leaves its `lag` and
    /// `lag_intensity` out-parameters untouched.
    pub lag: Option<(i32, f64)>,
}

/// The per-feature data [`MasstraceCorrelator::create_pseudo_spectra`] works
/// on (source `createConsensusMapCache`'s three out-vectors).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ConsensusMapCache {
    /// Each feature's mass trace, sorted by RT.
    pub feature_points: Vec<MasstracePoints>,
    /// Each feature's most intense handle as `(m/z, intensity)` (the source
    /// names the vector `max_intensities` and fills it in this order).
    pub max_intensities: Vec<(f64, f64)>,
    /// Each feature's retention time.
    pub rt_cache: Vec<f64>,
}

/// The scores with which one mass trace joined a pseudo-spectrum.
///
/// The source writes these as five float data arrays named `RT_apex`,
/// `RT_diff`, `lag`, `pearson_score` and `lag_intensity` on the spectrum, one
/// entry per joined trace: the seed peak has none, so the arrays are one
/// shorter than the peaks. A native spectrum's arrays must match its peaks,
/// so the port returns them typed, in `f32` as the source stores them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TraceCorrelation {
    /// `RT_apex`: the joined feature's retention time.
    pub rt_apex: f32,
    /// `RT_diff`: its distance from the seed's retention time.
    pub rt_diff: f32,
    /// `lag`: the cross-correlation lag.
    pub lag: f32,
    /// `pearson_score`.
    pub pearson_score: f32,
    /// `lag_intensity`: the cross-correlation at the lag.
    pub lag_intensity: f32,
}

/// One pseudo-spectrum: the seed trace and the traces that joined it.
#[derive(Clone, Debug, PartialEq)]
pub struct PseudoSpectrum {
    /// MS level 2 at the seed feature's RT; peak 0 is the seed, peak `k + 1`
    /// the `k`-th joined trace, each at its most intense handle's m/z and
    /// intensity.
    pub spectrum: MSSpectrum,
    /// One entry per joined trace (the source's float arrays; see
    /// [`TraceCorrelation`]).
    pub correlations: Vec<TraceCorrelation>,
}

/// Correlates individual mass traces found in mass spectrometric maps (source
/// `MasstraceCorrelator`), using the normalised cross-correlation and Pearson
/// scoring of the OpenSWATH module.
///
/// The source class is a `DefaultParamHandler` with three parameters that
/// nothing reads (it has no `updateMembers_`); they are kept for round trips.
/// Its `ProgressLogger` base is the `_with_progress` methods.
#[derive(Clone, Debug)]
pub struct MasstraceCorrelator {
    /// Source `sgolay_frame_length` (default 15); unused, as in the source.
    pub sgolay_frame_length: i64,
    /// Source `sgolay_polynomial_order` (default 3); unused, as in the source.
    pub sgolay_polynomial_order: i64,
    /// Source `gauss_width` (an integer parameter, default 50); unused, as in
    /// the source.
    pub gauss_width: i64,
    /// Native budget for [`create_pseudo_spectra`](Self::create_pseudo_spectra):
    /// matched points plus squared cross-correlation length, summed over the
    /// scored pairs; the source has none.
    pub max_work: u64,
}

impl Default for MasstraceCorrelator {
    fn default() -> Self {
        Self {
            sgolay_frame_length: 15,
            sgolay_polynomial_order: 3,
            gauss_width: 50,
            max_work: 20_000_000_000,
        }
    }
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

impl MasstraceCorrelator {
    /// The source defaults (`getDefaults`): `sgolay_frame_length` (15),
    /// `sgolay_polynomial_order` (3) and `gauss_width` (the integer 50), with
    /// the source descriptions.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The current parameters (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: vec![
                int_entry(
                    "sgolay_frame_length",
                    self.sgolay_frame_length,
                    "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added.",
                    false,
                    None,
                    None,
                ),
                int_entry(
                    "sgolay_polynomial_order",
                    self.sgolay_polynomial_order,
                    "Order or the polynomial that is fitted.",
                    false,
                    None,
                    None,
                ),
                int_entry(
                    "gauss_width",
                    self.gauss_width,
                    "Gaussian width.",
                    false,
                    None,
                    None,
                ),
            ],
            nodes: Vec::new(),
        })
    }

    /// Build a correlator from parameters (source `setParameters`); the values
    /// are stored and, as in the source, used for nothing.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a wrong value type.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`MasstraceCorrelator::from_param`], also returning the
    /// unknown-parameter warnings, under the source's handler name
    /// [`MASSTRACE_CORRELATOR_NAME`].
    ///
    /// # Errors
    ///
    /// As [`MasstraceCorrelator::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(MASSTRACE_CORRELATOR_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |param| {
            Ok(Self {
                sgolay_frame_length: integer(param, "sgolay_frame_length")?,
                sgolay_polynomial_order: integer(param, "sgolay_polynomial_order")?,
                gauss_width: integer(param, "gauss_width")?,
                ..Self::default()
            })
        })
    }

    /// Score two mass traces against each other (source `scoreHullpoints`).
    ///
    /// The traces are first matched in RT (points less than `mindiff` apart
    /// pair up; an unmatched point pairs with a zero, and the longer trace's
    /// ends are padded with zeros), then the matched intensities are scored:
    /// always with the Pearson correlation, and, when it is above
    /// `min_corr`, with the normalised cross-correlation, whose maximum gives
    /// the lag (co-elution) and the correlation there (shape). `max_lag` is
    /// unused, as in the source (the caller applies it). Both traces must be
    /// sorted by RT (the source's `@note`; unchecked here as there).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a non-finite RT (the source's matching loop
    /// then never advances) and [`Error::InvalidRange`] when both traces are
    /// empty (the source's Pearson throws `Exception::InvalidRange`).
    pub fn score_hullpoints(
        &self,
        hull_points1: &[(f64, f64)],
        hull_points2: &[(f64, f64)],
        min_corr: f64,
        _max_lag: i32,
        mindiff: f64,
    ) -> Result<HullpointScores> {
        let (mut first, mut second) = match_mass_traces(hull_points1, hull_points2, mindiff, true)?;
        let pearson_score = pearson_correlation_coefficient(&first, &second)?;
        if pearson_score <= min_corr {
            return Ok(HullpointScores {
                pearson_score,
                lag: None,
            });
        }
        let maxdelay = i32::try_from(first.len())
            .map_err(|_| invalid("mass trace too long for cross-correlation"))?;
        let xcorr = normalized_cross_correlation(&mut first, &mut second, maxdelay, 1)?;
        let best = xcorr.data[xcorr_array_get_max_peak(&xcorr)?];
        Ok(HullpointScores {
            pearson_score,
            lag: Some(best),
        })
    }

    /// Cache the mass traces of a consensus map (source
    /// `createConsensusMapCache`): each feature's handles as `(RT, intensity)`
    /// sorted by RT, its most intense handle's `(m/z, intensity)` (the first
    /// on ties; `(0, 0)` without a positive intensity) and its RT.
    ///
    /// Equal RTs keep the handle order (a stable sort); the source's
    /// `std::sort` is unstable, which for traces longer than 16 points may
    /// order equal RTs differently.
    ///
    /// # Errors
    ///
    /// None in practice; the signature matches
    /// [`MasstraceCorrelator::create_consensus_map_cache_with_progress`].
    pub fn create_consensus_map_cache(&self, map: &ConsensusMap) -> Result<ConsensusMapCache> {
        self.cache_reporting(map, &mut ProgressReporter::silent())
    }

    /// [`MasstraceCorrelator::create_consensus_map_cache`], reporting progress
    /// as the source's `ProgressLogger` base does: a section labelled
    /// [`CACHE_PROGRESS_LABEL`] from 0 to the feature count, `setProgress(i)`
    /// before feature `i`.
    ///
    /// # Errors
    ///
    /// The errors of `progress`.
    pub fn create_consensus_map_cache_with_progress(
        &self,
        map: &ConsensusMap,
        progress: &mut ProgressLogger,
    ) -> Result<ConsensusMapCache> {
        self.cache_reporting(map, &mut ProgressReporter::new(Some(progress)))
    }

    fn cache_reporting(
        &self,
        map: &ConsensusMap,
        reporter: &mut ProgressReporter<'_>,
    ) -> Result<ConsensusMapCache> {
        let count = progress_value(map.features.len())?;
        reporter.section(0, count, CACHE_PROGRESS_LABEL, |reporter| {
            let mut cache = ConsensusMapCache::default();
            for (i, feature) in map.features.iter().enumerate() {
                reporter.set_count(i)?;
                let mut points: MasstracePoints = feature
                    .handles()
                    .iter()
                    .map(|h| (h.rt, f64::from(h.intensity)))
                    .collect();
                points.sort_by(|a, b| {
                    if a.0 < b.0 {
                        std::cmp::Ordering::Less
                    } else if b.0 < a.0 {
                        std::cmp::Ordering::Greater
                    } else {
                        std::cmp::Ordering::Equal
                    }
                });
                cache.feature_points.push(points);
                let (mut max_int, mut max_mz) = (0.0, 0.0);
                for handle in feature.handles() {
                    if f64::from(handle.intensity) > max_int {
                        max_int = f64::from(handle.intensity);
                        max_mz = handle.mz;
                    }
                }
                cache.max_intensities.push((max_mz, max_int));
                cache.rt_cache.push(feature.base.rt);
            }
            Ok(cache)
        })
    }

    /// Compute pseudo-spectra from the mass traces of a consensus map (source
    /// `createPseudoSpectra`).
    ///
    /// All features are compared all-vs-all: each feature not yet used seeds a
    /// pseudo-spectrum, and every later feature whose RT lies within
    /// `max_rt_apex_difference` of the seed's and inside the seed's trace, and
    /// whose trace correlates above `min_correlation` with a lag within
    /// `max_lag`, joins it and can no longer seed one. A pseudo-spectrum with
    /// more than `min_peak_nr` peaks is kept. The map is assumed to hold the
    /// features of one (SWATH) map, sorted by intensity.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a seed feature has no handles (the source
    /// reads the first point of an empty trace) or the work budget
    /// [`max_work`](Self::max_work) is exhausted, and the errors of
    /// [`MasstraceCorrelator::score_hullpoints`]. Nothing is returned on
    /// error.
    pub fn create_pseudo_spectra(
        &self,
        map: &ConsensusMap,
        min_peak_nr: usize,
        min_correlation: f64,
        max_lag: i32,
        max_rt_apex_difference: f64,
    ) -> Result<Vec<PseudoSpectrum>> {
        self.pseudo_spectra_reporting(
            map,
            min_peak_nr,
            min_correlation,
            max_lag,
            max_rt_apex_difference,
            &mut ProgressReporter::silent(),
        )
    }

    /// [`MasstraceCorrelator::create_pseudo_spectra`], reporting progress as
    /// the source does: the cache section, then a section labelled
    /// [`CORRELATION_PROGRESS_LABEL`] from 0 to the feature count with
    /// `setProgress(i)` before seed candidate `i`.
    ///
    /// # Errors
    ///
    /// As [`MasstraceCorrelator::create_pseudo_spectra`], and the errors of
    /// `progress`.
    pub fn create_pseudo_spectra_with_progress(
        &self,
        map: &ConsensusMap,
        min_peak_nr: usize,
        min_correlation: f64,
        max_lag: i32,
        max_rt_apex_difference: f64,
        progress: &mut ProgressLogger,
    ) -> Result<Vec<PseudoSpectrum>> {
        self.pseudo_spectra_reporting(
            map,
            min_peak_nr,
            min_correlation,
            max_lag,
            max_rt_apex_difference,
            &mut ProgressReporter::new(Some(progress)),
        )
    }

    fn pseudo_spectra_reporting(
        &self,
        map: &ConsensusMap,
        min_peak_nr: usize,
        min_correlation: f64,
        max_lag: i32,
        max_rt_apex_difference: f64,
        reporter: &mut ProgressReporter<'_>,
    ) -> Result<Vec<PseudoSpectrum>> {
        let cache = self.cache_reporting(map, reporter)?;
        let count = progress_value(map.features.len())?;
        let mut remaining = self.max_work;
        reporter.section(0, count, CORRELATION_PROGRESS_LABEL, |reporter| {
            let mut used = vec![false; map.features.len()];
            let mut spectra = Vec::new();
            for i in 0..map.features.len() {
                reporter.set_count(i)?;
                if used[i] {
                    continue;
                }
                used[i] = true;
                let seed = &cache.feature_points[i];
                let (Some(&(first_point, _)), Some(&(last_point, _))) = (seed.first(), seed.last())
                else {
                    return Err(invalid("a seed consensus feature has no feature handles"));
                };
                let current_rt = cache.rt_cache[i];
                let mut spectrum = MSSpectrum {
                    rt: current_rt,
                    ms_level: 2,
                    ..MSSpectrum::default()
                };
                spectrum.peaks.push(Peak1D::new(
                    cache.max_intensities[i].0,
                    cache.max_intensities[i].1 as f32,
                ));
                let mut correlations = Vec::new();
                for j in i + 1..map.features.len() {
                    let rt = cache.rt_cache[j];
                    if (rt - current_rt).abs() > max_rt_apex_difference {
                        continue;
                    }
                    if rt < first_point || rt > last_point {
                        continue;
                    }
                    let other = &cache.feature_points[j];
                    let length = (seed.len() + other.len()) as u64;
                    let cost = length.saturating_add(length.saturating_mul(length));
                    remaining = remaining
                        .checked_sub(cost)
                        .ok_or_else(|| invalid("mass trace correlation work limit exceeded"))?;
                    let scores =
                        self.score_hullpoints(seed, other, min_correlation, max_lag, 0.1)?;
                    let Some((lag, lag_intensity)) = scores.lag else {
                        continue;
                    };
                    if scores.pearson_score > min_correlation && lag >= -max_lag && lag <= max_lag {
                        used[j] = true;
                        spectrum.peaks.push(Peak1D::new(
                            cache.max_intensities[j].0,
                            cache.max_intensities[j].1 as f32,
                        ));
                        correlations.push(TraceCorrelation {
                            rt_apex: map.features[j].base.rt as f32,
                            rt_diff: (map.features[i].base.rt - map.features[j].base.rt).abs()
                                as f32,
                            lag: lag as f32,
                            pearson_score: scores.pearson_score as f32,
                            lag_intensity: lag_intensity as f32,
                        });
                    }
                }
                if spectrum.peaks.len() > min_peak_nr {
                    spectra.push(PseudoSpectrum {
                        spectrum,
                        correlations,
                    });
                }
            }
            Ok(spectra)
        })
    }
}

/// Source `matchMassTraces_`: pair up two RT-sorted traces.
///
/// Points less than `mindiff` apart pair up; otherwise the earlier point pairs
/// with a zero. With `pad_ends`, the rest of the longer trace is paired with
/// zeros; without, the traces are first advanced to their first match and the
/// unmatched tails are dropped.
fn match_mass_traces(
    hull_points1: &[(f64, f64)],
    hull_points2: &[(f64, f64)],
    mindiff: f64,
    pad_ends: bool,
) -> Result<(Vec<f64>, Vec<f64>)> {
    if hull_points1
        .iter()
        .chain(hull_points2)
        .any(|p| !p.0.is_finite())
        || mindiff.is_nan()
    {
        return Err(invalid(
            "mass trace retention times and the matching distance must be finite",
        ));
    }
    let (mut k, mut m) = (0, 0);
    let (mut first, mut second) = (Vec::new(), Vec::new());
    if !pad_ends {
        while k < hull_points1.len() && m < hull_points2.len() {
            if (hull_points1[k].0 - hull_points2[m].0).abs() < mindiff {
                break;
            } else if hull_points1[k].0 > hull_points2[m].0 {
                m += 1;
            } else {
                k += 1;
            }
        }
    }
    while k < hull_points1.len() && m < hull_points2.len() {
        if (hull_points1[k].0 - hull_points2[m].0).abs() < mindiff {
            first.push(hull_points1[k].1);
            second.push(hull_points2[m].1);
            m += 1;
            k += 1;
        } else if hull_points1[k].0 > hull_points2[m].0 {
            first.push(0.0);
            second.push(hull_points2[m].1);
            m += 1;
        } else {
            first.push(hull_points1[k].1);
            second.push(0.0);
            k += 1;
        }
    }
    if pad_ends {
        for point in &hull_points1[k..] {
            first.push(point.1);
            second.push(0.0);
        }
        for point in &hull_points2[m..] {
            first.push(0.0);
            second.push(point.1);
        }
    }
    Ok((first, second))
}

#[cfg(test)]
mod tests {
    use super::match_mass_traces;

    #[test]
    fn matching_pads_ends_or_drops_unmatched_tails() {
        let a = [(1.0, 10.0), (2.0, 20.0), (3.0, 30.0)];
        let b = [(2.05, 5.0), (3.0, 6.0), (4.0, 7.0)];
        let (x, y) = match_mass_traces(&a, &b, 0.1, true).unwrap();
        assert_eq!(x, [10.0, 20.0, 30.0, 0.0]);
        assert_eq!(y, [0.0, 5.0, 6.0, 7.0]);
        let (x, y) = match_mass_traces(&a, &b, 0.1, false).unwrap();
        assert_eq!(x, [20.0, 30.0]);
        assert_eq!(y, [5.0, 6.0]);
        assert!(match_mass_traces(&[(f64::NAN, 1.0)], &b, 0.1, true).is_err());
    }
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `ANALYSIS/ID/IDDecoyProbability.h`: probabilities from a
//! target/decoy search, by modelling the score distributions.

use crate::identification::PeptideIdentification;
use crate::math::fitters::gamma::{GammaDistributionFitResult, GammaDistributionFitter};
use crate::math::fitters::gauss::{GaussFitResult, GaussFitter};
use crate::metadata::MetaValue;
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{Error, Result};

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// Largest `number_of_bins` accepted: the histograms are allocated up front.
pub const MAX_BINS: usize = 1_000_000;

/// A histogram normalization (source `Transformation_`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Transformation {
    max_intensity: f64,
    diff_score: f64,
    min_score: f64,
    max_score: f64,
    max_intensity_bin: usize,
}

/// Calculates probabilities with a target/decoy approach (source
/// `IDDecoyProbability`).
///
/// As in PeptideProphet, the distribution of forward (target) scores is
/// modelled with a Gaussian and that of reverse (decoy) scores with a Gamma
/// distribution; each hit's new score is `rho_fwd / (rho_fwd + rho_rev)` at
/// its score.
///
/// Scores where lower is better (e.g. E-values) are transformed to
/// `-log10(score)` first; a score below `10^-lower_score_better_default_value_if_zero`
/// becomes that parameter's value instead, since its logarithm would be
/// infinite. Every hit keeps its original score as the meta value
/// `<score type>_score`, and each identification's score type becomes
/// `<score type>_DecoyProbability`, higher is better. Identifications without
/// hits are dropped from the output.
///
/// # Notes
///
/// These points reproduce the source arithmetic:
///
/// - the final probability transform applies `-log10` without the
///   zero-score default, so a lower-is-better score of 0 is infinite there;
/// - the score range of a histogram starts from `f64::MIN_POSITIVE` (the
///   source's `numeric_limits<double>::min()`), not from the lowest double, so
///   when all scores are negative the range's maximum is that tiny positive
///   number;
/// - the forward histogram used for the Gaussian's data has
///   `number_of_bins - 1` as its scale, the reverse one `number_of_bins`.
///
/// These diverge from it:
///
/// - the reverse histogram puts a decoy score equal to the overall maximum
///   into bin `number_of_bins`, one past the end of the source's vector (an
///   out-of-bounds write); it goes to the last bin here;
/// - the Gaussian fit falls back to the initial guess on any fit error, the
///   source only on `Exception::UnableToFit` (its other failures are
///   undefined arithmetic).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IDDecoyProbability {
    /// `number_of_bins`: number of bins used for the fitting; if sparse
    /// datasets are used, this number should be smaller. Advanced.
    pub number_of_bins: usize,
    /// `lower_score_better_default_value_if_zero`: used when e.g. an E-value
    /// is 0 and cannot be transformed to a real number (the log of the
    /// E-value). Advanced.
    pub lower_score_better_default_value_if_zero: f64,
}

impl Default for IDDecoyProbability {
    fn default() -> Self {
        Self {
            number_of_bins: 40,
            lower_score_better_default_value_if_zero: 50.0,
        }
    }
}

impl IDDecoyProbability {
    /// The default parameters, in source declaration order. The source's
    /// debug-build parameters `rev_filename` and `fwd_filename` (which name
    /// gnuplot output) are not part of a release build and are not ported.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters in the layout of [`IDDecoyProbability::defaults`],
    /// carrying the current values (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `number_of_bins` does not fit the
    /// source 32-bit parameter.
    pub fn to_param(&self) -> Result<Param> {
        let bins = i32::try_from(self.number_of_bins)
            .map_err(|_| bad("'number_of_bins' does not fit the source 32-bit parameter"))?;
        let advanced = ["advanced".to_owned()];
        let mut param = Param::new();
        param.set_value(
            "number_of_bins",
            ParamValue::Integer(i64::from(bins)),
            "Number of bins used for the fitting, if sparse datasets are used, this number should be smaller",
            &advanced,
        )?;
        param.set_value(
            "lower_score_better_default_value_if_zero",
            ParamValue::Float(self.lower_score_better_default_value_if_zero),
            "This value is used if e.g. a E-value score is 0 and cannot be transformed in a real number (log of E-value)",
            &advanced,
        )?;
        Ok(param)
    }

    /// Build from parameters, as source `setParameters`. Missing parameters
    /// take their defaults; unknown ones are ignored (see
    /// [`IDDecoyProbability::from_param_with_warnings`]).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// `number_of_bins` is negative.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`IDDecoyProbability::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`IDDecoyProbability::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new("IDDecoyProbability")?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |param| {
            Ok(Self {
                number_of_bins: usize::try_from(param.value("number_of_bins")?.to_i64()?)
                    .map_err(|_| bad("'number_of_bins' must not be negative"))?,
                lower_score_better_default_value_if_zero: param
                    .value("lower_score_better_default_value_if_zero")?
                    .to_f64()?,
            })
        })
    }

    /// The score the histograms use: `-log10` for lower-is-better scores,
    /// with the zero default.
    fn transformed(&self, score: f64, higher_better: bool) -> f64 {
        if higher_better {
            return score;
        }
        let threshold = 10_f64.powf(-self.lower_score_better_default_value_if_zero);
        if score < threshold {
            self.lower_score_better_default_value_if_zero
        } else {
            -score.log10()
        }
    }

    /// Convert identifications of a combined target/decoy search into
    /// probabilities (source `apply(PeptideIdentificationList&)`). Decoy hits
    /// are those annotated `target_decoy` = `decoy`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a hit's target/decoy annotation is
    /// invalid, or as the modelling step: `number_of_bins` below 2 or above
    /// [`MAX_BINS`]; no forward or no reverse score, or all scores equal (the
    /// source divides by a zero range and indexes with the result); a Gamma
    /// fit failure (source `Exception::UnableToFit`); or a probability that is
    /// not finite (the source writes NaN or infinite scores). `ids` is
    /// unchanged on error.
    pub fn apply(&self, ids: &mut Vec<PeptideIdentification>) -> Result<()> {
        let mut next = ids.clone();
        let (mut rev, mut fwd, mut all) = (Vec::new(), Vec::new(), Vec::new());
        for pep in &mut next {
            let key = format!("{}_Score", pep.score_type);
            let higher_better = pep.higher_score_better;
            for hit in &mut pep.hits {
                hit.metadata
                    .insert(key.clone(), MetaValue::try_from(hit.score)?);
                let score = self.transformed(hit.score, higher_better);
                if hit.is_decoy()? {
                    rev.push(score);
                } else {
                    fwd.push(score);
                }
                all.push(score);
            }
        }
        self.model(&mut next, &rev, &fwd, &all)?;
        *ids = next;
        Ok(())
    }

    /// Convert the identifications of a forward search into probabilities,
    /// using a separate reversed (decoy) search (source
    /// `apply(prob_ids, fwd_ids, rev_ids)`); returns the converted forward
    /// identifications.
    ///
    /// # Errors
    ///
    /// As [`IDDecoyProbability::apply`].
    pub fn apply_forward_reverse(
        &self,
        fwd_ids: &[PeptideIdentification],
        rev_ids: &[PeptideIdentification],
    ) -> Result<Vec<PeptideIdentification>> {
        let mut next = fwd_ids.to_vec();
        let (mut rev, mut fwd, mut all) = (Vec::new(), Vec::new(), Vec::new());
        for pep in &mut next {
            let key = format!("{}_Score", pep.score_type);
            let higher_better = pep.higher_score_better;
            for hit in &mut pep.hits {
                hit.metadata
                    .insert(key.clone(), MetaValue::try_from(hit.score)?);
                let score = self.transformed(hit.score, higher_better);
                fwd.push(score);
                all.push(score);
            }
        }
        for pep in rev_ids {
            for hit in &pep.hits {
                let score = self.transformed(hit.score, pep.higher_score_better);
                rev.push(score);
                all.push(score);
            }
        }
        self.model(&mut next, &rev, &fwd, &all)?;
        Ok(next)
    }

    /// Source `normalizeBins_`: histogram `scores` into `bins` bins scaled so
    /// that the fullest bin is 4 ("the best value for the gamma
    /// distribution").
    fn normalize_bins(&self, scores: &[f64]) -> Result<(Vec<f64>, Transformation)> {
        let bins = self.number_of_bins;
        let mut max = f64::MIN_POSITIVE;
        let mut min = f64::MAX;
        for &score in scores {
            if score > max {
                max = score;
            }
            if score < min {
                min = score;
            }
        }
        let diff = max - min;
        if scores.is_empty() || !diff.is_finite() || diff <= 0.0 {
            return Err(bad(
                "IDDecoyProbability needs scores spanning a finite, non-empty range",
            ));
        }
        let mut binned = vec![0.0; bins];
        let mut max_bin = 0.0;
        let mut max_bin_number = 0;
        for &score in scores {
            let bin = bin_index((score - min) / diff * (bins - 1) as f64, bins);
            binned[bin] += 1.0;
            if binned[bin] > max_bin {
                max_bin = binned[bin];
                max_bin_number = bin;
            }
        }
        for value in &mut binned {
            *value /= max_bin / 4.0;
        }
        Ok((
            binned,
            Transformation {
                max_intensity: 4.0 / max_bin,
                diff_score: diff,
                min_score: min,
                max_score: max,
                max_intensity_bin: max_bin_number,
            },
        ))
    }

    /// Source `apply_`: fit the two models and rescore every hit.
    fn model(
        &self,
        ids: &mut Vec<PeptideIdentification>,
        rev_scores: &[f64],
        fwd_scores: &[f64],
        all_scores: &[f64],
    ) -> Result<()> {
        let bins = self.number_of_bins;
        if !(2..=MAX_BINS).contains(&bins) {
            return Err(bad(format!(
                "'number_of_bins' must be between 2 and {MAX_BINS}"
            )));
        }
        let (rev_normalized, rev_trafo) = self.normalize_bins(rev_scores)?;
        let (_, fwd_trafo) = self.normalize_bins(fwd_scores)?;
        let (_, all_trafo) = self.normalize_bins(all_scores)?;

        let rev_data: Vec<(f64, f64)> = (0..bins)
            .map(|i| (i as f64 / bins as f64 + 0.0001, rev_normalized[i]))
            .collect();
        let mut gamma_fitter = GammaDistributionFitter::new();
        gamma_fitter.set_initial_parameters(GammaDistributionFitResult::new(1.0, 3.0));
        let gamma = gamma_fitter.fit(&rev_data)?;

        let (min, diff) = (all_trafo.min_score, all_trafo.diff_score);
        let mut fwd_bins = vec![0_usize; bins];
        let mut rev_bins = vec![0_usize; bins];
        let mut max_bin = 0;
        for &score in fwd_scores {
            let bin = bin_index((score - min) / diff * (bins - 1) as f64, bins);
            fwd_bins[bin] += 1;
            max_bin = max_bin.max(fwd_bins[bin]);
        }
        let (mut max_reverse_bin, mut max_reverse_bin_value) = (0, 0);
        for &score in rev_scores {
            let bin = bin_index((score - min) / diff * bins as f64, bins);
            rev_bins[bin] += 1;
            max_bin = max_bin.max(rev_bins[bin]);
            if rev_bins[bin] > max_reverse_bin_value {
                max_reverse_bin = bin;
                max_reverse_bin_value = rev_bins[bin];
            }
        }
        let diff_scores: Vec<f64> = (0..bins)
            .map(|i| {
                let (fwd, rev) = (fwd_bins[i], rev_bins[i]);
                if fwd as f64 > 1.3 * rev as f64 && max_reverse_bin < i {
                    (fwd - rev) as f64 / max_bin as f64
                } else {
                    0.0
                }
            })
            .collect();

        let mut diff_data = Vec::with_capacity(bins);
        let (mut gauss_a, mut gauss_x0, mut norm_factor) = (0.0_f64, 0.0, 0.0);
        for (i, &y) in diff_scores.iter().enumerate() {
            let x = i as f64 / bins as f64;
            if y > gauss_a {
                gauss_a = y;
            }
            gauss_x0 += x * y;
            norm_factor += y;
            diff_data.push((x, y));
        }
        gauss_x0 /= diff_data.len() as f64;
        gauss_x0 /= norm_factor;
        let mut gauss_sigma = 0.0;
        for i in 0..=bins {
            gauss_sigma += (gauss_x0 - i as f64 / bins as f64).abs();
        }
        gauss_sigma /= diff_data.len() as f64;
        let initial = GaussFitResult::new(gauss_a, gauss_x0, gauss_sigma);
        let mut gauss_fitter = GaussFitter::new();
        gauss_fitter.set_initial_parameters(initial);
        let gauss = gauss_fitter.fit(&diff_data).unwrap_or(initial);

        let mut output = Vec::with_capacity(ids.len());
        for pep in ids.iter() {
            if pep.hits.is_empty() {
                continue;
            }
            let key = format!("{}_score", pep.score_type);
            let mut next = pep.clone();
            for hit in &mut next.hits {
                let score = if pep.higher_score_better {
                    hit.score
                } else {
                    -hit.score.log10()
                };
                hit.metadata
                    .insert(key.clone(), MetaValue::try_from(hit.score)?);
                let probability = self.probability(&gamma, &rev_trafo, &gauss, &fwd_trafo, score);
                if !probability.is_finite() {
                    return Err(bad(format!(
                        "decoy probability of score {} is not finite",
                        hit.score
                    )));
                }
                hit.score = probability;
            }
            next.higher_score_better = true;
            next.score_type = format!("{}_DecoyProbability", pep.score_type);
            output.push(next);
        }
        *ids = output;
        Ok(())
    }

    /// Source `getProbability_`.
    fn probability(
        &self,
        gamma: &GammaDistributionFitResult,
        gamma_trafo: &Transformation,
        gauss: &GaussFitResult,
        gauss_trafo: &Transformation,
        score: f64,
    ) -> f64 {
        let bins = self.number_of_bins as f64;
        let rev_trans = (score - gamma_trafo.min_score) / gamma_trafo.diff_score;
        let rho_rev = if rev_trans < gamma_trafo.max_intensity_bin as f64 / bins {
            1.0 / gamma_trafo.max_intensity
        } else {
            gamma.b.powf(gamma.p) / libm::tgamma(gamma.p)
                * rev_trans.powf(gamma.p - 1.0)
                * (-gamma.b * rev_trans).exp()
        };
        let fwd_trans = (score - gauss_trafo.min_score) / gauss_trafo.diff_score;
        let rho_fwd = if fwd_trans < gauss.x0 {
            gauss.a * (-(fwd_trans - gauss.x0).powi(2) / 2.0 / gauss.sigma.powi(2)).exp()
        } else {
            1.0
        };
        rho_fwd / (rho_fwd + rho_rev)
    }
}

/// The source's `(Size)` cast of a non-negative bin position, kept inside the
/// histogram (see the struct notes).
fn bin_index(position: f64, bins: usize) -> usize {
    (position as usize).min(bins - 1)
}

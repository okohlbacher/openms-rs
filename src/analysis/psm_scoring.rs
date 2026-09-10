// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Native Morpheus and HyperScore peptide-spectrum matching scores.

use crate::comparison::{Tolerance, finite_score, matched_alignment, validate_spectrum};
use crate::kernel::nearest;
use crate::{Error, MSSpectrum, Result};

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn tolerance_value(tolerance: Tolerance) -> Result<f64> {
    let (Tolerance::Absolute(value) | Tolerance::Ppm(value)) = tolerance;
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(bad("fragment tolerance must be finite and nonnegative"))
    }
}
fn window(tolerance: Tolerance, mz: f64) -> Result<f64> {
    finite_score(match tolerance {
        Tolerance::Absolute(v) => v,
        Tolerance::Ppm(v) => mz * v * 1e-6,
    })
}
fn f32_result(value: f64) -> Result<f32> {
    let result = value as f32;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(bad("PSM result exceeds finite f32 range"))
    }
}
fn validate_pair(exp: &MSSpectrum, theo: &MSSpectrum, tolerance: Tolerance) -> Result<()> {
    tolerance_value(tolerance)?;
    validate_spectrum(exp, true)?;
    validate_spectrum(theo, true)
}
fn validate_charges(
    exp: &MSSpectrum,
    theo: &MSSpectrum,
    charges: Option<(&[i32], &[i32])>,
) -> Result<()> {
    if charges.is_some_and(|(e, t)| e.len() != exp.len() || t.len() != theo.len()) {
        return Err(bad("charge array length differs from spectrum length"));
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MorpheusResult {
    pub matches: usize,
    pub n_peaks: usize,
    pub score: f32,
    pub matched_intensity: f32,
    pub total_intensity: f32,
    pub mean_error_da: f32,
    pub mean_error_ppm: f32,
}
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MorpheusScore {
    pub tolerance: Tolerance,
}
impl MorpheusScore {
    pub fn compute(
        &self,
        experimental: &MSSpectrum,
        theoretical: &MSSpectrum,
    ) -> Result<MorpheusResult> {
        self.compute_impl(experimental, theoretical, None)
    }
    pub fn compute_with_charges(
        &self,
        experimental: &MSSpectrum,
        experimental_charges: &[i32],
        theoretical: &MSSpectrum,
        theoretical_charges: &[i32],
    ) -> Result<MorpheusResult> {
        self.compute_impl(
            experimental,
            theoretical,
            Some((experimental_charges, theoretical_charges)),
        )
    }
    fn compute_impl(
        &self,
        exp: &MSSpectrum,
        theo: &MSSpectrum,
        charges: Option<(&[i32], &[i32])>,
    ) -> Result<MorpheusResult> {
        validate_pair(exp, theo, self.tolerance)?;
        validate_charges(exp, theo, charges)?;
        if exp.is_empty() || theo.is_empty() {
            return Ok(MorpheusResult::default());
        }
        let charge_matches = |e: usize, t: usize| charges.is_none_or(|(ez, tz)| ez[e] == tz[t]);
        let total = finite_score(exp.peaks.iter().map(|p| f64::from(p.intensity)).sum())?;
        if total == 0.0 {
            return Err(bad(
                "Morpheus score requires positive experimental total intensity",
            ));
        }
        let (mut t, mut e, mut matches) = (0, 0, 0);
        // The first pass counts theoretical peaks, allowing experimental reuse.
        while t < theo.len() && e < exp.len() {
            let delta = exp.peaks[e].mz - theo.peaks[t].mz;
            if delta.abs() <= window(self.tolerance, theo.peaks[t].mz)? {
                if charge_matches(e, t) {
                    matches += 1;
                }
                t += 1;
            } else if delta < 0.0 {
                e += 1;
            } else {
                t += 1;
            }
        }
        let (mut t, mut e) = (0, 0);
        let (mut matched, mut error, mut ppm) = (0.0, 0.0, 0.0);
        // The second pass counts experimental ion current, allowing theoretical
        // reuse. Charge mismatches still advance the source's chosen pointer.
        while t < theo.len() && e < exp.len() {
            let delta = exp.peaks[e].mz - theo.peaks[t].mz;
            if delta.abs() <= window(self.tolerance, theo.peaks[t].mz)? {
                if charge_matches(e, t) {
                    matched += f64::from(exp.peaks[e].intensity);
                    error += delta.abs();
                    ppm += finite_score(delta.abs() / theo.peaks[t].mz * 1e6)?;
                }
                e += 1;
            } else if delta < 0.0 {
                e += 1;
            } else {
                t += 1;
            }
        }
        Ok(MorpheusResult {
            matches,
            n_peaks: theo.len(),
            score: f32_result(matches as f64 + matched / total)?,
            matched_intensity: f32_result(matched)?,
            total_intensity: f32_result(total)?,
            mean_error_da: f32_result(if matches == 0 {
                1e10
            } else {
                error / matches as f64
            })?,
            mean_error_ppm: f32_result(if matches == 0 {
                1e10
            } else {
                ppm / matches as f64
            })?,
        })
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HyperScoreDetail {
    pub score: f64,
    pub matched_prefix_ions: usize,
    pub matched_suffix_ions: usize,
    /// Da for absolute tolerance, ppm for a ppm tolerance.
    pub mean_error: f64,
}
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HyperScore {
    pub tolerance: Tolerance,
}
fn ion_names(theoretical: &MSSpectrum) -> Result<&[String]> {
    let first = theoretical.string_data_arrays.first().ok_or_else(|| {
        bad("HyperScore requires theoretical ion annotations in the first string array")
    })?;
    if first.data.len() != theoretical.len() || first.data.iter().any(String::is_empty) {
        return Err(bad(
            "theoretical ion annotations must be aligned and nonempty",
        ));
    }
    Ok(&first.data)
}
fn by_series(name: &str) -> Option<bool> {
    if name.starts_with('y') || name.contains("$y") {
        Some(false)
    } else if name.starts_with('b') || name.contains("$b") {
        Some(true)
    } else {
        None
    }
}
fn terminal_series(name: &str) -> Option<bool> {
    let category = |c| match c {
        b'a' | b'b' | b'c' => Some(true),
        b'x' | b'y' | b'z' => Some(false),
        _ => None,
    };
    category(name.as_bytes()[0]).or_else(|| {
        name.find('$')
            .and_then(|i| name.as_bytes().get(i + 1))
            .and_then(|&c| category(c))
    })
}
fn log_factorial(n: usize) -> f64 {
    (2..=n).map(|i| (i as f64).ln()).sum()
}
fn hyperscore(dot: f64, prefix: usize, suffix: usize) -> Result<f64> {
    finite_score(dot.ln_1p() + log_factorial(prefix) + log_factorial(suffix))
}
impl HyperScore {
    /// Source basic HyperScore counts b/y ions, including cross-link names with
    /// $b/$y. Every matched peak contributes to the dot product.
    pub fn compute(&self, experimental: &MSSpectrum, theoretical: &MSSpectrum) -> Result<f64> {
        Ok(self
            .compute_detail_impl(experimental, theoretical, false)?
            .score)
    }
    /// Source detailed overload counts a/b/c and x/y/z rather than only b/y.
    pub fn compute_with_detail(
        &self,
        experimental: &MSSpectrum,
        theoretical: &MSSpectrum,
    ) -> Result<HyperScoreDetail> {
        self.compute_detail_impl(experimental, theoretical, true)
    }
    fn compute_detail_impl(
        &self,
        exp: &MSSpectrum,
        theo: &MSSpectrum,
        detail: bool,
    ) -> Result<HyperScoreDetail> {
        validate_pair(exp, theo, self.tolerance)?;
        if exp.is_empty() || theo.is_empty() {
            return Ok(HyperScoreDetail::default());
        }
        let names = ion_names(theo)?;
        let (mut dot, mut error, mut prefix, mut suffix) = (0.0, 0.0, 0, 0);
        for (t, e) in matched_alignment(theo, exp, self.tolerance)? {
            // The source multiplies its f32 intensities before adding to f64.
            dot += f64::from(theo.peaks[t].intensity * exp.peaks[e].intensity);
            if detail {
                let delta = (exp.peaks[e].mz - theo.peaks[t].mz).abs();
                error += match self.tolerance {
                    Tolerance::Absolute(_) => delta,
                    Tolerance::Ppm(_) => finite_score(delta / theo.peaks[t].mz * 1e6)?,
                };
            }
            match if detail {
                terminal_series(&names[t])
            } else {
                by_series(&names[t])
            } {
                Some(true) => prefix += 1,
                Some(false) => suffix += 1,
                None => (),
            }
        }
        finite_score(dot)?;
        Ok(HyperScoreDetail {
            score: hyperscore(dot, prefix, suffix)?,
            matched_prefix_ions: prefix,
            matched_suffix_ions: suffix,
            mean_error: if prefix + suffix == 0 {
                0.0
            } else {
                finite_score(error / (prefix + suffix) as f64)?
            },
        })
    }
    /// Source charge-aware overload uses exact f64 nearest peaks and a strict
    /// tolerance boundary. Only charge-matched b/y ions contribute intensity.
    pub fn compute_with_charges(
        &self,
        experimental: &MSSpectrum,
        experimental_charges: &[i32],
        theoretical: &MSSpectrum,
        theoretical_charges: &[i32],
    ) -> Result<f64> {
        let (score, _) = self.charge_score(
            experimental,
            experimental_charges,
            theoretical,
            theoretical_charges,
            None,
        )?;
        Ok(score)
    }
    /// Also accumulate source per-cleavage b/y experimental intensities into an
    /// existing peptide-length array. Invalid ion indices leave the array intact.
    pub fn compute_with_intensity_sum(
        &self,
        experimental: &MSSpectrum,
        experimental_charges: &[i32],
        theoretical: &MSSpectrum,
        theoretical_charges: &[i32],
        intensity_sum: &mut [f64],
    ) -> Result<f64> {
        if intensity_sum.is_empty()
            || intensity_sum.len() > 100_000
            || intensity_sum.iter().any(|v| !v.is_finite())
        {
            return Err(bad(
                "peptide intensity array requires 1..=100000 finite entries",
            ));
        }
        let (score, additions) = self.charge_score(
            experimental,
            experimental_charges,
            theoretical,
            theoretical_charges,
            Some(intensity_sum.len()),
        )?;
        let next = intensity_sum
            .iter()
            .zip(additions)
            .map(|(&a, b)| finite_score(a + b))
            .collect::<Result<Vec<_>>>()?;
        intensity_sum.copy_from_slice(&next);
        Ok(score)
    }
    fn charge_score(
        &self,
        exp: &MSSpectrum,
        ez: &[i32],
        theo: &MSSpectrum,
        tz: &[i32],
        length: Option<usize>,
    ) -> Result<(f64, Vec<f64>)> {
        validate_pair(exp, theo, self.tolerance)?;
        validate_charges(exp, theo, Some((ez, tz)))?;
        let n = length.unwrap_or(0);
        let mut b = vec![0.0; n];
        let mut y = vec![0.0; n];
        if exp.is_empty() || theo.is_empty() {
            return Ok((0.0, b));
        }
        let names = ion_names(theo)?;
        let (mut dot, mut prefix, mut suffix) = (0.0, 0, 0);
        for (t, peak) in theo.peaks.iter().enumerate() {
            let e = nearest(&exp.peaks, peak.mz, |p| p.mz).expect("nonempty validated spectrum");
            if (exp.peaks[e].mz - peak.mz).abs() >= window(self.tolerance, peak.mz)?
                || ez[e] != tz[t]
            {
                continue;
            }
            let series = if length.is_some() {
                match names[t].as_bytes()[0] {
                    b'b' => Some(true),
                    b'y' => Some(false),
                    _ => None,
                }
            } else {
                by_series(&names[t])
            };
            let Some(is_b) = series else {
                continue;
            };
            if let Some(n) = length {
                let digits: String = names[t][1..]
                    .chars()
                    .take_while(char::is_ascii_digit)
                    .collect();
                if digits.is_empty() {
                    continue;
                }
                let index: usize = digits.parse().map_err(|_| bad("ion ordinal overflows"))?;
                if index == 0 || index > n {
                    return Err(bad("ion ordinal falls outside peptide intensity array"));
                }
                let cell = if is_b {
                    &mut b[index - 1]
                } else {
                    &mut y[n - index]
                };
                *cell = finite_score(*cell + f64::from(exp.peaks[e].intensity))?;
            } else if is_b {
                prefix += 1;
            } else {
                suffix += 1;
            }
            // Here theo_intensity is explicitly f64 in the source overload.
            dot += f64::from(peak.intensity) * f64::from(exp.peaks[e].intensity);
        }
        if length.is_some() {
            prefix = b.iter().filter(|&&v| v > 0.0).count();
            suffix = y.iter().filter(|&&v| v > 0.0).count();
        }
        finite_score(dot)?;
        let additions = b
            .into_iter()
            .zip(y)
            .map(|(b, y)| finite_score(b + y))
            .collect::<Result<Vec<_>>>()?;
        Ok((
            if prefix + suffix == 0 {
                0.0
            } else {
                hyperscore(dot, prefix, suffix)?
            },
            additions,
        ))
    }
}

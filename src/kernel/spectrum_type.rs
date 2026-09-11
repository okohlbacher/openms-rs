// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Spectrum metadata precedence and the shared source shoulder estimator.
use super::{MSSpectrum, SpectrumType};
use crate::{Error, Result, metadata::ProcessingAction};

/// Limits for actual spectrum-type estimation and cumulative history work.
/// Explicit stored types and unrequested peak data are not traversed.
#[derive(Clone, Copy, Debug)]
pub struct SpectrumTypeQueryLimits {
    pub max_points: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for SpectrumTypeQueryLimits {
    fn default() -> Self {
        Self {
            max_points: 1_000_000,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}
fn resource() -> Error {
    Error::InvalidValue("spectrum type query resource limit exceeded".into())
}
fn spend(work: &mut usize, amount: usize) -> Result<()> {
    *work = work.checked_sub(amount).ok_or_else(resource)?;
    Ok(())
}
impl MSSpectrum {
    /// Stored type wins, then PeakPicking history, then optional peak estimation.
    /// This does not cache the answer, sort peaks or validate unrelated graphs.
    pub fn get_type(&self, query_data: bool) -> Result<SpectrumType> {
        self.get_type_with_limits(query_data, SpectrumTypeQueryLimits::default())
    }
    /// As [`Self::get_type`], with explicit resource ceilings.
    pub fn get_type_with_limits(
        &self,
        query_data: bool,
        mut limits: SpectrumTypeQueryLimits,
    ) -> Result<SpectrumType> {
        self.get_type_with_budget(
            query_data,
            limits.max_points,
            &mut limits.max_work,
            &mut limits.max_bytes,
        )
    }
    pub(crate) fn get_type_with_budget(
        &self,
        query_data: bool,
        max_points: usize,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<SpectrumType> {
        if self.spectrum_type != SpectrumType::Unknown {
            return Ok(self.spectrum_type);
        }
        for processing in &self.data_processing {
            let height = (usize::BITS - processing.actions.len().leading_zeros()) as usize;
            spend(work, 1 + 12 * height)?;
            if processing.actions.contains(&ProcessingAction::PeakPicking) {
                return Ok(SpectrumType::Centroid);
            }
        }
        if !query_data || self.peaks.len() < 5 {
            return Ok(SpectrumType::Unknown);
        }
        let n = self.peaks.len();
        if n > max_points {
            return Err(resource());
        }
        // Two scalar copies/finite checks, total sum and five complete maximum,
        // left/right shoulder scans. Conservative charge precedes allocation.
        spend(work, n.checked_mul(32).ok_or_else(resource)?)?;
        let storage = n
            .checked_mul(2 * std::mem::size_of::<f64>())
            .ok_or_else(resource)?;
        *bytes = bytes.checked_sub(storage).ok_or_else(resource)?;
        if self
            .peaks
            .iter()
            .any(|p| !p.mz.is_finite() || !p.intensity.is_finite())
        {
            return Err(Error::InvalidValue(
                "spectrum type estimation requires finite consumed peak values".into(),
            ));
        }
        let mut x = Vec::new();
        let mut y = Vec::new();
        x.try_reserve_exact(n).map_err(|_| resource())?;
        y.try_reserve_exact(n).map_err(|_| resource())?;
        for peak in &self.peaks {
            x.push(peak.mz);
            y.push(f64::from(peak.intensity));
        }
        Ok(estimate(&x, &mut y))
    }
}

// Both callers establish equal-length finite scalar arrays. The public picker
// retains its stricter nonnegative, sorted/distinct and full-record validation.
pub(crate) fn estimate(x: &[f64], y: &mut [f64]) -> SpectrumType {
    if y.len() < 5 {
        return SpectrumType::Unknown;
    }
    let total = y.iter().sum::<f64>();
    let mut explained = 0.0;
    let (mut profile, mut centroid) = (0, 0);
    for _ in 0..5 {
        if explained > 0.5 * total {
            break;
        }
        let mut maximum = 0.0;
        let mut index = None;
        for (i, &value) in y.iter().enumerate() {
            if value > maximum {
                maximum = value;
                index = Some(i);
            }
        }
        let Some(index) = index else {
            break;
        };
        let mut cursor = index;
        let mut last = maximum;
        while cursor > 0
            && y[cursor] <= last
            && y[cursor] > 0.0
            && y[cursor] / last > 0.1
            && x[cursor] + 1.0 > x[index]
        {
            last = y[cursor];
            explained += last;
            y[cursor] = 0.0;
            cursor -= 1;
        }
        if y[cursor] > last && cursor + 1 < y.len() {
            y[cursor + 1] = last;
        }
        let break_left = index - cursor < 3;
        y[index] = maximum;
        explained -= maximum;
        cursor = index;
        last = maximum;
        while cursor < y.len()
            && y[cursor] <= last
            && y[cursor] > 0.0
            && y[cursor] / last > 0.1
            && x[cursor] - 1.0 < x[index]
        {
            last = y[cursor];
            explained += last;
            y[cursor] = 0.0;
            cursor += 1;
        }
        if cursor < y.len() && y[cursor] > last && cursor > 0 {
            y[cursor - 1] = last;
        }
        if break_left || cursor - index < 3 {
            centroid += 1;
        } else {
            profile += 1;
        }
    }
    if profile as f32 / (profile + centroid) as f32 > 0.75 {
        SpectrumType::Profile
    } else {
        SpectrumType::Centroid
    }
}

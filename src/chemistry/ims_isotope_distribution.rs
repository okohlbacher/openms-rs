// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Nominal-bin IMS isotopes from OpenMS 54a232f, distinct from fine isotopes.
//! Explicit operation settings replace the source's two mutable global statics.

use crate::{Error, Result, param::ParamValue};

pub(super) const MAX_PEAKS: usize = 1_000_000;
const MAX_WORK: usize = 50_000_000;
const MAX_BYTES: usize = 64 * 1024 * 1024;
pub(super) const MAX_TEXT: usize = 8 * 1024 * 1024;

/// A stored mass defect and abundance. Neither quantity is implicitly normalized.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IMSIsotopePeak {
    pub mass: f64,
    pub abundance: f64,
}

/// Per-operation replacements for source `SIZE` and `ABUNDANCES_SUM_ERROR`.
/// Both source statics start at zero, as do these defaults.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IMSIsotopeOptions {
    pub size: usize,
    pub abundances_sum_error: f64,
}

impl IMSIsotopeOptions {
    fn checked_size(self) -> Result<usize> {
        if self.size > MAX_PEAKS {
            return Err(invalid("isotope size exceeds 1,000,000 peaks"));
        }
        Ok(self.size)
    }

    fn checked_error(self) -> Result<f64> {
        finite(self.abundances_sum_error, "abundance sum error")
    }
}

/// Owned nominal bins, retaining hidden trailing peaks in equality and averages.
/// Stored values are finite, but signed masses and abundances are supported.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IMSIsotopeDistribution {
    peaks: Vec<IMSIsotopePeak>,
    nominal_mass: u32,
}

impl IMSIsotopeDistribution {
    /// An empty distribution with the supplied nominal mass.
    pub fn new(nominal_mass: u32) -> Self {
        Self {
            peaks: Vec::new(),
            nominal_mass,
        }
    }

    /// A single unit-abundance isotope, with nominal mass zero.
    pub fn from_mass(mass: f64) -> Result<Self> {
        Self::from_peaks(
            vec![IMSIsotopePeak {
                mass,
                abundance: 1.0,
            }],
            0,
        )
    }

    /// Retains every input peak, without truncation or normalization.
    pub fn from_peaks(peaks: Vec<IMSIsotopePeak>, nominal_mass: u32) -> Result<Self> {
        check_count(peaks.len())?;
        for peak in &peaks {
            validate_peak(*peak)?;
        }
        Ok(Self {
            peaks,
            nominal_mass,
        })
    }

    pub fn peaks(&self) -> &[IMSIsotopePeak] {
        &self.peaks
    }
    pub fn stored_len(&self) -> usize {
        self.peaks.len()
    }
    pub fn is_empty(&self) -> bool {
        self.peaks.is_empty()
    }
    pub fn nominal_mass(&self) -> u32 {
        self.nominal_mass
    }
    pub fn set_nominal_mass(&mut self, nominal_mass: u32) {
        self.nominal_mass = nominal_mass;
    }

    /// Accessible length; unlike `is_empty`, this applies the requested size.
    pub fn size(&self, options: IMSIsotopeOptions) -> Result<usize> {
        Ok(self.peaks.len().min(options.checked_size()?))
    }

    /// Absolute mass `(stored defect + nominal mass) + index`, ignoring size.
    pub fn mass(&self, index: usize) -> Result<f64> {
        let peak = self
            .peaks
            .get(index)
            .ok_or_else(|| invalid("isotope index is out of range"))?;
        finite(
            (peak.mass + f64::from(self.nominal_mass)) + index as f64,
            "isotope mass",
        )
    }

    /// Stored abundance, ignoring the truncation setting.
    pub fn abundance(&self, index: usize) -> Result<f64> {
        self.peaks
            .get(index)
            .map(|peak| peak.abundance)
            .ok_or_else(|| invalid("isotope index is out of range"))
    }

    /// Source ordered weighted sum over all stored bins, without normalization.
    pub fn average_mass(&self) -> Result<f64> {
        let mut average = 0.0;
        for (index, peak) in self.peaks.iter().enumerate() {
            // Do not change the source addition/multiplication order.
            average += ((peak.mass + f64::from(self.nominal_mass)) + index as f64) * peak.abundance;
        }
        finite(average, "average isotope mass")
    }

    pub fn masses(&self, options: IMSIsotopeOptions) -> Result<Vec<f64>> {
        let len = self.size(options)?;
        let mut values = Work::default().vector(len)?;
        for index in 0..len {
            values.push(self.mass(index)?);
        }
        Ok(values)
    }

    pub fn abundances(&self, options: IMSIsotopeOptions) -> Result<Vec<f64>> {
        let len = self.size(options)?;
        let mut values = Work::default().vector(len)?;
        values.extend(self.peaks[..len].iter().map(|peak| peak.abundance));
        Ok(values)
    }

    /// Normalizes all stored bins, retaining the source's signed-sum behavior.
    /// A finite-input sum overflowing to positive infinity scales to signed zero.
    pub fn normalize(&mut self, options: IMSIsotopeOptions) -> Result<()> {
        let error = options.checked_error()?;
        let mut work = Work::default();
        work.consume(self.peaks.len().saturating_mul(3))?;
        let scale = normalization_scale(&self.peaks, error);
        if let Some(scale) = scale {
            // Validate every result before the first mutation; no clone is needed.
            for peak in &self.peaks {
                finite(peak.abundance * scale, "normalized abundance")?;
            }
            for peak in &mut self.peaks {
                peak.abundance *= scale;
            }
        }
        Ok(())
    }

    /// Ordered source convolution, with an immutable right-hand input.
    pub fn convolve(&self, rhs: &Self, options: IMSIsotopeOptions) -> Result<Self> {
        self.fold(rhs, options, &mut Work::default())
    }

    /// Atomically replaces this distribution; the borrowed right side is preserved.
    pub fn convolve_assign(&mut self, rhs: &Self, options: IMSIsotopeOptions) -> Result<()> {
        *self = self.convolve(rhs, options)?;
        Ok(())
    }

    /// Explicitly reproduces source right-hand padding on a successful nonempty fold.
    /// Both inputs are unchanged on error; an existing right-hand tail is retained.
    pub fn convolve_assign_with_padded_rhs(
        &mut self,
        rhs: &mut Self,
        options: IMSIsotopeOptions,
    ) -> Result<()> {
        let mut work = Work::default();
        let result = self.fold(rhs, options, &mut work)?;
        if !self.is_empty() && !rhs.is_empty() && rhs.peaks.len() < options.size {
            let mut padded = work.vector(options.size)?;
            padded.extend_from_slice(&rhs.peaks);
            padded.resize(options.size, IMSIsotopePeak::default());
            rhs.peaks = padded;
        }
        *self = result;
        Ok(())
    }

    /// Source binary folding order. In particular, powers zero and one are no-ops.
    pub fn pow(&self, power: u32, options: IMSIsotopeOptions) -> Result<Self> {
        let mut work = Work::default();
        let mut current = self.copy_with_work(&mut work)?;
        if power <= 1 {
            return Ok(current);
        }
        let mut result = Self::default();
        if power & 1 != 0 {
            result = current.copy_with_work(&mut work)?;
        }
        let mut remaining = power >> 1;
        while remaining != 0 {
            current = current.fold(&current, options, &mut work)?;
            // Once empty, every subsequent square and collection is the source's
            // empty-right no-op. Avoid materializing redundant result copies.
            if current.is_empty() {
                break;
            }
            if remaining & 1 != 0 {
                result = result.fold(&current, options, &mut work)?;
            }
            remaining >>= 1;
        }
        Ok(result)
    }

    pub fn pow_assign(&mut self, power: u32, options: IMSIsotopeOptions) -> Result<()> {
        if power <= 1 {
            return Ok(());
        }
        *self = self.pow(power, options)?;
        Ok(())
    }

    /// Classic-locale, six-significant-digit source stream representation.
    pub fn to_text(&self, options: IMSIsotopeOptions) -> Result<String> {
        self.text_with_work(options, &mut Work::default())
    }

    pub(super) fn text_with_work(
        &self,
        options: IMSIsotopeOptions,
        work: &mut Work,
    ) -> Result<String> {
        let len = self.size(options)?;
        // Each finite double needs at most 14 bytes in the shared six-digit formatter.
        let bytes = len
            .checked_mul(32)
            .ok_or_else(|| invalid("isotope text length overflow"))?;
        if bytes > MAX_TEXT {
            return Err(invalid("isotope text exceeds 8 MiB bound"));
        }
        // Two shared scalar renders each account for 256-byte scratch, 64-byte
        // output, and value measurement. Precharge 1 KiB per bin for both,
        // in addition to our destination, before invoking their local budgets.
        let formatting = len
            .checked_mul(1024)
            .ok_or_else(|| invalid("isotope formatting work overflow"))?;
        work.consume(formatting)?;
        work.allocate(formatting)?;
        work.allocate(bytes)?;
        let mut text = String::new();
        text.try_reserve_exact(bytes)
            .map_err(|_| invalid("cannot allocate isotope text"))?;
        for index in 0..len {
            text.push_str(&ParamValue::Float(self.mass(index)?).to_stream_text()?);
            text.push(' ');
            text.push_str(&ParamValue::Float(self.peaks[index].abundance).to_stream_text()?);
            text.push('\n');
        }
        Ok(text)
    }

    fn copy_with_work(&self, work: &mut Work) -> Result<Self> {
        let mut peaks = work.vector(self.peaks.len())?;
        peaks.extend_from_slice(&self.peaks);
        Ok(Self {
            peaks,
            nominal_mass: self.nominal_mass,
        })
    }

    fn fold(&self, rhs: &Self, options: IMSIsotopeOptions, work: &mut Work) -> Result<Self> {
        // These source shortcuts precede both global settings and nominal addition.
        if rhs.is_empty() {
            return self.copy_with_work(work);
        }
        if self.is_empty() {
            return rhs.copy_with_work(work);
        }
        let size = options.checked_size()?;
        let error = options.checked_error()?;
        let nominal_mass = self
            .nominal_mass
            .checked_add(rhs.nominal_mass)
            .ok_or_else(|| invalid("nominal isotope mass exceeds u32"))?;
        let pairs = size
            .checked_mul(size + 1)
            .and_then(|n| n.checked_div(2))
            .ok_or_else(|| invalid("isotope convolution work overflow"))?;
        work.consume(
            pairs
                .checked_mul(8)
                .ok_or_else(|| invalid("isotope convolution work overflow"))?,
        )?;
        let mut peaks = work.vector(size)?;
        // Virtual zero padding preserves arithmetic without copying either operand.
        for index in 0..size {
            let mut abundance = 0.0;
            let mut weighted_mass = 0.0;
            for left in 0..=index {
                let a = self.peaks.get(left).copied().unwrap_or_default();
                let b = rhs.peaks.get(index - left).copied().unwrap_or_default();
                abundance += a.abundance * b.abundance;
                weighted_mass += a.abundance * b.abundance * (a.mass + b.mass);
            }
            peaks.push(IMSIsotopePeak {
                mass: if abundance != 0.0 {
                    weighted_mass / abundance
                } else {
                    0.0
                },
                abundance,
            });
        }
        work.consume(size.saturating_mul(3))?;
        if let Some(scale) = normalization_scale(&peaks, error) {
            for peak in &mut peaks {
                peak.abundance *= scale;
            }
        }
        // A discarded zero-bin numerator need not be finite; every stored result must be.
        for peak in &peaks {
            validate_peak(*peak)?;
        }
        Ok(Self {
            peaks,
            nominal_mass,
        })
    }
}

fn normalization_scale(peaks: &[IMSIsotopePeak], error: f64) -> Option<f64> {
    let mut sum = 0.0;
    for peak in peaks {
        sum += peak.abundance;
    }
    if sum > 0.0 && (sum - 1.0).abs() > error {
        Some(1.0 / sum)
    } else {
        None
    }
}

fn validate_peak(peak: IMSIsotopePeak) -> Result<()> {
    finite(peak.mass, "stored isotope mass")?;
    finite(peak.abundance, "stored isotope abundance")?;
    Ok(())
}

fn check_count(len: usize) -> Result<()> {
    if len > MAX_PEAKS {
        return Err(invalid("isotope storage exceeds 1,000,000 peaks"));
    }
    Ok(())
}

pub(super) fn finite(value: f64, name: &str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(format!("{name} must be finite")))
    }
}

pub(super) fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

pub(super) struct Work {
    remaining: usize,
    bytes: usize,
}
impl Default for Work {
    fn default() -> Self {
        Self {
            remaining: MAX_WORK,
            bytes: MAX_BYTES,
        }
    }
}
impl Work {
    pub(super) fn consume(&mut self, amount: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(amount)
            .ok_or_else(|| invalid("IMS isotope operation exceeds 50,000,000 work units"))?;
        Ok(())
    }
    pub(super) fn allocate(&mut self, amount: usize) -> Result<()> {
        self.bytes = self.bytes.checked_sub(amount).ok_or_else(|| {
            invalid("IMS isotope operation exceeds 64 MiB logical allocation budget")
        })?;
        Ok(())
    }
    fn vector<T>(&mut self, len: usize) -> Result<Vec<T>> {
        self.consume(len)?;
        self.allocate(
            len.checked_mul(std::mem::size_of::<T>())
                .ok_or_else(|| invalid("isotope vector size overflow"))?,
        )?;
        let mut values = Vec::new();
        values
            .try_reserve_exact(len)
            .map_err(|_| invalid("cannot allocate isotope vector"))?;
        Ok(values)
    }
}

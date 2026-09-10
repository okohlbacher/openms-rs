// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source-ordered RNA fragmentation using independently declared nucleoside masses.
//! Single- and multiple-spectrum generation retain their distinct charge,
//! precursor, sulfur-linkage and metadata semantics from OpenMS 7c029e8.

use super::{EmpiricalFormula, NAFragmentType, NASequence, PROTON_MASS_U};
use crate::kernel::{DataArray, MSSpectrum, Peak1D};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::mem::size_of;

pub const MAX_RNA_SPECTRUM_RESIDUES: usize = 1_000_000;
/// Combined old/new single-spectrum peaks, or sum across all returned spectra.
pub const MAX_RNA_SPECTRUM_PEAKS: usize = 100_000;
pub const MAX_RNA_SPECTRUM_KEYS: usize = 4_096;
pub const MAX_RNA_SPECTRUM_ARRAYS: usize = 1_024;
pub const MAX_RNA_SPECTRUM_WORK: usize = 50_000_000;
pub const MAX_RNA_SPECTRUM_BYTES: usize = 256 * 1024 * 1024;
pub const MAX_RNA_SPECTRUM_LABEL_BYTES: usize = 16 * 1024 * 1024;
const ION_NAMES: &str = "IonNames";
const CHARGES: &str = "Charges";

/// All source settings. Only b/y ions default on; all intensities default to one.
/// Intensities stay f64 until each emitted peak is narrowed to f32. Unused
/// intensity settings are not validated; finite signed output is supported.
#[derive(Clone, Debug, PartialEq)]
pub struct NucleicAcidSpectrumGenerator {
    pub add_metainfo: bool,
    pub add_precursor_peaks: bool,
    pub add_all_precursor_charges: bool,
    pub add_first_prefix_ion: bool,
    pub add_a_ions: bool,
    pub add_b_ions: bool,
    pub add_c_ions: bool,
    pub add_d_ions: bool,
    pub add_w_ions: bool,
    pub add_x_ions: bool,
    pub add_y_ions: bool,
    pub add_z_ions: bool,
    pub add_a_minus_b_ions: bool,
    pub a_intensity: f64,
    pub b_intensity: f64,
    pub c_intensity: f64,
    pub d_intensity: f64,
    pub w_intensity: f64,
    pub x_intensity: f64,
    pub y_intensity: f64,
    pub z_intensity: f64,
    pub a_minus_b_intensity: f64,
    pub precursor_intensity: f64,
}
impl Default for NucleicAcidSpectrumGenerator {
    fn default() -> Self {
        Self {
            add_metainfo: false,
            add_precursor_peaks: false,
            add_all_precursor_charges: false,
            add_first_prefix_ion: false,
            add_a_ions: false,
            add_b_ions: true,
            add_c_ions: false,
            add_d_ions: false,
            add_w_ions: false,
            add_x_ions: false,
            add_y_ions: true,
            add_z_ions: false,
            add_a_minus_b_ions: false,
            a_intensity: 1.0,
            b_intensity: 1.0,
            c_intensity: 1.0,
            d_intensity: 1.0,
            w_intensity: 1.0,
            x_intensity: 1.0,
            y_intensity: 1.0,
            z_intensity: 1.0,
            a_minus_b_intensity: 1.0,
            precursor_intensity: 1.0,
        }
    }
}

impl NucleicAcidSpectrumGenerator {
    /// Creates a spectrum with source defaults (MS level one, unknown type).
    /// The M peak does not create an acquisition `Precursor` record.
    pub fn generate(
        &self,
        oligo: &NASequence,
        min_charge: i32,
        max_charge: i32,
    ) -> Result<MSSpectrum> {
        let mut spectrum = MSSpectrum::default();
        self.append_to(&mut spectrum, oligo, min_charge, max_charge)?;
        Ok(spectrum)
    }

    /// Appends and stably sorts all peaks, including when there are no new peaks.
    /// Source first-array names are retained. Missing first annotations for old
    /// peaks are padded with empty names/zero charge when additions require it.
    /// Other populated arrays cannot be extended without defined values and are
    /// rejected. Only peaks/arrays are staged; unrelated settings remain in place.
    pub fn append_to(
        &self,
        spectrum: &mut MSSpectrum,
        oligo: &NASequence,
        min_charge: i32,
        max_charge: i32,
    ) -> Result<()> {
        let mut work = Work::default();
        let (min, max, negative) = single_charges(min_charge, max_charge)?;
        let addition = self.single(oligo, min, max, negative, &mut work)?;
        let staged = stage_append(spectrum, &addition, self.add_metainfo, &mut work)?;
        staged.commit(spectrum);
        Ok(())
    }

    /// Source cumulative multi-charge generation. Charges need not be less than
    /// sequence length. The smallest key selects mode; skipped keys remain empty
    /// only with metadata enabled. Final-only M uses target charge +/- one, and
    /// the positive-mode final M deliberately has no absolute-value operation.
    /// Zero charges error only if a row would actually be divided by zero.
    pub fn generate_multiple(
        &self,
        oligo: &NASequence,
        charges: &BTreeSet<i32>,
        base_charge: i32,
    ) -> Result<BTreeMap<i32, MSSpectrum>> {
        let mut output = BTreeMap::new();
        let Some(&first) = charges.first() else {
            return Ok(output);
        };
        if charges.len() > MAX_RNA_SPECTRUM_KEYS {
            return Err(invalid("RNA spectrum key limit exceeded"));
        }
        let mut work = Work::default();
        work.allocate::<(i32, MSSpectrum)>(charges.len().saturating_mul(3).saturating_add(12))?;
        work.allocate::<u8>(1024)?;
        work.consume(charges.len().saturating_mul(128))?;
        if self.add_metainfo {
            for &charge in charges {
                output.insert(charge, empty_spectrum(true, &mut work)?);
            }
        }
        let uncharged = self.uncharged(oligo, &mut work)?;
        let negative = first < 0;
        let mut charge = if negative && base_charge > 0 {
            -base_charge
        } else {
            base_charge
        };
        let mut targets = collect_values(charges.iter().copied(), charges.len(), &mut work)?;
        if negative {
            targets.reverse();
        }
        let all_precursors = self.add_precursor_peaks && self.add_all_precursor_charges;
        let final_precursor = self.add_precursor_peaks && !self.add_all_precursor_charges;
        let mut cumulative = Rows::default();
        for target in targets {
            if (negative && target > charge) || (!negative && target < charge) {
                continue;
            }
            let next = if negative {
                target.checked_sub(1)
            } else {
                target.checked_add(1)
            }
            .ok_or_else(|| invalid("RNA charge counter overflows i32"))?;
            let span = usize::try_from((i64::from(target) - i64::from(charge)).unsigned_abs() + 1)
                .map_err(|_| invalid("RNA charge span overflows usize"))?;
            work.consume(span)?;
            let visible = self.visible_len(&uncharged, all_precursors);
            if visible != 0 {
                check_peaks(
                    cumulative
                        .peaks
                        .len()
                        .checked_add(
                            visible
                                .checked_mul(span)
                                .ok_or_else(|| invalid("RNA spectrum peak count overflows"))?,
                        )
                        .ok_or_else(|| invalid("RNA spectrum peak count overflows"))?,
                )?;
                while charge != next {
                    self.add_charged(
                        &mut cumulative,
                        &uncharged,
                        charge,
                        all_precursors,
                        &mut work,
                    )?;
                    charge = if negative { charge - 1 } else { charge + 1 };
                }
            } else {
                charge = next;
            }
            // The source copies to the next target before final M and before sort.
            // Keep an unsorted cumulative base and snapshot only the returned rows.
            let final_peak = if final_precursor {
                let precursor = *uncharged.peaks.last().ok_or_else(|| {
                    invalid("RNA final precursor requires a nonempty oligonucleotide")
                })?;
                let mz = charged_mass(precursor.mz, charge, negative)?;
                Some(Peak1D::new(mz, precursor.intensity))
            } else {
                None
            };
            let count = cumulative.peaks.len() + usize::from(final_peak.is_some());
            work.retain(count)?;
            let mut current = cumulative.copy(&mut work)?;
            if let Some(peak) = final_peak {
                current.push(
                    peak,
                    if self.add_metainfo { Some("M") } else { None },
                    charge,
                    &mut work,
                )?;
            }
            let mut spectrum = empty_spectrum(false, &mut work)?;
            let staged =
                stage_append_without_retain(&spectrum, &current, self.add_metainfo, &mut work)?;
            staged.commit(&mut spectrum);
            output.insert(target, spectrum);
        }
        Ok(output)
    }

    /// Atomic replacement form of the source map-output operation.
    pub fn replace_multiple(
        &self,
        spectra: &mut BTreeMap<i32, MSSpectrum>,
        oligo: &NASequence,
        charges: &BTreeSet<i32>,
        base_charge: i32,
    ) -> Result<()> {
        let result = self.generate_multiple(oligo, charges, base_charge)?;
        *spectra = result;
        Ok(())
    }

    fn single(
        &self,
        oligo: &NASequence,
        min: u32,
        max: u32,
        negative: bool,
        work: &mut Work,
    ) -> Result<Rows> {
        let uncharged = self.uncharged(oligo, work)?;
        let mut result = Rows::default();
        let stop = max.min(oligo.len().saturating_sub(1) as u32);
        if oligo.is_empty() || min > stop {
            return Ok(result);
        }
        work.consume((stop - min) as usize + 1)?;
        for z in min..=stop {
            let include = self.add_precursor_peaks && (self.add_all_precursor_charges || z == max);
            let charge = if negative { -(z as i32) } else { z as i32 };
            self.add_charged(&mut result, &uncharged, charge, include, work)?;
        }
        Ok(result)
    }
    fn visible_len(&self, rows: &Rows, include_precursor: bool) -> usize {
        if self.add_precursor_peaks && !include_precursor {
            rows.peaks.len().saturating_sub(1)
        } else {
            rows.peaks.len()
        }
    }
    fn add_charged(
        &self,
        output: &mut Rows,
        uncharged: &Rows,
        charge: i32,
        include_precursor: bool,
        work: &mut Work,
    ) -> Result<()> {
        let count = self.visible_len(uncharged, include_precursor);
        check_peaks(
            output
                .peaks
                .len()
                .checked_add(count)
                .ok_or_else(|| invalid("RNA peak count overflow"))?,
        )?;
        work.consume(count)?;
        for i in 0..count {
            let source = uncharged.peaks[i];
            let peak = Peak1D::new(charged_mass(source.mz, charge, true)?, source.intensity);
            output.push(
                peak,
                self.add_metainfo.then(|| uncharged.names[i].as_str()),
                charge,
                work,
            )?;
        }
        Ok(())
    }

    fn uncharged(&self, oligo: &NASequence, work: &mut Work) -> Result<Rows> {
        let mut output = Rows::default();
        if oligo.is_empty() {
            return Ok(output);
        }
        let n = oligo.len();
        if n > MAX_RNA_SPECTRUM_RESIDUES {
            return Err(invalid("RNA spectrum residue limit exceeded"));
        }
        // Source visits every record once to construct its mass/thiol tables.
        // Immutable record access lets us use indexed values without those copies.
        work.consume(n)?;
        let offsets = Offsets::new(work)?;
        let five = match oligo.five_prime_mod() {
            Some(end) => finite(end.mono_mass() - offsets.h)?,
            None => 0.0,
        };
        let three = match oligo.three_prime_mod() {
            Some(end) => finite(end.mono_mass() - offsets.h)?,
            None => 0.0,
        };
        let start = usize::from(!self.add_first_prefix_ion);
        let have_left = (self.add_a_ions
            || self.add_b_ions
            || self.add_c_ions
            || self.add_d_ions
            || self.add_a_minus_b_ions)
            && n > start + 1;
        let have_right =
            (self.add_w_ions || self.add_x_ions || self.add_y_ions || self.add_z_ions) && n > 1;
        let mut left = Vec::new();
        if have_left {
            work.allocate::<f64>(n - 1)?;
            left.try_reserve_exact(n - 1)
                .map_err(|_| invalid("RNA prefix allocation failed"))?;
            left.push(finite(oligo.residues()[0].mono_mass() + five)?);
            for i in 1..n - 1 {
                work.consume(4)?;
                left.push(finite(
                    finite(
                        finite(left[i - 1] + oligo.residues()[i].mono_mass())? + offsets.backbone,
                    )? + sulfur(oligo, i - 1, offsets.sulfur),
                )?);
            }
            for (enabled, letter, offset, intensity, extra_sulfur) in [
                (self.add_a_ions, "a", offsets.a, self.a_intensity, false),
                (self.add_b_ions, "b", 0.0, self.b_intensity, false),
                (
                    self.add_c_ions,
                    "c",
                    offsets.backbone,
                    self.c_intensity,
                    true,
                ),
                (
                    self.add_d_ions,
                    "d",
                    offsets.phosphate,
                    self.d_intensity,
                    true,
                ),
            ] {
                if enabled {
                    for (i, &mass) in left.iter().enumerate().skip(start) {
                        work.consume(2)?;
                        let mass = if extra_sulfur {
                            finite(mass + sulfur(oligo, i, offsets.sulfur))?
                        } else {
                            mass
                        };
                        self.emit(
                            &mut output,
                            finite(mass + offset)?,
                            intensity,
                            letter,
                            i + 1,
                            "",
                            work,
                        )?;
                    }
                }
            }
            if self.add_a_minus_b_ions {
                for i in start..left.len() {
                    let record = &oligo.residues()[i];
                    work.consume(record.baseloss_formula().atoms.len() + 4)?;
                    let mut mass = finite(record.baseloss_formula().mono_mass())?;
                    if i > 0 {
                        mass = finite(mass + finite(left[i - 1] + offsets.a_minus_b)?)?;
                        if record_before_is_star(oligo, i) {
                            mass = finite(mass + offsets.sulfur)?;
                        }
                    } else {
                        mass = finite(mass + offsets.first_a_minus_b)?;
                    }
                    let intensity = if record.is_ambiguous() {
                        self.a_minus_b_intensity * 0.5
                    } else {
                        self.a_minus_b_intensity
                    };
                    self.emit(&mut output, mass, intensity, "a", i + 1, "-B", work)?;
                    if record.is_ambiguous() {
                        self.emit(
                            &mut output,
                            finite(mass + offsets.methyl)?,
                            intensity,
                            "a",
                            i + 1,
                            "-B",
                            work,
                        )?;
                    }
                }
            }
        }
        let mut right = Vec::new();
        if have_right {
            work.allocate::<f64>(n - 1)?;
            right
                .try_reserve_exact(n - 1)
                .map_err(|_| invalid("RNA suffix allocation failed"))?;
            right.push(finite(oligo.residues()[n - 1].mono_mass() + three)?);
            for i in 1..n - 1 {
                work.consume(4)?;
                let r = n - i - 1;
                right.push(finite(
                    finite(
                        finite(right[i - 1] + oligo.residues()[r].mono_mass())? + offsets.backbone,
                    )? + sulfur(oligo, r, offsets.sulfur),
                )?);
            }
            for (enabled, letter, offset, intensity, extra_sulfur) in [
                (
                    self.add_w_ions,
                    "w",
                    offsets.phosphate,
                    self.w_intensity,
                    true,
                ),
                (
                    self.add_x_ions,
                    "x",
                    offsets.backbone,
                    self.x_intensity,
                    true,
                ),
                (self.add_y_ions, "y", 0.0, self.y_intensity, false),
                (self.add_z_ions, "z", offsets.a, self.z_intensity, false),
            ] {
                if enabled {
                    for (i, &mass) in right.iter().enumerate() {
                        work.consume(2)?;
                        let mass = if extra_sulfur {
                            finite(mass + sulfur(oligo, n - 2 - i, offsets.sulfur))?
                        } else {
                            mass
                        };
                        self.emit(
                            &mut output,
                            finite(mass + offset)?,
                            intensity,
                            letter,
                            i + 1,
                            "",
                            work,
                        )?;
                    }
                }
            }
        }
        if self.add_precursor_peaks {
            let mass = match (left.first(), right.last()) {
                (Some(&first), Some(&last)) => finite(finite(first + last)? + offsets.backbone)?,
                (Some(_), None) => finite(
                    finite(
                        finite(left[left.len() - 1] + oligo.residues()[n - 1].mono_mass())?
                            + offsets.backbone,
                    )? + three,
                )?,
                (None, Some(&last)) => finite(
                    finite(finite(last + oligo.residues()[0].mono_mass())? + offsets.backbone)?
                        + five,
                )?,
                (None, None) => oligo.mono_mass_with_budget(
                    NAFragmentType::Full,
                    0,
                    &mut work.remaining,
                    &mut work.bytes,
                )?,
            };
            output.push(
                Peak1D::new(mass, peak_intensity(self.precursor_intensity)?),
                self.add_metainfo.then_some("M"),
                0,
                work,
            )?;
        }
        Ok(output)
    }
    #[allow(clippy::too_many_arguments)]
    fn emit(
        &self,
        rows: &mut Rows,
        mass: f64,
        intensity: f64,
        letter: &str,
        ordinal: usize,
        suffix: &str,
        work: &mut Work,
    ) -> Result<()> {
        let name = if self.add_metainfo {
            let len = letter.len() + decimal_digits(ordinal) + suffix.len();
            work.text(len)?;
            Some(format!("{letter}{ordinal}{suffix}"))
        } else {
            None
        };
        rows.push_owned(Peak1D::new(mass, peak_intensity(intensity)?), name, 0, work)
    }
}

#[derive(Default)]
struct Rows {
    peaks: Vec<Peak1D>,
    names: Vec<String>,
    charges: Vec<i32>,
}
impl Rows {
    fn push(
        &mut self,
        peak: Peak1D,
        name: Option<&str>,
        charge: i32,
        work: &mut Work,
    ) -> Result<()> {
        let name = name.map(|text| work.copy_text(text)).transpose()?;
        self.push_owned(peak, name, charge, work)
    }
    fn push_owned(
        &mut self,
        peak: Peak1D,
        name: Option<String>,
        charge: i32,
        work: &mut Work,
    ) -> Result<()> {
        check_peaks(self.peaks.len() + 1)?;
        work.reserve(&mut self.peaks, 1)?;
        if name.is_some() {
            work.reserve(&mut self.names, 1)?;
            work.reserve(&mut self.charges, 1)?;
        }
        self.peaks.push(peak);
        if let Some(name) = name {
            self.names.push(name);
            self.charges.push(charge);
        }
        Ok(())
    }
    fn copy(&self, work: &mut Work) -> Result<Self> {
        work.allocate::<String>(self.names.len())?;
        let mut names = Vec::new();
        names
            .try_reserve_exact(self.names.len())
            .map_err(|_| invalid("RNA name vector allocation failed"))?;
        for name in &self.names {
            names.push(work.copy_text(name)?);
        }
        Ok(Self {
            peaks: copy_values(&self.peaks, work)?,
            names,
            charges: copy_values(&self.charges, work)?,
        })
    }
}

struct Offsets {
    h: f64,
    backbone: f64,
    a: f64,
    phosphate: f64,
    sulfur: f64,
    a_minus_b: f64,
    first_a_minus_b: f64,
    methyl: f64,
}
impl Offsets {
    fn new(work: &mut Work) -> Result<Self> {
        work.consume(128)?;
        work.allocate::<u8>(8192)?;
        let mass = |text: &str| -> Result<f64> { Ok(EmpiricalFormula::parse(text)?.mono_mass()) };
        Ok(Self {
            h: mass("H")?,
            backbone: mass("H-1PO2")?,
            a: -mass("H2O")?,
            phosphate: mass("HPO3")?,
            sulfur: mass("SO-1")?,
            a_minus_b: mass("H-5P")?,
            first_a_minus_b: -mass("H4O2")?,
            methyl: mass("CH2")?,
        })
    }
}
fn sulfur(oligo: &NASequence, index: usize, delta: f64) -> f64 {
    if oligo.residues()[index].code().ends_with('*') {
        delta
    } else {
        0.0
    }
}
fn record_before_is_star(oligo: &NASequence, index: usize) -> bool {
    oligo.residues()[index - 1].code().ends_with('*')
}
fn single_charges(min: i32, max: i32) -> Result<(u32, u32, bool)> {
    if min == i32::MIN || max == i32::MIN {
        return Err(invalid("RNA charge magnitude overflows i32"));
    }
    if (min < 0 && max > 0) || (min > 0 && max < 0) {
        return Err(invalid("RNA charge endpoints have different signs"));
    }
    let negative = min < 0 && max < 0;
    let a = min.unsigned_abs();
    let b = max.unsigned_abs();
    Ok((a.min(b), a.max(b), negative))
}
fn charged_mass(mass: f64, charge: i32, absolute: bool) -> Result<f64> {
    if charge == 0 {
        return Err(invalid("RNA peak charge cannot be zero"));
    }
    let value = finite(finite(mass / f64::from(charge))? + PROTON_MASS_U)?;
    Ok(if absolute { value.abs() } else { value })
}
fn peak_intensity(value: f64) -> Result<f32> {
    let narrowed = value as f32;
    if narrowed.is_finite() {
        Ok(narrowed)
    } else {
        Err(invalid("RNA peak intensity is nonfinite or exceeds f32"))
    }
}
fn decimal_digits(mut n: usize) -> usize {
    let mut digits = 1;
    while n >= 10 {
        n /= 10;
        digits += 1;
    }
    digits
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("RNA spectrum mass calculation is nonfinite"))
    }
}
fn check_peaks(n: usize) -> Result<()> {
    if n > MAX_RNA_SPECTRUM_PEAKS {
        Err(invalid("RNA spectrum peak limit exceeded"))
    } else {
        Ok(())
    }
}
fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

// Stage only the arrays' data. Existing names and unrelated spectrum metadata
// remain owned by the destination, even after a checked failure.
enum StagedArrays<T> {
    Existing(Vec<Vec<T>>),
    New(Vec<DataArray<T>>),
}
impl<T> StagedArrays<T> {
    fn commit(self, destination: &mut Vec<DataArray<T>>) {
        match self {
            Self::Existing(data) => {
                for (array, values) in destination.iter_mut().zip(data) {
                    array.data = values;
                }
            }
            Self::New(arrays) => *destination = arrays,
        }
    }
}
struct Staged {
    peaks: Vec<Peak1D>,
    floats: StagedArrays<f32>,
    integers: StagedArrays<i32>,
    strings: StagedArrays<String>,
}
impl Staged {
    fn commit(self, target: &mut MSSpectrum) {
        target.peaks = self.peaks;
        self.floats.commit(&mut target.float_data_arrays);
        self.integers.commit(&mut target.integer_data_arrays);
        self.strings.commit(&mut target.string_data_arrays);
    }
}
fn stage_append(old: &MSSpectrum, new: &Rows, metadata: bool, work: &mut Work) -> Result<Staged> {
    work.retain(
        old.len()
            .checked_add(new.peaks.len())
            .ok_or_else(|| invalid("RNA combined peak count overflow"))?,
    )?;
    stage_append_without_retain(old, new, metadata, work)
}
fn stage_append_without_retain(
    old: &MSSpectrum,
    new: &Rows,
    metadata: bool,
    work: &mut Work,
) -> Result<Staged> {
    let n = old
        .len()
        .checked_add(new.peaks.len())
        .ok_or_else(|| invalid("RNA combined peak count overflow"))?;
    check_peaks(n)?;
    let arrays = old
        .float_data_arrays
        .len()
        .saturating_add(old.integer_data_arrays.len())
        .saturating_add(old.string_data_arrays.len());
    let arrays = arrays
        .saturating_add(usize::from(metadata && old.integer_data_arrays.is_empty()))
        .saturating_add(usize::from(metadata && old.string_data_arrays.is_empty()));
    if arrays > MAX_RNA_SPECTRUM_ARRAYS {
        return Err(invalid("RNA spectrum array limit exceeded"));
    }
    work.consume(old.len() + arrays)?;
    for peak in &old.peaks {
        finite(peak.mz)?;
        if !peak.intensity.is_finite() {
            return Err(invalid("RNA existing peak intensity is nonfinite"));
        }
    }
    old.validate_data_arrays()?;
    let at = |i: usize| {
        if i < old.len() {
            old.peaks[i]
        } else {
            new.peaks[i - old.len()]
        }
    };
    let mut order = work.order(n)?;
    order.sort_by(|&a, &b| at(a).mz.partial_cmp(&at(b).mz).expect("finite peak m/z"));
    let peaks = collect_values(order.iter().map(|&i| at(i)), n, work)?;
    let floats = stage_numeric(
        &old.float_data_arrays,
        &[],
        None,
        old.len(),
        new.peaks.len(),
        &order,
        work,
    )?;
    let integers = stage_numeric(
        &old.integer_data_arrays,
        &new.charges,
        metadata.then_some(CHARGES),
        old.len(),
        new.peaks.len(),
        &order,
        work,
    )?;
    let strings = stage_strings(
        &old.string_data_arrays,
        &new.names,
        metadata,
        old.len(),
        new.peaks.len(),
        &order,
        work,
    )?;
    Ok(Staged {
        peaks,
        floats,
        integers,
        strings,
    })
}
#[allow(clippy::too_many_arguments)]
fn stage_numeric<T: Copy + Default>(
    arrays: &[DataArray<T>],
    added: &[T],
    first_name: Option<&str>,
    old_len: usize,
    added_len: usize,
    order: &[usize],
    work: &mut Work,
) -> Result<StagedArrays<T>> {
    let count = arrays.len().max(usize::from(first_name.is_some()));
    work.allocate::<Vec<T>>(count)?;
    let mut data = Vec::new();
    data.try_reserve_exact(count)
        .map_err(|_| invalid("RNA array allocation failed"))?;
    for i in 0..count {
        let previous = arrays.get(i).map_or(&[][..], |a| a.data.as_slice());
        let first = i == 0 && first_name.is_some();
        if added_len > 0 && !first && !previous.is_empty() {
            return Err(Error::Unsupported(
                "RNA append has no values for an existing populated array".into(),
            ));
        }
        if previous.is_empty() && (!first || added_len == 0) {
            data.push(Vec::new());
            continue;
        }
        data.push(collect_values(
            order.iter().map(|&j| {
                if j < old_len {
                    previous.get(j).copied().unwrap_or_default()
                } else {
                    added[j - old_len]
                }
            }),
            order.len(),
            work,
        )?);
    }
    if let Some(name) = first_name.filter(|_| arrays.is_empty()) {
        let mut result = Vec::new();
        work.allocate::<DataArray<T>>(1)?;
        result
            .try_reserve_exact(1)
            .map_err(|_| invalid("RNA array allocation failed"))?;
        result.push(DataArray::new(
            work.copy_text(name)?,
            data.pop().expect("first annotation data"),
        ));
        Ok(StagedArrays::New(result))
    } else {
        Ok(StagedArrays::Existing(data))
    }
}
fn stage_strings(
    arrays: &[DataArray<String>],
    added: &[String],
    metadata: bool,
    old_len: usize,
    added_len: usize,
    order: &[usize],
    work: &mut Work,
) -> Result<StagedArrays<String>> {
    let count = arrays.len().max(usize::from(metadata));
    work.allocate::<Vec<String>>(count)?;
    let mut data = Vec::new();
    data.try_reserve_exact(count)
        .map_err(|_| invalid("RNA string array allocation failed"))?;
    for i in 0..count {
        let previous = arrays.get(i).map_or(&[][..], |a| a.data.as_slice());
        let first = i == 0 && metadata;
        if added_len > 0 && !first && !previous.is_empty() {
            return Err(Error::Unsupported(
                "RNA append has no values for an existing populated array".into(),
            ));
        }
        if previous.is_empty() && (!first || added_len == 0) {
            data.push(Vec::new());
            continue;
        }
        work.allocate::<String>(order.len())?;
        let mut values = Vec::new();
        values
            .try_reserve_exact(order.len())
            .map_err(|_| invalid("RNA label allocation failed"))?;
        for &j in order {
            let text = if j < old_len {
                previous.get(j).map_or("", String::as_str)
            } else {
                added[j - old_len].as_str()
            };
            values.push(work.copy_text(text)?);
        }
        data.push(values);
    }
    if arrays.is_empty() && metadata {
        let mut result = Vec::new();
        work.allocate::<DataArray<String>>(1)?;
        result
            .try_reserve_exact(1)
            .map_err(|_| invalid("RNA array allocation failed"))?;
        result.push(DataArray::new(
            work.copy_text(ION_NAMES)?,
            data.pop().expect("first annotation data"),
        ));
        Ok(StagedArrays::New(result))
    } else {
        Ok(StagedArrays::Existing(data))
    }
}
fn copy_values<T: Copy>(source: &[T], work: &mut Work) -> Result<Vec<T>> {
    collect_values(source.iter().copied(), source.len(), work)
}
fn collect_values<T>(iter: impl Iterator<Item = T>, len: usize, work: &mut Work) -> Result<Vec<T>> {
    work.allocate::<T>(len)?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(len)
        .map_err(|_| invalid("RNA output allocation failed"))?;
    result.extend(iter);
    Ok(result)
}
fn empty_spectrum(metadata: bool, work: &mut Work) -> Result<MSSpectrum> {
    let mut result = MSSpectrum::default();
    if metadata {
        work.allocate::<DataArray<i32>>(1)?;
        work.allocate::<DataArray<String>>(1)?;
        result
            .integer_data_arrays
            .try_reserve_exact(1)
            .map_err(|_| invalid("RNA integer array allocation failed"))?;
        result
            .string_data_arrays
            .try_reserve_exact(1)
            .map_err(|_| invalid("RNA string array allocation failed"))?;
        result
            .integer_data_arrays
            .push(DataArray::new(work.copy_text(CHARGES)?, Vec::new()));
        result
            .string_data_arrays
            .push(DataArray::new(work.copy_text(ION_NAMES)?, Vec::new()));
    }
    Ok(result)
}

struct Work {
    remaining: usize,
    bytes: usize,
    label_bytes: usize,
    retained: usize,
}
impl Default for Work {
    fn default() -> Self {
        Self {
            remaining: MAX_RNA_SPECTRUM_WORK,
            bytes: MAX_RNA_SPECTRUM_BYTES,
            label_bytes: MAX_RNA_SPECTRUM_LABEL_BYTES,
            retained: 0,
        }
    }
}
impl Work {
    fn consume(&mut self, n: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(n)
            .ok_or_else(|| invalid("RNA spectrum work limit exceeded"))?;
        Ok(())
    }
    fn allocate<T>(&mut self, n: usize) -> Result<()> {
        self.consume(n)?;
        self.bytes = self
            .bytes
            .checked_sub(n.saturating_mul(size_of::<T>().max(1)))
            .ok_or_else(|| invalid("RNA spectrum allocation limit exceeded"))?;
        Ok(())
    }
    fn reserve<T>(&mut self, values: &mut Vec<T>, additional: usize) -> Result<()> {
        let needed = values
            .len()
            .checked_add(additional)
            .ok_or_else(|| invalid("RNA buffer size overflows"))?;
        if needed > values.capacity() {
            let capacity = needed.max(
                values
                    .capacity()
                    .saturating_mul(2)
                    .min(MAX_RNA_SPECTRUM_PEAKS),
            );
            self.allocate::<T>(capacity)?;
            values
                .try_reserve_exact(capacity - values.len())
                .map_err(|_| invalid("RNA buffer allocation failed"))?;
        }
        Ok(())
    }
    fn text(&mut self, n: usize) -> Result<()> {
        self.label_bytes = self
            .label_bytes
            .checked_sub(n)
            .ok_or_else(|| invalid("RNA spectrum label byte limit exceeded"))?;
        self.allocate::<u8>(n)
    }
    fn copy_text(&mut self, text: &str) -> Result<String> {
        self.text(text.len())?;
        Ok(text.to_owned())
    }
    fn retain(&mut self, n: usize) -> Result<()> {
        self.retained = self
            .retained
            .checked_add(n)
            .ok_or_else(|| invalid("RNA cumulative peak count overflows"))?;
        check_peaks(self.retained)
    }
    fn order(&mut self, n: usize) -> Result<Vec<usize>> {
        let depth = (usize::BITS - n.max(1).leading_zeros()) as usize;
        self.consume(n.saturating_mul(depth + 1).saturating_mul(4))?;
        // Index vector plus conservative full-length stable-sort scratch.
        self.allocate::<usize>(n.saturating_mul(2))?;
        let mut result = Vec::new();
        result
            .try_reserve_exact(n)
            .map_err(|_| invalid("RNA sort allocation failed"))?;
        result.extend(0..n);
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formula_only_precursor_consumes_the_existing_generator_budget() {
        let generator = NucleicAcidSpectrumGenerator {
            add_b_ions: false,
            add_y_ions: false,
            add_precursor_peaks: true,
            ..Default::default()
        };
        let oligo = NASequence::parse("AC").unwrap();
        assert!(generator.generate(&oligo, 1, 1).is_ok());
        // Enough for generator setup, but not the subsequent formula arithmetic.
        // A fresh budget inside the formula fallback would incorrectly succeed.
        let mut work = Work {
            remaining: 9_000,
            ..Default::default()
        };
        let error = match generator.uncharged(&oligo, &mut work) {
            Err(error) => error,
            Ok(_) => panic!("formula fallback reset the generator budget"),
        };
        assert!(
            error.to_string().contains("RNA operation work limit"),
            "{error}"
        );
        assert!(work.remaining < 9_000);
    }
}

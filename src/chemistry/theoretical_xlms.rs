// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source XLMS ladders, including the documented CPP-042/043 finite quirks.
//! See `docs/THEORETICAL_XLMS_SUPPORT.md`; this is not the ordinary TSG algorithm.
use super::{AASequence, PROTON_MASS_U, ProteinProteinCrossLink};
use crate::{Error, MSSpectrum, Peak1D, Result};
#[path = "theoretical_xlms_helpers.rs"]
mod helpers;
use helpers::{PreparedPeptide, publish};
const ISOTOPE: f64 = crate::concept::constants::C13C12_MASSDIFF_U;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct LossIndex {
    pub has_h2o_loss: bool,
    pub has_nh3_loss: bool,
}
impl LossIndex {
    fn union(self, other: Self) -> Self {
        Self {
            has_h2o_loss: self.has_h2o_loss || other.has_h2o_loss,
            has_nh3_loss: self.has_nh3_loss || other.has_nh3_loss,
        }
    }
}
/// All 25 source options. The two source TODO flags remain ineffective.
#[derive(Clone, Debug, PartialEq)]
pub struct XLMSOptions {
    pub add_isotopes: bool,
    pub max_isotope: i32,
    pub add_metainfo: bool,
    pub add_charges: bool,
    pub add_losses: bool,
    pub add_precursor_peaks: bool,
    pub add_abundant_immonium_ions: bool,
    pub add_k_linked_ions: bool,
    pub add_first_prefix_ion: bool,
    pub add_a_ions: bool,
    pub add_b_ions: bool,
    pub add_c_ions: bool,
    pub add_x_ions: bool,
    pub add_y_ions: bool,
    pub add_z_ions: bool,
    pub a_intensity: f64,
    pub b_intensity: f64,
    pub c_intensity: f64,
    pub x_intensity: f64,
    pub y_intensity: f64,
    pub z_intensity: f64,
    pub relative_loss_intensity: f64,
    pub precursor_intensity: f64,
    pub precursor_h2o_intensity: f64,
    pub precursor_nh3_intensity: f64,
}
impl Default for XLMSOptions {
    fn default() -> Self {
        Self {
            add_isotopes: false,
            max_isotope: 2,
            add_metainfo: true,
            add_charges: true,
            add_losses: false,
            add_precursor_peaks: true,
            add_abundant_immonium_ions: false,
            add_k_linked_ions: true,
            add_first_prefix_ion: true,
            add_a_ions: true,
            add_b_ions: true,
            add_c_ions: false,
            add_x_ions: false,
            add_y_ions: true,
            add_z_ions: false,
            a_intensity: 1.,
            b_intensity: 1.,
            c_intensity: 1.,
            x_intensity: 1.,
            y_intensity: 1.,
            z_intensity: 1.,
            relative_loss_intensity: 0.1,
            precursor_intensity: 1.,
            precursor_h2o_intensity: 1.,
            precursor_nh3_intensity: 1.,
        }
    }
}
/// Per-operation logical, cumulative allowances (not a physical RSS bound).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct XLMSLimits {
    pub max_residues: usize,
    pub max_peaks: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for XLMSLimits {
    fn default() -> Self {
        Self {
            max_residues: 4096,
            max_peaks: 100_000,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TheoreticalSpectrumGeneratorXLMS {
    pub options: XLMSOptions,
    pub limits: XLMSLimits,
}

struct Work {
    remaining: usize,
    bytes: usize,
}
impl Work {
    fn new(settings: &TheoreticalSpectrumGeneratorXLMS) -> Self {
        Self {
            remaining: settings.limits.max_work,
            bytes: settings.limits.max_bytes,
        }
    }
    fn charge(&mut self, work: usize, bytes: usize) -> Result<()> {
        self.remaining = self.remaining.checked_sub(work).ok_or_else(resource)?;
        self.bytes = self.bytes.checked_sub(bytes).ok_or_else(resource)?;
        Ok(())
    }
    fn slots<T>(&mut self, n: usize) -> Result<()> {
        self.charge(
            n,
            n.checked_mul(std::mem::size_of::<T>())
                .ok_or_else(resource)?,
        )
    }
}
fn invalid(s: &str) -> Error {
    Error::InvalidValue(format!("XLMS {s}"))
}
fn resource() -> Error {
    invalid("resource allowance exceeded")
}
fn finite(v: f64) -> Result<f64> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(invalid("nonfinite mass or intensity"))
    }
}
#[derive(Clone, Copy)]
enum Series {
    B,
    Y,
    A,
    X,
    C,
    Z,
}
impl Series {
    fn prefix(self) -> bool {
        matches!(self, Self::A | Self::B | Self::C)
    }
    fn letter(self) -> char {
        match self {
            Self::A => 'a',
            Self::B => 'b',
            Self::C => 'c',
            Self::X => 'x',
            Self::Y => 'y',
            Self::Z => 'z',
        }
    }
    fn correction(self) -> f64 {
        helpers::formula_mass(match self {
            Self::A => [-1, 0, 0, -1, 0, 0],
            Self::B => [0; 6],
            Self::C => [0, 3, 1, 0, 0, 0],
            Self::X => [1, 0, 0, 2, 0, 0],
            Self::Y => [0, 2, 0, 1, 0, 0],
            Self::Z => [0, -1, -1, 1, 0, 0],
        })
    }
}
struct Row {
    peak: Peak1D,
    charge: i32,
    name: String,
}
struct Emitter<'a> {
    options: &'a XLMSOptions,
    rows: Vec<Row>,
    limit: usize,
    work: &'a mut Work,
    water: f64,
    ammonia: f64,
}
impl<'a> Emitter<'a> {
    fn new(
        settings: &'a TheoreticalSpectrumGeneratorXLMS,
        old: usize,
        work: &'a mut Work,
    ) -> Result<Self> {
        let o = &settings.options;
        for value in [
            o.a_intensity,
            o.b_intensity,
            o.c_intensity,
            o.x_intensity,
            o.y_intensity,
            o.z_intensity,
            o.relative_loss_intensity,
            o.precursor_intensity,
            o.precursor_h2o_intensity,
            o.precursor_nh3_intensity,
        ] {
            finite(value)?;
        }
        work.charge(2048, 2048)?;
        Ok(Self {
            options: o,
            rows: Vec::new(),
            limit: settings
                .limits
                .max_peaks
                .checked_sub(old)
                .ok_or_else(resource)?,
            work,
            water: helpers::formula_mass([0, 2, 0, 1, 0, 0]),
            ammonia: helpers::formula_mass([0, 3, 1, 0, 0, 0]),
        })
    }
    fn isotope(&self) -> bool {
        self.options.add_isotopes && self.options.max_isotope >= 2
    }
    fn push(
        &mut self,
        mz: f64,
        intensity: f64,
        charge: i32,
        name: impl FnOnce() -> String,
    ) -> Result<()> {
        finite(mz)?;
        finite(intensity)?;
        let intensity = intensity as f32;
        finite(f64::from(intensity))?;
        if self.rows.len() >= self.limit {
            return Err(resource());
        }
        if self.rows.len() == self.rows.capacity() {
            let capacity = self
                .rows
                .capacity()
                .saturating_mul(2)
                .max(4)
                .min(self.limit);
            self.work.slots::<Row>(capacity)?;
            self.work.charge(self.rows.len(), 0)?;
            self.rows
                .try_reserve_exact(capacity - self.rows.len())
                .map_err(|_| resource())?;
        }
        self.work
            .charge(128, if self.options.add_metainfo { 128 } else { 0 })?;
        self.rows.push(Row {
            peak: Peak1D::new(mz, intensity),
            charge,
            name: if self.options.add_metainfo {
                name()
            } else {
                String::new()
            },
        });
        Ok(())
    }
    fn fragment(
        &mut self,
        mz: f64,
        intensity: f64,
        s: Series,
        index: usize,
        z: i32,
        kind: &str,
    ) -> Result<()> {
        finite(mz)?;
        if mz < 0. {
            return Ok(());
        }
        self.push(mz, intensity, z, || {
            format!("[{kind}${}{index}]", s.letter())
        })
    }
    #[allow(clippy::too_many_arguments)] // Explicit source ladder state; no generic builder.
    fn losses(
        &mut self,
        mass: f64,
        intensity: f64,
        z: i32,
        state: LossIndex,
        s: Series,
        index: usize,
        kind: &str,
    ) -> Result<()> {
        for (yes, loss, label) in [
            (state.has_h2o_loss, self.water, "H2O1"),
            (state.has_nh3_loss, self.ammonia, "H3N1"),
        ] {
            if yes {
                let m = finite(mass - loss)?;
                if m > 0. {
                    let p = divide(m, z)?;
                    self.push(
                        p,
                        finite(intensity * self.options.relative_loss_intensity)?,
                        z,
                        || format!("[{kind}${}{index}-{label}]", s.letter()),
                    )?;
                }
            }
        }
        Ok(())
    }
    #[allow(clippy::too_many_arguments)] // Explicit source ladder state; no generic builder.
    fn ladder(
        &mut self,
        p: &PreparedPeptide,
        first: usize,
        second: usize,
        alpha: bool,
        z: i32,
        s: Series,
        intensity: f64,
        precursor: Option<f64>,
        other: LossIndex,
    ) -> Result<()> {
        if p.masses.is_empty() {
            return Ok(());
        }
        if matches!(s, Series::C | Series::X) && p.masses.len() < 2 {
            return Err(invalid(
                "peptide must have at least two residues for c/x ions",
            ));
        }
        let kind = match (alpha, precursor.is_some()) {
            (true, false) => "alpha|ci",
            (false, false) => "beta|ci",
            (true, true) => "alpha|xi",
            (false, true) => "beta|xi",
        };
        let mut mass = finite(PROTON_MASS_U * f64::from(z))?;
        if let Some(precursor) = precursor {
            mass = finite(finite(mass + precursor)? - self.water)?;
            mass = finite(mass - if s.prefix() { p.c_term } else { p.n_term })?;
        } else {
            mass = finite(mass + if s.prefix() { p.n_term } else { p.c_term })?;
        }
        mass = finite(mass + p.correction(s))?;
        let n = p.masses.len();
        // A single explicit source-order loop keeps prefix/suffix subtraction
        // independent of a separately summed complementary peptide.
        let reverse = s.prefix() == precursor.is_some();
        let count = if reverse {
            n.saturating_sub(second.saturating_add(1))
        } else {
            first
        };
        for step in 0..count {
            let i = if reverse { n - 1 - step } else { step };
            self.work.charge(1, 0)?;
            let residue = *p
                .masses
                .get(i)
                .ok_or_else(|| invalid("link position exceeds peptide"))?;
            mass = finite(if precursor.is_some() {
                mass - residue
            } else {
                mass + residue
            })?;
            let pos = divide(mass, z)?;
            let index = match (s.prefix(), precursor.is_some()) {
                (true, false) => i + 1,
                (false, false) => n - i,
                (true, true) => i,
                (false, true) => n - 1 - i,
            };
            self.fragment(pos, intensity, s, index, z, kind)?;
            if self.options.add_losses {
                let loss_index = match (s.prefix(), precursor.is_some()) {
                    (true, false) => i,
                    (false, false) => i,
                    (true, true) => i - 1,
                    (false, true) => i + 1,
                };
                let states = if s.prefix() { &p.forward } else { &p.backward };
                // The source linked ladders guard the loss array independently
                // of their base/isotope peaks, including the length endpoint.
                if let Some(losses) = states.get(loss_index).copied() {
                    let losses = losses.union(other);
                    // CPP-042: source suffix linear losses divide an already divided m/z.
                    self.losses(
                        if !s.prefix() && precursor.is_none() {
                            pos
                        } else {
                            mass
                        },
                        intensity,
                        z,
                        losses,
                        s,
                        index,
                        kind,
                    )?;
                }
            }
            if self.isotope() {
                self.fragment(
                    finite(pos + divide(ISOTOPE, z)?)?,
                    intensity,
                    s,
                    index,
                    z,
                    kind,
                )?;
            }
        }
        Ok(())
    }
    fn precursor(&mut self, mass: f64, z: i32) -> Result<()> {
        for (loss, intensity, name) in [
            (0., self.options.precursor_intensity, "[M+H]"),
            (
                self.water,
                self.options.precursor_h2o_intensity,
                "[M+H]-H2O",
            ),
            (
                self.ammonia,
                self.options.precursor_nh3_intensity,
                "[M+H]-NH3",
            ),
        ] {
            let charged = finite(finite(mass + finite(PROTON_MASS_U * f64::from(z))?)? - loss)?;
            self.push(divide(charged, z)?, intensity, z, || name.into())?;
            // CPP-043: the source leaves the companion's charged mass undivided.
            if self.isotope() {
                self.push(finite(charged + divide(ISOTOPE, z)?)?, intensity, z, || {
                    name.into()
                })?;
            }
        }
        Ok(())
    }
    fn k_linked(
        &mut self,
        p: &PreparedPeptide,
        pos: usize,
        mass: f64,
        alpha: bool,
        z: i32,
    ) -> Result<()> {
        if pos == 0 {
            return Ok(());
        }
        let mut mass = finite(mass - p.span_mass(0, pos, Series::B)?)?;
        if pos >= p.masses.len() {
            return Ok(());
        }
        mass = finite(mass - p.span_mass(pos + 1, p.masses.len(), Series::X)?)?;
        mass = finite(mass + finite(PROTON_MASS_U * f64::from(z))?)?;
        if mass < 0. {
            return Ok(());
        }
        let mz = divide(mass, z)?;
        let letter = p.letters[pos] as char;
        self.push(mz, 1., z, || {
            format!("[{letter}-linked-{}]", if alpha { "beta" } else { "alpha" })
        })?;
        if self.isotope() {
            self.push(finite(mz + divide(ISOTOPE, z)?)?, 1., z, || {
                format!("[{letter}-linked-{}]", if alpha { "beta" } else { "alpha" })
            })?;
        }
        Ok(())
    }
    #[allow(clippy::too_many_arguments)] // Explicit source ladder state; no generic builder.
    fn generate(
        &mut self,
        p: &PreparedPeptide,
        positions: (usize, usize),
        alpha: bool,
        min: i32,
        max: i32,
        precursor: Option<f64>,
        other: LossIndex,
        k: bool,
        ladders: bool,
    ) -> Result<()> {
        let o = self.options;
        let series = [
            (o.add_b_ions, Series::B, o.b_intensity),
            (o.add_y_ions, Series::Y, o.y_intensity),
            (o.add_a_ions, Series::A, o.a_intensity),
            (o.add_x_ions, Series::X, o.x_intensity),
            (o.add_c_ions, Series::C, o.c_intensity),
            (o.add_z_ions, Series::Z, o.z_intensity),
        ];
        let charges = if min <= max {
            usize::try_from(i64::from(max) - i64::from(min) + 1).map_err(|_| resource())?
        } else {
            0
        };
        let per = p.masses.len().saturating_add(1).saturating_mul(16);
        self.work
            .charge(charges.checked_mul(per).ok_or_else(resource)?, 0)?;
        for z in min..=max {
            for (enabled, s, intensity) in series {
                if enabled && ladders {
                    self.ladder(
                        p,
                        positions.0,
                        positions.1,
                        alpha,
                        z,
                        s,
                        intensity,
                        precursor,
                        other,
                    )?;
                }
            }
            if k {
                self.k_linked(p, positions.0, precursor.unwrap(), alpha, z)?;
            }
        }
        Ok(())
    }
}
fn divide(mass: f64, z: i32) -> Result<f64> {
    if z == 0 {
        Err(invalid("zero charge division"))
    } else {
        finite(mass / f64::from(z))
    }
}
impl TheoreticalSpectrumGeneratorXLMS {
    pub fn new() -> Self {
        Self::default()
    }
    /// Append linear fragments. Source defaults: charge=1, link_pos_2=0.
    pub fn get_linear_ion_spectrum(
        &self,
        spectrum: &mut MSSpectrum,
        peptide: &AASequence,
        link_pos: usize,
        frag_alpha: bool,
        charge: i32,
        link_pos_2: usize,
    ) -> Result<()> {
        let mut work = Work::new(self);
        let mut e = Emitter::new(self, spectrum.len(), &mut work)?;
        let p = PreparedPeptide::new(peptide, self, e.work)?;
        let second = if link_pos_2 == 0 {
            link_pos
        } else {
            link_pos_2
        };
        e.generate(
            &p,
            (link_pos, second),
            frag_alpha,
            1,
            charge,
            None,
            LossIndex::default(),
            false,
            true,
        )?;
        publish(spectrum, e)
    }
    /// Append linked fragments using an explicit uncharged precursor mass.
    #[allow(clippy::too_many_arguments)] // Complete source overload, explicit arguments.
    pub fn get_xlink_ion_spectrum(
        &self,
        spectrum: &mut MSSpectrum,
        peptide: &AASequence,
        link_pos: usize,
        precursor_mass: f64,
        frag_alpha: bool,
        min_charge: i32,
        max_charge: i32,
        link_pos_2: usize,
    ) -> Result<()> {
        finite(precursor_mass)?;
        let mut work = Work::new(self);
        let mut e = Emitter::new(self, spectrum.len(), &mut work)?;
        let p = PreparedPeptide::new(peptide, self, e.work)?;
        let second = if link_pos_2 == 0 {
            link_pos
        } else {
            link_pos_2
        };
        e.generate(
            &p,
            (link_pos, second),
            frag_alpha,
            min_charge,
            max_charge,
            Some(precursor_mass),
            LossIndex::default(),
            self.options.add_k_linked_ions,
            true,
        )?;
        if self.options.add_precursor_peaks {
            e.precursor(precursor_mass, max_charge)?;
        }
        publish(spectrum, e)
    }
    /// Append one side of an owned crosslink. A missing alpha is a source no-op.
    pub fn get_crosslink_ion_spectrum(
        &self,
        spectrum: &mut MSSpectrum,
        crosslink: &ProteinProteinCrossLink,
        frag_alpha: bool,
        min_charge: i32,
        max_charge: i32,
    ) -> Result<()> {
        self.crosslink_with_work(
            spectrum,
            crosslink,
            frag_alpha,
            min_charge,
            max_charge,
            &mut Work::new(self),
        )
    }
    /// Composition-only adapter; both remaining counters propagate on errors.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn crosslink_with_budget(
        &self,
        spectrum: &mut MSSpectrum,
        crosslink: &ProteinProteinCrossLink,
        frag_alpha: bool,
        min_charge: i32,
        max_charge: i32,
        remaining_work: &mut usize,
        remaining_bytes: &mut usize,
    ) -> Result<()> {
        let initial_work = (*remaining_work).min(self.limits.max_work);
        let initial_bytes = (*remaining_bytes).min(self.limits.max_bytes);
        let mut work = Work {
            remaining: initial_work,
            bytes: initial_bytes,
        };
        let result = self.crosslink_with_work(
            spectrum, crosslink, frag_alpha, min_charge, max_charge, &mut work,
        );
        *remaining_work -= initial_work - work.remaining;
        *remaining_bytes -= initial_bytes - work.bytes;
        result
    }
    fn crosslink_with_work(
        &self,
        spectrum: &mut MSSpectrum,
        crosslink: &ProteinProteinCrossLink,
        frag_alpha: bool,
        min_charge: i32,
        max_charge: i32,
        work: &mut Work,
    ) -> Result<()> {
        let Some(alpha) = &crosslink.alpha else {
            return Ok(());
        };
        let mut e = Emitter::new(self, spectrum.len(), work)?;
        let alpha = PreparedPeptide::new(alpha, self, e.work)?;
        let empty = AASequence::default();
        let beta = PreparedPeptide::new(crosslink.beta.as_deref().unwrap_or(&empty), self, e.work)?;
        let series = self.options.add_a_ions
            || self.options.add_b_ions
            || self.options.add_c_ions
            || self.options.add_x_ions
            || self.options.add_y_ions
            || self.options.add_z_ions;
        let needed = self.options.add_precursor_peaks
            || min_charge <= max_charge
                && ((!alpha.masses.is_empty() && series)
                    || (self.options.add_k_linked_ions && !beta.masses.is_empty()));
        let mass = if needed {
            e.work.charge(
                alpha
                    .masses
                    .len()
                    .saturating_add(beta.masses.len())
                    .saturating_add(4),
                0,
            )?;
            let mass = finite(alpha.full_mass()? + crosslink.cross_linker_mass())?;
            if beta.masses.is_empty() {
                mass
            } else {
                finite(mass + beta.full_mass()?)?
            }
        } else {
            0.
        };
        let (p, other, pos) = if frag_alpha {
            (&alpha, &beta, crosslink.cross_link_position.0)
        } else {
            (&beta, &alpha, crosslink.cross_link_position.1)
        };
        // Source converts SignedSize to Size. Checked slice access below turns
        // actual out-of-range consumption into an error, preserving empty loops.
        let position = pos as usize;
        if !frag_alpha
            && p.masses.is_empty()
            && !alpha.masses.is_empty()
            && series
            && min_charge <= max_charge
        {
            return Err(invalid(
                "empty fragmented beta would index outside the source peptide",
            ));
        }
        let losses = if self.options.add_losses {
            *other
                .backward
                .first()
                .ok_or_else(|| invalid("empty peptide has no source loss index"))?
        } else {
            LossIndex::default()
        };
        e.generate(
            p,
            (position, position),
            frag_alpha,
            min_charge,
            max_charge,
            Some(mass),
            losses,
            self.options.add_k_linked_ions && !beta.masses.is_empty(),
            !alpha.masses.is_empty(),
        )?;
        if self.options.add_precursor_peaks {
            e.precursor(mass, max_charge)?;
        }
        publish(spectrum, e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mass_preparation_uses_one_shared_allowance() {
        let p = AASequence::parse("ACDE").unwrap();
        let g = TheoreticalSpectrumGeneratorXLMS::default();
        let mut probe = Work {
            remaining: 1_000_000,
            bytes: 1_000_000,
        };
        PreparedPeptide::new(&p, &g, &mut probe).unwrap();
        let used = 1_000_000 - probe.remaining;
        let mut shared = Work {
            remaining: used,
            bytes: 1_000_000,
        };
        PreparedPeptide::new(&p, &g, &mut shared).unwrap();
        assert_eq!(shared.remaining, 0);
        assert!(PreparedPeptide::new(&p, &g, &mut shared).is_err());
    }
    #[test]
    fn negative_filtered_rows_still_consume_charge_work() {
        let mut g = TheoreticalSpectrumGeneratorXLMS::default();
        g.options.add_a_ions = false;
        g.options.add_y_ions = false;
        g.options.add_precursor_peaks = false;
        g.options.add_k_linked_ions = false;
        let mut allowance = Work::new(&g);
        let mut e = Emitter::new(&g, 0, &mut allowance).unwrap();
        let p = PreparedPeptide::new(&AASequence::parse("AAAA").unwrap(), &g, e.work).unwrap();
        e.work.remaining = (p.masses.len() + 1) * 16;
        assert!(
            e.generate(
                &p,
                (0, 0),
                true,
                1,
                2,
                Some(-1e6),
                LossIndex::default(),
                false,
                true
            )
            .is_err()
        );
        assert!(e.rows.is_empty());
    }
    #[test]
    fn late_permutation_allocation_error_is_atomic() {
        let g = TheoreticalSpectrumGeneratorXLMS::default();
        let mut target = MSSpectrum::from_peaks(vec![Peak1D::new(100., 3.)]);
        let before = target.clone();
        let mut allowance = Work::new(&g);
        let mut e = Emitter::new(&g, target.len(), &mut allowance).unwrap();
        e.push(20., 1., 1, || "test".into()).unwrap();
        e.work.bytes = 1;
        assert!(publish(&mut target, e).is_err());
        assert_eq!(target, before);
    }
    #[test]
    fn second_side_failure_preserves_first_side_and_shared_counters() {
        use std::sync::Arc;
        let mut link = ProteinProteinCrossLink::default();
        link.alpha = Some(Arc::new(AASequence::parse("AMAA").unwrap()));
        link.beta = Some(Arc::new(AASequence::parse("AAMA").unwrap()));
        link.cross_link_position = (1, 2);
        let g = TheoreticalSpectrumGeneratorXLMS::default();
        let mut work = 50_000_000;
        let mut bytes = 64 * 1024 * 1024;
        let mut first = MSSpectrum::default();
        g.crosslink_with_budget(&mut first, &link, true, 1, 2, &mut work, &mut bytes)
            .unwrap();
        let work_after_first = work;
        let bytes_after_first = bytes;
        let mut complete = first.clone();
        g.crosslink_with_budget(&mut complete, &link, false, 1, 2, &mut work, &mut bytes)
            .unwrap();
        assert!(complete.len() > first.len());
        let used = [work_after_first - work, bytes_after_first - bytes];
        for index in 0..2 {
            let mut target = first.clone();
            let mut work = work_after_first;
            let mut bytes = bytes_after_first;
            if index == 0 {
                work = used[0] - 1
            } else {
                bytes = used[1] - 1
            }
            let before = (work, bytes);
            assert!(
                g.crosslink_with_budget(&mut target, &link, false, 1, 2, &mut work, &mut bytes)
                    .is_err()
            );
            assert_eq!(target, first);
            assert!(work < before.0 && bytes < before.1);
        }
    }
}

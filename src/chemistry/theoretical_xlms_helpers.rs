// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use super::*;
use crate::chemistry::{SequenceModification, composition_formula, residue_composition};
use crate::kernel::DataArray;

pub(super) fn formula_mass(counts: [i32; 6]) -> f64 {
    composition_formula(counts).mono_mass()
}
pub(super) struct PreparedPeptide {
    pub letters: Vec<u8>,
    pub masses: Vec<f64>,
    pub n_term: f64,
    pub c_term: f64,
    pub forward: Vec<LossIndex>,
    pub backward: Vec<LossIndex>,
    bare_x: Vec<bool>,
    corrections: [f64; 6],
    water: f64,
}
impl PreparedPeptide {
    pub fn new(
        p: &AASequence,
        settings: &TheoreticalSpectrumGeneratorXLMS,
        w: &mut Work,
    ) -> Result<Self> {
        if p.len() > settings.limits.max_residues {
            return Err(resource());
        }
        w.charge(p.len().saturating_add(8192), 8192)?;
        let payload = p.generation_payload_bytes()?;
        w.charge(payload, payload)?;
        w.slots::<(u8, f64, bool, LossIndex, LossIndex)>(p.len())?;
        let letters = p.as_str().as_bytes().to_vec();
        let water = composition_formula([0, 2, 0, 1, 0, 0]);
        let water_mass = water.mono_mass();
        let mut masses = Vec::new();
        masses.try_reserve_exact(p.len()).map_err(|_| resource())?;
        let mut bare_x = Vec::new();
        bare_x.try_reserve_exact(p.len()).map_err(|_| resource())?;
        for (i, &letter) in letters.iter().enumerate() {
            let modification = p.residue_modification(i)?;
            bare_x.push(letter == b'X' && modification.is_none());
            // Source Residue stores a free amino-acid formula/mass, then removes
            // water. General AASequence's internal formula sum is not identical.
            w.charge(4096, 4096)?;
            let base = match residue_composition(letter) {
                Some(c) => composition_formula(c).checked_add(&water)?,
                None if matches!(letter, b'B' | b'Z' | b'X') => Default::default(),
                _ => return Err(invalid("unsupported residue mass")),
            };
            let mut mass = base.mono_mass();
            match modification {
                Some(SequenceModification::MassTag(tag)) => {
                    mass = if let Some(delta) = tag.delta_mono_mass() {
                        finite(mass + delta)?
                    } else {
                        finite(
                            tag.residue_mono_mass()
                                .ok_or_else(|| invalid("unavailable residue mass tag"))?
                                + water_mass,
                        )?
                    };
                }
                Some(SequenceModification::Known(m)) => {
                    let entries = base
                        .atoms
                        .len()
                        .saturating_add(m.diff_formula().atoms.len())
                        .saturating_add(m.absolute_formula().map_or(0, |f| f.atoms.len()));
                    w.charge(
                        entries.saturating_mul(256),
                        entries.saturating_mul(512).saturating_add(2048),
                    )?;
                    if m.diff_mono_mass() != 0.
                        || m.diff_average_mass() != 0.
                        || !m.diff_formula().is_empty()
                    {
                        if !m.diff_formula().is_empty() {
                            mass = base.checked_add(m.diff_formula())?.mono_mass();
                        } else if let Some(f) = m.absolute_formula() {
                            mass = f.mono_mass();
                        } else if m.mono_mass() != 0. {
                            mass = m.mono_mass();
                        } else if m.diff_mono_mass() != 0. {
                            mass = finite(mass + m.diff_mono_mass())?;
                        }
                    }
                }
                None => {}
            }
            masses.push(finite(mass - water_mass)?);
        }
        let n_term = p
            .n_terminal_modification()
            .map(SequenceModification::diff_mono_mass)
            .transpose()?
            .unwrap_or(0.);
        let c_term = p
            .c_terminal_modification()
            .map(SequenceModification::diff_mono_mass)
            .transpose()?
            .unwrap_or(0.);
        let mut forward = Vec::new();
        let mut backward = Vec::new();
        if settings.options.add_losses {
            if p.is_empty() {
                return Err(invalid("empty peptide has no source loss index"));
            }
            forward.try_reserve_exact(p.len()).map_err(|_| resource())?;
            backward
                .try_reserve_exact(p.len())
                .map_err(|_| resource())?;
            let mut cumulative = LossIndex::default();
            for &c in &letters {
                if !b"RHKDESTNQCUGPAVILMFYW".contains(&c) {
                    return Err(invalid("residue missing from source loss alphabet"));
                }
                let state = LossIndex {
                    has_h2o_loss: b"DEST".contains(&c),
                    has_nh3_loss: b"KNQR".contains(&c),
                };
                cumulative = cumulative.union(state);
                forward.push(cumulative);
                backward.push(state);
            }
            for i in (0..p.len() - 1).rev() {
                backward[i] = backward[i].union(backward[i + 1]);
            }
        }
        Ok(Self {
            letters,
            masses,
            n_term,
            c_term,
            forward,
            backward,
            bare_x,
            water: water_mass,
            corrections: [
                Series::B,
                Series::Y,
                Series::A,
                Series::X,
                Series::C,
                Series::Z,
            ]
            .map(Series::correction),
        })
    }
    pub fn correction(&self, s: Series) -> f64 {
        self.corrections[match s {
            Series::B => 0,
            Series::Y => 1,
            Series::A => 2,
            Series::X => 3,
            Series::C => 4,
            Series::Z => 5,
        }]
    }
    pub fn full_mass(&self) -> Result<f64> {
        if self.masses.is_empty() {
            return Ok(0.);
        }
        let mut mass = finite(finite(0. + self.n_term)? + self.c_term)?;
        for (&m, &x) in self.masses.iter().zip(&self.bare_x) {
            if x {
                return Err(invalid("source AASequence mass rejects bare X"));
            }
            mass = finite(mass + m)?;
        }
        finite(mass + self.water)
    }
    pub fn span_mass(&self, start: usize, end: usize, s: Series) -> Result<f64> {
        let values = self
            .masses
            .get(start..end)
            .ok_or_else(|| invalid("link position exceeds peptide"))?;
        if values.is_empty() {
            return Ok(0.);
        }
        let mut mass = 0.;
        if start == 0 && s.prefix() {
            mass = finite(mass + self.n_term)?;
        }
        if end == self.masses.len() && !s.prefix() {
            mass = finite(mass + self.c_term)?;
        }
        for (&m, &x) in values.iter().zip(&self.bare_x[start..end]) {
            if x {
                return Err(invalid("source AASequence mass rejects bare X"));
            }
            mass = finite(mass + m)?;
        }
        finite(mass + self.correction(s))
    }
}

struct Arrays<T> {
    data: Vec<Vec<T>>,
    outer: Vec<DataArray<T>>,
    name: Option<String>,
}
impl<T: Clone + Default> Arrays<T> {
    #[allow(clippy::too_many_arguments)] // Atomic array staging uses one explicit shared budget.
    fn prepare(
        old: &[DataArray<T>],
        old_len: usize,
        indices: &[usize],
        rows: &[Row],
        name: Option<&str>,
        new: impl Fn(&Row) -> T,
        cost: impl Fn(&T) -> usize,
        w: &mut Work,
    ) -> Result<Self> {
        let count = old.len().max(usize::from(name.is_some()));
        w.slots::<Vec<T>>(count)?;
        w.slots::<DataArray<T>>(count)?;
        let mut data = Vec::new();
        data.try_reserve_exact(count).map_err(|_| resource())?;
        let mut outer = Vec::new();
        outer.try_reserve_exact(count).map_err(|_| resource())?;
        for i in 0..count {
            let previous = old.get(i).map_or(&[][..], |v| v.data.as_slice());
            if !previous.is_empty() && previous.len() != old_len {
                return Err(invalid("unaligned existing data array"));
            }
            let annotate = i == 0 && name.is_some();
            if !annotate && !previous.is_empty() && !rows.is_empty() {
                return Err(invalid("cannot extend unrelated populated data array"));
            }
            if previous.is_empty() && !annotate {
                data.push(Vec::new());
                continue;
            }
            w.slots::<T>(indices.len())?;
            for value in previous {
                let bytes = cost(value);
                w.charge(bytes, bytes)?;
            }
            // Each generated annotation is bounded before its original format
            // call. Reserve that same spelling bound for this final permutation.
            if annotate {
                w.charge(
                    rows.len().saturating_mul(128),
                    rows.len().saturating_mul(128),
                )?;
            }
            let mut values = Vec::new();
            values
                .try_reserve_exact(indices.len())
                .map_err(|_| resource())?;
            for &index in indices {
                values.push(if index < old_len {
                    previous.get(index).cloned().unwrap_or_default()
                } else {
                    new(&rows[index - old_len])
                });
            }
            data.push(values);
        }
        if let Some(name) = name {
            w.charge(
                name.len()
                    .saturating_add(old.first().map_or(0, |a| a.name.len())),
                name.len(),
            )?;
        }
        let name = name.map(|name| name.to_owned());
        Ok(Self { data, outer, name })
    }
    fn apply(mut self, target: &mut Vec<DataArray<T>>) {
        let mut old = std::mem::take(target).into_iter();
        for data in self.data {
            let mut array = old.next().unwrap_or_else(|| DataArray::new("", Vec::new()));
            array.data = data;
            if self.outer.is_empty() {
                if let Some(name) = self.name.take() {
                    array.name = name;
                }
            }
            self.outer.push(array);
        }
        *target = self.outer;
    }
}
pub(super) fn publish(target: &mut MSSpectrum, mut emitter: Emitter<'_>) -> Result<()> {
    let w = &mut emitter.work;
    let old_len = target.len();
    let total = old_len
        .checked_add(emitter.rows.len())
        .ok_or_else(resource)?;
    w.slots::<Peak1D>(old_len)?;
    for peak in &target.peaks {
        finite(peak.mz)?;
        finite(f64::from(peak.intensity))?;
    }
    let arrays = target
        .float_data_arrays
        .len()
        .saturating_add(target.integer_data_arrays.len())
        .saturating_add(target.string_data_arrays.len());
    w.charge(arrays, 0)?;
    // Description objects are moved, not traversed/cloned. Old array data and
    // renamed first-array text are metered before their replacement/destruction.
    let levels = (usize::BITS - total.max(1).leading_zeros()) as usize;
    w.charge(total.saturating_mul(levels.saturating_add(4)), 0)?;
    w.slots::<usize>(total.saturating_mul(2))?;
    w.slots::<Peak1D>(total)?;
    let mut indices = Vec::new();
    indices.try_reserve_exact(total).map_err(|_| resource())?;
    indices.extend(0..total);
    let peak = |i: usize| {
        if i < old_len {
            target.peaks[i]
        } else {
            emitter.rows[i - old_len].peak
        }
    };
    indices.sort_by(|&a, &b| peak(a).mz.partial_cmp(&peak(b).mz).unwrap());
    let mut peaks = Vec::new();
    peaks.try_reserve_exact(total).map_err(|_| resource())?;
    peaks.extend(indices.iter().map(|&i| peak(i)));
    let floats = Arrays::prepare(
        &target.float_data_arrays,
        old_len,
        &indices,
        &emitter.rows,
        None,
        |_| 0.,
        |_| 0,
        w,
    )?;
    let integers = Arrays::prepare(
        &target.integer_data_arrays,
        old_len,
        &indices,
        &emitter.rows,
        emitter.options.add_charges.then_some("charge"),
        |r| r.charge,
        |_| 0,
        w,
    )?;
    let strings = Arrays::prepare(
        &target.string_data_arrays,
        old_len,
        &indices,
        &emitter.rows,
        emitter.options.add_metainfo.then_some("IonNames"),
        |r| r.name.clone(),
        String::len,
        w,
    )?;
    // No fallible operations remain. All untouched spectrum fields and array
    // descriptions keep their original ownership and Arc identities.
    target.peaks = peaks;
    floats.apply(&mut target.float_data_arrays);
    integers.apply(&mut target.integer_data_arrays);
    strings.apply(&mut target.string_data_arrays);
    Ok(())
}

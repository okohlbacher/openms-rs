// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Conversions between OpenMS kernel records and the OpenSWATH light records.
//!
//! Ports the spectrum and chromatogram half of
//! `ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h` (class
//! `OpenSwathDataAccessHelper`, whose static members become free functions).
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md).
//!
//! `convertPeptideToAASequence` is [`convert_peptide_to_aa_sequence`](crate::analysis::openswath::data_access_helper::convert_peptide_to_aa_sequence), added
//! with the OpenSWATH scoring port. `convertTargetedExp` and the two
//! `convertTargetedCompound` overloads, which read the heavyweight
//! `ANALYSIS/TARGETED/TargetedExperiment`, are
//! [`convert_targeted_exp`](crate::analysis::openswath::data_access_helper::convert_targeted_exp),
//! [`convert_targeted_peptide`](crate::analysis::openswath::data_access_helper::convert_targeted_peptide)
//! and [`convert_targeted_compound`](crate::analysis::openswath::data_access_helper::convert_targeted_compound),
//! added with the MRM feature finder port.
//!
//! # Lengths
//!
//! The source walks the m/z (time) array and reads the intensity array at the
//! same index without checking its length. The port checks, and returns
//! [`Error::InvalidValue`](crate::Error::InvalidValue) before modifying the destination.

use super::data_structures::{Chromatogram, Spectrum};
use super::mrm_feature_support::data_value_to_string;
use super::reaction_monitoring_transition::DecoyTransitionType;
use super::targeted_experiment::TargetedExperiment;
use super::targeted_experiment_helper::{
    Compound, Peptide, PeptideCompound, RTUnit, RetentionTime, get_aa_sequence, set_modification,
};
use super::transition_experiment::{
    LightCompound, LightModification, LightProtein, LightTargetedExperiment, LightTransition,
};
use crate::chemistry::AASequence;
use crate::kernel::{ChromatogramPeak, MSChromatogram, MSSpectrum, Peak1D};
use crate::{Error, Result};

fn check_pair(position: usize, intensity: usize, what: &str) -> Result<()> {
    if intensity != position {
        return Err(Error::InvalidValue(format!(
            "the {what} array has {position} entries but the intensity array {intensity}"
        )));
    }
    Ok(())
}

/// Convert an OpenSWATH spectrum into the peaks of an OpenMS spectrum.
///
/// The peaks of `spectrum` are replaced by one peak per m/z value; its
/// metadata is kept. As the kernel's `clear(false)` does, the data arrays of
/// `spectrum` are cleared too, since they would no longer align with the new
/// peaks; further arrays of `sptr` (such as ion mobility) are not carried over,
/// as in the source. Intensities are narrowed to `f32`, the kernel's intensity
/// type, as in the source.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the intensity array is not as long as the m/z
/// array; `spectrum` is then unchanged.
pub fn convert_to_openms_spectrum(sptr: &Spectrum, spectrum: &mut MSSpectrum) -> Result<()> {
    let mz = &sptr.mz_array().data;
    let intensity = &sptr.intensity_array().data;
    check_pair(mz.len(), intensity.len(), "m/z")?;
    spectrum.clear(false);
    spectrum.peaks.reserve(mz.len());
    spectrum.peaks.extend(
        mz.iter()
            .zip(intensity)
            .map(|(&m, &i)| Peak1D::new(m, i as f32)),
    );
    Ok(())
}

/// Convert an OpenMS spectrum into an OpenSWATH spectrum holding only its m/z
/// and intensity arrays. Source `convertToSpectrumPtr`.
///
/// Data arrays are **not** carried over; the spectrum access view uses
/// [`spectrum_with_data_arrays`](super::spectrum_access_openms::spectrum_with_data_arrays),
/// which does carry them.
pub fn convert_to_spectrum(spectrum: &MSSpectrum) -> Spectrum {
    Spectrum::from_arrays(
        spectrum.peaks.iter().map(|p| p.mz).collect(),
        spectrum
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect(),
    )
}

/// Convert an OpenMS chromatogram into an OpenSWATH chromatogram holding only
/// its time and intensity arrays. Source `convertToChromatogramPtr`.
pub fn convert_to_chromatogram(chromatogram: &MSChromatogram) -> Chromatogram {
    Chromatogram::from_arrays(
        chromatogram.peaks.iter().map(|p| p.rt).collect(),
        chromatogram
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect(),
    )
}

/// Convert an OpenSWATH chromatogram into the peaks of an OpenMS chromatogram.
///
/// The rules of [`convert_to_openms_spectrum`] apply: peaks and data arrays are
/// replaced, metadata kept, intensities narrowed to `f32`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the intensity array is not as long as the time
/// array; `chromatogram` is then unchanged.
pub fn convert_to_openms_chromatogram(
    cptr: &Chromatogram,
    chromatogram: &mut MSChromatogram,
) -> Result<()> {
    convert_to_openms_chromatogram_filter(chromatogram, cptr, f64::NEG_INFINITY, f64::INFINITY)
        .map(|_| ())
}

/// Convert the points of an OpenSWATH chromatogram with `rt_min <= RT <=
/// rt_max` into the peaks of an OpenMS chromatogram.
///
/// Rules as [`convert_to_openms_chromatogram`]; points outside the window are
/// skipped. A NaN RT is never outside the window by the source's test
/// (`rt < rt_min || rt > rt_max` is false), so it is kept, as in the source.
/// Returns the number of points kept.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the intensity array is not as long as the time
/// array; `chromatogram` is then unchanged.
pub fn convert_to_openms_chromatogram_filter(
    chromatogram: &mut MSChromatogram,
    cptr: &Chromatogram,
    rt_min: f64,
    rt_max: f64,
) -> Result<usize> {
    let rt = &cptr.time_array().data;
    let intensity = &cptr.intensity_array().data;
    check_pair(rt.len(), intensity.len(), "time")?;
    chromatogram.clear(false);
    chromatogram.peaks.reserve(rt.len());
    for (&t, &i) in rt.iter().zip(intensity) {
        if t < rt_min || t > rt_max {
            continue;
        }
        chromatogram.peaks.push(ChromatogramPeak::new(t, i as f32));
    }
    Ok(chromatogram.peaks.len())
}

/// Convert a light peptide into an [`AASequence`] with its modifications
/// (source `convertPeptideToAASequence`).
///
/// The sequence is parsed, then every modification with a UniMod id other
/// than -1 is set as `UniMod:<id>` at its location (-1 the N-terminus, the
/// sequence length the C-terminus), through
/// [`set_modification`].
///
/// # Errors
///
/// [`Error::InvalidValue`] for a compound that is not a peptide (a debug-only
/// precondition in the source), for a sequence longer than `i32::MAX`, and the
/// errors of [`AASequence::parse`] and of `set_modification` (an unknown
/// modification or a location out of range).
pub fn convert_peptide_to_aa_sequence(peptide: &LightCompound) -> Result<AASequence> {
    if !peptide.is_peptide() {
        return Err(Error::InvalidValue(
            "Function needs peptide, not metabolite".into(),
        ));
    }
    let mut aa_sequence = AASequence::parse(&peptide.sequence)?;
    let max_size = i32::try_from(peptide.sequence.len())
        .map_err(|_| Error::InvalidValue("peptide sequence is too long".into()))?;
    for modification in &peptide.modifications {
        if modification.unimod_id != -1 {
            set_modification(
                modification.location,
                max_size,
                &format!("UniMod:{}", modification.unimod_id),
                &mut aa_sequence,
            )?;
        }
    }
    Ok(aa_sequence)
}

/// The retention time in seconds of the first entry of an analyte, and the
/// range of two entries, as both `convertTargetedCompound` overloads write
/// them into a [`LightCompound`].
fn convert_retention_times(analyte: &PeptideCompound, compound: &mut LightCompound) -> Result<()> {
    // Reset the RT fields up front so a reused compound does not keep a stale
    // value when the source has no RT information.
    compound.rt = f64::NAN;
    compound.rt_start = f64::NAN;
    compound.rt_end = f64::NAN;
    if analyte.has_retention_time() {
        compound.rt = analyte.retention_time()?;
        if analyte.retention_time_unit()? == RTUnit::Minute {
            compound.rt = 60.0 * analyte.retention_time()?;
        }
    }
    // A retention time range is encoded by exactly two entries; three or more
    // mean something else and are not truncated.
    if analyte.rts.len() == 2 {
        let seconds = |r: &RetentionTime| -> f64 {
            let value = r.rt().unwrap_or(0.0);
            if r.retention_time_unit == RTUnit::Minute {
                60.0 * value
            } else {
                value
            }
        };
        compound.rt_start = seconds(&analyte.rts[0]);
        compound.rt_end = seconds(&analyte.rts[1]);
    }
    compound.set_drift_time(analyte.drift_time());
    if analyte.has_charge() {
        compound.charge = analyte.charge_state().unwrap_or(0);
    }
    Ok(())
}

fn unimod_record_id(id: Option<u32>) -> i32 {
    id.and_then(|v| i32::try_from(v).ok()).unwrap_or(-1)
}

/// Convert a heavyweight peptide into a [`LightCompound`] (source
/// `convertTargetedCompound(const TargetedExperiment::Peptide&, LightCompound&)`).
///
/// Id, retention time (minutes are converted to seconds; `NaN` when unset;
/// `rt_start` and `rt_end` from exactly two entries), drift time, charge,
/// sequence, peptide group label and protein references are copied; the
/// `GeneName`, `SumFormula` and `CompoundName` metavalues, when present, become
/// `gene_name`, `sum_formula` and `compound_name` (a compound with a name is
/// not a peptide, so its modifications are not read). For a peptide the
/// modifications of the parsed sequence become [`LightModification`]s: the
/// N-terminal at location `-1`, the C-terminal at the sequence length, then the
/// modified residues, each with its UniMod record id (`-1` when the
/// modification has none).
///
/// # Errors
///
/// [`Error::InvalidValue`] when the sequence of a peptide, or one of its
/// modifications, cannot be resolved (the source's `Exception::InvalidValue`
/// from `AASequence::fromString`).
pub fn convert_targeted_peptide(pep: &Peptide) -> Result<LightCompound> {
    let mut p = LightCompound {
        id: pep.id.clone(),
        ..LightCompound::default()
    };
    convert_retention_times(&pep.base, &mut p)?;
    p.sequence = pep.sequence.clone();
    p.peptide_group_label = pep.peptide_group_label().to_owned();
    let meta = &pep.cv_terms.metadata;
    if meta.contains_key("GeneName") {
        p.gene_name = data_value_to_string(meta.get("GeneName"));
    }
    // Is it potentially a metabolomics compound
    if meta.contains_key("SumFormula") {
        p.sum_formula = data_value_to_string(meta.get("SumFormula"));
    }
    if meta.contains_key("CompoundName") {
        p.compound_name = data_value_to_string(meta.get("CompoundName"));
    }
    p.protein_refs = pep.protein_refs.clone();
    // Mapping of peptide modifications (not for metabolites)
    if p.is_peptide() {
        let aa_sequence = get_aa_sequence(pep)?;
        if let Some(modification) = aa_sequence.n_terminal_modification() {
            p.modifications.push(LightModification {
                location: -1,
                unimod_id: unimod_record_id(modification.record_id()),
            });
        }
        if let Some(modification) = aa_sequence.c_terminal_modification() {
            p.modifications.push(LightModification {
                location: i32::try_from(aa_sequence.len())
                    .map_err(|_| Error::InvalidValue("peptide sequence is too long".into()))?,
                unimod_id: unimod_record_id(modification.record_id()),
            });
        }
        for i in 0..aa_sequence.len() {
            if let Some(modification) = aa_sequence.residue_modification(i)? {
                p.modifications.push(LightModification {
                    location: i32::try_from(i)
                        .map_err(|_| Error::InvalidValue("peptide sequence is too long".into()))?,
                    unimod_id: unimod_record_id(modification.record_id()),
                });
            }
        }
    }
    Ok(p)
}

/// Convert a heavyweight small-molecule compound into a [`LightCompound`]
/// (source `convertTargetedCompound(const TargetedExperiment::Compound&, LightCompound&)`).
///
/// Id, retention time (as [`convert_targeted_peptide`]), drift time, charge,
/// the molecular formula as `sum_formula` and the `CompoundName` metavalue as
/// `compound_name` are copied.
///
/// # Errors
///
/// None in practice; the retention time accessors cannot fail after the
/// `has_retention_time` test.
pub fn convert_targeted_compound(compound: &Compound) -> Result<LightCompound> {
    let mut comp = LightCompound {
        id: compound.id.clone(),
        ..LightCompound::default()
    };
    convert_retention_times(&compound.base, &mut comp)?;
    comp.sum_formula = compound.molecular_formula.clone();
    let meta = &compound.cv_terms.metadata;
    if meta.contains_key("CompoundName") {
        comp.compound_name = data_value_to_string(meta.get("CompoundName"));
    }
    Ok(comp)
}

/// Convert a heavyweight targeted experiment into a [`LightTargetedExperiment`]
/// (source `convertTargetedExp`).
///
/// Proteins keep only their id; peptides then compounds become light
/// compounds; each transition becomes a [`LightTransition`] with its native
/// ID, m/z values, library intensity and peptide reference (the compound
/// reference when the peptide reference is empty), the precursor ion mobility
/// of the referenced peptide, the product charge when set, the three flags and
/// the decoy status. The decoy status is read first from a `decoy` CV term
/// whose first value is `1`, then from the accessions `MS:1002007` (target)
/// and `MS:1002008` (decoy), then from the decoy transition type.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a peptide reference matches no peptide or
/// compound (the source dereferences a null pointer), and the errors of
/// [`convert_targeted_peptide`].
pub fn convert_targeted_exp(
    transition_exp: &TargetedExperiment,
) -> Result<LightTargetedExperiment> {
    let mut light = LightTargetedExperiment::new();
    // copy proteins
    for protein in transition_exp.proteins() {
        light.proteins.push(LightProtein {
            id: protein.id.clone(),
            ..LightProtein::default()
        });
    }
    // copy peptides and store as compounds
    for peptide in transition_exp.peptides() {
        light.compounds.push(convert_targeted_peptide(peptide)?);
    }
    // copy compounds and store as compounds
    for compound in transition_exp.compounds() {
        light.compounds.push(convert_targeted_compound(compound)?);
    }
    // mapping of transitions
    for transition in transition_exp.transitions() {
        let mut t = LightTransition {
            transition_name: transition.native_id().to_owned(),
            product_mz: transition.product_mz(),
            precursor_mz: transition.precursor_mz(),
            library_intensity: transition.library_intensity(),
            peptide_ref: transition.peptide_ref().to_owned(),
            ..LightTransition::default()
        };
        // if the compound is a peptide, the ion mobility comes from the compound
        if !t.peptide_ref.is_empty() {
            let p = light.peptide_by_ref(&t.peptide_ref).ok_or_else(|| {
                Error::InvalidValue(format!(
                    "transition '{}' refers to the peptide '{}', which the experiment does not hold",
                    t.transition_name, t.peptide_ref
                ))
            })?;
            t.precursor_im = p.drift_time();
        } else {
            // try the compound ref
            t.peptide_ref = transition.compound_ref().to_owned();
        }
        if transition.is_product_charge_state_set() {
            t.fragment_charge = transition.product_charge_state().unwrap_or(0) as i8;
        }
        t.set_decoy(false);
        // legacy CV terms first, then the decoy transition type
        let cv_terms = &transition.cv_terms;
        let legacy_decoy = cv_terms
            .get("decoy")
            .and_then(|v| v.first())
            .map(|term| data_value_to_string(Some(&term.value)) == "1");
        if legacy_decoy == Some(true) {
            t.set_decoy(true);
        } else if cv_terms.contains("MS:1002007") {
            // target SRM transition
            t.set_decoy(false);
        } else if cv_terms.contains("MS:1002008") {
            // decoy SRM transition
            t.set_decoy(true);
        } else {
            match transition.decoy_transition_type() {
                DecoyTransitionType::Unknown | DecoyTransitionType::Target => t.set_decoy(false),
                DecoyTransitionType::Decoy => t.set_decoy(true),
            }
        }
        t.set_detecting_transition(transition.is_detecting_transition());
        t.set_identifying_transition(transition.is_identifying_transition());
        t.set_quantifying_transition(transition.is_quantifying_transition());
        light.transitions.push(t);
    }
    Ok(light)
}

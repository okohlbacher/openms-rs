// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Map targeted assays to chromatograms.
//!
//! Ports `ANALYSIS/TARGETED/MRMMapping.h` (with its `.cpp`), the algorithm of
//! the `MRMMapper` TOPP tool. [`MRMMapping::map_experiment`] matches each
//! chromatogram of an experiment against the transitions of a
//! [`TargetedExperiment`] by precursor and product m/z and returns annotated
//! copies of the mapped chromatograms.
//!
//! The source's second overload, taking an OpenSWATH
//! `LightTargetedExperiment`, converts it to a `TargetedExperiment` (peptides
//! with id and sequence, proteins with id and sequence, transitions with
//! native id, peptide reference, m/z, library intensity and flags) and calls
//! the first: [`MRMMapping::map_experiment_light`].
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/TARGETED_EXPERIMENT_SUPPORT.md).
//!
//! The source contains no `#pragma omp`; the port is serial.
//!
//! [`MRMMapping::map_experiment`]: crate::analysis::openswath::mrm_mapping::MRMMapping::map_experiment
//! [`TargetedExperiment`]: crate::analysis::openswath::targeted_experiment::TargetedExperiment

use crate::analysis::openswath::reaction_monitoring_transition::{
    DecoyTransitionType, ReactionMonitoringTransition,
};
use crate::analysis::openswath::targeted_experiment::TargetedExperiment;
use crate::analysis::openswath::targeted_experiment_helper::{Peptide, Protein};
use crate::analysis::openswath::transition_experiment::LightTargetedExperiment;
use crate::kernel::{MSChromatogram, MSExperiment};
use crate::metadata::MetaValue;
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use std::collections::HashMap;

/// The source `DefaultParamHandler` name.
pub const MRM_MAPPING_NAME: &str = "MRMMapping";

/// Default ceiling on the number of chromatogram data points
/// [`MRMMapping::map_experiment`] may copy into its output: 2^30 points.
pub const DEFAULT_MAX_OUTPUT_POINTS: usize = 1 << 30;

/// Maps chromatograms to the assays of a targeted experiment.
///
/// Parameters (source defaults): `precursor_tolerance` 0.1 and
/// `product_tolerance` 0.1 (in Th), `map_multiple_assays` false and
/// `error_on_unmapped` false.
#[derive(Clone, Debug, PartialEq)]
pub struct MRMMapping {
    /// Precursor tolerance when mapping, in Th. A chromatogram maps to a
    /// transition when the absolute precursor m/z difference is strictly
    /// below it.
    pub precursor_tolerance: f64,
    /// Product tolerance when mapping, in Th; strict, as the precursor one.
    pub product_tolerance: f64,
    /// Allow a chromatogram to map to several assays; each mapped assay then
    /// gets its own copy of the chromatogram.
    pub map_multiple_assays: bool,
    /// Treat unmapped chromatograms as an error.
    pub error_on_unmapped: bool,
    /// Native ceiling on the chromatogram data points copied into the output
    /// ([`DEFAULT_MAX_OUTPUT_POINTS`]); not a source parameter. With
    /// `map_multiple_assays`, `n` chromatograms and `m` transitions can
    /// produce `n * m` copies, so the total is checked before anything is
    /// copied.
    pub max_output_points: usize,
}

impl Default for MRMMapping {
    fn default() -> Self {
        Self {
            precursor_tolerance: 0.1,
            product_tolerance: 0.1,
            map_multiple_assays: false,
            error_on_unmapped: false,
            max_output_points: DEFAULT_MAX_OUTPUT_POINTS,
        }
    }
}

fn entry(name: &str, value: ParamValue, description: &str) -> ParamEntry {
    ParamEntry {
        name: name.into(),
        description: description.into(),
        value,
        ..ParamEntry::default()
    }
}

fn flag_entry(name: &str, value: bool, description: &str) -> ParamEntry {
    let mut e = entry(
        name,
        ParamValue::String(if value { "true" } else { "false" }.into()),
        description,
    );
    e.valid_strings = ["true", "false"].iter().map(|s| (*s).to_owned()).collect();
    e
}

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

impl MRMMapping {
    /// A mapper with the source defaults.
    pub fn new() -> Self {
        Self::default()
    }

    /// The source defaults as a [`Param`] (`getDefaults`).
    ///
    /// # Errors
    ///
    /// Propagates a [`Param`] construction failure, which does not occur for
    /// these four fixed entries.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The current parameters in the layout of [`defaults`](Self::defaults)
    /// (`getParameters`). The native
    /// [`max_output_points`](Self::max_output_points) is not a parameter.
    ///
    /// # Errors
    ///
    /// As [`defaults`](Self::defaults).
    pub fn to_param(&self) -> Result<Param> {
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: vec![
                entry(
                    "precursor_tolerance",
                    ParamValue::Float(self.precursor_tolerance),
                    "Precursor tolerance when mapping (in Th)",
                ),
                entry(
                    "product_tolerance",
                    ParamValue::Float(self.product_tolerance),
                    "Product tolerance when mapping (in Th)",
                ),
                flag_entry(
                    "map_multiple_assays",
                    self.map_multiple_assays,
                    "Allow to map multiple assays to chromatograms and duplicate these chromatograms in the output.",
                ),
                flag_entry(
                    "error_on_unmapped",
                    self.error_on_unmapped,
                    "Treat remaining, unmapped chromatograms as an error",
                ),
            ],
            nodes: Vec::new(),
        })
    }

    /// Build a mapper from parameters, as the source `setParameters`
    /// followed by `updateMembers_`. Missing parameters take their defaults;
    /// unknown ones are ignored (see
    /// [`from_param_with_warnings`](Self::from_param_with_warnings)).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its valid strings (the source `Exception::InvalidParameter`),
    /// and, as a native check, when a tolerance is NaN (which the source
    /// accepts and which then maps nothing).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`from_param`](Self::from_param), also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`from_param`](Self::from_param).
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(MRM_MAPPING_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        let tolerance = |key: &str| -> Result<f64> {
            let value = param.value(key)?.to_f64()?;
            if value.is_nan() {
                return Err(bad(format!("parameter '{key}' must not be NaN")));
            }
            Ok(value)
        };
        Ok(Self {
            precursor_tolerance: tolerance("precursor_tolerance")?,
            product_tolerance: tolerance("product_tolerance")?,
            map_multiple_assays: param.value("map_multiple_assays")?.to_bool()?,
            error_on_unmapped: param.value("error_on_unmapped")?.to_bool()?,
            ..Self::default()
        })
    }

    fn matches(&self, chromatogram: &MSChromatogram, targeted: &TargetedExperiment) -> Vec<usize> {
        targeted
            .transitions()
            .iter()
            .enumerate()
            .filter(|(_, t)| {
                (chromatogram.precursor.mz - t.precursor_mz()).abs() < self.precursor_tolerance
                    && (chromatogram.product.mz - t.product_mz()).abs() < self.product_tolerance
            })
            .map(|(j, _)| j)
            .collect()
    }

    /// Map the chromatograms of `input` to the assays of `targeted`.
    ///
    /// The output is a new experiment carrying `input`'s experimental
    /// settings (the source copies the `ExperimentalSettings` base and
    /// clears the rest), no spectra, and one annotated copy of each mapped
    /// chromatogram. For every chromatogram, in order:
    ///
    /// 1. When both its precursor and product m/z are below `1e-5` in
    ///    magnitude it is skipped, unless `map_multiple_assays` is set, in
    ///    which case it is still tried (the source only logs).
    /// 2. It matches every transition whose precursor and product m/z lie
    ///    strictly within the tolerances.
    /// 3. Each match yields a copy of the chromatogram whose native id is the
    ///    transition's native id, and whose precursor carries the metadata
    ///    `description` (the transition's native id) and `peptide_sequence`
    ///    (the transition's peptide reference, replaced by the referenced
    ///    peptide's sequence when the first peptide with that id has a
    ///    non-empty one).
    /// 4. No match counts the chromatogram as unmapped; one match is
    ///    appended; several are all appended with `map_multiple_assays` and
    ///    are an error without it.
    ///
    /// A chromatogram skipped in step 1 does not count as unmapped, as in the
    /// source.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] (the source's `Exception::IllegalArgument`,
    /// with its messages) when a chromatogram maps to several assays without
    /// `map_multiple_assays`, and, with `error_on_unmapped`, at the first
    /// unmapped chromatogram. Also [`Error::InvalidValue`] when the output
    /// would exceed [`max_output_points`](Self::max_output_points); this is
    /// checked before any chromatogram is copied. No partial output is
    /// returned on error (the source leaves its out-parameter half filled).
    pub fn map_experiment(
        &self,
        input: &MSExperiment,
        targeted: &TargetedExperiment,
    ) -> Result<MSExperiment> {
        // First occurrence of each peptide id, as the source's linear
        // search with `break`.
        let mut first_peptide: HashMap<&str, usize> = HashMap::new();
        for (index, peptide) in targeted.peptides().iter().enumerate() {
            first_peptide.entry(peptide.id.as_str()).or_insert(index);
        }

        // Pass 1: decide every chromatogram without copying anything.
        let mut plan: Vec<(usize, Vec<usize>)> = Vec::new();
        let mut points: usize = 0;
        for (i, chromatogram) in input.chromatograms.iter().enumerate() {
            let prec_product_set =
                !(chromatogram.precursor.mz.abs() < 1e-5 && chromatogram.product.mz.abs() < 1e-5);
            if !prec_product_set && !self.map_multiple_assays {
                continue;
            }
            let mapped = self.matches(chromatogram, targeted);
            if mapped.is_empty() {
                if self.error_on_unmapped {
                    return Err(bad(format!(
                        "Did not find a mapping for chromatogram {i}! Maybe try to increase your mapping tolerance."
                    )));
                }
                continue;
            }
            if mapped.len() > 1 && !self.map_multiple_assays {
                return Err(bad(format!(
                    "Chromatogram {} with {} -> {} maps to multiple assays! Either decrease your mapping tolerance or set map_multiple_assays to true.",
                    chromatogram.native_id,
                    crate::param::value::format_float(chromatogram.precursor.mz, true),
                    crate::param::value::format_float(chromatogram.product.mz, true),
                )));
            }
            points = points
                .checked_add(chromatogram.peaks.len().saturating_mul(mapped.len()))
                .filter(|&total| total <= self.max_output_points)
                .ok_or_else(|| {
                    bad(format!(
                        "MRMMapping output would exceed {} chromatogram points",
                        self.max_output_points
                    ))
                })?;
            plan.push((i, mapped));
        }

        // Pass 2: build the output.
        let mut output = MSExperiment {
            settings: input.settings.clone(),
            ..MSExperiment::default()
        };
        for (i, mapped) in plan {
            let chromatogram = &input.chromatograms[i];
            for j in mapped {
                let transition = &targeted.transitions()[j];
                let mut copy = chromatogram.clone();
                let mut sequence = transition.peptide_ref().to_owned();
                if let Some(&index) = first_peptide.get(transition.peptide_ref()) {
                    let peptide = &targeted.peptides()[index];
                    if !peptide.sequence.is_empty() {
                        sequence = peptide.sequence.clone();
                    }
                }
                let metadata = &mut copy.precursor.cv_terms.metadata;
                metadata.insert("peptide_sequence".into(), MetaValue::from(sequence));
                metadata.insert(
                    "description".into(),
                    MetaValue::from(transition.native_id()),
                );
                copy.native_id = transition.native_id().to_owned();
                output.chromatograms.push(copy);
            }
        }
        Ok(output)
    }

    /// Map chromatograms to the assays of a [`LightTargetedExperiment`]
    /// (the source overload `mapExperiment(PeakMap, LightTargetedExperiment,
    /// PeakMap&)`).
    ///
    /// The light experiment is converted to a [`TargetedExperiment`] as the
    /// source does and [`map_experiment`](Self::map_experiment) runs on it:
    /// each compound becomes a peptide carrying only its id and sequence,
    /// each protein a protein with its id and sequence, and each transition a
    /// transition with its native id (the light transition name), peptide
    /// reference, precursor and product m/z, library intensity, the detecting,
    /// quantifying and identifying flags, and the decoy type
    /// ([`DecoyTransitionType::Decoy`] when the light transition is a decoy,
    /// otherwise [`DecoyTransitionType::Target`]). Everything else the light
    /// types carry (retention times, charges, fragment annotations) is
    /// dropped, as in the source, and none of it is read by the mapping.
    ///
    /// # Errors
    ///
    /// As [`map_experiment`](Self::map_experiment).
    pub fn map_experiment_light(
        &self,
        input: &MSExperiment,
        targeted: &LightTargetedExperiment,
    ) -> Result<MSExperiment> {
        let mut converted = TargetedExperiment::new();
        for compound in targeted.compounds() {
            let mut peptide = Peptide::default();
            peptide.base.id = compound.id.clone();
            peptide.sequence = compound.sequence.clone();
            converted.add_peptide(peptide);
        }
        for protein in targeted.proteins() {
            converted.add_protein(Protein {
                id: protein.id.clone(),
                sequence: protein.sequence.clone(),
                ..Protein::default()
            });
        }
        for light in targeted.transitions() {
            let mut transition = ReactionMonitoringTransition::new();
            transition.set_native_id(light.transition_name.clone());
            transition.set_peptide_ref(light.peptide_ref.clone());
            transition.set_precursor_mz(light.precursor_mz);
            transition.set_product_mz(light.product_mz);
            transition.set_library_intensity(light.library_intensity);
            transition.set_detecting_transition(light.flags.detecting());
            transition.set_quantifying_transition(light.flags.quantifying());
            transition.set_identifying_transition(light.flags.identifying());
            transition.set_decoy_transition_type(if light.decoy() {
                DecoyTransitionType::Decoy
            } else {
                DecoyTransitionType::Target
            });
            converted.add_transition(transition);
        }
        self.map_experiment(input, &converted)
    }
}

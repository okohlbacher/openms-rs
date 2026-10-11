// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Target/decoy bookkeeping for experiment-specific decoys in targeted
//! metabolomics.
//!
//! Ports `ANALYSIS/TARGETED/MetaboTargetedTargetDecoy.h` (with its `.cpp`):
//! resolve overlapping target and decoy fragments, and generate the decoys
//! that fragmentation-tree re-rooting could not produce, for the assay
//! libraries that `AssayGeneratorMetaboSirius` writes. The source class is a
//! bag of static functions; here they are free functions over a
//! [`TargetedExperiment`].
//!
//! Targets and decoys are paired through the metadata of the compounds and
//! transitions: `m_ids_id` (the shared identifier), `decoy` (integer `0` for
//! a target compound, `1` for a decoy) and, for generated decoys,
//! `CompoundName` and `annotation`.
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/TARGETED_EXPERIMENT_SUPPORT.md).
//! The source contains no `#pragma omp`; the port is serial.
//!
//! [`TargetedExperiment`]: crate::analysis::openswath::targeted_experiment::TargetedExperiment

use crate::analysis::openswath::reaction_monitoring_transition::{
    DecoyTransitionType, ReactionMonitoringTransition,
};
use crate::analysis::openswath::targeted_experiment::TargetedExperiment;
use crate::analysis::openswath::targeted_experiment_helper::Compound;
use crate::metadata::{MetaValue, MetaValueData};
use crate::{Error, Result};
use std::collections::BTreeMap;

/// The masses of one target and its decoy, paired by a shared identifier
/// (source `MetaboTargetDecoyMassMapping`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MetaboTargetDecoyMassMapping {
    /// Unique identifier shared by target and decoy (the `m_ids_id`
    /// metadata).
    pub identifier: String,
    /// Compound reference of the target transitions.
    pub target_compound_ref: String,
    /// Compound reference of the decoy transitions; empty when there is no
    /// decoy.
    pub decoy_compound_ref: String,
    /// Product m/z of the target transitions, in transition order.
    pub target_product_masses: Vec<f64>,
    /// Product m/z of the decoy transitions, in transition order.
    pub decoy_product_masses: Vec<f64>,
}

/// Unit of the m/z tolerance of
/// [`resolve_overlapping_target_decoy_masses_by_decoy_mass_shift`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MzToleranceUnit {
    /// Parts per million of the target m/z (the source string `"ppm"`).
    Ppm,
    /// Dalton (the source treats every string other than `"ppm"` as Da).
    Da,
}

impl MzToleranceUnit {
    /// Parse the source's unit string: `"ppm"` is ppm, anything else Da, as
    /// the source's `mz_tol_unit == "ppm"` test.
    pub fn from_source_str(unit: &str) -> Self {
        if unit == "ppm" { Self::Ppm } else { Self::Da }
    }
}

/// `DataValue::toString()` for the metadata identifiers: strings verbatim,
/// integers in decimal, floats as `StringUtils::toStr(double)`, empty as
/// the empty string.
fn data_value_string(value: Option<&MetaValue>) -> String {
    match value.map(MetaValue::data) {
        None | Some(MetaValueData::Empty) => String::new(),
        Some(MetaValueData::Float(x)) => crate::param::value::format_float(*x, true),
        Some(_) => value.map(ToString::to_string).unwrap_or_default(),
    }
}

/// Group the transitions by compound reference, keeping transition order
/// within each group (source `constructTransitionsMap_`).
fn transitions_by_compound_ref(
    t_exp: &TargetedExperiment,
) -> BTreeMap<&str, Vec<&ReactionMonitoringTransition>> {
    let mut map: BTreeMap<&str, Vec<&ReactionMonitoringTransition>> = BTreeMap::new();
    for transition in t_exp.transitions() {
        map.entry(transition.compound_ref())
            .or_default()
            .push(transition);
    }
    map
}

/// Build the target/decoy mass mappings of an experiment.
///
/// One mapping is made per **target** compound, i.e. per compound whose
/// `decoy` metadata equals the integer `0` (a missing or differently typed
/// value is not a target, as the source's `DataValue` comparison). Its
/// identifier is the compound's `m_ids_id` rendered as text. The mapping then
/// collects, in transition order, the product m/z of every target transition
/// and every decoy transition whose `m_ids_id` metadata equals that identifier
/// **as a string** value; the compound references are those of the last such
/// target and decoy transitions.
pub fn construct_target_decoy_mass_mapping(
    t_exp: &TargetedExperiment,
) -> Vec<MetaboTargetDecoyMassMapping> {
    let zero = MetaValue::from(0);
    let identifiers: Vec<String> = t_exp
        .compounds()
        .iter()
        .filter(|compound| compound.cv_terms.metadata.get("decoy") == Some(&zero))
        .map(|compound| data_value_string(compound.cv_terms.metadata.get("m_ids_id")))
        .collect();
    identifiers
        .into_iter()
        .map(|identifier| {
            let wanted = MetaValue::from(identifier.as_str());
            let mut mapping = MetaboTargetDecoyMassMapping {
                identifier,
                ..MetaboTargetDecoyMassMapping::default()
            };
            for transition in t_exp.transitions() {
                if transition.cv_terms.metadata.get("m_ids_id") != Some(&wanted) {
                    continue;
                }
                match transition.decoy_transition_type() {
                    DecoyTransitionType::Target => {
                        mapping.target_product_masses.push(transition.product_mz());
                        mapping.target_compound_ref = transition.compound_ref().to_owned();
                    }
                    DecoyTransitionType::Decoy => {
                        mapping.decoy_product_masses.push(transition.product_mz());
                        mapping.decoy_compound_ref = transition.compound_ref().to_owned();
                    }
                    DecoyTransitionType::Unknown => {}
                }
            }
            mapping
        })
        .collect()
}

/// Shift decoy fragments that overlap a target fragment by `mass_to_add`
/// (e.g. CH2).
///
/// For each mapping, a decoy m/z overlaps when it lies within `mz_tol` of
/// any target m/z of the same mapping: `|decoy - target| <= mz_tol` in Da, or
/// `|decoy - target| / target * 1e6 <= mz_tol` in ppm. Then every transition
/// of the experiment whose compound reference is a mapping's decoy reference
/// and whose product m/z equals (exactly) one of that mapping's overlapping
/// decoy m/z gets `mass_to_add` added to its product m/z. When two mappings
/// share a decoy reference, the later mapping's overlap set wins, as the
/// source's map assignment. `mappings` is not modified (the source takes it
/// by non-const reference but only reads it).
pub fn resolve_overlapping_target_decoy_masses_by_decoy_mass_shift(
    t_exp: &mut TargetedExperiment,
    mappings: &[MetaboTargetDecoyMassMapping],
    mass_to_add: f64,
    mz_tol: f64,
    mz_tol_unit: MzToleranceUnit,
) {
    let mut matched_by_ref: BTreeMap<&str, Vec<f64>> = BTreeMap::new();
    for mapping in mappings {
        let mut matched = Vec::new();
        for &decoy_mz in &mapping.decoy_product_masses {
            let overlaps = mapping.target_product_masses.iter().any(|&target_mz| {
                let difference = match mz_tol_unit {
                    MzToleranceUnit::Ppm => (decoy_mz - target_mz).abs() / target_mz * 1e6,
                    MzToleranceUnit::Da => (decoy_mz - target_mz).abs(),
                };
                difference <= mz_tol
            });
            if overlaps {
                matched.push(decoy_mz);
            }
        }
        matched_by_ref.insert(mapping.decoy_compound_ref.as_str(), matched);
    }
    let transitions = t_exp
        .transitions()
        .iter()
        .map(|transition| {
            let shifted = matched_by_ref
                .get(transition.compound_ref())
                // std::set<double>::contains: equality, so -0.0 == 0.0.
                .is_some_and(|set| set.contains(&transition.product_mz()));
            let mut out = transition.clone();
            if shifted {
                out.set_product_mz(transition.product_mz() + mass_to_add);
            }
            out
        })
        .collect();
    t_exp.set_transitions(transitions);
}

/// Replace `_[` by `_decoy_[` everywhere, as the source
/// `std::regex_replace(s, std::regex(R"(_\[)"), "_decoy_[")`.
fn decoy_name(name: &str) -> String {
    name.replace("_[", "_decoy_[")
}

/// Generate a decoy by mass shift for each target whose fragmentation-tree
/// re-rooting produced none, then rebuild the experiment from the mappings.
///
/// First, every mapping whose decoy and target mass counts differ **and**
/// whose decoy reference is empty gets a decoy reference (the target
/// reference with `_[` replaced by `_decoy_[`) and the target masses plus
/// `mass_to_add` appended to its decoy masses.
///
/// Then the compounds and transitions of `t_exp` are **replaced** by, per
/// mapping in order: the target compound and its transitions (when a
/// compound has the target reference as id), and either the decoy compound
/// and its transitions (when one exists) or a generated decoy — a copy of
/// the target compound with id set to the decoy reference, `CompoundName`
/// suffixed `_decoy` and `decoy` set to `1`, and copies of the target's
/// transitions with `_[` replaced in the native id, decoy type, `annotation`
/// set to `NA`, the decoy compound reference and the mapping's decoy masses
/// as product m/z. Everything not reached through a mapping is dropped, as
/// in the source. A generated decoy of a target without transitions has no
/// transitions (the source logs a warning).
///
/// # Errors
///
/// [`Error::InvalidValue`] when a mapping needs a generated decoy but no
/// compound has its target reference as id (the source dereferences the end
/// iterator: undefined behaviour), and when a generated decoy has fewer decoy
/// masses than the target has transitions (the source indexes past the end
/// of `decoy_product_masses`). `t_exp` and `mappings` are unchanged on error.
pub fn generate_missing_decoys_by_mass_shift(
    t_exp: &mut TargetedExperiment,
    mappings: &mut [MetaboTargetDecoyMassMapping],
    mass_to_add: f64,
) -> Result<()> {
    let mut updated: Vec<MetaboTargetDecoyMassMapping> = mappings.to_vec();
    for mapping in &mut updated {
        if mapping.decoy_product_masses.len() != mapping.target_product_masses.len()
            && mapping.decoy_compound_ref.is_empty()
        {
            mapping.decoy_compound_ref = decoy_name(&mapping.target_compound_ref);
            let shifted: Vec<f64> = mapping
                .target_product_masses
                .iter()
                .map(|&d| d + mass_to_add)
                .collect();
            mapping.decoy_product_masses.extend(shifted);
        }
    }

    let by_ref = transitions_by_compound_ref(t_exp);
    let find = |id: &str| t_exp.compounds().iter().find(|c| c.id == id);
    let mut compounds: Vec<Compound> = Vec::new();
    let mut transitions: Vec<ReactionMonitoringTransition> = Vec::new();
    for mapping in &updated {
        let target = find(&mapping.target_compound_ref);
        let decoy = find(&mapping.decoy_compound_ref);
        if let Some(target) = target {
            compounds.push(target.clone());
            if let Some(list) = by_ref.get(target.id.as_str()) {
                transitions.extend(list.iter().map(|t| (*t).clone()));
            }
        }
        if let Some(decoy) = decoy {
            compounds.push(decoy.clone());
            if let Some(list) = by_ref.get(decoy.id.as_str()) {
                transitions.extend(list.iter().map(|t| (*t).clone()));
            }
            continue;
        }
        let target = target.ok_or_else(|| {
            Error::InvalidValue(format!(
                "no target compound '{}' to generate the decoy '{}' from",
                mapping.target_compound_ref, mapping.decoy_compound_ref
            ))
        })?;
        let mut decoy_compound = target.clone();
        let name = data_value_string(target.cv_terms.metadata.get("CompoundName"));
        decoy_compound.cv_terms.metadata.insert(
            "CompoundName".into(),
            MetaValue::from(format!("{name}_decoy")),
        );
        decoy_compound.id = mapping.decoy_compound_ref.clone();
        decoy_compound
            .cv_terms
            .metadata
            .insert("decoy".into(), MetaValue::from(1));
        if let Some(list) = by_ref.get(target.id.as_str()) {
            if mapping.decoy_product_masses.len() < list.len() {
                return Err(Error::InvalidValue(format!(
                    "decoy '{}' has {} masses for {} target transitions",
                    mapping.decoy_compound_ref,
                    mapping.decoy_product_masses.len(),
                    list.len()
                )));
            }
            for (transition, &mz) in list.iter().zip(&mapping.decoy_product_masses) {
                let mut decoy = (*transition).clone();
                decoy.set_native_id(decoy_name(transition.native_id()));
                decoy.set_decoy_transition_type(DecoyTransitionType::Decoy);
                decoy
                    .cv_terms
                    .metadata
                    .insert("annotation".into(), MetaValue::from("NA"));
                decoy.set_product_mz(mz);
                decoy.set_compound_ref(mapping.decoy_compound_ref.clone());
                transitions.push(decoy);
            }
        }
        compounds.push(decoy_compound);
    }
    t_exp.set_compounds(compounds);
    t_exp.set_transitions(transitions);
    mappings.clone_from_slice(&updated);
    Ok(())
}

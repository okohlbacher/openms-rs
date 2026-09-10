// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Hendrik Weisser, OpenMS Rust contributors $

//! Molecule references and observation-match payloads from Core SDK 6bfc0e4.
//! Reference identity and registration keys are distinct from full record equality.

use super::records::{Meter, charge_history_merge, merge_cost};
use super::{
    AdductId, CompoundId, GraphWork, MoleculeType, ObservationId, OligoId, PeptideId,
    ProcessingStepId, ScoredProcessingResult,
};
use crate::chemistry::{AdductInfo, NAFragmentType, PeptideFragmentType};
use crate::identification::PeakAnnotation;
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::BTreeMap;

/// A registered molecule. Variant order follows the source; IDs replace address order.
/// A reference is resolved and checked against its owning graph by graph operations.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IdentifiedMolecule {
    Peptide(PeptideId),
    Compound(CompoundId),
    Oligo(OligoId),
}
impl IdentifiedMolecule {
    pub fn molecule_type(self) -> MoleculeType {
        match self {
            Self::Peptide(_) => MoleculeType::Protein,
            Self::Compound(_) => MoleculeType::Compound,
            Self::Oligo(_) => MoleculeType::RNA,
        }
    }
    pub fn peptide(self) -> Result<PeptideId> {
        match self {
            Self::Peptide(id) => Ok(id),
            _ => Err(invalid("identified molecule is not a peptide")),
        }
    }
    pub fn compound(self) -> Result<CompoundId> {
        match self {
            Self::Compound(id) => Ok(id),
            _ => Err(invalid("identified molecule is not a compound")),
        }
    }
    pub fn oligo(self) -> Result<OligoId> {
        match self {
            Self::Oligo(id) => Ok(id),
            _ => Err(invalid("identified molecule is not an oligonucleotide")),
        }
    }
}
impl From<PeptideId> for IdentifiedMolecule {
    fn from(value: PeptideId) -> Self {
        Self::Peptide(value)
    }
}
impl From<CompoundId> for IdentifiedMolecule {
    fn from(value: CompoundId) -> Self {
        Self::Compound(value)
    }
}
impl From<OligoId> for IdentifiedMolecule {
    fn from(value: OligoId) -> Self {
        Self::Oligo(value)
    }
}

/// Typed formula selection avoids the source's overlapping peptide/RNA enum numbers.
/// Compound formula dispatch ignores this selection and the requested charge.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MoleculeFormulaKind {
    #[default]
    Full,
    Peptide(PeptideFragmentType),
    Rna(NAFragmentType),
}

/// Fragment annotation vectors by optional processing step. None means unassigned.
/// Vector order and duplicates are retained; equal-step merge keeps the first vector.
pub type PeakAnnotationSteps = BTreeMap<Option<ProcessingStepId>, Vec<PeakAnnotation>>;

/// A match between an observation and a registered peptide, compound or oligonucleotide.
/// Registration identity is observation, molecule and optional adduct; charge is payload.
#[derive(Clone, Debug, PartialEq)]
pub struct ObservationMatch {
    pub identified_molecule: IdentifiedMolecule,
    pub observation: ObservationId,
    pub charge: i32,
    pub adduct: Option<AdductId>,
    pub result: ScoredProcessingResult,
    pub peak_annotations: PeakAnnotationSteps,
}
impl ObservationMatch {
    pub fn new(
        identified_molecule: impl Into<IdentifiedMolecule>,
        observation: ObservationId,
    ) -> Self {
        Self {
            identified_molecule: identified_molecule.into(),
            observation,
            charge: 0,
            adduct: None,
            result: ScoredProcessingResult::default(),
            peak_annotations: PeakAnnotationSteps::new(),
        }
    }

    pub(crate) fn key_cmp(&self, other: &Self) -> Ordering {
        (self.observation, self.identified_molecule, self.adduct).cmp(&(
            other.observation,
            other.identified_molecule,
            other.adduct,
        ))
    }

    pub(crate) fn measure(&self, work: &mut GraphWork) -> Result<usize> {
        let mut meter = Meter::new::<Self>(work)?;
        meter.scored(&self.result)?;
        meter.slots::<(Option<ProcessingStepId>, Vec<PeakAnnotation>)>(
            self.peak_annotations.len(),
        )?;
        for annotations in self.peak_annotations.values() {
            meter.slots::<PeakAnnotation>(annotations.len())?;
            for annotation in annotations {
                meter.text(&annotation.annotation)?;
                annotation.validate()?;
            }
        }
        Ok(meter.bytes)
    }

    /// Merge payload atomically, keeping this record's molecule and observation.
    /// Existing nonzero charge must equal the incoming charge, including incoming zero.
    /// Existing Some(adduct) must equal the incoming adduct. Graph registration handles
    /// different optional adducts as distinct keys instead of invoking this merge.
    pub fn merge(&mut self, other: &Self) -> Result<()> {
        self.merge_with_work(other, &mut GraphWork::default())
    }

    pub(crate) fn merge_with_work(&mut self, other: &Self, work: &mut GraphWork) -> Result<()> {
        let bytes = self
            .measure(work)?
            .checked_add(other.measure(work)?)
            .ok_or_else(|| invalid("observation-match merge size overflows"))?;
        merge_cost(
            work,
            bytes,
            self.peak_annotations
                .len()
                .saturating_add(other.peak_annotations.len())
                .saturating_add(self.result.metadata.len())
                .saturating_add(other.result.metadata.len()),
        )?;
        charge_history_merge(&self.result, &other.result, work)?;
        if self.charge != 0 && self.charge != other.charge {
            return Err(invalid("conflicting observation-match charges"));
        }
        if self.adduct.is_some() && self.adduct != other.adduct {
            return Err(invalid("conflicting observation-match adducts"));
        }
        // The complete incoming payload is charged even when first-wins rules discard it.
        work.copy(bytes)?;
        let mut next = self.clone();
        next.result.merge_unchecked(&other.result);
        if next.charge == 0 {
            next.charge = other.charge;
        }
        if next.adduct.is_none() {
            next.adduct = other.adduct;
        }
        for (&step, annotations) in &other.peak_annotations {
            next.peak_annotations
                .entry(step)
                .or_insert_with(|| annotations.clone());
        }
        next.measure(work)?;
        *self = next;
        Ok(())
    }
}

/// Measure the existing immutable chemistry record without allocating a second representation.
pub(super) fn measure_adduct(value: &AdductInfo, work: &mut GraphWork) -> Result<usize> {
    let mut meter = Meter::new::<AdductInfo>(work)?;
    meter.text(value.name())?;
    meter.formula(value.empirical_formula())?;
    Ok(meter.bytes)
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::identification::graph::{GraphLimits, MAX_GRAPH_WORK};

    fn record() -> ObservationMatch {
        ObservationMatch::new(
            PeptideId { owner: 1, slot: 0 },
            ObservationId { owner: 1, slot: 0 },
        )
    }

    #[test]
    fn discarded_annotations_still_consume_work_and_failed_merge_is_atomic() {
        let mut left = record();
        left.peak_annotations.insert(None, Vec::new());
        let mut right = record();
        right.peak_annotations.insert(
            None,
            vec![PeakAnnotation {
                annotation: "x".repeat(10_000),
                ..PeakAnnotation::default()
            }],
        );
        let saved = left.clone();
        let mut work = GraphWork::new(GraphLimits {
            max_work: 1000,
            ..GraphLimits::default()
        });
        assert!(left.merge_with_work(&right, &mut work).is_err());
        assert_eq!(left, saved);
    }

    #[test]
    fn annotation_clone_is_precharged_against_shared_allocation_budget() {
        let mut left = record();
        let mut right = record();
        right
            .peak_annotations
            .insert(None, vec![PeakAnnotation::default()]);
        let saved = left.clone();
        let mut work = GraphWork::new(GraphLimits {
            max_work: MAX_GRAPH_WORK,
            max_bytes: 0,
            ..GraphLimits::default()
        });
        assert!(left.merge_with_work(&right, &mut work).is_err());
        assert_eq!(left, saved);
    }
}

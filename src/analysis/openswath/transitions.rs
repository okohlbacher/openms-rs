// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A nested (protein → peptide → transition) targeted experiment.
//!
//! Ports `OPENSWATHALGO/DATAACCESS/Transitions.h` (with its empty
//! `Transitions.cpp`). Unlike the flat
//! [`LightTargetedExperiment`](super::transition_experiment::LightTargetedExperiment),
//! these records nest: a [`TargetedExperiment`] owns proteins, a [`Protein`]
//! owns peptides, and a [`Peptide`] owns its transitions. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md).
//!
//! These types are not the heavyweight `ANALYSIS/TARGETED/TargetedExperiment`;
//! only the names coincide.

use super::transition_experiment::{LightModification, LightTransition};

/// A peptide with its modifications and transitions.
///
/// The source leaves `rt` and `charge` uninitialised; the port defaults them to 0.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Peptide {
    /// Retention time.
    pub rt: f64,
    /// Precursor charge.
    pub charge: i32,
    /// Amino-acid sequence.
    pub sequence: String,
    /// Peptide id.
    pub id: String,
    /// Modifications.
    pub modifications: Vec<LightModification>,
    /// Transitions of this peptide.
    pub transitions: Vec<LightTransition>,
}

impl Peptide {
    /// The charge state.
    pub fn charge_state(&self) -> i32 {
        self.charge
    }
}

/// A protein with its peptides.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Protein {
    /// Protein id.
    pub id: String,
    /// Protein sequence.
    pub sequence: String,
    /// Peptides of this protein.
    pub peptides: Vec<Peptide>,
}

/// A targeted experiment as a list of proteins.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TargetedExperiment {
    /// The proteins.
    pub proteins: Vec<Protein>,
}

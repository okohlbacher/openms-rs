// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Helpers over a [`LightTargetedExperiment`](crate::analysis::openswath::transition_experiment::LightTargetedExperiment).
//!
//! Ports `OPENSWATHALGO/DATAACCESS/TransitionHelper.h` with
//! `TransitionHelper.cpp` (the source struct `TransitionHelper` holds two static
//! functions, which become free functions here). See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md).

use super::transition_experiment::{LightCompound, LightTargetedExperiment, LightTransition};
use std::collections::BTreeMap;

/// Group the transitions of `lte` by compound reference into `transmap`.
///
/// Each transition is appended, in library order, to the list under its
/// [`LightTransition::peptide_ref`]. Entries already in `transmap` are kept and
/// appended to, as in the source, which inserts into the caller's map without
/// clearing it. The source takes `lte` by non-`const` reference but only reads
/// it (it copies the transition vector first).
pub fn convert(
    lte: &LightTargetedExperiment,
    transmap: &mut BTreeMap<String, Vec<LightTransition>>,
) {
    for transition in &lte.transitions {
        transmap
            .entry(transition.peptide_ref.clone())
            .or_default()
            .push(transition.clone());
    }
}

/// The **first** compound of `lte` whose id is `peptide_ref`, or `None`.
///
/// Source `findPeptide`, which copies the compound into an out-parameter and
/// returns whether one was found. Note that
/// [`LightTargetedExperiment::compound_by_ref`] resolves a duplicated id to the
/// last compound instead; both rules are the source's.
pub fn find_peptide<'a>(
    lte: &'a LightTargetedExperiment,
    peptide_ref: &str,
) -> Option<&'a LightCompound> {
    lte.compounds.iter().find(|c| c.id == peptide_ref)
}

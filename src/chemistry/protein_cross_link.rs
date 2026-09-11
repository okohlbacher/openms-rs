// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Owned source crosslink candidate. Sequence equality is allocation identity.
use super::{AASequence, TermSpecificity};
use crate::{Error, Result};
use std::{
    hash::{Hash, Hasher},
    sync::Arc,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ProteinProteinCrossLinkType {
    Cross,
    Mono,
    Loop,
}
impl ProteinProteinCrossLinkType {
    pub const COUNT: usize = 3;
}

/// Source pointer identity is retained through owned, immutable Arc handles.
/// Linker mass is finite; negative values are retained. Hash values themselves
/// are native and are not portable source pointer/hash encodings.
#[derive(Clone, Debug)]
pub struct ProteinProteinCrossLink {
    pub alpha: Option<Arc<AASequence>>,
    pub beta: Option<Arc<AASequence>>,
    pub cross_link_position: (isize, isize),
    cross_linker_mass: f64,
    pub cross_linker_name: String,
    pub term_spec_alpha: TermSpecificity,
    pub term_spec_beta: TermSpecificity,
    pub precursor_correction: i32,
}
impl Default for ProteinProteinCrossLink {
    fn default() -> Self {
        Self {
            alpha: None,
            beta: None,
            cross_link_position: (0, 0),
            cross_linker_mass: 0.,
            cross_linker_name: String::new(),
            term_spec_alpha: TermSpecificity::Anywhere,
            term_spec_beta: TermSpecificity::Anywhere,
            precursor_correction: 0,
        }
    }
}
impl ProteinProteinCrossLink {
    pub fn new(cross_linker_mass: f64) -> Result<Self> {
        let mut result = Self::default();
        result.set_cross_linker_mass(cross_linker_mass)?;
        Ok(result)
    }
    pub fn cross_linker_mass(&self) -> f64 {
        self.cross_linker_mass
    }
    pub fn set_cross_linker_mass(&mut self, mass: f64) -> Result<()> {
        if !mass.is_finite() {
            return Err(Error::InvalidValue("crosslink mass must be finite".into()));
        }
        self.cross_linker_mass = mass;
        Ok(())
    }
    pub fn get_type(&self) -> ProteinProteinCrossLinkType {
        if self.beta.as_ref().is_some_and(|b| !b.is_empty()) {
            ProteinProteinCrossLinkType::Cross
        } else if self.cross_link_position.1 == -1 {
            ProteinProteinCrossLinkType::Mono
        } else {
            ProteinProteinCrossLinkType::Loop
        }
    }
}
fn identity(value: &Option<Arc<AASequence>>) -> Option<usize> {
    value.as_ref().map(|v| Arc::as_ptr(v) as usize)
}
impl PartialEq for ProteinProteinCrossLink {
    fn eq(&self, other: &Self) -> bool {
        identity(&self.alpha) == identity(&other.alpha)
            && identity(&self.beta) == identity(&other.beta)
            && self.cross_link_position == other.cross_link_position
            && self.cross_linker_mass == other.cross_linker_mass
            && self.cross_linker_name == other.cross_linker_name
            && self.term_spec_alpha == other.term_spec_alpha
            && self.term_spec_beta == other.term_spec_beta
            && self.precursor_correction == other.precursor_correction
    }
}
impl Eq for ProteinProteinCrossLink {}
impl Hash for ProteinProteinCrossLink {
    fn hash<H: Hasher>(&self, state: &mut H) {
        identity(&self.alpha).hash(state);
        identity(&self.beta).hash(state);
        self.cross_link_position.hash(state);
        (if self.cross_linker_mass == 0. {
            0
        } else {
            self.cross_linker_mass.to_bits()
        })
        .hash(state);
        self.cross_linker_name.hash(state);
        (self.term_spec_alpha as u8).hash(state);
        (self.term_spec_beta as u8).hash(state);
        self.precursor_correction.hash(state);
    }
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! An inclusion or exclusion target of a targeted experiment.
//!
//! Ports `ANALYSIS/TARGETED/IncludeExcludeTarget.h` (with its `.cpp`): a
//! TraML `<Target>` from a `<TargetIncludeList>` or `<TargetExcludeList>`.
//! Structurally it is a reduced transition: precursor and product m/z with
//! their CV terms, interpretations, configurations, a prediction, a
//! retention time and a peptide or compound reference.
//!
//! The source header's brief calls this class "a SRM/MRM transition", copied
//! from `ReactionMonitoringTransition`; it stores a target, not a transition.
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/TARGETED_EXPERIMENT_SUPPORT.md)
//! for the member mapping and the evidence. The source contains no
//! `#pragma omp`; the port is serial.

use crate::Result;
use crate::analysis::openswath::targeted_experiment_helper::{
    Configuration, RetentionTime, hash_f64,
};
use crate::metadata::{CVTerm, CVTermList};
use std::hash::{Hash, Hasher};

/// An inclusion or exclusion target (TraML `<Target>`).
///
/// Precursor and product m/z default to `f64::MAX`, the source's
/// `std::numeric_limits<double>::max()`, meaning "not set". The source
/// header also promises a default charge of `numeric_limits<Int>::max()`,
/// but the class has no charge member; nothing is ported for it.
#[derive(Clone, Debug, PartialEq)]
pub struct IncludeExcludeTarget {
    /// CV terms of the target itself (the source `CVTermList` base).
    pub cv_terms: CVTermList,
    name: String,
    precursor_mz: f64,
    precursor_cv_terms: CVTermList,
    product_mz: f64,
    product_cv_terms: CVTermList,
    interpretation_list: Vec<CVTermList>,
    peptide_ref: String,
    compound_ref: String,
    configurations: Vec<Configuration>,
    prediction: CVTermList,
    rts: RetentionTime,
}

impl Default for IncludeExcludeTarget {
    fn default() -> Self {
        Self {
            cv_terms: CVTermList::default(),
            name: String::new(),
            precursor_mz: f64::MAX,
            precursor_cv_terms: CVTermList::default(),
            product_mz: f64::MAX,
            product_cv_terms: CVTermList::default(),
            interpretation_list: Vec::new(),
            peptide_ref: String::new(),
            compound_ref: String::new(),
            configurations: Vec::new(),
            prediction: CVTermList::default(),
            rts: RetentionTime::default(),
        }
    }
}

/// Equality compares every member, as the source `operator==`, except that
/// [`Configuration`] equality here includes its references and validations
/// (the source compares only a configuration's CV terms; see
/// [`Configuration`]). A value holding NaN is unequal to itself.
impl Eq for IncludeExcludeTarget {}

impl Hash for IncludeExcludeTarget {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.cv_terms.hash(state);
        self.name.hash(state);
        hash_f64(self.precursor_mz, state);
        self.precursor_cv_terms.hash(state);
        hash_f64(self.product_mz, state);
        self.product_cv_terms.hash(state);
        self.interpretation_list.hash(state);
        self.peptide_ref.hash(state);
        self.compound_ref.hash(state);
        self.configurations.hash(state);
        self.prediction.hash(state);
        self.rts.hash(state);
    }
}

impl IncludeExcludeTarget {
    /// A target with both m/z values unset (`f64::MAX`) and nothing else.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the target id.
    pub fn set_name(&mut self, name: impl Into<String>) {
        self.name = name.into();
    }

    /// The target id; empty when unset.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Set the reference to the targeted peptide.
    pub fn set_peptide_ref(&mut self, peptide_ref: impl Into<String>) {
        self.peptide_ref = peptide_ref.into();
    }

    /// Reference to the targeted peptide; empty when unset.
    pub fn peptide_ref(&self) -> &str {
        &self.peptide_ref
    }

    /// Set the reference to the targeted compound.
    pub fn set_compound_ref(&mut self, compound_ref: impl Into<String>) {
        self.compound_ref = compound_ref.into();
    }

    /// Reference to the targeted compound; empty when unset.
    pub fn compound_ref(&self) -> &str {
        &self.compound_ref
    }

    /// Set the precursor m/z (the Q1 value).
    pub fn set_precursor_mz(&mut self, mz: f64) {
        self.precursor_mz = mz;
    }

    /// The precursor m/z; `f64::MAX` when unset.
    pub fn precursor_mz(&self) -> f64 {
        self.precursor_mz
    }

    /// Replace the precursor CV terms.
    pub fn set_precursor_cv_term_list(&mut self, list: CVTermList) {
        self.precursor_cv_terms = list;
    }

    /// Add a precursor CV term.
    ///
    /// # Errors
    ///
    /// Whatever [`CVTermList::add`] returns for an invalid term; the list is
    /// then unchanged.
    pub fn add_precursor_cv_term(&mut self, cv_term: CVTerm) -> Result<()> {
        self.precursor_cv_terms.add(cv_term)
    }

    /// The precursor CV terms (always present, possibly empty).
    pub fn precursor_cv_term_list(&self) -> &CVTermList {
        &self.precursor_cv_terms
    }

    /// Set the product m/z (the Q3 value).
    pub fn set_product_mz(&mut self, mz: f64) {
        self.product_mz = mz;
    }

    /// The product m/z; `f64::MAX` when unset.
    pub fn product_mz(&self) -> f64 {
        self.product_mz
    }

    /// Replace the product CV terms.
    pub fn set_product_cv_term_list(&mut self, list: CVTermList) {
        self.product_cv_terms = list;
    }

    /// Add a product CV term.
    ///
    /// # Errors
    ///
    /// Whatever [`CVTermList::add`] returns for an invalid term; the list is
    /// then unchanged.
    pub fn add_product_cv_term(&mut self, cv_term: CVTerm) -> Result<()> {
        self.product_cv_terms.add(cv_term)
    }

    /// The product CV terms (always present, possibly empty).
    pub fn product_cv_term_list(&self) -> &CVTermList {
        &self.product_cv_terms
    }

    /// Replace the interpretations, each a CV term list.
    pub fn set_interpretations(&mut self, interpretations: Vec<CVTermList>) {
        self.interpretation_list = interpretations;
    }

    /// The interpretations.
    pub fn interpretations(&self) -> &[CVTermList] {
        &self.interpretation_list
    }

    /// Append an interpretation.
    pub fn add_interpretation(&mut self, interpretation: CVTermList) {
        self.interpretation_list.push(interpretation);
    }

    /// Replace the instrument configurations.
    pub fn set_configurations(&mut self, configurations: Vec<Configuration>) {
        self.configurations = configurations;
    }

    /// The instrument configurations.
    pub fn configurations(&self) -> &[Configuration] {
        &self.configurations
    }

    /// Append an instrument configuration.
    pub fn add_configuration(&mut self, configuration: Configuration) {
        self.configurations.push(configuration);
    }

    /// Replace the prediction, a CV term list here (the transition class uses
    /// a `Prediction` record instead).
    pub fn set_prediction(&mut self, prediction: CVTermList) {
        self.prediction = prediction;
    }

    /// Add a prediction CV term.
    ///
    /// # Errors
    ///
    /// Whatever [`CVTermList::add`] returns for an invalid term; the list is
    /// then unchanged.
    pub fn add_prediction_term(&mut self, term: CVTerm) -> Result<()> {
        self.prediction.add(term)
    }

    /// The prediction CV terms (always present, possibly empty).
    pub fn prediction(&self) -> &CVTermList {
        &self.prediction
    }

    /// Replace the retention time.
    pub fn set_retention_time(&mut self, rt: RetentionTime) {
        self.rts = rt;
    }

    /// The retention time.
    pub fn retention_time(&self) -> &RetentionTime {
        &self.rts
    }
}

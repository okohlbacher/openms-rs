// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! One SRM/MRM transition of a targeted experiment.
//!
//! Ports `ANALYSIS/MRM/ReactionMonitoringTransition.h` (with its `.cpp`). A
//! [`ReactionMonitoringTransition`] represents a TraML `<Transition>` element
//! completely: its precursor (Q1) m/z and CV terms, its product (Q3) ion and
//! any intermediate products, an optional retention time and prediction, a
//! reference to the peptide or compound it measures, its library intensity,
//! its target/decoy type and the three detecting/identifying/quantifying
//! flags.
//!
//! It implements [`Transition`], so an [`MRMTransitionGroup`] can now be
//! built from the transitions of a real assay library rather than only from
//! [`SimpleTransition`](crate::kernel::mrm::SimpleTransition).
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/TARGETED_EXPERIMENT_SUPPORT.md)
//! for the member mapping and the evidence. The source contains no
//! `#pragma omp`; the port is serial.
//!
//! [`MRMTransitionGroup`]: crate::kernel::mrm::MRMTransitionGroup
//! [`ReactionMonitoringTransition`]: crate::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition
//! [`Transition`]: crate::kernel::mrm::Transition

use crate::analysis::openswath::targeted_experiment_helper::{
    Prediction, RetentionTime, TraMLProduct, hash_f64,
};
use crate::kernel::mrm::Transition;
use crate::metadata::{CVTerm, CVTermList};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

/// The source typedef `ReactionMonitoringTransition::Product`. The other three source typedefs
/// (`Configuration`, `RetentionTime`, `Prediction`) name the helper types
/// directly and need no alias.
pub type Product = TraMLProduct;

/// Whether a transition belongs to a target or a decoy assay.
///
/// The variant order is the source enum order, which [`Ord`] follows; the
/// summary statistics of a targeted experiment count transitions per type in
/// that order.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DecoyTransitionType {
    /// Unknown type; the default.
    #[default]
    Unknown,
    /// Target transition.
    Target,
    /// Decoy transition.
    Decoy,
}

/// A SRM/MRM transition, as a TraML `<Transition>`.
///
/// The default precursor m/z of `0` means "uninitialised". The default
/// library intensity is `-101`, the source's sentinel, and the default flags
/// are detecting, not identifying, quantifying — the same default the TraML
/// reader assumes, which the source warns must not change.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionMonitoringTransition {
    /// CV terms of the transition itself (the source `CVTermList` base).
    pub cv_terms: CVTermList,
    name: String,
    peptide_ref: String,
    compound_ref: String,
    library_intensity: f64,
    decoy_type: DecoyTransitionType,
    precursor_mz: f64,
    precursor_cv_terms: Option<CVTermList>,
    product: TraMLProduct,
    intermediate_products: Vec<TraMLProduct>,
    rts: RetentionTime,
    prediction: Option<Prediction>,
    detecting: bool,
    identifying: bool,
    quantifying: bool,
}

impl Default for ReactionMonitoringTransition {
    fn default() -> Self {
        Self {
            cv_terms: CVTermList::default(),
            name: String::new(),
            peptide_ref: String::new(),
            compound_ref: String::new(),
            library_intensity: -101.0,
            decoy_type: DecoyTransitionType::Unknown,
            precursor_mz: 0.0,
            precursor_cv_terms: None,
            product: TraMLProduct::default(),
            intermediate_products: Vec::new(),
            rts: RetentionTime::default(),
            prediction: None,
            detecting: true,
            identifying: false,
            quantifying: true,
        }
    }
}

/// Equality compares every member, as the source `operator==`: an absent
/// precursor CV term list or prediction equals only another absent one
/// (`Helpers::cmpPtrSafe`). A value holding NaN is unequal to itself.
impl Eq for ReactionMonitoringTransition {}

/// Consistent with equality over every member. The source hash covers a
/// subset (CV-term accessions and counts but not values, for instance); a
/// subset is equally valid, the full set only collides less.
impl Hash for ReactionMonitoringTransition {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.cv_terms.hash(state);
        self.name.hash(state);
        self.peptide_ref.hash(state);
        self.compound_ref.hash(state);
        hash_f64(self.library_intensity, state);
        self.decoy_type.hash(state);
        hash_f64(self.precursor_mz, state);
        self.precursor_cv_terms.hash(state);
        self.product.hash(state);
        self.intermediate_products.hash(state);
        self.rts.hash(state);
        self.prediction.hash(state);
        self.detecting.hash(state);
        self.identifying.hash(state);
        self.quantifying.hash(state);
    }
}

impl ReactionMonitoringTransition {
    /// A default transition: see the type documentation for the defaults.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the transition id. Same member as [`set_native_id`](Self::set_native_id).
    pub fn set_name(&mut self, name: impl Into<String>) {
        self.name = name.into();
    }

    /// The transition id (a required TraML attribute). Same member as
    /// [`native_id`](Self::native_id).
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Set the transition id; the source `setNativeID` writes the same member
    /// as `setName`.
    pub fn set_native_id(&mut self, name: impl Into<String>) {
        self.name = name.into();
    }

    /// The transition id, under its native-id name.
    pub fn native_id(&self) -> &str {
        &self.name
    }

    /// Set the reference to the measured peptide.
    pub fn set_peptide_ref(&mut self, peptide_ref: impl Into<String>) {
        self.peptide_ref = peptide_ref.into();
    }

    /// Reference to the measured peptide; empty when unset.
    pub fn peptide_ref(&self) -> &str {
        &self.peptide_ref
    }

    /// Set the reference to the measured compound.
    pub fn set_compound_ref(&mut self, compound_ref: impl Into<String>) {
        self.compound_ref = compound_ref.into();
    }

    /// Reference to the measured compound; empty when unset.
    pub fn compound_ref(&self) -> &str {
        &self.compound_ref
    }

    /// Set the precursor m/z (the Q1 value).
    pub fn set_precursor_mz(&mut self, mz: f64) {
        self.precursor_mz = mz;
    }

    /// The precursor m/z (the Q1 value); `0` means uninitialised. A TraML
    /// transition has exactly one precursor, which must supply `MS:1000827`
    /// (isolation window target m/z).
    pub fn precursor_mz(&self) -> f64 {
        self.precursor_mz
    }

    /// Whether precursor CV terms exist, i.e. whether
    /// [`precursor_cv_term_list`](Self::precursor_cv_term_list) returns a list.
    pub fn has_precursor_cv_terms(&self) -> bool {
        self.precursor_cv_terms.is_some()
    }

    /// Replace the precursor CV term list. An empty list still counts as
    /// present, as in the source.
    pub fn set_precursor_cv_term_list(&mut self, list: CVTermList) {
        self.precursor_cv_terms = Some(list);
    }

    /// Add a precursor CV term, creating the list on first use.
    ///
    /// # Errors
    ///
    /// Whatever [`CVTermList::add`] returns for an invalid term; the
    /// transition is then unchanged (the list is not created either). The
    /// source validates nothing.
    pub fn add_precursor_cv_term(&mut self, cv_term: CVTerm) -> Result<()> {
        cv_term.validate()?;
        self.precursor_cv_terms
            .get_or_insert_with(CVTermList::default)
            .add(cv_term)
    }

    /// The precursor CV terms, or `None` when there are none.
    ///
    /// The source returns a reference guarded by an `OPENMS_PRECONDITION`
    /// that callers must satisfy with `hasPrecursorCVTerms()`; the absent case
    /// is `None` here.
    pub fn precursor_cv_term_list(&self) -> Option<&CVTermList> {
        self.precursor_cv_terms.as_ref()
    }

    /// Set the product (Q3) m/z.
    pub fn set_product_mz(&mut self, mz: f64) {
        self.product.set_mz(mz);
    }

    /// The product (Q3) m/z; `0` by default.
    pub fn product_mz(&self) -> f64 {
        self.product.mz()
    }

    /// The product charge, or `None` when unset (the source's unchecked
    /// precondition returns `0`).
    pub fn product_charge_state(&self) -> Option<i32> {
        self.product.charge_state()
    }

    /// Whether the product charge is set.
    pub fn is_product_charge_state_set(&self) -> bool {
        self.product.has_charge()
    }

    /// Add a CV term to the product.
    ///
    /// # Errors
    ///
    /// Whatever [`CVTermList::add`] returns for an invalid term; the product
    /// is then unchanged.
    pub fn add_product_cv_term(&mut self, cv_term: CVTerm) -> Result<()> {
        self.product.cv_terms.add(cv_term)
    }

    /// Intermediate product ions (MS3 and above); empty by default.
    pub fn intermediate_products(&self) -> &[TraMLProduct] {
        &self.intermediate_products
    }

    /// Append an intermediate product ion.
    pub fn add_intermediate_product(&mut self, product: TraMLProduct) {
        self.intermediate_products.push(product);
    }

    /// Replace the intermediate product ions.
    pub fn set_intermediate_products(&mut self, products: Vec<TraMLProduct>) {
        self.intermediate_products = products;
    }

    /// Replace the product (Q3) ion.
    pub fn set_product(&mut self, product: TraMLProduct) {
        self.product = product;
    }

    /// The product (Q3) ion.
    pub fn product(&self) -> &TraMLProduct {
        &self.product
    }

    /// Whether a prediction exists, i.e. whether
    /// [`prediction`](Self::prediction) succeeds.
    pub fn has_prediction(&self) -> bool {
        self.prediction.is_some()
    }

    /// Replace the prediction.
    pub fn set_prediction(&mut self, prediction: Prediction) {
        self.prediction = Some(prediction);
    }

    /// Add a CV term to the prediction, creating an empty prediction first
    /// when there is none.
    ///
    /// # Errors
    ///
    /// Whatever [`CVTermList::add`] returns for an invalid term; the
    /// transition is then unchanged (no prediction is created).
    pub fn add_prediction_term(&mut self, term: CVTerm) -> Result<()> {
        term.validate()?;
        self.prediction
            .get_or_insert_with(Prediction::default)
            .cv_terms
            .add(term)
    }

    /// The prediction.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when there is none, as the source's
    /// `Exception::MissingInformation`, which it throws in every build
    /// configuration.
    pub fn prediction(&self) -> Result<&Prediction> {
        self.prediction.as_ref().ok_or_else(|| {
            Error::MissingInformation(
                "ReactionMonitoringTransition has no Prediction object; check with hasPrediction() before calling getPrediction().".into(),
            )
        })
    }

    /// The target/decoy type of the transition.
    pub fn decoy_transition_type(&self) -> DecoyTransitionType {
        self.decoy_type
    }

    /// Set the target/decoy type of the transition.
    pub fn set_decoy_transition_type(&mut self, decoy_type: DecoyTransitionType) {
        self.decoy_type = decoy_type;
    }

    /// Library intensity: the ion count or normalised ion count from a
    /// spectral library (TraML stores it in CV term `MS:1001226`). `-101` when
    /// never set; negative values are not rejected.
    pub fn library_intensity(&self) -> f64 {
        self.library_intensity
    }

    /// Set the library intensity.
    pub fn set_library_intensity(&mut self, intensity: f64) {
        self.library_intensity = intensity;
    }

    /// Set a retention time on the level of this individual transition.
    ///
    /// This is very likely not what you want: retention time is normally set
    /// on the analyte (see
    /// [`PeptideCompound`](crate::analysis::openswath::targeted_experiment_helper::PeptideCompound)),
    /// and a transition-level value is usually ignored. Use it only with good
    /// reason.
    pub fn set_retention_time(&mut self, rt: RetentionTime) {
        self.rts = rt;
    }

    /// The transition-level retention time; see
    /// [`set_retention_time`](Self::set_retention_time) for why this is
    /// usually not the value you want.
    pub fn retention_time(&self) -> &RetentionTime {
        &self.rts
    }

    /// Whether the transition is a detecting transition: one of the set of an
    /// assay's transitions used for candidate peak-group detection. Ideally
    /// they were observed and validated in earlier experiments and carry a
    /// library intensity (Schubert et al., *Nat Protoc* 2015, 10(3):426-41,
    /// doi:10.1038/nprot.2015.015). `true` by default.
    pub fn is_detecting_transition(&self) -> bool {
        self.detecting
    }

    /// Set the detecting flag.
    pub fn set_detecting_transition(&mut self, val: bool) {
        self.detecting = val;
    }

    /// Whether the transition is an identifying transition: used to identify
    /// a candidate peak group independently, usually against another isomer
    /// or peptidoform (for example site-determining transitions in a
    /// phosphoproteomic DIA experiment), and scored separately from the
    /// detecting transitions. `false` by default.
    pub fn is_identifying_transition(&self) -> bool {
        self.identifying
    }

    /// Set the identifying flag.
    pub fn set_identifying_transition(&mut self, val: bool) {
        self.identifying = val;
    }

    /// Whether the transition is a quantifying transition: used to quantify
    /// the peptide. Excluding one (for example an interfered transition of a
    /// light/heavy pair isolated in the same SWATH window) keeps it out of the
    /// quantification. `true` by default.
    pub fn is_quantifying_transition(&self) -> bool {
        self.quantifying
    }

    /// Set the quantifying flag.
    pub fn set_quantifying_transition(&mut self, val: bool) {
        self.quantifying = val;
    }

    /// Order two transitions by product m/z, the source `ProductMZLess`.
    ///
    /// The source compares with `<`; this uses [`f64::total_cmp`] so that the
    /// order is total. The two agree on every pair of non-NaN values except
    /// `-0.0`, which sorts before `0.0` here and ties with it in the source.
    /// NaN sorts after every number (before them when negative).
    pub fn cmp_product_mz(left: &Self, right: &Self) -> Ordering {
        left.product_mz().total_cmp(&right.product_mz())
    }

    /// Order two transitions by name, the source `NameLess` (byte-wise, as
    /// `std::string` `<`).
    pub fn cmp_name(left: &Self, right: &Self) -> Ordering {
        left.name.as_bytes().cmp(right.name.as_bytes())
    }
}

impl Transition for ReactionMonitoringTransition {
    fn native_id(&self) -> &str {
        &self.name
    }
    fn precursor_mz(&self) -> f64 {
        self.precursor_mz
    }
    fn product_mz(&self) -> f64 {
        self.product.mz()
    }
    fn library_intensity(&self) -> f64 {
        self.library_intensity
    }
    fn is_detecting(&self) -> bool {
        self.detecting
    }
    fn is_identifying(&self) -> bool {
        self.identifying
    }
    fn is_quantifying(&self) -> bool {
        self.quantifying
    }
}

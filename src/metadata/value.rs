// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use crate::param::{MAX_PARAM_BYTES, ParamBudget, ParamValue, ParameterMetaSink};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::mem::size_of;

/// Full controlled-vocabulary unit identity. No ontology lookup is performed.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Unit {
    accession: String,
    name: String,
    cv_ref: String,
}

impl Unit {
    pub fn new(
        accession: impl Into<String>,
        name: impl Into<String>,
        cv_ref: impl Into<String>,
    ) -> Result<Self> {
        let unit = Self {
            accession: accession.into(),
            name: name.into(),
            cv_ref: cv_ref.into(),
        };
        if unit.accession.is_empty() || unit.accession.chars().any(char::is_whitespace) {
            return Err(invalid(
                "unit accession must be nonempty and contain no whitespace",
            ));
        }
        Ok(unit)
    }
    /// Maps a DataValue UO/MS numeric accession to an explicit unit identity.
    pub fn from_ontology(ontology: UnitOntology, accession: u32) -> Self {
        let prefix = match ontology {
            UnitOntology::Unit => "UO",
            UnitOntology::MassSpectrometry => "MS",
        };
        Self {
            accession: format!("{prefix}:{accession:07}"),
            name: String::new(),
            cv_ref: prefix.into(),
        }
    }
    pub fn accession(&self) -> &str {
        &self.accession
    }
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn cv_ref(&self) -> &str {
        &self.cv_ref
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UnitOntology {
    Unit,
    MassSpectrometry,
}

/// Source DataValue alternatives. Floating alternatives must be finite when
/// passed to [`MetaValue::new`]; stored values cannot subsequently be mutated.
/// A crate port that reproduces a source storing a non-finite value creates
/// it without that check, and [`MetaValue::validate`] reports it.
#[derive(Clone, Debug, Default, PartialEq)]
pub enum MetaValueData {
    #[default]
    Empty,
    String(String),
    Integer(i64),
    Float(f64),
    StringList(Vec<String>),
    IntegerList(Vec<i64>),
    FloatList(Vec<f64>),
}

impl Hash for MetaValueData {
    fn hash<H: Hasher>(&self, state: &mut H) {
        std::mem::discriminant(self).hash(state);
        match self {
            Self::Empty => {}
            Self::String(value) => value.hash(state),
            Self::Integer(value) => value.hash(state),
            Self::Float(value) => hash_float(*value, state),
            Self::StringList(values) => values.hash(state),
            Self::IntegerList(values) => values.hash(state),
            Self::FloatList(values) => {
                values.len().hash(state);
                for &value in values {
                    hash_float(value, state);
                }
            }
        }
    }
}

// Equality identifies both signs of zero. Other IEEE representations, including
// NaN payloads in unvalidated public records, are forwarded without alteration.
pub(super) fn hash_float<H: Hasher>(value: f64, state: &mut H) {
    (if value == 0.0 { 0 } else { value.to_bits() }).hash(state);
}

/// Typed owned value and optional unit. Equality is exact and includes units;
/// unlike C++ DataValue scalar epsilon equality it is transitive.
#[derive(Clone, Debug, Default, PartialEq, Hash)]
pub struct MetaValue {
    data: MetaValueData,
    unit: Option<Unit>,
}

impl MetaValue {
    pub fn new(data: MetaValueData) -> Result<Self> {
        let value = Self { data, unit: None };
        value.validate()?;
        Ok(value)
    }
    /// A floating scalar stored as a ported source algorithm stores it,
    /// without the finite check of [`MetaValue::new`] and `TryFrom<f64>`.
    ///
    /// Only for crate ports whose executed source stores a non-finite value
    /// in a `DataValue` and whose callers can observe it
    /// (`FeatureFinderAlgorithmPicked`'s `FWHM`, `score_fit`,
    /// `score_correlation` and `EGH_*` meta values). [`MetaValue::validate`]
    /// still reports such a value, so every checked consumer (writers,
    /// merges) refuses it as before.
    pub(crate) fn source_float(value: f64) -> Self {
        Self {
            data: MetaValueData::Float(value),
            unit: None,
        }
    }
    /// A floating list stored as a ported source reader stores it, without the
    /// finite check of [`MetaValue::new`] and `TryFrom<Vec<f64>>`.
    ///
    /// The list counterpart of [`MetaValue::source_float`], for the featureXML
    /// reader: `writeUserParam_` writes a `DOUBLE_LIST` element-wise through
    /// the same conversion as a scalar (`ListUtilsIO.h:29-44`), so a stored
    /// document can spell `[inf, -inf, NaN]`, and `ListUtils::create<double>`
    /// reads it back. [`MetaValue::validate`] still reports such a value.
    #[cfg(any(feature = "idxml", feature = "featurexml", feature = "consensusxml"))]
    pub(crate) fn source_float_list(values: Vec<f64>) -> Self {
        Self {
            data: MetaValueData::FloatList(values),
            unit: None,
        }
    }
    pub fn data(&self) -> &MetaValueData {
        &self.data
    }
    pub fn unit(&self) -> Option<&Unit> {
        self.unit.as_ref()
    }
    pub fn with_unit(mut self, unit: Unit) -> Result<Self> {
        self.unit = Some(unit);
        self.validate()?;
        Ok(self)
    }
    pub fn without_unit(mut self) -> Self {
        self.unit = None;
        self
    }
    /// Empty is distinct from a present empty string or empty list.
    pub fn is_empty(&self) -> bool {
        matches!(self.data, MetaValueData::Empty)
    }
    pub fn validate(&self) -> Result<()> {
        let valid = match &self.data {
            MetaValueData::Float(value) => value.is_finite(),
            MetaValueData::FloatList(values) => values.iter().all(|value| value.is_finite()),
            _ => true,
        };
        if valid {
            Ok(())
        } else {
            Err(invalid("metadata floating values must be finite"))
        }
    }
    pub fn as_str(&self) -> Result<&str> {
        match &self.data {
            MetaValueData::String(value) => Ok(value),
            _ => Err(invalid("metadata value is not a string")),
        }
    }
    pub fn as_i64(&self) -> Result<i64> {
        match self.data {
            MetaValueData::Integer(value) => Ok(value),
            _ => Err(invalid("metadata value is not an integer")),
        }
    }
    /// Accepts integers as well as floats, matching the numeric source conversion.
    /// Large i64 values can lose precision when represented as f64.
    pub fn as_f64(&self) -> Result<f64> {
        match self.data {
            MetaValueData::Integer(value) => Ok(value as f64),
            MetaValueData::Float(value) => Ok(value),
            _ => Err(invalid("metadata value is not numeric")),
        }
    }
    pub fn as_string_list(&self) -> Result<&[String]> {
        match &self.data {
            MetaValueData::StringList(value) => Ok(value),
            _ => Err(invalid("metadata value is not a string list")),
        }
    }
    pub fn as_integer_list(&self) -> Result<&[i64]> {
        match &self.data {
            MetaValueData::IntegerList(value) => Ok(value),
            _ => Err(invalid("metadata value is not an integer list")),
        }
    }
    pub fn as_float_list(&self) -> Result<&[f64]> {
        match &self.data {
            MetaValueData::FloatList(value) => Ok(value),
            _ => Err(invalid("metadata value is not a float list")),
        }
    }
    /// C++ DataValue's boolean convention is the exact strings true and false.
    pub fn to_bool(&self) -> Result<bool> {
        match self.as_str()? {
            "true" => Ok(true),
            "false" => Ok(false),
            _ => Err(invalid("boolean metadata must be the string true or false")),
        }
    }
}

impl From<String> for MetaValue {
    fn from(value: String) -> Self {
        Self {
            data: MetaValueData::String(value),
            unit: None,
        }
    }
}
impl From<&str> for MetaValue {
    fn from(value: &str) -> Self {
        value.to_owned().into()
    }
}
impl From<i64> for MetaValue {
    fn from(value: i64) -> Self {
        Self {
            data: MetaValueData::Integer(value),
            unit: None,
        }
    }
}
impl From<i32> for MetaValue {
    fn from(value: i32) -> Self {
        i64::from(value).into()
    }
}
impl From<u32> for MetaValue {
    fn from(value: u32) -> Self {
        i64::from(value).into()
    }
}
impl TryFrom<u64> for MetaValue {
    type Error = Error;
    fn try_from(value: u64) -> Result<Self> {
        Ok(i64::try_from(value)
            .map_err(|_| invalid("unsigned metadata integer exceeds i64"))?
            .into())
    }
}
impl TryFrom<f64> for MetaValue {
    type Error = Error;
    fn try_from(value: f64) -> Result<Self> {
        Self::new(MetaValueData::Float(value))
    }
}
impl TryFrom<f32> for MetaValue {
    type Error = Error;
    fn try_from(value: f32) -> Result<Self> {
        Self::try_from(f64::from(value))
    }
}
impl From<Vec<String>> for MetaValue {
    fn from(value: Vec<String>) -> Self {
        Self {
            data: MetaValueData::StringList(value),
            unit: None,
        }
    }
}
impl From<Vec<i64>> for MetaValue {
    fn from(value: Vec<i64>) -> Self {
        Self {
            data: MetaValueData::IntegerList(value),
            unit: None,
        }
    }
}
impl TryFrom<Vec<f64>> for MetaValue {
    type Error = Error;
    fn try_from(value: Vec<f64>) -> Result<Self> {
        Self::new(MetaValueData::FloatList(value))
    }
}

/// Lenient stringification, with source-style list delimiters and Rust's
/// round-trippable scalar formatting. Units are deliberately not appended.
impl fmt::Display for MetaValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fn list<T: fmt::Display>(values: &[T], f: &mut fmt::Formatter<'_>) -> fmt::Result {
            f.write_str("[")?;
            for (i, value) in values.iter().enumerate() {
                if i != 0 {
                    f.write_str(", ")?;
                }
                write!(f, "{value}")?;
            }
            f.write_str("]")
        }
        match &self.data {
            MetaValueData::Empty => Ok(()),
            MetaValueData::String(value) => f.write_str(value),
            MetaValueData::Integer(value) => write!(f, "{value}"),
            MetaValueData::Float(value) => write!(f, "{value}"),
            MetaValueData::StringList(values) => list(values, f),
            MetaValueData::IntegerList(values) => list(values, f),
            MetaValueData::FloatList(values) => list(values, f),
        }
    }
}

pub type MetaInfo = BTreeMap<String, MetaValue>;

pub fn validate_meta(metadata: &MetaInfo) -> Result<()> {
    for value in metadata.values() {
        value.validate()?;
    }
    Ok(())
}

/// Receive a parameter tree as metadata, for
/// `DefaultParamHandler::write_parameters_to_meta_values`.
///
/// The conversion is written here rather than in `param` because `param` is the
/// lower module: `metadata` reaches `chemistry`, which reaches `param`, so a
/// `param` that named `MetaInfo` would close a module cycle. `param` states
/// what it needs as [`ParameterMetaSink`] and this supplies it, which leaves
/// the public path on `DefaultParamHandler` where the source class puts it.
impl ParameterMetaSink for MetaInfo {
    type Value = MetaValue;

    fn measure_existing(&self, budget: &mut ParamBudget<'_>) -> Result<usize> {
        budget.consume(self.len())?;
        let mut total = ParamBudget::mul(self.len(), 128)?;
        for (key, value) in self {
            let mut bytes = ParamBudget::add(key.len(), size_of::<MetaValue>())?;
            budget.consume(bytes)?;
            match value.data() {
                MetaValueData::String(s) => {
                    budget.consume(s.len())?;
                    bytes = ParamBudget::add(bytes, s.len())?;
                }
                MetaValueData::StringList(values) => {
                    let slots = ParamBudget::mul(values.len(), size_of::<String>())?;
                    budget.consume(slots)?;
                    bytes = ParamBudget::add(bytes, slots)?;
                    for s in values {
                        budget.consume(s.len())?;
                        bytes = ParamBudget::add(bytes, s.len())?;
                    }
                }
                MetaValueData::IntegerList(values) => {
                    let slots = ParamBudget::mul(values.len(), size_of::<i64>())?;
                    budget.consume(slots)?;
                    bytes = ParamBudget::add(bytes, slots)?;
                }
                MetaValueData::FloatList(values) => {
                    let slots = ParamBudget::mul(values.len(), size_of::<f64>())?;
                    budget.consume(slots)?;
                    bytes = ParamBudget::add(bytes, slots)?;
                }
                _ => {}
            }
            if let Some(unit) = value.unit() {
                for s in [unit.accession(), unit.name(), unit.cv_ref()] {
                    budget.consume(s.len())?;
                    bytes = ParamBudget::add(bytes, s.len())?;
                }
            }
            total = ParamBudget::add(total, bytes)?;
            if total > MAX_PARAM_BYTES {
                return Err(Error::InvalidValue(
                    "metadata payload limit exceeded".into(),
                ));
            }
        }
        Ok(total)
    }

    fn value_of(value: &ParamValue, budget: &mut ParamBudget<'_>) -> Result<MetaValue> {
        let data = match value {
            ParamValue::Empty => MetaValueData::Empty,
            ParamValue::String(value) => MetaValueData::String(value.clone()),
            ParamValue::Integer(value) => MetaValueData::Integer(*value),
            ParamValue::Float(value) => MetaValueData::Float(*value),
            ParamValue::StringList(value) => MetaValueData::StringList(value.clone()),
            ParamValue::IntegerList(value) => {
                budget.copy(ParamBudget::mul(value.len(), size_of::<i64>())?)?;
                MetaValueData::IntegerList(value.iter().copied().map(i64::from).collect())
            }
            ParamValue::FloatList(value) => MetaValueData::FloatList(value.clone()),
        };
        MetaValue::new(data)
    }

    fn stage(&mut self, key: String, value: MetaValue) {
        self.insert(key, value);
    }

    fn staged(&self) -> usize {
        self.len()
    }

    fn absorb(&mut self, staged: &mut Self) {
        self.append(staged);
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MetaMergePolicy {
    KeepExisting,
    Overwrite,
    RejectConflicts,
}

/// Merges metadata transactionally. RejectConflicts permits equal shared values.
pub fn merge_meta(target: &mut MetaInfo, source: &MetaInfo, policy: MetaMergePolicy) -> Result<()> {
    validate_meta(target)?;
    validate_meta(source)?;
    if policy == MetaMergePolicy::RejectConflicts
        && source
            .iter()
            .any(|(key, value)| target.get(key).is_some_and(|existing| existing != value))
    {
        return Err(invalid("conflicting metadata values"));
    }
    for (key, value) in source {
        if policy == MetaMergePolicy::Overwrite {
            target.insert(key.clone(), value.clone());
        } else {
            target.entry(key.clone()).or_insert_with(|| value.clone());
        }
    }
    Ok(())
}

/// Preserves legacy strings without guessing numbers or units.
pub fn meta_from_strings(values: &BTreeMap<String, String>) -> MetaInfo {
    values
        .iter()
        .map(|(key, value)| (key.clone(), value.clone().into()))
        .collect()
}

/// Only strings with no units can be exported without losing type information.
pub fn meta_to_strings(values: &MetaInfo) -> Result<BTreeMap<String, String>> {
    values
        .iter()
        .map(|(key, value)| {
            if value.unit().is_some() {
                return Err(invalid("string metadata export would lose a unit"));
            }
            Ok((key.clone(), value.as_str()?.to_owned()))
        })
        .collect()
}

/// Explicitly discards units and type information, using [`MetaValue`]'s Display.
pub fn meta_to_strings_lossy(values: &MetaInfo) -> BTreeMap<String, String> {
    values
        .iter()
        .map(|(key, value)| (key.clone(), value.to_string()))
        .collect()
}

/// Controlled vocabulary term. A value's optional unit is the term's unit.
#[derive(Clone, Debug, Default, PartialEq, Hash)]
pub struct CVTerm {
    pub accession: String,
    pub name: String,
    pub cv_ref: String,
    pub value: MetaValue,
}

impl CVTerm {
    pub fn new(
        accession: impl Into<String>,
        name: impl Into<String>,
        cv_ref: impl Into<String>,
    ) -> Self {
        Self {
            accession: accession.into(),
            name: name.into(),
            cv_ref: cv_ref.into(),
            value: MetaValue::default(),
        }
    }
    pub fn has_value(&self) -> bool {
        !self.value.is_empty()
    }
    pub fn has_unit(&self) -> bool {
        self.value.unit().is_some()
    }
    pub fn validate(&self) -> Result<()> {
        if self.accession.is_empty() || self.accession.chars().any(char::is_whitespace) {
            return Err(invalid(
                "CV term accession must be nonempty and contain no whitespace",
            ));
        }
        self.value.validate()
    }
}

/// Accession-indexed terms. Duplicates within each accession are retained in order.
#[derive(Clone, Debug, Default, PartialEq, Hash)]
pub struct CVTermList {
    terms: BTreeMap<String, Vec<CVTerm>>,
    pub metadata: MetaInfo,
}

impl CVTermList {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn terms(&self) -> &BTreeMap<String, Vec<CVTerm>> {
        &self.terms
    }
    pub fn get(&self, accession: &str) -> Option<&[CVTerm]> {
        self.terms.get(accession).map(Vec::as_slice)
    }
    pub fn contains(&self, accession: &str) -> bool {
        self.terms.contains_key(accession)
    }
    /// Refers to CV terms only; ordinary metadata does not make this list nonempty.
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
    pub fn add(&mut self, term: CVTerm) -> Result<()> {
        term.validate()?;
        self.terms
            .entry(term.accession.clone())
            .or_default()
            .push(term);
        Ok(())
    }
    /// Appends terms, matching the surprisingly additive C++ setCVTerms method.
    pub fn add_terms(&mut self, terms: &[CVTerm]) -> Result<()> {
        for term in terms {
            term.validate()?;
        }
        for term in terms {
            self.terms
                .entry(term.accession.clone())
                .or_default()
                .push(term.clone());
        }
        Ok(())
    }
    pub fn replace(&mut self, term: CVTerm) -> Result<()> {
        term.validate()?;
        self.terms.insert(term.accession.clone(), vec![term]);
        Ok(())
    }
    /// Replaces one accession; an empty vector retains an explicitly present key.
    pub fn replace_accession(&mut self, accession: &str, terms: Vec<CVTerm>) -> Result<()> {
        if accession.is_empty() || accession.chars().any(char::is_whitespace) {
            return Err(invalid("invalid CV accession"));
        }
        for term in &terms {
            term.validate()?;
            if term.accession != accession {
                return Err(invalid("CV term accession does not match map key"));
            }
        }
        self.terms.insert(accession.into(), terms);
        Ok(())
    }
    /// Appends all incoming terms without duplicate checking; metadata is retained.
    pub fn consume(&mut self, other: &Self) -> Result<()> {
        other.validate()?;
        for (key, terms) in &other.terms {
            self.terms
                .entry(key.clone())
                .or_default()
                .extend(terms.iter().cloned());
        }
        Ok(())
    }
    /// Replaces the complete CV map while retaining ordinary metadata.
    pub fn replace_all(&mut self, terms: BTreeMap<String, Vec<CVTerm>>) -> Result<()> {
        let mut replacement = Self::new();
        for (accession, terms) in terms {
            replacement.replace_accession(&accession, terms)?;
        }
        self.terms = replacement.terms;
        Ok(())
    }
    pub fn remove(&mut self, accession: &str) -> Option<Vec<CVTerm>> {
        self.terms.remove(accession)
    }
    pub fn validate(&self) -> Result<()> {
        validate_meta(&self.metadata)?;
        for terms in self.terms.values() {
            for term in terms {
                term.validate()?;
            }
        }
        Ok(())
    }
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

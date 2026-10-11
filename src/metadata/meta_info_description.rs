// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The description of a meta-data array: a name, meta values and data
//! processing, with the source's ordering.
//!
//! Ports `src/openms/include/OpenMS/METADATA/MetaInfoDescription.h` and
//! `src/openms/source/METADATA/MetaInfoDescription.cpp` (core SDK `bc9cc12`).
//! See `docs/META_INFO_DESCRIPTION_SUPPORT.md`.
//!
//! The source class is the base of the three data arrays of
//! `METADATA/DataArrays.h`; this crate's [`DataArray`](crate::kernel::DataArray)
//! carries the same three fields inline and reaches the ordering through
//! [`DescriptionView`](crate::metadata::meta_info_description::DescriptionView). [`MetaInfoDescription`](crate::metadata::meta_info_description::MetaInfoDescription) is the standalone value.
//!
//! # Ordering
//!
//! The source's `operator<` is not the order its `operator==` implies, and
//! this port keeps it exactly ([`MetaInfoDescription::source_lt`](crate::metadata::meta_info_description::MetaInfoDescription::source_lt)):
//!
//! 1. fewer data-processing entries sort first;
//! 2. only when the meta values differ: a description **with** meta values
//!    sorts before one without; then the names decide; then the meta-value
//!    keys as registry indices, lexicographically; then each value's
//!    `DataValue::toString()` text in key order;
//! 3. otherwise, the number of processing actions of each data-processing
//!    entry in turn;
//! 4. otherwise neither is less.
//!
//! So two descriptions that differ **only in their names** are unequal and
//! yet neither is less than the other, because the name is consulted only
//! inside step 2 (recorded as a source defect candidate). The comparison needs
//! the registry because the source orders keys by their registry index.
//!
//! # Example
//!
//! ```
//! use openms::metadata::meta_info_description::MetaInfoDescription;
//! use openms::metadata::meta_info_indexed::IndexedMetaInfo;
//! use openms::metadata::meta_info_registry::MetaInfoRegistry;
//!
//! let mut registry = MetaInfoRegistry::default();
//! let mut a = MetaInfoDescription::new();
//! a.name = "A".into();
//! let mut b = a.clone();
//! b.name = "B".into();
//! // Names alone do not order two descriptions.
//! assert!(a != b);
//! assert!(!a.source_lt(&b, &registry)? && !b.source_lt(&a, &registry)?);
//! // With differing meta values, the name decides.
//! a.metadata.set_value_registered(&mut registry, "label", "x".into())?;
//! b.metadata.set_value_registered(&mut registry, "label", "y".into())?;
//! assert!(a.source_lt(&b, &registry)?);
//! # Ok::<(), openms::Error>(())
//! ```

use std::sync::Arc;

use crate::Result;
use crate::metadata::meta_info_indexed::IndexedMetaInfo;
use crate::metadata::meta_info_registry::MetaInfoRegistry;
use crate::metadata::{DataProcessing, MetaInfo, MetaValue, MetaValueData};

/// Description of a meta-data array (`MetaInfoDescription`).
///
/// | C++ | Rust |
/// |---|---|
/// | default, copy and move constructors and assignments, destructor | [`Default`], [`Clone`], moves, [`Drop`] |
/// | `MetaInfoInterface` base | the public [`metadata`](Self::metadata) field |
/// | `getName()`, `setName()` | the public [`name`](Self::name) field |
/// | `getDataProcessing()` (const and mutable), `setDataProcessing()` | the public [`data_processing`](Self::data_processing) field |
/// | `operator==`, `operator!=` | [`PartialEq`] (data processing compared by pointee, as `Helpers::cmpPtrSafe`) |
/// | `operator<`, `<=`, `>`, `>=` | [`source_lt`](Self::source_lt), [`source_le`](Self::source_le), [`source_gt`](Self::source_gt), [`source_ge`](Self::source_ge) |
///
/// The source's shared pointers may be null (`getDataProcessing().resize(1)`
/// in its class test); an `Arc` cannot, so the null branches of its ordering
/// have no counterpart.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MetaInfoDescription {
    /// The name of the peak annotations.
    pub name: String,
    /// Meta values, the source's `MetaInfoInterface` base.
    pub metadata: MetaInfo,
    /// Applied data-processing steps, shared as in the source.
    pub data_processing: Vec<Arc<DataProcessing>>,
}

impl MetaInfoDescription {
    /// An empty description, the source's default constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// A borrowed view of this description, for [`description_lt`].
    pub fn view(&self) -> DescriptionView<'_> {
        DescriptionView {
            name: &self.name,
            metadata: &self.metadata,
            data_processing: &self.data_processing,
        }
    }

    /// The source's `operator<` (see the module documentation).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`](crate::Error::InvalidValue) when the meta
    /// values differ and one of them is stored under a name `registry` does
    /// not hold: the source orders by registry index, and its setters always
    /// register the name.
    pub fn source_lt(&self, rhs: &Self, registry: &MetaInfoRegistry) -> Result<bool> {
        description_lt(self.view(), rhs.view(), registry)
    }

    /// The source's `operator<=`: `*this < rhs || *this == rhs`.
    ///
    /// # Errors
    ///
    /// As [`source_lt`](Self::source_lt).
    pub fn source_le(&self, rhs: &Self, registry: &MetaInfoRegistry) -> Result<bool> {
        Ok(self.source_lt(rhs, registry)? || self == rhs)
    }

    /// The source's `operator>`: `!(*this <= rhs)`.
    ///
    /// # Errors
    ///
    /// As [`source_lt`](Self::source_lt).
    pub fn source_gt(&self, rhs: &Self, registry: &MetaInfoRegistry) -> Result<bool> {
        Ok(!self.source_le(rhs, registry)?)
    }

    /// The source's `operator>=`: `!(*this < rhs)`.
    ///
    /// # Errors
    ///
    /// As [`source_lt`](Self::source_lt).
    pub fn source_ge(&self, rhs: &Self, registry: &MetaInfoRegistry) -> Result<bool> {
        Ok(!self.source_lt(rhs, registry)?)
    }
}

/// The three description fields of any type that carries them, borrowed.
///
/// [`MetaInfoDescription::view`](crate::metadata::meta_info_description::MetaInfoDescription::view) and
/// [`DataArray::description_view`](crate::kernel::DataArray::description_view)
/// produce one.
#[derive(Clone, Copy, Debug)]
pub struct DescriptionView<'a> {
    /// The name.
    pub name: &'a str,
    /// The meta values.
    pub metadata: &'a MetaInfo,
    /// The data processing.
    pub data_processing: &'a [Arc<DataProcessing>],
}

impl PartialEq for DescriptionView<'_> {
    /// The source's `MetaInfoDescription::operator==`: equal meta values,
    /// names and data processing (compared by pointee).
    fn eq(&self, rhs: &Self) -> bool {
        self.metadata == rhs.metadata
            && self.name == rhs.name
            && self.data_processing.len() == rhs.data_processing.len()
            && self
                .data_processing
                .iter()
                .zip(rhs.data_processing)
                .all(|(a, b)| a == b)
    }
}

/// The source's `MetaInfoDescription::operator<` over two views (see the
/// module documentation).
///
/// # Errors
///
/// As [`MetaInfoDescription::source_lt`](crate::metadata::meta_info_description::MetaInfoDescription::source_lt).
pub fn description_lt(
    lhs: DescriptionView<'_>,
    rhs: DescriptionView<'_>,
    registry: &MetaInfoRegistry,
) -> Result<bool> {
    if lhs.data_processing.len() != rhs.data_processing.len() {
        return Ok(lhs.data_processing.len() < rhs.data_processing.len());
    }
    if lhs.metadata != rhs.metadata {
        let lhs_empty = lhs.metadata.is_empty();
        let rhs_empty = rhs.metadata.is_empty();
        if lhs_empty != rhs_empty {
            // `bool < bool`: true only for a non-empty lhs and an empty rhs.
            return Ok(!lhs_empty & rhs_empty);
        }
        if lhs.name != rhs.name {
            return Ok(lhs.name < rhs.name);
        }
        let lhs_keys = lhs.metadata.keys_as_indices(registry)?;
        let rhs_keys = rhs.metadata.keys_as_indices(registry)?;
        if lhs_keys != rhs_keys {
            return Ok(lhs_keys < rhs_keys);
        }
        for key in lhs_keys {
            let lhs_text = lhs
                .metadata
                .value_at(registry, key)
                .map(data_value_text)
                .unwrap_or_default();
            let rhs_text = rhs
                .metadata
                .value_at(registry, key)
                .map(data_value_text)
                .unwrap_or_default();
            if lhs_text != rhs_text {
                return Ok(lhs_text < rhs_text);
            }
        }
    }
    for (a, b) in lhs.data_processing.iter().zip(rhs.data_processing) {
        let (lhs_size, rhs_size) = (a.actions.len(), b.actions.len());
        if lhs_size != rhs_size {
            return Ok(lhs_size < rhs_size);
        }
    }
    Ok(false)
}

/// `DataValue::toString()` with its default full precision: strings as they
/// are, integers in decimal, doubles as `StringUtils::toStr(double)`, lists
/// as `ListUtilsIO`'s `[a, b, c]`; the unit is not part of the text.
pub fn data_value_text(value: &MetaValue) -> String {
    fn list<T>(items: &[T], text: impl Fn(&T) -> String) -> String {
        let parts: Vec<String> = items.iter().map(text).collect();
        format!("[{}]", parts.join(", "))
    }
    let float = |v: &f64| crate::param::value::format_float(*v, true);
    match value.data() {
        MetaValueData::Empty => String::new(),
        MetaValueData::String(text) => text.clone(),
        MetaValueData::Integer(number) => number.to_string(),
        MetaValueData::Float(number) => float(number),
        MetaValueData::StringList(items) => list(items, String::clone),
        MetaValueData::IntegerList(items) => list(items, i64::to_string),
        MetaValueData::FloatList(items) => list(items, float),
    }
}

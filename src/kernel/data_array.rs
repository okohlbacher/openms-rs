// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::DataArray;
use crate::metadata::meta_info_description::{
    DescriptionView, MetaInfoDescription, description_lt,
};
use crate::metadata::meta_info_registry::MetaInfoRegistry;
use crate::metadata::{CVTermList, MetaInfo, MetaValue, MetaValueData};
use crate::{Error, Result};
use std::mem::size_of;

impl<T> DataArray<T> {
    /// Validate owned metadata and processing records. Algorithm callers meter
    /// descriptions first when recursive traversal must share a resource budget.
    pub fn validate_description(&self) -> Result<()> {
        crate::metadata::validate_meta(&self.metadata)?;
        for processing in &self.data_processing {
            processing.validate()?;
        }
        Ok(())
    }
    /// Charge the complete annotation description before cloning, validating or
    /// destroying it. Shared processing payload is conservatively counted even
    /// when a clone would only increment its Arc reference count.
    pub(crate) fn description_with_budget(
        &self,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<()> {
        let mut meter = Meter { work, bytes };
        meter.text(&self.name)?;
        meter.meta(&self.metadata)?;
        meter
            .slots::<std::sync::Arc<crate::metadata::DataProcessing>>(self.data_processing.len())?;
        for processing in &self.data_processing {
            meter.slots::<crate::metadata::DataProcessing>(1)?;
            meter.text(&processing.software.name)?;
            meter.text(&processing.software.version)?;
            meter.cv(&processing.software.cv_terms)?;
            meter.tree::<crate::metadata::ProcessingAction>(processing.actions.len())?;
            meter.meta(&processing.metadata)?;
        }
        Ok(())
    }

    /// Copy the source description onto a newly generated array. Callers with
    /// resource budgets must charge description_with_budget first.
    pub(crate) fn copy_description_to<U>(&self, target: &mut DataArray<U>) {
        target.metadata = self.metadata.clone();
        target.data_processing = self.data_processing.clone();
    }

    /// Whether a name-and-values-only transport would lose information.
    pub fn has_description_metadata(&self) -> bool {
        !self.metadata.is_empty() || !self.data_processing.is_empty()
    }

    /// The array's `MetaInfoDescription` base, borrowed: its name, meta
    /// values and data processing.
    pub fn description_view(&self) -> DescriptionView<'_> {
        DescriptionView {
            name: &self.name,
            metadata: &self.metadata,
            data_processing: &self.data_processing,
        }
    }

    /// The array's `MetaInfoDescription` base as an owned value.
    pub fn description(&self) -> MetaInfoDescription {
        MetaInfoDescription {
            name: self.name.clone(),
            metadata: self.metadata.clone(),
            data_processing: self.data_processing.clone(),
        }
    }
}

impl<T: PartialOrd> DataArray<T> {
    /// The source's `FloatDataArray`/`IntegerDataArray`/`StringDataArray`
    /// `operator<` (`METADATA/DataArrays.h`): when the descriptions differ,
    /// their order ([`description_lt`]); otherwise `std::vector`'s `<`, which
    /// is `std::lexicographical_compare` and is reproduced element by
    /// element, so an unordered pair (a NaN) is neither less nor greater and
    /// the comparison moves on. See `docs/META_INFO_DESCRIPTION_SUPPORT.md`.
    ///
    /// # Errors
    ///
    /// As [`MetaInfoDescription::source_lt`](crate::metadata::meta_info_description::MetaInfoDescription::source_lt).
    pub fn source_lt(&self, rhs: &Self, registry: &MetaInfoRegistry) -> Result<bool> {
        let (lhs_description, rhs_description) = (self.description_view(), rhs.description_view());
        if lhs_description != rhs_description {
            return description_lt(lhs_description, rhs_description, registry);
        }
        for (a, b) in self.data.iter().zip(&rhs.data) {
            if a < b {
                return Ok(true);
            }
            if b < a {
                return Ok(false);
            }
        }
        Ok(self.data.len() < rhs.data.len())
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

fn limit() -> Error {
    Error::InvalidValue("data array description resource limit exceeded".into())
}
pub(crate) struct Meter<'a> {
    pub(crate) work: &'a mut usize,
    pub(crate) bytes: &'a mut usize,
}
impl Meter<'_> {
    pub(crate) fn charge(&mut self, work: usize, bytes: usize) -> Result<()> {
        *self.work = self.work.checked_sub(work).ok_or_else(limit)?;
        *self.bytes = self.bytes.checked_sub(bytes).ok_or_else(limit)?;
        Ok(())
    }
    pub(crate) fn slots<T>(&mut self, count: usize) -> Result<()> {
        self.charge(count, count.checked_mul(size_of::<T>()).ok_or_else(limit)?)
    }
    pub(crate) fn text(&mut self, text: &str) -> Result<()> {
        self.charge(text.len(), text.len())
    }
    pub(crate) fn tree<T>(&mut self, count: usize) -> Result<()> {
        if count != 0 {
            // BTree nodes can be sparsely occupied. Cover their root, spare
            // element slots and child pointers, not just logical element bytes.
            self.charge(count, 512)?;
            self.slots::<(T, [usize; 4])>(count.checked_mul(3).ok_or_else(limit)?.max(11))?;
        }
        Ok(())
    }
    pub(crate) fn meta(&mut self, meta: &MetaInfo) -> Result<()> {
        self.tree::<(String, MetaValue)>(meta.len())?;
        for (name, value) in meta {
            self.text(name)?;
            self.value(value)?;
        }
        Ok(())
    }
    pub(crate) fn meta_update(&mut self, meta: &MetaInfo, key: &str) -> Result<()> {
        let height = (usize::BITS - meta.len().max(1).leading_zeros()) as usize;
        self.charge(
            height
                .saturating_mul(key.len().saturating_add(1))
                .saturating_mul(16),
            0,
        )?;
        self.tree::<(String, MetaValue)>(height.saturating_add(1))?;
        self.text(key)?;
        if let Some(old) = meta.get(key) {
            self.value(old)?;
        }
        Ok(())
    }
    pub(crate) fn value(&mut self, value: &MetaValue) -> Result<()> {
        if let Some(unit) = value.unit() {
            self.text(unit.accession())?;
            self.text(unit.name())?;
            self.text(unit.cv_ref())?;
        }
        match value.data() {
            MetaValueData::String(value) => self.text(value)?,
            MetaValueData::StringList(values) => {
                self.slots::<String>(values.len())?;
                for value in values {
                    self.text(value)?;
                }
            }
            MetaValueData::IntegerList(values) => self.slots::<i64>(values.len())?,
            MetaValueData::FloatList(values) => self.slots::<f64>(values.len())?,
            _ => self.charge(1, 0)?,
        }
        Ok(())
    }
    pub(crate) fn cv(&mut self, terms: &CVTermList) -> Result<()> {
        self.meta(&terms.metadata)?;
        self.tree::<(String, Vec<crate::metadata::CVTerm>)>(terms.terms().len())?;
        for (key, values) in terms.terms() {
            self.text(key)?;
            self.slots::<crate::metadata::CVTerm>(values.len())?;
            for value in values {
                self.text(&value.accession)?;
                self.text(&value.name)?;
                self.text(&value.cv_ref)?;
                self.value(&value.value)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod sparse_tree_tests {
    use super::*;
    #[test]
    fn singleton_tree_charges_the_full_large_inline_root_before_allocation() {
        type Large = (String, crate::metadata::Instrument);
        let required = 512 + 11 * size_of::<(Large, [usize; 4])>();
        let mut work = usize::MAX;
        let mut bytes = required - 1;
        assert!(
            Meter {
                work: &mut work,
                bytes: &mut bytes
            }
            .tree::<Large>(1)
            .is_err()
        );
        let mut bytes = required;
        Meter {
            work: &mut work,
            bytes: &mut bytes,
        }
        .tree::<Large>(1)
        .unwrap();
        assert_eq!(bytes, 0);
    }
    #[test]
    fn replacing_typed_metadata_precharges_old_list_and_sparse_map_atomically() {
        let meta = MetaInfo::from([("key".into(), vec!["x".repeat(100); 4].into())]);
        let original = meta.clone();
        let mut work = 100_000;
        let mut bytes = 100_000;
        Meter {
            work: &mut work,
            bytes: &mut bytes,
        }
        .meta_update(&meta, "key")
        .unwrap();
        let cost = 100_000 - bytes;
        let mut work = 100_000;
        let mut bytes = cost - 1;
        assert!(
            Meter {
                work: &mut work,
                bytes: &mut bytes
            }
            .meta_update(&meta, "key")
            .is_err()
        );
        assert_eq!(meta, original);
        assert!(work < 100_000);
    }
}

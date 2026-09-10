// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The source's typed protein-group quantity transport and ownership guard.

use crate::identification::ProteinGroup;
use crate::kernel::DataArray;
use crate::metadata::{MetaInfo, MetaValue, MetaValueData};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

const FLOATS: &[&str] = &[
    "psm_count",
    "distinct_peptides",
    "file_channel_level_abundance",
];
const LEGACY_FLOATS: &[&str] = &[
    "abundances",
    "psm_count",
    "distinct_peptides",
    "file_channel_level_abundance",
];
const INTEGERS: &[&str] = &["file_channel_level_channel", "file_level_psm_count"];
const STRINGS: &[&str] = &["file_channel_level_filename", "file_level_filename"];
fn bad(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}
fn owned_suffix<'a>(key: &'a str, group: &str) -> Option<&'a str> {
    let rest = key.strip_prefix(group)?.strip_prefix('_')?;
    let digits = rest.split('_').next()?;
    (!digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit())).then_some(rest)
}

pub(super) fn read(
    meta: &mut MetaInfo,
    name: &str,
    refs: &BTreeMap<String, String>,
    warnings: &mut Vec<String>,
    max_items: usize,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<Vec<ProteinGroup>> {
    let mut owned = MetaInfo::new();
    let mut kept = MetaInfo::new();
    for (key, value) in std::mem::take(meta) {
        if owned_suffix(&key, name).is_some() {
            owned.insert(key, value);
        } else {
            kept.insert(key, value);
        }
    }
    *meta = kept;
    let mut result = Vec::new();
    for index in 0..max_items {
        let key = format!("{name}_{index}");
        let Some(base) = owned.remove(&key) else {
            break;
        };
        let text = base.as_str()?;
        let mut parts = text.split(',');
        let probability: f64 = parts
            .next()
            .ok_or_else(|| bad("missing group probability"))?
            .parse()
            .map_err(|_| bad("invalid group probability"))?;
        let mut group = ProteinGroup {
            probability,
            ..Default::default()
        };
        for reference in parts {
            if group.accessions.len() >= max_items {
                return Err(bad("protein group exceeds list limit"));
            }
            let accession = refs
                .get(reference)
                .ok_or_else(|| bad(format!("unknown protein group reference {reference}")))?;
            // Short PH identifiers can each resolve to a very long accession.
            // Charge resolved storage and comparisons before expanding references.
            super::charge(bytes, accession.len().saturating_add(32))?;
            super::charge(work, accession.len().saturating_mul(4).saturating_add(1))?;
            group.accessions.push(accession.clone());
        }
        if group.accessions.is_empty() {
            return Err(bad("protein group needs at least one reference"));
        }
        let prefix = format!("{key}_");
        let guard = owned.remove(&format!("{key}_quantified_proteins"));
        let has_guard = guard.is_some();
        let matches = if let Some(guard) = guard {
            guard.as_string_list()? == group.accessions
        } else {
            false
        };
        let mut floats = BTreeMap::new();
        let mut integers = BTreeMap::new();
        let mut strings = BTreeMap::new();
        // Remove only this group's contiguous entries: never rebuild the
        // remaining map once per group (quadratic on quantified experiments).
        let keys: Vec<_> = owned
            .range(prefix.clone()..)
            .take_while(|(key, _)| key.starts_with(&prefix))
            .map(|(key, _)| key.clone())
            .collect();
        for key in keys {
            let value = owned.remove(&key).expect("collected existing quantity key");
            if !matches {
                continue;
            }
            if value.unit().is_some() {
                return Err(bad("protein group quantities cannot have units"));
            }
            let array = key[prefix.len()..].to_owned();
            match value.data() {
                MetaValueData::FloatList(values) => {
                    if values.len() > max_items {
                        return Err(bad("protein group quantity list exceeds limit"));
                    }
                    let values = values
                        .iter()
                        .map(|v| {
                            let value = *v as f32;
                            if value.is_finite() {
                                Ok(value)
                            } else {
                                Err(bad("protein group float exceeds f32"))
                            }
                        })
                        .collect::<Result<Vec<_>>>()?;
                    floats.insert(array, values);
                }
                MetaValueData::IntegerList(values) => {
                    if values.len() > max_items {
                        return Err(bad("protein group quantity list exceeds limit"));
                    }
                    integers.insert(
                        array,
                        values
                            .iter()
                            .map(|&v| {
                                i32::try_from(v)
                                    .map_err(|_| bad("protein quantity integer exceeds i32"))
                            })
                            .collect::<Result<Vec<_>>>()?,
                    );
                }
                MetaValueData::StringList(values) => {
                    if values.len() > max_items {
                        return Err(bad("protein group quantity list exceeds limit"));
                    }
                    strings.insert(array, values.clone());
                }
                _ => {} // Source-owned non-list leftovers are not quantity arrays.
            }
        }
        if has_guard && !matches {
            warnings.push(format!(
                "discarded quantities for {key}: protein ownership differs"
            ));
        }
        if !floats.is_empty() || !integers.is_empty() || !strings.is_empty() {
            let legacy = floats.contains_key("abundances");
            group.float_data_arrays =
                canonical(if legacy { LEGACY_FLOATS } else { FLOATS }, floats);
            group.integer_data_arrays = canonical(INTEGERS, integers);
            group.string_data_arrays = canonical(STRINGS, strings);
            if legacy {
                let len = group.float_data_arrays[0].data.len();
                for i in 1..=2 {
                    if group.float_data_arrays[i].data.is_empty() {
                        group.float_data_arrays[i].data.resize(len, 0.0);
                    }
                }
            }
        }
        group.validate()?;
        result.push(group);
    }
    if owned.contains_key(&format!("{name}_{max_items}")) {
        return Err(bad("protein group count exceeds limit"));
    }
    Ok(result)
}

fn canonical<T>(names: &[&str], mut parsed: BTreeMap<String, Vec<T>>) -> Vec<DataArray<T>> {
    let mut arrays = Vec::with_capacity(names.len() + parsed.len());
    for &name in names {
        arrays.push(DataArray::new(
            name,
            parsed.remove(name).unwrap_or_default(),
        ));
    }
    arrays.extend(
        parsed
            .into_iter()
            .map(|(name, data)| DataArray::new(name, data)),
    );
    arrays
}

pub(super) fn write(
    meta: &mut MetaInfo,
    name: &str,
    groups: &[ProteinGroup],
    refs: &BTreeMap<String, String>,
) -> Result<()> {
    // Regenerate all owned entries. Unlike the source, also remove stale base
    // entries beyond the new group count so filtering cannot resurrect a group.
    meta.retain(|key, _| owned_suffix(key, name).is_none());
    for (index, group) in groups.iter().enumerate() {
        group.validate()?;
        if group.accessions.is_empty() {
            return Err(bad("protein group needs at least one accession"));
        }
        let key = format!("{name}_{index}");
        let mut entry = group.probability.to_string();
        for accession in &group.accessions {
            entry.push(',');
            entry.push_str(
                refs.get(accession)
                    .ok_or_else(|| bad(format!("unknown protein accession {accession}")))?,
            );
        }
        meta.insert(key.clone(), entry.into());
        let mut names = BTreeSet::new();
        let mut wrote = false;
        for array in &group.float_data_arrays {
            if array.name.is_empty()
                || array.data.is_empty()
                || (matches!(array.name.as_str(), "psm_count" | "distinct_peptides")
                    && array.data.iter().all(|v| *v == 0.0))
            {
                continue;
            }
            check_name(&mut names, &array.name)?;
            meta.insert(
                format!("{key}_{}", array.name),
                MetaValue::try_from(array.data.iter().map(|v| f64::from(*v)).collect::<Vec<_>>())?,
            );
            wrote = true;
        }
        for array in &group.integer_data_arrays {
            if array.name.is_empty()
                || array.data.is_empty()
                || (matches!(array.name.as_str(), "psm_count" | "distinct_peptides")
                    && array.data.iter().all(|v| *v == 0))
            {
                continue;
            }
            check_name(&mut names, &array.name)?;
            meta.insert(
                format!("{key}_{}", array.name),
                array
                    .data
                    .iter()
                    .map(|v| i64::from(*v))
                    .collect::<Vec<_>>()
                    .into(),
            );
            wrote = true;
        }
        for array in &group.string_data_arrays {
            if array.name.is_empty() || array.data.is_empty() {
                continue;
            }
            check_name(&mut names, &array.name)?;
            meta.insert(format!("{key}_{}", array.name), array.data.clone().into());
            wrote = true;
        }
        if wrote {
            meta.insert(
                format!("{key}_quantified_proteins"),
                group.accessions.clone().into(),
            );
        }
    }
    Ok(())
}
fn check_name<'a>(names: &mut BTreeSet<&'a str>, name: &'a str) -> Result<()> {
    if name == "quantified_proteins" || !names.insert(name) {
        return Err(bad("duplicate or reserved protein quantity array name"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn short_protein_references_cannot_expand_past_shared_payload_limit() {
        let refs = BTreeMap::from([("PH_0".into(), "A".repeat(512))]);
        let mut metadata = MetaInfo::from([("protein_group_0".into(), "0.9,PH_0,PH_0".into())]);
        let mut remaining = 700;
        assert!(
            read(
                &mut metadata,
                "protein_group",
                &refs,
                &mut Vec::new(),
                100,
                &mut 100_000,
                &mut remaining
            )
            .is_err()
        );
        assert_eq!(remaining, 156);
    }
}

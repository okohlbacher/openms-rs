// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Helpers for the Ion Identity Molecular Networking (IIMN) workflow on GNPS.
//!
//! Port of `ANALYSIS/ID/IonIdentityMolecularNetworking.h`
//! (`IonIdentityMolecularNetworking.cpp`). See `docs/ID_MAPPER_SUPPORT.md`
//! (section "IonIdentityMolecularNetworking") for the API mapping.
//!
//! IIMN extends GNPS Feature-Based Molecular Networking by grouping consensus
//! features that represent the same compound seen as different adducts (for
//! example `[M+H]+` and `[M+Na]+`). The grouping is usually produced upstream
//! by an adduct decharger (`MetaboliteAdductDecharger`), which records a list of
//! adduct-group IDs in the `LinkedGroups` meta value of each consensus
//! feature.
//!
//! Two stages live here:
//!
//! 1. [`annotate_consensus_map`](crate::analysis::id::ion_identity_molecular_networking::annotate_consensus_map) consumes `LinkedGroups`, builds a
//!    bipartite graph (consensus feature ↔ adduct group), computes its
//!    connected components and writes back three IIMN meta values: a 1-based
//!    row ID, the network number and the list of partner-feature row IDs. The
//!    linked-groups meta value is removed.
//! 2. [`write_supplementary_pair_table`](crate::analysis::id::ion_identity_molecular_networking::write_supplementary_pair_table) turns those annotations into the CSV
//!    file GNPS IIMN takes as its "supplementary pairs" upload.
//!
//! The output of stage 1 is also a precondition of the GNPS meta-value and
//! quantification tables written for the `GNPSExport` tool.
//!
//! The source builds the graph with Boost.Graph (`adjacency_list` plus
//! `connected_components`). The component numbering that function produces —
//! components numbered in order of their lowest vertex index, vertices in the
//! order they were added — is reproduced by a breadth-first labelling over the
//! same vertex order, so no graph crate is needed.

use crate::concept::constants::user_param::{
    IIMN_ADDUCT_PARTNERS, IIMN_ANNOTATION_NETWORK_NUMBER, IIMN_BEST_ION, IIMN_LINKED_GROUPS,
    IIMN_ROW_ID,
};
use crate::format::sv_out_stream::{QuotingMethod, SVOutStream};
use crate::kernel::features::ConsensusMap;
use crate::metadata::MetaValue;
use crate::param::value::format_float;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::io::Write;
use std::path::Path;

/// Most graph vertices plus edges [`annotate_consensus_map`] builds before
/// refusing. Native ceiling; the source has none.
pub const MAX_GRAPH_ELEMENTS: usize = 50_000_000;

/// Source `std::stoull` acceptance: whether the C library's conversion
/// ([`crate::concept::string_conversion::stoul`]) converts. Only whether it
/// throws matters here, because the group vertex's numeric label is never read
/// again.
fn stoull_accepts(text: &str) -> bool {
    crate::concept::string_conversion::stoul(text).is_ok()
}

/// Populate the IIMN meta values of `consensus_map` in place, the source
/// static `annotateConsensusMap`.
///
/// Reads [`IIMN_LINKED_GROUPS`] (a string list of adduct-group IDs), builds a
/// bipartite graph whose vertices are consensus features and adduct groups,
/// computes its connected components, and writes per feature:
///
/// - [`IIMN_ROW_ID`] — the 1-based index of the feature, on *every* feature,
///   including those without a group annotation;
/// - [`IIMN_ANNOTATION_NETWORK_NUMBER`] — the 1-based connected-component
///   number, shared by direct and transitive adduct partners; only on features
///   that carried [`IIMN_LINKED_GROUPS`];
/// - [`IIMN_ADDUCT_PARTNERS`] — the `;`-separated row IDs of the partners in
///   ascending feature order; only on features sharing at least one group with
///   another feature.
///
/// [`IIMN_LINKED_GROUPS`] is removed from every feature. Groups are compared
/// as text, so `"01"` and `"1"` are different groups, as in the source.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a linked-groups value is not a string
/// list (the source `DataValue::toStringList` throws), when a group ID is not
/// a number `std::stoull` accepts (the source throws `std::invalid_argument` or
/// `std::out_of_range`), or when the graph would exceed
/// [`MAX_GRAPH_ELEMENTS`]. The map is unchanged on error.
pub fn annotate_consensus_map(consensus_map: &mut ConsensusMap) -> Result<()> {
    // Vertices in insertion order: Some(feature index) or None for a group.
    let mut vertex_feature: Vec<Option<usize>> = Vec::new();
    let mut adjacency: Vec<Vec<usize>> = Vec::new();
    let mut group_vertex: HashMap<String, usize> = HashMap::new();
    let mut elements = 0usize;
    for (i, feature) in consensus_map.features.iter().enumerate() {
        let Some(value) = feature.metadata.get(IIMN_LINKED_GROUPS) else {
            continue;
        };
        let groups = value.as_string_list().map_err(|_| {
            Error::InvalidValue(format!(
                "IonIdentityMolecularNetworking: '{IIMN_LINKED_GROUPS}' must be a string list"
            ))
        })?;
        let feature_vertex = vertex_feature.len();
        vertex_feature.push(Some(i));
        adjacency.push(Vec::new());
        elements = elements.saturating_add(1 + 2 * groups.len());
        if elements > MAX_GRAPH_ELEMENTS {
            return Err(Error::InvalidValue(
                "IonIdentityMolecularNetworking: graph exceeds its size limit".into(),
            ));
        }
        for group in groups {
            let target = match group_vertex.get(group.as_str()) {
                Some(&vertex) => vertex,
                None => {
                    if !stoull_accepts(group) {
                        return Err(Error::InvalidValue(format!(
                            "IonIdentityMolecularNetworking: adduct group '{group}' is not an unsigned integer"
                        )));
                    }
                    let vertex = vertex_feature.len();
                    vertex_feature.push(None);
                    adjacency.push(Vec::new());
                    group_vertex.insert(group.clone(), vertex);
                    vertex
                }
            };
            adjacency[feature_vertex].push(target);
            adjacency[target].push(feature_vertex);
        }
    }

    // Boost connected_components numbering: components in order of their
    // first vertex.
    let mut component = vec![usize::MAX; vertex_feature.len()];
    let mut next = 0usize;
    for start in 0..vertex_feature.len() {
        if component[start] != usize::MAX {
            continue;
        }
        component[start] = next;
        let mut queue = VecDeque::from([start]);
        while let Some(v) = queue.pop_front() {
            for &w in &adjacency[v] {
                if component[w] == usize::MAX {
                    component[w] = next;
                    queue.push_back(w);
                }
            }
        }
        next += 1;
    }

    // Plan every meta-value change before touching the map.
    let mut network: Vec<(usize, i64)> = Vec::new();
    let mut partners: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    for (v, feature) in vertex_feature.iter().enumerate() {
        let Some(feature) = *feature else { continue };
        let number = i64::try_from(component[v] + 1)
            .map_err(|_| Error::InvalidValue("network number overflow".into()))?;
        network.push((feature, number));
        for &group in &adjacency[v] {
            for &partner in &adjacency[group] {
                if partner == v {
                    continue;
                }
                if let Some(partner_feature) = vertex_feature[partner] {
                    partners.entry(feature).or_default().insert(partner_feature);
                }
            }
        }
    }
    let row_id = |index: usize| -> Result<MetaValue> {
        i64::try_from(index + 1)
            .map(MetaValue::from)
            .map_err(|_| Error::InvalidValue("row ID overflow".into()))
    };

    // Commit.
    for i in 0..consensus_map.features.len() {
        let id = row_id(i)?;
        consensus_map.features[i]
            .metadata
            .insert(IIMN_ROW_ID.into(), id);
    }
    for (feature, number) in network {
        consensus_map.features[feature].metadata.insert(
            IIMN_ANNOTATION_NETWORK_NUMBER.into(),
            MetaValue::from(number),
        );
    }
    for (feature, set) in partners {
        let text = set
            .iter()
            .map(|&j| (j + 1).to_string())
            .collect::<Vec<_>>()
            .join(";");
        consensus_map.features[feature]
            .metadata
            .insert(IIMN_ADDUCT_PARTNERS.into(), MetaValue::from(text));
    }
    for feature in &mut consensus_map.features {
        feature.metadata.remove(IIMN_LINKED_GROUPS);
    }
    Ok(())
}

/// Write the IIMN supplementary-pairs CSV consumed by GNPS, the source static
/// `writeSupplementaryPairTable`, to `writer`.
///
/// The table has the header `ID1,ID2,EdgeType,Score,Annotation` and one row per
/// unordered pair of features that share at least one adduct group (the
/// reverse edge is dropped as the file is written), in ascending order of the
/// first feature and then its partners:
///
/// - `ID1`, `ID2` — the [`IIMN_ROW_ID`] of the two features;
/// - `EdgeType` — the fixed text `MS1 annotation`;
/// - `Score` — the number of direct partners of both features minus 2;
/// - `Annotation` — the [`IIMN_BEST_ION`] of each feature (`default` when
///   unset) and `dm/z=` followed by the absolute m/z difference, formatted as
///   `StringUtils::toStr(double)`.
///
/// Fields are separated by commas; a comma inside a field is replaced by `_`,
/// without quoting (source `SVOutStream(…, ",", "_", QuotingMethod::NONE)`).
/// The schema is GNPS's
/// <https://ccms-ucsd.github.io/GNPSDocumentation/fbmn-iin/#supplementary-pairs>.
///
/// Returns `Ok(false)` without writing anything when the first feature has no
/// [`IIMN_ROW_ID`], i.e. [`annotate_consensus_map`] has not been run; the
/// source returns silently there. An empty map is treated the same way — the
/// source reads its first element without checking.
///
/// The score is computed in unsigned 64-bit arithmetic with wrap-around, as
/// the source's `size_t` expression is: a partner list that names a feature
/// without a partner list of its own therefore scores `2^64 - 1`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when an [`IIMN_ADDUCT_PARTNERS`] value is
/// not a string, or one of its entries is not a row ID of this map (the source
/// throws from `std::stoi` or indexes out of range), and propagates write
/// failures.
pub fn write_supplementary_pair_table(
    consensus_map: &ConsensusMap,
    writer: impl Write,
) -> Result<bool> {
    let Some(first) = consensus_map.features.first() else {
        return Ok(false);
    };
    if !first.metadata.contains_key(IIMN_ROW_ID) {
        return Ok(false);
    }
    let count = consensus_map.features.len();
    let mut sorted: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    for (i, feature) in consensus_map.features.iter().enumerate() {
        let Some(value) = feature.metadata.get(IIMN_ADDUCT_PARTNERS) else {
            continue;
        };
        let text = value.as_str().map_err(|_| {
            Error::InvalidValue(format!(
                "IonIdentityMolecularNetworking: '{IIMN_ADDUCT_PARTNERS}' must be a string"
            ))
        })?;
        for entry in text.split(';') {
            let index = entry
                .trim()
                .parse::<usize>()
                .ok()
                .and_then(|row| row.checked_sub(1))
                .filter(|&index| index < count)
                .ok_or_else(|| {
                    Error::InvalidValue(format!(
                        "IonIdentityMolecularNetworking: partner '{entry}' is not a row ID of this map"
                    ))
                })?;
            sorted.entry(i).or_default().insert(index);
        }
    }
    let num_partners: HashMap<usize, u64> = sorted
        .iter()
        .map(|(&key, set)| (key, set.len() as u64))
        .collect();

    let mut out = SVOutStream::with_options(writer, ",", "_", QuotingMethod::None)?;
    for header in ["ID1", "ID2", "EdgeType", "Score", "Annotation"] {
        out.write_field(header)?;
    }
    out.newline()?;
    let best_ion = |index: usize| -> String {
        consensus_map.features[index]
            .metadata
            .get(IIMN_BEST_ION)
            .map_or_else(|| "default".to_owned(), ToString::to_string)
    };
    let row = |index: usize| -> String {
        consensus_map.features[index]
            .metadata
            .get(IIMN_ROW_ID)
            .map(ToString::to_string)
            .unwrap_or_default()
    };
    let keys: Vec<usize> = sorted.keys().copied().collect();
    for key in keys {
        let current: Vec<usize> = sorted
            .get(&key)
            .map(|set| set.iter().copied().collect())
            .unwrap_or_default();
        for partner in current {
            out.write_display(row(key))?;
            out.write_display(row(partner))?;
            out.write_field("MS1 annotation")?;
            let score = num_partners
                .get(&key)
                .copied()
                .unwrap_or(0)
                .wrapping_add(num_partners.get(&partner).copied().unwrap_or(0))
                .wrapping_sub(2);
            out.write_number(score)?;
            let delta = (consensus_map.features[key].mz - consensus_map.features[partner].mz).abs();
            out.write_field(&format!(
                "{} {} dm/z={}",
                best_ion(key),
                best_ion(partner),
                format_float(delta, true)
            ))?;
            out.newline()?;
            if let Some(set) = sorted.get_mut(&partner) {
                set.remove(&key);
            }
        }
    }
    out.flush()?;
    Ok(true)
}

/// [`write_supplementary_pair_table`] into a file at `path`, created or
/// truncated. As in the source, no file is created when the map carries no
/// IIMN annotation.
///
/// # Errors
///
/// As [`write_supplementary_pair_table`], plus [`Error::Io`] when the file
/// cannot be created.
pub fn store_supplementary_pair_table(
    consensus_map: &ConsensusMap,
    path: impl AsRef<Path>,
) -> Result<bool> {
    let annotated = consensus_map
        .features
        .first()
        .is_some_and(|f| f.metadata.contains_key(IIMN_ROW_ID));
    if !annotated {
        return Ok(false);
    }
    let file = std::fs::File::create(path)?;
    write_supplementary_pair_table(consensus_map, std::io::BufWriter::new(file))
}

#[cfg(test)]
mod tests {
    use super::stoull_accepts;

    #[test]
    fn stoull_acceptance_follows_the_c_library() {
        assert!(stoull_accepts("12"));
        assert!(stoull_accepts(" +7"));
        assert!(stoull_accepts("-1"));
        assert!(stoull_accepts("3abc"));
        assert!(!stoull_accepts("abc"));
        assert!(!stoull_accepts(""));
        assert!(!stoull_accepts("99999999999999999999"));
    }
}

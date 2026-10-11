// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Renaming the tools of an old INI file: the native form of the
//! `OpenMS4-cli` package's `APPLICATIONS/INIUpdater.h` and `INIUpdater.cpp`
//! (cli `c19e494`).
//!
//! Earlier OpenMS releases had one tool per task with a `-type` choosing the
//! algorithm (`FeatureFinder -type centroided`); later releases split those
//! into one tool each (`FeatureFinderCentroided`). [`INIUpdater`] holds the
//! fixed old-name-to-new-name table and answers which tool of this release an
//! INI section of an older one belongs to. The `INIUpdater` utility is its one
//! consumer.
//!
//! See `docs/TOPP_CLI_SUPPORT.md` (*Tool base classes*) for the member mapping.

use super::tool_handler::ToolHandler;
use crate::Result;
use crate::data_structures::ToolDescriptionInternal;
use crate::data_structures::list;
use crate::param::Param;
use std::collections::BTreeMap;
use std::sync::OnceLock;

/// One entry of the old-to-new tool table: the old tool with its `-type`, and
/// the tool that replaces it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ToolMappingEntry {
    /// The old tool, with the `-type` value that selected the algorithm.
    pub old: ToolDescriptionInternal,
    /// The tool of this release that replaces it.
    pub new: ToolDescriptionInternal,
}

/// Each old tool mapped to its new name (source typedef `ToolMapping`,
/// `std::map<Internal::ToolDescriptionInternal,
/// Internal::ToolDescriptionInternal>`).
///
/// The source map is ordered, and looked up, by `ToolDescriptionInternal`'s
/// `operator<`, which compares only the text `<name>.<types joined by ','>`
/// ([`ToolDescriptionInternal::sort_key`]): two descriptions with the same
/// name and types are the same key, whatever their category or internal
/// flag. This map is keyed by that text, so it orders and finds exactly as
/// the source's does.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ToolMapping {
    entries: BTreeMap<String, ToolMappingEntry>,
}

impl ToolMapping {
    /// Map `old` to `new`, replacing an entry with the same key (source
    /// `map_[old] = new`).
    pub fn insert(&mut self, old: ToolDescriptionInternal, new: ToolDescriptionInternal) {
        self.entries
            .insert(old.sort_key(), ToolMappingEntry { old, new });
    }

    /// The tool that replaces `old`, found by the source map's key
    /// equivalence (source `map_.find(old)`).
    pub fn get(&self, old: &ToolDescriptionInternal) -> Option<&ToolDescriptionInternal> {
        self.entries.get(&old.sort_key()).map(|entry| &entry.new)
    }

    /// Whether `old` has an entry (source `map_.contains(old)`).
    pub fn contains(&self, old: &ToolDescriptionInternal) -> bool {
        self.entries.contains_key(&old.sort_key())
    }

    /// Number of entries (source `size()`).
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the table is empty (source `empty()`).
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// The entries in the source map's order.
    pub fn iter(&self) -> impl Iterator<Item = &ToolMappingEntry> {
        self.entries.values()
    }
}

/// A description with `name` and the types of the comma-separated `types`,
/// as the source's `TDE(name, ListUtils::create<std::string>(types))`.
fn tool(name: &str, types: &str) -> ToolDescriptionInternal {
    ToolDescriptionInternal::with_types(name, &split_types(types))
}

/// `ListUtils::create<std::string>(text)`: the text split at every comma,
/// with nothing for an empty text.
///
/// The split cannot fail for strings (a `String` item is taken verbatim, and a
/// comma is one byte, so it never splits a UTF-8 character); should the list
/// helper ever refuse, the text is kept whole, which is what the one-item case
/// gives.
fn split_types(types: &str) -> Vec<String> {
    list::create::<String>(types, b',').unwrap_or_else(|_| vec![types.to_owned()])
}

/// The fixed table, in the source's order of insertion
/// (`INIUpdater.cpp:43-63`).
fn build_mapping() -> ToolMapping {
    let mut map = ToolMapping::default();
    let rows: [(&str, &str, &str); 13] = [
        ("FeatureFinder", "centroided", "FeatureFinderCentroided"),
        (
            "FeatureFinder",
            "isotope_wavelet",
            "FeatureFinderIsotopeWavelet",
        ),
        ("FeatureFinder", "mrm", "FeatureFinderMRM"),
        ("FeatureLinker", "labeled", "FeatureLinkerLabeled"),
        ("FeatureLinker", "unlabeled", "FeatureLinkerUnlabeled"),
        ("FeatureLinker", "unlabeled_qt", "FeatureLinkerUnlabeledQT"),
        ("NoiseFilter", "gaussian", "NoiseFilterGaussian"),
        ("NoiseFilter", "sgolay", "NoiseFilterSGolay"),
        ("MapAligner", "apply_given_trafo", "MapRTTransformer"),
        ("MapAligner", "identification", "MapAlignerIdentification"),
        ("MapAligner", "pose_clustering", "MapAlignerPoseClustering"),
        ("PeakPicker", "wavelet", "PeakPickerWavelet"),
        ("PeakPicker", "high_res", "PeakPickerHiRes"),
    ];
    for (old, old_type, new) in rows {
        map.insert(tool(old, old_type), tool(new, ""));
    }
    map
}

/// Updates an INI file of an older release (source class `INIUpdater`).
///
/// The source keeps the table in a static member that the first
/// `getNameMapping` call fills and that every instance shares; here it is an
/// immutable table built once per process ([`OnceLock`]), which every
/// instance reads.
#[derive(Clone, Copy, Debug, Default)]
pub struct INIUpdater;

impl INIUpdater {
    /// An updater (source constructor, which fills the shared table).
    pub fn new() -> Self {
        let _ = Self::mapping();
        Self
    }

    fn mapping() -> &'static ToolMapping {
        static MAPPING: OnceLock<ToolMapping> = OnceLock::new();
        MAPPING.get_or_init(build_mapping)
    }

    /// The tools an INI file configures (source `getToolNamesFromINI`,
    /// `INIUpdater.cpp:26-39`): for every entry whose name holds exactly one
    /// `:` and ends in `:version`, the text before the `:`, in the
    /// parameter tree's traversal order.
    ///
    /// A tool section of an INI file is `<tool>:<instance>:…`, so its
    /// `version` entry, `<tool>:version`, is the one top-level entry with a
    /// single colon; the source relies on that and so does this.
    ///
    /// # Errors
    ///
    /// Returns the parameter tree's traversal error, when the tree exceeds
    /// its bounds.
    pub fn get_tool_names_from_ini(&self, ini: &Param) -> Result<Vec<String>> {
        let mut tool_names = Vec::new();
        for item in ini.iter()? {
            let name = item.key.as_str();
            if name.bytes().filter(|&byte| byte == b':').count() == 1 && name.ends_with(":version")
            {
                if let Some((prefix, _)) = name.split_once(':') {
                    tool_names.push(prefix.to_owned());
                }
            }
        }
        Ok(tool_names)
    }

    /// The old-to-new tool table (source `getNameMapping`,
    /// `INIUpdater.cpp:41-67`): thirteen entries mapping `FeatureFinder`,
    /// `FeatureLinker`, `NoiseFilter`, `MapAligner` and `PeakPicker` with a
    /// `-type` to the tool that replaced them, each new tool with an empty
    /// type list.
    pub fn get_name_mapping(&self) -> &'static ToolMapping {
        Self::mapping()
    }

    /// The name of the tool that replaces `old_name` (source
    /// `getNewToolName`, `INIUpdater.cpp:69-97`), or `None` when there is
    /// none.
    ///
    /// `tools_type` is the old tool's `-type`, and empty when it had none.
    /// The lookups, in order:
    ///
    /// 1. the table entry for `old_name` with the types of `tools_type`
    ///    (split at commas, as `ListUtils::create`);
    /// 2. the table entry for `old_name` without types, so a type is ignored
    ///    when a mapping without a type exists;
    /// 3. `old_name` itself when it is a TOPP tool of this release
    ///    ([`ToolHandler::get_topp_tool_list`]), whatever the type.
    ///
    /// The source returns `bool` and sets an out-parameter, which it empties
    /// first; `None` stands for `false` and its empty name.
    ///
    /// # Errors
    ///
    /// As [`ToolHandler::get_topp_tool_list`], which the third lookup asks:
    /// a malformed or duplicate manifest.
    pub fn get_new_tool_name(
        &self,
        old_name: &str,
        tools_type: &str,
        registry: &ToolHandler,
    ) -> Result<Option<String>> {
        let map = Self::mapping();
        let old_with_type = tool(old_name, tools_type);
        if let Some(new) = map.get(&old_with_type) {
            return Ok(Some(new.name.clone()));
        }
        let old_without_type = ToolDescriptionInternal::with_types(old_name, &[]);
        if let Some(new) = map.get(&old_without_type) {
            return Ok(Some(new.name.clone()));
        }
        if registry.get_topp_tool_list()?.contains_key(old_name) {
            return Ok(Some(old_name.to_owned()));
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_table_is_keyed_by_the_source_sort_key() {
        let map = INIUpdater::new().get_name_mapping();
        assert_eq!(map.len(), 13);
        let keys: Vec<String> = map.iter().map(|entry| entry.old.sort_key()).collect();
        let mut sorted = keys.clone();
        sorted.sort();
        assert_eq!(keys, sorted);
        // Category and the internal flag are not part of the key.
        let mut other = tool("PeakPicker", "wavelet");
        other.category = "Signal processing".into();
        other.is_internal = true;
        assert_eq!(
            map.get(&other).map(|new| new.name.as_str()),
            Some("PeakPickerWavelet")
        );
    }

    #[test]
    fn an_empty_type_text_is_no_type() {
        assert!(split_types("").is_empty());
        assert_eq!(split_types("a,b"), vec!["a".to_owned(), "b".to_owned()]);
    }
}

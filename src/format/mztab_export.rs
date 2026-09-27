// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! MzTab exporters: feature maps, consensus maps and identification runs to an
//! [`MzTab`](crate::format::mztab::MzTab) document — the export surface of `FORMAT/MzTab.h`.
//!
//! The source keeps the data model and the exporters in one class. This port
//! splits them: [`crate::format::mztab`] is the data model and this module is
//! the export layer, adding four associated functions to
//! [`MzTab`](crate::format::mztab::MzTab) (`export_feature_map`,
//! `export_identifications`, `export_consensus_map` and
//! `extract_modification_list`)
//! plus the two row streams the source nests in the class,
//! [`IdMzTabStream`](crate::format::mztab_export::IdMzTabStream) (`MzTab::IDMzTabStream`) and
//! [`ConsensusMapMzTabStream`](crate::format::mztab_export::ConsensusMapMzTabStream)
//! (`MzTab::CMMzTabStream`). The streams are what
//! [`MzTabFile::store_identifications`](crate::format::mztab_file::MzTabFile::store_identifications)
//! and
//! [`MzTabFile::store_consensus_map`](crate::format::mztab_file::MzTabFile::store_consensus_map)
//! write from, one row at a time.
//!
//! Every exporter reproduces the source's row construction, including its
//! defects, because the retained MzTabExporter outputs depend on them; each
//! preserved defect is named at the item that carries it and listed in
//! `docs/MZTAB_SUPPORT.md`. The recurring ones:
//!
//! * the target/decoy annotation is renamed to a PRIDE column
//!   (`opt_global_cv_MS:1002217_decoy_peptide` in `PEP`/`PSM`,
//!   `opt_global_cv_PRIDE:0000303_decoy_hit` in `PRT`) and its values mapped to
//!   `0`/`1`;
//! * an optional column is written from the *first* row entry of its name, so
//!   two metadata keys that differ only by spaces versus underscores share one
//!   column;
//! * fixed modifications are compared by the modification's short identifier
//!   against the search parameters' *full* identifiers, so they are never
//!   recognised as fixed and are always reported.
//!
//! ```
//! use openms::format::mztab::MzTab;
//! use openms::format::mztab_export::IdentificationExportOptions;
//! use openms::identification::{PeptideIdentification, ProteinIdentification};
//!
//! let mut run = ProteinIdentification::new();
//! run.identifier = "run_0".into();
//! run.primary_ms_run_paths = vec!["file.mzML".into()];
//! let mut spectrum = PeptideIdentification::new();
//! spectrum.identifier = "run_0".into();
//!
//! // A spectrum without hits yields one placeholder PSM row whose sequence is
//! // null, exactly as the source's stream does.
//! let document = MzTab::export_identifications(
//!     &[run],
//!     &[spectrum],
//!     "ids.idXML",
//!     &IdentificationExportOptions::default(),
//! )?;
//! assert_eq!(document.meta_data.mz_tab_type.get(), "Identification");
//! assert_eq!(document.psm_data.len(), 1);
//! assert!(document.psm_data[0].sequence.is_null());
//! # Ok::<(), openms::Error>(())
//! ```

use crate::chemistry::SequenceModification;
use crate::concept::constants::user_param::{
    ID_MERGE_INDEX, LOCALIZED_MODIFICATIONS_USERPARAM, SIGNIFICANCE_THRESHOLD,
};
use crate::format::file_types::strip_extension;
use crate::format::mztab::{
    MzTab, MzTabBoolean, MzTabDouble, MzTabDoubleList, MzTabInteger, MzTabMSRunMetaData,
    MzTabMetaData, MzTabModification, MzTabModificationList, MzTabModificationMetaData,
    MzTabOptionalColumnEntry, MzTabPSMSectionRow, MzTabParameter, MzTabParameterList,
    MzTabPeptideSectionRow, MzTabProteinSectionRow, MzTabSoftwareMetaData, MzTabSpectraRef,
    MzTabString, add_meta_info_to_optional_columns, meta_value_text, modification_metadata,
};
use crate::identification::{
    PeptideHit, PeptideIdentification, ProteinGroup, ProteinHit, ProteinIdentification,
};
use crate::kernel::{ConsensusFeature, ConsensusMap, Feature, FeatureMap};
use crate::metadata::{ExperimentalDesign, MetaInfo};
use crate::param::value::format_float;
use crate::system::file::basename;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

/// Title the identification exporters write when the caller gives none.
pub const DEFAULT_IDENTIFICATION_TITLE: &str = "ID export from OpenMS";
/// Title the consensus-map exporters write when the caller gives none.
pub const DEFAULT_CONSENSUS_MAP_TITLE: &str = "ConsensusMap export from OpenMS";

/// Optional column that carries the OpenMS peptidoform string in `PEP`/`PSM`.
const PEPTIDOFORM_COLUMN: &str = "opt_global_cv_MS:1000889_peptidoform_sequence";
/// The target/decoy optional column as the metadata key spells it.
const TARGET_DECOY_COLUMN: &str = "opt_global_target_decoy";
/// PRIDE's decoy column for `PEP` and `PSM` rows.
const DECOY_PEPTIDE_COLUMN: &str = "opt_global_cv_MS:1002217_decoy_peptide";
/// PRIDE's decoy column for `PRT` rows.
const DECOY_HIT_COLUMN: &str = "opt_global_cv_PRIDE:0000303_decoy_hit";
/// The `PRT` column that says which kind of protein record a row is.
const RESULT_TYPE_COLUMN: &str = "opt_global_result_type";

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

fn missing(message: impl Into<String>) -> Error {
    Error::MissingInformation(message.into())
}

fn parameter(text: &str) -> Result<MzTabParameter> {
    let mut value = MzTabParameter::default();
    value.from_cell_string(text)?;
    Ok(value)
}

fn string(text: &str) -> MzTabString {
    MzTabString::from_text(text)
}

fn one_value_list(value: f64) -> MzTabDoubleList {
    let mut list = MzTabDoubleList::default();
    list.set(vec![MzTabDouble::new(value)]);
    list
}

fn int_cell(value: i64, what: &str) -> Result<MzTabInteger> {
    i32::try_from(value)
        .map(MzTabInteger::new)
        .map_err(|_| bad(format!("{what} does not fit an MzTab Integer cell")))
}

/// The `file://`-prefixed location the source writes for an MS run path.
fn file_location(path: &str) -> MzTabString {
    if path.starts_with("file://") {
        string(path)
    } else {
        string(&format!("file://{path}"))
    }
}

/// `MzTabSpectraRef::setMSFile` as a release build of the source executes it:
/// the assertion is compiled out, so a zero index is ignored rather than
/// rejected and the reference keeps its previous (null) run.
fn set_ms_file_lenient(reference: &mut MzTabSpectraRef, index: usize) {
    if index >= 1 {
        // Cannot fail for a nonzero index.
        let _ = reference.set_ms_file(index);
    }
}

/// `MzTabSpectraRef::setSpecRef` as a release build executes it: empty text is
/// logged and ignored.
fn set_spec_ref_lenient(reference: &mut MzTabSpectraRef, spec_ref: &str) {
    if !spec_ref.is_empty() {
        // Cannot fail for nonempty text.
        let _ = reference.set_spec_ref(spec_ref);
    }
}

/// A key with every space replaced by an underscore, as the source's
/// `StringUtils::substitute(s, ' ', '_')`.
fn underscored(key: &str) -> String {
    key.replace(' ', "_")
}

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------

/// Options of [`MzTab::export_identifications`], [`IdMzTabStream::new`] and
/// [`MzTabFile::store_identifications`](crate::format::mztab_file::MzTabFile::store_identifications).
///
/// The source passes these as positional arguments; the three with defaults
/// there have the same defaults here, and `first_run_inference_only`, which the
/// source requires, defaults to `false`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct IdentificationExportOptions {
    /// Whether all protein inference information (groups and scores) is stored
    /// in the first run. Only honoured when that run carries inference data.
    pub first_run_inference_only: bool,
    /// Export spectra without PSMs as well. Source default `false`.
    pub export_empty_pep_ids: bool,
    /// Export every PSM of a spectrum instead of only the best one. Source
    /// default `false`.
    pub export_all_psms: bool,
    /// The metadata section's `title`. Source default
    /// [`DEFAULT_IDENTIFICATION_TITLE`].
    pub title: String,
}

impl Default for IdentificationExportOptions {
    fn default() -> Self {
        Self {
            first_run_inference_only: false,
            export_empty_pep_ids: false,
            export_all_psms: false,
            title: DEFAULT_IDENTIFICATION_TITLE.into(),
        }
    }
}

/// Options of [`MzTab::export_consensus_map`], [`ConsensusMapMzTabStream::new`]
/// and
/// [`MzTabFile::store_consensus_map`](crate::format::mztab_file::MzTabFile::store_consensus_map).
///
/// The first four are required positional arguments in the source and default
/// to `false` here; the remaining three have the source's defaults.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ConsensusMapExportOptions {
    /// Whether all protein inference information (groups and scores) is stored
    /// in the first run. Only honoured when that run carries inference data.
    pub first_run_inference_only: bool,
    /// Export consensus features that carry no peptide hit.
    pub export_unidentified_features: bool,
    /// Also export the map's unassigned peptide identifications as PSMs.
    pub export_unassigned_ids: bool,
    /// Add per-subfeature m/z and RT entries to every `PEP` row.
    ///
    /// The entries are added to the rows, but the stream's peptide column list
    /// never names them, so the streaming writer
    /// [`MzTabFile::store_consensus_map`](crate::format::mztab_file::MzTabFile::store_consensus_map)
    /// drops them — as in the source, where `-opt_columns subfeatures` of
    /// MzTabExporter therefore has no visible effect. The in-memory
    /// [`MzTab::export_consensus_map`] keeps them, because a document's column
    /// list is collected from its rows.
    pub export_subfeatures: bool,
    /// Export spectra without PSMs as well. Source default `false`.
    pub export_empty_pep_ids: bool,
    /// Export every PSM of a spectrum instead of only the best one. Source
    /// default `false`.
    pub export_all_psms: bool,
    /// The metadata section's `title`. Source default
    /// [`DEFAULT_CONSENSUS_MAP_TITLE`].
    pub title: String,
}

impl Default for ConsensusMapExportOptions {
    fn default() -> Self {
        Self {
            first_run_inference_only: false,
            export_unidentified_features: false,
            export_unassigned_ids: false,
            export_subfeatures: false,
            export_empty_pep_ids: false,
            export_all_psms: false,
            title: DEFAULT_CONSENSUS_MAP_TITLE.into(),
        }
    }
}

// ---------------------------------------------------------------------------
// Target/decoy remapping and optional columns
// ---------------------------------------------------------------------------

/// Rename the target/decoy column to `new_header` and map its values.
///
/// Source `remapTargetDecoyPSMAndPeptideSection_` and
/// `remapTargetDecoyProteinSection_`, which differ only in the new name:
/// `target` and `target+decoy` become `0`, `decoy` becomes `1`, anything else
/// is kept verbatim. An entry already carrying the new name is remapped too.
fn remap_target_decoy(opt: &mut [MzTabOptionalColumnEntry], new_header: &str) {
    for entry in opt.iter_mut() {
        if entry.name == TARGET_DECOY_COLUMN || entry.name == new_header {
            entry.name = new_header.to_owned();
            let replacement = match entry.value.get() {
                "target" | "target+decoy" => Some("0"),
                "decoy" => Some("1"),
                _ => None,
            };
            if let Some(value) = replacement {
                entry.value = string(value);
            }
        }
    }
}

fn keys_of(meta: &MetaInfo) -> impl Iterator<Item = &String> {
    meta.keys()
}

/// Append one `opt_global_<key>` entry per key, filled by `fill`.
fn push_entries<F>(opt: &mut Vec<MzTabOptionalColumnEntry>, keys: &BTreeSet<String>, fill: F)
where
    F: Fn(&str) -> MzTabString,
{
    for key in keys {
        opt.push(MzTabOptionalColumnEntry::new(
            format!("opt_global_{key}"),
            fill(key),
        ));
    }
}

/// Set every entry named `opt_global_<key with spaces as underscores>` to the
/// value of `key` in `meta`, for every key of `meta`.
///
/// The source loops over `getKeys()` and, for each, over the whole entry list
/// without breaking, so every entry of a matching name is overwritten. Where
/// two keys map to one name the later key wins; the source visits keys in the
/// global `MetaInfoRegistry` order, this port in lexical order.
fn overwrite_entries(opt: &mut [MzTabOptionalColumnEntry], meta: &MetaInfo) {
    for (key, value) in meta {
        let name = format!("opt_global_{}", underscored(key));
        let text = meta_value_text(value);
        for entry in opt.iter_mut() {
            if entry.name == name {
                entry.value = string(&text);
            }
        }
    }
}

fn set_entry(opt: &mut [MzTabOptionalColumnEntry], name: &str, value: &str) {
    for entry in opt.iter_mut() {
        if entry.name == name {
            entry.value = string(value);
        }
    }
}

// ---------------------------------------------------------------------------
// Modifications
// ---------------------------------------------------------------------------

/// The MzTab identifier of `modification`.
///
/// Source `getModificationIdentifier_`: the UniMod accession upper-cased
/// (`UNIMOD:21`) when the record has one, otherwise `CHEMMOD:` followed by the
/// monoisotopic mass difference in the source's full-precision number format.
fn modification_identifier(modification: &SequenceModification) -> Result<MzTabString> {
    if let Some(record) = modification.known() {
        if let Some(accession) = record.unimod_accession() {
            if !accession.is_empty() {
                return Ok(string(&accession.to_uppercase()));
            }
        }
    }
    let mass = modification.diff_mono_mass()?;
    Ok(string(&format!("CHEMMOD:{}", format_float(mass, true))))
}

/// The source's `ResidueModification::getId()`: the short identifier, which
/// is what the fixed-modification test compares.
fn short_id(modification: &SequenceModification) -> &str {
    modification.name()
}

impl MzTab {
    /// The MzTab modification list of a peptide hit.
    ///
    /// Source `MzTab::extractModificationList`. Every modification of the hit's
    /// sequence is reported — N-terminal at position `0`, residues at their
    /// one-based position, C-terminal at `length + 1` — except those whose
    /// short identifier is listed in `fixed_mods`. When `localization_mods` is
    /// not empty, a residue modification whose full identifier (`Phospho (S)`)
    /// it lists carries the hit's `Luciphor_global_flr` as the parameter
    /// `[MS, MS:1002380, false localization rate, <value>]`; without that
    /// metadata the parameter is null.
    ///
    /// The source documents that fixed modifications are not reported in the
    /// `PRT` and `PEP` sections, but it compares
    /// `ResidueModification::getId()` (`Carbamidomethyl`) against the search
    /// parameters' full identifiers (`Carbamidomethyl (C)`), so the comparison
    /// never matches when the list comes from search parameters and fixed
    /// modifications are reported after all. The comparison is preserved: the
    /// retained MzTabExporter outputs show those modifications.
    ///
    /// # Errors
    ///
    /// [`Error::Unsupported`] (from the sequence chemistry) when an anonymous
    /// mass tag without a UniMod record has no known mass difference, and
    /// [`Error::Parse`] when the `Luciphor_global_flr` value turns the
    /// localization parameter into unparseable text.
    pub fn extract_modification_list(
        pep_hit: &PeptideHit,
        fixed_mods: &[String],
        localization_mods: &[String],
    ) -> Result<MzTabModificationList> {
        let sequence = &pep_hit.sequence;
        let has_loc_mods = !localization_mods.is_empty();
        let mut localization_score = MzTabParameter::default();
        if has_loc_mods {
            if let Some(flr) = pep_hit.metadata.get("Luciphor_global_flr") {
                localization_score = parameter(&format!(
                    "[MS,MS:1002380,false localization rate,{}]",
                    meta_value_text(flr)
                ))?;
            }
        }
        let is_fixed = |modification: &SequenceModification| {
            fixed_mods.iter().any(|name| name == short_id(modification))
        };
        let mut mods = Vec::new();
        if sequence.is_modified() {
            if let Some(modification) = sequence.n_terminal_modification() {
                if !is_fixed(modification) {
                    let mut entry = MzTabModification::default();
                    entry.set_modification_identifier(modification_identifier(modification)?);
                    entry.set_positions_and_parameters(vec![(0, MzTabParameter::default())]);
                    mods.push(entry);
                }
            }
            for index in 0..sequence.len() {
                if let Some(modification) = sequence.residue_modification(index)? {
                    if !is_fixed(modification) {
                        let localized = has_loc_mods
                            && localization_mods
                                .iter()
                                .any(|name| name == modification.full_id());
                        let score = if localized {
                            localization_score.clone()
                        } else {
                            MzTabParameter::default()
                        };
                        let mut entry = MzTabModification::default();
                        entry.set_positions_and_parameters(vec![(index + 1, score)]);
                        entry.set_modification_identifier(modification_identifier(modification)?);
                        mods.push(entry);
                    }
                }
            }
            if let Some(modification) = sequence.c_terminal_modification() {
                if !is_fixed(modification) {
                    let mut entry = MzTabModification::default();
                    entry.set_positions_and_parameters(vec![(
                        sequence.len() + 1,
                        MzTabParameter::default(),
                    )]);
                    entry.set_modification_identifier(modification_identifier(modification)?);
                    mods.push(entry);
                }
            }
        }
        let mut list = MzTabModificationList::default();
        list.set(mods);
        Ok(list)
    }
}

/// `fixed_mod`/`variable_mod` metadata without the empty-list placeholder.
///
/// Source `generateMzTabStringFromModifications`, which the identification and
/// consensus exporters call directly: an empty list yields no entry at all. A
/// name the registry does not know leaves a gap in the numbering.
fn modifications_metadata(names: &[String]) -> Result<BTreeMap<usize, MzTabModificationMetaData>> {
    Ok(modification_metadata(names)?.metadata)
}

// ---------------------------------------------------------------------------
// Run and search-engine bookkeeping (the source's protected statics)
// ---------------------------------------------------------------------------

/// Source `mapIDRunIdentifier2IDRunIndex_`: identifier to run index, the last
/// run of a repeated identifier winning.
fn map_run_identifier_to_index(prot_ids: &[&ProteinIdentification]) -> BTreeMap<String, usize> {
    let mut map = BTreeMap::new();
    for (index, run) in prot_ids.iter().enumerate() {
        map.insert(run.identifier.clone(), index);
    }
    map
}

/// Source `mapBetweenMSFileNameAndMSRunIndex_`.
///
/// Every distinct primary MS run path gets the next one-based `ms_run` index.
/// A run without paths is registered under its index as a string key, the
/// source's "hack" for mixed annotated and unannotated runs.
fn map_ms_file_names(
    prot_ids: &[&ProteinIdentification],
    skip_first: bool,
) -> (BTreeMap<String, usize>, BTreeMap<usize, String>) {
    let mut name_to_index = BTreeMap::new();
    let mut index_to_name = BTreeMap::new();
    let mut current = 1usize;
    for (position, run) in prot_ids.iter().enumerate() {
        if skip_first && position == 0 {
            continue;
        }
        if run.primary_ms_run_paths.is_empty() {
            let key = current.to_string();
            index_to_name.entry(current).or_insert_with(|| key.clone());
            name_to_index.entry(key).or_insert(current);
            current += 1;
        } else {
            for path in &run.primary_ms_run_paths {
                index_to_name.entry(current).or_insert_with(|| path.clone());
                if !name_to_index.contains_key(path) {
                    name_to_index.insert(path.clone(), current);
                    current += 1;
                }
            }
        }
    }
    (name_to_index, index_to_name)
}

/// Source `addMSRunMetaData_`: one `ms_run[n]-location` per registered file.
fn add_ms_run_meta_data(index_to_name: &BTreeMap<usize, String>, meta: &mut MzTabMetaData) {
    for (index, name) in index_to_name {
        let run = MzTabMSRunMetaData {
            location: file_location(name),
            ..MzTabMSRunMetaData::default()
        };
        meta.ms_run.insert(*index, run);
    }
}

/// Source `mapIDRunFileIndex2MSFileIndex_`: `(run index, file index within
/// the run)` to `ms_run` index. A run without paths maps file `0` to its *run*
/// index, not to an `ms_run` index, as the source does.
fn map_run_file_to_ms_file(
    prot_ids: &[&ProteinIdentification],
    name_to_index: &BTreeMap<String, usize>,
    skip_first: bool,
) -> Result<BTreeMap<(usize, usize), usize>> {
    let mut map = BTreeMap::new();
    for (run_index, run) in prot_ids.iter().enumerate() {
        if skip_first && run_index == 0 {
            continue;
        }
        if run.primary_ms_run_paths.is_empty() {
            map.insert((run_index, 0), run_index);
        } else {
            for (file_index, path) in run.primary_ms_run_paths.iter().enumerate() {
                let ms_file = name_to_index
                    .get(path)
                    .copied()
                    .ok_or_else(|| missing(format!("MS run path {path:?} is not registered")))?;
                map.insert((run_index, file_index), ms_file);
            }
        }
    }
    Ok(map)
}

/// Search engine name, version and PSM score type.
type EngineKey = (String, String, String);
/// Ordered `(key, value)` settings of one search engine.
type EngineSettings = Vec<(String, String)>;

/// What [`map_runs_and_search_engines`] collects.
#[derive(Default)]
struct SearchEngines {
    to_runs: BTreeMap<EngineKey, BTreeSet<usize>>,
    run_to_engines: BTreeMap<usize, Vec<(String, String)>>,
    to_settings: BTreeMap<String, EngineSettings>,
}

/// Source `mapBetweenRunAndSearchEngines_`.
///
/// The main engine of each run is its search engine; `SE:<name>` keys of the
/// search parameters (left by PercolatorAdapter or ConsensusID) add the
/// engines they replaced. The PSM score type of a run is that of the first
/// peptide identification carrying the run's identifier, or `unknown_score`.
/// Only the first settings seen for an engine name are kept.
fn map_runs_and_search_engines(
    prot_ids: &[&ProteinIdentification],
    pep_ids: &[&PeptideIdentification],
    skip_first: bool,
) -> Result<SearchEngines> {
    let mut engines = SearchEngines::default();
    for (run_index, run) in prot_ids.iter().enumerate() {
        if skip_first && run_index == 0 {
            continue;
        }
        let score_type = pep_ids
            .iter()
            .find(|pep| pep.identifier == run.identifier)
            .map_or_else(|| "unknown_score".to_owned(), |pep| pep.score_type.clone());
        engines
            .to_runs
            .entry((
                run.search_engine.clone(),
                run.search_engine_version.clone(),
                score_type,
            ))
            .or_default()
            .insert(run_index);
        let list = engines.run_to_engines.entry(run_index).or_default();
        list.push((run.search_engine.clone(), run.search_engine_version.clone()));
        for (key, value) in &run.search_parameters.metadata {
            if let Some(name) = key.strip_prefix("SE:") {
                list.push((name.to_owned(), meta_value_text(value)));
            }
        }
        for (name, _) in list.iter() {
            let settings = run.search_engine_settings_as_pairs(name)?;
            engines.to_settings.entry(name.clone()).or_insert(settings);
        }
    }
    Ok(engines)
}

/// Source `getSearchModifications_`: variable and fixed modification names of
/// every run, each list sorted and deduplicated.
fn search_modifications(prot_ids: &[&ProteinIdentification]) -> (Vec<String>, Vec<String>) {
    let mut variable = Vec::new();
    let mut fixed = Vec::new();
    for run in prot_ids {
        variable.extend(run.search_parameters.variable_modifications.iter().cloned());
        fixed.extend(run.search_parameters.fixed_modifications.iter().cloned());
    }
    variable.sort();
    variable.dedup();
    fixed.sort();
    fixed.dedup();
    (variable, fixed)
}

fn named_float<'a>(group: &'a ProteinGroup, name: &str) -> Option<&'a [f32]> {
    group
        .float_data_arrays
        .iter()
        .rfind(|array| array.name == name)
        .map(|array| array.data.as_slice())
}

fn named_integer<'a>(group: &'a ProteinGroup, name: &str) -> Option<&'a [i32]> {
    group
        .integer_data_arrays
        .iter()
        .rfind(|array| array.name == name)
        .map(|array| array.data.as_slice())
}

const ASSAY_ABUNDANCE: &str = "fraction_group_level_abundance";
const ASSAY_FRACTION_GROUP: &str = "fraction_group_level_fraction_group";
const ASSAY_LABEL: &str = "fraction_group_level_label";

/// Source `getQuantStudyVariables_`: the length of the assay abundance arrays
/// when every indistinguishable group carries the three parallel assay arrays
/// with equal lengths, otherwise zero. Only its non-zero-ness is used.
fn quant_study_variables(run: &ProteinIdentification) -> usize {
    let mut count = 0;
    for group in &run.indistinguishable_groups {
        match (
            named_float(group, ASSAY_ABUNDANCE),
            named_integer(group, ASSAY_FRACTION_GROUP),
            named_integer(group, ASSAY_LABEL),
        ) {
            (Some(abundance), Some(groups), Some(labels))
                if abundance.len() == groups.len() && abundance.len() == labels.len() =>
            {
                count = abundance.len();
            }
            _ => return 0,
        }
    }
    count
}

/// Source `getProteinScoreType_`.
fn protein_score_type(run: &ProteinIdentification) -> Result<MzTabParameter> {
    if run.has_inference_data()? {
        parameter(&format!(
            "[,,{} {},]",
            run.inference_engine()?,
            run.score_type
        ))
    } else {
        parameter(&format!("[,,one-peptide-rule {},]", run.score_type))
    }
}

/// Source `addSearchMetaData_`: protein score types, inference engines, search
/// engines with their settings, and the PSM/peptide score types.
fn add_search_meta_data(
    prot_ids: &[&ProteinIdentification],
    engines: &SearchEngines,
    meta: &mut MzTabMetaData,
    first_run_inference_only: bool,
) -> Result<()> {
    let mut protein_score_types = BTreeSet::new();
    let mut inference_engines = BTreeSet::new();
    for run in prot_ids {
        protein_score_types.insert(protein_score_type(run)?.to_cell_string());
        if run.has_inference_data()? {
            inference_engines.insert((
                run.inference_engine()?.to_owned(),
                run.inference_engine_version()?.to_owned(),
            ));
        }
        if first_run_inference_only {
            break;
        }
    }
    for (offset, text) in protein_score_types.iter().enumerate() {
        let term = if text.contains("q-value") {
            parameter("[MS,MS:1003117,OpenMS:Target-decoy protein q-value, ]")?
        } else if text.contains("Epifany") {
            parameter("[MS,MS:1003119,EPIFANY:Protein posterior probability,]")?
        } else {
            parameter(text)?
        };
        meta.protein_search_engine_score.insert(offset + 1, term);
    }
    for (engine, version) in &inference_engines {
        let software = match engine.as_str() {
            "Epifany" => parameter(&format!("[MS,MS:1003118,EPIFANY,{version}]"))?,
            "TOPPProteinInference" => {
                parameter(&format!("[MS,MS:1002203,TOPP ProteinInference,{version}]"))?
            }
            _ => parameter(&format!("[,,{engine},{version}]"))?,
        };
        let index = meta.software.len() + 1;
        meta.software.insert(
            index,
            MzTabSoftwareMetaData {
                software,
                ..MzTabSoftwareMetaData::default()
            },
        );
    }
    let first_index = meta.software.len() + 1;
    for (software_index, (name, version, _)) in (first_index..).zip(engines.to_runs.keys()) {
        let software = if name == "Comet" {
            parameter(&format!("[MS,MS:1002251,Comet,{version}]"))?
        } else if name == "Percolator" {
            parameter(&format!("[MS,MS:1001490,Percolator,{version}]"))?
        } else if name == "MS-GF+" || name == "MSGFPlus" {
            parameter(&format!("[MS,MS:1002048,MS-GF+,{version}]"))?
        } else if name == "XTandem" {
            parameter(&format!("[MS,MS:1001476,X!Tandem,{version}]"))?
        } else if name.contains("ConsensusID") {
            parameter(&format!("[MS,MS:1002188,TOPP ConsensusID,{version}]"))?
        } else {
            parameter(&format!("[,,{name},{version}]"))?
        };
        let mut entry = MzTabSoftwareMetaData {
            software,
            ..MzTabSoftwareMetaData::default()
        };
        let settings = engines
            .to_settings
            .get(name)
            .ok_or_else(|| missing(format!("no settings recorded for search engine {name:?}")))?;
        for (offset, (key, value)) in settings.iter().enumerate() {
            let value = if value.is_empty() { "null" } else { value };
            entry
                .setting
                .insert(offset + 1, string(&format!("{key}:{value}")));
        }
        meta.software.insert(software_index, entry);
    }
    for (offset, (name, _, score)) in engines.to_runs.keys().enumerate() {
        let mut psm = parameter(&format!("[,,{name} {score},]"))?;
        let mut peptide = parameter("[MS,MS:1003114,OpenMS:Best PSM Score,]")?;
        if score == "peptide-level q-value" {
            peptide = if name == "Percolator" {
                parameter("[MS,MS:1002360,distinct peptide-level FDRScore,Percolator]")?
            } else {
                parameter("[MS,MS:1003116,OpenMS:Target-decoy peptide q-value,]")?
            };
            psm = peptide.clone();
        } else if score.contains("q-value") {
            psm = if name == "Percolator" {
                parameter("[MS,MS:1001491,percolator:Q value,]")?
            } else {
                parameter("[MS,MS:1003115,OpenMS:Target-decoy PSM q-value,]")?
            };
        } else if score == "Posterior Error Probability" || score == "pep" {
            if name.contains("ConsensusID") {
                // StringUtils::suffix(name, '_'): the text after the last '_',
                // or the whole name when there is none.
                let algorithm = name.rsplit('_').next().unwrap_or(name);
                psm = parameter(&format!(
                    "[MS,MS:1003113,OpenMS:ConsensusID PEP,{algorithm}]"
                ))?;
            } else if name.contains("Percolator") {
                psm = parameter("[MS,MS:1001493,percolator:PEP,]")?;
            }
        } else if name == "Comet" {
            if score == "expect" {
                psm = parameter("[MS,MS:1002257,Comet:expectation value,]")?;
            }
        } else if (name == "MSGFPlus" || name == "MS-GF+") && score == "SpecEValue" {
            psm = parameter("[MS,MS:1002052,MS-GF:SpecEValue,]")?;
        }
        meta.psm_search_engine_score.insert(offset + 1, psm);
        meta.peptide_search_engine_score.insert(offset + 1, peptide);
    }
    Ok(())
}

/// Source `getMSRunSpectrumIdentifierType_`.
///
/// The source returns from inside its loop on the first iteration whatever it
/// finds, so only the *first* peptide identification's spectrum reference is
/// inspected; the default is `mzML unique identifier`. Preserved.
fn ms_run_spectrum_identifier_type(
    peptide_ids: &[&PeptideIdentification],
) -> Result<MzTabParameter> {
    let default = "[MS,MS:1001530,mzML unique identifier,]";
    let Some(first) = peptide_ids.first() else {
        return parameter(default);
    };
    let reference = first.spectrum_reference();
    let term = if reference.contains("controllerNumber=") {
        "[MS,MS:1000768,Thermo nativeID format,]"
    } else if reference.contains("process=") {
        "[MS,MS:1000769,Waters nativeID format,]"
    } else if reference.contains("cycle=") {
        "[MS,MS:1000770,WIFF nativeID format,]"
    } else if reference.contains("scan=") {
        "[MS,MS:1000776,scan number only nativeID format,]"
    } else if reference.contains("spectrum=") {
        "[MS,MS:1000777,spectrum identifier nativeID format,]"
    } else {
        default
    };
    parameter(term)
}

/// Source `mapGroupsToProteins_`: group index to the indices of the protein
/// hits its accessions name (first hit of an accession), groups whose
/// accessions all miss being absent from the map.
fn map_groups_to_proteins(
    groups: &[ProteinGroup],
    proteins: &[ProteinHit],
) -> BTreeMap<usize, BTreeSet<usize>> {
    let mut first_index: BTreeMap<&str, usize> = BTreeMap::new();
    for (index, protein) in proteins.iter().enumerate() {
        first_index
            .entry(protein.accession.as_str())
            .or_insert(index);
    }
    let mut map: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    for (group_index, group) in groups.iter().enumerate() {
        for accession in &group.accessions {
            if let Some(index) = first_index.get(accession.as_str()) {
                map.entry(group_index).or_default().insert(*index);
            }
        }
    }
    map
}

// ---------------------------------------------------------------------------
// Metadata-key collection
// ---------------------------------------------------------------------------

/// Source `getIdentificationMetaValues_`, with the underscore substitution the
/// source applies to protein-hit, peptide-identification and peptide-hit keys.
fn identification_meta_keys(
    prot_ids: &[&ProteinIdentification],
    peptide_ids: &[&PeptideIdentification],
) -> (BTreeSet<String>, BTreeSet<String>, BTreeSet<String>) {
    let mut protein_hit = BTreeSet::new();
    for run in prot_ids {
        for hit in &run.hits {
            protein_hit.extend(keys_of(&hit.metadata).map(|key| underscored(key)));
        }
    }
    let mut peptide_id = BTreeSet::new();
    let mut peptide_hit = BTreeSet::new();
    for identification in peptide_ids {
        peptide_id.extend(keys_of(&identification.metadata).map(|key| underscored(key)));
        for hit in &identification.hits {
            peptide_hit.extend(keys_of(&hit.metadata).map(|key| underscored(key)));
        }
    }
    peptide_id.remove(SIGNIFICANCE_THRESHOLD);
    (protein_hit, peptide_id, peptide_hit)
}

// ---------------------------------------------------------------------------
// Protein rows
// ---------------------------------------------------------------------------

/// Source `proteinSectionRowFromProteinHit_`.
///
/// The source builds a modification entry for each modification of the hit and
/// then never adds it to the list, so the `modifications` cell of a protein
/// row is always null. Preserved.
fn protein_row_from_hit(
    hit: &ProteinHit,
    db: &MzTabString,
    db_version: &MzTabString,
    user_keys: &BTreeSet<String>,
) -> Result<MzTabProteinSectionRow> {
    let mut row = MzTabProteinSectionRow {
        accession: string(&hit.accession),
        description: string(&hit.description()),
        database: db.clone(),
        database_version: db_version.clone(),
        ..MzTabProteinSectionRow::default()
    };
    row.best_search_engine_score
        .insert(1, MzTabDouble::new(hit.score));
    let coverage = hit.coverage.unwrap_or(-1.0) / 100.0;
    row.coverage = if coverage >= 0.0 {
        MzTabDouble::new(coverage)
    } else {
        MzTabDouble::null()
    };
    add_meta_info_to_optional_columns(user_keys, &mut row.opt, "global", &hit.metadata)?;
    row.opt.push(MzTabOptionalColumnEntry::new(
        RESULT_TYPE_COLUMN,
        string("protein_details"),
    ));
    remap_target_decoy(&mut row.opt, DECOY_HIT_COLUMN);
    Ok(row)
}

/// Source `nextProteinSectionRowFromProteinGroup_`: a general protein group.
fn protein_row_from_group(
    group: &ProteinGroup,
    db: &MzTabString,
    db_version: &MzTabString,
) -> MzTabProteinSectionRow {
    let mut row = MzTabProteinSectionRow {
        database: db.clone(),
        database_version: db_version.clone(),
        ..MzTabProteinSectionRow::default()
    };
    if let Some(first) = group.accessions.first() {
        row.accession = string(first);
    }
    row.ambiguity_members
        .set(group.accessions.iter().map(|a| string(a)).collect());
    row.best_search_engine_score
        .insert(1, MzTabDouble::new(group.probability));
    row.coverage = MzTabDouble::null();
    row.opt.push(MzTabOptionalColumnEntry::new(
        RESULT_TYPE_COLUMN,
        string("general_protein_group"),
    ));
    remap_target_decoy(&mut row.opt, DECOY_HIT_COLUMN);
    row
}

fn int_list(hit: &ProteinHit, key: &str) -> Result<BTreeMap<usize, MzTabInteger>> {
    let mut map = BTreeMap::new();
    if let Some(value) = hit.metadata.get(key) {
        for (offset, number) in value.as_integer_list()?.iter().enumerate() {
            map.insert(offset + 1, int_cell(*number, key)?);
        }
    }
    Ok(map)
}

/// Source `nextProteinSectionRowFromIndistinguishableGroup_`.
///
/// The leader is the protein hit with the smallest index among those the group
/// names; accession, description, taxonomy, species, URI, GO terms,
/// modifications and the per-run counts come from it, the coverage is the mean
/// of the member coverages, and a leading float array named `abundances`
/// fills both the assay and the study-variable abundances.
///
/// The source replaces the modification list on every iteration of its loop,
/// so only the leader's *last* modification (in position order) is reported.
/// Preserved; the source orders equal positions by `ResidueModification`'s
/// comparison, this port by full identifier.
///
/// # Errors
///
/// [`Error::MissingInformation`] when no protein hit carries any of the
/// group's accessions (the source's `std::map::at` throws
/// `std::out_of_range`), and [`Error::InvalidValue`] when `TaxID`, `Species`,
/// `URI`, `GO` or a `num_*_ms_run` value has the wrong type.
fn protein_row_from_indistinguishable_group(
    protein_hits: &[ProteinHit],
    group: &ProteinGroup,
    group_index: usize,
    ind2prot: &BTreeMap<usize, BTreeSet<usize>>,
    db: &MzTabString,
    db_version: &MzTabString,
) -> Result<MzTabProteinSectionRow> {
    let members = ind2prot.get(&group_index).ok_or_else(|| {
        missing("indistinguishable protein group names no protein hit of its run")
    })?;
    let leader_index = *members
        .iter()
        .next()
        .ok_or_else(|| missing("indistinguishable protein group has no member"))?;
    let leader = &protein_hits[leader_index];
    let mut row = MzTabProteinSectionRow {
        database: db.clone(),
        database_version: db_version.clone(),
        accession: string(&leader.accession),
        ..MzTabProteinSectionRow::default()
    };
    let entries: Vec<MzTabString> = group.accessions.iter().map(|a| string(a)).collect();
    let single = entries.len() == 1;
    row.ambiguity_members.set(entries);
    row.opt.push(MzTabOptionalColumnEntry::new(
        RESULT_TYPE_COLUMN,
        string(if single {
            "single_protein"
        } else {
            "indistinguishable_protein_group"
        }),
    ));
    let mut coverage = 0.0_f64;
    for index in members {
        coverage +=
            (1.0 / members.len() as f64) * 0.01 * protein_hits[*index].coverage.unwrap_or(-1.0);
    }
    if coverage >= 0.0 {
        row.coverage = MzTabDouble::new(coverage);
    }
    if let Some(first) = group.float_data_arrays.first() {
        if first.name == "abundances" {
            for (offset, value) in first.data.iter().enumerate() {
                let index = offset + 1;
                row.protein_abundance_assay
                    .insert(index, MzTabDouble::new(f64::from(*value)));
                row.protein_abundance_study_variable
                    .insert(index, MzTabDouble::new(f64::from(*value)));
                row.protein_abundance_stdev_study_variable
                    .insert(index, MzTabDouble::null());
                row.protein_abundance_std_error_study_variable
                    .insert(index, MzTabDouble::null());
            }
        }
    }
    row.description = string(&leader.description());
    if let Some(value) = leader.metadata.get("TaxID") {
        row.taxid = int_cell(value.as_i64()?, "TaxID")?;
    }
    if let Some(value) = leader.metadata.get("Species") {
        row.species = string(value.as_str()?);
    }
    if let Some(value) = leader.metadata.get("URI") {
        row.uri = string(value.as_str()?);
    }
    if let Some(value) = leader.metadata.get("GO") {
        row.go_terms
            .from_cell_string(&value.as_string_list()?.join(","))?;
    }
    row.best_search_engine_score
        .insert(1, MzTabDouble::new(group.probability));
    row.reliability = MzTabInteger::null();
    row.search_engine = MzTabParameterList::default();
    let mut ordered: Vec<_> = leader.modifications.iter().collect();
    ordered.sort_by(|a, b| {
        a.position
            .cmp(&b.position)
            .then_with(|| a.modification.full_id().cmp(b.modification.full_id()))
    });
    if let Some(last) = ordered.last() {
        let mut entry = MzTabModification::default();
        entry.set_modification_identifier(modification_identifier(&last.modification)?);
        entry.set_positions_and_parameters(vec![(last.position + 1, MzTabParameter::default())]);
        row.modifications.set(vec![entry]);
    }
    row.num_psms_ms_run = int_list(leader, "num_psms_ms_run")?;
    row.num_peptides_distinct_ms_run = int_list(leader, "num_peptides_distinct_ms_run")?;
    row.num_peptides_unique_ms_run = int_list(leader, "num_peptides_unique_ms_run")?;
    remap_target_decoy(&mut row.opt, DECOY_HIT_COLUMN);
    Ok(row)
}

/// Source `applyAssayProteinQuantification` (a file-local function of
/// `MzTab.cpp`), used only by the consensus-map stream.
///
/// When the group carries none of the three parallel assay arrays the row is
/// left as it is; when it carries an inconsistent set (missing arrays, unequal
/// lengths, a key below one, a repeated `(fraction group, label)` key) the
/// source logs a warning and also leaves the row; otherwise the abundance
/// columns are rebuilt: one assay column per design assay, filled through the
/// parallel keys, and a study variable copying its assay when it owns exactly
/// one, null otherwise.
fn apply_assay_protein_quantification(
    row: &mut MzTabProteinSectionRow,
    group: &ProteinGroup,
    assay_keys: &[(u32, u32)],
    study_variable_to_assays: &[Vec<usize>],
) -> Result<bool> {
    let abundance = named_float(group, ASSAY_ABUNDANCE);
    let groups = named_integer(group, ASSAY_FRACTION_GROUP);
    let labels = named_integer(group, ASSAY_LABEL);
    if abundance.is_none() && groups.is_none() && labels.is_none() {
        return Ok(false);
    }
    let (Some(abundance), Some(groups), Some(labels)) = (abundance, groups, labels) else {
        return Ok(false);
    };
    if abundance.len() != groups.len() || abundance.len() != labels.len() {
        return Ok(false);
    }
    let mut by_key: BTreeMap<(u32, u32), f32> = BTreeMap::new();
    for index in 0..abundance.len() {
        if groups[index] <= 0 || labels[index] <= 0 {
            return Ok(false);
        }
        let key = (groups[index] as u32, labels[index] as u32);
        if by_key.insert(key, abundance[index]).is_some() {
            return Ok(false);
        }
    }
    row.protein_abundance_assay.clear();
    row.protein_abundance_study_variable.clear();
    row.protein_abundance_stdev_study_variable.clear();
    row.protein_abundance_std_error_study_variable.clear();
    for (offset, key) in assay_keys.iter().enumerate() {
        let value = by_key.get(key).map_or_else(MzTabDouble::null, |value| {
            MzTabDouble::new(f64::from(*value))
        });
        row.protein_abundance_assay.insert(offset + 1, value);
    }
    for (study_variable, assays) in study_variable_to_assays.iter().enumerate().skip(1) {
        let value = if assays.len() == 1 {
            *row.protein_abundance_assay
                .get(&assays[0])
                .ok_or_else(|| missing("study variable references an assay without a column"))?
        } else {
            MzTabDouble::null()
        };
        row.protein_abundance_study_variable
            .insert(study_variable, value);
        row.protein_abundance_stdev_study_variable
            .insert(study_variable, MzTabDouble::null());
        row.protein_abundance_std_error_study_variable
            .insert(study_variable, MzTabDouble::null());
    }
    Ok(true)
}

// ---------------------------------------------------------------------------
// The protein-row state machine shared by both streams
// ---------------------------------------------------------------------------

/// Which kind of protein record the `PRT` state machine emits next.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum ProteinState {
    /// Protein hits of the current run.
    #[default]
    Hits,
    /// General protein groups of the current run.
    Groups,
    /// Indistinguishable groups of the current run.
    Indistinguishable,
}

/// The source's `PRT_STATE_` machine with its four counters.
#[derive(Clone, Debug, Default)]
struct ProteinCursor {
    state: ProteinState,
    run: usize,
    hit: usize,
    group: usize,
    indistinguishable: usize,
    ind2prot: BTreeMap<usize, BTreeSet<usize>>,
}

/// `(assay keys, study variable to assays)` of the consensus stream.
type AssayQuantification<'s> = (&'s [(u32, u32)], &'s [Vec<usize>]);

/// Everything the `PRT` rows need from a stream.
struct ProteinContext<'s, 'a> {
    prot_ids: &'s [&'a ProteinIdentification],
    first_run_inference: bool,
    quant_study_variables: usize,
    db: &'s MzTabString,
    db_version: &'s MzTabString,
    user_keys: &'s BTreeSet<String>,
    /// `(assay keys, study variable to assays)` for the consensus stream.
    assay_quantification: Option<AssayQuantification<'s>>,
}

impl ProteinCursor {
    /// One step of the source's `nextPRTRow`: all protein hits of a run, then
    /// its general groups (only when no quantified groups exist), then its
    /// indistinguishable groups, then the next run — or only the first run when
    /// it holds the inference result.
    fn next(&mut self, context: &ProteinContext<'_, '_>) -> Result<Option<MzTabProteinSectionRow>> {
        loop {
            if context.prot_ids.is_empty()
                || (context.first_run_inference && self.run > 0)
                || self.run >= context.prot_ids.len()
            {
                return Ok(None);
            }
            let run = context.prot_ids[self.run];
            let hits = &run.hits;
            let groups: &[ProteinGroup] = if context.quant_study_variables == 0 {
                &run.protein_groups
            } else {
                &[]
            };
            let indistinguishable = &run.indistinguishable_groups;
            if self.hit == 0 && self.state == ProteinState::Hits {
                self.ind2prot = map_groups_to_proteins(indistinguishable, hits);
            }
            if self.state == ProteinState::Hits {
                if self.hit >= hits.len() {
                    self.hit = 0;
                    self.state = ProteinState::Groups;
                } else {
                    let row = protein_row_from_hit(
                        &hits[self.hit],
                        context.db,
                        context.db_version,
                        context.user_keys,
                    )?;
                    self.hit += 1;
                    return Ok(Some(row));
                }
            }
            if self.state == ProteinState::Groups {
                if self.group >= groups.len() {
                    self.group = 0;
                    self.state = ProteinState::Indistinguishable;
                } else {
                    let row =
                        protein_row_from_group(&groups[self.group], context.db, context.db_version);
                    self.group += 1;
                    return Ok(Some(row));
                }
            }
            if self.indistinguishable >= indistinguishable.len() {
                self.indistinguishable = 0;
                self.hit = 0;
                self.state = ProteinState::Hits;
                self.run += 1;
                continue;
            }
            let group = &indistinguishable[self.indistinguishable];
            let mut row = protein_row_from_indistinguishable_group(
                hits,
                group,
                self.indistinguishable,
                &self.ind2prot,
                context.db,
                context.db_version,
            )?;
            if let Some((assay_keys, study_variables)) = context.assay_quantification {
                apply_assay_protein_quantification(&mut row, group, assay_keys, study_variables)?;
            }
            self.indistinguishable += 1;
            return Ok(Some(row));
        }
    }
}

// ---------------------------------------------------------------------------
// PSM rows
// ---------------------------------------------------------------------------

/// The best hit of one identification: the first hit, replaced by any later
/// hit with a strictly better score. Source `IDFilter::getBestHit` over a
/// one-element list; NaN scores never replace.
fn best_hit(identification: &PeptideIdentification) -> Option<&PeptideHit> {
    let mut best: Option<&PeptideHit> = None;
    for hit in &identification.hits {
        match best {
            None => best = Some(hit),
            Some(current) => {
                let better = if identification.higher_score_better {
                    hit.score > current.score
                } else {
                    hit.score < current.score
                };
                if better {
                    best = Some(hit);
                }
            }
        }
    }
    best
}

/// The shared lookups a `PSM` row needs.
struct PsmContext<'s, 'a> {
    prot_ids: &'s [&'a ProteinIdentification],
    run_index: &'s BTreeMap<String, usize>,
    run_file_to_ms_file: &'s BTreeMap<(usize, usize), usize>,
    run_to_engines: &'s BTreeMap<usize, Vec<(String, String)>>,
    db: &'s MzTabString,
    db_version: &'s MzTabString,
    export_empty_pep_ids: bool,
    export_all_psms: bool,
}

fn merge_index(identification: &PeptideIdentification) -> Result<Option<usize>> {
    match identification.metadata.get(ID_MERGE_INDEX) {
        None => Ok(None),
        Some(value) => usize::try_from(value.as_i64()?)
            .map(Some)
            .map_err(|_| bad("id_merge_index must be a nonnegative integer")),
    }
}

/// Source `PSMSectionRowFromPeptideID_`.
///
/// `None` for an identification without hits unless empty identifications are
/// exported; such a row carries the spectrum information only.
///
/// # Errors
///
/// [`Error::MissingInformation`] when the identification's run identifier is
/// unknown (the source's `std::map::at` throws `std::out_of_range`), when a run
/// with several files has no `id_merge_index` (the source's message), or when
/// the run has no registered search engine (undefined behaviour in the
/// source, which dereferences an empty vector's `begin()`).
fn psm_row_from_peptide_id(
    identification: &PeptideIdentification,
    context: &PsmContext<'_, '_>,
    current_psm_index: usize,
    psm_id: usize,
) -> Result<Option<MzTabPSMSectionRow>> {
    if identification.hits.is_empty() && !context.export_empty_pep_ids {
        return Ok(None);
    }
    let mut row = MzTabPSMSectionRow {
        psm_id: int_cell(
            i64::try_from(psm_id).map_err(|_| bad("PSM_ID overflows"))?,
            "PSM_ID",
        )?,
        database: context.db.clone(),
        database_version: context.db_version.clone(),
        retention_time: one_value_list(identification.rt.unwrap_or(f64::NAN)),
        exp_mass_to_charge: MzTabDouble::new(identification.mz.unwrap_or(f64::NAN)),
        ..MzTabPSMSectionRow::default()
    };
    let mut id_keys: BTreeSet<String> = identification.metadata.keys().cloned().collect();
    id_keys.remove(SIGNIFICANCE_THRESHOLD);
    add_meta_info_to_optional_columns(&id_keys, &mut row.opt, "global", &identification.metadata)?;

    let native_id = identification.spectrum_reference();
    let run_index = *context
        .run_index
        .get(&identification.identifier)
        .ok_or_else(|| {
            missing(format!(
                "peptide identification references unknown run {:?}",
                identification.identifier
            ))
        })?;
    let run = context.prot_ids[run_index];
    let localization_mods: Vec<String> = match run
        .search_parameters
        .metadata
        .get(LOCALIZED_MODIFICATIONS_USERPARAM)
    {
        Some(value) => value.as_string_list()?.to_vec(),
        None => Vec::new(),
    };
    let ms_file = if run.primary_ms_run_paths.len() <= 1 {
        context
            .run_file_to_ms_file
            .get(&(run_index, 0))
            .copied()
            .unwrap_or(0)
    } else {
        match merge_index(identification)? {
            Some(file) => context
                .run_file_to_ms_file
                .get(&(run_index, file))
                .copied()
                .unwrap_or(0),
            None => {
                return Err(missing(
                    "Multiple files in a run, but no id_merge_index in PeptideIdentification found.",
                ));
            }
        }
    };
    set_ms_file_lenient(&mut row.spectra_ref, ms_file);
    set_spec_ref_lenient(&mut row.spectra_ref, &native_id);
    if identification.hits.is_empty() {
        return Ok(Some(row));
    }

    let hit = if context.export_all_psms {
        identification
            .hits
            .get(current_psm_index)
            .ok_or_else(|| bad("PSM index outside the identification's hits"))?
    } else {
        best_hit(identification).ok_or_else(|| missing("identification has no best hit"))?
    };
    let sequence = &hit.sequence;
    row.sequence = string(sequence.as_str());
    row.modifications = MzTab::extract_modification_list(hit, &[], &localization_mods)?;
    let (engine, version) = context
        .run_to_engines
        .get(&run_index)
        .and_then(|engines| engines.first())
        .ok_or_else(|| missing("the identification's run has no registered search engine"))?;
    row.search_engine
        .from_cell_string(&format!("[,,{engine},{version}]"))?;
    row.search_engine_score
        .insert(1, MzTabDouble::new(hit.score));
    row.charge = MzTabInteger::new(hit.charge);
    row.calc_mass_to_charge = if hit.charge != 0 {
        MzTabDouble::new(sequence.mz(hit.charge)?)
    } else {
        MzTabDouble::null()
    };
    row.opt.push(MzTabOptionalColumnEntry::new(
        PEPTIDOFORM_COLUMN,
        string(&sequence.to_string()),
    ));
    let hit_keys: BTreeSet<String> = hit.metadata.keys().cloned().collect();
    add_meta_info_to_optional_columns(&hit_keys, &mut row.opt, "global", &hit.metadata)?;
    row.unique = MzTabBoolean::new(hit.protein_accessions().len() == 1);
    row.add_pep_evidence_to_rows(&hit.evidences)?;
    remap_target_decoy(&mut row.opt, DECOY_PEPTIDE_COLUMN);
    Ok(Some(row))
}

/// The source's `nextPSMRow` advance: past the last hit (or always, when only
/// the best PSM is exported) move to the next identification.
fn advance_psm(
    export_all_psms: bool,
    hits: usize,
    identification: &mut usize,
    current_psm_index: &mut usize,
) {
    if !export_all_psms || *current_psm_index + 1 >= hits {
        *identification += 1;
        *current_psm_index = 0;
    } else {
        *current_psm_index += 1;
    }
}

// ---------------------------------------------------------------------------
// IDMzTabStream
// ---------------------------------------------------------------------------

/// Row-at-a-time MzTab export of identification runs. Source
/// `MzTab::IDMzTabStream`.
///
/// Construction computes the complete metadata section and the optional column
/// names; [`IdMzTabStream::next_prt_row`] and [`IdMzTabStream::next_psm_row`]
/// then produce the rows one at a time, which is how
/// [`MzTabFile::store_identifications`](crate::format::mztab_file::MzTabFile::store_identifications)
/// writes without holding the document. There is no peptide section: the
/// stream has no feature information, and [`IdMzTabStream::next_pep_row`]
/// always ends immediately.
///
/// The stream borrows the runs and identifications; the source holds raw
/// pointers to them.
#[derive(Debug)]
pub struct IdMzTabStream<'a> {
    prot_ids: Vec<&'a ProteinIdentification>,
    peptide_ids: Vec<&'a PeptideIdentification>,
    protein_hit_user_value_keys: BTreeSet<String>,
    run_index: BTreeMap<String, usize>,
    run_to_engines: BTreeMap<usize, Vec<(String, String)>>,
    run_file_to_ms_file: BTreeMap<(usize, usize), usize>,
    first_run_inference: bool,
    export_empty_pep_ids: bool,
    export_all_psms: bool,
    quant_study_variables: usize,
    protein_cursor: ProteinCursorDebug,
    pep_id: usize,
    psm_id: usize,
    current_psm_index: usize,
    db: MzTabString,
    db_version: MzTabString,
    prt_optional_column_names: Vec<String>,
    pep_optional_column_names: Vec<String>,
    psm_optional_column_names: Vec<String>,
    meta_data: MzTabMetaData,
}

/// [`ProteinCursor`] wrapped so the public streams can derive `Debug` without
/// exposing it.
#[derive(Clone, Default)]
struct ProteinCursorDebug(ProteinCursor);

impl std::fmt::Debug for ProteinCursorDebug {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ProteinCursor")
            .field("state", &self.0.state)
            .field("run", &self.0.run)
            .finish_non_exhaustive()
    }
}

/// The `[MS, MS:1000752, TOPP software, <version>]` entry both streams write
/// first. The source writes `VersionInfo::getVersion()`; this port writes the
/// crate version. The upstream tests exclude this line from comparison.
fn topp_software() -> Result<MzTabSoftwareMetaData> {
    Ok(MzTabSoftwareMetaData {
        software: parameter(&format!(
            "[MS,MS:1000752,TOPP software,{}]",
            env!("CARGO_PKG_VERSION")
        ))?,
        ..MzTabSoftwareMetaData::default()
    })
}

/// The `database` and `database_version` cells of every row: the database
/// file's base name without extension, and its version or null.
fn database_cells(run: &ProteinIdentification) -> (MzTabString, MzTabString) {
    let parameters = &run.search_parameters;
    let db = string(strip_extension(basename(&parameters.database)));
    let version = if parameters.database_version.is_empty() {
        MzTabString::null()
    } else {
        string(&parameters.database_version)
    };
    (db, version)
}

impl<'a> IdMzTabStream<'a> {
    /// Prepare the export of `prot_ids` and `peptide_ids`.
    ///
    /// `filename` is accepted for parity with the source, which stores it and
    /// never reads it. The metadata section is complete when this returns.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the input exceeds
    /// [`MzTab::MAX_ROWS`] identifications or runs, or when a run's inference
    /// or search-setting metadata has the wrong type;
    /// [`Error::MissingInformation`] when a run's MS file cannot be registered;
    /// [`Error::Parse`] when an engine name or version turns a metadata term
    /// into unparseable text.
    pub fn new(
        prot_ids: &[&'a ProteinIdentification],
        peptide_ids: &[&'a PeptideIdentification],
        filename: &str,
        options: &IdentificationExportOptions,
    ) -> Result<Self> {
        let _ = filename;
        if prot_ids.len() > MzTab::MAX_ROWS || peptide_ids.len() > MzTab::MAX_ROWS {
            return Err(bad("MzTab export input exceeds its row limit"));
        }
        let prot_ids: Vec<&'a ProteinIdentification> = prot_ids.to_vec();
        let peptide_ids: Vec<&'a PeptideIdentification> = peptide_ids.to_vec();
        let run_index = map_run_identifier_to_index(&prot_ids);
        let has_inference = match prot_ids.first() {
            Some(run) => run.has_inference_data()?,
            None => false,
        };
        let first_run_inference = has_inference && options.first_run_inference_only;
        let (name_to_index, index_to_name) = map_ms_file_names(&prot_ids, first_run_inference);
        let run_file_to_ms_file =
            map_run_file_to_ms_file(&prot_ids, &name_to_index, first_run_inference)?;
        let (variable_mods, fixed_mods) = search_modifications(&prot_ids);
        let engines = map_runs_and_search_engines(&prot_ids, &peptide_ids, first_run_inference)?;
        let (mut protein_hit_keys, peptide_id_keys, peptide_hit_keys) =
            identification_meta_keys(&prot_ids, &peptide_ids);
        protein_hit_keys.remove("Description");

        let mut prt_names: Vec<String> = protein_hit_keys
            .iter()
            .map(|key| format!("opt_global_{key}"))
            .collect();
        let mut psm_names: Vec<String> = peptide_id_keys
            .iter()
            .chain(&peptide_hit_keys)
            .map(|key| format!("opt_global_{key}"))
            .collect();
        for name in &mut prt_names {
            if name == TARGET_DECOY_COLUMN {
                *name = DECOY_HIT_COLUMN.into();
            }
        }
        prt_names.push(RESULT_TYPE_COLUMN.into());
        for name in &mut psm_names {
            if name == TARGET_DECOY_COLUMN {
                *name = DECOY_PEPTIDE_COLUMN.into();
            }
        }
        psm_names.push(PEPTIDOFORM_COLUMN.into());

        let quant_study_variables = prot_ids.first().map_or(0, |run| quant_study_variables(run));

        let mut meta = MzTabMetaData {
            mz_tab_type: string("Identification"),
            mz_tab_mode: string("Summary"),
            description: string("OpenMS export from ID data"),
            title: string(&options.title),
            variable_mod: modifications_metadata(&variable_mods)?,
            fixed_mod: modifications_metadata(&fixed_mods)?,
            ..MzTabMetaData::default()
        };
        let software_index = (meta.software.len() + 1).max(1);
        meta.software.insert(software_index, topp_software()?);
        let (mut db, mut db_version) = (MzTabString::null(), MzTabString::null());
        if let Some(first) = prot_ids.first() {
            add_ms_run_meta_data(&index_to_name, &mut meta);
            add_search_meta_data(&prot_ids, &engines, &mut meta, first_run_inference)?;
            (db, db_version) = database_cells(first);
        }

        Ok(Self {
            prot_ids,
            peptide_ids,
            protein_hit_user_value_keys: protein_hit_keys,
            run_index,
            run_to_engines: engines.run_to_engines,
            run_file_to_ms_file,
            first_run_inference,
            export_empty_pep_ids: options.export_empty_pep_ids,
            export_all_psms: options.export_all_psms,
            quant_study_variables,
            protein_cursor: ProteinCursorDebug::default(),
            pep_id: 0,
            psm_id: 0,
            current_psm_index: 0,
            db,
            db_version,
            prt_optional_column_names: prt_names,
            pep_optional_column_names: Vec::new(),
            psm_optional_column_names: psm_names,
            meta_data: meta,
        })
    }

    /// The complete metadata section.
    pub fn meta_data(&self) -> &MzTabMetaData {
        &self.meta_data
    }

    /// The `PRT` optional column names: `opt_global_<key>` for every protein
    /// hit metadata key except `Description` (spaces as underscores, the
    /// target/decoy key renamed), then `opt_global_result_type`.
    pub fn protein_optional_column_names(&self) -> &[String] {
        &self.prt_optional_column_names
    }

    /// The `PEP` optional column names, always empty for this stream.
    pub fn peptide_optional_column_names(&self) -> &[String] {
        &self.pep_optional_column_names
    }

    /// The `PSM` optional column names: the identification keys (without
    /// `significance_threshold`), then the hit keys — a key present at both
    /// levels appears twice — with the target/decoy key renamed, then the
    /// peptidoform column.
    pub fn psm_optional_column_names(&self) -> &[String] {
        &self.psm_optional_column_names
    }

    /// The next `PRT` row, or `None` when the section is complete.
    ///
    /// Source `nextPRTRow`, which fills an out-parameter and returns `false` at
    /// the end.
    ///
    /// # Errors
    ///
    /// As the protein row builders: an indistinguishable group that names no
    /// protein hit, or protein metadata of the wrong type.
    pub fn next_prt_row(&mut self) -> Result<Option<MzTabProteinSectionRow>> {
        let context = ProteinContext {
            prot_ids: &self.prot_ids,
            first_run_inference: self.first_run_inference,
            quant_study_variables: self.quant_study_variables,
            db: &self.db,
            db_version: &self.db_version,
            user_keys: &self.protein_hit_user_value_keys,
            assay_quantification: None,
        };
        self.protein_cursor.0.next(&context)
    }

    /// Always `None`: identification data has no peptide section. Source
    /// `nextPEPRow`, which returns `false`.
    pub fn next_pep_row(&mut self) -> Option<MzTabPeptideSectionRow> {
        None
    }

    /// The next `PSM` row, or `None` when every identification was visited.
    ///
    /// An identification that yields no row (no hits, and empty
    /// identifications not exported) produces a *default* row whose sequence
    /// is null, exactly as the source's `nextPSMRow`; writers skip such rows.
    /// The `PSM_ID` counter advances for every call that returns a row,
    /// placeholder or not. With `export_all_psms` an identification without
    /// hits still advances, the fix the source's regression test pins.
    ///
    /// # Errors
    ///
    /// As `psm_row_from_peptide_id`: an unknown run identifier, a multi-file
    /// run without `id_merge_index`, a run without search engine, metadata of
    /// the wrong type, or sequence chemistry that has no mass.
    pub fn next_psm_row(&mut self) -> Result<Option<MzTabPSMSectionRow>> {
        let Some(identification) = self.peptide_ids.get(self.pep_id).copied() else {
            return Ok(None);
        };
        let context = PsmContext {
            prot_ids: &self.prot_ids,
            run_index: &self.run_index,
            run_file_to_ms_file: &self.run_file_to_ms_file,
            run_to_engines: &self.run_to_engines,
            db: &self.db,
            db_version: &self.db_version,
            export_empty_pep_ids: self.export_empty_pep_ids,
            export_all_psms: self.export_all_psms,
        };
        let row = psm_row_from_peptide_id(
            identification,
            &context,
            self.current_psm_index,
            self.psm_id,
        )?;
        advance_psm(
            self.export_all_psms,
            identification.hits.len(),
            &mut self.pep_id,
            &mut self.current_psm_index,
        );
        self.psm_id += 1;
        Ok(Some(row.unwrap_or_default()))
    }
}

// ---------------------------------------------------------------------------
// Feature map export
// ---------------------------------------------------------------------------

/// Source `getFeatureMapMetaValues_`: every feature, identification and hit
/// key, without `spectrum_reference` at hit level and
/// `significance_threshold` at identification level. Keys keep their spaces.
fn feature_map_meta_keys(
    map: &FeatureMap,
) -> (BTreeSet<String>, BTreeSet<String>, BTreeSet<String>) {
    let mut feature = BTreeSet::new();
    let mut identification = BTreeSet::new();
    let mut hit = BTreeSet::new();
    for f in &map.features {
        feature.extend(f.metadata.keys().cloned());
        for id in &f.peptide_identifications {
            identification.extend(id.metadata.keys().cloned());
            for h in &id.hits {
                hit.extend(h.metadata.keys().cloned());
            }
        }
    }
    hit.remove("spectrum_reference");
    identification.remove(SIGNIFICANCE_THRESHOLD);
    (feature, identification, hit)
}

/// Source `peptideSectionRowFromFeature_`.
fn peptide_row_from_feature(
    feature: &Feature,
    feature_keys: &BTreeSet<String>,
    identification_keys: &BTreeSet<String>,
    hit_keys: &BTreeSet<String>,
    fixed_mods: &[String],
) -> Result<MzTabPeptideSectionRow> {
    let mut row = MzTabPeptideSectionRow {
        mass_to_charge: MzTabDouble::new(feature.mz),
        retention_time: one_value_list(feature.rt),
        charge: MzTabInteger::new(feature.charge),
        ..MzTabPeptideSectionRow::default()
    };
    let mut window = Vec::new();
    if let Some(bounds) = feature.convex_hull().bounding_box() {
        window.push(MzTabDouble::new(bounds.rt_range().min));
        window.push(MzTabDouble::new(bounds.rt_range().max));
    }
    row.retention_time_window.set(window);
    row.peptide_abundance_stdev_study_variable
        .insert(1, MzTabDouble::null());
    row.peptide_abundance_std_error_study_variable
        .insert(1, MzTabDouble::null());
    row.peptide_abundance_study_variable
        .insert(1, MzTabDouble::new(f64::from(feature.intensity)));
    row.best_search_engine_score.insert(1, MzTabDouble::null());
    row.search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(1, MzTabDouble::null());
    row.opt.push(MzTabOptionalColumnEntry::new(
        PEPTIDOFORM_COLUMN,
        MzTabString::null(),
    ));
    add_meta_info_to_optional_columns(feature_keys, &mut row.opt, "global", &feature.metadata)?;
    let empty = MetaInfo::new();
    let ids = &feature.peptide_identifications;
    let Some(first) = ids.first() else {
        add_meta_info_to_optional_columns(identification_keys, &mut row.opt, "global", &empty)?;
        add_meta_info_to_optional_columns(hit_keys, &mut row.opt, "global", &empty)?;
        return Ok(row);
    };
    add_meta_info_to_optional_columns(
        identification_keys,
        &mut row.opt,
        "global",
        &first.metadata,
    )?;
    let all_hits: Vec<PeptideHit> = ids.iter().flat_map(|id| id.hits.iter().cloned()).collect();
    if all_hits.is_empty() {
        // The source appends the identification columns a second time here,
        // not the hit columns; preserved.
        add_meta_info_to_optional_columns(identification_keys, &mut row.opt, "global", &empty)?;
        return Ok(row);
    }
    let mut combined = first.clone();
    combined.hits = all_hits;
    let higher = combined.higher_score_better;
    // PeptideIdentification::sort: a stable sort by score, NaN-tolerant.
    combined.hits.sort_by(|a, b| {
        let order = a
            .score
            .partial_cmp(&b.score)
            .unwrap_or(std::cmp::Ordering::Equal);
        if higher { order.reverse() } else { order }
    });
    let best = &combined.hits[0];
    let sequence = &best.sequence;
    row.sequence = string(sequence.as_str());
    row.modifications = MzTab::extract_modification_list(best, fixed_mods, &[])?;
    row.unique = MzTabBoolean::new(best.protein_accessions().len() == 1);
    row.accession = best
        .evidences
        .first()
        .map_or_else(MzTabString::null, |e| string(&e.protein_accession));
    row.best_search_engine_score
        .insert(1, MzTabDouble::new(best.score));
    row.search_engine_score_ms_run
        .entry(1)
        .or_default()
        .insert(1, MzTabDouble::new(best.score));
    set_entry(&mut row.opt, PEPTIDOFORM_COLUMN, &sequence.to_string());
    add_meta_info_to_optional_columns(hit_keys, &mut row.opt, "global", &best.metadata)?;
    remap_target_decoy(&mut row.opt, DECOY_PEPTIDE_COLUMN);
    Ok(row)
}

impl MzTab {
    /// Export a feature map as a `Quantification`/`Summary` document with one
    /// `PEP` row per feature. Source `MzTab::exportFeatureMapToMzTab`.
    ///
    /// Modifications and the database come from the first protein
    /// identification run's search parameters; the variable and fixed
    /// modification metadata use the placeholder terms when the lists are
    /// empty. The single `ms_run[1]` location is the map's first primary MS run
    /// path (the `UNKNOWN` placeholder when it has none) with `file://`
    /// prepended, and `uri[1]` is `filename`.
    ///
    /// A row's retention time window is the RT extent of the feature's convex
    /// hull, empty when it has none. Its identification is the best hit over
    /// all of the feature's identifications, taking the score orientation of
    /// the first one.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the map exceeds [`MzTab::MAX_ROWS`] features
    /// or a row its optional-column limit, or when the map's `spectra_data` is
    /// not a string list; errors of [`MzTab::extract_modification_list`].
    pub fn export_feature_map(feature_map: &FeatureMap, filename: &str) -> Result<MzTab> {
        if feature_map.features.len() > MzTab::MAX_ROWS {
            return Err(bad("MzTab export input exceeds its row limit"));
        }
        let mut document = MzTab::default();
        let mut meta = MzTabMetaData::default();
        let (mut variable_mods, mut fixed_mods) = (Vec::new(), Vec::new());
        if let Some(first) = feature_map.protein_identifications.first() {
            variable_mods = first.search_parameters.variable_modifications.clone();
            fixed_mods = first.search_parameters.fixed_modifications.clone();
        }
        meta.variable_mod =
            crate::format::mztab::variable_modification_metadata(&variable_mods)?.metadata;
        meta.fixed_mod = crate::format::mztab::fixed_modification_metadata(&fixed_mods)?.metadata;
        meta.mz_tab_type = string("Quantification");
        meta.mz_tab_mode = string("Summary");
        meta.description = string("OpenMS export from featureXML");
        let paths = feature_map.primary_ms_run_path()?;
        let run = MzTabMSRunMetaData {
            location: paths
                .first()
                .map_or_else(MzTabString::null, |path| file_location(path)),
            ..MzTabMSRunMetaData::default()
        };
        meta.ms_run.insert(1, run);
        meta.uri.insert(1, string(filename));
        meta.psm_search_engine_score
            .insert(1, MzTabParameter::default());
        meta.peptide_search_engine_score
            .insert(1, MzTabParameter::default());
        document.meta_data = meta;

        let (feature_keys, identification_keys, hit_keys) = feature_map_meta_keys(feature_map);
        let mut rows = Vec::with_capacity(feature_map.features.len());
        for feature in &feature_map.features {
            rows.push(peptide_row_from_feature(
                feature,
                &feature_keys,
                &identification_keys,
                &hit_keys,
                &fixed_mods,
            )?);
        }
        document.peptide_data = rows;
        Ok(document)
    }

    /// Export identification runs and peptide identifications as an
    /// `Identification`/`Summary` document with `PRT` and `PSM` sections.
    /// Source `MzTab::exportIdentificationsToMzTab`, which drains an
    /// [`IdMzTabStream`].
    ///
    /// Every row the stream yields is kept, including the placeholder rows with
    /// a null sequence that the stream emits for skipped identifications — the
    /// source keeps them here, unlike its file writer.
    ///
    /// # Errors
    ///
    /// As [`IdMzTabStream::new`], [`IdMzTabStream::next_prt_row`] and
    /// [`IdMzTabStream::next_psm_row`], plus [`Error::InvalidValue`] when a
    /// section would exceed [`MzTab::MAX_ROWS`].
    pub fn export_identifications(
        prot_ids: &[ProteinIdentification],
        peptide_ids: &[PeptideIdentification],
        filename: &str,
        options: &IdentificationExportOptions,
    ) -> Result<MzTab> {
        let protein_refs: Vec<&ProteinIdentification> = prot_ids.iter().collect();
        let peptide_refs: Vec<&PeptideIdentification> = peptide_ids.iter().collect();
        let mut stream = IdMzTabStream::new(&protein_refs, &peptide_refs, filename, options)?;
        let mut document = MzTab {
            meta_data: stream.meta_data().clone(),
            ..MzTab::default()
        };
        while let Some(row) = stream.next_prt_row()? {
            push_bounded(&mut document.protein_data, row)?;
        }
        while let Some(row) = stream.next_psm_row()? {
            push_bounded(&mut document.psm_data, row)?;
        }
        Ok(document)
    }
}

fn push_bounded<T>(rows: &mut Vec<T>, row: T) -> Result<()> {
    if rows.len() >= MzTab::MAX_ROWS {
        return Err(bad("MzTab section exceeds its row limit"));
    }
    rows.push(row);
    Ok(())
}

// ---------------------------------------------------------------------------
// Consensus map export
// ---------------------------------------------------------------------------

/// Source `getConsensusMapMetaValues_` (with its local helper
/// `extractMetaValuesFromIDs`): consensus-feature, identification and hit keys
/// with spaces as underscores, from the unassigned identifications and from
/// every feature; `spectrum_reference` and `significance_threshold` removed at
/// identification level.
fn consensus_map_meta_keys(
    map: &ConsensusMap,
) -> (BTreeSet<String>, BTreeSet<String>, BTreeSet<String>) {
    let mut feature = BTreeSet::new();
    let mut identification = BTreeSet::new();
    let mut hit = BTreeSet::new();
    let mut from_ids = |ids: &[PeptideIdentification]| {
        for id in ids {
            identification.extend(id.metadata.keys().map(|key| underscored(key)));
            for h in &id.hits {
                hit.extend(h.metadata.keys().map(|key| underscored(key)));
            }
        }
    };
    from_ids(&map.unassigned_peptide_identifications);
    for f in &map.features {
        feature.extend(f.metadata.keys().map(|key| underscored(key)));
        from_ids(&f.peptide_identifications);
    }
    identification.remove("spectrum_reference");
    identification.remove(SIGNIFICANCE_THRESHOLD);
    (feature, identification, hit)
}

/// Source `checkSequenceUniqueness_`.
///
/// # Errors
///
/// [`Error::InvalidValue`] with the source's message when two identifications
/// of one consensus feature disagree on their first hit's sequence, and when
/// an identification has no hit (undefined behaviour in the source, which
/// indexes the empty hit list).
fn check_sequence_uniqueness(ids: &[PeptideIdentification]) -> Result<()> {
    let first = first_hit(&ids[0])?;
    for id in ids {
        if first_hit(id)?.sequence != first.sequence {
            return Err(bad(
                "Consensus features may contain at most one identification. Run IDConflictResolver first to remove ambiguities!",
            ));
        }
    }
    Ok(())
}

fn first_hit(identification: &PeptideIdentification) -> Result<&PeptideHit> {
    identification.hits.first().ok_or_else(|| {
        bad("a consensus feature's peptide identification has no hit; remove empty identifications first")
    })
}

/// The design-derived assay structure of a consensus-map export.
#[derive(Clone, Debug, Default)]
struct AssayDesign {
    /// `(fraction group, label)` of assay `i + 1`, ascending.
    keys: Vec<(u32, u32)>,
    /// One-based study variable of each one-based assay; index 0 unused.
    assay_to_study_variable: Vec<usize>,
    /// Assays of each one-based study variable; index 0 unused.
    study_variable_to_assays: Vec<Vec<usize>>,
    /// `(path, label)` to the zero-based assay index.
    path_label_to_assay: BTreeMap<(String, u32), usize>,
    study_variables: usize,
}

impl AssayDesign {
    /// The source's assay construction from `ExperimentalDesign::fromConsensusMap`:
    /// one assay per distinct `(fraction group, label)` quantification unit, in
    /// ascending key order; one study variable per sample.
    fn from_map(map: &ConsensusMap) -> Result<Self> {
        let design = ExperimentalDesign::from_consensus_map(map)?;
        let mut key_to_sample: BTreeMap<(u32, u32), u32> = BTreeMap::new();
        for entry in design.ms_file_section() {
            let key = (entry.fraction_group, entry.label);
            match key_to_sample.get(&key) {
                Some(sample) if *sample != entry.sample => {
                    return Err(missing(
                        "The same (fraction_group, label) quantification unit resolves to several samples.",
                    ));
                }
                Some(_) => {}
                None => {
                    key_to_sample.insert(key, entry.sample);
                }
            }
        }
        let keys: Vec<(u32, u32)> = key_to_sample.keys().copied().collect();
        let study_variables = design.number_of_samples() as usize;
        let mut assay_to_study_variable = vec![0; keys.len() + 1];
        let mut study_variable_to_assays = vec![Vec::new(); study_variables + 1];
        let mut index_by_key = BTreeMap::new();
        for (offset, key) in keys.iter().enumerate() {
            let assay = offset + 1;
            let study_variable = key_to_sample[key] as usize + 1;
            if study_variable > study_variables {
                return Err(missing(
                    "A quantification unit resolves to a sample outside the experimental design's sample section.",
                ));
            }
            index_by_key.insert(*key, offset);
            assay_to_study_variable[assay] = study_variable;
            study_variable_to_assays[study_variable].push(assay);
        }
        let mut path_label_to_assay = BTreeMap::new();
        for entry in design.ms_file_section() {
            path_label_to_assay.insert(
                (entry.path.clone(), entry.label),
                index_by_key[&(entry.fraction_group, entry.label)],
            );
        }
        Ok(Self {
            keys,
            assay_to_study_variable,
            study_variable_to_assays,
            path_label_to_assay,
            study_variables,
        })
    }

    fn assay_of(&self, filename: &str, label: u32) -> Result<usize> {
        self.path_label_to_assay
            .get(&(filename.to_owned(), label))
            .map(|index| index + 1)
            .ok_or_else(|| {
                missing(format!(
                    "The consensus map column ('{filename}', label {label}) has no quantification unit in the experimental design."
                ))
            })
    }
}

/// The lookups a consensus `PEP` row needs.
struct ConsensusRowContext<'s, 'a> {
    map: &'a ConsensusMap,
    ms_runs: usize,
    design: &'s AssayDesign,
    feature_keys: &'s BTreeSet<String>,
    identification_keys: &'s BTreeSet<String>,
    hit_keys: &'s BTreeSet<String>,
    run_index: &'s BTreeMap<String, usize>,
    run_file_to_ms_file: &'s BTreeMap<(usize, usize), usize>,
    fixed_mods: &'s [String],
    export_subfeatures: bool,
}

/// Source `peptideSectionRowFromConsensusFeature_`.
///
/// # Errors
///
/// [`Error::MissingInformation`] when a handle's map index has no column
/// header, a column has no quantification unit, an identification's run is
/// unknown, a multi-file run lacks `id_merge_index`, or a PSM maps to an MS
/// file the quantitative metadata does not register — each with the source's
/// message where it has one; [`Error::InvalidValue`] from
/// [`check_sequence_uniqueness`].
fn peptide_row_from_consensus_feature(
    feature: &ConsensusFeature,
    context: &ConsensusRowContext<'_, '_>,
) -> Result<MzTabPeptideSectionRow> {
    let map = context.map;
    let design = context.design;
    let mut row = MzTabPeptideSectionRow::default();
    row.opt.push(MzTabOptionalColumnEntry::new(
        PEPTIDOFORM_COLUMN,
        MzTabString::null(),
    ));
    push_entries(&mut row.opt, context.feature_keys, |key| {
        feature
            .metadata
            .get(key)
            .map_or_else(MzTabString::null, |value| string(&meta_value_text(value)))
    });
    // The source reads these from the consensus feature, not from its
    // identifications; an absent key yields an empty, i.e. null, cell.
    push_entries(&mut row.opt, context.identification_keys, |key| {
        feature
            .metadata
            .get(key)
            .map_or_else(MzTabString::null, |value| string(&meta_value_text(value)))
    });
    push_entries(&mut row.opt, context.hit_keys, |_| MzTabString::null());

    row.mass_to_charge = MzTabDouble::new(feature.mz);
    row.retention_time = one_value_list(feature.rt);
    row.charge = MzTabInteger::new(feature.charge);
    row.best_search_engine_score.insert(1, MzTabDouble::null());

    let distinct_assay_grain = design
        .study_variable_to_assays
        .iter()
        .skip(1)
        .any(|assays| assays.len() != 1);
    if distinct_assay_grain {
        for assay in 1..=design.keys.len() {
            row.peptide_abundance_assay
                .insert(assay, MzTabDouble::null());
        }
    }
    for study_variable in 1..=design.study_variables {
        row.peptide_abundance_stdev_study_variable
            .insert(study_variable, MzTabDouble::null());
        row.peptide_abundance_std_error_study_variable
            .insert(study_variable, MzTabDouble::null());
        row.peptide_abundance_study_variable
            .insert(study_variable, MzTabDouble::null());
    }
    let scores = row.search_engine_score_ms_run.entry(1).or_default();
    for ms_run in 1..=context.ms_runs {
        scores.insert(ms_run, MzTabDouble::null());
    }

    let mut assay_abundances: BTreeMap<usize, f64> = BTreeMap::new();
    for handle in feature.handles() {
        let header = map.column_headers.get(&handle.map_index).ok_or_else(|| {
            missing(format!(
                "consensus feature handle references map index {} without a column header",
                handle.map_index
            ))
        })?;
        let label = header.label_as_uint(&map.experiment_type)?;
        let assay = design.assay_of(&header.filename, label)?;
        let study_variable = design.assay_to_study_variable[assay];
        *assay_abundances.entry(assay).or_insert(0.0) += f64::from(handle.intensity);
        if context.export_subfeatures {
            let suffix = if design.study_variable_to_assays[study_variable].len() == 1 {
                format!("study_variable[{study_variable}]")
            } else {
                format!("assay[{assay}]")
            };
            row.opt.push(MzTabOptionalColumnEntry::new(
                format!("opt_global_mass_to_charge_{suffix}"),
                string(&format_float(handle.mz, true)),
            ));
            row.opt.push(MzTabOptionalColumnEntry::new(
                format!("opt_global_retention_time_{suffix}"),
                string(&format_float(handle.rt, true)),
            ));
        }
    }
    if distinct_assay_grain {
        for (assay, abundance) in &assay_abundances {
            row.peptide_abundance_assay
                .insert(*assay, MzTabDouble::new(*abundance));
        }
    }
    for study_variable in 1..=design.study_variables {
        let assays = &design.study_variable_to_assays[study_variable];
        if assays.len() == 1 {
            if let Some(abundance) = assay_abundances.get(&assays[0]) {
                row.peptide_abundance_study_variable
                    .insert(study_variable, MzTabDouble::new(*abundance));
            }
        }
    }

    let ids = &feature.peptide_identifications;
    if !ids.is_empty() {
        check_sequence_uniqueness(ids)?;
        let best_id = &ids[0];
        let best = first_hit(best_id)?;
        let sequence = &best.sequence;
        row.sequence = string(sequence.as_str());
        row.modifications = MzTab::extract_modification_list(best, context.fixed_mods, &[])?;
        row.unique = MzTabBoolean::new(best.protein_accessions().len() == 1);
        row.accession = best
            .evidences
            .first()
            .map_or_else(MzTabString::null, |e| string(&e.protein_accession));
        overwrite_entries(&mut row.opt, &best_id.metadata);
        set_entry(&mut row.opt, PEPTIDOFORM_COLUMN, &sequence.to_string());
        overwrite_entries(&mut row.opt, &best.metadata);

        let mut best_score = best.score;
        for id in ids {
            let run_index = *context.run_index.get(&id.identifier).ok_or_else(|| {
                missing(format!(
                    "peptide identification references unknown run {:?}",
                    id.identifier
                ))
            })?;
            let run = map
                .protein_identifications
                .get(run_index)
                .ok_or_else(|| missing("identification run index outside the map's runs"))?;
            let key = if run.primary_ms_run_paths.len() <= 1 {
                (run_index, 0)
            } else {
                match merge_index(id)? {
                    Some(file) => (run_index, file),
                    None => {
                        return Err(missing(
                            "Multiple files in a run, but no id_merge_index in PeptideIdentification found.",
                        ));
                    }
                }
            };
            let ms_file = *context.run_file_to_ms_file.get(&key).ok_or_else(|| {
                missing("identification run and file index have no registered MS file")
            })?;
            let current = first_hit(id)?.score;
            let cell = row
                .search_engine_score_ms_run
                .entry(1)
                .or_default()
                .get_mut(&ms_file)
                .ok_or_else(|| {
                    missing(format!(
                        "PSM {} does not map to an MS file registered in the quantitative metadata. Check your merging and filtering steps and/or report the issue, please.",
                        id.metadata
                            .get("spectrum_reference")
                            .map(meta_value_text)
                            .unwrap_or_default()
                    ))
                })?;
            *cell = MzTabDouble::new(current);
            let improves = if id.higher_score_better {
                current >= best_score
            } else {
                current <= best_score
            };
            if improves {
                best_score = current;
                if id.metadata.contains_key("spectrum_reference") {
                    set_spec_ref_lenient(&mut row.spectra_ref, &id.spectrum_reference());
                    set_ms_file_lenient(&mut row.spectra_ref, ms_file);
                }
            }
        }
        row.best_search_engine_score
            .insert(1, MzTabDouble::new(best_score));
    }
    remap_target_decoy(&mut row.opt, DECOY_PEPTIDE_COLUMN);
    Ok(row)
}

/// Row-at-a-time MzTab export of a consensus map. Source
/// `MzTab::CMMzTabStream`.
///
/// Construction computes the complete metadata section — including the assays
/// and study variables derived from the map's experimental design — and the
/// optional column names; the three `next_*_row` methods then produce the
/// `PRT`, `PEP` and `PSM` rows one at a time. The stream borrows the map.
///
/// An OpenMS sample becomes an MzTab study variable and each of its
/// `(fraction group, label)` quantification units an assay; an assay lists
/// every `ms_run` (fraction) of its unit.
#[derive(Debug)]
pub struct ConsensusMapMzTabStream<'a> {
    map: &'a ConsensusMap,
    prot_ids: Vec<&'a ProteinIdentification>,
    peptide_ids: Vec<&'a PeptideIdentification>,
    protein_hit_user_value_keys: BTreeSet<String>,
    feature_keys: BTreeSet<String>,
    identification_keys: BTreeSet<String>,
    hit_keys: BTreeSet<String>,
    run_index: BTreeMap<String, usize>,
    run_to_engines: BTreeMap<usize, Vec<(String, String)>>,
    run_file_to_ms_file: BTreeMap<(usize, usize), usize>,
    design: AssayDesignDebug,
    ms_runs: usize,
    first_run_inference: bool,
    fixed_mods: Vec<String>,
    export_unidentified_features: bool,
    export_subfeatures: bool,
    export_empty_pep_ids: bool,
    export_all_psms: bool,
    quant_study_variables: usize,
    protein_cursor: ProteinCursorDebug,
    pep_id: usize,
    pep_counter: usize,
    psm_id: usize,
    current_psm_index: usize,
    db: MzTabString,
    db_version: MzTabString,
    prt_optional_column_names: Vec<String>,
    pep_optional_column_names: Vec<String>,
    psm_optional_column_names: Vec<String>,
    meta_data: MzTabMetaData,
}

/// [`AssayDesign`] behind an opaque `Debug`.
#[derive(Clone, Default)]
struct AssayDesignDebug(AssayDesign);

impl std::fmt::Debug for AssayDesignDebug {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AssayDesign")
            .field("assays", &self.0.keys.len())
            .field("study_variables", &self.0.study_variables)
            .finish()
    }
}

/// Keys present in every protein hit of `hits` (source
/// `MetaInfoInterfaceUtils::findCommonMetaKeys` at 100 %). No hits, no keys.
fn common_protein_hit_keys(hits: &[ProteinHit]) -> BTreeSet<String> {
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for hit in hits {
        for key in hit.metadata.keys() {
            *counts.entry(key.as_str()).or_insert(0) += 1;
        }
    }
    counts
        .into_iter()
        .filter(|(_, count)| *count >= hits.len())
        .map(|(key, _)| key.to_owned())
        .collect()
}

impl<'a> ConsensusMapMzTabStream<'a> {
    /// Prepare the export of `consensus_map`.
    ///
    /// The identifications exported as PSMs are those of every consensus
    /// feature, in feature order, then — with `export_unassigned_ids` — the
    /// unassigned ones. `filename` is accepted for parity with the source,
    /// which stores it and never reads it.
    ///
    /// The protein-hit optional columns are the metadata keys common to *all*
    /// hits of a run, united over the runs, without `Description`, with spaces
    /// as underscores. The `PEP` optional columns are the consensus-feature
    /// keys, the PRIDE decoy column when some hit carries `target_decoy`, and
    /// the peptidoform column; identification-level keys are deliberately not
    /// written there, as in the source.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] with the source's messages when the
    /// experimental design maps one quantification unit to several samples or
    /// to a sample outside its sample section, when a column header has no
    /// quantification unit, or when a column's file has no `ms_run` entry;
    /// otherwise as [`IdMzTabStream::new`].
    pub fn new(
        consensus_map: &'a ConsensusMap,
        filename: &str,
        options: &ConsensusMapExportOptions,
    ) -> Result<Self> {
        let _ = filename;
        let map = consensus_map;
        if map.features.len() > MzTab::MAX_ROWS
            || map.protein_identifications.len() > MzTab::MAX_ROWS
        {
            return Err(bad("MzTab export input exceeds its row limit"));
        }
        let prot_ids: Vec<&'a ProteinIdentification> = map.protein_identifications.iter().collect();
        let mut peptide_ids: Vec<&'a PeptideIdentification> = Vec::new();
        for feature in &map.features {
            peptide_ids.extend(feature.peptide_identifications.iter());
        }
        if options.export_unassigned_ids {
            peptide_ids.extend(map.unassigned_peptide_identifications.iter());
        }
        if peptide_ids.len() > MzTab::MAX_ROWS {
            return Err(bad("MzTab export input exceeds its row limit"));
        }
        let run_index = map_run_identifier_to_index(&prot_ids);
        let has_inference = match prot_ids.first() {
            Some(run) => run.has_inference_data()?,
            None => false,
        };
        let first_run_inference = has_inference && options.first_run_inference_only;
        let (name_to_index, index_to_name) = map_ms_file_names(&prot_ids, first_run_inference);
        let run_file_to_ms_file =
            map_run_file_to_ms_file(&prot_ids, &name_to_index, first_run_inference)?;
        let (variable_mods, fixed_mods) = search_modifications(&prot_ids);
        let id_format = ms_run_spectrum_identifier_type(&peptide_ids)?;
        let engines = map_runs_and_search_engines(&prot_ids, &peptide_ids, first_run_inference)?;
        let (feature_keys, identification_keys, hit_keys) = consensus_map_meta_keys(map);

        let mut pep_names: Vec<String> = feature_keys
            .iter()
            .map(|key| format!("opt_global_{key}"))
            .collect();
        if hit_keys.contains("target_decoy") {
            pep_names.push(DECOY_PEPTIDE_COLUMN.into());
        }
        pep_names.push(PEPTIDOFORM_COLUMN.into());
        let mut psm_names: Vec<String> = hit_keys
            .iter()
            .map(|key| format!("opt_global_{key}"))
            .collect();
        for name in &mut psm_names {
            if name == TARGET_DECOY_COLUMN {
                *name = DECOY_PEPTIDE_COLUMN.into();
            }
        }
        psm_names.push(PEPTIDOFORM_COLUMN.into());

        let quant_study_variables = prot_ids.first().map_or(0, |run| quant_study_variables(run));

        let mut meta = MzTabMetaData {
            mz_tab_type: string("Quantification"),
            mz_tab_mode: string("Summary"),
            description: string("OpenMS export from consensusXML"),
            variable_mod: modifications_metadata(&variable_mods)?,
            fixed_mod: modifications_metadata(&fixed_mods)?,
            ..MzTabMetaData::default()
        };
        let software_index = (meta.software.len() + 1).max(1);
        meta.software.insert(software_index, topp_software()?);

        let mut protein_keys = BTreeSet::new();
        let (mut db, mut db_version) = (MzTabString::null(), MzTabString::null());
        if let Some(first) = prot_ids.first() {
            add_ms_run_meta_data(&index_to_name, &mut meta);
            add_search_meta_data(&prot_ids, &engines, &mut meta, first_run_inference)?;
            (db, db_version) = database_cells(first);
            for run in &prot_ids {
                let mut common = common_protein_hit_keys(&run.hits);
                common.remove("Description");
                protein_keys.extend(common);
            }
        }
        let protein_keys: BTreeSet<String> =
            protein_keys.iter().map(|key| underscored(key)).collect();
        let mut prt_names: Vec<String> = protein_keys
            .iter()
            .map(|key| format!("opt_global_{key}"))
            .collect();
        for name in &mut prt_names {
            if name == TARGET_DECOY_COLUMN {
                *name = DECOY_HIT_COLUMN.into();
            }
        }
        prt_names.push(RESULT_TYPE_COLUMN.into());

        let design = AssayDesign::from_map(map)?;

        meta.title = string(&options.title);
        let experiment_type = map.experiment_type.as_str();
        meta.quantification_method = match experiment_type {
            "label-free" => parameter("[MS,MS:1001834,LC-MS label-free quantitation analysis,]")?,
            "labeled_MS1" => parameter("[PRIDE,PRIDE_0000316,MS1 based isotope labeling,]")?,
            "labeled_MS2" => parameter("[PRIDE,PRIDE_0000317,MS2 based isotope labeling,]")?,
            _ => MzTabParameter::default(),
        };
        meta.protein_quantification_unit = parameter("[,,Abundance,]")?;
        meta.peptide_quantification_unit = parameter("[,,Abundance,]")?;

        let mut ms_run_paths = map.primary_ms_run_path();
        ms_run_paths.dedup();
        let mut path_to_run: BTreeMap<String, usize> = BTreeMap::new();
        for (offset, path) in ms_run_paths.iter().enumerate() {
            let run_number = offset + 1;
            path_to_run.entry(path.clone()).or_insert(run_number);
            let run = MzTabMSRunMetaData {
                format: parameter("[MS,MS:1000584,mzML file,]")?,
                id_format: id_format.clone(),
                location: file_location(path),
                ..MzTabMSRunMetaData::default()
            };
            meta.ms_run.insert(run_number, run);
        }
        for header in map.column_headers.values() {
            let label = header.label_as_uint(experiment_type)?;
            let assay = design.assay_of(&header.filename, label)?;
            let reagent = match experiment_type {
                "label-free" => parameter("[MS,MS:1002038,unlabeled sample,]")?,
                "labeled_MS1" => parameter(&format!(
                    "[PRIDE,PRIDE:0000316,MS1 based isotope labeling,{}]",
                    header.label
                ))?,
                "labeled_MS2" => parameter(&format!(
                    "[PRIDE,PRIDE:0000317,MS2 based isotope labeling,{}]",
                    header.label
                ))?,
                _ => MzTabParameter::default(),
            };
            let run_number = *path_to_run.get(&header.filename).ok_or_else(|| {
                missing(format!(
                    "The consensus map column '{}' has no ms_run entry; the primary MS run paths and the column headers disagree.",
                    header.filename
                ))
            })?;
            let entry = meta.assay.entry(assay).or_default();
            entry.quantification_reagent = reagent;
            let reference = i32::try_from(run_number).map_err(|_| bad("ms_run index overflows"))?;
            if !entry.ms_run_ref.contains(&reference) {
                entry.ms_run_ref.push(reference);
            }
        }
        for study_variable in 1..=design.study_variables {
            let entry = meta.study_variable.entry(study_variable).or_default();
            entry.description.from_cell_string("no description given");
            for assay in &design.study_variable_to_assays[study_variable] {
                entry
                    .assay_refs
                    .push(i32::try_from(*assay).map_err(|_| bad("assay index overflows"))?);
            }
        }

        Ok(Self {
            map,
            prot_ids,
            peptide_ids,
            protein_hit_user_value_keys: protein_keys,
            feature_keys,
            identification_keys,
            hit_keys,
            run_index,
            run_to_engines: engines.run_to_engines,
            run_file_to_ms_file,
            design: AssayDesignDebug(design),
            ms_runs: ms_run_paths.len(),
            first_run_inference,
            fixed_mods,
            export_unidentified_features: options.export_unidentified_features,
            export_subfeatures: options.export_subfeatures,
            export_empty_pep_ids: options.export_empty_pep_ids,
            export_all_psms: options.export_all_psms,
            quant_study_variables,
            protein_cursor: ProteinCursorDebug::default(),
            pep_id: 0,
            pep_counter: 0,
            psm_id: 0,
            current_psm_index: 0,
            db,
            db_version,
            prt_optional_column_names: prt_names,
            pep_optional_column_names: pep_names,
            psm_optional_column_names: psm_names,
            meta_data: meta,
        })
    }

    /// The complete metadata section.
    pub fn meta_data(&self) -> &MzTabMetaData {
        &self.meta_data
    }

    /// The `PRT` optional column names.
    pub fn protein_optional_column_names(&self) -> &[String] {
        &self.prt_optional_column_names
    }

    /// The `PEP` optional column names. Subfeature entries are never named
    /// here; see [`ConsensusMapExportOptions::export_subfeatures`].
    pub fn peptide_optional_column_names(&self) -> &[String] {
        &self.pep_optional_column_names
    }

    /// The `PSM` optional column names: the hit keys, the target/decoy key
    /// renamed, then the peptidoform column.
    pub fn psm_optional_column_names(&self) -> &[String] {
        &self.psm_optional_column_names
    }

    /// The next `PRT` row, or `None` when the section is complete.
    ///
    /// As [`IdMzTabStream::next_prt_row`], except that an indistinguishable
    /// group's abundances are rebuilt from its parallel
    /// `fraction_group_level_*` arrays over the design's assays when it carries
    /// a consistent set of them.
    ///
    /// # Errors
    ///
    /// As [`IdMzTabStream::next_prt_row`].
    pub fn next_prt_row(&mut self) -> Result<Option<MzTabProteinSectionRow>> {
        let context = ProteinContext {
            prot_ids: &self.prot_ids,
            first_run_inference: self.first_run_inference,
            quant_study_variables: self.quant_study_variables,
            db: &self.db,
            db_version: &self.db_version,
            user_keys: &self.protein_hit_user_value_keys,
            assay_quantification: Some((
                &self.design.0.keys,
                &self.design.0.study_variable_to_assays,
            )),
        };
        self.protein_cursor.0.next(&context)
    }

    /// The next `PEP` row, one per consensus feature, or `None` at the end.
    ///
    /// Without `export_unidentified_features`, features none of whose
    /// identifications has a hit are skipped.
    ///
    /// # Errors
    ///
    /// As `peptide_row_from_consensus_feature`.
    pub fn next_pep_row(&mut self) -> Result<Option<MzTabPeptideSectionRow>> {
        let features = &self.map.features;
        let has_hits = |feature: &ConsensusFeature| {
            feature
                .peptide_identifications
                .iter()
                .any(|id| !id.hits.is_empty())
        };
        loop {
            let Some(feature) = features.get(self.pep_id) else {
                return Ok(None);
            };
            if self.export_unidentified_features || has_hits(feature) {
                break;
            }
            self.pep_id += 1;
        }
        let feature = &features[self.pep_id];
        let context = ConsensusRowContext {
            map: self.map,
            ms_runs: self.ms_runs,
            design: &self.design.0,
            feature_keys: &self.feature_keys,
            identification_keys: &self.identification_keys,
            hit_keys: &self.hit_keys,
            run_index: &self.run_index,
            run_file_to_ms_file: &self.run_file_to_ms_file,
            fixed_mods: &self.fixed_mods,
            export_subfeatures: self.export_subfeatures,
        };
        let row = peptide_row_from_consensus_feature(feature, &context)?;
        self.pep_id += 1;
        Ok(Some(row))
    }

    /// The next `PSM` row, or `None` when every identification was visited.
    ///
    /// As [`IdMzTabStream::next_psm_row`], including the placeholder rows.
    ///
    /// # Errors
    ///
    /// As [`IdMzTabStream::next_psm_row`].
    pub fn next_psm_row(&mut self) -> Result<Option<MzTabPSMSectionRow>> {
        let Some(identification) = self.peptide_ids.get(self.pep_counter).copied() else {
            return Ok(None);
        };
        let context = PsmContext {
            prot_ids: &self.prot_ids,
            run_index: &self.run_index,
            run_file_to_ms_file: &self.run_file_to_ms_file,
            run_to_engines: &self.run_to_engines,
            db: &self.db,
            db_version: &self.db_version,
            export_empty_pep_ids: self.export_empty_pep_ids,
            export_all_psms: self.export_all_psms,
        };
        let row = psm_row_from_peptide_id(
            identification,
            &context,
            self.current_psm_index,
            self.psm_id,
        )?;
        advance_psm(
            self.export_all_psms,
            identification.hits.len(),
            &mut self.pep_counter,
            &mut self.current_psm_index,
        );
        self.psm_id += 1;
        Ok(Some(row.unwrap_or_default()))
    }
}

impl MzTab {
    /// Export a consensus map as a `Quantification`/`Summary` document with
    /// `PRT`, `PEP` and `PSM` sections. Source `MzTab::exportConsensusMapToMzTab`,
    /// which drains a [`ConsensusMapMzTabStream`].
    ///
    /// PSM placeholder rows (null sequence) are dropped, as the source does
    /// here.
    ///
    /// # Errors
    ///
    /// As [`ConsensusMapMzTabStream::new`] and its three `next_*_row` methods,
    /// plus [`Error::InvalidValue`] when a section would exceed
    /// [`MzTab::MAX_ROWS`].
    pub fn export_consensus_map(
        consensus_map: &ConsensusMap,
        filename: &str,
        options: &ConsensusMapExportOptions,
    ) -> Result<MzTab> {
        let mut stream = ConsensusMapMzTabStream::new(consensus_map, filename, options)?;
        let mut document = MzTab {
            meta_data: stream.meta_data().clone(),
            ..MzTab::default()
        };
        while let Some(row) = stream.next_prt_row()? {
            push_bounded(&mut document.protein_data, row)?;
        }
        while let Some(row) = stream.next_pep_row()? {
            push_bounded(&mut document.peptide_data, row)?;
        }
        while let Some(row) = stream.next_psm_row()? {
            if !row.sequence.is_null() {
                push_bounded(&mut document.psm_data, row)?;
            }
        }
        Ok(document)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn remap_renames_and_maps_values() {
        let mut opt = vec![
            MzTabOptionalColumnEntry::new(TARGET_DECOY_COLUMN, string("target+decoy")),
            MzTabOptionalColumnEntry::new(DECOY_PEPTIDE_COLUMN, string("decoy")),
            MzTabOptionalColumnEntry::new(TARGET_DECOY_COLUMN, string("other")),
        ];
        remap_target_decoy(&mut opt, DECOY_PEPTIDE_COLUMN);
        assert!(opt.iter().all(|entry| entry.name == DECOY_PEPTIDE_COLUMN));
        assert_eq!(opt[0].value.get(), "0");
        assert_eq!(opt[1].value.get(), "1");
        assert_eq!(opt[2].value.get(), "other");
    }

    #[test]
    fn ms_file_names_register_each_path_once() {
        let mut a = ProteinIdentification::new();
        a.primary_ms_run_paths = vec!["x.mzML".into(), "y.mzML".into()];
        let mut b = ProteinIdentification::new();
        b.primary_ms_run_paths = vec!["y.mzML".into()];
        let c = ProteinIdentification::new();
        let runs = [&a, &b, &c];
        let (names, indices) = map_ms_file_names(&runs, false);
        assert_eq!(names["x.mzML"], 1);
        assert_eq!(names["y.mzML"], 2);
        assert_eq!(names["3"], 3);
        assert_eq!(indices.len(), 3);
        let files = map_run_file_to_ms_file(&runs, &names, false).unwrap();
        assert_eq!(files[&(1, 0)], 2);
        // A run without paths maps to its run index, not an ms_run index.
        assert_eq!(files[&(2, 0)], 2);
    }

    #[test]
    fn best_hit_keeps_the_first_of_equal_scores() {
        let mut id = PeptideIdentification::new();
        for (score, sequence) in [(1.0, "AAA"), (2.0, "CCC"), (2.0, "DDD")] {
            id.hits.push(
                PeptideHit::new(
                    score,
                    1,
                    2,
                    crate::chemistry::AASequence::parse(sequence).unwrap(),
                )
                .unwrap(),
            );
        }
        assert_eq!(best_hit(&id).unwrap().sequence.as_str(), "CCC");
        id.higher_score_better = false;
        assert_eq!(best_hit(&id).unwrap().sequence.as_str(), "AAA");
    }
}

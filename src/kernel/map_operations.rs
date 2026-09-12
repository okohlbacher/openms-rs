// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl, Timo Sachsenberg, OpenMS Rust contributors $

//! Container operations of `KERNEL/FeatureMap.h` and `KERNEL/ConsensusMap.h`.
//!
//! The containers themselves live in
//! [`kernel::features`](crate::kernel::features); this module adds the members
//! that the member-by-member review of the two headers found missing: the
//! annotation-state summary
//! ([`AnnotationStatistics`](crate::kernel::map_operations::AnnotationStatistics)),
//! map arithmetic, swapping, protein-identification lookup, the primary MS run
//! path, the container-level unique-ID walk, the identification-graph queries,
//! splitting a consensus map back into feature maps
//! ([`SplitMeta`](crate::kernel::map_operations::SplitMeta)) and the two stream
//! layouts. See `docs/MAP_OPERATIONS_SUPPORT.md`.
//!
//! Operations live here rather than in `kernel/features.rs` for the same reason
//! as in [`kernel::feature_identification`](crate::kernel::feature_identification)
//! and [`kernel::gap_closures`](crate::kernel::gap_closures): that file is
//! frozen for this work package, and Rust permits inherent and trait
//! implementations for a crate-local type from any module of the same crate, so
//! callers see one type with one API.
//!
//! # Identification records attach by reference
//!
//! The source maps own an `IdentificationData id_data_` member
//! (`FeatureMap.h:294`, `ConsensusMap.h:376`). This port deliberately does not
//! embed the graph: it is not `Clone` by design, so embedding it would strip
//! `Clone` and `PartialEq` from both maps. Operations that must resolve a
//! reference take the graph, or the
//! [`ReferenceTranslator`](crate::identification::graph::ReferenceTranslator),
//! as a parameter, exactly as
//! [`BaseFeature::update_id_references`](crate::kernel::features::BaseFeature::update_id_references)
//! already does. The source's `getIdentificationData()` accessors therefore
//! have no counterpart: the graph is the caller's.
//!
//! # Serial by design
//!
//! Neither source file uses OpenMP, so no parallel section is lost here. The
//! port is serial throughout, as `docs/REPOSITORY_ANALYSIS.md` requires.

use super::feature_identification::AnnotationState;
use super::features::{
    BaseFeature, ColumnHeader, ConsensusFeature, ConsensusMap, Feature, FeatureHandle, FeatureMap,
};
use crate::concept::{HasUniqueId, UniqueIdGenerator};
use crate::format::FileType;
use crate::identification::graph::{IdentificationData, ObservationMatchId, ReferenceTranslator};
use crate::identification::{PeptideIdentification, ProteinIdentification};
use crate::metadata::{DataProcessing, MetaValue};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

/// Source metadata key recording the map a peptide identification came from.
const MAP_INDEX: &str = "map_index";
/// Source metadata key holding the primary MS run paths of a feature map.
const SPECTRA_DATA: &str = "spectra_data";
/// Source placeholder written when no MS run is annotated.
const UNKNOWN_RUN: &str = "UNKNOWN";
/// Source filename written into every merged column header by `appendRows`.
const MERGED_FILENAME: &str = "mergedConsensusXMLFile";
/// Software name that marks a consensus map as isobarically quantified.
const ISOBARIC_ANALYZER: &str = "IsobaricAnalyzer";
/// Redraws a single element may consume while conflicting unique IDs are
/// resolved. The source `while` loop is unbounded.
const MAX_UNIQUE_ID_REDRAWS: usize = 64;

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn limit() -> Error {
    Error::InvalidValue("map operation work limit exceeded".into())
}

// ---------------------------------------------------------------------------
// AnnotationStatistics
// ---------------------------------------------------------------------------

/// Number of annotation states, the source `SIZE_OF_ANNOTATIONSTATE`.
const STATE_COUNT: usize = AnnotationState::NAMES.len();

/// Summary of the peptide identifications assigned to the features of a map.
///
/// Ports the `AnnotationStatistics` struct declared in `KERNEL/FeatureMap.h`.
/// Each feature contributes one vote, namely its
/// [`AnnotationState`];
/// subordinate features are not counted, matching
/// `FeatureMap::getAnnotationStatistics`, which iterates top-level features
/// only.
///
/// The source member is `std::vector<Size> states`, sized at construction by
/// the `SIZE_OF_ANNOTATIONSTATE` sentinel and indexed by the state's
/// discriminant. Here it is a fixed-length array, so the vector can never have
/// the wrong length and the sentinel has no counterpart. The source's copy
/// constructor, copy assignment and `operator==` are `Clone` and `PartialEq`;
/// `operator+=` is [`AnnotationStatistics::add`] and the
/// [`std::ops::AddAssign`] implementation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct AnnotationStatistics {
    states: [usize; STATE_COUNT],
}

impl AnnotationStatistics {
    /// All counts zero, as the source default constructor.
    pub const fn new() -> Self {
        Self {
            states: [0; STATE_COUNT],
        }
    }

    /// The counts, in [`AnnotationState`] discriminant order.
    ///
    /// The source exposes the `states` vector as a public member; this returns
    /// the array by reference, so a caller cannot resize it.
    pub const fn states(&self) -> &[usize; STATE_COUNT] {
        &self.states
    }

    /// How many features reported `state`.
    pub const fn count(&self, state: AnnotationState) -> usize {
        self.states[state as usize]
    }

    /// Record one more feature in `state`.
    ///
    /// Ports `AnnotationStatistics::operator+=(BaseFeature::AnnotationState)`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the count for `state` would
    /// overflow `usize`; the statistics are unchanged. The source increments an
    /// unchecked `Size`.
    pub fn add(&mut self, state: AnnotationState) -> Result<()> {
        let slot = &mut self.states[state as usize];
        *slot = slot.checked_add(1).ok_or_else(limit)?;
        Ok(())
    }
}

/// Saturating form of [`AnnotationStatistics::add`] for the `+=` spelling the
/// source uses. A count at `usize::MAX` stays there rather than wrapping.
impl std::ops::AddAssign<AnnotationState> for AnnotationStatistics {
    fn add_assign(&mut self, state: AnnotationState) {
        let slot = &mut self.states[state as usize];
        *slot = slot.saturating_add(1);
    }
}

/// Source `operator<<(std::ostream&, const AnnotationStatistics&)`: a heading,
/// one indented `name: count` line per state in discriminant order, and a
/// trailing blank line from the closing `std::endl`.
impl fmt::Display for AnnotationStatistics {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Feature annotation with identifications:")?;
        for (index, name) in AnnotationState::NAMES.iter().enumerate() {
            writeln!(f, "    {}: {}", name, self.states[index])?;
        }
        writeln!(f)
    }
}

// ---------------------------------------------------------------------------
// SplitMeta
// ---------------------------------------------------------------------------

/// What [`ConsensusMap::split`] does with the meta values of a consensus feature.
///
/// Ports `ConsensusMap::SplitMeta`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SplitMeta {
    /// Do not copy any meta values. The source default.
    #[default]
    Discard,
    /// Copy all meta values to every feature map the consensus feature reaches.
    CopyAll,
    /// Copy all meta values to the first feature map only.
    CopyFirst,
}

/// The column headers of a consensus map, keyed by map index.
///
/// Ports `ConsensusMap::ColumnHeaders`, the source `std::map<UInt64,
/// ColumnHeader>` type alias, as a deterministic `BTreeMap`. It is the type of
/// the public
/// [`ConsensusMap::column_headers`](crate::kernel::features::ConsensusMap)
/// field.
pub type ColumnHeaders = BTreeMap<u64, ColumnHeader>;

// ---------------------------------------------------------------------------
// FeatureMap
// ---------------------------------------------------------------------------

impl FeatureMap {
    /// Features, identifications and processing records one checked operation
    /// will combine or traverse.
    pub const MAX_ITEMS: usize = 10_000_000;

    /// Summarise the annotation state of every top-level feature.
    ///
    /// Ports `FeatureMap::getAnnotationStatistics`. Subordinate features are
    /// not visited, as in the source.
    ///
    /// # Arguments
    ///
    /// * `graph` — the `IdentificationData` in which any observation matches
    ///   referenced by the features were registered. Pass `None` when no
    ///   feature references a match; see
    ///   [`BaseFeature::annotation_state`](crate::kernel::features::BaseFeature::annotation_state)
    ///   for why the graph is a parameter here.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the map holds more than
    /// [`FeatureMap::MAX_ITEMS`] features, and any error of
    /// `BaseFeature::annotation_state` — in particular
    /// [`Error::MissingInformation`] when a feature references matches and
    /// `graph` is `None`. The map is only read.
    pub fn annotation_statistics(
        &self,
        graph: Option<&IdentificationData>,
    ) -> Result<AnnotationStatistics> {
        if self.features.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        let mut result = AnnotationStatistics::new();
        for feature in &self.features {
            result.add(feature.annotation_state(graph)?)?;
        }
        Ok(result)
    }

    /// Append every feature and record of `rhs` to this map.
    ///
    /// Ports `FeatureMap::operator+=`. Features, protein identifications,
    /// unassigned peptide identifications and data-processing records are
    /// concatenated. The document identifier and the container's own unique ID
    /// are reset to their defaults, and the source's note that conflicting
    /// unique IDs receive new ones is preserved: after appending, any feature
    /// whose non-zero unique ID now occurs twice is redrawn from `generator`,
    /// and the number of redraws is returned.
    ///
    /// The source additionally resets the cached range information and merges
    /// the two `IdentificationData` graphs. Ranges are computed on demand by
    /// [`FeatureMap::ranges`](crate::kernel::features::FeatureMap::ranges), so
    /// there is no cache to reset; the graphs are the caller's, so merge them
    /// with `IdentificationData::merge_from` and pass the resulting translator
    /// to [`FeatureMap::update_id_references`] for the appended features. The
    /// source logs that document identifiers are lost; this port does not log,
    /// because no module of this crate does.
    ///
    /// # Arguments
    ///
    /// * `rhs` — the map to append. It is only read.
    /// * `generator` — caller-owned replacement for the source's process-wide
    ///   `UniqueIdGenerator` singleton, used only when unique IDs collide.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the combined number of features,
    /// identifications or processing records would exceed
    /// [`FeatureMap::MAX_ITEMS`], or when one element needs more than 64
    /// redraws to obtain a free unique ID — the source loops until the
    /// generator happens to produce one. The result is built into a temporary
    /// and committed only on success, so a failure leaves this map exactly as
    /// it was.
    pub fn append(&mut self, rhs: &Self, generator: &mut UniqueIdGenerator) -> Result<usize> {
        combine_limit(self.features.len(), rhs.features.len())?;
        combine_limit(
            self.protein_identifications.len(),
            rhs.protein_identifications.len(),
        )?;
        combine_limit(
            self.unassigned_peptide_identifications.len(),
            rhs.unassigned_peptide_identifications.len(),
        )?;
        combine_limit(self.data_processing.len(), rhs.data_processing.len())?;
        let mut merged = self.clone();
        merged.identifier = String::new();
        merged.loaded_file_path = String::new();
        merged.loaded_file_type = FileType::Unknown;
        merged.unique_id = 0;
        merged
            .protein_identifications
            .extend(rhs.protein_identifications.iter().cloned());
        merged
            .unassigned_peptide_identifications
            .extend(rhs.unassigned_peptide_identifications.iter().cloned());
        merged
            .data_processing
            .extend(rhs.data_processing.iter().cloned());
        merged.features.extend(rhs.features.iter().cloned());
        let replaced =
            if has_duplicate_unique_ids(merged.features.iter().map(|feature| feature.unique_id)) {
                resolve_unique_id_conflicts(
                    &mut merged.features,
                    |feature| &mut feature.base.unique_id,
                    generator,
                )?
            } else {
                0
            };
        *self = merged;
        Ok(replaced)
    }

    /// This map with `rhs` appended, leaving both operands unchanged.
    ///
    /// Ports `FeatureMap::operator+`, which the source implements as a copy
    /// followed by `operator+=`. The redraw count [`FeatureMap::append`]
    /// returns is not reported here, matching the source's `FeatureMap` return
    /// type; call `append` on a clone when it is wanted.
    ///
    /// # Errors
    ///
    /// As [`FeatureMap::append`].
    pub fn merged(&self, rhs: &Self, generator: &mut UniqueIdGenerator) -> Result<Self> {
        let mut result = self.clone();
        result.append(rhs, generator)?;
        Ok(result)
    }

    /// Exchange the features of this map with those of `from`.
    ///
    /// Ports `FeatureMap::swapFeaturesOnly`, which swaps the feature vector and
    /// the cached range information so that neither map is left with a stale
    /// cache. Ranges are computed on demand here, so only the features move;
    /// every other member, the document identifier and unique ID included,
    /// stays with its map.
    pub fn swap_features_only(&mut self, from: &mut Self) {
        std::mem::swap(&mut self.features, &mut from.features);
    }

    /// Exchange the contents of this map with those of `from`.
    ///
    /// Ports `FeatureMap::swap`, which swaps the features and ranges, the
    /// document identifier, the unique ID and its index, the protein and
    /// unassigned peptide identifications, the data-processing records and the
    /// identification data.
    ///
    /// The source does **not** swap the inherited `MetaInfoInterface`, so meta
    /// values stay with their map while everything else moves. That asymmetry
    /// is preserved here: `metadata` is the one field this method leaves alone.
    /// It is almost certainly a source oversight — `clear(true)` and
    /// `operator==` both treat the meta values as part of the map — so a caller
    /// who wants a complete exchange should use [`std::mem::swap`] instead.
    ///
    /// Identification references are unaffected: they are owner-tagged IDs into
    /// a caller-held graph, not iterators into a member, so nothing needs
    /// translating after a swap.
    pub fn swap(&mut self, from: &mut Self) {
        self.swap_features_only(from);
        std::mem::swap(&mut self.identifier, &mut from.identifier);
        std::mem::swap(&mut self.loaded_file_path, &mut from.loaded_file_path);
        std::mem::swap(&mut self.loaded_file_type, &mut from.loaded_file_type);
        std::mem::swap(&mut self.unique_id, &mut from.unique_id);
        std::mem::swap(
            &mut self.protein_identifications,
            &mut from.protein_identifications,
        );
        std::mem::swap(
            &mut self.unassigned_peptide_identifications,
            &mut from.unassigned_peptide_identifications,
        );
        std::mem::swap(&mut self.data_processing, &mut from.data_processing);
    }

    /// The first protein identification whose identifier is `identifier`.
    ///
    /// Ports the `const` `FeatureMap::findProteinIdentification`, which returns
    /// `nullptr` when there is no match; this returns `None`. The search is
    /// linear and stops at the first match, as in the source, so a repeated
    /// identifier hides the later records.
    pub fn find_protein_identification(&self, identifier: &str) -> Option<&ProteinIdentification> {
        self.protein_identifications
            .iter()
            .find(|record| record.identifier == identifier)
    }

    /// Mutable form of [`FeatureMap::find_protein_identification`].
    ///
    /// Ports the non-`const` `FeatureMap::findProteinIdentification`.
    pub fn find_protein_identification_mut(
        &mut self,
        identifier: &str,
    ) -> Option<&mut ProteinIdentification> {
        self.protein_identifications
            .iter_mut()
            .find(|record| record.identifier == identifier)
    }

    /// Annotate the file paths of the primary MS runs behind this map.
    ///
    /// Ports `FeatureMap::setPrimaryMSRunPath(const StringList&)`, which stores
    /// the list in the `spectra_data` meta value. An empty list is stored as an
    /// empty list, as in the source.
    ///
    /// The source logs a warning for an empty list and for every path that does
    /// not end in `mzML` or `mzml`, because only mzML keeps results traceable.
    /// No module of this crate logs, so the advisory is documented rather than
    /// emitted; nothing about the stored value depends on it.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the list holds more than
    /// [`FeatureMap::MAX_ITEMS`] paths. The meta value is replaced only on
    /// success.
    pub fn set_primary_ms_run_path(&mut self, paths: &[String]) -> Result<()> {
        if paths.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        self.metadata
            .insert(SPECTRA_DATA.into(), MetaValue::from(paths.to_vec()));
        Ok(())
    }

    /// Annotate the primary MS run paths from `experiment`, falling back to
    /// `paths`.
    ///
    /// Ports `FeatureMap::setPrimaryMSRunPath(const StringList&, MSExperiment&)`.
    /// The experiment's own primary run path is used when it names exactly one
    /// file, that file ends in `mzML`, and it exists on disk; otherwise `paths`
    /// is stored. The existence test touches the filesystem, as the source
    /// `File::exists` does, and a path that cannot be examined counts as
    /// absent.
    ///
    /// The source takes the experiment by non-`const` reference although it
    /// only reads it; this takes it by shared reference.
    ///
    /// # Errors
    ///
    /// As [`FeatureMap::set_primary_ms_run_path`].
    pub fn set_primary_ms_run_path_from_experiment(
        &mut self,
        paths: &[String],
        experiment: &super::MSExperiment,
    ) -> Result<()> {
        match usable_experiment_run_path(experiment) {
            Some(path) => self.set_primary_ms_run_path(&[path]),
            None => self.set_primary_ms_run_path(paths),
        }
    }

    /// The annotated primary MS run paths.
    ///
    /// Ports `FeatureMap::getPrimaryMSRunPath(StringList&)`. When nothing is
    /// annotated, or the annotation is an empty list, the source pushes the
    /// literal `UNKNOWN` and logs a warning; the placeholder is preserved here
    /// and the warning is not logged.
    ///
    /// The source appends into the caller's list without clearing it, so a
    /// non-empty output argument suppresses the placeholder and keeps stale
    /// entries. This returns a fresh vector, which is the evident intent.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `spectra_data` is present but is
    /// not a string list.
    pub fn primary_ms_run_path(&self) -> Result<Vec<String>> {
        let mut result = match self.metadata.get(SPECTRA_DATA) {
            Some(value) => value.as_string_list()?.to_vec(),
            None => Vec::new(),
        };
        if result.is_empty() {
            result.push(UNKNOWN_RUN.into());
        }
        Ok(result)
    }

    /// Apply `visit` to the unique ID of the map itself, of every feature and
    /// of every subordinate feature, returning the sum of its results.
    ///
    /// Ports the mutating `FeatureMap::applyMemberFunction`, whose C++
    /// member-function pointer is a closure here. The container itself is
    /// visited first, then each feature's own subtree in order, depth first,
    /// exactly as the source accumulates.
    ///
    /// ```
    /// use openms::concept::HasUniqueId;
    /// use openms::kernel::features::{Feature, FeatureMap};
    ///
    /// let mut map = FeatureMap::from_features(vec![Feature::default(), Feature::default()]);
    /// map.features[1].subordinates.push(Feature::default());
    /// // Source: `fm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId)`.
    /// let invalid = map.for_each_unique_id(|id| usize::from(id.has_invalid_unique_id()))?;
    /// assert_eq!(invalid, 4); // the container, two features and one subordinate
    /// # Ok::<(), openms::Error>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the map holds more than
    /// [`FeatureMap::MAX_ITEMS`] features, when a feature tree exceeds the
    /// limits of
    /// [`Feature::for_each_unique_id`](crate::kernel::features::Feature::for_each_unique_id),
    /// or when the accumulated count overflows. The limits are checked while
    /// walking, so `visit` may already have run on earlier features; the source
    /// is unbounded.
    pub fn for_each_unique_id(
        &mut self,
        mut visit: impl FnMut(&mut u64) -> usize,
    ) -> Result<usize> {
        if self.features.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        let mut total = visit(&mut self.unique_id);
        for feature in &mut self.features {
            total = total
                .checked_add(feature.for_each_unique_id(&mut visit)?)
                .ok_or_else(limit)?;
        }
        Ok(total)
    }

    /// Apply `visit` to the unique ID of the map, of every feature and of every
    /// subordinate feature without mutating them.
    ///
    /// Ports the `const` `FeatureMap::applyMemberFunction`, which exists in the
    /// source only because a member-function pointer is either `const` or not.
    /// Order, accumulation and limits match [`FeatureMap::for_each_unique_id`].
    ///
    /// # Errors
    ///
    /// As [`FeatureMap::for_each_unique_id`].
    pub fn count_unique_ids(&self, mut visit: impl FnMut(u64) -> usize) -> Result<usize> {
        if self.features.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        let mut total = visit(self.unique_id);
        for feature in &self.features {
            total = total
                .checked_add(feature.count_unique_ids(&mut visit)?)
                .ok_or_else(limit)?;
        }
        Ok(total)
    }

    /// Observation matches in `graph` that no top-level feature references.
    ///
    /// Ports `FeatureMap::getUnassignedIDMatches`, the set difference between
    /// every match registered in the graph and the matches the features claim.
    /// Only top-level features are considered; the source's `@TODO` about
    /// subordinates is answered the same way, by not descending.
    ///
    /// The source reaches the graph through its own `id_data_` member; here it
    /// is a parameter, so a caller can ask the question of any graph the
    /// references belong to.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the map holds more than
    /// [`FeatureMap::MAX_ITEMS`] features or a feature references more than
    /// [`BaseFeature::MAX_ID_MATCHES`](crate::kernel::features::BaseFeature::MAX_ID_MATCHES)
    /// matches.
    pub fn unassigned_id_matches(
        &self,
        graph: &IdentificationData,
    ) -> Result<BTreeSet<ObservationMatchId>> {
        if self.features.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        unassigned_matches(
            graph,
            self.features.iter().map(|feature| &feature.base.id_matches),
        )
    }

    /// Translate the identification references of every feature, subordinates
    /// included, into the graph the translator describes.
    ///
    /// The source performs this inside its copy constructor, copy assignment
    /// and `operator+=`, because those merge the embedded `IdentificationData`
    /// and must repoint the features afterwards. Here the graph is the
    /// caller's, so the step is an explicit call made after
    /// `IdentificationData::merge_from` or
    /// `IdentificationData::try_clone_with_translation`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the translator has no entry for one
    /// of the referenced records, and the traversal limits of
    /// [`Feature::update_all_id_references`](crate::kernel::features::Feature::update_all_id_references).
    /// Every reference in the whole map is verified in a first, read-only pass,
    /// so a missing translation anywhere leaves the entire map unchanged; the
    /// source updates feature by feature.
    pub fn update_id_references(&mut self, translator: &ReferenceTranslator) -> Result<()> {
        if self.features.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        for feature in &self.features {
            verify_feature_tree(feature, translator)?;
        }
        for feature in &mut self.features {
            feature.update_all_id_references(translator)?;
        }
        Ok(())
    }
}

/// Source `operator<<(std::ostream&, const FeatureMap&)`: a banner, a column
/// header line, then one tab-separated line per feature holding the position
/// (retention time and m/z separated by a space, as `DPosition`'s own stream
/// operator writes it), intensity, overall quality, charge and unique ID, and
/// finally a closing banner. Subordinate features are not printed.
///
/// Numbers use Rust formatting rather than the C++ stream's
/// six-significant-digit default, matching every other kernel `Display`.
impl fmt::Display for FeatureMap {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "# -- DFEATUREMAP BEGIN --")?;
        writeln!(f, "# POS \tINTENS\tOVALLQ\tCHARGE\tUniqueID")?;
        for feature in &self.features {
            writeln!(
                f,
                "{} {}\t{}\t{}\t{}\t{}",
                feature.rt,
                feature.mz,
                feature.intensity,
                feature.quality,
                feature.charge,
                feature.unique_id
            )?;
        }
        writeln!(f, "# -- DFEATUREMAP END --")
    }
}

// ---------------------------------------------------------------------------
// ConsensusMap
// ---------------------------------------------------------------------------

impl ConsensusMap {
    /// Consensus features, identifications and processing records one checked
    /// operation will combine or traverse.
    pub const MAX_ITEMS: usize = 10_000_000;
    /// Column headers one checked operation will combine.
    pub const MAX_COLUMNS: usize = 1_000_000;

    /// A map of `n` default consensus features.
    ///
    /// Ports `explicit ConsensusMap(size_type n)`. Every other member takes its
    /// default, so the experiment type is `label-free`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `n` exceeds
    /// [`ConsensusMap::MAX_ITEMS`]; the source allocates unchecked.
    pub fn with_size(n: usize) -> Result<Self> {
        if n > Self::MAX_ITEMS {
            return Err(limit());
        }
        Ok(Self::from_features(vec![ConsensusFeature::new(); n]))
    }

    /// Set the type of experiment this map quantifies.
    ///
    /// Ports `ConsensusMap::setExperimentType`, which accepts only
    /// `label-free`, `labeled_MS1` and `labeled_MS2`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for any other string, matching the
    /// source `Exception::IllegalArgument`. The stored type is unchanged. The
    /// public `experiment_type` field bypasses this check, as the source's
    /// private member cannot;
    /// [`ConsensusMap::validate`](crate::kernel::features::ConsensusMap::validate)
    /// rejects an invalid value written that way.
    pub fn set_experiment_type(&mut self, experiment_type: &str) -> Result<()> {
        if !matches!(
            experiment_type,
            "label-free" | "labeled_MS1" | "labeled_MS2"
        ) {
            return Err(invalid(
                "unknown experiment type; must be one of label-free, labeled_MS1, labeled_MS2",
            ));
        }
        self.experiment_type = experiment_type.to_owned();
        Ok(())
    }

    /// Append the consensus features of `rhs` as new rows.
    ///
    /// Ports `ConsensusMap::appendRows`. The number of columns does not grow:
    /// `rhs`'s column headers are inserted only for indices this map does not
    /// already use, matching `std::map::insert`, which never overwrites.
    ///
    /// Three source behaviours are preserved exactly, because callers can
    /// observe them:
    ///
    /// * The document identifier and the container's unique ID are reset.
    /// * After the headers are merged, the source walks the **merged** header
    ///   map and `rhs`'s header map in parallel and, for as many entries as the
    ///   shorter of the two has, rewrites the merged entry's filename to
    ///   `mergedConsensusXMLFile` and sets its size to the sum of the two
    ///   entries at that *position*. The pairing is positional, not by column
    ///   index, so which sizes are added depends on iteration order rather than
    ///   on which columns describe the same run. Headers beyond that prefix
    ///   keep their filename and size.
    /// * Every protein identification in the result — including the ones this
    ///   map already held — has its fixed and variable modification lists
    ///   sorted and deduplicated.
    ///
    /// Conflicting unique IDs are redrawn from `generator` and the number of
    /// redraws is returned, as for [`FeatureMap::append`]. The source also
    /// merges the embedded identification graphs; pass the translator from
    /// `IdentificationData::merge_from` to
    /// [`ConsensusMap::update_id_references`] instead.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the combined number of consensus
    /// features, identifications or processing records exceeds
    /// [`ConsensusMap::MAX_ITEMS`], when the combined header count exceeds
    /// [`ConsensusMap::MAX_COLUMNS`], when a header size sum overflows, or when
    /// an element needs more than 64 redraws for a free unique ID. The result
    /// is built into a temporary and committed only on success.
    pub fn append_rows(&mut self, rhs: &Self, generator: &mut UniqueIdGenerator) -> Result<usize> {
        let mut merged = self.preflight_append(rhs)?;
        let mut headers = merged.column_headers.clone();
        for (index, header) in &rhs.column_headers {
            headers.entry(*index).or_insert_with(|| header.clone());
        }
        let keys: Vec<u64> = headers.keys().copied().collect();
        for (key, other) in keys.iter().zip(rhs.column_headers.values()) {
            let header = headers
                .get_mut(key)
                .ok_or_else(|| invalid("merged consensus column disappeared"))?;
            header.filename = MERGED_FILENAME.to_owned();
            header.size = header.size.checked_add(other.size).ok_or_else(limit)?;
        }
        merged.column_headers = headers;
        merged
            .protein_identifications
            .extend(rhs.protein_identifications.iter().cloned());
        for record in &mut merged.protein_identifications {
            dedup_modifications(record);
        }
        merged
            .unassigned_peptide_identifications
            .extend(rhs.unassigned_peptide_identifications.iter().cloned());
        merged.features.extend(rhs.features.iter().cloned());
        let replaced = resolve_map_unique_ids(&mut merged.features, generator)?;
        *self = merged;
        Ok(replaced)
    }

    /// Append the columns of `rhs`, shifting its map indices past this map's.
    ///
    /// Ports `ConsensusMap::appendColumns`. The number of columns becomes the
    /// sum of both maps': every column header of `rhs` keyed `k` is inserted at
    /// `k + n`, where `n` is this map's header *count* before the merge — not
    /// its largest index, so a map whose headers are not keyed `0..n-1` can
    /// still collide, and a collision keeps the existing header because
    /// `std::map::insert` does not overwrite. Each appended consensus feature's
    /// handles and each `map_index` meta value, on assigned and unassigned
    /// peptide identifications alike, are shifted by the same `n`.
    ///
    /// Conflicting unique IDs are redrawn from `generator` and the number of
    /// redraws is returned. The source also merges the embedded identification
    /// graphs; use [`ConsensusMap::update_id_references`] instead.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a combined count exceeds
    /// [`ConsensusMap::MAX_ITEMS`] or [`ConsensusMap::MAX_COLUMNS`], when a
    /// shifted column index or map index overflows `u64`, when a `map_index`
    /// meta value is not a nonnegative integer, or when an element needs more
    /// than 64 redraws for a free unique ID. The result is built into a
    /// temporary and committed only on success.
    pub fn append_columns(
        &mut self,
        rhs: &Self,
        generator: &mut UniqueIdGenerator,
    ) -> Result<usize> {
        let mut merged = self.preflight_append(rhs)?;
        let shift = u64::try_from(merged.column_headers.len())
            .map_err(|_| invalid("consensus column count exceeds the map index range"))?;
        for (index, header) in &rhs.column_headers {
            let shifted = index.checked_add(shift).ok_or_else(|| {
                invalid("shifted consensus column index exceeds the map index range")
            })?;
            merged
                .column_headers
                .entry(shifted)
                .or_insert_with(|| header.clone());
        }
        merged
            .protein_identifications
            .extend(rhs.protein_identifications.iter().cloned());
        for record in &mut merged.protein_identifications {
            dedup_modifications(record);
        }
        for record in &rhs.unassigned_peptide_identifications {
            let mut copy = record.clone();
            shift_map_index(&mut copy, shift)?;
            merged.unassigned_peptide_identifications.push(copy);
        }
        for feature in &rhs.features {
            let mut copy = feature.clone();
            for record in &mut copy.base.peptide_identifications {
                shift_map_index(record, shift)?;
            }
            let mut handles = copy.handles().to_vec();
            for handle in &mut handles {
                handle.map_index = handle.map_index.checked_add(shift).ok_or_else(|| {
                    invalid("shifted consensus map index exceeds the map index range")
                })?;
            }
            copy.set_handles(handles)?;
            merged.features.push(copy);
        }
        let replaced = resolve_map_unique_ids(&mut merged.features, generator)?;
        *self = merged;
        Ok(replaced)
    }

    /// Sort each consensus feature's peptide identifications by map index.
    ///
    /// Ports `ConsensusMap::sortPeptideIdentificationsByMapIndex`. The sort is
    /// stable and ascending in the `map_index` meta value; identifications that
    /// carry no `map_index` move to the end and keep their relative order,
    /// exactly as the source comparator places them. Unassigned peptide
    /// identifications are not touched, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the map holds more than
    /// [`ConsensusMap::MAX_ITEMS`] consensus features, or when a `map_index`
    /// meta value is not an integer — the source compares the meta values
    /// themselves, so a non-integer value there yields an order that depends on
    /// the value's type. Every key is read before anything is reordered, so a
    /// failure leaves the map unchanged.
    pub fn sort_peptide_identifications_by_map_index(&mut self) -> Result<()> {
        if self.features.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        let mut sorted: Vec<Vec<PeptideIdentification>> = Vec::with_capacity(self.features.len());
        for feature in &self.features {
            let mut keyed: Vec<(u8, i64, PeptideIdentification)> =
                Vec::with_capacity(feature.base.peptide_identifications.len());
            for record in &feature.base.peptide_identifications {
                match record.metadata.get(MAP_INDEX) {
                    Some(value) => keyed.push((0, value.as_i64()?, record.clone())),
                    None => keyed.push((1, 0, record.clone())),
                }
            }
            keyed.sort_by_key(|entry| (entry.0, entry.1));
            sorted.push(keyed.into_iter().map(|entry| entry.2).collect());
        }
        for (feature, records) in self.features.iter_mut().zip(sorted) {
            feature.base.peptide_identifications = records;
        }
        Ok(())
    }

    /// Exchange the contents of this map with those of `from`.
    ///
    /// Ports `ConsensusMap::swap`, which swaps the consensus features and
    /// ranges, the document identifier, the unique ID and its index, the column
    /// headers, the experiment type, the protein and unassigned peptide
    /// identifications, the data-processing records and the identification
    /// data.
    ///
    /// As for [`FeatureMap::swap`], the source does **not** swap the inherited
    /// `MetaInfoInterface`; `metadata` therefore stays with its map here too.
    /// Use [`std::mem::swap`] for a complete exchange.
    pub fn swap(&mut self, from: &mut Self) {
        std::mem::swap(&mut self.features, &mut from.features);
        std::mem::swap(&mut self.identifier, &mut from.identifier);
        std::mem::swap(&mut self.loaded_file_path, &mut from.loaded_file_path);
        std::mem::swap(&mut self.loaded_file_type, &mut from.loaded_file_type);
        std::mem::swap(&mut self.unique_id, &mut from.unique_id);
        std::mem::swap(&mut self.column_headers, &mut from.column_headers);
        std::mem::swap(&mut self.experiment_type, &mut from.experiment_type);
        std::mem::swap(
            &mut self.protein_identifications,
            &mut from.protein_identifications,
        );
        std::mem::swap(
            &mut self.unassigned_peptide_identifications,
            &mut from.unassigned_peptide_identifications,
        );
        std::mem::swap(&mut self.data_processing, &mut from.data_processing);
    }

    /// The first protein identification whose identifier is `identifier`.
    ///
    /// Ports the `const` `ConsensusMap::findProteinIdentification`; `None`
    /// replaces the source `nullptr`.
    pub fn find_protein_identification(&self, identifier: &str) -> Option<&ProteinIdentification> {
        self.protein_identifications
            .iter()
            .find(|record| record.identifier == identifier)
    }

    /// Mutable form of [`ConsensusMap::find_protein_identification`].
    ///
    /// Ports the non-`const` `ConsensusMap::findProteinIdentification`.
    pub fn find_protein_identification_mut(
        &mut self,
        identifier: &str,
    ) -> Option<&mut ProteinIdentification> {
        self.protein_identifications
            .iter_mut()
            .find(|record| record.identifier == identifier)
    }

    /// Annotate the file path of each column's primary MS run.
    ///
    /// Ports `ConsensusMap::setPrimaryMSRunPath(const StringList&)`, which
    /// stores the paths in the column headers rather than in a meta value — the
    /// one place where the consensus and feature containers disagree about
    /// where this belongs.
    ///
    /// The source behaviour is preserved in full:
    ///
    /// * An empty list sets *every* existing header's filename to `UNKNOWN`.
    /// * Otherwise, when headers already exist, the number of paths must equal
    ///   the number of headers.
    /// * The paths are then written to the headers keyed `0..paths.len()-1`,
    ///   **by position, not by existing key**. `std::map::operator[]` inserts a
    ///   default header for a key that does not exist, so a map whose headers
    ///   are keyed otherwise gains new columns instead of being renamed. That
    ///   also means the count check above can pass while every write lands on a
    ///   fresh header.
    ///
    /// The source's advisory warnings — for an empty list, and for a path that
    /// does not end in `mzML` or `mzml` — are not logged, as no module of this
    /// crate logs.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the list is non-empty, headers
    /// exist and the two counts differ (the source
    /// `Exception::InvalidParameter`), or when more than
    /// [`ConsensusMap::MAX_COLUMNS`] paths are given. The headers are rebuilt
    /// into a temporary, so a failure leaves them untouched.
    pub fn set_primary_ms_run_path(&mut self, paths: &[String]) -> Result<()> {
        if paths.len() > Self::MAX_COLUMNS {
            return Err(limit());
        }
        let mut headers = self.column_headers.clone();
        if paths.is_empty() {
            for header in headers.values_mut() {
                header.filename = UNKNOWN_RUN.to_owned();
            }
            self.column_headers = headers;
            return Ok(());
        }
        if !headers.is_empty() && paths.len() != headers.len() {
            return Err(invalid(
                "number of MS run paths must match the number of consensus columns",
            ));
        }
        for (index, path) in paths.iter().enumerate() {
            let key = u64::try_from(index)
                .map_err(|_| invalid("consensus column index exceeds the map index range"))?;
            headers.entry(key).or_default().filename = path.clone();
        }
        self.column_headers = headers;
        Ok(())
    }

    /// Annotate the primary MS run paths from `experiment`, falling back to
    /// `paths`.
    ///
    /// Ports `ConsensusMap::setPrimaryMSRunPath(const StringList&, MSExperiment&)`,
    /// which selects between the two lists exactly as the `FeatureMap` overload
    /// does; see [`FeatureMap::set_primary_ms_run_path_from_experiment`] for
    /// the rule and for the filesystem access it implies.
    ///
    /// # Errors
    ///
    /// As [`ConsensusMap::set_primary_ms_run_path`]. Note that a single
    /// experiment path is stored into a map with more than one column only when
    /// that map has no headers yet, because the counts must otherwise match.
    pub fn set_primary_ms_run_path_from_experiment(
        &mut self,
        paths: &[String],
        experiment: &super::MSExperiment,
    ) -> Result<()> {
        match usable_experiment_run_path(experiment) {
            Some(path) => self.set_primary_ms_run_path(&[path]),
            None => self.set_primary_ms_run_path(paths),
        }
    }

    /// The filename of every column header, in column-index order.
    ///
    /// Ports `ConsensusMap::getPrimaryMSRunPath(StringList&)`. Unlike the
    /// `FeatureMap` counterpart there is no `UNKNOWN` placeholder: a map with
    /// no columns yields an empty list. The source appends into the caller's
    /// list without clearing it; this returns a fresh vector.
    pub fn primary_ms_run_path(&self) -> Vec<String> {
        self.column_headers
            .values()
            .map(|header| header.filename.clone())
            .collect()
    }

    /// Apply `visit` to the unique ID of the map and of every consensus
    /// feature, returning the sum of its results.
    ///
    /// Ports the mutating `ConsensusMap::applyMemberFunction`. The container is
    /// visited first, then each consensus feature. Unlike the `FeatureMap`
    /// overload there is no recursion, because a consensus feature has no
    /// subordinates; it holds handles, which are not visited.
    ///
    /// ```
    /// use openms::concept::HasUniqueId;
    /// use openms::kernel::features::{ConsensusFeature, ConsensusMap};
    ///
    /// let mut map = ConsensusMap::from_features(vec![ConsensusFeature::new(); 3]);
    /// // Source: `cm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId)`.
    /// let invalid = map.for_each_unique_id(|id| usize::from(id.has_invalid_unique_id()))?;
    /// assert_eq!(invalid, 4); // the container and three consensus features
    /// # Ok::<(), openms::Error>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the map holds more than
    /// [`ConsensusMap::MAX_ITEMS`] consensus features or the accumulated count
    /// overflows; the source is unbounded.
    pub fn for_each_unique_id(
        &mut self,
        mut visit: impl FnMut(&mut u64) -> usize,
    ) -> Result<usize> {
        if self.features.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        let mut total = visit(&mut self.unique_id);
        for feature in &mut self.features {
            total = total
                .checked_add(visit(&mut feature.base.unique_id))
                .ok_or_else(limit)?;
        }
        Ok(total)
    }

    /// Apply `visit` to the unique ID of the map and of every consensus feature
    /// without mutating them.
    ///
    /// Ports the `const` `ConsensusMap::applyMemberFunction`. Order,
    /// accumulation and limits match [`ConsensusMap::for_each_unique_id`].
    ///
    /// # Errors
    ///
    /// As [`ConsensusMap::for_each_unique_id`].
    pub fn count_unique_ids(&self, mut visit: impl FnMut(u64) -> usize) -> Result<usize> {
        if self.features.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        let mut total = visit(self.unique_id);
        for feature in &self.features {
            total = total
                .checked_add(visit(feature.base.unique_id))
                .ok_or_else(limit)?;
        }
        Ok(total)
    }

    /// Observation matches in `graph` that no consensus feature references.
    ///
    /// Ports `ConsensusMap::getUnassignedIDMatches`; see
    /// [`FeatureMap::unassigned_id_matches`] for why the graph is a parameter.
    ///
    /// # Errors
    ///
    /// As [`FeatureMap::unassigned_id_matches`].
    pub fn unassigned_id_matches(
        &self,
        graph: &IdentificationData,
    ) -> Result<BTreeSet<ObservationMatchId>> {
        if self.features.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        unassigned_matches(
            graph,
            self.features.iter().map(|feature| &feature.base.id_matches),
        )
    }

    /// Translate the identification references of every consensus feature into
    /// the graph the translator describes.
    ///
    /// The source performs this inside its copy constructor, `appendRows` and
    /// `appendColumns`, after merging the embedded `IdentificationData`; here
    /// it is an explicit call, as for [`FeatureMap::update_id_references`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the translator has no entry for one
    /// of the referenced records. Every reference is verified in a first,
    /// read-only pass, so a missing translation leaves the map unchanged.
    pub fn update_id_references(&mut self, translator: &ReferenceTranslator) -> Result<()> {
        if self.features.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        for feature in &self.features {
            verify_base(&feature.base, translator)?;
        }
        for feature in &mut self.features {
            feature.base.update_id_references(translator)?;
        }
        Ok(())
    }

    /// Split this consensus map back into the feature maps it was linked from.
    ///
    /// Ports `ConsensusMap::split`. One [`FeatureMap`] is produced per column
    /// header, in column-index order, and every feature handle becomes a
    /// [`Feature`] in the map its `map_index` names. A handle contributes its
    /// retention time, m/z, intensity, charge and width; its unique ID is
    /// **not** carried, because the source constructs the feature through
    /// `BaseFeature(const FeatureHandle&)`, which slices the handle to its
    /// `Peak2D` base. Every produced map receives a copy of this map's
    /// data-processing records.
    ///
    /// Where the peptide identifications go depends on whether this map has
    /// been through `IsobaricAnalyzer`, which the source decides by looking for
    /// that software name in the data-processing records:
    ///
    /// * **Isobaric.** Every assigned identification of a consensus feature
    ///   goes to the feature with the smallest map index, and all unassigned
    ///   identifications, together with all protein identifications, go to the
    ///   first feature map. The source comment marks copying every protein
    ///   identification to the first map as wrong but deliberate; it is
    ///   preserved.
    /// * **Not isobaric.** Every identification, assigned and unassigned, must
    ///   carry a `map_index` meta value and is routed by it. Protein
    ///   identifications are not distributed at all, as in the source.
    ///
    /// `mode` decides what happens to the consensus feature's own meta values;
    /// see [`SplitMeta`].
    ///
    /// One further source quirk is preserved deliberately: an identification
    /// whose `map_index` names a column this consensus feature has no handle
    /// for still creates a feature there, because the source routes through
    /// `std::map::operator[]`, which default-inserts. That feature is a default
    /// `BaseFeature` carrying only the identification.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when a non-isobaric map has an
    /// identification without `map_index`, matching the source exception of the
    /// same name, and [`Error::InvalidValue`] when
    ///
    /// * an isobaric map has a consensus feature whose smallest map index is
    ///   not zero, or which has no handle at all (the source
    ///   `Exception::ElementNotFound`, and an unguarded dereference of an empty
    ///   map, respectively);
    /// * [`SplitMeta::CopyFirst`] is requested and the smallest *handle* map
    ///   index of a consensus feature is not zero (the source
    ///   `Exception::ElementNotFound`). The source takes that index before it
    ///   routes the identifications, so a feature created at index 0 by an
    ///   identification alone does not satisfy the check, and this port matches;
    /// * a map index does not name a column of this map. The source indexes its
    ///   result vector with the map index and reads out of bounds when the
    ///   headers are not keyed `0..n-1`;
    /// * this map has no columns but an isobaric map's records would need a
    ///   first feature map to go to;
    /// * the map exceeds [`ConsensusMap::MAX_ITEMS`] or
    ///   [`ConsensusMap::MAX_COLUMNS`].
    ///
    /// The feature maps are built independently of this map, which is only
    /// read, so a failure changes nothing.
    pub fn split(&self, mode: SplitMeta) -> Result<Vec<FeatureMap>> {
        if self.features.len() > Self::MAX_ITEMS {
            return Err(limit());
        }
        if self.column_headers.len() > Self::MAX_COLUMNS {
            return Err(limit());
        }
        let columns = self.column_headers.len();
        let mut maps = vec![FeatureMap::new(); columns];
        let isobaric = has_isobaric_analyzer(&self.data_processing);
        for consensus in &self.features {
            let mut features: BTreeMap<u64, BaseFeature> = BTreeMap::new();
            for handle in consensus.handles() {
                features.insert(handle.map_index, base_feature_from_handle(handle));
            }
            let min_index = features.keys().next().copied();
            if isobaric && min_index != Some(0) {
                return Err(invalid(
                    "map seems to have gone through IsobaricAnalyzer, but a consensus feature has \
                     no handle with map index 0",
                ));
            }
            for record in &consensus.base.peptide_identifications {
                let target = if isobaric {
                    min_index.ok_or_else(|| {
                        invalid("isobaric consensus feature has no handle to attach an ID to")
                    })?
                } else {
                    map_index_of(record)?.ok_or_else(|| {
                        Error::MissingInformation(
                            "map did not undergo IsobaricAnalyzer, but a peptide identification \
                             carries no map_index"
                                .into(),
                        )
                    })?
                };
                check_column(target, columns)?;
                features
                    .entry(target)
                    .or_default()
                    .peptide_identifications
                    .push(record.clone());
            }
            match mode {
                SplitMeta::Discard => {}
                SplitMeta::CopyAll => {
                    for feature in features.values_mut() {
                        feature.metadata = consensus.base.metadata.clone();
                    }
                }
                SplitMeta::CopyFirst => {
                    // The source checks the smallest map index of the HANDLES,
                    // taken before the identifications were routed, and then
                    // writes to the smallest key of the map that routing may
                    // have extended. Both halves are reproduced.
                    if min_index != Some(0) {
                        return Err(invalid(
                            "no feature with map index 0 to copy meta values to",
                        ));
                    }
                    if let Some(feature) = features.values_mut().next() {
                        feature.metadata = consensus.base.metadata.clone();
                    }
                }
            }
            for (index, feature) in features {
                let slot = check_column(index, columns)?;
                maps[slot].features.push(Feature::from(feature));
            }
        }
        if isobaric {
            if columns == 0 {
                return Err(invalid(
                    "isobaric consensus map has no columns to attach unassigned identifications to",
                ));
            }
            maps[0]
                .unassigned_peptide_identifications
                .clone_from(&self.unassigned_peptide_identifications);
            maps[0]
                .protein_identifications
                .clone_from(&self.protein_identifications);
        } else {
            for record in &self.unassigned_peptide_identifications {
                let index = map_index_of(record)?.ok_or_else(|| {
                    Error::MissingInformation(
                        "map did not undergo IsobaricAnalyzer, but an unassigned peptide \
                         identification carries no map_index"
                            .into(),
                    )
                })?;
                let slot = check_column(index, columns)?;
                maps[slot]
                    .unassigned_peptide_identifications
                    .push(record.clone());
            }
        }
        for map in &mut maps {
            map.data_processing.clone_from(&self.data_processing);
        }
        Ok(maps)
    }

    /// Shared preflight and reset for the two append operations.
    fn preflight_append(&self, rhs: &Self) -> Result<Self> {
        combine_limit(self.features.len(), rhs.features.len())?;
        combine_limit(
            self.protein_identifications.len(),
            rhs.protein_identifications.len(),
        )?;
        combine_limit(
            self.unassigned_peptide_identifications.len(),
            rhs.unassigned_peptide_identifications.len(),
        )?;
        combine_limit(self.data_processing.len(), rhs.data_processing.len())?;
        let columns = self
            .column_headers
            .len()
            .checked_add(rhs.column_headers.len())
            .ok_or_else(limit)?;
        if columns > Self::MAX_COLUMNS {
            return Err(limit());
        }
        let mut merged = self.clone();
        merged.identifier = String::new();
        merged.loaded_file_path = String::new();
        merged.loaded_file_type = FileType::Unknown;
        merged.unique_id = 0;
        merged
            .data_processing
            .extend(rhs.data_processing.iter().cloned());
        Ok(merged)
    }
}

/// Source `operator<<(std::ostream&, const ConsensusMap&)`: one `Map <index>:
/// <filename> - <label> - <size>` line per column header, then every consensus
/// feature followed by a blank line.
///
/// The per-feature block is the source `operator<<(std::ostream&, const
/// ConsensusFeature&)` layout, written here rather than through a `Display`
/// implementation on
/// [`ConsensusFeature`], which
/// `KERNEL/ConsensusFeature.h` — a different work package — has not yet gained.
/// Numbers use Rust formatting rather than the C++ stream's
/// six-significant-digit default, matching every other kernel `Display`.
impl fmt::Display for ConsensusMap {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (index, header) in &self.column_headers {
            writeln!(
                f,
                "Map {}: {} - {} - {}",
                index, header.filename, header.label, header.size
            )?;
        }
        for feature in &self.features {
            write_consensus_feature(f, feature)?;
            writeln!(f)?;
        }
        Ok(())
    }
}

/// Source `operator<<(std::ostream&, const ConsensusFeature&)`: a banner, the
/// position, intensity and quality, one indented block per grouped handle in
/// handle order, the meta values in key order, then a closing banner whose line
/// ends with a space before the newline.
fn write_consensus_feature(f: &mut fmt::Formatter<'_>, feature: &ConsensusFeature) -> fmt::Result {
    writeln!(f, "---------- CONSENSUS ELEMENT BEGIN -----------------")?;
    writeln!(f, "Position: {} {}", feature.rt, feature.mz)?;
    writeln!(f, "Intensity {}", feature.intensity)?;
    writeln!(f, "Quality {}", feature.quality)?;
    writeln!(f, "Grouped features: ")?;
    for handle in feature.handles() {
        writeln!(f, " - Map index: {}", handle.map_index)?;
        writeln!(f, "   Feature id: {}", handle.unique_id)?;
        writeln!(f, "   RT: {}", handle.rt)?;
        writeln!(f, "   m/z: {}", handle.mz)?;
        writeln!(f, "   Intensity: {}", handle.intensity)?;
    }
    writeln!(f, "Meta information: ")?;
    for (key, value) in &feature.base.metadata {
        writeln!(f, "   {key}: {value}")?;
    }
    writeln!(f, "---------- CONSENSUS ELEMENT END ----------------- ")
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

/// Reject a combination that would exceed the shared container ceiling.
fn combine_limit(left: usize, right: usize) -> Result<()> {
    let total = left.checked_add(right).ok_or_else(limit)?;
    if total > FeatureMap::MAX_ITEMS {
        return Err(limit());
    }
    Ok(())
}

/// Whether any non-zero unique ID occurs more than once. Zero is the source
/// `UniqueIdInterface::INVALID` and is not counted, so several unassigned
/// elements are not a conflict — which is why the source's
/// `updateUniqueIdToIndex` postcondition accepts a map of unassigned features.
fn has_duplicate_unique_ids(ids: impl Iterator<Item = u64>) -> bool {
    let mut seen = BTreeSet::new();
    ids.filter(|id| *id != 0).any(|id| !seen.insert(id))
}

/// Source `UniqueIdIndexer::resolveUniqueIdConflicts`, reachable only through
/// the container merges of the two headers ported here: an element without a
/// unique ID is given one, and an element whose ID is already taken is redrawn
/// until it is free. `CONCEPT/UniqueIdIndexer.h` itself is not part of this
/// work package, so this reproduces only the behaviour those merges depend on.
///
/// The returned count is the number of redraws caused by a collision, which is
/// what the source counts: IDs freshly assigned to previously unassigned
/// elements are not counted.
fn resolve_unique_id_conflicts<T>(
    items: &mut [T],
    id_of: impl Fn(&mut T) -> &mut u64,
    generator: &mut UniqueIdGenerator,
) -> Result<usize> {
    let mut seen: BTreeSet<u64> = BTreeSet::new();
    let mut replaced = 0usize;
    for item in items.iter_mut() {
        let id = id_of(item);
        id.ensure_unique_id(generator);
        let mut redraws = 0usize;
        while seen.contains(id) {
            id.assign_new_unique_id(generator);
            replaced = replaced.checked_add(1).ok_or_else(limit)?;
            redraws += 1;
            if redraws > MAX_UNIQUE_ID_REDRAWS {
                return Err(invalid(
                    "could not find a free unique ID; the generator keeps producing collisions",
                ));
            }
        }
        seen.insert(*id);
    }
    Ok(replaced)
}

/// Resolve unique-ID conflicts among consensus features, but only when there
/// are any, matching the source's `updateUniqueIdToIndex` postcondition check
/// before it falls back to `resolveUniqueIdConflicts`.
fn resolve_map_unique_ids(
    features: &mut [ConsensusFeature],
    generator: &mut UniqueIdGenerator,
) -> Result<usize> {
    if !has_duplicate_unique_ids(features.iter().map(|feature| feature.base.unique_id)) {
        return Ok(0);
    }
    resolve_unique_id_conflicts(features, |feature| &mut feature.base.unique_id, generator)
}

/// Source `appendRows`/`appendColumns` modification hygiene: sort and
/// deduplicate the fixed and variable modification lists of a search parameter
/// set. The source applies this to every protein identification of the merged
/// map, not only to the appended ones.
fn dedup_modifications(record: &mut ProteinIdentification) {
    for list in [
        &mut record.search_parameters.variable_modifications,
        &mut record.search_parameters.fixed_modifications,
    ] {
        list.sort();
        list.dedup();
    }
}

/// The `map_index` meta value of an identification, if it carries one.
fn map_index_of(record: &PeptideIdentification) -> Result<Option<u64>> {
    match record.metadata.get(MAP_INDEX) {
        None => Ok(None),
        Some(value) => u64::try_from(value.as_i64()?)
            .map(Some)
            .map_err(|_| invalid("map_index must be a nonnegative index")),
    }
}

/// Add `shift` to an identification's `map_index`, if it carries one.
fn shift_map_index(record: &mut PeptideIdentification, shift: u64) -> Result<()> {
    let Some(index) = map_index_of(record)? else {
        return Ok(());
    };
    let shifted = index
        .checked_add(shift)
        .ok_or_else(|| invalid("shifted map_index exceeds the map index range"))?;
    record
        .metadata
        .insert(MAP_INDEX.into(), MetaValue::try_from(shifted)?);
    Ok(())
}

/// Turn a map index into a position in the split result, rejecting the
/// out-of-range index the source reads past the end of its vector.
fn check_column(index: u64, columns: usize) -> Result<usize> {
    let slot = usize::try_from(index)
        .map_err(|_| invalid("map index does not name a consensus column"))?;
    if slot >= columns {
        return Err(invalid("map index does not name a consensus column"));
    }
    Ok(slot)
}

/// Source `BaseFeature(const FeatureHandle&)`: position, intensity, charge and
/// width are carried; quality starts at zero and the handle's unique ID is not
/// carried, because the source slices the handle to its `Peak2D` base.
fn base_feature_from_handle(handle: &FeatureHandle) -> BaseFeature {
    BaseFeature {
        rt: handle.rt,
        mz: handle.mz,
        intensity: handle.intensity,
        charge: handle.charge,
        width: handle.width,
        ..BaseFeature::default()
    }
}

/// Source `DataProcessingUtils::hasIsobaricAnalyzer`, the query
/// `ConsensusMap::split` dispatches on. `METADATA/DataProcessingUtils.h` is not
/// part of this work package, so the one query the split depends on is
/// reproduced privately rather than exposed.
fn has_isobaric_analyzer(processing: &[DataProcessing]) -> bool {
    processing
        .iter()
        .any(|record| record.software.name == ISOBARIC_ANALYZER)
}

/// Source `MSExperiment::getPrimaryMSRunPath` reduced to the question the two
/// `setPrimaryMSRunPath` overloads ask: is there exactly one annotated run, is
/// it an mzML, and does it exist on disk?
///
/// The source assembles each location from a source file's path and name and
/// chooses the separator from the path with any `file:///` prefix stripped, but
/// then joins using the *unstripped* path; that spelling is reproduced so the
/// existence test sees what the source tests. `KERNEL/MSExperiment.h` is a
/// different work package, so this stays private.
fn usable_experiment_run_path(experiment: &super::MSExperiment) -> Option<String> {
    let mut paths = Vec::new();
    for file in &experiment.settings.source_files {
        if file.path.is_empty() || file.name.is_empty() {
            continue;
        }
        let actual = file.path.strip_prefix("file:///").unwrap_or(&file.path);
        let separator = if actual.contains('\\') && !actual.contains('/') {
            "\\"
        } else {
            "/"
        };
        paths.push(format!("{}{}{}", file.path, separator, file.name));
    }
    if paths.len() != 1 {
        return None;
    }
    let path = paths.remove(0);
    if !path.ends_with("mzML") || !std::path::Path::new(&path).exists() {
        return None;
    }
    Some(path)
}

/// Source `getUnassignedIDMatches`: every match registered in `graph` that none
/// of the given match sets claims.
fn unassigned_matches<'a>(
    graph: &IdentificationData,
    claimed: impl Iterator<Item = &'a BTreeSet<ObservationMatchId>>,
) -> Result<BTreeSet<ObservationMatchId>> {
    let mut assigned: BTreeSet<ObservationMatchId> = BTreeSet::new();
    for matches in claimed {
        if matches.len() > BaseFeature::MAX_ID_MATCHES {
            return Err(limit());
        }
        assigned.extend(matches.iter().copied());
    }
    Ok(graph
        .observation_matches()
        .map(|(id, _)| id)
        .filter(|id| !assigned.contains(id))
        .collect())
}

/// Read-only check that every identification reference of one feature can be
/// translated, so the map-wide update can be atomic.
fn verify_base(feature: &BaseFeature, translator: &ReferenceTranslator) -> Result<()> {
    if feature.id_matches.len() > BaseFeature::MAX_ID_MATCHES {
        return Err(limit());
    }
    if let Some(id) = feature.primary_id {
        translator.molecule(id)?;
    }
    for &id in &feature.id_matches {
        translator.observation_match(id)?;
    }
    Ok(())
}

/// [`verify_base`] over a feature and its subordinates, under the same
/// traversal ceilings as `Feature::update_all_id_references`.
fn verify_feature_tree(feature: &Feature, translator: &ReferenceTranslator) -> Result<()> {
    let mut visited = 0usize;
    let mut pending: Vec<(&Feature, usize)> = vec![(feature, 0)];
    while let Some((current, depth)) = pending.pop() {
        if depth > Feature::MAX_SUBORDINATE_DEPTH {
            return Err(invalid("subordinate feature depth exceeds 128"));
        }
        visited = visited.checked_add(1).ok_or_else(limit)?;
        if visited > Feature::MAX_TRAVERSED_FEATURES {
            return Err(limit());
        }
        verify_base(&current.base, translator)?;
        pending.extend(current.subordinates.iter().map(|sub| (sub, depth + 1)));
    }
    Ok(())
}

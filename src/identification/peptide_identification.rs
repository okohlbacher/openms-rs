// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The `METADATA/PeptideIdentification.h` members beyond the record itself:
//! unique identifiers across merged runs, Universal Spectrum Identifiers, the
//! score comparator factory, the source `empty()` and equality that lets an
//! identification key a hash container (its `Hash` is in [`super::hits`]).
//!
//! The record and its accessors live in the parent module; this file holds
//! the members that need run mappings, consensus maps, ProForma or USI. Every
//! behaviour here was checked against the executed Linux x86_64 Release build
//! (oracle `B8b-peptide-id`); see `docs/IDENTIFICATION_SUPPORT.md`.

use super::run_mapping::IdentifierMSRunMapper;
use super::{PeptideHit, PeptideIdentification};
use crate::chemistry::proforma::{ChargeState, Peptidoform, PeptidoformIon, WriteMode};
use crate::data_structures::list::ListFormat;
use crate::kernel::ConsensusMap;
use crate::metadata::usi::{IndexType, Usi};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

/// Where one identification of a [`ConsensusMap`] sits: the consensus feature
/// index (`None` for an unassigned identification) and the position in that
/// feature's (or the unassigned) identification list.
///
/// The source multimap stores `std::pair<Size, Size>` and marks unassigned
/// identifications with the feature index `Size(-1)`; `None` replaces that
/// sentinel.
pub type ConsensusIdLocation = (Option<usize>, usize);

/// The source `std::multimap<std::string, std::pair<Size, Size>>` returned by
/// `buildUIDsFromAllPepIDs`: every UID in ascending byte order, each with its
/// locations in insertion order (map order, then the unassigned identifications).
pub type ConsensusUidMap = BTreeMap<String, Vec<ConsensusIdLocation>>;

/// Most identifications [`PeptideIdentification::build_uids_from_all`] indexes.
///
/// The source has no ceiling. This bounds the key strings the index retains,
/// one per identification, and is checked before any UID is built.
pub const MAX_UID_IDENTIFICATIONS: usize = 10_000_000;

impl PeptideIdentification {
    /// Source static `getScoreComparator(bool higher_score_better)`: the strict
    /// "ranks before" predicate [`PeptideIdentification::sort`] orders by.
    ///
    /// With `higher_score_better` it is the source `PeptideHit::ScoreMore`
    /// (`a.score > b.score`), otherwise `PeptideHit::ScoreLess`
    /// (`a.score < b.score`). Equal scores compare false both ways, so a stable
    /// sort keeps their input order; a NaN score also compares false both ways,
    /// which is why [`PeptideIdentification::sort`] refuses non-finite scores.
    pub fn score_comparator(higher_score_better: bool) -> fn(&PeptideHit, &PeptideHit) -> bool {
        fn more(a: &PeptideHit, b: &PeptideHit) -> bool {
            a.score > b.score
        }
        fn less(a: &PeptideHit, b: &PeptideHit) -> bool {
            a.score < b.score
        }
        if higher_score_better { more } else { less }
    }

    /// Source `empty()`: no identifier, no hits, significance threshold 0, no
    /// score type and the default higher-is-better orientation.
    ///
    /// RT, m/z and other metadata are not consulted, exactly as in the source.
    /// This is not [`PeptideIdentification::is_empty`], which (as the Rust
    /// counterpart of `len`) looks at the candidate hits only.
    pub fn source_empty(&self) -> bool {
        self.identifier.is_empty()
            && self.hits.is_empty()
            && self.significance_threshold == 0.0
            && self.score_type.is_empty()
            && self.higher_score_better
    }

    /// Source static `getReferencingHits(hits, accessions)`: the hits of
    /// `hits` with at least one peptide evidence whose protein accession is in
    /// `accessions`, in input order. The source returns copies; these are
    /// borrows. [`PeptideIdentification::referencing_hits`] is the same query
    /// over this identification's own hits.
    pub fn hits_referencing<'a>(
        hits: &'a [PeptideHit],
        accessions: &BTreeSet<String>,
    ) -> Vec<&'a PeptideHit> {
        hits.iter()
            .filter(|hit| {
                hit.evidences
                    .iter()
                    .any(|evidence| accessions.contains(&evidence.protein_accession))
            })
            .collect()
    }

    /// Source static `buildUIDFromPepID(pep_id, mapping)`: a key that tells
    /// identifications of merged runs apart.
    ///
    /// When the run this identification belongs to ([`Self::identifier`]) has
    /// exactly one MS run path in `mapping`, the UID is that path, `|` and the
    /// spectrum reference. Otherwise, when a `map_index` metadata value
    /// exists, it is that value's source text form (`3`, `x7`, `2.5`), `|` and
    /// the spectrum reference; an identifier missing from `mapping` takes this
    /// branch as well, since it has no paths. The spectrum reference is
    /// [`Self::spectrum_reference`], so an integer reference is written as a
    /// number and an empty one leaves nothing after the `|`.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] with the source messages when no
    /// `spectrum_reference` metadata value exists ("Spectrum reference missing
    /// at PeptideIdentification.") or when the run does not have exactly one
    /// path and there is no `map_index` ("Multiple files in a run, but no
    /// map_index in PeptideIdentification found."). The source throws
    /// `Exception::MissingInformation` in both cases.
    pub fn build_uid(&self, mapping: &IdentifierMSRunMapper) -> Result<String> {
        if !self.metadata.contains_key("spectrum_reference") {
            return Err(Error::MissingInformation(
                "Spectrum reference missing at PeptideIdentification.".into(),
            ));
        }
        let paths = mapping.ms_run_paths(&self.identifier);
        let head = if paths.len() == 1 {
            paths[0].clone()
        } else if let Some(map_index) = self.metadata.get("map_index") {
            map_index.to_list_text()?.into_owned()
        } else {
            return Err(Error::MissingInformation(
                "Multiple files in a run, but no map_index in PeptideIdentification found.".into(),
            ));
        };
        Ok(format!("{head}|{}", self.spectrum_reference()))
    }

    /// Source static `buildUIDsFromAllPepIDs(cmap)`: the UID of every
    /// identification in a consensus map, mapped to where it sits.
    ///
    /// The run mapping is built from the map's protein identifications, then
    /// every consensus feature's identifications are visited in map order and
    /// the unassigned identifications last, as in the source. Each UID is
    /// [`PeptideIdentification::build_uid`]; equal UIDs keep their locations in
    /// that visiting order, which is the source multimap's insertion order.
    ///
    /// # Errors
    ///
    /// Whatever [`IdentifierMSRunMapper::from_runs`] or
    /// [`PeptideIdentification::build_uid`] reports, and
    /// [`Error::InvalidValue`] when the map holds more than
    /// [`MAX_UID_IDENTIFICATIONS`] identifications. Nothing is returned on
    /// error; the source would have thrown with a partly built multimap that
    /// its caller never sees either.
    pub fn build_uids_from_all(map: &ConsensusMap) -> Result<ConsensusUidMap> {
        let total = map
            .features
            .iter()
            .map(|feature| feature.base.peptide_identifications.len())
            .try_fold(map.unassigned_peptide_identifications.len(), |sum, n| {
                sum.checked_add(n)
            })
            .filter(|&total| total <= MAX_UID_IDENTIFICATIONS)
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "consensus map holds more than {MAX_UID_IDENTIFICATIONS} peptide identifications"
                ))
            })?;
        let mapping = IdentifierMSRunMapper::from_runs(&map.protein_identifications)?;
        let mut uids = ConsensusUidMap::new();
        let mut visited = 0_usize;
        let features = map
            .features
            .iter()
            .enumerate()
            .map(|(index, feature)| (Some(index), &feature.base.peptide_identifications));
        let unassigned = std::iter::once((None, &map.unassigned_peptide_identifications));
        for (feature, identifications) in features.chain(unassigned) {
            for (position, identification) in identifications.iter().enumerate() {
                let uid = identification.build_uid(&mapping)?;
                uids.entry(uid).or_default().push((feature, position));
                visited += 1;
            }
        }
        debug_assert_eq!(visited, total);
        Ok(uids)
    }

    /// Source `buildUSI(ms_run_name, dataset_id = "local",
    /// include_interpretation = false)`: the Universal Spectrum Identifier of
    /// this identification in a single-run context.
    ///
    /// The scan number comes from [`Self::spectrum_reference`] through
    /// [`Usi::extract_scan_number_from_native_id`]. With a scan number the USI
    /// has index type `scan`; without one but with a non-empty reference it is
    /// `nativeId` carrying the whole reference; with an empty reference the
    /// result is the default, invalid [`Usi`]. Pass [`Usi::LOCAL_COLLECTION`]
    /// for the source's default `dataset_id` of unpublished data.
    ///
    /// With `include_interpretation` and at least one hit, the **first** hit
    /// (call [`Self::sort`] first to make it the best one) becomes the
    /// interpretation: its sequence as a lossless ProForma peptidoform ion,
    /// with `/charge` appended unless the charge is 0 (`EM[UNIMOD:35]K/2`,
    /// `[UNIMOD:1]-PEPTIDEK`, `PEPTIDEK-[UNIMOD:2]/-2`).
    ///
    /// The source notes that this overload assumes one MS run; for merged
    /// files use [`PeptideIdentification::build_usi_with_mapping`].
    ///
    /// # Errors
    ///
    /// Whatever ProForma conversion and serialization of the first hit's
    /// sequence reports; the source would throw from the same step.
    pub fn build_usi(
        &self,
        ms_run_name: &str,
        dataset_id: &str,
        include_interpretation: bool,
    ) -> Result<Usi> {
        let reference = self.spectrum_reference();
        let scan = Usi::extract_scan_number_from_native_id(&reference)?;
        let mut interpretation = String::new();
        if include_interpretation {
            if let Some(first) = self.hits.first() {
                let ion = PeptidoformIon {
                    chains: vec![Peptidoform::from_aa_sequence(&first.sequence)?],
                    charge: (first.charge != 0).then_some(ChargeState::Simple(first.charge)),
                    ..PeptidoformIon::default()
                };
                interpretation = ion.to_text(WriteMode::Lossless)?;
            }
        }
        Ok(match scan {
            Some(scan) => Usi::new(
                dataset_id,
                ms_run_name,
                IndexType::Scan,
                scan.to_string(),
                interpretation,
            ),
            None if !reference.is_empty() => Usi::new(
                dataset_id,
                ms_run_name,
                IndexType::NativeId,
                reference,
                interpretation,
            ),
            None => Usi::default(),
        })
    }

    /// Source `buildUSI(mapping, dataset_id = "local",
    /// include_interpretation = false)`: [`PeptideIdentification::build_usi`]
    /// with the MS run resolved for merged workflows.
    ///
    /// The run file is [`IdentifierMSRunMapper::primary_ms_run_path`]: the
    /// path at this identification's `id_merge_index` (default 0) in its run,
    /// else the legacy `base_name` metadata value. An empty result gives the
    /// default, invalid [`Usi`], as does an index beyond the run's paths
    /// without a `base_name`. The USI carries the path's basename
    /// ([`Usi::extract_basename`]), so `file:///C:/data/fileB.mzML` becomes
    /// `fileB.mzML`.
    ///
    /// # Errors
    ///
    /// Whatever [`IdentifierMSRunMapper::primary_ms_run_path`] reports - a
    /// negative or non-integer `id_merge_index` in a mapped run, where the
    /// source throws `Exception::ConversionError` - and whatever
    /// [`PeptideIdentification::build_usi`] reports.
    pub fn build_usi_with_mapping(
        &self,
        mapping: &IdentifierMSRunMapper,
        dataset_id: &str,
        include_interpretation: bool,
    ) -> Result<Usi> {
        let path = mapping.primary_ms_run_path(self)?;
        if path.is_empty() {
            return Ok(Usi::default());
        }
        self.build_usi(
            &Usi::extract_basename(&path),
            dataset_id,
            include_interpretation,
        )
    }
}

/// Identifications can key a `HashSet` or `HashMap`, as the source's
/// `std::unordered_set<PeptideIdentification>` does; the `Hash` impl is in
/// [`super::hits`].
///
/// Equality is reflexive for every value [`PeptideIdentification::validate`]
/// accepts, because it refuses NaN scores, coordinates and thresholds. A NaN
/// placed in a public field by hand makes the value unequal to itself, which a
/// hash container then treats as a distinct key; nothing unsafe follows.
impl Eq for PeptideIdentification {}

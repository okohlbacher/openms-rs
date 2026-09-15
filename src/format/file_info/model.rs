// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The FileInfo options and structured result (`FORMAT/FileInfo.h`).
//!
//! These are the nested value types of the source class `OpenMS::FileInfo`:
//! the options a run takes (the command-line flags of the FileInfo tool) and
//! the aggregate it returns, which holds the extracted file-level information
//! and caches the text and TSV reports. The computation and the report live in
//! [`crate::format::file_info::report`]; the per-branch work in
//! [`crate::format::file_info::peaks`] and [`crate::format::file_info::features`].
//! `docs/FILE_INFO_SUPPORT.md` holds the API mapping, the preserved source
//! conventions, the native differences and the evidence.
//!
//! The source declares every type, including those of branches this port does
//! not run yet (consensus maps, identifications, FASTA, mzTab, validation,
//! corrupt-data checks and detailed listings). They are ported here with their
//! fields so the result keeps the source shape; a run that would fill them
//! returns [`crate::Error::Unsupported`] instead of an incompletely filled
//! result.
//!
//! Several declared fields are never written by the source `run` at core
//! `bc9cc12` either:
//! [`FileInfoResult::experiment_meta`](crate::format::file_info::model::FileInfoResult::experiment_meta),
//! [`FileInfoResult::statistics`](crate::format::file_info::model::FileInfoResult::statistics),
//! [`FileInfoResult::corruption`](crate::format::file_info::model::FileInfoResult::corruption),
//! [`FileInfoResult::detail`](crate::format::file_info::model::FileInfoResult::detail),
//! the validation `schema_version` and `detail`, and, for a feature map, every
//! consensus-only field. The port leaves them at their defaults in the same
//! places.

use crate::concept::progress_logger::ProgressLogType;
use crate::format::FileType;
use crate::math::statistic_functions::SummaryStatistics;
use std::collections::BTreeMap;

/// One `[min, max]` interval of a single dimension, the source
/// `FileInfo::Range` when its `present` flag is set.
///
/// The source struct carries `bool present` next to `min` and `max`, both zero
/// when absent. Here an absent range is `None` in the containing
/// [`RangeSet`], so a present range and a zero-valued range cannot be
/// confused.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Range {
    /// Smallest value of the dimension.
    pub min: f64,
    /// Largest value of the dimension.
    pub max: f64,
}

/// The four range dimensions FileInfo reports for a map or a spectrum group,
/// the source `FileInfo::RangeSet`.
///
/// A dimension without any value is `None`; the report prints it as `<none>`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RangeSet {
    /// Retention time in seconds.
    pub rt: Option<Range>,
    /// Mass-to-charge ratio.
    pub mz: Option<Range>,
    /// Ion mobility; only filled when the source range manager has a mobility
    /// dimension and a value in it.
    pub mobility: Option<Range>,
    /// Intensity.
    pub intensity: Option<Range>,
    /// Whether the mobility dimension applies at all: `true` for the combined,
    /// overall-spectrum and per-MS-level ranges of an experiment, `false` for
    /// chromatogram ranges and feature maps, whose source range managers have
    /// no mobility dimension.
    pub has_mobility: bool,
}

/// The range categories of a report, the source `FileInfo::Ranges`.
///
/// An experiment fills all four categories; a feature map fills only
/// [`Ranges::combined`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Ranges {
    /// Spectra and chromatograms of an experiment, or the whole map.
    pub combined: RangeSet,
    /// All spectra of an experiment; experiments only.
    pub spectra_overall: RangeSet,
    /// Spectra of one MS level, keyed by that level; experiments only.
    pub per_ms_level: BTreeMap<u32, RangeSet>,
    /// All chromatograms of an experiment; experiments only.
    pub chromatograms: RangeSet,
    /// `true` when the per-level and chromatogram categories apply.
    pub is_experiment: bool,
}

/// The general header block, always populated, the source `FileInfo::FileMeta`.
#[derive(Clone, Debug, PartialEq)]
pub struct FileMeta {
    /// The file name exactly as passed to the run.
    pub file_name: String,
    /// The forced or detected file type; [`FileType::Unknown`] when neither
    /// gave one.
    pub file_type: FileType,
    /// [`FileType::name`] of [`FileMeta::file_type`], the source
    /// `FileTypes::typeToName`.
    pub file_type_name: String,
}

impl Default for FileMeta {
    /// An empty name, [`FileType::Unknown`] and an empty type name, as the
    /// source member initialisers.
    fn default() -> Self {
        Self {
            file_name: String::new(),
            file_type: FileType::Unknown,
            file_type_name: String::new(),
        }
    }
}

/// A contact person of [`ExperimentMeta`], the source
/// `FileInfo::ExperimentMeta::Contact`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Contact {
    /// First name.
    pub first_name: String,
    /// Last name.
    pub last_name: String,
    /// Email address.
    pub email: String,
}

/// Experiment, instrument, sample and contact metadata of a peak file (the
/// `-m` block), the source `FileInfo::ExperimentMeta`.
///
/// The source declares this aggregate and the `Result::experiment_meta` slot,
/// but its `run` at core `bc9cc12` never fills them: the `-m` block is written
/// only to the text and TSV reports. The port keeps that, so
/// [`FileInfoResult::experiment_meta`] is always `None`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ExperimentMeta {
    /// Whether the block was filled.
    pub present: bool,
    /// Document identifier.
    pub document_id: String,
    /// Date as `DateTime::get` text.
    pub date: String,
    /// Sample name.
    pub sample_name: String,
    /// Sample organism.
    pub sample_organism: String,
    /// Sample comment.
    pub sample_comment: String,
    /// Instrument name.
    pub instrument_name: String,
    /// Instrument model.
    pub instrument_model: String,
    /// Instrument vendor.
    pub instrument_vendor: String,
    /// Ionization-method names of the ion sources.
    pub ion_sources: Vec<String>,
    /// Analyzer-type names of the mass analyzers.
    pub mass_analyzers: Vec<String>,
    /// Detector-type names of the ion detectors.
    pub detectors: Vec<String>,
    /// Contact persons.
    pub contacts: Vec<Contact>,
}

/// One data-processing step (the `-p` block), the source
/// `FileInfo::ProcessingStep`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ProcessingStep {
    /// Software name.
    pub software_name: String,
    /// Software version.
    pub software_version: String,
    /// Completion time as `DateTime::get` text, `0000-00-00 00:00:00` when
    /// unset.
    pub completion_time: String,
    /// Processing-action names in enum order.
    pub actions: Vec<String>,
}

/// One named summary-statistics block, the source `FileInfo::NamedStats`;
/// `title` is the label used in both renderers.
///
/// Declared by the source, never filled by its `run` at core `bc9cc12`: the
/// `-s` blocks are written only to the reports. [`FileInfoResult::statistics`]
/// therefore stays empty here too.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NamedStats {
    /// Block label.
    pub title: String,
    /// The statistics, `Math::SummaryStatistics<std::vector<double>>`.
    pub stats: SummaryStatistics,
}

/// Peak-file (`MSExperiment`) specifics, the source `FileInfo::PeakInfo`.
///
/// Integer widths follow the source: `Int` keys are `i32`, `UInt64` counts
/// are `u64`, and `std::map` keys iterate in the same order as the
/// [`BTreeMap`] keys, strings bytewise.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PeakInfo {
    /// Instrument name of the experiment settings.
    pub instrument_name: String,
    /// `(analyzer-type name, resolution)` per mass analyzer, in instrument
    /// order.
    pub mass_analyzers: Vec<(String, f64)>,
    /// Distinct MS levels, ascending.
    pub ms_levels: Vec<i32>,
    /// Peaks of all spectra plus points of all chromatograms, the source
    /// `MSExperiment::getSize`.
    pub total_peaks: u64,
    /// Number of spectra.
    pub num_spectra: u64,
    /// Spectra per MS level.
    pub spectra_per_ms_level: BTreeMap<i32, u64>,
    /// `"<annotated> (<estimated>)"` spectrum-type names per MS level, as the
    /// report prints them.
    pub peak_type_per_ms_level: BTreeMap<i32, String>,
    /// Precursor activation methods: `(MS level, full method name)` to count,
    /// over every precursor of every spectrum.
    pub activation_methods: BTreeMap<(i32, String), u64>,
    /// Charge of the first precursor of each spectrum that has one, to count.
    pub precursor_charges: BTreeMap<i32, u64>,
    /// Float data-array name to the number of arrays with that name.
    pub float_arrays: BTreeMap<String, u64>,
    /// Integer data-array name to the number of arrays with that name.
    pub int_arrays: BTreeMap<String, u64>,
    /// String data-array name to the number of arrays with that name.
    pub string_arrays: BTreeMap<String, u64>,
    /// Distinct FAIMS compensation voltages in volts, ascending.
    pub faims_cvs: Vec<f64>,
    /// Number of chromatograms.
    pub num_chromatograms: u64,
    /// Points of all chromatograms.
    pub num_chrom_peaks: u64,
    /// Chromatogram-type name to count, keyed bytewise by name.
    pub chromatogram_types: BTreeMap<String, u64>,
}

impl PeakInfo {
    /// The activation methods as `(MS level, full method name, count)` rows in
    /// map order, the source `activationMethodsFlat` (a flattened view that
    /// mirrors the TSV columns).
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::format::file_info::model::PeakInfo;
    ///
    /// let mut info = PeakInfo::default();
    /// info.activation_methods
    ///     .insert((2, "Collision-induced dissociation".to_owned()), 3);
    /// assert_eq!(
    ///     info.activation_methods_flat(),
    ///     [(2, "Collision-induced dissociation".to_owned(), 3)]
    /// );
    /// ```
    pub fn activation_methods_flat(&self) -> Vec<(i32, String, u64)> {
        self.activation_methods
            .iter()
            .map(|((level, method), count)| (*level, method.clone(), *count))
            .collect()
    }
}

/// One map column of a consensus map, the source
/// `FileInfo::FeatureInfo::MapColumn`. Filled by the consensusXML branch,
/// which this port does not run yet.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MapColumn {
    /// File name of the input map.
    pub filename: String,
    /// Identifier of the column.
    pub identifier: String,
    /// Label of the column.
    pub label: String,
    /// Number of features of the input map.
    pub size: u64,
}

/// Feature and consensus-map specifics, the source `FileInfo::FeatureInfo`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FeatureInfo {
    /// `true` for a consensus map.
    pub is_consensus: bool,
    /// Number of (consensus) features.
    pub num_features: u64,
    /// Total ion current of a featureXML map: the feature intensities summed
    /// as `double` in file order.
    pub tic: f64,
    /// Charge to number of features.
    pub charges: BTreeMap<i32, u64>,
    /// Number of peptide identifications to number of features carrying that
    /// many.
    pub ids_per_element: BTreeMap<u64, u64>,
    /// Peptide identifications attached to features.
    pub assigned_ids: u64,
    /// Peptide identifications not attached to any feature.
    pub unassigned_ids: u64,
    /// Consensus maps only: consensus-feature size to count.
    pub size_distribution: BTreeMap<u64, u64>,
    /// Consensus maps only: the map columns.
    pub map_columns: Vec<MapColumn>,
}

/// Identification specifics of idXML and mzIdentML files, the source
/// `FileInfo::IdentInfo`. Filled by the identification branch, which this
/// port does not run yet.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IdentInfo {
    /// Database name.
    pub db_name: String,
    /// Database version.
    pub db_version: String,
    /// Taxonomy.
    pub taxonomy: String,
    /// `"engine (version)"` per search engine.
    pub search_engines: Vec<String>,
    /// Number of identification runs.
    pub num_runs: u64,
    /// Protein hits.
    pub protein_hits: u64,
    /// Non-redundant protein hits.
    pub non_redundant_protein_hits: u64,
    /// Spectra with at least one peptide hit.
    pub matched_spectra: u64,
    /// Peptide hits.
    pub peptide_hits: u64,
    /// Peptide-spectrum matches per spectrum.
    pub psms_per_spectrum: f64,
    /// Average peptide length.
    pub avg_peptide_length: f64,
    /// Non-redundant peptides.
    pub non_redundant_peptides: u64,
    /// Modified top hits.
    pub modified_tophits: u64,
    /// Modification name to count.
    pub modification_counts: BTreeMap<String, u64>,
}

/// FASTA specifics, the source `FileInfo::FastaInfo`. Filled by the FASTA
/// branch, which this port does not run yet.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FastaInfo {
    /// Number of sequences.
    pub num_sequences: u64,
    /// Residues over all sequences.
    pub total_residues: u64,
    /// Whether the sequences were classified as nucleic acids.
    pub is_nucleic_acid: bool,
    /// Sequence-length statistics.
    pub length_stats: SummaryStatistics,
    /// Residue byte to count. The source key is `char`, whose signedness is
    /// platform-defined, so bytes of `0x80` and above may iterate first there.
    pub residue_counts: BTreeMap<u8, u64>,
    /// Sequences containing an ambiguous residue.
    pub seq_with_ambiguous: u64,
    /// Duplicate headers.
    pub dup_headers: u64,
    /// Duplicate sequences.
    pub dup_sequences: u64,
    /// Ambiguity bucket to count; the buckets depend on
    /// [`FastaInfo::is_nucleic_acid`].
    pub ambiguity_counts: BTreeMap<String, u64>,
}

/// mzTab specifics, the source `FileInfo::MzTabInfo`. Filled by the mzTab
/// branch, which this port does not run yet.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MzTabInfo {
    /// mzTab version.
    pub version: String,
    /// mzTab mode.
    pub mode: String,
    /// mzTab type, the source member `type` (a Rust keyword).
    pub kind: String,
    /// Peptide-spectrum matches.
    pub psms: u64,
    /// Peptides.
    pub peptides: u64,
    /// Proteins.
    pub proteins: u64,
    /// Oligonucleotides.
    pub oligonucleotides: u64,
    /// Oligonucleotide-spectrum matches.
    pub osms: u64,
    /// Small molecules.
    pub small_molecules: u64,
    /// Nucleic acids.
    pub nucleic_acids: u64,
}

/// The `-v` schema and semantic validation and `-i` indexed-mzML blocks, the
/// source `FileInfo::ValidationInfo`. Neither flag runs in this port yet.
///
/// The source `run` fills `performed`, `supported`, `valid`, `warnings`,
/// `errors` and the four index fields; it never writes `schema_version` or
/// `detail` at core `bc9cc12`, whose validator output goes only into the text.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ValidationInfo {
    /// Whether validation ran.
    pub performed: bool,
    /// Whether validation of this file type is supported.
    pub supported: bool,
    /// Whether the file was valid.
    pub valid: bool,
    /// Schema version validated against; never written by the source run.
    pub schema_version: String,
    /// Captured validator output, re-emitted verbatim; never written by the
    /// source run.
    pub detail: String,
    /// Semantic-validation warnings.
    pub warnings: Vec<String>,
    /// Semantic-validation errors.
    pub errors: Vec<String>,
    /// Whether the mzML index was checked.
    pub index_checked: bool,
    /// Whether a valid index was found.
    pub index_valid: bool,
    /// Spectra in the index.
    pub indexed_spectra: u64,
    /// Chromatograms in the index.
    pub indexed_chromatograms: u64,
}

impl Default for ValidationInfo {
    /// Nothing performed, `supported` set, as the source member initialisers.
    fn default() -> Self {
        Self {
            performed: false,
            supported: true,
            valid: false,
            schema_version: String::new(),
            detail: String::new(),
            warnings: Vec::new(),
            errors: Vec::new(),
            index_checked: false,
            index_valid: false,
            indexed_spectra: 0,
            indexed_chromatograms: 0,
        }
    }
}

/// The `-c` corrupt-data block, the source `FileInfo::CorruptionInfo`, with
/// already formatted message lines.
///
/// Declared by the source, never filled by its `run` at core `bc9cc12`: the
/// `-c` messages go only into the text report. [`FileInfoResult::corruption`]
/// therefore stays at its default here too.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CorruptionInfo {
    /// Whether the check ran.
    pub performed: bool,
    /// Error lines.
    pub errors: Vec<String>,
    /// Warning lines.
    pub warnings: Vec<String>,
}

/// The `-d` detailed per-spectrum listing, the source `FileInfo::DetailInfo`,
/// kept as rendered lines.
///
/// Declared by the source, never filled by its `run` at core `bc9cc12`: the
/// listing goes only into the text report. [`FileInfoResult::detail`]
/// therefore stays at its default here too.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DetailInfo {
    /// Whether the listing ran.
    pub performed: bool,
    /// Rendered lines.
    pub lines: Vec<String>,
}

/// The master result aggregate, the source `FileInfo::Result`.
///
/// Named `FileInfoResult` because `Result` is the crate's error-carrying
/// alias. Exactly one branch aggregate among [`FileInfoResult::peak`] and
/// [`FileInfoResult::feature`] is set after a successful run on a peak file or
/// featureXML map; none is set for an unknown type.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FileInfoResult {
    /// General header information, always populated.
    pub meta: FileMeta,
    /// Ranges of peak files and feature maps.
    pub ranges: Ranges,
    /// Peak-file specifics.
    pub peak: Option<PeakInfo>,
    /// Feature-map specifics.
    pub feature: Option<FeatureInfo>,
    /// Identification specifics; never set by this port yet.
    pub ident: Option<IdentInfo>,
    /// FASTA specifics; never set by this port yet.
    pub fasta: Option<FastaInfo>,
    /// mzTab specifics; never set by this port yet.
    pub mztab: Option<MzTabInfo>,
    /// Always `None`, as the source leaves it; see [`ExperimentMeta`].
    pub experiment_meta: Option<ExperimentMeta>,
    /// The data-processing steps of the `-p` block, filled only when `-p` was
    /// requested.
    pub processing: Vec<ProcessingStep>,
    /// Always empty, as the source leaves it; see [`NamedStats`].
    pub statistics: Vec<NamedStats>,
    /// Validation and index-check results; `-v` and `-i` do not run in this
    /// port yet, so it stays at its default.
    pub validation: ValidationInfo,
    /// Always at its default, as the source leaves it; see [`CorruptionInfo`].
    pub corruption: CorruptionInfo,
    /// Always at its default, as the source leaves it; see [`DetailInfo`].
    pub detail: DetailInfo,
    /// trafoXML model and summary; never set by this port yet.
    pub transformation_summary: String,
    /// PQP summary; never set by this port yet.
    pub targeted_summary: String,
    /// The human-readable report, identical to the FileInfo tool's `-out`
    /// output. Empty for an unknown file type.
    pub text: String,
    /// The TSV report, identical to the FileInfo tool's `-out_tsv` output.
    /// Empty for an unknown file type.
    pub tsv: String,
    /// Messages the source writes with `OPENMS_LOG_WARN` during the run
    /// instead of into the reports, in order: currently only
    /// `FAIMSHelper`'s missing-voltage warning, which the source logs twice
    /// (it asks for the voltages twice) and this field holds once. Native
    /// field: the library prints nothing, and a caller such as the FileInfo
    /// tool decides where the messages go.
    pub warnings: Vec<String>,
}

/// What to compute, mirroring the command-line flags, the source
/// `FileInfo::Options`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Options {
    /// A forced input type, the tool's `-in_type`; [`FileType::Unknown`]
    /// detects the type from the file name and then the content.
    pub forced_type: FileType,
    /// `-m`: metadata.
    pub meta: bool,
    /// `-p`: data processing.
    pub processing: bool,
    /// `-s`: summary statistics.
    pub statistics: bool,
    /// `-d`: detailed listing. Not ported yet.
    pub detailed: bool,
    /// `-c`: corrupt-data check. Not ported yet.
    pub check_corrupt: bool,
    /// `-v`: schema and semantic validation. Not ported yet.
    pub validate: bool,
    /// `-i`: indexed-mzML check. Not ported yet.
    pub check_index: bool,
    /// Progress-logging verbosity forwarded to the loaders, the tool's
    /// `-no_progress`; [`ProgressLogType::None`] keeps the library silent.
    /// The native loaders on this path report no progress, so the value has no
    /// effect.
    pub log_type: ProgressLogType,
    /// Read a dangling mzML header reference the way the source loader does.
    ///
    /// Native field; the source `Options` has none, because its mzML reader is
    /// always lenient. `false`, the default, keeps the strict reader, which
    /// refuses an mzML file whose `softwareRef`, `dataProcessingRef` or
    /// `defaultDataProcessingRef` names no definition with [`crate::Error::Parse`].
    /// `true` passes `crate::format::mzml::ReadOptions::source_dangling_references`
    /// to the mzML reader of the peak-file branch, which drops the reference as
    /// the source does and warns once per dangling ID on the crate's warning log
    /// stream. The FileInfo tool sets it, so its mzML loading matches the
    /// source's `FileHandler::loadExperiment` (decision D10 of the early TOPP
    /// bundle). The other readers ignore it, and so does a build without the
    /// `mzml` feature, which reads no mzML.
    pub source_dangling_references: bool,
}

impl Default for Options {
    /// Every flag off, no forced type and no progress logging, as the source
    /// member initialisers, and the strict mzML reader.
    fn default() -> Self {
        Self {
            forced_type: FileType::Unknown,
            meta: false,
            processing: false,
            statistics: false,
            detailed: false,
            check_corrupt: false,
            validate: false,
            check_index: false,
            log_type: ProgressLogType::None,
            source_dangling_references: false,
        }
    }
}

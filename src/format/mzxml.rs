// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded mzXML 3.1 reader and writer: `FORMAT/MzXMLFile.h` with
//! `FORMAT/HANDLERS/MzXMLHandler.h`. See `docs/MZXML_SUPPORT.md`.
//!
//! mzXML is mzML's predecessor and differs from it in five ways this module has
//! to handle. Its peak arrays interleave m/z and intensity as **pairs** in one
//! base64 payload instead of two separate arrays; the element width comes from a
//! `precision` attribute ("32"/"64") rather than a CV term; the payload is
//! **big-endian** (`byteOrder="network"`); compression is named by a
//! `compressionType` attribute; and MS2 scans are XML **children** of their MS1
//! parent rather than siblings. Nesting carries no data: a nested `<scan>`
//! becomes an ordinary spectrum in document order, exactly as upstream does.
//!
//! Entry points are
//! [`load`](crate::format::mzxml::load) /
//! [`store`](crate::format::mzxml::store) for paths,
//! [`read`](crate::format::mzxml::read) /
//! [`write`](crate::format::mzxml::write) for streams, and
//! [`MzXMLFile`](crate::format::mzxml::MzXMLFile) for the source class's
//! options-carrying adapter shape.
//!
//! The source file derives from `ProgressLogger` and hands itself to its
//! handler; [`load_with_progress`](crate::format::mzxml::load_with_progress)
//! and [`store_with_progress`](crate::format::mzxml::store_with_progress) make
//! the handler's calls on a caller's logger, and every other entry point runs
//! the same code and reports nothing.

use super::PeakFileOptions;
use super::path_io;
use crate::concept::progress_logger::{ProgressLogger, ProgressReporter, progress_value};
use crate::interfaces::MSDataConsumer;
use crate::kernel::{MSExperiment, MSSpectrum, NumericRange, Peak1D, Precursor};
use crate::metadata::{
    ActivationMethod, AnalyzerType, ChecksumType, ContactPerson, DataProcessing, DetectorType,
    ExperimentalSettings, Instrument, IonDetector, IonSource, IonizationMethod, MassAnalyzer,
    MetaInfo, MetaValue, MetaValueData, Polarity, ProcessingAction, ResolutionMethod, ScanMode,
    ScanWindow, SourceFile,
};
use crate::{Error, Result};
use base64::{
    Engine, alphabet,
    engine::{DecodePaddingMode, GeneralPurpose, GeneralPurposeConfig, general_purpose::STANDARD},
};
use flate2::{Compression, Decompress, FlushDecompress, Status, write::ZlibEncoder};
use quick_xml::{
    Reader,
    events::{BytesStart, Event},
};
use std::io::{BufRead, Write};
use std::ops::ControlFlow;
use std::path::Path;
use std::sync::Arc;

// ---------------------------------------------------------------------------
// Source constants
// ---------------------------------------------------------------------------

/// Schema resource the source `MzXMLFile` constructor registers
/// (`MzXMLFile.cpp:18`). The crate carries it, and the three schemas it
/// includes, unchanged, and validates against them with the `xml-schema`
/// feature (`MzXMLFile::is_valid`); see `docs/XML_SCHEMA_SUPPORT.md`.
pub const SCHEMA: &str = "/SCHEMAS/mzXML_idx_3.1.xsd";
/// Schema version the source `MzXMLFile` constructor registers
/// (`MzXMLFile.cpp:18`), and the version the writer emits.
pub const SCHEMA_VERSION: &str = "3.1";
/// Default namespace the writer emits, `MzXMLHandler.cpp:644`.
pub const NAMESPACE: &str = "http://sashimi.sourceforge.net/schema_revision/mzXML_3.1";

/// Metadata key holding a `<comment>` child of `<scan>` or `<msInstrument>`.
///
/// The source stores a scan comment in `SpectrumSettings::comment_` and an
/// instrument comment in a `#comment` meta value (`MzXMLHandler.cpp:597`,
/// `:605`). This port has no comment field on
/// [`MSSpectrum`], so both use the same `#comment`
/// key. Keys starting with `#` are internal and are not written as
/// `<nameValue>`, matching `MzXMLHandler::writeUserParam_`
/// (`MzXMLHandler.cpp:1141`).
pub const COMMENT_KEY: &str = "#comment";
/// Metadata key on a contact holding the `<operator phone="...">` attribute,
/// the source `#phone` (`MzXMLHandler.cpp:404`).
pub const PHONE_KEY: &str = "#phone";
/// Metadata key on a [`DataProcessing`]
/// holding the `<software type="...">` attribute, the source `#type`
/// (`MzXMLHandler.cpp:155`).
pub const PROCESSING_TYPE_KEY: &str = "#type";
/// Metadata key on a [`DataProcessing`]
/// holding `<dataProcessing intensityCutoff="...">`, the source
/// `#intensity_cutoff` (`MzXMLHandler.cpp:460`).
pub const INTENSITY_CUTOFF_KEY: &str = "#intensity_cutoff";
/// Spectrum metadata key holding `<scan filterLine="...">`, the source
/// `filter string` (`MzXMLHandler.cpp:325`), which is mzML's `MS:1000512`.
pub const FILTER_STRING_KEY: &str = "filter string";

/// Spectrum metadata keys the writer emits as `<scan>` attributes rather than
/// as `<nameValue>` children, so they are not written twice. The source writes
/// both (`MzXMLHandler.cpp:964-993` and `:1076`); see `docs/MZXML_SUPPORT.md`.
pub const ATTRIBUTE_METADATA_KEYS: [&str; 6] = [
    "lowest observed m/z",
    "highest observed m/z",
    "base peak m/z",
    "base peak intensity",
    "total ion current",
    FILTER_STRING_KEY,
];

// mzXML controlled-vocabulary tables, verbatim from `MzXMLHandler::init_`
// (`MzXMLHandler.cpp:1271-1294`). The index of a term in its table is the
// numeric value of the corresponding OpenMS enum, so the empty leading entry is
// the enum's `Unknown`. Source `cvStringToEnum_` takes the FIRST match, which is
// why the many later empty slots never win.
const POLARITY_TERMS: &str = "any;+;-";
const IONIZATION_TERMS: &str = ";ESI;EI;CI;FAB;;;;;;;;;;;;;APCI;;;NSI;;SELDI;;;MALDI";
const ANALYZER_TERMS: &str =
    ";Quadrupole;Quadrupole Ion Trap;;;TOF;Magnetic Sector;FT-ICR;;;;;;FTMS";
const DETECTOR_TERMS: &str = ";EMT;;;Faraday Cup;;;;;Channeltron;Daly;Microchannel plate";
const RESOLUTION_TERMS: &str = ";FWHM;TenPercentValley;Baseline";

/// Decoder for `<peaks>` payloads.
///
/// Upstream's `Base64::decode` keeps only as many whole bytes as the symbol
/// count yields and never inspects the surplus bits of the final symbol, so
/// mzXML in the wild carries payloads a canonical decoder refuses. The upstream
/// fixture `MzXMLFile_3_64bit.mzXML` is one: its first scan ends `AA1=`, whose
/// two surplus bits are set, and `MzXMLFile_test.cpp:313` asserts the resulting
/// intensity 100.0000991821289 against a 0.01 tolerance. Padding is likewise
/// treated as optional, because some converters omit it. Accepting these is
/// leniency about representation, not about values: the decoded bytes are
/// identical to what the source produces.
static PEAKS_BASE64: GeneralPurpose = GeneralPurpose::new(
    &alphabet::STANDARD,
    GeneralPurposeConfig::new()
        .with_decode_allow_trailing_bits(true)
        .with_decode_padding_mode(DecodePaddingMode::Indifferent),
);

fn term_index(table: &str, term: &str) -> Option<usize> {
    table.split(';').position(|entry| entry == term)
}
fn term_name(table: &str, index: usize) -> &str {
    table.split(';').nth(index).unwrap_or("")
}
fn enum_index<T: PartialEq + Copy>(all: &[T], value: T) -> usize {
    all.iter().position(|item| *item == value).unwrap_or(0)
}

// ---------------------------------------------------------------------------
// Errors and small helpers
// ---------------------------------------------------------------------------

fn invalid(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}
fn budget(message: &str) -> Error {
    Error::InvalidValue(format!("mzXML {message} limit exceeded"))
}
/// Valid XML 1.0 character content, as `crate::format::mzml` checks it.
fn xml_string(value: &str) -> Result<()> {
    if value.chars().all(|c| {
        matches!(c, '\t' | '\n' | '\r')
            || ('\u{20}'..='\u{d7ff}').contains(&c)
            || ('\u{e000}'..='\u{fffd}').contains(&c)
            || c >= '\u{10000}'
    }) {
        Ok(())
    } else {
        Err(Error::InvalidValue(
            "mzXML text contains invalid XML 1.0 characters".into(),
        ))
    }
}
/// Source `DRange::encloses` (`DRange.h:152-161`): closed below, open above.
/// `PeakFileOptions` ranges therefore exclude their own maximum.
fn encloses(range: NumericRange, value: f64) -> bool {
    !(value < range.min || value >= range.max)
}
fn finite(value: f64, label: &str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(format!("nonfinite {label}")))
    }
}

/// Source `operator<<(std::ostream&, double)` with `std::setprecision(digits)`:
/// `%g` semantics, trailing zeros removed, exponent at least two digits.
///
/// Used only when [`WriteOptions::significant_digits`] is `Some`. Nonfinite
/// input never reaches here; the writer's preflight rejects it.
fn general_format(value: f64, digits: usize) -> String {
    let digits = digits.max(1);
    if value == 0.0 {
        return if value.is_sign_negative() {
            "-0".into()
        } else {
            "0".into()
        };
    }
    let scientific = format!("{:.*e}", digits - 1, value);
    let (mantissa, exponent) = match scientific.split_once('e') {
        Some(parts) => parts,
        // `{:e}` always emits an exponent; this branch cannot be reached.
        None => return scientific,
    };
    let exponent: i32 = exponent.parse().unwrap_or(0);
    if exponent < -4 || exponent >= digits as i32 {
        let sign = if exponent < 0 { '-' } else { '+' };
        return format!("{}e{sign}{:02}", trim_zeros(mantissa), exponent.abs());
    }
    let decimals = usize::try_from(digits as i32 - 1 - exponent).unwrap_or(0);
    trim_zeros(&format!("{value:.decimals$}")).to_owned()
}
fn trim_zeros(text: &str) -> &str {
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.')
    } else {
        text
    }
}

// ---------------------------------------------------------------------------
// Read limits and options
// ---------------------------------------------------------------------------

/// Resource ceilings checked before any allocation from a file-derived length.
///
/// Every field bounds one quantity the file itself declares. A limit failure
/// returns [`Error::InvalidValue`] and leaves any caller-supplied destination
/// untouched, because a destination is replaced only after a complete parse.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Maximum decompressed XML bytes consumed, including base64 payloads.
    pub max_xml_bytes: u64,
    /// Maximum `<scan>` elements, counted before any filter is applied.
    pub max_scans: usize,
    /// Maximum `peaksCount` accepted for one scan, checked before reserving.
    pub max_peaks_per_scan: usize,
    /// Maximum declared peaks summed over the file, filtered scans included.
    pub max_total_peaks: usize,
    /// Maximum retained base64 characters for one `<peaks>` element.
    pub max_encoded_bytes: usize,
    /// Maximum decoded bytes for one `<peaks>` element, compression included.
    ///
    /// Charged three times: against the declared `peaksCount`, against the
    /// symbol count of the payload *before* the base64 decode allocates, and
    /// against the inflated length of a compressed payload.
    pub max_decoded_bytes: usize,
    /// Maximum open-element depth, which bounds nested `<scan>` recursion.
    pub max_depth: usize,
    /// Maximum `<parentFile>` elements.
    pub max_source_files: usize,
    /// Maximum `<dataProcessing>` elements.
    pub max_data_processing: usize,
    /// Maximum `<precursorMz>` elements per scan.
    pub max_precursors_per_scan: usize,
    /// Maximum metadata entries stored across the experiment and its spectra.
    pub max_metadata_entries: usize,
    /// Maximum characters retained for one non-`<peaks>` text node.
    pub max_text_bytes: usize,
    /// Maximum non-fatal diagnostics retained in a [`ReadReport`].
    pub max_diagnostics: usize,
}
impl Default for ReadLimits {
    fn default() -> Self {
        Self {
            max_xml_bytes: 512 * 1024 * 1024,
            max_scans: 1_000_000,
            max_peaks_per_scan: 10_000_000,
            max_total_peaks: 20_000_000,
            max_encoded_bytes: 128 * 1024 * 1024,
            max_decoded_bytes: 64 * 1024 * 1024,
            max_depth: 64,
            max_source_files: 100_000,
            max_data_processing: 100_000,
            max_precursors_per_scan: 10_000,
            max_metadata_entries: 1_000_000,
            max_text_bytes: 1024 * 1024,
            max_diagnostics: 1_000,
        }
    }
}

/// Scientific loading choices plus the independent resource ceilings.
///
/// `peaks` is the source `PeakFileOptions` the handler consults: `metadata_only`,
/// the RT/MS-level/precursor-m/z whole-scan filters, the m/z and intensity
/// per-peak filters, `fill_data`, `sort_spectra_by_mz` and `max_data_pool_size`.
/// The remaining `PeakFileOptions` fields are mzML- or writer-only and are
/// listed in `docs/MZXML_SUPPORT.md`.
#[derive(Clone, Debug, Default)]
pub struct ReadOptions {
    /// Source scientific loading options.
    pub peaks: PeakFileOptions,
    /// Native resource ceilings, with no source analogue.
    pub limits: ReadLimits,
}

/// Non-fatal observations and the raw scan count of a completed read.
///
/// The source logs these through `XMLHandler::error`/`warning`, which write to
/// the global log stream and never abort the load (`XMLHandler.cpp:71-109`).
/// This port returns them instead, so a caller can route or assert on them.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReadReport {
    /// `<scan>` elements seen, before filtering — the source `getScanCount()`.
    pub scan_count: usize,
    /// Messages, oldest first, truncated at [`ReadLimits::max_diagnostics`].
    pub diagnostics: Vec<String>,
    /// Whether [`ReadLimits::max_diagnostics`] discarded later messages.
    pub diagnostics_truncated: bool,
}
impl ReadReport {
    fn note(&mut self, limits: &ReadLimits, message: impl Into<String>) {
        if self.diagnostics.len() >= limits.max_diagnostics {
            self.diagnostics_truncated = true;
            return;
        }
        self.diagnostics.push(message.into());
    }
}

/// How much of a document a read materializes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Detail {
    /// Header and every selected scan.
    Full,
    /// Header only; parsing stops at the first `<scan>`, as source
    /// `MzXMLHandler` does by throwing `EndParsingSoftly`
    /// (`MzXMLHandler.cpp:242-245`).
    MetadataOnly,
    /// Header plus a scan count, storing no spectra — source
    /// `XMLHandler::LD_RAWCOUNTS` (`MzXMLHandler.cpp:286`).
    RawCounts,
}

// ---------------------------------------------------------------------------
// Write limits and options
// ---------------------------------------------------------------------------

/// Element width of the interleaved m/z–intensity array the writer emits.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PeakPrecision {
    /// `precision="64"`, which stores every `f64` m/z without loss.
    #[default]
    Float64,
    /// `precision="32"`, what the source writer always emits
    /// (`MzXMLHandler.cpp:1049`). m/z coordinates are narrowed to `f32`.
    Float32,
}
impl PeakPrecision {
    const fn bits(self) -> u32 {
        match self {
            Self::Float64 => 64,
            Self::Float32 => 32,
        }
    }
    const fn width(self) -> usize {
        match self {
            Self::Float64 => 8,
            Self::Float32 => 4,
        }
    }
}

/// Resource ceilings for one write, all checked before a byte is emitted.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WriteLimits {
    /// Maximum spectra written.
    pub max_spectra: usize,
    /// Maximum `<parentFile>` elements written.
    pub max_source_files: usize,
    /// Maximum peaks summed over all spectra.
    pub max_total_peaks: usize,
    /// Maximum precursors per spectrum.
    pub max_precursors_per_spectrum: usize,
    /// Maximum MS level, which bounds the indentation the source derives from
    /// it (`std::string(ms_level + 1, '\t')`, `MzXMLHandler.cpp:886`).
    pub max_ms_level: u32,
    /// Maximum metadata entries written for one record.
    pub max_metadata_entries: usize,
}
impl Default for WriteLimits {
    fn default() -> Self {
        Self {
            max_spectra: 1_000_000,
            max_source_files: 100_000,
            max_total_peaks: 100_000_000,
            max_precursors_per_spectrum: 10_000,
            max_ms_level: 64,
            max_metadata_entries: 1_000_000,
        }
    }
}

/// Writer configuration. [`WriteOptions::default`] is lossless;
/// [`WriteOptions::source`] reproduces `MzXMLHandler::writeTo`.
#[derive(Clone, Copy, Debug)]
pub struct WriteOptions {
    /// Stored element width. Defaults to [`PeakPrecision::Float64`] because the
    /// source's fixed `precision="32"` narrows every m/z coordinate.
    pub precision: PeakPrecision,
    /// Significant digits for numbers written as attributes or element text.
    /// `None`, the default, writes Rust's shortest round-tripping form.
    /// `Some(6)` reproduces the source's default `std::ostream` precision, which
    /// is where mzXML written by OpenMS loses m/z and retention-time digits.
    pub significant_digits: Option<usize>,
    /// Truncate `precursorIntensity` to an integer, as the source cast
    /// `(int)precursor.getIntensity()` does (`MzXMLHandler.cpp:1008`).
    pub integer_precursor_intensity: bool,
    /// Source `PeakFileOptions::force_mq_compatibility`: skip empty spectra,
    /// infer a Thermo `<msManufacturer>`, force `scanType="Full"` and a `CID`
    /// activation method, add `lowMz`/`highMz`/`basePeakIntensity`/
    /// `totIonCurrent`, break the `<peaks>` tag across lines and force an index.
    pub force_mq_compatibility: bool,
    /// Emit the `<index>`/`<indexOffset>` trailer, the source
    /// `PeakFileOptions::write_index` (default `true`).
    pub write_index: bool,
    /// zlib-compress each peak array and emit `compressionType="zlib"`. The
    /// source writer always emits `"none"`; its reader accepts both, so this is
    /// a native extension. Ignored under `force_mq_compatibility`.
    pub zlib_compression: bool,
    /// Native resource ceilings, with no source analogue.
    pub limits: WriteLimits,
}
impl Default for WriteOptions {
    fn default() -> Self {
        Self {
            precision: PeakPrecision::default(),
            significant_digits: None,
            integer_precursor_intensity: false,
            force_mq_compatibility: false,
            write_index: true,
            zlib_compression: false,
            limits: WriteLimits::default(),
        }
    }
}
impl WriteOptions {
    /// Exactly what `MzXMLHandler::writeTo` emits: 32-bit peaks, six
    /// significant digits, integer precursor intensities and no compression.
    ///
    /// This is lossy by construction; it exists so a caller can reproduce
    /// upstream output for tools that require it.
    pub fn source() -> Self {
        Self {
            precision: PeakPrecision::Float32,
            significant_digits: Some(6),
            integer_precursor_intensity: true,
            ..Self::default()
        }
    }
    /// Derive writer settings from a source `PeakFileOptions`, which carries the
    /// two mzXML writer switches `force_mq_compatibility` and `write_index`.
    pub fn from_peak_options(options: &PeakFileOptions) -> Self {
        Self {
            force_mq_compatibility: options.force_mq_compatibility,
            write_index: options.write_index,
            zlib_compression: options.zlib_compression,
            ..Self::default()
        }
    }
    fn number(&self, value: f64) -> String {
        match self.significant_digits {
            Some(digits) => general_format(value, digits),
            None => {
                let mut text = value.to_string();
                if text == "-0" {
                    text = "0".into();
                }
                text
            }
        }
    }
    fn meta_text(&self, value: &MetaValue) -> String {
        match value.data() {
            MetaValueData::Float(number) => self.number(*number),
            _ => value.to_string(),
        }
    }
}

// ---------------------------------------------------------------------------
// Transform options and report
// ---------------------------------------------------------------------------

/// Settings for a consumer-driven two-pass read, the source
/// `MzXMLFile::transform` (`MzXMLFile.cpp:59-89`).
#[derive(Clone, Debug, Default)]
pub struct TransformOptions {
    /// Source `skip_full_count`. When `true` the first pass stops at the first
    /// `<scan>` and reports an expected size of zero; when `false` it counts
    /// every scan without storing one.
    pub skip_full_count: bool,
    /// Loading options for both passes.
    pub read: ReadOptions,
}

/// Outcome of a consumer-driven read.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct TransformReport {
    /// Scans the first pass reported to the consumer through
    /// [`MSDataConsumer::set_expected_size`].
    pub expected_spectra: usize,
    /// Spectra delivered to the consumer, including one requesting a stop.
    pub delivered: usize,
    /// Whether a consumer returned [`ControlFlow::Break`] and ended the read.
    pub stopped: bool,
    /// Second-pass diagnostics and raw scan count.
    pub read: ReadReport,
}

// ---------------------------------------------------------------------------
// The file adapter
// ---------------------------------------------------------------------------

/// File adapter for mzXML 3.1, the source `MzXMLFile`.
///
/// The source class exposes its `PeakFileOptions` through
/// `getOptions`/`setOptions`; here the options are a public field, so
/// `file.options.peaks.add_ms_level(1)?` replaces
/// `file.getOptions().addMSLevel(1)`. Writer settings are a separate field
/// because this port's writer has lossless defaults the source does not.
///
/// The source also derives from `ProgressLogger`; progress reporting is not
/// ported (see `docs/MZXML_SUPPORT.md`).
#[derive(Clone, Debug, Default)]
pub struct MzXMLFile {
    /// Loading options, the source `options_` member.
    pub options: ReadOptions,
    /// Writing options, which the source derives from `options_` alone.
    pub write_options: WriteOptions,
}

impl MzXMLFile {
    /// A file adapter with source-default loading and lossless writing.
    pub fn new() -> Self {
        Self::default()
    }
    /// Schema version the source constructor registers, [`SCHEMA_VERSION`].
    pub fn version(&self) -> &'static str {
        SCHEMA_VERSION
    }
    /// Load a map from an mzXML file.
    ///
    /// Plain, gzip and bzip2 inputs are detected by content, independently of
    /// the suffix. The document path and the content-detected file type are
    /// recorded on the result, as source `MzXMLFile::load` does
    /// (`MzXMLFile.cpp:44-45`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when the file cannot be opened — the source
    /// `Exception::FileNotFound` — [`Error::Parse`] for malformed XML or peak
    /// data, [`Error::Unsupported`] for a feature this port does not represent,
    /// and [`Error::InvalidValue`] when a [`ReadLimits`] ceiling is exceeded.
    pub fn load(&self, path: impl AsRef<Path>) -> Result<MSExperiment> {
        load_with_options(path, &self.options)
    }
    /// [`MzXMLFile::load`], reporting progress to `logger` as the source's
    /// `load` does; see the free [`load_with_progress`].
    ///
    /// # Errors
    ///
    /// As [`load_with_progress`].
    pub fn load_with_progress(
        &self,
        path: impl AsRef<Path>,
        logger: &mut ProgressLogger,
    ) -> Result<MSExperiment> {
        load_with_progress(path, &self.options, logger)
    }
    /// Load into an existing map, replacing it only on success.
    ///
    /// # Errors
    ///
    /// As [`MzXMLFile::load`]. The destination is untouched on failure.
    pub fn load_into(&self, path: impl AsRef<Path>, destination: &mut MSExperiment) -> Result<()> {
        load_into_with_options(path, destination, &self.options)
    }
    /// Store a map in an mzXML file, publishing it atomically.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when the file cannot be created — the source
    /// `Exception::UnableToCreateFile` — [`Error::InvalidValue`] for a record
    /// this format cannot represent or a [`WriteLimits`] ceiling, and
    /// [`Error::UnsortedData`] when `force_mq_compatibility` needs sorted m/z.
    pub fn store(&self, path: impl AsRef<Path>, experiment: &MSExperiment) -> Result<()> {
        store_with_options(path, experiment, &self.write_options)
    }
    /// [`MzXMLFile::store`], reporting progress to `logger` as the source's
    /// `store` does; see the free [`store_with_progress`].
    ///
    /// # Errors
    ///
    /// As [`store_with_progress`].
    pub fn store_with_progress(
        &self,
        path: impl AsRef<Path>,
        experiment: &MSExperiment,
        logger: &mut ProgressLogger,
    ) -> Result<()> {
        store_with_progress(path, experiment, &self.write_options, logger)
    }
    /// Validate a file against the bundled `mzXML_idx_3.1.xsd`, the
    /// [`SCHEMA`] the source constructor registers, with the three schemas it
    /// includes.
    ///
    /// Source `MzXMLFile::isValid(filename, os)`, inherited from
    /// `Internal::XMLFile`: the messages the source writes to `os` are the
    /// report's diagnostics, and the source's `bool` is
    /// [`is_valid`](crate::format::xml_schema::SchemaValidationReport::is_valid).
    /// Available with the `xml-schema` feature, which brings in the libxml2
    /// validator; the source always has Xerces. The indexed schema is used for
    /// every document, indexed or not, as in the source; its `index` element
    /// is optional. An mzXML 2.1 document, whose namespace differs, is invalid
    /// against it, as it is in the source.
    ///
    /// # Errors
    ///
    /// As [`xml_schema::validate`](crate::format::xml_schema::validate): an
    /// I/O failure, where the source throws `Exception::FileNotFound`, and
    /// input that is not well-formed XML, where the source returns `false`.
    #[cfg(feature = "xml-schema")]
    pub fn is_valid(
        &self,
        path: impl AsRef<Path>,
    ) -> Result<crate::format::xml_schema::SchemaValidationReport> {
        crate::format::xml_schema::validate(crate::format::xml_schema::SchemaKind::MzXML, path)
    }
    /// Transform a file while loading, handing every scan to `consumer` and
    /// storing nothing, the source `MzXMLFile::transform`
    /// (`MzXMLFile.cpp:59-72`).
    ///
    /// `skip_full_count` suppresses the counting first pass; the consumer then
    /// receives an expected size of zero.
    ///
    /// # Errors
    ///
    /// As [`MzXMLFile::load`], plus any error the consumer returns.
    pub fn transform(
        &self,
        path: impl AsRef<Path>,
        consumer: &mut dyn MSDataConsumer,
        skip_full_count: bool,
    ) -> Result<TransformReport> {
        transform_with_options(
            path,
            consumer,
            &TransformOptions {
                skip_full_count,
                read: self.options.clone(),
            },
        )
    }
    /// Transform a file while loading and also store the result, the source
    /// `MzXMLFile::transform` overload that takes a map
    /// (`MzXMLFile.cpp:74-89`).
    ///
    /// The source forces `always_append_data` for this pass; so does this port,
    /// so `destination` receives every consumed spectrum whatever the option
    /// says.
    ///
    /// # Errors
    ///
    /// As [`MzXMLFile::transform`]. The destination is untouched on failure.
    pub fn transform_into(
        &self,
        path: impl AsRef<Path>,
        consumer: &mut dyn MSDataConsumer,
        destination: &mut MSExperiment,
        skip_full_count: bool,
    ) -> Result<TransformReport> {
        transform_into_with_options(
            path,
            consumer,
            destination,
            &TransformOptions {
                skip_full_count,
                read: self.options.clone(),
            },
        )
    }
}

// ---------------------------------------------------------------------------
// Free-function entry points
// ---------------------------------------------------------------------------

/// Read an mzXML document from a stream with source-default options.
///
/// # Errors
///
/// As [`MzXMLFile::load`], except that no document path is recorded.
pub fn read(input: impl BufRead) -> Result<MSExperiment> {
    read_with_options(input, &ReadOptions::default())
}

/// Read an mzXML document from a stream.
///
/// # Errors
///
/// As [`MzXMLFile::load`], except that no document path is recorded.
pub fn read_with_options(input: impl BufRead, options: &ReadOptions) -> Result<MSExperiment> {
    let mut report = ReadReport::default();
    read_with_report(input, options, &mut report)
}

/// Read an mzXML document and keep the non-fatal diagnostics the source logs.
///
/// # Errors
///
/// As [`MzXMLFile::load`]. `report` may already carry diagnostics when an error
/// is returned; they describe input seen before the failure.
pub fn read_with_report(
    input: impl BufRead,
    options: &ReadOptions,
    report: &mut ReadReport,
) -> Result<MSExperiment> {
    read_reporting(input, options, report, ProgressReporter::silent())
}

fn read_reporting(
    input: impl BufRead,
    options: &ReadOptions,
    report: &mut ReadReport,
    progress: ProgressReporter<'_>,
) -> Result<MSExperiment> {
    let detail = if options.peaks.metadata_only {
        Detail::MetadataOnly
    } else {
        Detail::Full
    };
    let mut run = Run::new(options, detail, report);
    run.progress = progress;
    run.parse(input, None)
}

/// Read only the experiment-wide metadata, stopping at the first `<scan>`.
///
/// This is the source `PeakFileOptions::metadata_only` path, which throws
/// `EndParsingSoftly` on the first scan (`MzXMLHandler.cpp:242-245`).
///
/// # Errors
///
/// As [`MzXMLFile::load`].
pub fn read_metadata(input: impl BufRead, options: &ReadOptions) -> Result<ExperimentalSettings> {
    let mut report = ReadReport::default();
    let mut run = Run::new(options, Detail::MetadataOnly, &mut report);
    Ok(run.parse(input, None)?.settings)
}

/// Count `<scan>` elements without storing spectra, the source
/// `XMLHandler::LD_RAWCOUNTS` load detail used by `transformFirstPass_`.
///
/// Returns the count and the header metadata the same pass collects.
///
/// # Errors
///
/// As [`MzXMLFile::load`].
pub fn read_scan_count(
    input: impl BufRead,
    options: &ReadOptions,
) -> Result<(usize, ExperimentalSettings)> {
    let mut report = ReadReport::default();
    let settings = {
        let mut run = Run::new(options, Detail::RawCounts, &mut report);
        run.parse(input, None)?.settings
    };
    Ok((report.scan_count, settings))
}

/// Load a map from a path with source-default options.
///
/// # Errors
///
/// As [`MzXMLFile::load`].
pub fn load(path: impl AsRef<Path>) -> Result<MSExperiment> {
    load_with_options(path, &ReadOptions::default())
}

/// Load a map from a path, recording its document identity.
///
/// # Errors
///
/// As [`MzXMLFile::load`].
pub fn load_with_options(path: impl AsRef<Path>, options: &ReadOptions) -> Result<MSExperiment> {
    load_reporting(path, options, ProgressReporter::silent())
}

/// Load a map from a path, reporting progress to `logger` as source
/// `MzXMLFile::load` does through its handler.
///
/// The calls are the handler's: `startProgress(0, scanCount, "loading mzXML
/// file")` at `<msRun>` (`MzXMLHandler.cpp:130-136`, with 0 when the attribute
/// is absent), `setProgress(n)` as the `n`-th `<scan>` begins, nested or not
/// and whether or not a filter drops it (`:281-282`), and `endProgress()` at
/// `</mzXML>`, once the spectra are complete (`:524-531`). A metadata-only load
/// stops at the first `<scan>` and makes no set there. The result is the one
/// [`load_with_options`] returns, and so is every error: both run the same
/// code, whose calls go nowhere for [`load_with_options`].
///
/// A failure after the start leaves the section open, as in the source, where
/// the exception bypasses `endProgress`: no `-- done` line is printed, the
/// nesting depth stays one level deeper, and a command backend of `logger`
/// refuses its next start.
///
/// # Errors
///
/// As [`MzXMLFile::load`], plus the errors of the progress calls
/// ([`ProgressLogger::start_progress`] and its siblings).
pub fn load_with_progress(
    path: impl AsRef<Path>,
    options: &ReadOptions,
    logger: &mut ProgressLogger,
) -> Result<MSExperiment> {
    load_reporting(path, options, ProgressReporter::new(Some(logger)))
}

fn load_reporting(
    path: impl AsRef<Path>,
    options: &ReadOptions,
    progress: ProgressReporter<'_>,
) -> Result<MSExperiment> {
    let path = path.as_ref();
    let mut document = crate::metadata::DocumentIdentifier::new();
    let text = path
        .to_str()
        .ok_or_else(|| Error::InvalidValue("mzXML filename is not UTF-8".into()))?;
    document.set_loaded_file_path(text)?;
    document.set_loaded_file_type(path)?;
    let mut report = ReadReport::default();
    let mut result = read_reporting(path_io::open(path)?, options, &mut report, progress)?;
    result.settings.document.loaded_file_path = document.loaded_file_path;
    result.settings.document.loaded_file_type = document.loaded_file_type;
    Ok(result)
}

/// Load into an existing map, replacing it only after a complete parse.
///
/// # Errors
///
/// As [`MzXMLFile::load`]. The destination is untouched on failure.
pub fn load_into(path: impl AsRef<Path>, destination: &mut MSExperiment) -> Result<()> {
    load_into_with_options(path, destination, &ReadOptions::default())
}

/// Load into an existing map with explicit options.
///
/// # Errors
///
/// As [`MzXMLFile::load`]. The destination is untouched on failure.
pub fn load_into_with_options(
    path: impl AsRef<Path>,
    destination: &mut MSExperiment,
    options: &ReadOptions,
) -> Result<()> {
    *destination = load_with_options(path, options)?;
    Ok(())
}

/// Write an mzXML document with lossless defaults.
///
/// # Errors
///
/// As [`MzXMLFile::store`].
pub fn write(output: &mut dyn Write, experiment: &MSExperiment) -> Result<()> {
    write_with_options(output, experiment, &WriteOptions::default())
}

/// Write an mzXML document.
///
/// Every record is checked before a byte is emitted, so a rejected experiment
/// produces no partial output on a fresh stream.
///
/// # Errors
///
/// As [`MzXMLFile::store`].
pub fn write_with_options(
    output: &mut dyn Write,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    write_engine(output, experiment, options, &mut ProgressReporter::silent())
}

/// Store a map at a path, publishing the file atomically. `.gz` and `.bz2`
/// suffixes select outer compression, as elsewhere in this crate.
///
/// # Errors
///
/// As [`MzXMLFile::store`].
pub fn store(path: impl AsRef<Path>, experiment: &MSExperiment) -> Result<()> {
    store_with_options(path, experiment, &WriteOptions::default())
}

/// Store a map at a path with explicit writer settings.
///
/// # Errors
///
/// As [`MzXMLFile::store`].
pub fn store_with_options(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    store_reporting(path, experiment, options, &mut ProgressReporter::silent())
}

/// Store a map at a path, reporting progress to `logger` as source
/// `MzXMLFile::store` does through its handler's `writeTo`.
///
/// The calls are the handler's: `startProgress(0, spectra, "storing mzXML
/// file")` before the document's first byte (`MzXMLHandler.cpp:636`),
/// `setProgress(s)` as spectrum `s` is reached, including an empty one that
/// MaxQuant compatibility skips (`:864`), and `endProgress()` after `</mzXML>`
/// (`:1119`). The destination is opened before the checks and the first
/// call, as the source's `XMLFile::save_` opens it before `writeTo`. The bytes
/// and every error are those of [`store_with_options`], which runs the same
/// code with the calls going nowhere; its checks precede the start.
///
/// A failure after the start leaves the section open, as described at
/// [`load_with_progress`].
///
/// # Errors
///
/// As [`MzXMLFile::store`], plus the errors of the progress calls.
pub fn store_with_progress(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
    logger: &mut ProgressLogger,
) -> Result<()> {
    store_reporting(
        path,
        experiment,
        options,
        &mut ProgressReporter::new(Some(logger)),
    )
}

fn store_reporting(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
    progress: &mut ProgressReporter<'_>,
) -> Result<()> {
    path_io::write(path.as_ref(), |writer| {
        write_engine(writer, experiment, options, progress)
    })
}

/// Transform a file while loading, storing nothing.
///
/// # Errors
///
/// As [`MzXMLFile::transform`].
pub fn transform_with_options(
    path: impl AsRef<Path>,
    consumer: &mut dyn MSDataConsumer,
    options: &TransformOptions,
) -> Result<TransformReport> {
    run_transform(path.as_ref(), consumer, None, options)
}

/// Transform a file while loading and put every consumed spectrum into
/// `destination`, replacing its previous contents only on success.
///
/// Source `MzXMLFile::transform` does not reset the map it is handed, so
/// repeated calls accumulate into it; this port replaces it, following the
/// `_into` convention the rest of this crate uses.
///
/// # Errors
///
/// As [`MzXMLFile::transform`]. The destination is untouched on failure.
pub fn transform_into_with_options(
    path: impl AsRef<Path>,
    consumer: &mut dyn MSDataConsumer,
    destination: &mut MSExperiment,
    options: &TransformOptions,
) -> Result<TransformReport> {
    run_transform(path.as_ref(), consumer, Some(destination), options)
}

fn run_transform(
    path: &Path,
    consumer: &mut dyn MSDataConsumer,
    destination: Option<&mut MSExperiment>,
    options: &TransformOptions,
) -> Result<TransformReport> {
    let mut report = TransformReport::default();
    // First pass: source transformFirstPass_ (MzXMLFile.cpp:91-109).
    {
        let detail = if options.skip_full_count {
            Detail::MetadataOnly
        } else {
            Detail::RawCounts
        };
        let mut first = ReadReport::default();
        let settings = {
            let mut run = Run::new(&options.read, detail, &mut first);
            run.parse(path_io::open(path)?, None)?.settings
        };
        report.expected_spectra = first.scan_count;
        consumer.set_expected_size(first.scan_count, 0)?;
        consumer.set_experimental_settings(&settings)?;
    }
    // Second pass: the source forces always_append_data for the map overload.
    let mut read_options = options.read.clone();
    if destination.is_some() {
        read_options.peaks.always_append_data = true;
    }
    let retain = destination.is_some() || read_options.peaks.always_append_data;
    let mut sink = Sink {
        consumer,
        retain,
        delivered: 0,
        stopped: false,
    };
    let parsed = {
        let mut run = Run::new(&read_options, Detail::Full, &mut report.read);
        run.parse(path_io::open(path)?, Some(&mut sink))?
    };
    report.delivered = sink.delivered;
    report.stopped = sink.stopped;
    if let Some(destination) = destination {
        *destination = parsed;
    }
    Ok(report)
}

struct Sink<'a> {
    consumer: &'a mut dyn MSDataConsumer,
    retain: bool,
    delivered: usize,
    stopped: bool,
}

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

struct PeaksBlock {
    slot: usize,
    precision: u32,
    compressed: bool,
    encoded: String,
}
struct Pending {
    spectrum: MSSpectrum,
    declared_peaks: usize,
    peaks_seen: bool,
}
struct Frame {
    /// Index into the pending pool, or `None` for a filtered or counted scan.
    slot: Option<usize>,
}

struct Run<'a> {
    /// The source handler's `logger_`.
    progress: ProgressReporter<'a>,
    options: &'a ReadOptions,
    detail: Detail,
    report: &'a mut ReadReport,
    experiment: MSExperiment,
    data_processing: Vec<Arc<DataProcessing>>,
    pending: Vec<Pending>,
    scans: Vec<Frame>,
    tags: Vec<String>,
    peaks: Option<PeaksBlock>,
    text: String,
    total_peaks: usize,
    metadata_entries: usize,
    stop: bool,
}

impl<'a> Run<'a> {
    fn new(options: &'a ReadOptions, detail: Detail, report: &'a mut ReadReport) -> Self {
        Self {
            progress: ProgressReporter::silent(),
            options,
            detail,
            report,
            experiment: MSExperiment::new(),
            data_processing: Vec::new(),
            pending: Vec::new(),
            scans: Vec::new(),
            tags: Vec::new(),
            peaks: None,
            text: String::new(),
            total_peaks: 0,
            metadata_entries: 0,
            stop: false,
        }
    }

    fn note(&mut self, message: impl Into<String>) {
        self.report.note(&self.options.limits, message);
    }

    fn parse(
        &mut self,
        input: impl BufRead,
        mut sink: Option<&mut Sink<'_>>,
    ) -> Result<MSExperiment> {
        let limits = self.options.limits;
        let cap = limits
            .max_xml_bytes
            .checked_add(1)
            .ok_or_else(|| Error::InvalidValue("mzXML byte limit must be below u64::MAX".into()))?;
        let mut reader = Reader::from_reader(input.take(cap));
        reader.config_mut().expand_empty_elements = true;
        reader.config_mut().enable_all_checks(true);
        let mut buffer = Vec::new();
        let mut seen_root = false;
        let mut seen_declaration = false;
        let mut ascii_only = false;
        let mut latin1 = false;
        loop {
            // The input is capped at the ceiling plus one byte, so an oversized
            // document is cut off mid-token and quick-xml reports the truncation
            // first. Report the ceiling in that case, so the cause is named.
            let event = match reader.read_event_into(&mut buffer) {
                Ok(event) => event,
                Err(error) => {
                    if reader.buffer_position() >= limits.max_xml_bytes {
                        return Err(budget("XML byte"));
                    }
                    return Err(invalid(error.to_string()));
                }
            };
            if ascii_only && !event.is_ascii() {
                return Err(Error::Unsupported(if latin1 {
                    "non-ASCII ISO-8859-1 mzXML requires transcoding".into()
                } else {
                    "non-ASCII bytes in US-ASCII mzXML".to_owned()
                }));
            }
            if reader.buffer_position() > limits.max_xml_bytes {
                return Err(budget("XML byte"));
            }
            match event {
                Event::Start(element) => {
                    let tag = local_name(&element)?;
                    if self.tags.is_empty() {
                        if seen_root {
                            return Err(invalid("more than one mzXML root element"));
                        }
                        seen_root = true;
                        if tag != "mzXML" {
                            return Err(invalid("root element is not mzXML"));
                        }
                    }
                    if self.tags.len() >= limits.max_depth {
                        return Err(budget("XML depth"));
                    }
                    self.tags.push(tag);
                    self.text.clear();
                    self.start(&element)?;
                    if self.stop {
                        break;
                    }
                }
                Event::End(_) => {
                    let tag = match self.tags.pop() {
                        Some(tag) => tag,
                        None => return Err(invalid("unbalanced mzXML end element")),
                    };
                    self.end(&tag, sink.as_deref_mut())?;
                    self.text.clear();
                    if self.stop {
                        break;
                    }
                }
                Event::Text(text) => {
                    let decoded = text.decode().map_err(|e| invalid(e.to_string()))?;
                    xml_string(&decoded)?;
                    self.characters(&decoded)?;
                }
                // A CDATA section carries no entity expansion, so it is accepted
                // as plain character content. Xerces splits long payloads into
                // several callbacks, which the upstream test calls "CDATA
                // splitting"; both spellings end up here.
                Event::CData(text) => {
                    let decoded = std::str::from_utf8(text.as_ref())
                        .map_err(|e| invalid(e.to_string()))?
                        .to_owned();
                    xml_string(&decoded)?;
                    self.characters(&decoded)?;
                }
                Event::Decl(declaration) => {
                    if seen_root || seen_declaration {
                        return Err(invalid("misplaced or duplicate XML declaration"));
                    }
                    seen_declaration = true;
                    if declaration
                        .version()
                        .map_err(|e| invalid(e.to_string()))?
                        .as_ref()
                        != b"1.0"
                    {
                        return Err(Error::Unsupported("only XML 1.0 is supported".into()));
                    }
                    if let Some(encoding) = declaration.encoding() {
                        let encoding = encoding.map_err(|e| invalid(e.to_string()))?;
                        latin1 = encoding.eq_ignore_ascii_case(b"ISO-8859-1");
                        ascii_only = latin1 || encoding.eq_ignore_ascii_case(b"US-ASCII");
                        if !ascii_only && !encoding.eq_ignore_ascii_case(b"UTF-8") {
                            return Err(Error::Unsupported(
                                "only UTF-8, US-ASCII or ISO-8859-1 mzXML XML is supported".into(),
                            ));
                        }
                    }
                }
                Event::Comment(text) => {
                    xml_string(&text.decode().map_err(|e| invalid(e.to_string()))?)?;
                }
                Event::PI(text) => {
                    xml_string(
                        std::str::from_utf8(text.as_ref()).map_err(|e| invalid(e.to_string()))?,
                    )?;
                }
                Event::DocType(_) => {
                    return Err(Error::Unsupported("XML DTDs are not supported".into()));
                }
                Event::GeneralRef(_) => {
                    return Err(Error::Unsupported(
                        "XML entity references are not supported".into(),
                    ));
                }
                // `expand_empty_elements` is set, so quick-xml turns every
                // self-closing element into a start/end pair. Refusing this
                // variant keeps a parser change from silently dropping the
                // `<peaks xsi:nil="true"/>` and `<nameValue/>` elements the
                // format is full of, rather than ignoring them.
                Event::Empty(_) => {
                    return Err(invalid("unexpanded empty XML element"));
                }
                Event::Eof => break,
            }
            buffer.clear();
        }
        if !seen_root {
            return Err(invalid("document has no mzXML root element"));
        }
        if self.peaks.is_some() {
            return Err(invalid("truncated mzXML peaks element"));
        }
        if !self.stop {
            self.flush(sink)?;
        } else {
            // MetadataOnly and a consumer stop both end the document early; a
            // partial scan tree is discarded rather than half-committed.
            self.pending.clear();
        }
        Ok(std::mem::take(&mut self.experiment))
    }

    /// Innermost open `<scan>` slot, or `None` inside a filtered scan.
    fn slot(&self) -> Option<usize> {
        self.scans.last().and_then(|frame| frame.slot)
    }
    /// Enclosing element while a start tag is on the stack, the source
    /// `*(open_tags_.end() - 2)`.
    fn parent_tag(&self) -> &str {
        self.tags
            .len()
            .checked_sub(2)
            .and_then(|index| self.tags.get(index))
            .map(String::as_str)
            .unwrap_or("")
    }
    /// Enclosing element in an end-tag handler, where the element's own name has
    /// already been popped. The source pops in `onEndElement` too, so its
    /// `onCharacters` sees the same parent through a still-pushed own tag.
    fn enclosing_tag(&self) -> &str {
        self.tags.last().map(String::as_str).unwrap_or("")
    }
    fn meta(&mut self, entries: usize) -> Result<()> {
        self.metadata_entries = self
            .metadata_entries
            .checked_add(entries)
            .filter(|count| *count <= self.options.limits.max_metadata_entries)
            .ok_or_else(|| budget("metadata entry"))?;
        Ok(())
    }

    fn start(&mut self, element: &BytesStart<'_>) -> Result<()> {
        let tag = self.tags.last().cloned().unwrap_or_default();
        // Source: every element inside a filtered scan is ignored except a
        // nested <scan>, which resets the filter (MzXMLHandler.cpp:126-129).
        if tag != "scan" && !self.scans.is_empty() && self.slot().is_none() {
            return Ok(());
        }
        match tag.as_str() {
            "msRun" => self.start_ms_run(element),
            "parentFile" => self.start_parent_file(element),
            "software" => self.start_software(element),
            "peaks" => self.start_peaks(element),
            "precursorMz" => self.start_precursor(element),
            "scan" => self.start_scan(element),
            "operator" => self.start_operator(element),
            "msManufacturer" => {
                let value = attribute(element, "value")?.unwrap_or_default();
                self.experiment.settings.instrument.vendor = value;
                Ok(())
            }
            "msModel" => {
                let value = attribute(element, "value")?.unwrap_or_default();
                self.experiment.settings.instrument.model = value;
                Ok(())
            }
            "msIonisation" => {
                let value = required(element, "value")?;
                let method = self.term::<IonizationMethod>(
                    IONIZATION_TERMS,
                    IonizationMethod::ALL,
                    &value,
                    "msIonization",
                );
                let instrument = &mut self.experiment.settings.instrument;
                if instrument.ion_sources.is_empty() {
                    instrument.ion_sources.push(IonSource::default());
                }
                instrument.ion_sources[0].ionization_method = method;
                Ok(())
            }
            "msMassAnalyzer" => {
                let value = required(element, "value")?;
                let kind = self.term::<AnalyzerType>(
                    ANALYZER_TERMS,
                    AnalyzerType::ALL,
                    &value,
                    "msMassAnalyzer",
                );
                let instrument = &mut self.experiment.settings.instrument;
                if instrument.mass_analyzers.is_empty() {
                    instrument.mass_analyzers.push(MassAnalyzer::default());
                }
                instrument.mass_analyzers[0].analyzer_type = kind;
                Ok(())
            }
            "msDetector" => {
                let value = required(element, "value")?;
                let kind = self.term::<DetectorType>(
                    DETECTOR_TERMS,
                    DetectorType::ALL,
                    &value,
                    "msDetector",
                );
                let instrument = &mut self.experiment.settings.instrument;
                if instrument.ion_detectors.is_empty() {
                    instrument.ion_detectors.push(IonDetector::default());
                }
                instrument.ion_detectors[0].detector_type = kind;
                Ok(())
            }
            "msResolution" => {
                let value = required(element, "value")?;
                let method = self.term::<ResolutionMethod>(
                    RESOLUTION_TERMS,
                    ResolutionMethod::ALL,
                    &value,
                    "msResolution",
                );
                let instrument = &mut self.experiment.settings.instrument;
                // Source indexes getMassAnalyzers()[0] unconditionally, which is
                // undefined behaviour when <msResolution> precedes
                // <msMassAnalyzer> (MzXMLHandler.cpp:436).
                if instrument.mass_analyzers.is_empty() {
                    instrument.mass_analyzers.push(MassAnalyzer::default());
                    self.note("msResolution appeared before msMassAnalyzer");
                }
                self.experiment.settings.instrument.mass_analyzers[0].resolution_method = method;
                Ok(())
            }
            "dataProcessing" => self.start_data_processing(element),
            "nameValue" => self.start_name_value(element),
            "processingOperation" => self.start_processing_operation(element),
            _ => Ok(()),
        }
    }

    fn term<T: PartialEq + Copy + Default>(
        &mut self,
        table: &str,
        all: &[T],
        value: &str,
        label: &str,
    ) -> T {
        match term_index(table, value).and_then(|index| all.get(index).copied()) {
            Some(found) => found,
            None => {
                self.note(format!("unexpected CV entry '{label}'='{value}'"));
                T::default()
            }
        }
    }

    fn start_ms_run(&mut self, element: &BytesStart<'_>) -> Result<()> {
        // scanCount is a reservation hint upstream; here it is only validated
        // against the scan ceiling so a hostile value allocates nothing.
        let mut declared = 0;
        if let Some(text) = attribute(element, "scanCount")? {
            match text.trim().parse::<i64>() {
                Ok(count) if count >= 0 => {
                    if usize::try_from(count).is_ok_and(|n| n > self.options.limits.max_scans) {
                        return Err(budget("declared scan count"));
                    }
                    declared = count;
                }
                _ => self.note(format!("invalid msRun scanCount '{text}'")),
            }
        }
        // `MzXMLHandler.cpp:135`: the declared count, 0 when absent.
        self.progress.start(0, declared, "loading mzXML file")?;
        self.data_processing.clear();
        self.report.scan_count = 0;
        Ok(())
    }

    fn start_parent_file(&mut self, element: &BytesStart<'_>) -> Result<()> {
        if self.experiment.settings.source_files.len() >= self.options.limits.max_source_files {
            return Err(budget("parentFile"));
        }
        let checksum = required(element, "fileSha1")?;
        // A SHA-1 checksum record must be 40 hexadecimal characters for
        // SourceFile::validate to accept it. Source stores any string with the
        // SHA1 tag; this port keeps the text and downgrades the algorithm so the
        // returned record stays valid.
        let checksum_type = if checksum.is_empty() {
            ChecksumType::Unknown
        } else if checksum.len() == 40 && checksum.bytes().all(|b| b.is_ascii_hexdigit()) {
            ChecksumType::Sha1
        } else {
            self.note(format!(
                "parentFile fileSha1 '{checksum}' is not a SHA-1 digest; checksum type recorded as unknown"
            ));
            ChecksumType::Unknown
        };
        self.experiment.settings.source_files.push(SourceFile {
            name: required(element, "fileName")?,
            file_type: required(element, "fileType")?,
            checksum,
            checksum_type,
            ..SourceFile::default()
        });
        Ok(())
    }

    fn start_software(&mut self, element: &BytesStart<'_>) -> Result<()> {
        match self.parent_tag() {
            "dataProcessing" => {
                let version = required(element, "version")?;
                let name = required(element, "name")?;
                let kind = required(element, "type")?;
                let time = attribute(element, "completionTime")?;
                // Source data_processing_.back() is undefined behaviour when
                // <software type="processing"> appears without an enclosing
                // <dataProcessing> (MzXMLHandler.cpp:153).
                let Some(processing) = self.data_processing.last().cloned() else {
                    self.note("software inside dataProcessing with no open dataProcessing");
                    return Ok(());
                };
                let mut processing = (*processing).clone();
                processing.software.version = version;
                processing.software.name = name;
                self.meta(1)?;
                processing
                    .metadata
                    .insert(PROCESSING_TYPE_KEY.into(), MetaValue::from(kind));
                processing.completion_time = self.completion_time(time.as_deref())?;
                if let Some(last) = self.data_processing.last_mut() {
                    *last = Arc::new(processing);
                }
                Ok(())
            }
            "msInstrument" => {
                let software = &mut self.experiment.settings.instrument.software;
                software.version = required(element, "version")?;
                software.name = required(element, "name")?;
                Ok(())
            }
            other => {
                self.note(format!("unexpected tag 'software' in tag '{other}'"));
                Ok(())
            }
        }
    }

    /// Source `asDateTime_` trims to 19 characters and logs a conversion error
    /// instead of failing the load (`XMLHandler.h:359-377`).
    fn completion_time(
        &mut self,
        text: Option<&str>,
    ) -> Result<Option<crate::data_structures::DateTime>> {
        let Some(text) = text else { return Ok(None) };
        let trimmed = text.trim();
        if trimmed.is_empty() {
            return Ok(None);
        }
        // Never byte-slice text from a file: take 19 characters, not 19 bytes.
        let clipped: String = trimmed.chars().take(19).collect();
        match crate::data_structures::DateTime::parse(&clipped) {
            Ok(value) => Ok(Some(value)),
            Err(_) => {
                self.note(format!("DateTime conversion error of \"{clipped}\""));
                Ok(None)
            }
        }
    }

    fn start_peaks(&mut self, element: &BytesStart<'_>) -> Result<()> {
        let Some(slot) = self.slot() else {
            // Source indexes spectrum_data_.back() unconditionally, which is
            // undefined behaviour for <peaks> outside a <scan>.
            return Err(invalid("peaks element outside a scan"));
        };
        if self.pending[slot].peaks_seen {
            return Err(invalid("more than one peaks element in one scan"));
        }
        let precision = attribute(element, "precision")?.unwrap_or_else(|| "32".into());
        let precision = match precision.as_str() {
            "32" => 32,
            "64" => 64,
            other => {
                self.note(format!("Invalid precision '{other}' in element 'peaks'"));
                return Err(Error::Unsupported(format!(
                    "mzXML peaks precision '{other}'; only 32 and 64 are defined"
                )));
            }
        };
        let byte_order = attribute(element, "byteOrder")?.unwrap_or_else(|| "network".into());
        if byte_order != "network" {
            self.note(format!(
                "Invalid or missing byte order '{byte_order}' in element 'peaks'. Must be 'network'!"
            ));
            return Err(Error::Unsupported(format!(
                "mzXML peaks byteOrder '{byte_order}'; only 'network' is defined"
            )));
        }
        let content = attribute(element, "contentType")?.unwrap_or_else(|| "m/z-int".into());
        if content != "m/z-int" {
            self.note(format!(
                "Invalid or missing pair order '{content}' in element 'peaks'. Must be 'm/z-int'!"
            ));
            return Err(Error::Unsupported(format!(
                "mzXML peaks contentType '{content}'; only 'm/z-int' is defined"
            )));
        }
        let compression = attribute(element, "compressionType")?.unwrap_or_else(|| "none".into());
        let compressed = match compression.as_str() {
            "none" => false,
            "zlib" => true,
            other => {
                self.note(format!(
                    "Invalid compression type {other}in elements 'peaks'. Must be 'none' or 'zlib'! "
                ));
                return Err(Error::Unsupported(format!(
                    "mzXML peaks compressionType '{other}'; only 'none' and 'zlib' are defined"
                )));
            }
        };
        if self.options.peaks.fill_data {
            self.peaks = Some(PeaksBlock {
                slot,
                precision,
                compressed,
                encoded: String::new(),
            });
        }
        self.pending[slot].peaks_seen = true;
        Ok(())
    }

    fn start_precursor(&mut self, element: &BytesStart<'_>) -> Result<()> {
        let Some(slot) = self.slot() else {
            return Err(invalid("precursorMz element outside a scan"));
        };
        if self.pending[slot].spectrum.precursors.len()
            >= self.options.limits.max_precursors_per_scan
        {
            return Err(budget("precursorMz"));
        }
        let mut precursor = Precursor::default();
        match attribute(element, "precursorIntensity")? {
            Some(text) => match text.trim().parse::<f64>() {
                Ok(value) if value.is_finite() => precursor.intensity = value as f32,
                _ => self.note(format!("invalid precursorIntensity '{text}'; using zero")),
            },
            None => self.note(
                "Mandatory attribute 'precursorIntensity' of tag 'precursorMz' not found! \
                 Setting precursor intensity to zero!",
            ),
        }
        if !precursor.intensity.is_finite() {
            precursor.intensity = 0.0;
        }
        if let Some(text) = attribute(element, "precursorCharge")? {
            match text.trim().parse::<i32>() {
                Ok(charge) => precursor.charge = charge,
                Err(_) => self.note(format!("invalid precursorCharge '{text}'")),
            }
        }
        // Source stores the full width in the lower offset and halves both
        // offsets when the m/z text arrives (MzXMLHandler.cpp:219-222, :574-579).
        if let Some(text) = attribute(element, "windowWideness")? {
            match text.trim().parse::<f64>() {
                Ok(width) if width.is_finite() && width >= 0.0 => {
                    precursor.isolation_window_lower_offset = width;
                }
                Ok(width) => {
                    self.note(format!("windowWideness '{width}' is negative or nonfinite"));
                }
                Err(_) => self.note(format!("invalid windowWideness '{text}'")),
            }
        }
        if let Some(text) = attribute(element, "activationMethod")? {
            if !text.is_empty() {
                match ActivationMethod::ALL
                    .iter()
                    .find(|method| method.short_name() == text)
                {
                    Some(method) => {
                        precursor.activation_methods.insert(*method);
                    }
                    None => self.note(format!("unknown activationMethod '{text}'")),
                }
            }
        }
        self.pending[slot].spectrum.precursors.push(precursor);
        Ok(())
    }

    fn start_scan(&mut self, element: &BytesStart<'_>) -> Result<()> {
        if self.detail == Detail::MetadataOnly {
            self.stop = true;
            self.tags.pop();
            return Ok(());
        }
        if self.report.scan_count >= self.options.limits.max_scans {
            return Err(budget("scan"));
        }
        self.report.scan_count += 1;
        let ms_level_text = required(element, "msLevel")?;
        let mut ms_level = match ms_level_text.trim().parse::<i64>() {
            Ok(level) if (0..=i64::from(u32::MAX)).contains(&level) => level as u32,
            _ => {
                return Err(invalid(format!("invalid scan msLevel '{ms_level_text}'")));
            }
        };
        if ms_level == 0 {
            self.note(
                "Invalid 'msLevel' attribute with value '0' in 'scan' element found. \
                 Assuming ms level 1!",
            );
            ms_level = 1;
        }
        let retention_time = match attribute(element, "retentionTime")? {
            Some(text) => {
                let (seconds, bad) = duration_seconds(&text);
                if bad {
                    self.note(format!(
                        "Double conversion error in retentionTime \"{text}\""
                    ));
                }
                finite(seconds, "scan retentionTime")?
            }
            None => 0.0,
        };
        // `MzXMLHandler.cpp:281-282`: the scans begun before this one, before
        // any filter decides whether it is kept.
        self.progress.set_count(self.report.scan_count - 1)?;
        let filtered = (self.options.peaks.has_rt_range()
            && !encloses(self.options.peaks.rt_range(), retention_time))
            || (self.options.peaks.has_ms_levels()
                && !i32::try_from(ms_level)
                    .is_ok_and(|level| self.options.peaks.contains_ms_level(level)))
            || self.detail == Detail::RawCounts;
        if filtered {
            self.scans.push(Frame { slot: None });
            return Ok(());
        }
        let declared_peaks = {
            let text = required(element, "peaksCount")?;
            match text.trim().parse::<i64>() {
                Ok(count) if count >= 0 => usize::try_from(count)
                    .ok()
                    .filter(|count| *count <= self.options.limits.max_peaks_per_scan)
                    .ok_or_else(|| budget("scan peaksCount"))?,
                _ => return Err(invalid(format!("invalid scan peaksCount '{text}'"))),
            }
        };
        self.total_peaks = self
            .total_peaks
            .checked_add(declared_peaks)
            .filter(|total| *total <= self.options.limits.max_total_peaks)
            .ok_or_else(|| budget("total peak"))?;
        let mut spectrum = MSSpectrum {
            ms_level,
            rt: retention_time,
            native_id: format!("scan={}", required(element, "num")?),
            data_processing: self.data_processing.clone(),
            ..MSSpectrum::default()
        };
        spectrum
            .peaks
            .try_reserve_exact(declared_peaks)
            .map_err(|_| budget("peak reservation"))?;
        // Source reads startMz/endMz into a default ScanWindow and stores it only
        // when either bound is nonzero (MzXMLHandler.cpp:308-314).
        let mut window = ScanWindow::default();
        if let Some(text) = attribute(element, "startMz")? {
            window.begin = self.number_attribute(&text, "startMz")?;
        }
        if let Some(text) = attribute(element, "endMz")? {
            window.end = self.number_attribute(&text, "endMz")?;
        }
        if window.begin != 0.0 || window.end != 0.0 {
            if window.begin > window.end {
                return Err(Error::InvalidRange(
                    "mzXML scan startMz exceeds endMz".into(),
                ));
            }
            spectrum.instrument_settings.scan_windows.push(window);
        }
        let polarity = attribute(element, "polarity")?.unwrap_or_else(|| "any".into());
        spectrum.instrument_settings.polarity =
            self.term::<Polarity>(POLARITY_TERMS, Polarity::ALL, &polarity, "polarity");
        if let Some(filter_line) = attribute(element, "filterLine")? {
            if !filter_line.is_empty() {
                self.meta(1)?;
                spectrum
                    .metadata
                    .insert(FILTER_STRING_KEY.into(), MetaValue::from(filter_line));
            }
        }
        match attribute(element, "scanType")?.unwrap_or_default().as_str() {
            // Unknown/unset scan type leaves the mode alone and warns about nothing.
            "" => {}
            "zoom" => {
                spectrum.instrument_settings.zoom_scan = true;
                spectrum.instrument_settings.scan_mode = ScanMode::MassSpectrum;
            }
            "Full" => {
                spectrum.instrument_settings.scan_mode = if ms_level > 1 {
                    ScanMode::MsnSpectrum
                } else {
                    ScanMode::MassSpectrum
                };
            }
            "SIM" => spectrum.instrument_settings.scan_mode = ScanMode::SelectedIonMonitoring,
            "SRM" | "MRM" => {
                spectrum.instrument_settings.scan_mode = ScanMode::SelectedReactionMonitoring;
            }
            "CRM" => {
                spectrum.instrument_settings.scan_mode = ScanMode::ConsecutiveReactionMonitoring;
            }
            // Q1, Q3 and the three non-standard ABI Sashimi types.
            "Q1" | "Q3" | "EMS" => {
                spectrum.instrument_settings.scan_mode = ScanMode::MassSpectrum;
            }
            "EPI" => {
                spectrum.instrument_settings.scan_mode = ScanMode::MassSpectrum;
                spectrum.ms_level = 2;
            }
            "ER" => {
                spectrum.instrument_settings.zoom_scan = true;
                spectrum.instrument_settings.scan_mode = ScanMode::MassSpectrum;
            }
            other => {
                spectrum.instrument_settings.scan_mode = ScanMode::MassSpectrum;
                self.note(format!("Unknown scan mode '{other}'. Assuming full scan"));
            }
        }
        let slot = self.pending.len();
        self.pending.push(Pending {
            spectrum,
            declared_peaks,
            peaks_seen: false,
        });
        self.scans.push(Frame { slot: Some(slot) });
        Ok(())
    }

    fn number_attribute(&mut self, text: &str, label: &str) -> Result<f64> {
        match text.trim().parse::<f64>() {
            Ok(value) if value.is_finite() => Ok(value),
            _ => {
                self.note(format!("Double conversion error of \"{text}\" in {label}"));
                Ok(0.0)
            }
        }
    }

    fn start_operator(&mut self, element: &BytesStart<'_>) -> Result<()> {
        let first = required(element, "first")?;
        let last = required(element, "last")?;
        let email = attribute(element, "email")?.unwrap_or_default();
        let phone = attribute(element, "phone")?.unwrap_or_default();
        let url = attribute(element, "URI")?.unwrap_or_default();
        let contacts = &mut self.experiment.settings.contacts;
        // Source resize(1) then back(): a second <operator> overwrites the first.
        if contacts.len() != 1 {
            contacts.clear();
            contacts.push(ContactPerson::default());
        }
        let contact = &mut contacts[0];
        contact.first_name = first;
        contact.last_name = last;
        contact.email = email;
        contact.url = url;
        if !phone.is_empty() {
            self.meta(1)?;
            self.experiment.settings.contacts[0]
                .metadata
                .insert(PHONE_KEY.into(), MetaValue::from(phone));
        }
        Ok(())
    }

    fn start_data_processing(&mut self, element: &BytesStart<'_>) -> Result<()> {
        if self.data_processing.len() >= self.options.limits.max_data_processing {
            return Err(budget("dataProcessing"));
        }
        let mut processing = DataProcessing::default();
        for (name, action) in [
            ("deisotoped", ProcessingAction::Deisotoping),
            ("chargeDeconvoluted", ProcessingAction::ChargeDeconvolution),
            ("centroided", ProcessingAction::PeakPicking),
        ] {
            if let Some(text) = attribute(element, name)? {
                if text == "true" || text == "1" {
                    processing.actions.insert(action);
                }
            }
        }
        if let Some(text) = attribute(element, "intensityCutoff")? {
            let cutoff = self.number_attribute(&text, "intensityCutoff")?;
            if cutoff != 0.0 {
                self.meta(1)?;
                processing
                    .metadata
                    .insert(INTENSITY_CUTOFF_KEY.into(), MetaValue::try_from(cutoff)?);
            }
        }
        self.data_processing.push(Arc::new(processing));
        Ok(())
    }

    fn start_name_value(&mut self, element: &BytesStart<'_>) -> Result<()> {
        let Some(name) = attribute(element, "name")?.filter(|name| !name.is_empty()) else {
            return Ok(());
        };
        let value = MetaValue::from(attribute(element, "value")?.unwrap_or_default());
        self.meta(1)?;
        match self.parent_tag() {
            "msInstrument" => {
                self.experiment
                    .settings
                    .instrument
                    .metadata
                    .insert(name, value);
            }
            "scan" => match self.slot() {
                Some(slot) => {
                    self.pending[slot].spectrum.metadata.insert(name, value);
                }
                None => return Err(invalid("nameValue in a scan with no spectrum")),
            },
            other => {
                self.note(format!("Unexpected tag 'nameValue' in tag '{other}'"));
            }
        }
        Ok(())
    }

    fn start_processing_operation(&mut self, element: &BytesStart<'_>) -> Result<()> {
        let Some(name) = attribute(element, "name")?.filter(|name| !name.is_empty()) else {
            return Ok(());
        };
        let value = MetaValue::from(attribute(element, "value")?.unwrap_or_default());
        let Some(processing) = self.data_processing.last().cloned() else {
            self.note("processingOperation with no open dataProcessing");
            return Ok(());
        };
        self.meta(1)?;
        let mut processing = (*processing).clone();
        processing.metadata.insert(name, value);
        if let Some(last) = self.data_processing.last_mut() {
            *last = Arc::new(processing);
        }
        Ok(())
    }

    fn characters(&mut self, text: &str) -> Result<()> {
        let Some(tag) = self.tags.last() else {
            if text.trim().is_empty() {
                return Ok(());
            }
            return Err(invalid("text outside the mzXML root element"));
        };
        if tag == "peaks" {
            // Source concatenates chunks and strips whitespace before decoding,
            // because "line breaks inside the base64 data are unfortunately no
            // exception" (MzXMLHandler.cpp:1158-1160).
            let limit = self.options.limits.max_encoded_bytes;
            if let Some(block) = &mut self.peaks {
                for byte in text.bytes() {
                    if byte.is_ascii_whitespace() {
                        continue;
                    }
                    if !byte.is_ascii() {
                        return Err(invalid("non-ASCII base64 text in peaks"));
                    }
                    if block.encoded.len() >= limit {
                        return Err(budget("base64 payload"));
                    }
                    block.encoded.push(char::from(byte));
                }
            }
            return Ok(());
        }
        if matches!(tag.as_str(), "offset" | "indexOffset" | "sha1") {
            return Ok(());
        }
        if self.text.len().saturating_add(text.len()) > self.options.limits.max_text_bytes {
            return Err(budget("text node"));
        }
        self.text.push_str(text);
        Ok(())
    }

    fn end(&mut self, tag: &str, sink: Option<&mut Sink<'_>>) -> Result<()> {
        match tag {
            "peaks" => {
                if let Some(block) = self.peaks.take() {
                    let peaks = self.decode(&block)?;
                    let spectrum = &mut self.pending[block.slot].spectrum;
                    spectrum.peaks = peaks;
                    if self.options.peaks.sort_spectra_by_mz && !spectrum.is_sorted() {
                        spectrum.sort_by_position()?;
                    }
                }
                Ok(())
            }
            "precursorMz" => self.end_precursor(),
            "comment" => self.end_comment(),
            "scan" => {
                if self.scans.pop().is_none() {
                    return Err(invalid("unbalanced mzXML scan element"));
                }
                // Source flushes only at nesting level zero, so a whole nested
                // scan tree is committed together (MzXMLHandler.cpp:539-542).
                if self.scans.is_empty()
                    && self.pending.len() >= self.options.peaks.max_data_pool_size.max(1)
                {
                    self.flush(sink)?;
                }
                Ok(())
            }
            "mzXML" => {
                self.flush(sink)?;
                // `MzXMLHandler.cpp:530`.
                self.progress.end()
            }
            _ => {
                if !self.text.trim().is_empty() && !matches!(tag, "offset" | "indexOffset" | "sha1")
                {
                    let text = self.text.clone();
                    self.note(format!(
                        "Unhandled character content '{text}' in element '{tag}'"
                    ));
                }
                Ok(())
            }
        }
    }

    fn end_precursor(&mut self) -> Result<()> {
        let Some(slot) = self.slot() else {
            return Ok(());
        };
        let text = self.text.trim().to_owned();
        if text.is_empty() {
            self.note("precursorMz element has no m/z text");
            return Ok(());
        }
        let mz = self.number_attribute(&text, "precursorMz")?;
        let Some(precursor) = self.pending[slot].spectrum.precursors.last_mut() else {
            return Ok(());
        };
        precursor.mz = mz;
        let width = precursor.isolation_window_lower_offset;
        if width != 0.0 {
            precursor.isolation_window_lower_offset = 0.5 * width;
            precursor.isolation_window_upper_offset = 0.5 * width;
        }
        if self.options.peaks.has_precursor_mz_range()
            && !encloses(self.options.peaks.precursor_mz_range(), mz)
        {
            // Source drops the spectrum already pushed and skips the rest of the
            // scan (MzXMLHandler.cpp:581-587). Only the innermost scan is
            // dropped here; upstream pops the last pushed entry, which for a
            // nested scan can be a different spectrum.
            self.pending.truncate(slot);
            if let Some(frame) = self.scans.last_mut() {
                frame.slot = None;
            }
        }
        Ok(())
    }

    fn end_comment(&mut self) -> Result<()> {
        let text = self.text.trim().to_owned();
        match self.enclosing_tag() {
            "msInstrument" => {
                self.meta(1)?;
                self.experiment
                    .settings
                    .instrument
                    .metadata
                    .insert(COMMENT_KEY.into(), MetaValue::from(text));
            }
            // Source ignores a <dataProcessing> comment.
            "dataProcessing" => {}
            "scan" => {
                if let Some(slot) = self.slot() {
                    self.meta(1)?;
                    self.pending[slot]
                        .spectrum
                        .metadata
                        .insert(COMMENT_KEY.into(), MetaValue::from(text));
                }
            }
            other => {
                if !text.is_empty() {
                    self.note(format!("Unhandled comment '{text}' in element '{other}'"));
                }
            }
        }
        Ok(())
    }

    fn decode(&mut self, block: &PeaksBlock) -> Result<Vec<Peak1D>> {
        let declared = self.pending[block.slot].declared_peaks;
        let width = if block.precision == 64 { 8 } else { 4 };
        let expected = declared
            .checked_mul(2 * width)
            .filter(|bytes| *bytes <= self.options.limits.max_decoded_bytes)
            .ok_or_else(|| budget("decoded array byte"))?;
        if block.encoded.is_empty() {
            // Source returns early for an empty payload, so peaksCount is not
            // compared against anything (MzXMLHandler.cpp:1153-1156).
            if declared != 0 {
                self.note(format!(
                    "peaksCount {declared} declared but the peaks element is empty"
                ));
            }
            return Ok(Vec::new());
        }
        // Charge the decoded ceiling against the symbol count *before* the
        // decode allocates. Whitespace was stripped as the payload
        // accumulated, so this is the exact decoded length of any well-formed
        // payload, padded or not: three bytes per four symbols, the trailing
        // partial group contributing one byte less than its symbol count, less
        // the padding. A malformed payload fails in `decode` below anyway.
        let symbols = block.encoded.as_bytes();
        let padding = symbols
            .iter()
            .rev()
            .take_while(|byte| **byte == b'=')
            .count()
            .min(2);
        let decoded_len = (symbols.len() / 4 * 3)
            .saturating_add((symbols.len() % 4).saturating_sub(1))
            .saturating_sub(padding);
        if decoded_len > self.options.limits.max_decoded_bytes {
            return Err(budget("decoded array byte"));
        }
        let raw = PEAKS_BASE64
            .decode(symbols)
            .map_err(|e| invalid(format!("invalid base64 in peaks: {e}")))?;
        if raw.len() > self.options.limits.max_decoded_bytes {
            return Err(budget("decoded array byte"));
        }
        let bytes = if block.compressed {
            inflate(&raw, expected, self.options.limits.max_decoded_bytes)?
        } else {
            raw
        };
        // Source asserts data.size() == 2 * peak_count_ and then iterates to
        // 2 * peak_count_ regardless; the assert is compiled out in release
        // builds, so a short payload reads past the decoded buffer
        // (MzXMLHandler.cpp:1175-1187, :1202-1213).
        if bytes.len() != expected {
            return Err(invalid(format!(
                "peaks payload decodes to {} bytes but peaksCount {declared} needs {expected}",
                bytes.len()
            )));
        }
        let mz_range = self
            .options
            .peaks
            .has_mz_range()
            .then(|| self.options.peaks.mz_range());
        let intensity_range = self
            .options
            .peaks
            .has_intensity_range()
            .then(|| self.options.peaks.intensity_range());
        let mut peaks = Vec::new();
        peaks
            .try_reserve_exact(declared)
            .map_err(|_| budget("peak reservation"))?;
        for pair in bytes.chunks_exact(2 * width) {
            let (first, second) = pair
                .split_at_checked(width)
                .ok_or_else(|| invalid("truncated mzXML peaks pair"))?;
            let (mz, intensity) = if width == 8 {
                (big_endian_f64(first)?, big_endian_f64(second)?)
            } else {
                (
                    f64::from(big_endian_f32(first)?),
                    f64::from(big_endian_f32(second)?),
                )
            };
            if !mz.is_finite() || !intensity.is_finite() {
                return Err(invalid("nonfinite value in mzXML peaks"));
            }
            // A 64-bit payload can hold an intensity beyond the f32 the kernel
            // stores. The source narrows it with an implicit conversion and
            // keeps the resulting infinity; this refuses instead.
            let stored = intensity as f32;
            if !stored.is_finite() {
                return Err(invalid("mzXML peak intensity overflows f32"));
            }
            if mz_range.is_none_or(|range| encloses(range, mz))
                && intensity_range.is_none_or(|range| encloses(range, intensity))
            {
                peaks.push(Peak1D {
                    mz,
                    intensity: stored,
                });
            }
        }
        Ok(peaks)
    }

    fn flush(&mut self, sink: Option<&mut Sink<'_>>) -> Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        let drained: Vec<MSSpectrum> = self
            .pending
            .drain(..)
            .map(|pending| pending.spectrum)
            .collect();
        match sink {
            Some(sink) => {
                for mut spectrum in drained {
                    let control = sink.consumer.consume_spectrum(&mut spectrum)?;
                    sink.delivered += 1;
                    if sink.retain {
                        self.experiment.spectra.push(spectrum);
                    }
                    if control == ControlFlow::Break(()) {
                        sink.stopped = true;
                        self.stop = true;
                        return Ok(());
                    }
                }
            }
            None => self.experiment.spectra.extend(drained),
        }
        Ok(())
    }
}

fn local_name(element: &BytesStart<'_>) -> Result<String> {
    std::str::from_utf8(element.local_name().as_ref())
        .map(str::to_owned)
        .map_err(|e| invalid(e.to_string()))
}
fn attribute(element: &BytesStart<'_>, name: &str) -> Result<Option<String>> {
    let mut found = None;
    for attribute in element.attributes().with_checks(false) {
        let attribute = attribute.map_err(|e| invalid(e.to_string()))?;
        let key = attribute.key;
        let matches =
            key.as_ref() == name.as_bytes() || key.local_name().as_ref() == name.as_bytes();
        if !matches {
            continue;
        }
        let value = attribute
            .unescape_value()
            .map_err(|e| invalid(e.to_string()))?
            .into_owned();
        xml_string(&value)?;
        if found.replace(value).is_some() {
            return Err(invalid(format!("duplicate attribute '{name}'")));
        }
    }
    Ok(found)
}
fn required(element: &BytesStart<'_>, name: &str) -> Result<String> {
    attribute(element, name)?.ok_or_else(|| {
        invalid(format!(
            "Required attribute '{name}' not present in mzXML element"
        ))
    })
}

/// `xs:duration` in seconds, following `MzXMLHandler.cpp:254-279`.
///
/// The source drops everything before the last `T`, so the date part of
/// `P1DT2H` is ignored, and it never inspects the leading sign, so `-PT1S`
/// reads as `+1`. This port removes a leading `-` and applies it to the total,
/// because the writer emits that spelling for a negative retention time and
/// losing it makes a store/load cycle change the data. That is the *only*
/// divergence: each component's text goes to `str::parse` exactly as the source
/// hands it to `XMLHandler::asDouble_`, with nothing stripped, so `PT-1S` reads
/// as `-1` in both and a component the source cannot convert — the `P1` of a
/// month-only `P1M`, say — contributes zero in both. Returns the seconds and
/// whether a component failed to parse; the source logs that and contributes
/// zero.
fn duration_seconds(text: &str) -> (f64, bool) {
    let trimmed = text.trim();
    let (negative, body) = match trimmed.strip_prefix('-') {
        Some(unsigned) => (true, unsigned),
        None => (false, trimmed),
    };
    let mut rest = body.rsplit('T').next().unwrap_or("");
    if !body.contains('T') {
        rest = body;
    }
    let mut seconds = 0.0;
    let mut bad = false;
    for (marker, scale) in [('H', 3600.0), ('M', 60.0), ('S', 1.0)] {
        if !rest.contains(marker) {
            continue;
        }
        let head = rest.split(marker).next().unwrap_or("");
        match head.trim().parse::<f64>() {
            Ok(value) if value.is_finite() => seconds += scale * value,
            _ => bad = true,
        }
        rest = rest.rsplit(marker).next().unwrap_or("");
    }
    (if negative { -seconds } else { seconds }, bad)
}

// mzXML declares byteOrder="network", so every stored value is big-endian
// (MzXMLHandler.cpp:177-182). A length that is not exactly the element width is
// an error rather than a panic, even though `chunks_exact` already guarantees it.
fn big_endian_f32(bytes: &[u8]) -> Result<f32> {
    bytes
        .try_into()
        .map(f32::from_be_bytes)
        .map_err(|_| invalid("truncated 32-bit mzXML peak value"))
}
fn big_endian_f64(bytes: &[u8]) -> Result<f64> {
    bytes
        .try_into()
        .map(f64::from_be_bytes)
        .map_err(|_| invalid("truncated 64-bit mzXML peak value"))
}

fn inflate(input: &[u8], expected: usize, ceiling: usize) -> Result<Vec<u8>> {
    let limit = expected.min(ceiling);
    let mut decoder = Decompress::new(true);
    let mut output = Vec::new();
    output
        .try_reserve_exact(limit)
        .map_err(|_| budget("decoded array byte"))?;
    let mut chunk = [0u8; 8192];
    loop {
        let before_in = decoder.total_in() as usize;
        let before_out = decoder.total_out();
        let capacity = limit
            .saturating_sub(output.len())
            .saturating_add(1)
            .min(chunk.len());
        let tail = input.get(before_in..).unwrap_or(&[]);
        let status = decoder
            .decompress(tail, &mut chunk[..capacity], FlushDecompress::None)
            .map_err(|e| invalid(format!("invalid zlib peaks array: {e}")))?;
        let written = usize::try_from(decoder.total_out() - before_out).unwrap_or(usize::MAX);
        if output.len().checked_add(written).is_none_or(|n| n > limit) {
            return Err(invalid("decompressed peaks array exceeds declared length"));
        }
        output.extend_from_slice(chunk.get(..written).unwrap_or(&[]));
        if status == Status::StreamEnd {
            if decoder.total_in() as usize != input.len() {
                return Err(invalid("trailing data after zlib peaks array"));
            }
            break;
        }
        if decoder.total_in() as usize == before_in && written == 0 {
            return Err(invalid("incomplete zlib peaks array"));
        }
    }
    Ok(output)
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

struct Counting<'a> {
    inner: &'a mut dyn Write,
    written: u64,
}
impl Write for Counting<'_> {
    fn write(&mut self, data: &[u8]) -> std::io::Result<usize> {
        let n = self.inner.write(data)?;
        self.written = self.written.saturating_add(n as u64);
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}
impl Counting<'_> {
    fn text(&mut self, value: &str) -> Result<()> {
        self.write_all(value.as_bytes())?;
        Ok(())
    }
    fn tabs(&mut self, count: u32) -> Result<()> {
        for _ in 0..count {
            self.write_all(b"\t")?;
        }
        Ok(())
    }
}

fn escape(value: &str) -> String {
    quick_xml::escape::escape(value)
        .replace('\n', "&#10;")
        .replace('\r', "&#13;")
        .replace('\t', "&#9;")
}

/// Preflight: everything the writer cannot represent is rejected here, before a
/// byte is emitted, so a failure leaves a fresh stream empty and the atomic
/// path writer leaves the destination file untouched.
fn preflight(experiment: &MSExperiment, options: &WriteOptions) -> Result<()> {
    let limits = options.limits;
    if experiment.spectra.len() > limits.max_spectra {
        return Err(budget("spectra written"));
    }
    let mut total = 0usize;
    for spectrum in &experiment.spectra {
        total = total
            .checked_add(spectrum.peaks.len())
            .filter(|count| *count <= limits.max_total_peaks)
            .ok_or_else(|| budget("peaks written"))?;
        if spectrum.ms_level == 0 || spectrum.ms_level > limits.max_ms_level {
            return Err(Error::InvalidValue(format!(
                "mzXML needs an MS level in 1..={}, found {}",
                limits.max_ms_level, spectrum.ms_level
            )));
        }
        if spectrum.precursors.len() > limits.max_precursors_per_spectrum {
            return Err(budget("precursors written"));
        }
        if spectrum.metadata.len() > limits.max_metadata_entries {
            return Err(budget("metadata entries written"));
        }
        finite(spectrum.rt, "spectrum retention time")?;
        xml_string(&spectrum.native_id)?;
        for peak in &spectrum.peaks {
            finite(peak.mz, "peak m/z")?;
            finite(f64::from(peak.intensity), "peak intensity")?;
        }
        for precursor in &spectrum.precursors {
            precursor.validate()?;
        }
        spectrum.instrument_settings.validate()?;
        for key in spectrum.metadata.keys() {
            xml_string(key)?;
        }
        if options.force_mq_compatibility && !spectrum.is_sorted() {
            // Source logs a non-fatal error and then writes begin()/rbegin()
            // m/z as lowMz/highMz anyway, producing wrong attributes
            // (MzXMLHandler.cpp:967, :973).
            return Err(Error::UnsortedData);
        }
    }
    let settings = &experiment.settings;
    if settings.source_files.len() > limits.max_source_files {
        return Err(budget("source files written"));
    }
    for source in &settings.source_files {
        source.validate()?;
        xml_string(&source.name)?;
        xml_string(&source.file_type)?;
    }
    for contact in &settings.contacts {
        xml_string(&contact.first_name)?;
        xml_string(&contact.last_name)?;
    }
    settings.instrument.software.validate()?;
    for analyzer in &settings.instrument.mass_analyzers {
        finite(analyzer.resolution, "mass analyzer resolution")?;
    }
    for processing in experiment
        .spectra
        .first()
        .map(|spectrum| spectrum.data_processing.as_slice())
        .unwrap_or(&[])
    {
        processing.validate()?;
    }
    Ok(())
}

/// Whether every string the writer will emit is ASCII.
///
/// The source writer hard-codes `encoding="ISO-8859-1"` (`MzXMLHandler.cpp:643`)
/// and streams `std::string` bytes unchanged, so a UTF-8 experiment produces a
/// document whose declaration contradicts its bytes — which this crate's own
/// readers then refuse. This port keeps the source declaration whenever it is
/// true and switches to `UTF-8` when it is not.
fn document_is_ascii(experiment: &MSExperiment, options: &WriteOptions) -> bool {
    fn meta_ascii(metadata: &MetaInfo, options: &WriteOptions) -> bool {
        metadata
            .iter()
            .all(|(key, value)| key.is_ascii() && options.meta_text(value).is_ascii())
    }
    let settings = &experiment.settings;
    let processing = experiment
        .spectra
        .first()
        .map(|spectrum| spectrum.data_processing.as_slice())
        .unwrap_or(&[]);
    settings.source_files.iter().all(|source| {
        source.name.is_ascii() && source.file_type.is_ascii() && source.checksum.is_ascii()
    }) && settings.contacts.iter().take(1).all(|contact| {
        contact.first_name.is_ascii()
            && contact.last_name.is_ascii()
            && contact.email.is_ascii()
            && contact.url.is_ascii()
            && meta_ascii(&contact.metadata, options)
    }) && settings.instrument.vendor.is_ascii()
        && settings.instrument.model.is_ascii()
        && settings.instrument.software.name.is_ascii()
        && settings.instrument.software.version.is_ascii()
        && meta_ascii(&settings.instrument.metadata, options)
        && experiment.spectra.iter().all(|spectrum| {
            spectrum.native_id.is_ascii() && meta_ascii(&spectrum.metadata, options)
        })
        && processing.iter().all(|step| {
            step.software.name.is_ascii()
                && step.software.version.is_ascii()
                && meta_ascii(&step.metadata, options)
        })
}

fn write_engine(
    output: &mut dyn Write,
    experiment: &MSExperiment,
    options: &WriteOptions,
    progress: &mut ProgressReporter<'_>,
) -> Result<()> {
    preflight(experiment, options)?;
    // `MzXMLHandler.cpp:636`.
    progress.start(
        0,
        progress_value(experiment.spectra.len())?,
        "storing mzXML file",
    )?;
    let mut out = Counting {
        inner: output,
        written: 0,
    };
    write_header(&mut out, experiment, options)?;
    write_instrument(&mut out, experiment, options)?;
    write_data_processing(&mut out, experiment, options)?;
    let index = write_scans(&mut out, experiment, options, progress)?;
    out.text("\t</msRun>\n")?;
    if options.write_index || options.force_mq_compatibility {
        let offset = out.written;
        out.text("<index name = \"scan\" >\n")?;
        for (id, position) in &index {
            out.text(&format!("<offset id = \"{id}\" >{position}</offset>\n"))?;
        }
        out.text("</index>\n")?;
        out.text(&format!("<indexOffset>{offset}</indexOffset>\n"))?;
    }
    out.text("</mzXML>\n")?;
    out.flush()?;
    // `MzXMLHandler.cpp:1119`.
    progress.end()
}

fn write_header(
    out: &mut Counting<'_>,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    // Source counts only nonempty spectra and bumps a zero to one, because an
    // empty mzXML is not schema-valid (MzXMLHandler.cpp:625-634).
    let mut count = experiment
        .spectra
        .iter()
        .filter(|spectrum| !spectrum.peaks.is_empty())
        .count();
    if count == 0 {
        count = 1;
    }
    // Source uses the FIRST and LAST spectrum's retention time, not the minimum
    // and maximum, so a descending file writes endTime before startTime
    // (MzXMLHandler.cpp:637-642).
    let first = experiment.spectra.first().map_or(0.0, |s| s.rt);
    let last = experiment.spectra.last().map_or(0.0, |s| s.rt);
    let encoding = if document_is_ascii(experiment, options) {
        "ISO-8859-1"
    } else {
        "UTF-8"
    };
    out.text(&format!(
        "<?xml version=\"1.0\" encoding=\"{encoding}\"?>\n\
         <mzXML xmlns=\"http://sashimi.sourceforge.net/schema_revision/mzXML_3.1\" \n\
        \x20xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \n\
        \x20xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/mzXML_3.1\
        \x20http://sashimi.sourceforge.net/schema_revision/mzXML_3.1/mzXML_idx_3.1.xsd\">\n"
    ))?;
    out.text(&format!(
        "\t<msRun scanCount=\"{count}\" startTime=\"{}\" endTime=\"{}\" >\n",
        duration_text(first, options),
        duration_text(last, options),
    ))?;
    if experiment.settings.source_files.is_empty() {
        out.text(
            "\t\t<parentFile fileName=\"\" fileType=\"processedData\" \
             fileSha1=\"0000000000000000000000000000000000000000\"/>\n",
        )?;
    } else {
        for source in &experiment.settings.source_files {
            // mzXML's fileType is an enumeration, so the source searches the
            // OpenMS free-text type for "raw" (MzXMLHandler.cpp:662-671).
            let kind = if source.file_type.to_ascii_lowercase().contains("raw") {
                "RAWData"
            } else {
                "processedData"
            };
            let checksum =
                if source.checksum.len() == 40 && source.checksum_type == ChecksumType::Sha1 {
                    source.checksum.clone()
                } else {
                    "0".repeat(40)
                };
            out.text(&format!(
                "\t\t<parentFile fileName=\"{}\" fileType=\"{kind}\" fileSha1=\"{}\"/>\n",
                escape(&source.name),
                escape(&checksum),
            ))?;
        }
    }
    Ok(())
}

fn write_instrument(
    out: &mut Counting<'_>,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    let settings = &experiment.settings;
    if settings.instrument == Instrument::default() && settings.contacts.is_empty() {
        return Ok(());
    }
    let instrument = &settings.instrument;
    // MaxQuant's parameter defaults need a Thermo manufacturer; the source infers
    // one from an Xcalibur acquisition software name (MzXMLHandler.cpp:694-700).
    let mut manufacturer = instrument.vendor.clone();
    if options.force_mq_compatibility
        || (manufacturer.is_empty()
            && instrument
                .software
                .name
                .to_ascii_lowercase()
                .contains("xcalibur"))
    {
        manufacturer = "Thermo Scientific".into();
    }
    out.text("\t\t<msInstrument>\n")?;
    out.text(&format!(
        "\t\t\t<msManufacturer category=\"msManufacturer\" value=\"{}\"/>\n",
        escape(&manufacturer)
    ))?;
    out.text(&format!(
        "\t\t\t<msModel category=\"msModel\" value=\"{}\"/>\n",
        escape(&instrument.model)
    ))?;
    let ionization = instrument
        .ion_sources
        .first()
        .map(|source| {
            term_name(
                IONIZATION_TERMS,
                enum_index(IonizationMethod::ALL, source.ionization_method),
            )
        })
        .unwrap_or("");
    out.text(&format!(
        "\t\t\t<msIonisation category=\"msIonisation\" value=\"{}\"/>\n",
        escape(ionization)
    ))?;
    let analyzer = instrument
        .mass_analyzers
        .first()
        .map(|analyzer| {
            term_name(
                ANALYZER_TERMS,
                enum_index(AnalyzerType::ALL, analyzer.analyzer_type),
            )
        })
        .unwrap_or("");
    out.text(&format!(
        "\t\t\t<msMassAnalyzer category=\"msMassAnalyzer\" value=\"{}\"/>\n",
        escape(analyzer)
    ))?;
    let detector = instrument
        .ion_detectors
        .first()
        .map(|detector| {
            term_name(
                DETECTOR_TERMS,
                enum_index(DetectorType::ALL, detector.detector_type),
            )
        })
        .unwrap_or("");
    out.text(&format!(
        "\t\t\t<msDetector category=\"msDetector\" value=\"{}\"/>\n",
        escape(detector)
    ))?;
    out.text(&format!(
        "\t\t\t<software type=\"acquisition\" name=\"{}\" version=\"{}\"/>\n",
        escape(&instrument.software.name),
        escape(&instrument.software.version)
    ))?;
    let resolution = instrument
        .mass_analyzers
        .first()
        .map(|analyzer| {
            term_name(
                RESOLUTION_TERMS,
                enum_index(ResolutionMethod::ALL, analyzer.resolution_method),
            )
        })
        .unwrap_or("");
    // Must not be empty, or MaxQuant crashes on loading (MzXMLHandler.cpp:733).
    if !resolution.is_empty() {
        out.text(&format!(
            "\t\t\t<msResolution category=\"msResolution\" value=\"{}\"/>\n",
            escape(resolution)
        ))?;
    }
    if let Some(contact) = settings.contacts.first() {
        out.text(&format!(
            "\t\t\t<operator first=\"{}\" last=\"{}\"",
            escape(&contact.first_name),
            escape(&contact.last_name)
        ))?;
        if !contact.email.is_empty() {
            out.text(&format!(" email=\"{}\"", escape(&contact.email)))?;
        }
        if !contact.url.is_empty() {
            out.text(&format!(" URI=\"{}\"", escape(&contact.url)))?;
        }
        if let Some(phone) = contact.metadata.get(PHONE_KEY) {
            out.text(&format!(" phone=\"{}\"", escape(&options.meta_text(phone))))?;
        }
        out.text("/>\n")?;
    }
    write_user_param(out, &instrument.metadata, 3, "nameValue", options, &[])?;
    if let Some(comment) = instrument.metadata.get(COMMENT_KEY) {
        out.text(&format!(
            "\t\t\t<comment>{}</comment>\n",
            escape(&options.meta_text(comment))
        ))?;
    }
    out.text("\t\t</msInstrument>\n")?;
    Ok(())
}

fn write_data_processing(
    out: &mut Counting<'_>,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    // Source assigns the first spectrum's processing to the whole file
    // (MzXMLHandler.cpp:774-819).
    let processing = experiment
        .spectra
        .first()
        .map(|spectrum| spectrum.data_processing.as_slice())
        .unwrap_or(&[]);
    if processing.is_empty() {
        out.text(
            "\t\t<dataProcessing>\n\t\t\t<software type=\"processing\" name=\"\" version=\"\"/>\n\
             \t\t</dataProcessing>\n",
        )?;
        return Ok(());
    }
    for step in processing {
        out.text(&format!(
            "\t\t<dataProcessing deisotoped=\"{}\" chargeDeconvoluted=\"{}\" centroided=\"{}\"",
            usize::from(step.actions.contains(&ProcessingAction::Deisotoping)),
            usize::from(
                step.actions
                    .contains(&ProcessingAction::ChargeDeconvolution)
            ),
            usize::from(step.actions.contains(&ProcessingAction::PeakPicking)),
        ))?;
        if let Some(cutoff) = step.metadata.get(INTENSITY_CUTOFF_KEY) {
            out.text(&format!(
                " intensityCutoff=\"{}\"",
                escape(&options.meta_text(cutoff))
            ))?;
        }
        out.text(">\n\t\t\t<software type=\"")?;
        match step.metadata.get(PROCESSING_TYPE_KEY) {
            Some(kind) => out.text(&escape(&options.meta_text(kind)))?,
            None => out.text("processing")?,
        }
        out.text(&format!(
            "\" name=\"{}\" version=\"{}",
            escape(&step.software.name),
            escape(&step.software.version)
        ))?;
        if let Some(time) = step.completion_time.filter(|time| !time.is_null()) {
            out.text(&format!(
                "\" completionTime=\"{}",
                escape(&time.get().replace(' ', "T"))
            ))?;
        }
        out.text("\"/>\n")?;
        write_user_param(out, &step.metadata, 3, "processingOperation", options, &[])?;
        out.text("\t\t</dataProcessing>\n")?;
    }
    Ok(())
}

/// Whether every native ID is a plain number or a number after a `scan=` prefix,
/// as the source decides before renumbering (`MzXMLHandler.cpp:821-855`).
///
/// The source additionally warns when at least one nonempty ID was unparsable,
/// because renumbering then loses it. That warning is not reproduced: this port
/// has no log stream and the writer returns no report. A caller that cares can
/// classify its own native IDs, which is what this function does.
struct NativeIds {
    all_numbers: bool,
    all_prefixed: bool,
}
fn classify_native_ids(experiment: &MSExperiment) -> NativeIds {
    let mut all_numbers = true;
    let mut all_prefixed = true;
    for spectrum in &experiment.spectra {
        let id = spectrum.native_id.as_str();
        let stripped = match id.strip_prefix("scan=") {
            Some(rest) => rest,
            None => {
                all_prefixed = false;
                id
            }
        };
        if stripped.parse::<i32>().is_err() {
            all_numbers = false;
            all_prefixed = false;
        }
    }
    NativeIds {
        all_numbers,
        all_prefixed,
    }
}

fn write_scans(
    out: &mut Counting<'_>,
    experiment: &MSExperiment,
    options: &WriteOptions,
    progress: &mut ProgressReporter<'_>,
) -> Result<Vec<(i64, u64)>> {
    let ids = classify_native_ids(experiment);
    let mut index = Vec::new();
    let mut open_scans = 0u32;
    let mut written = 0i64;
    let spectra = &experiment.spectra;
    for (position, spectrum) in spectra.iter().enumerate() {
        // `MzXMLHandler.cpp:864`, before the MaxQuant skip.
        progress.set_count(position)?;
        if spectrum.peaks.is_empty() && options.force_mq_compatibility {
            // MaxQuant's XML parser cannot deal with empty spectra.
            continue;
        }
        written += 1;
        let ms_level = spectrum.ms_level;
        open_scans += 1;
        let id = if ids.all_prefixed {
            spectrum
                .native_id
                .strip_prefix("scan=")
                .and_then(|rest| rest.parse::<i32>().ok())
                .map_or(written, i64::from)
        } else if ids.all_numbers {
            spectrum.native_id.parse::<i32>().map_or(written, i64::from)
        } else {
            written
        };
        out.tabs(ms_level + 1)?;
        index.push((id, out.written));
        out.text(&format!(
            "<scan num=\"{id}\" msLevel=\"{ms_level}\" peaksCount=\"{}\" polarity=\"{}\"",
            spectrum.peaks.len(),
            match spectrum.instrument_settings.polarity {
                Polarity::Positive => "+",
                Polarity::Negative => "-",
                Polarity::Unknown => "any",
            }
        ))?;
        let mut scan_type = match spectrum.instrument_settings.scan_mode {
            ScanMode::Unknown => String::new(),
            ScanMode::MassSpectrum | ScanMode::Ms1Spectrum | ScanMode::MsnSpectrum => {
                if spectrum.instrument_settings.zoom_scan {
                    "zoom".into()
                } else {
                    "Full".into()
                }
            }
            ScanMode::SelectedIonMonitoring => "SIM".into(),
            ScanMode::SelectedReactionMonitoring => "SRM".into(),
            ScanMode::ConsecutiveReactionMonitoring => "CRM".into(),
            // Source warns and falls back to Full for every other mode.
            _ => "Full".into(),
        };
        if scan_type.is_empty() && options.force_mq_compatibility {
            scan_type = "Full".into();
        }
        if !scan_type.is_empty() {
            out.text(&format!(" scanType=\"{scan_type}\""))?;
        }
        if let Some(filter) = spectrum.metadata.get(FILTER_STRING_KEY) {
            out.text(&format!(
                " filterLine=\"{}\"",
                escape(&options.meta_text(filter))
            ))?;
        }
        out.text(&format!(
            " retentionTime=\"{}\"",
            duration_text(spectrum.rt, options)
        ))?;
        if let Some(window) = spectrum.instrument_settings.scan_windows.first() {
            out.text(&format!(
                " startMz=\"{}\" endMz=\"{}\"",
                options.number(window.begin),
                options.number(window.end)
            ))?;
            // The mzXML format can store only one scan window for each scan.
        }
        write_scan_statistics(out, spectrum, options)?;
        if ms_level == 2 {
            if let Some(precursor) = spectrum.precursors.first() {
                if precursor.activation_energy != 0.0 {
                    out.text(&format!(
                        " collisionEnergy=\"{}\" ",
                        options.number(precursor.activation_energy)
                    ))?;
                }
            }
        }
        out.text(">\n")?;
        for precursor in &spectrum.precursors {
            write_precursor(out, spectrum, precursor, ms_level, options)?;
        }
        write_peaks(out, spectrum, ms_level, options)?;
        write_user_param(
            out,
            &spectrum.metadata,
            ms_level + 2,
            "nameValue",
            options,
            &ATTRIBUTE_METADATA_KEYS,
        )?;
        if let Some(comment) = spectrum.metadata.get(COMMENT_KEY) {
            out.tabs(ms_level + 2)?;
            out.text(&format!(
                "<comment>{}</comment>\n",
                escape(&options.meta_text(comment))
            ))?;
        }
        // Close as many scans as the next MS level allows, so an MS2 scan stays
        // nested inside its MS1 parent (MzXMLHandler.cpp:1082-1096).
        //
        // A next spectrum that MaxQuant mode will skip for being empty counts
        // as no next spectrum, which closes the open scans instead of leaving
        // one open for a child that is never written; upstream reads its MS
        // level regardless. This is a one-step lookahead in both: neither scans
        // forward to the next spectrum that will actually be written.
        let next_ms_level = spectra
            .get(position + 1)
            .filter(|next| !(next.peaks.is_empty() && options.force_mq_compatibility))
            .map_or(0, |next| next.ms_level);
        if next_ms_level <= ms_level {
            for step in 0..=(ms_level - next_ms_level) {
                if open_scans == 0 {
                    break;
                }
                out.tabs(ms_level - step + 1)?;
                out.text("</scan>\n")?;
                open_scans -= 1;
            }
        }
    }
    for _ in 0..open_scans {
        out.text("\t\t</scan>\n")?;
    }
    Ok(index)
}

fn write_scan_statistics(
    out: &mut Counting<'_>,
    spectrum: &MSSpectrum,
    options: &WriteOptions,
) -> Result<()> {
    let mq = options.force_mq_compatibility;
    for (key, attribute, forced) in [
        ("lowest observed m/z", "lowMz", mq),
        ("highest observed m/z", "highMz", mq),
        ("base peak m/z", "basePeakMz", true),
        ("base peak intensity", "basePeakIntensity", mq),
        ("total ion current", "totIonCurrent", mq),
    ] {
        if let Some(value) = spectrum.metadata.get(key) {
            out.text(&format!(
                " {attribute}=\"{}\"",
                escape(&options.meta_text(value))
            ))?;
            continue;
        }
        if !forced {
            continue;
        }
        let computed = match attribute {
            "lowMz" => spectrum.peaks.first().map_or(0.0, |peak| peak.mz),
            "highMz" => spectrum.peaks.last().map_or(0.0, |peak| peak.mz),
            "basePeakMz" => spectrum.base_peak().map_or(0.0, |peak| peak.mz),
            "basePeakIntensity" => spectrum
                .base_peak()
                .map_or(0.0, |peak| f64::from(peak.intensity)),
            _ => f64::from(spectrum.calculate_tic()),
        };
        out.text(&format!(" {attribute}=\"{}\"", options.number(computed)))?;
    }
    Ok(())
}

fn write_precursor(
    out: &mut Counting<'_>,
    spectrum: &MSSpectrum,
    precursor: &Precursor,
    ms_level: u32,
    options: &WriteOptions,
) -> Result<()> {
    out.tabs(ms_level + 2)?;
    let intensity = f64::from(precursor.intensity);
    let intensity = if options.integer_precursor_intensity {
        format!("{}", intensity.trunc() as i64)
    } else {
        options.number(intensity)
    };
    out.text(&format!("<precursorMz precursorIntensity=\"{intensity}\""))?;
    if precursor.charge != 0 {
        out.text(&format!(" precursorCharge=\"{}\"", precursor.charge))?;
    }
    let width = precursor.isolation_window_lower_offset + precursor.isolation_window_upper_offset;
    if width > 0.0 {
        out.text(&format!(" windowWideness=\"{}\"", options.number(width)))?;
    }
    match precursor.activation_methods.iter().next() {
        Some(method) => {
            out.text(&format!(" activationMethod=\"{}\" ", method.short_name()))?;
        }
        None if options.force_mq_compatibility => {
            // A missing activation makes old MaxQuant versions crash.
            out.text(&format!(
                " activationMethod=\"{}\" ",
                ActivationMethod::Cid.short_name()
            ))?;
        }
        None => {}
    }
    // The ReAdW converter's more accurate monoisotopic m/z wins when present.
    let mz = spectrum
        .acquisition_info
        .acquisitions
        .first()
        .and_then(|acquisition| {
            acquisition
                .metadata
                .get("[Thermo Trailer Extra]Monoisotopic M/Z:")
        })
        .and_then(|value| value.as_f64().ok())
        .filter(|value| value.is_finite())
        .unwrap_or(precursor.mz);
    out.text(&format!(">{}</precursorMz>\n", options.number(mz)))?;
    Ok(())
}

fn write_peaks(
    out: &mut Counting<'_>,
    spectrum: &MSSpectrum,
    ms_level: u32,
    options: &WriteOptions,
) -> Result<()> {
    let bits = options.precision.bits();
    let compression = if options.zlib_compression && !options.force_mq_compatibility {
        "zlib"
    } else {
        "none"
    };
    let mut payload = Vec::new();
    let width = options.precision.width();
    payload
        .try_reserve_exact(spectrum.peaks.len().saturating_mul(2 * width))
        .map_err(|_| budget("peaks written"))?;
    for peak in &spectrum.peaks {
        match options.precision {
            PeakPrecision::Float32 => {
                payload.extend_from_slice(&(peak.mz as f32).to_be_bytes());
                payload.extend_from_slice(&peak.intensity.to_be_bytes());
            }
            PeakPrecision::Float64 => {
                payload.extend_from_slice(&peak.mz.to_be_bytes());
                payload.extend_from_slice(&f64::from(peak.intensity).to_be_bytes());
            }
        }
    }
    let payload = if compression == "zlib" && !payload.is_empty() {
        let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(&payload)?;
        encoder.finish()?
    } else {
        payload
    };
    let compressed_len = if compression == "zlib" {
        payload.len()
    } else {
        0
    };
    out.tabs(ms_level + 2)?;
    // Some parsers require these line breaks (MaxQuant's mzXML reader fails
    // otherwise) while others, mostly TPP tools such as SpectraST, cannot deal
    // with them (MzXMLHandler.cpp:1039-1050).
    if options.force_mq_compatibility {
        out.text(&format!(
            "<peaks precision=\"{bits}\"\n byteOrder=\"network\"\n contentType=\"m/z-int\"\n \
             compressionType=\"{compression}\"\n compressedLen=\"{compressed_len}\" "
        ))?;
    } else {
        out.text(&format!(
            "<peaks precision=\"{bits}\" byteOrder=\"network\" contentType=\"m/z-int\" \
             compressionType=\"{compression}\" compressedLen=\"{compressed_len}\" "
        ))?;
    }
    if spectrum.peaks.is_empty() {
        out.text(" xsi:nil=\"true\" />\n")?;
    } else {
        out.text(&format!(">{}</peaks>\n", STANDARD.encode(&payload)))?;
    }
    Ok(())
}

fn write_user_param(
    out: &mut Counting<'_>,
    metadata: &MetaInfo,
    indent: u32,
    tag: &str,
    options: &WriteOptions,
    skip: &[&str],
) -> Result<()> {
    for (key, value) in metadata {
        // Internally used meta info starts with '#' (MzXMLHandler.cpp:1141).
        if key.starts_with('#') || skip.contains(&key.as_str()) {
            continue;
        }
        out.tabs(indent)?;
        out.text(&format!(
            "<{tag} name=\"{}\" value=\"{}\"/>\n",
            escape(key),
            escape(&options.meta_text(value))
        ))?;
    }
    Ok(())
}

/// `xs:duration` text for a retention time in seconds.
///
/// The source writes the sign before `PT` for a scan's `retentionTime`
/// (`MzXMLHandler.cpp:948-953`) but not for `msRun`'s `startTime`/`endTime`
/// (`:648`), which yields the invalid `PT-1S` there; this port signs both.
fn duration_text(seconds: f64, options: &WriteOptions) -> String {
    let sign = if seconds < 0.0 { "-" } else { "" };
    format!("{sign}PT{}S", options.number(seconds.abs()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn general_format_matches_ostream_defaults() {
        assert_eq!(general_format(60.0, 6), "60");
        assert_eq!(general_format(3661.0, 6), "3661");
        assert_eq!(general_format(97.3951, 6), "97.3951");
        assert_eq!(general_format(1999.93, 6), "1999.93");
        assert_eq!(general_format(100.083336, 6), "100.083");
        assert_eq!(general_format(0.0, 6), "0");
        assert_eq!(general_format(1234567.0, 6), "1.23457e+06");
        assert_eq!(general_format(0.0000123, 6), "1.23e-05");
        assert_eq!(general_format(9.999999, 6), "10");
    }

    #[test]
    fn durations_follow_the_source_algorithm() {
        assert_eq!(duration_seconds("PT60S").0, 60.0);
        assert_eq!(duration_seconds("PT2M1S").0, 121.0);
        assert_eq!(duration_seconds("PT1H61S").0, 3661.0);
        // The source drops the day component with everything before 'T'.
        assert_eq!(duration_seconds("P1DT2H").0, 7200.0);
        // No 'T' and no H/M/S contributes nothing, as upstream.
        assert_eq!(duration_seconds("60").0, 0.0);
        assert!(duration_seconds("PTxS").1);
        // Native: the sign is honoured; upstream reads +1 here.
        assert_eq!(duration_seconds("-PT1S").0, -1.0);
    }

    #[test]
    fn cv_tables_line_up_with_the_metadata_enums() {
        assert_eq!(term_index(POLARITY_TERMS, "+"), Some(1));
        assert_eq!(
            term_index(IONIZATION_TERMS, "MALDI"),
            Some(enum_index(IonizationMethod::ALL, IonizationMethod::Maldi))
        );
        assert_eq!(
            term_index(ANALYZER_TERMS, "Quadrupole Ion Trap"),
            Some(enum_index(AnalyzerType::ALL, AnalyzerType::PaulIonTrap))
        );
        assert_eq!(
            term_index(DETECTOR_TERMS, "Faraday Cup"),
            Some(enum_index(DetectorType::ALL, DetectorType::FaradayCup))
        );
        assert_eq!(
            term_index(RESOLUTION_TERMS, "FWHM"),
            Some(enum_index(ResolutionMethod::ALL, ResolutionMethod::Fwhm))
        );
        // An unknown term resolves to the enum's leading empty slot.
        assert_eq!(term_index(IONIZATION_TERMS, ""), Some(0));
    }

    #[test]
    fn drange_encloses_is_open_above() {
        let range = NumericRange {
            min: 115.0,
            max: 135.0,
        };
        assert!(encloses(range, 115.0));
        assert!(encloses(range, 130.0));
        assert!(!encloses(range, 135.0));
        assert!(!encloses(range, 114.9));
    }
}

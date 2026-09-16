// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded mzML 1.1 subset I/O. See `docs/MZML_SUPPORT.md` for metadata limitations.
//! Numeric peak arrays are little-endian f32/f64, optionally zlib compressed.
//! This is an event parser, but the returned experiment is held in memory.

pub use super::indexed_mzml::has_index;
#[cfg(feature = "mzml-schema")]
pub use super::mzml_schema::{
    SchemaDiagnostic, SchemaDiagnosticLevel, SchemaKind, SchemaValidationLimits,
    SchemaValidationOptions, SchemaValidationReport, validate_schema, validate_schema_reader,
    validate_schema_with_options,
};
#[path = "mzml_write_options.rs"]
mod peak_writer;
pub use peak_writer::{
    PeakWriteLimits, PeakWriteReport, store_with_peak_options, store_with_peak_options_and_limits,
    write_with_peak_options, write_with_peak_options_and_limits,
};
#[path = "mzml_counts.rs"]
mod counts;
pub use counts::{
    MAX_COUNT_EVENT_BYTES, MAX_COUNT_WORK, MAX_COUNT_XML_DEPTH, MzMLCounts, load_size,
    load_size_with_options, read_size, read_size_with_options,
};
#[path = "mzml_scaling.rs"]
mod scaling;
pub use scaling::{Allowance, InputScaling};
#[path = "mzml_centroid.rs"]
mod centroid;
pub use centroid::{CentroidInfoLimits, SpecInfo, centroid_info, centroid_info_with_options};
#[path = "mzml_consumer.rs"]
mod consumer;
#[path = "mzml_header.rs"]
mod header;
pub use crate::interfaces::MSDataConsumer;
pub use consumer::{
    TransformOptions, TransformReport, transform, transform_from, transform_from_into,
    transform_into, transform_into_with_options, transform_with_options,
};
#[path = "mzml_load.rs"]
mod load;
#[cfg(feature = "mzml-validation")]
pub use super::mzml_validator::{validate_semantics, validate_semantics_with_options};
use crate::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, Precursor,
    SpectrumType,
};
use crate::metadata::{DriftTimeUnit, MetaValue, MetaValueData, Product, Unit};
use crate::{Error, Result};
pub use load::LoadOptions;
#[path = "mzml_acquisition.rs"]
mod acquisition_metadata;
pub use acquisition_metadata::AcquisitionMode;
#[path = "mzml_numpress.rs"]
mod numpress_transport;
use super::numpress_coder::{self as coder, NumpressCompression, NumpressConfig};
pub use numpress_transport::{NumpressWriteOptions, NumpressWriteReport, write_with_numpress};

#[path = "mzml_paths.rs"]
mod paths;
pub use paths::{
    load, load_into, load_into_with_options, load_with_options, store, store_with_options,
};
#[path = "mzml_precursor.rs"]
mod precursor_metadata;
#[path = "mzml_record.rs"]
mod record_transport;
#[path = "mzml_settings.rs"]
mod settings_metadata;
use base64::{Engine, engine::general_purpose::STANDARD};
use flate2::{Compression, Decompress, FlushDecompress, Status, write::ZlibEncoder};
use quick_xml::{
    NsReader,
    events::{BytesStart, Event},
    name::ResolveResult,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::{BufRead, Write},
};

const NS: &[u8] = b"http://psi.hupo.org/ms/mzml";
const NAME_KEY: &str = "openms-rust:name";
/// FAIMS volt unit as spelled by the earlier OpenMS writer that produced the
/// upstream `FAIMS_CV-60C_V-45_Interleaved.mzML:324` and
/// `FAIMS_test_data.mzML:171`; the pinned writer emits `UO:0000218`
/// (`MzMLHandler.cpp:5415`). The source reader ignores unit attributes on
/// mobility terms, so it reads both spellings.
const LEGACY_FAIMS_VOLT: &str = "UO:000218";

/// Resource limits, acquisition normalization and source-compatibility switches
/// for mzML loading.
///
/// Limits apply before allocation from declared lengths and during decoding.
/// Every cumulative quantity is bounded twice: by the absolute `max_*` ceiling
/// here and by the size-derived allowance in [`ReadOptions::scaling`], which
/// grows with the XML bytes consumed so far. A charge fails when it exceeds
/// either. The absolute ceilings default to "unbounded" (`usize::MAX`, or 1 TiB
/// of XML), so by default the size-derived allowances decide, and a document of
/// any realistic size reads while a small one that declares huge arrays or
/// amplifies parameter-group references is refused. An explicit absolute
/// ceiling is always honoured as given. See
/// `docs/MZML_READER_SCALE_SUPPORT.md`.
///
/// Counts that each need their own start tag (`max_records`,
/// `max_total_arrays`, `max_param_groups`) are bounded by the XML size too, but
/// a record is the most memory-amplifying thing a document can declare per XML
/// byte, so their allowances grow once per group of consumed bytes
/// ([`InputScaling::records`], [`InputScaling::arrays`],
/// [`InputScaling::param_groups`]) rather than per byte.
#[derive(Clone, Copy, Debug)]
pub struct ReadOptions {
    /// Materialize source dummy scans, or preserve the canonical empty native form.
    pub acquisition_mode: AcquisitionMode,
    /// Maximum XML input bytes, including metadata and encoded arrays. The
    /// default, 1 TiB, only stops an unbounded stream; a file ends on its own.
    pub max_xml_bytes: u64,
    /// Maximum compressed or decoded bytes for each binary array. Default
    /// 64 MiB, 8 million `f64` values. The largest single array in any
    /// benchmark input declares 53,824 elements (the `UK222.mzML` TIC
    /// chromatogram, 431 KB decoded), so the default is 155 times the largest
    /// measured array and far beyond any instrument's single spectrum. This
    /// ceiling is per array and does not scale with the input: it also bounds
    /// the one transient value vector a decoded array allocates, which is what
    /// stops a small document whose few arrays each inflate to hundreds of
    /// megabytes.
    pub max_array_bytes: usize,
    /// Maximum raw declared spectrum/chromatogram points, including filtered
    /// records. Default unbounded; see [`InputScaling::peaks`].
    pub max_total_peaks: usize,
    /// Maximum decoded bytes across all primary and auxiliary arrays. Default
    /// unbounded; see [`InputScaling::array_bytes`].
    pub max_total_array_bytes: usize,
    /// Maximum elements across all arrays, including empty string elements.
    /// Default unbounded; see [`InputScaling::array_elements`].
    pub max_total_array_elements: usize,
    /// Maximum binary arrays, including empty placeholders. Default unbounded;
    /// see [`InputScaling::arrays`].
    pub max_total_arrays: usize,
    /// Maximum total spectra plus chromatograms. Default unbounded; see
    /// [`InputScaling::records`].
    pub max_records: usize,
    /// Maximum referenceable parameter groups (including unused and empty
    /// groups). Default unbounded; see [`InputScaling::param_groups`].
    pub max_param_groups: usize,
    /// Maximum groups, parameters/ref uses and supported acquisition
    /// descriptors combined. Default unbounded; see [`InputScaling::params`].
    pub max_total_params: usize,
    /// Cumulative parameter, acquisition descriptor and resolved-reference
    /// storage bytes. Default unbounded; see [`InputScaling::param_bytes`].
    pub max_param_bytes: usize,
    /// Size-derived allowances applied together with the absolute ceilings
    /// above. [`InputScaling::fixed`] restores the former fixed ceilings.
    pub scaling: InputScaling,
    /// Read an unparseable `run/@startTimeStamp` or processing completion time
    /// the way source `MzMLHandler` does.
    ///
    /// ProteoWizard writes `startTimeStamp="-infinity"` when the vendor file
    /// has no acquisition date (PXD001819 `UPS1_50amol_R1.mzML`), and Boost
    /// spells the other special values `infinity` and `not-a-date-time`.
    /// `DateTime::set` rejects all of them, as source `DateTime::set` throws
    /// `Exception::ParseError`.
    ///
    /// `true`, the default, is what source `MzMLHandler` does, which has no
    /// strict mode here. For `startTimeStamp`, `XMLHandler::asDateTime_`
    /// (`XMLHandler.h:359-377`) catches the error, logs `DateTime conversion
    /// error` as a non-fatal error and leaves the run date-time unset, so the
    /// writer omits the attribute. For a `processingMethod` completion time
    /// (`MS:1000747`), `XMLHandler::cvParamToValue` (`XMLHandler.cpp:232-243`)
    /// warns and drops the term, leaving the completion time unset. Both were
    /// executed on the C++ Release build, where `FileInfo` and `FileConverter`
    /// read such a file and exit 0. Each dropped value writes one line to the
    /// crate's warning log stream.
    ///
    /// `false` rejects such a document with [`Error::InvalidValue`]
    /// (`invalid DateTime input or calendar fields`) instead, for a caller that
    /// would rather not lose the value silently. The default is lenient because
    /// the alternative is refusing files that every OpenMS tool reads: it made
    /// all 60 Rust tool runs of the OpenMS4 smoke benchmark exit 6.
    pub source_invalid_timestamps: bool,
    /// Read a dangling header reference the way source `MzMLHandler` does.
    ///
    /// Covers a `softwareRef` on an `instrumentConfiguration` or
    /// `processingMethod`, and a `dataProcessingRef` or
    /// `defaultDataProcessingRef` on a record list, record or binary array,
    /// whose ID names no preceding definition. `false`, the default, rejects
    /// such a document with [`Error::Parse`] (`unresolved softwareRef`,
    /// `unresolved dataProcessingRef`), because the reference cannot be kept.
    ///
    /// `true` selects the source behaviour, which `std::map::operator[]`
    /// produces there (`MzMLHandler.cpp:920-952`, `:1034`, `:1264`, `:1288`):
    /// the software becomes `Software::default()` and the processing history
    /// becomes empty, so the reference is dropped. Each distinct dangling ID is
    /// reported once per read on the crate's warning log stream, where the
    /// source is silent. Malformed IDs, and unresolved `sourceFileRef`,
    /// `sampleRef`, instrument configuration and parameter group references,
    /// remain errors. Tool paths that reproduce source loading enable it; see
    /// `docs/MZML_HEADER_SUPPORT.md`.
    pub source_dangling_references: bool,
    /// Round a unit-converted **32-bit** `time array` back to `f32`, as source
    /// `MzMLHandlerHelper::decodeBase64Arrays` does.
    ///
    /// Applies only to an `MS:1000595` time array that is both 32-bit
    /// (`MS:1000521`) and carries `unitAccession="UO:0000031"` (minute), which
    /// is how ProteoWizard writes the TIC chromatogram of a Thermo run. The
    /// source decodes such an array into a `std::vector<float>` and then
    /// applies the minute multiplier in place with `it = it * unit_multiplier`
    /// (`MzMLHandlerHelper.cpp:217-222`), where `it` binds to `float&`: the
    /// product is computed in `double` and narrowed back to `float` on
    /// assignment, so the seconds value is kept at 32-bit precision. The
    /// 64-bit branch above it (`:210-216`) keeps full precision, and a
    /// Numpress array is forced to 64-bit before the multiplier runs
    /// (`:183-191`), so the loss is specific to this one combination.
    ///
    /// `false`, the default, computes `f64::from(value) * 60.0` and keeps the
    /// `f64` result, because narrowing back discards about seven decimal
    /// digits that the port has already recovered and nothing in the format
    /// asks for. `true` selects the source behaviour, and is what a tool path
    /// that reproduces source loading passes. Measured on the 40,856-point TIC
    /// chromatogram of `profile_hr_qe_silac_uk222/UK222.mzML`, the two differ
    /// on 38,107 of 40,856 times by at most 2.44e-4 s, which moves 21.6% of
    /// the point spacings by more than 1e-3 relative and at most 6.2e-3
    /// relative. That is enough to move the chromatogram `PeakPickerHiRes`
    /// picks from that file by up to 3.19e-3 s in apex position and 1.75e-3
    /// relative in intensity, which is what this switch exists for; see
    /// `docs/MZML_SUPPORT.md`, section "The minute conversion of a 32-bit time
    /// array", for the measurement and the executed evidence.
    ///
    /// The narrowing can overflow where the source silently stores an
    /// infinity: a finite 32-bit time above 5.67e36 minutes has no finite
    /// 32-bit product with 60. This port rejects that document with
    /// [`Error::Parse`] (`nonfinite binary value`) instead, as it already does
    /// for a decoded non-finite value.
    pub source_time_array_precision: bool,
    /// Accept NaN and infinite values in auxiliary float arrays (named data
    /// arrays other than m/z, intensity and time), as source `MzMLHandler`
    /// does, which decodes them without a check.
    ///
    /// The port rejects a non-finite decoded value with [`Error::Parse`]
    /// (`nonfinite binary value`) by default. Coordinates and intensities stay
    /// rejected either way, because the kernel requires them finite. The source
    /// writes such auxiliary values itself: the `debug/input.mzML` of
    /// `FeatureFinderAlgorithmPicked`'s debug mode holds NaN trace scores when
    /// `mass_trace:min_spectra` is 1. [`write_source_float_arrays`] is the
    /// matching writer. Not part of [`ReadOptions::source`], so the tools that
    /// use that preset keep refusing such values at load time.
    pub source_nonfinite_float_arrays: bool,
}
impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            acquisition_mode: AcquisitionMode::default(),
            max_xml_bytes: 1 << 40,
            max_array_bytes: 1 << 26,
            max_total_peaks: usize::MAX,
            max_records: usize::MAX,
            max_total_array_bytes: usize::MAX,
            max_total_array_elements: usize::MAX,
            max_total_arrays: usize::MAX,
            max_param_groups: usize::MAX,
            max_total_params: usize::MAX,
            max_param_bytes: usize::MAX,
            scaling: InputScaling::default(),
            source_invalid_timestamps: true,
            source_dangling_references: false,
            source_time_array_precision: false,
            source_nonfinite_float_arrays: false,
        }
    }
}
impl ReadOptions {
    /// The default limits with every source-compatibility switch enabled.
    ///
    /// Sets [`ReadOptions::source_dangling_references`] and
    /// [`ReadOptions::source_time_array_precision`] on top of the default
    /// [`ReadOptions::source_invalid_timestamps`], so a document reads as
    /// source `MzMLHandler` reads it wherever this port otherwise refuses a
    /// loss or keeps precision the source drops. Limits and acquisition
    /// normalization stay at their defaults. This is what a TOPP tool that
    /// reproduces source loading passes.
    pub fn source() -> Self {
        Self {
            source_dangling_references: true,
            source_invalid_timestamps: true,
            source_time_array_precision: true,
            ..Self::default()
        }
    }
}

/// Writer stores f64 coordinates and f32 intensities; compression changes only encoding.
#[derive(Clone, Copy, Debug, Default)]
pub struct WriteOptions {
    /// Compress binary arrays with zlib; false writes uncompressed arrays.
    pub zlib_compression: bool,
}

fn invalid(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}
fn number<T: std::str::FromStr>(value: &str, label: &str) -> Result<T> {
    value
        .parse()
        .map_err(|_| invalid(format!("invalid {label}: {value}")))
}
fn finite(value: &str, label: &str) -> Result<f64> {
    let value: f64 = number(value, label)?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(format!("nonfinite {label}")))
    }
}
// Narrowing to the stored `f32`, refusing a value that does not survive it.
// Load-bearing beyond the narrowing itself: this is what makes every intensity
// a record hands to the kernel validator finite, which is why `Record::finish`
// can skip the per-peak value loop there.
fn intensity(value: f64) -> Result<f32> {
    let result = value as f32;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(invalid("intensity overflows f32"))
    }
}
fn attributes(
    element: &BytesStart<'_>,
    decoder: quick_xml::encoding::Decoder,
    mut budget: Option<&mut ParameterBudget>,
) -> Result<BTreeMap<String, String>> {
    if let Some(budget) = &mut budget {
        budget.begin()?;
    }
    let mut result = BTreeMap::new();
    // The map detects duplicates logarithmically; quick-xml's optional duplicate
    // scan would repeatedly compare all preceding attribute names.
    for attribute in element.attributes().with_checks(false) {
        let attribute = attribute.map_err(|e| invalid(e.to_string()))?;
        if attribute.value.contains(&b'<') {
            return Err(invalid("raw less-than sign in XML attribute"));
        }
        if let Some(budget) = &mut budget {
            // Unescaping XML character references cannot enlarge the UTF-8 data.
            budget.attribute(attribute.key.as_ref().len(), attribute.value.len())?;
        }
        let key = std::str::from_utf8(attribute.key.as_ref())
            .map_err(|e| invalid(e.to_string()))?
            .to_owned();
        let value = attribute
            .decode_and_unescape_value(decoder)
            .map_err(|e| invalid(e.to_string()))?
            .into_owned();
        xml_string(&value)?;
        if result.insert(key, value).is_some() {
            return Err(invalid("duplicate XML attribute"));
        }
    }
    Ok(result)
}
fn required<'a>(attrs: &'a BTreeMap<String, String>, key: &str) -> Result<&'a str> {
    attrs
        .get(key)
        .map(String::as_str)
        .ok_or_else(|| invalid(format!("missing {key} attribute")))
}
fn seconds(unit: Option<&String>) -> Result<f64> {
    match unit.map(String::as_str) {
        Some("UO:0000010") => Ok(1.0),
        Some("UO:0000031") => Ok(60.0),
        other => Err(Error::Unsupported(format!(
            "time unit {other:?}; explicit seconds/minutes required"
        ))),
    }
}

#[derive(Clone, PartialEq, Eq)]
enum Kind {
    Mz,
    Time,
    Intensity,
    Auxiliary(String),
    Role(record_transport::Role),
}
#[derive(Clone, Copy, PartialEq, Eq)]
enum Encoding {
    Float32,
    Float64,
    Int32,
    Int64,
    Ascii,
}
impl Encoding {
    fn width(self) -> Option<usize> {
        match self {
            Self::Float32 | Self::Int32 => Some(4),
            Self::Float64 | Self::Int64 => Some(8),
            Self::Ascii => None,
        }
    }
}
enum Values {
    Floats(Vec<f64>),
    Integers(Vec<i32>),
    Strings(Vec<String>),
}
#[derive(Default)]
struct Binary {
    spectrum: bool,
    metadata: crate::metadata::MetaInfo,
    data_processing: Vec<std::sync::Arc<crate::metadata::DataProcessing>>,
    kind: Option<Kind>,
    encoding: Option<Encoding>,
    compressed: Option<bool>,
    numpress: Option<NumpressCompression>,
    // Base64 characters, held as bytes: they are ASCII by construction, and a
    // byte buffer is what both the strict decoder and the Numpress coder read.
    encoded: Vec<u8>,
    encoded_length: usize,
    array_length: Option<usize>,
    time_scale: f64,
    has_binary: bool,
}
// Scratch buffers reused across binary data arrays. One profile spectrum's
// base64 text and its decoded bytes are each a few hundred kilobytes, well over
// glibc's 128 KiB `mmap` threshold, so allocating them per array costs a fresh
// `mmap`, the matching `munmap` and a page fault per touched page, for every
// array of every record. Reuse holds one buffer at the high-water mark of the
// arrays actually read. Nothing here is sized from a declared length: an array
// can only enlarge a buffer to what its own payload already occupies, or to the
// upper bound the base64 decoder derives from that payload.
#[derive(Default)]
struct Buffers {
    /// Base64 characters of the array being read, lent to its `Binary`.
    encoded: Vec<u8>,
    /// What the base64 decoder wrote, of which only a prefix is meaningful.
    decoded: Vec<u8>,
    /// What zlib expanded out of the decoded bytes.
    inflated: Vec<u8>,
}
impl Buffers {
    /// An empty buffer for the next array's base64 text, keeping whatever
    /// capacity earlier arrays grew it to.
    fn encoded(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.encoded)
    }
    /// Takes a finished array's text buffer back, keeping the larger of the two
    /// so that an interleaved small array cannot shrink the pool.
    fn recycle(&mut self, mut buffer: Vec<u8>) {
        if buffer.capacity() >= self.encoded.capacity() {
            buffer.clear();
            self.encoded = buffer;
        }
    }
}
// Canonical non-primary descendants of binary data array in pinned PSI-MS.
// A zero type mask means this term declares no binary-data-type restriction.
const CANONICAL_ARRAYS: &[(&str, &str, u8)] = &[
    ("MS:1000516", "charge array", 1),
    ("MS:1000517", "signal to noise array", 10),
    ("MS:1000617", "wavelength array", 10),
    ("MS:1000820", "flow rate array", 10),
    ("MS:1000821", "pressure array", 10),
    ("MS:1000822", "temperature array", 10),
    ("MS:1002477", "mean ion mobility drift time array", 0),
    ("MS:1002478", "mean charge array", 2),
    ("MS:1002529", "resolution array", 10),
    ("MS:1002530", "baseline array", 10),
    ("MS:1002742", "noise array", 10),
    ("MS:1002743", "sampled noise m/z array", 10),
    ("MS:1002744", "sampled noise intensity array", 10),
    ("MS:1002745", "sampled noise baseline array", 10),
    ("MS:1002816", "mean ion mobility array", 0),
    ("MS:1002893", "ion mobility array", 0),
    ("MS:1003006", "mean inverse reduced ion mobility array", 0),
    ("MS:1003007", "raw ion mobility array", 0),
    ("MS:1003008", "raw inverse reduced ion mobility array", 0),
    ("MS:1003143", "mass array", 10),
    ("MS:1003153", "raw ion mobility drift time array", 0),
    ("MS:1003154", "deconvoluted ion mobility array", 0),
    (
        "MS:1003155",
        "deconvoluted inverse reduced ion mobility array",
        0,
    ),
    (
        "MS:1003156",
        "deconvoluted ion mobility drift time array",
        0,
    ),
    (
        "MS:1003157",
        "scanning quadrupole position lower bound m/z array",
        0,
    ),
    (
        "MS:1003158",
        "scanning quadrupole position upper bound m/z array",
        0,
    ),
];
fn canonical_array_name(accession: &str) -> Option<&'static str> {
    CANONICAL_ARRAYS
        .iter()
        .find(|t| t.0 == accession)
        .map(|t| t.1)
}
fn canonical_array_accession(name: &str) -> Option<&'static str> {
    CANONICAL_ARRAYS.iter().find(|t| t.1 == name).map(|t| t.0)
}
fn check_canonical_encoding(name: &str, encoding: Encoding) -> Result<()> {
    if let Some((_, _, mask)) = CANONICAL_ARRAYS.iter().find(|t| t.1 == name) {
        let bit = match encoding {
            Encoding::Int32 => 1,
            Encoding::Float32 => 2,
            Encoding::Int64 => 4,
            Encoding::Float64 => 8,
            Encoding::Ascii => 16,
        };
        if *mask != 0 && mask & bit == 0 {
            return Err(Error::Unsupported(
                "canonical auxiliary array binary type".into(),
            ));
        }
    }
    Ok(())
}
/// Whether a binary-array parameter carries any unit attribute.
fn has_unit_attributes(attrs: &BTreeMap<String, String>) -> bool {
    ["unitAccession", "unitCvRef", "unitName"]
        .iter()
        .any(|key| attrs.contains_key(*key))
}
fn integer_array_encoding(name: &str) -> Encoding {
    if name == "charge array" {
        Encoding::Int32
    } else {
        Encoding::Int64
    }
}

impl Binary {
    fn cv(&mut self, attrs: &BTreeMap<String, String>, budget: &mut ParameterBudget) -> Result<()> {
        let accession = required(attrs, "accession")?;
        let detected_role = record_transport::array_role(self, attrs, budget)?;
        let kind = if detected_role.is_some() {
            detected_role
        } else {
            match accession {
                "MS:1000514" => Some(Kind::Mz),
                "MS:1000515" => Some(Kind::Intensity),
                "MS:1000595" => {
                    self.time_scale = seconds(attrs.get("unitAccession"))?;
                    Some(Kind::Time)
                }
                "MS:1000786" => {
                    let name = required(attrs, "value")?;
                    if canonical_array_accession(name).is_some() {
                        return Err(Error::Unsupported(
                            "nonstandard array name collides with a canonical array".into(),
                        ));
                    }
                    if name.is_empty() {
                        return Err(invalid("auxiliary array name is empty"));
                    }
                    Some(Kind::Auxiliary(name.to_owned()))
                }
                _ => canonical_array_name(accession).map(|name| Kind::Auxiliary(name.into())),
            }
        };
        if let Some(Kind::Auxiliary(name)) = &kind {
            if has_unit_attributes(attrs) {
                self.mobility_array_unit(name, attrs, budget)?;
            }
        }
        if let Some(kind) = kind {
            if self.kind.replace(kind).is_some() {
                return Err(invalid("multiple binary array types"));
            }
        } else {
            if ["unitAccession", "unitCvRef", "unitName"]
                .iter()
                .any(|key| attrs.contains_key(*key))
            {
                return Err(Error::Unsupported(
                    "units on binary encoding/compression terms are not represented".into(),
                ));
            }
            let encoding = match accession {
                "MS:1000521" => Some(Encoding::Float32),
                "MS:1000523" => Some(Encoding::Float64),
                "MS:1000519" => Some(Encoding::Int32),
                "MS:1000522" => Some(Encoding::Int64),
                "MS:1001479" => Some(Encoding::Ascii),
                _ => None,
            };
            if let Some(encoding) = encoding {
                if self.encoding.replace(encoding).is_some() {
                    return Err(invalid("multiple binary precisions"));
                }
            } else if matches!(accession, "MS:1000574" | "MS:1000576") {
                if self.compressed.replace(accession == "MS:1000574").is_some()
                    || (accession == "MS:1000576" && self.numpress.is_some())
                {
                    return Err(invalid("conflicting or duplicate binary compression terms"));
                }
            } else if let Some((mode, zlib)) = numpress_transport::compression(accession) {
                if self.numpress.replace(mode).is_some()
                    || self.compressed == Some(false)
                    || (zlib && self.compressed.replace(true).is_some())
                {
                    return Err(invalid(
                        "conflicting or duplicate Numpress compression terms",
                    ));
                }
            } else {
                return Err(Error::Unsupported(format!("binary array CV {accession}")));
            }
        }
        Ok(())
    }
    /// Retain the unit of an ion-mobility array as the source's `unit_accession`
    /// array metadata (`MzMLHandlerHelper.cpp:291-295`), which the writer
    /// emits again and the reader restores.
    ///
    /// Only arrays that `MSSpectrum::contains_im_data` classifies as ion
    /// mobility keep a unit. The source stores a unit on every non-default
    /// array; other auxiliary arrays stay an explicit `Unsupported` error in
    /// this port, as before.
    fn mobility_array_unit(
        &mut self,
        name: &str,
        attrs: &BTreeMap<String, String>,
        budget: &mut ParameterBudget,
    ) -> Result<()> {
        let probe = MSSpectrum {
            float_data_arrays: vec![DataArray::new(name, Vec::new())],
            ..Default::default()
        };
        if !probe.contains_im_data() {
            return Err(Error::Unsupported(
                "units on auxiliary arrays are not represented".into(),
            ));
        }
        let accession = attrs
            .get("unitAccession")
            .ok_or_else(|| invalid("array unit attributes require an accession"))?;
        let prefix = accession
            .split_once(':')
            .map(|(prefix, _)| prefix)
            .ok_or_else(|| invalid("invalid array unit accession"))?;
        if !matches!(prefix, "MS" | "UO") {
            return Err(Error::Unsupported(
                "ion mobility array units require an MS or UO accession".into(),
            ));
        }
        if attrs
            .get("unitCvRef")
            .is_some_and(|cv_ref| cv_ref != prefix)
        {
            return Err(invalid("array unit CV identity mismatch"));
        }
        record_transport::slot(budget, "unit_accession")?;
        if self
            .metadata
            .insert("unit_accession".into(), accession.clone().into())
            .is_some()
        {
            return Err(invalid("duplicate array unit metadata"));
        }
        Ok(())
    }
    // Appends one `<binary>` text node's base64 characters. The source helper
    // `MzMLHandlerHelper::decodeBase64Arrays` removes the four XML whitespace
    // characters unless `skip_xml_checks` bypasses that step; every other
    // character is kept verbatim and left for the strict decoder to reject.
    //
    // A profile spectrum arrives as a few hundred kilobytes in one text node,
    // so the leading run of ordinary base64 characters is located with a single
    // vectorised scan and copied with one `extend_from_slice`, instead of being
    // tested and pushed one character at a time. The checks keep the order they
    // had per character: a run is charged against `encodedLength` before it is
    // copied, so the first byte outside the run is still the one that decides
    // which error a malformed node reports.
    fn append_encoded(&mut self, text: &[u8], normalize: bool, fill_data: bool) -> Result<()> {
        let mut rest = text;
        while !rest.is_empty() {
            let (run, tail) = rest.split_at(run_inside(rest, 0x21, 0x5e));
            if !run.is_empty() {
                if fill_data {
                    if self.encoded.len().saturating_add(run.len()) > self.encoded_length {
                        return Err(invalid("binary text exceeds encodedLength"));
                    }
                    self.encoded.extend_from_slice(run);
                }
                rest = tail;
                continue;
            }
            let Some((&byte, tail)) = tail.split_first() else {
                break;
            };
            if !byte.is_ascii() {
                return Err(invalid("non-ASCII base64 text"));
            }
            let skipped = normalize && matches!(byte, b' ' | b'\t' | b'\n' | b'\r');
            if fill_data && !skipped {
                if self.encoded.len() >= self.encoded_length {
                    return Err(invalid("binary text exceeds encodedLength"));
                }
                self.encoded.push(byte);
            }
            rest = tail;
        }
        Ok(())
    }
    // No payload access or allocation: shared by decoding and fill_data=false.
    fn descriptor(&self, default_count: usize) -> Result<(Encoding, bool, usize)> {
        if !self.has_binary {
            return Err(invalid("missing binary element"));
        }
        let kind = self
            .kind
            .as_ref()
            .ok_or_else(|| invalid("missing binary array type"))?;
        let encoding = match (self.numpress, self.encoding) {
            (Some(_), None | Some(Encoding::Float32 | Encoding::Float64))
            | (Some(NumpressCompression::Pic), Some(Encoding::Int32 | Encoding::Int64)) => {
                Encoding::Float64
            }
            (Some(_), _) => {
                return Err(invalid(
                    "Numpress requires floating data (or PIC integer repair)",
                ));
            }
            (None, Some(encoding)) => encoding,
            (None, None) => return Err(invalid("missing binary precision")),
        };
        let compressed = self
            .compressed
            .or(self.numpress.map(|_| false))
            .ok_or_else(|| invalid("missing binary compression term"))?;
        if let Kind::Auxiliary(name) = kind {
            check_canonical_encoding(name, encoding)?;
        }
        let auxiliary = matches!(kind, Kind::Auxiliary(_));
        let empty_wavelength = matches!(kind, Kind::Role(record_transport::Role::Wavelength))
            && self.array_length == Some(0);
        if !auxiliary && !matches!(encoding, Encoding::Float32 | Encoding::Float64) {
            return Err(Error::Unsupported(
                "primary peak arrays must use floating-point encoding".into(),
            ));
        }
        let count = self.array_length.unwrap_or(default_count);
        let valid_length = count == default_count
            || (auxiliary && count == 0)
            || empty_wavelength
            || matches!(kind, Kind::Role(role) if role.noise());
        if !valid_length {
            return Err(invalid(
                "nonempty arrayLength differs from defaultArrayLength",
            ));
        }
        Ok((encoding, compressed, count))
    }
    fn decode(
        &mut self,
        default_count: usize,
        options: &ReadOptions,
        remaining_bytes: &mut usize,
        remaining_elements: &mut usize,
        numpress_work: &mut coder::Work,
        buffers: &mut Buffers,
    ) -> Result<(Kind, Values)> {
        let (encoding, compressed, count) = self.descriptor(default_count)?;
        let kind = self.kind.take().unwrap(); // validated above
        *remaining_elements = remaining_elements
            .checked_sub(count)
            .ok_or_else(|| invalid("total binary element limit exceeded"))?;
        let byte_limit = options.max_array_bytes.min(*remaining_bytes);
        let expected = encoding
            .width()
            .map(|width| {
                count
                    .checked_mul(width)
                    .filter(|&n| n <= byte_limit)
                    .ok_or_else(|| invalid("binary array exceeds configured byte limit"))
            })
            .transpose()?;
        // Every string has at least its NUL terminator. Check before allocating.
        if expected.is_none() && count > byte_limit {
            return Err(invalid("string array exceeds configured byte limit"));
        }
        if self.encoded.len() != self.encoded_length {
            return Err(invalid(
                "encodedLength does not match base64 character count",
            ));
        }
        if let Some(mode) = self.numpress {
            // Unlike the standalone wrapper's source short-text no-op, mzML
            // requires a valid base64 payload and its exact declared point count.
            if (compressed && self.encoded.is_empty())
                || (!self.encoded.is_empty() && self.encoded.len() < 4)
            {
                return Err(invalid("truncated Numpress base64"));
            }
            numpress_work.limits.raw.max_values = count;
            let config = NumpressConfig {
                compression: mode,
                ..NumpressConfig::default()
            };
            let text =
                std::str::from_utf8(&self.encoded).map_err(|_| invalid("non-ASCII base64 text"))?;
            let mut values = coder::decode_text(text, compressed, &config, numpress_work)?;
            if values.len() != count {
                return Err(invalid("Numpress point count differs from declared length"));
            }
            let expected = expected.expect("Numpress effective type is f64");
            *remaining_bytes = remaining_bytes
                .checked_sub(expected)
                .ok_or_else(|| invalid("total binary byte limit exceeded"))?;
            numpress_work.spend(count)?;
            for value in &mut values {
                if kind == Kind::Time {
                    *value *= self.time_scale;
                }
                // Load-bearing beyond this function: `Record::finish` skips the
                // per-peak value loop of the record validation because this
                // guard and its `Float32`/`Float64` twin below have already
                // refused every nonfinite coordinate. Do not weaken either
                // without restoring the full `validate` there.
                if !value.is_finite() {
                    return Err(invalid("nonfinite Numpress value"));
                }
            }
            return Ok((kind, Values::Floats(values)));
        }
        // Decoding into a reused buffer sized from `decoded_len_estimate`, an
        // upper bound, rather than into a fresh `Vec`: the estimate makes the
        // too-small arm unreachable, and the buffer skips both the per-array
        // allocation and the zero fill `decode_vec` would pay for it.
        let estimate = base64::decoded_len_estimate(self.encoded.len());
        if buffers.decoded.len() < estimate {
            buffers.decoded.resize(estimate, 0);
        }
        let length = STANDARD
            .decode_slice(&self.encoded, &mut buffers.decoded)
            .map_err(|e| match e {
                base64::DecodeSliceError::DecodeError(e) => invalid(format!("invalid base64: {e}")),
                base64::DecodeSliceError::OutputSliceTooSmall => {
                    invalid("invalid base64: decoded array exceeds its buffer")
                }
            })?;
        let bytes = buffers.decoded.get(..length).unwrap_or_default();
        if bytes.len() > options.max_array_bytes {
            return Err(invalid("encoded array exceeds configured byte limit"));
        }
        let limit = expected.unwrap_or(byte_limit);
        let decoded: &[u8] = if compressed {
            // Bounded chunks also handle string arrays, whose byte lengths are
            // variable. Reject concatenated/truncated streams and expansion past
            // either the declared numeric length or cumulative storage budget.
            let mut decoder = Decompress::new(true);
            let output = &mut buffers.inflated;
            output.clear();
            let mut chunk = [0u8; 8192];
            loop {
                let before_in = decoder.total_in();
                let before_out = decoder.total_out();
                let capacity = limit
                    .saturating_sub(output.len())
                    .saturating_add(1)
                    .min(chunk.len());
                let status = decoder
                    .decompress(
                        &bytes[before_in as usize..],
                        &mut chunk[..capacity],
                        FlushDecompress::None,
                    )
                    .map_err(|e| invalid(format!("invalid zlib array: {e}")))?;
                let written = (decoder.total_out() - before_out) as usize;
                if output.len().checked_add(written).is_none_or(|n| n > limit) {
                    return Err(invalid("decoded binary array exceeds byte limit"));
                }
                output.extend_from_slice(&chunk[..written]);
                if status == Status::StreamEnd {
                    if decoder.total_in() != bytes.len() as u64 {
                        return Err(invalid("trailing compressed data"));
                    }
                    break;
                }
                if decoder.total_in() == before_in && written == 0 {
                    return Err(invalid("incomplete zlib stream"));
                }
            }
            output.as_slice()
        } else {
            bytes
        };
        if decoded.len() > limit || expected.is_some_and(|n| decoded.len() != n) {
            return Err(invalid(
                "binary byte length differs from declared length or exceeds limit",
            ));
        }
        *remaining_bytes = remaining_bytes
            .checked_sub(decoded.len())
            .ok_or_else(|| invalid("total binary byte limit exceeded"))?;
        let values = match encoding {
            Encoding::Ascii => {
                if !decoded.is_ascii() || (!decoded.is_empty() && decoded.last() != Some(&0)) {
                    return Err(invalid(
                        "string array must contain NUL-terminated ASCII strings",
                    ));
                }
                if decoded.iter().filter(|&&b| b == 0).count() != count {
                    return Err(invalid(
                        "string array element count differs from declared length",
                    ));
                }
                let strings = if decoded.is_empty() {
                    Vec::new()
                } else {
                    decoded[..decoded.len() - 1]
                        .split(|&b| b == 0)
                        .map(|value| String::from_utf8(value.to_vec()).expect("validated ASCII"))
                        .collect()
                };
                Values::Strings(strings)
            }
            Encoding::Int32 | Encoding::Int64 => {
                let width = encoding.width().unwrap();
                let mut values = Vec::with_capacity(count);
                for chunk in decoded.chunks_exact(width) {
                    let value = if width == 4 {
                        i64::from(i32::from_le_bytes(chunk.try_into().unwrap()))
                    } else {
                        i64::from_le_bytes(chunk.try_into().unwrap())
                    };
                    values.push(
                        i32::try_from(value)
                            .map_err(|_| invalid("integer array value overflows native i32"))?,
                    );
                }
                Values::Integers(values)
            }
            Encoding::Float32 | Encoding::Float64 => {
                let width = encoding.width().unwrap();
                let scale = if kind == Kind::Time {
                    self.time_scale
                } else {
                    1.0
                };
                // Source `MzMLHandlerHelper::decodeBase64Arrays` applies the
                // minute multiplier in place on a `std::vector<float>`, so a
                // converted 32-bit time array keeps only 32-bit precision;
                // see `ReadOptions::source_time_array_precision`.
                let narrow = options.source_time_array_precision
                    && kind == Kind::Time
                    && width == 4
                    && scale != 1.0;
                let mut values = Vec::with_capacity(count);
                for chunk in decoded.chunks_exact(width) {
                    let mut value = if width == 4 {
                        f64::from(f32::from_le_bytes(chunk.try_into().unwrap()))
                    } else {
                        f64::from_le_bytes(chunk.try_into().unwrap())
                    } * scale;
                    if narrow {
                        value = f64::from(value as f32);
                    }
                    // Load-bearing beyond this function: see the Numpress guard
                    // above. Every coordinate a record hands to the kernel
                    // validator passes through here or through that one, which
                    // is why `Record::finish` can skip the per-peak loop.
                    if !value.is_finite()
                        && !(options.source_nonfinite_float_arrays
                            && matches!(kind, Kind::Auxiliary(_)))
                    {
                        return Err(invalid("nonfinite binary value"));
                    }
                    values.push(value);
                }
                Values::Floats(values)
            }
        };
        Ok((kind, values))
    }
}

struct Record {
    spectrum: Option<MSSpectrum>,
    chromatogram: Option<MSChromatogram>,
    count: usize,
    array_names: BTreeSet<String>,
    primary_seen: [bool; 2],
    coordinates: Option<Vec<f64>>,
    intensities: Option<Vec<f64>>,
    precursor: Option<Precursor>,
    selected_ion: bool,
    ignore_selected_ion: bool,
    precursor_mz_selected_ion: bool,
    selected_mz: bool,
    rt_seen: bool,
    precursors_seen: usize,
    name_seen: bool,
    seen_fields: BTreeSet<&'static str>,
    precursor_fields: BTreeSet<&'static str>,
    precursor_list_seen: bool,
    selected_list_seen: bool,
    scan_list_seen: bool,
    product_active: bool,
    product_seen: bool,
    product_fields: BTreeSet<&'static str>,
    product_list_seen: bool,
    product: Option<Product>,
    scan_active: bool,
    scan_window_list_seen: bool,
    scan_window: Option<crate::metadata::ScanWindow>,
    scan_window_fields: BTreeSet<&'static str>,
    scan_window_unit: Option<String>,
    precursor_filter: Option<crate::kernel::NumericRange>,
    precursor_outside: bool,
    raw_float_arrays: Vec<DataArray<f64>>,
    /// [`ReadOptions::source_nonfinite_float_arrays`].
    nonfinite_float_arrays: bool,
    primary_metadata: [crate::metadata::MetaInfo; 2],
    intensity_first: bool,
    wavelength: Option<(usize, usize)>,
}
impl Record {
    fn check_array_kind(&mut self, kind: &Kind) -> Result<()> {
        let index = match kind {
            Kind::Auxiliary(name) => {
                if !self.array_names.insert(name.clone()) {
                    return Err(invalid("duplicate auxiliary array name"));
                }
                return Ok(());
            }
            Kind::Role(role) if *role == record_transport::Role::Wavelength || role.noise() => {
                if !self.array_names.insert(role.terms().1.into()) {
                    return Err(invalid("duplicate wavelength/noise array"));
                }
                return Ok(());
            }
            Kind::Role(_) | Kind::Intensity => 1,
            Kind::Mz if self.spectrum.is_some() => 0,
            Kind::Time if self.chromatogram.is_some() => 0,
            _ => return Err(invalid("coordinate array has wrong type for record")),
        };
        if std::mem::replace(&mut self.primary_seen[index], true) {
            return Err(invalid("duplicate coordinate/intensity array"));
        }
        Ok(())
    }

    fn metadata(&mut self) -> &mut crate::metadata::MetaInfo {
        if let Some(s) = &mut self.spectrum {
            &mut s.metadata
        } else {
            &mut self.chromatogram.as_mut().unwrap().metadata
        }
    }
    /// Spectrum- and scan-level ion mobility: `MS:1001581` directly below the
    /// spectrum (`MzMLHandler.cpp:1731-1738`) and the four mobility terms below
    /// a scan (2279-2311) set the spectrum's drift time and unit.
    ///
    /// The accession/unit table is the selected-ion one. A unit attribute must
    /// name that quantity's unit, except that the legacy FAIMS spelling
    /// [`LEGACY_FAIMS_VOLT`] is read as volts; the source ignores unit
    /// attributes entirely. An identical repeat is accepted. A conflicting
    /// repeat is `Unsupported`, where the source silently keeps the last value.
    fn spectrum_mobility(
        &mut self,
        accession: &str,
        value: &str,
        attrs: &BTreeMap<String, String>,
    ) -> Result<()> {
        let Some((unit, expected)) = precursor_metadata::mobility_term(accession) else {
            return Ok(());
        };
        if let Some(given) = attrs.get("unitAccession") {
            let legacy =
                unit == DriftTimeUnit::FaimsCompensationVoltage && given == LEGACY_FAIMS_VOLT;
            if given != expected && !legacy {
                return Err(Error::Unsupported(
                    "spectrum mobility unit does not match its typed quantity".into(),
                ));
            }
        }
        let drift_time = finite(value, "spectrum mobility")?;
        let spectrum = self
            .spectrum
            .as_mut()
            .ok_or_else(|| invalid("mobility outside spectrum"))?;
        if !self.seen_fields.insert("spectrum_mobility") {
            if spectrum.drift_time != drift_time || spectrum.drift_time_unit != unit {
                return Err(Error::Unsupported(
                    "conflicting spectrum ion mobility values".into(),
                ));
            }
            return Ok(());
        }
        spectrum.drift_time = drift_time;
        spectrum.drift_time_unit = unit;
        Ok(())
    }
    fn cv(
        &mut self,
        parent: &str,
        attrs: &BTreeMap<String, String>,
        budget: &mut ParameterBudget,
    ) -> Result<()> {
        let accession = required(attrs, "accession")?;
        let value = attrs.get("value").map(String::as_str).unwrap_or("");
        if self.product_active {
            if parent != "isolationWindow" {
                return Err(invalid("product CV outside isolation window"));
            }
            let field = match accession {
                "MS:1000827" => "target",
                "MS:1000828" => "lower",
                "MS:1000829" => "upper",
                _ => {
                    return Err(Error::Unsupported(format!(
                        "product isolation CV {accession}"
                    )));
                }
            };
            if !self.product_fields.insert(field) {
                return Err(invalid("duplicate product isolation quantity"));
            }
            if attrs
                .get("unitAccession")
                .is_some_and(|u| u != "MS:1000040")
                || attrs.get("unitCvRef").is_some_and(|u| u != "MS")
            {
                return Err(Error::Unsupported(
                    "product isolation quantity requires m/z units".into(),
                ));
            }
            let quantity = finite(value, "product isolation quantity")?;
            let product = self.product.as_mut().expect("active product");
            match field {
                "target" => product.mz = quantity,
                "lower" => product.isolation_window_lower_offset = quantity,
                _ => product.isolation_window_upper_offset = quantity,
            }
            if field != "target" && quantity < 0.0 {
                return Err(invalid("negative product isolation offset"));
            }
            return Ok(());
        }
        if record_transport::read_cv(self, parent, attrs, budget)? {
            return Ok(());
        }
        if acquisition_metadata::read_cv(self, parent, attrs)? {
            return Ok(());
        }
        if parent == "chromatogram"
            && matches!(
                required(attrs, "accession")?,
                "MS:1003019" | "MS:1003020" | "MS:1000626"
            )
        {
            record_transport::slot(budget, "chromatogram type accession")?;
        }
        if settings_metadata::read_cv(self, parent, attrs)? {
            return Ok(());
        }
        let field = match (parent, accession) {
            ("spectrum", "MS:1000511") => Some("ms_level"),
            ("spectrum", "MS:1000127" | "MS:1000128") => Some("spectrum_type"),
            ("selectedIon", "MS:1000744" | "MS:1000040") => Some("precursor_mz"),
            ("selectedIon", "MS:1000041") => Some("precursor_charge"),
            ("selectedIon", "MS:1000042") => Some("precursor_intensity"),
            ("isolationWindow", "MS:1000827") if self.precursor.is_some() => {
                Some("isolation_target")
            }
            _ if self.precursor.is_some() => precursor_metadata::field(parent, accession),
            _ => None,
        };
        if let Some(field) = field {
            let seen = if matches!(parent, "selectedIon" | "isolationWindow" | "activation") {
                &mut self.precursor_fields
            } else {
                &mut self.seen_fields
            };
            let repeated_mz = !self.precursor_mz_selected_ion
                && matches!(field, "precursor_mz" | "isolation_target");
            if !seen.insert(field) && !repeated_mz {
                return Err(invalid(format!("duplicate scientific field {field}")));
            }
        }
        if let Some(p) = self.precursor.as_mut() {
            if precursor_metadata::read_cv(
                p,
                parent,
                accession,
                value,
                attrs,
                self.selected_mz && self.precursor_mz_selected_ion,
            )? {
                if parent == "isolationWindow"
                    && accession == "MS:1000827"
                    && !self.precursor_mz_selected_ion
                    && self.spectrum.is_some()
                    && self
                        .precursor_filter
                        .is_some_and(|range| !load::contains(range, p.mz))
                {
                    self.precursor_outside = true;
                }
                return Ok(());
            }
        }
        match (parent, accession) {
            ("spectrum", "MS:1000511") => {
                self.spectrum.as_mut().unwrap().ms_level = number(value, "MS level")?
            }
            ("spectrum", "MS:1000127") => {
                self.spectrum.as_mut().unwrap().spectrum_type = SpectrumType::Centroid
            }
            ("spectrum", "MS:1000128") => {
                self.spectrum.as_mut().unwrap().spectrum_type = SpectrumType::Profile
            }
            ("scan", "MS:1000016") => {
                if self.rt_seen {
                    return Err(Error::Unsupported(
                        "multiple scan start times in one spectrum".into(),
                    ));
                }
                let rt = finite(value, "scan start time")? * seconds(attrs.get("unitAccession"))?;
                if !rt.is_finite() {
                    return Err(invalid("scan start time overflow"));
                }
                self.spectrum
                    .as_mut()
                    .ok_or_else(|| invalid("scan outside spectrum"))?
                    .rt = rt;
                self.rt_seen = true;
            }
            ("selectedIon", "MS:1000744" | "MS:1000040") => {
                let mz = finite(value, "precursor m/z")?;
                if self.precursor_mz_selected_ion
                    && self.spectrum.is_some()
                    && self.precursor.as_ref().is_some_and(|p| p.mz != mz)
                    && self
                        .precursor_filter
                        .is_some_and(|range| !load::contains(range, mz))
                {
                    self.precursor_outside = true;
                }
                let p = self
                    .precursor
                    .as_mut()
                    .ok_or_else(|| invalid("selected ion outside precursor"))?;
                if self.precursor_mz_selected_ion {
                    p.mz = mz;
                    self.selected_mz = true;
                } else if p.mz != mz {
                    // Source stores the first selected ion separately in target
                    // mode. Equal later events do not erase an earlier value.
                    let slot = std::mem::size_of::<(String, MetaValue)>() + 32;
                    budget.spend(512usize.saturating_add(11 * slot).saturating_add(17))?;
                    p.cv_terms
                        .metadata
                        .insert("selected ion m/z".into(), MetaValue::try_from(mz)?);
                }
            }
            ("selectedIon", "MS:1000041") => {
                self.precursor
                    .as_mut()
                    .ok_or_else(|| invalid("selected ion outside precursor"))?
                    .charge = number(value, "precursor charge")?
            }
            ("selectedIon", "MS:1000042") => {
                self.precursor
                    .as_mut()
                    .ok_or_else(|| invalid("selected ion outside precursor"))?
                    .intensity = intensity(finite(value, "precursor intensity")?)?
            }
            ("spectrum", "MS:1000525") => {
                // MzMLHandler.cpp:1642-1645: the generic representation term
                // resets an earlier centroid or profile term to unknown.
                self.spectrum
                    .as_mut()
                    .ok_or_else(|| invalid("spectrum representation outside spectrum"))?
                    .spectrum_type = SpectrumType::Unknown
            }
            ("spectrum", "MS:1001581")
            | ("scan", "MS:1002476" | "MS:1002815" | "MS:1001581" | "MS:1002954") => {
                self.spectrum_mobility(accession, value, attrs)?
            }
            _ => {} // Acquisition CVs outside the supported model are intentionally not retained.
        }
        Ok(())
    }
    fn finish(
        mut self,
        selection: Option<&mut load::State<'_>>,
        fill_data: bool,
    ) -> Result<Option<consumer::Completed>> {
        let (mut positions, mut intensities) = if !fill_data {
            (Vec::new(), Vec::new())
        } else {
            match (self.coordinates.take(), self.intensities.take()) {
                (Some(p), Some(i)) => (p, i),
                (None, None) if self.count == 0 => (Vec::new(), Vec::new()),
                _ => return Err(invalid("missing coordinate or intensity array")),
            }
        };
        let keep = if let Some(selection) = selection {
            selection.apply(&mut self, &mut positions, &mut intensities)?
        } else {
            true
        };
        let floats = if let Some(spectrum) = &mut self.spectrum {
            &mut spectrum.float_data_arrays
        } else {
            &mut self.chromatogram.as_mut().unwrap().float_data_arrays
        };
        let lenient = self.nonfinite_float_arrays;
        for array in self.raw_float_arrays {
            floats.push(DataArray {
                name: array.name,
                data: array
                    .data
                    .into_iter()
                    // The source narrows with `static_cast<float>`.
                    .map(|value| {
                        if lenient {
                            Ok(value as f32)
                        } else {
                            intensity(value)
                        }
                    })
                    .collect::<Result<_>>()?,
                metadata: array.metadata,
                data_processing: array.data_processing,
            });
        }
        if let Some(mut spectrum) = self.spectrum {
            // Written as a loop rather than a `collect::<Result<_>>`: the
            // fallible collect carries its error through every step of the
            // iterator pipeline, where one reserved buffer and a plain push
            // pair each coordinate with its intensity in the same order and
            // stop at the same first bad value.
            let mut peaks = Vec::with_capacity(positions.len().min(intensities.len()));
            for (mz, i) in positions.into_iter().zip(intensities) {
                peaks.push(Peak1D::new(mz, intensity(i)?));
            }
            spectrum.peaks = peaks;
            // Every peak value of this record is already known finite, so the
            // per-peak loop of `MSSpectrum::validate` cannot fail here and is
            // skipped; everything else that `validate` checks still runs. The
            // argument, for both fields:
            //
            // * `peak.intensity` is the value `intensity` returned three lines
            //   above, and `intensity` (this module, `fn intensity`) refuses
            //   every nonfinite `f32` as it narrows.
            // * `peak.mz` comes from `positions`, which is `self.coordinates`.
            //   That field is only ever assigned a decoded `Values::Floats` --
            //   from the primary-array slot in `apply_binary_array`, or from
            //   the wavelength array promoted to the coordinate slot -- and
            //   `Binary::decode` rejects every nonfinite decoded float before
            //   it is pushed, in both of its float paths: the `Float32`/
            //   `Float64` branch after the time scale and the optional 32-bit
            //   narrowing (`nonfinite binary value`), and the Numpress branch
            //   after its own time scale (`nonfinite Numpress value`).
            //   `load::State::apply` then only filters, permutes and truncates
            //   `positions`; it never produces a value.
            //
            // So this skips a full scan of the record's peaks, on every record
            // of every mzML read, for a verdict the decoder has already given
            // per value. `tests/mzml.rs::nonfinite_peak_values_are_refused_by_
            // the_decoder_not_the_validator` pins the premise: if a decode path
            // ever stops rejecting nonfinite values, that test fails here
            // rather than this skipped check failing to catch it downstream.
            spectrum.validate_given_finite_peaks()?;
            Ok(keep.then_some(consumer::Completed::Spectrum(spectrum)))
        } else {
            let mut chromatogram = self.chromatogram.unwrap();
            let mut peaks = Vec::with_capacity(positions.len().min(intensities.len()));
            for (rt, i) in positions.into_iter().zip(intensities) {
                peaks.push(ChromatogramPeak::new(rt, intensity(i)?));
            }
            chromatogram.peaks = peaks;
            // The same argument as for the spectrum above, with `peak.rt` in
            // place of `peak.mz`: chromatogram coordinates reach `positions`
            // through the same `Kind::Time` primary-array slot, so they are
            // decoded floats that passed the same guard.
            chromatogram.validate_given_finite_peaks()?;
            Ok(keep.then_some(consumer::Completed::Chromatogram(chromatogram)))
        }
    }
}

// Inline and referenced parameters take exactly the same scientific path.
#[allow(clippy::too_many_arguments)]
fn apply_parameter(
    tag: &str,
    attrs: &BTreeMap<String, String>,
    parent: &str,
    record: &mut Option<Record>,
    binary: &mut Option<Binary>,
    experiment: &mut MSExperiment,
    selection: Option<&mut load::State<'_>>,
    budget: &mut ParameterBudget,
) -> Result<()> {
    required(
        attrs,
        if tag == "cvParam" {
            "accession"
        } else {
            "name"
        },
    )?;
    // Later selected ions are scientifically ignored by the source. Attribute,
    // XML, group-reference and expansion checks have already run for this event.
    if parent == "selectedIon" && record.as_ref().is_some_and(|r| r.ignore_selected_ion) {
        return Ok(());
    }
    if tag == "cvParam" {
        if let Some(selection) = selection {
            let selected = selection.options.scientific.precursor_mz_selected_ion;
            let accession = attrs.get("accession").map(String::as_str);
            let target_event = !selected
                && parent == "isolationWindow"
                && accession == Some("MS:1000827")
                && record
                    .as_ref()
                    .is_some_and(|r| r.precursor.is_some() && !r.product_active);
            let selected_event = selected
                && parent == "selectedIon"
                && matches!(accession, Some("MS:1000744" | "MS:1000040"));
            if selection.options.scientific.has_precursor_mz_range()
                && (target_event || selected_event)
            {
                selection.spend(1)?;
            }
        }
    }
    match tag {
        "cvParam" => {
            if parent == "binaryDataArray" {
                binary
                    .as_mut()
                    .ok_or_else(|| invalid("CV outside binary array"))?
                    .cv(attrs, budget)?;
            } else if let Some(r) = record {
                r.cv(parent, attrs, budget)?;
            } else if parent == "run" && required(attrs, "accession")? == "MS:1000858" {
                experiment.settings.fraction_identifier =
                    attrs.get("value").cloned().unwrap_or_default();
            }
        }
        "userParam" if parent == "binaryDataArray" => {
            let name = required(attrs, "name")?;
            let value = product_user_value(attrs)?;
            let b = binary
                .as_mut()
                .ok_or_else(|| invalid("array userParam outside array"))?;
            if b.metadata.insert(name.into(), value).is_some() {
                return Err(invalid("duplicate array metadata key"));
            }
        }
        "userParam" if record.as_ref().is_some_and(|r| r.product_active) => {
            if parent != "isolationWindow" {
                return Err(invalid("product userParam outside isolation window"));
            }
            let name = required(attrs, "name")?;
            let value = product_user_value(attrs)?;
            let metadata = &mut record
                .as_mut()
                .unwrap()
                .product
                .as_mut()
                .unwrap()
                .cv_terms
                .metadata;
            if metadata.insert(name.into(), value).is_some() {
                return Err(invalid("duplicate product userParam name"));
            }
        }
        "userParam" if matches!(parent, "scan" | "scanList") => {
            acquisition_metadata::insert(
                record
                    .as_mut()
                    .ok_or_else(|| invalid("acquisition metadata outside record"))?,
                parent,
                required(attrs, "name")?,
                product_user_value(attrs)?,
            )?;
        }
        "userParam" if parent == "scanWindow" => {
            let name = required(attrs, "name")?;
            // The dedicated unit field is derived from endpoint CV parameters.
            // A competing userParam would silently overwrite that identity.
            if name == "unit_accession" {
                return Err(Error::Unsupported(
                    "scan window unit_accession must use endpoint units".into(),
                ));
            }
            let value = product_user_value(attrs)?;
            let window = record
                .as_mut()
                .and_then(|r| r.scan_window.as_mut())
                .ok_or_else(|| invalid("userParam outside scan window"))?;
            if window.metadata.insert(name.into(), value).is_some() {
                return Err(invalid("duplicate scan window userParam name"));
            }
        }
        "userParam"
            if matches!(parent, "activation" | "isolationWindow" | "selectedIon")
                && record.as_ref().is_some_and(|r| r.precursor.is_some()) =>
        {
            let name = required(attrs, "name")?;
            // The source writer's schema fallback carries no scientific data.
            if parent == "activation"
                && name == "activation information unavailable"
                && attrs.len() == 1
            {
                return Ok(());
            }
            let value = product_user_value(attrs)?;
            let r = record.as_mut().unwrap();
            if name == "peak intensity unit accession"
                && r.precursor_fields.contains("intensity_explicit_unit")
            {
                return Err(invalid("duplicate precursor intensity unit metadata"));
            }
            let p = r.precursor.as_mut().unwrap();
            if p.cv_terms.metadata.insert(name.into(), value).is_some() {
                return Err(invalid("duplicate precursor metadata key"));
            }
        }
        "userParam" if parent == "run" => {
            let name = required(attrs, "name")?;
            if name == NAME_KEY {
                return Err(invalid("reserved record name userParam at run level"));
            }
            let value = product_user_value(attrs)?;
            if experiment
                .settings
                .metadata
                .insert(name.into(), value)
                .is_some()
            {
                return Err(invalid("duplicate run userParam name"));
            }
        }
        "userParam" if matches!(parent, "spectrum" | "chromatogram") => {
            let name = required(attrs, "name")?.to_owned();
            let value = product_user_value(attrs)?;
            record_transport::slot(budget, &name)?;
            if name == NAME_KEY {
                if value.unit().is_some() {
                    return Err(invalid("record name cannot have units"));
                }
                let value = value.as_str()?.to_owned();
                let r = record
                    .as_mut()
                    .ok_or_else(|| invalid("name outside record"))?;
                if r.name_seen {
                    return Err(invalid("duplicate record name userParam"));
                }
                r.name_seen = true;
                if let Some(s) = &mut r.spectrum {
                    s.name = value;
                } else {
                    r.chromatogram.as_mut().unwrap().name = value;
                }
            } else {
                let metadata = record
                    .as_mut()
                    .ok_or_else(|| invalid("metadata outside record"))?
                    .metadata();
                if metadata.insert(name.clone(), value).is_some() {
                    return Err(Error::Unsupported(format!(
                        "duplicate userParam name {name}"
                    )));
                }
            }
        }
        _ => {}
    }
    Ok(())
}

#[derive(Debug)]
struct Parameter {
    tag: &'static str,
    attrs: BTreeMap<String, String>,
}

struct ParameterBudget {
    remaining: usize,
    bytes: usize,
}
impl ParameterBudget {
    fn begin(&mut self) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(1)
            .ok_or_else(|| invalid("parameter count exceeds configured limit"))?;
        // Covers even a sparsely occupied map node and Vec/String slots.
        self.spend(1024)
    }
    fn spend(&mut self, bytes: usize) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| invalid("parameter bytes exceed configured limit"))?;
        Ok(())
    }
    fn attribute(&mut self, key_bytes: usize, value_bytes: usize) -> Result<()> {
        let cost = key_bytes
            .checked_add(value_bytes)
            .and_then(|n| n.checked_mul(2))
            .and_then(|n| n.checked_add(256))
            .ok_or_else(|| invalid("parameter byte size overflow"))?;
        // Includes parsed strings plus any simultaneous retained metadata/ID copy.
        self.spend(cost)
    }
    fn charge(&mut self, attrs: &BTreeMap<String, String>) -> Result<()> {
        self.begin()?;
        for (key, value) in attrs {
            self.attribute(key.len(), value.len())?;
        }
        Ok(())
    }
}

// XML Schema ID/IDREF use NCName, including the XML 1.0 Unicode name ranges.
fn parameter_id(value: &str) -> Result<&str> {
    let value = value.trim_matches([' ', '\t', '\r', '\n']);
    fn start(c: char) -> bool {
        matches!(c, 'A'..='Z' | '_' | 'a'..='z' | '\u{c0}'..='\u{d6}' |
            '\u{d8}'..='\u{f6}' | '\u{f8}'..='\u{2ff}' | '\u{370}'..='\u{37d}' |
            '\u{37f}'..='\u{1fff}' | '\u{200c}'..='\u{200d}' | '\u{2070}'..='\u{218f}' |
            '\u{2c00}'..='\u{2fef}' | '\u{3001}'..='\u{d7ff}' | '\u{f900}'..='\u{fdcf}' |
            '\u{fdf0}'..='\u{fffd}' | '\u{10000}'..='\u{effff}')
    }
    let mut chars = value.chars();
    if !chars.next().is_some_and(start) || !chars.all(|c| start(c) ||
        matches!(c, '-' | '.' | '0'..='9' | '\u{b7}' | '\u{300}'..='\u{36f}' | '\u{203f}'..='\u{2040}')) {
        return Err(invalid("invalid parameter group ID/IDREF"));
    }
    Ok(value)
}

fn parameter_context(parent: &str) -> bool {
    // Every ParamGroupType or extension in mzML 1.1.0, including header contexts
    // whose metadata is outside the current native experiment representation.
    matches!(
        parent,
        "fileContent"
            | "sourceFile"
            | "contact"
            | "sample"
            | "source"
            | "analyzer"
            | "detector"
            | "instrumentConfiguration"
            | "software"
            | "processingMethod"
            | "scanSettings"
            | "target"
            | "run"
            | "scanList"
            | "scan"
            | "scanWindow"
            | "binaryDataArray"
            | "spectrum"
            | "chromatogram"
            | "isolationWindow"
            | "activation"
            | "selectedIon"
    )
}

fn apply_group(
    parameters: &[Parameter],
    parent: &str,
    budget: &mut ParameterBudget,
    record: &mut Option<Record>,
    binary: &mut Option<Binary>,
    experiment: &mut MSExperiment,
    mut selection: Option<&mut load::State<'_>>,
) -> Result<()> {
    for parameter in parameters {
        budget.charge(&parameter.attrs)?;
        apply_parameter(
            parameter.tag,
            &parameter.attrs,
            parent,
            record,
            binary,
            experiment,
            selection.as_deref_mut(),
            budget,
        )?;
    }
    Ok(())
}

/// Read a complete experiment with default resource limits.
pub fn read(reader: impl BufRead) -> Result<MSExperiment> {
    read_with_options(reader, &ReadOptions::default())
}

/// Read the documented mzML subset, rejecting unsupported binary encodings.
/// Errors have line zero when the XML parser cannot supply a line number.
pub fn read_with_options(reader: impl BufRead, options: &ReadOptions) -> Result<MSExperiment> {
    read_impl(reader, options, None, false)
}

/// Execute supported scientific loading choices without changing legacy reads.
/// Unimplemented read behavior is rejected before input. Write-only settings do
/// not affect loading. Excluded records remain fully decoded and validated.
pub fn read_with_load_options(
    reader: impl BufRead,
    load: &LoadOptions,
    limits: &ReadOptions,
) -> Result<MSExperiment> {
    read_impl(reader, limits, Some(load), load.scientific.metadata_only)
}

/// Read complete experimental settings and stop at the first record-list
/// opening tag, before parsing its count or consuming any binary payload.
/// The required default processing reference must resolve. With no record list,
/// the bounded header-only document is read through its ordinary closing tags.
pub fn read_metadata(reader: impl BufRead) -> Result<crate::metadata::ExperimentalSettings> {
    read_metadata_with_options(reader, &ReadOptions::default())
}
/// Read run metadata with caller-supplied XML and parameter limits.
///
/// Like [`read_metadata`], stops at the first record-list opening tag. Invalid
/// header metadata, unresolved references, and exceeded limits return an error.
pub fn read_metadata_with_options(
    reader: impl BufRead,
    options: &ReadOptions,
) -> Result<crate::metadata::ExperimentalSettings> {
    Ok(read_impl(reader, options, None, true)?.settings)
}

fn read_impl(
    reader: impl BufRead,
    options: &ReadOptions,
    load: Option<&LoadOptions>,
    metadata_only: bool,
) -> Result<MSExperiment> {
    read_engine(
        reader,
        options,
        load,
        metadata_only,
        MSExperiment::new(),
        None,
    )
}

fn read_engine(
    reader: impl BufRead,
    options: &ReadOptions,
    load: Option<&LoadOptions>,
    metadata_only: bool,
    mut experiment: MSExperiment,
    mut consumer: Option<&mut consumer::Sink<'_>>,
) -> Result<MSExperiment> {
    let fill_data = load.is_none_or(|o| o.scientific.fill_data);
    // The source option bypasses only the helper's four-character whitespace
    // removal, never XML syntax checks or strict native Base64 validation.
    let normalize_binary = load.is_none_or(|o| !o.scientific.skip_xml_checks);
    let mut selection = if metadata_only {
        None
    } else {
        load.map(load::State::new).transpose()?
    };
    let mut header_draft = header::Draft::default();
    let mut header_registry = header::Registry::default();
    // No absolute ceiling: `options.scaling.metadata_*` bound this allowance.
    let mut header_work = header::Work {
        remaining: usize::MAX,
        bytes: usize::MAX,
    };
    let mut default_processing = Vec::new();
    let limit = options
        .max_xml_bytes
        .checked_add(1)
        .ok_or_else(|| Error::InvalidValue("XML byte limit must be below u64::MAX".into()))?;
    let mut reader = NsReader::from_reader(reader.take(limit));
    reader.config_mut().expand_empty_elements = true;
    reader.config_mut().enable_all_checks(true);
    let mut buffer = Vec::new();
    let mut stack = Vec::<String>::new();
    let mut record: Option<Record> = None;
    let mut binary: Option<Binary> = None;
    let mut buffers = Buffers::default();
    let mut seen_root = false;
    let mut seen_mzml = false;
    let mut seen_run = false;
    // Declared points still available, under both `max_total_peaks` and
    // `scaling.peaks`.
    let mut remaining_peaks = options.max_total_peaks;
    // Records, binary arrays and parameter groups still available, under both
    // their absolute ceiling and `scaling.records`/`arrays`/`param_groups`.
    let mut remaining_records = options.max_records;
    let mut spectrum_list_seen = false;
    let mut chromatogram_list_seen = false;
    let mut array_list_seen = false;
    let mut remaining_arrays = options.max_total_arrays;
    let mut numpress_limits = coder::NumpressCoderLimits::default();
    numpress_limits.raw.max_encoded_bytes = options.max_array_bytes;
    numpress_limits.max_text_bytes = usize::try_from(options.max_xml_bytes).unwrap_or(usize::MAX);
    let mut remaining_array_bytes = options.max_total_array_bytes;
    let mut remaining_array_elements = options.max_total_array_elements;
    let mut ids = BTreeSet::new();
    // Record, binary-array, precursor/selectedIon/scan and product/scanWindow
    // list `count` attributes must be present and numeric, but a value
    // disagreeing with the actual number of children is **advisory on reading**:
    // the source reader never compares the two. `spectrumList` and
    // `chromatogramList` spend the count on a progress range and
    // `reserveSpaceSpectra`/`reserveSpaceChromatograms`
    // (MzMLHandler.cpp:965-979 and :996-1013), `binaryDataArrayList` on
    // `bin_data_.reserve(...)` (MzMLHandler.cpp:1015-1017), and
    // `selectedIonList` only warns when the count exceeds one
    // (MzMLHandler.cpp:1371-1375); `precursorList`, `productList` and
    // `scanWindowList` have no open-tag handler at all. The count attribute is
    // read in exactly those four places (MzMLHandler.cpp:965, :996, :1017,
    // :1374) and compared against nothing.
    //
    // Real files carry wrong counts. The upstream class-test fixture
    // `MzMLFile_1.mzML` declares `<binaryDataArrayList count="2">` with four
    // `binaryDataArray` children, and `MzMLFile_test.cpp` loads it; rejecting
    // the mismatch made this port unable to read its own reference data, the
    // same defect already fixed for header lists in
    // `mzml_header::Node::children` (`DTAExtractor_1_input.mzML` declares
    // `softwareList count="5"` with four entries). Declared counts are still
    // spent as resource ceilings before any allocation, and writing still emits
    // the true count.
    //
    // `referenceableParamGroupList` follows the same rule. An earlier pass kept
    // it strict so that `read` and `mzml_counts::read_size` would agree on a
    // file; both are advisory now, which keeps them agreeing and matches the
    // source, since upstream has a handler for `referenceableParamGroup`
    // (MzMLHandler.cpp:1055) and none for the enclosing list, so it never reads
    // that count at all.
    // (depth of the group children, declared count, groups seen)
    let mut group_list: Option<(usize, usize, usize)> = None;
    let mut groups = BTreeMap::<String, Vec<Parameter>>::new();
    let mut group: Option<(String, Vec<Parameter>)> = None;
    let mut group_list_seen = false;
    let mut remaining_groups = options.max_param_groups;
    let mut parameter_budget = ParameterBudget {
        remaining: options.max_total_params,
        bytes: options.max_param_bytes,
    };
    // Every cumulative counter holds the room left under its absolute ceiling
    // and, reconciled after each event, under its size-derived allowance.
    let scaling = &options.scaling;
    let mut ledger = scaling::Ledger::attach([
        (&mut remaining_peaks, scaling.peaks),
        (&mut remaining_array_bytes, scaling.array_bytes),
        (&mut remaining_array_elements, scaling.array_elements),
        (&mut remaining_records, scaling.records),
        (&mut remaining_arrays, scaling.arrays),
        (&mut remaining_groups, scaling.param_groups),
        (&mut parameter_budget.remaining, scaling.params),
        (&mut parameter_budget.bytes, scaling.param_bytes),
        (&mut header_work.remaining, scaling.metadata_work),
        (&mut header_work.bytes, scaling.metadata_bytes),
    ]);
    let mut selection_ledger = selection.as_mut().map(|state| {
        let [work, bytes] = state.allowances();
        scaling::Ledger::attach([
            (work, scaling.selection_work),
            (bytes, scaling.selection_bytes),
        ])
    });
    let mut seen_declaration = false;
    let mut ascii_only = false;
    let mut latin1_subset = false;
    loop {
        let decoder = reader.decoder();
        let (namespace, event) = reader
            .read_resolved_event_into(&mut buffer)
            .map_err(|e| invalid(e.to_string()))?;
        let namespace_ok = matches!(namespace, ResolveResult::Bound(ns) if ns.as_ref() == NS);
        // Credit the bytes of this event before any of its charges.
        let position = reader.buffer_position();
        ledger.sync(
            position,
            [
                &mut remaining_peaks,
                &mut remaining_array_bytes,
                &mut remaining_array_elements,
                &mut remaining_records,
                &mut remaining_arrays,
                &mut remaining_groups,
                &mut parameter_budget.remaining,
                &mut parameter_budget.bytes,
                &mut header_work.remaining,
                &mut header_work.bytes,
            ],
        );
        if let (Some(ledger), Some(state)) = (selection_ledger.as_mut(), selection.as_mut()) {
            ledger.sync(position, state.allowances());
        }
        if ascii_only && !event.is_ascii() {
            if latin1_subset {
                return Err(Error::Unsupported(
                    "non-ASCII ISO-8859-1 mzML requires transcoding".into(),
                ));
            }
            return Err(invalid("non-ASCII bytes in US-ASCII XML"));
        }
        if reader.buffer_position() > options.max_xml_bytes {
            return Err(invalid("XML exceeds configured byte limit"));
        }
        match event {
            Event::Start(element) => {
                if !namespace_ok {
                    return Err(invalid("element is outside the mzML namespace"));
                }
                let tag = std::str::from_utf8(element.local_name().as_ref())
                    .map_err(|e| invalid(e.to_string()))?
                    .to_owned();
                let is_parameter = !seen_run
                    || matches!(
                        tag.as_str(),
                        "cvParam"
                            | "userParam"
                            | "referenceableParamGroup"
                            | "referenceableParamGroupRef"
                            | "scan"
                            | "sourceFile"
                            | "instrumentConfiguration"
                            | "run"
                    );
                let attrs = attributes(
                    &element,
                    decoder,
                    is_parameter.then_some(&mut parameter_budget),
                )?;
                let parent = stack.last().map(String::as_str).unwrap_or("");
                if stack.is_empty() {
                    if seen_root || !matches!(tag.as_str(), "mzML" | "indexedmzML") {
                        return Err(invalid("expected a single mzML or indexedmzML root"));
                    }
                    seen_root = true;
                }
                if matches!(
                    parent,
                    "cvParam" | "userParam" | "referenceableParamGroupRef"
                ) {
                    return Err(invalid("parameter/ref elements cannot have children"));
                }
                if (parent == "referenceableParamGroupList" && tag != "referenceableParamGroup")
                    || (parent == "referenceableParamGroup"
                        && !matches!(tag.as_str(), "cvParam" | "userParam"))
                {
                    return Err(invalid(
                        "invalid child in referenceable parameter group/list",
                    ));
                }
                if stack.len() >= 128 {
                    return Err(invalid("XML nesting exceeds 128 levels"));
                }
                if parent == "product" && tag != "isolationWindow" {
                    return Err(invalid("product only permits an isolation window"));
                }
                if (parent == "scanList"
                    && !matches!(
                        tag.as_str(),
                        "cvParam" | "userParam" | "referenceableParamGroupRef" | "scan"
                    ))
                    || (parent == "scan"
                        && !matches!(
                            tag.as_str(),
                            "cvParam"
                                | "userParam"
                                | "referenceableParamGroupRef"
                                | "scanWindowList"
                        ))
                    || (parent == "scanWindowList" && tag != "scanWindow")
                    || (parent == "scanWindow"
                        && !matches!(
                            tag.as_str(),
                            "cvParam" | "userParam" | "referenceableParamGroupRef"
                        ))
                    || (parent == "productList" && tag != "product")
                {
                    return Err(invalid("invalid scan window/product list child"));
                }
                if header_draft.captures(&tag, parent) {
                    if seen_run {
                        return Err(invalid("header element after run start"));
                    }
                    header_draft.start(&tag, attrs, &mut header_work)?;
                    stack.push(tag);
                    buffer.clear();
                    continue;
                }
                if let Some((depth, _, actual)) = group_list.as_mut() {
                    if *depth == stack.len() && tag == "referenceableParamGroup" {
                        *actual += 1;
                    }
                }
                match tag.as_str() {
                    "indexedmzML" if !parent.is_empty() => {
                        return Err(invalid("nested indexedmzML wrapper"));
                    }
                    "isolationWindow" | "activation" => {
                        if parent == "product" && tag == "isolationWindow" {
                            let r = record
                                .as_mut()
                                .filter(|r| r.product_active)
                                .ok_or_else(|| invalid("isolation window outside product"))?;
                            if !r.product_fields.insert("isolation_element") {
                                return Err(invalid("duplicate product isolation window"));
                            }
                        } else {
                            if parent != "precursor" {
                                return Err(invalid("misplaced precursor acquisition element"));
                            }
                            let r = record
                                .as_mut()
                                .filter(|r| r.precursor.is_some())
                                .ok_or_else(|| invalid("acquisition element outside precursor"))?;
                            let key = if tag == "activation" {
                                "activation_element"
                            } else {
                                "isolation_element"
                            };
                            if !r.precursor_fields.insert(key) {
                                return Err(invalid("duplicate precursor acquisition element"));
                            }
                        }
                    }
                    "productList" => {
                        if parent != "spectrum" {
                            return Err(invalid("product list outside spectrum"));
                        }
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("product list outside record"))?;
                        if r.product_list_seen {
                            return Err(invalid("duplicate product list"));
                        }
                        r.product_list_seen = true;
                        // Advisory declared count; still a parameter ceiling.
                        let count = number::<usize>(required(&attrs, "count")?, "product count")?;
                        if count > parameter_budget.remaining {
                            return Err(invalid("product count exceeds parameter limit"));
                        }
                    }
                    "product" => {
                        if !matches!(parent, "productList" | "chromatogram") {
                            return Err(invalid(
                                "product outside spectrum productList/chromatogram",
                            ));
                        }
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("product outside record"))?;
                        if (parent == "productList") != r.spectrum.is_some() {
                            return Err(invalid("product owner mismatch"));
                        }
                        if (r.chromatogram.is_some() && r.product_seen)
                            || r.product_active
                            || r.precursor.is_some()
                        {
                            return Err(invalid("duplicate or nested chromatogram product"));
                        }
                        parameter_budget.begin()?;
                        r.product_seen = true;
                        r.product_active = true;
                        r.product_fields.clear();
                        r.product = Some(Product::default());
                    }
                    "scan" => {
                        if parent != "scanList" {
                            return Err(invalid("scan outside scanList"));
                        }
                        let r = record
                            .as_mut()
                            .filter(|r| r.spectrum.is_some())
                            .ok_or_else(|| invalid("scan outside spectrum"))?;
                        if r.scan_active {
                            return Err(invalid("nested scan"));
                        }
                        r.spectrum
                            .as_mut()
                            .unwrap()
                            .acquisition_info
                            .acquisitions
                            .push(header_registry.scan(
                                &attrs,
                                &mut parameter_budget,
                                &mut header_work,
                            )?);
                        r.scan_active = true;
                        r.scan_window_list_seen = false;
                    }
                    "scanWindowList" => {
                        if parent != "scan" {
                            return Err(invalid("scanWindowList outside scan"));
                        }
                        let r = record
                            .as_mut()
                            .filter(|r| r.scan_active)
                            .ok_or_else(|| invalid("scan window list outside active scan"))?;
                        if r.scan_window_list_seen {
                            return Err(invalid("duplicate scanWindowList"));
                        }
                        r.scan_window_list_seen = true;
                        // Advisory declared count; still a parameter ceiling.
                        let count =
                            number::<usize>(required(&attrs, "count")?, "scan window count")?;
                        if count > parameter_budget.remaining {
                            return Err(invalid("scan window count exceeds parameter limit"));
                        }
                    }
                    "scanWindow" => {
                        if parent != "scanWindowList" {
                            return Err(invalid("scanWindow outside scanWindowList"));
                        }
                        let r = record
                            .as_mut()
                            .filter(|r| r.scan_active)
                            .ok_or_else(|| invalid("scan window outside active scan"))?;
                        parameter_budget.begin()?;
                        if r.scan_window
                            .replace(crate::metadata::ScanWindow::default())
                            .is_some()
                        {
                            return Err(invalid("nested scan window"));
                        }
                        r.scan_window_fields.clear();
                        r.scan_window_unit = None;
                    }
                    "precursorList" | "selectedIonList" | "scanList" => {
                        let expected_parent = match tag.as_str() {
                            "precursorList" | "scanList" => "spectrum",
                            _ => "precursor",
                        };
                        if parent != expected_parent {
                            return Err(invalid("misplaced precursor/scan list"));
                        }
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("list outside spectrum"))?;
                        let seen = match tag.as_str() {
                            "precursorList" => &mut r.precursor_list_seen,
                            "selectedIonList" => &mut r.selected_list_seen,
                            _ => &mut r.scan_list_seen,
                        };
                        if *seen {
                            return Err(invalid("duplicate precursor/scan list"));
                        }
                        *seen = true;
                        // Advisory declared count; still a parameter ceiling for scans.
                        let declared: usize = number(required(&attrs, "count")?, "list count")?;
                        if tag == "scanList" && declared > parameter_budget.remaining {
                            return Err(invalid("scan count exceeds parameter limit"));
                        }
                    }
                    "mzML" => {
                        if seen_mzml || !matches!(parent, "" | "indexedmzML") {
                            return Err(invalid("unexpected mzML element"));
                        }
                        if !required(&attrs, "version")?.starts_with("1.1.") {
                            return Err(Error::Unsupported("only mzML 1.1 is supported".into()));
                        }
                        if let Some(accession) = attrs.get("accession") {
                            experiment.settings.document.identifier =
                                header_work.copy(accession)?;
                        }
                        if let Some(id) = attrs.get("id") {
                            header_work
                                .meter()
                                .tree::<(String, crate::metadata::MetaValue)>(1)?;
                            experiment
                                .settings
                                .metadata
                                .insert("mzml_id".into(), header_work.copy(id)?.into());
                        }
                        seen_mzml = true;
                    }
                    "fileDescription"
                    | "sourceFileList"
                    | "sourceFile"
                    | "contact"
                    | "fileContent"
                    | "sampleList"
                    | "sample"
                    | "softwareList"
                    | "software"
                    | "instrumentConfigurationList"
                    | "instrumentConfiguration"
                    | "componentList"
                    | "source"
                    | "analyzer"
                    | "detector"
                    | "softwareRef"
                    | "dataProcessingList"
                    | "dataProcessing"
                    | "processingMethod" => {
                        return Err(invalid("misplaced mzML header element"));
                    }
                    "run" => {
                        if parent != "mzML" || seen_run {
                            return Err(invalid("expected exactly one mzML run"));
                        }
                        header_registry = std::mem::take(&mut header_draft).finish(
                            &attrs,
                            &groups,
                            &mut parameter_budget,
                            &mut header_work,
                            &mut experiment.settings,
                            options,
                        )?;
                        seen_run = true;
                    }
                    "referenceableParamGroupList" => {
                        if parent != "mzML" || group_list_seen || seen_run {
                            return Err(invalid("misplaced/duplicate parameter group list"));
                        }
                        let expected =
                            number::<usize>(required(&attrs, "count")?, "parameter group count")?;
                        if expected == 0 || expected > remaining_groups {
                            return Err(invalid(
                                "parameter group count is zero or exceeds configured limit",
                            ));
                        }
                        group_list_seen = true;
                        group_list = Some((stack.len() + 1, expected, 0));
                    }
                    "referenceableParamGroup" => {
                        if parent != "referenceableParamGroupList" || group.is_some() {
                            return Err(invalid("misplaced parameter group"));
                        }
                        let id = parameter_id(required(&attrs, "id")?)?;
                        if groups.contains_key(id) {
                            return Err(invalid("duplicate parameter group ID"));
                        }
                        remaining_groups = remaining_groups.checked_sub(1).ok_or_else(|| {
                            invalid("parameter group count exceeds configured limit")
                        })?;
                        group = Some((id.to_owned(), Vec::new()));
                    }
                    "referenceableParamGroupRef" => {
                        if !parameter_context(parent) {
                            return Err(invalid("misplaced parameter group reference"));
                        }
                        let id = parameter_id(required(&attrs, "ref")?)?;
                        if let Some(parameters) = groups.get(id) {
                            apply_group(
                                parameters,
                                parent,
                                &mut parameter_budget,
                                &mut record,
                                &mut binary,
                                &mut experiment,
                                selection.as_mut(),
                            )?;
                        } else {
                            return Err(invalid(format!("unknown parameter group {id}")));
                        }
                    }
                    "spectrumList" | "chromatogramList" => {
                        if parent != "run" {
                            return Err(invalid("record list outside run"));
                        }
                        if metadata_only {
                            let id = required(&attrs, "defaultDataProcessingRef")?;
                            header_registry.processing(id, &mut header_work)?;
                            return Ok(experiment);
                        }
                        let slot = if tag == "spectrumList" {
                            &mut spectrum_list_seen
                        } else {
                            &mut chromatogram_list_seen
                        };
                        if *slot {
                            return Err(invalid("duplicate record list"));
                        }
                        *slot = true;
                        // Advisory declared count; records are bounded by
                        // `options.max_records` and `scaling.records` as each
                        // one opens.
                        let _declared: usize = number(required(&attrs, "count")?, "record count")?;
                        default_processing = attrs
                            .get("defaultDataProcessingRef")
                            .map(|id| header_registry.processing(id, &mut header_work))
                            .transpose()?
                            .unwrap_or_default();
                    }
                    "spectrum" | "chromatogram" => {
                        let expected_parent = if tag == "spectrum" {
                            "spectrumList"
                        } else {
                            "chromatogramList"
                        };
                        if parent != expected_parent || record.is_some() {
                            return Err(invalid("misplaced spectrum/chromatogram"));
                        }
                        remaining_records = remaining_records
                            .checked_sub(1)
                            .ok_or_else(|| invalid("record count exceeds configured limit"))?;
                        let count: usize = number(
                            required(&attrs, "defaultArrayLength")?,
                            "defaultArrayLength",
                        )?;
                        if fill_data {
                            remaining_peaks = remaining_peaks
                                .checked_sub(count)
                                .ok_or_else(|| invalid("peak count exceeds configured limit"))?;
                        }
                        let id = required(&attrs, "id")?.to_owned();
                        if id.is_empty() || !ids.insert((tag.clone(), id.clone())) {
                            return Err(invalid("empty or duplicate record id"));
                        }
                        record = Some(Record {
                            spectrum: (tag == "spectrum").then(|| MSSpectrum {
                                native_id: id.clone(),
                                ..Default::default()
                            }),
                            chromatogram: (tag == "chromatogram").then(|| MSChromatogram {
                                native_id: id,
                                ..Default::default()
                            }),
                            count,
                            array_names: BTreeSet::new(),
                            primary_seen: [false; 2],
                            coordinates: None,
                            intensities: None,
                            precursor: None,
                            selected_ion: false,
                            ignore_selected_ion: false,
                            precursor_mz_selected_ion: load
                                .is_none_or(|o| o.scientific.precursor_mz_selected_ion),
                            selected_mz: false,
                            rt_seen: false,
                            precursors_seen: 0,
                            name_seen: false,
                            seen_fields: BTreeSet::new(),
                            precursor_fields: BTreeSet::new(),
                            precursor_list_seen: false,
                            selected_list_seen: false,
                            scan_list_seen: false,
                            product_active: false,
                            product_seen: false,
                            product_fields: BTreeSet::new(),
                            product_list_seen: false,
                            product: None,
                            scan_active: false,
                            scan_window_list_seen: false,
                            scan_window: None,
                            scan_window_fields: BTreeSet::new(),
                            scan_window_unit: None,
                            precursor_filter: load.and_then(|o| {
                                o.scientific
                                    .has_precursor_mz_range()
                                    .then(|| o.scientific.precursor_mz_range())
                            }),
                            precursor_outside: false,
                            raw_float_arrays: Vec::new(),
                            nonfinite_float_arrays: options.source_nonfinite_float_arrays,
                            primary_metadata: Default::default(),
                            wavelength: None,
                            intensity_first: false,
                        });
                        let processing = if let Some(id) = attrs.get("dataProcessingRef") {
                            header_registry.processing(id, &mut header_work)?
                        } else {
                            header_work.slots::<std::sync::Arc<crate::metadata::DataProcessing>>(
                                default_processing.len(),
                            )?;
                            default_processing.clone()
                        };
                        if tag == "chromatogram" && attrs.contains_key("sourceFileRef") {
                            return Err(Error::Unsupported(
                                "chromatogram sourceFileRef is not permitted by mzML".into(),
                            ));
                        }
                        let source = attrs
                            .get("sourceFileRef")
                            .map(|id| header_registry.source(id, &mut header_work))
                            .transpose()?
                            .unwrap_or_default();
                        let r = record.as_mut().unwrap();
                        if let Some(s) = &mut r.spectrum {
                            s.data_processing = processing;
                            s.source_file = source;
                            if let Some(spot) = attrs.get("spotID") {
                                if !spot.is_empty() {
                                    header_work
                                        .meter()
                                        .meta_update(&s.metadata, "maldi_spot_id")?;
                                    s.metadata.insert(
                                        "maldi_spot_id".into(),
                                        header_work.copy(spot)?.into(),
                                    );
                                }
                            }
                        } else {
                            let c = r.chromatogram.as_mut().unwrap();
                            c.data_processing = processing;
                            c.source_file = source;
                        }
                        array_list_seen = false;
                    }
                    "binaryDataArrayList" => {
                        if !matches!(parent, "spectrum" | "chromatogram")
                            || record.is_none()
                            || array_list_seen
                        {
                            return Err(invalid("misplaced/duplicate binary array list"));
                        }
                        array_list_seen = true;
                        // Advisory declared count; arrays are bounded by
                        // `options.max_total_arrays` and `scaling.arrays` as
                        // each one opens.
                        let _declared: usize =
                            number(required(&attrs, "count")?, "binary array count")?;
                    }
                    "binaryDataArray" => {
                        if parent != "binaryDataArrayList" || binary.is_some() {
                            return Err(invalid("misplaced binary array"));
                        }
                        let encoded_length: usize =
                            number(required(&attrs, "encodedLength")?, "encodedLength")?;
                        let max_encoded = options
                            .max_array_bytes
                            .saturating_add(2)
                            .saturating_div(3)
                            .saturating_mul(4);
                        if fill_data && encoded_length > max_encoded {
                            return Err(invalid("encodedLength exceeds configured byte limit"));
                        }
                        binary = Some(Binary {
                            spectrum: record.as_ref().is_some_and(|r| r.spectrum.is_some()),
                            encoded_length,
                            encoded: buffers.encoded(),
                            data_processing: attrs
                                .get("dataProcessingRef")
                                .map(|id| header_registry.processing(id, &mut header_work))
                                .transpose()?
                                .unwrap_or_default(),
                            array_length: attrs
                                .get("arrayLength")
                                .map(|n| number(n, "arrayLength"))
                                .transpose()?,
                            ..Default::default()
                        });
                        remaining_arrays = remaining_arrays
                            .checked_sub(1)
                            .ok_or_else(|| invalid("total binary array count limit exceeded"))?;
                    }
                    "binary" => {
                        if parent != "binaryDataArray" {
                            return Err(invalid("misplaced binary element"));
                        }
                        let b = binary
                            .as_mut()
                            .ok_or_else(|| invalid("binary element outside array"))?;
                        if b.has_binary {
                            return Err(invalid("duplicate binary element"));
                        }
                        b.has_binary = true;
                    }
                    "precursor" => {
                        if !matches!(parent, "precursorList" | "chromatogram") {
                            return Err(invalid("misplaced precursor"));
                        }
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("precursor outside record"))?;
                        r.precursors_seen += 1;
                        if r.chromatogram.is_some() && r.precursors_seen > 1 {
                            return Err(invalid("multiple chromatogram precursors"));
                        }
                        if r.precursor.replace(Precursor::default()).is_some() {
                            return Err(invalid("nested precursor"));
                        }
                        let p = r.precursor.as_mut().unwrap();
                        p.spectrum_reference = attrs
                            .get("spectrumRef")
                            .map(|s| header_work.copy(s))
                            .transpose()?;
                        p.cv_terms.metadata = header_registry.source_metadata(
                            &attrs,
                            &mut parameter_budget,
                            &mut header_work,
                        )?;
                        if let Some(id) = attrs.get("externalSpectrumID") {
                            header_work
                                .meter()
                                .tree::<(String, crate::metadata::MetaValue)>(1)?;
                            header_work.charge(20, 20)?;
                            p.cv_terms.metadata.insert(
                                "external_spectrum_id".into(),
                                header_work.copy(id)?.into(),
                            );
                        }
                        r.selected_ion = false;
                        r.ignore_selected_ion = false;
                        r.selected_mz = false;
                        r.selected_list_seen = false;
                        r.precursor_fields.clear();
                    }
                    "selectedIon" => {
                        if parent != "selectedIonList" {
                            return Err(invalid("misplaced selected ion"));
                        }
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("selected ion outside record"))?;
                        if r.precursor.is_none() {
                            return Err(invalid("selected ion outside precursor"));
                        }
                        if r.selected_ion && r.precursor_mz_selected_ion {
                            return Err(Error::Unsupported(
                                "multiple selected ions in one precursor".into(),
                            ));
                        }
                        r.ignore_selected_ion = r.selected_ion;
                        r.selected_ion = true;
                    }
                    "cvParam" | "userParam" => {
                        required(
                            &attrs,
                            if tag == "cvParam" {
                                "accession"
                            } else {
                                "name"
                            },
                        )?;
                        if parent == "referenceableParamGroup" {
                            let (_, parameters) = group.as_mut().unwrap();
                            if tag == "cvParam"
                                && parameters.last().is_some_and(|p| p.tag == "userParam")
                            {
                                return Err(invalid(
                                    "cvParam follows userParam in parameter group",
                                ));
                            }
                            parameters.push(Parameter {
                                tag: if tag == "cvParam" {
                                    "cvParam"
                                } else {
                                    "userParam"
                                },
                                attrs,
                            });
                        } else {
                            apply_parameter(
                                &tag,
                                &attrs,
                                parent,
                                &mut record,
                                &mut binary,
                                &mut experiment,
                                selection.as_mut(),
                                &mut parameter_budget,
                            )?;
                        }
                    }
                    _ => {}
                }
                stack.push(tag);
            }
            Event::End(_) => {
                if let Some((depth, _declared, _actual)) = group_list {
                    if depth == stack.len() {
                        // The declared count stays a resource ceiling above, but a
                        // mismatch with the actual number of children is advisory on
                        // reading: upstream has a handler for `referenceableParamGroup`
                        // (MzMLHandler.cpp:1055) and none for the enclosing list, so it
                        // never reads this count at all. Writing still emits the true
                        // count. Same rule as the record, binary-array and
                        // precursor/scan list counts.
                        group_list = None;
                    }
                }
                let tag = stack
                    .pop()
                    .ok_or_else(|| invalid("unmatched closing tag"))?;
                if header_draft.end(&tag)? {
                    buffer.clear();
                    continue;
                }
                match tag.as_str() {
                    "referenceableParamGroup" => {
                        let (id, parameters) = group
                            .take()
                            .ok_or_else(|| invalid("missing parameter group state"))?;
                        groups.insert(id, parameters);
                    }
                    "binaryDataArray" => {
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("array outside record"))?;
                        let mut b = binary
                            .take()
                            .ok_or_else(|| invalid("missing binary array state"))?;
                        if !fill_data {
                            let (_, _, count) = b.descriptor(r.count)?;
                            r.check_array_kind(b.kind.as_ref().unwrap())?;
                            if matches!(
                                b.kind,
                                Some(Kind::Role(record_transport::Role::Wavelength))
                            ) {
                                r.wavelength = Some((0, count));
                            }
                            // Source creates no primary/auxiliary data or descriptions
                            // when data population is disabled. XML descriptors were
                            // parsed, but their encoded payload is never decoded.
                            buffer.clear();
                            continue;
                        }
                        let metadata = std::mem::take(&mut b.metadata);
                        let data_processing = std::mem::take(&mut b.data_processing);
                        // The Numpress coder's work and allocation allowances
                        // are per array: its defaults plus a multiple of this
                        // array's encoded text and declared values. One
                        // allowance for the whole document rejected any
                        // Numpress file beyond a few megabytes; the cumulative
                        // element and byte counters charged in `decode` bound
                        // the total.
                        let mut numpress_work = {
                            let mut limits = numpress_limits;
                            let text = b.encoded.len();
                            let values = b.array_length.unwrap_or(r.count);
                            limits.raw.max_work = limits
                                .raw
                                .max_work
                                .saturating_add(text.saturating_mul(512))
                                .saturating_add(values.saturating_mul(64));
                            limits.max_total_bytes = limits
                                .max_total_bytes
                                .saturating_add(text.saturating_mul(64))
                                .saturating_add(values.saturating_mul(32));
                            coder::Work::new(limits)
                        };
                        let decoded = b.decode(
                            r.count,
                            options,
                            &mut remaining_array_bytes,
                            &mut remaining_array_elements,
                            &mut numpress_work,
                            &mut buffers,
                        );
                        buffers.recycle(std::mem::take(&mut b.encoded));
                        let (kind, values) = decoded?;
                        r.check_array_kind(&kind)?;
                        if let Kind::Role(role) = kind {
                            if role.noise() {
                                if !metadata.is_empty() || !data_processing.is_empty() {
                                    return Err(Error::Unsupported("independent noise array descriptions/history have no owner".into()));
                                }
                                let Values::Floats(values) = values else {
                                    return Err(invalid("noise array must be floating point"));
                                };
                                let key = role.terms().1;
                                header_work.meter().meta_update(r.metadata(), key)?;
                                header_work.charge(values.len(), 0)?;
                                if r.metadata()
                                    .insert(key.into(), MetaValue::try_from(values)?)
                                    .is_some()
                                {
                                    return Err(invalid("duplicate noise metadata owner"));
                                }
                                buffer.clear();
                                continue;
                            }
                            if role == record_transport::Role::Wavelength {
                                let Values::Floats(values) = values else {
                                    return Err(invalid("wavelength array must be floating point"));
                                };
                                r.wavelength = Some((r.raw_float_arrays.len(), values.len()));
                                r.raw_float_arrays.push(DataArray {
                                    name: "wavelength array".into(),
                                    data: values,
                                    metadata,
                                    data_processing,
                                });
                                buffer.clear();
                                continue;
                            }
                        }
                        if let Kind::Auxiliary(name) = kind {
                            let (_floats, integers, strings) =
                                if let Some(spectrum) = &mut r.spectrum {
                                    (
                                        &mut spectrum.float_data_arrays,
                                        &mut spectrum.integer_data_arrays,
                                        &mut spectrum.string_data_arrays,
                                    )
                                } else {
                                    let chromatogram = r.chromatogram.as_mut().unwrap();
                                    (
                                        &mut chromatogram.float_data_arrays,
                                        &mut chromatogram.integer_data_arrays,
                                        &mut chromatogram.string_data_arrays,
                                    )
                                };
                            match values {
                                Values::Floats(values) => {
                                    r.raw_float_arrays.push(DataArray {
                                        name,
                                        data: values,
                                        metadata,
                                        data_processing,
                                    });
                                }
                                Values::Integers(values) => integers.push(DataArray {
                                    name,
                                    data: values,
                                    metadata,
                                    data_processing,
                                }),
                                Values::Strings(values) => strings.push(DataArray {
                                    name,
                                    data: values,
                                    metadata,
                                    data_processing,
                                }),
                            }
                        } else {
                            if !data_processing.is_empty() {
                                return Err(Error::Unsupported(
                                    "primary-array processing has no independent native owner"
                                        .into(),
                                ));
                            }
                            // Source spectra merge m/z metadata before intensity, regardless
                            // of XML order; chromatograms merge in encounter order.
                            let index =
                                usize::from(matches!(kind, Kind::Intensity | Kind::Role(_)));
                            if r.coordinates.is_none() && r.intensities.is_none() {
                                r.intensity_first = index == 1;
                            }
                            r.primary_metadata[index] = metadata;
                            let Values::Floats(values) = values else {
                                return Err(invalid("non-floating primary binary array"));
                            };
                            let slot = match kind {
                                Kind::Intensity | Kind::Role(_) => &mut r.intensities,
                                Kind::Mz if r.spectrum.is_some() => &mut r.coordinates,
                                Kind::Time if r.chromatogram.is_some() => &mut r.coordinates,
                                _ => {
                                    return Err(invalid(
                                        "coordinate array has wrong type for record",
                                    ));
                                }
                            };
                            if slot.replace(values).is_some() {
                                return Err(invalid("duplicate coordinate/intensity array"));
                            }
                        }
                    }
                    "precursor" => {
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("precursor outside record"))?;
                        let mut p = r
                            .precursor
                            .take()
                            .ok_or_else(|| invalid("missing precursor state"))?;
                        if p.isolation_target_mz == Some(p.mz) {
                            p.isolation_target_mz = None;
                        }
                        p.validate()?;
                        if let Some(s) = &mut r.spectrum {
                            s.precursors.push(p);
                        } else {
                            r.chromatogram.as_mut().unwrap().precursor = p;
                        }
                    }
                    "product" => {
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("product outside record"))?;
                        if !r.product_active {
                            return Err(invalid("missing product state"));
                        }
                        r.product_active = false;
                        let product = r.product.take().ok_or_else(|| invalid("missing product"))?;
                        product.validate()?;
                        if let Some(spectrum) = &mut r.spectrum {
                            spectrum.products.push(product);
                        } else {
                            r.chromatogram.as_mut().unwrap().product = product;
                        }
                    }
                    "scanWindow" => {
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("scan window outside record"))?;
                        let window = r
                            .scan_window
                            .take()
                            .ok_or_else(|| invalid("missing scan window"))?;
                        window.validate()?;
                        r.spectrum
                            .as_mut()
                            .ok_or_else(|| invalid("scan window outside spectrum"))?
                            .instrument_settings
                            .scan_windows
                            .push(window);
                    }
                    "scan" => {
                        record
                            .as_mut()
                            .ok_or_else(|| invalid("scan outside record"))?
                            .scan_active = false;
                    }
                    "spectrum" | "chromatogram" => {
                        let mut record = record
                            .take()
                            .ok_or_else(|| invalid("missing record state"))?;
                        if let Some(s) = &mut record.spectrum {
                            acquisition_metadata::normalize(
                                &mut s.acquisition_info,
                                options.acquisition_mode,
                            );
                            // Source fallback precedes primary-array metadata merging
                            // and does not manufacture a scan-time filtering event.
                            if !record.rt_seen {
                                let height = (usize::BITS - s.metadata.len().max(1).leading_zeros())
                                    as usize;
                                header_work.charge(height.saturating_mul(16 * 23), 0)?;
                                if let Some(value) = s.metadata.get("elution time (seconds)") {
                                    s.rt = value.as_f64()?;
                                }
                            }
                        }
                        if let Some((index, count)) = record.wavelength.take() {
                            if !record.primary_seen[0] && count != record.count {
                                return Err(invalid(
                                    "primary wavelength length differs from defaultArrayLength",
                                ));
                            }
                            if fill_data && record.coordinates.is_none() {
                                if !record.raw_float_arrays[index].data_processing.is_empty() {
                                    return Err(Error::Unsupported(
                                        "wavelength processing has no primary owner".into(),
                                    ));
                                }
                                header_work.charge(record.raw_float_arrays.len(), 0)?;
                                let array = record.raw_float_arrays.remove(index);
                                record.coordinates = Some(array.data);
                                record.primary_metadata[0] = array.metadata;
                            } else if fill_data {
                                record.raw_float_arrays[index]
                                    .metadata
                                    .remove(record_transport::COORDINATE);
                            }
                        }
                        let order = if record.chromatogram.is_some() && record.intensity_first {
                            [1, 0]
                        } else {
                            [0, 1]
                        };
                        for index in order {
                            let metadata = std::mem::take(&mut record.primary_metadata[index]);
                            for (name, value) in metadata {
                                if record.spectrum.is_some()
                                    && record_transport::NOISE.contains(&name.as_str())
                                {
                                    return Err(invalid(
                                        "primary metadata collides with independent noise ownership",
                                    ));
                                }
                                header_work.meter().meta_update(record.metadata(), &name)?;
                                record.metadata().insert(name, value);
                            }
                        }
                        if let Some(completed) = record.finish(selection.as_mut(), fill_data)? {
                            if let Some(sink) = consumer.as_deref_mut() {
                                if !sink.accept(completed)? {
                                    return Ok(experiment);
                                }
                            } else {
                                match completed {
                                    consumer::Completed::Spectrum(s) => experiment.spectra.push(s),
                                    consumer::Completed::Chromatogram(c) => {
                                        experiment.chromatograms.push(c)
                                    }
                                }
                            }
                        }
                    }
                    "mzML" => {
                        if let Some(sink) = consumer.as_deref_mut() {
                            if !sink.finish()? {
                                return Ok(experiment);
                            }
                        }
                    }
                    _ => {}
                }
            }
            Event::Text(text) => {
                let text = text.decode().map_err(|e| invalid(e.to_string()))?;
                xml_string(&text)?;
                if stack.last().is_some_and(|tag| tag == "binary") {
                    binary
                        .as_mut()
                        .ok_or_else(|| invalid("text outside binary array"))?
                        .append_encoded(text.as_bytes(), normalize_binary, fill_data)?;
                } else if !seen_run && !text.trim().is_empty() {
                    return Err(invalid("non-whitespace text in mzML header"));
                } else if stack.last().is_some_and(|tag| {
                    matches!(
                        tag.as_str(),
                        "referenceableParamGroupList"
                            | "referenceableParamGroup"
                            | "referenceableParamGroupRef"
                            | "cvParam"
                            | "userParam"
                    )
                }) && !text.trim().is_empty()
                {
                    return Err(invalid("text in parameter group/list/parameter/ref"));
                } else if stack.is_empty() && !text.trim().is_empty() {
                    return Err(invalid("text outside XML root"));
                }
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
                    latin1_subset = encoding.eq_ignore_ascii_case(b"ISO-8859-1");
                    ascii_only = encoding.eq_ignore_ascii_case(b"US-ASCII") || latin1_subset;
                    if !encoding.eq_ignore_ascii_case(b"UTF-8") && !ascii_only {
                        return Err(Error::Unsupported(
                            "only UTF-8 or ASCII-compatible US-ASCII/ISO-8859-1 mzML XML is supported".into(),
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
            Event::GeneralRef(_) | Event::CData(_) => {
                return Err(Error::Unsupported(
                    "XML entity references in text and CDATA are not supported".into(),
                ));
            }
            Event::Eof => break,
            _ => {}
        }
        buffer.clear();
    }
    if !seen_mzml || !seen_run || !stack.is_empty() || record.is_some() {
        return Err(invalid("incomplete mzML document"));
    }
    Ok(experiment)
}

fn escape(value: &str) -> String {
    quick_xml::escape::escape(value)
        .replace('\n', "&#10;")
        .replace('\r', "&#13;")
        .replace('\t', "&#9;")
}
// Length of the leading run of `bytes` inside `first ..= first + span`,
// wrapping, so a run is expressed as one subtraction and one comparison.
//
// Text nodes carry the base64 payload of a profile spectrum, tens of megabytes
// per record, so the scan runs a fixed-size block at a time: the block reduces
// to a branchless `or` of saturating subtractions, which the optimiser turns
// into one `psubb`/`psubusb`/`por` chain per vector instead of one compare and
// one branch per byte. Only a block that holds something outside the run is
// walked byte by byte.
fn run_inside(bytes: &[u8], first: u8, span: u8) -> usize {
    const BLOCK: usize = 64;
    let mut index = 0;
    while let Some(block) = bytes.get(index..).and_then(<[u8]>::first_chunk::<BLOCK>) {
        let mut outside = 0u8;
        for &byte in block {
            outside |= byte.wrapping_sub(first).saturating_sub(span);
        }
        if outside != 0 {
            break;
        }
        index += BLOCK;
    }
    while index < bytes.len() && bytes[index].wrapping_sub(first) <= span {
        index += 1;
    }
    index
}
// The XML 1.0 `Char` production, tested on the UTF-8 bytes rather than on
// decoded scalars: `value` is already valid UTF-8, so the only excluded
// scalars it can hold are the C0 controls other than tab/newline/return and
// the two noncharacters U+FFFE and U+FFFF. Surrogates cannot appear in a
// `str`, and every scalar at or above U+10000 is admitted, so no other range
// needs testing. `0xef` only ever starts a three-byte sequence in valid UTF-8,
// which makes the two noncharacters pure byte patterns.
fn xml_string(value: &str) -> Result<()> {
    let bytes = value.as_bytes();
    let mut index = 0;
    while index < bytes.len() {
        index += run_inside(&bytes[index..], 0x20, 0x5f);
        let Some(&byte) = bytes.get(index) else { break };
        if matches!(byte, b'\t' | b'\n' | b'\r') {
            index += 1;
            continue;
        }
        let rest = &bytes[index..];
        if byte < 0x80 || rest.starts_with(b"\xef\xbf\xbe") || rest.starts_with(b"\xef\xbf\xbf") {
            return Err(Error::InvalidValue(
                "text contains invalid XML 1.0 characters".into(),
            ));
        }
        index += 1;
    }
    Ok(())
}

// XMLHandler::fromXSDString retains scalar numeric types, and treats every
// other XSD type as text. Lists and Empty have no reversible source encoding.
fn product_user_value(attrs: &BTreeMap<String, String>) -> Result<MetaValue> {
    scalar_user_value(attrs, attrs.get("type").map(String::as_str).unwrap_or(""))
}
fn scalar_user_value(attrs: &BTreeMap<String, String>, kind: &str) -> Result<MetaValue> {
    use crate::data_structures::list::ListParse;
    let text = attrs.get("value").map(String::as_str).unwrap_or("");
    let data = match kind {
        "xsd:double" | "xsd:float" | "xsd:decimal" => {
            MetaValueData::Float(f64::from_list_item(text)?)
        }
        "xsd:byte" | "xsd:int" | "xsd:unsignedShort" | "xsd:short" | "xsd:unsignedByte"
        | "xsd:unsignedInt" => MetaValueData::Integer(i64::from(i32::from_list_item(text)?)),
        "xsd:long"
        | "xsd:unsignedLong"
        | "xsd:integer"
        | "xsd:negativeInteger"
        | "xsd:nonNegativeInteger"
        | "xsd:nonPositiveInteger"
        | "xsd:positiveInteger" => {
            let token = text.trim_matches([' ', '\t', '\n', '\r']);
            let token = token.strip_prefix('+').unwrap_or(token);
            if token.starts_with('+') {
                return Err(invalid("invalid product metadata integer"));
            }
            MetaValueData::Integer(number(token, "product metadata integer")?)
        }
        _ => MetaValueData::String(text.into()),
    };
    let mut result = MetaValue::new(data)?;
    if let Some(accession) = attrs.get("unitAccession") {
        let cv_ref = attrs
            .get("unitCvRef")
            .map(String::as_str)
            .unwrap_or_else(|| accession.split_once(':').map_or("", |(prefix, _)| prefix));
        let unit = Unit::new(
            accession,
            attrs.get("unitName").map(String::as_str).unwrap_or(""),
            cv_ref,
        )?;
        product_unit(&unit)?;
        result = result.with_unit(unit)?;
    } else if attrs.contains_key("unitCvRef") || attrs.contains_key("unitName") {
        return Err(invalid("product metadata unit has no accession"));
    }
    Ok(result)
}

fn product_unit(unit: &Unit) -> Result<()> {
    if !matches!(unit.cv_ref(), "MS" | "UO")
        || unit.accession().split_once(':').map(|(prefix, _)| prefix) != Some(unit.cv_ref())
    {
        return Err(Error::Unsupported(
            "product metadata units require MS or UO identity".into(),
        ));
    }
    for text in [unit.accession(), unit.name(), unit.cv_ref()] {
        xml_string(text)?;
    }
    Ok(())
}

fn validate_product_write(product: &Product) -> Result<()> {
    product.validate()?;
    // Source MzMLHandler does not preserve arbitrary Product CVTermList entries
    // with their typed identity through this isolation-window representation.
    if !product.cv_terms.is_empty() {
        return Err(Error::Unsupported(
            "arbitrary product CV terms are not represented in mzML".into(),
        ));
    }
    validate_scalar_metadata(&product.cv_terms.metadata)
}
/// Writer preflight for the spectrum-level mobility the scan writer emits.
///
/// The source writes a FAIMS voltage whenever the unit is FAIMS, including the
/// unset `-1` sentinel, and any other unit only for a set drift time; a drift
/// time without a unit is written as milliseconds with a warning
/// (`MzMLHandler.cpp:5412-5440`). This port writes the same two lossless cases
/// and refuses the lossy ones: a drift time without a unit, and a non-FAIMS
/// unit without a drift time, which the source drops. A non-finite drift time
/// could not be read back and is refused too.
fn validate_spectrum_mobility(s: &MSSpectrum) -> Result<()> {
    if !s.drift_time.is_finite() {
        return Err(Error::InvalidValue(
            "spectrum drift time must be finite".into(),
        ));
    }
    let faims = s.drift_time_unit == DriftTimeUnit::FaimsCompensationVoltage;
    if (s.has_drift_time() && s.drift_time_unit == DriftTimeUnit::None)
        || (!s.has_drift_time() && !faims && s.drift_time_unit != DriftTimeUnit::None)
    {
        return Err(Error::Unsupported(
            "spectrum mobility requires both a value and an explicit unit".into(),
        ));
    }
    Ok(())
}
fn validate_scalar_metadata(metadata: &crate::metadata::MetaInfo) -> Result<()> {
    for (key, value) in metadata {
        validate_scalar_value(key, value)?;
    }
    Ok(())
}
fn validate_scalar_value(key: &str, value: &MetaValue) -> Result<()> {
    xml_string(key)?;
    match value.data() {
        MetaValueData::String(text) => xml_string(text)?,
        MetaValueData::Integer(_) | MetaValueData::Float(_) => {}
        _ => {
            return Err(Error::Unsupported(
                "Empty/list metadata has no lossless mzML scalar encoding".into(),
            ));
        }
    }
    if let Some(unit) = value.unit() {
        product_unit(unit)?;
    }
    Ok(())
}

fn write_product(w: &mut impl Write, product: &Product) -> Result<()> {
    const MZ_UNIT: &str = " unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"";
    writeln!(w, "<product><isolationWindow>")?;
    cv(
        w,
        "MS:1000827",
        "isolation window target m/z",
        &product.mz.to_string(),
        MZ_UNIT,
    )?;
    for (quantity, accession, name) in [
        (
            product.isolation_window_lower_offset,
            "MS:1000828",
            "isolation window lower offset",
        ),
        (
            product.isolation_window_upper_offset,
            "MS:1000829",
            "isolation window upper offset",
        ),
    ] {
        if quantity > 0.0 || quantity.is_sign_negative() {
            cv(w, accession, name, &quantity.to_string(), MZ_UNIT)?;
        }
    }
    write_scalar_metadata(w, &product.cv_terms.metadata, None)?;
    writeln!(w, "</isolationWindow></product>")?;
    Ok(())
}
fn write_scalar_metadata(
    w: &mut impl Write,
    metadata: &crate::metadata::MetaInfo,
    skip: Option<&str>,
) -> Result<()> {
    write_scalar_metadata_skipping(w, metadata, skip.as_slice())
}
fn write_scalar_metadata_skipping(
    w: &mut impl Write,
    metadata: &crate::metadata::MetaInfo,
    skip: &[&str],
) -> Result<()> {
    for (name, value) in metadata {
        if skip.contains(&name.as_str()) {
            continue;
        }
        let (kind, text) = match value.data() {
            MetaValueData::String(text) => {
                ("xsd:string", std::borrow::Cow::Borrowed(text.as_str()))
            }
            MetaValueData::Integer(n) => ("xsd:integer", std::borrow::Cow::Owned(n.to_string())),
            MetaValueData::Float(n) => ("xsd:double", std::borrow::Cow::Owned(float_text(*n))),
            _ => unreachable!("product preflight checked scalar metadata"),
        };
        write!(
            w,
            "<userParam name=\"{}\" type=\"{kind}\" value=\"{}\"",
            escape(name),
            escape(&text)
        )?;
        if let Some(unit) = value.unit() {
            write!(
                w,
                " unitAccession=\"{}\" unitCvRef=\"{}\" unitName=\"{}\"",
                escape(unit.accession()),
                escape(unit.cv_ref()),
                escape(unit.name())
            )?;
        }
        writeln!(w, "/>")?;
    }
    Ok(())
}

/// Write one `cvParam`, omitting `value` when it is empty.
///
/// Source `MzMLHandler::writeCV_` (3600-3606) writes the attribute only for a
/// non-empty `DataValue`, and its literal valueless terms carry no `value`
/// either; a reader cannot distinguish an absent value from an empty one, so
/// nothing is lost.
fn cv(w: &mut impl Write, accession: &str, name: &str, value: &str, unit: &str) -> Result<()> {
    write!(
        w,
        "<cvParam cvRef=\"MS\" accession=\"{accession}\" name=\"{name}\""
    )?;
    if !value.is_empty() {
        write!(w, " value=\"{}\"", escape(value))?;
    }
    writeln!(w, "{unit}/>")?;
    Ok(())
}
/// C++ `StringUtils::toStr(double)` text, when it reads back as the same value.
///
/// Every number the source writes into a `cvParam` or `userParam` goes through
/// `DataValue::toString` and thus `NumericFormatting::appendNumeric`
/// (`StringUtils.cpp:384`): 15 fraction digits for magnitudes in `[1e-2, 1e4)`
/// and zero, the shortest round-tripping scientific text otherwise. An
/// inherited `3.0`, `1.0e20` or `2.027586375e06` is therefore written back
/// unchanged instead of being reformatted. The fixed branch keeps only 15
/// fraction digits, which loses precision for some values; those keep Rust's
/// shortest round-tripping text, because this port does not discard data it
/// was given.
fn float_text(value: f64) -> String {
    let text = crate::format::file_info::text_format::to_str(value);
    let exact = crate::data_structures::list::ListParse::from_list_item(&text)
        .is_ok_and(|parsed: f64| parsed.to_bits() == value.to_bits());
    if exact { text } else { value.to_string() }
}
const SECOND: &str = " unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"";
/// The intensity array's source unit (`MzMLHandler.cpp:5688`).
const COUNTS: &str =
    " unitCvRef=\"MS\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\"";
fn write_precursor(w: &mut impl Write, precursor: &Precursor, tpp: bool) -> Result<()> {
    precursor_metadata::write_start(w, precursor, tpp)?;
    cv(
        w,
        "MS:1000744",
        "selected ion m/z",
        &precursor_metadata::selected_mz(precursor)?.to_string(),
        " unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"",
    )?;
    if tpp || precursor.charge != 0 {
        cv(
            w,
            "MS:1000041",
            "charge state",
            &precursor.charge.to_string(),
            "",
        )?;
    }
    precursor_metadata::write_intensity(w, precursor)?;
    precursor_metadata::write_end(w, precursor)
}
// Keep the established scalar codec arguments plus its precomputed header.
#[allow(clippy::too_many_arguments)]
fn write_array(
    w: &mut impl Write,
    bytes: impl FnOnce() -> Vec<u8>,
    kind: Kind,
    encoding: Encoding,
    array_length: Option<usize>,
    options: &WriteOptions,
    prepared: &mut Option<std::slice::Iter<'_, numpress_transport::PreparedArray>>,
    header: Option<&header::ArrayHeader>,
) -> Result<()> {
    let (encoded, encoding, mode) = if let Some(arrays) = prepared {
        let array = arrays.next().expect("preflight array order matches writer");
        (
            std::borrow::Cow::Borrowed(array.encoded.as_str()),
            array.encoding,
            array.mode,
        )
    } else {
        let bytes = bytes();
        let bytes = if options.zlib_compression {
            let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
            encoder.write_all(&bytes)?;
            encoder.finish()?
        } else {
            bytes
        };
        let encoded = STANDARD.encode(bytes);
        (
            std::borrow::Cow::Owned(encoded),
            encoding,
            NumpressCompression::None,
        )
    };
    write!(w, "<binaryDataArray encodedLength=\"{}\"", encoded.len())?;
    if let Some(length) = array_length {
        write!(w, " arrayLength=\"{length}\"")?;
    }
    if let Some(header) = header {
        w.write_all(header.attrs.as_bytes())?;
    }
    writeln!(w, ">")?;
    let (accession, name) = match encoding {
        Encoding::Float32 => ("MS:1000521", "32-bit float"),
        Encoding::Float64 => ("MS:1000523", "64-bit float"),
        Encoding::Int32 => ("MS:1000519", "32-bit integer"),
        Encoding::Int64 => ("MS:1000522", "64-bit integer"),
        Encoding::Ascii => ("MS:1001479", "null-terminated ASCII string"),
    };
    cv(w, accession, name, "", "")?;
    if let Some((accession, name)) =
        numpress_transport::compression_term(mode, options.zlib_compression)
    {
        cv(w, accession, name, "", "")?;
    } else if options.zlib_compression {
        cv(w, "MS:1000574", "zlib compression", "", "")?;
    } else {
        cv(w, "MS:1000576", "no compression", "", "")?;
    }
    match kind {
        Kind::Role(role) => {
            let (accession, name, unit, unit_name) = role.terms();
            let unit = if unit.is_empty() {
                String::new()
            } else {
                format!(
                    " unitCvRef=\"{}\" unitAccession=\"{unit}\" unitName=\"{unit_name}\"",
                    unit.split_once(':').unwrap().0
                )
            };
            cv(
                w,
                accession,
                name,
                if role == record_transport::Role::Detector {
                    "detector signal"
                } else {
                    ""
                },
                &unit,
            )?;
        }
        Kind::Mz => cv(
            w,
            "MS:1000514",
            "m/z array",
            "",
            " unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"",
        )?,
        Kind::Time => cv(w, "MS:1000595", "time array", "", SECOND)?,
        // Source `MzMLHandler.cpp:5688` always writes the counts unit here.
        Kind::Intensity => cv(w, "MS:1000515", "intensity array", "", COUNTS)?,
        Kind::Auxiliary(name) => {
            if let Some(accession) = canonical_array_accession(&name) {
                cv(w, accession, &name, "", "")?;
            } else {
                cv(w, "MS:1000786", "non-standard data array", &name, "")?;
            }
        }
    }
    if let Some(header) = header {
        w.write_all(header.params.as_bytes())?;
    }
    writeln!(w, "<binary>{encoded}</binary></binaryDataArray>")?;
    Ok(())
}

fn write_auxiliary_arrays(
    w: &mut impl Write,
    floats: &[DataArray<f32>],
    integers: &[DataArray<i32>],
    strings: &[DataArray<String>],
    options: &WriteOptions,
    prepared: &mut Option<std::slice::Iter<'_, numpress_transport::PreparedArray>>,
    headers: &mut std::slice::Iter<'_, header::ArrayHeader>,
) -> Result<()> {
    for array in floats {
        write_array(
            w,
            || array.data.iter().flat_map(|v| v.to_le_bytes()).collect(),
            Kind::Auxiliary(array.name.clone()),
            Encoding::Float32,
            Some(array.data.len()),
            options,
            prepared,
            headers.next(),
        )?;
    }
    for array in integers {
        // Source ordinary annotations use i64; canonical charge requires i32.
        let encoding = integer_array_encoding(&array.name);
        let bytes = || {
            if matches!(encoding, Encoding::Int32) {
                array.data.iter().flat_map(|v| v.to_le_bytes()).collect()
            } else {
                array
                    .data
                    .iter()
                    .flat_map(|&v| i64::from(v).to_le_bytes())
                    .collect()
            }
        };
        write_array(
            w,
            bytes,
            Kind::Auxiliary(array.name.clone()),
            encoding,
            Some(array.data.len()),
            options,
            prepared,
            headers.next(),
        )?;
    }
    for array in strings {
        write_array(
            w,
            || {
                array
                    .data
                    .iter()
                    .flat_map(|v| v.bytes().chain(std::iter::once(0)))
                    .collect()
            },
            Kind::Auxiliary(array.name.clone()),
            Encoding::Ascii,
            Some(array.data.len()),
            options,
            prepared,
            headers.next(),
        )?;
    }
    Ok(())
}

fn check_array_descriptions(
    floats: &[DataArray<f32>],
    integers: &[DataArray<i32>],
    strings: &[DataArray<String>],
) -> Result<()> {
    for meta in floats
        .iter()
        .map(|a| &a.metadata)
        .chain(integers.iter().map(|a| &a.metadata))
        .chain(strings.iter().map(|a| &a.metadata))
    {
        validate_scalar_metadata(meta)?;
    }
    Ok(())
}
fn check_auxiliary_arrays(
    floats: &[DataArray<f32>],
    integers: &[DataArray<i32>],
    strings: &[DataArray<String>],
    nonfinite_floats: bool,
) -> Result<()> {
    check_array_descriptions(floats, integers, strings)?;
    for array in floats {
        check_canonical_encoding(&array.name, Encoding::Float32)?;
    }
    for array in integers {
        check_canonical_encoding(&array.name, integer_array_encoding(&array.name))?;
    }
    for array in strings {
        check_canonical_encoding(&array.name, Encoding::Ascii)?;
    }
    let mut names = BTreeSet::new();
    for name in floats
        .iter()
        .map(|a| &a.name)
        .chain(integers.iter().map(|a| &a.name))
        .chain(strings.iter().map(|a| &a.name))
    {
        xml_string(name)?;
        if name.is_empty() || !names.insert(name) {
            return Err(Error::InvalidValue(
                "empty or duplicate auxiliary array name".into(),
            ));
        }
    }
    if !nonfinite_floats && floats.iter().flat_map(|a| &a.data).any(|v| !v.is_finite()) {
        return Err(Error::InvalidValue(
            "nonfinite auxiliary float value".into(),
        ));
    }
    for value in strings.iter().flat_map(|a| &a.data) {
        if !value.is_ascii() || value.as_bytes().contains(&0) {
            return Err(Error::Unsupported(
                "mzML string arrays require ASCII without embedded NUL".into(),
            ));
        }
    }
    Ok(())
}

/// Write indexed mzML 1.1 (`indexedmzML`) with uncompressed binary arrays, as
/// source `MzMLFile::store` does with its default `PeakFileOptions`
/// (`write_index_ = true`, `PeakFileOptions.h:244`).
///
/// This is the writer behind `FileHandler::store_experiment` and therefore
/// behind every TOPP tool that stores mzML. It makes one pass: each record's
/// byte offset is taken as its `<spectrum` or `<chromatogram` tag starts, the
/// index follows `</mzML>`, and `fileChecksum` is the SHA-1 of every byte from
/// the start of the document through the opening `<fileChecksum>` tag, as the
/// indexed mzML schema specifies. The source writes the constant `0` there
/// (CPP-049). Offsets count bytes of the XML text written to `writer`.
///
/// An experiment with neither spectra nor chromatograms has no record to
/// index and is written as plain mzML; the source emits an index with a dummy
/// `-1` offset instead (CPP-050).
///
/// Validation and the header plan complete before the first byte is written,
/// so a rejected experiment leaves `writer` untouched. For binary encoding
/// options and the prepared two-pass writer use
/// [`write_with_peak_options`]; for plain mzML use [`write_with_options`].
///
/// # Errors
///
/// Returns the validation errors of [`write_with_options`] and any I/O error.
pub fn write(writer: impl Write, experiment: &MSExperiment) -> Result<()> {
    write_indexed(writer, experiment, false)
}

/// [`write`], accepting NaN and infinite values in auxiliary float arrays, as
/// source `MzMLFile::store` does.
///
/// The default writer rejects them (`nonfinite auxiliary float value`). The
/// source writes whatever the arrays hold; `FeatureFinderAlgorithmPicked`'s
/// debug mode stores NaN trace scores this way. Coordinates and intensities
/// must still be finite. [`ReadOptions::source_nonfinite_float_arrays`] reads
/// such a file back.
///
/// # Errors
///
/// As [`write`], without the non-finite auxiliary check.
pub fn write_source_float_arrays(writer: impl Write, experiment: &MSExperiment) -> Result<()> {
    write_indexed(writer, experiment, true)
}

fn write_indexed(
    writer: impl Write,
    experiment: &MSExperiment,
    nonfinite_floats: bool,
) -> Result<()> {
    if experiment.spectra.is_empty() && experiment.chromatograms.is_empty() {
        return write_with_options(writer, experiment, &WriteOptions::default());
    }
    let header = header::prepare(experiment)?;
    validate_write_arrays(experiment, nonfinite_floats)?;
    let mut output = peak_writer::Output::streamed(writer, experiment)?;
    write_document(
        &mut output,
        experiment,
        &WriteOptions::default(),
        &mut None,
        &header,
        false,
    )
}

/// Write plain (unindexed) mzML 1.1 after preflight validation.
///
/// Named float, integer and ASCII string arrays are preserved. Unsupported
/// metadata or unrepresentable array values are rejected before output. The
/// streaming `MSDataWritingConsumer` splits this layout per record, which is
/// why it stays unindexed; [`write()`] is the indexed default.
pub fn write_with_options(
    mut w: impl Write,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    let header = header::prepare(experiment)?;
    validate_write(experiment)?;
    write_impl(&mut w, experiment, options, &mut None, &header)
}
fn experiment_header_guard(experiment: &MSExperiment) -> Result<()> {
    header::guard(experiment)
}
/// How many whole-document allowances an mzML write of `experiment` receives:
/// one for the experiment-level header and one more for every spectrum and
/// chromatogram.
///
/// The writer preflights (header plan, settings validation and the prepared
/// writers' markup, index and binary budgets) used to share one fixed
/// allowance across the whole document. That refused realistic runs after
/// about 650 records: the 2026-09-14 smoke benchmark measured 647 passing and
/// 648 failing spectra on `UK222_picked`, while the C++ Release writer stores
/// the complete 44k-spectrum runs. Every such allowance is now multiplied by
/// this share count, so a ceiling grows linearly with the records it has to
/// cover and still bounds amplification within the document. The source
/// enforces no ceilings at all.
fn writer_shares(experiment: &MSExperiment) -> usize {
    experiment
        .spectra
        .len()
        .saturating_add(experiment.chromatograms.len())
        .saturating_add(1)
}
fn validate_write(experiment: &MSExperiment) -> Result<()> {
    validate_write_arrays(experiment, false)
}
fn validate_write_arrays(experiment: &MSExperiment, nonfinite_floats: bool) -> Result<()> {
    experiment_header_guard(experiment)?;
    // O(1) loss guards precede validation of newly supported owned settings.
    for spectrum in &experiment.spectra {
        settings_metadata::spectrum_guard(spectrum)?;
        record_transport::array_owners(
            &spectrum.metadata,
            false,
            &spectrum.float_data_arrays,
            &spectrum.integer_data_arrays,
            &spectrum.string_data_arrays,
        )?;
        check_array_descriptions(
            &spectrum.float_data_arrays,
            &spectrum.integer_data_arrays,
            &spectrum.string_data_arrays,
        )?;
    }
    for chromatogram in &experiment.chromatograms {
        settings_metadata::chromatogram_guard(chromatogram)?;
        record_transport::array_owners(
            &chromatogram.metadata,
            true,
            &chromatogram.float_data_arrays,
            &chromatogram.integer_data_arrays,
            &chromatogram.string_data_arrays,
        )?;
        check_array_descriptions(
            &chromatogram.float_data_arrays,
            &chromatogram.integer_data_arrays,
            &chromatogram.string_data_arrays,
        )?;
    }
    // Cumulative settings preflight, independent of binary encoding: one fixed
    // allowance per record share (`writer_shares`). Cover owned scalar
    // metadata before validation/rendering traverses it.
    let shares = writer_shares(experiment);
    let mut settings_work = 50_000_000usize.saturating_mul(shares);
    let mut settings_bytes = (256usize * 1024 * 1024).saturating_mul(shares);
    let initial_work = settings_work;
    let initial_bytes = settings_bytes;
    experiment
        .settings
        .with_budget(&mut settings_work, &mut settings_bytes)?;
    // Cover XML validation, scalar rendering and worst-case entity escaping
    // before any typed run value is formatted. Separate from binary budgets.
    settings_work = settings_work
        .checked_sub(
            (initial_work - settings_work)
                .checked_mul(7)
                .ok_or_else(|| invalid("mzML run metadata work overflow"))?,
        )
        .ok_or_else(|| invalid("mzML run metadata work limit"))?;
    settings_bytes = settings_bytes
        .checked_sub(
            (initial_bytes - settings_bytes)
                .checked_mul(7)
                .ok_or_else(|| invalid("mzML run metadata byte overflow"))?,
        )
        .ok_or_else(|| invalid("mzML run metadata byte limit"))?;
    settings_work = settings_work
        .checked_sub(
            experiment
                .spectra
                .len()
                .checked_mul(2)
                .and_then(|n| n.checked_add(256))
                .ok_or_else(|| invalid("mzML settings work overflow"))?,
        )
        .ok_or_else(|| invalid("mzML settings record limit exceeded"))?;
    for spectrum in &experiment.spectra {
        spectrum.record_metadata_with_budget(&mut settings_work, &mut settings_bytes)?;
        spectrum.acquisition_with_budget(&mut settings_work, &mut settings_bytes)?;
        settings_metadata::validate_spectrum(spectrum)?;
    }
    for chromatogram in &experiment.chromatograms {
        chromatogram.record_metadata_with_budget(&mut settings_work, &mut settings_bytes)?;
    }
    // New independent grids are never bounded by the aligned peak count. Cover
    // ordinary binary bytes, compression scratch and Base64 before emission.
    for spectrum in &experiment.spectra {
        for index in 0..3 {
            if let Some(values) = record_transport::noise_values(&spectrum.metadata, index)? {
                let cost = values
                    .len()
                    .checked_mul(64)
                    .and_then(|n| n.checked_add(4096))
                    .ok_or_else(|| invalid("noise output budget overflow"))?;
                settings_work = settings_work
                    .checked_sub(cost)
                    .ok_or_else(|| invalid("noise output work limit"))?;
                settings_bytes = settings_bytes
                    .checked_sub(cost)
                    .ok_or_else(|| invalid("noise output byte limit"))?;
            }
        }
    }
    experiment.validate()?;
    if experiment.settings.metadata.contains_key(NAME_KEY) {
        return Err(invalid("reserved record name userParam at run level"));
    }
    validate_scalar_metadata(&experiment.settings.metadata)?;
    let mut spectrum_ids = BTreeSet::new();
    let mut chromatogram_ids = BTreeSet::new();
    for (i, s) in experiment.spectra.iter().enumerate() {
        for precursor in &s.precursors {
            precursor_metadata::validate_write(precursor)?;
        }
        if !s.peptide_identifications.is_empty() {
            return Err(Error::Unsupported(
                "mzML writer cannot store peptide identifications".into(),
            ));
        }
        validate_spectrum_mobility(s)?;
        xml_string(&s.name)?;
        record_transport::validate(&s.metadata, false)?;
        let id = if s.native_id.is_empty() {
            format!("index={i}")
        } else {
            s.native_id.clone()
        };
        xml_string(&id)?;
        if id.split(' ').any(|part| {
            part.split_once('=').is_none_or(|(key, value)| {
                key.is_empty()
                    || value.is_empty()
                    || key.chars().any(char::is_whitespace)
                    || value.chars().any(char::is_whitespace)
            })
        }) {
            return Err(Error::InvalidValue(
                "spectrum native_id must have mzML key=value form".into(),
            ));
        }
        if !spectrum_ids.insert(id) {
            return Err(Error::InvalidValue("duplicate spectrum native_id".into()));
        }
        check_auxiliary_arrays(
            &s.float_data_arrays,
            &s.integer_data_arrays,
            &s.string_data_arrays,
            nonfinite_floats,
        )?;
    }
    for (i, c) in experiment.chromatograms.iter().enumerate() {
        precursor_metadata::validate_write(&c.precursor)?;
        validate_product_write(&c.product)?;
        xml_string(&c.name)?;
        record_transport::validate(&c.metadata, true)?;
        let id = if c.native_id.is_empty() {
            format!("chromatogram={i}")
        } else {
            c.native_id.clone()
        };
        xml_string(&id)?;
        if !chromatogram_ids.insert(id) {
            return Err(Error::InvalidValue(
                "duplicate chromatogram native_id".into(),
            ));
        }
        check_auxiliary_arrays(
            &c.float_data_arrays,
            &c.integer_data_arrays,
            &c.string_data_arrays,
            nonfinite_floats,
        )?;
    }
    // A precursor `spectrumRef` that names no spectrum of this document is
    // written as it stands, as source `MzMLHandler::writePrecursor_` does: it
    // emits the `spectrum_ref` meta value verbatim and resolves nothing
    // (`MzMLHandler.cpp:4535-4547` at core bc9cc12). Requiring resolution here
    // made the port unable to split a real run: `MzMLSplitter` moves whole
    // spectra into parts, so every MS2 whose precursor stayed in an earlier
    // part carries a reference out of the part. The C++ tool writes those
    // parts, and the benchmark's 1.2 GB Velos run splits into four parts of
    // which part 2 alone carries four such references (to `scan=10929`, which
    // part 1 holds). The reference is still checked to be a well-formed XML
    // string by `precursor_metadata::validate_write` above, which is the
    // property the writer can actually establish about it.
    Ok(())
}
fn write_impl(
    w: impl Write,
    experiment: &MSExperiment,
    options: &WriteOptions,
    prepared: &mut Option<std::slice::Iter<'_, numpress_transport::PreparedArray>>,
    header: &header::Plan,
) -> Result<()> {
    let mut w = peak_writer::Output::legacy(w);
    write_document(&mut w, experiment, options, prepared, header, false)
}
fn write_document<W: Write>(
    mut w: &mut peak_writer::Output<'_, W>,
    experiment: &MSExperiment,
    options: &WriteOptions,
    prepared: &mut Option<std::slice::Iter<'_, numpress_transport::PreparedArray>>,
    header: &header::Plan,
    tpp: bool,
) -> Result<()> {
    w.header(&header.prefix)?;
    let mut array_headers = header.arrays.iter();
    if !experiment.spectra.is_empty() {
        writeln!(
            w,
            "<spectrumList count=\"{}\" defaultDataProcessingRef=\"dp_00000000000000000000\">",
            experiment.spectra.len()
        )?;
        for (i, spectrum) in experiment.spectra.iter().enumerate() {
            let id = peak_writer::native_id(&spectrum.native_id, i, false);
            w.record(false)?;
            writeln!(
                w,
                "<spectrum id=\"{}\" index=\"{i}\" defaultArrayLength=\"{}\"{}>",
                escape(&id),
                spectrum.len(),
                header.spectra[i]
            )?;
            cv(
                &mut w,
                "MS:1000511",
                "ms level",
                &spectrum.ms_level.to_string(),
                "",
            )?;
            match spectrum.spectrum_type {
                SpectrumType::Centroid => cv(&mut w, "MS:1000127", "centroid spectrum", "", "")?,
                SpectrumType::Profile => cv(&mut w, "MS:1000128", "profile spectrum", "", "")?,
                SpectrumType::Unknown => {}
            }
            settings_metadata::write_spectrum(&mut w, spectrum)?;
            w.write_all(header.spectrum_metadata[i].as_bytes())?;
            settings_metadata::write_scan(&mut w, spectrum)?;
            if !spectrum.precursors.is_empty() {
                writeln!(w, "<precursorList count=\"{}\">", spectrum.precursors.len())?;
                for precursor in &spectrum.precursors {
                    write_precursor(&mut w, precursor, tpp)?;
                }
                writeln!(w, "</precursorList>")?;
            }
            if !spectrum.products.is_empty() {
                writeln!(w, "<productList count=\"{}\">", spectrum.products.len())?;
                for product in &spectrum.products {
                    write_product(&mut w, product)?;
                }
                writeln!(w, "</productList>")?;
            }
            writeln!(
                w,
                "<binaryDataArrayList count=\"{}\">",
                2 + record_transport::NOISE
                    .iter()
                    .filter(|key| spectrum.metadata.contains_key(**key))
                    .count()
                    + spectrum.float_data_arrays.len()
                    + spectrum.integer_data_arrays.len()
                    + spectrum.string_data_arrays.len()
            )?;
            write_array(
                &mut w,
                || {
                    spectrum
                        .peaks
                        .iter()
                        .flat_map(|p| p.mz.to_le_bytes())
                        .collect()
                },
                record_transport::primary_kind(&spectrum.metadata, false, true)?,
                Encoding::Float64,
                None,
                options,
                prepared,
                None,
            )?;
            write_array(
                &mut w,
                || {
                    spectrum
                        .peaks
                        .iter()
                        .flat_map(|p| p.intensity.to_le_bytes())
                        .collect()
                },
                record_transport::primary_kind(&spectrum.metadata, false, false)?,
                Encoding::Float32,
                None,
                options,
                prepared,
                None,
            )?;
            for index in 0..3 {
                if let Some(values) = record_transport::noise_values(&spectrum.metadata, index)? {
                    write_array(
                        &mut w,
                        || values.iter().flat_map(|v| v.to_le_bytes()).collect(),
                        Kind::Role(record_transport::noise_role(index)),
                        Encoding::Float64,
                        Some(values.len()),
                        options,
                        prepared,
                        None,
                    )?;
                }
            }
            write_auxiliary_arrays(
                &mut w,
                &spectrum.float_data_arrays,
                &spectrum.integer_data_arrays,
                &spectrum.string_data_arrays,
                options,
                prepared,
                &mut array_headers,
            )?;
            writeln!(w, "</binaryDataArrayList></spectrum>")?;
        }
        writeln!(w, "</spectrumList>")?;
    }
    if !experiment.chromatograms.is_empty() {
        writeln!(
            w,
            "<chromatogramList count=\"{}\" defaultDataProcessingRef=\"dp_00000000000000000000\">",
            experiment.chromatograms.len()
        )?;
        for (i, chromatogram) in experiment.chromatograms.iter().enumerate() {
            let id = peak_writer::native_id(&chromatogram.native_id, i, true);
            w.record(true)?;
            writeln!(
                w,
                "<chromatogram id=\"{}\" index=\"{i}\" defaultArrayLength=\"{}\"{}>",
                escape(&id),
                chromatogram.len(),
                header.chromatograms[i]
            )?;
            settings_metadata::write_chromatogram(&mut w, chromatogram)?;
            w.write_all(header.chromatogram_metadata[i].as_bytes())?;
            if tpp || chromatogram.precursor != Precursor::default() {
                write_precursor(&mut w, &chromatogram.precursor, tpp)?;
            }
            write_product(&mut w, &chromatogram.product)?;
            writeln!(
                w,
                "<binaryDataArrayList count=\"{}\">",
                2 + chromatogram.float_data_arrays.len()
                    + chromatogram.integer_data_arrays.len()
                    + chromatogram.string_data_arrays.len()
            )?;
            write_array(
                &mut w,
                || {
                    chromatogram
                        .peaks
                        .iter()
                        .flat_map(|p| p.rt.to_le_bytes())
                        .collect()
                },
                record_transport::primary_kind(&chromatogram.metadata, true, true)?,
                Encoding::Float64,
                None,
                options,
                prepared,
                None,
            )?;
            write_array(
                &mut w,
                || {
                    chromatogram
                        .peaks
                        .iter()
                        .flat_map(|p| p.intensity.to_le_bytes())
                        .collect()
                },
                record_transport::primary_kind(&chromatogram.metadata, true, false)?,
                Encoding::Float32,
                None,
                options,
                prepared,
                None,
            )?;
            write_auxiliary_arrays(
                &mut w,
                &chromatogram.float_data_arrays,
                &chromatogram.integer_data_arrays,
                &chromatogram.string_data_arrays,
                options,
                prepared,
                &mut array_headers,
            )?;
            writeln!(w, "</binaryDataArrayList></chromatogram>")?;
        }
        writeln!(w, "</chromatogramList>")?;
    }
    writeln!(w, "</run></mzML>")?;
    w.footer(experiment)?;
    w.flush()?;
    Ok(())
}

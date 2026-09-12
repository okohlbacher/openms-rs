// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! imzML imaging index, dataset geometry and companion `.ibd` array reads.
//!
//! Ports `FORMAT/HANDLERS/ImzMLHandlerHelper.h` and
//! `FORMAT/HANDLERS/ImzMLHandler.h`. See `docs/IMZML_HANDLER_SUPPORT.md`.
//!
//! An imzML dataset (HUPO-PSI imzML 1.1.0) is two files. The `.imzML` is mzML
//! XML that carries metadata plus, per spectrum, the IMS CV params giving the
//! byte offset (`IMS:1000102`), element count (`IMS:1000103`) and encoded byte
//! length (`IMS:1000104`) of that spectrum's arrays inside the companion `.ibd`
//! binary file. The `.ibd` begins with a 16-byte UUID that must equal the
//! `IMS:1000080` identifier in the XML. Two storage modes exist: *continuous*
//! (`IMS:1000030`, one shared m/z array stored once) and *processed*
//! (`IMS:1000031`, a private m/z array per spectrum);
//! [`ImagingMode`](crate::format::imzml_handler::ImagingMode) records which.
//!
//! The division of labour matches the source. C++ `ImzMLHandler` derives from
//! `MzMLHandler` and intercepts only the IMS terms its base does not know;
//! everything else — instrument, data processing, retention time, MS level — is
//! the base class's work. Here the base class is
//! [`mzml`](crate::format::mzml) and this module is the interception layer:
//! [`read_index`](crate::format::imzml_handler::read_index) scans an `.imzML`
//! for IMS terms only, and
//! [`ImzMLHandler`](crate::format::imzml_handler::ImzMLHandler) pairs the
//! resulting index with the `.ibd` to decode one pixel's arrays per call. There
//! is no second mzML parser here: nothing in this module interprets
//! referenceable groups or header metadata beyond the IMS vocabulary.
//!
//! One binary payload is the exception. A `<binary>` element belonging to an
//! m/z or intensity array that carries **no** `IMS:1000101` holds that array
//! inline, and the source fills such an array from the peaks its `MzMLHandler`
//! base decoded rather than from the `.ibd` (`ImzMLHandler.cpp:214-232`). This
//! module has no base class to borrow those peaks from, so it keeps that one
//! payload — whitespace removed, charged against
//! [`max_text_bytes`](crate::format::imzml_handler::ImzMLReadLimits::max_text_bytes)
//! as it accumulates — and decodes it at the same point the `.ibd` would have
//! been read. Nothing else inline is kept: an auxiliary array without
//! `IMS:1000101` is reported as
//! [`AuxSkipReason::Inline`](crate::format::imzml_handler::AuxSkipReason::Inline)
//! and dropped, as the source drops it.
//!
//! Every offset and length in the index comes from the XML and indexes into a
//! second file, so this is an attacker-shaped format. Each read is preflighted
//! against [`ImzMLReadLimits`](crate::format::imzml_handler::ImzMLReadLimits)
//! **and** against the actual `.ibd` length before anything is allocated, so a
//! malformed or hostile index costs a bounded amount of memory and nothing
//! else. The source checks the element count against a ceiling and lets a short
//! `fread` fail afterwards, which allocates first and discovers the truncation
//! second.
//!
//! Nothing here starts a thread. Every decode takes `&mut self`, so the
//! compiler enforces the exclusive access that the source's shared `FILE*`
//! only documents.

use crate::format::controlled_vocabulary::ControlledVocabulary;
use crate::kernel::{DataArray, MSSpectrum, Peak1D};
use crate::metadata::MetaValue;
use crate::{Error, Result};
use base64::{Engine, engine::general_purpose::STANDARD};
use quick_xml::{Reader, events::Event};
use sha1::{Digest, Sha1};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

/// Length in bytes of the UUID header at the start of a `.ibd` file.
///
/// Source `ImzMLFile.cpp` reads `unsigned char actual[16]` from offset 0 and
/// compares it with the `IMS:1000080` identifier declared in the XML.
pub const IBD_UUID_BYTES: usize = 16;

/// Scalar type of the elements of one `.ibd` binary array.
///
/// Source `ImzMLSpectrumIndex::DataType`. `Unknown` is the source's `UNKNOWN`:
/// the array carried none of the four PSI-MS binary-data-type terms, which is
/// not an error while the array is only indexed and becomes one when it is
/// decoded.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum ImzMLDataType {
    /// `MS:1000521`, 32-bit IEEE float.
    Float32,
    /// `MS:1000523`, 64-bit IEEE float.
    Float64,
    /// `MS:1000519`, signed 32-bit integer.
    Int32,
    /// `MS:1000522`, signed 64-bit integer.
    Int64,
    /// No supported binary-data-type term was seen on the array.
    #[default]
    Unknown,
}

impl ImzMLDataType {
    /// Stored width of one element in bytes, or `None` for
    /// [`Unknown`](Self::Unknown).
    pub const fn width(self) -> Option<usize> {
        match self {
            Self::Float32 | Self::Int32 => Some(4),
            Self::Float64 | Self::Int64 => Some(8),
            Self::Unknown => None,
        }
    }

    /// The dataset-level spelling the source records in
    /// `ImzMLMeta::mz_data_type` and `ImzMLMeta::int_data_type`.
    ///
    /// These are the exact strings of source `ImzMLInterceptConsumer::dtStr_`,
    /// including `"unknown"`, which that function returns for `UNKNOWN`.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Float32 => "float32",
            Self::Float64 => "float64",
            Self::Int32 => "int32",
            Self::Int64 => "int64",
            Self::Unknown => "unknown",
        }
    }

    /// The type a PSI-MS binary-data-type accession selects, or `None` when the
    /// accession is not one of the four the source recognises.
    pub fn from_accession(accession: &str) -> Option<Self> {
        match accession {
            "MS:1000521" => Some(Self::Float32),
            "MS:1000523" => Some(Self::Float64),
            "MS:1000519" => Some(Self::Int32),
            "MS:1000522" => Some(Self::Int64),
            _ => None,
        }
    }
}

/// Which of the two imzML binary-layout modes a dataset declares.
///
/// The source stores this as a `std::string` that is `"continuous"`,
/// `"processed"` or empty; `Option<ImagingMode>` is that same three-way state,
/// and [`as_str`](Self::as_str) reproduces the source spelling exactly so a
/// writer or a `MetaValue` keeps the on-disk vocabulary.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum ImagingMode {
    /// `IMS:1000030`: every spectrum shares one m/z array, stored once.
    Continuous,
    /// `IMS:1000031`: every spectrum carries its own m/z array.
    Processed,
}

impl ImagingMode {
    /// The source's string for this mode, `"continuous"` or `"processed"`.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Continuous => "continuous",
            Self::Processed => "processed",
        }
    }
}

/// Dataset-level imaging metadata of one imzML file.
///
/// Source `ImzMLMeta`. imzML extends mzML with a companion binary file and the
/// IMS ontology terms that describe imaging geometry, binary layout and file
/// checksum; this is the whole of that vocabulary as the source models it.
///
/// The source's separate `IMAGING/MSImagingGeometry.h` pixel grid is a
/// different header and is not part of this package; see
/// `docs/IMZML_HANDLER_SUPPORT.md`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ImzMLMeta {
    /// Max pixel column, `IMS:1000042`, raised to the largest observed
    /// [`ImzMLSpectrumIndex::x`].
    pub max_count_x: u32,
    /// Max pixel row, `IMS:1000043`, raised to the largest observed
    /// [`ImzMLSpectrumIndex::y`].
    pub max_count_y: u32,
    /// Max depth slice; no IMS term sets it, so it is the largest observed
    /// [`ImzMLSpectrumIndex::z`] and 1 for a 2-D dataset. Source default is 1
    /// and [`Default`] keeps 0, because an index with no spectra has observed
    /// nothing; [`read_index`] raises it to at least 1 for any dataset that has
    /// one.
    pub max_count_z: u32,
    /// Physical pixel width in µm, `IMS:1000046`.
    pub pixel_size_x: f64,
    /// Physical pixel height in µm, `IMS:1000047`.
    pub pixel_size_y: f64,
    /// Physical x extent in µm, `IMS:1000044`.
    pub max_dim_x: f64,
    /// Physical y extent in µm, `IMS:1000045`.
    pub max_dim_y: f64,
    /// Binary layout mode, or `None` when the file declares neither
    /// `IMS:1000030` nor `IMS:1000031`.
    pub imaging_mode: Option<ImagingMode>,
    /// Path of the companion `.ibd`, as source `ibd_file_path`. The source has
    /// the loader assign this after construction; here
    /// [`ImzMLHandler::open`] fills it in with the file it actually opened.
    pub ibd_file_path: PathBuf,
    /// SHA-1 checksum of the `.ibd` as declared by `IMS:1000091`, empty when
    /// absent. Parsed, and verifiable on demand with
    /// [`ImzMLHandler::verify_ibd_sha1`].
    pub ibd_sha1: String,
    /// MD5 checksum of the `.ibd` as declared by `IMS:1000090`, empty when
    /// absent. Parsed but never verified; this crate has no MD5 implementation.
    pub ibd_md5: String,
    /// Dataset UUID, `IMS:1000080`, as written in the XML including any dashes
    /// or braces. The source ignores an empty value and keeps the previous one.
    pub uuid: String,
    /// Data type of the first m/z array in document order that declares one.
    pub mz_data_type: ImzMLDataType,
    /// Data type of the first intensity array in document order that declares one.
    pub int_data_type: ImzMLDataType,
    /// `"top down"` (`IMS:1000401`) or `"bottom up"` (`IMS:1000402`), empty when
    /// absent.
    pub scan_pattern: String,
    /// `"flyback"` (`IMS:1000413`), `"meander"` (`IMS:1000412`),
    /// `"horizontal"` (`IMS:1000480`) or `"vertical"` (`IMS:1000481`), empty
    /// when absent.
    pub scan_direction: String,
    /// `"left-right"` (`IMS:1000491`) or `"right-left"` (`IMS:1000492`), empty
    /// when absent.
    pub line_scan_direction: String,
    /// `"positive"` (`MS:1000130`) or `"negative"` (`MS:1000129`), empty when
    /// absent.
    pub polarity: String,
}

/// Index entry for one auxiliary — neither m/z nor intensity — external array.
///
/// Source `ImzMLSpectrumIndex::AuxArray`. imzML permits extra
/// `binaryDataArray` entries per spectrum beyond the mandatory m/z and
/// intensity pair, most notably ion mobility (`MS:1003006`). The source keeps
/// them in the index so that an on-disc reader and an in-memory load expose the
/// same viewer contract: one equal-length float array per auxiliary array.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ImzMLAuxArray {
    /// Ontology or free-text array name, the name a decoded
    /// [`DataArray`] carries. For a child of `MS:1000513` it is the CV term's
    /// own name, not the XML `name` attribute; for `MS:1000786` it is the
    /// param's `value`.
    pub name: String,
    /// PSI-MS accession of the array-identity param, empty when none was seen.
    pub accession: String,
    /// `unitAccession` of the array-identity param, empty when absent.
    pub unit_accession: String,
    /// Byte offset in the `.ibd`, `IMS:1000102`.
    pub offset: u64,
    /// Element count, `IMS:1000103`.
    pub length: u64,
    /// Stored byte length, `IMS:1000104`. Only a compressed array needs it, and
    /// compressed arrays are rejected at decode time.
    pub encoded_bytes: u64,
    /// Scalar type of the stored elements.
    pub data_type: ImzMLDataType,
    /// True for any child of `MS:1000572` other than uncompressed
    /// `MS:1000576`. Such an array is rejected when decoded.
    pub compressed: bool,
}

/// Per-spectrum binary index entry for an imzML dataset.
///
/// Source `ImzMLSpectrumIndex`. It records where one pixel's arrays live in the
/// companion `.ibd` together with that pixel's image coordinates, so random
/// access needs one seek and one read per array and no XML re-parse.
///
/// Coordinates are 1-based, as imzML writes them.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ImzMLSpectrumIndex {
    /// 0-based document order. Source field is `int32_t index`; a count cannot
    /// be negative, so this is unsigned.
    pub index: u32,
    /// The `id` attribute of the `<spectrum>` element. Native addition: the
    /// source leaves native identifiers to its `MzMLHandler` base, and a
    /// stand-alone index is unusable for stage-3 record matching without it.
    pub native_id: String,
    /// Pixel column, 1-based, `IMS:1000050`.
    pub x: u32,
    /// Pixel row, 1-based, `IMS:1000051`.
    pub y: u32,
    /// Depth slice, 1-based, `IMS:1000052`; 1 when the param is absent.
    pub z: u32,
    /// Byte offset of the m/z array, `IMS:1000102`.
    pub mz_offset: u64,
    /// Element count of the m/z array, `IMS:1000103`.
    pub mz_length: u64,
    /// Stored byte length of the m/z array, `IMS:1000104`. Native addition: the
    /// source keeps `IMS:1000104` only for auxiliary arrays.
    pub mz_encoded_bytes: u64,
    /// Scalar type of the m/z array.
    pub mz_type: ImzMLDataType,
    /// True when the m/z array carries a compression term other than
    /// `MS:1000576`.
    pub mz_compressed: bool,
    /// True when the m/z array declares `IMS:1000101`, so its payload is in the
    /// `.ibd` rather than inline base64. Native addition: the source keeps this
    /// in its private `ArrayMeta` and decides there whether to read the `.ibd`
    /// or keep the base class's inline peaks.
    pub mz_external: bool,
    /// The m/z array's inline base64 payload with ASCII whitespace removed, or
    /// empty when [`Self::mz_external`] is set and the payload therefore lives
    /// in the `.ibd`.
    ///
    /// Native addition. The source reads the inline peaks back from the
    /// `MSSpectrum` its `MzMLHandler` base already populated
    /// (`ImzMLHandler.cpp:214-218`); this module parses the `.imzML` alone and
    /// so keeps the encoded text instead, decoding it in
    /// [`ImzMLHandler::spectrum`] and [`ImzMLHandler::mz_array`] where the
    /// `.ibd` read would otherwise happen. It is charged against
    /// [`ImzMLReadLimits::max_text_bytes`] while it accumulates.
    pub mz_inline: String,
    /// Byte offset of the intensity array, `IMS:1000102`.
    pub int_offset: u64,
    /// Element count of the intensity array, `IMS:1000103`.
    pub int_length: u64,
    /// Stored byte length of the intensity array, `IMS:1000104`. Native
    /// addition, as [`Self::mz_encoded_bytes`].
    pub int_encoded_bytes: u64,
    /// Scalar type of the intensity array.
    pub int_type: ImzMLDataType,
    /// True when the intensity array carries a compression term other than
    /// `MS:1000576`.
    pub int_compressed: bool,
    /// True when the intensity array declares `IMS:1000101`. Native addition,
    /// as [`Self::mz_external`].
    pub int_external: bool,
    /// The intensity array's inline base64 payload. Native addition, as
    /// [`Self::mz_inline`].
    pub int_inline: String,
    /// Extra external arrays, in document order.
    ///
    /// Only named arrays appear here, as in the source: its index builder skips
    /// an external array that no `MS:1000513` child and no `MS:1000786` named,
    /// because there is nothing to attach the values to. See
    /// [`Self::unnamed_aux`].
    pub aux: Vec<ImzMLAuxArray>,
    /// How many external auxiliary arrays were dropped for having no name.
    ///
    /// The source drops them from the index entry exactly as this does, and
    /// warns about each one from a separate per-spectrum snapshot that the
    /// index does not carry. This count is native, and it is what lets
    /// [`ImzMLHandler::spectrum`] report the same omission without keeping a
    /// nameless entry that would make `aux.len()` disagree with the source's.
    pub unnamed_aux: u32,
    /// Names of this spectrum's non-external auxiliary arrays. The source warns
    /// and drops these when the peaks come from the `.ibd`; it keeps them in its
    /// private `SpecIMS`, and they are public here so the drop is reportable
    /// rather than only loggable.
    pub inline_aux_names: Vec<String>,
}

/// A parsed `.imzML`: the dataset metadata plus one index entry per spectrum.
///
/// Source `ImzMLFile::loadSpectraIndex` hands back the same pair through two
/// out-parameters, `ImzMLMeta&` and `std::vector<ImzMLSpectrumIndex>&`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ImzMLIndex {
    /// Dataset-level imaging metadata.
    pub meta: ImzMLMeta,
    /// One entry per spectrum, in document order.
    pub spectra: Vec<ImzMLSpectrumIndex>,
}

impl ImzMLIndex {
    /// The number of indexed spectra, as source `getNrSpectra()`.
    pub fn len(&self) -> usize {
        self.spectra.len()
    }

    /// Whether the file indexed no spectra at all.
    pub fn is_empty(&self) -> bool {
        self.spectra.is_empty()
    }

    /// The entry at `index`, or `None` when `index` is out of range.
    pub fn get(&self, index: usize) -> Option<&ImzMLSpectrumIndex> {
        self.spectra.get(index)
    }

    /// Position of the spectrum acquired at 1-based pixel `(x, y, z)`, or
    /// `None` when no spectrum claims it.
    ///
    /// The first spectrum in document order wins a duplicated coordinate, which
    /// is the rule the source's geometry builder uses and the upstream suite
    /// asserts (`buildImagingGeometry tolerates duplicate pixels by default`:
    /// the first spectrum keeps the shared pixel). Duplicates are tolerated,
    /// not rejected, because the source's readers accept them.
    ///
    /// This lookup is a linear scan over the index, not a stored grid. The
    /// source builds an `MSImagingGeometry` during `open()` and answers in O(1)
    /// from it; that type belongs to a different header and this package does
    /// not own it.
    pub fn index_at_coord(&self, x: u32, y: u32, z: u32) -> Option<usize> {
        self.spectra
            .iter()
            .position(|entry| entry.x == x && entry.y == y && entry.z == z)
    }
}

/// Explicit ceilings checked before any allocation, seek or read.
///
/// The source guards only the element count, against a fixed
/// `MAX_IBD_ARRAY_ELEMENTS` of 100,000,000. Nothing compares the requested
/// range against the length of the `.ibd`, so the output vector is sized from
/// the declared count and a short `fread` reports the truncation once the
/// memory is already committed — 800 MB for a float64 array at that ceiling, and
/// 1.2 GB for a float32 one, which `readMzArray` stages through a second
/// `std::vector<float>` before widening. Every field below turns one of those
/// into a checked error raised before the allocation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ImzMLReadLimits {
    /// Maximum elements in one array, as source `MAX_IBD_ARRAY_ELEMENTS`.
    pub max_array_elements: u64,
    /// Maximum stored bytes of one array, checked against the declared count
    /// times the element width.
    pub max_array_bytes: usize,
    /// Maximum indexed spectra in one file.
    pub max_spectra: usize,
    /// Maximum auxiliary arrays on one spectrum, including those later skipped.
    pub max_aux_arrays: usize,
    /// Maximum auxiliary arrays across the whole file.
    pub max_total_aux_arrays: usize,
    /// Maximum `.imzML` bytes consumed by the index scan.
    pub max_xml_bytes: u64,
    /// Maximum referenceable parameter groups captured for replay.
    pub max_param_groups: usize,
    /// Maximum parameters captured across all groups.
    pub max_group_params: usize,
    /// Cumulative bytes of every string the index stores: identifiers, array
    /// names, accessions, checksums, captured group parameters and the inline
    /// base64 of a peak array that declares no `IMS:1000101`
    /// ([`ImzMLSpectrumIndex::mz_inline`]).
    ///
    /// The inline payload is charged before it is appended, so this is the
    /// ceiling on how much encoded peak data one scan retains, and on the peak
    /// allocation that retention costs. A conformant imzML 1.1.0 stores both
    /// peak arrays externally and charges nothing here.
    pub max_text_bytes: usize,
    /// Maximum distinct accessions resolved against the PSI-MS vocabulary.
    /// Results are memoised, so this bounds the vocabulary work of one scan.
    pub max_cv_lookups: usize,
    /// Maximum `.ibd` bytes hashed by [`ImzMLHandler::verify_ibd_sha1`].
    pub max_checksum_bytes: u64,
}

impl Default for ImzMLReadLimits {
    fn default() -> Self {
        Self {
            max_array_elements: 100_000_000,
            max_array_bytes: 256 << 20,
            max_spectra: 5_000_000,
            max_aux_arrays: 256,
            max_total_aux_arrays: 10_000_000,
            max_xml_bytes: 512 << 20,
            max_param_groups: 100_000,
            max_group_params: 1_000_000,
            max_text_bytes: 64 << 20,
            max_cv_lookups: 100_000,
            max_checksum_bytes: 16 << 30,
        }
    }
}

/// Outcome of comparing the `.ibd` UUID header with the XML's `IMS:1000080`.
///
/// Source `verifyIbdUuid_` logs a warning for every state but `Match` and lets
/// the load continue, deliberately, so that a non-conformant dataset still
/// opens. This port returns the verdict instead of writing to a log, so a
/// caller that wants strict conformance can enforce it and a caller that does
/// not is unaffected. [`ImzMLHandler::open`] does not reject any of these.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum UuidStatus {
    /// The 16 header bytes equal the declared identifier.
    Match,
    /// Both are present and differ. Usually the `.imzML` and `.ibd` do not
    /// belong together.
    Mismatch {
        /// Lower-case hex of the 16 bytes found at the start of the `.ibd`.
        found: String,
        /// Lower-case hex of the identifier the XML declares.
        declared: String,
    },
    /// The XML carries no `IMS:1000080`, or its value is not 32 hex digits
    /// after dashes and braces are removed.
    NotDeclared,
    /// The `.ibd` is shorter than [`IBD_UUID_BYTES`].
    IbdTooShort,
}

/// Outcome of recomputing a declared `.ibd` checksum.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ChecksumStatus {
    /// The recomputed digest equals the declared one, compared case-insensitively.
    Match,
    /// The recomputed digest differs from the declared one.
    Mismatch {
        /// Lower-case hex of the recomputed digest.
        found: String,
        /// The digest as the XML declares it.
        declared: String,
    },
    /// The XML declares no checksum of this kind.
    NotDeclared,
}

/// Why one auxiliary array of a spectrum produced no data array.
///
/// The source logs a warning for each of these and continues, so that one bad
/// auxiliary array cannot abort a whole dataset load — a policy the upstream
/// suite asserts (`load skips one bad aux length and still returns all
/// spectra`). This port keeps the policy and reports the skips in
/// [`DecodedSpectrum::skipped_aux`] rather than only logging them.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AuxSkipReason {
    /// No `MS:1000513` child and no `MS:1000786` gave the array a name, so it
    /// never entered the index. Counted by
    /// [`ImzMLSpectrumIndex::unnamed_aux`].
    Unnamed,
    /// `IMS:1000103` is zero, so there is nothing to read.
    ZeroLength,
    /// The declared element count is not the spectrum's peak count. Viewers
    /// require an auxiliary array to be as long as the peak array.
    LengthMismatch {
        /// The array's declared `IMS:1000103`.
        length: u64,
        /// Peaks actually decoded for this spectrum.
        peaks: usize,
    },
    /// The array carries none of the four supported binary-data-type terms.
    UnknownDataType,
    /// The array has no `IMS:1000101`, so its payload is inline base64 while
    /// the peaks come from the `.ibd`. The source does not decode these.
    Inline,
}

/// One auxiliary array that a decode skipped, with the reason.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SkippedAux {
    /// Array name, empty for an [`AuxSkipReason::Unnamed`] array, since a name
    /// is what it lacked.
    pub name: String,
    /// Why it was skipped.
    pub reason: AuxSkipReason,
}

/// One decoded pixel: the spectrum plus everything the decode chose not to
/// include.
///
/// The source's equivalent path returns only the `MSSpectrum` and writes the
/// rest to the log. Carrying the omissions in the return value is what lets a
/// caller notice them.
#[derive(Clone, Debug, PartialEq)]
pub struct DecodedSpectrum {
    /// Peaks, pixel coordinates as `imzml:x`/`imzml:y`/`imzml:z` meta values,
    /// and one float data array per decoded auxiliary array.
    pub spectrum: MSSpectrum,
    /// True when the m/z or intensity array is not external, so that side's
    /// values came from the array's inline base64 rather than from the `.ibd`.
    ///
    /// The source fills the non-external side from the peaks its
    /// `MzMLHandler` base decoded (`ImzMLHandler.cpp:214-232`); this module
    /// decodes [`ImzMLSpectrumIndex::mz_inline`] /
    /// [`ImzMLSpectrumIndex::int_inline`] instead. Either way the two sides
    /// end up the same length for a well-formed file, which is why the
    /// mismatch below it is a genuine-corruption guard and not the normal
    /// outcome for a mixed spectrum.
    ///
    /// A conformant imzML 1.1.0 stores both arrays externally, so this is
    /// `false` for every conformant dataset.
    pub inline_peaks: bool,
    /// Auxiliary arrays that produced no data array: the unnamed ones the
    /// index already dropped, then the named ones in document order, then the
    /// inline ones.
    pub skipped_aux: Vec<SkippedAux>,
}

/// Random-access reader for one `.ibd` companion file.
///
/// Source `ImzMLBinaryIO` is a class of static methods over a `FILE*` that
/// `ImzMLHandler` owns and closes in its destructor. Binding the handle, its
/// length and the ceilings into one value is what makes the preflight possible:
/// every read knows the file it is reading and can refuse a range that leaves
/// it.
///
/// imzML stores `.ibd` arrays little-endian. The source byte-swaps under
/// `OPENMS_IS_BIG_ENDIAN`; `f64::from_le_bytes` and friends do the same work on
/// every host with no conditional compilation, so there is no big-endian code
/// path to get wrong here.
#[derive(Debug)]
pub struct ImzMLBinaryIO {
    path: PathBuf,
    file: File,
    length: u64,
    limits: ImzMLReadLimits,
}

impl ImzMLBinaryIO {
    /// Open `path` for reading with default limits.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened or its length cannot be
    /// determined. Source `ImzMLHandler::openIBD` throws
    /// `Exception::FileNotFound` for the same condition.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Self::open_with_limits(path, ImzMLReadLimits::default())
    }

    /// Open `path` for reading with explicit ceilings.
    ///
    /// # Errors
    ///
    /// As [`open`](Self::open).
    pub fn open_with_limits(path: impl AsRef<Path>, limits: ImzMLReadLimits) -> Result<Self> {
        let path = path.as_ref();
        let mut file = File::open(path)?;
        let length = file.seek(SeekFrom::End(0))?;
        Ok(Self {
            path: path.to_path_buf(),
            file,
            length,
            limits,
        })
    }

    /// The path this reader was opened on.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// The length of the `.ibd` in bytes, as measured when it was opened.
    ///
    /// Every range check uses this value. A file that grows or shrinks
    /// underneath the reader is not re-measured, and a range that has become
    /// short fails on the read rather than the preflight.
    pub fn len(&self) -> u64 {
        self.length
    }

    /// Whether the `.ibd` is empty.
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    /// The ceilings this reader enforces.
    pub fn limits(&self) -> ImzMLReadLimits {
        self.limits
    }

    /// The 16-byte UUID header, or `None` when the file is shorter than
    /// [`IBD_UUID_BYTES`].
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the seek or read fails. A file too short to hold the
    /// header is `Ok(None)`, because the source treats it as a reason to skip
    /// verification rather than as an I/O failure.
    pub fn uuid(&mut self) -> Result<Option<[u8; IBD_UUID_BYTES]>> {
        if self.length < IBD_UUID_BYTES as u64 {
            return Ok(None);
        }
        let mut bytes = [0u8; IBD_UUID_BYTES];
        self.file.seek(SeekFrom::Start(0))?;
        self.file.read_exact(&mut bytes)?;
        Ok(Some(bytes))
    }

    /// Lower-case hex SHA-1 of the whole `.ibd`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the file is longer than
    /// [`ImzMLReadLimits::max_checksum_bytes`], and [`Error::Io`] when the read
    /// fails. The file is hashed in fixed-size chunks, so the digest costs a
    /// constant amount of memory whatever the file's size.
    pub fn sha1_hex(&mut self) -> Result<String> {
        if self.length > self.limits.max_checksum_bytes {
            return Err(Error::InvalidValue(
                ".ibd exceeds the configured checksum byte limit".into(),
            ));
        }
        self.file.seek(SeekFrom::Start(0))?;
        let mut digest = Sha1::new();
        let mut buffer = [0u8; 64 << 10];
        loop {
            let read = self.file.read(&mut buffer)?;
            if read == 0 {
                break;
            }
            digest.update(&buffer[..read]);
        }
        Ok(hex(&digest.finalize()))
    }

    /// Read `count` m/z values stored at `offset` with scalar type `data_type`.
    ///
    /// Source `ImzMLBinaryIO::readMzArray`. Integer types are widened to `f64`
    /// exactly as the source does; a 64-bit integer above 2^53 therefore loses
    /// precision in both.
    ///
    /// A `count` of zero yields an empty vector without touching the file, as
    /// the source's early return does, so a garbage offset on a zero-length
    /// array is not an error.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `count` exceeds
    /// [`ImzMLReadLimits::max_array_elements`], when the stored size exceeds
    /// [`ImzMLReadLimits::max_array_bytes`] or `usize`, or when the allocation
    /// fails; [`Error::Unsupported`] when `data_type` is
    /// [`ImzMLDataType::Unknown`]; [`Error::Parse`] when the range leaves the
    /// `.ibd`; [`Error::Io`] when the seek or read fails. The source raises
    /// `Exception::ParseError` for all of these and discovers the overlong
    /// range only when the read comes up short.
    pub fn read_mz_array(
        &mut self,
        offset: u64,
        count: u64,
        data_type: ImzMLDataType,
    ) -> Result<Vec<f64>> {
        let Some((count, bytes)) = self.preflight(offset, count, data_type, "m/z array")? else {
            return Ok(Vec::new());
        };
        let raw = self.read_exact_at(offset, bytes, "m/z array")?;
        widen_to_f64(&raw, data_type, count, "m/z array")
    }

    /// Read `count` intensity values stored at `offset` with scalar type
    /// `data_type`.
    ///
    /// Source `ImzMLBinaryIO::readIntArray`, which shares its body with
    /// `readAuxArray` through a private helper. Values become `f32` because
    /// that is how OpenMS stores intensities, so a `float64` or 64-bit integer
    /// array is narrowed exactly as the source narrows it.
    ///
    /// # Errors
    ///
    /// As [`read_mz_array`](Self::read_mz_array).
    pub fn read_intensity_array(
        &mut self,
        offset: u64,
        count: u64,
        data_type: ImzMLDataType,
    ) -> Result<Vec<f32>> {
        self.read_f32_array(offset, count, data_type, "intensity array")
    }

    /// Read `count` values of one auxiliary array stored at `offset`.
    ///
    /// Source `ImzMLBinaryIO::readAuxArray`. `name` appears in error messages,
    /// as the source's `array_name` does; an empty `name` becomes
    /// `"auxiliary array"`, again as the source does. Values are `f32` because
    /// OpenMS stores auxiliary spectrum data in a float data array.
    ///
    /// A compressed array — any child of `MS:1000572` other than uncompressed
    /// `MS:1000576` — must be rejected by the caller before this is invoked;
    /// nothing here inflates a payload. [`ImzMLHandler::spectrum`] performs that
    /// rejection.
    ///
    /// # Errors
    ///
    /// As [`read_mz_array`](Self::read_mz_array).
    pub fn read_aux_array(
        &mut self,
        offset: u64,
        count: u64,
        data_type: ImzMLDataType,
        name: &str,
    ) -> Result<Vec<f32>> {
        let label = if name.is_empty() {
            "auxiliary array"
        } else {
            name
        };
        self.read_f32_array(offset, count, data_type, label)
    }

    fn read_f32_array(
        &mut self,
        offset: u64,
        count: u64,
        data_type: ImzMLDataType,
        what: &str,
    ) -> Result<Vec<f32>> {
        let Some((count, bytes)) = self.preflight(offset, count, data_type, what)? else {
            return Ok(Vec::new());
        };
        let raw = self.read_exact_at(offset, bytes, what)?;
        narrow_to_f32(&raw, data_type, count, what)
    }

    /// `Ok(None)` for an empty array, otherwise the element count and the
    /// stored byte length, both already checked against the ceilings and
    /// against the length of the open file.
    fn preflight(
        &self,
        offset: u64,
        count: u64,
        data_type: ImzMLDataType,
        what: &str,
    ) -> Result<Option<(usize, usize)>> {
        if count == 0 {
            return Ok(None);
        }
        let width = data_type
            .width()
            .ok_or_else(|| Error::Unsupported(format!("unsupported {what} data type in .ibd")))?;
        if count > self.limits.max_array_elements {
            return Err(Error::InvalidValue(format!(
                "{what} element count {count} exceeds the configured limit of {}",
                self.limits.max_array_elements
            )));
        }
        let count = usize::try_from(count)
            .map_err(|_| Error::InvalidValue(format!("{what} element count exceeds usize")))?;
        let bytes = count
            .checked_mul(width)
            .filter(|&bytes| bytes <= self.limits.max_array_bytes)
            .ok_or_else(|| {
                Error::InvalidValue(format!("{what} exceeds the configured byte limit"))
            })?;
        let end = offset
            .checked_add(bytes as u64)
            .ok_or_else(|| Error::InvalidValue(format!("{what} byte range overflows u64")))?;
        if end > self.length {
            return Err(parse(format!(
                "{what} byte range {offset}..{end} extends past the {} byte .ibd",
                self.length
            )));
        }
        Ok(Some((count, bytes)))
    }

    fn read_exact_at(&mut self, offset: u64, bytes: usize, what: &str) -> Result<Vec<u8>> {
        let mut raw = try_vec::<u8>(bytes, what)?;
        self.file.seek(SeekFrom::Start(offset))?;
        Read::by_ref(&mut self.file)
            .take(bytes as u64)
            .read_to_end(&mut raw)?;
        if raw.len() != bytes {
            return Err(parse(format!(
                "the .ibd ended inside the requested {what} byte range"
            )));
        }
        Ok(raw)
    }
}

/// Write a float32 array to a `.ibd` stream, little-endian.
///
/// Source `ImzMLBinaryIO::writeFloat32Array`. The source guards the count
/// against `MAX_IBD_ARRAY_ELEMENTS` and byte-swaps on a big-endian host;
/// `f32::to_le_bytes` needs no host test. An empty slice writes nothing, as the
/// source's early return does.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the slice is longer than
/// [`ImzMLReadLimits::max_array_elements`], [`Error::Io`] when the write fails.
/// The source raises `Exception::ParseError` in both cases, and also for a null
/// `FILE*`, a state `impl Write` cannot be in.
pub fn write_float32_array(
    mut out: impl Write,
    data: &[f32],
    limits: &ImzMLReadLimits,
) -> Result<()> {
    check_write_count(data.len(), limits, "float32 array write")?;
    for &value in data {
        out.write_all(&value.to_le_bytes())?;
    }
    Ok(())
}

/// Write m/z values to a `.ibd` stream as float32, little-endian.
///
/// Source `ImzMLBinaryIO::writeMzAsFloat32`, which narrows every value with a
/// `static_cast<float>`; this does the same, so both lose the same precision.
///
/// # Errors
///
/// As [`write_float32_array`].
pub fn write_mz_as_float32(
    mut out: impl Write,
    mz: &[f64],
    limits: &ImzMLReadLimits,
) -> Result<()> {
    check_write_count(mz.len(), limits, "float32 array write")?;
    for &value in mz {
        out.write_all(&(value as f32).to_le_bytes())?;
    }
    Ok(())
}

/// Write a float64 array to a `.ibd` stream, little-endian.
///
/// Source `ImzMLBinaryIO::writeFloat64Array`.
///
/// # Errors
///
/// As [`write_float32_array`].
pub fn write_float64_array(
    mut out: impl Write,
    data: &[f64],
    limits: &ImzMLReadLimits,
) -> Result<()> {
    check_write_count(data.len(), limits, "float64 array write")?;
    for &value in data {
        out.write_all(&value.to_le_bytes())?;
    }
    Ok(())
}

/// Write m/z values to a `.ibd` stream as float64, little-endian.
///
/// Source `ImzMLBinaryIO::writeMzAsFloat64`, which forwards the vector
/// unchanged to `writeFloat64Array`.
///
/// # Errors
///
/// As [`write_float32_array`].
pub fn write_mz_as_float64(out: impl Write, mz: &[f64], limits: &ImzMLReadLimits) -> Result<()> {
    write_float64_array(out, mz, limits)
}

/// The 16 raw bytes of an imzML UUID string, or `None` when the string does not
/// hold exactly 32 hex digits.
///
/// Source `ImzMLFile.cpp`'s file-static `uuidStringToBytes_`, which the writer
/// mirrors so that reader and writer agree on the byte order stored in the
/// `.ibd` header. Dashes and the braces some writers add are removed first;
/// nothing else is.
pub fn uuid_bytes(uuid: &str) -> Option<[u8; IBD_UUID_BYTES]> {
    let mut digits = [0u8; IBD_UUID_BYTES * 2];
    let mut seen = 0;
    for byte in uuid.bytes() {
        if matches!(byte, b'-' | b'{' | b'}') {
            continue;
        }
        if seen == digits.len() || !byte.is_ascii_hexdigit() {
            return None;
        }
        digits[seen] = byte;
        seen += 1;
    }
    if seen != digits.len() {
        return None;
    }
    let mut out = [0u8; IBD_UUID_BYTES];
    for (target, pair) in out.iter_mut().zip(digits.chunks_exact(2)) {
        *target = nibble(pair[0]) * 16 + nibble(pair[1]);
    }
    Some(out)
}

/// The `.ibd` path that belongs to an `.imzML` path.
///
/// Source `ImzMLFile::inferIbdPath_`: a case-insensitive `.imzML` suffix is
/// replaced by `.ibd`, and any other name simply gains `.ibd`. That function is
/// private to a header this package does not own; it is reproduced here so a
/// handler can open a dataset from the `.imzML` path alone, which is what the
/// source's loaders do on the caller's behalf.
///
/// The suffix is matched on character boundaries with [`str::get`], never by
/// byte slicing: `path[len - 6..]` panics whenever the sixth-from-last byte is a
/// UTF-8 continuation byte, and every public load and store path in this family
/// reaches this function first, so `dir/日本語.txt` aborted the process before a
/// single validation ran.
///
/// The suffix is replaced by truncation, as the source's
/// `p.substr(0, p.size() - 6) + ".ibd"` does, not with
/// [`PathBuf::set_extension`]. Rust treats a name that is entirely `.imzML` as
/// having no extension, so `set_extension` appended and produced
/// `.imzML.ibd` where the source yields `.ibd`.
///
/// For a path that is not valid UTF-8 the suffix is tested against a lossy view.
/// The suffix itself is ASCII, so only a path whose final six bytes are
/// themselves malformed could be judged differently, and such a path names no
/// `.imzML` file.
pub fn infer_ibd_path(imzml_path: impl AsRef<Path>) -> PathBuf {
    let path = imzml_path.as_ref();
    let text = path.as_os_str().to_string_lossy();
    if let Some(cut) = text.len().checked_sub(6) {
        if text
            .get(cut..)
            .is_some_and(|suffix| suffix.eq_ignore_ascii_case(".imzml"))
        {
            return PathBuf::from(format!("{}.ibd", &text[..cut]));
        }
    }
    let mut appended = path.as_os_str().to_os_string();
    appended.push(".ibd");
    PathBuf::from(appended)
}

/// A reader for one imzML dataset: the parsed `.imzML` index plus its `.ibd`.
///
/// Source `Internal::ImzMLHandler`. The source's own documentation says not to
/// use the class directly but to go through `ImzMLFile::load()` or
/// `OnDiscImzMLExperiment::open()`; the same holds here once those layers exist.
/// What this type owns is the imzML-specific half: the IMS vocabulary, the
/// per-spectrum index, the dataset geometry and the `.ibd` reads. Retention
/// time, MS level, instrument and data processing come from the mzML half,
/// which is [`mzml`](crate::format::mzml).
///
/// The source parses the XML through a SAX handler that must intercept terms
/// before its base class sees them, and it re-associates spectra with their IMS
/// state through a counter because `MzMLHandler` delivers the whole spectrum
/// list at once. None of that machinery is needed here: the index is built in
/// one pass and handed over as a value, so there is no delivery order to
/// re-synchronise and no `skip_spectrum_` interaction to get wrong.
///
/// # Examples
///
/// ```
/// use openms::format::imzml_handler::{ImagingMode, ImzMLHandler, UuidStatus};
///
/// let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
///     .join("tests/data/ImzMLFile_1_Example_Continuous.imzML");
/// let mut handler = ImzMLHandler::open(path)?;
///
/// // Continuous mode: nine pixels on a 3x3 grid sharing one m/z array.
/// assert_eq!(handler.len(), 9);
/// assert_eq!(handler.meta().imaging_mode, Some(ImagingMode::Continuous));
/// assert_eq!((handler.meta().max_count_x, handler.meta().max_count_y), (3, 3));
/// assert_eq!(handler.uuid_status()?, UuidStatus::Match);
///
/// let decoded = handler.spectrum(0)?;
/// assert_eq!(decoded.spectrum.peaks.len(), 8399);
/// assert!((decoded.spectrum.peaks[0].mz - 100.0).abs() < 1e-4);
///
/// // Every pixel reads the same m/z array offset.
/// let shared = handler.index()[0].mz_offset;
/// assert!(handler.index().iter().all(|entry| entry.mz_offset == shared));
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Debug)]
pub struct ImzMLHandler {
    imzml_path: PathBuf,
    index: ImzMLIndex,
    ibd: ImzMLBinaryIO,
}

impl ImzMLHandler {
    /// Open `imzml_path` and the `.ibd` sibling [`infer_ibd_path`] names.
    ///
    /// # Errors
    ///
    /// As [`open_with_limits`](Self::open_with_limits).
    pub fn open(imzml_path: impl AsRef<Path>) -> Result<Self> {
        let imzml_path = imzml_path.as_ref();
        let ibd_path = infer_ibd_path(imzml_path);
        Self::open_with_limits(imzml_path, ibd_path, ImzMLReadLimits::default())
    }

    /// Open `imzml_path` with an explicit `.ibd` path.
    ///
    /// The source's `OnDiscImzMLExperiment::open(imzml, ibd)` threads the
    /// override through the index load *and* the UUID check, so that both
    /// target the file that will actually be read rather than an inferred
    /// sibling that may be missing or stale. This does the same.
    ///
    /// # Errors
    ///
    /// As [`open_with_limits`](Self::open_with_limits).
    pub fn open_with_ibd(imzml_path: impl AsRef<Path>, ibd_path: impl AsRef<Path>) -> Result<Self> {
        Self::open_with_limits(imzml_path, ibd_path, ImzMLReadLimits::default())
    }

    /// Open both files with explicit ceilings.
    ///
    /// The UUID header is not compared here. The source reports a mismatch as a
    /// warning and loads the dataset anyway, deliberately, so that a
    /// non-conformant `.ibd` still opens; ask [`uuid_status`](Self::uuid_status)
    /// for the verdict.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when either file cannot be opened or read;
    /// [`Error::Parse`] when the `.imzML` is not well formed, when an IMS value
    /// does not parse, or when a spectrum declares no `IMS:1000050` /
    /// `IMS:1000051` pixel coordinate; [`Error::InvalidValue`] when the file
    /// exceeds one of the `limits`.
    pub fn open_with_limits(
        imzml_path: impl AsRef<Path>,
        ibd_path: impl AsRef<Path>,
        limits: ImzMLReadLimits,
    ) -> Result<Self> {
        let imzml_path = imzml_path.as_ref();
        let ibd = ImzMLBinaryIO::open_with_limits(ibd_path, limits)?;
        let file = BufReader::new(File::open(imzml_path)?);
        let mut index = read_index_with_limits(file, &limits)?;
        index.meta.ibd_file_path = ibd.path().to_path_buf();
        Ok(Self {
            imzml_path: imzml_path.to_path_buf(),
            index,
            ibd,
        })
    }

    /// The `.imzML` path this handler parsed.
    pub fn imzml_path(&self) -> &Path {
        &self.imzml_path
    }

    /// The `.ibd` reader, for a caller that needs a raw range.
    pub fn ibd(&mut self) -> &mut ImzMLBinaryIO {
        &mut self.ibd
    }

    /// The `.ibd` path this handler opened, as source
    /// `ImzMLMeta::ibd_file_path`.
    pub fn ibd_path(&self) -> &Path {
        self.ibd.path()
    }

    /// Dataset-level imaging metadata, as source `getImzMLMeta()`.
    pub fn meta(&self) -> &ImzMLMeta {
        &self.index.meta
    }

    /// The whole per-spectrum index, as source `getIndex()`.
    pub fn index(&self) -> &[ImzMLSpectrumIndex] {
        &self.index.spectra
    }

    /// The parsed index and metadata together.
    pub fn parsed(&self) -> &ImzMLIndex {
        &self.index
    }

    /// The index entry for spectrum `index`, as source
    /// `OnDiscImzMLExperiment::getIndex(i)`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `index` is out of range, matching the
    /// source's `Exception::IndexOverflow`.
    pub fn entry(&self, index: usize) -> Result<&ImzMLSpectrumIndex> {
        self.index
            .get(index)
            .ok_or_else(|| out_of_range(index, self.len()))
    }

    /// The number of indexed spectra, as source `getNrSpectra()`.
    pub fn len(&self) -> usize {
        self.index.len()
    }

    /// Whether the dataset indexed no spectra.
    pub fn is_empty(&self) -> bool {
        self.index.is_empty()
    }

    /// Position of the spectrum at 1-based pixel `(x, y, z)`, or `None`.
    ///
    /// See [`ImzMLIndex::index_at_coord`] for the duplicate-coordinate rule.
    pub fn index_at_coord(&self, x: u32, y: u32, z: u32) -> Option<usize> {
        self.index.index_at_coord(x, y, z)
    }

    /// The ceilings this handler enforces.
    pub fn limits(&self) -> ImzMLReadLimits {
        self.ibd.limits()
    }

    /// Compare the `.ibd` UUID header with the XML's `IMS:1000080`.
    ///
    /// Source `verifyIbdUuid_`, which runs this check after the parse on every
    /// read path and logs a warning for each non-matching outcome. Returning
    /// the verdict leaves the policy to the caller.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the header cannot be read.
    pub fn uuid_status(&mut self) -> Result<UuidStatus> {
        let Some(declared) = uuid_bytes(&self.index.meta.uuid) else {
            return Ok(UuidStatus::NotDeclared);
        };
        let Some(found) = self.ibd.uuid()? else {
            return Ok(UuidStatus::IbdTooShort);
        };
        if found == declared {
            Ok(UuidStatus::Match)
        } else {
            Ok(UuidStatus::Mismatch {
                found: hex(&found),
                declared: hex(&declared),
            })
        }
    }

    /// Recompute the `.ibd` SHA-1 and compare it with `IMS:1000091`.
    ///
    /// The source parses `IMS:1000091` into `ImzMLMeta::ibd_sha1` and mirrors it
    /// onto the loaded experiment, but never recomputes it; no OpenMS read path
    /// verifies either declared checksum. This port can verify SHA-1 because
    /// the crate already depends on a SHA-1 implementation for indexed mzML.
    /// `IMS:1000090` MD5 is parsed only: there is no MD5 implementation in this
    /// crate's dependency set and adding one is out of this package's scope.
    ///
    /// Verification is never automatic — hashing the `.ibd` costs a full pass
    /// over a file that can be tens of gigabytes.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the `.ibd` exceeds
    /// [`ImzMLReadLimits::max_checksum_bytes`]; [`Error::Io`] when the read
    /// fails.
    pub fn verify_ibd_sha1(&mut self) -> Result<ChecksumStatus> {
        if self.index.meta.ibd_sha1.is_empty() {
            return Ok(ChecksumStatus::NotDeclared);
        }
        let found = self.ibd.sha1_hex()?;
        if found.eq_ignore_ascii_case(&self.index.meta.ibd_sha1) {
            Ok(ChecksumStatus::Match)
        } else {
            Ok(ChecksumStatus::Mismatch {
                found,
                declared: self.index.meta.ibd_sha1.clone(),
            })
        }
    }

    /// The m/z array of spectrum `index`.
    ///
    /// The values come from the `.ibd` at [`ImzMLSpectrumIndex::mz_offset`]
    /// when the array declares `IMS:1000101`, and from its inline base64
    /// otherwise. That is the one rule [`spectrum`](Self::spectrum) applies,
    /// and applying it here too is what makes the two agree: an earlier
    /// revision read the `.ibd` unconditionally, so a non-conformant file with
    /// an inline m/z array *and* an `IMS:1000102` gave this accessor — and
    /// through it `OnDiscImzMLExperiment::extract_ion_image` — peaks that
    /// `spectrum` did not return for the same pixel.
    ///
    /// In continuous mode every spectrum names the same offset, so this returns
    /// the shared axis and re-reads it per call; the source's on-disc reader
    /// does the same. A caller that walks a whole image should read it once for
    /// pixel 0 and reuse it, which the equality of
    /// [`ImzMLSpectrumIndex::mz_offset`] across the index makes checkable.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `index` is out of range;
    /// [`Error::Unsupported`] when the array is compressed or has no supported
    /// data type; [`Error::Parse`] when an inline payload is not base64 or does
    /// not hold a whole number of elements; otherwise as
    /// [`ImzMLBinaryIO::read_mz_array`].
    pub fn mz_array(&mut self, index: usize) -> Result<Vec<f64>> {
        let Self {
            index: parsed, ibd, ..
        } = self;
        let entry = parsed
            .get(index)
            .ok_or_else(|| out_of_range(index, parsed.len()))?;
        if entry.mz_compressed {
            return Err(compressed_error("m/z array"));
        }
        if entry.mz_external {
            return ibd.read_mz_array(entry.mz_offset, entry.mz_length, entry.mz_type);
        }
        inline_mz_array(&entry.mz_inline, entry.mz_type, &ibd.limits())
    }

    /// The intensity array of spectrum `index`, by the same rule as
    /// [`mz_array`](Self::mz_array).
    ///
    /// # Errors
    ///
    /// As [`mz_array`](Self::mz_array).
    pub fn intensity_array(&mut self, index: usize) -> Result<Vec<f32>> {
        let Self {
            index: parsed, ibd, ..
        } = self;
        let entry = parsed
            .get(index)
            .ok_or_else(|| out_of_range(index, parsed.len()))?;
        if entry.int_compressed {
            return Err(compressed_error("intensity array"));
        }
        if entry.int_external {
            return ibd.read_intensity_array(entry.int_offset, entry.int_length, entry.int_type);
        }
        inline_intensity_array(&entry.int_inline, entry.int_type, &ibd.limits())
    }

    /// Read auxiliary array `aux` of spectrum `index` from the `.ibd`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when either index is out of range;
    /// [`Error::Unsupported`] when the array is compressed, matching the
    /// source's refusal to inflate an external array; otherwise as
    /// [`ImzMLBinaryIO::read_aux_array`].
    pub fn aux_array(&mut self, index: usize, aux: usize) -> Result<Vec<f32>> {
        let entry = self.entry(index)?;
        let array = entry.aux.get(aux).ok_or_else(|| {
            Error::InvalidValue(format!(
                "auxiliary array {aux} is not below the {} arrays of spectrum {index}",
                entry.aux.len()
            ))
        })?;
        let (offset, length, data_type, compressed, name) = (
            array.offset,
            array.length,
            array.data_type,
            array.compressed,
            array.name.clone(),
        );
        if compressed {
            return Err(compressed_error(&format!("auxiliary array '{name}'")));
        }
        self.ibd.read_aux_array(offset, length, data_type, &name)
    }

    /// Decode spectrum `index`: its peaks, its pixel coordinates and its
    /// auxiliary arrays.
    ///
    /// The returned spectrum carries `imzml:x`, `imzml:y` and `imzml:z` meta
    /// values and one float data array per decoded auxiliary array, named after
    /// the array's ontology term and carrying its `unit_accession` where the
    /// XML declared one. That is the contract source
    /// `OnDiscImzMLExperiment::getSpectrum` fulfils; mzML scan metadata is not
    /// loaded on this path, in the source either.
    ///
    /// Peaks are returned in stored order. The source re-sorts them when
    /// `PeakFileOptions::getSortSpectraByMZ()` is set, because the base class's
    /// sort ran before the external arrays overwrote the peaks; this handler
    /// holds no `PeakFileOptions`, so a caller that needs the guarantee calls
    /// [`MSSpectrum::sort_by_position`].
    ///
    /// Auxiliary arrays the source warns about and drops are reported in
    /// [`DecodedSpectrum::skipped_aux`] instead of being logged. A compressed
    /// array is the one auxiliary condition that is an error rather than a skip,
    /// as in the source.
    ///
    /// Each of the two peak arrays is read from the `.ibd` when it declares
    /// `IMS:1000101` and decoded from its own inline base64 when it does not,
    /// which is the rule source `ImzMLInterceptConsumer::consumeSpectrum`
    /// applies (`ImzMLHandler.cpp:207-232`, where the non-external side is
    /// filled from the peaks `MzMLHandler` decoded). A spectrum with exactly
    /// one external array therefore decodes rather than failing the length
    /// check; a conformant imzML 1.1.0 never has one, so reaching that path
    /// takes a non-conformant file. The length check remains, as the guard
    /// against a file whose two arrays genuinely disagree.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `index` is out of range, or when an inline
    /// payload exceeds [`ImzMLReadLimits::max_array_bytes`] or
    /// [`ImzMLReadLimits::max_array_elements`];
    /// [`Error::Unsupported`] when the m/z, intensity or any auxiliary array is
    /// compressed, where the source raises `Exception::ParseError` with the
    /// advice to re-export without compression, or when a non-external peak
    /// array has no supported binary data type; [`Error::Parse`] when the
    /// decoded m/z and intensity arrays have different lengths, naming the
    /// pixel as the source's message does, or when an inline payload is not
    /// valid base64 or does not hold a whole number of elements; otherwise as
    /// [`ImzMLBinaryIO::read_mz_array`].
    pub fn spectrum(&mut self, index: usize) -> Result<DecodedSpectrum> {
        let Self {
            index: parsed, ibd, ..
        } = self;
        let entry = parsed
            .get(index)
            .ok_or_else(|| out_of_range(index, parsed.len()))?;
        if entry.mz_compressed || entry.int_compressed {
            return Err(compressed_error("m/z or intensity arrays"));
        }
        let limits = ibd.limits();
        // Source `ImzMLInterceptConsumer::consumeSpectrum` reads the external
        // side out of the .ibd and fills the other side from the inline peaks
        // its `MzMLHandler` base decoded (ImzMLHandler.cpp:207-232). The two
        // therefore always end up the same length for a well-formed file, which
        // is what makes the mismatch below a corruption guard.
        let mz = if entry.mz_external {
            ibd.read_mz_array(entry.mz_offset, entry.mz_length, entry.mz_type)?
        } else {
            inline_mz_array(&entry.mz_inline, entry.mz_type, &limits)?
        };
        let intensity = if entry.int_external {
            ibd.read_intensity_array(entry.int_offset, entry.int_length, entry.int_type)?
        } else {
            inline_intensity_array(&entry.int_inline, entry.int_type, &limits)?
        };
        if mz.len() != intensity.len() {
            return Err(parse(format!(
                "m/z and intensity array length mismatch at pixel ({},{},{}): mz={} intensity={}",
                entry.x,
                entry.y,
                entry.z,
                mz.len(),
                intensity.len()
            )));
        }

        let mut spectrum = MSSpectrum {
            native_id: entry.native_id.clone(),
            ..MSSpectrum::default()
        };
        spectrum.peaks = mz
            .iter()
            .zip(&intensity)
            .map(|(&mz, &intensity)| Peak1D::new(mz, intensity))
            .collect();
        spectrum
            .metadata
            .insert("imzml:x".into(), MetaValue::from(entry.x));
        spectrum
            .metadata
            .insert("imzml:y".into(), MetaValue::from(entry.y));
        spectrum
            .metadata
            .insert("imzml:z".into(), MetaValue::from(entry.z));

        let peaks = spectrum.peaks.len();
        let mut skipped_aux = Vec::new();
        for _ in 0..entry.unnamed_aux {
            skipped_aux.push(SkippedAux {
                name: String::new(),
                reason: AuxSkipReason::Unnamed,
            });
        }
        for array in &entry.aux {
            if array.compressed {
                return Err(compressed_error(&format!(
                    "auxiliary array '{}' (IMS:1000104 encoded length={})",
                    array.name, array.encoded_bytes
                )));
            }
            if array.length == 0 {
                skipped_aux.push(skipped(array, AuxSkipReason::ZeroLength));
                continue;
            }
            if array.length != peaks as u64 {
                skipped_aux.push(skipped(
                    array,
                    AuxSkipReason::LengthMismatch {
                        length: array.length,
                        peaks,
                    },
                ));
                continue;
            }
            if array.data_type == ImzMLDataType::Unknown {
                skipped_aux.push(skipped(array, AuxSkipReason::UnknownDataType));
                continue;
            }
            let values =
                ibd.read_aux_array(array.offset, array.length, array.data_type, &array.name)?;
            let mut decoded = DataArray::new(array.name.clone(), values);
            if !array.unit_accession.is_empty() {
                decoded.metadata.insert(
                    "unit_accession".into(),
                    MetaValue::from(array.unit_accession.as_str()),
                );
            }
            spectrum.float_data_arrays.push(decoded);
        }
        for name in &entry.inline_aux_names {
            skipped_aux.push(SkippedAux {
                name: name.clone(),
                reason: AuxSkipReason::Inline,
            });
        }

        Ok(DecodedSpectrum {
            spectrum,
            inline_peaks: !entry.mz_external || !entry.int_external,
            skipped_aux,
        })
    }

    /// Decode the spectrum acquired at 1-based pixel `(x, y, z)`.
    ///
    /// Source `OnDiscImzMLExperiment::getSpectrumAtCoord(x, y, z = 1)`, whose
    /// own documentation notes that only the `z == 1` plane is addressable
    /// because the lookup goes through a 2-D geometry. This searches the index
    /// itself, so every `z` present in the file is addressable.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no spectrum was acquired at that pixel,
    /// matching the source's `Exception::ElementNotFound`; otherwise as
    /// [`spectrum`](Self::spectrum).
    pub fn spectrum_at_coord(&mut self, x: u32, y: u32, z: u32) -> Result<DecodedSpectrum> {
        let index = self.index_at_coord(x, y, z).ok_or_else(|| {
            Error::InvalidValue(format!(
                "no imzML spectrum was acquired at pixel ({x},{y},{z})"
            ))
        })?;
        self.spectrum(index)
    }
}

/// Parse the imaging metadata and the per-spectrum index out of an `.imzML`
/// document, with default limits.
///
/// # Errors
///
/// As [`read_index_with_limits`].
pub fn read_index(reader: impl BufRead) -> Result<ImzMLIndex> {
    read_index_with_limits(reader, &ImzMLReadLimits::default())
}

/// Parse the imaging metadata and the per-spectrum index out of an `.imzML`
/// document, with explicit ceilings.
///
/// This is the whole of source `ImzMLHandler`'s own parsing work: the IMS
/// vocabulary and nothing else. mzML metadata belongs to the
/// [`mzml`](crate::format::mzml) reader, as it belongs to `MzMLHandler` in the
/// source.
///
/// A `cvParam` whose accession no imzML rule can act on is skipped without its
/// `value` being decoded at all. The source reaches the same conclusion one
/// step later — its dispatch ignores every unlisted accession, and explicitly
/// so inside a `binaryDataArray`, "so they cannot overwrite the array name" —
/// but it has already transcoded the value by then. Not decoding it is why an
/// `.imzML` that declares a non-UTF-8 encoding, as the upstream processed
/// fixture declares ISO-8859-1, still indexes: its non-ASCII bytes are all in
/// contact and instrument params. A non-UTF-8 value on a param this parser does
/// read is an error, not a silent replacement, because the source relies on
/// Xerces to transcode from the declared encoding and this port implements no
/// transcoder.
///
/// The one binary payload this scan keeps is the inline base64 of an m/z or
/// intensity array that declares no `IMS:1000101`, which the source takes from
/// its `MzMLHandler` base instead. It is kept as encoded ASCII, whitespace
/// removed, charged against [`ImzMLReadLimits::max_text_bytes`] before it is
/// appended, and decoded only when a peak array is asked for. Every other
/// `<binary>` is ignored, as before.
///
/// The reader is capped at one byte past [`ImzMLReadLimits::max_xml_bytes`], so
/// that ceiling bounds the parser's peak buffer and not only its cumulative
/// progress.
///
/// # Errors
///
/// [`Error::Parse`] when the document is not well formed, when an IMS numeric
/// value is empty, negative, non-numeric, non-finite or out of range, when a
/// spectrum carries neither `IMS:1000050` nor `IMS:1000051`, or when a kept
/// inline payload contains non-ASCII bytes; the source raises
/// `Exception::ParseError` for all of these, with the same reasons.
/// [`Error::InvalidValue`] when the document exceeds one of `limits`.
/// [`Error::Io`] when the reader fails.
pub fn read_index_with_limits(
    reader: impl BufRead,
    limits: &ImzMLReadLimits,
) -> Result<ImzMLIndex> {
    let mut parser = Parser::new(limits);
    // `buffer_position()` is cumulative progress, so testing it after the read
    // bounds how far the scan may get but not what one event may allocate:
    // `read_event_into` grows `buffer` to hold a whole start tag or text node
    // before it returns, so a single oversized event was fully committed before
    // the ceiling could fire. Capping the input at one byte past the ceiling
    // caps the buffer with it, because the parser cannot buffer bytes it was
    // never handed. This is the pattern `mzml::read` already uses.
    let limit = limits
        .max_xml_bytes
        .checked_add(1)
        .ok_or_else(|| Error::InvalidValue("imzML XML byte limit must be below u64::MAX".into()))?;
    let mut reader = Reader::from_reader(reader.take(limit));
    reader.config_mut().expand_empty_elements = true;
    reader.config_mut().enable_all_checks(true);
    let mut buffer = Vec::new();
    loop {
        let event = match reader.read_event_into(&mut buffer) {
            Ok(event) => event,
            // The cap truncates an oversized document mid-markup, which the
            // parser reports as a syntax error; the ceiling is the real reason,
            // so say so rather than blaming the document's shape.
            Err(e) => {
                if reader.buffer_position() > limits.max_xml_bytes {
                    return Err(Error::InvalidValue(
                        "imzML XML exceeds the configured byte limit".into(),
                    ));
                }
                return Err(parse(e.to_string()));
            }
        };
        if reader.buffer_position() > limits.max_xml_bytes {
            return Err(Error::InvalidValue(
                "imzML XML exceeds the configured byte limit".into(),
            ));
        }
        match event {
            Event::Start(element) => {
                let name = local_name(element.name().as_ref()).to_vec();
                parser.start(&name, &element)?;
            }
            Event::End(element) => {
                let name = local_name(element.name().as_ref()).to_vec();
                parser.end(&name)?;
            }
            Event::Text(text) => parser.text(text.as_ref())?,
            Event::Eof => break,
            _ => {}
        }
    }
    parser.finish()
}

/// One `cvParam` captured inside a `referenceableParamGroup` for later replay.
///
/// Source `ImzMLHandler::CvEntry`.
#[derive(Clone, Debug)]
struct CvEntry {
    accession: String,
    value: String,
    unit_accession: String,
}

/// What an accession means to an imzML parser, resolved once per accession.
#[derive(Clone, Debug, Eq, PartialEq)]
enum Role {
    PositionX,
    PositionY,
    PositionZ,
    DataType(ImzMLDataType),
    MzArray,
    IntensityArray,
    External,
    Offset,
    ArrayLength,
    EncodedLength,
    NoCompression,
    MaxCountX,
    MaxCountY,
    Mode(ImagingMode),
    Sha1,
    Md5,
    Uuid,
    PixelSizeX,
    PixelSizeY,
    MaxDimX,
    MaxDimY,
    Polarity(&'static str),
    ScanPattern(&'static str),
    ScanDirection(&'static str),
    LineScanDirection(&'static str),
    /// `MS:1000786`, whose `value` is the array's free-text name.
    FreeTextArray,
    /// A child of `MS:1000572` other than `MS:1000576`.
    Compression,
    /// A child of `MS:1000513`, whose CV term name names the array.
    ArrayIdentity(String),
}

/// Per-`binaryDataArray` IMS state, as source `ImzMLHandler::ArrayMeta`.
#[derive(Clone, Debug, Default)]
struct ArrayMeta {
    data_type: ImzMLDataType,
    is_mz: bool,
    is_int: bool,
    is_external: bool,
    offset: u64,
    count: u64,
    encoded_bytes: u64,
    compressed: bool,
    accession: String,
    name: String,
    unit_accession: String,
    /// The array's inline base64 payload, whitespace already removed. Captured
    /// only for a non-external m/z or intensity array, which is the only array
    /// whose inline payload a decode can use.
    inline: String,
}

struct Parser<'a> {
    limits: &'a ImzMLReadLimits,
    meta: ImzMLMeta,
    spectra: Vec<ImzMLSpectrumIndex>,
    in_spectrum: bool,
    in_scan: bool,
    in_bda: bool,
    in_binary: bool,
    binary: String,
    in_ref_group: bool,
    ref_id: String,
    ref_groups: BTreeMap<String, Vec<CvEntry>>,
    group_params: usize,
    native_id: String,
    x: u32,
    y: u32,
    z: u32,
    x_seen: bool,
    y_seen: bool,
    array: ArrayMeta,
    mz: ArrayMeta,
    int: ArrayMeta,
    aux: Vec<ArrayMeta>,
    inline_aux: Vec<String>,
    total_aux: usize,
    text_bytes: usize,
    roles: BTreeMap<String, Option<Role>>,
}

impl<'a> Parser<'a> {
    fn new(limits: &'a ImzMLReadLimits) -> Self {
        Self {
            limits,
            meta: ImzMLMeta::default(),
            spectra: Vec::new(),
            in_spectrum: false,
            in_scan: false,
            in_bda: false,
            in_binary: false,
            binary: String::new(),
            in_ref_group: false,
            ref_id: String::new(),
            ref_groups: BTreeMap::new(),
            group_params: 0,
            native_id: String::new(),
            x: 0,
            y: 0,
            z: 1,
            x_seen: false,
            y_seen: false,
            array: ArrayMeta::default(),
            mz: ArrayMeta::default(),
            int: ArrayMeta::default(),
            aux: Vec::new(),
            inline_aux: Vec::new(),
            total_aux: 0,
            text_bytes: 0,
            roles: BTreeMap::new(),
        }
    }

    /// Charge stored string bytes against the text budget.
    fn charge(&mut self, bytes: usize) -> Result<()> {
        self.text_bytes = self.text_bytes.saturating_add(bytes);
        if self.text_bytes > self.limits.max_text_bytes {
            return Err(Error::InvalidValue(
                "imzML index text exceeds the configured byte limit".into(),
            ));
        }
        Ok(())
    }

    fn start(&mut self, name: &[u8], element: &quick_xml::events::BytesStart<'_>) -> Result<()> {
        // The source tracks its IMS state before delegating to MzMLHandler, and
        // resets the per-spectrum fields on the spectrum start tag.
        match name {
            b"spectrum" => {
                self.in_spectrum = true;
                self.x = 0;
                self.y = 0;
                self.z = 1;
                self.x_seen = false;
                self.y_seen = false;
                self.mz = ArrayMeta::default();
                self.int = ArrayMeta::default();
                self.aux.clear();
                self.inline_aux.clear();
                self.native_id = attribute(element, b"id")?.unwrap_or_default();
                let bytes = self.native_id.len();
                self.charge(bytes)?;
            }
            b"scan" if self.in_spectrum => self.in_scan = true,
            b"binaryDataArray" if self.in_spectrum => {
                self.in_bda = true;
                self.array = ArrayMeta::default();
            }
            b"binary" if self.in_bda => {
                self.in_binary = true;
                self.binary.clear();
            }
            b"referenceableParamGroup" => {
                self.in_ref_group = true;
                self.ref_id = attribute(element, b"id")?.unwrap_or_default();
                if self.ref_groups.len() >= self.limits.max_param_groups
                    && !self.ref_groups.contains_key(&self.ref_id)
                {
                    return Err(Error::InvalidValue(
                        "imzML referenceable parameter groups exceed the configured limit".into(),
                    ));
                }
                let bytes = self.ref_id.len();
                self.charge(bytes)?;
                let id = self.ref_id.clone();
                self.ref_groups.entry(id).or_default();
            }
            b"referenceableParamGroupRef" => {
                if let Some(id) = attribute(element, b"ref")? {
                    self.apply_ref_group(&id)?;
                }
            }
            b"cvParam" => self.cv_param(element)?,
            _ => {}
        }
        Ok(())
    }

    /// Accumulate one `<binary>` text chunk of a non-external peak array.
    ///
    /// `raw` is the event's bytes, undecoded: base64 is ASCII in every charset
    /// an `.imzML` can declare, so no transcoder is needed and the module's
    /// tolerance of an ISO-8859-1 document is preserved. ASCII whitespace is
    /// removed, as source `MzMLHandlerHelper::decodeBase64Arrays` removes it
    /// ("line breaks inside the base64 data are unfortunately no exception"),
    /// and anything else non-ASCII is an error rather than a silent
    /// replacement.
    ///
    /// Only a non-external m/z or intensity array is captured. The mzML 1.1
    /// schema makes `<binary>` the last child of `<binaryDataArray>`, after
    /// every `cvParam`, so the array's identity and `IMS:1000101` are already
    /// known here; an out-of-order document leaves the payload uncaptured and
    /// decodes as a zero-length array.
    ///
    /// The retained length is charged against
    /// [`ImzMLReadLimits::max_text_bytes`] *before* the bytes are appended, so
    /// the ceiling bounds the allocation rather than discovering it afterwards.
    fn text(&mut self, raw: &[u8]) -> Result<()> {
        if !self.in_binary || self.array.is_external || !(self.array.is_mz || self.array.is_int) {
            return Ok(());
        }
        let mut kept = 0usize;
        for &byte in raw {
            if byte.is_ascii_whitespace() {
                continue;
            }
            if !byte.is_ascii() {
                return Err(parse(
                    "non-ASCII bytes in the inline base64 of an imzML peak array",
                ));
            }
            kept = kept.saturating_add(1);
        }
        if kept == 0 {
            return Ok(());
        }
        self.charge(kept)?;
        self.binary
            .try_reserve(kept)
            .map_err(|_| Error::InvalidValue("cannot allocate imzML inline peak array".into()))?;
        self.binary.extend(
            raw.iter()
                .copied()
                .filter(|byte| !byte.is_ascii_whitespace())
                .map(char::from),
        );
        Ok(())
    }

    fn end(&mut self, name: &[u8]) -> Result<()> {
        match name {
            b"binary" if self.in_bda => {
                self.in_binary = false;
                self.array.inline = std::mem::take(&mut self.binary);
            }
            b"binaryDataArray" if self.in_spectrum => {
                self.in_binary = false;
                self.binary.clear();
                let array = std::mem::take(&mut self.array);
                if array.is_mz {
                    self.mz = array;
                } else if array.is_int {
                    self.int = array;
                } else if array.is_external {
                    if self.aux.len() >= self.limits.max_aux_arrays {
                        return Err(Error::InvalidValue(
                            "imzML auxiliary arrays on one spectrum exceed the configured limit"
                                .into(),
                        ));
                    }
                    self.total_aux += 1;
                    if self.total_aux > self.limits.max_total_aux_arrays {
                        return Err(Error::InvalidValue(
                            "imzML auxiliary arrays exceed the configured limit".into(),
                        ));
                    }
                    self.aux.push(array);
                } else {
                    // Source: a non-external auxiliary array is recorded by name
                    // (or accession when unnamed) and later warned about.
                    let name = if array.name.is_empty() {
                        array.accession
                    } else {
                        array.name
                    };
                    self.charge(name.len())?;
                    self.inline_aux.push(name);
                }
                self.in_bda = false;
            }
            b"scan" => self.in_scan = false,
            b"spectrum" => {
                if !self.x_seen || !self.y_seen {
                    return Err(parse(
                        "imzML spectrum missing required IMS pixel coordinate (x and/or y)",
                    ));
                }
                if self.spectra.len() >= self.limits.max_spectra {
                    return Err(Error::InvalidValue(
                        "imzML spectra exceed the configured limit".into(),
                    ));
                }
                self.push_spectrum()?;
                self.in_spectrum = false;
                self.x = 0;
                self.y = 0;
                self.z = 1;
                self.mz = ArrayMeta::default();
                self.int = ArrayMeta::default();
                self.aux.clear();
                self.inline_aux.clear();
            }
            b"referenceableParamGroup" => {
                self.in_ref_group = false;
                self.ref_id.clear();
            }
            _ => {}
        }
        Ok(())
    }

    fn push_spectrum(&mut self) -> Result<()> {
        let index = u32::try_from(self.spectra.len())
            .map_err(|_| Error::InvalidValue("imzML spectrum count exceeds u32".into()))?;
        let mz = std::mem::take(&mut self.mz);
        let int = std::mem::take(&mut self.int);
        // Source records the first decoded data type at dataset level. It does
        // so inside the decode branch, so an index-only load leaves both empty;
        // taking the first declared type instead keeps a metadata-only parse
        // informative and agrees on every file whose first spectrum declares a
        // type, which includes both upstream fixtures.
        if self.meta.mz_data_type == ImzMLDataType::Unknown {
            self.meta.mz_data_type = mz.data_type;
        }
        if self.meta.int_data_type == ImzMLDataType::Unknown {
            self.meta.int_data_type = int.data_type;
        }
        let aux = std::mem::take(&mut self.aux);
        let mut entries = Vec::new();
        entries
            .try_reserve_exact(aux.len())
            .map_err(|_| Error::InvalidValue("cannot allocate imzML auxiliary index".into()))?;
        let mut unnamed_aux = 0u32;
        for array in aux {
            // Source drops an unnamed auxiliary array from the index entry but
            // keeps it in the per-spectrum snapshot, where it is warned about.
            if array.name.is_empty() {
                unnamed_aux = unnamed_aux.saturating_add(1);
                continue;
            }
            self.charge(array.name.len() + array.accession.len() + array.unit_accession.len())?;
            entries.push(ImzMLAuxArray {
                name: array.name,
                accession: array.accession,
                unit_accession: array.unit_accession,
                offset: array.offset,
                length: array.count,
                encoded_bytes: array.encoded_bytes,
                data_type: array.data_type,
                compressed: array.compressed,
            });
        }
        self.spectra.push(ImzMLSpectrumIndex {
            index,
            native_id: std::mem::take(&mut self.native_id),
            x: self.x,
            y: self.y,
            z: self.z,
            mz_offset: mz.offset,
            mz_length: mz.count,
            mz_encoded_bytes: mz.encoded_bytes,
            mz_type: mz.data_type,
            mz_compressed: mz.compressed,
            mz_external: mz.is_external,
            mz_inline: mz.inline,
            int_offset: int.offset,
            int_length: int.count,
            int_encoded_bytes: int.encoded_bytes,
            int_type: int.data_type,
            int_compressed: int.compressed,
            int_external: int.is_external,
            int_inline: int.inline,
            aux: entries,
            unnamed_aux,
            inline_aux_names: std::mem::take(&mut self.inline_aux),
        });
        Ok(())
    }

    /// Replay a captured group, as source `applyRefGroup_`. An unknown
    /// reference is ignored, exactly as the source's failed map lookup is; the
    /// group must therefore already have been declared, which imzML's element
    /// order guarantees.
    fn apply_ref_group(&mut self, id: &str) -> Result<()> {
        let Some(entries) = self.ref_groups.get(id) else {
            return Ok(());
        };
        for entry in entries.clone() {
            let Some(role) = self.role(&entry.accession)? else {
                continue;
            };
            self.handle(&role, &entry.accession, &entry.value, &entry.unit_accession)?;
        }
        Ok(())
    }

    fn cv_param(&mut self, element: &quick_xml::events::BytesStart<'_>) -> Result<()> {
        let Some(accession) = raw_attribute(element, b"accession") else {
            return Ok(());
        };
        let Ok(accession) = std::str::from_utf8(&accession) else {
            // An accession is ASCII in every controlled vocabulary; one that is
            // not cannot match any imzML rule, so it is ignored rather than
            // failing a document this parser is otherwise able to index.
            return Ok(());
        };
        let accession = accession.to_owned();
        let Some(role) = self.role(&accession)? else {
            return Ok(());
        };
        let value = attribute(element, b"value")?.unwrap_or_default();
        let unit = attribute(element, b"unitAccession")?.unwrap_or_default();
        if self.in_ref_group {
            self.group_params += 1;
            if self.group_params > self.limits.max_group_params {
                return Err(Error::InvalidValue(
                    "imzML referenceable parameters exceed the configured limit".into(),
                ));
            }
            self.charge(accession.len() + value.len() + unit.len())?;
            let entry = CvEntry {
                accession,
                value,
                unit_accession: unit,
            };
            let id = self.ref_id.clone();
            self.ref_groups.entry(id).or_default().push(entry);
            return Ok(());
        }
        self.handle(&role, &accession, &value, &unit)
    }

    /// Resolve an accession's meaning, memoised. `None` means no imzML rule can
    /// act on it in any context.
    fn role(&mut self, accession: &str) -> Result<Option<Role>> {
        if let Some(role) = self.roles.get(accession) {
            return Ok(role.clone());
        }
        let role = literal_role(accession);
        let role = match role {
            Some(role) => Some(role),
            None => {
                if self.roles.len() >= self.limits.max_cv_lookups {
                    return Err(Error::InvalidValue(
                        "imzML vocabulary lookups exceed the configured limit".into(),
                    ));
                }
                let cv = ControlledVocabulary::psi_ms()?;
                if !cv.exists(accession) {
                    None
                } else if cv.is_child_of(accession, "MS:1000572")? {
                    // Source: any child of MS:1000572 that is not MS:1000576.
                    Some(Role::Compression)
                } else if cv.is_child_of(accession, "MS:1000513")? {
                    // Source uses the CV term's own name, matching how
                    // MzMLHandler names the float data array it creates, rather
                    // than the XML name attribute.
                    Some(Role::ArrayIdentity(cv.get_term(accession)?.name.clone()))
                } else {
                    None
                }
            }
        };
        self.charge(accession.len())?;
        if let Some(Role::ArrayIdentity(name)) = &role {
            self.charge(name.len())?;
        }
        self.roles.insert(accession.to_owned(), role.clone());
        Ok(role)
    }

    /// Source `handleIMSCvParam_`: three context-ordered dispatch blocks.
    fn handle(&mut self, role: &Role, accession: &str, value: &str, unit: &str) -> Result<()> {
        // 1. Pixel coordinates, only inside a scan inside a spectrum. Only the
        //    three coordinate terms return here; anything else falls through,
        //    which is how a polarity term inside a scan still reaches the
        //    dataset block.
        if self.in_scan && self.in_spectrum {
            match role {
                Role::PositionX => {
                    self.x = ims_u32(accession, value)?;
                    self.x_seen = true;
                    return Ok(());
                }
                Role::PositionY => {
                    self.y = ims_u32(accession, value)?;
                    self.y_seen = true;
                    return Ok(());
                }
                Role::PositionZ => {
                    self.z = ims_u32(accession, value)?;
                    return Ok(());
                }
                _ => {}
            }
        }

        // 2. Inside a binaryDataArray. The source returns unconditionally at
        //    the end of this block, so a dataset-level term inside an array is
        //    ignored and cannot overwrite the array's name.
        if self.in_bda {
            match role {
                Role::DataType(data_type) => self.array.data_type = *data_type,
                Role::MzArray => self.array.is_mz = true,
                Role::IntensityArray => self.array.is_int = true,
                Role::External => self.array.is_external = true,
                Role::Offset => self.array.offset = ims_u64(accession, value)?,
                Role::ArrayLength => self.array.count = ims_u64(accession, value)?,
                Role::EncodedLength => self.array.encoded_bytes = ims_u64(accession, value)?,
                Role::NoCompression => self.array.compressed = false,
                Role::Compression => self.array.compressed = true,
                Role::FreeTextArray => {
                    self.array.accession = accession.to_owned();
                    if !value.is_empty() {
                        self.array.name = value.to_owned();
                    }
                    if !unit.is_empty() {
                        self.array.unit_accession = unit.to_owned();
                    }
                }
                Role::ArrayIdentity(name) => {
                    self.array.accession = accession.to_owned();
                    self.array.name = name.clone();
                    if !unit.is_empty() {
                        self.array.unit_accession = unit.to_owned();
                    }
                }
                _ => {}
            }
            return Ok(());
        }

        // 3. Dataset-level imaging metadata.
        match role {
            Role::MaxCountX => self.meta.max_count_x = ims_u32(accession, value)?,
            Role::MaxCountY => self.meta.max_count_y = ims_u32(accession, value)?,
            Role::Mode(mode) => self.meta.imaging_mode = Some(*mode),
            Role::Sha1 => {
                self.charge(value.len())?;
                self.meta.ibd_sha1 = value.to_owned();
            }
            Role::Md5 => {
                self.charge(value.len())?;
                self.meta.ibd_md5 = value.to_owned();
            }
            Role::Uuid => {
                // Source keeps the previous value for an empty IMS:1000080.
                if !value.is_empty() {
                    self.charge(value.len())?;
                    self.meta.uuid = value.to_owned();
                }
            }
            Role::PixelSizeX => self.meta.pixel_size_x = ims_f64(accession, value)?,
            Role::PixelSizeY => self.meta.pixel_size_y = ims_f64(accession, value)?,
            Role::MaxDimX => self.meta.max_dim_x = ims_f64(accession, value)?,
            Role::MaxDimY => self.meta.max_dim_y = ims_f64(accession, value)?,
            Role::Polarity(text) => self.meta.polarity = (*text).to_owned(),
            Role::ScanPattern(text) => self.meta.scan_pattern = (*text).to_owned(),
            Role::ScanDirection(text) => self.meta.scan_direction = (*text).to_owned(),
            Role::LineScanDirection(text) => self.meta.line_scan_direction = (*text).to_owned(),
            _ => {}
        }
        Ok(())
    }

    fn finish(mut self) -> Result<ImzMLIndex> {
        // Source raises the dataset bounding box to the largest coordinate it
        // delivers, so a file whose IMS:1000042/43 understate the grid still
        // reports the real extent. max_count_z has no CV term at all.
        for entry in &self.spectra {
            self.meta.max_count_x = self.meta.max_count_x.max(entry.x);
            self.meta.max_count_y = self.meta.max_count_y.max(entry.y);
            self.meta.max_count_z = self.meta.max_count_z.max(entry.z);
        }
        Ok(ImzMLIndex {
            meta: self.meta,
            spectra: self.spectra,
        })
    }
}

/// Accessions with a fixed imzML meaning, in the source's dispatch order.
fn literal_role(accession: &str) -> Option<Role> {
    Some(match accession {
        "IMS:1000050" => Role::PositionX,
        "IMS:1000051" => Role::PositionY,
        "IMS:1000052" => Role::PositionZ,
        "MS:1000521" => Role::DataType(ImzMLDataType::Float32),
        "MS:1000523" => Role::DataType(ImzMLDataType::Float64),
        "MS:1000519" => Role::DataType(ImzMLDataType::Int32),
        "MS:1000522" => Role::DataType(ImzMLDataType::Int64),
        "MS:1000514" => Role::MzArray,
        "MS:1000515" => Role::IntensityArray,
        "IMS:1000101" => Role::External,
        "IMS:1000102" => Role::Offset,
        "IMS:1000103" => Role::ArrayLength,
        "IMS:1000104" => Role::EncodedLength,
        "MS:1000576" => Role::NoCompression,
        "MS:1000786" => Role::FreeTextArray,
        "IMS:1000042" => Role::MaxCountX,
        "IMS:1000043" => Role::MaxCountY,
        "IMS:1000030" => Role::Mode(ImagingMode::Continuous),
        "IMS:1000031" => Role::Mode(ImagingMode::Processed),
        "IMS:1000091" => Role::Sha1,
        "IMS:1000090" => Role::Md5,
        "IMS:1000080" => Role::Uuid,
        "IMS:1000046" => Role::PixelSizeX,
        "IMS:1000047" => Role::PixelSizeY,
        "IMS:1000044" => Role::MaxDimX,
        "IMS:1000045" => Role::MaxDimY,
        "MS:1000129" => Role::Polarity("negative"),
        "MS:1000130" => Role::Polarity("positive"),
        "IMS:1000401" => Role::ScanPattern("top down"),
        "IMS:1000402" => Role::ScanPattern("bottom up"),
        "IMS:1000413" => Role::ScanDirection("flyback"),
        "IMS:1000412" => Role::ScanDirection("meander"),
        "IMS:1000480" => Role::ScanDirection("horizontal"),
        "IMS:1000481" => Role::ScanDirection("vertical"),
        "IMS:1000491" => Role::LineScanDirection("left-right"),
        "IMS:1000492" => Role::LineScanDirection("right-left"),
        _ => return None,
    })
}

fn local_name(name: &[u8]) -> &[u8] {
    match name.iter().position(|&byte| byte == b':') {
        Some(colon) => &name[colon + 1..],
        None => name,
    }
}

fn raw_attribute<'b>(
    element: &'b quick_xml::events::BytesStart<'_>,
    key: &[u8],
) -> Option<std::borrow::Cow<'b, [u8]>> {
    element
        .attributes()
        .with_checks(false)
        .flatten()
        .find(|attribute| local_name(attribute.key.as_ref()) == key)
        .map(|attribute| attribute.value)
}

/// The unescaped text of one attribute, or `None` when it is absent.
fn attribute(element: &quick_xml::events::BytesStart<'_>, key: &[u8]) -> Result<Option<String>> {
    let Some(raw) = raw_attribute(element, key) else {
        return Ok(None);
    };
    let text = std::str::from_utf8(&raw).map_err(|_| {
        parse("imzML attribute this parser must read is not valid UTF-8; the declared encoding is not transcoded")
    })?;
    let text = quick_xml::escape::unescape(text).map_err(|e| parse(e.to_string()))?;
    Ok(Some(text.into_owned()))
}

/// Source `parseImsUInt32_`.
fn ims_u32(accession: &str, value: &str) -> Result<u32> {
    digits(accession, value)?
        .parse::<u32>()
        .map_err(|_| ims_error(accession, value, "out of range for uint32"))
}

/// Source `parseImsUInt64_`, which rejects a value containing `-` outright
/// because `std::stoull` would otherwise wrap it to a huge unsigned number.
fn ims_u64(accession: &str, value: &str) -> Result<u64> {
    digits(accession, value)?
        .parse::<u64>()
        .map_err(|_| ims_error(accession, value, "out of range for uint64"))
}

/// The source's numeric helpers accept the leading whitespace `std::stoull` and
/// `std::stod` skip, then reject any trailing character. Requiring digits only
/// keeps that contract and additionally rejects a `+` sign, which no imzML
/// writer emits and no upstream fixture carries.
fn digits<'b>(accession: &str, value: &'b str) -> Result<&'b str> {
    let trimmed = value.trim_matches(|c: char| c.is_ascii_whitespace());
    if trimmed.is_empty() {
        return Err(ims_error(accession, value, "empty value"));
    }
    if trimmed.starts_with('-') {
        return Err(ims_error(accession, value, "negative value not allowed"));
    }
    if !trimmed.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(ims_error(accession, value, "not a valid unsigned integer"));
    }
    Ok(trimmed)
}

/// Source `parseImsDouble_`. The source's `std::stod` accepts `inf` and `nan`;
/// this rejects them, because a non-finite pixel size or image extent cannot be
/// used by any consumer and would propagate silently.
fn ims_f64(accession: &str, value: &str) -> Result<f64> {
    let trimmed = value.trim_matches(|c: char| c.is_ascii_whitespace());
    if trimmed.is_empty() {
        return Err(ims_error(accession, value, "empty value"));
    }
    trimmed
        .parse::<f64>()
        .ok()
        .filter(|number| number.is_finite())
        .ok_or_else(|| ims_error(accession, value, "not a finite floating-point number"))
}

fn ims_error(accession: &str, value: &str, reason: &str) -> Error {
    parse(format!(
        "Invalid IMS CV value for {accession} ('{value}'): {reason}"
    ))
}

fn compressed_error(what: &str) -> Error {
    Error::Unsupported(format!(
        "Compressed external {what} are not supported (need uncompressed MS:1000576). \
         Re-export without compression."
    ))
}

fn out_of_range(index: usize, count: usize) -> Error {
    Error::InvalidValue(format!(
        "imzML spectrum index {index} is not below the indexed count {count}"
    ))
}

fn skipped(array: &ImzMLAuxArray, reason: AuxSkipReason) -> SkippedAux {
    SkippedAux {
        name: array.name.clone(),
        reason,
    }
}

/// Widen `count` little-endian elements of `data_type` out of `raw` into `f64`.
///
/// Shared by the `.ibd` read and the inline-base64 decode so the two cannot
/// widen differently. Integer types are widened exactly as source
/// `ImzMLBinaryIO::readMzArray` widens them, so a 64-bit integer above 2^53
/// loses precision in both.
fn widen_to_f64(
    raw: &[u8],
    data_type: ImzMLDataType,
    count: usize,
    what: &str,
) -> Result<Vec<f64>> {
    let mut out = try_vec::<f64>(count, what)?;
    match data_type {
        ImzMLDataType::Float32 => {
            for chunk in raw.chunks_exact(4) {
                out.push(f64::from(f32::from_le_bytes(four(chunk))));
            }
        }
        ImzMLDataType::Float64 => {
            for chunk in raw.chunks_exact(8) {
                out.push(f64::from_le_bytes(eight(chunk)));
            }
        }
        ImzMLDataType::Int32 => {
            for chunk in raw.chunks_exact(4) {
                out.push(f64::from(i32::from_le_bytes(four(chunk))));
            }
        }
        ImzMLDataType::Int64 => {
            for chunk in raw.chunks_exact(8) {
                out.push(i64::from_le_bytes(eight(chunk)) as f64);
            }
        }
        // Both callers resolve the element width first, which fails for
        // `Unknown`, so this arm is unreachable; it is an error rather than a
        // panic because the input that would reach it is file-derived.
        ImzMLDataType::Unknown => {
            return Err(Error::Unsupported(format!("unsupported {what} data type")));
        }
    }
    Ok(out)
}

/// Narrow `count` little-endian elements of `data_type` out of `raw` into
/// `f32`, as source `readIntArray` / `readAuxArray` narrow them.
fn narrow_to_f32(
    raw: &[u8],
    data_type: ImzMLDataType,
    count: usize,
    what: &str,
) -> Result<Vec<f32>> {
    let mut out = try_vec::<f32>(count, what)?;
    match data_type {
        ImzMLDataType::Float32 => {
            for chunk in raw.chunks_exact(4) {
                out.push(f32::from_le_bytes(four(chunk)));
            }
        }
        ImzMLDataType::Float64 => {
            for chunk in raw.chunks_exact(8) {
                out.push(f64::from_le_bytes(eight(chunk)) as f32);
            }
        }
        ImzMLDataType::Int32 => {
            for chunk in raw.chunks_exact(4) {
                out.push(i32::from_le_bytes(four(chunk)) as f32);
            }
        }
        ImzMLDataType::Int64 => {
            for chunk in raw.chunks_exact(8) {
                out.push(i64::from_le_bytes(eight(chunk)) as f32);
            }
        }
        ImzMLDataType::Unknown => {
            return Err(Error::Unsupported(format!("unsupported {what} data type")));
        }
    }
    Ok(out)
}

/// Decode one array's inline base64 payload into raw little-endian bytes and
/// the element count they hold.
///
/// `Ok(None)` for an empty payload, which is what an external array's
/// placeholder `<binary/>` leaves behind and what an absent `<binary>` leaves
/// behind, so neither is an error — the same tolerance
/// [`ImzMLBinaryIO::read_mz_array`] gives a zero `IMS:1000103`.
///
/// Every ceiling the `.ibd` preflight applies is applied here too, and the
/// decoded-size ceiling is checked against the encoded character count
/// *before* the decode allocates: four base64 characters carry at most three
/// bytes, so `encoded.len() / 4 * 3` is an upper bound on the output.
fn inline_preflight(
    encoded: &str,
    data_type: ImzMLDataType,
    limits: &ImzMLReadLimits,
    what: &str,
) -> Result<Option<(Vec<u8>, usize)>> {
    if encoded.is_empty() {
        return Ok(None);
    }
    let width = data_type.width().ok_or_else(|| {
        Error::Unsupported(format!("unsupported inline {what} data type in .imzML"))
    })?;
    // Checked before the decoder allocates its output buffer.
    if encoded.len() / 4 * 3 > limits.max_array_bytes {
        return Err(Error::InvalidValue(format!(
            "inline {what} exceeds the configured byte limit"
        )));
    }
    let raw = STANDARD
        .decode(encoded.as_bytes())
        .map_err(|e| parse(format!("invalid inline {what} base64: {e}")))?;
    if raw.len() % width != 0 {
        return Err(parse(format!(
            "inline {what} holds {} bytes, which is not a whole number of {width}-byte elements",
            raw.len()
        )));
    }
    let count = raw.len() / width;
    if count as u64 > limits.max_array_elements {
        return Err(Error::InvalidValue(format!(
            "inline {what} element count {count} exceeds the configured limit of {}",
            limits.max_array_elements
        )));
    }
    Ok(Some((raw, count)))
}

/// The m/z values an array's inline base64 payload holds.
fn inline_mz_array(
    encoded: &str,
    data_type: ImzMLDataType,
    limits: &ImzMLReadLimits,
) -> Result<Vec<f64>> {
    let Some((raw, count)) = inline_preflight(encoded, data_type, limits, "m/z array")? else {
        return Ok(Vec::new());
    };
    widen_to_f64(&raw, data_type, count, "m/z array")
}

/// The intensity values an array's inline base64 payload holds.
fn inline_intensity_array(
    encoded: &str,
    data_type: ImzMLDataType,
    limits: &ImzMLReadLimits,
) -> Result<Vec<f32>> {
    let Some((raw, count)) = inline_preflight(encoded, data_type, limits, "intensity array")?
    else {
        return Ok(Vec::new());
    };
    narrow_to_f32(&raw, data_type, count, "intensity array")
}

fn check_write_count(count: usize, limits: &ImzMLReadLimits, what: &str) -> Result<()> {
    if count as u64 > limits.max_array_elements {
        return Err(Error::InvalidValue(format!(
            "{what}: element count {count} exceeds the configured limit of {}",
            limits.max_array_elements
        )));
    }
    Ok(())
}

fn parse(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

fn try_vec<T>(capacity: usize, what: &str) -> Result<Vec<T>> {
    let mut out = Vec::new();
    out.try_reserve_exact(capacity)
        .map_err(|_| Error::InvalidValue(format!("cannot allocate imzML {what}")))?;
    Ok(out)
}

fn four(chunk: &[u8]) -> [u8; 4] {
    [chunk[0], chunk[1], chunk[2], chunk[3]]
}

fn eight(chunk: &[u8]) -> [u8; 8] {
    [
        chunk[0], chunk[1], chunk[2], chunk[3], chunk[4], chunk[5], chunk[6], chunk[7],
    ]
}

fn nibble(digit: u8) -> u8 {
    match digit {
        b'0'..=b'9' => digit - b'0',
        b'a'..=b'f' => digit - b'a' + 10,
        _ => digit - b'A' + 10,
    }
}

/// Source `bytesToHex_`: lower-case, two digits per byte, no separators.
fn hex(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len() * 2);
    for &byte in bytes {
        out.push(char::from_digit(u32::from(byte >> 4), 16).unwrap_or('0'));
        out.push(char::from_digit(u32::from(byte & 0x0f), 16).unwrap_or('0'));
    }
    out
}

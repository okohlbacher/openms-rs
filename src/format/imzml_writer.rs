// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Writer for an imzML 1.1.0 dataset: the `.imzML` XML and its companion `.ibd`.
//!
//! Ports `FORMAT/HANDLERS/ImzMLWriter.h` with
//! `FORMAT/HANDLERS/ImzMLWriter.cpp`. See `docs/IMZML_WRITER_SUPPORT.md`.
//!
//! An imzML dataset is two files. The `.ibd` opens with a 16-byte UUID header
//! and then carries every binary array back to back, uncompressed and
//! little-endian. The `.imzML` is mzML XML carrying the metadata and, per
//! spectrum, the IMS params naming the byte offset (`IMS:1000102`), element
//! count (`IMS:1000103`) and stored byte length (`IMS:1000104`) of that
//! spectrum's arrays inside the `.ibd`; the payload elements themselves are
//! empty `<binary/>` tags. Both storage modes are written:
//! [`ImagingMode::Continuous`](crate::format::imzml_handler::ImagingMode::Continuous)
//! stores one shared m/z array once and points every spectrum's m/z offset at
//! it, and
//! [`ImagingMode::Processed`](crate::format::imzml_handler::ImagingMode::Processed)
//! gives each spectrum its own m/z array.
//!
//! # Relationship to the mzML writer
//!
//! The source writes this XML by hand, with raw `std::ostream` inserters in a
//! file-static `writeImzMLXml_`, rather than delegating to `MzMLHandler`. It
//! has to: every `binaryDataArray` here is an *external* array whose
//! `<binary/>` element is empty and whose payload lives in a second file, which
//! is the one thing an mzML serialiser cannot emit. This port keeps that
//! division. What it reuses is the binary half — the `.ibd` array writers, the
//! UUID codec and the `.ibd` path rule already ported next door in
//! [`imzml_handler`](crate::format::imzml_handler) — plus `quick_xml`'s escape,
//! the same escape [`mzml`](crate::format::mzml) uses. Nothing here encodes
//! base64, because imzML external arrays are not base64.
//!
//! # Bounded and atomic
//!
//! Every offset written into the XML is an absolute byte position in a second
//! file, so a writer that miscounts produces a document that makes a reader
//! read out of bounds. Offsets are therefore accumulated with checked `u64`
//! arithmetic against
//! [`ImzMLWriteLimits`](crate::format::imzml_writer::ImzMLWriteLimits), and the
//! complete plan — pixel coordinates, every offset, every length, every
//! resolved ontology term — is built and checked **before either file is
//! created**. A rejected experiment therefore leaves no file behind at all.
//! The source opens the `.ibd` first and validates as it goes, so a failure
//! part-way through leaves a truncated `.ibd` and no `.imzML` at all.
//!
//! Nothing here starts a thread, and nothing buffers a payload: arrays are
//! streamed to the `.ibd` one at a time.

use crate::concept::progress_logger::ProgressLogger;
use crate::format::PeakFileOptions;
use crate::format::controlled_vocabulary::ControlledVocabulary;
use crate::format::imzml_handler::{
    IBD_UUID_BYTES, ImagingMode, ImzMLDataType, ImzMLMeta, ImzMLReadLimits, infer_ibd_path,
    uuid_bytes, write_float32_array, write_float64_array, write_mz_as_float32, write_mz_as_float64,
};
use crate::kernel::{DataArray, MSExperiment, MSSpectrum, NumericRange};
use crate::metadata::{MetaValue, MetaValueData};
use crate::{Error, Result};
use sha1::{Digest, Sha1};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

/// The absolute m/z agreement the source accepts before it calls a dataset
/// continuous, from `spectraShareMz_`'s default argument.
///
/// Continuous mode stores one m/z array and every spectrum reads it back, so
/// any disagreement below this bound is silently replaced by the reference
/// spectrum's value. [`ImzMLWriteOptions::source`] selects it;
/// [`ImzMLWriteOptions::default`] uses `0.0` and so discards nothing.
pub const SOURCE_SHARED_MZ_TOLERANCE: f64 = 1e-5;

/// Namespace of the RFC 4122 version-5 identifier [`derive_uuid`] computes.
///
/// Fixed for the life of this crate, so a dataset written twice from the same
/// arrays carries the same UUID. It is itself a version-5 identifier,
/// `b29f06ef-b005-5b88-b246-0bc4d30466d8`, obtained by stamping the first 16
/// bytes of `SHA-1("openms-rs:imzml:ibd:v1")`. The source has no analogue: it
/// draws 16 bytes from `std::random_device` instead.
pub const IBD_UUID_NAMESPACE: [u8; IBD_UUID_BYTES] = [
    0xb2, 0x9f, 0x06, 0xef, 0xb0, 0x05, 0x5b, 0x88, 0xb2, 0x46, 0x0b, 0xc4, 0xd3, 0x04, 0x66, 0xd8,
];

/// Explicit ceilings checked before either output file is created.
///
/// The source has no write-side ceiling other than the shared
/// `MAX_IBD_ARRAY_ELEMENTS` its array writers enforce, and it accumulates
/// `.ibd` offsets in a plain `uint64_t` that wraps on overflow. Each field here
/// turns one of those into an error raised during the preflight, before a byte
/// is written.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ImzMLWriteLimits {
    /// Maximum spectra in one dataset.
    pub max_spectra: usize,
    /// Maximum peaks in one spectrum. Also handed to the `.ibd` array writers
    /// as [`ImzMLReadLimits::max_array_elements`], the only field of that
    /// struct they consult.
    pub max_peaks_per_spectrum: usize,
    /// Maximum peaks summed over the whole dataset.
    pub max_total_peaks: u64,
    /// Maximum auxiliary float arrays exported for one spectrum.
    pub max_aux_arrays_per_spectrum: usize,
    /// Maximum auxiliary float arrays exported across the dataset.
    pub max_total_aux_arrays: usize,
    /// Maximum size of the `.ibd` the plan may describe, UUID header included.
    pub max_ibd_bytes: u64,
    /// Cumulative bytes of every string the XML carries over from the
    /// experiment: native identifiers, array names, the UUID, the instrument
    /// model and the geometry vocabulary. Counted per occurrence, so this
    /// bounds the experiment-derived text of the document rather than the set
    /// of distinct strings.
    pub max_text_bytes: usize,
    /// Maximum bytes of one such string.
    pub max_name_bytes: usize,
    /// Maximum distinct array names resolved against the PSI-MS vocabulary.
    /// Results are memoised, so this bounds the vocabulary work of one store.
    pub max_cv_lookups: usize,
    /// Maximum `.ibd` bytes re-read to compute the declared SHA-1 and MD5.
    pub max_checksum_bytes: u64,
    /// Maximum entries kept in each list of [`StoreReport`]. The source caps
    /// its own warning listings at 20 for the same reason: a dataset where
    /// every spectrum is affected must not make the report grow with it.
    pub max_reported_items: usize,
}

impl Default for ImzMLWriteLimits {
    fn default() -> Self {
        Self {
            max_spectra: 5_000_000,
            max_peaks_per_spectrum: 100_000_000,
            max_total_peaks: 4_000_000_000,
            max_aux_arrays_per_spectrum: 256,
            max_total_aux_arrays: 10_000_000,
            max_ibd_bytes: 1 << 40,
            max_text_bytes: 1 << 30,
            max_name_bytes: 64 << 10,
            max_cv_lookups: 100_000,
            max_checksum_bytes: 16 << 30,
            max_reported_items: 20,
        }
    }
}

impl ImzMLWriteLimits {
    /// The ceilings the stage-1 `.ibd` array writers enforce, derived from
    /// [`Self::max_peaks_per_spectrum`]. Only
    /// [`ImzMLReadLimits::max_array_elements`] is read by them, so every other
    /// field of the returned value is its default and is never consulted.
    fn array_limits(self) -> ImzMLReadLimits {
        ImzMLReadLimits {
            max_array_elements: self.max_peaks_per_spectrum as u64,
            ..ImzMLReadLimits::default()
        }
    }
}

/// How much m/z disagreement a continuous dataset may hide, plus the ceilings.
///
/// Continuous mode is the only lossy choice this writer can make: it stores the
/// reference spectrum's m/z array once and every other spectrum reads that
/// array back, so m/z differences within [`Self::shared_mz_tolerance`] are
/// discarded. The default is `0.0`, which discards nothing — a dataset whose
/// m/z axes are not bit-identical is written in processed mode instead, and an
/// explicit `imzml:imaging_mode` of `"continuous"` over such a dataset is
/// refused. [`Self::source`] selects the source's
/// [`SOURCE_SHARED_MZ_TOLERANCE`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImzMLWriteOptions {
    /// Resource ceilings for the preflight.
    pub limits: ImzMLWriteLimits,
    /// Largest absolute m/z difference two spectra may have and still be
    /// treated as sharing one axis. Must be finite and non-negative.
    pub shared_mz_tolerance: f64,
}

impl Default for ImzMLWriteOptions {
    fn default() -> Self {
        Self {
            limits: ImzMLWriteLimits::default(),
            shared_mz_tolerance: 0.0,
        }
    }
}

impl ImzMLWriteOptions {
    /// The source `ImzMLWriter::store` behaviour: accept up to
    /// [`SOURCE_SHARED_MZ_TOLERANCE`] of m/z disagreement when choosing or
    /// validating continuous mode, and discard it.
    pub fn source() -> Self {
        Self {
            limits: ImzMLWriteLimits::default(),
            shared_mz_tolerance: SOURCE_SHARED_MZ_TOLERANCE,
        }
    }
}

/// Why one `FloatDataArray` of a spectrum produced no external array.
///
/// The source logs a warning for each of these and carries on, so that one
/// unexportable auxiliary array cannot abort a whole store. This port keeps the
/// policy and reports the omissions in [`StoreReport::skipped_float_arrays`]
/// rather than only logging them.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum FloatArraySkipReason {
    /// The array holds no values. The source skips it before any other check
    /// and without a warning; it is reported here because the array's name and
    /// existence are still lost.
    Empty,
    /// The array has no name, so no ontology term and no `MS:1000786` value
    /// could identify it on read.
    Unnamed,
    /// The array is not as long as the peak array. Viewers require an auxiliary
    /// array to hold one value per peak.
    LengthMismatch {
        /// Values the array holds.
        length: usize,
        /// Peaks the spectrum holds.
        peaks: usize,
    },
    /// The name resolves to `MS:1000514` "m/z array" or `MS:1000515`
    /// "intensity array". Written out, such an array would become a second m/z
    /// or intensity `binaryDataArray` whose offset and type replace the real
    /// peak metadata on read, because the last `cvParam` wins.
    ReservedPeakArrayName {
        /// The reserved accession the name resolved to.
        accession: String,
    },
}

/// One `FloatDataArray` that the store did not export, with the reason.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SkippedFloatArray {
    /// 0-based position of the spectrum in the filtered experiment.
    pub spectrum: usize,
    /// The array's name, empty for a [`FloatArraySkipReason::Unnamed`] array
    /// since a name is what it lacked.
    pub name: String,
    /// Why it was not exported.
    pub reason: FloatArraySkipReason,
}

/// One spectrum whose pixel was already claimed by an earlier spectrum.
///
/// Source `validatePixelMetadataForStore_` warns about these and writes them
/// out unchanged, because the readers accept duplicated coordinates and a
/// dataset that loads has to be storable again. Readers map only the first
/// spectrum per pixel into the imaging geometry, so the later ones stay
/// reachable by index but not by coordinate.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct DuplicatePixel {
    /// 0-based position of the spectrum in the filtered experiment.
    pub spectrum: usize,
    /// Pixel column, 1-based.
    pub x: u32,
    /// Pixel row, 1-based.
    pub y: u32,
    /// Depth slice, 1-based.
    pub z: u32,
}

/// What one [`store`] wrote, and what it chose not to write.
///
/// The source returns `void` and sends all of this to `OPENMS_LOG_WARN`.
/// Returning it is what lets a caller notice that a store silently dropped a
/// charge array or collapsed nine pixels into one.
#[derive(Clone, Debug, PartialEq)]
pub struct StoreReport {
    /// The mode actually written. Mirrors
    /// [`ImzMLMeta::imaging_mode`](crate::format::imzml_handler::ImzMLMeta::imaging_mode)
    /// of [`Self::meta`], which is always `Some` in a report.
    pub mode: ImagingMode,
    /// The dataset metadata as written: the resolved UUID, the grid raised to
    /// the largest pixel, the recomputed image extents, the binary precisions,
    /// the `.ibd` path and the `.ibd` SHA-1 and MD5 digests.
    pub meta: ImzMLMeta,
    /// Spectra written to the `.imzML`.
    pub spectra_written: usize,
    /// Spectra that [`apply_store_options`] removed before the write.
    pub spectra_filtered_out: usize,
    /// Size of the `.ibd` in bytes, UUID header included.
    pub ibd_bytes: u64,
    /// Auxiliary float arrays written across the dataset.
    pub aux_arrays_written: usize,
    /// How many spectra reused a pixel an earlier spectrum had claimed.
    pub duplicate_pixel_count: usize,
    /// The first [`ImzMLWriteLimits::max_reported_items`] of them.
    pub duplicate_pixels: Vec<DuplicatePixel>,
    /// How many non-empty integer and string data arrays were dropped. imzML
    /// carries per-peak data only as float arrays, so these have no external
    /// representation at all.
    pub dropped_data_array_count: usize,
    /// Their distinct names in first-seen order, at most
    /// [`ImzMLWriteLimits::max_reported_items`] of them. Integer and string
    /// array names share one list, as they share one warning in the source.
    pub dropped_data_array_names: Vec<String>,
    /// How many `FloatDataArray`s were not exported.
    pub skipped_float_array_count: usize,
    /// The first [`ImzMLWriteLimits::max_reported_items`] of them.
    pub skipped_float_arrays: Vec<SkippedFloatArray>,
}

/// The PSI-MS identity of a `FloatDataArray` name, as written to the XML.
///
/// Source file-static `ResolvedArrayCv`.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ResolvedArrayCv {
    /// Accession of the matching child of `MS:1000513`, or `MS:1000786` when
    /// no child carries that name.
    pub accession: String,
    /// The ontology term's own name, or `"non-standard data array"` for
    /// `MS:1000786`.
    pub name: String,
    /// `unitAccession` to write, empty when the term declares no unit or
    /// declares more than one.
    pub unit_accession: String,
    /// `unitName` to write, empty under the same conditions.
    pub unit_name: String,
    /// `unitCvRef` to write: the first two characters of
    /// [`Self::unit_accession`], so `MS:1002814` gives `MS` and `UO:0000110`
    /// gives `UO`. Empty when there is no unit.
    pub unit_cv_ref: String,
    /// True when the name matched no ontology term, so the array is written as
    /// `MS:1000786` with the name as the param's `value`.
    pub non_standard: bool,
}

/// Map a `FloatDataArray` name to a PSI-MS binary-array term.
///
/// Source `resolveFloatArrayCv_`, which walks the whole `MS:1000513` subtree
/// looking for a term whose *name* equals `array_name`. Here that search is
/// [`ControlledVocabulary::first_child_with_name`], the same depth-first
/// descendant walk with the vocabulary's own work and byte budget applied, so
/// one lookup cannot become unbounded.
///
/// A unit is written only when the term declares exactly one allowed unit. A
/// term with several — `MS:1003007` "raw ion mobility array" allows both
/// milliseconds and seconds — gets no unit attributes rather than an
/// arbitrarily chosen one, and the upstream suite asserts exactly that.
///
/// # Arguments
///
/// * `array_name` — the name as `DataArray::name` carries it. An empty name
///   never reaches this function from [`store`]; passed here it resolves to
///   `MS:1000786`.
/// * `cv` — the vocabulary to search, normally
///   [`ControlledVocabulary::psi_ms`].
///
/// # Errors
///
/// [`Error::InvalidValue`] when the vocabulary walk exceeds its own budget.
/// A term that declares a unit identifier the vocabulary does not define
/// yields no unit attributes rather than an error, because the source's
/// `cv.exists(unit_id)` guard does the same.
pub fn resolve_float_array_cv(
    array_name: &str,
    cv: &ControlledVocabulary,
) -> Result<ResolvedArrayCv> {
    let Some(term) = cv.first_child_with_name("MS:1000513", array_name)? else {
        return Ok(ResolvedArrayCv {
            accession: "MS:1000786".to_owned(),
            name: "non-standard data array".to_owned(),
            non_standard: true,
            ..ResolvedArrayCv::default()
        });
    };
    let mut resolved = ResolvedArrayCv {
        accession: term.id.clone(),
        name: term.name.clone(),
        ..ResolvedArrayCv::default()
    };
    // Source: write a unit only when the term allows exactly one.
    if term.units.len() == 1 {
        if let Some(unit_id) = term.units.iter().next() {
            if cv.exists(unit_id) {
                let unit = cv.get_term(unit_id)?;
                resolved.unit_cv_ref = unit_id.chars().take(2).collect();
                resolved.unit_accession = unit_id.clone();
                resolved.unit_name = unit.name.clone();
            }
        }
    }
    Ok(resolved)
}

/// Read the dataset-level imaging metadata off an experiment's meta values.
///
/// Source `extractMeta_`. Every field is optional: an absent key leaves the
/// [`ImzMLMeta`] default in place, which is how an experiment that never came
/// from an imzML file can still be stored. The keys are the ones the readers in
/// [`imzml_handler`](crate::format::imzml_handler) mirror onto an experiment,
/// so a load followed by a store round-trips them: `imzml:imaging_mode`,
/// `imzml:max_count_x`, `imzml:max_count_y`, `imzml:max_count_z`,
/// `imzml:pixel_size_x`, `imzml:pixel_size_y`, `imzml:max_dim_x`,
/// `imzml:max_dim_y`, `imzml:uuid`, `imzml:scan_pattern`,
/// `imzml:scan_direction`, `imzml:line_scan_direction` and `imzml:polarity`.
///
/// `imzml:imaging_mode` is only read here; whether it is honoured is
/// [`storage_mode`]'s decision. A value that is neither `"continuous"` nor
/// `"processed"` becomes `None`, which is how the source's later string
/// comparisons treat it: neither branch matches and the mode is auto-detected.
///
/// # Notes
///
/// The six vocabulary keys — `imzml:imaging_mode`, `imzml:uuid`,
/// `imzml:scan_pattern`, `imzml:scan_direction`, `imzml:line_scan_direction`
/// and `imzml:polarity` — are read with the source's *lenient*
/// stringification. `extractMeta_` reads each of them through
/// `DataValue::toString()`, which `StringUtils.h` documents in terms as the
/// stringification that never throws: a number becomes its decimal text, a
/// list is joined as `[a, b]`, an empty value becomes `""`, and a genuine
/// string is returned verbatim. A private `meta_text` reproduces that, so an
/// experiment whose `imzml:polarity` is an integer — which stores fine
/// upstream — stores here too rather than being refused. None of the six is a
/// number in practice; the leniency matters because none of them is *required*
/// to be a string for the store to succeed.
///
/// The seven numeric keys keep the source's strictness, because the source is
/// strict there: `imzml:max_count_x`, `imzml:max_count_y` and
/// `imzml:max_count_z` go through `static_cast<UInt>`, which throws for
/// anything but an integer and for a negative one, and `imzml:pixel_size_x`,
/// `imzml:pixel_size_y`, `imzml:max_dim_x` and `imzml:max_dim_y` go through
/// `static_cast<double>`, which accepts an integer and throws for an empty
/// value. That `static_cast<double>` reads the `DataValue` union's `double`
/// member unconditionally for every other type, so a *string* pixel size is
/// undefined behaviour upstream rather than an error; refusing it here is the
/// one place this function is deliberately stricter than the source, alongside
/// the out-of-range count below.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a count is not an integer, is negative or is
/// above `2^32`, or when a size is neither an integer nor a float. The source
/// reaches the same outcome for the wrong type through `DataValue`'s
/// `static_cast`, which throws `Exception::ConversionError`; for an
/// out-of-range count its `static_cast<UInt>` wraps silently instead. No
/// vocabulary key can fail.
pub fn dataset_meta(exp: &MSExperiment) -> Result<ImzMLMeta> {
    let metadata = &exp.settings.metadata;
    let text = |key: &str| -> String {
        match metadata.get(key) {
            Some(value) => meta_text(value),
            None => String::new(),
        }
    };
    let count = |key: &str| -> Result<u32> {
        match metadata.get(key) {
            Some(value) => {
                let raw = value.as_i64().map_err(|_| meta_type(key, "an integer"))?;
                u32::try_from(raw).map_err(|_| {
                    Error::InvalidValue(format!(
                        "imzML meta value '{key}' ({raw}) is outside the uint32 range an imzML \
                         pixel count occupies"
                    ))
                })
            }
            None => Ok(0),
        }
    };
    let size = |key: &str| -> Result<f64> {
        match metadata.get(key) {
            Some(value) => value.as_f64().map_err(|_| meta_type(key, "numeric")),
            None => Ok(0.0),
        }
    };

    Ok(ImzMLMeta {
        imaging_mode: match text("imzml:imaging_mode").as_str() {
            "continuous" => Some(ImagingMode::Continuous),
            "processed" => Some(ImagingMode::Processed),
            _ => None,
        },
        max_count_x: count("imzml:max_count_x")?,
        max_count_y: count("imzml:max_count_y")?,
        max_count_z: count("imzml:max_count_z")?,
        pixel_size_x: size("imzml:pixel_size_x")?,
        pixel_size_y: size("imzml:pixel_size_y")?,
        max_dim_x: size("imzml:max_dim_x")?,
        max_dim_y: size("imzml:max_dim_y")?,
        uuid: text("imzml:uuid"),
        scan_pattern: text("imzml:scan_pattern"),
        scan_direction: text("imzml:scan_direction"),
        line_scan_direction: text("imzml:line_scan_direction"),
        polarity: text("imzml:polarity"),
        ..ImzMLMeta::default()
    })
}

/// The source's lenient `DataValue::toString()` over a [`MetaValue`].
///
/// `DataValue.cpp`'s `toString(full_precision = true)`, which is what
/// `extractMeta_` reaches for and what `StringUtils::toStr(const DataValue&)`
/// documents as the stringification that never throws:
///
/// * an empty value is the empty string,
/// * a string is returned verbatim, with no quoting,
/// * an integer is its decimal text,
/// * a float is [`float_text`] — the same 15-digit rule every float this
///   writer emits goes through, so `5.0` is `"5.0"` and not `"5"`,
/// * a list is `[`, its elements joined with `", "`, then `]`, each element
///   stringified as the corresponding scalar. That is the source's
///   `operator<<(std::ostream&, const std::vector<T>&)`, which also writes the
///   brackets and the `", "` separator, and stringifies each element with
///   `StringUtils::toStr`.
///
/// A unit is not appended, matching `DataValue::toString`, which has no access
/// to one.
///
/// This is deliberately not [`MetaValue`]'s own `Display`: that renders a float
/// with Rust's shortest round-trip formatting, and the source renders it with
/// `Internal::NumericFormatting::appendNumeric`.
fn meta_text(value: &MetaValue) -> String {
    fn joined<T>(values: &[T], element: impl Fn(&T) -> String) -> String {
        let mut out = String::from("[");
        for (index, value) in values.iter().enumerate() {
            if index != 0 {
                out.push_str(", ");
            }
            out.push_str(&element(value));
        }
        out.push(']');
        out
    }
    match value.data() {
        MetaValueData::Empty => String::new(),
        MetaValueData::String(text) => text.clone(),
        MetaValueData::Integer(number) => number.to_string(),
        MetaValueData::Float(number) => float_text(*number),
        MetaValueData::StringList(values) => joined(values, Clone::clone),
        MetaValueData::IntegerList(values) => joined(values, i64::to_string),
        MetaValueData::FloatList(values) => joined(values, |&number| float_text(number)),
    }
}

/// Whether every spectrum carries the same m/z axis, within `tolerance`.
///
/// Source `spectraShareMz_`. The first non-empty spectrum in document order is
/// the reference; every other spectrum must have the same peak count and every
/// m/z within `tolerance` of the reference's at the same position. An
/// experiment with no spectra, or one where a spectrum is empty while another
/// is not, does not share an axis. A single non-empty spectrum trivially does,
/// which is why a one-pixel dataset is written continuous unless
/// `imzml:imaging_mode` says otherwise.
///
/// # Arguments
///
/// * `tolerance` — absolute m/z agreement required, in the m/z unit. `0.0`
///   demands bit-identical axes, up to `-0.0 == 0.0`. A negative or NaN
///   tolerance makes every comparison fail, so the answer is `false` for any
///   experiment with more than one non-empty spectrum; [`store`] rejects such a
///   tolerance before it gets here.
///
/// # Notes
///
/// This is the one comparison in the writer whose cost is the product of the
/// spectrum count and the peak count. [`store`] checks
/// [`ImzMLWriteLimits::max_spectra`] and [`ImzMLWriteLimits::max_total_peaks`]
/// before calling it, which bounds it; called directly it costs what the
/// experiment in hand costs. A non-finite m/z compares unequal to itself, so a
/// spectrum carrying one never shares an axis.
pub fn spectra_share_mz(exp: &MSExperiment, tolerance: f64) -> bool {
    let Some(reference) = exp.spectra.iter().find(|s| !s.peaks.is_empty()) else {
        return false;
    };
    exp.spectra.iter().all(|spectrum| {
        std::ptr::eq(spectrum, reference)
            || (spectrum.peaks.len() == reference.peaks.len()
                && spectrum
                    .peaks
                    .iter()
                    .zip(&reference.peaks)
                    .all(|(a, b)| (a.mz - b.mz).abs() <= tolerance))
    })
}

/// Decide which of the two imzML binary layouts to write.
///
/// Source `isContinuousMode_`. An explicit `"continuous"` in `meta` is honoured
/// only if the spectra really do share an m/z axis, an explicit `"processed"`
/// is honoured unconditionally, and no declaration at all leaves
/// [`spectra_share_mz`] to decide.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `meta.imaging_mode` is
/// [`ImagingMode::Continuous`] but the spectra do not share an m/z axis within
/// `tolerance`. The source raises `Exception::InvalidParameter`, for which this
/// crate has no variant. With the default `tolerance` of `0.0` this is also the
/// outcome for axes that agree only approximately; use
/// [`ImzMLWriteOptions::source`] to accept and discard that disagreement.
pub fn storage_mode(exp: &MSExperiment, meta: &ImzMLMeta, tolerance: f64) -> Result<ImagingMode> {
    match meta.imaging_mode {
        Some(ImagingMode::Continuous) if spectra_share_mz(exp, tolerance) => {
            Ok(ImagingMode::Continuous)
        }
        Some(ImagingMode::Continuous) => Err(Error::InvalidValue(format!(
            "continuous imzML requires all spectra to share the same m/z axis and array length \
             (compared with an absolute tolerance of {tolerance})"
        ))),
        Some(ImagingMode::Processed) => Ok(ImagingMode::Processed),
        None if spectra_share_mz(exp, tolerance) => Ok(ImagingMode::Continuous),
        None => Ok(ImagingMode::Processed),
    }
}

/// Apply the spectrum and peak filters of `options` to `exp` in place, and
/// return how many spectra were removed.
///
/// Source `applyStoreOptions_`, in its order: MS level, retention time,
/// precursor m/z, then the per-peak m/z and intensity filters, then
/// `metadata_only`. Every range comparison reproduces `DRange<1>::encloses`,
/// which is closed at the minimum and **open at the maximum**.
///
/// The precursor filter looks only at spectra above MS level 1 that carry a
/// precursor, and only at the first of them; a spectrum with none is kept.
///
/// `sort_spectra_by_mz` sorts each non-empty, not-already-sorted spectrum by
/// m/z. The source runs that sort inside the peak-filter branch and again in an
/// `else` branch, so it happens either way; this does it once, before the peak
/// filters, which is the same observable order.
///
/// `metadata_only` clears peaks and every data array but keeps the metadata, so
/// the pixel coordinates survive and the dataset is still storable — with
/// zero-length arrays.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a spectrum carries a non-finite m/z or
/// intensity. That is checked by [`MSSpectrum::sort_by_position`] before
/// anything is mutated, so a rejected spectrum is left untouched. The source
/// does not sort safely and would reorder peaks away from their annotations.
///
/// # Notes
///
/// A data array whose length is neither zero nor the peak count is **not** an
/// error here, because it is not one in the writer: `store` skips such an array
/// and reports it in [`StoreReport::skipped_float_arrays`], whatever `options`
/// says. [`MSSpectrum::sort_by_position`] and [`MSSpectrum::select`] do reject
/// one — they must, since reordering peaks under a mis-sized annotation array
/// is exactly the corruption they exist to prevent — so a misaligned array is
/// lifted off the spectrum for the duration of the sort and the peak filters
/// and put back, at its original position, before this function returns. Both
/// the array and its original length therefore survive into the report.
///
/// Without that, a store's treatment of one misaligned array turned on an
/// unrelated option: with default `options` neither the sort nor
/// [`MSSpectrum::select`] runs and the array is skipped, while any sort or
/// trimming filter made the same experiment fail outright. The source diverges
/// the same way — `applyStoreOptions_` reaches `MSSpectrum::sortByPosition` and
/// `MSSpectrum::select`, both of which run `checkDataArraySizes_` and throw
/// `Exception::Precondition` — but its writer's own explicit and commented
/// policy, in `appendAndWriteFloatDataArrays_`, is to skip and warn. That is
/// the behaviour made unconditional.
///
/// The source ends with `exp.updateRanges()`. This crate computes ranges on
/// demand in [`MSExperiment::ranges`], so there is no cache to refresh and the
/// call has no counterpart.
///
/// # Warning
///
/// `metadata_only` and a declared continuous mode do not combine. Clearing
/// every peak leaves [`spectra_share_mz`] with no non-empty reference, so
/// [`storage_mode`] refuses the `"continuous"` the experiment declares — which
/// is exactly what loading any continuous imzML puts there. The source behaves
/// identically and raises `Exception::InvalidParameter` from
/// `isContinuousMode_`, so a metadata-only store of a continuous dataset fails
/// there too. Remove `imzml:imaging_mode` to store the geometry alone; the
/// result is a processed dataset with zero-length arrays.
pub fn apply_store_options(exp: &mut MSExperiment, options: &PeakFileOptions) -> Result<usize> {
    let before = exp.spectra.len();
    if options.has_ms_levels() {
        exp.spectra.retain(|spectrum| {
            i32::try_from(spectrum.ms_level).is_ok_and(|level| options.contains_ms_level(level))
        });
    }
    if options.has_rt_range() {
        let range = options.rt_range();
        exp.spectra.retain(|spectrum| encloses(range, spectrum.rt));
    }
    if options.has_precursor_mz_range() {
        let range = options.precursor_mz_range();
        exp.spectra.retain(|spectrum| {
            if spectrum.ms_level <= 1 {
                return true;
            }
            match spectrum.precursors.first() {
                Some(precursor) => encloses(range, precursor.mz),
                None => true,
            }
        });
    }

    let mz_range = options.has_mz_range().then(|| options.mz_range());
    let intensity_range = options
        .has_intensity_range()
        .then(|| options.intensity_range());
    for spectrum in &mut exp.spectra {
        // Neither the sort nor the peak filter may see a misaligned annotation
        // array: both refuse one, and the writer's policy for one is to skip
        // it, not to fail. Lifted off here, restored below whether the filters
        // succeeded or not, so a rejected spectrum is left as it was found.
        let misaligned = detach_misaligned_arrays(spectrum);
        let filtered = filter_peaks(spectrum, options, mz_range, intensity_range);
        restore_misaligned_arrays(spectrum, misaligned);
        filtered?;
    }

    if options.metadata_only {
        for spectrum in &mut exp.spectra {
            spectrum.clear(false);
        }
    }
    Ok(before - exp.spectra.len())
}

/// One spectrum's m/z sort and peak filters, source `applyStoreOptions_`'s
/// per-spectrum body.
///
/// Split out of [`apply_store_options`] only so that its failure can be
/// returned past the restoration of the arrays
/// [`detach_misaligned_arrays`] lifted off.
fn filter_peaks(
    spectrum: &mut MSSpectrum,
    options: &PeakFileOptions,
    mz_range: Option<NumericRange>,
    intensity_range: Option<NumericRange>,
) -> Result<()> {
    if options.sort_spectra_by_mz && !spectrum.peaks.is_empty() && !spectrum.is_sorted() {
        spectrum.sort_by_position()?;
    }
    if mz_range.is_none() && intensity_range.is_none() {
        return Ok(());
    }
    let kept: Vec<usize> = spectrum
        .peaks
        .iter()
        .enumerate()
        .filter(|(_, peak)| {
            mz_range.is_none_or(|range| encloses(range, peak.mz))
                && intensity_range.is_none_or(|range| encloses(range, f64::from(peak.intensity)))
        })
        .map(|(index, _)| index)
        .collect();
    if kept.len() != spectrum.peaks.len() {
        spectrum.select(&kept)?;
    }
    Ok(())
}

/// The float, integer and string data arrays of one spectrum whose length is
/// neither zero nor the peak count, with the position each held.
///
/// Only [`detach_misaligned_arrays`] builds one and only
/// [`restore_misaligned_arrays`] consumes it, so it never outlives one
/// spectrum's pass through the peak filters.
struct MisalignedArrays {
    float: Vec<(usize, DataArray<f32>)>,
    integer: Vec<(usize, DataArray<i32>)>,
    string: Vec<(usize, DataArray<String>)>,
}

/// Lift every data array that does not hold one value per peak off `spectrum`.
///
/// A zero-length array is left in place: that is the length
/// [`MSSpectrum::validate_data_arrays`] accepts for an unpopulated array, and
/// the one [`Plan::aux`] reports as [`FloatArraySkipReason::Empty`].
fn detach_misaligned_arrays(spectrum: &mut MSSpectrum) -> MisalignedArrays {
    fn detach<T>(arrays: &mut Vec<DataArray<T>>, peaks: usize) -> Vec<(usize, DataArray<T>)> {
        let mut detached = Vec::new();
        let mut aligned = Vec::with_capacity(arrays.len());
        for (position, array) in std::mem::take(arrays).into_iter().enumerate() {
            if array.data.is_empty() || array.data.len() == peaks {
                aligned.push(array);
            } else {
                detached.push((position, array));
            }
        }
        *arrays = aligned;
        detached
    }
    let peaks = spectrum.peaks.len();
    MisalignedArrays {
        float: detach(&mut spectrum.float_data_arrays, peaks),
        integer: detach(&mut spectrum.integer_data_arrays, peaks),
        string: detach(&mut spectrum.string_data_arrays, peaks),
    }
}

/// Put back what [`detach_misaligned_arrays`] lifted off, each array at the
/// position it held.
///
/// `detached` is in ascending position order, so inserting in that order places
/// each array back at its original index: every array that was before it is
/// already present again. The insertion point is clamped to the current length,
/// which no caller can reach — `restore` is only ever handed this function's own
/// output — so that a future one cannot panic here either.
fn restore_misaligned_arrays(spectrum: &mut MSSpectrum, detached: MisalignedArrays) {
    fn restore<T>(arrays: &mut Vec<DataArray<T>>, detached: Vec<(usize, DataArray<T>)>) {
        for (position, array) in detached {
            let at = position.min(arrays.len());
            arrays.insert(at, array);
        }
    }
    restore(&mut spectrum.float_data_arrays, detached.float);
    restore(&mut spectrum.integer_data_arrays, detached.integer);
    restore(&mut spectrum.string_data_arrays, detached.string);
}

/// Lower-case hex RFC 1321 MD5 of `data`.
///
/// Ported from the source's own inline `md5ProcessBlock_` and `md5Hex_`, which
/// exist because imzML 1.1.0 defines `IMS:1000090` "ibd MD5" and OpenMS has no
/// MD5 elsewhere. Neither does this crate, and adding a dependency is a
/// manifest change outside this package, so the algorithm is ported rather than
/// imported.
///
/// # Warning
///
/// MD5 is not collision resistant and must not be used to authenticate
/// anything. It is here only because imzML specifies it as a file-integrity
/// fingerprint, alongside the `IMS:1000091` SHA-1 that [`store`] also writes.
///
/// # Notes
///
/// [`store`] does not call this: it streams the `.ibd` through the same block
/// function in fixed-size chunks, so the digest costs a constant amount of
/// memory whatever the file's size. This entry point takes a slice because that
/// is what a caller checking it against a published test vector has.
pub fn md5_hex(data: &[u8]) -> String {
    let mut digest = Md5::new();
    digest.update(data);
    digest.finish_hex()
}

/// The RFC 4122 version-5 identifier of `seed` under [`IBD_UUID_NAMESPACE`].
///
/// `SHA-1(namespace || seed)` truncated to 16 bytes, with the version nibble
/// set to 5 and the variant bits to `10`. [`store`] calls this with the `.ibd`
/// payload — everything after the 16-byte header — when the experiment declares
/// no usable `imzml:uuid`, so the same arrays always produce the same
/// identifier and a written dataset is reproducible byte for byte.
///
/// The source instead draws 16 bytes from `std::random_device` and stamps them
/// as version 4. That is a reasonable choice for a "universally unique"
/// identifier and an impossible one here: this crate has no random-number
/// dependency. The consequence to know is that two datasets with identical
/// binary payloads receive identical identifiers. The identifier's job in imzML
/// is to bind one `.imzML` to one `.ibd`, which it still does exactly; under
/// this derivation it is not a dataset-level provenance key.
pub fn derive_uuid(seed: &[u8]) -> [u8; IBD_UUID_BYTES] {
    let mut digest = Sha1::new();
    digest.update(IBD_UUID_NAMESPACE);
    digest.update(seed);
    stamp_uuid_v5(&digest.finalize())
}

/// Write `exp` as an imzML dataset with default ceilings and no progress
/// output.
///
/// Equivalent to [`store_with_options`] with [`ImzMLWriteOptions::default`] and
/// a fresh [`ProgressLogger`], whose default log type is
/// [`ProgressLogType::None`](crate::concept::progress_logger::ProgressLogType::None)
/// and therefore silent.
///
/// # Errors
///
/// As [`store_with_options`].
///
/// # Examples
///
/// ```
/// use openms::format::PeakFileOptions;
/// use openms::format::imzml_handler::{ImagingMode, ImzMLHandler};
/// use openms::format::imzml_writer::store;
/// use openms::kernel::{MSExperiment, MSSpectrum, Peak1D};
/// use openms::metadata::MetaValue;
///
/// let mut experiment = MSExperiment::new();
/// for (x, y) in [(1u32, 1u32), (2, 1)] {
///     let mut spectrum =
///         MSSpectrum::from(vec![Peak1D::new(100.0, 10.0), Peak1D::new(200.0, 20.0)]);
///     spectrum.metadata.insert("imzml:x".into(), MetaValue::from(x));
///     spectrum.metadata.insert("imzml:y".into(), MetaValue::from(y));
///     experiment.spectra.push(spectrum);
/// }
///
/// let directory = std::env::temp_dir().join("openms-imzml-writer-doctest");
/// std::fs::create_dir_all(&directory)?;
/// let path = directory.join("example.imzML");
/// let report = store(&path, &experiment, &PeakFileOptions::default())?;
///
/// // Both spectra share an m/z axis, so the shared array is stored once.
/// assert_eq!(report.mode, ImagingMode::Continuous);
/// assert_eq!(report.spectra_written, 2);
///
/// let mut handler = ImzMLHandler::open(&path)?;
/// assert_eq!(handler.len(), 2);
/// assert_eq!(handler.index()[0].mz_offset, handler.index()[1].mz_offset);
/// let decoded = handler.spectrum_at_coord(2, 1, 1)?;
/// assert_eq!(decoded.spectrum.peaks.len(), 2);
/// assert!((decoded.spectrum.peaks[1].mz - 200.0).abs() < 1e-9);
/// # std::fs::remove_dir_all(&directory).ok();
/// # Ok::<(), openms::Error>(())
/// ```
pub fn store(
    imzml_path: impl AsRef<Path>,
    exp: &MSExperiment,
    options: &PeakFileOptions,
) -> Result<StoreReport> {
    store_with_options(
        imzml_path,
        exp,
        options,
        &ImzMLWriteOptions::default(),
        &mut ProgressLogger::new(),
    )
}

/// Write `exp` as an imzML dataset: `imzml_path` and the `.ibd` sibling
/// [`infer_ibd_path`] names.
///
/// Source `ImzMLWriter::store`. `exp` is never modified: the filters of
/// `options` are applied to a clone, as the source applies them to its
/// `MSExperiment work = exp;`.
///
/// What is written, in order: the `.ibd` — its 16-byte UUID header, then in
/// continuous mode one shared m/z array followed by every spectrum's intensity
/// and auxiliary arrays, or in processed mode every spectrum's m/z, intensity
/// and auxiliary arrays — and then the `.imzML`, whose per-spectrum
/// `IMS:1000102`, `IMS:1000103` and `IMS:1000104` params point into it. m/z and
/// intensity precision follow [`PeakFileOptions::mz_32_bit`] and
/// [`PeakFileOptions::intensity_32_bit`]; auxiliary arrays are always 32-bit
/// float. Nothing is compressed: `MS:1000576` "no compression" is written on
/// every array, which the readers require of an external array.
///
/// Pixel coordinates come from each spectrum's `imzml:x`, `imzml:y` and
/// optional `imzml:z` meta values, and the dataset metadata from the
/// experiment's own through [`dataset_meta`]. Auxiliary arrays come from
/// `MSSpectrum::float_data_arrays` and are named through
/// [`resolve_float_array_cv`]; integer and string data arrays have no imzML
/// representation and are dropped. Everything dropped or skipped appears in the
/// returned [`StoreReport`], where the source puts it in the log.
///
/// # Arguments
///
/// * `imzml_path` — the output `.imzML`. The `.ibd` is its sibling: a
///   case-insensitive `.imzML` suffix is replaced, any other name gains `.ibd`.
/// * `exp` — the experiment to store. It must have at least one spectrum left
///   after filtering, and every spectrum must carry `imzml:x` and `imzml:y`.
/// * `options` — filtering, m/z sorting and binary precision.
/// * `write` — resource ceilings and the shared-m/z tolerance.
/// * `logger` — progress output, reported over `0..=spectra + 2`. The range is
///   always ended, including on the error paths, so the logger's nesting depth
///   stays balanced; the source leaks its own nesting counter when `store`
///   throws.
///
/// # Errors
///
/// * [`Error::MissingInformation`] when the filtered experiment has no spectrum
///   left, or when a spectrum has no `imzml:x` or no `imzml:y`.
/// * [`Error::InvalidValue`] when a pixel coordinate is below 1 or above
///   `u32`, when a coordinate or dataset meta value has the wrong type, when
///   `write.shared_mz_tolerance` is not finite and non-negative, when
///   `meta.imaging_mode` demands a continuous layout the spectra cannot
///   support, when a string the XML must carry holds characters XML 1.0 cannot
///   represent, or when the dataset exceeds one of `write.limits`. The source
///   raises `Exception::InvalidValue` for the coordinates,
///   `Exception::InvalidParameter` for the mode and nothing at all for the
///   rest.
/// * [`Error::Io`] when either file cannot be created, written, re-read for its
///   checksums or flushed, where the source raises
///   `Exception::UnableToCreateFile` and `Exception::ParseError`.
///
/// Every check that does not need the filesystem runs before either file is
/// created, so a rejected experiment leaves nothing behind. Once writing has
/// begun an I/O failure leaves the partial files in place, as it does in the
/// source; the `.imzML` is written last, so a dataset that has an `.imzML` has
/// a complete `.ibd`.
///
/// # Notes
///
/// Spectra sharing a pixel coordinate are written out as-is and counted in
/// [`StoreReport::duplicate_pixel_count`], matching what the readers accept: a
/// dataset that loads must be storable again. Readers map only the first
/// spectrum per pixel into the imaging geometry.
///
/// `imzml:max_dim_x` and `imzml:max_dim_y` are recomputed as pixel size times
/// pixel count whenever the corresponding pixel size is positive, so a value
/// carried on the experiment is overwritten rather than trusted. The grid
/// counts are only ever raised, never lowered, so a declared
/// `imzml:max_count_x` larger than any pixel survives.
///
/// The retention time is written as `MS:1000016` "scan start time" in seconds
/// for every spectrum, including the OpenMS unset default of `-1`, as the
/// source does. Every float `cvParam` this writer emits — that retention time,
/// the pixel sizes and the image extents — is rendered by a private
/// `float_text`, the source's 15-digit `NumericFormatting` rule, so the unset
/// default reads `value="-1.0"` rather than `value="-1"`.
///
/// `PeakFileOptions::metadata_only` over an experiment that declares
/// `imzml:imaging_mode = "continuous"` is refused; see the warning on
/// [`apply_store_options`].
pub fn store_with_options(
    imzml_path: impl AsRef<Path>,
    exp: &MSExperiment,
    options: &PeakFileOptions,
    write: &ImzMLWriteOptions,
    logger: &mut ProgressLogger,
) -> Result<StoreReport> {
    let imzml_path = imzml_path.as_ref();
    let ibd_path = infer_ibd_path(imzml_path);
    if !write.shared_mz_tolerance.is_finite() || write.shared_mz_tolerance < 0.0 {
        return Err(Error::InvalidValue(
            "imzML shared m/z tolerance must be finite and non-negative".into(),
        ));
    }

    let mut work = exp.clone();
    let filtered_out = apply_store_options(&mut work, options)?;
    if work.spectra.is_empty() {
        return Err(Error::MissingInformation(
            "cannot store an empty MSExperiment as imzML (all spectra removed by the peak file \
             options?)"
                .into(),
        ));
    }

    let plan = Plan::build(&work, options, write, &ibd_path)?;
    let steps = i64::try_from(plan.spectra.len().saturating_add(2))
        .map_err(|_| Error::InvalidValue("imzML progress range exceeds i64".into()))?;
    logger.start_progress(0, steps, "storing imzML file")?;
    let written = plan.write(&work, imzml_path, filtered_out, logger);
    let ended = logger.end_progress(written.as_ref().map_or(0, |report| report.ibd_bytes));
    let report = written?;
    ended?;
    Ok(report)
}

// ---------------------------------------------------------------------------
// Planning — everything that can fail without touching the filesystem.
// ---------------------------------------------------------------------------

/// One auxiliary array's place in the `.ibd` and identity in the XML.
#[derive(Clone, Debug)]
struct AuxPlan {
    /// Position in `MSSpectrum::float_data_arrays`.
    array: usize,
    cv: ResolvedArrayCv,
    /// `DataArray::name`, the `value` written for a non-standard array.
    array_name: String,
    offset: u64,
    count: u64,
    encoded: u64,
}

/// One spectrum's pixel and array placement.
#[derive(Clone, Debug)]
struct SpectrumPlan {
    x: u32,
    y: u32,
    z: u32,
    mz_offset: u64,
    mz_count: u64,
    mz_encoded: u64,
    int_offset: u64,
    int_count: u64,
    int_encoded: u64,
    aux: Vec<AuxPlan>,
}

/// The whole store, resolved and checked, with no file yet created.
struct Plan {
    limits: ImzMLWriteLimits,
    mode: ImagingMode,
    meta: ImzMLMeta,
    instrument_model: String,
    /// Stored m/z width, from `PeakFileOptions::mz_32_bit`.
    mz_32_bit: bool,
    /// Stored intensity width, from `PeakFileOptions::intensity_32_bit`.
    int_32_bit: bool,
    /// Index of the spectrum whose m/z array is the shared axis, continuous
    /// mode only.
    shared_mz: Option<usize>,
    /// Total `.ibd` size, UUID header included.
    ibd_bytes: u64,
    spectra: Vec<SpectrumPlan>,
    diagnostics: Diagnostics,
}

/// The report's variable-length halves, accumulated during planning.
#[derive(Default)]
struct Diagnostics {
    duplicate_pixel_count: usize,
    duplicate_pixels: Vec<DuplicatePixel>,
    dropped_data_array_count: usize,
    dropped_data_array_names: Vec<String>,
    skipped_float_array_count: usize,
    skipped_float_arrays: Vec<SkippedFloatArray>,
    aux_arrays_written: usize,
}

impl Diagnostics {
    fn duplicate(&mut self, pixel: DuplicatePixel, limits: &ImzMLWriteLimits) {
        self.duplicate_pixel_count += 1;
        if self.duplicate_pixels.len() < limits.max_reported_items {
            self.duplicate_pixels.push(pixel);
        }
    }

    fn dropped(&mut self, name: &str, limits: &ImzMLWriteLimits) {
        self.dropped_data_array_count += 1;
        if self
            .dropped_data_array_names
            .iter()
            .any(|seen| seen == name)
        {
            return;
        }
        if self.dropped_data_array_names.len() < limits.max_reported_items {
            self.dropped_data_array_names.push(name.to_owned());
        }
    }

    fn skipped(
        &mut self,
        spectrum: usize,
        name: &str,
        reason: FloatArraySkipReason,
        limits: &ImzMLWriteLimits,
    ) {
        self.skipped_float_array_count += 1;
        if self.skipped_float_arrays.len() < limits.max_reported_items {
            self.skipped_float_arrays.push(SkippedFloatArray {
                spectrum,
                name: name.to_owned(),
                reason,
            });
        }
    }
}

/// Cumulative byte budget and XML 1.0 check for experiment-derived strings.
struct TextBudget {
    spent: usize,
    limits: ImzMLWriteLimits,
}

impl TextBudget {
    fn new(limits: ImzMLWriteLimits) -> Self {
        Self { spent: 0, limits }
    }

    fn charge(&mut self, what: &str, value: &str) -> Result<()> {
        if value.len() > self.limits.max_name_bytes {
            return Err(Error::InvalidValue(format!(
                "imzML {what} is longer than the configured {} byte limit",
                self.limits.max_name_bytes
            )));
        }
        self.spent = self.spent.saturating_add(value.len());
        if self.spent > self.limits.max_text_bytes {
            return Err(Error::InvalidValue(
                "imzML document text exceeds the configured byte limit".into(),
            ));
        }
        // XMLHandler::writeXMLEscape escapes the five entity characters and
        // nothing else, so a value holding a character XML 1.0 cannot encode
        // would produce a document no parser can read back.
        if value.chars().any(|c| {
            !(matches!(c, '\t' | '\n' | '\r')
                || ('\u{20}'..='\u{d7ff}').contains(&c)
                || ('\u{e000}'..='\u{fffd}').contains(&c)
                || c >= '\u{10000}')
        }) {
            return Err(Error::InvalidValue(format!(
                "imzML {what} contains characters XML 1.0 cannot represent"
            )));
        }
        Ok(())
    }
}

/// Mutable state of one planning pass.
struct Planner<'a> {
    limits: ImzMLWriteLimits,
    text: TextBudget,
    diagnostics: Diagnostics,
    cache: BTreeMap<String, ResolvedArrayCv>,
    cv: &'a ControlledVocabulary,
    offset: u64,
}

impl<'a> Planner<'a> {
    fn new(limits: ImzMLWriteLimits, cv: &'a ControlledVocabulary) -> Self {
        Self {
            limits,
            text: TextBudget::new(limits),
            diagnostics: Diagnostics::default(),
            cache: BTreeMap::new(),
            cv,
            offset: IBD_UUID_BYTES as u64,
        }
    }

    /// Reserve `bytes` at the current `.ibd` position and return where they
    /// start. This is the one place an offset is produced, so the ceiling and
    /// the overflow check cannot be bypassed.
    fn reserve(&mut self, bytes: u64) -> Result<u64> {
        let start = self.offset;
        let end = start
            .checked_add(bytes)
            .ok_or_else(|| Error::InvalidValue("imzML .ibd byte offset overflows u64".into()))?;
        if end > self.limits.max_ibd_bytes {
            return Err(Error::InvalidValue(format!(
                "the imzML .ibd would reach {end} bytes, past the configured limit of {}",
                self.limits.max_ibd_bytes
            )));
        }
        self.offset = end;
        Ok(start)
    }

    /// Source `resolveFloatArrayCv_`'s `ArrayCvCache`: the same names repeat on
    /// every spectrum, while resolving one walks the whole `MS:1000513`
    /// subtree.
    fn resolve(&mut self, name: &str) -> Result<ResolvedArrayCv> {
        if let Some(resolved) = self.cache.get(name) {
            return Ok(resolved.clone());
        }
        if self.cache.len() >= self.limits.max_cv_lookups {
            return Err(Error::InvalidValue(
                "imzML array-name vocabulary lookups exceed the configured limit".into(),
            ));
        }
        let resolved = resolve_float_array_cv(name, self.cv)?;
        self.cache.insert(name.to_owned(), resolved.clone());
        Ok(resolved)
    }

    /// Plan the exportable auxiliary arrays of one spectrum, in document order.
    fn aux(&mut self, index: usize, spectrum: &MSSpectrum) -> Result<Vec<AuxPlan>> {
        let peaks = spectrum.peaks.len();
        let mut planned = Vec::new();
        for (position, array) in spectrum.float_data_arrays.iter().enumerate() {
            // Source order: empty first (silently), then unnamed, then the
            // per-peak length contract, then the reserved peak-array names.
            if array.data.is_empty() {
                self.diagnostics.skipped(
                    index,
                    &array.name,
                    FloatArraySkipReason::Empty,
                    &self.limits,
                );
                continue;
            }
            if array.name.is_empty() {
                self.diagnostics
                    .skipped(index, "", FloatArraySkipReason::Unnamed, &self.limits);
                continue;
            }
            if array.data.len() != peaks {
                self.diagnostics.skipped(
                    index,
                    &array.name,
                    FloatArraySkipReason::LengthMismatch {
                        length: array.data.len(),
                        peaks,
                    },
                    &self.limits,
                );
                continue;
            }
            let cv = self.resolve(&array.name)?;
            if matches!(cv.accession.as_str(), "MS:1000514" | "MS:1000515") {
                self.diagnostics.skipped(
                    index,
                    &array.name,
                    FloatArraySkipReason::ReservedPeakArrayName {
                        accession: cv.accession,
                    },
                    &self.limits,
                );
                continue;
            }
            if planned.len() >= self.limits.max_aux_arrays_per_spectrum {
                return Err(Error::InvalidValue(format!(
                    "imzML auxiliary arrays on spectrum {index} exceed the configured limit of {}",
                    self.limits.max_aux_arrays_per_spectrum
                )));
            }
            self.diagnostics.aux_arrays_written += 1;
            if self.diagnostics.aux_arrays_written > self.limits.max_total_aux_arrays {
                return Err(Error::InvalidValue(format!(
                    "imzML auxiliary arrays exceed the configured limit of {}",
                    self.limits.max_total_aux_arrays
                )));
            }
            self.text.charge("float data array name", &array.name)?;
            // FloatDataArray is always written as float32.
            let encoded = array_bytes(array.data.len(), 4)?;
            planned.push(AuxPlan {
                array: position,
                cv,
                array_name: array.name.clone(),
                offset: self.reserve(encoded)?,
                count: array.data.len() as u64,
                encoded,
            });
        }
        Ok(planned)
    }
}

impl Plan {
    fn build(
        work: &MSExperiment,
        options: &PeakFileOptions,
        write: &ImzMLWriteOptions,
        ibd_path: &Path,
    ) -> Result<Self> {
        let limits = write.limits;
        if work.spectra.len() > limits.max_spectra {
            return Err(Error::InvalidValue(format!(
                "imzML spectra ({}) exceed the configured limit of {}",
                work.spectra.len(),
                limits.max_spectra
            )));
        }
        let mut total_peaks = 0u64;
        for spectrum in &work.spectra {
            if spectrum.peaks.len() > limits.max_peaks_per_spectrum {
                return Err(Error::InvalidValue(format!(
                    "an imzML spectrum's peak count ({}) exceeds the configured limit of {}",
                    spectrum.peaks.len(),
                    limits.max_peaks_per_spectrum
                )));
            }
            total_peaks = total_peaks.saturating_add(spectrum.peaks.len() as u64);
        }
        if total_peaks > limits.max_total_peaks {
            return Err(Error::InvalidValue(format!(
                "the imzML dataset's peak count ({total_peaks}) exceeds the configured limit of {}",
                limits.max_total_peaks
            )));
        }

        let mut meta = dataset_meta(work)?;
        let mode = storage_mode(work, &meta, write.shared_mz_tolerance)?;
        meta.imaging_mode = Some(mode);
        meta.ibd_file_path = ibd_path.to_path_buf();
        meta.mz_data_type = precision(options.mz_32_bit);
        meta.int_data_type = precision(options.intensity_32_bit);
        let mz_width = element_bytes(options.mz_32_bit);
        let int_width = element_bytes(options.intensity_32_bit);

        let mut planner = Planner::new(limits, ControlledVocabulary::psi_ms()?);
        for value in [
            &meta.uuid,
            &meta.scan_pattern,
            &meta.scan_direction,
            &meta.line_scan_direction,
            &meta.polarity,
        ] {
            planner.text.charge("dataset meta value", value)?;
        }
        let instrument_model = instrument_model(work);
        planner.text.charge("instrument model", &instrument_model)?;

        // Source validatePixelMetadataForStore_ and warnOnDroppedDataArrays_,
        // merged into the pass that reads the coordinates: all three walk the
        // same spectra.
        let mut seen: BTreeSet<(u32, u32, u32)> = BTreeSet::new();
        let mut pixels = Vec::new();
        pixels
            .try_reserve_exact(work.spectra.len())
            .map_err(|_| Error::InvalidValue("cannot allocate the imzML pixel plan".into()))?;
        for (index, spectrum) in work.spectra.iter().enumerate() {
            let pixel = pixel_coord(spectrum, index)?;
            if !seen.insert(pixel) {
                planner.diagnostics.duplicate(
                    DuplicatePixel {
                        spectrum: index,
                        x: pixel.0,
                        y: pixel.1,
                        z: pixel.2,
                    },
                    &limits,
                );
            }
            planner
                .text
                .charge("spectrum native identifier", &spectrum.native_id)?;
            for name in spectrum
                .integer_data_arrays
                .iter()
                .filter(|array| !array.data.is_empty())
                .map(|array| array.name.as_str())
                .chain(
                    spectrum
                        .string_data_arrays
                        .iter()
                        .filter(|array| !array.data.is_empty())
                        .map(|array| array.name.as_str()),
                )
            {
                planner.diagnostics.dropped(name, &limits);
            }
            pixels.push(pixel);
        }

        // Lay the .ibd out. Continuous mode puts the one shared m/z array
        // immediately after the UUID header and points every spectrum at it.
        let shared_mz = if mode == ImagingMode::Continuous {
            Some(
                work.spectra
                    .iter()
                    .position(|spectrum| !spectrum.peaks.is_empty())
                    .ok_or_else(|| {
                        Error::InvalidValue(
                            "continuous imzML needs at least one non-empty spectrum to hold the \
                             shared m/z array"
                                .into(),
                        )
                    })?,
            )
        } else {
            None
        };
        let shared = match shared_mz {
            Some(index) => {
                let peaks = work.spectra[index].peaks.len();
                let encoded = array_bytes(peaks, mz_width)?;
                Some((planner.reserve(encoded)?, peaks as u64, encoded))
            }
            None => None,
        };

        let mut spectra = Vec::new();
        spectra
            .try_reserve_exact(work.spectra.len())
            .map_err(|_| Error::InvalidValue("cannot allocate the imzML spectrum plan".into()))?;
        for (index, spectrum) in work.spectra.iter().enumerate() {
            let peaks = spectrum.peaks.len();
            let (mz_offset, mz_count, mz_encoded) = match shared {
                Some(shared) => shared,
                None => {
                    let encoded = array_bytes(peaks, mz_width)?;
                    (planner.reserve(encoded)?, peaks as u64, encoded)
                }
            };
            let int_encoded = array_bytes(peaks, int_width)?;
            let int_offset = planner.reserve(int_encoded)?;
            let (x, y, z) = pixels[index];
            spectra.push(SpectrumPlan {
                x,
                y,
                z,
                mz_offset,
                mz_count,
                mz_encoded,
                int_offset,
                int_count: peaks as u64,
                int_encoded,
                aux: planner.aux(index, spectrum)?,
            });
        }

        // Source updateGridFromPixels_: the declared grid is only ever raised,
        // and the physical extents are recomputed from the pixel size.
        for plan in &spectra {
            meta.max_count_x = meta.max_count_x.max(plan.x);
            meta.max_count_y = meta.max_count_y.max(plan.y);
            meta.max_count_z = meta.max_count_z.max(plan.z);
        }
        if meta.max_count_z == 0 {
            meta.max_count_z = 1;
        }
        if meta.pixel_size_x > 0.0 && meta.max_count_x > 0 {
            meta.max_dim_x = meta.pixel_size_x * f64::from(meta.max_count_x);
        }
        if meta.pixel_size_y > 0.0 && meta.max_count_y > 0 {
            meta.max_dim_y = meta.pixel_size_y * f64::from(meta.max_count_y);
        }

        if planner.offset > limits.max_checksum_bytes {
            return Err(Error::InvalidValue(format!(
                "the imzML .ibd would reach {} bytes, past the configured checksum limit of {}",
                planner.offset, limits.max_checksum_bytes
            )));
        }

        Ok(Self {
            limits,
            mode,
            meta,
            instrument_model,
            mz_32_bit: options.mz_32_bit,
            int_32_bit: options.intensity_32_bit,
            shared_mz,
            ibd_bytes: planner.offset,
            spectra,
            diagnostics: planner.diagnostics,
        })
    }

    // ---------------------------------------------------------------------
    // Writing.
    // ---------------------------------------------------------------------

    fn write(
        mut self,
        work: &MSExperiment,
        imzml_path: &Path,
        spectra_filtered_out: usize,
        logger: &mut ProgressLogger,
    ) -> Result<StoreReport> {
        let ibd_path = self.meta.ibd_file_path.clone();
        let payload_digest = self.write_ibd(work, &ibd_path, logger)?;
        let uuid = match uuid_bytes(&self.meta.uuid) {
            // Source keeps the declared string verbatim, dashes and braces and
            // all, and writes only its 16 decoded bytes to the .ibd.
            Some(bytes) => bytes,
            None => {
                let bytes = stamp_uuid_v5(&payload_digest);
                self.meta.uuid = uuid_to_string(&bytes);
                bytes
            }
        };
        write_ibd_uuid_header(&ibd_path, &uuid)?;

        let (sha1, md5) = ibd_checksums(&ibd_path, self.ibd_bytes, self.limits)?;
        self.meta.ibd_sha1 = sha1;
        self.meta.ibd_md5 = md5;
        logger.set_progress(
            i64::try_from(self.spectra.len().saturating_add(1)).unwrap_or(i64::MAX),
        )?;

        self.write_xml(work, imzml_path)?;
        Ok(StoreReport {
            mode: self.mode,
            meta: self.meta,
            spectra_written: self.spectra.len(),
            spectra_filtered_out,
            ibd_bytes: self.ibd_bytes,
            aux_arrays_written: self.diagnostics.aux_arrays_written,
            duplicate_pixel_count: self.diagnostics.duplicate_pixel_count,
            duplicate_pixels: self.diagnostics.duplicate_pixels,
            dropped_data_array_count: self.diagnostics.dropped_data_array_count,
            dropped_data_array_names: self.diagnostics.dropped_data_array_names,
            skipped_float_array_count: self.diagnostics.skipped_float_array_count,
            skipped_float_arrays: self.diagnostics.skipped_float_arrays,
        })
    }

    /// Stream the payload, returning the SHA-1 of the namespace followed by it.
    ///
    /// The 16 header bytes are written as zeros and filled in afterwards,
    /// because the identifier may be derived from the payload the header
    /// precedes.
    fn write_ibd(
        &self,
        work: &MSExperiment,
        ibd_path: &Path,
        logger: &mut ProgressLogger,
    ) -> Result<[u8; 20]> {
        let array_limits = self.limits.array_limits();
        let mut buffered = BufWriter::new(File::create(ibd_path)?);
        buffered.write_all(&[0u8; IBD_UUID_BYTES])?;
        // Only the payload seeds the identifier, so the digest starts after the
        // placeholder header and with the namespace [`derive_uuid`] prepends.
        let mut out = Hashed::new(buffered);
        out.digest.update(IBD_UUID_NAMESPACE);

        if let Some(index) = self.shared_mz {
            write_mz(
                &mut out,
                &work.spectra[index],
                self.mz_32_bit,
                &array_limits,
            )?;
        }
        for (index, (plan, spectrum)) in self.spectra.iter().zip(&work.spectra).enumerate() {
            logger.set_progress(i64::try_from(index + 1).unwrap_or(i64::MAX))?;
            if self.shared_mz.is_none() {
                write_mz(&mut out, spectrum, self.mz_32_bit, &array_limits)?;
            }
            let intensities: Vec<f32> = spectrum.peaks.iter().map(|peak| peak.intensity).collect();
            if self.int_32_bit {
                write_float32_array(&mut out, &intensities, &array_limits)?;
            } else {
                let widened: Vec<f64> = intensities.iter().map(|&v| f64::from(v)).collect();
                write_float64_array(&mut out, &widened, &array_limits)?;
            }
            for aux in &plan.aux {
                write_float32_array(
                    &mut out,
                    &spectrum.float_data_arrays[aux.array].data,
                    &array_limits,
                )?;
            }
        }

        let (mut buffered, digest) = out.into_parts();
        buffered.flush()?;
        let mut file = buffered.into_inner().map_err(|e| Error::Io(e.into()))?;
        file.flush()?;
        let mut seed = [0u8; 20];
        seed.copy_from_slice(&digest.finalize());
        Ok(seed)
    }

    fn write_xml(&self, work: &MSExperiment, imzml_path: &Path) -> Result<()> {
        let mut w = BufWriter::new(File::create(imzml_path)?);
        let mode = match self.mode {
            ImagingMode::Continuous => ("IMS:1000030", "continuous"),
            ImagingMode::Processed => ("IMS:1000031", "processed"),
        };
        let mz_precision = precision_term(self.meta.mz_data_type);
        let int_precision = precision_term(self.meta.int_data_type);

        writeln!(w, "<?xml version=\"1.0\" encoding=\"utf-8\"?>")?;
        writeln!(
            w,
            "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">"
        )?;
        writeln!(w, "\t<cvList count=\"3\">")?;
        writeln!(
            w,
            "\t\t<cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry \
             Ontology\" version=\"4.1.30\" \
             URI=\"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo\"/>"
        )?;
        writeln!(
            w,
            "\t\t<cv id=\"IMS\" fullName=\"Imaging MS Ontology\" version=\"1.1.0\" \
             URI=\"https://raw.githubusercontent.com/imzML/imzML/master/imagingMS.obo\"/>"
        )?;
        writeln!(
            w,
            "\t\t<cv id=\"UO\" fullName=\"Unit Ontology\" version=\"09:04:2014\" \
             URI=\"https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo\"/>"
        )?;
        writeln!(w, "\t</cvList>")?;

        writeln!(w, "\t<fileDescription>")?;
        writeln!(w, "\t\t<fileContent>")?;
        cv_param(
            &mut w,
            3,
            "IMS",
            "IMS:1000080",
            "universally unique identifier",
            &self.meta.uuid,
            "",
        )?;
        if !self.meta.ibd_sha1.is_empty() {
            cv_param(
                &mut w,
                3,
                "IMS",
                "IMS:1000091",
                "ibd SHA-1",
                &self.meta.ibd_sha1,
                "",
            )?;
        }
        if !self.meta.ibd_md5.is_empty() {
            cv_param(
                &mut w,
                3,
                "IMS",
                "IMS:1000090",
                "ibd MD5",
                &self.meta.ibd_md5,
                "",
            )?;
        }
        cv_param(&mut w, 3, "IMS", mode.0, mode.1, "", "")?;
        match self.meta.polarity.as_str() {
            "positive" => cv_param(&mut w, 3, "MS", "MS:1000130", "positive scan", "", "")?,
            "negative" => cv_param(&mut w, 3, "MS", "MS:1000129", "negative scan", "", "")?,
            _ => {}
        }
        cv_param(&mut w, 3, "MS", "MS:1000294", "mass spectrum", "", "")?;
        writeln!(w, "\t\t</fileContent>")?;
        writeln!(w, "\t</fileDescription>")?;

        writeln!(w, "\t<referenceableParamGroupList count=\"2\">")?;
        for (id, accession, name, term) in [
            ("mzArray", "MS:1000514", "m/z array", mz_precision),
            (
                "intensityArray",
                "MS:1000515",
                "intensity array",
                int_precision,
            ),
        ] {
            writeln!(w, "\t\t<referenceableParamGroup id=\"{id}\">")?;
            cv_param(&mut w, 3, "MS", accession, name, "", "")?;
            cv_param(&mut w, 3, "MS", term.0, term.1, "", "")?;
            cv_param(&mut w, 3, "MS", "MS:1000576", "no compression", "", "")?;
            cv_param(&mut w, 3, "IMS", "IMS:1000101", "external data", "true", "")?;
            writeln!(w, "\t\t</referenceableParamGroup>")?;
        }
        writeln!(w, "\t</referenceableParamGroupList>")?;

        writeln!(w, "\t<softwareList count=\"1\">")?;
        writeln!(
            w,
            "\t\t<software id=\"sw1\" version=\"{}\">",
            crate::CORE_SDK_VERSION
        )?;
        cv_param(
            &mut w,
            3,
            "MS",
            "MS:1000799",
            "custom unreleased software tool",
            "OpenMS",
            "",
        )?;
        writeln!(w, "\t\t</software>")?;
        writeln!(w, "\t</softwareList>")?;

        writeln!(w, "\t<scanSettingsList count=\"1\">")?;
        writeln!(w, "\t\t<scanSettings id=\"scanSettings1\">")?;
        if self.meta.max_count_x > 0 {
            let value = self.meta.max_count_x.to_string();
            cv_param(
                &mut w,
                3,
                "IMS",
                "IMS:1000042",
                "max count of pixels x",
                &value,
                "",
            )?;
        }
        if self.meta.max_count_y > 0 {
            let value = self.meta.max_count_y.to_string();
            cv_param(
                &mut w,
                3,
                "IMS",
                "IMS:1000043",
                "max count of pixels y",
                &value,
                "",
            )?;
        }
        for (accession, name, value) in [
            ("IMS:1000044", "max dimension x", self.meta.max_dim_x),
            ("IMS:1000045", "max dimension y", self.meta.max_dim_y),
            ("IMS:1000046", "pixel size x", self.meta.pixel_size_x),
            ("IMS:1000047", "pixel size y", self.meta.pixel_size_y),
        ] {
            if value > 0.0 {
                cv_param(
                    &mut w,
                    3,
                    "IMS",
                    accession,
                    name,
                    &float_text(value),
                    MICROMETER,
                )?;
            }
        }
        for (accession, name) in geometry_terms(&self.meta) {
            cv_param(&mut w, 3, "IMS", accession, name, "", "")?;
        }
        writeln!(w, "\t\t</scanSettings>")?;
        writeln!(w, "\t</scanSettingsList>")?;

        writeln!(w, "\t<instrumentConfigurationList count=\"1\">")?;
        writeln!(w, "\t\t<instrumentConfiguration id=\"IC1\">")?;
        cv_param(
            &mut w,
            3,
            "MS",
            "MS:1000031",
            "instrument model",
            &self.instrument_model,
            "",
        )?;
        writeln!(w, "\t\t</instrumentConfiguration>")?;
        writeln!(w, "\t</instrumentConfigurationList>")?;

        writeln!(w, "\t<dataProcessingList count=\"1\">")?;
        writeln!(w, "\t\t<dataProcessing id=\"dp1\">")?;
        writeln!(
            w,
            "\t\t\t<processingMethod order=\"1\" softwareRef=\"sw1\">"
        )?;
        cv_param(&mut w, 4, "MS", "MS:1000544", "Conversion to imzML", "", "")?;
        writeln!(w, "\t\t\t</processingMethod>")?;
        writeln!(w, "\t\t</dataProcessing>")?;
        writeln!(w, "\t</dataProcessingList>")?;

        writeln!(
            w,
            "\t<run id=\"ru_0\" defaultInstrumentConfigurationRef=\"IC1\">"
        )?;
        writeln!(
            w,
            "\t\t<spectrumList count=\"{}\" defaultDataProcessingRef=\"dp1\">",
            self.spectra.len()
        )?;
        for (index, (plan, spectrum)) in self.spectra.iter().zip(&work.spectra).enumerate() {
            self.write_spectrum(&mut w, index, plan, spectrum)?;
        }
        writeln!(w, "\t\t</spectrumList>")?;
        writeln!(w, "\t</run>")?;
        writeln!(w, "</mzML>")?;
        w.flush()?;
        Ok(())
    }

    fn write_spectrum(
        &self,
        w: &mut impl Write,
        index: usize,
        plan: &SpectrumPlan,
        spectrum: &MSSpectrum,
    ) -> Result<()> {
        let native_id = if spectrum.native_id.is_empty() {
            format!("spectrum={}", index + 1)
        } else {
            spectrum.native_id.clone()
        };
        writeln!(
            w,
            "\t\t\t<spectrum index=\"{index}\" id=\"{}\" defaultArrayLength=\"{}\">",
            escape(&native_id),
            spectrum.peaks.len()
        )?;
        if self.mode == ImagingMode::Continuous {
            writeln!(w, "\t\t\t\t<referenceableParamGroupRef ref=\"mzArray\"/>")?;
        }
        if spectrum.ms_level != 0 {
            let value = spectrum.ms_level.to_string();
            cv_param(w, 4, "MS", "MS:1000511", "ms level", &value, "")?;
        }
        writeln!(w, "\t\t\t\t<scanList count=\"1\">")?;
        cv_param(w, 5, "MS", "MS:1000795", "no combination", "", "")?;
        writeln!(w, "\t\t\t\t\t<scan instrumentConfigurationRef=\"IC1\">")?;
        cv_param(
            w,
            6,
            "MS",
            "MS:1000016",
            "scan start time",
            &float_text(spectrum.rt),
            SECOND,
        )?;
        for (accession, name, value) in [
            ("IMS:1000050", "position x", plan.x),
            ("IMS:1000051", "position y", plan.y),
            ("IMS:1000052", "position z", plan.z),
        ] {
            cv_param(w, 6, "IMS", accession, name, &value.to_string(), "")?;
        }
        writeln!(w, "\t\t\t\t\t</scan>")?;
        writeln!(w, "\t\t\t\t</scanList>")?;
        writeln!(
            w,
            "\t\t\t\t<binaryDataArrayList count=\"{}\">",
            2 + plan.aux.len()
        )?;
        external_array(
            w,
            4,
            "mzArray",
            plan.mz_offset,
            plan.mz_count,
            plan.mz_encoded,
        )?;
        external_array(
            w,
            4,
            "intensityArray",
            plan.int_offset,
            plan.int_count,
            plan.int_encoded,
        )?;
        for aux in &plan.aux {
            external_aux_array(w, 4, aux)?;
        }
        writeln!(w, "\t\t\t\t</binaryDataArrayList>")?;
        writeln!(w, "\t\t\t</spectrum>")?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Output helpers.
// ---------------------------------------------------------------------------

const MICROMETER: &str = "unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"";
const SECOND: &str = "unitAccession=\"UO:0000010\" unitName=\"second\" unitCvRef=\"UO\"";

/// Source `writeCvParam_`: the `value` attribute is omitted when the value is
/// empty, so a term whose meaning is its presence carries no `value=""`.
///
/// The `name` attribute is escaped here and is not in the source, which writes
/// it raw. The names this writer emits come from the PSI-MS ontology, and that
/// ontology does contain term names holding `<` and `!` — the cleavage-rule
/// regexes, none of which is a binary-array term. Escaping a name that needs no
/// escape is the identity, so this cannot diverge on any array term while it
/// removes the chance of emitting a document no parser can read.
fn cv_param(
    w: &mut impl Write,
    indent: usize,
    cv_ref: &str,
    accession: &str,
    name: &str,
    value: &str,
    unit: &str,
) -> Result<()> {
    for _ in 0..indent {
        w.write_all(b"\t")?;
    }
    write!(
        w,
        "<cvParam cvRef=\"{cv_ref}\" accession=\"{accession}\" name=\"{}\"",
        escape(name)
    )?;
    if !value.is_empty() {
        write!(w, " value=\"{}\"", escape(value))?;
    }
    if !unit.is_empty() {
        write!(w, " {unit}")?;
    }
    writeln!(w, "/>")?;
    Ok(())
}

/// Source `writeExternalBinaryArray_`: an m/z or intensity array, identified
/// through its `referenceableParamGroup`.
fn external_array(
    w: &mut impl Write,
    indent: usize,
    group: &str,
    offset: u64,
    count: u64,
    encoded: u64,
) -> Result<()> {
    let pad = "\t".repeat(indent);
    writeln!(w, "{pad}<binaryDataArray encodedLength=\"0\">")?;
    writeln!(w, "{pad}\t<referenceableParamGroupRef ref=\"{group}\"/>")?;
    external_extents(w, indent + 1, offset, count, encoded)?;
    writeln!(w, "{pad}\t<binary/>")?;
    writeln!(w, "{pad}</binaryDataArray>")?;
    Ok(())
}

/// Source `writeExternalAuxBinaryArray_`: an auxiliary array, identified
/// through inline `cvParam`s so that an arbitrary array name round-trips.
fn external_aux_array(w: &mut impl Write, indent: usize, aux: &AuxPlan) -> Result<()> {
    let pad = "\t".repeat(indent);
    writeln!(w, "{pad}<binaryDataArray encodedLength=\"0\">")?;
    cv_param(w, indent + 1, "MS", "MS:1000576", "no compression", "", "")?;
    cv_param(w, indent + 1, "MS", "MS:1000521", "32-bit float", "", "")?;
    cv_param(
        w,
        indent + 1,
        "IMS",
        "IMS:1000101",
        "external data",
        "true",
        "",
    )?;
    let unit = if aux.cv.unit_accession.is_empty() {
        String::new()
    } else {
        format!(
            "unitCvRef=\"{}\" unitAccession=\"{}\" unitName=\"{}\"",
            aux.cv.unit_cv_ref,
            aux.cv.unit_accession,
            escape(&aux.cv.unit_name)
        )
    };
    cv_param(
        w,
        indent + 1,
        "MS",
        &aux.cv.accession,
        &aux.cv.name,
        if aux.cv.non_standard {
            aux.array_name.as_str()
        } else {
            ""
        },
        &unit,
    )?;
    external_extents(w, indent + 1, aux.offset, aux.count, aux.encoded)?;
    writeln!(w, "{pad}\t<binary/>")?;
    writeln!(w, "{pad}</binaryDataArray>")?;
    Ok(())
}

/// The three IMS params that place one array inside the `.ibd`.
fn external_extents(
    w: &mut impl Write,
    indent: usize,
    offset: u64,
    count: u64,
    encoded: u64,
) -> Result<()> {
    for (accession, name, value) in [
        ("IMS:1000102", "external offset", offset),
        ("IMS:1000103", "external array length", count),
        ("IMS:1000104", "external encoded length", encoded),
    ] {
        cv_param(w, indent, "IMS", accession, name, &value.to_string(), "")?;
    }
    Ok(())
}

/// Source `writeImsGeometryCvParams_`: the recognised scan-pattern, scan-
/// direction and line-scan-direction terms, in that order. A value the source
/// does not recognise writes nothing at all.
fn geometry_terms(meta: &ImzMLMeta) -> Vec<(&'static str, &'static str)> {
    let mut terms = Vec::new();
    match meta.scan_pattern.as_str() {
        "top down" => terms.push(("IMS:1000401", "top down")),
        "bottom up" => terms.push(("IMS:1000402", "bottom up")),
        _ => {}
    }
    match meta.scan_direction.as_str() {
        "flyback" => terms.push(("IMS:1000413", "flyback")),
        "meander" => terms.push(("IMS:1000412", "meander")),
        "horizontal" => terms.push(("IMS:1000480", "horizontal")),
        "vertical" => terms.push(("IMS:1000481", "vertical")),
        _ => {}
    }
    match meta.line_scan_direction.as_str() {
        "left-right" => terms.push(("IMS:1000491", "left-right")),
        "right-left" => terms.push(("IMS:1000492", "right-left")),
        _ => {}
    }
    terms
}

/// Source `instrumentModelForExport_`: the model, else the name, else a
/// placeholder, because `MS:1000031` has to carry something.
fn instrument_model(exp: &MSExperiment) -> String {
    let instrument = &exp.settings.instrument;
    if !instrument.model.is_empty() {
        instrument.model.clone()
    } else if !instrument.name.is_empty() {
        instrument.name.clone()
    } else {
        "OpenMS export".to_owned()
    }
}

/// Source `writeMzArray_`: `PeakFileOptions::getMz32Bit` chooses the narrowing
/// writer, which loses precision on every value.
fn write_mz(
    out: impl Write,
    spectrum: &MSSpectrum,
    use_32_bit: bool,
    limits: &ImzMLReadLimits,
) -> Result<()> {
    let mz: Vec<f64> = spectrum.peaks.iter().map(|peak| peak.mz).collect();
    if use_32_bit {
        write_mz_as_float32(out, &mz, limits)
    } else {
        write_mz_as_float64(out, &mz, limits)
    }
}

/// Fill in the 16 header bytes the payload pass left as zeros.
fn write_ibd_uuid_header(ibd_path: &Path, uuid: &[u8; IBD_UUID_BYTES]) -> Result<()> {
    let mut file = std::fs::OpenOptions::new().write(true).open(ibd_path)?;
    file.seek(SeekFrom::Start(0))?;
    file.write_all(uuid)?;
    file.flush()?;
    Ok(())
}

/// Recompute both declared `.ibd` digests in one pass, as lower-case hex.
///
/// The source reads the file twice, once per digest. `expected` is the size the
/// plan reserved; a file that does not have it means an earlier write silently
/// came up short.
fn ibd_checksums(
    ibd_path: &Path,
    expected: u64,
    limits: ImzMLWriteLimits,
) -> Result<(String, String)> {
    if expected > limits.max_checksum_bytes {
        return Err(Error::InvalidValue(format!(
            "the imzML .ibd ({expected} bytes) exceeds the configured checksum limit of {}",
            limits.max_checksum_bytes
        )));
    }
    let mut file = File::open(ibd_path)?;
    let mut sha1 = Sha1::new();
    let mut md5 = Md5::new();
    let mut buffer = [0u8; 64 << 10];
    let mut total = 0u64;
    loop {
        let read = file.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        sha1.update(&buffer[..read]);
        md5.update(&buffer[..read]);
        total = total.saturating_add(read as u64);
    }
    if total != expected {
        return Err(Error::Io(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            format!("the imzML .ibd is {total} bytes where the plan reserved {expected}"),
        )));
    }
    Ok((hex(&sha1.finalize()), md5.finish_hex()))
}

/// A writer that also feeds a SHA-1, so the payload is hashed in the pass that
/// writes it rather than in a second read.
struct Hashed<W: Write> {
    inner: W,
    digest: Sha1,
}

impl<W: Write> Hashed<W> {
    fn new(inner: W) -> Self {
        Self {
            inner,
            digest: Sha1::new(),
        }
    }

    fn into_parts(self) -> (W, Sha1) {
        (self.inner, self.digest)
    }
}

impl<W: Write> Write for Hashed<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(buf)?;
        self.digest.update(&buf[..written]);
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

// ---------------------------------------------------------------------------
// Small shared helpers.
// ---------------------------------------------------------------------------

/// Read and validate one spectrum's pixel coordinate.
///
/// Source `validatePixelMetadataForStore_` and `readPixelCoord_`, merged: the
/// source validates the three meta values as `Int` and then re-reads them as
/// `uint32_t` in a second pass.
fn pixel_coord(spectrum: &MSSpectrum, index: usize) -> Result<(u32, u32, u32)> {
    let required = |key: &str| -> Result<i64> {
        let value = spectrum.metadata.get(key).ok_or_else(|| {
            Error::MissingInformation(format!(
                "spectrum {index} is missing the '{key}' meta value imzML export requires"
            ))
        })?;
        value.as_i64().map_err(|_| meta_type(key, "an integer"))
    };
    let x = required("imzml:x")?;
    let y = required("imzml:y")?;
    let z = match spectrum.metadata.get("imzml:z") {
        Some(value) => value
            .as_i64()
            .map_err(|_| meta_type("imzml:z", "an integer"))?,
        None => 1,
    };
    let coordinate = |key: &str, value: i64| -> Result<u32> {
        if value < 1 {
            return Err(Error::InvalidValue(format!(
                "imzML pixel coordinates must be >= 1: '{key}' is {value} at spectrum {index}"
            )));
        }
        u32::try_from(value).map_err(|_| {
            Error::InvalidValue(format!(
                "imzML pixel coordinate '{key}' ({value}) exceeds uint32 at spectrum {index}"
            ))
        })
    };
    Ok((
        coordinate("imzml:x", x)?,
        coordinate("imzml:y", y)?,
        coordinate("imzml:z", z)?,
    ))
}

/// Stored byte length of `count` elements of the given width.
fn array_bytes(count: usize, width: u64) -> Result<u64> {
    (count as u64)
        .checked_mul(width)
        .ok_or_else(|| Error::InvalidValue("imzML array byte length overflows u64".into()))
}

/// Source `DRange<1>::encloses`: closed at the minimum, open at the maximum.
fn encloses(range: NumericRange, value: f64) -> bool {
    !(value < range.min || value >= range.max)
}

/// Source `elementByteSize_`.
fn element_bytes(use_32_bit: bool) -> u64 {
    if use_32_bit { 4 } else { 8 }
}

fn precision(use_32_bit: bool) -> ImzMLDataType {
    if use_32_bit {
        ImzMLDataType::Float32
    } else {
        ImzMLDataType::Float64
    }
}

/// The PSI-MS binary-data-type term for a stored precision. Only the two float
/// widths can occur: [`store`] never plans an integer array.
fn precision_term(data_type: ImzMLDataType) -> (&'static str, &'static str) {
    match data_type {
        ImzMLDataType::Float32 => ("MS:1000521", "32-bit float"),
        _ => ("MS:1000523", "64-bit float"),
    }
}

fn meta_type(key: &str, expected: &str) -> Error {
    Error::InvalidValue(format!(
        "imzML meta value '{key}' is not {expected}, which imzML export requires"
    ))
}

/// The text the source writes for a `double` `cvParam` value.
///
/// `StringConversions::toString(double)` → `StringUtils::appendToStr(double,
/// std::string&)` → `Internal::NumericFormatting::appendNumeric(value, target,
/// writtenDigits<double>() == 15, fixed_format = false)`. That rule is:
///
/// * `NaN` is `"NaN"`, an infinity is `"inf"` or `"-inf"`;
/// * a non-zero `|v|` at or above `1e4`, or below `1e-2`, is written in
///   scientific notation with the shortest round-tripping mantissa, and its
///   exponent is rewritten from `std::to_chars`' `printf`-style `e+05` into
///   the historical `e05` — the `'+'` dropped, the `'-'` and the two-digit
///   zero padding kept;
/// * everything else is written fixed with 15 fractional digits;
/// * in both forms trailing zeros after the decimal point are trimmed, but at
///   least one digit is kept, so `5.0` stays `"5.0"` and never becomes `"5"`.
///
/// The implementation is the port of that header already in this crate,
/// [`crate::param::value`]'s `format_float` with `full_precision = true`, which
/// is the same `appendNumeric` call with the same 15-digit precision. Sharing
/// it is what keeps a pixel size written here byte-identical to the same value
/// written through a `Param`.
///
/// The writer's own self-audit claimed the source used `std::ostringstream`
/// here, i.e. six significant digits and no forced `.0`. It does not, and the
/// difference is visible in every float `cvParam`: a retention time of `-1`
/// is `"-1.0"`, a pixel size of `1e5` is `"1.0e05"`, and a pixel size of
/// `9999.9` carries all 15 fractional digits of the nearest double.
fn float_text(value: f64) -> String {
    crate::param::value::format_float(value, true)
}

/// Source `XMLHandler::writeXMLEscape` plus the control-character entities
/// [`mzml`](crate::format::mzml) writes, so an attribute survives a round trip
/// through a parser that normalises whitespace.
fn escape(value: &str) -> String {
    quick_xml::escape::escape(value)
        .replace('\n', "&#10;")
        .replace('\r', "&#13;")
        .replace('\t', "&#9;")
}

/// Stamp 20 SHA-1 bytes into an RFC 4122 version-5 identifier.
fn stamp_uuid_v5(digest: &[u8]) -> [u8; IBD_UUID_BYTES] {
    let mut out = [0u8; IBD_UUID_BYTES];
    out.copy_from_slice(&digest[..IBD_UUID_BYTES]);
    out[6] = (out[6] & 0x0f) | 0x50;
    out[8] = (out[8] & 0x3f) | 0x80;
    out
}

/// Source `uuidBytesToString_`: lower-case hex with dashes after bytes 4, 6, 8
/// and 10.
fn uuid_to_string(bytes: &[u8; IBD_UUID_BYTES]) -> String {
    let mut out = String::with_capacity(36);
    for (index, byte) in bytes.iter().enumerate() {
        if matches!(index, 4 | 6 | 8 | 10) {
            out.push('-');
        }
        out.push(hex_digit(byte >> 4));
        out.push(hex_digit(byte & 0x0f));
    }
    out
}

/// Source `bytesToHex_` and `sha1DigestToHex_`: lower-case, two digits per
/// byte, no separators.
fn hex(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        out.push(hex_digit(byte >> 4));
        out.push(hex_digit(byte & 0x0f));
    }
    out
}

fn hex_digit(nibble: u8) -> char {
    char::from_digit(u32::from(nibble), 16).unwrap_or('0')
}

// ---------------------------------------------------------------------------
// RFC 1321 MD5, ported from the source's inline implementation.
// ---------------------------------------------------------------------------

struct Md5 {
    state: [u32; 4],
    block: [u8; 64],
    filled: usize,
    bytes: u64,
}

impl Md5 {
    fn new() -> Self {
        Self {
            state: [0x6745_2301, 0xefcd_ab89, 0x98ba_dcfe, 0x1032_5476],
            block: [0; 64],
            filled: 0,
            bytes: 0,
        }
    }

    fn update(&mut self, data: &[u8]) {
        self.bytes = self.bytes.wrapping_add(data.len() as u64);
        let mut rest = data;
        while !rest.is_empty() {
            let take = (64 - self.filled).min(rest.len());
            self.block[self.filled..self.filled + take].copy_from_slice(&rest[..take]);
            self.filled += take;
            rest = &rest[take..];
            if self.filled == 64 {
                self.compress();
                self.filled = 0;
            }
        }
    }

    fn finish_hex(mut self) -> String {
        let bits = self.bytes.wrapping_mul(8);
        self.block[self.filled] = 0x80;
        self.filled += 1;
        if self.filled > 56 {
            self.block[self.filled..].fill(0);
            self.compress();
            self.filled = 0;
        }
        self.block[self.filled..56].fill(0);
        self.block[56..].copy_from_slice(&bits.to_le_bytes());
        self.compress();
        let mut out = String::with_capacity(32);
        for word in self.state {
            out.push_str(&hex(&word.to_le_bytes()));
        }
        out
    }

    fn compress(&mut self) {
        let mut m = [0u32; 16];
        for (word, chunk) in m.iter_mut().zip(self.block.chunks_exact(4)) {
            *word = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
        }
        let [mut a, mut b, mut c, mut d] = self.state;
        for i in 0..64 {
            let (f, g) = match i / 16 {
                0 => ((b & c) | (!b & d), i),
                1 => ((b & d) | (c & !d), (5 * i + 1) % 16),
                2 => (b ^ c ^ d, (3 * i + 5) % 16),
                _ => (c ^ (b | !d), (7 * i) % 16),
            };
            let rotated = a
                .wrapping_add(f)
                .wrapping_add(MD5_SINE[i])
                .wrapping_add(m[g])
                .rotate_left(MD5_SHIFT[i]);
            a = d;
            d = c;
            c = b;
            b = b.wrapping_add(rotated);
        }
        self.state[0] = self.state[0].wrapping_add(a);
        self.state[1] = self.state[1].wrapping_add(b);
        self.state[2] = self.state[2].wrapping_add(c);
        self.state[3] = self.state[3].wrapping_add(d);
    }
}

/// `floor(abs(sin(i + 1)) * 2^32)`, the RFC 1321 table the source spells out as
/// literal round constants.
const MD5_SINE: [u32; 64] = [
    0xd76a_a478,
    0xe8c7_b756,
    0x2420_70db,
    0xc1bd_ceee,
    0xf57c_0faf,
    0x4787_c62a,
    0xa830_4613,
    0xfd46_9501,
    0x6980_98d8,
    0x8b44_f7af,
    0xffff_5bb1,
    0x895c_d7be,
    0x6b90_1122,
    0xfd98_7193,
    0xa679_438e,
    0x49b4_0821,
    0xf61e_2562,
    0xc040_b340,
    0x265e_5a51,
    0xe9b6_c7aa,
    0xd62f_105d,
    0x0244_1453,
    0xd8a1_e681,
    0xe7d3_fbc8,
    0x21e1_cde6,
    0xc337_07d6,
    0xf4d5_0d87,
    0x455a_14ed,
    0xa9e3_e905,
    0xfcef_a3f8,
    0x676f_02d9,
    0x8d2a_4c8a,
    0xfffa_3942,
    0x8771_f681,
    0x6d9d_6122,
    0xfde5_380c,
    0xa4be_ea44,
    0x4bde_cfa9,
    0xf6bb_4b60,
    0xbebf_bc70,
    0x289b_7ec6,
    0xeaa1_27fa,
    0xd4ef_3085,
    0x0488_1d05,
    0xd9d4_d039,
    0xe6db_99e5,
    0x1fa2_7cf8,
    0xc4ac_5665,
    0xf429_2244,
    0x432a_ff97,
    0xab94_23a7,
    0xfc93_a039,
    0x655b_59c3,
    0x8f0c_cc92,
    0xffef_f47d,
    0x8584_5dd1,
    0x6fa8_7e4f,
    0xfe2c_e6e0,
    0xa301_4314,
    0x4e08_11a1,
    0xf753_7e82,
    0xbd3a_f235,
    0x2ad7_d2bb,
    0xeb86_d391,
];

/// The RFC 1321 per-round left-rotation amounts.
const MD5_SHIFT: [u32; 64] = [
    7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9,
    14, 20, 5, 9, 14, 20, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 6, 10, 15,
    21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21,
];

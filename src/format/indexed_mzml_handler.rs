// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Random access to a single spectrum or chromatogram of an indexed mzML file.
//!
//! Ports `FORMAT/HANDLERS/IndexedMzMLHandler.h` and the record-decoding half of
//! `FORMAT/HANDLERS/MzMLSpectrumDecoder.h`. See `docs/INDEXED_MZML_HANDLER_SUPPORT.md`.
//!
//! [`IndexedMzMLHandler`](crate::format::indexed_mzml_handler::IndexedMzMLHandler)
//! parses the footer index with
//! [`IndexedMzMLDecoder`](crate::format::indexed_mzml::IndexedMzMLDecoder), then
//! reads the byte range of one record, wraps it in the document's own cached mzML
//! header and hands the result to the existing reader in
//! [`mzml`](crate::format::mzml). There is no second XML parser: the record is
//! decoded by exactly the code that reads a whole file, so RT, MS level,
//! precursors, products, scan settings and auxiliary arrays are all available.
//! The source decodes only `binaryDataArray` payloads plus the `id` attribute and
//! leaves every other field of the returned record default-constructed.
//!
//! The source keeps one file handle and moves it per access, and its class
//! documentation warns that this is not thread-safe. Every fetch here takes
//! `&mut self`, so the compiler enforces exclusive access; a second reader is a
//! second [`IndexedMzMLHandler`](crate::format::indexed_mzml_handler::IndexedMzMLHandler)
//! on the same path. Nothing in this module starts a thread.
//!
//! A malformed or hostile index cannot cause an unbounded read: every byte range
//! is checked against the file length and against
//! [`RecordReadLimits`](crate::format::indexed_mzml_handler::RecordReadLimits)
//! before any buffer is allocated.

use crate::format::indexed_mzml::{IndexReadLimits, IndexedMzMLDecoder};
use crate::format::mzml::{self, LoadOptions, ReadOptions};
use crate::format::peak_options::PeakFileOptions;
use crate::kernel::{MSChromatogram, MSSpectrum};
use crate::{Error, Result};
use quick_xml::{Reader, events::Event};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

/// Which of the two record lists of an mzML run a byte offset belongs to.
///
/// The source keeps the distinction in two parallel offset vectors and two
/// `getXById` overload families; one enumeration replaces both.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum RecordKind {
    /// A `<spectrum>` element, indexed by the `name="spectrum"` index section.
    Spectrum,
    /// A `<chromatogram>` element, indexed by the `name="chromatogram"` section.
    Chromatogram,
}

impl RecordKind {
    /// The mzML element name of one record of this kind.
    pub const fn element(self) -> &'static str {
        match self {
            Self::Spectrum => "spectrum",
            Self::Chromatogram => "chromatogram",
        }
    }

    /// The mzML element name of the list holding records of this kind.
    pub const fn list_element(self) -> &'static str {
        match self {
            Self::Spectrum => "spectrumList",
            Self::Chromatogram => "chromatogramList",
        }
    }

    const fn open_bytes(self) -> &'static [u8] {
        match self {
            Self::Spectrum => b"<spectrum",
            Self::Chromatogram => b"<chromatogram",
        }
    }

    const fn close_bytes(self) -> &'static [u8] {
        match self {
            Self::Spectrum => b"</spectrum",
            Self::Chromatogram => b"</chromatogram",
        }
    }

    const fn list_open_bytes(self) -> &'static [u8] {
        match self {
            Self::Spectrum => b"<spectrumList",
            Self::Chromatogram => b"<chromatogramList",
        }
    }
}

/// Explicit ceilings checked before any allocation or file read.
///
/// The source computes `endidx - startidx` from two index entries and passes the
/// difference straight to `new char[]`. A decreasing pair makes that length
/// negative and an out-of-range pair reads past the end of the file into an
/// uninitialised buffer. Each field below turns one of those into a checked
/// error, so an untrusted index costs a bounded amount of memory and nothing else.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct RecordReadLimits {
    /// Maximum bytes in the range spanned by one record, before it is read.
    pub max_record_bytes: usize,
    /// Maximum bytes of the cached document header preceding the first record list.
    pub max_header_bytes: usize,
    /// Maximum bytes scanned backwards from a first record to find its list tag.
    pub max_list_tag_bytes: usize,
    /// Maximum records per kind accepted from the index.
    pub max_records: usize,
    /// Cumulative bytes of stored native identifiers across both index sections,
    /// counting the ordered vector and the lookup map separately.
    pub max_native_id_bytes: usize,
    /// Maximum open elements in the cached header, bounding the synthesised
    /// closing-tag suffix.
    pub max_depth: usize,
    /// Limits handed to the footer decoder that produces the offsets.
    pub index: IndexReadLimits,
}

impl Default for RecordReadLimits {
    fn default() -> Self {
        Self {
            max_record_bytes: 256 << 20,
            max_header_bytes: 16 << 20,
            max_list_tag_bytes: 64 << 10,
            max_records: 1_000_000,
            max_native_id_bytes: 64 << 20,
            max_depth: 64,
            index: IndexReadLimits::default(),
        }
    }
}

/// The opening tag of one record list, as far as re-emitting it requires.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct ListTag {
    /// Raw, still-escaped `defaultDataProcessingRef` attribute value, if present.
    default_processing: Option<String>,
}

/// A low-level reader for one indexed mzML file.
///
/// Source `Internal::IndexedMzMLHandler`. It gives random access to the spectra
/// and chromatograms of an indexed mzML file without holding the whole file in
/// memory; for a container-shaped API use the facade that will sit on top of it,
/// as the source directs callers to `IndexedMzMLFileLoader` and
/// `OnDiscMSExperiment`.
///
/// Construction parses the footer index. The source instead reports failure
/// through `getParsingSuccess()` and documents that calling `getSpectrumById` or
/// `getChromatogramById` before checking it is invalid; here a handler only
/// exists when its index parsed, so that unchecked state cannot be reached.
///
/// Repeated `openFile` calls on one source object append to the offset vectors
/// instead of replacing them, because `parseFooter_` never clears them. There is
/// no reopen here: [`open`](IndexedMzMLHandler::open) always yields a fresh
/// handler, so the accumulation cannot happen.
///
/// # Examples
///
/// ```
/// use openms::format::indexed_mzml_handler::IndexedMzMLHandler;
///
/// let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
///     .join("tests/data/indexed_mzml/IndexedmzMLFile_1.mzML");
/// let mut handler = IndexedMzMLHandler::open(path)?;
/// assert_eq!(handler.spectrum_count(), 2);
/// assert_eq!(handler.chromatogram_count(), 1);
///
/// // Only this one record's byte range is read and decoded.
/// let spectrum = handler.spectrum(0)?.expect("no filter excludes it");
/// assert_eq!(spectrum.native_id, "controllerType=0 controllerNumber=1 scan=1");
/// assert_eq!(spectrum.peaks.len(), 19914);
/// assert_eq!(spectrum.ms_level, 1);
///
/// let chromatogram = handler.chromatogram_by_native_id("TIC")?.unwrap();
/// assert_eq!(chromatogram.peaks.len(), 48);
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Debug)]
pub struct IndexedMzMLHandler {
    path: PathBuf,
    file: File,
    file_length: u64,
    index_offset: u64,
    spectra: Vec<(String, u64)>,
    spectra_ids: BTreeMap<String, usize>,
    chromatograms: Vec<(String, u64)>,
    chromatogram_ids: BTreeMap<String, usize>,
    spectra_before_chromatograms: bool,
    header: Vec<u8>,
    closers: String,
    spectrum_list: Option<ListTag>,
    chromatogram_list: Option<ListTag>,
    limits: RecordReadLimits,
    read: ReadOptions,
    options: PeakFileOptions,
}

impl IndexedMzMLHandler {
    /// Open `path` and parse its footer index, with default limits.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when the file cannot be opened or read, and
    /// [`Error::Parse`] when it carries no usable `indexListOffset` footer, when
    /// the index does not decode, or when an index entry lies outside the file.
    /// The source constructor swallows all of these and records `false` in
    /// `parsing_success_`; only the numeric conversion of the footer offset
    /// escapes it, as `Exception::ConversionError`.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Self::open_with_limits(path, RecordReadLimits::default(), ReadOptions::default())
    }

    /// Open `path` with explicit resource ceilings and mzML decoding limits.
    ///
    /// # Errors
    ///
    /// As [`open`](Self::open), plus [`Error::InvalidValue`] when the index
    /// declares more records than `limits.max_records`, when the identifiers
    /// exceed `limits.max_native_id_bytes`, or when the document header before
    /// the first record list exceeds `limits.max_header_bytes`.
    pub fn open_with_limits(
        path: impl AsRef<Path>,
        limits: RecordReadLimits,
        read: ReadOptions,
    ) -> Result<Self> {
        let path = path.as_ref();
        let decoder = IndexedMzMLDecoder {
            limits: limits.index,
        };
        let index_offset = decoder.find_index_list_offset(path)?.ok_or_else(|| {
            parse("file has no indexListOffset footer; it is not an indexed mzML")
        })?;
        let offsets = decoder.parse_offsets(path, index_offset)?;
        let mut file = File::open(path)?;
        let file_length = file.seek(SeekFrom::End(0))?;
        if index_offset > file_length {
            return Err(parse("indexListOffset points past the end of the file"));
        }

        let spectra = offsets.spectra;
        let chromatograms = offsets.chromatograms;
        if spectra.len() > limits.max_records || chromatograms.len() > limits.max_records {
            return Err(Error::InvalidValue(
                "indexed mzML record count exceeds the configured limit".into(),
            ));
        }
        let mut budget = limits.max_native_id_bytes;
        for (id, _) in spectra.iter().chain(chromatograms.iter()) {
            budget = budget
                .checked_sub(id.len().saturating_mul(2))
                .ok_or_else(|| {
                    Error::InvalidValue(
                        "indexed mzML native identifiers exceed the configured byte limit".into(),
                    )
                })?;
        }
        let spectra_ids = identifier_map(&spectra);
        let chromatogram_ids = identifier_map(&chromatograms);

        // Source `parseFooter_` defaults to spectra-first and only compares when
        // both sections are populated.
        let spectra_before_chromatograms = match (spectra.first(), chromatograms.first()) {
            (Some((_, s)), Some((_, c))) => s < c,
            _ => true,
        };

        let mut spectrum_list = None;
        let mut chromatogram_list = None;
        let mut header_end = None;
        for (kind, first) in [
            (RecordKind::Spectrum, spectra.first()),
            (RecordKind::Chromatogram, chromatograms.first()),
        ] {
            let Some(&(_, offset)) = first else { continue };
            if offset > file_length {
                return Err(parse("index entry points past the end of the file"));
            }
            let (start, tag) = locate_list_tag(&mut file, kind, offset, &limits)?;
            match kind {
                RecordKind::Spectrum => spectrum_list = Some(tag),
                RecordKind::Chromatogram => chromatogram_list = Some(tag),
            }
            header_end = Some(header_end.map_or(start, |previous: u64| previous.min(start)));
        }

        let (header, closers) = match header_end {
            Some(end) => {
                let bytes = read_range(&mut file, 0, end, limits.max_header_bytes, "header")?;
                let closers = closing_suffix(&bytes, &limits)?;
                (bytes, closers)
            }
            None => (Vec::new(), String::new()),
        };

        Ok(Self {
            path: path.to_path_buf(),
            file,
            file_length,
            index_offset,
            spectra,
            spectra_ids,
            chromatograms,
            chromatogram_ids,
            spectra_before_chromatograms,
            header,
            closers,
            spectrum_list,
            chromatogram_list,
            limits,
            read,
            options: PeakFileOptions::default(),
        })
    }

    /// The path this handler was opened on.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// The ceilings this handler enforces.
    pub fn limits(&self) -> RecordReadLimits {
        self.limits
    }

    /// The mzML decoding limits applied to each assembled single-record document.
    pub fn read_options(&self) -> ReadOptions {
        self.read
    }

    /// The scientific filtering options applied to each fetched record.
    ///
    /// Source `IndexedMzMLHandler` has no options at all; `OnDiscMSExperiment`
    /// owns the `PeakFileOptions` and applies them around the handler. They live
    /// here so that a filtered-out record costs no peak decoding.
    pub fn options(&self) -> &PeakFileOptions {
        &self.options
    }

    /// Mutable access to the filtering options, as source `getOptions()`.
    pub fn options_mut(&mut self) -> &mut PeakFileOptions {
        &mut self.options
    }

    /// Replace the filtering options, as source `setOptions()`.
    pub fn set_options(&mut self, options: PeakFileOptions) {
        self.options = options;
    }

    /// The byte offset of the `<indexList>` element, as source `index_offset_`.
    pub fn index_list_offset(&self) -> u64 {
        self.index_offset
    }

    /// Whether the spectrum index precedes the chromatogram index in the file.
    ///
    /// Source `spectra_before_chroms_`: it decides where the final record of each
    /// kind ends. With either section empty the source assumes `true`, and so
    /// does this.
    pub fn spectra_before_chromatograms(&self) -> bool {
        self.spectra_before_chromatograms
    }

    /// The number of spectra in the index, as source `getNrSpectra()`.
    pub fn spectrum_count(&self) -> usize {
        self.spectra.len()
    }

    /// The number of chromatograms in the index, as source `getNrChromatograms()`.
    pub fn chromatogram_count(&self) -> usize {
        self.chromatograms.len()
    }

    /// Whether the index lists no records of either kind.
    pub fn is_empty(&self) -> bool {
        self.spectra.is_empty() && self.chromatograms.is_empty()
    }

    /// The number of records of `kind` in the index.
    pub fn count(&self, kind: RecordKind) -> usize {
        self.offsets(kind).len()
    }

    /// The native identifier the index records for `index`, or `None` when the
    /// index is out of range.
    pub fn native_id(&self, kind: RecordKind, index: usize) -> Option<&str> {
        self.offsets(kind).get(index).map(|(id, _)| id.as_str())
    }

    /// The byte offset the index records for `index`, or `None` when the index is
    /// out of range.
    pub fn offset(&self, kind: RecordKind, index: usize) -> Option<u64> {
        self.offsets(kind).get(index).map(|&(_, offset)| offset)
    }

    /// The position of `native_id` in the index of `kind`, or `None` when no
    /// entry carries it.
    ///
    /// The source builds this map with `unordered_map::emplace`, which keeps the
    /// first entry when a native identifier repeats; so does this. Lookups here
    /// are ordered and therefore reproducible run to run.
    pub fn index_of(&self, kind: RecordKind, native_id: &str) -> Option<usize> {
        match kind {
            RecordKind::Spectrum => self.spectra_ids.get(native_id),
            RecordKind::Chromatogram => self.chromatogram_ids.get(native_id),
        }
        .copied()
    }

    /// The raw XML of one record, from its `<spectrum`/`<chromatogram` start tag
    /// through its matching closing tag.
    ///
    /// Source `getSpectrumById_helper_` and `getChromatogramById_helper_` return
    /// everything between this record's offset and the next one, which for the
    /// last record of a kind runs to the start of the other list or to
    /// `<indexList>`. That tail carries `</spectrumList>` and the next list's
    /// opening tag, and the source's DOM parser is constructed without an error
    /// handler, so the resulting well-formedness errors are discarded silently.
    /// This trims the range at the record's own closing tag instead, so trailing
    /// content is never fed to a parser.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `index` is out of range or the range exceeds
    /// [`RecordReadLimits::max_record_bytes`]; [`Error::Parse`] when the index
    /// entries do not increase, when the range leaves the file, when the range
    /// does not start with the expected element, or when the element is not
    /// closed inside it; [`Error::Io`] when the read fails.
    pub fn record_xml(&mut self, kind: RecordKind, index: usize) -> Result<Vec<u8>> {
        let (start, end) = self.record_range(kind, index)?;
        let raw = read_range(
            &mut self.file,
            start,
            end,
            self.limits.max_record_bytes,
            "record",
        )?;
        let length = record_length(&raw, kind)?;
        let mut bytes = raw;
        bytes.truncate(length);
        Ok(bytes)
    }

    /// Read the spectrum at `index`, as source `getMSSpectrumById(int)`.
    ///
    /// Returns `Ok(None)` when [`options`](Self::options) exclude the record as a
    /// whole — an MS level outside the configured set, a retention time outside
    /// the RT range, a precursor m/z outside the precursor range, or
    /// `metadata_only`. m/z and intensity ranges select peaks, so a record that
    /// keeps none of them is returned as `Some` with no peaks.
    ///
    /// Source `OnDiscMSExperiment::getSpectrum` instead returns the metadata-only
    /// spectrum for an excluded record, explicitly so that the caller's index
    /// mapping survives; `None` says the same thing without inventing a record,
    /// and the index mapping here is the caller's `index` argument.
    ///
    /// # Errors
    ///
    /// As [`record_xml`](Self::record_xml), plus any error the mzML reader raises
    /// on the assembled single-record document, and [`Error::Parse`] when the
    /// decoded record's `id` attribute disagrees with the identifier the index
    /// recorded for this offset. The source checks neither.
    pub fn spectrum(&mut self, index: usize) -> Result<Option<MSSpectrum>> {
        let expected = self.expected_id(RecordKind::Spectrum, index)?;
        let mut experiment = self.decode(RecordKind::Spectrum, index)?;
        let Some(spectrum) = experiment.spectra.pop() else {
            return Ok(None);
        };
        check_identity(&expected, &spectrum.native_id)?;
        Ok(Some(spectrum))
    }

    /// Read the chromatogram at `index`, as source `getMSChromatogramById(int)`.
    ///
    /// Returns `Ok(None)` when [`options`](Self::options) exclude chromatograms
    /// or request metadata only; the RT and intensity ranges select points within
    /// the chromatogram rather than dropping it.
    ///
    /// # Errors
    ///
    /// As [`spectrum`](Self::spectrum).
    pub fn chromatogram(&mut self, index: usize) -> Result<Option<MSChromatogram>> {
        let expected = self.expected_id(RecordKind::Chromatogram, index)?;
        let mut experiment = self.decode(RecordKind::Chromatogram, index)?;
        let Some(chromatogram) = experiment.chromatograms.pop() else {
            return Ok(None);
        };
        check_identity(&expected, &chromatogram.native_id)?;
        Ok(Some(chromatogram))
    }

    /// Read the spectrum whose native identifier is `native_id`, as source
    /// `getMSSpectrumByNativeId`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no index entry carries `native_id`, matching
    /// the source's `Exception::IllegalArgument`; otherwise as
    /// [`spectrum`](Self::spectrum).
    ///
    /// The source's copy constructor does not copy either native-id map, so every
    /// lookup on a copied handler fails with that same exception. This port has
    /// no copy constructor and the maps are ordinary owned state.
    pub fn spectrum_by_native_id(&mut self, native_id: &str) -> Result<Option<MSSpectrum>> {
        let index = self
            .index_of(RecordKind::Spectrum, native_id)
            .ok_or_else(|| missing(RecordKind::Spectrum, native_id))?;
        self.spectrum(index)
    }

    /// Read the chromatogram whose native identifier is `native_id`, as source
    /// `getMSChromatogramByNativeId`.
    ///
    /// # Errors
    ///
    /// As [`spectrum_by_native_id`](Self::spectrum_by_native_id).
    pub fn chromatogram_by_native_id(&mut self, native_id: &str) -> Result<Option<MSChromatogram>> {
        let index = self
            .index_of(RecordKind::Chromatogram, native_id)
            .ok_or_else(|| missing(RecordKind::Chromatogram, native_id))?;
        self.chromatogram(index)
    }

    fn offsets(&self, kind: RecordKind) -> &[(String, u64)] {
        match kind {
            RecordKind::Spectrum => &self.spectra,
            RecordKind::Chromatogram => &self.chromatograms,
        }
    }

    fn expected_id(&self, kind: RecordKind, index: usize) -> Result<String> {
        self.native_id(kind, index)
            .map(str::to_owned)
            .ok_or_else(|| out_of_range(kind, index, self.count(kind)))
    }

    /// Byte range `[start, end)` of one record, as the source's two helpers.
    fn record_range(&self, kind: RecordKind, index: usize) -> Result<(u64, u64)> {
        let offsets = self.offsets(kind);
        let &(_, start) = offsets
            .get(index)
            .ok_or_else(|| out_of_range(kind, index, offsets.len()))?;
        let end = match offsets.get(index + 1) {
            Some(&(_, next)) => next,
            None => self.tail_end(kind),
        };
        if end < start {
            return Err(parse("index offsets do not increase across the record"));
        }
        if end > self.file_length {
            return Err(parse("record byte range extends past the end of the file"));
        }
        if end - start > self.limits.max_record_bytes as u64 {
            return Err(Error::InvalidValue(
                "record byte range exceeds the configured limit".into(),
            ));
        }
        Ok((start, end))
    }

    /// Where the last record of `kind` ends: at the other list when it follows,
    /// otherwise at `<indexList>`. This mirrors the source branch exactly.
    fn tail_end(&self, kind: RecordKind) -> u64 {
        let (other, other_follows) = match kind {
            RecordKind::Spectrum => (&self.chromatograms, self.spectra_before_chromatograms),
            RecordKind::Chromatogram => (&self.spectra, !self.spectra_before_chromatograms),
        };
        match other.first() {
            Some(&(_, offset)) if other_follows => offset,
            _ => self.index_offset,
        }
    }

    /// Assemble a one-record mzML document and hand it to the whole-file reader.
    fn decode(&mut self, kind: RecordKind, index: usize) -> Result<crate::kernel::MSExperiment> {
        let record = self.record_xml(kind, index)?;
        let list = match kind {
            RecordKind::Spectrum => self.spectrum_list.as_ref(),
            RecordKind::Chromatogram => self.chromatogram_list.as_ref(),
        }
        .ok_or_else(|| parse("no record list opening tag was located for this record kind"))?;

        let mut document = Vec::new();
        let size = self.header.len()
            + record.len()
            + self.closers.len()
            + list.default_processing.as_ref().map_or(0, String::len)
            + 128;
        document
            .try_reserve(size)
            .map_err(|_| Error::InvalidValue("cannot allocate single-record document".into()))?;
        document.extend_from_slice(&self.header);
        document.extend_from_slice(b"<");
        document.extend_from_slice(kind.list_element().as_bytes());
        document.extend_from_slice(b" count=\"1\"");
        if let Some(reference) = &list.default_processing {
            document.extend_from_slice(b" defaultDataProcessingRef=\"");
            document.extend_from_slice(reference.as_bytes());
            document.extend_from_slice(b"\"");
        }
        document.extend_from_slice(b">");
        document.extend_from_slice(&record);
        document.extend_from_slice(b"</");
        document.extend_from_slice(kind.list_element().as_bytes());
        document.extend_from_slice(b">");
        document.extend_from_slice(self.closers.as_bytes());

        let load = LoadOptions {
            scientific: self.options.clone(),
            skip_spectra: false,
            ..LoadOptions::default()
        };
        mzml::read_with_load_options(document.as_slice(), &load, &self.read)
    }
}

fn parse(message: &str) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

fn out_of_range(kind: RecordKind, index: usize, count: usize) -> Error {
    Error::InvalidValue(format!(
        "{} index {index} is not below the indexed count {count}",
        kind.element()
    ))
}

fn missing(kind: RecordKind, native_id: &str) -> Error {
    Error::InvalidValue(format!(
        "no {} in the index carries native id {native_id:?}",
        kind.element()
    ))
}

fn check_identity(expected: &str, decoded: &str) -> Result<()> {
    if expected == decoded {
        Ok(())
    } else {
        Err(parse(
            "index native id does not match the record found at its offset",
        ))
    }
}

/// First entry wins, as the source's `unordered_map::emplace`.
fn identifier_map(offsets: &[(String, u64)]) -> BTreeMap<String, usize> {
    let mut map = BTreeMap::new();
    for (index, (id, _)) in offsets.iter().enumerate() {
        map.entry(id.clone()).or_insert(index);
    }
    map
}

/// Read `[start, end)` after checking it against `limit` and the open file.
fn read_range(file: &mut File, start: u64, end: u64, limit: usize, what: &str) -> Result<Vec<u8>> {
    let length = end
        .checked_sub(start)
        .ok_or_else(|| parse("byte range does not increase"))?;
    let length = usize::try_from(length)
        .ok()
        .filter(|&length| length <= limit)
        .ok_or_else(|| {
            Error::InvalidValue(format!(
                "indexed mzML {what} exceeds the configured byte limit"
            ))
        })?;
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(length)
        .map_err(|_| Error::InvalidValue(format!("cannot allocate indexed mzML {what}")))?;
    file.seek(SeekFrom::Start(start))?;
    file.take(length as u64).read_to_end(&mut bytes)?;
    if bytes.len() != length {
        return Err(parse("file ended inside the requested byte range"));
    }
    Ok(bytes)
}

fn whitespace(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\r' | b'\n')
}

/// Length of the record element at the start of `bytes`, including its closing
/// tag. Leading whitespace is rejected rather than skipped, because an index
/// offset is defined to point at the element itself.
fn record_length(bytes: &[u8], kind: RecordKind) -> Result<usize> {
    let open = kind.open_bytes();
    let delimited = bytes.len() > open.len()
        && bytes.starts_with(open)
        && (matches!(bytes[open.len()], b'>' | b'/') || whitespace(bytes[open.len()]));
    if !delimited {
        return Err(parse("index offset does not point at the expected element"));
    }
    let mut quote = None;
    let mut slash = false;
    for (position, &byte) in bytes.iter().enumerate().skip(open.len()) {
        match quote {
            Some(delimiter) if byte == delimiter => quote = None,
            Some(_) => {}
            None if matches!(byte, b'\'' | b'"') => quote = Some(byte),
            None if byte == b'>' => {
                if slash {
                    return Ok(position + 1);
                }
                return closing_tag(bytes, position + 1, kind);
            }
            None => slash = byte == b'/',
        }
    }
    Err(parse(
        "record start tag is not terminated inside its byte range",
    ))
}

/// Position just past the record's closing tag, found by a literal scan.
///
/// mzML forbids an unescaped `<` in attribute values and character data, so the
/// literal can only be the real closing tag outside a comment or CDATA section.
/// Inside one it would truncate early and the assembled document would then fail
/// to parse, which is a loud failure rather than a silent partial record.
fn closing_tag(bytes: &[u8], from: usize, kind: RecordKind) -> Result<usize> {
    let close = kind.close_bytes();
    let mut position = from;
    while position + close.len() <= bytes.len() {
        if !bytes[position..].starts_with(close) {
            position += 1;
            continue;
        }
        let mut end = position + close.len();
        while end < bytes.len() && whitespace(bytes[end]) {
            end += 1;
        }
        if bytes.get(end) == Some(&b'>') {
            return Ok(end + 1);
        }
        position += 1;
    }
    Err(parse("record is not closed inside its byte range"))
}

/// Find the `<spectrumList`/`<chromatogramList` tag that opens the list holding
/// the record at `first`, scanning backwards within a bounded window.
fn locate_list_tag(
    file: &mut File,
    kind: RecordKind,
    first: u64,
    limits: &RecordReadLimits,
) -> Result<(u64, ListTag)> {
    let window = limits.max_list_tag_bytes as u64;
    let start = first.saturating_sub(window);
    let bytes = read_range(file, start, first, limits.max_list_tag_bytes, "list tag")?;
    let open = kind.list_open_bytes();
    let mut found = None;
    let mut position = 0;
    while position + open.len() < bytes.len() {
        let next = bytes[position + open.len()];
        if bytes[position..].starts_with(open) && (matches!(next, b'>' | b'/') || whitespace(next))
        {
            found = Some(position);
        }
        position += 1;
    }
    let found = found.ok_or_else(|| {
        parse("no record list opening tag precedes the first record within the search window")
    })?;
    let end = tag_end(&bytes, found)?;
    Ok((start + found as u64, list_tag(&bytes[found..end])?))
}

/// Position just past the `>` that ends the start tag beginning at `from`.
fn tag_end(bytes: &[u8], from: usize) -> Result<usize> {
    let mut quote = None;
    for (position, &byte) in bytes.iter().enumerate().skip(from) {
        match quote {
            Some(delimiter) if byte == delimiter => quote = None,
            Some(_) => {}
            None if matches!(byte, b'\'' | b'"') => quote = Some(byte),
            None if byte == b'>' => return Ok(position + 1),
            None => {}
        }
    }
    Err(parse("record list opening tag is not terminated"))
}

/// Extract the one attribute of a record list that the assembled document needs.
///
/// The value is re-emitted verbatim between double quotes, so it must not contain
/// a double quote or a `<`. mzML declares it `xsd:IDREF`, which permits neither.
fn list_tag(bytes: &[u8]) -> Result<ListTag> {
    let mut reader = Reader::from_reader(bytes);
    reader.config_mut().check_end_names = false;
    let mut buffer = Vec::new();
    let element = match reader.read_event_into(&mut buffer) {
        Ok(Event::Start(element)) => element,
        Ok(Event::Empty(element)) => element,
        Ok(_) => return Err(parse("record list opening tag is not an element")),
        Err(error) => return Err(parse(&error.to_string())),
    };
    let mut default_processing = None;
    for attribute in element.attributes().with_checks(false) {
        let attribute = attribute.map_err(|error| parse(&error.to_string()))?;
        if attribute.key.as_ref() != b"defaultDataProcessingRef" {
            continue;
        }
        let value = std::str::from_utf8(&attribute.value)
            .map_err(|_| parse("defaultDataProcessingRef is not UTF-8"))?;
        if value.is_empty()
            || value
                .bytes()
                .any(|byte| matches!(byte, b'<' | b'"' | b'&') || byte < 0x20)
        {
            return Err(parse(
                "defaultDataProcessingRef is empty or not a plain IDREF",
            ));
        }
        default_processing = Some(value.to_owned());
    }
    Ok(ListTag { default_processing })
}

fn element_name(raw: &[u8]) -> Result<String> {
    std::str::from_utf8(raw)
        .map(str::to_owned)
        .map_err(|_| parse("mzML header element name is not UTF-8"))
}

/// Closing tags for every element still open at the end of the cached header,
/// innermost first, so that header plus one record plus this suffix is a
/// complete document.
fn closing_suffix(header: &[u8], limits: &RecordReadLimits) -> Result<String> {
    let mut reader = Reader::from_reader(header);
    reader.config_mut().check_end_names = false;
    let mut buffer = Vec::new();
    let mut stack: Vec<String> = Vec::new();
    loop {
        let event = reader
            .read_event_into(&mut buffer)
            .map_err(|error| parse(&error.to_string()))?;
        match event {
            Event::Start(element) => {
                if stack.len() >= limits.max_depth {
                    return Err(parse("mzML header nesting exceeds the configured depth"));
                }
                stack.push(element_name(element.name().as_ref())?);
            }
            Event::End(element) => {
                let name = element_name(element.name().as_ref())?;
                if stack.pop() != Some(name) {
                    return Err(parse("mzML header has a mismatched closing tag"));
                }
            }
            Event::Eof => break,
            _ => {}
        }
        buffer.clear();
    }
    if stack.is_empty() {
        return Err(parse(
            "mzML header closes every element before the record list",
        ));
    }
    let mut closers = String::new();
    for name in stack.iter().rev() {
        closers.push_str("</");
        closers.push_str(name);
        closers.push('>');
    }
    Ok(closers)
}

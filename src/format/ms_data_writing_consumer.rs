// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Streaming mzML writer that consumes spectra and chromatograms one at a time.
//!
//! Ports `FORMAT/DATAACCESS/MSDataWritingConsumer.h` and its `.cpp`. See
//! `docs/MS_DATA_WRITING_CONSUMER_SUPPORT.md` for the full API mapping.
//!
//! The source class is an abstract `MzMLHandler` subclass that also implements
//! `Interfaces::IMSDataConsumer`: it writes the mzML header when the first
//! record arrives and appends each following record immediately, so an
//! experiment never has to be held in memory. Its two pure virtual hooks let a
//! derived class transform a record before it is written.
//!
//! [`MSDataWritingConsumer`](crate::format::ms_data_writing_consumer::MSDataWritingConsumer)
//! is that driver. The template-method pair becomes the
//! [`MSDataWritingProcessor`](crate::format::ms_data_writing_consumer::MSDataWritingProcessor)
//! trait,
//! [`PlainProcessor`](crate::format::ms_data_writing_consumer::PlainProcessor)
//! reproduces `PlainMSDataWritingConsumer`, and
//! [`NoopMSDataWritingConsumer`](crate::format::ms_data_writing_consumer::NoopMSDataWritingConsumer)
//! reproduces the do-nothing variant. The consumer implements
//! [`MSDataConsumer`](crate::interfaces::MSDataConsumer), so
//! [`transform`](crate::format::mzml::transform) can drive it directly.
//!
//! Record XML comes from [`crate::format::mzml`]'s writer, one record at a
//! time, exactly as the source delegates to `MzMLHandler::writeSpectrum_` and
//! `writeChromatogram_`. That keeps one definition of the mzML encoding
//! instead of a second, drifting one; the cost and the guards are described at
//! [`consume_spectrum`](crate::format::ms_data_writing_consumer::MSDataWritingConsumer::consume_spectrum).

use crate::format::mzml::{self, IndexedOutput, WriteOptions};
use crate::interfaces::MSDataConsumer;
use crate::kernel::{MSChromatogram, MSExperiment, MSSpectrum};
use crate::metadata::{DataProcessing, ExperimentalSettings};
use crate::{Error, Result};
use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::ops::ControlFlow;
use std::path::Path;
use std::sync::Arc;

/// The `spectrumList` opening tag the whole-document writer emits for exactly
/// one spectrum, and the marker a rendered record is split on.
const SPECTRUM_LIST_OPEN: &str =
    "<spectrumList count=\"1\" defaultDataProcessingRef=\"dp_00000000000000000000\">\n";
/// The matching closing tag.
const SPECTRUM_LIST_CLOSE: &str = "</spectrumList>\n";
/// The `chromatogramList` opening tag for exactly one chromatogram.
const CHROMATOGRAM_LIST_OPEN: &str =
    "<chromatogramList count=\"1\" defaultDataProcessingRef=\"dp_00000000000000000000\">\n";
/// The matching closing tag.
const CHROMATOGRAM_LIST_CLOSE: &str = "</chromatogramList>\n";
/// The document's final tags, written once by the equivalent of `doCleanup_`.
const DOCUMENT_CLOSE: &str = "</run></mzML>\n";
/// The XML declaration every rendered document opens with. The streaming
/// consumer strips it from the rendered header and lets
/// [`mzml::IndexedOutput::header`] write it, so that the `indexedmzML` opening
/// tag lands between the declaration and `<mzML`, where the whole-document
/// indexed writer puts it.
const XML_DECLARATION: &str = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
/// The default processing reference both list tags carry.
const DEFAULT_PROCESSING: &str = "dp_00000000000000000000";
/// The `index` attribute of a record rendered on its own, always zero, which
/// the running record number replaces.
const INDEX_ZERO: &str = " index=\"0\" defaultArrayLength=\"";

/// Attribute prefixes with which a record element references a header element.
const RECORD_REFERENCES: [&str; 3] = [
    " dataProcessingRef=\"",
    " sourceFileRef=\"",
    " instrumentConfigurationRef=\"",
];
/// Header element openings whose `id` a record element may reference.
const HEADER_DECLARATIONS: [&str; 3] = [
    "<dataProcessing id=\"",
    "<sourceFile id=\"",
    "<instrumentConfiguration id=\"",
];
/// The header lists whose entries a record's references are numbered against.
///
/// An identifier a record emits is an index into one of these, so two records
/// may only share one header when both lists come out identical. Comparing the
/// rendered text is exact and needs no knowledge of how the indices are
/// assigned. `fileContent` and `instrumentConfigurationList` are deliberately
/// not compared: the first is descriptive and legitimately differs between an
/// MS1 and an MS2 record, and the second comes from the frozen experimental
/// settings rather than from the record.
const DECLARATION_BLOCKS: [(&str, &str); 2] = [
    ("<sourceFileList", "</sourceFileList>"),
    ("<dataProcessingList", "</dataProcessingList>"),
];

fn limit(what: &str) -> Error {
    Error::InvalidValue(format!("mzML writing consumer {what} limit exceeded"))
}

fn layout(what: &str) -> Error {
    Error::Unsupported(format!(
        "the mzML writer's {what} is not in the layout this streaming consumer splits on"
    ))
}

/// Whether the list `count` attributes must agree with the records written.
///
/// The source's class `@note` states that the expected size "will not be
/// enforced but it will lead to an inconsistent mzML if the count attribute of
/// spectrumList or chromatogramList is incorrect". The crate's policy for a
/// source behaviour that loses or corrupts information is to refuse by default
/// and to offer the source behaviour explicitly, as `dta::WriteOptions`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum CountPolicy {
    /// [`MSDataWritingConsumer::finish`] fails when the announced counts and
    /// the written counts differ. The document is still closed first, so the
    /// file left behind is well-formed XML with wrong `count` attributes.
    #[default]
    Checked,
    /// Accept the mismatch, as the source does.
    SourceInconsistent,
}

/// Administrative ceilings for one consumer, charged before anything is written.
///
/// The source enforces nothing: `consumeSpectrum` copies, processes and writes
/// whatever arrives. These ceilings exist because every value here can come
/// from a file the caller did not produce.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WritingLimits {
    /// Maximum spectra one consumer may write.
    pub max_spectra: usize,
    /// Maximum chromatograms one consumer may write.
    pub max_chromatograms: usize,
    /// Maximum bytes of the intermediate document rendered for one record.
    pub max_record_bytes: usize,
    /// Maximum bytes of native identifiers retained for duplicate detection.
    pub max_native_id_bytes: usize,
}

impl WritingLimits {
    /// Default record count for each kind, one million.
    pub const MAX_RECORDS: usize = 1_000_000;
    /// Default intermediate document size for one record: 256 MiB.
    pub const MAX_RECORD_BYTES: usize = 256 * 1024 * 1024;
    /// Default retained identifier text: 64 MiB.
    pub const MAX_NATIVE_ID_BYTES: usize = 64 * 1024 * 1024;
}

impl Default for WritingLimits {
    fn default() -> Self {
        Self {
            max_spectra: Self::MAX_RECORDS,
            max_chromatograms: Self::MAX_RECORDS,
            max_record_bytes: Self::MAX_RECORD_BYTES,
            max_native_id_bytes: Self::MAX_NATIVE_ID_BYTES,
        }
    }
}

/// Per-record transformation applied before a record is written.
///
/// Ports the source's private pure virtual pair `processSpectrum_` and
/// `processChromatogram_`, the template-method hooks whose documentation says
/// "Redefine this function to determine spectra processing before writing to
/// disk". Both methods are required, so an implementation cannot silently
/// inherit a no-op it did not intend; [`PlainProcessor`] is the explicit
/// no-op.
///
/// The record handed over is the consumer's own copy, so a processor may
/// rewrite it freely without the caller observing the change - the source
/// makes the same copy in `SpectrumType scpy = s`.
pub trait MSDataWritingProcessor {
    /// Transform a spectrum before it is written.
    ///
    /// # Errors
    ///
    /// Any error aborts the record before a byte of it is written. The
    /// source's hook returns `void` and can only report a problem by throwing.
    fn process_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<()>;
    /// Transform a chromatogram before it is written.
    ///
    /// # Errors
    ///
    /// As [`Self::process_spectrum`].
    fn process_chromatogram(&mut self, chromatogram: &mut MSChromatogram) -> Result<()>;
}

/// A processor that changes nothing, as `PlainMSDataWritingConsumer` does.
///
/// The source class documentation calls it "probably the class you want if you
/// want to write mzML files to disk by providing spectra and chromatograms
/// sequentially"; [`PlainMSDataWritingConsumer`] is that combination.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct PlainProcessor;

impl MSDataWritingProcessor for PlainProcessor {
    fn process_spectrum(&mut self, _spectrum: &mut MSSpectrum) -> Result<()> {
        Ok(())
    }
    fn process_chromatogram(&mut self, _chromatogram: &mut MSChromatogram) -> Result<()> {
        Ok(())
    }
}

/// A consumer that writes records as they arrive, without holding them.
///
/// Ports `MSDataWritingConsumer`, including its state machine: nothing is
/// written until the first record arrives, the header is then rendered from
/// the experimental settings plus that record, each record is appended
/// immediately, and [`finish`](Self::finish) closes the open list and the
/// document.
///
/// # Examples
///
/// ```
/// use openms::format::ms_data_writing_consumer::PlainMSDataWritingConsumer;
/// use openms::kernel::MSSpectrum;
/// use openms::metadata::ExperimentalSettings;
///
/// let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
/// consumer.set_experimental_settings(&ExperimentalSettings::default())?;
/// consumer.set_expected_size(1, 0)?;
/// let mut spectrum = MSSpectrum {
///     native_id: "scan=1".into(),
///     ..Default::default()
/// };
/// consumer.consume_spectrum(&mut spectrum)?;
/// assert_eq!(consumer.spectra_written(), 1);
/// let bytes = consumer.finish()?;
/// let text = String::from_utf8(bytes).unwrap();
/// assert!(text.contains("<spectrumList count=\"1\""));
/// assert!(text.ends_with("</run></mzML>\n"));
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Notes
///
/// The source's three class-level `@note` paragraphs all hold here. The first
/// record starts the header. Spectra may not follow chromatograms, because an
/// mzML file cannot carry two `spectrumList` elements. And the expected sizes
/// are written into the list `count` attributes without being enforced, which
/// [`CountPolicy`] turns from a silent inconsistency into a reported one.
#[derive(Debug)]
pub struct MSDataWritingConsumer<W: Write, P: MSDataWritingProcessor = PlainProcessor> {
    writer: IndexedOutput<'static, W>,
    processor: P,
    /// Identifier and byte offset of every spectrum written, for the index.
    spectrum_ids: Vec<(String, u64)>,
    /// Identifier and byte offset of every chromatogram written.
    chromatogram_ids: Vec<(String, u64)>,
    options: WriteOptions,
    limits: WritingLimits,
    counts: CountPolicy,
    settings: ExperimentalSettings,
    additional_data_processing: Option<Arc<DataProcessing>>,
    declared: BTreeSet<String>,
    declarations: String,
    native_ids: BTreeSet<String>,
    native_id_bytes: usize,
    started_writing: bool,
    writing_spectra: bool,
    writing_chromatograms: bool,
    spectra_written: usize,
    chromatograms_written: usize,
    spectra_expected: usize,
    chromatograms_expected: usize,
}

/// A consumer that writes records unchanged, the source's
/// `PlainMSDataWritingConsumer`.
pub type PlainMSDataWritingConsumer<W> = MSDataWritingConsumer<W, PlainProcessor>;

impl<W: Write, P: MSDataWritingProcessor> MSDataWritingConsumer<W, P> {
    /// A consumer writing to `writer` and transforming records with `processor`.
    ///
    /// The source constructor takes only a filename because it owns its
    /// `std::ofstream`; it opens the file in binary mode so that no line
    /// endings are rewritten, which a Rust writer does not do in the first
    /// place. [`create`](Self::create) is the filename form.
    pub fn new(writer: W, processor: P) -> Self {
        Self {
            // The offset tables of the wrapped output stay empty: this consumer
            // learns its records one at a time and keeps its own tables, filled
            // from `IndexedOutput::position`.
            writer: IndexedOutput::streamed_with_capacity(writer, 0, 0)
                .expect("two empty offset tables need no allocation"),
            processor,
            spectrum_ids: Vec::new(),
            chromatogram_ids: Vec::new(),
            options: WriteOptions::default(),
            limits: WritingLimits::default(),
            counts: CountPolicy::default(),
            settings: ExperimentalSettings::default(),
            additional_data_processing: None,
            declared: BTreeSet::new(),
            declarations: String::new(),
            native_ids: BTreeSet::new(),
            native_id_bytes: 0,
            started_writing: false,
            writing_spectra: false,
            writing_chromatograms: false,
            spectra_written: 0,
            chromatograms_written: 0,
            spectra_expected: 0,
            chromatograms_expected: 0,
        }
    }

    /// Replace the mzML writer options, which select binary array compression.
    ///
    /// Native: the source consumer inherits `MzMLHandler`'s `PeakFileOptions`
    /// and no TOPP consumer reaches them through this class.
    pub fn with_write_options(mut self, options: WriteOptions) -> Self {
        self.options = options;
        self
    }

    /// Replace the administrative ceilings of [`WritingLimits`].
    pub fn with_limits(mut self, limits: WritingLimits) -> Self {
        self.limits = limits;
        self
    }

    /// Replace the list-count policy; see [`CountPolicy`].
    pub fn with_count_policy(mut self, counts: CountPolicy) -> Self {
        self.counts = counts;
        self
    }

    /// The mzML writer options in force.
    pub fn write_options(&self) -> WriteOptions {
        self.options
    }

    /// The administrative ceilings in force.
    pub fn limits(&self) -> WritingLimits {
        self.limits
    }

    /// The list-count policy in force.
    pub fn count_policy(&self) -> CountPolicy {
        self.counts
    }

    /// The per-record processor, for a caller that wants to read its state.
    pub fn processor(&self) -> &P {
        &self.processor
    }

    /// Set the experimental settings used for the whole file.
    ///
    /// Ports `setExperimentalSettings`. The settings and the first record
    /// together determine almost the whole mzML header, as the source's
    /// `@param` note says.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] once the header has been written, so a
    /// caller cannot believe settings took effect when they did not. The
    /// source overwrites `settings_` silently at any time and a late value is
    /// simply never used.
    pub fn set_experimental_settings(&mut self, settings: &ExperimentalSettings) -> Result<()> {
        self.apply_settings(settings)
    }

    /// The experimental settings currently held, the source's `settings_`.
    pub fn settings(&self) -> &ExperimentalSettings {
        &self.settings
    }

    /// Set the record counts written into the two list `count` attributes.
    ///
    /// Ports `setExpectedSize`. The source's own documentation warns that
    /// these "will contain wrong numbers if the expected size is not set
    /// correctly", and its class `@note` repeats that nothing enforces them.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when either count exceeds its ceiling
    /// in [`WritingLimits`], or once the corresponding list tag has already
    /// been written, where a later value could not take effect.
    pub fn set_expected_size(&mut self, spectra: usize, chromatograms: usize) -> Result<()> {
        self.expect_size(spectra, chromatograms)
    }

    /// The expected spectrum and chromatogram counts, in that order.
    pub fn expected_size(&self) -> (usize, usize) {
        (self.spectra_expected, self.chromatograms_expected)
    }

    /// Add a data-processing entry to every record written from now on.
    ///
    /// Ports `addDataProcessing`. The source wraps the argument in a fresh
    /// `DataProcessingPtr` and sets its `add_dataprocessing_` flag; every
    /// later record gets that pointer appended to its own history. Calling it
    /// twice replaces the entry rather than adding a second one, as upstream.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] once the header has been written. The
    /// header declares the `dataProcessingList`, so an entry added afterwards
    /// would make every following record reference an element that does not
    /// exist - see the note at [`Self::consume_spectrum`].
    pub fn add_data_processing(&mut self, processing: DataProcessing) -> Result<()> {
        if self.started_writing {
            return Err(Error::InvalidValue(
                "data processing cannot be added after the mzML header is written".into(),
            ));
        }
        self.additional_data_processing = Some(Arc::new(processing));
        Ok(())
    }

    /// The extra data-processing entry, when one was added.
    pub fn additional_data_processing(&self) -> Option<&Arc<DataProcessing>> {
        self.additional_data_processing.as_ref()
    }

    /// Spectra written so far, the source's `getNrSpectraWritten`.
    pub fn spectra_written(&self) -> usize {
        self.spectra_written
    }

    /// Chromatograms written so far, the source's `getNrChromatogramsWritten`.
    pub fn chromatograms_written(&self) -> usize {
        self.chromatograms_written
    }

    /// Whether the mzML header has been written, the source's `started_writing_`.
    pub fn started_writing(&self) -> bool {
        self.started_writing
    }

    /// Whether a `spectrumList` is open, the source's `writing_spectra_`.
    pub fn writing_spectra(&self) -> bool {
        self.writing_spectra
    }

    /// Whether a `chromatogramList` is open, the source's
    /// `writing_chromatograms_`.
    pub fn writing_chromatograms(&self) -> bool {
        self.writing_chromatograms
    }

    /// Process and write one spectrum.
    ///
    /// Ports `consumeSpectrum`. The spectrum is copied, the copy is handed to
    /// [`MSDataWritingProcessor::process_spectrum`], the optional extra
    /// data-processing entry is appended, the header and the `spectrumList`
    /// opening tag are written if this is the first record of their kind, and
    /// the spectrum element is appended. The caller's spectrum is never
    /// modified, matching the source's `SpectrumType scpy = s`.
    ///
    /// An empty `native_id` is filled in with `index=N` for the record's
    /// position, which is what the whole-document writer would have produced
    /// there. The source writes the empty identifier through, and a file whose
    /// records all carry `id=""` has duplicate identifiers.
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidValue`] when a chromatogram has already been written.
    ///   The source throws `Exception::IllegalArgument` with "Cannot write
    ///   spectra after writing chromatograms.", because two `spectrumList`
    ///   elements cannot appear in one mzML file.
    /// * [`Error::InvalidValue`] when [`WritingLimits::max_spectra`] would be
    ///   exceeded, or when the native identifier repeats an earlier one. The
    ///   source detects neither; the whole-document writer rejects duplicate
    ///   identifiers, and the streaming path must not be the weaker one.
    /// * [`Error::Unsupported`] when the record would need a `sourceFileList`
    ///   or `dataProcessingList` other than the one already written, or when
    ///   it references an element the header did not declare at all. Only the
    ///   first record contributes to the header, so a later record with a
    ///   different data-processing history or source file - its own, or one on
    ///   an auxiliary array - either has nothing to point at or points at the
    ///   wrong entry. The source emits the dangling reference instead:
    ///   `MzMLHandler::writeSpectrum_` builds `sourceFileRef="sf_sp_<n>"` and
    ///   `dataProcessingRef="dp_sp_<n>"` from the running record number while
    ///   its `dps` vector holds only the first record's history, and the
    ///   `// TODO ... assert this here` comment at
    ///   `MSDataWritingConsumer.cpp:93` marks the same gap.
    /// * [`Error::Io`] from a failed write, including an intermediate document
    ///   larger than [`WritingLimits::max_record_bytes`], which is refused
    ///   inside the rendering sink.
    /// * Any error from the processor or from the mzML writer.
    ///
    /// A rejected record leaves the file exactly as it was: the record is
    /// rendered into memory and every check runs before a byte reaches the
    /// writer. Bytes already written for earlier records are not rolled back,
    /// which no streaming writer can do.
    pub fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<()> {
        self.write_spectrum(spectrum)
    }

    /// Process and write one chromatogram.
    ///
    /// Ports `consumeChromatogram`, which first closes an open `spectrumList`
    /// and then follows the same path as [`Self::consume_spectrum`]. An empty
    /// `native_id` becomes `chromatogram=N`.
    ///
    /// # Errors
    ///
    /// As [`Self::consume_spectrum`], except that no ordering rule forbids a
    /// chromatogram: [`WritingLimits::max_chromatograms`] bounds the count,
    /// duplicate identifiers are rejected, an undeclared reference is
    /// [`Error::Unsupported`], and I/O, processor and writer failures
    /// propagate.
    pub fn consume_chromatogram(&mut self, chromatogram: &mut MSChromatogram) -> Result<()> {
        self.write_chromatogram(chromatogram)
    }

    /// Close the open list and the document, then return the writer.
    ///
    /// Ports `doCleanup_`, which the source destructor calls: it closes an
    /// open `spectrumList` or `chromatogramList`, writes the document footer
    /// only when writing actually started, and closes the file.
    ///
    /// Rust cannot report an error from a drop, so this replaces the
    /// destructor. Dropping a consumer without calling it leaves a file whose
    /// `run` and `mzML` elements are unclosed - deliberately, because the
    /// alternative is a silently swallowed I/O failure on a file the caller
    /// believes is complete. A consumer that received no record writes nothing
    /// at all, as upstream; [`Self::started_writing`] reports that case.
    ///
    /// # Errors
    ///
    /// Under [`CountPolicy::Checked`], returns [`Error::InvalidValue`] when
    /// the records written do not match the expected sizes - after the
    /// document has been closed, so the file is complete but its `count`
    /// attributes disagree with its contents, which the source produces
    /// silently. [`Error::Io`] propagates a failed write or flush.
    pub fn finish(mut self) -> Result<W> {
        if self.writing_spectra {
            self.writer.write_all(SPECTRUM_LIST_CLOSE.as_bytes())?;
            self.writing_spectra = false;
        } else if self.writing_chromatograms {
            self.writer.write_all(CHROMATOGRAM_LIST_CLOSE.as_bytes())?;
            self.writing_chromatograms = false;
        }
        if self.started_writing {
            self.writer.write_all(DOCUMENT_CLOSE.as_bytes())?;
            let (spectra, chromatograms) = (self.spectrum_ids, self.chromatogram_ids);
            self.writer.footer_ids(&spectra, &chromatograms)?;
        }
        self.writer.flush()?;
        let expected = (self.spectra_expected, self.chromatograms_expected);
        let written = (self.spectra_written, self.chromatograms_written);
        if self.counts == CountPolicy::Checked && self.started_writing && expected != written {
            return Err(Error::InvalidValue(format!(
                "mzML list counts announce {expected:?} records but {written:?} were written"
            )));
        }
        Ok(self.writer.into_inner())
    }

    fn apply_settings(&mut self, settings: &ExperimentalSettings) -> Result<()> {
        if self.started_writing {
            return Err(Error::InvalidValue(
                "experimental settings cannot change after the mzML header is written".into(),
            ));
        }
        self.settings = settings.clone();
        Ok(())
    }

    fn expect_size(&mut self, spectra: usize, chromatograms: usize) -> Result<()> {
        if spectra > self.limits.max_spectra || chromatograms > self.limits.max_chromatograms {
            return Err(limit("expected record count"));
        }
        if self.writing_spectra || self.writing_chromatograms {
            return Err(Error::InvalidValue(
                "expected size cannot change after a list tag is written".into(),
            ));
        }
        self.spectra_expected = spectra;
        self.chromatograms_expected = chromatograms;
        Ok(())
    }

    fn write_spectrum(&mut self, spectrum: &MSSpectrum) -> Result<()> {
        if self.writing_chromatograms {
            return Err(Error::InvalidValue(
                "cannot write spectra after writing chromatograms".into(),
            ));
        }
        if self.spectra_written >= self.limits.max_spectra {
            return Err(limit("spectrum count"));
        }
        let index = self.spectra_written;
        let mut copy = spectrum.clone();
        self.processor.process_spectrum(&mut copy)?;
        if let Some(processing) = &self.additional_data_processing {
            copy.data_processing.push(Arc::clone(processing));
        }
        if copy.native_id.is_empty() {
            copy.native_id = format!("index={index}");
        }
        self.remember_native_id(&copy.native_id)?;
        let id = copy.native_id.clone();
        let mut document = self.document();
        document.spectra.push(copy);
        let (head, block) = self.prepare(
            &document,
            index,
            SPECTRUM_LIST_OPEN,
            SPECTRUM_LIST_CLOSE,
            "spectrum",
        )?;
        if let Some(head) = head {
            self.write_head(&head)?;
        }
        if !self.writing_spectra {
            let count = self.spectra_expected;
            self.writer.write_all(
                format!(
                    "<spectrumList count=\"{count}\" defaultDataProcessingRef=\"{DEFAULT_PROCESSING}\">\n"
                )
                .as_bytes(),
            )?;
            self.writing_spectra = true;
        }
        self.spectrum_ids.push((id, self.writer.position()));
        self.writer.write_all(block.as_bytes())?;
        self.spectra_written = index.saturating_add(1);
        Ok(())
    }

    fn write_chromatogram(&mut self, chromatogram: &MSChromatogram) -> Result<()> {
        if self.chromatograms_written >= self.limits.max_chromatograms {
            return Err(limit("chromatogram count"));
        }
        let index = self.chromatograms_written;
        let mut copy = chromatogram.clone();
        self.processor.process_chromatogram(&mut copy)?;
        if let Some(processing) = &self.additional_data_processing {
            copy.data_processing.push(Arc::clone(processing));
        }
        if copy.native_id.is_empty() {
            copy.native_id = format!("chromatogram={index}");
        }
        self.remember_native_id(&copy.native_id)?;
        let id = copy.native_id.clone();
        let mut document = self.document();
        document.chromatograms.push(copy);
        let (head, block) = self.prepare(
            &document,
            index,
            CHROMATOGRAM_LIST_OPEN,
            CHROMATOGRAM_LIST_CLOSE,
            "chromatogram",
        )?;
        if let Some(head) = head {
            self.write_head(&head)?;
        }
        if self.writing_spectra {
            self.writer.write_all(SPECTRUM_LIST_CLOSE.as_bytes())?;
            self.writing_spectra = false;
        }
        if !self.writing_chromatograms {
            let count = self.chromatograms_expected;
            self.writer.write_all(
                format!(
                    "<chromatogramList count=\"{count}\" defaultDataProcessingRef=\"{DEFAULT_PROCESSING}\">\n"
                )
                .as_bytes(),
            )?;
            self.writing_chromatograms = true;
        }
        self.chromatogram_ids.push((id, self.writer.position()));
        self.writer.write_all(block.as_bytes())?;
        self.chromatograms_written = index.saturating_add(1);
        Ok(())
    }

    /// Write the document header, letting the indexed output put the
    /// `indexedmzML` opening tag between the XML declaration and `<mzML`,
    /// where the whole-document indexed writer puts it.
    ///
    /// The rendered header always begins with the declaration
    /// [`XML_DECLARATION`]; anything else means the writer's layout changed
    /// under this module and is refused rather than guessed at.
    fn write_head(&mut self, head: &str) -> Result<()> {
        let body = head
            .strip_prefix(XML_DECLARATION)
            .ok_or_else(|| layout("XML declaration"))?;
        self.writer.header(body)?;
        self.started_writing = true;
        Ok(())
    }


    /// A one-record document carrying the experimental settings.
    ///
    /// Every record is rendered against the same settings, not only the first.
    /// A record's `sourceFileRef` is an index into a `sourceFileList` that
    /// starts with the settings' own source files, so rendering a later record
    /// without them would number its reference from a different base and point
    /// it at the wrong entry. The header text this produces for records after
    /// the first is discarded, which is the cost of reusing one mzML encoder
    /// rather than maintaining a second; the source renders the header once.
    fn document(&self) -> MSExperiment {
        MSExperiment {
            settings: self.settings.clone(),
            ..MSExperiment::new()
        }
    }

    /// Render one record, split it out, renumber it and check its references.
    ///
    /// Returns the header text when this render is the one that establishes
    /// the header, and the record element either way.
    fn prepare(
        &mut self,
        document: &MSExperiment,
        index: usize,
        open: &str,
        close: &str,
        what: &str,
    ) -> Result<(Option<String>, String)> {
        let rendered = self.render(document)?;
        let (head, body) = split(&rendered, open, close, what)?;
        let block = stamp_index(body, index)?;
        let declarations = declaration_text(head);
        if self.started_writing {
            if declarations != self.declarations {
                return Err(Error::Unsupported(
                    "record needs a different mzML sourceFileList or dataProcessingList than the \
                     header written for the first record"
                        .into(),
                ));
            }
            check_references(&block, &self.declared)?;
            return Ok((None, block));
        }
        let declared = declared_ids(head);
        check_references(&block, &declared)?;
        let head = head.to_owned();
        self.declared = declared;
        self.declarations = declarations;
        Ok((Some(head), block))
    }

    /// Render a one-record document with the crate's whole-document writer.
    fn render(&self, document: &MSExperiment) -> Result<String> {
        let mut buffer = Bounded {
            data: Vec::new(),
            remaining: self.limits.max_record_bytes,
        };
        mzml::write_with_options(&mut buffer, document, &self.options)?;
        String::from_utf8(buffer.data)
            .map_err(|_| Error::Unsupported("the mzML writer produced non-UTF-8 output".into()))
    }

    /// Reject a repeated native identifier and charge the retained text.
    fn remember_native_id(&mut self, id: &str) -> Result<()> {
        let bytes = self
            .native_id_bytes
            .checked_add(id.len())
            .ok_or_else(|| limit("native identifier"))?;
        if bytes > self.limits.max_native_id_bytes {
            return Err(limit("native identifier"));
        }
        if !self.native_ids.insert(id.to_owned()) {
            return Err(Error::InvalidValue(format!(
                "duplicate mzML native identifier {id:?}"
            )));
        }
        self.native_id_bytes = bytes;
        Ok(())
    }
}

impl<W: Write, P: MSDataWritingProcessor> MSDataConsumer for MSDataWritingConsumer<W, P> {
    fn set_expected_size(&mut self, spectra: usize, chromatograms: usize) -> Result<()> {
        self.expect_size(spectra, chromatograms)
    }
    fn set_experimental_settings(&mut self, settings: &ExperimentalSettings) -> Result<()> {
        self.apply_settings(settings)
    }
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        self.write_spectrum(spectrum)?;
        Ok(ControlFlow::Continue(()))
    }
    fn consume_chromatogram(
        &mut self,
        chromatogram: &mut MSChromatogram,
    ) -> Result<ControlFlow<()>> {
        self.write_chromatogram(chromatogram)?;
        Ok(ControlFlow::Continue(()))
    }
}

impl<P: MSDataWritingProcessor> MSDataWritingConsumer<BufWriter<File>, P> {
    /// Create or truncate `path` and write mzML to it.
    ///
    /// Ports the source's only constructor,
    /// `MSDataWritingConsumer(const std::string& filename)`.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be created. The source opens its
    /// stream without checking `is_open`, so a consumer on an unwritable path
    /// reports nothing at all.
    pub fn create(path: impl AsRef<Path>, processor: P) -> Result<Self> {
        Ok(Self::new(
            BufWriter::new(File::create(path.as_ref())?),
            processor,
        ))
    }
}

impl<W: Write> PlainMSDataWritingConsumer<W> {
    /// A consumer writing records unchanged to `writer`.
    ///
    /// Ports `PlainMSDataWritingConsumer`, which is `MSDataWritingConsumer`
    /// with both processing hooks defined as empty bodies.
    pub fn plain(writer: W) -> Self {
        Self::new(writer, PlainProcessor)
    }
}

impl PlainMSDataWritingConsumer<BufWriter<File>> {
    /// A consumer writing records unchanged to a newly created `path`.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be created.
    pub fn create_plain(path: impl AsRef<Path>) -> Result<Self> {
        Self::create(path, PlainProcessor)
    }
}

/// A consumer that accepts records and does nothing with them.
///
/// Ports `NoopMSDataWritingConsumer`, whose class documentation explains it is
/// "sometimes necessary to fulfill the requirement of passing an valid
/// MSDataWritingConsumer object or pointer but no operation is required". The
/// source overrides every interface method with an empty body and overrides
/// `doCleanup_` so that even the destructor writes nothing - but it still
/// takes a filename, and its base constructor still creates and truncates that
/// file. This port takes no path and touches no file, so asking for a consumer
/// that does nothing cannot destroy an existing output.
///
/// The counters report what arrived, which the source's overrides discard.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct NoopMSDataWritingConsumer {
    spectra: usize,
    chromatograms: usize,
}

impl NoopMSDataWritingConsumer {
    /// A consumer that discards everything.
    pub fn new() -> Self {
        Self::default()
    }
    /// Spectra offered so far.
    pub fn spectra_written(&self) -> usize {
        self.spectra
    }
    /// Chromatograms offered so far.
    pub fn chromatograms_written(&self) -> usize {
        self.chromatograms
    }
}

impl MSDataConsumer for NoopMSDataWritingConsumer {
    fn set_expected_size(&mut self, _spectra: usize, _chromatograms: usize) -> Result<()> {
        Ok(())
    }
    fn set_experimental_settings(&mut self, _settings: &ExperimentalSettings) -> Result<()> {
        Ok(())
    }
    fn consume_spectrum(&mut self, _spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        self.spectra = self.spectra.saturating_add(1);
        Ok(ControlFlow::Continue(()))
    }
    fn consume_chromatogram(
        &mut self,
        _chromatogram: &mut MSChromatogram,
    ) -> Result<ControlFlow<()>> {
        self.chromatograms = self.chromatograms.saturating_add(1);
        Ok(ControlFlow::Continue(()))
    }
}

/// A `Write` sink that refuses to grow past a ceiling.
///
/// The rendered intermediate document is the one allocation whose size follows
/// the record, so it is bounded while it is produced rather than measured
/// afterwards.
struct Bounded {
    data: Vec<u8>,
    remaining: usize,
}

impl Write for Bounded {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self.remaining.checked_sub(buf.len()) {
            None => Err(std::io::Error::other(
                "mzML writing consumer record byte limit exceeded",
            )),
            Some(remaining) => {
                self.data
                    .try_reserve(buf.len())
                    .map_err(|_| std::io::Error::other("mzML record buffer allocation failed"))?;
                self.data.extend_from_slice(buf);
                self.remaining = remaining;
                Ok(buf.len())
            }
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Split a one-record document into its header text and its record element.
///
/// `open` and `close` are the exact list tags the whole-document writer emits
/// for a single record. Both are ASCII and both are located with
/// [`str::find`], so every index used here is a character boundary and no
/// slice can split a multi-byte character in a caller-supplied identifier.
fn split<'a>(document: &'a str, open: &str, close: &str, what: &str) -> Result<(&'a str, &'a str)> {
    let open_at = document.find(open).ok_or_else(|| layout(what))?;
    let prefix = document.get(..open_at).ok_or_else(|| layout(what))?;
    let body_at = open_at
        .checked_add(open.len())
        .ok_or_else(|| layout(what))?;
    let rest = document.get(body_at..).ok_or_else(|| layout(what))?;
    let close_at = rest.find(close).ok_or_else(|| layout(what))?;
    let block = rest.get(..close_at).ok_or_else(|| layout(what))?;
    Ok((prefix, block))
}

/// Replace a lone record's `index="0"` with its position in the stream.
///
/// The source passes the running record number into `writeSpectrum_`; a record
/// rendered on its own is always at index zero, so the attribute is restamped.
/// The marker must occur exactly once - it cannot appear inside an attribute
/// value, because the writer escapes `"` as `&quot;` - and a document where it
/// does not is refused rather than silently mis-numbered.
fn stamp_index(block: &str, index: usize) -> Result<String> {
    if block.matches(INDEX_ZERO).count() != 1 {
        return Err(layout("record index attribute"));
    }
    Ok(block.replacen(
        INDEX_ZERO,
        &format!(" index=\"{index}\" defaultArrayLength=\""),
        1,
    ))
}

/// The rendered text of the header lists a record's references index into.
///
/// A list that is absent contributes nothing, so two headers that both omit it
/// still compare equal. Both markers are ASCII and located with [`str::find`],
/// so every index used here is a character boundary.
fn declaration_text(head: &str) -> String {
    let mut text = String::new();
    for (open, close) in DECLARATION_BLOCKS {
        let Some(open_at) = head.find(open) else {
            continue;
        };
        let Some(rest) = head.get(open_at..) else {
            continue;
        };
        let Some(close_at) = rest.find(close) else {
            continue;
        };
        if let Some(block) = close_at
            .checked_add(close.len())
            .and_then(|end| rest.get(..end))
        {
            text.push_str(block);
        }
    }
    text
}

/// Collect the `id` values of header elements a record may reference.
fn declared_ids(head: &str) -> BTreeSet<String> {
    let mut ids = BTreeSet::new();
    for pattern in HEADER_DECLARATIONS {
        for value in attribute_values(head, pattern) {
            ids.insert(value.to_owned());
        }
    }
    ids
}

/// Refuse a record element referencing an element the header did not declare.
fn check_references(block: &str, declared: &BTreeSet<String>) -> Result<()> {
    for pattern in RECORD_REFERENCES {
        for value in attribute_values(block, pattern) {
            if !declared.contains(value) {
                return Err(Error::Unsupported(format!(
                    "record references undeclared mzML element {value:?}; only the first record \
                     contributes to the header of a streamed file"
                )));
            }
        }
    }
    Ok(())
}

/// Every attribute value introduced by `pattern`, which must end in `="`.
///
/// The scan matches ASCII markers and each value ends at the next `"`; an
/// unterminated value yields nothing rather than a slice past the end.
fn attribute_values<'a>(text: &'a str, pattern: &str) -> Vec<&'a str> {
    let mut found = Vec::new();
    for (at, _) in text.match_indices(pattern) {
        let Some(start) = at.checked_add(pattern.len()) else {
            continue;
        };
        let Some(rest) = text.get(start..) else {
            continue;
        };
        let Some(end) = rest.find('"') else {
            continue;
        };
        if let Some(value) = rest.get(..end) {
            found.push(value);
        }
    }
    found
}

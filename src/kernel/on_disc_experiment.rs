// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Random access to an indexed mzML file without loading it into memory.
//!
//! Ports `KERNEL/OnDiscMSExperiment.h` and its implementation. See
//! `docs/ON_DISC_EXPERIMENT_SUPPORT.md`.
//!
//! [`OnDiscMSExperiment`](crate::kernel::on_disc_experiment::OnDiscMSExperiment)
//! is a thin facade: the index parsing, the per-record byte ranges and the record
//! decoding all live in
//! [`IndexedMzMLHandler`](crate::format::indexed_mzml_handler::IndexedMzMLHandler),
//! exactly as the source class body is almost pure delegation to
//! `Internal::IndexedMzMLHandler`. What the facade adds is the source's own two
//! contributions: an in-memory metadata experiment loaded once with
//! `fillData=false`, and the `PeakFileOptions` filtering applied around each
//! fetch.
//!
//! The metadata copy is what makes a filtered fetch cheap. Retention time, MS
//! level and precursor m/z are known before any peak data is read, so a record
//! the options exclude costs no I/O and is returned as its metadata record with
//! no peaks — which keeps every index valid, unlike in-memory loading where a
//! filtered spectrum disappears from the container.
//!
//! The source documents this class as not thread-safe, because the handler holds
//! one file position, and recommends `#pragma omp parallel for
//! firstprivate(ondisc_map)`. Every fetch here takes `&mut self`, so the
//! compiler enforces exclusive access; a second reader is a second
//! [`open`](crate::kernel::on_disc_experiment::OnDiscMSExperiment::open) or a
//! [`try_clone`](crate::kernel::on_disc_experiment::OnDiscMSExperiment::try_clone),
//! which reopens the file as the source's copy constructor does. Nothing here
//! starts a thread.

use crate::format::indexed_mzml_handler::{IndexedMzMLHandler, RecordReadLimits};
use crate::format::mzml::{self, LoadOptions, ReadOptions};
use crate::format::peak_options::PeakFileOptions;
use crate::kernel::{MSChromatogram, MSExperiment, MSSpectrum, NumericRange};
use crate::metadata::ExperimentalSettings;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// The source `typedef OpenMS::OnDiscMSExperiment OnDiscPeakMap`, which is the
/// name its own class test uses throughout.
pub type OnDiscPeakMap = OnDiscMSExperiment;

/// Explicit ceilings checked before the facade allocates or mutates anything.
///
/// The source has no ceilings of its own: `loadMetaData_` hands the file to
/// `FileHandler` and the per-record reads are the handler's. These fields add
/// the two the facade itself needs — the native-identifier caches it builds from
/// the metadata, and the whole-file materialisation — and carry the handler's
/// and the mzML reader's own limits through unchanged.
#[derive(Clone, Copy, Debug)]
pub struct OnDiscLimits {
    /// Ceilings for one record's byte range, handed to the handler.
    pub record: RecordReadLimits,
    /// XML and binary-array ceilings, used for the metadata load and each record.
    pub read: ReadOptions,
    /// Maximum metadata spectra plus chromatograms accepted from one file.
    pub max_records: usize,
    /// Cumulative bytes of the two native-identifier caches, counting each
    /// identifier once for the key and once for the metadata record it names.
    pub max_native_id_bytes: usize,
    /// Maximum peaks plus chromatogram points
    /// [`load_experiment`](OnDiscMSExperiment::load_experiment) may hold at once.
    pub max_materialized_points: usize,
}

impl Default for OnDiscLimits {
    fn default() -> Self {
        Self {
            record: RecordReadLimits::default(),
            read: ReadOptions::default(),
            max_records: 1_000_000,
            max_native_id_bytes: 64 << 20,
            max_materialized_points: 10_000_000,
        }
    }
}

/// Representation of a mass spectrometry experiment on disc.
///
/// Random access to the spectra and chromatograms of an indexed mzML file
/// without loading the whole file into memory. A default-constructed value holds
/// no file; [`open_file`](Self::open_file) opens one, as the source's
/// documentation for its default constructor directs.
///
/// # Filtering with `PeakFileOptions`
///
/// [`options`](Self::options) filter the data a fetch returns, in two groups:
///
/// - RT range, MS level and precursor m/z range are checked **before** peak data
///   is loaded, so the I/O is skipped entirely;
/// - m/z range and intensity range are applied **after** loading, because they
///   select peaks inside a record.
///
/// Unlike in-memory loading, where a filtered spectrum is removed from the
/// container, this facade preserves every index: a spectrum that fails the RT,
/// MS-level or precursor filter comes back with its metadata and no peaks. That
/// is the source's `@note`, and it is why [`spectrum`](Self::spectrum) returns a
/// spectrum rather than an `Option`.
///
/// The by-native-identifier fetches are the exception: the source applies no
/// filtering at all on those paths, and neither does this port. See
/// [`spectrum_by_native_id`](Self::spectrum_by_native_id).
///
/// # Examples
///
/// ```
/// use openms::kernel::on_disc_experiment::OnDiscMSExperiment;
///
/// let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
///     .join("tests/data/indexed_mzml/IndexedmzMLFile_1.mzML");
/// let mut experiment = OnDiscMSExperiment::open(path)?;
/// assert_eq!(experiment.len(), 2);
/// assert_eq!(experiment.chromatogram_count(), 1);
/// assert_eq!(
///     experiment.experimental_settings().unwrap().instrument.name,
///     "LTQ FT"
/// );
///
/// // Only this record's byte range is read.
/// let spectrum = experiment.spectrum(0)?;
/// assert_eq!(spectrum.peaks.len(), 19914);
///
/// // An MS level the options exclude costs no peak I/O, and keeps its metadata.
/// experiment.options_mut().set_ms_levels(&[2])?;
/// let filtered = experiment.spectrum(0)?;
/// assert!(filtered.peaks.is_empty());
/// assert_eq!(filtered.ms_level, 1);
///
/// // The source's own example walks every index and skips the empty records,
/// // which is how a filtered spectrum is recognised without losing its slot.
/// let mut kept = 0;
/// for index in 0..experiment.len() {
///     if experiment.spectrum(index)?.peaks.is_empty() {
///         continue; // filtered out; no peak data was read
///     }
///     kept += 1;
/// }
/// assert_eq!(kept, 0); // this file has no MS2 spectrum
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Debug, Default)]
pub struct OnDiscMSExperiment {
    path: PathBuf,
    handler: Option<IndexedMzMLHandler>,
    meta: Option<MSExperiment>,
    spectra_native_ids: BTreeMap<String, usize>,
    chromatograms_native_ids: BTreeMap<String, usize>,
    options: PeakFileOptions,
    skip_xml_checks: bool,
    limits: OnDiscLimits,
}

impl OnDiscMSExperiment {
    /// An experiment with no file open, as the source default constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// An experiment with no file open and explicit ceilings.
    pub fn with_limits(limits: OnDiscLimits) -> Self {
        Self {
            limits,
            ..Self::default()
        }
    }

    /// Open `path`, requiring a usable index.
    ///
    /// Native convenience over [`open_file`](Self::open_file): the source's
    /// `false` return — the file exists but carries no usable mzML index —
    /// becomes [`Error::Parse`] here, so a value returned by this function can
    /// always be fetched from.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when the file is not an indexed mzML, plus anything
    /// [`open_file`](Self::open_file) reports.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let mut experiment = Self::new();
        let path = path.as_ref();
        if !experiment.open_file(path, false)? {
            return Err(Error::Parse {
                line: 0,
                message: format!("{} is not an indexed mzML file", path.display()),
            });
        }
        Ok(experiment)
    }

    /// Open a specific file on disc.
    ///
    /// Parses the index and, unless `skip_metadata` is set, reads the meta
    /// information into memory. Returns whether the index parsed: `false` means
    /// the file most likely was not an indexed mzML, and every count is then
    /// zero while the metadata — if it was requested and the file is readable
    /// mzML — is still available.
    ///
    /// `skip_metadata` is the source's `skipMetaData` parameter and defaults to
    /// `false` there. Skipping it means [`experimental_settings`](Self::experimental_settings)
    /// and [`metadata`](Self::metadata) stay `None`, [`is_sorted_by_rt`](Self::is_sorted_by_rt)
    /// reports `false`, and the RT, MS-level and precursor filters are not
    /// applied at all, because the values they test are only known from the
    /// metadata. The source has the same three consequences and documents none
    /// of them.
    ///
    /// The source skips the metadata load for an empty filename as well; so does
    /// this, and an empty path additionally never reaches the index parser.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be read, and any error the mzML reader
    /// raises on the metadata load — the source lets the corresponding
    /// `FileHandler::loadExperiment` exception escape too. Resource ceilings
    /// report [`Error::InvalidValue`]: more metadata records than
    /// [`OnDiscLimits::max_records`], or native identifiers exceeding
    /// [`OnDiscLimits::max_native_id_bytes`].
    ///
    /// Everything is built into temporaries and committed at the end, so a
    /// failed open leaves the previous file, metadata, identifier caches and
    /// options in place. The source assigns `filename_` first and replaces
    /// `meta_ms_experiment_` inside `loadMetaData_`, so a throwing metadata load
    /// leaves it naming the new file with the old — or a half-built — metadata.
    pub fn open_file(&mut self, path: impl AsRef<Path>, skip_metadata: bool) -> Result<bool> {
        let path = path.as_ref();
        let empty = path.as_os_str().is_empty();

        // Source `openFile` always calls the handler, which records failure in
        // `parsing_success_` rather than raising it.
        let handler = if empty {
            None
        } else {
            IndexedMzMLHandler::open_with_limits(path, self.limits.record, self.limits.read).ok()
        };

        let meta = if empty || skip_metadata {
            None
        } else {
            Some(self.load_meta_data(path)?)
        };
        let (spectra_native_ids, chromatograms_native_ids) = match &meta {
            Some(meta) => self.identifier_caches(meta)?,
            None => (BTreeMap::new(), BTreeMap::new()),
        };

        let parsed = handler.is_some();
        let mut handler = handler;
        if let Some(handler) = handler.as_mut() {
            // Source `setSkipXMLChecks` is remembered by the handler across
            // `openFile`, because it is a handler member that `openFile` never
            // resets. A fresh handler here is given the remembered flag.
            handler.options_mut().skip_xml_checks = self.skip_xml_checks;
        }

        self.path = path.to_path_buf();
        self.handler = handler;
        self.meta = meta;
        self.spectra_native_ids = spectra_native_ids;
        self.chromatograms_native_ids = chromatograms_native_ids;
        Ok(parsed)
    }

    /// The path the experiment was opened on, as the source's protected
    /// `filename_`. Empty until [`open_file`](Self::open_file) succeeds.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// The ceilings this experiment enforces.
    pub fn limits(&self) -> OnDiscLimits {
        self.limits
    }

    /// Whether the index parsed, i.e. the last [`open_file`](Self::open_file)
    /// return value.
    ///
    /// Native: the source keeps this only in the handler's `parsing_success_`,
    /// which the facade does not expose. Fetching from an experiment for which
    /// this is `false` is an error rather than the source's exception from deep
    /// inside the handler.
    pub fn is_indexed(&self) -> bool {
        self.handler.is_some()
    }

    /// Reopen the same file as an independent reader, as the source copy
    /// constructor.
    ///
    /// The source copies `filename_`, the handler (whose own copy constructor
    /// deliberately reopens the file, because "this is critical for parallel
    /// access to the same file"), the metadata pointer and the options. This
    /// does the same, except that the metadata experiment is cloned rather than
    /// shared and the two identifier caches are rebuilt: the source's copy
    /// constructor omits them from its initialiser list, and so does the
    /// handler's, which is why native-identifier lookups on a source copy always
    /// throw.
    ///
    /// `Clone` is not implemented because reopening the file can fail.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file can no longer be opened, and
    /// [`Error::Parse`] when its index no longer decodes. An experiment whose
    /// index never parsed clones without a handler and without error.
    pub fn try_clone(&self) -> Result<Self> {
        let handler = match &self.handler {
            Some(_) => {
                let mut handler = IndexedMzMLHandler::open_with_limits(
                    &self.path,
                    self.limits.record,
                    self.limits.read,
                )?;
                handler.options_mut().skip_xml_checks = self.skip_xml_checks;
                Some(handler)
            }
            None => None,
        };
        Ok(Self {
            path: self.path.clone(),
            handler,
            meta: self.meta.clone(),
            spectra_native_ids: self.spectra_native_ids.clone(),
            chromatograms_native_ids: self.chromatograms_native_ids.clone(),
            options: self.options.clone(),
            skip_xml_checks: self.skip_xml_checks,
            limits: self.limits,
        })
    }

    /// Whether all spectra are sorted by ascending retention time.
    ///
    /// Only the metadata is consulted, so this says nothing about the m/z order
    /// inside each spectrum — the source notes that checking that would mean
    /// loading every spectrum. Without metadata the answer is `false`, as in the
    /// source, which cannot distinguish "unsorted" from "unknown" either.
    pub fn is_sorted_by_rt(&self) -> bool {
        match &self.meta {
            Some(meta) => meta.is_sorted(false),
            None => false,
        }
    }

    /// The number of spectra available, as source `size()` and `getNrSpectra()`.
    pub fn len(&self) -> usize {
        self.spectrum_count()
    }

    /// Whether no spectra are available, as source `empty()`.
    ///
    /// A file with chromatograms but no spectra is empty here, exactly as in the
    /// source, whose `empty()` also asks only about spectra.
    pub fn is_empty(&self) -> bool {
        self.spectrum_count() == 0
    }

    /// The total number of spectra available, as source `getNrSpectra()`.
    ///
    /// Read from the index, so zero whenever the index did not parse, even when
    /// the metadata describes spectra.
    pub fn spectrum_count(&self) -> usize {
        self.handler
            .as_ref()
            .map_or(0, IndexedMzMLHandler::spectrum_count)
    }

    /// The total number of chromatograms available, as source
    /// `getNrChromatograms()`.
    pub fn chromatogram_count(&self) -> usize {
        self.handler
            .as_ref()
            .map_or(0, IndexedMzMLHandler::chromatogram_count)
    }

    /// The meta information of this experiment, as source
    /// `getExperimentalSettings()`.
    ///
    /// `None` when the file was opened with `skip_metadata`, matching the null
    /// `shared_ptr` the source returns; its own class test asserts that null.
    pub fn experimental_settings(&self) -> Option<&ExperimentalSettings> {
        self.meta.as_ref().map(|meta| &meta.settings)
    }

    /// The metadata experiment, as source `getMetaData()`.
    ///
    /// Its spectra and chromatograms carry every field except peak data: the
    /// metadata load sets `fillData=false`. Indices correspond one-to-one with
    /// the index, which is why the source loads the metadata unfiltered and
    /// filters only at retrieval time.
    ///
    /// The source hands out a non-const `std::shared_ptr<PeakMap>`, so a caller
    /// can mutate the metadata the facade is using — and, because the source
    /// builds its native-identifier maps lazily and only when they are empty,
    /// such a mutation silently desynchronises them. This returns a shared
    /// reference instead; clone it if an owned copy is wanted.
    pub fn metadata(&self) -> Option<&MSExperiment> {
        self.meta.as_ref()
    }

    /// Mutable access to the options for loading/storing, as source
    /// `getOptions()`.
    ///
    /// The source calls them loading/storing options because the type is shared
    /// with the writers; nothing in this facade stores.
    ///
    /// Note that `skip_xml_checks` on these options is inert, as in the source:
    /// the facade's options are never handed to the decoder. Use
    /// [`set_skip_xml_checks`](Self::set_skip_xml_checks).
    pub fn options_mut(&mut self) -> &mut PeakFileOptions {
        &mut self.options
    }

    /// Non-mutable access to the options for loading/storing, as source
    /// `getOptions() const`.
    pub fn options(&self) -> &PeakFileOptions {
        &self.options
    }

    /// Replace the options for loading/storing, as source `setOptions()`.
    ///
    /// Only the RT, MS-level, precursor-m/z, m/z and intensity selections are
    /// consulted by the fetches. `fill_data`, `metadata_only`,
    /// `skip_chromatograms` and the writing options are carried but unused, as
    /// in the source, and `skip_xml_checks` here does not reach the decoder —
    /// see [`options_mut`](Self::options_mut).
    pub fn set_options(&mut self, options: PeakFileOptions) {
        self.options = options;
    }

    /// Whether XML checks are skipped, the value of
    /// [`set_skip_xml_checks`](Self::set_skip_xml_checks). Native: the source
    /// has no getter.
    pub fn skip_xml_checks(&self) -> bool {
        self.skip_xml_checks
    }

    /// Set whether to skip some XML checks and be fast instead, as source
    /// `setSkipXMLChecks()`.
    ///
    /// Forwarded to the handler, where it suppresses the four-character Base64
    /// whitespace strip (`src/format/mzml.rs:2386`). That is the whole effect in
    /// the source too: `skip_xml_checks_` is passed to
    /// `MzMLHandlerHelper::decodeBase64Arrays` and never to XML syntax checking.
    /// The flag is remembered across [`open_file`](Self::open_file), because in
    /// the source it is a handler member that `openFile` does not reset.
    pub fn set_skip_xml_checks(&mut self, skip: bool) {
        self.skip_xml_checks = skip;
        if let Some(handler) = self.handler.as_mut() {
            handler.options_mut().skip_xml_checks = skip;
        }
    }

    /// A single spectrum by index, as source `getSpectrum()` and `operator[]`.
    ///
    /// With metadata available, the RT range, MS level and precursor m/z range
    /// of [`options`](Self::options) are tested first and an excluded record is
    /// returned as its metadata spectrum with no peaks, performing no I/O.
    /// Otherwise the record is read and the m/z and intensity ranges then select
    /// peaks within it. Without metadata only those two ranges apply.
    ///
    /// A record's own metadata is preserved: the source merges the decoded peak
    /// arrays into the metadata spectrum, and so does this.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `index` is not below
    /// [`spectrum_count`](Self::spectrum_count), when the metadata holds no such
    /// spectrum, or when the index did not parse and peaks would be needed;
    /// plus anything the handler reports for this record. The source indexes
    /// both the metadata vector and the offset vector without a bound check on
    /// the former.
    pub fn spectrum(&mut self, index: usize) -> Result<MSSpectrum> {
        let mut spectrum = match self.meta.as_ref() {
            Some(meta) => {
                let base = meta.spectra.get(index).ok_or_else(|| {
                    Error::InvalidValue(format!(
                        "spectrum index {index} is not below the {} metadata spectra",
                        meta.spectra.len()
                    ))
                })?;
                // Source order: each check returns the metadata spectrum before
                // any peak data is read.
                if self.options.has_rt_range() && !encloses(self.options.rt_range(), base.rt) {
                    return Ok(base.clone());
                }
                if self.options.has_ms_levels() && !self.contains_ms_level(base.ms_level) {
                    return Ok(base.clone());
                }
                if self.options.has_precursor_mz_range() && !base.precursors.is_empty() {
                    let precursor_mz = base.precursors[0].mz;
                    if !encloses(self.options.precursor_mz_range(), precursor_mz) {
                        return Ok(base.clone());
                    }
                }
                let mut base = base.clone();
                let fetched = self.fetch_spectrum(index)?;
                overlay_spectrum(&mut base, fetched);
                base
            }
            None => self.fetch_spectrum(index)?,
        };

        if self.options.has_mz_range() || self.options.has_intensity_range() {
            let (mz, intensity) = (self.options.mz_range(), self.options.intensity_range());
            let (has_mz, has_intensity) = (
                self.options.has_mz_range(),
                self.options.has_intensity_range(),
            );
            spectrum.retain_peaks(|peak| {
                (!has_mz || encloses(mz, peak.mz))
                    && (!has_intensity || encloses(intensity, f64::from(peak.intensity)))
            })?;
        }
        Ok(spectrum)
    }

    /// A single chromatogram by index, as source `getChromatogram()`.
    ///
    /// The RT range and intensity range of [`options`](Self::options) select
    /// points within the chromatogram; unlike for spectra there is no
    /// before-I/O check, so the record is always read. That asymmetry is the
    /// source's: RT decides a whole spectrum but only individual chromatogram
    /// points.
    ///
    /// # Errors
    ///
    /// As [`spectrum`](Self::spectrum), for the chromatogram index.
    pub fn chromatogram(&mut self, index: usize) -> Result<MSChromatogram> {
        let mut chromatogram = match self.meta.as_ref() {
            Some(meta) => {
                let mut base = meta
                    .chromatograms
                    .get(index)
                    .ok_or_else(|| {
                        Error::InvalidValue(format!(
                            "chromatogram index {index} is not below the {} metadata chromatograms",
                            meta.chromatograms.len()
                        ))
                    })?
                    .clone();
                let fetched = self.fetch_chromatogram(index)?;
                overlay_chromatogram(&mut base, fetched);
                base
            }
            None => self.fetch_chromatogram(index)?,
        };

        if self.options.has_rt_range() || self.options.has_intensity_range() {
            let (rt, intensity) = (self.options.rt_range(), self.options.intensity_range());
            let (has_rt, has_intensity) = (
                self.options.has_rt_range(),
                self.options.has_intensity_range(),
            );
            chromatogram.retain_peaks(|peak| {
                (!has_rt || encloses(rt, peak.rt))
                    && (!has_intensity || encloses(intensity, f64::from(peak.intensity)))
            })?;
        }
        Ok(chromatogram)
    }

    /// A single spectrum by its native identifier, as source
    /// `getSpectrumByNativeId()`.
    ///
    /// No option is consulted on this path: the source resolves the identifier,
    /// merges the peaks into the metadata spectrum and returns it, never
    /// touching `options_`. A caller who has configured an m/z range therefore
    /// gets unfiltered peaks here and filtered peaks from
    /// [`spectrum`](Self::spectrum). This port preserves that; the difference is
    /// documented in `docs/ON_DISC_EXPERIMENT_SUPPORT.md` rather than silently
    /// repaired.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no spectrum carries the identifier, matching
    /// the source's `Exception::IllegalArgument`, or when the index did not
    /// parse; plus anything the handler reports for this record.
    pub fn spectrum_by_native_id(&mut self, native_id: &str) -> Result<MSSpectrum> {
        let index = match self.meta.as_ref() {
            Some(_) => Some(native_index(
                &self.spectra_native_ids,
                native_id,
                "spectrum",
            )?),
            None => None,
        };
        let fetched = self.handler_mut()?.spectrum_by_native_id(native_id)?;
        let mut fetched = fetched.ok_or_else(|| excluded("spectrum"))?;
        if let Some(index) = index {
            let meta = self.meta.as_ref().expect("metadata present");
            let mut base = meta.spectra[index].clone();
            overlay_spectrum(&mut base, fetched);
            fetched = base;
        }
        Ok(fetched)
    }

    /// A single chromatogram by its native identifier, as source
    /// `getChromatogramByNativeId()`.
    ///
    /// Unfiltered, as [`spectrum_by_native_id`](Self::spectrum_by_native_id).
    ///
    /// # Errors
    ///
    /// As [`spectrum_by_native_id`](Self::spectrum_by_native_id).
    pub fn chromatogram_by_native_id(&mut self, native_id: &str) -> Result<MSChromatogram> {
        let index = match self.meta.as_ref() {
            Some(_) => Some(native_index(
                &self.chromatograms_native_ids,
                native_id,
                "chromatogram",
            )?),
            None => None,
        };
        let fetched = self.handler_mut()?.chromatogram_by_native_id(native_id)?;
        let mut fetched = fetched.ok_or_else(|| excluded("chromatogram"))?;
        if let Some(index) = index {
            let meta = self.meta.as_ref().expect("metadata present");
            let mut base = meta.chromatograms[index].clone();
            overlay_chromatogram(&mut base, fetched);
            fetched = base;
        }
        Ok(fetched)
    }

    /// Materialise the whole file into an in-memory experiment.
    ///
    /// Walks every index entry through [`spectrum`](Self::spectrum) and
    /// [`chromatogram`](Self::chromatogram), so [`options`](Self::options) apply
    /// as they do for a single fetch and an excluded record arrives as its
    /// metadata with no peaks, keeping the indices aligned with the file. The
    /// experimental settings are copied from the metadata, or left default when
    /// the file was opened with `skip_metadata`.
    ///
    /// There is no such member in the source; it is the loop of
    /// `IndexedMzMLFileLoader::store`
    /// (`src/openms/source/FORMAT/IndexedMzMLFileLoader.cpp:50`), which feeds
    /// the same two walks to a writing consumer, collected into an
    /// [`MSExperiment`] instead. That loop also dereferences
    /// `getExperimentalSettings()` without a null check, which this cannot do.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the index lists more than
    /// [`OnDiscLimits::max_records`] records, or when the peaks and points
    /// collected would exceed [`OnDiscLimits::max_materialized_points`]; plus
    /// anything a single fetch reports. The experiment is built in a temporary,
    /// so a failure leaves this object and any destination untouched.
    pub fn load_experiment(&mut self) -> Result<MSExperiment> {
        let spectra = self.spectrum_count();
        let chromatograms = self.chromatogram_count();
        let records = spectra
            .checked_add(chromatograms)
            .ok_or_else(|| limit("record count"))?;
        if records > self.limits.max_records {
            return Err(limit("record count"));
        }

        let mut experiment = MSExperiment::new();
        experiment
            .spectra
            .try_reserve_exact(spectra)
            .map_err(|_| limit("spectrum vector"))?;
        experiment
            .chromatograms
            .try_reserve_exact(chromatograms)
            .map_err(|_| limit("chromatogram vector"))?;

        let mut points = 0usize;
        for index in 0..spectra {
            let spectrum = self.spectrum(index)?;
            points = points
                .checked_add(spectrum.peaks.len())
                .filter(|&total| total <= self.limits.max_materialized_points)
                .ok_or_else(|| limit("materialised point count"))?;
            experiment.spectra.push(spectrum);
        }
        for index in 0..chromatograms {
            let chromatogram = self.chromatogram(index)?;
            points = points
                .checked_add(chromatogram.peaks.len())
                .filter(|&total| total <= self.limits.max_materialized_points)
                .ok_or_else(|| limit("materialised point count"))?;
            experiment.chromatograms.push(chromatogram);
        }
        if let Some(meta) = self.meta.as_ref() {
            experiment.settings = meta.settings.clone();
        }
        Ok(experiment)
    }

    /// Source `loadMetaData_`: always the complete metadata, unfiltered, with
    /// `fillData=false`, so the indices correspond to the index sections.
    fn load_meta_data(&self, path: &Path) -> Result<MSExperiment> {
        let mut load = LoadOptions::default();
        load.scientific.fill_data = false;
        mzml::load_with_options(path, &load, &self.limits.read)
    }

    /// Source `getMetaSpectrumById_`/`getMetaChromatogramById_` build these
    /// lazily with `unordered_map::emplace`, so the first entry wins for a
    /// duplicate identifier. Built eagerly here and replaced together with the
    /// metadata, because the source never clears them on a second `openFile`.
    fn identifier_caches(
        &self,
        meta: &MSExperiment,
    ) -> Result<(BTreeMap<String, usize>, BTreeMap<String, usize>)> {
        let records = meta
            .spectra
            .len()
            .checked_add(meta.chromatograms.len())
            .ok_or_else(|| limit("metadata record count"))?;
        if records > self.limits.max_records {
            return Err(limit("metadata record count"));
        }
        let mut budget = self.limits.max_native_id_bytes;
        let lengths = meta
            .spectra
            .iter()
            .map(|spectrum| spectrum.native_id.len())
            .chain(
                meta.chromatograms
                    .iter()
                    .map(|chromatogram| chromatogram.native_id.len()),
            );
        for length in lengths {
            budget = budget
                .checked_sub(length.saturating_mul(2))
                .ok_or_else(|| limit("native identifier bytes"))?;
        }

        let mut spectra = BTreeMap::new();
        for (index, spectrum) in meta.spectra.iter().enumerate() {
            spectra
                .entry(spectrum.native_id.clone())
                .or_insert_with(|| index);
        }
        let mut chromatograms = BTreeMap::new();
        for (index, chromatogram) in meta.chromatograms.iter().enumerate() {
            chromatograms
                .entry(chromatogram.native_id.clone())
                .or_insert_with(|| index);
        }
        Ok((spectra, chromatograms))
    }

    fn handler_mut(&mut self) -> Result<&mut IndexedMzMLHandler> {
        self.handler.as_mut().ok_or_else(|| {
            Error::InvalidValue(
                "no mzML index was parsed for this experiment; peak data cannot be read".into(),
            )
        })
    }

    fn fetch_spectrum(&mut self, index: usize) -> Result<MSSpectrum> {
        self.handler_mut()?
            .spectrum(index)?
            .ok_or_else(|| excluded("spectrum"))
    }

    fn fetch_chromatogram(&mut self, index: usize) -> Result<MSChromatogram> {
        self.handler_mut()?
            .chromatogram(index)?
            .ok_or_else(|| excluded("chromatogram"))
    }

    /// Source `PeakFileOptions::containsMSLevel(Int)` against a `UInt` level, an
    /// implicit narrowing conversion. A level no `i32` can hold matches nothing.
    fn contains_ms_level(&self, level: u32) -> bool {
        i32::try_from(level).is_ok_and(|level| self.options.contains_ms_level(level))
    }
}

/// Equality of the underlying file and the parsed meta information only.
///
/// The source's `operator==` compares `filename_` and the metadata, and
/// deliberately does not compare the index: two experiments on the same file
/// have the same index by construction, and the file reader may be at different
/// positions. When either side has no metadata the source falls back to
/// comparing the two `shared_ptr`s, so an experiment opened with
/// `skip_metadata` equals only another one that also has no metadata and the
/// same path — never one that loaded its metadata, even from the same file.
/// That is reproduced here, including the reflexive `true`.
impl PartialEq for OnDiscMSExperiment {
    fn eq(&self, other: &Self) -> bool {
        match (&self.meta, &other.meta) {
            (Some(mine), Some(theirs)) => self.path == other.path && mine == theirs,
            _ => self.path == other.path && self.meta.is_none() && other.meta.is_none(),
        }
    }
}

/// Source `DRange<1>::encloses` (`DATASTRUCTURES/DRange.h:152`): the interval is
/// half-open, `[min, max)`, so a value equal to the maximum is outside it. This
/// is the same rule `src/format/mzml_load.rs` applies to a whole-file load.
fn encloses(range: NumericRange, value: f64) -> bool {
    !(value < range.min || value >= range.max)
}

/// Source `MzMLSpectrumDecoder::domParseSpectrum(const std::string&, MSSpectrum&)`
/// appends the decoded peaks and auxiliary arrays to the caller's spectrum and
/// overwrites its native identifier, leaving every other field alone. The
/// metadata load fills no arrays, so appending and replacing coincide;
/// replacing is used so that a metadata record which did carry arrays cannot
/// produce duplicates.
fn overlay_spectrum(base: &mut MSSpectrum, fetched: MSSpectrum) {
    base.peaks = fetched.peaks;
    base.float_data_arrays = fetched.float_data_arrays;
    base.integer_data_arrays = fetched.integer_data_arrays;
    base.string_data_arrays = fetched.string_data_arrays;
    base.native_id = fetched.native_id;
}

/// As [`overlay_spectrum`], for source `domParseChromatogram`.
fn overlay_chromatogram(base: &mut MSChromatogram, fetched: MSChromatogram) {
    base.peaks = fetched.peaks;
    base.float_data_arrays = fetched.float_data_arrays;
    base.integer_data_arrays = fetched.integer_data_arrays;
    base.string_data_arrays = fetched.string_data_arrays;
    base.native_id = fetched.native_id;
}

/// Source `getMetaSpectrumById_`/`getMetaChromatogramById_` raise
/// `Exception::IllegalArgument` with exactly this wording for an identifier the
/// metadata does not carry.
fn native_index(cache: &BTreeMap<String, usize>, native_id: &str, what: &str) -> Result<usize> {
    cache
        .get(native_id)
        .copied()
        .ok_or_else(|| Error::InvalidValue(format!("could not find {what} with id '{native_id}'.")))
}

fn limit(what: &str) -> Error {
    Error::InvalidValue(format!(
        "on-disc experiment {what} exceeds its configured limit"
    ))
}

/// The handler reports a record its own options exclude as `None`. The facade
/// keeps those options unfiltered, so this cannot be reached through the public
/// API; it stays an error rather than a silent empty record.
fn excluded(what: &str) -> Error {
    Error::InvalidValue(format!(
        "the indexed mzML reader excluded this {what} before the on-disc filters ran"
    ))
}

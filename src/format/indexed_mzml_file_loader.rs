// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Loads an indexed mzML file into an on-disc experiment, and stores on-disc or
//! in-memory experiments as indexed mzML.
//!
//! Ports `FORMAT/IndexedMzMLFileLoader.h` and its `.cpp`. See
//! `docs/INDEXED_MZML_FILE_LOADER_SUPPORT.md` for the API mapping and the
//! evidence.
//!
//! The `OnDiscPeakMap` typedef the header declares is
//! [`OnDiscPeakMap`](crate::kernel::on_disc_experiment::OnDiscPeakMap).

use crate::format::ms_data_writing_consumer::{PlainMSDataWritingConsumer, ReferencePolicy};
use crate::format::mzml;
use crate::format::path_io;
use crate::format::peak_options::PeakFileOptions;
use crate::kernel::MSExperiment;
use crate::kernel::on_disc_experiment::OnDiscPeakMap;
use crate::{Error, Result};
use std::path::Path;

/// A class to load an indexed mzML file.
///
/// Provides the same interface as the other file classes (mzML, mzXML and so
/// on) to load and store a file. Reading a file from disk loads it into an
/// [`OnDiscPeakMap`], while both an [`OnDiscPeakMap`] and an in-memory
/// [`MSExperiment`] can be written to disk. Both stores write `indexedmzML`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IndexedMzMLFileLoader {
    /// Options for storing; source `options_`.
    options: PeakFileOptions,
}

impl IndexedMzMLFileLoader {
    /// A loader with default options; source `IndexedMzMLFileLoader()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Mutable access to the options for loading/storing; source
    /// `PeakFileOptions& getOptions()`.
    pub fn options_mut(&mut self) -> &mut PeakFileOptions {
        &mut self.options
    }

    /// Non-mutable access to the options for loading/storing; source
    /// `const PeakFileOptions& getOptions() const`.
    pub fn options(&self) -> &PeakFileOptions {
        &self.options
    }

    /// Sets the options for loading/storing; source `setOptions`.
    pub fn set_options(&mut self, options: PeakFileOptions) {
        self.options = options;
    }

    /// Loads `path` into `experiment`; source `bool load(filename, OnDiscPeakMap&)`.
    ///
    /// Tries to parse the file; success is the return value and must be
    /// checked. `Ok(true)` means an indexed mzML was opened. `Ok(false)` means
    /// the file exists but most likely is not an mzML or is not indexed; the
    /// experiment then has no records to fetch. The call is exactly the
    /// source's `exp.openFile(filename)`, which loads the metadata as well.
    ///
    /// The loader's own [`options`](Self::options) are **not** passed to the
    /// experiment, as in the source, whose `load` never reads `options_`; an
    /// [`OnDiscPeakMap`] carries its own options
    /// ([`OnDiscMSExperiment::options_mut`](crate::kernel::on_disc_experiment::OnDiscMSExperiment::options_mut)).
    ///
    /// # Errors
    ///
    /// As [`OnDiscMSExperiment::open_file`](crate::kernel::on_disc_experiment::OnDiscMSExperiment::open_file):
    /// [`Error::Io`] when the file does not exist or cannot be read (the
    /// source's `Exception::FileNotFound` and `Exception::FileNotReadable`),
    /// and any error of the metadata load. The source's `@note` adds that a
    /// corrupt `indexListOffset` throws `Exception::ConversionError`; here the
    /// on-disc experiment's index parser treats every unusable index alike, as
    /// "not indexed", so that case is `Ok(false)`.
    pub fn load(&self, path: impl AsRef<Path>, experiment: &mut OnDiscPeakMap) -> Result<bool> {
        experiment.open_file(path, false)
    }

    /// Stores an on-disc experiment at `path`; source `store(filename, OnDiscPeakMap&)`.
    ///
    /// Streams the experiment record by record, as the source does: every
    /// spectrum is fetched from disk and handed to a
    /// [`PlainMSDataWritingConsumer`], then every chromatogram, so only one
    /// record is held in memory at a time. The expected sizes are the
    /// experiment's record counts, the settings are its metadata, and the
    /// consumer writes with the loader's options (precision, zlib, Numpress,
    /// TPP compatibility), as the source's `consumer.setOptions(options_)`
    /// does. The output is byte-identical to [`store`](Self::store) of the
    /// same records with the same options.
    ///
    /// Like the source, this first sets the loader's own
    /// [`PeakFileOptions::write_index`] to `true` ("ensure that we write the
    /// index"), and the change persists in [`options`](Self::options) after
    /// the call. The file is written plain even when `path` ends in `.gz` or
    /// `.bz2`, as the source's consumer writes through an `std::ofstream`; it
    /// is published atomically, so a failed store leaves no partial file.
    ///
    /// Records whose data-processing history or source file differs from the
    /// first record's are refused ([`ReferencePolicy::Checked`]); the source
    /// writes them with a dangling reference. Use
    /// [`store_on_disc_with_policy`](Self::store_on_disc_with_policy) with
    /// [`ReferencePolicy::SourceDangling`] for the source's output.
    ///
    /// # Errors
    ///
    /// * [`Error::MissingInformation`] when the experiment holds no metadata,
    ///   because it was opened with `skip_metadata` or never opened. The
    ///   source dereferences the null settings pointer
    ///   (`IndexedMzMLFileLoader.cpp:47`).
    /// * Every error a record fetch or the writing consumer reports, including
    ///   [`Error::Io`].
    pub fn store_on_disc(
        &mut self,
        path: impl AsRef<Path>,
        experiment: &mut OnDiscPeakMap,
    ) -> Result<()> {
        self.store_on_disc_with_policy(path, experiment, ReferencePolicy::Checked)
    }

    /// [`store_on_disc`](Self::store_on_disc) with an explicit
    /// [`ReferencePolicy`] for records whose header references differ from the
    /// first record's.
    ///
    /// No source counterpart; the source always writes the dangling reference,
    /// which [`ReferencePolicy::SourceDangling`] reproduces.
    ///
    /// # Errors
    ///
    /// As [`store_on_disc`](Self::store_on_disc).
    pub fn store_on_disc_with_policy(
        &mut self,
        path: impl AsRef<Path>,
        experiment: &mut OnDiscPeakMap,
        references: ReferencePolicy,
    ) -> Result<()> {
        // Source order: the consumer is created and given the sizes and the
        // settings, then `options_.setWriteIndex(true)` and `setOptions`, all
        // before the first record; nothing is written before that record.
        self.options.write_index = true;
        let options = self.options.clone();
        let settings = experiment
            .experimental_settings()
            .ok_or_else(|| {
                Error::MissingInformation(
                    "the on-disc experiment holds no metadata to store; it was opened \
                     without metadata or not opened"
                        .into(),
                )
            })?
            .clone();
        let spectra = experiment.spectrum_count();
        let chromatograms = experiment.chromatogram_count();
        path_io::write_plain(path.as_ref(), |writer| {
            let mut consumer =
                PlainMSDataWritingConsumer::plain(writer).with_reference_policy(references);
            consumer.set_expected_size(spectra, chromatograms)?;
            consumer.set_experimental_settings(&settings)?;
            consumer.set_peak_file_options(options)?;
            for index in 0..spectra {
                let mut spectrum = experiment.spectrum(index)?;
                consumer.consume_spectrum(&mut spectrum)?;
            }
            for index in 0..chromatograms {
                let mut chromatogram = experiment.chromatogram(index)?;
                consumer.consume_chromatogram(&mut chromatogram)?;
            }
            consumer.finish()?;
            Ok(())
        })
    }

    /// Stores an in-memory experiment at `path`; source `store(filename, PeakMap&)`.
    ///
    /// Sets the loader's own [`PeakFileOptions::write_index`] to `true`, which
    /// persists after the call as in the source, and writes the experiment
    /// with the loader's options through the mzML writer, as the source's
    /// `MzMLFile::store` with those options does. A `.gz` or `.bz2` suffix
    /// compresses the output, as the source's `XMLFile::save_` does.
    ///
    /// # Errors
    ///
    /// Every error of
    /// [`mzml::store_with_peak_options`](crate::format::mzml::store_with_peak_options).
    /// In particular an experiment with neither spectra nor chromatograms is
    /// refused with [`Error::InvalidValue`], because an indexed mzML needs a
    /// record to index; the source writes an index with a dummy `-1` offset
    /// (CPP-050).
    pub fn store(&mut self, path: impl AsRef<Path>, experiment: &MSExperiment) -> Result<()> {
        self.options.write_index = true;
        mzml::store_with_peak_options(path, experiment, &self.options)?;
        Ok(())
    }
}

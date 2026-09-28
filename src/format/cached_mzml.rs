// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Random-access reader and writer for the OpenMS on-disk cache format, ported
//! from `FORMAT/CachedMzML.h` (class `CachedmzML`).
//!
//! Spectra and chromatograms are split across two files: an `.mzML` file
//! holding the experiment-level metadata and a side-car `.mzML.cached` binary
//! file holding the per-spectrum and per-chromatogram payload, in the layout
//! of [`crate::format::cached_mzml_handler`]. After [`CachedmzML::load`](crate::format::cached_mzml::CachedmzML::load) (or
//! [`CachedmzML::open`](crate::format::cached_mzml::CachedmzML::open)) the side-car's byte offsets are kept in an in-memory
//! index, so single spectra or chromatograms can be fetched in any order
//! without scanning the whole file. See `docs/CACHED_MZML_SUPPORT.md`.

use crate::format::cached_mzml_handler::CachedMzMLHandler;
use crate::kernel::{MSChromatogram, MSExperiment, MSSpectrum};
use crate::{Error, Result};
use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};

/// Random-access reader and writer for a cached mzML pair (source
/// `CachedmzML`).
#[derive(Debug, Default)]
pub struct CachedmzML {
    /// The metadata loaded from the `.mzML` file.
    meta_ms_experiment: MSExperiment,
    /// The open side-car stream; `None` before [`Self::load`].
    ifs: Option<BufReader<File>>,
    /// Name of the mzML file.
    filename: PathBuf,
    /// Name of the cached mzML file.
    filename_cached: PathBuf,
    spectra_index: Vec<u64>,
    chrom_index: Vec<u64>,
}

/// The side-car name: the mzML path with `.cached` appended.
fn cached_name(filename: &Path) -> PathBuf {
    let mut name = filename.as_os_str().to_owned();
    name.push(".cached");
    PathBuf::from(name)
}

impl Clone for CachedmzML {
    /// Each instance keeps its own handle on the side-car, so the copy opens a
    /// fresh stream on the same path; the original's stream is not shared.
    /// When the side-car cannot be reopened the copy has no stream, and its
    /// reads fail with [`Error::InvalidValue`]; the source's copy has a failed
    /// stream then.
    fn clone(&self) -> Self {
        let ifs = self
            .ifs
            .as_ref()
            .and_then(|_| File::open(&self.filename_cached).ok())
            .map(BufReader::new);
        Self {
            meta_ms_experiment: self.meta_ms_experiment.clone(),
            ifs,
            filename: self.filename.clone(),
            filename_cached: self.filename_cached.clone(),
            spectra_index: self.spectra_index.clone(),
            chrom_index: self.chrom_index.clone(),
        }
    }
}

impl CachedmzML {
    /// An empty instance with no associated file (source default constructor).
    pub fn new() -> Self {
        Self::default()
    }

    /// Construct and load metadata and side-car index from `filename` in one
    /// step: default construction followed by [`Self::load`] (source
    /// `CachedmzML(const std::string&)`).
    ///
    /// # Errors
    ///
    /// As [`Self::load`].
    pub fn open(filename: impl AsRef<Path>) -> Result<Self> {
        let mut map = Self::new();
        Self::load(filename, &mut map)?;
        Ok(map)
    }

    /// Read one spectrum from the side-car file (source `getSpectrum`).
    ///
    /// The metadata of the spectrum (instrument settings, native ID, ...) is
    /// taken from the experiment loaded by [`Self::load`]; the MS level,
    /// retention time, peaks and float data arrays are read from the side-car
    /// and replace or extend it as
    /// [`CachedMzMLHandler::read_spectrum`] describes.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `id` is not below
    /// [`Self::nr_spectra`] (a precondition the source checks in debug builds
    /// only) or nothing is loaded; [`Error::Io`] when seeking fails, where the
    /// source throws `Exception::ParseError` ("Error while changing position of
    /// input stream pointer."); and the errors of
    /// [`CachedMzMLHandler::read_spectrum`].
    pub fn spectrum(&mut self, id: usize) -> Result<MSSpectrum> {
        let offset = *self
            .spectra_index
            .get(id)
            .filter(|_| id < self.meta_ms_experiment.spectra.len())
            .ok_or_else(|| {
                Error::InvalidValue(format!("cached mzML spectrum {id} does not exist"))
            })?;
        let ifs = self.stream()?;
        ifs.seek(SeekFrom::Start(offset))?;
        let mut spectrum = self.meta_ms_experiment.spectra[id].clone();
        CachedMzMLHandler::read_spectrum(&mut spectrum, self.stream()?)?;
        Ok(spectrum)
    }

    /// Read one chromatogram from the side-car file (source
    /// `getChromatogram`), completing its metadata from the loaded experiment
    /// as [`CachedMzMLHandler::read_chromatogram`] describes.
    ///
    /// # Errors
    ///
    /// As [`Self::spectrum`], for [`Self::nr_chromatograms`].
    pub fn chromatogram(&mut self, id: usize) -> Result<MSChromatogram> {
        let offset = *self
            .chrom_index
            .get(id)
            .filter(|_| id < self.meta_ms_experiment.chromatograms.len())
            .ok_or_else(|| {
                Error::InvalidValue(format!("cached mzML chromatogram {id} does not exist"))
            })?;
        let ifs = self.stream()?;
        ifs.seek(SeekFrom::Start(offset))?;
        let mut chromatogram = self.meta_ms_experiment.chromatograms[id].clone();
        CachedMzMLHandler::read_chromatogram(&mut chromatogram, self.stream()?)?;
        Ok(chromatogram)
    }

    fn stream(&mut self) -> Result<&mut BufReader<File>> {
        self.ifs
            .as_mut()
            .ok_or_else(|| Error::InvalidValue("no cached mzML file is open".into()))
    }

    /// Path of the metadata mzML file (the source's protected `filename_`);
    /// empty before [`Self::load`].
    pub fn filename(&self) -> &Path {
        &self.filename
    }

    /// Path of the binary side-car (the source's protected
    /// `filename_cached_`); empty before [`Self::load`].
    pub fn cached_filename(&self) -> &Path {
        &self.filename_cached
    }

    /// Number of spectra in the loaded experiment (source `getNrSpectra`).
    pub fn nr_spectra(&self) -> usize {
        self.meta_ms_experiment.spectra.len()
    }

    /// Number of chromatograms in the loaded experiment (source
    /// `getNrChromatograms`).
    pub fn nr_chromatograms(&self) -> usize {
        self.meta_ms_experiment.chromatograms.len()
    }

    /// Experiment-level metadata loaded from the `.mzML` file (source
    /// `getMetaData`): instrument and acquisition settings and the
    /// per-spectrum and per-chromatogram headers. Empty until [`Self::load`].
    pub fn meta_data(&self) -> &MSExperiment {
        &self.meta_ms_experiment
    }

    /// Write `map` to disk as a cached mzML pair (source `store`): the binary
    /// side-car at `filename` + `.cached` holding the payload, and the mzML
    /// file at `filename` holding the metadata, each spectrum and chromatogram
    /// marked with the `cached_data` data-processing entry of
    /// [`CachedMzMLHandler::write_metadata`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when either output file cannot be created (the
    /// source's `Exception::UnableToCreateFile`), and the errors of the mzML
    /// writer.
    pub fn store(filename: impl AsRef<Path>, map: &MSExperiment) -> Result<()> {
        let filename = filename.as_ref();
        let handler = CachedMzMLHandler::new();
        handler.write_memdump(map, cached_name(filename))?;
        handler.write_metadata_x(map, filename, true)
    }

    /// Load a cached mzML pair into `map` (source `load`): the metadata from
    /// `filename` and the side-car index from `filename` + `.cached`;
    /// subsequent [`Self::spectrum`] and [`Self::chromatogram`] calls fetch the
    /// payload on demand. Previous contents of `map` are replaced.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when either file cannot be opened (the source's
    /// `Exception::FileNotFound`), and [`Error::Parse`] or the mzML reader's
    /// errors when parsing the metadata or building the side-car index fails.
    /// Unlike the source, `map` is unchanged on an error; the source may have
    /// replaced the index and names by then.
    pub fn load(filename: impl AsRef<Path>, map: &mut CachedmzML) -> Result<()> {
        let filename = filename.as_ref();
        let filename_cached = cached_name(filename);
        // Create the index from the given file.
        let mut cache = CachedMzMLHandler::new();
        cache.create_memdump_index(&filename_cached)?;
        // Open the file stream.
        let ifs = BufReader::new(File::open(&filename_cached)?);
        // Load the metadata from disk.
        let meta = crate::format::mzml::load(filename)?;
        *map = Self {
            meta_ms_experiment: meta,
            ifs: Some(ifs),
            filename: filename.to_path_buf(),
            filename_cached,
            spectra_index: cache.spectra_index().to_vec(),
            chrom_index: cache.chromatogram_index().to_vec(),
        };
        Ok(())
    }
}

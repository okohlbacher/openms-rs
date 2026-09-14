// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The imzML file adapter: `load`, `store` and the imaging geometry builders.
//!
//! Ports `FORMAT/ImzMLFile.h` with `FORMAT/ImzMLFile.cpp`. See
//! `docs/IMZML_FILE_SUPPORT.md`.
//!
//! This is the user-facing half of the imzML family. An imzML dataset is two
//! files: an `.imzML` carrying mzML 1.1.0 XML with IMS ontology terms, and a
//! companion `.ibd` carrying every binary array, prefixed by a 16-byte UUID
//! that must equal the XML's `IMS:1000080`. Per spectrum the XML declares the
//! byte offset (`IMS:1000102`), element count (`IMS:1000103`) and stored byte
//! length (`IMS:1000104`) of that pixel's arrays inside the `.ibd`, plus the
//! pixel's image coordinates (`IMS:1000050` / `IMS:1000051` / `IMS:1000052`).
//!
//! The three loading modes of the header are all here:
//!
//! * **In memory** —
//!   [`ImzMLFile::load`](crate::format::imzml_file::ImzMLFile::load) decodes
//!   every pixel and derives the pixel grid from the parsed coordinate index,
//!   giving zero-based `(x, y)` random access in RAM.
//! * **Streaming** —
//!   [`ImzMLFile::load_into_consumer`](crate::format::imzml_file::ImzMLFile::load_into_consumer)
//!   hands each decoded spectrum to an
//!   [`MSDataConsumer`](crate::interfaces::MSDataConsumer) and retains nothing.
//! * **Index only** —
//!   [`ImzMLFile::load_spectra_index`](crate::format::imzml_file::ImzMLFile::load_spectra_index)
//!   parses the XML and returns the byte-offset index for
//!   [`OnDiscImzMLExperiment`](crate::kernel::on_disc_imzml_experiment::OnDiscImzMLExperiment).
//!
//! # What this module adds, and what it reuses
//!
//! Nothing here parses XML or touches the `.ibd`. The index scan, the `.ibd`
//! reads and their preflight live in
//! [`imzml_handler`](crate::format::imzml_handler); the writer and its `.ibd`
//! layout live in [`imzml_writer`](crate::format::imzml_writer); the pixel
//! grid, the regions and the ion image live in
//! [`on_disc_imzml_experiment`](crate::kernel::on_disc_imzml_experiment),
//! whose
//! [`build_imaging_geometry`](crate::kernel::on_disc_imzml_experiment::build_imaging_geometry)
//! is re-exported here because it is a public static member of *this* header.
//! What this module owns is the composition: which ceilings apply, which
//! `PeakFileOptions` are honoured, how the dataset metadata is mirrored onto an
//! [`MSExperiment`](crate::kernel::MSExperiment), and the two geometry-driven
//! conversions the header declares.
//!
//! # Bounded work
//!
//! Every offset and length that reaches the `.ibd` came out of the XML, so a
//! hostile `.imzML` is an index into a second file that the attacker also
//! controls. The scan and every read are preflighted one layer down against
//! [`ImzMLReadLimits`](crate::format::imzml_handler::ImzMLReadLimits) and the
//! measured `.ibd` length. On top of that, a whole image load sums the declared
//! element counts of the parsed index and refuses the dataset **before decoding
//! the first array** when the total exceeds
//! [`max_loaded_peaks`](crate::format::imzml_file::ImzMLLoadLimits::max_loaded_peaks)
//! or
//! [`max_loaded_float_values`](crate::format::imzml_file::ImzMLLoadLimits::max_loaded_float_values),
//! so a file that declares nine spectra of two billion peaks each costs one XML
//! scan and no allocation.
//!
//! Nothing here starts a thread. The source parallelises its array decode with
//! `#pragma omp parallel for` inside `ImzMLHandler`; this port is serial
//! throughout, so a large image decodes on one core.

use crate::concept::progress_logger::{ProgressLogType, ProgressLogger};
use crate::format::PeakFileOptions;
use crate::format::imzml_handler::{
    ImzMLHandler, ImzMLIndex, ImzMLMeta, ImzMLReadLimits, ImzMLSpectrumIndex, SkippedAux,
    UuidStatus, infer_ibd_path,
};
use crate::format::imzml_writer::{
    ImzMLWriteOptions, StoreReport, apply_store_options, store_with_options,
};
use crate::interfaces::MSDataConsumer;
use crate::kernel::on_disc_imzml_experiment::{
    GeometryReport, ImagingGeometry, IonImage, MAX_LISTED_PROBLEMS, build_imaging_geometry,
};
use crate::kernel::ranges::RangeBase;
use crate::kernel::{MSExperiment, MSSpectrum};
use crate::metadata::{MetaInfo, MetaValue};
use crate::{Error, Result};
use std::path::Path;

/// Build the pixel grid straight from a parsed imzML index.
///
/// Source `static buildImagingGeometry(const std::vector<ImzMLSpectrumIndex>&,
/// const ImzMLMeta&, MSImagingGeometry&)`, a public static member of *this*
/// header whose Rust implementation lives next to the grid type it fills, in
/// [`on_disc_imzml_experiment`](crate::kernel::on_disc_imzml_experiment),
/// because the on-disc reader needs it too. This alias is the name the header
/// declares it under; the two paths are the same function.
pub use crate::kernel::on_disc_imzml_experiment::build_imaging_geometry as build_imaging_geometry_from_index;

/// Resource ceilings for a whole-dataset imzML load.
///
/// The source has none: `ImzMLFile::loadImpl_` parses and decodes whatever the
/// file declares. These bound the work an `.imzML` can ask for before any
/// array is read.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ImzMLLoadLimits {
    /// Ceilings for the XML scan and each individual `.ibd` read, passed
    /// straight to [`ImzMLHandler::open_with_limits`].
    pub index: ImzMLReadLimits,
    /// Maximum peaks summed over every spectrum of one in-memory load,
    /// checked against the index's declared `IMS:1000103` counts before the
    /// first array is read.
    pub max_loaded_peaks: usize,
    /// Maximum auxiliary float values summed over every spectrum of one
    /// in-memory load, checked in the same preflight.
    pub max_loaded_float_values: usize,
}

impl Default for ImzMLLoadLimits {
    fn default() -> Self {
        Self {
            index: ImzMLReadLimits::default(),
            max_loaded_peaks: 500_000_000,
            max_loaded_float_values: 500_000_000,
        }
    }
}

/// An experiment paired with the pixel grid its spectra were acquired on.
///
/// This is the return type of [`ImzMLFile::load`] and the argument of
/// [`ImzMLFile::store_imaging`], standing in for source
/// `MSImagingExperiment`. That class lives in `IMAGING/MSImagingExperiment.h`,
/// which is a **different header with its own class test** and is not ported;
/// this type reproduces exactly the surface `ImzMLFile`'s two
/// `MSImagingExperiment`-typed members need, so that neither member has to be
/// left unported. `docs/IMZML_FILE_SUPPORT.md` records which members of that
/// header are covered here and which are not, and a later IMAGING package owns
/// the header itself.
///
/// The geometry is the source of truth for coordinates: a pixel binds a
/// zero-based `(x, y)` to a *linear spectrum index*, so every spectrum stays
/// reachable by index even when it never reached the grid.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ImagingExperiment {
    experiment: MSExperiment,
    geometry: ImagingGeometry,
}

impl ImagingExperiment {
    /// An experiment with no spectra and an empty grid.
    ///
    /// Source default constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// Bind `experiment` to `geometry` without checking that the pixels
    /// reference existing spectra.
    ///
    /// Source has no such constructor; it pairs `setMSExperiment` with
    /// `setGeometry` and leaves `validate()` to the caller, which is what this
    /// reproduces. Call [`validate`](Self::validate) when the two halves come
    /// from different places.
    pub fn from_parts(experiment: MSExperiment, geometry: ImagingGeometry) -> Self {
        Self {
            experiment,
            geometry,
        }
    }

    /// Read access to the spectra.
    ///
    /// Source `getMSExperiment() const`.
    pub fn ms_experiment(&self) -> &MSExperiment {
        &self.experiment
    }

    /// Mutable access to the spectra.
    ///
    /// Source `getMSExperiment()`. The grid is not revalidated, so removing a
    /// spectrum through this reference can leave a pixel dangling; ask
    /// [`validate`](Self::validate) afterwards.
    pub fn ms_experiment_mut(&mut self) -> &mut MSExperiment {
        &mut self.experiment
    }

    /// Replace the spectra, keeping the grid.
    ///
    /// Source `setMSExperiment(MSExperiment)`, which likewise keeps the
    /// geometry. Its `operator=(MSExperiment)` overload instead *clears* the
    /// geometry, because assigning a new experiment invalidates every pixel's
    /// spectrum index; that overload is [`set_ms_experiment_and_clear`] here.
    ///
    /// [`set_ms_experiment_and_clear`]: Self::set_ms_experiment_and_clear
    pub fn set_ms_experiment(&mut self, experiment: MSExperiment) {
        self.experiment = experiment;
    }

    /// Replace the spectra and clear the grid.
    ///
    /// Source `operator=(MSExperiment exp)`, whose documentation says the
    /// pixel-to-spectrum bindings cannot survive a wholesale replacement, so
    /// the geometry is cleared. Rust has no assignment operator to overload,
    /// so this is a named method.
    pub fn set_ms_experiment_and_clear(&mut self, experiment: MSExperiment) {
        self.experiment = experiment;
        self.geometry.clear();
    }

    /// Read access to the pixel grid.
    ///
    /// Source `getGeometry() const`.
    pub fn geometry(&self) -> &ImagingGeometry {
        &self.geometry
    }

    /// Mutable access to the pixel grid, for adding regions.
    ///
    /// Source `getGeometry()`. The upstream suite uses exactly this to attach
    /// the same region to an in-memory and an on-disc experiment before
    /// comparing their ion images.
    pub fn geometry_mut(&mut self) -> &mut ImagingGeometry {
        &mut self.geometry
    }

    /// Replace the pixel grid.
    ///
    /// Source `setGeometry(MSImagingGeometry)`.
    pub fn set_geometry(&mut self, geometry: ImagingGeometry) {
        self.geometry = geometry;
    }

    /// Pixels in the grid, as source `getNumberOfPixels()`.
    ///
    /// This is at most [`number_of_spectra`](Self::number_of_spectra): a
    /// spectrum whose coordinates were out of grid, below 1, duplicated or on
    /// another plane is not in the grid.
    pub fn number_of_pixels(&self) -> usize {
        self.geometry.number_of_pixels()
    }

    /// Spectra in the underlying experiment, as source
    /// `getNumberOfSpectra()`.
    pub fn number_of_spectra(&self) -> usize {
        self.experiment.spectra.len()
    }

    /// Whether a spectrum was acquired at zero-based pixel `(x, y)`.
    ///
    /// Source `hasPixel(UInt, UInt)`. Coordinates here are zero-based grid
    /// indices, not the 1-based coordinates the `.imzML` stores.
    pub fn has_pixel(&self, x: u32, y: u32) -> bool {
        self.geometry.has_pixel(x, y)
    }

    /// The spectrum bound to zero-based pixel `(x, y)`.
    ///
    /// Source `getSpectrum(UInt x, UInt y) const`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no pixel exists at that coordinate, where
    /// the source throws `Exception::ElementNotFound`, and when the pixel
    /// references a spectrum index at or above
    /// [`number_of_spectra`](Self::number_of_spectra), where the source throws
    /// `Exception::InvalidValue`.
    pub fn spectrum(&self, x: u32, y: u32) -> Result<&MSSpectrum> {
        let index = self.pixel_spectrum_index(x, y)?;
        self.experiment
            .spectra
            .get(index)
            .ok_or_else(|| dangling_pixel(x, y, index))
    }

    /// Mutable access to the spectrum bound to zero-based pixel `(x, y)`.
    ///
    /// Source `getSpectrum(UInt x, UInt y)`.
    ///
    /// # Errors
    ///
    /// As [`spectrum`](Self::spectrum).
    pub fn spectrum_mut(&mut self, x: u32, y: u32) -> Result<&mut MSSpectrum> {
        let index = self.pixel_spectrum_index(x, y)?;
        let spectra = self.experiment.spectra.len();
        if index >= spectra {
            return Err(dangling_pixel(x, y, index));
        }
        Ok(&mut self.experiment.spectra[index])
    }

    /// Check that every pixel references an existing spectrum.
    ///
    /// Source `validate() const`, which throws on the first dangling
    /// reference. [`ImzMLFile::load`] calls it, as the source's `load` does.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] naming the first pixel whose spectrum index is
    /// at or above [`number_of_spectra`](Self::number_of_spectra), matching
    /// the source's `Exception::InvalidValue` "Pixel references missing
    /// spectrum".
    pub fn validate(&self) -> Result<()> {
        let spectra = self.experiment.spectra.len();
        for pixel in self.geometry.pixels() {
            if pixel.spectrum_index >= spectra {
                return Err(dangling_pixel(pixel.x, pixel.y, pixel.spectrum_index));
            }
        }
        Ok(())
    }

    /// Spectrum indices of the acquired pixels belonging to one region.
    ///
    /// Source `getRegionSpectrumIndices(Size)`, which delegates to the
    /// geometry. The values are indices into
    /// [`ms_experiment`](Self::ms_experiment), not pixel positions.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `region_id` is unknown, where the source
    /// throws `Exception::ElementNotFound`.
    pub fn region_spectrum_indices(&self, region_id: usize) -> Result<Vec<usize>> {
        self.geometry.region_spectrum_indices(region_id)
    }

    /// Sum the intensities inside `[mz - dm, mz + dm]`, `dm = mz *
    /// tolerance_ppm * 1e-6`, at every acquired pixel.
    ///
    /// Source `extractIonImage(double, double) const`. A pixel with no peak in
    /// the window is valid with intensity 0; a pixel absent from the grid stays
    /// invalid. The returned image's m/z range is the window.
    ///
    /// # Arguments
    ///
    /// * `mz` — window centre; must be finite and `>= 0`.
    /// * `tolerance_ppm` — half-window width in ppm; must be finite and
    ///   `>= 0`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `mz` or `tolerance_ppm` is negative or not
    /// finite, when a pixel references a missing spectrum, or when the grid
    /// exceeds the image ceiling; the source throws `Exception::InvalidValue`
    /// for the first two and `Exception::IndexOverflow` for a pixel outside the
    /// declared dimensions. [`Error::UnsortedData`] when a referenced spectrum
    /// is not sorted by m/z. The source states sortedness as an unchecked
    /// precondition — "Phase 1 callers must ensure it manually" — and would
    /// silently return a wrong sum; the window search here checks, because the
    /// cost is one pass and the alternative is an ion image that is quietly
    /// wrong. Spectra that [`ImzMLFile::load`] produced are sorted whenever
    /// [`PeakFileOptions::sort_spectra_by_mz`] is set, which is its default.
    pub fn extract_ion_image(&self, mz: f64, tolerance_ppm: f64) -> Result<IonImage> {
        self.extract(mz, tolerance_ppm, None)
    }

    /// The same extraction restricted to the acquired pixels of one region.
    ///
    /// Source `extractIonImage(double, double, Size) const`. Pixels outside
    /// the region stay invalid, so the image keeps the grid's full dimensions.
    ///
    /// # Arguments
    ///
    /// * `mz` — window centre; finite and `>= 0`.
    /// * `tolerance_ppm` — half-window width in ppm; finite and `>= 0`.
    /// * `region_id` — the region, as added through
    ///   [`geometry_mut`](Self::geometry_mut).
    ///
    /// # Errors
    ///
    /// As [`extract_ion_image`](Self::extract_ion_image), plus
    /// [`Error::InvalidValue`] when `region_id` is unknown, where the source
    /// throws `Exception::ElementNotFound`.
    pub fn extract_ion_image_in_region(
        &self,
        mz: f64,
        tolerance_ppm: f64,
        region_id: usize,
    ) -> Result<IonImage> {
        let pixels = self.geometry.region_pixels(region_id)?;
        self.extract(mz, tolerance_ppm, Some(&pixels))
    }

    /// Source `Internal::extractIonImage`, whose only caller-specific part is
    /// how a spectrum is obtained; `positions` is `None` for the whole grid.
    fn extract(
        &self,
        mz: f64,
        tolerance_ppm: f64,
        positions: Option<&[usize]>,
    ) -> Result<IonImage> {
        if !mz.is_finite() || !tolerance_ppm.is_finite() || mz < 0.0 || tolerance_ppm < 0.0 {
            return Err(Error::InvalidValue(format!(
                "ion image mz and tolerance_ppm must be finite and non-negative: mz={mz}, \
                 tolerance_ppm={tolerance_ppm}"
            )));
        }
        let dm = mz * tolerance_ppm * 1e-6;
        let lo = mz - dm;
        let hi = mz + dm;

        let mut image = IonImage::new(self.geometry.width(), self.geometry.height())?;
        image.set_mz_range(RangeBase::from_min_max(lo, hi)?);

        let spectra = self.experiment.spectra.len();
        let count = positions.map_or_else(|| self.geometry.number_of_pixels(), <[usize]>::len);
        for step in 0..count {
            let position = positions.map_or(step, |positions| positions[step]);
            let pixel = *self.geometry.pixels().get(position).ok_or_else(|| {
                Error::InvalidValue(format!(
                    "imaging pixel {position} is not below the {} pixels of the grid",
                    self.geometry.number_of_pixels()
                ))
            })?;
            if pixel.spectrum_index >= spectra {
                return Err(dangling_pixel(pixel.x, pixel.y, pixel.spectrum_index));
            }
            let spectrum = &self.experiment.spectra[pixel.spectrum_index];
            let begin = spectrum.mz_begin(lo)?;
            let end = spectrum.mz_end(hi)?;
            let sum: f64 = spectrum.peaks[begin..end]
                .iter()
                .map(|peak| f64::from(peak.intensity))
                .sum();
            image.set_intensity(pixel.x, pixel.y, sum)?;
        }
        Ok(image)
    }

    /// The spectrum index bound to a pixel, or the source's
    /// `ElementNotFound` equivalent.
    fn pixel_spectrum_index(&self, x: u32, y: u32) -> Result<usize> {
        self.geometry
            .spectrum_index(x, y)
            .ok_or_else(|| Error::InvalidValue(format!("no imaging pixel exists at ({x},{y})")))
    }
}

impl From<MSExperiment> for ImagingExperiment {
    /// Source `explicit MSImagingExperiment(MSExperiment exp)`, which takes the
    /// spectra and leaves the geometry empty.
    fn from(experiment: MSExperiment) -> Self {
        Self {
            experiment,
            geometry: ImagingGeometry::new(),
        }
    }
}

/// What one [`ImzMLFile`] load read, and what it could not place.
///
/// The source returns `void` from every load and sends all of this to
/// `OPENMS_LOG_WARN`. Returning it is what lets a caller notice that a dataset
/// loaded with a `.ibd` whose UUID does not match its `.imzML`, or that three
/// pixels never reached the grid.
#[derive(Clone, Debug, PartialEq)]
pub struct LoadReport {
    /// Dataset-level imaging metadata as parsed, with
    /// [`ImzMLMeta::ibd_file_path`] set to the `.ibd` actually opened.
    pub meta: ImzMLMeta,
    /// Verdict of the `.ibd` UUID header check that source `verifyIbdUuid_`
    /// runs on every read path. Every value but [`UuidStatus::Match`] is a
    /// warning in the source and does not stop the load; it does not stop this
    /// load either.
    pub uuid: UuidStatus,
    /// What the geometry build could not place. Empty for
    /// [`ImzMLFile::load_into_consumer`] and
    /// [`ImzMLFile::load_spectra_index`], which build no geometry.
    pub geometry: GeometryReport,
    /// Spectra delivered: retained for [`ImzMLFile::load`], passed to the
    /// consumer for [`ImzMLFile::load_into_consumer`].
    pub spectra_loaded: usize,
    /// Spectra removed by the `PeakFileOptions` filters after decoding.
    pub spectra_filtered_out: usize,
    /// Spectra whose m/z or intensity array is not external, so this port
    /// contributed no peaks for them. The source takes those peaks from its
    /// `MzMLHandler` base's inline base64; see the module documentation.
    pub inline_peak_spectra: usize,
    /// Auxiliary arrays that produced no data array, paired with the
    /// zero-based spectrum position that carried them. At most
    /// [`MAX_LISTED_PROBLEMS`] entries, while
    /// [`Self::skipped_aux_count`] counts every one.
    pub skipped_aux: Vec<(usize, SkippedAux)>,
    /// How many auxiliary arrays were skipped in total.
    pub skipped_aux_count: usize,
    /// True when the consumer asked to stop before the last spectrum, so
    /// [`Self::spectra_loaded`] is short of the dataset's spectrum count.
    pub stopped_early: bool,
}

impl LoadReport {
    /// Whether the dataset loaded with no omission and a matching `.ibd`
    /// UUID.
    ///
    /// Native convenience: the source has no report to summarise.
    pub fn is_clean(&self) -> bool {
        self.uuid == UuidStatus::Match
            && self.geometry.is_clean()
            && self.spectra_filtered_out == 0
            && self.inline_peak_spectra == 0
            && self.skipped_aux_count == 0
            && !self.stopped_early
    }
}

/// What [`build_imaging_geometry_from_experiment`] could not place.
///
/// The index-based builder's [`GeometryReport`] covers four of the five
/// outcomes. The fifth is unique to the meta-value path: a spectrum can simply
/// carry no `imzml:x` or `imzml:y`, which the source skips **silently** — not
/// even a warning — because an `MSExperiment` may hold spectra that never came
/// from an imzML file.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct MetaGeometryReport {
    /// The four outcomes shared with the index-based builder.
    pub placement: GeometryReport,
    /// Spectra carrying neither `imzml:x` nor `imzml:y`, at most
    /// [`MAX_LISTED_PROBLEMS`] of them.
    pub without_coordinates: Vec<usize>,
    /// How many spectra were skipped for that reason.
    pub without_coordinates_count: usize,
}

impl MetaGeometryReport {
    /// Whether every spectrum of the experiment reached the grid.
    pub fn is_clean(&self) -> bool {
        self.placement.is_clean() && self.without_coordinates_count == 0
    }
}

/// File adapter for imzML 1.1.0 mass spectrometry imaging data.
///
/// Source `ImzMLFile`, which derives from `Internal::XMLFile` — contributing
/// the mzML 1.1.0 schema and `is_valid` (available with the `mzml-schema` feature) — and from
/// `ProgressLogger`, contributing [`log_type`](Self::log_type) and
/// [`set_log_type`](Self::set_log_type). Rust has no implementation
/// inheritance, so both are members here.
///
/// The `.ibd` path is derived from the `.imzML` path: a case-insensitive
/// `.imzML` suffix is replaced with `.ibd`, any other name gains `.ibd`. Every
/// load takes an explicit-`.ibd` variant for the case where the sibling is
/// missing or stale.
///
/// # Imaging metadata on the experiment
///
/// After a load the dataset fields are meta values on the experiment —
/// `imzml:imaging_mode`, `imzml:ibd_path`, `imzml:uuid`, the checksums,
/// `imzml:max_count_x` / `_y` / `_z`, the pixel sizes, the scan geometry, the
/// polarity and the array data types when the file declares them — and each
/// spectrum carries `imzml:x`, `imzml:y` and `imzml:z`. That is what makes a
/// loaded experiment storable again, and what
/// [`build_imaging_geometry_from_experiment`] reads.
///
/// # `.ibd` integrity
///
/// Every load path compares the first 16 bytes of the `.ibd` with the XML's
/// `IMS:1000080`. A mismatch, a missing or unparsable identifier, or an `.ibd`
/// shorter than 16 bytes is reported in [`LoadReport::uuid`] and does **not**
/// fail the load, so legacy and non-conformant datasets still open — the
/// source's documented choice, which it implements by logging a warning.
/// Callers that require strict conformance check that field.
///
/// # Examples
///
/// ```
/// use openms::format::imzml_file::ImzMLFile;
/// use openms::format::imzml_handler::ImagingMode;
/// use openms::format::imzml_writer::store;
/// use openms::format::PeakFileOptions;
/// use openms::kernel::{MSExperiment, MSSpectrum, Peak1D};
/// use openms::metadata::MetaValue;
/// use openms::system::file::TempDir;
///
/// // A 2 x 1 image with a shared m/z axis.
/// let mut experiment = MSExperiment::new();
/// for (x, y) in [(1u32, 1u32), (2, 1)] {
///     let mut spectrum =
///         MSSpectrum::from(vec![Peak1D::new(100.0, 10.0), Peak1D::new(200.0, 20.0)]);
///     spectrum.metadata.insert("imzml:x".into(), MetaValue::from(x));
///     spectrum.metadata.insert("imzml:y".into(), MetaValue::from(y));
///     experiment.spectra.push(spectrum);
/// }
///
/// // A directory of its own, removed when `directory` is dropped.
/// let directory = TempDir::new_in(std::env::temp_dir(), false)?;
/// let path = directory.path().join("example.imzML");
///
/// let file = ImzMLFile::new();
/// file.store(&path, &experiment)?;
///
/// let (imaging, report) = file.load(&path)?;
/// assert_eq!(report.meta.imaging_mode, Some(ImagingMode::Continuous));
/// assert_eq!(imaging.number_of_spectra(), 2);
/// assert_eq!(imaging.geometry().width(), 2);
/// assert_eq!(imaging.geometry().height(), 1);
/// // imzML coordinates are 1-based; the grid is 0-based.
/// assert!(imaging.has_pixel(1, 0));
/// assert_eq!(imaging.spectrum(1, 0)?.peaks.len(), 2);
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Clone, Debug, Default)]
pub struct ImzMLFile {
    options: PeakFileOptions,
    limits: ImzMLLoadLimits,
    write: ImzMLWriteOptions,
    log_type: ProgressLogType,
}

impl ImzMLFile {
    /// A file adapter with default options, ceilings and no progress output.
    ///
    /// Source default constructor, which registers the mzML 1.1.0 schema with
    /// its `Internal::XMLFile` base.
    pub fn new() -> Self {
        Self::default()
    }

    /// Read access to the load/store options, as source `getOptions() const`.
    pub fn options(&self) -> &PeakFileOptions {
        &self.options
    }

    /// Mutable access to the load/store options, as source `getOptions()`.
    pub fn options_mut(&mut self) -> &mut PeakFileOptions {
        &mut self.options
    }

    /// Replace the load/store options, as source `setOptions`.
    pub fn set_options(&mut self, options: PeakFileOptions) {
        self.options = options;
    }

    /// The ceilings this adapter enforces. Native; the source has none.
    pub fn limits(&self) -> ImzMLLoadLimits {
        self.limits
    }

    /// Replace the ceilings. Native; the source has none.
    pub fn set_limits(&mut self, limits: ImzMLLoadLimits) {
        self.limits = limits;
    }

    /// The writer options used by [`store`](Self::store).
    ///
    /// Native. They carry the writer's own ceilings and the shared-m/z
    /// tolerance, which the source hard-codes at `1e-5`; see
    /// [`ImzMLWriteOptions`].
    pub fn write_options(&self) -> &ImzMLWriteOptions {
        &self.write
    }

    /// Replace the writer options. Native, as
    /// [`write_options`](Self::write_options).
    pub fn set_write_options(&mut self, write: ImzMLWriteOptions) {
        self.write = write;
    }

    /// Progress reporting mode, from the source's `ProgressLogger` base.
    pub fn log_type(&self) -> ProgressLogType {
        self.log_type
    }

    /// Set the progress reporting mode, from the source's `ProgressLogger`
    /// base. `ImzMLFile::loadImpl_` copies it onto the logger it hands the
    /// handler, and `store` onto the writer's.
    pub fn set_log_type(&mut self, log_type: ProgressLogType) {
        self.log_type = log_type;
    }

    /// Load an imzML dataset into memory with zero-based pixel random access.
    ///
    /// Source `load(const std::string&, MSImagingExperiment&)`. Every pixel is
    /// decoded, the grid is derived **from the parsed coordinate index** — not
    /// by reading `imzml:x` / `imzml:y` back off the spectra — and the two are
    /// bound together and validated. `index[i]` corresponds to spectrum `i`,
    /// because both follow document order.
    ///
    /// Coordinates are converted from imzML's 1-based convention to zero-based
    /// grid indices, and only spectra with `z == 1` are placed: the grid is
    /// two-dimensional. A coordinate below 1, a pixel outside the declared
    /// `IMS:1000042` x `IMS:1000043` grid and a coordinate an earlier spectrum
    /// already claimed are each skipped and reported in
    /// [`LoadReport::geometry`]; the source warns about each and loads anyway.
    /// Every spectrum stays reachable through
    /// [`ImagingExperiment::ms_experiment`] whether or not it reached the grid.
    ///
    /// # Arguments
    ///
    /// * `imzml_path` — path to the `.imzML`. The `.ibd` is its sibling; use
    ///   [`load_with_ibd`](Self::load_with_ibd) to name it explicitly.
    ///
    /// # Errors
    ///
    /// * [`Error::Io`] when the `.imzML` or the `.ibd` cannot be opened, where
    ///   the source throws `Exception::FileNotFound`.
    /// * [`Error::Parse`] when the XML is malformed, when an IMS value does not
    ///   parse, when a spectrum declares no pixel coordinate, or when a
    ///   declared array range leaves the `.ibd`; the source throws
    ///   `Exception::ParseError` for the first three.
    /// * [`Error::Unsupported`] when an array is compressed — imzML external
    ///   arrays must carry `MS:1000576` — or when `options` carry a retention
    ///   time, MS level or precursor m/z filter; see the notes.
    /// * [`Error::InvalidValue`] when the dataset exceeds one of
    ///   [`limits`](Self::limits), or when the derived grid exceeds the image
    ///   ceiling.
    /// * [`Error::UnsortedData`] never from here:
    ///   [`PeakFileOptions::sort_spectra_by_mz`] sorts on the way in.
    ///
    /// # Notes
    ///
    /// Which `PeakFileOptions` are honoured differs from the source, because
    /// the source applies them inside `MzMLHandler` and this port does not
    /// reach that parser for an imzML file. `sort_spectra_by_mz`, the m/z and
    /// intensity peak filters, `metadata_only` and `fill_data` are honoured. A
    /// retention-time, MS-level or precursor-m/z filter is **refused** with
    /// [`Error::Unsupported`]: this port does not parse `MS:1000016` scan start
    /// time, `MS:1000511` ms level or the precursor list from an `.imzML`, so
    /// every spectrum would carry the unset defaults and such a filter would
    /// silently discard the whole dataset. Refusing is the crate's rule for a
    /// lossy operation; `docs/IMZML_FILE_SUPPORT.md` records the reader gap
    /// this comes from.
    pub fn load(&self, imzml_path: impl AsRef<Path>) -> Result<(ImagingExperiment, LoadReport)> {
        let imzml_path = imzml_path.as_ref();
        self.load_with_ibd(imzml_path, infer_ibd_path(imzml_path))
    }

    /// Load an imzML dataset into memory, reading the arrays from `ibd_path`.
    ///
    /// Source threads an explicit `.ibd` through `loadImpl_`'s
    /// `ibd_path_override`, so that the index load *and* the UUID check target
    /// the file that will actually be read rather than an inferred sibling that
    /// may be missing or stale.
    ///
    /// # Errors
    ///
    /// As [`load`](Self::load).
    pub fn load_with_ibd(
        &self,
        imzml_path: impl AsRef<Path>,
        ibd_path: impl AsRef<Path>,
    ) -> Result<(ImagingExperiment, LoadReport)> {
        let (experiment, mut report) = self.load_experiment_with_ibd(imzml_path, ibd_path)?;
        // The grid binds each pixel to a *linear spectrum index*, so `index[i]`
        // must still be spectrum `i`. Only the three refused filters can drop a
        // whole spectrum, so this cannot trigger today; it is here so that a
        // later filter which does drop one fails loudly instead of silently
        // shifting every pixel's binding.
        if report.filtered_out != 0 || experiment.spectra.len() != report.index_entries.len() {
            return Err(Error::InvalidValue(format!(
                "imzML load produced {} spectra for {} indexed pixels, so the imaging geometry \
                 cannot be bound by linear index",
                experiment.spectra.len(),
                report.index_entries.len()
            )));
        }
        // Source: the geometry comes from the parsed index, the source of
        // truth, so the in-memory and on-disc paths cannot diverge.
        let (geometry, placement) = build_imaging_geometry(&report.index_entries, &report.meta)?;
        report.geometry = placement;
        let imaging = ImagingExperiment::from_parts(experiment, geometry);
        imaging.validate()?;
        Ok((imaging, report.into_report()))
    }

    /// Load an imzML dataset as a flat [`MSExperiment`], with no grid.
    ///
    /// The upstream suite's own helper: it calls
    /// `load(path, MSImagingExperiment&)` and returns `getMSExperiment()`. The
    /// header declares no `MSExperiment` overload, and the section named
    /// `void load(const std::string& filename, MSExperiment& exp)` goes through
    /// that helper. Building no geometry is the only difference from
    /// [`load`](Self::load), so an experiment loaded here can be handed to
    /// [`build_imaging_geometry_from_experiment`].
    ///
    /// # Errors
    ///
    /// As [`load`](Self::load), minus the grid ceiling.
    pub fn load_experiment(
        &self,
        imzml_path: impl AsRef<Path>,
    ) -> Result<(MSExperiment, LoadReport)> {
        let imzml_path = imzml_path.as_ref();
        self.load_experiment_with_ibd(imzml_path, infer_ibd_path(imzml_path))
            .map(|(experiment, staged)| (experiment, staged.into_report()))
    }

    /// Stream an imzML dataset to a consumer, retaining nothing.
    ///
    /// Source `load(const std::string&, Interfaces::IMSDataConsumer&)`.
    /// `set_expected_size` is called with the indexed spectrum count and zero
    /// chromatograms, `set_experimental_settings` with the settings carrying the
    /// mirrored dataset metadata, and then each decoded spectrum is passed to
    /// `consume_spectrum`. `consume_chromatogram` is never called: imzML has no
    /// chromatograms.
    ///
    /// Delivery order is document order. The source delivers in a batch once
    /// the `spectrumList` section has been parsed rather than per
    /// `startElement`; here the whole XML is scanned for the index first and the
    /// arrays are read afterwards, which is the same guarantee — no spectrum
    /// reaches the consumer before the list is fully parsed.
    ///
    /// A consumer returning [`ControlFlow::Break`](std::ops::ControlFlow::Break) stops the load; that is
    /// recorded in [`LoadReport::stopped_early`]. The source's interface returns
    /// `void` and cannot ask to stop.
    ///
    /// # Arguments
    ///
    /// * `imzml_path` — path to the `.imzML`; the `.ibd` is its sibling.
    /// * `consumer` — receiver for the settings and the spectra.
    ///
    /// # Errors
    ///
    /// As [`load`](Self::load), plus whatever `consumer` returns.
    pub fn load_into_consumer(
        &self,
        imzml_path: impl AsRef<Path>,
        consumer: &mut dyn MSDataConsumer,
    ) -> Result<LoadReport> {
        let imzml_path = imzml_path.as_ref();
        self.load_into_consumer_with_ibd(imzml_path, infer_ibd_path(imzml_path), consumer)
    }

    /// Stream an imzML dataset to a consumer, reading the arrays from
    /// `ibd_path`.
    ///
    /// # Errors
    ///
    /// As [`load_into_consumer`](Self::load_into_consumer).
    pub fn load_into_consumer_with_ibd(
        &self,
        imzml_path: impl AsRef<Path>,
        ibd_path: impl AsRef<Path>,
        consumer: &mut dyn MSDataConsumer,
    ) -> Result<LoadReport> {
        self.refuse_unparsed_filters()?;
        let mut handler = self.open(imzml_path.as_ref(), ibd_path.as_ref())?;
        let uuid = handler.uuid_status()?;
        let meta = handler.meta().clone();
        let entries: Vec<ImzMLSpectrumIndex> = handler.index().to_vec();
        preflight(&entries, &self.limits)?;

        let mut settings = MSExperiment::new();
        attach_dataset_meta(&mut settings, &meta)?;
        consumer.set_expected_size(entries.len(), 0)?;
        consumer.set_experimental_settings(&settings.settings)?;

        let mut staged = Staged::new(meta, uuid, entries);
        let mut logger = self.logger();
        let steps = i64::try_from(staged.index_entries.len())
            .map_err(|_| Error::InvalidValue("imzML progress range exceeds i64".into()))?;
        logger.start_progress(0, steps, "loading imzML file")?;
        let streamed = self.stream(&mut handler, &mut staged, &mut logger, consumer);
        let ended = logger.end_progress(0);
        streamed?;
        ended?;
        Ok(staged.into_report())
    }

    /// Parse the XML and build the per-spectrum `.ibd` index without reading
    /// any array.
    ///
    /// Source `loadSpectraIndex(filename, meta, index, ibd_path = "")`, whose
    /// two out-parameters are the two halves of the returned [`ImzMLIndex`].
    /// The companion `.ibd` is opened — so a missing one is an error here, as
    /// in the source — and its path is recorded in
    /// [`ImzMLMeta::ibd_file_path`].
    ///
    /// The source reaches this by setting `PeakFileOptions::setFillData(false)`
    /// on the handler; this port reaches it by not decoding, which is why the
    /// returned index is unaffected by [`options`](Self::options).
    ///
    /// # Arguments
    ///
    /// * `imzml_path` — path to the `.imzML`; the `.ibd` is its sibling.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when either file cannot be opened, [`Error::Parse`] when
    /// the XML or an IMS value is malformed or a spectrum declares no pixel
    /// coordinate, and [`Error::InvalidValue`] when the file exceeds one of
    /// [`limits`](Self::limits). The source throws `Exception::FileNotFound`
    /// and `Exception::ParseError` for the first two.
    pub fn load_spectra_index(&self, imzml_path: impl AsRef<Path>) -> Result<ImzMLIndex> {
        let imzml_path = imzml_path.as_ref();
        self.load_spectra_index_with_ibd(imzml_path, infer_ibd_path(imzml_path))
    }

    /// Parse the XML and build the index, opening `ibd_path` instead of the
    /// inferred sibling.
    ///
    /// Source `loadSpectraIndex`'s fourth parameter, which
    /// `OnDiscImzMLExperiment::open(imzml, ibd)` uses.
    ///
    /// # Errors
    ///
    /// As [`load_spectra_index`](Self::load_spectra_index).
    pub fn load_spectra_index_with_ibd(
        &self,
        imzml_path: impl AsRef<Path>,
        ibd_path: impl AsRef<Path>,
    ) -> Result<ImzMLIndex> {
        let handler = self.open(imzml_path.as_ref(), ibd_path.as_ref())?;
        Ok(handler.parsed().clone())
    }

    /// Parse the XML, build the index and report the `.ibd` UUID verdict.
    ///
    /// Source `loadSpectraIndex` runs `verifyIbdUuid_` too — the check is in
    /// `loadImpl_`, which every read path goes through — and logs the result.
    /// This returns it, so the index-only path is not the one place a caller
    /// cannot ask.
    ///
    /// # Errors
    ///
    /// As [`load_spectra_index`](Self::load_spectra_index), plus
    /// [`Error::Io`] when the 16-byte header cannot be read.
    pub fn load_spectra_index_checked(
        &self,
        imzml_path: impl AsRef<Path>,
        ibd_path: impl AsRef<Path>,
    ) -> Result<(ImzMLIndex, UuidStatus)> {
        let mut handler = self.open(imzml_path.as_ref(), ibd_path.as_ref())?;
        let uuid = handler.uuid_status()?;
        Ok((handler.parsed().clone(), uuid))
    }

    /// Validate an `.imzML` against the mzML 1.1.0 XML schema.
    ///
    /// Source `isValid(filename, os)`, which delegates to
    /// `Internal::XMLFile::isValid` with the `mzML_1_10.xsd` its constructor
    /// registered, and writes the engine's messages to `os`. Here the messages
    /// are the report's [`diagnostics`], and
    /// [`is_valid`] is the source's `bool` return.
    ///
    /// imzML extends mzML with IMS CV terms only, never with new elements, so
    /// the unmodified mzML schema is the right one — which is why the source
    /// registers no imzML schema of its own.
    ///
    /// Available only with the `mzml-schema` feature, which brings in the
    /// libxml2 validator; the source always has Xerces.
    ///
    /// # Arguments
    ///
    /// * `imzml_path` — path to the `.imzML` to validate.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be read and [`Error::Parse`] when it
    /// is not well-formed XML. A well-formed file that violates the schema is
    /// not an error: it is a report whose `is_valid` is false, as the source's
    /// `false` return.
    ///
    /// [`diagnostics`]: crate::format::mzml_schema::SchemaValidationReport::diagnostics
    /// [`is_valid`]: crate::format::mzml_schema::SchemaValidationReport::is_valid
    #[cfg(feature = "mzml-schema")]
    pub fn is_valid(
        &self,
        imzml_path: impl AsRef<Path>,
    ) -> Result<crate::format::mzml_schema::SchemaValidationReport> {
        crate::format::mzml_schema::validate_schema(imzml_path)
    }

    /// Store an experiment as an imzML dataset: `imzml_path` and its `.ibd`
    /// sibling.
    ///
    /// Source `store(const std::string&, const MSExperiment&) const`, which
    /// forwards to `Internal::ImzMLWriter::store` with its own options. This
    /// forwards to [`store_with_options`] the same way, so every guarantee of
    /// that function holds: external float32 or float64 arrays selected by
    /// [`PeakFileOptions::mz_32_bit`] and
    /// [`PeakFileOptions::intensity_32_bit`], a 16-byte UUID header in the
    /// `.ibd` linked to `IMS:1000080`, continuous mode when
    /// `imzml:imaging_mode` says so or all spectra share an m/z axis and
    /// processed mode otherwise, and nothing written at all when the experiment
    /// is rejected.
    ///
    /// Each spectrum must carry `imzml:x` and `imzml:y` as 1-based imzML pixel
    /// coordinates. Use [`store_imaging`](Self::store_imaging) for an
    /// experiment whose coordinates live only in a grid.
    ///
    /// # Arguments
    ///
    /// * `imzml_path` — path to the output `.imzML`.
    /// * `exp` — experiment with spectra and optional `imzml:*` meta values.
    ///   It is never modified; the filters are applied to a clone.
    ///
    /// # Errors
    ///
    /// * [`Error::MissingInformation`] when `exp` has no spectrum left after
    ///   filtering, or a spectrum lacks `imzml:x` or `imzml:y`, as the source's
    ///   `Exception::MissingInformation`.
    /// * [`Error::InvalidValue`] when a pixel coordinate is invalid, as the
    ///   source's `Exception::InvalidValue`; when a continuous export is
    ///   requested that the spectra cannot support, where the source throws
    ///   `Exception::InvalidParameter`, for which this crate has no variant;
    ///   and when the dataset exceeds one of [`write_options`].
    /// * [`Error::Io`] when either file cannot be created or written, as the
    ///   source's `Exception::UnableToCreateFile` and, for a failed array
    ///   serialisation, `Exception::ParseError`.
    ///
    /// # Notes
    ///
    /// Spectra sharing a pixel coordinate are written out as-is and counted in
    /// [`StoreReport::duplicate_pixel_count`], matching what the loaders
    /// accept for the same dataset: a dataset that loads must be storable
    /// again.
    ///
    /// [`store_with_options`]: crate::format::imzml_writer::store_with_options
    /// [`write_options`]: Self::write_options
    pub fn store(&self, imzml_path: impl AsRef<Path>, exp: &MSExperiment) -> Result<StoreReport> {
        let mut logger = self.logger();
        store_with_options(imzml_path, exp, &self.options, &self.write, &mut logger)
    }

    /// Store an [`ImagingExperiment`] as imzML, taking the pixel coordinates
    /// and grid dimensions from its geometry.
    ///
    /// Source `store(const std::string&, const MSImagingExperiment&) const`.
    /// Unlike the [`MSExperiment`] overload this does not require `imzml:x` /
    /// `imzml:y` on the spectra: the grid is the source of truth, so any
    /// imaging experiment can be written — including one built without those
    /// meta values, as `BrukerTimsImagingFile` produces. They are synthesised
    /// from the grid onto a clone and the writer is reused, which is exactly
    /// what the source does.
    ///
    /// Dataset-level metadata already on the wrapped experiment — imaging mode,
    /// data types, scan geometry — is preserved; `imzml:max_count_x` /
    /// `_y` and the pixel sizes are overwritten from the grid whenever the grid
    /// declares them as positive.
    ///
    /// A spectrum that no pixel references keeps whatever `imzml:x` / `imzml:y`
    /// it already carried and is written with those, or is refused when it
    /// carries none. The source behaves identically: it only ever writes
    /// coordinates for the spectra the grid names.
    ///
    /// # Arguments
    ///
    /// * `imzml_path` — path to the output `.imzML`.
    /// * `exp` — imaging experiment to store; never modified.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a pixel references a spectrum index outside
    /// the experiment — the source's `Exception::InvalidValue` "imaging
    /// geometry references a spectrum index outside the experiment" — or when a
    /// 1-based coordinate derived from the grid overflows `u32`; otherwise as
    /// [`store`](Self::store).
    pub fn store_imaging(
        &self,
        imzml_path: impl AsRef<Path>,
        exp: &ImagingExperiment,
    ) -> Result<StoreReport> {
        let geometry = exp.geometry();
        let mut out = exp.ms_experiment().clone();
        let spectra = out.spectra.len();
        for pixel in geometry.pixels() {
            if pixel.spectrum_index >= spectra {
                return Err(Error::InvalidValue(format!(
                    "imaging geometry references a spectrum index outside the experiment: {}",
                    pixel.spectrum_index
                )));
            }
            // 0-based geometry -> 1-based imzML.
            let x = pixel
                .x
                .checked_add(1)
                .ok_or_else(|| coordinate_overflow(pixel.x))?;
            let y = pixel
                .y
                .checked_add(1)
                .ok_or_else(|| coordinate_overflow(pixel.y))?;
            let metadata = &mut out.spectra[pixel.spectrum_index].metadata;
            metadata.insert("imzml:x".into(), MetaValue::from(x));
            metadata.insert("imzml:y".into(), MetaValue::from(y));
            metadata
                .entry("imzml:z".into())
                .or_insert_with(|| MetaValue::from(1_u32));
        }

        let settings = &mut out.settings.metadata;
        if geometry.width() > 0 {
            settings.insert(
                "imzml:max_count_x".into(),
                MetaValue::from(geometry.width()),
            );
        }
        if geometry.height() > 0 {
            settings.insert(
                "imzml:max_count_y".into(),
                MetaValue::from(geometry.height()),
            );
        }
        if geometry.pixel_size_x() > 0.0 {
            settings.insert(
                "imzml:pixel_size_x".into(),
                MetaValue::try_from(geometry.pixel_size_x())?,
            );
        }
        if geometry.pixel_size_y() > 0.0 {
            settings.insert(
                "imzml:pixel_size_y".into(),
                MetaValue::try_from(geometry.pixel_size_y())?,
            );
        }
        self.store(imzml_path, &out)
    }

    /// Decode every pixel into a flat experiment and stage the report.
    fn load_experiment_with_ibd(
        &self,
        imzml_path: impl AsRef<Path>,
        ibd_path: impl AsRef<Path>,
    ) -> Result<(MSExperiment, Staged)> {
        self.refuse_unparsed_filters()?;
        let mut handler = self.open(imzml_path.as_ref(), ibd_path.as_ref())?;
        let uuid = handler.uuid_status()?;
        let meta = handler.meta().clone();
        let entries: Vec<ImzMLSpectrumIndex> = handler.index().to_vec();
        preflight(&entries, &self.limits)?;

        let mut staged = Staged::new(meta, uuid, entries);
        let mut experiment = MSExperiment::new();
        experiment.spectra.reserve(staged.index_entries.len());

        let mut logger = self.logger();
        let steps = i64::try_from(staged.index_entries.len())
            .map_err(|_| Error::InvalidValue("imzML progress range exceeds i64".into()))?;
        logger.start_progress(0, steps, "loading imzML file")?;
        let decoded = self.decode_all(&mut handler, &mut staged, &mut logger, &mut experiment);
        let ended = logger.end_progress(0);
        decoded?;
        ended?;

        attach_dataset_meta(&mut experiment, &staged.meta)?;
        // Source applies the PeakFileOptions filters inside MzMLHandler while
        // parsing; this port applies the ones it can honour afterwards, which
        // is the same observable result for a format whose peaks are all
        // external.
        staged.filtered_out = apply_store_options(&mut experiment, &self.options)?;
        staged.spectra_loaded = experiment.spectra.len();
        Ok((experiment, staged))
    }

    /// Decode each indexed pixel into `experiment`, recording the omissions.
    fn decode_all(
        &self,
        handler: &mut ImzMLHandler,
        staged: &mut Staged,
        logger: &mut ProgressLogger,
        experiment: &mut MSExperiment,
    ) -> Result<()> {
        let mut budget = PeakBudget::new(&self.limits);
        for position in 0..staged.index_entries.len() {
            let spectrum = self.decode_one(handler, staged, position)?;
            budget.charge(&spectrum)?;
            experiment.spectra.push(spectrum);
            logger.next_progress()?;
        }
        Ok(())
    }

    /// Decode each indexed pixel and hand it to `consumer`, retaining nothing.
    fn stream(
        &self,
        handler: &mut ImzMLHandler,
        staged: &mut Staged,
        logger: &mut ProgressLogger,
        consumer: &mut dyn MSDataConsumer,
    ) -> Result<()> {
        let mut budget = PeakBudget::new(&self.limits);
        for position in 0..staged.index_entries.len() {
            let mut spectrum = self.decode_one(handler, staged, position)?;
            budget.charge(&spectrum)?;
            // The peak filters are per-spectrum, so a streaming load can honour
            // them without materialising the dataset. A spectrum the filters
            // empty is still delivered, as the source delivers it.
            let mut one = MSExperiment::new();
            one.spectra.push(std::mem::take(&mut spectrum));
            let removed = apply_store_options(&mut one, &self.options)?;
            staged.filtered_out += removed;
            let Some(mut kept) = one.spectra.pop() else {
                logger.next_progress()?;
                continue;
            };
            staged.spectra_loaded += 1;
            let flow = consumer.consume_spectrum(&mut kept)?;
            logger.next_progress()?;
            if flow.is_break() {
                staged.stopped_early = true;
                return Ok(());
            }
        }
        Ok(())
    }

    /// One pixel, honouring `fill_data` and recording its skipped auxiliary
    /// arrays.
    fn decode_one(
        &self,
        handler: &mut ImzMLHandler,
        staged: &mut Staged,
        position: usize,
    ) -> Result<MSSpectrum> {
        if !self.options.fill_data {
            // Source sets PeakFileOptions::setFillData(false) on the handler
            // for its index-only path; honoured here for any load, so a caller
            // can have the coordinates without the arrays.
            let entry = handler.entry(position)?;
            return Ok(coordinate_only_spectrum(entry));
        }
        let decoded = handler.spectrum(position)?;
        if decoded.inline_peaks {
            staged.inline_peak_spectra += 1;
        }
        for skipped in decoded.skipped_aux {
            staged.skipped_aux_count += 1;
            if staged.skipped_aux.len() < MAX_LISTED_PROBLEMS {
                staged.skipped_aux.push((position, skipped));
            }
        }
        Ok(decoded.spectrum)
    }

    /// Open both files with this adapter's ceilings.
    fn open(&self, imzml_path: &Path, ibd_path: &Path) -> Result<ImzMLHandler> {
        ImzMLHandler::open_with_limits(imzml_path, ibd_path, self.limits.index)
    }

    /// A logger carrying this adapter's `ProgressLogger` log type.
    fn logger(&self) -> ProgressLogger {
        let mut logger = ProgressLogger::new();
        logger.set_log_type(self.log_type);
        logger
    }

    /// Refuse the filters whose inputs an `.imzML` scan does not parse.
    fn refuse_unparsed_filters(&self) -> Result<()> {
        let unsupported = if self.options.has_rt_range() {
            "a retention time range (no MS:1000016 scan start time is parsed)"
        } else if self.options.has_ms_levels() {
            "an MS level selection (no MS:1000511 ms level is parsed)"
        } else if self.options.has_precursor_mz_range() {
            "a precursor m/z range (no precursor list is parsed)"
        } else {
            return Ok(());
        };
        Err(Error::Unsupported(format!(
            "imzML load cannot honour {unsupported}; the source applies it inside MzMLHandler, \
             which this port does not reach for an .imzML. Filtering on an unparsed field would \
             discard every spectrum, so it is refused instead"
        )))
    }
}

/// Mirror the dataset-level imaging metadata onto an experiment's meta values.
///
/// Source file-static `attachImzMLMeta_`, which `loadImpl_` calls on every
/// non-index load. The grid counts and the `.ibd` path are always written; the
/// pixel sizes and image extents only when positive, and the identifier, the
/// two checksums, the two array data types and the four acquisition-geometry
/// terms only when the file declared them. That is what makes the keys a
/// faithful record: a key present means the file said so.
///
/// The keys are the ones
/// [`dataset_meta`](crate::format::imzml_writer::dataset_meta) reads back, so a
/// load followed by a store round-trips them.
///
/// # Arguments
///
/// * `exp` — experiment to annotate; its existing meta values are kept unless
///   a key collides.
/// * `meta` — the parsed dataset metadata.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a pixel size or image extent is not finite.
/// The source writes whatever `double` it parsed; the index scan one layer down
/// already refuses a non-finite dimension, so this cannot trigger for a parsed
/// dataset and exists for a hand-built [`ImzMLMeta`].
///
/// # Notes
///
/// Source writes `imzml:mz_data_type` and `imzml:int_data_type` only when the
/// parsed string is non-empty. Here the corresponding
/// [`ImzMLDataType`](crate::format::imzml_handler::ImzMLDataType) has an
/// `Unknown` variant for "no type term seen", and that is the case that writes
/// no key — the same condition under a typed spelling.
pub fn attach_dataset_meta(exp: &mut MSExperiment, meta: &ImzMLMeta) -> Result<()> {
    use crate::format::imzml_handler::ImzMLDataType;

    let settings = &mut exp.settings.metadata;
    settings.insert(
        "imzml:imaging_mode".into(),
        MetaValue::from(meta.imaging_mode.map_or("", |mode| mode.as_str())),
    );
    settings.insert(
        "imzml:ibd_path".into(),
        MetaValue::from(meta.ibd_file_path.to_string_lossy().as_ref()),
    );
    settings.insert(
        "imzml:max_count_x".into(),
        MetaValue::from(meta.max_count_x),
    );
    settings.insert(
        "imzml:max_count_y".into(),
        MetaValue::from(meta.max_count_y),
    );
    settings.insert(
        "imzml:max_count_z".into(),
        MetaValue::from(meta.max_count_z),
    );
    for (key, value) in [
        ("imzml:pixel_size_x", meta.pixel_size_x),
        ("imzml:pixel_size_y", meta.pixel_size_y),
        ("imzml:max_dim_x", meta.max_dim_x),
        ("imzml:max_dim_y", meta.max_dim_y),
    ] {
        if value > 0.0 {
            settings.insert(key.into(), MetaValue::try_from(value)?);
        }
    }
    for (key, value) in [
        ("imzml:uuid", meta.uuid.as_str()),
        ("imzml:ibd_sha1", meta.ibd_sha1.as_str()),
        ("imzml:ibd_md5", meta.ibd_md5.as_str()),
        ("imzml:scan_pattern", meta.scan_pattern.as_str()),
        ("imzml:scan_direction", meta.scan_direction.as_str()),
        (
            "imzml:line_scan_direction",
            meta.line_scan_direction.as_str(),
        ),
        ("imzml:polarity", meta.polarity.as_str()),
    ] {
        if !value.is_empty() {
            settings.insert(key.into(), MetaValue::from(value));
        }
    }
    for (key, value) in [
        ("imzml:mz_data_type", meta.mz_data_type),
        ("imzml:int_data_type", meta.int_data_type),
    ] {
        if value != ImzMLDataType::Unknown {
            settings.insert(key.into(), MetaValue::from(value.name()));
        }
    }
    Ok(())
}

/// Build the pixel grid from an experiment's `imzml:*` meta values.
///
/// Source `static buildImagingGeometry(const MSExperiment&,
/// MSImagingGeometry&)`. The grid extent comes from `imzml:max_count_x` and
/// `imzml:max_count_y` when present, and every spectrum carrying `imzml:x`,
/// `imzml:y` and an `imzml:z` of 1 is registered at the zero-based pixel its
/// 1-based coordinates name.
///
/// Prefer [`build_imaging_geometry`] when a parsed index is at hand: that is
/// the path the loaders use and it does not depend on meta values surviving a
/// round trip. This one is for an experiment already in memory — one loaded
/// through [`ImzMLFile::load_experiment`] or a `FileHandler`.
///
/// # Arguments
///
/// * `exp` — an experiment previously loaded from imzML, or one carrying the
///   same keys.
///
/// # Returns
///
/// The grid and a [`MetaGeometryReport`] naming every spectrum it could not
/// place.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `imzml:max_count_x` / `_y`, `imzml:x` /
/// `imzml:y` / `imzml:z` or `imzml:pixel_size_x` / `_y` is present with the
/// wrong type or out of range, where the source's `DataValue` conversions throw
/// `Exception::ConversionError`; and when the declared or derived grid exceeds
/// the image ceiling, which the source has no equivalent for.
///
/// # Notes
///
/// The two builders disagree on one ordering, and this reproduces it. The
/// index-based builder tests `z != 1` **first**, so a spectrum at `(0, 0, 2)`
/// counts as another plane; this meta-value builder tests the coordinates
/// first, so the same spectrum counts as a non-positive coordinate
/// (`ImzMLFile.cpp:277` versus `:397`). Both skip it, so only the report
/// differs.
///
/// A spectrum carrying only one of `imzml:x` and `imzml:y` is skipped like one
/// carrying neither: the source requires both.
///
/// The pixel size is copied only when **both** `imzml:pixel_size_x` and
/// `imzml:pixel_size_y` are present, as the source's single `&&` condition
/// requires, and it is copied even when zero or negative — where the
/// index-based builder copies it only when both are strictly positive
/// (`ImzMLFile.cpp:363` versus `:472`). That difference is reproduced, so the
/// two builders can still disagree about the pixel size of one dataset; a
/// non-finite value is the only one either refuses, through
/// [`ImagingGeometry::set_pixel_size`].
pub fn build_imaging_geometry_from_experiment(
    exp: &MSExperiment,
) -> Result<(ImagingGeometry, MetaGeometryReport)> {
    let mut geometry = ImagingGeometry::new();
    let mut report = MetaGeometryReport::default();

    let settings = &exp.settings.metadata;
    let mut width = optional_count(settings, "imzml:max_count_x")?;
    let mut height = optional_count(settings, "imzml:max_count_y")?;
    if width > 0 && height > 0 {
        geometry.set_dimensions(width, height)?;
    }

    let mut max_x = 0_u32;
    let mut max_y = 0_u32;
    for (position, spectrum) in exp.spectra.iter().enumerate() {
        let metadata = &spectrum.metadata;
        let (Some(x_raw), Some(y_raw)) = (metadata.get("imzml:x"), metadata.get("imzml:y")) else {
            // Source skips these silently: an MSExperiment may hold spectra
            // that never came from an imzML file.
            report.without_coordinates_count += 1;
            push_capped(&mut report.without_coordinates, position);
            continue;
        };
        let x_imz = coordinate(x_raw, "imzml:x")?;
        let y_imz = coordinate(y_raw, "imzml:y")?;
        if x_imz < 1 || y_imz < 1 {
            report.placement.non_positive_count += 1;
            push_capped(&mut report.placement.non_positive_coordinates, position);
            continue;
        }
        let z_imz = match metadata.get("imzml:z") {
            Some(value) => coordinate(value, "imzml:z")?,
            None => 1,
        };
        if z_imz != 1 {
            report.placement.other_plane_count += 1;
            continue;
        }
        // Both are >= 1, so the subtraction and the cast cannot wrap.
        let x = u32::try_from(x_imz - 1).map_err(|_| coordinate_range("imzml:x", x_imz))?;
        let y = u32::try_from(y_imz - 1).map_err(|_| coordinate_range("imzml:y", y_imz))?;
        if width > 0 && height > 0 && (x >= width || y >= height) {
            report.placement.out_of_grid_count += 1;
            push_capped(&mut report.placement.out_of_grid, position);
            continue;
        }
        max_x = max_x.max(x);
        max_y = max_y.max(y);
        if geometry.has_pixel(x, y) {
            report.placement.duplicate_count += 1;
            push_capped(&mut report.placement.duplicate_pixels, position);
            continue;
        }
        geometry.add_pixel(x, y, position)?;
    }

    if width == 0 && (max_x > 0 || geometry.number_of_pixels() > 0) {
        width = max_x
            .checked_add(1)
            .ok_or_else(|| Error::InvalidValue("derived image width overflows u32".into()))?;
    }
    if height == 0 && (max_y > 0 || geometry.number_of_pixels() > 0) {
        height = max_y
            .checked_add(1)
            .ok_or_else(|| Error::InvalidValue("derived image height overflows u32".into()))?;
    }
    if width > 0 && height > 0 && (geometry.width() != width || geometry.height() != height) {
        geometry.set_dimensions(width, height)?;
    }

    if let (Some(x), Some(y)) = (
        settings.get("imzml:pixel_size_x"),
        settings.get("imzml:pixel_size_y"),
    ) {
        let x = x
            .as_f64()
            .map_err(|_| meta_type("imzml:pixel_size_x", "numeric"))?;
        let y = y
            .as_f64()
            .map_err(|_| meta_type("imzml:pixel_size_y", "numeric"))?;
        geometry.set_pixel_size(x, y, "micrometer")?;
    }

    Ok((geometry, report))
}

// ---------------------------------------------------------------------------
// Internals
// ---------------------------------------------------------------------------

/// The report under construction, plus the index the geometry build needs.
struct Staged {
    meta: ImzMLMeta,
    uuid: UuidStatus,
    index_entries: Vec<ImzMLSpectrumIndex>,
    geometry: GeometryReport,
    spectra_loaded: usize,
    filtered_out: usize,
    inline_peak_spectra: usize,
    skipped_aux: Vec<(usize, SkippedAux)>,
    skipped_aux_count: usize,
    stopped_early: bool,
}

impl Staged {
    fn new(meta: ImzMLMeta, uuid: UuidStatus, index_entries: Vec<ImzMLSpectrumIndex>) -> Self {
        Self {
            meta,
            uuid,
            index_entries,
            geometry: GeometryReport::default(),
            spectra_loaded: 0,
            filtered_out: 0,
            inline_peak_spectra: 0,
            skipped_aux: Vec::new(),
            skipped_aux_count: 0,
            stopped_early: false,
        }
    }

    fn into_report(self) -> LoadReport {
        LoadReport {
            meta: self.meta,
            uuid: self.uuid,
            geometry: self.geometry,
            spectra_loaded: self.spectra_loaded,
            spectra_filtered_out: self.filtered_out,
            inline_peak_spectra: self.inline_peak_spectra,
            skipped_aux: self.skipped_aux,
            skipped_aux_count: self.skipped_aux_count,
            stopped_early: self.stopped_early,
        }
    }
}

/// Running peak total for one load, so an inline peak array counts against the
/// same ceiling as an external one.
///
/// [`preflight`] can only see the index, whose `mz_length` comes from the
/// external-array CV params. A spectrum storing its peaks inline carries the
/// length in the XML instead, so those peaks reached the caller uncounted. The
/// gap was unreachable while the reader refused inline arrays; once it decodes
/// them the only remaining bound is `max_text_bytes`, and 512 MiB of Base64 is
/// roughly 96 million `f32` peaks. Charged per spectrum, so a load fails as
/// soon as the ceiling is crossed rather than after materialising everything.
struct PeakBudget {
    seen: u64,
    ceiling: usize,
}

impl PeakBudget {
    fn new(limits: &ImzMLLoadLimits) -> Self {
        Self {
            seen: 0,
            ceiling: limits.max_loaded_peaks,
        }
    }

    fn charge(&mut self, spectrum: &MSSpectrum) -> Result<()> {
        self.seen = self
            .seen
            .checked_add(spectrum.peaks.len() as u64)
            .ok_or_else(|| load_ceiling("peaks", self.ceiling))?;
        if self.seen > self.ceiling as u64 {
            return Err(load_ceiling("peaks", self.ceiling));
        }
        Ok(())
    }
}

/// Refuse a dataset whose declared element counts exceed the load ceilings,
/// before the first array is read.
fn preflight(entries: &[ImzMLSpectrumIndex], limits: &ImzMLLoadLimits) -> Result<()> {
    let mut peaks = 0_u64;
    let mut values = 0_u64;
    for entry in entries {
        peaks = peaks
            .checked_add(entry.mz_length)
            .ok_or_else(|| load_ceiling("peaks", limits.max_loaded_peaks))?;
        for aux in &entry.aux {
            values = values
                .checked_add(aux.length)
                .ok_or_else(|| load_ceiling("auxiliary values", limits.max_loaded_float_values))?;
        }
    }
    if peaks > limits.max_loaded_peaks as u64 {
        return Err(load_ceiling("peaks", limits.max_loaded_peaks));
    }
    if values > limits.max_loaded_float_values as u64 {
        return Err(load_ceiling(
            "auxiliary values",
            limits.max_loaded_float_values,
        ));
    }
    Ok(())
}

fn load_ceiling(what: &str, ceiling: usize) -> Error {
    Error::InvalidValue(format!(
        "imzML dataset declares more {what} than the configured ceiling of {ceiling}"
    ))
}

/// A spectrum with its pixel coordinates and identifier but no arrays, as a
/// `fill_data = false` load produces.
fn coordinate_only_spectrum(entry: &ImzMLSpectrumIndex) -> MSSpectrum {
    let mut spectrum = MSSpectrum {
        native_id: entry.native_id.clone(),
        ..MSSpectrum::default()
    };
    spectrum
        .metadata
        .insert("imzml:x".into(), MetaValue::from(entry.x));
    spectrum
        .metadata
        .insert("imzml:y".into(), MetaValue::from(entry.y));
    spectrum
        .metadata
        .insert("imzml:z".into(), MetaValue::from(entry.z));
    spectrum
}

fn optional_count(settings: &MetaInfo, key: &str) -> Result<u32> {
    match settings.get(key) {
        Some(value) => {
            let raw = value.as_i64().map_err(|_| meta_type(key, "an integer"))?;
            u32::try_from(raw).map_err(|_| coordinate_range(key, raw))
        }
        None => Ok(0),
    }
}

fn coordinate(value: &MetaValue, key: &str) -> Result<i64> {
    value.as_i64().map_err(|_| meta_type(key, "an integer"))
}

fn meta_type(key: &str, expected: &str) -> Error {
    Error::InvalidValue(format!("imzML meta value '{key}' must be {expected}"))
}

fn coordinate_range(key: &str, value: i64) -> Error {
    Error::InvalidValue(format!(
        "imzML meta value '{key}' ({value}) is outside the uint32 range an imzML pixel \
         coordinate occupies"
    ))
}

fn coordinate_overflow(value: u32) -> Error {
    Error::InvalidValue(format!(
        "imaging pixel coordinate {value} has no 1-based imzML representation in uint32"
    ))
}

fn dangling_pixel(x: u32, y: u32, index: usize) -> Error {
    Error::InvalidValue(format!(
        "imaging pixel ({x},{y}) references the missing spectrum {index}"
    ))
}

fn push_capped(listed: &mut Vec<usize>, position: usize) {
    if listed.len() < MAX_LISTED_PROBLEMS {
        listed.push(position);
    }
}

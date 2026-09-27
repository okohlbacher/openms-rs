// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! File adapter for Bruker "XMass Analysis (fid)" spectra, the source
//! `XMassFile` (`FORMAT/XMassFile.h`).
//!
//! XMass Analysis is the native format of the Bruker Flex series (MALDI-TOF)
//! spectrometers. Each spectrum is stored in its own directory holding several
//! files, of which two are read:
//!
//! - **`acqus`**: the calibration that converts a time-of-flight point into an
//!   m/z, the instrument description and the acquisition method
//!   ([`AcqusHandler`]);
//! - **`fid`**: the intensities, one 4-byte integer per point
//!   ([`FidHandler`]).
//!
//! The m/z formula follows Titulaer *et al.*, "A database application for
//! pre-processing, storage and comparison of mass spectra derived from patients
//! and controls", BMC Bioinformatics 7:403 (2006); see
//! [`AcqusHandler::position`].
//!
//! Writing is not available, in the source or here.
//!
//! # Locating `acqus`
//!
//! The source derives the `acqus` path from the `fid` path by replacing its
//! last three *bytes* with `acqus` (`StringUtils::prefix(filename,
//! filename.length() - 3) + "acqus"`), whatever those bytes are. So
//! `…/1/fid` reads `…/1/acqus`, and a file that is not called `fid` still
//! loses three characters. The port keeps that rule.
//!
//! See `docs/XMASS_FILE_SUPPORT.md` for the member mapping, the native
//! differences and the evidence.

use crate::concept::progress_logger::ProgressLogType;
use crate::data_structures::DateTime;
use crate::format::acqus_handler::AcqusHandler;
use crate::format::fid_handler::FidHandler;
use crate::kernel::SpectrumType;
use crate::metadata::{
    AcquisitionInfo, AnalyzerType, DataProcessing, InletType, InstrumentSettings, IonSource,
    IonizationMethod, MassAnalyzer, Polarity, ProcessingAction, ScanMode, Software, SourceFile,
};
use crate::{Error, MSExperiment, MSSpectrum, Peak1D, Result};
use std::path::Path;
use std::sync::Arc;

/// Default ceiling on the points of one spectrum, [`LoadOptions::max_points`].
///
/// A Flex series spectrum has about 10^5 points; 2^24 (16.7 million, 256 MiB
/// of peaks) is two orders of magnitude above that.
pub const DEFAULT_MAX_POINTS: u64 = 1 << 24;

/// What [`XMassFile::load_with_options`] does with a `fid` file that holds
/// fewer points than `acqus` declares in `$TD`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum TruncatedFid {
    /// Refuse the file with [`Error::Parse`] before any peak is allocated.
    #[default]
    Refuse,
    /// Reproduce the source: the missing points are read anyway, as
    /// [`FidHandler::intensity_source`] reads them — a partial last value, then
    /// intensity 0 — so the spectrum still has `$TD` points.
    SourcePadding,
}

/// Options of [`XMassFile::load_with_options`]; all native, the source has
/// none.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LoadOptions {
    /// Handling of a `fid` shorter than `$TD` points. Default:
    /// [`TruncatedFid::Refuse`].
    pub truncated_fid: TruncatedFid,
    /// Largest `$TD` accepted, checked before any peak is allocated. Default:
    /// [`DEFAULT_MAX_POINTS`]. A negative `$TD`, which the source's `Size`
    /// cast wraps to nearly 2^64, always exceeds it.
    pub max_points: u64,
}

impl Default for LoadOptions {
    fn default() -> Self {
        Self {
            truncated_fid: TruncatedFid::Refuse,
            max_points: DEFAULT_MAX_POINTS,
        }
    }
}

impl LoadOptions {
    /// The source behaviour for a truncated `fid`
    /// ([`TruncatedFid::SourcePadding`]), with the default point ceiling.
    pub fn source() -> Self {
        Self {
            truncated_fid: TruncatedFid::SourcePadding,
            ..Self::default()
        }
    }
}

/// File adapter for XMass Analysis `fid` spectra.
///
/// The source class `XMassFile`, a `ProgressLogger`. The source reports no
/// progress from any member; the log type is kept for the inherited surface.
#[derive(Clone, Debug, Default)]
pub struct XMassFile {
    log_type: ProgressLogType,
}

impl XMassFile {
    /// A file adapter with no progress output.
    ///
    /// The source default constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// Progress reporting mode, from the source's `ProgressLogger` base.
    pub fn log_type(&self) -> ProgressLogType {
        self.log_type
    }

    /// Set the progress reporting mode, from the source's `ProgressLogger`
    /// base. No member reports progress, as in the source.
    pub fn set_log_type(&mut self, log_type: ProgressLogType) {
        self.log_type = log_type;
    }

    /// Load the spectrum of an XMass `fid` file with the default
    /// [`LoadOptions`].
    ///
    /// The source `load(const std::string& filename, MSSpectrum& spectrum)`.
    /// See [`Self::load_with_options`] for what is read and the errors.
    ///
    /// # Errors
    ///
    /// As [`Self::load_with_options`].
    pub fn load(&self, path: impl AsRef<Path>) -> Result<MSSpectrum> {
        self.load_with_options(path, &LoadOptions::default())
    }

    /// Load the spectrum of an XMass `fid` file.
    ///
    /// The source `load`, which fills a caller's spectrum after
    /// `spectrum.clear(true)`; here a new spectrum is returned, so a failure
    /// leaves nothing half-filled. The spectrum has `$TD` peaks: peak `i` has
    /// the m/z [`AcqusHandler::position`]`(i)` and the `fid` intensity of point
    /// `i`, a non-negative 32-bit integer stored as `f32`. The source converts
    /// the same way, so an intensity above 2^24 is rounded to the nearest
    /// `f32` in both.
    ///
    /// The metadata is the source's, field by field:
    ///
    /// | Field | Value |
    /// |---|---|
    /// | `rt` | 0 |
    /// | `ms_level` | 1 |
    /// | `name` | `Xmass analysis file ` + `$ID_raw`, angle brackets kept |
    /// | `spectrum_type` | profile |
    /// | `native_id` | `spectrum=xsd:` + `$ID_raw` without `<` and `>` |
    /// | instrument settings | scan mode mass spectrum, no zoom scan, polarity positive for `.IONIZATION MODE` `LD+`, negative for `LD-`, otherwise unknown |
    /// | acquisition info | method of combination `Sum of ` + `$NoSHOTS` + ` raw spectrum` |
    /// | source file | name `fid`, path the `fid` path without its last three bytes, size `4 * $TD` bytes in MiB, type `Xmass analysis file (fid)` |
    /// | data processing | one step by software `FlexControl`, version from `$FCVer` (`<flexControl 2.4.41.0>` gives `2.4.41.0`), software meta value `Acquisition method` = `$ACQMETH` without `<` and `>`, actions smoothing, baseline reduction and m/z calibration, completed now |
    ///
    /// The source also sets the spectrum comment to `no comment`. The kernel
    /// [`MSSpectrum`] has no comment field (the source `SpectrumSettings`
    /// member), so that one assignment has nowhere to go; it is the only
    /// source assignment not reproduced.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] for a path that is not UTF-8, is shorter than
    ///   three bytes, or whose last three bytes do not end on a character
    ///   boundary. The source's `StringUtils::prefix` throws
    ///   `Exception::IndexOverflow` for the short path; the other two cannot
    ///   occur for its byte strings.
    /// - Every error of [`AcqusHandler::open`] for the derived `acqus` path,
    ///   which the source opens first: a missing `acqus` is an
    ///   [`Error::Io`] of kind `NotFound`, the source's
    ///   `Exception::FileNotFound`.
    /// - The I/O error of opening the `fid` file, whose kind distinguishes the
    ///   source's `FileNotFound`, `FileNotReadable` and `IOException`.
    /// - [`Error::InvalidValue`] when `$TD` exceeds
    ///   [`LoadOptions::max_points`]; native, checked before allocation.
    /// - [`Error::Parse`] when the `fid` holds fewer than `$TD` points and
    ///   [`LoadOptions::truncated_fid`] is [`TruncatedFid::Refuse`]; checked
    ///   from the file length before allocation for a regular file, and at the
    ///   first short read otherwise.
    pub fn load_with_options(
        &self,
        path: impl AsRef<Path>,
        options: &LoadOptions,
    ) -> Result<MSSpectrum> {
        let fid_path = path.as_ref();
        let directory = directory_prefix(fid_path)?;
        let acqus = AcqusHandler::open(format!("{directory}acqus"))?;
        let mut fid = FidHandler::open(fid_path)?;

        let points = acqus.size();
        if points > options.max_points {
            return Err(Error::InvalidValue(format!(
                "acqus $TD of {points} points exceeds the ceiling of {} points",
                options.max_points
            )));
        }
        if options.truncated_fid == TruncatedFid::Refuse {
            if let Ok(metadata) = std::fs::metadata(fid_path) {
                if metadata.is_file() && metadata.len() / 4 < points {
                    return Err(Error::Parse {
                        line: 0,
                        message: format!(
                            "fid file holds {} complete points, acqus $TD declares {points}",
                            metadata.len() / 4
                        ),
                    });
                }
            }
        }
        // `points <= max_points`, which the caller chose to allocate.
        let capacity = usize::try_from(points)
            .map_err(|_| Error::InvalidValue("acqus $TD exceeds the address space".into()))?;

        let mut spectrum = MSSpectrum::default();
        spectrum.peaks.reserve_exact(capacity);
        while (spectrum.peaks.len() as u64) < points {
            let mz = acqus.position(fid.index());
            let intensity = match options.truncated_fid {
                TruncatedFid::Refuse => fid.intensity()?,
                TruncatedFid::SourcePadding => fid.intensity_source()?,
            };
            spectrum.peaks.push(Peak1D {
                mz,
                intensity: intensity as f32,
            });
        }

        spectrum.rt = 0.0;
        spectrum.ms_level = 1;
        let id_raw = acqus.param("$ID_raw");
        spectrum.name = format!("Xmass analysis file {id_raw}");
        spectrum.spectrum_type = SpectrumType::Profile;
        spectrum.native_id = format!("spectrum=xsd:{}", strip_angle_brackets(&id_raw));

        spectrum.instrument_settings = InstrumentSettings {
            scan_mode: ScanMode::MassSpectrum,
            zoom_scan: false,
            polarity: polarity(&acqus),
            ..InstrumentSettings::default()
        };
        spectrum.acquisition_info = AcquisitionInfo {
            method_of_combination: format!("Sum of {} raw spectrum", acqus.param("$NoSHOTS")),
            ..AcquisitionInfo::default()
        };
        spectrum.source_file = SourceFile {
            name: "fid".into(),
            path: directory.to_owned(),
            // The source computes `4.0 * getSize() / 1024 / 1024` in double and
            // `setFileSize(float)` narrows it.
            size_mb: (4.0 * points as f64 / 1024.0 / 1024.0) as f32,
            file_type: "Xmass analysis file (fid)".into(),
            ..SourceFile::default()
        };

        let mut software = Software {
            name: "FlexControl".into(),
            version: flex_control_version(&acqus.param("$FCVer")),
            ..Software::default()
        };
        software.cv_terms.metadata.insert(
            "Acquisition method".into(),
            strip_angle_brackets(&acqus.param("$ACQMETH")).into(),
        );
        let processing = DataProcessing {
            software,
            actions: [
                ProcessingAction::Smoothing,
                ProcessingAction::BaselineReduction,
                ProcessingAction::MzCalibration,
            ]
            .into_iter()
            .collect(),
            completion_time: Some(DateTime::now()),
            ..DataProcessing::default()
        };
        spectrum.data_processing = vec![Arc::new(processing)];
        Ok(spectrum)
    }

    /// Import the experiment-wide settings of an XMass spectrum into `exp`.
    ///
    /// The source `importExperimentalSettings(const std::string& filename,
    /// PeakMap& exp)`, which reads only the `acqus` file next to `path` and
    /// sets these fields of `exp`'s experimental settings, leaving every other
    /// field as it was:
    ///
    /// - instrument name `SPECTROMETER/DATASYSTEM`, vendor `ORIGIN`, model
    ///   `$InstrID` without `<` and `>`;
    /// - exactly one ion source (the list is cleared first): inlet direct when
    ///   `.INLET` is `DIRECT`, otherwise inlet unknown with MALDI ionization
    ///   (so a direct inlet leaves the ionization method unknown, as in the
    ///   source); polarity from `.IONIZATION MODE` as in [`Self::load`]; meta
    ///   value `MALDI target reference` = `$TgIDS` without `<` and `>`; order 0;
    /// - exactly one mass analyzer (the list is cleared first): time-of-flight
    ///   when `.SPECTROMETER TYPE` is `TOF`, otherwise unknown;
    /// - the date from `$AQ_DATE` without `<` and `>`, parsed by
    ///   [`DateTime::set`].
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] for a path the `acqus` rule cannot shorten,
    ///   as in [`Self::load_with_options`], and every error of
    ///   [`AcqusHandler::open`].
    /// - The [`DateTime::set`] error when `$AQ_DATE` is missing or does not
    ///   parse, the source's `Exception::ParseError`. The source has already
    ///   written the instrument by then; here `exp` is left unchanged on any
    ///   error.
    pub fn import_experimental_settings(
        &self,
        path: impl AsRef<Path>,
        exp: &mut MSExperiment,
    ) -> Result<()> {
        let directory = directory_prefix(path.as_ref())?;
        let acqus = AcqusHandler::open(format!("{directory}acqus"))?;

        let mut date = DateTime::default();
        date.set(&strip_angle_brackets(&acqus.param("$AQ_DATE")))?;

        let mut ion_source = IonSource::default();
        if acqus.param(".INLET") == "DIRECT" {
            ion_source.inlet_type = InletType::Direct;
        } else {
            ion_source.inlet_type = InletType::Unknown;
            ion_source.ionization_method = IonizationMethod::Maldi;
        }
        ion_source.polarity = polarity(&acqus);
        ion_source.metadata.insert(
            "MALDI target reference".into(),
            strip_angle_brackets(&acqus.param("$TgIDS")).into(),
        );
        ion_source.order = 0;

        let analyzer = MassAnalyzer {
            analyzer_type: if acqus.param(".SPECTROMETER TYPE") == "TOF" {
                AnalyzerType::Tof
            } else {
                AnalyzerType::Unknown
            },
            ..MassAnalyzer::default()
        };

        let settings = &mut exp.settings;
        let instrument = &mut settings.instrument;
        instrument.name = acqus.param("SPECTROMETER/DATASYSTEM");
        instrument.vendor = acqus.param("ORIGIN");
        instrument.model = strip_angle_brackets(&acqus.param("$InstrID"));
        instrument.ion_sources = vec![ion_source];
        instrument.mass_analyzers = vec![analyzer];
        settings.date_time = date;
        Ok(())
    }

    /// Store a spectrum as an XMass file: not available.
    ///
    /// The source `store(const std::string&, const MSSpectrum&)`, which always
    /// throws `Exception::NotImplemented`. Its documentation names
    /// `Exception::FileNotWritable` instead; the code is followed.
    ///
    /// # Errors
    ///
    /// Always [`Error::Unsupported`].
    pub fn store(&self, _path: impl AsRef<Path>, _spectrum: &MSSpectrum) -> Result<()> {
        Err(Error::Unsupported(
            "XMassFile::store is not implemented (XMass files are read-only)".into(),
        ))
    }
}

/// The `fid` path without its last three bytes, as `StringUtils::prefix(filename,
/// filename.length() - 3)`.
fn directory_prefix(path: &Path) -> Result<&str> {
    let text = path
        .to_str()
        .ok_or_else(|| Error::InvalidValue("XMass file path is not valid UTF-8".into()))?;
    let end = text.len().checked_sub(3).ok_or_else(|| {
        Error::InvalidValue(format!(
            "XMass file path '{text}' is shorter than the three bytes the acqus rule removes"
        ))
    })?;
    text.get(..end).ok_or_else(|| {
        Error::InvalidValue(format!(
            "XMass file path '{text}' does not end in three single-byte characters"
        ))
    })
}

/// The source lambda that removes every `<` and `>`.
fn strip_angle_brackets(text: &str) -> String {
    text.chars().filter(|c| !matches!(c, '<' | '>')).collect()
}

/// Polarity from `.IONIZATION MODE`: `LD+` positive, `LD-` negative, else unknown.
fn polarity(acqus: &AcqusHandler) -> Polarity {
    match acqus.param(".IONIZATION MODE").as_str() {
        "LD+" => Polarity::Positive,
        "LD-" => Polarity::Negative,
        _ => Polarity::Unknown,
    }
}

/// The source's `$FCVer` rule: after the last space when the value starts with
/// `<flexControl `, then before the first `>` when it ends with `>`.
fn flex_control_version(value: &str) -> String {
    let mut version = value;
    if version.starts_with("<flexControl ") {
        version = version.rsplit(' ').next().unwrap_or(version);
    }
    if version.ends_with('>') {
        version = version.split('>').next().unwrap_or(version);
    }
    version.to_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn directory_prefix_removes_three_bytes() {
        assert_eq!(directory_prefix(Path::new("a/1/fid")).unwrap(), "a/1/");
        assert_eq!(directory_prefix(Path::new("abc")).unwrap(), "");
        assert!(matches!(
            directory_prefix(Path::new("ab")),
            Err(Error::InvalidValue(_))
        ));
        assert!(matches!(
            directory_prefix(Path::new("a\u{e9}d")),
            Err(Error::InvalidValue(_))
        ));
    }

    #[test]
    fn flex_control_version_follows_the_source() {
        assert_eq!(flex_control_version("<flexControl 2.4.41.0>"), "2.4.41.0");
        assert_eq!(flex_control_version("<flexControl a b>"), "b");
        assert_eq!(flex_control_version("3.0>x>"), "3.0");
        assert_eq!(flex_control_version("<other 1.0>"), "<other 1.0");
        assert_eq!(flex_control_version(""), "");
    }

    #[test]
    fn store_is_unsupported() {
        assert!(matches!(
            XMassFile::new().store("x", &MSSpectrum::default()),
            Err(Error::Unsupported(_))
        ));
    }
}

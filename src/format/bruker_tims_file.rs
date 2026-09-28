// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reader for Bruker timsTOF `.d` directories (TDF), the source
//! `BrukerTimsFile` (`FORMAT/BrukerTimsFile.h`), behind the `bruker-tdf`
//! feature.
//!
//! Ion mobility is reported as inverse reduced mobility (1/K0, V·s/cm², the
//! source's `DriftTimeUnit::VSSC`), per peak in a float data array named
//! `raw inverse reduced ion mobility array` (`MS:1003008`, the name the
//! source's `IMDataConverter::setIMUnit` gives a VSSC array).
//!
//! # What is ported
//!
//! The frame codec and the linear calibrations live in
//! [`crate::format::tims_data_handle`]; the per-frame rational 1/K0 model in
//! [`crate::format::rational_scan2im_converter`]; the IM-axis hill
//! centroider in [`crate::format::pasef_hill_centroider`]. On top of them this
//! module ports:
//!
//! - [`Config`](crate::format::bruker_tims_file::Config), every field of the source struct with its default;
//! - the calibration strategy of the source's `openTimsDataHandle`;
//! - the frame-range / RT-range resolution (`resolveEffectiveConfig`);
//! - `FRAME` export (`loadFrames_`) with MS1 centroiding `Off` or
//!   `HillBased`, the isotopic prefilter and the hill-bounds arrays;
//! - `AUTO` export of DIA-PASEF data without frame aggregation (`loadDIA_`,
//!   raw path): MS1 frames as in `FRAME` export, MS2 frames split by
//!   diaPASEF window, with `Off` or `HillBased` MS2 centroiding;
//! - `transform` onto an [`MSDataConsumer`](crate::interfaces::MSDataConsumer).
//!
//! # What is not ported (each refused with [`Error::Unsupported`])
//!
//! - DDA-PASEF export (`loadDDA_`: `SPECTRUM` mode, and `AUTO` on DDA data);
//! - frame aggregation (`ms1_n_neighbors > 0` outside `FRAME` export,
//!   `dia_ms2_n_neighbors > 0`), which needs the source's `FrameAggregator`;
//! - `Greedy2D` MS1 centroiding, which needs the source's `FrameCentroider`;
//! - the Bruker SDK calibration (`TimsCalibrationStrategy::BrukerSdk`, and the
//!   SDK step of `Auto`): the SDK is a closed-source shared library loaded at
//!   run time, which this `unsafe`-free crate cannot load. `Auto` falls
//!   through to the rational model exactly as the source does when the SDK
//!   fails to load;
//! - `readDIAMetadata` and `loadDIAStreaming`, which are typed on the OpenSWATH
//!   `SwathMap` and `FullSwathFileConsumer`.
//!
//! See `docs/BRUKER_TIMS_SUPPORT.md` for the member mapping and the evidence.

use crate::concept::constants::C13C12_MASSDIFF_U;
use crate::concept::progress_logger::{ProgressLogType, ProgressLogger};
use crate::format::pasef_hill_centroider::{Centroid, Params as HillParams, centroid_frame};
use crate::format::tims_data_handle::{TdfFrame, TdfLimits, TimsDataHandle, sql_error};
use crate::interfaces::MSDataConsumer;
use crate::kernel::{DataArray, Precursor, SpectrumType};
use crate::metadata::{DriftTimeUnit, ExperimentalSettings, MetaValue, SourceFile};
use crate::{Error, MSExperiment, MSSpectrum, Peak1D, Result};
use rusqlite::types::ValueRef;
use rusqlite::{Connection, OpenFlags};
use std::collections::BTreeMap;
use std::path::Path;

/// Name of the per-peak 1/K0 float data array, `MS:1003008`.
pub const IM_ARRAY_NAME: &str = "raw inverse reduced ion mobility array";

/// Centroiding algorithm for IM-PASEF frames, the source
/// `Config::CentroidAlgo`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CentroidAlgo {
    /// No IM-axis centroiding: raw detector peaks pass through.
    #[default]
    Off,
    /// 2D greedy clustering in (m/z, IM) (the source's legacy path, adapted
    /// from Sage). Not ported.
    Greedy2d,
    /// IM-axis hill detection ([`centroid_frame`]), modelled on Biosaur2.
    HillBased,
}

/// Export mode, the source `Config::ExportMode`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ExportMode {
    /// Detect DDA vs DIA by querying for diaPASEF windows.
    #[default]
    Auto,
    /// Force the DDA (per-precursor) path.
    Spectrum,
    /// Raw 4D frames as spectra with per-peak IM, for MS1 and MS2 alike.
    Frame,
}

/// Scan-to-1/K0 calibration strategy, the source
/// `Config::TimsCalibrationStrategy`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum TimsCalibrationStrategy {
    /// Bruker SDK, then the rational `TimsCalibration` model, then linear.
    #[default]
    Auto,
    /// The Bruker SDK only. Not available in this port.
    BrukerSdk,
    /// The rational `TimsCalibration` model only.
    Rational,
    /// The linear `GlobalMetadata` model.
    Linear,
}

/// Pressure compensation, the source `Config::PressureCompensation`. Only the
/// Bruker SDK applies it, so it has no effect in this port beyond the
/// source's warning.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PressureCompensation {
    /// None.
    #[default]
    None,
    /// One correction for the whole analysis.
    Global,
    /// A correction per frame.
    PerFrame,
}

/// Processing and export configuration, the source `BrukerTimsFile::Config`.
///
/// Every field of the source struct, with its default. [`Self::limits`] is
/// native.
#[derive(Clone, Debug, PartialEq)]
pub struct Config {
    /// m/z recalibration tolerance in Da (0 selects the default 0.1 Da). DDA
    /// only.
    pub calibration_tolerance: f64,
    /// Enable DDA m/z recalibration. Default off.
    pub calibrate: bool,
    /// Load MS1 spectra; disable for MS2-only workflows. All export modes.
    pub load_ms1: bool,
    /// Frame Id range filter, inclusive. Frame Ids are 1-based and dense.
    /// Defaults (0, `u32::MAX`) load every frame. Intersected with the RT
    /// range.
    pub frame_id_min: u32,
    /// Upper end of the frame Id range, inclusive.
    pub frame_id_max: u32,
    /// Retention time range filter in seconds, inclusive, resolved to a frame
    /// Id range through `Frames.Time`. Defaults (0, +inf) impose no filter.
    pub rt_min_sec: f64,
    /// Upper end of the RT range, inclusive.
    pub rt_max_sec: f64,
    /// MS1 centroiding algorithm. Left at `Off`, both
    /// [`Self::ms1_centroid_mz_ppm`] and [`Self::ms1_centroid_im_pct`] above 0
    /// select `Greedy2d` implicitly (the source's back-compatibility rule).
    pub ms1_centroid_algo: CentroidAlgo,
    /// MS2 centroiding algorithm for DIA-PASEF and DDA-PASEF. Takes precedence
    /// over [`Self::dia_ms2_centroid`].
    pub ms2_centroid_algo: CentroidAlgo,
    /// MS1 centroiding m/z tolerance in ppm (Greedy2D and HillBased); 0
    /// disables MS1 centroiding. Default 10.
    pub ms1_centroid_mz_ppm: f32,
    /// Greedy2D only: MS1 IM tolerance in percent (0 disables). Default 0.
    pub ms1_centroid_im_pct: f32,
    /// Greedy2D only: cap on centroided peaks per MS1 spectrum. Default
    /// 100000.
    pub ms1_centroid_max_peaks: i32,
    /// HillBased MS2 m/z tolerance in ppm; 0 disables HillBased MS2 with a
    /// warning. Default 20.
    pub ms2_centroid_mz_ppm: f32,
    /// HillBased valley factor, shared by MS1 and MS2. Default 1.3.
    pub centroid_valley_factor: f64,
    /// HillBased MS1: minimum IM scans per hill. Default 1.
    pub ms1_centroid_min_hill_length: usize,
    /// HillBased MS2: minimum IM scans per hill. Default 2 (DIA-tuned; use 1
    /// for DDA).
    pub ms2_centroid_min_hill_length: usize,
    /// HillBased: largest run of empty IM scans a hill may bridge. Default 0.
    pub centroid_max_scan_gap: usize,
    /// HillBased (MS1 and DIA-MS2): attach the four hill bounding-box arrays
    /// (`im lower bound`, `im upper bound`, `m/z lower bound`,
    /// `m/z upper bound`). Default off.
    pub expose_hill_bounds: bool,
    /// MS1 and DIA-MS2 isotopic-partner prefilter. Default off.
    pub isotopic_prefilter: bool,
    /// Isotopic prefilter tolerance in ppm. Default 50.
    pub isotopic_prefilter_tol_ppm: f64,
    /// DIA MS2 frame aggregation: neighbours on each side (0 disables).
    /// Default 0. Aggregation is not ported.
    pub dia_ms2_n_neighbors: i32,
    /// DIA MS2 denoising support after aggregation. Default 1.
    pub dia_ms2_min_support: i32,
    /// Legacy DIA MS2 2D peak picking (implicit Greedy2D). Default off.
    pub dia_ms2_centroid: bool,
    /// MS1 frame aggregation: neighbours on each side (0 disables); ignored in
    /// `Frame` export. Default 0. Aggregation is not ported.
    pub ms1_n_neighbors: i32,
    /// MS1 denoising support after aggregation. Default 0.
    pub ms1_min_support: i32,
    /// MS1 aggregation RT cap in seconds (0 = none). Default 0.
    pub ms1_max_rt_distance_sec: f64,
    /// Export mode. Default `Auto`.
    pub export_mode: ExportMode,
    /// Scan-to-1/K0 calibration strategy. Default `Auto`.
    pub tims_calibration_strategy: TimsCalibrationStrategy,
    /// Pressure compensation (Bruker SDK only). Default `None`.
    pub pressure_compensation: PressureCompensation,
    /// Path of the Bruker SDK library; empty reads `OPENMS_BRUKER_SDK_PATH`.
    pub bruker_sdk_path: String,
    /// Native resource ceilings of the frame reader.
    pub limits: TdfLimits,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            calibration_tolerance: 0.0,
            calibrate: false,
            load_ms1: true,
            frame_id_min: 0,
            frame_id_max: u32::MAX,
            rt_min_sec: 0.0,
            rt_max_sec: f64::INFINITY,
            ms1_centroid_algo: CentroidAlgo::Off,
            ms2_centroid_algo: CentroidAlgo::Off,
            ms1_centroid_mz_ppm: 10.0,
            ms1_centroid_im_pct: 0.0,
            ms1_centroid_max_peaks: 100_000,
            ms2_centroid_mz_ppm: 20.0,
            centroid_valley_factor: 1.3,
            ms1_centroid_min_hill_length: 1,
            ms2_centroid_min_hill_length: 2,
            centroid_max_scan_gap: 0,
            expose_hill_bounds: false,
            isotopic_prefilter: false,
            isotopic_prefilter_tol_ppm: 50.0,
            dia_ms2_n_neighbors: 0,
            dia_ms2_min_support: 1,
            dia_ms2_centroid: false,
            ms1_n_neighbors: 0,
            ms1_min_support: 0,
            ms1_max_rt_distance_sec: 0.0,
            export_mode: ExportMode::Auto,
            tims_calibration_strategy: TimsCalibrationStrategy::Auto,
            pressure_compensation: PressureCompensation::None,
            bruker_sdk_path: String::new(),
            limits: TdfLimits::default(),
        }
    }
}

impl Config {
    /// Frame Id inside `[frame_id_min, frame_id_max]`, the source's
    /// `frameInRange`.
    pub fn frame_in_range(&self, id: u32) -> bool {
        id >= self.frame_id_min && id <= self.frame_id_max
    }

    /// The MS1 algorithm actually used, the source's `effectiveMS1Algo`:
    /// `Off` when an algorithm is chosen but `ms1_centroid_mz_ppm <= 0`, or
    /// Greedy2D without `ms1_centroid_im_pct > 0`; with `Off` chosen, both
    /// tolerances above 0 select Greedy2D. `emit_warnings` logs the source's
    /// warnings for the inconsistent combinations.
    pub fn effective_ms1_algo(&self, emit_warnings: bool) -> CentroidAlgo {
        if self.ms1_centroid_algo != CentroidAlgo::Off {
            if self.ms1_centroid_mz_ppm <= 0.0 {
                if emit_warnings {
                    let _ = crate::openms_log_warn!(
                        "Warning: ms1_centroid_algo selected but ms1_centroid_mz_ppm <= 0; MS1 centroiding disabled."
                    );
                }
                return CentroidAlgo::Off;
            }
            if self.ms1_centroid_algo == CentroidAlgo::Greedy2d && self.ms1_centroid_im_pct <= 0.0 {
                if emit_warnings {
                    let _ = crate::openms_log_warn!(
                        "Warning: Greedy2D MS1 centroiding requires ms1_centroid_im_pct > 0; centroiding disabled."
                    );
                }
                return CentroidAlgo::Off;
            }
            return self.ms1_centroid_algo;
        }
        let has_mz = self.ms1_centroid_mz_ppm > 0.0;
        let has_im = self.ms1_centroid_im_pct > 0.0;
        if has_mz != has_im {
            if emit_warnings {
                let _ = crate::openms_log_warn!(
                    "Warning: ms1_centroid_mz_ppm and ms1_centroid_im_pct must both be > 0 to enable MS1 centroiding. Centroiding disabled."
                );
            }
            return CentroidAlgo::Off;
        }
        if has_mz && has_im {
            CentroidAlgo::Greedy2d
        } else {
            CentroidAlgo::Off
        }
    }

    /// The MS2 algorithm actually used, the source's `effectiveMS2Algo`:
    /// HillBased needs `ms2_centroid_mz_ppm > 0`; with `Off` chosen, the legacy
    /// `dia_ms2_centroid` selects Greedy2D.
    pub fn effective_ms2_algo(&self, emit_warnings: bool) -> CentroidAlgo {
        if self.ms2_centroid_algo != CentroidAlgo::Off {
            if self.ms2_centroid_algo == CentroidAlgo::HillBased && self.ms2_centroid_mz_ppm <= 0.0
            {
                if emit_warnings {
                    let _ = crate::openms_log_warn!(
                        "Warning: HillBased DIA-MS2 centroiding requires ms2_centroid_mz_ppm > 0; centroiding disabled."
                    );
                }
                return CentroidAlgo::Off;
            }
            return self.ms2_centroid_algo;
        }
        if self.dia_ms2_centroid {
            CentroidAlgo::Greedy2d
        } else {
            CentroidAlgo::Off
        }
    }
}

/// One diaPASEF isolation window, the source's `DIAWindow`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiaWindow {
    /// `WindowGroup` (or a synthetic group per distinct window of the older
    /// `DiaFrameMsMsInfo` layout).
    pub window_group: i32,
    /// Isolation window centre, `IsolationMz`.
    pub mz_center: f64,
    /// Isolation window width, `IsolationWidth`.
    pub mz_width: f64,
    /// Lower 1/K0 bound, converted from the scan bounds with frame 1's
    /// calibration.
    pub im_lower: f64,
    /// Upper 1/K0 bound.
    pub im_upper: f64,
    /// `ScanNumBegin`.
    pub scan_begin: u32,
    /// `ScanNumEnd`.
    pub scan_end: u32,
}

/// Reader for Bruker timsTOF `.d` directories, the source class
/// `BrukerTimsFile` (a `ProgressLogger`).
#[derive(Clone, Debug, Default)]
pub struct BrukerTimsFile {
    log_type: ProgressLogType,
}

impl BrukerTimsFile {
    /// A reader with no progress output, the source default constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// Progress reporting mode, from the source's `ProgressLogger` base.
    pub fn log_type(&self) -> ProgressLogType {
        self.log_type
    }

    /// Set the progress reporting mode.
    pub fn set_log_type(&mut self, log_type: ProgressLogType) {
        self.log_type = log_type;
    }

    /// Load a `.d` directory, the source `load(path, exp, config)`; the source
    /// overload without a config is `load(path, &Config::default())`.
    ///
    /// The experiment holds, sorted by RT with every spectrum sorted by m/z
    /// (the source's closing `sortSpectra(true)`), the spectra of the chosen
    /// export mode, and one source file (`Bruker TDF`, native Id format
    /// `MS:1002818`).
    ///
    /// # Errors
    ///
    /// - [`Error::Io`] (kind `NotFound` for a missing `analysis.tdf`) and the
    ///   other errors of [`TimsDataHandle::open`], where the source throws
    ///   `Exception::FileNotReadable`;
    /// - [`Error::Unsupported`] for the unported paths listed in the module
    ///   documentation, and for `BrukerSdk` calibration;
    /// - [`Error::Io`] for `Rational` calibration when the `TimsCalibration`
    ///   table is missing or unusable (the source's `FileNotReadable`);
    /// - the frame errors of [`TimsDataHandle::read_frame`];
    /// - [`Error::InvalidValue`] for a non-finite `Frames.Time`, which the
    ///   final sort refuses.
    pub fn load(&self, path: impl AsRef<Path>, config: &Config) -> Result<MSExperiment> {
        let path = path.as_ref();
        let handle = open_handle(path, config)?;
        let db = open_db(path)?;
        let effective =
            resolve_effective_config(&db, config, handle.min_frame_id(), handle.max_frame_id());
        let is_dia = is_dia(&db);
        let mut exp = MSExperiment::new();
        self.load_mode(&handle, &db, &effective, is_dia, &mut exp)?;
        exp.settings.source_files.push(source_file(path));
        exp.sort_spectra(true)?;
        Ok(exp)
    }

    /// Feed the spectra of a `.d` directory to `consumer`, the source
    /// `transform(path, consumer, config)`.
    ///
    /// As in the source, the whole dataset is loaded first (the source notes
    /// the streaming optimisation as pending), the consumer is told the
    /// expected spectrum count the source computes from SQL (which counts
    /// empty frames that the loaders skip), the settings carry the one source
    /// file, and the spectra follow sorted as [`Self::load`] sorts them.
    ///
    /// # Errors
    ///
    /// As [`Self::load`], plus the consumer's errors. A consumer's
    /// `ControlFlow::Break` stops the feed successfully.
    pub fn transform(
        &self,
        path: impl AsRef<Path>,
        consumer: &mut dyn MSDataConsumer,
        config: &Config,
    ) -> Result<()> {
        let path = path.as_ref();
        let handle = open_handle(path, config)?;
        let db = open_db(path)?;
        let effective =
            resolve_effective_config(&db, config, handle.min_frame_id(), handle.max_frame_id());
        let is_dia = is_dia(&db);
        let expected = expected_spectra(&db, &handle, &effective, is_dia)?;
        consumer.set_expected_size(expected, 0)?;
        let mut settings = ExperimentalSettings::default();
        settings.source_files.push(source_file(path));
        consumer.set_experimental_settings(&settings)?;
        let mut exp = MSExperiment::new();
        self.load_mode(&handle, &db, &effective, is_dia, &mut exp)?;
        exp.sort_spectra(true)?;
        for spectrum in &mut exp.spectra {
            if consumer.consume_spectrum(spectrum)?.is_break() {
                break;
            }
        }
        Ok(())
    }

    fn load_mode(
        &self,
        handle: &TimsDataHandle,
        db: &Connection,
        config: &Config,
        is_dia: bool,
        exp: &mut MSExperiment,
    ) -> Result<()> {
        match config.export_mode {
            ExportMode::Frame => self.load_frames(handle, config, exp),
            ExportMode::Spectrum => Err(unsupported_dda()),
            ExportMode::Auto if !is_dia => Err(unsupported_dda()),
            ExportMode::Auto => self.load_dia(handle, db, config, exp),
        }
    }

    /// The source `loadFrames_`: every in-range frame as one spectrum, MS1
    /// (`MsMsType` 0) first, then MS2.
    fn load_frames(
        &self,
        handle: &TimsDataHandle,
        config: &Config,
        exp: &mut MSExperiment,
    ) -> Result<()> {
        if config.ms1_n_neighbors > 0 {
            let _ = crate::openms_log_warn!(
                "Warning: ms1_n_neighbors (={}) is ignored when export_mode=FRAME because FRAME mode returns raw frames. Set export_mode=AUTO or SPECTRUM to enable MS1 frame aggregation.",
                config.ms1_n_neighbors
            );
        }
        for level in 1..=2u32 {
            if level == 1 && !config.load_ms1 {
                continue;
            }
            let frames: Vec<&TdfFrame> = handle
                .frames()
                .values()
                .filter(|frame| config.frame_in_range(frame.id))
                .filter(|frame| (level == 1) == (frame.msms_type == 0))
                .collect();
            if frames.is_empty() {
                if level == 2 {
                    let _ = crate::openms_log_warn!("Warning: no MS2 frames found (skipping)");
                }
                continue;
            }
            let mut progress = self.progress();
            progress.start_progress(
                0,
                frames.len() as i64,
                &format!("Loading MS{level} frames"),
            )?;
            for (index, frame) in frames.iter().enumerate() {
                let spectrum = if level == 1 {
                    ms1_spectrum(handle, frame, config)?
                } else {
                    frame_to_spectrum(handle, frame, level)?
                };
                exp.spectra.push(spectrum);
                progress.set_progress(index as i64)?;
            }
            progress.end_progress(0)?;
        }
        Ok(())
    }

    /// The source `loadDIA_` without frame aggregation: MS1 frames, then per
    /// window group, per window, per frame the peaks inside the window's 1/K0
    /// range.
    fn load_dia(
        &self,
        handle: &TimsDataHandle,
        db: &Connection,
        config: &Config,
        exp: &mut MSExperiment,
    ) -> Result<()> {
        if config.ms1_n_neighbors > 0 && config.load_ms1 {
            return Err(Error::Unsupported(
                "BrukerTimsFile: MS1 frame aggregation (ms1_n_neighbors > 0) is not ported".into(),
            ));
        }
        if config.dia_ms2_n_neighbors > 0 {
            return Err(Error::Unsupported(
                "BrukerTimsFile: DIA MS2 frame aggregation (dia_ms2_n_neighbors > 0) is not ported"
                    .into(),
            ));
        }
        let ms1_frames: Vec<&TdfFrame> = handle
            .frames()
            .values()
            .filter(|frame| config.frame_in_range(frame.id) && frame.msms_type == 0)
            .collect();
        warn_partial_ms1_aggregation_config(config, ms1_frames.len());
        if config.load_ms1 {
            let mut progress = self.progress();
            progress.start_progress(0, ms1_frames.len() as i64, "Loading DIA-PASEF MS1 frames")?;
            for (index, frame) in ms1_frames.iter().enumerate() {
                exp.spectra.push(ms1_spectrum(handle, frame, config)?);
                progress.set_progress(index as i64)?;
            }
            progress.end_progress(0)?;
        }

        let windows = read_dia_windows(db, handle)?;
        if windows.is_empty() {
            let _ =
                crate::openms_log_warn!("Warning: DIA dataset detected but no SWATH windows found");
            return Ok(());
        }
        let group_to_frames = in_range_groups(read_frame_to_window_group(db, &windows)?, config);
        let mut group_to_windows: BTreeMap<i32, Vec<&DiaWindow>> = BTreeMap::new();
        for window in &windows {
            group_to_windows
                .entry(window.window_group)
                .or_default()
                .push(window);
        }
        let total: usize = group_to_frames
            .iter()
            .filter_map(|(group, frames)| {
                group_to_windows
                    .get(group)
                    .map(|windows| frames.len() * windows.len())
            })
            .sum();
        let algo = config.effective_ms2_algo(false);
        let hill = algo == CentroidAlgo::HillBased;
        let mut progress = self.progress();
        progress.start_progress(
            0,
            total as i64,
            if hill {
                "Loading DIA-PASEF MS2 frames (hill, no RT aggregation)"
            } else {
                "Loading DIA-PASEF MS2 frames"
            },
        )?;
        let mut done = 0i64;
        for (group, frame_ids) in &group_to_frames {
            let Some(dia_windows) = group_to_windows.get(group) else {
                continue;
            };
            for window in dia_windows {
                for &frame_id in frame_ids {
                    progress.set_progress(done)?;
                    done += 1;
                    let frame = handle.frame(frame_id).ok_or_else(|| {
                        Error::InvalidValue(format!(
                            "DiaFrameMsMsInfo names frame {frame_id}, which Frames does not hold"
                        ))
                    })?;
                    if frame.num_peaks == 0 {
                        continue;
                    }
                    if let Some(spectrum) = dia_ms2_spectrum(handle, frame, window, config, hill)? {
                        exp.spectra.push(spectrum);
                    }
                }
            }
        }
        progress.end_progress(0)?;
        Ok(())
    }

    fn progress(&self) -> ProgressLogger {
        let mut logger = ProgressLogger::new();
        logger.set_log_type(self.log_type);
        logger
    }
}

fn unsupported_dda() -> Error {
    Error::Unsupported(
        "BrukerTimsFile: DDA-PASEF export (SPECTRUM mode, or AUTO on DDA data) is not ported; \
         use export_mode FRAME"
            .into(),
    )
}

/// The source's `loadExperimentalSettings_`.
fn source_file(path: &Path) -> SourceFile {
    let text = path.to_string_lossy();
    SourceFile {
        name: crate::system::file::basename(&text).to_string(),
        path: crate::system::file::path(&text).to_string(),
        file_type: "Bruker TDF".into(),
        native_id_type: "Bruker TDF nativeID format".into(),
        native_id_type_accession: "MS:1002818".into(),
        ..SourceFile::default()
    }
}

fn open_db(path: &Path) -> Result<Connection> {
    Connection::open_with_flags(
        path.join("analysis.tdf"),
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )
    .map_err(sql_error)
}

/// The source's `openTimsDataHandle`: the linear handle, then the calibration
/// strategy.
fn open_handle(path: &Path, config: &Config) -> Result<TimsDataHandle> {
    let mut handle = TimsDataHandle::open(path, config.limits).map_err(|error| match error {
        Error::Io(io) => Error::Io(std::io::Error::new(
            io.kind(),
            format!("{} (opentims: {io})", path.display()),
        )),
        other => other,
    })?;
    let strategy = config.tims_calibration_strategy;
    if config.pressure_compensation != PressureCompensation::None
        && strategy != TimsCalibrationStrategy::BrukerSdk
        && strategy != TimsCalibrationStrategy::Auto
    {
        let _ =
            crate::openms_log_warn!("Pressure compensation requires BRUKER_SDK strategy, ignoring");
    }
    if matches!(
        strategy,
        TimsCalibrationStrategy::Auto | TimsCalibrationStrategy::BrukerSdk
    ) {
        let sdk_path = if config.bruker_sdk_path.is_empty() {
            std::env::var("OPENMS_BRUKER_SDK_PATH").unwrap_or_default()
        } else {
            config.bruker_sdk_path.clone()
        };
        if strategy == TimsCalibrationStrategy::BrukerSdk {
            return Err(Error::Unsupported(if sdk_path.is_empty() {
                "Bruker SDK path not set (use OPENMS_BRUKER_SDK_PATH or Config::bruker_sdk_path)"
                    .into()
            } else {
                "Bruker SDK failed: the closed-source Bruker SDK cannot be loaded by this port"
                    .into()
            }));
        }
        if !sdk_path.is_empty() {
            let _ = crate::openms_log_debug_nofile!(
                "Bruker SDK not available (not loadable by this port), trying rational"
            );
        }
    }
    if matches!(
        strategy,
        TimsCalibrationStrategy::Auto | TimsCalibrationStrategy::Rational
    ) {
        if handle.use_rational_calibration() {
            return Ok(handle);
        }
        if strategy == TimsCalibrationStrategy::Rational {
            return Err(Error::Io(std::io::Error::other(format!(
                "TimsCalibration table not found or unsupported in {}",
                path.display()
            ))));
        }
    }
    let _ = crate::openms_log_info!("TIMS calibration: linear (GlobalMetadata)");
    Ok(handle)
}

/// SQLite's integer conversion of a column (`getInt`): NULL and blobs 0,
/// reals truncated, text by its leading number; the low 32 bits.
pub(crate) fn sql_int(value: ValueRef<'_>) -> i32 {
    match value {
        ValueRef::Integer(i) => i as i32,
        ValueRef::Real(r) => (r as i64) as i32,
        ValueRef::Text(t) => {
            crate::format::tims_data_handle::c_atof(&String::from_utf8_lossy(t)) as i64 as i32
        }
        ValueRef::Null | ValueRef::Blob(_) => 0,
    }
}

/// SQLite's real conversion of a column (`getDouble`).
pub(crate) fn sql_double(value: ValueRef<'_>) -> f64 {
    match value {
        ValueRef::Integer(i) => i as f64,
        ValueRef::Real(r) => r,
        ValueRef::Text(t) => crate::format::tims_data_handle::c_atof(&String::from_utf8_lossy(t)),
        ValueRef::Null | ValueRef::Blob(_) => 0.0,
    }
}

/// The source's `resolveEffectiveConfig`: warn about empty frame ranges and
/// narrow the frame range to the frames whose `Time` lies in the RT range.
pub fn resolve_effective_config(
    db: &Connection,
    config: &Config,
    file_min: u32,
    file_max: u32,
) -> Config {
    let mut effective = config.clone();
    let frame_range_set = config.frame_id_min != 0 || config.frame_id_max != u32::MAX;
    if frame_range_set {
        if config.frame_id_min > config.frame_id_max {
            let _ = crate::openms_log_warn!(
                "Warning: BrukerTimsFile::Config::frame_id_min ({}) > frame_id_max ({}). No frames will be loaded.",
                config.frame_id_min,
                config.frame_id_max
            );
        } else if config.frame_id_max < file_min || config.frame_id_min > file_max {
            let _ = crate::openms_log_warn!(
                "Warning: BrukerTimsFile::Config frame range [{}, {}] does not intersect file frame range [{file_min}, {file_max}]. No frames will be loaded.",
                config.frame_id_min,
                config.frame_id_max
            );
        }
    }
    let rt_set = config.rt_min_sec > 0.0 || config.rt_max_sec.is_finite();
    if !rt_set {
        return effective;
    }
    if config.rt_min_sec > config.rt_max_sec {
        let _ = crate::openms_log_warn!(
            "Warning: BrukerTimsFile::Config::rt_min_sec ({}) > rt_max_sec ({}). No frames will be loaded.",
            config.rt_min_sec,
            config.rt_max_sec
        );
        effective.frame_id_min = 1;
        effective.frame_id_max = 0;
        return effective;
    }
    let upper = if config.rt_max_sec.is_finite() {
        config.rt_max_sec
    } else {
        f64::MAX
    };
    let bounds = db.query_row(
        "SELECT MIN(Id), MAX(Id) FROM Frames WHERE Time >= ? AND Time <= ?",
        rusqlite::params![config.rt_min_sec, upper],
        |row| {
            let min = row.get_ref(0)?;
            let max = row.get_ref(1)?;
            Ok(if matches!(min, ValueRef::Null) {
                None
            } else {
                Some((sql_int(min) as u32, sql_int(max) as u32))
            })
        },
    );
    match bounds {
        Ok(Some((rt_min_id, rt_max_id))) => {
            let (before_min, before_max) = (effective.frame_id_min, effective.frame_id_max);
            effective.frame_id_min = effective.frame_id_min.max(rt_min_id);
            effective.frame_id_max = effective.frame_id_max.min(rt_max_id);
            if before_min <= before_max && effective.frame_id_min > effective.frame_id_max {
                let _ = crate::openms_log_warn!(
                    "Warning: BrukerTimsFile::Config frame_id range [{before_min}, {before_max}] and RT range [{}, {}] s do not intersect. No frames will be loaded.",
                    config.rt_min_sec,
                    config.rt_max_sec
                );
            }
        }
        Ok(None) => {
            let _ = crate::openms_log_warn!(
                "Warning: BrukerTimsFile::Config RT range [{}, {}] s does not intersect any frame's acquisition time. No frames will be loaded.",
                config.rt_min_sec,
                config.rt_max_sec
            );
            effective.frame_id_min = 1;
            effective.frame_id_max = 0;
        }
        Err(error) => {
            let _ = crate::openms_log_warn!(
                "Warning: failed to resolve RT range via SQL: {error}. RT filter will be ignored."
            );
        }
    }
    effective
}

/// The source's `warnPartialMS1AggregationConfig`, including its one-shot
/// centroiding-configuration warnings.
fn warn_partial_ms1_aggregation_config(config: &Config, ms1_frames: usize) {
    if config.ms1_min_support > 0 && config.ms1_n_neighbors == 0 {
        let _ = crate::openms_log_warn!(
            "Warning: ms1_min_support (={}) ignored: ms1_n_neighbors=0 disables aggregation.",
            config.ms1_min_support
        );
    }
    if config.ms1_max_rt_distance_sec > 0.0 && config.ms1_n_neighbors == 0 {
        let _ = crate::openms_log_warn!(
            "Warning: ms1_max_rt_distance_sec (={}) ignored: ms1_n_neighbors=0 disables aggregation.",
            config.ms1_max_rt_distance_sec
        );
    }
    if config.ms1_n_neighbors > 0
        && ms1_frames > 0
        && (2 * i64::from(config.ms1_n_neighbors) + 1) as u64 > ms1_frames as u64
    {
        let _ = crate::openms_log_warn!(
            "Warning: ms1_n_neighbors (={}) exceeds half the run's MS1 frame count (={ms1_frames}); effective window will be clamped at run edges.",
            config.ms1_n_neighbors
        );
    }
    if !config.load_ms1 && config.ms1_n_neighbors > 0 {
        let _ = crate::openms_log_warn!(
            "Warning: ms1_n_neighbors (={}) ignored: load_ms1=false disables MS1 loading entirely.",
            config.ms1_n_neighbors
        );
    }
    let _ = config.effective_ms1_algo(true);
    let _ = config.effective_ms2_algo(true);
}

/// The source's `isDIA`: a non-empty `DiaFrameMsMsWindows`, else a non-empty
/// `DiaFrameMsMsInfo`. A missing table counts as empty.
pub fn is_dia(db: &Connection) -> bool {
    let non_empty = |table: &str| {
        db.query_row(&format!("SELECT COUNT(*) FROM {table}"), [], |row| {
            Ok(sql_int(row.get_ref(0)?))
        })
        .is_ok_and(|count| count > 0)
    };
    non_empty("DiaFrameMsMsWindows") || non_empty("DiaFrameMsMsInfo")
}

/// The source's `readDIAWindows`: `DiaFrameMsMsWindows` when non-empty, else
/// the distinct windows of `DiaFrameMsMsInfo` numbered 0, 1, ... in
/// `(IsolationMz, ScanNumBegin)` order. Scan bounds are converted to 1/K0 with
/// frame 1's calibration.
///
/// # Errors
///
/// The SQLite error as [`Error::Io`] when the chosen table cannot be read.
pub fn read_dia_windows(db: &Connection, handle: &TimsDataHandle) -> Result<Vec<DiaWindow>> {
    let converter = handle.scan2im().as_converter();
    let use_new = db
        .query_row("SELECT COUNT(*) FROM DiaFrameMsMsWindows", [], |row| {
            Ok(sql_int(row.get_ref(0)?))
        })
        .is_ok_and(|count| count > 0);
    let (sql, offset) = if use_new {
        (
            "SELECT WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth \
             FROM DiaFrameMsMsWindows ORDER BY WindowGroup, ScanNumBegin",
            1,
        )
    } else {
        (
            "SELECT DISTINCT ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth \
             FROM DiaFrameMsMsInfo ORDER BY IsolationMz, ScanNumBegin",
            0,
        )
    };
    let mut statement = db.prepare(sql).map_err(sql_error)?;
    let mut rows = statement.query([]).map_err(sql_error)?;
    let mut windows = Vec::new();
    let mut next_group = 0i32;
    while let Some(row) = rows.next().map_err(sql_error)? {
        let column = |index: usize| row.get_ref(index).map_err(sql_error);
        let window_group = if use_new {
            sql_int(column(0)?)
        } else {
            let group = next_group;
            next_group = next_group.wrapping_add(1);
            group
        };
        let scan_begin = sql_int(column(offset)?) as u32;
        let scan_end = sql_int(column(offset + 1)?) as u32;
        let ims = converter.convert_indices(1, &[scan_begin, scan_end]);
        windows.push(DiaWindow {
            window_group,
            mz_center: sql_double(column(offset + 2)?),
            mz_width: sql_double(column(offset + 3)?),
            im_lower: ims[0].min(ims[1]),
            im_upper: ims[0].max(ims[1]),
            scan_begin,
            scan_end,
        });
    }
    Ok(windows)
}

/// The source's `readFrameToWindowGroupMapping`: window group to its MS2
/// frame Ids in frame order.
///
/// # Errors
///
/// The SQLite error as [`Error::Io`].
pub fn read_frame_to_window_group(
    db: &Connection,
    windows: &[DiaWindow],
) -> Result<BTreeMap<i32, Vec<u32>>> {
    let mut mapping: BTreeMap<i32, Vec<u32>> = BTreeMap::new();
    let has_group_column = db
        .prepare("SELECT WindowGroup FROM DiaFrameMsMsInfo LIMIT 1")
        .and_then(|mut statement| statement.query([])?.next().map(|row| row.is_some()))
        .unwrap_or(false);
    if has_group_column {
        let mut statement = db
            .prepare("SELECT DISTINCT Frame, WindowGroup FROM DiaFrameMsMsInfo ORDER BY Frame")
            .map_err(sql_error)?;
        let mut rows = statement.query([]).map_err(sql_error)?;
        while let Some(row) = rows.next().map_err(sql_error)? {
            let frame = sql_int(row.get_ref(0).map_err(sql_error)?) as u32;
            let group = sql_int(row.get_ref(1).map_err(sql_error)?);
            mapping.entry(group).or_default().push(frame);
        }
        return Ok(mapping);
    }
    // Older layout: match each frame's window geometry to the synthetic groups.
    let mut statement = db
        .prepare(
            "SELECT Frame, IsolationMz, IsolationWidth, ScanNumBegin, ScanNumEnd \
             FROM DiaFrameMsMsInfo ORDER BY Frame",
        )
        .map_err(sql_error)?;
    let mut rows = statement.query([]).map_err(sql_error)?;
    while let Some(row) = rows.next().map_err(sql_error)? {
        let column = |index: usize| row.get_ref(index).map_err(sql_error);
        let frame = sql_int(column(0)?) as u32;
        let key = (
            sql_double(column(1)?),
            sql_double(column(2)?),
            sql_int(column(3)?) as u32,
            sql_int(column(4)?) as u32,
        );
        // The source keys a std::map on (mz, width, begin, end); with duplicate
        // keys its later window wins. The last matching window here is that one.
        let group = windows
            .iter()
            .rev()
            .find(|w| (w.mz_center, w.mz_width, w.scan_begin, w.scan_end) == key)
            .map(|w| w.window_group);
        if let Some(group) = group {
            let frames = mapping.entry(group).or_default();
            if frames.last() != Some(&frame) {
                frames.push(frame);
            }
        }
    }
    Ok(mapping)
}

fn in_range_groups(
    mut groups: BTreeMap<i32, Vec<u32>>,
    config: &Config,
) -> BTreeMap<i32, Vec<u32>> {
    for frames in groups.values_mut() {
        frames.retain(|&frame| config.frame_in_range(frame));
    }
    groups
}

/// The source's `transform` expected-size computation.
fn expected_spectra(
    db: &Connection,
    handle: &TimsDataHandle,
    config: &Config,
    is_dia: bool,
) -> Result<usize> {
    let count = |extra: &str| -> Result<usize> {
        db.query_row(
            &format!(
                "SELECT COUNT(*) FROM Frames WHERE Id >= {} AND Id <= {}{extra}",
                config.frame_id_min, config.frame_id_max
            ),
            [],
            |row| Ok(sql_int(row.get_ref(0)?)),
        )
        .map(|count| count as u32 as usize)
        .map_err(sql_error)
    };
    let ms1 = if config.load_ms1 {
        count(" AND MsMsType = 0")?
    } else {
        0
    };
    match config.export_mode {
        ExportMode::Frame => Ok(ms1 + count(" AND MsMsType != 0")?),
        ExportMode::Auto if is_dia => {
            let windows = read_dia_windows(db, handle)?;
            let groups = in_range_groups(read_frame_to_window_group(db, &windows)?, config);
            let mut per_group: BTreeMap<i32, usize> = BTreeMap::new();
            for window in &windows {
                *per_group.entry(window.window_group).or_default() += 1;
            }
            Ok(ms1
                + groups
                    .iter()
                    .filter_map(|(group, frames)| per_group.get(group).map(|n| frames.len() * n))
                    .sum::<usize>())
        }
        _ => Err(unsupported_dda()),
    }
}

/// The decoded peaks of one frame with m/z and 1/K0.
struct FrameBuffers {
    mz: Vec<f64>,
    intensity: Vec<f64>,
    im: Vec<f64>,
    scan_ids: Vec<u32>,
}

impl FrameBuffers {
    fn read(handle: &TimsDataHandle, frame: &TdfFrame) -> Result<Self> {
        let peaks = handle.read_frame(frame.id)?;
        Ok(Self {
            mz: handle.tof2mz().convert(&peaks.tofs),
            intensity: peaks.intensities.iter().map(|&i| f64::from(i)).collect(),
            im: handle
                .scan2im()
                .as_converter()
                .convert_indices(frame.id, &peaks.scan_ids),
            scan_ids: peaks.scan_ids,
        })
    }

    fn len(&self) -> usize {
        self.mz.len()
    }

    fn is_empty(&self) -> bool {
        self.mz.is_empty()
    }

    /// Keep the peaks whose `keep` entry is true, in order.
    fn retain(&mut self, keep: &[bool]) {
        fn apply<T: Copy>(values: &mut Vec<T>, keep: &[bool]) {
            *values = values
                .iter()
                .zip(keep)
                .filter(|(_, kept)| **kept)
                .map(|(value, _)| *value)
                .collect();
        }
        apply(&mut self.mz, keep);
        apply(&mut self.intensity, keep);
        apply(&mut self.im, keep);
        apply(&mut self.scan_ids, keep);
    }

    /// The source's `sortFrameBuffersByMz`, as a stable sort (the source's
    /// `std::sort` leaves equal m/z in unspecified order).
    fn sort_by_mz(&mut self) {
        let mut order: Vec<usize> = (0..self.len()).collect();
        order.sort_by(|&a, &b| self.mz[a].total_cmp(&self.mz[b]));
        self.mz = order.iter().map(|&i| self.mz[i]).collect();
        self.intensity = order.iter().map(|&i| self.intensity[i]).collect();
        self.im = order.iter().map(|&i| self.im[i]).collect();
        self.scan_ids = order.iter().map(|&i| self.scan_ids[i]).collect();
    }
}

/// The source's `isotopicPrefilter`: keep a peak only if some peak at `m/z ±
/// C13C12_MASSDIFF_U / q` (q = 1..5) within `tol_ppm` of that target lies at
/// most one scan away. Leaves the buffers sorted by m/z; returns the number
/// dropped.
fn isotopic_prefilter(buf: &mut FrameBuffers, tol_ppm: f64) -> usize {
    if buf.is_empty() {
        return 0;
    }
    buf.sort_by_mz();
    let n = buf.len();
    let mut keep = vec![false; n];
    for i in 0..n {
        if keep[i] {
            continue;
        }
        'charges: for q in 1..=5 {
            let delta = C13C12_MASSDIFF_U / f64::from(q);
            for sign in [-1.0, 1.0] {
                let target = buf.mz[i] + sign * delta;
                let tolerance = target * tol_ppm * 1e-6;
                let low = buf.mz.partition_point(|&mz| mz < target - tolerance);
                let high = buf.mz.partition_point(|&mz| mz <= target + tolerance);
                for j in low..high {
                    if j == i {
                        continue;
                    }
                    let scan_delta = i64::from(buf.scan_ids[j]) - i64::from(buf.scan_ids[i]);
                    if scan_delta.abs() <= 1 {
                        keep[i] = true;
                        keep[j] = true;
                        break 'charges;
                    }
                }
            }
        }
    }
    buf.retain(&keep);
    n - buf.len()
}

/// A spectrum with the metadata every frame spectrum carries.
fn frame_spectrum(frame: &TdfFrame, ms_level: u32, native_id: String) -> MSSpectrum {
    MSSpectrum {
        rt: frame.time,
        ms_level,
        drift_time_unit: DriftTimeUnit::InverseReducedMobility,
        spectrum_type: SpectrumType::Centroid,
        native_id,
        ..MSSpectrum::default()
    }
}

/// The source's `frameToSpectrum`: every peak of the frame, in storage order,
/// with its 1/K0 in the float data array. An empty frame gives a spectrum
/// without the array, as in the source.
fn frame_to_spectrum(
    handle: &TimsDataHandle,
    frame: &TdfFrame,
    ms_level: u32,
) -> Result<MSSpectrum> {
    let mut spectrum = frame_spectrum(frame, ms_level, format!("frame={}", frame.id));
    if frame.num_peaks == 0 {
        return Ok(spectrum);
    }
    let buf = FrameBuffers::read(handle, frame)?;
    spectrum.peaks = buf
        .mz
        .iter()
        .zip(&buf.intensity)
        .map(|(&mz, &intensity)| Peak1D {
            mz,
            intensity: intensity as f32,
        })
        .collect();
    spectrum.float_data_arrays.push(DataArray::new(
        IM_ARRAY_NAME,
        buf.im.iter().map(|&im| im as f32).collect(),
    ));
    Ok(spectrum)
}

/// The source's `loadMS1Spectrum`: [`frame_to_spectrum`], or the MS1
/// centroiding path (`centroidMS1Frame`) when MS1 centroiding is on.
fn ms1_spectrum(handle: &TimsDataHandle, frame: &TdfFrame, config: &Config) -> Result<MSSpectrum> {
    let algo = config.effective_ms1_algo(false);
    match algo {
        CentroidAlgo::Off => frame_to_spectrum(handle, frame, 1),
        CentroidAlgo::Greedy2d => Err(Error::Unsupported(
            "BrukerTimsFile: Greedy2D MS1 centroiding (FrameCentroider) is not ported".into(),
        )),
        CentroidAlgo::HillBased => {
            let mut spectrum = frame_spectrum(frame, 1, format!("frame={}", frame.id));
            if frame.num_peaks == 0 {
                return Ok(spectrum);
            }
            let mut buf = FrameBuffers::read(handle, frame)?;
            if config.isotopic_prefilter {
                let before = buf.len();
                let dropped = isotopic_prefilter(&mut buf, config.isotopic_prefilter_tol_ppm);
                log_prefilter("MS1", &format!("frame={}", frame.id), dropped, before);
                if buf.is_empty() {
                    return Ok(spectrum);
                }
            }
            let params = hill_params(
                config.ms1_centroid_mz_ppm,
                config.ms1_centroid_min_hill_length,
                config,
            );
            let centroids = centroid_frame(
                &buf.mz,
                &buf.intensity,
                &buf.im,
                &params,
                Some(&buf.scan_ids),
            )?;
            fill_centroided(&mut spectrum, &centroids, config.expose_hill_bounds);
            Ok(spectrum)
        }
    }
}

fn log_prefilter(stage: &str, context: &str, dropped: usize, before: usize) {
    let _ = crate::openms_log_info!(
        "Isotopic prefilter ({stage}, {context}): dropped {dropped} / {before} peaks ({}%)",
        100.0 * dropped as f64 / before.max(1) as f64
    );
}

fn hill_params(mz_ppm: f32, min_hill_length: usize, config: &Config) -> HillParams {
    HillParams {
        mz_tol_ppm: f64::from(mz_ppm),
        valley_factor: config.centroid_valley_factor,
        min_hill_length,
        max_scan_gap: config.centroid_max_scan_gap,
        ..HillParams::default()
    }
}

/// Centroids as peaks, the apex 1/K0 as the IM array, and optionally the
/// source's four hill-bounds arrays.
fn fill_centroided(spectrum: &mut MSSpectrum, centroids: &[Centroid], bounds: bool) {
    spectrum.peaks = centroids
        .iter()
        .map(|c| Peak1D {
            mz: c.mz,
            intensity: c.intensity as f32,
        })
        .collect();
    spectrum.float_data_arrays.push(DataArray::new(
        IM_ARRAY_NAME,
        centroids.iter().map(|c| c.im_apex as f32).collect(),
    ));
    if bounds {
        let array = |name: &str, value: fn(&Centroid) -> f64| {
            DataArray::new(name, centroids.iter().map(|c| value(c) as f32).collect())
        };
        spectrum
            .float_data_arrays
            .push(array("im lower bound", |c| c.im_lower));
        spectrum
            .float_data_arrays
            .push(array("im upper bound", |c| c.im_upper));
        spectrum
            .float_data_arrays
            .push(array("m/z lower bound", |c| c.mz_lower));
        spectrum
            .float_data_arrays
            .push(array("m/z upper bound", |c| c.mz_upper));
    }
}

/// One DIA MS2 spectrum of the source's raw `loadDIA_` path: the frame's
/// peaks whose 1/K0 lies inside the window, `None` when none (or, for
/// HillBased, no centroid) survives.
fn dia_ms2_spectrum(
    handle: &TimsDataHandle,
    frame: &TdfFrame,
    window: &DiaWindow,
    config: &Config,
    hill: bool,
) -> Result<Option<MSSpectrum>> {
    let all = FrameBuffers::read(handle, frame)?;
    let keep: Vec<bool> = all
        .im
        .iter()
        .map(|&im| im >= window.im_lower && im <= window.im_upper)
        .collect();
    let mut buf = all;
    buf.retain(&keep);
    if buf.is_empty() {
        return Ok(None);
    }
    if config.isotopic_prefilter {
        let before = buf.len();
        let dropped = isotopic_prefilter(&mut buf, config.isotopic_prefilter_tol_ppm);
        log_prefilter(
            "DIA-MS2 raw",
            &format!("frame={} window={}", frame.id, window.window_group),
            dropped,
            before,
        );
        if buf.is_empty() {
            return Ok(None);
        }
    }
    let mut spectrum = frame_spectrum(
        frame,
        2,
        format!(
            "frame={} windowGroup={} scan={}",
            frame.id, window.window_group, window.scan_begin
        ),
    );
    let mut precursor = Precursor::new(window.mz_center, 0);
    precursor.isolation_window_lower_offset = window.mz_width / 2.0;
    precursor.isolation_window_upper_offset = window.mz_width / 2.0;
    spectrum.precursors = vec![precursor];
    spectrum.metadata.insert(
        "ion mobility lower limit".into(),
        MetaValue::try_from(window.im_lower)?,
    );
    spectrum.metadata.insert(
        "ion mobility upper limit".into(),
        MetaValue::try_from(window.im_upper)?,
    );
    if hill {
        let params = hill_params(
            config.ms2_centroid_mz_ppm,
            config.ms2_centroid_min_hill_length,
            config,
        );
        let centroids = centroid_frame(
            &buf.mz,
            &buf.intensity,
            &buf.im,
            &params,
            Some(&buf.scan_ids),
        )?;
        if centroids.is_empty() {
            return Ok(None);
        }
        fill_centroided(&mut spectrum, &centroids, config.expose_hill_bounds);
    } else {
        spectrum.peaks = buf
            .mz
            .iter()
            .zip(&buf.intensity)
            .map(|(&mz, &intensity)| Peak1D {
                mz,
                intensity: intensity as f32,
            })
            .collect();
        spectrum.float_data_arrays.push(DataArray::new(
            IM_ARRAY_NAME,
            buf.im.iter().map(|&im| im as f32).collect(),
        ));
    }
    Ok(Some(spectrum))
}

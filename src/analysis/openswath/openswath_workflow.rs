// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The OpenSWATH workflow: extraction and scoring of an assay library in
//! SWATH (DIA), PRM, diaPASEF and SRM/MRM data.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathWorkflow.h` with its `.cpp`.
//! [`OpenSwathWorkflow::perform_extraction`] runs the complete analysis on a
//! vector of [`SwathMap`]s:
//!
//! 1. obtain the precursor (MS1) chromatograms when enabled
//!    ([`OpenSwathWorkflowBase::ms1_extraction`]) and score them alone for
//!    MS1-only data;
//! 2. select the transitions of every SWATH window with
//!    [`select_swath_transitions`](crate::analysis::openswath::openswath_helper::select_swath_transitions)
//!    (or the best-matching window for PRM and diaPASEF data);
//! 3. extract the chromatograms of each batch of compounds
//!    ([`prepare_extraction_coordinates`](OpenSwathWorkflowBase::prepare_extraction_coordinates),
//!    [`ChromatogramExtractor`], [`return_chromatogram`]);
//! 4. score them ([`OpenSwathWorkflow::score_all_chromatograms`]) with the
//!    transition group picker and [`MRMFeatureFinderScoring`];
//! 5. write the chromatograms and the features to the consumer, the feature
//!    map and, as rows, the OSW writer.
//!
//! # Execution order and parallelism
//!
//! The source schedules this with OpenMP in three ways: a wave scheduler for
//! in-memory data (`batch_size <= 0`), per-SWATH batches otherwise, and a
//! nested-parallelism variant. A single-threaded run of each has one result,
//! and this port reproduces it: the SWATH windows of a wave are extracted in
//! parallel, the scoring jobs (sorted by their transition count, descending)
//! are scored in parallel, and their results are *consumed in job order* on
//! one thread. The output of the per-SWATH path is therefore bit-identical for
//! any thread count (set with [`OpenSwathWorkflow::set_threads`]). The wave
//! scheduler sorts the features of a wave by position when it scored several
//! jobs, as the source does, and which SWATHs share a wave follows the thread
//! count (memory and `scoring_threads` cap it), so there the *set* of features
//! with every value is what does not depend on the thread count, not their
//! order. A mobilogram consumer, which cannot be shared, serializes the
//! scoring.
//!
//! The unique identifiers of the features are keyed by the SWATH and the
//! assay they belong to (a generator seeded from the workflow's generator,
//! once per call), never by the schedule: seed the workflow generator
//! ([`OpenSwathWorkflow::set_unique_id_generator`]) for reproducible
//! identifiers.
//!
//! # Not ported
//!
//! - the nested OpenMP variant (`threads_outer_loop >= 0` in a build with
//!   `MT_ENABLE_NESTED_OPENMP`): the source falls back to the wave scheduler
//!   without it, and so does this port;
//! - the source's `OSWBufferedWriter` thread: rows are written, in job order,
//!   as soon as a chunk of jobs is scored; the memory bound is that chunk;
//! - the informational log lines of the source are kept in
//!   [`OpenSwathWorkflow::take_log`] instead of a global log stream;
//! - `ScopedResamplingWarningSuppression`, whose warnings this port's
//!   resampler does not emit;
//! - SRM/MRM mode needs a [`ChromatogramHandler`]
//!   (`ANALYSIS/TARGETED/IChromatogramHandler.h` with `DefaultChromHandler.h`
//!   and `MRMChromHandler.h`); without one it returns
//!   [`Error::Unsupported`](crate::Error::Unsupported).
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_WORKFLOW_SUPPORT.md).
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::openswath_workflow::{ChromExtractParams, OpenSwathWorkflow};
//!
//! let workflow = OpenSwathWorkflow::new(false, false, false, false, false, -1);
//! assert!(!workflow.mrm());
//! let params = ChromExtractParams::default();
//! assert_eq!(params.rt_extraction_window, 0.0);
//! ```

use super::chromatogram_extractor::{
    ChromatogramExtractor, prepare_coordinates, return_chromatogram,
};
use super::chromatogram_extractor_algorithm::ExtractionCoordinates;
use super::data_structures::Chromatogram;
use super::ion_mobility_scoring::MobilogramSink;
use super::mrm_feature_finder_scoring::MRMFeatureFinderScoring;
use super::mrm_transition_group_picker::MRMTransitionGroupPicker;
use super::openswath_helper::{
    compute_precursor_id, select_swath_transitions, select_swath_transitions_pasef,
};
use super::openswath_osw_writer::{OpenSwathOSWWriter, OswData};
use super::openswath_workflow_scheduler::{
    ConcurrencyEstimate, ConcurrencyLimiter, OpenSwathWorkflowScheduler, SchedulerOptions,
    ScopedSlot, free_memory_kib,
};
use super::spectrum_access::{SpectrumAccess, SpectrumAccessPtr};
use super::spectrum_access_openms::SpectrumAccessOpenMSInMemory;
use super::swath_map::SwathMap;
use super::transition_experiment::{LightCompound, LightTargetedExperiment, LightTransition};
use crate::analysis::transformations::TransformationDescription;
use crate::concept::parallel::{Threads, map_collect};
use crate::concept::progress_logger::ProgressLogger;
use crate::concept::unique_id::UniqueIdGenerator;
use crate::interfaces::MSDataConsumer;
use crate::kernel::mrm::MRMTransitionGroup;
use crate::kernel::{FeatureMap, MSChromatogram, MSSpectrum};
use crate::math::source_sort::source_sort_by;
use crate::param::Param;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::ops::ControlFlow;
use std::sync::{Arc, Mutex};

/// The parameters of the chromatogram extraction (source `ChromExtractParams`).
///
/// The source leaves every member uninitialized; the default here is all zero
/// and an empty extraction function.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChromExtractParams {
    /// Do not extract anything closer than this (in Da) to the upper edge of
    /// the isolation window.
    pub min_upper_edge_dist: f64,
    /// The m/z extraction window in Da or ppm (50 ppm extracts +/- 25 ppm).
    pub mz_extraction_window: f64,
    /// The ion mobility extraction window.
    pub im_extraction_window: f64,
    /// Whether `mz_extraction_window` is in ppm instead of Da.
    pub ppm: bool,
    /// The extraction function in mass space (`tophat`).
    pub extraction_function: String,
    /// The retention time extraction window.
    pub rt_extraction_window: f64,
    /// Extra retention time to extract (to look at the chromatogram outside
    /// the window).
    pub extra_rt_extract: f64,
}

/// A provider of chromatograms for the workflow's SRM/MRM mode (source
/// `IChromatogramHandler`, `ANALYSIS/TARGETED/IChromatogramHandler.h`).
///
/// The interface isolates the workflow from the source and mapping strategy
/// of the chromatograms. It is defined here, next to its only caller, until
/// the `ANALYSIS/TARGETED` handlers are ported; the factory
/// `IChromatogramHandler::createDefault` is part of that port.
pub trait ChromatogramHandler: Send + Sync {
    /// Collect iRT chromatograms from the SWATH maps and try to map them to
    /// the iRT transitions (source `collectIrtChromatogramsForIrt`).
    ///
    /// # Errors
    ///
    /// Implementation-specific.
    #[allow(clippy::too_many_arguments)]
    fn collect_irt_chromatograms_for_irt(
        &self,
        swath_maps: &[SwathMap],
        irt_transitions: &LightTargetedExperiment,
        mrm_mapping_param: &Param,
        cp: &ChromExtractParams,
        trafo: &TransformationDescription,
        pasef: bool,
        load_into_memory: bool,
    ) -> Result<Vec<MSChromatogram>>;

    /// Extract or select chromatograms for the transitions and return the
    /// mapped and filtered ones (source
    /// `extractAndMapChromatogramsForTransitions`).
    ///
    /// # Errors
    ///
    /// Implementation-specific.
    fn extract_and_map_chromatograms_for_transitions(
        &self,
        swath_maps: &[SwathMap],
        transition_exp: &LightTargetedExperiment,
        cp: &ChromExtractParams,
        mrm_mapping_param: &Param,
    ) -> Result<Vec<MSChromatogram>>;
}

/// The severity of an entry of the workflow log.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LogLevel {
    /// `OPENMS_LOG_INFO`.
    Info,
    /// `OPENMS_LOG_WARN`.
    Warn,
    /// `OPENMS_LOG_ERROR`.
    Error,
}

/// One message the workflow would have written to the OpenMS log.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LogEntry {
    /// The severity.
    pub level: LogLevel,
    /// The text, as the source writes it.
    pub message: String,
}

/// The state shared by the OpenSWATH workflows (source `OpenSwathWorkflowBase`,
/// the base class of [`OpenSwathWorkflow`] and `OpenSwathCalibrationWorkflow`).
///
/// The source members are `protected`; here they are public so that a
/// calibration workflow can hold a `OpenSwathWorkflowBase` by composition.
pub struct OpenSwathWorkflowBase {
    /// Access to the MS1 map (source `ms1_map_`); `None` unless set by
    /// [`OpenSwathWorkflow::load_ms1_map`] use. Whoever uses it concurrently
    /// calls [`SpectrumAccess::light_clone`] first.
    pub ms1_map: Option<SpectrumAccessPtr>,
    /// Whether to use the MS1 traces.
    pub use_ms1_traces: bool,
    /// Whether to use ion mobility extraction on MS1 traces.
    pub use_ms1_ion_mobility: bool,
    /// Whether the data is targeted DIA (PRM) with potentially overlapping
    /// windows: a precursor is extracted from the single matching window whose
    /// m/z centre is closest to its library m/z.
    pub prm: bool,
    /// Whether the data is diaPASEF: a precursor is extracted from the single
    /// window matching in m/z and ion mobility whose ion mobility centre is
    /// closest.
    pub pasef: bool,
    /// Whether the data is chromatogram-only SRM/MRM data.
    pub mrm: bool,
    /// The legacy nested outer-loop thread count; -1 leaves scheduling
    /// automatic. Not used by this port, see the module documentation.
    pub threads_outer_loop: i32,
    /// The progress reporter (the source class is a `ProgressLogger`), behind a
    /// lock because the parallel stages report through a shared reference.
    pub logger: Mutex<ProgressLogger>,
}

impl std::fmt::Debug for OpenSwathWorkflowBase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OpenSwathWorkflowBase")
            .field("ms1_map", &self.ms1_map.as_ref().map(|_| "SpectrumAccess"))
            .field("use_ms1_traces", &self.use_ms1_traces)
            .field("use_ms1_ion_mobility", &self.use_ms1_ion_mobility)
            .field("prm", &self.prm)
            .field("pasef", &self.pasef)
            .field("mrm", &self.mrm)
            .field("threads_outer_loop", &self.threads_outer_loop)
            .finish()
    }
}

impl Default for OpenSwathWorkflowBase {
    /// The source default constructor: no MS1 traces, no PRM, no diaPASEF, no
    /// MRM, all threads in the outer loop.
    fn default() -> Self {
        Self::new(false, false, false, false, false, -1)
    }
}

impl OpenSwathWorkflowBase {
    /// Construct a workflow base (source constructor).
    ///
    /// `threads_outer_loop >= 0` requests the legacy nested per-SWATH
    /// parallelism of the source, `-1` automatic scheduling.
    pub fn new(
        use_ms1_traces: bool,
        use_ms1_ion_mobility: bool,
        prm: bool,
        pasef: bool,
        mrm: bool,
        threads_outer_loop: i32,
    ) -> Self {
        Self {
            ms1_map: None,
            use_ms1_traces,
            use_ms1_ion_mobility,
            prm,
            pasef,
            mrm,
            threads_outer_loop,
            logger: Mutex::new(ProgressLogger::new()),
        }
    }

    /// Start a progress section (source `startProgress`).
    ///
    /// # Errors
    ///
    /// The errors of the progress backend.
    pub fn start_progress(&self, begin: i64, end: i64, label: &str) -> Result<()> {
        lock(&self.logger).start_progress(begin, end, label)
    }

    /// Report progress (source `setProgress`).
    ///
    /// # Errors
    ///
    /// The errors of the progress backend.
    pub fn set_progress(&self, value: i64) -> Result<()> {
        lock(&self.logger).set_progress(value)
    }

    /// End the progress section (source `endProgress`).
    ///
    /// # Errors
    ///
    /// The errors of the progress backend.
    pub fn end_progress(&self, bytes_processed: u64) -> Result<()> {
        lock(&self.logger).end_progress(bytes_processed)
    }

    /// Prepare the extraction coordinates and empty chromatograms, handling
    /// the retention time transformation (source `prepareExtractionCoordinates_`).
    ///
    /// Creates one empty chromatogram and one coordinate with the right id,
    /// m/z and RT start and end per transition (or per precursor isotope with
    /// `ms1`). With a negative `cp.rt_extraction_window` the extraction is
    /// not restricted in RT. Otherwise the RT of the assay is transformed
    /// back to the run's RT space with `trafo_inverse` and the window
    /// `(rt_extraction_window + extra_rt_extract) / 2` is added on either side.
    ///
    /// # Errors
    ///
    /// The errors of
    /// [`prepare_coordinates`](super::chromatogram_extractor::prepare_coordinates)
    /// and of the transformation.
    pub fn prepare_extraction_coordinates(
        &self,
        transition_exp_used: &LightTargetedExperiment,
        trafo_inverse: &TransformationDescription,
        cp: &ChromExtractParams,
        ms1: bool,
        ms1_isotopes: i32,
    ) -> Result<(Vec<Chromatogram>, Vec<ExtractionCoordinates>)> {
        let mut chrom_list = Vec::new();
        let mut coordinates = Vec::new();
        if cp.rt_extraction_window < 0.0 {
            prepare_coordinates(
                &mut chrom_list,
                &mut coordinates,
                transition_exp_used,
                cp.rt_extraction_window,
                ms1,
                ms1_isotopes,
            )?;
        } else {
            // an RT window of 0 just writes the RT into start and end; those
            // are then transformed and widened by the extra window
            prepare_coordinates(
                &mut chrom_list,
                &mut coordinates,
                transition_exp_used,
                0.0,
                ms1,
                ms1_isotopes,
            )?;
            let half = (cp.rt_extraction_window + cp.extra_rt_extract) / 2.0;
            for coordinate in &mut coordinates {
                coordinate.rt_start = trafo_inverse.apply(coordinate.rt_start)? - half;
                coordinate.rt_end = trafo_inverse.apply(coordinate.rt_end)? + half;
            }
        }
        Ok((chrom_list, coordinates))
    }

    /// Extract the precursor (MS1) chromatograms of the assays (source
    /// `MS1Extraction_`); one per compound and isotope, named
    /// `<group id>_Precursor_i<n>`.
    ///
    /// `swath_maps` and `ms1_only` are unused by the source as well.
    ///
    /// # Errors
    ///
    /// The errors of the coordinate preparation, the extraction and the
    /// conversion back to chromatograms.
    #[allow(clippy::too_many_arguments)]
    pub fn ms1_extraction(
        &self,
        ms1_map: &dyn SpectrumAccess,
        _swath_maps: &[SwathMap],
        cp: &ChromExtractParams,
        transition_exp: &LightTargetedExperiment,
        trafo_inverse: &TransformationDescription,
        _ms1_only: bool,
        ms1_isotopes: i32,
    ) -> Result<Vec<MSChromatogram>> {
        let (mut chrom_list, coordinates) = self.prepare_extraction_coordinates(
            transition_exp,
            trafo_inverse,
            cp,
            true,
            ms1_isotopes,
        )?;
        ChromatogramExtractor::new().extract_chromatograms_im(
            ms1_map,
            &mut chrom_list,
            &coordinates,
            cp.mz_extraction_window,
            cp.ppm,
            cp.im_extraction_window,
            &cp.extraction_function,
        )?;
        let mut out = Vec::new();
        return_chromatogram(
            &chrom_list,
            &coordinates,
            transition_exp,
            &MSSpectrum::default(),
            &mut out,
            true,
            cp.im_extraction_window,
        )?;
        Ok(out)
    }
}

// --- helpers --------------------------------------------------------------

const OSW_MIN_WRITE_BUFFER_BYTES: u64 = 64 * 1024 * 1024;
const OSW_FALLBACK_WRITE_BUFFER_BYTES: u64 = 512 * 1024 * 1024;
const OSW_MAX_WRITE_BUFFER_BYTES: u64 = 2 * 1024 * 1024 * 1024;
const OSW_SYSTEM_MEMORY_RESERVE_BYTES: u64 = 4 * 1024 * 1024 * 1024;

/// The OSW write buffer budget (source `determineOSWWriteBufferBytes`): a
/// quarter of the free memory beyond a 4 GiB reserve (an eighth of the free
/// memory when that is not above the reserve), clamped to 64 MiB .. 2 GiB;
/// 512 MiB when the free memory is unknown.
pub fn determine_osw_write_buffer_bytes() -> u64 {
    osw_write_buffer_bytes_from(free_memory_kib())
}

/// [`determine_osw_write_buffer_bytes`] for a given amount of free memory in
/// KiB (`None`: not reported).
pub fn osw_write_buffer_bytes_from(free_kib: Option<u64>) -> u64 {
    let Some(kib) = free_kib else {
        return OSW_FALLBACK_WRITE_BUFFER_BYTES;
    };
    let available = kib.saturating_mul(1024);
    let usable = if available > OSW_SYSTEM_MEMORY_RESERVE_BYTES {
        available - OSW_SYSTEM_MEMORY_RESERVE_BYTES
    } else {
        available / 8
    };
    (usable / 4).clamp(OSW_MIN_WRITE_BUFFER_BYTES, OSW_MAX_WRITE_BUFFER_BYTES)
}

/// The inner scoring batch size of the per-SWATH path (source
/// `calculateInnerBatchSize`): a positive user value capped at the number of
/// compounds; otherwise 5% of the memory available for scoring over 2 KiB
/// per compound, clamped to 2000 .. 10000 and capped at the number of
/// compounds (5000 when the memory is unknown).
pub fn calculate_inner_batch_size(total_compounds: usize, user_inner_batch_size: i32) -> usize {
    calculate_inner_batch_size_from(
        total_compounds,
        user_inner_batch_size,
        OpenSwathWorkflowScheduler::estimate_available_memory_for_scoring(
            &SchedulerOptions::default(),
        ),
    )
}

/// [`calculate_inner_batch_size`] for a given memory budget in bytes (0:
/// unknown).
pub fn calculate_inner_batch_size_from(
    total_compounds: usize,
    user_inner_batch_size: i32,
    safe_memory: u64,
) -> usize {
    if user_inner_batch_size > 0 {
        return (user_inner_batch_size as usize).min(total_compounds);
    }
    if safe_memory == 0 {
        return 5000.min(total_compounds);
    }
    let feature_buffer_safe = safe_memory / 20;
    let bytes_per_compound = 2 * 1024;
    let inner = (feature_buffer_safe / bytes_per_compound).min(usize::MAX as u64) as usize;
    inner.clamp(2000, 10000).min(total_compounds)
}

/// One feature of a job, its chromatograms and its OSW rows.
struct ScoredJob {
    features: FeatureMap,
    osw: OswData,
}

/// One extracted and scored batch of the per-SWATH path.
struct ScoredBatch {
    chromatograms: Vec<MSChromatogram>,
    ms1_chromatograms: Vec<MSChromatogram>,
    features: FeatureMap,
    osw: OswData,
}

/// Per-SWATH state of the wave scheduler (source `SwathSchedulerContext`).
#[derive(Default)]
struct SwathContext {
    active: bool,
    score_job_size: usize,
    nr_score_jobs: usize,
    remaining_score_jobs: usize,
    finalized: bool,
    transition_exp_used_all: LightTargetedExperiment,
    current_swath_map: Option<SpectrumAccessPtr>,
    ms1_chromatograms: Vec<MSChromatogram>,
    chrom_exp: Vec<MSChromatogram>,
    feature_file: FeatureMap,
}

#[derive(Clone, Debug)]
struct InnerBatchScoringJob {
    context_index: usize,
    batch_index: usize,
    estimated_compounds: usize,
    estimated_transitions: usize,
}

/// The OpenSWATH workflow (source `OpenSwathWorkflow`).
pub struct OpenSwathWorkflow {
    base: OpenSwathWorkflowBase,
    threads: Threads,
    ids: Mutex<UniqueIdGenerator>,
    chrom_handler: Option<Arc<dyn ChromatogramHandler>>,
    log: Mutex<Vec<LogEntry>>,
}

impl std::fmt::Debug for OpenSwathWorkflow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OpenSwathWorkflow")
            .field("base", &self.base)
            .field("threads", &self.threads)
            .field("chrom_handler", &self.chrom_handler.is_some())
            .finish()
    }
}

/// A well-mixed combination of two 64-bit keys (SplitMix64 finalizer).
fn mix_seed(a: u64, b: u64) -> u64 {
    let mut z = a ^ b
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add(0x1234_5678_9ABC_DEF1);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// The generator key of one assay of one swath: `seed` combined with the FNV-1a
/// hash of the assay identifier.
fn assay_seed(seed: u64, assay: &str) -> u64 {
    let mut hash: u64 = 0xCBF2_9CE4_8422_2325;
    for byte in assay.as_bytes() {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01B3);
    }
    mix_seed(seed, hash)
}

fn lock<T>(mutex: &Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    mutex.lock().unwrap_or_else(|e| e.into_inner())
}

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// The in-memory copy of a swath map access (source
/// `SpectrumAccessOpenMSInMemory(*sptr)`).
fn in_memory(sptr: &SpectrumAccessPtr) -> Result<SpectrumAccessPtr> {
    Ok(Arc::new(SpectrumAccessOpenMSInMemory::new(&**sptr)?))
}

/// Select the compounds of batch `batch_idx` (source `selectCompoundsForBatch_`):
/// all proteins, the compounds `batch_idx * batch_size ..` (at most
/// `batch_size`) and the transitions of those compounds.
pub fn select_compounds_for_batch(
    transition_exp_used_all: &LightTargetedExperiment,
    batch_size: usize,
    batch_idx: usize,
) -> LightTargetedExperiment {
    let total = transition_exp_used_all.compounds.len();
    let start = batch_idx.saturating_mul(batch_size).min(total);
    let end = batch_idx
        .saturating_mul(batch_size)
        .saturating_add(batch_size)
        .min(total);
    let mut used = LightTargetedExperiment::new();
    used.proteins = transition_exp_used_all.proteins.clone();
    used.compounds = transition_exp_used_all.compounds[start..end].to_vec();
    used.transitions =
        copy_batch_transitions(&used.compounds, &transition_exp_used_all.transitions);
    used
}

/// The transitions of `all_transitions` that reference one of the compounds
/// (source `copyBatchTransitions_`), in library order.
pub fn copy_batch_transitions(
    used_compounds: &[LightCompound],
    all_transitions: &[LightTransition],
) -> Vec<LightTransition> {
    let selected: HashSet<&str> = used_compounds.iter().map(|c| c.id.as_str()).collect();
    all_transitions
        .iter()
        .filter(|t| selected.contains(t.peptide_ref.as_str()))
        .cloned()
        .collect()
}

/// The window of every transition for PRM data (source `prm_` branch): of
/// the windows that contain the precursor (strictly, and not closer than
/// `min_upper_edge_dist` to the upper edge) the one whose m/z centre is
/// closest; the earlier window on a tie.
fn select_prm_windows(
    transition_exp: &LightTargetedExperiment,
    min_upper_edge_dist: f64,
    swath_maps: &[SwathMap],
) -> Vec<Option<usize>> {
    let mut map: Vec<Option<usize>> = vec![None; transition_exp.transitions.len()];
    for (i, swath) in swath_maps.iter().enumerate() {
        for (k, tr) in transition_exp.transitions.iter().enumerate() {
            let mz = tr.precursor_mz;
            if swath.lower < mz
                && mz < swath.upper
                && (swath.upper - mz).abs() >= min_upper_edge_dist
            {
                if map[k].is_none() {
                    map[k] = Some(i);
                }
                if let Some(current) = map[k] {
                    if (swath_maps[current].center - mz).abs() > (swath.center - mz).abs() {
                        map[k] = Some(i);
                    }
                }
            }
        }
    }
    map
}

/// The transitions, compounds and proteins whose transition was assigned to
/// window `index` by a PRM or diaPASEF window map.
fn select_assigned_transitions(
    transition_exp: &LightTargetedExperiment,
    tr_win_map: &[Option<usize>],
    index: usize,
) -> LightTargetedExperiment {
    let mut used = LightTargetedExperiment::new();
    let mut matching_compounds: BTreeSet<&str> = BTreeSet::new();
    for (k, window) in tr_win_map.iter().enumerate() {
        if *window == Some(index) {
            let tr = &transition_exp.transitions[k];
            used.transitions.push(tr.clone());
            matching_compounds.insert(tr.peptide_ref.as_str());
        }
    }
    let mut matching_proteins: BTreeSet<&str> = BTreeSet::new();
    for compound in &transition_exp.compounds {
        if matching_compounds.contains(compound.id.as_str()) {
            used.compounds.push(compound.clone());
            for reference in &compound.protein_refs {
                matching_proteins.insert(reference.as_str());
            }
        }
    }
    for protein in &transition_exp.proteins {
        if matching_proteins.contains(protein.id.as_str()) {
            used.proteins.push(protein.clone());
        }
    }
    used
}

/// Reborrow an optional mobilogram sink for one call, with the object lifetime
/// shortened to the borrow (the callee's signature ties the two together).
fn reborrow_sink<'s>(
    sink: &'s mut Option<&mut (dyn MobilogramSink + '_)>,
) -> Option<&'s mut dyn MobilogramSink> {
    match sink {
        Some(s) => Some(&mut **s),
        None => None,
    }
}

fn sort_features_by_position(map: &mut FeatureMap) -> Result<()> {
    // `Peak2D::PositionLess` on the feature map, ordered like the Release
    // build's `std::sort`
    source_sort_by(&mut map.features, |a, b| {
        if a.rt < b.rt {
            return true;
        }
        if a.rt > b.rt {
            return false;
        }
        a.mz < b.mz
    })
}

impl OpenSwathWorkflow {
    /// Construct the workflow (source constructor).
    ///
    /// `threads_outer_loop` is kept for the nested variant of the source,
    /// which this port does not implement; use
    /// [`set_threads`](Self::set_threads) to choose the thread count.
    pub fn new(
        use_ms1_traces: bool,
        use_ms1_ion_mobility: bool,
        prm: bool,
        pasef: bool,
        mrm: bool,
        threads_outer_loop: i32,
    ) -> Self {
        Self {
            base: OpenSwathWorkflowBase::new(
                use_ms1_traces,
                use_ms1_ion_mobility,
                prm,
                pasef,
                mrm,
                threads_outer_loop,
            ),
            threads: Threads::all(),
            ids: Mutex::new(UniqueIdGenerator::new()),
            chrom_handler: None,
            log: Mutex::new(Vec::new()),
        }
    }

    /// The shared workflow state.
    pub fn base(&self) -> &OpenSwathWorkflowBase {
        &self.base
    }

    /// The shared workflow state, mutable (progress reporting, MS1 map).
    pub fn base_mut(&mut self) -> &mut OpenSwathWorkflowBase {
        &mut self.base
    }

    /// Whether the data is SRM/MRM data.
    pub fn mrm(&self) -> bool {
        self.base.mrm
    }

    /// Set the number of threads of the extraction and scoring stages (the
    /// source uses OpenMP's maximum). The result does not depend on it.
    pub fn set_threads(&mut self, threads: Threads) {
        self.threads = threads;
    }

    /// Replace the generator of the unique feature ids (the source draws from
    /// the global `UniqueIdGenerator`).
    pub fn set_unique_id_generator(&mut self, generator: UniqueIdGenerator) {
        self.ids = Mutex::new(generator);
    }

    /// Provide the chromatogram handler of the SRM/MRM mode (source
    /// `IChromatogramHandler::createDefault`).
    pub fn set_chromatogram_handler(&mut self, handler: Option<Arc<dyn ChromatogramHandler>>) {
        self.chrom_handler = handler;
    }

    /// The messages the source writes to the OpenMS log, in order, emptying
    /// the buffer.
    pub fn take_log(&self) -> Vec<LogEntry> {
        std::mem::take(&mut *lock(&self.log))
    }

    fn log(&self, level: LogLevel, message: impl Into<String>) {
        lock(&self.log).push(LogEntry {
            level,
            message: message.into(),
        });
    }

    fn draw_seed(&self) -> u64 {
        lock(&self.ids).get_unique_id()
    }

    /// The first MS1 map of `swath_maps` (source `loadMS1Map`), or `None`.
    ///
    /// With `load_into_memory` the map is copied into a
    /// [`SpectrumAccessOpenMSInMemory`]. The source keeps the *last* map
    /// flagged MS1 when there are several; so does this.
    ///
    /// # Errors
    ///
    /// The errors of reading the spectra into memory.
    pub fn load_ms1_map(
        &self,
        swath_maps: &[SwathMap],
        load_into_memory: bool,
    ) -> Result<Option<SpectrumAccessPtr>> {
        let mut ms1_map: Option<SpectrumAccessPtr> = None;
        for map in swath_maps {
            if map.ms1 {
                ms1_map = map.sptr.clone();
            }
        }
        if load_into_memory {
            if let Some(map) = &ms1_map {
                return Ok(Some(in_memory(map)?));
            }
        }
        Ok(ms1_map)
    }

    /// Execute the OpenSWATH analysis on a set of SWATH maps and an assay
    /// library (source `performExtraction`).
    ///
    /// # Arguments
    ///
    /// - `swath_maps`: the raw data; an MS1 map has `ms1` set.
    /// - `rt_trafo`: translates the RT of this run to the normalized RT space
    ///   of the library.
    /// - `chromatogram_extraction_params`, `ms1_chromatogram_extraction_params`:
    ///   the extraction parameters of the fragment and precursor traces.
    /// - `feature_finder_param`: the parameters of the feature finding in
    ///   the chromatographic dimension ([`MRMFeatureFinderScoring`] with its
    ///   `TransitionGroupPicker:` section).
    /// - `assay_library`: the assays to extract and score.
    /// - `result_feature_file`, `store_features_in_feature_file`: the features
    ///   are appended to the map when the flag is set (`result_osw` then
    ///   normally is inactive; the two outputs are exclusive in the source).
    /// - `result_osw`: receives the features as rows when it is active.
    /// - `result_chromatograms`: receives the extracted chromatograms.
    /// - `batch_size`: values `<= 0` select the automatic scheduling, which
    ///   with `load_into_memory` is the SWATH wave scheduler; a positive
    ///   value forces the per-SWATH path with that many compounds per batch.
    /// - `ms1_isotopes`: the number of MS1 isotopes to extract besides the
    ///   monoisotopic peak.
    /// - `load_into_memory`: cache each SWATH map in memory while it is used.
    /// - `mrm_mapping_param`: the parameters of the chromatogram mapping in
    ///   SRM/MRM mode.
    /// - `mobilogram_consumer`: receives the extracted ion mobilograms.
    /// - `inner_batch_size`: the scoring batch size of the wave scheduler;
    ///   `<= 0` sizes it from the thread count and the number of SWATHs.
    /// - `max_concurrent_swaths`: the SWATHs kept resident by the wave
    ///   scheduler; `<= 0` plans by memory.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when MS1-only data is given without
    /// `use_ms1_traces`, when a SWATH map has no data, when a transition has
    /// no chromatogram (not in SRM/MRM mode, where it is skipped), when an
    /// assay has no detecting chromatogram or transition (same), and the
    /// errors of picking and scoring in DIA mode; [`Error::Unsupported`]
    /// when `prm` and `pasef` are both set (source `NotImplemented`), and in
    /// SRM/MRM mode without a [`ChromatogramHandler`]; the errors of the
    /// extraction, the consumers and the OSW writer.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_extraction(
        &mut self,
        swath_maps: &[SwathMap],
        rt_trafo: &TransformationDescription,
        chromatogram_extraction_params: &ChromExtractParams,
        ms1_chromatogram_extraction_params: &ChromExtractParams,
        feature_finder_param: &Param,
        assay_library: &LightTargetedExperiment,
        result_feature_file: &mut FeatureMap,
        store_features_in_feature_file: bool,
        result_osw: &OpenSwathOSWWriter,
        mut result_chromatograms: Option<&mut (dyn MSDataConsumer + '_)>,
        batch_size: i32,
        ms1_isotopes: i32,
        load_into_memory: bool,
        mrm_mapping_param: &Param,
        mut mobilogram_consumer: Option<&mut (dyn MobilogramSink + '_)>,
        inner_batch_size: i32,
        max_concurrent_swaths: i32,
    ) -> Result<()> {
        let cp = chromatogram_extraction_params;
        let user_inner_batch_size = inner_batch_size;
        let user_max_concurrent_swaths = max_concurrent_swaths;
        // every identifier is keyed by this and the swath and assay it belongs
        // to, so none depends on the thread count or the scheduling
        let base_seed = self.draw_seed();
        let ms1_only = swath_maps.len() == 1 && swath_maps[0].ms1;
        let store_features = store_features_in_feature_file;

        if self.base.mrm {
            self.base.start_progress(0, 1, "Extraction and Scoring")?;
            let Some(handler) = self.chrom_handler.clone() else {
                return Err(Error::Unsupported(
                    "SRM/MRM mode needs a ChromatogramHandler (ANALYSIS/TARGETED/IChromatogramHandler, \
                     DefaultChromHandler, MRMChromHandler)"
                        .into(),
                ));
            };
            let mut filtered = handler.extract_and_map_chromatograms_for_transitions(
                swath_maps,
                assay_library,
                cp,
                mrm_mapping_param,
            )?;
            let mut feature_file = FeatureMap::new();
            let seed = mix_seed(base_seed, u64::MAX);
            self.score_all_chromatograms_seeded(
                &filtered,
                &[],
                swath_maps,
                assay_library,
                feature_finder_param,
                rt_trafo,
                cp.rt_extraction_window,
                &mut feature_file,
                result_osw,
                ms1_isotopes,
                false,
                mobilogram_consumer.as_deref_mut(),
                None,
                seed,
            )?;
            let mut empty_ms1: Vec<MSChromatogram> = Vec::new();
            write_out_features_and_chroms(
                &mut filtered,
                &mut empty_ms1,
                &feature_file,
                result_feature_file,
                store_features,
                result_chromatograms.as_deref_mut(),
            )?;
            self.base.end_progress(0)?;
            return Ok(());
        }

        // compute the inversion of the transformation
        let mut trafo_inverse = rt_trafo.clone();
        trafo_inverse.invert()?;

        self.log(
            LogLevel::Info,
            format!(
                "Will analyze {} transitions in total.",
                assay_library.transitions.len()
            ),
        );
        let mut progress: i64 = 0;
        self.base.start_progress(
            0,
            swath_maps.len() as i64,
            "Extracting and scoring transitions",
        )?;

        // (i) obtain precursor chromatograms (MS1) if precursor extraction is enabled
        let ms1_cp = ms1_chromatogram_extraction_params.clone();
        if ms1_only && !self.base.use_ms1_traces {
            return Err(invalid(
                "Error, you need to enable use_ms1_traces when run in MS1 mode.",
            ));
        }
        if self.base.use_ms1_traces {
            self.base.ms1_map = self.load_ms1_map(swath_maps, load_into_memory)?;
        }

        // (ii) precursor extraction only
        if ms1_only {
            let ms1_map = self
                .base
                .ms1_map
                .clone()
                .ok_or_else(|| invalid("the MS1 map has no data"))?;
            let mut ms1_chromatograms = self.base.ms1_extraction(
                &*ms1_map,
                swath_maps,
                &ms1_cp,
                assay_library,
                &trafo_inverse,
                ms1_only,
                ms1_isotopes,
            )?;
            let mut feature_file = FeatureMap::new();
            let seed = mix_seed(base_seed, u64::MAX - 1);
            self.score_all_chromatograms_seeded(
                &[],
                &ms1_chromatograms,
                swath_maps,
                assay_library,
                feature_finder_param,
                rt_trafo,
                cp.rt_extraction_window,
                &mut feature_file,
                result_osw,
                ms1_isotopes,
                true,
                mobilogram_consumer.as_deref_mut(),
                None,
                seed,
            )?;
            let mut chromatograms: Vec<MSChromatogram> = Vec::new();
            write_out_features_and_chroms(
                &mut chromatograms,
                &mut ms1_chromatograms,
                &feature_file,
                result_feature_file,
                store_features,
                result_chromatograms.as_deref_mut(),
            )?;
        }

        // (iii) map transitions to individual DIA windows for the cases where
        // this is non-trivial (overlapping windows)
        let mut tr_win_map: Vec<Option<usize>> = Vec::new();
        if self.base.prm && self.base.pasef {
            self.log(
                LogLevel::Error,
                "Setting -pasef and -matching_window_only flags simultaneously is not currently supported.",
            );
            return Err(Error::Unsupported(
                "PRM and diaPASEF matching cannot be combined".into(),
            ));
        } else if self.base.prm {
            tr_win_map = select_prm_windows(assay_library, cp.min_upper_edge_dist, swath_maps);
        } else if self.base.pasef {
            tr_win_map =
                select_swath_transitions_pasef(assay_library, cp.min_upper_edge_dist, swath_maps);
        }

        let scoring_threads = self.threads.get().max(1);
        let mut scheduler_options = SchedulerOptions {
            max_concurrent_swaths: user_max_concurrent_swaths,
            scoring_threads,
            ..SchedulerOptions::default()
        };
        let non_ms1_swath_count = swath_maps.iter().filter(|m| !m.ms1).count();
        if non_ms1_swath_count > 0 {
            // chromatogram storage dominates the memory for large libraries:
            // use the transition density in the wave estimate
            scheduler_options.avg_transitions_per_swath = assay_library
                .transitions
                .len()
                .div_ceil(non_ms1_swath_count);
        }
        scheduler_options.osw_buffer_bytes = if result_osw.is_active() {
            determine_osw_write_buffer_bytes()
        } else {
            0
        };
        let estimate =
            OpenSwathWorkflowScheduler::estimate_concurrency(swath_maps, &scheduler_options);

        // the nested OpenMP variant is not part of this port, so it never
        // vetoes the wave scheduler
        let use_swath_range_scheduler = batch_size <= 0 && load_into_memory && !ms1_only;

        if use_swath_range_scheduler {
            self.run_wave_scheduler(
                swath_maps,
                &trafo_inverse,
                rt_trafo,
                cp,
                &ms1_cp,
                feature_finder_param,
                assay_library,
                result_feature_file,
                store_features,
                result_osw,
                result_chromatograms.as_deref_mut(),
                ms1_isotopes,
                mobilogram_consumer.as_deref_mut(),
                user_inner_batch_size,
                &scheduler_options,
                &estimate,
                &tr_win_map,
                &mut progress,
                base_seed,
            )?;
            self.base.end_progress(0)?;
            return Ok(());
        }

        // (iv) extraction and scoring of the fragment ion chromatograms, one
        // SWATH after the other, in batches of compounds
        self.run_per_swath(
            swath_maps,
            &trafo_inverse,
            rt_trafo,
            cp,
            &ms1_cp,
            feature_finder_param,
            assay_library,
            result_feature_file,
            store_features,
            result_osw,
            result_chromatograms,
            batch_size,
            ms1_isotopes,
            load_into_memory,
            mobilogram_consumer,
            user_inner_batch_size,
            ms1_only,
            &estimate,
            &tr_win_map,
            &mut progress,
            base_seed,
        )?;
        self.base.end_progress(0)?;
        Ok(())
    }

    /// Extract, for one SWATH, the transitions to use (source step 1).
    fn swath_transitions(
        &self,
        transition_exp: &LightTargetedExperiment,
        swath_maps: &[SwathMap],
        index: usize,
        cp: &ChromExtractParams,
        tr_win_map: &[Option<usize>],
    ) -> LightTargetedExperiment {
        if !(self.base.prm || self.base.pasef) {
            let mut used = LightTargetedExperiment::new();
            select_swath_transitions(
                transition_exp,
                &mut used,
                cp.min_upper_edge_dist,
                swath_maps[index].lower,
                swath_maps[index].upper,
            );
            used
        } else {
            select_assigned_transitions(transition_exp, tr_win_map, index)
        }
    }

    /// Extract the MS1 and MS2 chromatograms of one set of assays.
    #[allow(clippy::too_many_arguments)]
    fn extract_batch(
        &self,
        swath_access: &dyn SpectrumAccess,
        swath_maps: &[SwathMap],
        transition_exp_used: &LightTargetedExperiment,
        trafo_inverse: &TransformationDescription,
        cp: &ChromExtractParams,
        ms1_cp: &ChromExtractParams,
        ms1_isotopes: i32,
        ms1_only: bool,
    ) -> Result<(Vec<MSChromatogram>, Vec<MSChromatogram>)> {
        let mut ms1_chromatograms = Vec::new();
        if let Some(ms1_map) = &self.base.ms1_map {
            let threadsafe_ms1 = ms1_map.light_clone()?;
            ms1_chromatograms = self.base.ms1_extraction(
                &*threadsafe_ms1,
                swath_maps,
                ms1_cp,
                transition_exp_used,
                trafo_inverse,
                ms1_only,
                ms1_isotopes,
            )?;
        }
        let (mut chrom_list, coordinates) = self.base.prepare_extraction_coordinates(
            transition_exp_used,
            trafo_inverse,
            cp,
            false,
            -1,
        )?;
        ChromatogramExtractor::new().extract_chromatograms_im(
            swath_access,
            &mut chrom_list,
            &coordinates,
            cp.mz_extraction_window,
            cp.ppm,
            cp.im_extraction_window,
            &cp.extraction_function,
        )?;
        let mut chrom_exp = Vec::new();
        return_chromatogram(
            &chrom_list,
            &coordinates,
            transition_exp_used,
            &MSSpectrum::default(),
            &mut chrom_exp,
            false,
            cp.im_extraction_window,
        )?;
        Ok((chrom_exp, ms1_chromatograms))
    }

    /// The SWATH wave scheduler path: waves of SWATHs are extracted in
    /// parallel, then their scoring jobs.
    #[allow(clippy::too_many_arguments)]
    fn run_wave_scheduler(
        &self,
        swath_maps: &[SwathMap],
        trafo_inverse: &TransformationDescription,
        trafo: &TransformationDescription,
        cp: &ChromExtractParams,
        ms1_cp: &ChromExtractParams,
        feature_finder_param: &Param,
        transition_exp: &LightTargetedExperiment,
        out_feature_file: &mut FeatureMap,
        store_features: bool,
        osw_writer: &OpenSwathOSWWriter,
        mut chrom_consumer: Option<&mut (dyn MSDataConsumer + '_)>,
        ms1_isotopes: i32,
        mut mobilogram_consumer: Option<&mut (dyn MobilogramSink + '_)>,
        user_inner_batch_size: i32,
        scheduler_options: &SchedulerOptions,
        estimate: &ConcurrencyEstimate,
        tr_win_map: &[Option<usize>],
        progress: &mut i64,
        base_seed: u64,
    ) -> Result<()> {
        let scoring_threads = scheduler_options.scoring_threads;
        let waves = OpenSwathWorkflowScheduler::plan_waves(swath_maps, scheduler_options, estimate);
        self.log(
            LogLevel::Info,
            format!(
                "Use SWATH wave scheduler with {} waves, max_concurrent_swaths={}, memory_budget={}, avg_spectra={}, avg_transitions={}, estimated_swath={}.",
                waves.len(),
                estimate.max_concurrent_swaths,
                estimate.memory_budget_bytes,
                estimate.avg_spectra_per_swath,
                scheduler_options.avg_transitions_per_swath,
                estimate.estimated_bytes_per_swath
            ),
        );
        for map in swath_maps {
            if map.ms1 {
                *progress += 1;
                self.base.set_progress(*progress)?;
            }
        }
        self.log(
            LogLevel::Info,
            format!("Use SWATH wave scheduler with {scoring_threads} scoring threads."),
        );

        let mut contexts: Vec<SwathContext> = (0..swath_maps.len())
            .map(|_| SwathContext::default())
            .collect();
        let threads = self.threads;

        for wave in &waves {
            // stage 1: select and extract every SWATH of the wave
            let prepared: Vec<Result<SwathContext>> =
                map_collect(&wave.swath_indices, threads, |&swath_index| {
                    let mut context = SwathContext {
                        transition_exp_used_all: self.swath_transitions(
                            transition_exp,
                            swath_maps,
                            swath_index,
                            cp,
                            tr_win_map,
                        ),
                        ..SwathContext::default()
                    };
                    if context.transition_exp_used_all.transitions.is_empty() {
                        return Ok(context);
                    }
                    let sptr = swath_maps[swath_index]
                        .sptr
                        .as_ref()
                        .ok_or_else(|| invalid(format!("SWATH map {swath_index} has no data")))?;
                    let current = in_memory(sptr)?;
                    let (chrom_exp, ms1_chromatograms) = self.extract_batch(
                        &*current,
                        swath_maps,
                        &context.transition_exp_used_all,
                        trafo_inverse,
                        cp,
                        ms1_cp,
                        ms1_isotopes,
                        false,
                    )?;
                    context.current_swath_map = Some(current);
                    context.chrom_exp = chrom_exp;
                    context.ms1_chromatograms = ms1_chromatograms;
                    context.active = true;
                    Ok(context)
                });
            for (position, result) in prepared.into_iter().enumerate() {
                contexts[wave.swath_indices[position]] = result?;
            }

            // create the inner scoring jobs of all active SWATHs
            let mut wave_compounds = 0;
            let mut active_swaths = 0;
            for &swath_index in &wave.swath_indices {
                let context = &contexts[swath_index];
                if !context.active {
                    continue;
                }
                wave_compounds += context.transition_exp_used_all.compounds.len();
                active_swaths += 1;
            }
            let wave_inner_batch = OpenSwathWorkflowScheduler::choose_inner_batch_size(
                wave_compounds,
                active_swaths,
                scoring_threads,
                user_inner_batch_size,
                scheduler_options,
            );

            let mut queue: Vec<InnerBatchScoringJob> = Vec::new();
            for &swath_index in &wave.swath_indices {
                let context = &mut contexts[swath_index];
                if !context.active {
                    continue;
                }
                let n_compounds = context.transition_exp_used_all.compounds.len();
                let inner_batch = wave_inner_batch.min(n_compounds);
                context.score_job_size = inner_batch;
                if context.score_job_size > 0 {
                    context.nr_score_jobs = n_compounds.div_ceil(context.score_job_size);
                    context.remaining_score_jobs = context.nr_score_jobs;
                }
                let n_transitions = context.transition_exp_used_all.transitions.len();
                for batch_index in 0..context.nr_score_jobs {
                    let batch_start = batch_index * context.score_job_size;
                    let estimated_compounds = if batch_start < n_compounds {
                        context.score_job_size.min(n_compounds - batch_start)
                    } else {
                        0
                    };
                    // ceil(estimated_compounds * n_transitions / n_compounds); no
                    // compounds, no estimated compounds
                    let estimated_transitions =
                        (estimated_compounds * n_transitions).div_ceil(n_compounds.max(1));
                    queue.push(InnerBatchScoringJob {
                        context_index: swath_index,
                        batch_index,
                        estimated_compounds,
                        estimated_transitions,
                    });
                }
            }
            let wave_has_parallel_scoring = queue.len() > 1;
            // largest jobs first, then by batch and context
            source_sort_by(&mut queue, |lhs, rhs| {
                if lhs.estimated_transitions != rhs.estimated_transitions {
                    return lhs.estimated_transitions > rhs.estimated_transitions;
                }
                if lhs.estimated_compounds != rhs.estimated_compounds {
                    return lhs.estimated_compounds > rhs.estimated_compounds;
                }
                if lhs.batch_index != rhs.batch_index {
                    return lhs.batch_index < rhs.batch_index;
                }
                lhs.context_index < rhs.context_index
            })?;

            // stage 2: score the jobs in chunks, consume in job order
            let chunk_size = if mobilogram_consumer.is_some() {
                1
            } else {
                (threads.get().max(1) * 4).max(1)
            };
            let mut next_job = 0;
            while next_job < queue.len() {
                let chunk_end = (next_job + chunk_size).min(queue.len());
                let chunk: Vec<(InnerBatchScoringJob, u64)> = queue[next_job..chunk_end]
                    .iter()
                    .map(|job| (job.clone(), mix_seed(base_seed, job.context_index as u64)))
                    .collect();
                next_job = chunk_end;

                let run = |entry: &(InnerBatchScoringJob, u64),
                           mobilogram: Option<&mut (dyn MobilogramSink + '_)>|
                 -> Result<ScoredJob> {
                    let (job, seed) = entry;
                    let context = &contexts[job.context_index];
                    let n_compounds = context.transition_exp_used_all.compounds.len();
                    let inner_batch = context.score_job_size.min(n_compounds);
                    let transition_exp_used = select_compounds_for_batch(
                        &context.transition_exp_used_all,
                        inner_batch,
                        job.batch_index,
                    );
                    let mut tmp = vec![swath_maps[job.context_index].clone()];
                    let current = context
                        .current_swath_map
                        .as_ref()
                        .ok_or_else(|| invalid("active SWATH without data"))?;
                    tmp[0].sptr = Some(current.light_clone()?);
                    let mut features = FeatureMap::new();
                    let mut osw = OswData::default();
                    self.score_all_chromatograms_seeded(
                        &context.chrom_exp,
                        &context.ms1_chromatograms,
                        &tmp,
                        &transition_exp_used,
                        feature_finder_param,
                        trafo,
                        cp.rt_extraction_window,
                        &mut features,
                        osw_writer,
                        ms1_isotopes,
                        false,
                        mobilogram,
                        Some(&mut osw),
                        *seed,
                    )?;
                    Ok(ScoredJob { features, osw })
                };
                let results: Vec<Result<ScoredJob>> =
                    if let Some(consumer) = mobilogram_consumer.as_deref_mut() {
                        let mut results = Vec::with_capacity(chunk.len());
                        for entry in &chunk {
                            results.push(run(entry, Some(&mut *consumer)));
                        }
                        results
                    } else {
                        map_collect(&chunk, threads, |entry| run(entry, None))
                    };

                for (entry, result) in chunk.iter().zip(results) {
                    let scored = result?;
                    let job = &entry.0;
                    if osw_writer.is_active() && !scored.osw.is_empty() {
                        osw_writer.write_rows(&scored.osw)?;
                    }
                    let context = &mut contexts[job.context_index];
                    if !osw_writer.is_active() && store_features {
                        context
                            .feature_file
                            .features
                            .extend(scored.features.features);
                        context
                            .feature_file
                            .protein_identifications
                            .extend(scored.features.protein_identifications);
                    }
                    context.remaining_score_jobs -= 1;
                    if context.remaining_score_jobs == 0 {
                        self.finalize_context(
                            &mut contexts[job.context_index],
                            out_feature_file,
                            store_features,
                            chrom_consumer.as_deref_mut(),
                            progress,
                        )?;
                    }
                }
            }

            for &swath_index in &wave.swath_indices {
                let context = &mut contexts[swath_index];
                if context.active && !context.finalized {
                    self.finalize_context(
                        context,
                        out_feature_file,
                        store_features,
                        chrom_consumer.as_deref_mut(),
                        progress,
                    )?;
                } else if !context.active {
                    *progress += 1;
                    self.base.set_progress(*progress)?;
                }
            }

            for &swath_index in &wave.swath_indices {
                contexts[swath_index] = SwathContext {
                    finalized: true,
                    ..SwathContext::default()
                };
            }

            if store_features && wave_has_parallel_scoring {
                // parallel scoring jobs append their peak groups in
                // completion order: sort once after the wave
                sort_features_by_position(out_feature_file)?;
            }
        }
        Ok(())
    }

    fn finalize_context(
        &self,
        context: &mut SwathContext,
        out_feature_file: &mut FeatureMap,
        store_features: bool,
        chrom_consumer: Option<&mut (dyn MSDataConsumer + '_)>,
        progress: &mut i64,
    ) -> Result<()> {
        if !context.active || context.finalized {
            return Ok(());
        }
        write_out_features_and_chroms(
            &mut context.chrom_exp,
            &mut context.ms1_chromatograms,
            &context.feature_file,
            out_feature_file,
            store_features,
            chrom_consumer,
        )?;
        context.finalized = true;
        *progress += 1;
        self.base.set_progress(*progress)?;
        Ok(())
    }

    /// The per-SWATH path: for each SWATH window, batches of compounds are
    /// extracted, scored and written.
    #[allow(clippy::too_many_arguments)]
    fn run_per_swath(
        &self,
        swath_maps: &[SwathMap],
        trafo_inverse: &TransformationDescription,
        trafo: &TransformationDescription,
        cp: &ChromExtractParams,
        ms1_cp: &ChromExtractParams,
        feature_finder_param: &Param,
        transition_exp: &LightTargetedExperiment,
        out_feature_file: &mut FeatureMap,
        store_features: bool,
        osw_writer: &OpenSwathOSWWriter,
        mut chrom_consumer: Option<&mut (dyn MSDataConsumer + '_)>,
        batch_size_arg: i32,
        ms1_isotopes: i32,
        load_into_memory: bool,
        mut mobilogram_consumer: Option<&mut (dyn MobilogramSink + '_)>,
        user_inner_batch_size: i32,
        ms1_only: bool,
        estimate: &ConcurrencyEstimate,
        tr_win_map: &[Option<usize>],
        progress: &mut i64,
        base_seed: u64,
    ) -> Result<()> {
        let limiter = if load_into_memory
            && !ms1_only
            && estimate.non_ms1_swath_count > estimate.max_concurrent_swaths
        {
            self.log(
                LogLevel::Info,
                format!(
                    "Limit concurrently loaded SWATHs to {} in the streaming scheduler.",
                    estimate.max_concurrent_swaths
                ),
            );
            Some(ConcurrencyLimiter::new(estimate.max_concurrent_swaths))
        } else {
            None
        };
        let threads = if mobilogram_consumer.is_some() {
            Threads::serial()
        } else {
            self.threads
        };

        // the SWATHs are processed in chunks, the batches of a chunk in
        // parallel; the outputs are consumed in SWATH and batch order
        let indices: Vec<usize> = (0..swath_maps.len()).collect();
        let chunk_size = if mobilogram_consumer.is_some() {
            1
        } else {
            threads.get().max(1)
        };
        for chunk in indices.chunks(chunk_size) {
            let seeds: Vec<u64> = chunk
                .iter()
                .map(|&i| mix_seed(base_seed, i as u64))
                .collect();
            let entries: Vec<(usize, u64)> = chunk.iter().copied().zip(seeds).collect();
            let process = |entry: &(usize, u64),
                           mut mobilogram: Option<&mut (dyn MobilogramSink + '_)>|
             -> Result<Vec<ScoredBatch>> {
                let (i, seed) = *entry;
                if swath_maps[i].ms1 {
                    return Ok(Vec::new());
                }
                let transition_exp_used_all =
                    self.swath_transitions(transition_exp, swath_maps, i, cp, tr_win_map);
                if transition_exp_used_all.transitions.is_empty() {
                    return Ok(Vec::new());
                }
                let _slot = ScopedSlot::new(limiter.as_ref());
                let sptr = swath_maps[i]
                    .sptr
                    .as_ref()
                    .ok_or_else(|| invalid(format!("SWATH map {i} has no data")))?;
                let current_swath_map: SpectrumAccessPtr = if load_into_memory {
                    in_memory(sptr)?
                } else {
                    sptr.light_clone()?
                };

                let n_compounds = transition_exp_used_all.compounds.len();
                let batch_size = if batch_size_arg <= 0 {
                    // automatic: bound the memory of the extracted batch
                    calculate_inner_batch_size(n_compounds, user_inner_batch_size).max(1)
                } else if batch_size_arg as usize >= n_compounds {
                    n_compounds
                } else {
                    batch_size_arg as usize
                };
                let nr_batches = if batch_size > 0 {
                    n_compounds.div_ceil(batch_size)
                } else {
                    0
                };
                let mut batches = Vec::with_capacity(nr_batches);
                for pep_idx in 0..nr_batches {
                    let transition_exp_used =
                        select_compounds_for_batch(&transition_exp_used_all, batch_size, pep_idx);
                    let (chromatograms, ms1_chromatograms) = self.extract_batch(
                        &*current_swath_map,
                        swath_maps,
                        &transition_exp_used,
                        trafo_inverse,
                        cp,
                        ms1_cp,
                        ms1_isotopes,
                        ms1_only,
                    )?;
                    let mut tmp = vec![swath_maps[i].clone()];
                    tmp[0].sptr = Some(Arc::clone(&current_swath_map));
                    let mut features = FeatureMap::new();
                    let mut osw = OswData::default();
                    self.score_all_chromatograms_seeded(
                        &chromatograms,
                        &ms1_chromatograms,
                        &tmp,
                        &transition_exp_used,
                        feature_finder_param,
                        trafo,
                        cp.rt_extraction_window,
                        &mut features,
                        osw_writer,
                        ms1_isotopes,
                        false,
                        mobilogram.as_deref_mut(),
                        Some(&mut osw),
                        seed,
                    )?;
                    batches.push(ScoredBatch {
                        chromatograms,
                        ms1_chromatograms,
                        features,
                        osw,
                    });
                }
                Ok(batches)
            };
            let results: Vec<Result<Vec<ScoredBatch>>> =
                if let Some(consumer) = mobilogram_consumer.as_deref_mut() {
                    let mut results = Vec::with_capacity(entries.len());
                    for entry in &entries {
                        results.push(process(entry, Some(&mut *consumer)));
                    }
                    results
                } else {
                    map_collect(&entries, threads, |entry| process(entry, None))
                };
            for result in results {
                for mut batch in result? {
                    if osw_writer.is_active() && !batch.osw.is_empty() {
                        osw_writer.write_rows(&batch.osw)?;
                    }
                    write_out_features_and_chroms(
                        &mut batch.chromatograms,
                        &mut batch.ms1_chromatograms,
                        &batch.features,
                        out_feature_file,
                        store_features,
                        chrom_consumer.as_deref_mut(),
                    )?;
                }
                *progress += 1;
                self.base.set_progress(*progress)?;
            }
        }
        Ok(())
    }

    /// Score a set of chromatograms (source `scoreAllChromatograms_`).
    ///
    /// For each assay (compound or peptide) of `transition_exp` a transition
    /// group is built from its transitions and the extracted chromatograms
    /// (cut to the expected RT of the assay `+/-` `rt_extraction_window` when
    /// that is positive), the MS1 chromatograms of the precursor isotopes are
    /// added, the peak groups are picked and scored with
    /// [`MRMFeatureFinderScoring`], and the scored features are appended to
    /// `output` and, as rows, to `deferred_osw_output` (or written at once
    /// with `osw_writer`).
    ///
    /// With an active OSW writer `output` is cleared per assay: the OSW and
    /// the feature map outputs are exclusive.
    ///
    /// In SRM/MRM mode (`mrm`) a transition without chromatogram, an assay
    /// without detecting chromatogram and an assay whose picking or scoring
    /// fails with [`Error::InvalidValue`] or [`Error::InvalidRange`] are
    /// skipped; in DIA mode these errors are returned. An assay that fails
    /// with any other error is skipped in both modes, as the source catches
    /// every `std::exception` there.
    ///
    /// The identifiers of the new features are drawn from the generator of
    /// the workflow.
    ///
    /// # Errors
    ///
    /// As above; also the errors of the OSW row preparation and writing, and
    /// [`Error::InvalidValue`] when the library has a detecting-less feature
    /// ("did not find any detection transition").
    #[allow(clippy::too_many_arguments)]
    pub fn score_all_chromatograms(
        &self,
        ms2_chromatograms: &[MSChromatogram],
        ms1_chromatograms: &[MSChromatogram],
        swath_maps: &[SwathMap],
        transition_exp: &LightTargetedExperiment,
        feature_finder_param: &Param,
        trafo: &TransformationDescription,
        rt_extraction_window: f64,
        output: &mut FeatureMap,
        osw_writer: &OpenSwathOSWWriter,
        nr_ms1_isotopes: i32,
        ms1only: bool,
        mobilogram_consumer: Option<&mut (dyn MobilogramSink + '_)>,
        deferred_osw_output: Option<&mut OswData>,
    ) -> Result<()> {
        let seed = self.draw_seed();
        self.score_all_chromatograms_seeded(
            ms2_chromatograms,
            ms1_chromatograms,
            swath_maps,
            transition_exp,
            feature_finder_param,
            trafo,
            rt_extraction_window,
            output,
            osw_writer,
            nr_ms1_isotopes,
            ms1only,
            mobilogram_consumer,
            deferred_osw_output,
            seed,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn score_all_chromatograms_seeded(
        &self,
        ms2_chromatograms: &[MSChromatogram],
        ms1_chromatograms: &[MSChromatogram],
        swath_maps: &[SwathMap],
        transition_exp: &LightTargetedExperiment,
        feature_finder_param: &Param,
        trafo: &TransformationDescription,
        rt_extraction_window: f64,
        output: &mut FeatureMap,
        osw_writer: &OpenSwathOSWWriter,
        nr_ms1_isotopes: i32,
        ms1only: bool,
        mut mobilogram_consumer: Option<&mut (dyn MobilogramSink + '_)>,
        deferred_osw_output: Option<&mut OswData>,
        seed: u64,
    ) -> Result<()> {
        let mut trafo_inv = trafo.clone();
        trafo_inv.invert()?;

        let mut feature_finder = MRMFeatureFinderScoring::from_param(feature_finder_param)?;

        // to ensure thread-safe access to the spectra, a light clone of the
        // MS1 map is used
        if self.base.use_ms1_traces {
            if let Some(ms1_map) = &self.base.ms1_map {
                feature_finder.set_ms1_map(Some(ms1_map.light_clone()?));
            } else {
                self.log(
                    LogLevel::Warn,
                    "WARNING: Attempted to use MS1 traces but no MS1 map was provided: Will not use MS1 signal!",
                );
            }
        }

        // if use_total_mi_score is defined, the picker computes the score
        let mut picker_param = feature_finder_param.copy("TransitionGroupPicker:", true)?;
        if feature_finder_param
            .value("Scores:use_total_mi_score")?
            .to_bool()?
        {
            picker_param.set_value(
                "compute_total_mi",
                crate::param::ParamValue::String("true".into()),
                "",
                &[],
            )?;
        }
        let picker = MRMTransitionGroupPicker::from_param(&picker_param)?;
        feature_finder.prepare_protein_peptide_maps(transition_exp);

        let ms1_chromatogram_map: HashMap<&str, usize> = ms1_chromatograms
            .iter()
            .enumerate()
            .map(|(i, c)| (c.native_id.as_str(), i))
            .collect();
        let chromatogram_map: HashMap<&str, usize> = ms2_chromatograms
            .iter()
            .enumerate()
            .map(|(i, c)| (c.native_id.as_str(), i))
            .collect();
        let mut assay_peptide_map: HashMap<&str, usize> = HashMap::new();
        for (i, compound) in transition_exp.compounds.iter().enumerate() {
            assay_peptide_map.insert(compound.id.as_str(), i);
        }

        // map peptide id to its transitions; every compound gets an entry
        // even without transitions (e.g. MS1 only)
        let mut assay_map: BTreeMap<&str, Vec<&LightTransition>> = BTreeMap::new();
        for compound in &transition_exp.compounds {
            assay_map.insert(compound.id.as_str(), Vec::new());
        }
        for transition in &transition_exp.transitions {
            assay_map
                .entry(transition.peptide_ref.as_str())
                .or_default()
                .push(transition);
        }

        let mut to_osw_output = OswData::default();
        if osw_writer.is_active() {
            let stop_report = feature_finder_param
                .value("stop_report_after_feature")?
                .to_i64()?;
            let expected = if stop_report > 0 {
                stop_report as usize
            } else {
                5
            };
            to_osw_output.reserve_rows(
                transition_exp.compounds.len().saturating_mul(expected),
                transition_exp.transitions.len().saturating_mul(expected),
            );
        }

        for (&id, transitions) in &assay_map {
            let mut transition_group = MRMTransitionGroup::<MSChromatogram, LightTransition>::new();
            transition_group.set_transition_group_id(id);
            // an unknown assay reads the first compound, as the source's
            // `operator[]` default-inserts index 0
            let compound_index = assay_peptide_map.get(id).copied().unwrap_or(0);
            let expected_rt = transition_exp
                .compounds
                .get(compound_index)
                .ok_or_else(|| invalid(format!("no compound for assay '{id}'")))?
                .rt;

            // 1. add every transition and its chromatogram to the group
            let mut detection_assay: Option<&LightTransition> = None;
            for &transition in transitions {
                if transition.is_detecting_transition() {
                    detection_assay = Some(transition);
                }
                // MS1-only: there are no chromatograms for the transitions
                if ms1only {
                    continue;
                }
                let Some(&chrom_index) = chromatogram_map.get(transition.native_id()) else {
                    if self.base.mrm {
                        // in SRM/MRM mode missing chromatograms are skipped
                        continue;
                    }
                    return Err(invalid(format!(
                        "Error, did not find chromatogram for transition {}",
                        transition.native_id()
                    )));
                };
                let mut chromatogram = ms2_chromatograms[chrom_index].clone();
                chromatogram.native_id = transition.native_id().to_owned();
                if rt_extraction_window > 0.0 {
                    let de_normalized_rt = trafo_inv.apply(expected_rt)?;
                    let rt_max = de_normalized_rt + rt_extraction_window;
                    let rt_min = de_normalized_rt - rt_extraction_window;
                    let keep: Vec<usize> = chromatogram
                        .peaks
                        .iter()
                        .enumerate()
                        .filter(|(_, p)| p.rt >= rt_min && p.rt <= rt_max)
                        .map(|(i, _)| i)
                        .collect();
                    if keep.len() != chromatogram.peaks.len() {
                        chromatogram.select(&keep)?;
                    }
                }
                let native_id = transition.native_id().to_owned();
                transition_group.add_transition(transition.clone(), &native_id)?;
                transition_group.add_chromatogram(chromatogram, &native_id)?;
            }

            // .osw and .featureXML are mutually exclusive
            if osw_writer.is_active() {
                output.clear(true);
            }

            // 2. the MS1 chromatograms of the precursor isotopes
            for iso in 0..=nr_ms1_isotopes {
                let prec_id = compute_precursor_id(transition_group.transition_group_id(), iso);
                if !ms1_chromatograms.is_empty() {
                    if let Some(&index) = ms1_chromatogram_map.get(prec_id.as_str()) {
                        let chromatogram = ms1_chromatograms[index].clone();
                        let native_id = chromatogram.native_id.clone();
                        transition_group.add_precursor_chromatogram(chromatogram, &native_id)?;
                    }
                }
            }

            // 3. / 4. find and score the peak groups; with no chromatograms
            // at all (every transition skipped) the picker would fail
            if transition_group.chromatograms().is_empty()
                && transition_group.precursor_chromatograms().is_empty()
            {
                continue;
            }

            // the picker needs at least one chromatogram of a detecting
            // transition
            let mut has_detecting_chrom = false;
            for chrom in transition_group.chromatograms() {
                let native_id = chrom.native_id.as_str();
                if transition_group.transitions().is_empty() {
                    has_detecting_chrom = true;
                    break;
                }
                if transition_group.has_transition(native_id)
                    && transition_group
                        .transition(native_id)?
                        .is_detecting_transition()
                {
                    has_detecting_chrom = true;
                    break;
                }
            }
            if !has_detecting_chrom && !transition_group.precursor_chromatograms().is_empty() {
                has_detecting_chrom = true;
            }
            if !has_detecting_chrom {
                if self.base.mrm {
                    continue;
                }
                return Err(invalid(format!(
                    "Error, did not find any detecting chromatogram for assay {id}"
                )));
            }

            // the identifiers of this assay's features: the picker's and the
            // finder's (subordinates) streams are keyed by the swath (`seed`)
            // and the assay, never by the schedule
            let assay_key = assay_seed(seed, id);
            let mut ids = UniqueIdGenerator::from_seed(assay_key);
            feature_finder
                .set_unique_id_generator(UniqueIdGenerator::from_seed(mix_seed(assay_key, 1)));
            match picker.pick_transition_group(&mut transition_group, &mut ids) {
                Ok(()) => {}
                Err(error) => {
                    if self.skippable(&error, id, "picking") {
                        continue;
                    }
                    return Err(error);
                }
            }
            match feature_finder.score_peakgroups(
                &mut transition_group,
                trafo,
                swath_maps,
                output,
                ms1only,
                reborrow_sink(&mut mobilogram_consumer),
            ) {
                Ok(()) => {}
                Err(error) => {
                    if self.skippable(&error, id, "scoring") {
                        continue;
                    }
                    return Err(error);
                }
            }

            // a detecting transition is needed to derive features
            if detection_assay.is_none() && !output.features.is_empty() {
                return Err(invalid(format!(
                    "Error, did not find any detection transition for feature {id}"
                )));
            }

            // 5. add to the OSW output
            if osw_writer.is_active() && !output.features.is_empty() {
                osw_writer.prepare_rows_into(
                    &mut to_osw_output,
                    &LightCompound::default(),
                    None,
                    output,
                    id,
                )?;
            }
        }

        // write only at the very end: one barrier per call
        if osw_writer.is_active() {
            if let Some(deferred) = deferred_osw_output {
                deferred.append(to_osw_output);
                return Ok(());
            }
            osw_writer.write_rows(&to_osw_output)?;
        }
        Ok(())
    }

    /// Whether an error of picking or scoring skips the assay: always for
    /// errors other than [`Error::InvalidValue`] and [`Error::InvalidRange`]
    /// (the source's catch of `std::exception`), and in SRM/MRM mode for
    /// those two.
    fn skippable(&self, error: &Error, id: &str, stage: &str) -> bool {
        let mode_specific = matches!(error, Error::InvalidValue(_) | Error::InvalidRange(_));
        if mode_specific && !self.base.mrm {
            return false;
        }
        self.log(
            LogLevel::Error,
            format!(
                "Exception while {stage} transition group for assay {id}: {error} - skipping assay."
            ),
        );
        true
    }
}

/// Write the chromatograms and the features of one batch (source
/// `writeOutFeaturesAndChroms_`): the non-empty MS1 chromatograms, then the
/// non-empty chromatograms to the consumer, then, with `store_features`, the
/// features and the protein identifications of `feature_file` to
/// `out_feature_file`.
///
/// A consumer that answers [`ControlFlow::Break`] receives no further
/// chromatograms of this call. Without a consumer the chromatograms are not
/// written (the source dereferences a null consumer).
///
/// # Errors
///
/// The errors of the consumer.
pub fn write_out_features_and_chroms(
    chromatograms: &mut [MSChromatogram],
    ms1_chromatograms: &mut [MSChromatogram],
    feature_file: &FeatureMap,
    out_feature_file: &mut FeatureMap,
    store_features: bool,
    chrom_consumer: Option<&mut (dyn MSDataConsumer + '_)>,
) -> Result<()> {
    if let Some(consumer) = chrom_consumer {
        'consume: {
            for chromatogram in ms1_chromatograms.iter_mut() {
                if chromatogram.peaks.is_empty() {
                    continue;
                }
                if consumer.consume_chromatogram(chromatogram)? == ControlFlow::Break(()) {
                    break 'consume;
                }
            }
            for chromatogram in chromatograms.iter_mut() {
                if chromatogram.peaks.is_empty() {
                    continue;
                }
                if consumer.consume_chromatogram(chromatogram)? == ControlFlow::Break(()) {
                    break 'consume;
                }
            }
        }
    }
    if store_features {
        out_feature_file
            .features
            .extend(feature_file.features.iter().cloned());
        out_feature_file
            .protein_identifications
            .extend(feature_file.protein_identifications.iter().cloned());
    }
    Ok(())
}

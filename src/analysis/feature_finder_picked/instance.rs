// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The picked feature finder as a reusable instance
//! (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`, class
//! `FeatureFinderAlgorithmPicked`).
//!
//! [`FeatureFinderAlgorithmPicked`](crate::analysis::feature_finder_picked::instance::FeatureFinderAlgorithmPicked)
//! is the source object with its state: the `DefaultParamHandler` base and the
//! typed members of `updateMembers_`, the user seeds of `setSeeds`, the abort
//! counts `aborts_`, the debug abort map `abort_reasons_`, the debug stream
//! `log_` and the `ProgressLogger` base. Its
//! [`run`](crate::analysis::feature_finder_picked::instance::FeatureFinderAlgorithmPicked::run)
//! is source `run` followed by `run_`, in the source's order: the steps before
//! seed selection, then charge by charge the pattern scores, the seeds and
//! their extension, then the overlap resolution. It writes into the caller's
//! feature map as the source does.
//!
//! # Reusing an instance
//!
//! Source `run` stores a pointer to the caller's map (`setData_`,
//! `FeatureFinderAlgorithmPicked.cpp:134-138`) and never clears it unless the
//! input is empty (`:1059-1063`). Step 3.3 appends to it (`:844`), and step 4
//! then sorts, resolves, filters, re-sorts and annotates the *whole* map
//! (`:863-1011`). A caller's features therefore take part in the overlap
//! resolution, can be moved into a new feature's subordinates or have their
//! intensity set to zero and be removed, and get `spectrum_index` and
//! `spectrum_native_id` of the current input. Feature labels restart at 0 in
//! every run (`:437-438`). `aborts_` is only ever incremented (`:1135`) and
//! `abort_reasons_` only ever assigned (`:1138`), so both accumulate over the
//! runs of one instance, and the "reasons for not finalizing" block of every
//! run prints the accumulated counts. This port does all of that.
//!
//! See `docs/FEATURE_FINDER_PICKED_SUPPORT.md`, sections "Reusing an instance"
//! and "Debug mode".

use std::collections::BTreeMap;
use std::ops::Range;

use crate::analysis::feature_finder_picked::algorithm::{
    HANDLER_NAME, Limits, Options, PseudoRtShiftKey, RejectedParameters, Settings,
    check_parameters, check_run_conversions, default_parameters, validate_input,
};
use crate::analysis::feature_finder_picked::debug::{
    AbortReasons, DebugLogFile, DebugOutput, DebugTermination, FEATURE_DEBUG_PATH,
    FeatureDebugFiles, FeatureDebugInput, LogFragment, LogSink, NoLog, ReportLine, SeedMap,
    TerminationKind, TerminationPoint, abort_map_source, debug_experiment, g, g32, put_all,
    read_pseudo_rt_shift, seed_map, write_feature_debug_info,
};
use crate::analysis::feature_finder_picked::extension::{
    EMPTY_PATTERN_WHAT, OverallScores, extend_mass_traces_logged, find_best_isotope_fit_logged,
};
use crate::analysis::feature_finder_picked::fitting::{
    ABORT_COULD_NOT_EXTEND, ABORT_NO_ISOTOPE_PATTERN, FeatureInput, FittedModel, QualityOutcome,
    build_feature_checked, check_feature_quality_logged, crop_feature_logged,
};
use crate::analysis::feature_finder_picked::helper_structs::{MassTraces, NAN_RT_MERGE_WHAT, Seed};
use crate::analysis::feature_finder_picked::resolution::{
    annotate_apex, invalid_apex_warning, resolve_overlaps_logged,
};
use crate::analysis::feature_finder_picked::seeds::{IsotopeWindows, SeedStage, sort_user_seeds};
use crate::analysis::feature_finder_picked::source_sort::source_sort_by;
use crate::analysis::feature_finder_picked::trace_fitter::TraceFitterParams;
use crate::concept::parallel::{Threads, map_collect};
use crate::concept::progress_logger::{ProgressLogType, ProgressLogger};
use crate::kernel::{Feature, FeatureMap, MSExperiment, Point2D};
use crate::metadata::MetaValue;
use crate::param::{DefaultParamHandler, Param};
use crate::{Error, Result};

/// The `OPENMS_LOG_INFO` heading of the abort block, verbatim.
pub const ABORT_BLOCK_HEADING: &str =
    "Info: reasons for not finalizing a feature during its construction:";

/// The source's progress calls, sent to an optional progress logger.
///
/// Without a logger every call is a no-op, which is the source's default
/// `ProgressLogger` type `NONE` (`ProgressLogger.cpp:127-128`) without its
/// cost. The source passes some ranges whose begin exceeds their end (step 2
/// and step 3.2 on an input of fewer than `2 * min_spectra` scans); like the
/// Release build, [`ProgressLogger::start_progress`] accepts them, so they are
/// passed unchanged (`S 5 0` for a custom backend, as the Release build
/// passes it).
pub struct Progress<'a> {
    logger: Option<&'a mut ProgressLogger>,
}

impl Progress<'static> {
    /// Calls that go nowhere.
    pub(crate) fn silent() -> Self {
        Self { logger: None }
    }
}

impl<'a> Progress<'a> {
    /// Calls that go to `logger`, if any.
    pub(crate) fn new(logger: Option<&'a mut ProgressLogger>) -> Self {
        Self { logger }
    }

    /// `startProgress(begin, end, label)`.
    pub(crate) fn start(&mut self, begin: i64, end: i64, label: &str) -> Result<()> {
        match self.logger.as_deref_mut() {
            Some(logger) => logger.start_progress(begin, end, label),
            None => Ok(()),
        }
    }

    /// `setProgress(value)`.
    pub(crate) fn set(&mut self, value: i64) -> Result<()> {
        match self.logger.as_deref_mut() {
            Some(logger) => logger.set_progress(value),
            None => Ok(()),
        }
    }

    /// `setProgress` for every value of `values`, in order.
    pub(crate) fn set_each(&mut self, values: Range<i64>) -> Result<()> {
        if self.logger.is_some() {
            for value in values {
                self.set(value)?;
            }
        }
        Ok(())
    }

    /// `endProgress()`.
    pub(crate) fn end(&mut self) -> Result<()> {
        match self.logger.as_deref_mut() {
            Some(logger) => logger.end_progress(0),
            None => Ok(()),
        }
    }
}

/// A callback that sees each console line of a run as it arises.
pub type ConsoleSink = Box<dyn FnMut(&ReportLine) + Send>;

/// The run's console lines, kept in order and forwarded to an optional sink.
struct Console<'a> {
    lines: &'a mut Vec<ReportLine>,
    sink: Option<&'a mut ConsoleSink>,
}

impl Console<'_> {
    fn push(&mut self, line: ReportLine) {
        if let Some(sink) = self.sink.as_deref_mut() {
            sink(&line);
        }
        self.lines.push(line);
    }
}

/// The state of the source's `std::ofstream log_` member.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum LogState {
    /// Never opened: the next debug run opens (and truncates) the file.
    #[default]
    Closed,
    /// Opened by a debug run and still open, as the member is never closed:
    /// what that run wrote into the stream and how much of it the file buffer
    /// has handed to the file.
    Open { written: usize, flushed: usize },
    /// A later open failed: `failbit` is set and every write is dropped; the
    /// buffer still holds the opening run's unflushed tail.
    Failed { written: usize, flushed: usize },
}

impl LogState {
    /// The stream as [`DebugLogFile`], once a debug run opened it.
    fn file(self) -> Option<DebugLogFile> {
        match self {
            Self::Closed => None,
            Self::Open { written, flushed } => Some(DebugLogFile {
                written,
                flushed,
                failed: false,
            }),
            Self::Failed { written, flushed } => Some(DebugLogFile {
                written,
                flushed,
                failed: true,
            }),
        }
    }
}

/// The source algorithm object with its state.
///
/// Create it with [`Self::new`], configure it like the source through
/// [`Self::set_parameters`], [`Self::set_seeds`], [`Self::set_log_type`] and
/// [`Self::set_options`], and call [`Self::run`] as often as needed. After
/// every run, successful or not, [`Self::report`] holds the run's console
/// lines and [`Self::debug_output`] its debug output.
pub struct FeatureFinderAlgorithmPicked {
    handler: DefaultParamHandler,
    /// The refused set `getParameters()` shows after a refused
    /// `setParameters` ([`RejectedParameters::Shown`]).
    rejected: Option<Param>,
    settings: Settings,
    seeds: FeatureMap,
    aborts: BTreeMap<String, u32>,
    abort_reasons: AbortReasons,
    log_state: LogState,
    progress: Option<ProgressLogger>,
    options: Options,
    report: Vec<ReportLine>,
    debug: Option<DebugOutput>,
    /// Where the source process ends in the last run, if it does.
    termination: Option<DebugTermination>,
    console: Option<ConsoleSink>,
    /// Source member `isotope_distributions_`, which `run_` extends and never
    /// clears ([`IsotopeWindows::precalculate_onto`]).
    windows: Option<IsotopeWindows>,
}

impl std::fmt::Debug for FeatureFinderAlgorithmPicked {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FeatureFinderAlgorithmPicked")
            .field("parameters", self.handler.parameters())
            .field("settings", &self.settings)
            .field("seeds", &self.seeds.len())
            .field("aborts", &self.aborts)
            .field("abort_reasons", &self.abort_reasons.len())
            .field("log_state", &self.log_state)
            .field("log_type", &self.log_type())
            .field("options", &self.options)
            .finish_non_exhaustive()
    }
}

impl FeatureFinderAlgorithmPicked {
    /// The source default constructor: the 29 defaults
    /// ([`default_parameters`]) become the current parameters
    /// (`defaultsToParam_`) and the typed members are read from them
    /// (`updateMembers_`). No seeds, no aborts, progress type `NONE`, default
    /// [`Options`].
    ///
    /// # Errors
    ///
    /// Returns an error only if the fixed defaults cannot be built, which does
    /// not happen.
    pub fn new() -> Result<Self> {
        Self::with_options(Options::default())
    }

    /// [`Self::new`] with explicit native options.
    ///
    /// # Errors
    ///
    /// As [`Self::new`].
    pub fn with_options(options: Options) -> Result<Self> {
        let mut handler = DefaultParamHandler::new(HANDLER_NAME)?;
        let defaults = default_parameters()?;
        handler.set_defaults(defaults.clone())?;
        // Every default has a description, so `defaultsToParam_` warns about
        // nothing.
        let (settings, _) =
            handler.defaults_to_parameters_with(|merged| Settings::read(merged, &defaults))?;
        Ok(Self {
            handler,
            rejected: None,
            settings,
            seeds: FeatureMap::new(),
            aborts: BTreeMap::new(),
            abort_reasons: AbortReasons::new(),
            log_state: LogState::Closed,
            progress: None,
            options,
            report: Vec::new(),
            debug: None,
            termination: None,
            console: None,
            windows: None,
        })
    }

    /// `getName()`: `FeatureFinderAlgorithmPicked` unless [`Self::set_name`]
    /// changed it.
    pub fn name(&self) -> &str {
        self.handler.name()
    }

    /// `setName(name)` of the `DefaultParamHandler` base: the name the
    /// unknown-parameter warnings and the parameter errors of later
    /// [`Self::set_parameters`] calls and runs carry (`error_name_`,
    /// `DefaultParamHandler.cpp:65`, `Param.cpp:1085`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a name beyond the parameter text
    /// bounds; the name stays.
    pub fn set_name(&mut self, name: &str) -> Result<()> {
        self.handler.set_name(name)
    }

    /// `getSubsections()` (`DefaultParamHandler.cpp:120-123`): the subsections
    /// `setParameters` does not check.
    /// The source constructor registers none, so this is empty.
    pub fn subsections(&self) -> &[String] {
        self.handler.subsections()
    }

    /// `operator==` of the `DefaultParamHandler` base, which the source object
    /// inherits: the current parameters as [`Self::parameters`] shows them (a
    /// shown refused set included), the defaults, the subsections, the name
    /// and the two check flags, the parameter trees with the source's `Param`
    /// predicate (`DefaultParamHandler.cpp:32-40`). Like the source, it
    /// compares nothing else of the object.
    ///
    /// # Errors
    ///
    /// Returns an error only if a comparison exceeds the parameter work bound.
    pub fn handler_equal(&self, other: &Self) -> Result<bool> {
        Ok(self.handler.name() == other.handler.name()
            && self.handler.subsections() == other.handler.subsections()
            && self.handler.check_defaults() == other.handler.check_defaults()
            && self.handler.warn_empty_defaults() == other.handler.warn_empty_defaults()
            && self.parameters().source_equal(other.parameters())?
            && self
                .handler
                .defaults()
                .source_equal(other.handler.defaults())?)
    }

    /// `getDefaults()`.
    pub fn defaults(&self) -> &Param {
        self.handler.defaults()
    }

    /// `getDefaultParameters()`: a copy of the defaults.
    pub fn default_parameters(&self) -> Param {
        self.handler.defaults().clone()
    }

    /// `getParameters()`: the current parameters, which every [`Self::run`]
    /// replaces by its argument merged with the defaults.
    ///
    /// After a refused parameter set this is, as in the source, the refused set
    /// merged with the defaults, unless [`Options::rejected_parameters`] is
    /// [`RejectedParameters::Discarded`] ([`Self::set_parameters`]).
    pub fn parameters(&self) -> &Param {
        self.rejected
            .as_ref()
            .unwrap_or_else(|| self.handler.parameters())
    }

    /// The typed members `updateMembers_` set, and the values `run_` reads.
    pub fn settings(&self) -> &Settings {
        &self.settings
    }

    /// `setParameters(param)`: `param` merged with the defaults becomes the
    /// current parameter set, which is checked against the defaults, and the
    /// typed members are updated.
    ///
    /// Returns the warnings the source logs, one per parameter the defaults do
    /// not know (such an entry is kept), in the source's wording and in the
    /// order `Param::checkDefaults` visits the merged set.
    ///
    /// The checks are the source's ([`check_parameters`]): an integer is
    /// narrowed to its low 32 bits before its restriction is checked, and the
    /// typed members are converted as the Linux x86_64 Release build converts
    /// them ([`Settings`]), so a 64-bit value such as `intensity:bins = 2^32 +
    /// 10` is accepted and read as 10 (executed:
    /// `../oracle/ffap-complete-fix2`, case `bigint`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a value of the wrong type or outside
    /// its restriction, with the text of the source's
    /// `Exception::InvalidParameter`; the typed members keep their values. For
    /// a negative `mass_trace:max_missing` or `intensity:bins` whose low 32
    /// bits pass the restriction, `updateMembers_` throws
    /// `Exception::ConversionError` half way
    /// ([`NEGATIVE_UNSIGNED_WHAT`](crate::analysis::feature_finder_picked::algorithm::NEGATIVE_UNSIGNED_WHAT)): the members before that one have their new
    /// values, as in the source. The source has assigned the merged set to
    /// `param_` before it checks it, so [`Self::parameters`] then shows the
    /// refused set ([`RejectedParameters::Shown`], the default; executed:
    /// `params_after_failed_set.txt`), or the last accepted one under
    /// [`RejectedParameters::Discarded`]. The warnings the source logged before
    /// it threw are dropped here; [`Self::set_parameters_logged`] keeps them.
    pub fn set_parameters(&mut self, parameters: &Param) -> Result<Vec<String>> {
        let mut warnings = Vec::new();
        self.set_parameters_logged(parameters, &mut warnings)?;
        Ok(warnings)
    }

    /// [`Self::set_parameters`], appending the unknown-parameter warnings to
    /// `warnings` whether or not the set is refused: on a refusal, those of the
    /// entries `Param::checkDefaults` visits before the refused one, which the
    /// source logs before it throws (executed: `rejected_stderr.txt`).
    ///
    /// # Errors
    ///
    /// As [`Self::set_parameters`].
    pub fn set_parameters_logged(
        &mut self,
        parameters: &Param,
        warnings: &mut Vec<String>,
    ) -> Result<()> {
        let defaults = self.handler.defaults().clone();
        let name = self.handler.name().to_owned();
        // `param_ = tmp` with the defaults filled in, before any check: the
        // merged set, which the shared handler builds without its own checks.
        let check = self.handler.check_defaults();
        let mut merged_handler = self.handler.clone();
        merged_handler.set_check_defaults(false);
        merged_handler.set_parameters(parameters)?;
        merged_handler.set_check_defaults(check);
        let merged = merged_handler.parameters();
        // `checkDefaults` (when enabled), then `updateMembers_`.
        let mut unknown = Vec::new();
        let checked = if check {
            check_parameters(merged, &defaults, &name, &mut unknown)
        } else {
            Ok(())
        };
        warnings.extend(
            unknown
                .iter()
                .map(|key| format!("Warning: {name} received the unknown parameter '{key}'!")),
        );
        let mut members = self.settings.clone();
        let result = checked
            .and_then(|()| members.update_members(merged))
            .and_then(|()| members.read_run_values(merged, &defaults));
        match result {
            Ok(()) => {
                self.handler = merged_handler;
                self.settings = members;
                self.rejected = None;
                Ok(())
            }
            Err(error) => {
                // A conversion that throws half way through `updateMembers_`
                // leaves the members before it updated.
                self.settings = members;
                if self.options.rejected_parameters == RejectedParameters::Shown {
                    self.rejected = Some(merged.clone());
                }
                Err(error)
            }
        }
    }

    /// Source member `isotope_distributions_`: the isotope windows the last
    /// run precalculated, which the next run extends.
    pub fn isotope_windows(&self) -> Option<&IsotopeWindows> {
        self.windows.as_ref()
    }

    /// `setSeeds(seeds)`: the user-specified seeds of the next run. [`Self::run`]
    /// replaces them with its own argument, as the source's `run` does.
    pub fn set_seeds(&mut self, seeds: &FeatureMap) {
        self.seeds = seeds.clone();
    }

    /// The user-specified seeds (source member `seeds_`). A run sorts them by
    /// m/z in place, as the Release build's `seeds_.sortByMZ()` leaves them.
    pub fn seeds(&self) -> &FeatureMap {
        &self.seeds
    }

    /// Source member `aborts_`: how often each abort reason kept a seed from
    /// becoming a feature, over every run of this instance.
    ///
    /// The source counts in `UInt`, and so does this, wrapping at `2^32`. The
    /// source increments the counts from inside its parallel region without
    /// synchronisation, which is a data race with more than one thread; the
    /// port counts serially in seed order, so the counts are the source's
    /// single-thread counts at any thread count.
    pub fn aborts(&self) -> &BTreeMap<String, u32> {
        &self.aborts
    }

    /// Source member `abort_reasons_`, filled by debug runs only and never
    /// cleared.
    pub fn abort_reasons(&self) -> &AbortReasons {
        &self.abort_reasons
    }

    /// The native options.
    pub fn options(&self) -> &Options {
        &self.options
    }

    /// Replace the native options.
    pub fn set_options(&mut self, options: Options) {
        self.options = options;
    }

    /// `setLogType(type)` of the `ProgressLogger` base: `NONE` (the default)
    /// removes the logger, any other type installs a
    /// [`ProgressLogger`] of that type.
    pub fn set_log_type(&mut self, log_type: ProgressLogType) {
        if log_type == ProgressLogType::None {
            self.progress = None;
        } else {
            let mut logger = ProgressLogger::new();
            logger.set_log_type(log_type);
            self.progress = Some(logger);
        }
    }

    /// `getLogType()` of the `ProgressLogger` base.
    pub fn log_type(&self) -> ProgressLogType {
        self.progress
            .as_ref()
            .map_or(ProgressLogType::None, ProgressLogger::log_type)
    }

    /// Install a configured progress logger (for example one with a custom
    /// backend or clock), or remove it.
    pub fn set_progress_logger(&mut self, logger: Option<ProgressLogger>) {
        self.progress = logger;
    }

    /// The installed progress logger, if any.
    pub fn progress_logger_mut(&mut self) -> Option<&mut ProgressLogger> {
        self.progress.as_mut()
    }

    /// Forward every console line of later runs to `sink` as it arises, in
    /// addition to [`Self::report`]; a progress logger writing to the same
    /// place then interleaves with the lines as the source's output does.
    pub fn set_console(&mut self, sink: Option<ConsoleSink>) {
        self.console = sink;
    }

    /// The console lines of the last [`Self::run`], in source order, up to
    /// where it ended.
    pub fn report(&self) -> &[ReportLine] {
        &self.report
    }

    /// The debug output of the last [`Self::run`], when `write_debug` was set;
    /// kept after an error, up to where the run ended.
    pub fn debug_output(&self) -> Option<&DebugOutput> {
        self.debug.as_ref()
    }

    /// Take the debug output of the last run.
    pub fn take_debug_output(&mut self) -> Option<DebugOutput> {
        self.debug.take()
    }

    /// Where the source process ends in the last [`Self::run`], with or
    /// without `write_debug`: the run was refused there, and the executed
    /// process leaves `debug/log.txt` at
    /// [`DebugTermination::log_file_bytes`]. `None` when the source returns
    /// or throws to its caller, or the run failed on one of the port's own
    /// ceilings.
    pub fn termination(&self) -> Option<&DebugTermination> {
        self.termination.as_ref()
    }

    /// The instance's never-closed `debug/log.txt` stream (source member
    /// `log_`), once a debug run opened it: what that run wrote, how much of
    /// it is in the file while the instance lives, and whether a later debug
    /// run failed to reopen it. Dropping the instance writes the rest, as the
    /// source's destructor does; a process that ends first leaves the file at
    /// [`DebugLogFile::flushed`].
    pub fn debug_log_file(&self) -> Option<DebugLogFile> {
        self.log_state.file()
    }

    /// Find features: source `run(PeakMap&&, FeatureMap&, const Param&, const FeatureMap& seeds)`.
    ///
    /// `experiment` holds centroided MS1 spectra and is consumed. The features
    /// found are *added* to `features`, whose existing features take part in
    /// the overlap resolution and the final sort (module documentation).
    ///
    /// 1. An experiment without spectra clears `features` completely
    ///    (`FeatureMap::clear(true)`: features, identifications, processing,
    ///    metadata and unique id) and returns; the parameters are not touched.
    /// 2. The input checks of [`validate_input`]; their warnings are
    ///    [`ReportLine::Warn`].
    /// 3. [`Self::set_parameters`]; its warnings follow.
    /// 4. `seeds` become the user seeds.
    /// 5. Source `run_`: steps 1 to 2.5, then per charge the pattern scores,
    ///    the seeds and their extension in parallel on [`Options::threads`],
    ///    then step 4 and the apex annotation.
    ///
    /// # Errors
    ///
    /// Every error of the steps above, and in particular:
    ///
    /// - [`Error::Unsupported`] where the source process terminates:
    ///   `write_debug` with [`PseudoRtShiftKey::Source`] and no usable
    ///   `debug:pseudo_rt_shift`, at the first seed that reaches the fit
    ///   ([`DebugOutput::termination`] records it);
    /// - [`Error::InvalidValue`] where the source terminates without a debug
    ///   key issue: a feature m/z without an isotope window in step 3.3.5
    ///   ([`DebugOutput::termination`] records it in a debug run, after that
    ///   seed's debug files);
    /// - [`Error::InvalidValue`] where the source is undefined on a caller's
    ///   map: a feature of charge 0 in an overlapping pair of different
    ///   charges, and a stale abort-reason seed outside the current input;
    /// - [`Error::InvalidValue`] with the source's `what()` text where the
    ///   source throws a catchable exception: [`NEGATIVE_UNSIGNED_WHAT`](crate::analysis::feature_finder_picked::algorithm::NEGATIVE_UNSIGNED_WHAT) for a
    ///   negative `fit:max_iterations` at the start of `run_`, and
    ///   `std::length_error`'s text for more isotope windows than a vector can
    ///   hold in step 2.5
    ///   ([`LENGTH_ERROR_WHAT`](crate::analysis::feature_finder_picked::seeds::LENGTH_ERROR_WHAT)).
    ///
    /// The instance and `features` are left as the source leaves them at that
    /// point: features of the charges already processed are in `features`, and
    /// their aborts are counted. [`Self::report`] and [`Self::debug_output`]
    /// hold what the run produced up to there: a debug run that fails after
    /// the point where the source opens `debug/log.txt` (after the score
    /// arrays, before step 1) has opened it and written what the source wrote.
    /// Where the source process ends instead of returning,
    /// [`Self::termination`] records the point, whether or not the run writes
    /// debug output: a wrapping score-array count that writes out of bounds, a
    /// seed-loop crash, hang or exception, a charge-0 remainder in step 4, and
    /// a stale abort seed outside the input.
    pub fn run(
        &mut self,
        experiment: MSExperiment,
        features: &mut FeatureMap,
        parameters: &Param,
        seeds: &FeatureMap,
    ) -> Result<()> {
        let mut report = Vec::new();
        self.debug = None;
        self.termination = None;
        let mut sink = self.console.take();
        let result = {
            let mut console = Console {
                lines: &mut report,
                sink: sink.as_mut(),
            };
            self.run_inner(experiment, features, parameters, seeds, &mut console)
        };
        self.console = sink;
        self.report = report;
        // A run that opened the stream leaves it open with what it wrote; the
        // stream keeps those counts for every later run (`LogState`).
        if let Some(out) = self.debug.as_ref().filter(|out| out.log_opened) {
            self.log_state = LogState::Open {
                written: out.log.len(),
                flushed: out.log.flushed_bytes(),
            };
        }
        if let Some(termination) = self.termination.as_mut() {
            termination.log_file_bytes = self.log_state.file().map(|file| file.flushed);
            if let Some(out) = self.debug.as_mut() {
                out.termination = Some(termination.clone());
            }
        }
        result
    }

    fn run_inner(
        &mut self,
        mut experiment: MSExperiment,
        features: &mut FeatureMap,
        parameters: &Param,
        seeds: &FeatureMap,
        report: &mut Console<'_>,
    ) -> Result<()> {
        if experiment.spectra.is_empty() {
            features.clear(true);
            return Ok(());
        }
        let mut warnings = Vec::new();
        let checked = validate_input(&mut experiment, &mut warnings);
        for warning in warnings {
            report.push(ReportLine::Warn(warning));
        }
        if !checked? {
            return Ok(());
        }
        let mut warnings = Vec::new();
        let applied = self.set_parameters_logged(parameters, &mut warnings);
        for warning in warnings {
            report.push(ReportLine::Warn(warning));
        }
        applied?;
        self.seeds = seeds.clone();
        self.run_core(experiment, features, report)
    }

    /// Source `run_` (`FeatureFinderAlgorithmPicked.cpp:140-1055`).
    fn run_core(
        &mut self,
        experiment: MSExperiment,
        features: &mut FeatureMap,
        report: &mut Console<'_>,
    ) -> Result<()> {
        let settings = self.settings.clone();
        let options = self.options;
        let debug = settings.write_debug;
        let limits = options.limits;
        // `UInt max_iterations = param_.getValue(...)` (`:152`) throws before
        // anything else happens.
        check_run_conversions(self.handler.parameters())?;
        let mut progress = Progress::new(self.progress.as_mut());

        // `seeds_.sortByMZ()` (`FeatureFinderAlgorithmPicked.cpp:190`) sorts the
        // member in place: every run replaces `seeds_` with the caller's map
        // (`run_inner`, as the source's `run` calls `setSeeds`) and sorts it
        // here, and `seeds()` returns the sorted copy.
        sort_user_seeds(&mut self.seeds)?;

        // Steps 0 to 2.5.
        let mut prefix = LogFragment::new();
        let stage_log = Vec::new();
        let mut opened = false;
        let prepared = if debug {
            SeedStage::prepare(
                experiment,
                &self.seeds,
                settings.clone(),
                &options,
                stage_log,
                self.windows.as_ref(),
                &mut prefix,
                &mut progress,
                &mut opened,
            )
        } else {
            SeedStage::prepare(
                experiment,
                &self.seeds,
                settings.clone(),
                &options,
                stage_log,
                self.windows.as_ref(),
                &mut NoLog,
                &mut progress,
                &mut opened,
            )
        };
        let mut stage = match prepared {
            Ok(stage) => stage,
            Err(error) => {
                // The source has opened the stream and created debug/features
                // before step 1 (`:226-232`) and written its first line; a
                // failure in steps 1 to 2.5 (the step-2.5 `std::length_error`,
                // or one of the port's ceilings there) leaves both, and the
                // stream stays open for the next run.
                if debug && opened {
                    let mut out = Some(open_debug_log(&mut self.log_state));
                    append_log(&mut out, &prefix, &limits)?;
                    self.debug = out;
                } else if !opened
                    && settings.score_arrays_overrun(&limits)
                    && settings
                        .charge_count()
                        .is_err_and(|refusal| refusal.to_string() == error.to_string())
                {
                    // The wrapped score-array count: the source writes past
                    // the arrays before it opens the stream (`:196-221`).
                    self.termination = Some(DebugTermination {
                        point: TerminationPoint::ScoreArrays,
                        kind: TerminationKind::OutOfBounds,
                        exception: "SIGSEGV",
                        message: error.to_string(),
                        log_file_bytes: None,
                    });
                }
                return Err(error);
            }
        };
        // Step 2.5 has replaced the member.
        self.windows = Some(stage.windows().clone());
        // `debug_` and the `log_.open` of `:226-232`, which follow the score
        // arrays.
        let mut out = debug.then(|| open_debug_log(&mut self.log_state));
        let result = (|| -> Result<()> {
            append_log(&mut out, &prefix, &limits)?;

            let fitter_parameters = TraceFitterParams {
                max_iteration: i64::from(settings.max_iterations),
                weighted: false,
            };
            let mut plot_nr_global: i64 = -1;
            let mut feature_nr_global: i64 = 0;
            let mut seed_work = 0u64;
            loop {
                // Steps 3.1 and 3.2.
                let mut pattern_log = LogFragment::new();
                let selected = if debug {
                    stage.select_next_charge(&mut pattern_log, &mut progress)?
                } else {
                    stage.select_next_charge(&mut NoLog, &mut progress)?
                };
                if !selected {
                    break;
                }
                let charge_index = stage.charges().len() - 1;
                append_log(&mut out, &pattern_log, &limits)?;
                if let Some(out) = out.as_mut() {
                    let charge_seeds = &stage.charges()[charge_index];
                    out.seed_maps.push(SeedMap {
                        charge: charge_seeds.charge,
                        map: seed_map(
                            stage.experiment(),
                            stage.scores(),
                            charge_index,
                            &charge_seeds.seeds,
                        )?,
                    });
                    report.push(ReportLine::StoreSeedMap(out.seed_maps.len() - 1));
                }
                progress.end()?;
                stage.log_seed_count();
                if let Some(line) = stage.log().last() {
                    report.push(ReportLine::Out(line.clone()));
                }

                // Step 3.3.
                preflight_charge(&stage, charge_index, &limits, &mut seed_work)?;
                let charge_seeds = &stage.charges()[charge_index];
                let charge = charge_seeds.charge;
                let count = charge_seeds.seeds.len();
                progress.start(
                    0,
                    i64::try_from(count)
                        .map_err(|_| Error::InvalidValue("seed count overflow".into()))?,
                    &format!("Extending seeds for charge {charge}"),
                )?;
                let outcomes = extend_charge(
                    &stage,
                    charge_index,
                    &fitter_parameters,
                    options.threads,
                    debug,
                );
                let mut book = Bookkeeping {
                    aborts: &mut self.aborts,
                    abort_reasons: &mut self.abort_reasons,
                    out: &mut out,
                    termination: &mut self.termination,
                    features: &mut features.features,
                    plot_nr_global: &mut plot_nr_global,
                    feature_nr_global: &mut feature_nr_global,
                    progress: &mut progress,
                    limits: &limits,
                };
                let debug_key = DebugKey {
                    policy: options.pseudo_rt_shift,
                    parameters: self.handler.parameters(),
                };
                let candidates =
                    settle_charge(&stage, charge_index, outcomes, &mut book, &debug_key)?;
                progress.end()?;
                report.push(ReportLine::Out(format!(
                    "Found {candidates} feature candidates for charge {charge}."
                )));
            }

            // Step 4.
            let all = &mut features.features;
            let n = all.len() as u64;
            progress.start(
                0,
                n.wrapping_mul(n) as i64,
                "Resolving overlapping features",
            )?;
            let mut resolution_log = LogFragment::new();
            if debug {
                put_all(
                    &mut resolution_log,
                    &[
                        "Resolving intersecting features (",
                        &all.len().to_string(),
                        " candidates)\n",
                    ],
                );
            }
            // `sortByMZ()`: `std::sort` with `Feature::MZLess`.
            source_sort_by(all, |a, b| a.mz < b.mz)?;
            let mut trap = None;
            let removed = if debug {
                resolve_overlaps_logged(
                    all,
                    settings.max_feature_intersection,
                    &mut resolution_log,
                    // `size_t` to `SignedSize`: two's complement.
                    &mut |value| progress.set(value as i64),
                    &mut trap,
                )
            } else {
                resolve_overlaps_logged(
                    all,
                    settings.max_feature_intersection,
                    &mut NoLog,
                    &mut |value| progress.set(value as i64),
                    &mut trap,
                )
            };
            append_log(&mut out, &resolution_log, &limits)?;
            if let (Some((first, second)), Err(error)) = (trap, &removed) {
                // `f2.getCharge() % f1.getCharge()` traps (`:936`, `:945`):
                // the process dies with the pair's `Intersection` line in the
                // stream's buffer.
                self.termination = Some(DebugTermination {
                    point: TerminationPoint::OverlapResolution { first, second },
                    kind: TerminationKind::ArithmeticTrap,
                    exception: "SIGFPE",
                    message: error.to_string(),
                    log_file_bytes: None,
                });
            }
            let removed = removed?;
            report.push(ReportLine::Info(format!(
                "Removed {removed} overlapping features."
            )));
            all.retain(|feature| feature.intensity != 0.0);
            // `sortByIntensity(true)`: `std::sort` with the arguments of
            // `Feature::IntensityLess` swapped.
            source_sort_by(all, |a, b| b.intensity < a.intensity)?;
            progress.end()?;
            let invalid_apex = annotate_apex(all, stage.experiment())?;
            if invalid_apex > 0 {
                report.push(ReportLine::Warn(invalid_apex_warning(invalid_apex)));
            }
            report.push(ReportLine::Info(String::new()));
            report.push(ReportLine::Info(ABORT_BLOCK_HEADING.into()));
            for (reason, count) in &self.aborts {
                report.push(ReportLine::Info(format!(" - {reason}: {count} times")));
            }
            report.push(ReportLine::Info(String::new()));
            report.push(ReportLine::Info(format!("{} features found.", all.len())));
            Ok(())
        })();
        if let Err(error) = result {
            self.debug = out;
            return Err(error);
        }
        if let Some(mut debug_out) = out {
            // The abort map (`:1028-1045`) may refuse; what came before stays.
            let map = match abort_map_source(&self.abort_reasons, stage.experiment()) {
                Ok(map) => map,
                Err(error) => {
                    self.debug = Some(debug_out);
                    return Err(error);
                }
            };
            match map {
                Ok(map) => {
                    debug_out.abort_reasons = Some(map);
                    report.push(ReportLine::StoreAbortReasons);
                }
                Err(stale) => {
                    // `map_[it2->first.spectrum]` reads out of bounds: the
                    // process dies before it stores the map (executed:
                    // SIGSEGV).
                    self.termination = Some(DebugTermination {
                        point: TerminationPoint::AbortMap { entry: stale.entry },
                        kind: TerminationKind::OutOfBounds,
                        exception: "SIGSEGV",
                        message: stale.error.to_string(),
                        log_file_bytes: None,
                    });
                    self.debug = Some(debug_out);
                    return Err(stale.error);
                }
            }
            let scores = stage.scores().clone();
            let input = debug_experiment(stage.into_experiment(), &scores);
            match input {
                Ok(input) => {
                    debug_out.input = Some(input);
                    report.push(ReportLine::StoreInput);
                    self.debug = Some(debug_out);
                }
                Err(error) => {
                    self.debug = Some(debug_out);
                    return Err(error);
                }
            }
        }
        Ok(())
    }
}

/// `log_.open("debug/log.txt")` (`FeatureFinderAlgorithmPicked.cpp:231`) on
/// the instance's never-closed stream: the first open succeeds and truncates
/// the file, every later one fails and leaves the stream failed.
fn open_debug_log(state: &mut LogState) -> DebugOutput {
    let opened = match *state {
        LogState::Closed => {
            *state = LogState::Open {
                written: 0,
                flushed: 0,
            };
            true
        }
        LogState::Open { written, flushed } | LogState::Failed { written, flushed } => {
            *state = LogState::Failed { written, flushed };
            false
        }
    };
    DebugOutput {
        log_opened: opened,
        ..DebugOutput::default()
    }
}

/// Append a fragment to the run's debug log, unless the stream is not open.
fn append_log(
    out: &mut Option<DebugOutput>,
    fragment: &LogFragment,
    limits: &Limits,
) -> Result<()> {
    if let Some(out) = out.as_mut() {
        if out.log_opened {
            out.log.append(fragment, limits.max_debug_bytes)?;
        }
    }
    Ok(())
}

/// The ceilings of one charge's seed loop, checked before it starts: its seed
/// count, and the seed work accumulated over the charges so far.
pub(crate) fn preflight_charge(
    stage: &SeedStage,
    charge_index: usize,
    limits: &Limits,
    work: &mut u64,
) -> Result<()> {
    let spectra = stage.experiment().spectra.len() as u64;
    let isotopes = stage.settings().max_isotopes() as u64;
    let per_seed = isotopes.saturating_mul(isotopes.saturating_add(spectra));
    let charge = &stage.charges()[charge_index];
    if charge.seeds.len() > limits.max_seeds {
        return Err(Error::InvalidValue(format!(
            "{} seeds for charge {} exceed the limit of {}",
            charge.seeds.len(),
            charge.charge,
            limits.max_seeds
        )));
    }
    *work = work.saturating_add((charge.seeds.len() as u64).saturating_mul(per_seed));
    if *work > limits.max_seed_work {
        return Err(Error::InvalidValue(format!(
            "the seed loop may take {} work units, exceeding the limit of {}",
            *work, limits.max_seed_work
        )));
    }
    Ok(())
}

/// What `writeFeatureDebugInfo_` needs of one seed, kept until its `plot_nr`
/// is known.
pub(crate) struct DebugWrite {
    model: FittedModel,
    traces: MassTraces,
    new_traces: MassTraces,
    feature_ok: bool,
    error_msg: String,
    final_score: f64,
    seed_mz: f64,
}

/// What one seed produced in step 3.3.
pub(crate) struct SeedOutcome {
    /// Whether the seed reached the fit and therefore consumed a `plot_nr`.
    plot_nr_used: bool,
    /// The candidate, or the source abort reason that dropped the seed.
    result: std::result::Result<SeedCandidate, String>,
    /// The seed's debug lines up to the debug write (or its abort).
    log: Option<LogFragment>,
    /// The debug write, when the seed reached `:714` in a debug run.
    debug_write: Option<DebugWrite>,
    /// Where the source's process ends with this seed.
    terminated: Option<SeedTermination>,
}

/// How the source's process ends at one seed.
pub(crate) enum SeedTermination {
    /// After the debug write, step 3.3.5 throws `Exception::InvalidValue`
    /// inside the OpenMP region (`:790`); the `what()` text.
    Exception(String),
    /// The port refuses where the executed process dies (an out-of-bounds
    /// read) or never returns (the NaN profile merge), after the seed's log
    /// lines so far.
    Refused(Error, TerminationKind),
}

/// One accepted candidate and the later seeds it swallows.
pub(crate) struct SeedCandidate {
    feature: Feature,
    /// Indices of the seeds after this one that lie inside the feature: the
    /// source's `seeds_in_features[i]`.
    contained: Vec<usize>,
}

/// The parallel part of step 3.3 for one charge: every seed's outcome, in seed
/// order, on `threads` workers.
pub(crate) fn extend_charge(
    stage: &SeedStage,
    charge_index: usize,
    fitter_parameters: &TraceFitterParams,
    threads: Threads,
    debug: bool,
) -> Vec<Result<SeedOutcome>> {
    let charge_seeds = &stage.charges()[charge_index];
    let charge = charge_seeds.charge;
    let seeds = &charge_seeds.seeds;
    let indices: Vec<usize> = (0..seeds.len()).collect();
    let overall = OverallScores::new(stage.scores(), charge_index);
    map_collect(&indices, threads, |&index| {
        if debug {
            let mut log = LogFragment::new();
            let outcome = extend_seed(
                stage,
                overall,
                fitter_parameters,
                charge,
                seeds,
                index,
                &mut log,
            );
            outcome.map(|mut outcome| {
                outcome.log = Some(log);
                outcome
            })
        } else {
            extend_seed(
                stage,
                overall,
                fitter_parameters,
                charge,
                seeds,
                index,
                &mut NoLog,
            )
        }
    })
}

/// One seed of step 3.3: isotope fit, extension, fit, cropping, quality checks
/// and feature creation (`FeatureFinderAlgorithmPicked.cpp:596-817`).
fn extend_seed<L: LogSink>(
    stage: &SeedStage,
    overall: OverallScores<'_>,
    fitter_parameters: &TraceFitterParams,
    charge: i32,
    seeds: &[Seed],
    index: usize,
    log: &mut L,
) -> Result<SeedOutcome> {
    let settings = stage.settings();
    let spectra = &stage.experiment().spectra;
    let seed = seeds[index];
    let aborted = |plot_nr_used: bool, reason: &str| SeedOutcome {
        plot_nr_used,
        result: Err(reason.to_owned()),
        log: None,
        debug_write: None,
        terminated: None,
    };
    let seed_spectrum = &spectra[seed.spectrum];
    let seed_peak = seed_spectrum.peaks[seed.peak];
    if log.enabled() {
        // The source writes these lines on the master thread only; with one
        // thread that is every seed.
        put_all(log, &["\n", "Seed ", &index.to_string(), ":\n"]);
        put_all(log, &[" - Int: ", &g32(seed_peak.intensity), "\n"]);
        put_all(log, &[" - RT: ", &g(seed_spectrum.rt), "\n"]);
        put_all(log, &[" - MZ: ", &g(seed_peak.mz), "\n"]);
    }

    // Refusals where the executed process dies or never returns keep the
    // seed's log lines so far, as the process had written them to its buffer.
    let refused = |plot_nr_used: bool, error: Error, kind: TerminationKind| SeedOutcome {
        plot_nr_used,
        result: Err(String::new()),
        log: None,
        debug_write: None,
        terminated: Some(SeedTermination::Refused(error, kind)),
    };

    let (isotope_fit_quality, pattern) =
        find_best_isotope_fit_logged(spectra, stage.windows(), settings, seed, charge, log)?;
    if isotope_fit_quality < settings.min_isotope_fit {
        return Ok(aborted(false, ABORT_NO_ISOTOPE_PATTERN));
    }
    let mut traces = match extend_mass_traces_logged(spectra, overall, settings, &pattern, log) {
        Ok(traces) => traces,
        // `extendMassTraces_` reads the first entry of an empty best pattern
        // (`feature:min_isotope_fit` 0 lets a seed without a placement through):
        // executed SIGSEGV.
        Err(Error::InvalidValue(what)) if what == EMPTY_PATTERN_WHAT => {
            return Ok(refused(
                false,
                Error::InvalidValue(what),
                TerminationKind::OutOfBounds,
            ));
        }
        Err(error) => return Err(error),
    };
    let seed_mz = seed_peak.mz;
    if !traces.is_valid(seed_mz, settings.trace_tolerance) {
        return Ok(aborted(false, ABORT_COULD_NOT_EXTEND));
    }

    // Source: the baseline estimate is three quarters of the lowest peak.
    traces.update_baseline();
    traces.baseline *= 0.75;
    traces
        .get_mut(traces.max_trace)
        .ok_or_else(|| {
            Error::InvalidValue(
                "FeatureFinderAlgorithmPicked seed extension: the maximum trace index is out of \
                 range; the source dereferences it here"
                    .into(),
            )
        })?
        .update_maximum();

    let mut model = FittedModel::new(settings.rt_shape, *fitter_parameters);
    // The source's fit can throw `Exception::UnableToFit` (`TraceFitter.cpp:111`,
    // `:129`), which would escape its parallel region and end the process, but
    // no input reaches either throw from here (see `FittedModel::fit`). An error
    // here is therefore one of the port's own ceilings or the refused merge of a
    // NaN retention time into the intensity profile, where the source never
    // returns; the run fails with it rather than turning it into an abort
    // reason the source never records.
    match model.fit(&traces) {
        Ok(()) => {}
        // `computeIntensityProfile` never ends (`CPP-242`).
        Err(Error::InvalidValue(what)) if what == NAN_RT_MERGE_WHAT => {
            return Ok(refused(
                true,
                Error::InvalidValue(what),
                TerminationKind::NeverReturns,
            ));
        }
        Err(error) => return Err(error),
    }
    let new_traces =
        crop_feature_logged(model.as_fitter(), &traces, settings.min_trace_score, log)?;
    let outcome =
        check_feature_quality_logged(model.as_fitter(), &new_traces, seed_mz, settings, log)?;
    let debug_write = |model: &FittedModel, outcome: &QualityOutcome| DebugWrite {
        model: model.clone(),
        traces: traces.clone(),
        new_traces: new_traces.clone(),
        feature_ok: matches!(outcome, QualityOutcome::Accepted(_)),
        error_msg: match outcome {
            QualityOutcome::Rejected(reason) => (*reason).to_owned(),
            QualityOutcome::Accepted(_) => String::new(),
        },
        final_score: match outcome {
            QualityOutcome::Accepted(quality) => quality.final_score,
            QualityOutcome::Rejected(_) => 0.0,
        },
        seed_mz,
    };
    let write = log.enabled().then(|| debug_write(&model, &outcome));
    let quality = match outcome {
        QualityOutcome::Rejected(reason) => {
            let mut rejected = aborted(true, reason);
            rejected.debug_write = write;
            return Ok(rejected);
        }
        QualityOutcome::Accepted(quality) => quality,
    };
    let traces = new_traces;
    let feature = match build_feature_checked(FeatureInput {
        model: &model,
        traces: &traces,
        pattern: &pattern,
        windows: stage.windows(),
        settings,
        charge,
        // Overwritten serially; see `settle_charge`.
        plot_nr: -1,
        quality,
    })? {
        Ok(feature) => feature,
        Err(what) => {
            // The source has written this seed's log lines and debug files
            // (`:717`) before `:790` throws.
            return Ok(SeedOutcome {
                plot_nr_used: true,
                result: Err(String::new()),
                log: None,
                debug_write: write,
                terminated: Some(SeedTermination::Exception(what)),
            });
        }
    };

    // Source: every later seed inside both the overall bounding box and one of
    // the mass-trace hulls.
    let mut contained = Vec::new();
    if let Some(bounds) = feature.hull_bounding_box() {
        for (offset, later) in seeds.iter().enumerate().skip(index + 1) {
            let rt = spectra[later.spectrum].rt;
            let mz = spectra[later.spectrum].peaks[later.peak].mz;
            if bounds.encloses(Point2D::new(rt, mz))? && feature.encloses(rt, mz)? {
                contained.push(offset);
            }
        }
    }
    Ok(SeedOutcome {
        plot_nr_used: true,
        result: Ok(SeedCandidate { feature, contained }),
        log: None,
        debug_write: write,
        terminated: None,
    })
}

/// Where the serial part of step 3.3 records its results.
pub(crate) struct Bookkeeping<'a, 'p> {
    pub(crate) aborts: &'a mut BTreeMap<String, u32>,
    pub(crate) abort_reasons: &'a mut AbortReasons,
    pub(crate) out: &'a mut Option<DebugOutput>,
    /// Where the source process ends, debug run or not; the run completes
    /// the stream length afterwards.
    pub(crate) termination: &'a mut Option<DebugTermination>,
    pub(crate) features: &'a mut Vec<Feature>,
    pub(crate) plot_nr_global: &'a mut i64,
    pub(crate) feature_nr_global: &'a mut i64,
    pub(crate) progress: &'a mut Progress<'p>,
    pub(crate) limits: &'a Limits,
}

/// Which parameter `writeFeatureDebugInfo_` reads, and from where.
pub(crate) struct DebugKey<'a> {
    pub(crate) policy: PseudoRtShiftKey,
    pub(crate) parameters: &'a Param,
}

/// The serial part of step 3.3 for one charge, in seed order: progress, the
/// debug lines and files, the aborts, and the candidates that survive the
/// containment pass, appended to the output map with fresh labels. Returns
/// the number of candidates.
pub(crate) fn settle_charge(
    stage: &SeedStage,
    charge_index: usize,
    outcomes: Vec<Result<SeedOutcome>>,
    book: &mut Bookkeeping<'_, '_>,
    key: &DebugKey<'_>,
) -> Result<usize> {
    let charge_seeds = &stage.charges()[charge_index];
    let charge = charge_seeds.charge;
    let seeds = &charge_seeds.seeds;
    let features_before = book.features.len();
    let mut accepted: Vec<(usize, SeedCandidate)> = Vec::new();
    for (index, outcome) in outcomes.into_iter().enumerate() {
        // `IF_MASTERTHREAD setProgress(gl_progress++)`: every seed with one
        // thread.
        book.progress.set(index as i64)?;
        let outcome = outcome?;
        let plot_nr = if outcome.plot_nr_used {
            *book.plot_nr_global += 1;
            *book.plot_nr_global
        } else {
            -1
        };
        if let Some(fragment) = &outcome.log {
            append_log(book.out, fragment, book.limits)?;
        }
        if let (Some(out), Some(write)) = (book.out.as_mut(), &outcome.debug_write) {
            let key_name = match key.policy {
                PseudoRtShiftKey::Source => "debug:pseudo_rt_shift",
                PseudoRtShiftKey::Declared => "advanced:pseudo_rt_shift",
            };
            match read_pseudo_rt_shift(key.parameters, key_name)? {
                Ok(shift) => {
                    // A heap-address shift may refuse here; the files of the
                    // earlier seeds stay in the output.
                    let files: FeatureDebugFiles = write_feature_debug_info(FeatureDebugInput {
                        fitter: write.model.as_fitter(),
                        traces: &write.traces,
                        new_traces: &write.new_traces,
                        feature_ok: write.feature_ok,
                        error_msg: &write.error_msg,
                        final_score: write.final_score,
                        plot_nr,
                        seed_mz: write.seed_mz,
                        features_len: features_before,
                        pseudo_rt_shift: shift,
                        path: FEATURE_DEBUG_PATH,
                    })?;
                    out.feature_files.push(files);
                }
                Err((exception, message)) => {
                    *book.termination = Some(DebugTermination {
                        point: TerminationPoint::Seed {
                            charge,
                            seed_index: index,
                            plot_nr,
                        },
                        kind: TerminationKind::Exception,
                        exception,
                        message: message.clone(),
                        log_file_bytes: None,
                    });
                    return Err(Error::Unsupported(format!(
                        "write_debug: writeFeatureDebugInfo_ reads '{key_name}', and the source \
                         throws {exception} ({message}) inside its OpenMP region for seed \
                         {index} of charge {charge}, which terminates the process; select \
                         PseudoRtShiftKey::Declared to read advanced:pseudo_rt_shift instead"
                    )));
                }
            }
        }
        match outcome.terminated {
            None => {}
            Some(SeedTermination::Exception(what)) => {
                // `getIsotopeDistribution_` throws at `:790`, inside the OpenMP
                // region: the process ends after this seed's debug files.
                *book.termination = Some(DebugTermination {
                    point: TerminationPoint::Seed {
                        charge,
                        seed_index: index,
                        plot_nr,
                    },
                    kind: TerminationKind::Exception,
                    exception: "InvalidValue",
                    message: what.clone(),
                    log_file_bytes: None,
                });
                return Err(Error::InvalidValue(format!(
                    "FeatureFinderAlgorithmPicked step 3.3.5: {what}; the source throws this \
                     inside its OpenMP region for seed {index} of charge {charge}, where \
                     std::terminate ends the process"
                )));
            }
            Some(SeedTermination::Refused(error, kind)) => {
                // The executed process dies here (or never returns): its log
                // file holds what its buffer had flushed, this seed's lines
                // included in the buffer.
                *book.termination = Some(DebugTermination {
                    point: TerminationPoint::Seed {
                        charge,
                        seed_index: index,
                        plot_nr,
                    },
                    kind,
                    exception: match kind {
                        TerminationKind::OutOfBounds => "SIGSEGV",
                        TerminationKind::ArithmeticTrap => "SIGFPE",
                        TerminationKind::Exception | TerminationKind::NeverReturns => "",
                    },
                    message: error.to_string(),
                    log_file_bytes: None,
                });
                return Err(error);
            }
        }
        match outcome.result {
            Err(reason) => {
                // `abort_` (`:1129-1140`).
                if let Some(out) = book.out.as_mut() {
                    if out.log_opened {
                        let mut line = LogFragment::new();
                        put_all(&mut line, &["Abort: ", &reason, "\n"]);
                        out.log.append(&line, book.limits.max_debug_bytes)?;
                    }
                }
                let count = book.aborts.entry(reason.clone()).or_insert(0);
                *count = count.wrapping_add(1);
                if book.out.is_some() {
                    book.abort_reasons.insert(seeds[index], &reason)?;
                }
            }
            Ok(mut candidate) => {
                // The source assigns `plot_nr` inside a critical section, so
                // its value depends on the schedule with several threads; with
                // one it is this seed-order number, which the survivors'
                // labels overwrite below.
                candidate
                    .feature
                    .metadata
                    .insert("label".into(), MetaValue::from(plot_nr));
                accepted.push((index, candidate));
            }
        }
    }

    let mut contained_seeds: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    let mut feature_candidates = 0usize;
    for (seed_nr, candidate) in accepted {
        if contained_seeds.contains(&seed_nr) {
            continue;
        }
        feature_candidates += 1;
        let mut feature = candidate.feature;
        feature
            .metadata
            .insert("label".into(), MetaValue::from(*book.feature_nr_global));
        *book.feature_nr_global += 1;
        book.features
            .try_reserve(1)
            .map_err(|_| Error::InvalidValue("cannot allocate a feature".into()))?;
        book.features.push(feature);
        contained_seeds.extend(candidate.contained);
    }
    Ok(feature_candidates)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::feature_finder_picked::algorithm::RtShape;
    use crate::analysis::feature_finder_picked::debug::LogSink;
    use crate::kernel::{MSSpectrum, Peak1D};

    /// Source review of `.cpp:710-790`: a seed whose feature creation throws
    /// (a feature m/z without an isotope window, `getIsotopeDistribution_` at
    /// `:790`, inside the OpenMP region) has already written its log lines and
    /// its `writeFeatureDebugInfo_` files (`:717`) when the process ends. The
    /// serial settlement keeps both, records the termination and fails the
    /// run at that seed; no abort is counted for it. No executed input reaches
    /// this path (the executed and generated inputs of the oracles all
    /// returned), so the outcome is built here.
    #[test]
    fn a_step_3_3_5_termination_keeps_the_seed_debug_output() {
        let spectra = (0..3)
            .map(|index| MSSpectrum {
                rt: f64::from(index),
                ms_level: 1,
                peaks: vec![Peak1D::new(500.0, 10.0)],
                ..MSSpectrum::default()
            })
            .collect();
        let experiment = MSExperiment {
            spectra,
            ..MSExperiment::default()
        };
        let defaults = default_parameters().unwrap();
        let (settings, _) = Settings::from_parameters(&Param::new()).unwrap();
        let stage = SeedStage::compute(
            experiment,
            &FeatureMap::new(),
            settings,
            &Options::default(),
            Vec::new(),
        )
        .unwrap();
        assert_eq!(stage.charges()[0].charge, 1);

        let mut fragment = LogFragment::new();
        fragment.put("\nSeed 0:\n");
        let write = DebugWrite {
            model: FittedModel::new(
                RtShape::Symmetric,
                TraceFitterParams {
                    max_iteration: 500,
                    weighted: false,
                },
            ),
            traces: MassTraces::default(),
            new_traces: MassTraces::default(),
            feature_ok: true,
            error_msg: String::new(),
            final_score: 0.5,
            seed_mz: 500.0,
        };
        let what = "the value '9223372036854775808' was used but is not valid; \
                    IsotopeDistribution not precalculated. Maximum allowed index is 43";
        let outcome = SeedOutcome {
            plot_nr_used: true,
            result: Err(String::new()),
            log: Some(fragment),
            debug_write: Some(write),
            terminated: Some(SeedTermination::Exception(what.to_owned())),
        };
        let mut aborts = BTreeMap::new();
        let mut abort_reasons = AbortReasons::new();
        let mut out = Some(DebugOutput {
            log_opened: true,
            ..DebugOutput::default()
        });
        let mut termination = None;
        let mut features = Vec::new();
        let (mut plot_nr_global, mut feature_nr_global) = (-1, 0);
        let mut progress = Progress::silent();
        let limits = Limits::default();
        let mut book = Bookkeeping {
            aborts: &mut aborts,
            abort_reasons: &mut abort_reasons,
            out: &mut out,
            termination: &mut termination,
            features: &mut features,
            plot_nr_global: &mut plot_nr_global,
            feature_nr_global: &mut feature_nr_global,
            progress: &mut progress,
            limits: &limits,
        };
        let key = DebugKey {
            policy: PseudoRtShiftKey::Declared,
            parameters: &defaults,
        };
        let error = settle_charge(&stage, 0, vec![Ok(outcome)], &mut book, &key).unwrap_err();
        assert!(error.to_string().contains(what), "{error}");
        assert!(aborts.is_empty() && abort_reasons.is_empty() && features.is_empty());
        assert_eq!(plot_nr_global, 0);
        let out = out.unwrap();
        assert_eq!(out.log.text(), "\nSeed 0:\n");
        assert_eq!(out.feature_files.len(), 1);
        assert_eq!(out.feature_files[0].dta_name(), "debug/features/0.dta");
        // The run completes the stream length (`FeatureFinderAlgorithmPicked::run`).
        assert_eq!(out.termination, None);
        assert_eq!(
            termination,
            Some(DebugTermination {
                point: TerminationPoint::Seed {
                    charge: 1,
                    seed_index: 0,
                    plot_nr: 0,
                },
                kind: TerminationKind::Exception,
                exception: "InvalidValue",
                message: what.to_owned(),
                log_file_bytes: None,
            })
        );
    }
}

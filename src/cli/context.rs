// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Resolved tool parameters and run services: the native form of the source
//! `TOPPBase` `get*_` accessors, `parseRange_`, `inputFileReadable_`,
//! `outputFileWritable_`, `setMaxNumberOfThreads`, the progress log type and
//! the test-mode unique-id seed.
//!
//! See `docs/TOPP_CLI_SUPPORT.md` for the supported source subset.

use super::logging::ToolLog;
use super::parameter::ExitCode;
use super::processing::{self, AddDataProcessing};
use crate::concept::UniqueIdGenerator;
use crate::concept::parallel::Threads;
use crate::concept::progress_logger::ProgressLogType;
use crate::metadata::{DataProcessing, ProcessingAction};
use crate::param::{Param, ParamValue};
use crate::system::file;
use crate::{Error, Result};
use std::io::Write;
use std::sync::Arc;

/// Seed of the unique-id generator under `-test`, as `TOPPBase::main` sets it
/// (`TOPPBase.cpp:369-376`). Its first two raw draws are
/// 5233264595117471314 and 4835329514588776807.
pub const TEST_MODE_UNIQUE_ID_SEED: u64 = 19_991_231_235_959;

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// Parameters and services of one tool run, after defaults, INI file and
/// command line are merged and validated.
///
/// The source exposes these through protected `getStringOption_` style methods
/// and members on the tool itself. Here they are a borrowed context passed to
/// the tool body, so a tool cannot mutate its own resolved parameters mid-run,
/// and every service a later tool needs — progress logging, the thread policy,
/// the unique-id generator and processing annotations — derives from the same
/// resolved values.
#[derive(Clone, Debug)]
pub struct ToolContext {
    tool_name: String,
    version: String,
    ini_location: String,
    param: Param,
    debug_level: i64,
    test_mode: bool,
    no_progress: bool,
    force: bool,
    threads: i64,
    log: Arc<ToolLog>,
}

impl ToolContext {
    /// A context over `param`, the resolved tree without the `<tool>:<n>:`
    /// prefix, reading and writing the run's log.
    pub(crate) fn new(
        tool_name: &str,
        version: &str,
        ini_location: &str,
        param: Param,
        log: Arc<ToolLog>,
    ) -> Self {
        let int = |key: &str, fallback: i64| match param.value(key) {
            Ok(ParamValue::Integer(value)) => *value,
            _ => fallback,
        };
        let flag =
            |key: &str| matches!(param.value(key), Ok(ParamValue::String(text)) if text == "true");
        let debug_level = int("debug", 0);
        let threads = int("threads", 1);
        let test_mode = flag("test");
        let no_progress = flag("no_progress");
        let force = flag("force");
        Self {
            tool_name: tool_name.to_owned(),
            version: version.to_owned(),
            ini_location: ini_location.to_owned(),
            param,
            debug_level,
            test_mode,
            no_progress,
            force,
            threads,
            log,
        }
    }

    /// The executable and INI section name, as `toolName_`.
    pub fn tool_name(&self) -> &str {
        &self.tool_name
    }
    /// The product version the tool reports, as the source `version_`.
    pub fn version(&self) -> &str {
        &self.version
    }
    /// The INI section this run reads, as `getIniLocation_`: `<tool>:1:` for
    /// every run that reaches the tool body, because `-instance` is rejected
    /// like in the source, which leaves it out of the defaults.
    pub fn ini_location(&self) -> &str {
        &self.ini_location
    }
    /// The complete resolved parameter tree, as `getParam_`.
    ///
    /// Keys carry no `<tool>:1:` prefix, so a subsection parameter reads as
    /// `algorithm:peakcount`, exactly as in the source.
    pub fn param(&self) -> &Param {
        &self.param
    }
    /// Source `-debug`, the debug level; zero disables debug output.
    pub fn debug_level(&self) -> i64 {
        self.debug_level
    }
    /// Source `-test`: outputs must not depend on the clock, the machine or
    /// absolute paths.
    pub fn test_mode(&self) -> bool {
        self.test_mode
    }
    /// Source `-no_progress`: progress logging is disabled.
    pub fn no_progress(&self) -> bool {
        self.no_progress
    }
    /// Source `-force`: the tool may override its own safety checks.
    pub fn force(&self) -> bool {
        self.force
    }
    /// Source `-threads` as resolved from the defaults, the INI file and the
    /// command line; zero or a negative count requests every available
    /// processor. [`thread_policy`](Self::thread_policy) turns it into a worker
    /// count.
    pub fn threads(&self) -> i64 {
        self.threads
    }

    /// The progress logger type the tool hands to loaders and algorithms.
    ///
    /// Source `log_type_` (`TOPPBase.cpp:400-403`): `CMD` unless `-no_progress`
    /// is given, when it stays `NONE`.
    pub fn progress_log_type(&self) -> ProgressLogType {
        if self.no_progress {
            ProgressLogType::None
        } else {
            ProgressLogType::Cmd
        }
    }
    /// The worker-thread policy of this run: how many workers the tool body
    /// and the parallel algorithms it calls may use.
    ///
    /// Source `setMaxNumberOfThreads(getParamAsInt_("threads", 1))`
    /// (`TOPPBase.cpp:84-98, 408` at cli c19e494), which runs after the INI
    /// file and the command line are merged, so `-threads` and an INI
    /// `threads` value act alike:
    ///
    /// * a positive count `n` gives `n` workers;
    /// * zero **or a negative count** gives every available processor, as the
    ///   source's `if (num_threads <= 0) num_threads = omp_get_num_procs()`.
    ///   [`Threads::from_cli`] maps a negative count to one worker; this
    ///   method does not use it for non-positive counts, because the source
    ///   does not. Executed on the C++ Release build: `-threads -1` and
    ///   `-threads -7` start the same 128-thread team as `-threads 0` on a
    ///   128-processor node (`oracle/tool-threads/cpp_results_ibminode06.jsonl`).
    ///
    /// "Available" is [`std::thread::available_parallelism`]. Like
    /// `omp_get_num_procs`, it counts the processors in the affinity mask
    /// (`taskset -c 0-3` gives 4 in both); on Linux it also honours a cgroup
    /// CPU quota, which `omp_get_num_procs` ignores, so under a quota the port
    /// may start fewer workers than the source. Results do not change, by the
    /// determinism contract.
    ///
    /// `OMP_NUM_THREADS` does not enter the policy, as it does not enter the
    /// source's: `omp_set_num_threads` overrides it for every parallel region
    /// of the tool body. Executed: `-threads 4` gives the C++ tool body a
    /// four-thread team under `OMP_NUM_THREADS=1` and `=16` alike. The extra
    /// threads the C++ tools start without `OMP_NUM_THREADS` (129 at
    /// `-threads 1`, also for `-write_ini`) belong to the OpenBLAS library
    /// linked into that build, which sizes its server pool from the variable
    /// at library load (`blas_thread_init`); they do no tool work, and this
    /// port links no BLAS. `RAYON_NUM_THREADS` is ignored too, because the
    /// pool size is always explicit.
    ///
    /// A count above both [`THREAD_CEILING`](Self::THREAD_CEILING) and the
    /// available processors is clamped to the larger of the two, so a typing
    /// slip such as `-threads 2147483647` cannot ask the operating system for
    /// two billion threads; libgomp would try and abort. Clamping never changes
    /// a result, by the determinism contract.
    ///
    /// The source sets a process-wide OpenMP limit. Here the policy is passed
    /// to each computation, and [`in_thread_pool`](Self::in_thread_pool) runs
    /// the tool body on a pool of exactly this size, which keeps the parallel
    /// result bit-identical to the serial one (`src/concept/parallel.rs`). The
    /// source `@note` that the setting only works when OpenMS is compiled with
    /// OpenMP carries over: without the `parallel` feature every computation is
    /// serial, whatever this returns.
    pub fn thread_policy(&self) -> Threads {
        let Ok(requested) = usize::try_from(self.threads) else {
            return Threads::all();
        };
        if requested == 0 {
            return Threads::all();
        }
        if requested <= Self::THREAD_CEILING {
            return Threads::from_cli(self.threads);
        }
        let ceiling = Self::THREAD_CEILING.max(Threads::all().get());
        Threads::from_cli(i64::try_from(requested.min(ceiling)).unwrap_or(1))
    }

    /// Worker count above which [`thread_policy`](Self::thread_policy) clamps
    /// a request, unless the machine reports more available processors.
    ///
    /// Native bound; the source has none. 1024 exceeds the logical processors
    /// of the benchmark nodes (128) and of the largest gate host (384), and it
    /// never reduces a request on a machine with more processors, because the
    /// ceiling then rises to the available count.
    pub const THREAD_CEILING: usize = 1024;

    /// Stack size of each pool worker, in bytes: 8 MiB.
    ///
    /// [`in_thread_pool`](Self::in_thread_pool) moves the tool body off the
    /// main thread, whose stack on Linux and macOS is 8 MiB, onto a worker,
    /// whose Rust default is 2 MiB. Matching the main thread keeps recursion
    /// depth what it was before the pool existed. The source body runs on the
    /// main thread and its OpenMP workers use the runtime's default.
    pub const WORKER_STACK_BYTES: usize = 8 << 20;

    /// Run `work`, normally the tool body, on a scoped worker pool sized by
    /// [`thread_policy`](Self::thread_policy), and return its result.
    ///
    /// This is the native form of the source applying `-threads` before
    /// `main_` (`TOPPBase.cpp:408-415`). A tool calls it from
    /// [`Tool::run`](crate::cli::Tool::run) with its body, so the setting
    /// reaches the run phase instead of being ignored; the five ported tools
    /// do, and a new tool is expected to. Inside `work`, rayon reports the policy's
    /// worker count (`rayon::current_num_threads`), and rayon parallel
    /// iterators run on this pool rather than on rayon's global one, so no
    /// computation of the run can exceed the requested count by accident. The
    /// pool is built for this call and its threads end when it returns; they
    /// are named `openms-<index>`, which is what `/proc/<pid>/task/*/comm`
    /// shows.
    ///
    /// A pool is built for one worker too. `-threads 1` therefore runs `work`
    /// on one pool thread rather than on the calling thread, which is what
    /// bounds rayon to one worker inside it.
    ///
    /// `work` must be [`Send`] because it runs on a pool thread. A tool that
    /// writes a report to the `out` and `err` streams of
    /// [`Tool::run_io`](crate::cli::Tool::run_io), which are not `Send`,
    /// computes inside `work` and writes after this returns.
    ///
    /// The determinism contract applies: a tool's output must not depend on
    /// the worker count. Outside `-test` the processing record still lists the
    /// `threads` parameter as given, as the source's does.
    ///
    /// Without the `parallel` feature, `work` runs on the calling thread.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when the operating system refuses to start the
    /// worker threads; nothing of `work` has run then. The source has no such
    /// path: libgomp aborts the process when it cannot create a thread. A
    /// panic inside `work` propagates to the caller unchanged.
    pub fn in_thread_pool<R, F>(&self, work: F) -> Result<R>
    where
        F: FnOnce() -> R + Send,
        R: Send,
    {
        #[cfg(feature = "parallel")]
        {
            let workers = self.thread_policy().get();
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .stack_size(Self::WORKER_STACK_BYTES)
                .thread_name(|index| format!("openms-{index}"))
                .build()
                .map_err(|error| {
                    Error::Io(std::io::Error::other(format!(
                        "cannot start {workers} worker threads: {error}"
                    )))
                })?;
            Ok(pool.install(work))
        }
        #[cfg(not(feature = "parallel"))]
        {
            Ok(work())
        }
    }
    /// A unique-id generator for this run's outputs.
    ///
    /// Under `-test` the generator starts from [`TEST_MODE_UNIQUE_ID_SEED`], as
    /// `UniqueIdGenerator::setSeed(19991231235959)` in `TOPPBase::main`, so the
    /// same draws are made on every machine. Otherwise it is seeded from the
    /// clock and process id. The source seeds one process-wide generator; each
    /// call here returns an independent generator starting at the first draw,
    /// so a tool should create one and reuse it.
    pub fn unique_id_generator(&self) -> UniqueIdGenerator {
        if self.test_mode {
            UniqueIdGenerator::from_seed(TEST_MODE_UNIQUE_ID_SEED)
        } else {
            UniqueIdGenerator::new()
        }
    }
    /// The processing record of this run, as `getProcessingInfo_`.
    ///
    /// See [`TEST_MODE_VERSION`](super::TEST_MODE_VERSION) for the values
    /// recorded under `-test`; otherwise the product version, the current
    /// local time and every resolved parameter as `parameter: <name>` are
    /// recorded.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a floating-point parameter cannot
    /// be represented as metadata.
    pub fn processing_info(&self, actions: &[ProcessingAction]) -> Result<DataProcessing> {
        processing::processing_info(
            &self.tool_name,
            &self.version,
            self.test_mode,
            &self.param,
            actions,
        )
    }
    /// Attach `processing` to an output map, as `addDataProcessing_`.
    ///
    /// A feature map and a consensus map gain one entry; every spectrum and
    /// chromatogram of an experiment gains one shared entry. Under `-test` a
    /// consensus map's column-header file names are reduced to base names.
    pub fn add_data_processing<T: AddDataProcessing + ?Sized>(
        &self,
        target: &mut T,
        processing: &DataProcessing,
    ) {
        target.add_data_processing(processing, self.test_mode);
    }

    fn value(&self, name: &str) -> Result<&ParamValue> {
        self.param
            .value(name)
            .map_err(|_| bad(format!("parameter '{name}' was not registered")))
    }

    /// A string, input-file, output-file or output-prefix option, as
    /// `getStringOption_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered or does
    /// not hold a string; the source throws `UnregisteredParameter` or
    /// `WrongParameterType`.
    pub fn string(&self, name: &str) -> Result<&str> {
        match self.value(name)? {
            ParamValue::String(text) => {
                self.write_debug(&format!("Value of string option '{name}': {text}"), 1);
                Ok(text)
            }
            _ => Err(bad(format!("parameter '{name}' is not a string"))),
        }
    }
    /// An output-directory option, as `getOutputDirOption`: the directory,
    /// created with its missing parents when it does not exist yet
    /// (`TOPPBase.cpp:1419-1439`). An empty value is returned as is, and
    /// nothing is created for it.
    ///
    /// # Errors
    ///
    /// As [`string`](Self::string), and [`Error::Io`] when the directory
    /// cannot be created. The source ignores `File::makeDir`'s result; the
    /// port reports the failure rather than hand the tool a directory that is
    /// not there.
    pub fn output_dir(&self, name: &str) -> Result<&str> {
        match self.value(name)? {
            ParamValue::String(text) => {
                self.write_debug(
                    &format!("Value of string(outdir) option '{name}': {text}"),
                    1,
                );
                if !text.is_empty() {
                    file::make_dir(text)?;
                }
                Ok(text)
            }
            _ => Err(bad(format!("parameter '{name}' is not a string"))),
        }
    }
    /// An integer option, as `getIntOption_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered or does
    /// not hold an integer.
    pub fn int(&self, name: &str) -> Result<i64> {
        match self.value(name)? {
            ParamValue::Integer(value) => {
                self.write_debug(&format!("Value of int option '{name}': {value}"), 1);
                Ok(*value)
            }
            _ => Err(bad(format!("parameter '{name}' is not an integer"))),
        }
    }
    /// A floating-point option, as `getDoubleOption_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered or does
    /// not hold a floating-point value. An integer is not widened: the source
    /// throws `WrongParameterType` for an integer option, and the strict update
    /// already refuses an INI value whose type differs from the registered one.
    pub fn double(&self, name: &str) -> Result<f64> {
        match self.value(name)? {
            ParamValue::Float(value) => {
                self.write_debug(
                    &format!("Value of double option '{name}': {}", double_text(*value)),
                    1,
                );
                Ok(*value)
            }
            _ => Err(bad(format!(
                "parameter '{name}' is not a floating-point number"
            ))),
        }
    }
    /// A string, input-file or output-file list, as `getStringList_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered or does
    /// not hold a string list.
    pub fn string_list(&self, name: &str) -> Result<&[String]> {
        match self.value(name)? {
            ParamValue::StringList(values) => {
                for value in values {
                    self.write_debug(&format!("Value of string option '{name}': {value}"), 1);
                }
                Ok(values)
            }
            _ => Err(bad(format!("parameter '{name}' is not a string list"))),
        }
    }
    /// An integer list, as `getIntList_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered or does
    /// not hold an integer list.
    pub fn int_list(&self, name: &str) -> Result<&[i32]> {
        match self.value(name)? {
            ParamValue::IntegerList(values) => {
                for value in values {
                    self.write_debug(&format!("Value of string option '{name}': {value}"), 1);
                }
                Ok(values)
            }
            _ => Err(bad(format!("parameter '{name}' is not an integer list"))),
        }
    }
    /// A floating-point list, as `getDoubleList_`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered or does
    /// not hold a floating-point list.
    pub fn double_list(&self, name: &str) -> Result<&[f64]> {
        match self.value(name)? {
            ParamValue::FloatList(values) => {
                for value in values {
                    self.write_debug(
                        &format!("Value of string option '{name}': {}", double_text(*value)),
                        1,
                    );
                }
                Ok(values)
            }
            _ => Err(bad(format!("parameter '{name}' is not a float list"))),
        }
    }
    /// Source `getFlag_`: a registered flag is true only when it was given.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not registered, does not
    /// hold a string, or holds a string other than `true` or `false`; the source
    /// throws `WrongParameterType` and `InvalidParameter` respectively.
    pub fn flag(&self, name: &str) -> Result<bool> {
        match self.value(name)? {
            ParamValue::String(text) if text == "true" => {
                self.write_debug(&format!("Value of string option '{name}': 1"), 1);
                Ok(true)
            }
            ParamValue::String(text) if text == "false" => {
                self.write_debug(&format!("Value of string option '{name}': 0"), 1);
                Ok(false)
            }
            ParamValue::String(text) => Err(bad(format!(
                "Invalid value '{text}' for flag parameter '{name}'. Valid values are 'true' and 'false' only."
            ))),
            _ => Err(bad(format!("parameter '{name}' is not a flag"))),
        }
    }
    /// Values of a registered subsection with the subsection prefix removed,
    /// as `getParam_().copy("<name>:", true)`.
    ///
    /// # Errors
    ///
    /// Propagates parameter-tree failures; an unknown subsection yields an
    /// empty tree, as the source copy does.
    pub fn subsection(&self, name: &str) -> Result<Param> {
        self.param.copy(&format!("{name}:"), true)
    }

    /// Source `writeLogInfo_`: `text` on `out` (the source's
    /// `OPENMS_LOG_INFO`) and a line in the `-log` file.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when `out` fails; the log file is best effort, as the
    /// source's unchecked stream is.
    pub fn write_log_info(&self, out: &mut dyn Write, text: &str) -> Result<()> {
        writeln!(out, "{text}")?;
        self.log.line(text);
        Ok(())
    }
    /// Source `writeLogWarn_`: `text` on `err` and a line in the `-log` file.
    ///
    /// # Errors
    ///
    /// As [`write_log_info`](Self::write_log_info).
    pub fn write_log_warn(&self, err: &mut dyn Write, text: &str) -> Result<()> {
        writeln!(err, "{text}")?;
        self.log.line(text);
        Ok(())
    }
    /// Source `writeLogError_`: `text` on `err` and a line in the `-log` file.
    ///
    /// # Errors
    ///
    /// As [`write_log_info`](Self::write_log_info).
    pub fn write_log_error(&self, err: &mut dyn Write, text: &str) -> Result<()> {
        self.write_log_warn(err, text)
    }
    /// Source `writeDebug_(text, min_level)`: a line in the `-log` file when
    /// the debug level is at least `min_level`. As in the Release build,
    /// nothing reaches the console.
    pub fn write_debug(&self, text: &str, min_level: u32) {
        self.log.debug(text, min_level);
    }
    /// Source `writeDebug_(text, param, min_level)`: `text` and `param`
    /// between separator lines in the `-log` file.
    pub fn write_debug_param(&self, text: &str, param: &Param, min_level: u32) {
        self.log.debug_param(text, param, min_level);
    }
}

/// `StringUtils::toStr(double)`, as the source's debug lines print a value.
fn double_text(value: f64) -> String {
    ParamValue::Float(value).to_text(true).unwrap_or_default()
}

/// Source `parseRange_` for floating-point bounds (`TOPPBase.cpp:2016-2052`).
///
/// `":8"`, `"2:"`, `"2:8"` and `":"` are accepted, and an absent side leaves
/// that bound untouched. The part before the first colon sets `low` and the
/// part after the last colon sets `high`, each converted like
/// `StringUtils::toDouble`. Returns whether any bound was set.
///
/// As in the source, `low > high` is not an error, and both bounds are returned
/// as given. The source's consumers may still swap them: DTAExtractor passes
/// its retention-time bounds to `DRange<1>(rt_l, rt_u)`, whose constructor
/// normalises reversed bounds (`DTAExtractor.cpp:144` at topp 174b576,
/// `DIntervalBase.h:85-90` at core bc9cc12), so the C++ tool reads `-rt 70:50`
/// as 50 to 70. A caller that needs that must order the bounds itself.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the colon is missing, where the source
/// throws `ConversionError` rather than reading `400` as `400:400`, or when a
/// bound is not a number. Neither bound is modified on error; the source may
/// already have assigned `low` when `high` fails to convert.
pub fn parse_range(text: &str, low: &mut f64, high: &mut f64) -> Result<bool> {
    if !text.contains(':') {
        return Err(bad(format!(
            "Invalid range '{text}': expected format '[min]:[max]' (the ':' separator is missing)"
        )));
    }
    let conversion = || {
        bad(format!(
            "Could not convert string '{text}' to a range of floating point values"
        ))
    };
    let start = text.split(':').next().unwrap_or("");
    let end = text.rsplit(':').next().unwrap_or("");
    let new_low = if start.is_empty() {
        None
    } else {
        Some(to_double(start).map_err(|_| conversion())?)
    };
    let new_high = if end.is_empty() {
        None
    } else {
        Some(to_double(end).map_err(|_| conversion())?)
    };
    if let Some(value) = new_low {
        *low = value;
    }
    if let Some(value) = new_high {
        *high = value;
    }
    Ok(new_low.is_some() || new_high.is_some())
}

/// Source `parseRange_` for integer bounds (`TOPPBase.cpp:2054-2090`).
///
/// The integer overload of [`parse_range`]: the same `[min]:[max]` grammar,
/// each bound converted like `StringUtils::toInt32`. Returns whether any bound
/// was set.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the colon is missing, with the
/// source's message, or when a bound is not a 32-bit integer
/// (`Could not convert string '<text>' to a range of integer values`).
/// Neither bound is modified on error.
pub fn parse_range_int(text: &str, low: &mut i32, high: &mut i32) -> Result<bool> {
    if !text.contains(':') {
        return Err(bad(format!(
            "Invalid range '{text}': expected format '[min]:[max]' (the ':' separator is missing)"
        )));
    }
    let conversion = || {
        bad(format!(
            "Could not convert string '{text}' to a range of integer values"
        ))
    };
    let start = text.split(':').next().unwrap_or("");
    let end = text.rsplit(':').next().unwrap_or("");
    let new_low = if start.is_empty() {
        None
    } else {
        Some(to_int32(start).map_err(|_| conversion())?)
    };
    let new_high = if end.is_empty() {
        None
    } else {
        Some(to_int32(end).map_err(|_| conversion())?)
    };
    if let Some(value) = new_low {
        *low = value;
    }
    if let Some(value) = new_high {
        *high = value;
    }
    Ok(new_low.is_some() || new_high.is_some())
}

/// Check that `filename` can be read, as `inputFileReadable_`
/// (`TOPPBase.cpp:1968-1995`).
///
/// Returns `None` when the file exists, is readable and — unless it is a
/// directory — holds at least one byte. Otherwise writes the source's two
/// diagnostic lines to `err` and returns the exit code the source exception
/// maps to: [`ExitCode::InputFileNotFound`] (`FileNotFound`),
/// [`ExitCode::InputFileNotReadable`] (`FileNotReadable`) or
/// [`ExitCode::InputFileEmpty`] (`FileEmpty`). `param_name` names the option
/// in the first line; an empty name gives the source's generic wording.
/// Diagnostics are best effort: a failing `err` does not change the result.
pub fn input_file_readable(
    filename: &str,
    param_name: &str,
    err: &mut dyn Write,
) -> Option<ExitCode> {
    let (code, heading, detail) = input_file_problem(filename, param_name)?;
    let _ = writeln!(err, "{heading}");
    let _ = writeln!(err, "{detail}");
    Some(code)
}

/// What [`input_file_readable`] reports, without writing it: the exit code,
/// the heading the source writes to its error log only and the `Error: …`
/// line its catch block writes to the error log and the `-log` file.
pub(crate) fn input_file_problem(
    filename: &str,
    param_name: &str,
) -> Option<(ExitCode, String, String)> {
    let (code, detail) = if !file::exists(filename) {
        (
            ExitCode::InputFileNotFound,
            format!("Error: File not found (the file '{filename}' could not be found)"),
        )
    } else if !file::readable(filename) {
        (
            ExitCode::InputFileNotReadable,
            format!(
                "Error: File not readable (the file '{filename}' is not readable for the current user)"
            ),
        )
    } else if !file::is_directory(filename) && file::empty(filename) {
        (
            ExitCode::InputFileEmpty,
            format!("Error: File empty (the file '{filename}' is empty)"),
        )
    } else {
        return None;
    };
    let heading = if param_name.is_empty() {
        "Cannot read input file!".to_owned()
    } else {
        format!("Cannot read input file given from parameter '-{param_name}'!")
    };
    Some((code, heading, detail))
}

/// Check that `filename` can be written, as `outputFileWritable_`
/// (`TOPPBase.cpp:1997-2013`).
///
/// Returns `None` when [`file::writable`] answers yes; that query never creates
/// or removes a file under the caller's name. Otherwise writes the source's two
/// diagnostic lines to `err` and returns [`ExitCode::CannotWriteOutputFile`],
/// the code of the source's `UnableToCreateFile`. Diagnostics are best effort.
pub fn output_file_writable(
    filename: &str,
    param_name: &str,
    err: &mut dyn Write,
) -> Option<ExitCode> {
    let (code, heading, detail) = output_file_problem(filename, param_name)?;
    let _ = writeln!(err, "{heading}");
    let _ = writeln!(err, "{detail}");
    Some(code)
}

/// What [`output_file_writable`] reports, without writing it, split as
/// [`input_file_problem`] splits it.
pub(crate) fn output_file_problem(
    filename: &str,
    param_name: &str,
) -> Option<(ExitCode, String, String)> {
    if file::writable(filename) {
        return None;
    }
    let heading = if param_name.is_empty() {
        "Cannot write output file!".to_owned()
    } else {
        format!("Cannot write output file given from parameter '-{param_name}'!")
    };
    Some((
        ExitCode::CannotWriteOutputFile,
        heading,
        format!("Error: Unable to write file (the file '{filename}' could not be created. )"),
    ))
}

/// The source's four whitespace characters, skipped from `index` on.
fn skip_whitespace(bytes: &[u8], mut index: usize) -> usize {
    while bytes
        .get(index)
        .is_some_and(|b| matches!(*b, b' ' | b'\t' | b'\n' | b'\r'))
    {
        index += 1;
    }
    index
}

/// Source `StringUtils::toInt32` (`StringUtils.cpp:136-166`), returning the
/// source's `ConversionError` message as the error.
///
/// Leading and trailing space, tab, newline and carriage return are skipped,
/// one `+` may precede the number, and the rest must be a complete decimal
/// `i32`. Like `std::from_chars` after the source strips `+`, a `-` may still
/// follow it.
pub(crate) fn to_int32(text: &str) -> std::result::Result<i32, String> {
    let bytes = text.as_bytes();
    let not_converted = || format!("Could not convert string '{text}' to an integer value");
    let mut cursor = skip_whitespace(bytes, 0);
    if cursor == bytes.len() {
        return Err(not_converted());
    }
    if bytes[cursor] == b'+' {
        cursor += 1;
    }
    let number_start = cursor;
    if bytes.get(cursor) == Some(&b'-') {
        cursor += 1;
    }
    let digits_start = cursor;
    while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
        cursor += 1;
    }
    if cursor == digits_start {
        return Err(not_converted());
    }
    let value = text
        .get(number_start..cursor)
        .and_then(|digits| digits.parse::<i32>().ok())
        .ok_or_else(not_converted)?;
    let after = skip_whitespace(bytes, cursor);
    if after != bytes.len() {
        return Err(format!(
            "Prefix of string '{text}' successfully converted to an int32 value. Additional characters found at position {}",
            after + 1
        ));
    }
    Ok(value)
}

/// Length of a case-insensitive ASCII `word` at `bytes[index..]`, if present.
fn word_at(bytes: &[u8], index: usize, word: &str) -> Option<usize> {
    let candidate = bytes.get(index..index.checked_add(word.len())?)?;
    candidate
        .eq_ignore_ascii_case(word.as_bytes())
        .then_some(word.len())
}

/// Source `StringUtils::toDouble` (`StringUtils.cpp:239-276`), returning the
/// source's `ConversionError` message as the error.
///
/// Whitespace is skipped on both sides and `nan`, optionally followed by a
/// parenthesised payload, is accepted before anything else. Then one `+` may
/// precede a `std::from_chars` general-format number: an optional `-`, digits
/// with an optional fraction and exponent, or `inf`/`infinity`. A finite
/// literal that overflows is an error, as `result_out_of_range` is; underflow
/// rounds, as the oracle's libc++ fallback accepts it. Hexadecimal floats are
/// rejected, following `std::from_chars` rather than the `strtod` fallback the
/// source uses on libc++.
pub(crate) fn to_double(text: &str) -> std::result::Result<f64, String> {
    let bytes = text.as_bytes();
    let not_converted = || format!("Could not convert string '{text}' to a double value");
    let first = skip_whitespace(bytes, 0);
    if first == bytes.len() {
        return Err(not_converted());
    }
    if let Some(length) = word_at(bytes, first, "nan") {
        let mut end = first + length;
        if bytes.get(end) == Some(&b'(') {
            match bytes[end..].iter().position(|b| *b == b')') {
                Some(close) => end += close + 1,
                None => end = first,
            }
        }
        if end != first && skip_whitespace(bytes, end) == bytes.len() {
            return Ok(f64::NAN);
        }
    }
    let mut start = first;
    if bytes[start] == b'+' {
        start += 1;
    }
    let mut cursor = start;
    if bytes.get(cursor) == Some(&b'-') {
        cursor += 1;
    }
    let infinity = word_at(bytes, cursor, "infinity").or_else(|| word_at(bytes, cursor, "inf"));
    let nan = word_at(bytes, cursor, "nan");
    if let Some(length) = infinity.or(nan) {
        cursor += length;
    } else {
        let integer_start = cursor;
        while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
            cursor += 1;
        }
        let mut digits = cursor - integer_start;
        if bytes.get(cursor) == Some(&b'.') {
            let fraction_start = cursor + 1;
            let mut fraction_end = fraction_start;
            while bytes.get(fraction_end).is_some_and(u8::is_ascii_digit) {
                fraction_end += 1;
            }
            if digits > 0 || fraction_end > fraction_start {
                digits += fraction_end - fraction_start;
                cursor = fraction_end;
            }
        }
        if digits == 0 {
            return Err(not_converted());
        }
        if matches!(bytes.get(cursor), Some(b'e' | b'E')) {
            let mut exponent = cursor + 1;
            if matches!(bytes.get(exponent), Some(b'+' | b'-')) {
                exponent += 1;
            }
            let exponent_digits = exponent;
            while bytes.get(exponent).is_some_and(u8::is_ascii_digit) {
                exponent += 1;
            }
            if exponent > exponent_digits {
                cursor = exponent;
            }
        }
    }
    let literal = text.get(start..cursor).ok_or_else(not_converted)?;
    let value: f64 = if nan.is_some() && infinity.is_none() {
        f64::NAN
    } else {
        literal.parse().map_err(|_| not_converted())?
    };
    if value.is_infinite() && infinity.is_none() {
        return Err(not_converted());
    }
    let after = skip_whitespace(bytes, cursor);
    if after != bytes.len() {
        return Err(format!(
            "Prefix of string '{text}' successfully converted to a double value. Additional characters found at position {}",
            after + 1
        ));
    }
    Ok(value)
}

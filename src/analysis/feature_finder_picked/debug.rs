// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The debug mode of the picked feature finder (`write_debug`,
//! `FEATUREFINDER/FeatureFinderAlgorithmPicked.h`), returned as data.
//!
//! With `write_debug = true` the source writes intermediate results below the
//! working directory (`FeatureFinderAlgorithmPicked.cpp:226-232`, `550-571`,
//! `714-718`, `1028-1053`, `2130-2203`):
//!
//! | File | Source | Here |
//! |---|---|---|
//! | `debug/` and `debug/features/` | `File::makeDir` at `:230` | the caller creates them ([`DebugOutput`]) |
//! | `debug/log.txt` | the member `log_`, opened at `:231`, 58 write statements | [`DebugOutput::log`] |
//! | `debug/seeds_<charge>.featureXML` | `:550-571`, for every charge | [`DebugOutput::seed_maps`] ([`seed_map`]) |
//! | `debug/features/<n>.dta`, `<n>_cropped.dta`, `<n>.plot` | `writeFeatureDebugInfo_` at `:714-718` | [`DebugOutput::feature_files`] ([`write_feature_debug_info`]) |
//! | `debug/abort_reasons.featureXML` | the member `abort_reasons_`, `:1028-1045` | [`DebugOutput::abort_reasons`] ([`abort_map`]) |
//! | `debug/input.mzML` | `:1047-1053` | [`DebugOutput::input`] ([`debug_experiment`]) |
//!
//! Library code does not write files. The algorithm instance
//! ([`crate::analysis::feature_finder_picked::instance::FeatureFinderAlgorithmPicked`])
//! fills a [`DebugOutput`] instead, and the `FeatureFinderCentroided` tool
//! writes it under the names above, in the source's order.
//!
//! # Where the source is undefined
//!
//! - **`writeFeatureDebugInfo_` terminates the process.** It reads the
//!   parameter `debug:pseudo_rt_shift` (`:2137`), which is not declared; the
//!   declared one is `advanced:pseudo_rt_shift` (`:124`). `Param::getValue`
//!   throws `ElementNotFound`, and because the call sits in an `omp critical`
//!   section inside the `omp parallel for` of `:595`, the exception escapes an
//!   OpenMP region and the runtime calls `std::terminate`: the executed
//!   FeatureFinderCentroided prints OpenMS's "FATAL: uncaught exception!"
//!   block and is killed by `SIGABRT` at the first seed that reaches the fit.
//!   A safe port cannot end the process. The source-following
//!   [`PseudoRtShiftKey::Source`](crate::analysis::feature_finder_picked::algorithm::PseudoRtShiftKey::Source)
//!   returns [`Error::Unsupported`] at exactly that seed and records the point
//!   in [`DebugOutput::termination`];
//!   [`PseudoRtShiftKey::Declared`](crate::analysis::feature_finder_picked::algorithm::PseudoRtShiftKey::Declared)
//!   reads the declared parameter and produces the member's files.
//! - **A string or list under `debug:pseudo_rt_shift`.** The source's
//!   `double` conversion then returns the bits of a heap pointer, a tiny
//!   positive number that changes from process to process
//!   ([`PseudoRtShift::HeapAddress`]). Only a shifted value that small an
//!   addend can change shows it (for the first few traces, a retention time
//!   below about `1e-293` or a fitted centre below about `1e-303`); everywhere
//!   else the port writes the source's bytes, and [`write_feature_debug_info`]
//!   refuses at the first value that the address can change.
//! - **`abort_` races.** `abort_` (`:1129-1140`) writes `log_` and
//!   `abort_reasons_` from inside the parallel region without synchronisation;
//!   the parameter text itself says "do not use in parallel mode". With more
//!   than one thread the source is undefined (executed: the log differs from
//!   run to run). This port collects every seed's lines and abort in seed
//!   order, so its output at any thread count is the source's single-thread
//!   output: a superset of the source's defined domain.
//! - **Stale seeds of an earlier run.** `abort_reasons_` is never cleared, and
//!   its seeds hold spectrum and peak indices of the run that stored them.
//!   `:1037-1039` reads them from the current map without a bounds check; an
//!   index outside the current map is an out-of-bounds read, and [`abort_map`]
//!   refuses exactly there.
//!
//! See `docs/FEATURE_FINDER_PICKED_SUPPORT.md`, section "Debug mode".
//!
//! [`DebugOutput`]: crate::analysis::feature_finder_picked::debug::DebugOutput
//! [`DebugOutput::log`]: crate::analysis::feature_finder_picked::debug::DebugOutput::log
//! [`DebugOutput::seed_maps`]: crate::analysis::feature_finder_picked::debug::DebugOutput::seed_maps
//! [`seed_map`]: crate::analysis::feature_finder_picked::debug::seed_map
//! [`DebugOutput::feature_files`]: crate::analysis::feature_finder_picked::debug::DebugOutput::feature_files
//! [`write_feature_debug_info`]: crate::analysis::feature_finder_picked::debug::write_feature_debug_info
//! [`DebugOutput::abort_reasons`]: crate::analysis::feature_finder_picked::debug::DebugOutput::abort_reasons
//! [`abort_map`]: crate::analysis::feature_finder_picked::debug::abort_map
//! [`DebugOutput::input`]: crate::analysis::feature_finder_picked::debug::DebugOutput::input
//! [`debug_experiment`]: crate::analysis::feature_finder_picked::debug::debug_experiment
//! [`Error::Unsupported`]: crate::Error::Unsupported
//! [`DebugOutput::termination`]: crate::analysis::feature_finder_picked::debug::DebugOutput::termination
//! [`PseudoRtShift::HeapAddress`]: crate::analysis::feature_finder_picked::debug::PseudoRtShift::HeapAddress

use std::collections::BTreeMap;

use crate::analysis::feature_finder_picked::helper_structs::{MassTraces, Seed};
use crate::analysis::feature_finder_picked::scoring::ScoreArrays;
use crate::analysis::feature_finder_picked::trace_fitter::TraceFitter;
use crate::format::file_info::text_format::{fixed_truncated, ostream_g, to_str, to_str_f32};
use crate::kernel::{DataArray, Feature, FeatureMap, MSExperiment};
use crate::metadata::MetaValue;
use crate::{Error, Result};

/// The directory `writeFeatureDebugInfo_` writes into by default: its `path`
/// argument, `"debug/features/"`.
pub const FEATURE_DEBUG_PATH: &str = "debug/features/";

/// The source's `log_` line of step 1 (`FeatureFinderAlgorithmPicked.cpp:238`).
pub const LOG_PRECALCULATING: &str = "Precalculating intensity thresholds ...\n";

/// Capacity of the put area of libstdc++'s `basic_filebuf<char>`:
/// `_M_buf_size - 1` with `_M_buf_size = BUFSIZ = 8192`.
const FILEBUF_CAPACITY: usize = 8191;

/// A destination of the source's `log_` stream.
///
/// Each call of [`Self::put`] is one `<<` insertion of the source, with its
/// formatted text; the boundaries matter to [`DebugLog::flushed_bytes`].
pub(crate) trait LogSink {
    /// Whether the stream is written at all (the source's `debug_`).
    fn enabled(&self) -> bool;
    /// One insertion.
    fn put(&mut self, text: &str);
}

/// The stream of a run without `write_debug`: nothing is written.
pub(crate) struct NoLog;

impl LogSink for NoLog {
    #[inline(always)]
    fn enabled(&self) -> bool {
        false
    }
    #[inline(always)]
    fn put(&mut self, _: &str) {}
}

/// Text written to the stream by one part of a run (one seed, one charge's
/// pattern scores), with the length of every insertion.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct LogFragment {
    text: String,
    insertions: Vec<u32>,
}

impl LogFragment {
    pub(crate) fn new() -> Self {
        Self::default()
    }
}

/// The insertion length [`LogFragment`] records for a single `char`.
const CHAR_INSERTION: u32 = 0;

impl LogSink for LogFragment {
    fn enabled(&self) -> bool {
        true
    }
    fn put(&mut self, text: &str) {
        self.text.push_str(text);
        if text == "\n" {
            // Every lone newline the source writes to `log_` is the `char`
            // '\n', which libstdc++ inserts with `sputc`, not `sputn`.
            self.insertions.push(CHAR_INSERTION);
        } else if !text.is_empty() {
            // An empty `sputn` writes only when the put area is full
            // (`0 >= 0` in `xsputn`), and only a `sputc` leaves it full; every
            // empty insertion of the source follows a non-empty string
            // (`" peaks (abort: "`, `"Abort: "`), so it changes nothing.
            self.insertions
                .push(u32::try_from(text.len()).unwrap_or(u32::MAX));
        }
    }
}

/// glibc's spelling of a NaN in `%g` and `%f`: `-nan` when the sign bit is
/// set, `nan` otherwise.
///
/// The shared formatters print `nan` whatever the sign, as Apple libc does; the
/// reference build is Linux x86_64, where glibc keeps the sign, and a NaN made
/// by `0.0 / 0.0` there carries it. IEEE 754 leaves the sign of a computed NaN
/// unspecified, so this matches the C++ wherever both sides produce the NaN by
/// the same operations.
fn glibc_nan(value: f64) -> Option<String> {
    value.is_nan().then(|| {
        if value.is_sign_negative() {
            "-nan".to_owned()
        } else {
            "nan".to_owned()
        }
    })
}

/// `std::ostream << double` with the default precision 6.
pub(crate) fn g(value: f64) -> String {
    glibc_nan(value).unwrap_or_else(|| ostream_g(value, 6))
}

/// `std::ostream << float`: the promoted value at precision 6.
pub(crate) fn g32(value: f32) -> String {
    g(f64::from(value))
}

/// `StringUtils::number(value, digits)`: `snprintf("%.*f")` into 64 bytes.
pub(crate) fn number(value: f64, digits: u32) -> String {
    glibc_nan(value).unwrap_or_else(|| fixed_truncated(value, digits))
}

/// Write `parts` as consecutive insertions, when the stream is enabled. A
/// part that is exactly `"\n"` stands for the source's `'\n'` character.
pub(crate) fn put_all<L: LogSink + ?Sized>(log: &mut L, parts: &[&str]) {
    if log.enabled() {
        for part in parts {
            log.put(part);
        }
    }
}

/// The text of the source's `debug/log.txt` for one run, and where the
/// libstdc++ file buffer stood.
///
/// The source writes through an `std::ofstream` with the default 8192-byte
/// buffer, whose put area holds 8191 bytes. When the process ends normally the
/// destructor flushes everything; when it is killed by `std::terminate` the
/// buffered tail is lost, and the file holds only what the buffer had already
/// handed to the operating system. [`Self::flushed_bytes`] tracks that length
/// exactly, following libstdc++: a string or number goes through
/// `basic_filebuf::xsputn`, which copies it into the buffer when it is
/// shorter than the free space and otherwise writes the buffered bytes and the
/// insertion directly (`fstream.tcc`); a single `char`, which the source
/// writes for every lone `'\n'`, goes through `sputc`, which fills the buffer
/// to the last byte and then writes all 8192. The executed files confirm the
/// model: it reproduces the lengths of all three log files the Release build
/// left while buffered (two terminated runs and one open stream). See
/// `docs/FEATURE_FINDER_PICKED_SUPPORT.md`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DebugLog {
    text: String,
    flushed: usize,
    buffered: usize,
}

impl DebugLog {
    /// Everything written in this run, in source order.
    pub fn text(&self) -> &str {
        &self.text
    }

    /// The lines of [`Self::text`], each without its `'\n'`. A final
    /// unterminated line, which the source never writes, is included.
    pub fn lines(&self) -> impl Iterator<Item = &str> {
        self.text.lines()
    }

    /// How many leading bytes of [`Self::text`] the source's file buffer had
    /// written to the file at the end of this run's last insertion.
    pub fn flushed_bytes(&self) -> usize {
        self.flushed
    }

    /// Bytes in total.
    pub fn len(&self) -> usize {
        self.text.len()
    }

    /// Whether nothing was written.
    pub fn is_empty(&self) -> bool {
        self.text.is_empty()
    }

    /// Append a fragment, advancing the file-buffer model insertion by
    /// insertion.
    pub(crate) fn append(&mut self, fragment: &LogFragment, limit: usize) -> Result<()> {
        if self.text.len().saturating_add(fragment.text.len()) > limit {
            return Err(Error::InvalidValue(format!(
                "the debug log exceeds the limit of {limit} bytes"
            )));
        }
        self.text
            .try_reserve(fragment.text.len())
            .map_err(|_| Error::InvalidValue("cannot allocate the debug log".into()))?;
        self.text.push_str(&fragment.text);
        for &length in &fragment.insertions {
            if length == CHAR_INSERTION {
                // `sputc`: into the buffer while there is room, otherwise
                // `overflow` writes the full buffer and the character.
                if self.buffered < FILEBUF_CAPACITY {
                    self.buffered += 1;
                } else {
                    self.flushed += self.buffered + 1;
                    self.buffered = 0;
                }
                continue;
            }
            // `xsputn`: an insertion at least as long as the free space is
            // written directly together with the buffer.
            let length = length as usize;
            let available = FILEBUF_CAPACITY - self.buffered;
            if length >= available {
                self.flushed += self.buffered + length;
                self.buffered = 0;
            } else {
                self.buffered += length;
            }
        }
        Ok(())
    }
}

/// One line of the run's console output, tagged with the source stream that
/// prints it, or a point where the source stores a debug file.
///
/// The executed tool routes `std::cout` straight to standard output, while
/// `OPENMS_LOG_INFO` and `OPENMS_LOG_WARN` each hold back a line that repeats
/// one of the two they printed last and report the repeats later as
/// `<line> occurred N times` (`LogStreamBuf`, `LogStream.cpp:180-300`). A
/// caller that reproduces the console keeps the streams apart for that
/// reason. The store points matter because `FeatureXMLFile::store` itself
/// prints on `OPENMS_LOG_INFO`.
#[derive(Clone, Debug, PartialEq)]
pub enum ReportLine {
    /// A line the source writes to `std::cout`.
    Out(String),
    /// A line the source writes to `OPENMS_LOG_INFO`; may be empty.
    Info(String),
    /// A line the source writes to `OPENMS_LOG_WARN`.
    Warn(String),
    /// The source stores [`DebugOutput::seed_maps`]`[index]` here.
    StoreSeedMap(usize),
    /// The source stores [`DebugOutput::abort_reasons`] here.
    StoreAbortReasons,
    /// The source stores [`DebugOutput::input`] here.
    StoreInput,
}

impl ReportLine {
    /// The text of a console line, or `None` for a store point.
    pub fn text(&self) -> Option<&str> {
        match self {
            Self::Out(text) | Self::Info(text) | Self::Warn(text) => Some(text),
            Self::StoreSeedMap(_) | Self::StoreAbortReasons | Self::StoreInput => None,
        }
    }
}

/// The seeds of one charge as the source stores them in debug mode
/// (`debug/seeds_<charge>.featureXML`).
#[derive(Clone, Debug, PartialEq)]
pub struct SeedMap {
    /// The charge; the file name is `seeds_<charge>.featureXML`.
    pub charge: i32,
    /// The seeds as features; the map and the features keep unique id 0.
    pub map: FeatureMap,
}

/// The three files `writeFeatureDebugInfo_` writes for one seed that reached
/// the fit.
#[derive(Clone, Debug, PartialEq)]
pub struct FeatureDebugFiles {
    /// The seed's `plot_nr`, which names the files.
    pub plot_nr: i64,
    /// The directory prefix the file names and the gnuplot script use.
    pub path: String,
    /// `<path><plot_nr>.dta`: every peak of the extended traces, as `pseudo RT`
    /// and intensity separated by a tab, one per line.
    pub dta: String,
    /// `<path><plot_nr>_cropped.dta`, the same for the cropped traces; `None`
    /// when they hold no peak, where the source writes no file.
    pub cropped_dta: Option<String>,
    /// `<path><plot_nr>.plot`, the gnuplot script. Bytes, because the source
    /// names the trace functions `'f' + k` in a C++ `char`, which leaves ASCII
    /// from the 27th trace on.
    pub plot: Vec<u8>,
}

impl FeatureDebugFiles {
    /// The file name of [`Self::dta`].
    pub fn dta_name(&self) -> String {
        format!("{}{}.dta", self.path, self.plot_nr)
    }

    /// The file name of [`Self::cropped_dta`].
    pub fn cropped_dta_name(&self) -> String {
        format!("{}{}_cropped.dta", self.path, self.plot_nr)
    }

    /// The file name of [`Self::plot`].
    pub fn plot_name(&self) -> String {
        format!("{}{}.plot", self.path, self.plot_nr)
    }
}

/// Where the source process ends in `writeFeatureDebugInfo_`, and why.
#[derive(Clone, Debug, PartialEq)]
pub struct DebugTermination {
    /// The charge being extended.
    pub charge: i32,
    /// The seed's position in that charge's seed list.
    pub seed_index: usize,
    /// The `plot_nr` the seed received.
    pub plot_nr: i64,
    /// The C++ exception class that escapes the OpenMP region.
    pub exception: &'static str,
    /// Its `what()` text, as the executed build prints it.
    pub message: String,
}

/// Everything `write_debug` produces in one run, in the order the source
/// writes it.
///
/// A caller that writes it follows the source: create `debug/features`;
/// create (truncate) `debug/log.txt` with [`Self::log`] when
/// [`Self::log_opened`] is set, and leave the file alone otherwise; store each
/// seed map, the feature files, the abort map and the input as named in the
/// module documentation. After a [`Self::termination`] the executed build has
/// written only [`DebugLog::flushed_bytes`] of the log and nothing after the
/// last seed map.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DebugOutput {
    /// Whether this run opened `debug/log.txt`, truncating it.
    ///
    /// The source's `log_` is an instance member that is opened in every debug
    /// run and never closed. A second open of an open `std::ofstream` fails,
    /// which sets `failbit`: every later write of the instance is dropped, and
    /// the file keeps what the first debug run wrote. `false` here, with an
    /// empty [`Self::log`], is that case.
    pub log_opened: bool,
    /// The run's `log_` output.
    pub log: DebugLog,
    /// One seed map per charge the run reached, from `charge_low` upwards.
    pub seed_maps: Vec<SeedMap>,
    /// The `writeFeatureDebugInfo_` output, by `plot_nr`.
    pub feature_files: Vec<FeatureDebugFiles>,
    /// The abort reasons as features, when the run reached `:1028`.
    pub abort_reasons: Option<FeatureMap>,
    /// The input with its score arrays, when the run reached `:1047`.
    pub input: Option<MSExperiment>,
    /// Where the source process terminates, when it does.
    pub termination: Option<DebugTermination>,
}

/// A key of `std::map<Seed, std::string>`: `Seed::operator<` compares only
/// the `float` intensities, so seeds of equal intensity, `-0.0` and `+0.0`
/// included, are one key.
#[derive(Clone, Copy, Debug, PartialEq)]
struct IntensityKey(f32);

impl Eq for IntensityKey {}

impl PartialOrd for IntensityKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for IntensityKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Keys are not NaN (checked on insertion) and their zeros are
        // normalised, so the total order is the `<` order.
        self.0.total_cmp(&other.0)
    }
}

/// The source member `abort_reasons_`: `std::map<Seed, std::string>`.
///
/// Keyed by the seed's intensity alone ([`Seed::is_less_intense_than`]). The
/// first seed stored under an intensity keeps its spectrum and peak, and every
/// later abort of an equally intense seed replaces the reason only. Entries
/// iterate by ascending intensity. The source fills it only in debug runs and
/// never clears it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AbortReasons {
    entries: BTreeMap<IntensityKey, (Seed, String)>,
}

impl AbortReasons {
    /// An empty map.
    pub fn new() -> Self {
        Self::default()
    }

    /// `abort_reasons_[seed] = reason`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a NaN intensity, which the ordering
    /// of the source map cannot hold (it is not a strict weak ordering there);
    /// validated input has no such peak.
    pub fn insert(&mut self, seed: Seed, reason: &str) -> Result<()> {
        if seed.intensity.is_nan() {
            return Err(Error::InvalidValue(
                "a seed with a NaN intensity cannot key the abort-reason map".into(),
            ));
        }
        // `-0.0 < 0.0` is false both ways: one key.
        let key = IntensityKey(if seed.intensity == 0.0 {
            0.0
        } else {
            seed.intensity
        });
        match self.entries.get_mut(&key) {
            Some(entry) => {
                entry.1.clear();
                entry.1.push_str(reason);
            }
            None => {
                self.entries.insert(key, (seed, reason.to_owned()));
            }
        }
        Ok(())
    }

    /// The entries by ascending intensity: the stored seed and its last
    /// reason.
    pub fn iter(&self) -> impl Iterator<Item = (&Seed, &str)> {
        self.entries
            .values()
            .map(|(seed, reason)| (seed, reason.as_str()))
    }

    /// Number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the map is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// The debug seed map of one charge (`FeatureFinderAlgorithmPicked.cpp:550-571`).
///
/// One feature per seed, in seed order: intensity the seed's, overall quality
/// the peak's overall score of this charge, RT the scan's, m/z the peak's, and
/// the meta values `intensity_score`, `pattern_score` and `trace_score`, each
/// the stored `float` score promoted to `double` as `DataValue(float)` stores
/// it. Unique ids stay 0, as in the source, which never assigns them here.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a seed does not address a peak or a
/// score is missing, which a seed of the run cannot cause, and when a score is
/// not finite, which the native meta values cannot hold; a seed's scores are
/// finite because its overall score reached the seed threshold.
pub fn seed_map(
    experiment: &MSExperiment,
    scores: &ScoreArrays,
    charge_index: usize,
    seeds: &[Seed],
) -> Result<FeatureMap> {
    let missing = || Error::InvalidValue("a debug seed does not address a scored peak".into());
    let mut features = Vec::new();
    features
        .try_reserve_exact(seeds.len())
        .map_err(|_| Error::InvalidValue("cannot allocate the debug seed map".into()))?;
    for seed in seeds {
        let spectrum = experiment.spectra.get(seed.spectrum).ok_or_else(missing)?;
        let peak = spectrum.peaks.get(seed.peak).ok_or_else(missing)?;
        let score = |values: Option<&[f32]>| -> Result<f32> {
            values
                .and_then(|values| values.get(seed.peak).copied())
                .ok_or_else(missing)
        };
        let overall = score(scores.overall(charge_index, seed.spectrum))?;
        let mut feature = Feature::new(spectrum.rt, peak.mz, seed.intensity);
        feature.quality = overall;
        for (key, value) in [
            ("intensity_score", score(scores.intensity(seed.spectrum))?),
            (
                "pattern_score",
                score(scores.pattern(charge_index, seed.spectrum))?,
            ),
            ("trace_score", score(scores.trace(seed.spectrum))?),
        ] {
            feature
                .metadata
                .insert(key.into(), MetaValue::try_from(f64::from(value))?);
        }
        features.push(feature);
    }
    Ok(FeatureMap::from_features(features))
}

/// The debug abort map (`FeatureFinderAlgorithmPicked.cpp:1028-1045`).
///
/// One feature per entry of `reasons`, by ascending intensity: RT, m/z and
/// intensity read from `experiment` at the stored seed's spectrum and peak,
/// the meta value `label` holding the reason, and unique id `k` for the `k`-th
/// entry, so the first one has the invalid id 0. The source then gives the map
/// a fresh unique id from the process-wide generator (`setUniqueId()`); the map
/// here keeps id 0 and the caller draws it, as the tool does.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a stored seed does not address a peak
/// of `experiment`. The source reads `map_[spectrum][peak]` without a check;
/// that happens for seeds an earlier run of the same instance stored, whose
/// indices refer to that run's input, and it is an out-of-bounds read with no
/// reproducible result. Entries inside the current input are read from it, as
/// the source does, whatever run stored them.
pub fn abort_map(reasons: &AbortReasons, experiment: &MSExperiment) -> Result<FeatureMap> {
    let mut features = Vec::new();
    features
        .try_reserve_exact(reasons.len())
        .map_err(|_| Error::InvalidValue("cannot allocate the debug abort map".into()))?;
    for (counter, (seed, reason)) in reasons.iter().enumerate() {
        let peak = experiment.spectra.get(seed.spectrum).and_then(|spectrum| {
            spectrum
                .peaks
                .get(seed.peak)
                .map(|peak| (spectrum.rt, peak))
        });
        let Some((rt, peak)) = peak else {
            return Err(Error::InvalidValue(format!(
                "abort_reasons_ holds a seed at spectrum {} peak {} from an earlier run, outside \
                 the current input of {} spectra; the source reads it without a bounds check, \
                 which is undefined behaviour",
                seed.spectrum,
                seed.peak,
                experiment.spectra.len()
            )));
        };
        let mut feature = Feature::new(rt, peak.mz, peak.intensity);
        feature
            .metadata
            .insert("label".into(), MetaValue::from(reason.to_owned()));
        feature.unique_id = counter as u64;
        features.push(feature);
    }
    Ok(FeatureMap::from_features(features))
}

/// The input as the source stores it in `debug/input.mzML`
/// (`FeatureFinderAlgorithmPicked.cpp:198-224` and `1047-1052`).
///
/// Each spectrum's float data arrays are resized to `3 + 2 * charges`
/// (`std::vector::resize`: arrays the input already had keep their position
/// and description, a longer list is cut, a shorter one is padded), then named
/// and filled: `trace_score`, `intensity_score`, `local_max`, then
/// `pattern_score_<c>` and `overall_score_<c>` for every charge. Array 2,
/// `local_max`, is then erased. The source's comment there says "without
/// overall score", but the code removes the local-maximum flags; the port
/// follows the code. Integer and string arrays are untouched.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `scores` does not match `experiment`,
/// or an array cannot be allocated.
pub fn debug_experiment(
    mut experiment: MSExperiment,
    scores: &ScoreArrays,
) -> Result<MSExperiment> {
    if scores.spectrum_count() != experiment.spectra.len() {
        return Err(Error::InvalidValue(
            "the score arrays do not belong to the experiment".into(),
        ));
    }
    let names = scores.array_names();
    let charges = scores.charge_count();
    for (index, spectrum) in experiment.spectra.iter_mut().enumerate() {
        let missing = || Error::InvalidValue("a spectrum has no score arrays".into());
        let mut values: Vec<&[f32]> = vec![
            scores.trace(index).ok_or_else(missing)?,
            scores.intensity(index).ok_or_else(missing)?,
            scores.local_max(index).ok_or_else(missing)?,
        ];
        for charge in 0..charges {
            values.push(scores.pattern(charge, index).ok_or_else(missing)?);
        }
        for charge in 0..charges {
            values.push(scores.overall(charge, index).ok_or_else(missing)?);
        }
        let arrays = &mut spectrum.float_data_arrays;
        arrays.truncate(names.len());
        arrays
            .try_reserve(names.len() - arrays.len())
            .map_err(|_| Error::InvalidValue("cannot allocate the debug arrays".into()))?;
        while arrays.len() < names.len() {
            arrays.push(DataArray::default());
        }
        for ((array, name), data) in arrays.iter_mut().zip(&names).zip(values) {
            array.name.clone_from(name);
            array.data.clear();
            array
                .data
                .try_reserve_exact(data.len())
                .map_err(|_| Error::InvalidValue("cannot allocate the debug arrays".into()))?;
            array.data.extend_from_slice(data);
        }
        arrays.remove(2);
    }
    Ok(experiment)
}

/// The pseudo-RT shift `writeFeatureDebugInfo_` works with.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PseudoRtShift {
    /// A numeric parameter value, converted to `double` as the source converts
    /// it (an integer by `double(value)`).
    Value(f64),
    /// A string or list parameter value. `ParamValue::operator double()`
    /// (`ParamValue.cpp:397-408`) returns the union member `dou_` for every
    /// type but `EMPTY` and `INT`; the Release build's code is one `movsd
    /// 0x8(%rdi),%xmm0` (`libOpenMS.so`, `_ZNK6OpenMS10ParamValuecvdEv`), so
    /// the shift is the bit pattern of the pointer the union holds, the
    /// `std::string` or `std::vector` the constructor allocated with `new`
    /// (`ParamValue.cpp:104-116`).
    ///
    /// That pointer lies below [`HEAP_ADDRESS_END`], so read as a `double` it
    /// is a subnormal number below `2^-1027`, and it changes from process to
    /// process (executed: nine values in three processes, all below `2^47` as
    /// integers). The port cannot know it. Trace `k` of the `.dta` files is
    /// written as `pseudo_rt_shift * k + rt`, and the formula of trace `k` in
    /// the `.plot` file prints `pseudo_rt_shift * k + centre` with six
    /// significant digits; rounded addition, rounded multiplication and that
    /// rounding are monotone, so a value that the largest possible shift,
    /// `HEAP_ADDRESS_END * k` as a `double`, leaves unchanged is unchanged by
    /// every possible address. There the source writes the text of shift `0`,
    /// and so does the port (executed: FeatureFinderCentroided_1 with a string
    /// and a string-list value, and the same input moved so that one scan sits
    /// at RT `5e-275` or `1e-289`; byte-identical files in three processes
    /// each, with a different address in each). Where that shift changes a
    /// written value, a peak of trace `k >= 1` (extended or cropped) or the
    /// fitted centre in such a trace's formula, the text is the address's and
    /// [`write_feature_debug_info`] refuses (executed: the scan moved to RT
    /// `0`, `1e-295` or `1e-300` gave a different `0.dta` in each of three
    /// processes). The refusal is exact for that address range: every value it
    /// refuses is changed by an address just below [`HEAP_ADDRESS_END`] (the
    /// product `HEAP_ADDRESS_END * k` is never a power of two, so no rounding
    /// tie lies between the two), while the lowest mappable address
    /// (`mmap_min_addr`, 65536 on the reference host) changes it less or not
    /// at all, so the written text depends on the address the process got. A
    /// NaN or an infinity absorbs every shift.
    HeapAddress,
}

/// The end of the address range a heap pointer can take in a Linux x86_64
/// process: `DEFAULT_MAP_WINDOW = (1 << 47) - PAGE_SIZE`
/// (`arch/x86/include/asm/page_64_types.h`, kernel 6.8 on the reference
/// host).
///
/// With four-level paging (the reference host: 48-bit virtual addresses, no
/// `la57`) it is the end of user space, `TASK_SIZE_MAX`. With five-level
/// paging the kernel maps above it only for an `mmap` hint above it; the
/// program break and the default `mmap` area both lie below it
/// (`ELF_ET_DYN_BASE = DEFAULT_MAP_WINDOW / 3 * 2`, `TASK_UNMAPPED_BASE` from
/// `TASK_SIZE_LOW = DEFAULT_MAP_WINDOW`), and `malloc` passes no hint. The
/// bound of [`PseudoRtShift::HeapAddress`].
pub const HEAP_ADDRESS_END: u64 = (1 << 47) - 4096;

/// The largest value `pseudo_rt_shift * k` takes for a heap-address shift: the
/// source's product, `double` times `double(k)`, with the address at
/// [`HEAP_ADDRESS_END`]. Every address below it gives a product no larger.
fn largest_heap_shift(k: usize) -> f64 {
    f64::from_bits(HEAP_ADDRESS_END) * k as f64
}

/// Whether every heap-address shift leaves the `.dta` value `rt` of trace `k`
/// as it is: `toStr(pseudo_rt_shift * k + rt)`. [`to_str`] writes a distinct
/// text for every distinct value that a shift can reach, so the text is the
/// shift-`0` text exactly when the value is.
fn dta_absorbs(rt: f64, k: usize) -> bool {
    !rt.is_finite() || rt + largest_heap_shift(k) == rt
}

/// Whether every heap-address shift leaves the formula text of trace `k` as
/// it is: `operator<<(pseudo_rt_shift * k + centre)`, six significant digits.
fn formula_absorbs(center: f64, k: usize) -> bool {
    !center.is_finite() || ostream_g(center + largest_heap_shift(k), 6) == ostream_g(center, 6)
}

/// The value `writeFeatureDebugInfo_` reads for `pseudo_rt_shift`: its
/// `ParamValue` converted to `double` (`ParamValue::operator double`,
/// `ParamValue.cpp:397-408`).
///
/// # Errors
///
/// `Ok(Err(termination))` when the source throws inside the OpenMP region:
/// the key is absent (`ElementNotFound`) or holds no value
/// (`ConversionError`). The caller turns that into the terminated run. A
/// string or list value is [`PseudoRtShift::HeapAddress`].
pub(crate) fn read_pseudo_rt_shift(
    parameters: &crate::param::Param,
    key: &str,
) -> Result<std::result::Result<PseudoRtShift, (&'static str, String)>> {
    use crate::param::ParamValue;
    if !parameters.exists(key)? {
        return Ok(Err((
            "ElementNotFound",
            format!("the element '{key}' could not be found"),
        )));
    }
    Ok(match parameters.value(key)? {
        ParamValue::Empty => Err((
            "ConversionError",
            "Could not convert ParamValue::EMPTY to double".to_owned(),
        )),
        // `double(data_.ssize_)`: `cvtsi2sdq`, round to nearest.
        ParamValue::Integer(value) => Ok(PseudoRtShift::Value(*value as f64)),
        ParamValue::Float(value) => Ok(PseudoRtShift::Value(*value)),
        ParamValue::String(_)
        | ParamValue::StringList(_)
        | ParamValue::IntegerList(_)
        | ParamValue::FloatList(_) => Ok(PseudoRtShift::HeapAddress),
    })
}

/// The inputs of `writeFeatureDebugInfo_` besides the algorithm's state.
#[derive(Clone, Copy)]
pub struct FeatureDebugInput<'a> {
    /// The fitted model.
    pub fitter: &'a dyn TraceFitter,
    /// The extended traces, after the baseline and maximum update.
    pub traces: &'a MassTraces,
    /// The cropped traces.
    pub new_traces: &'a MassTraces,
    /// Whether the candidate passed the quality checks.
    pub feature_ok: bool,
    /// The rejection reason when it did not.
    pub error_msg: &'a str,
    /// The final score when it did.
    pub final_score: f64,
    /// The seed's `plot_nr`.
    pub plot_nr: i64,
    /// The seed peak's m/z.
    pub seed_mz: f64,
    /// The size of the output map at this point, `(*features_).size()`: the
    /// caller's features plus those of earlier charges.
    pub features_len: usize,
    /// The shift that places trace `k` at `rt + k * pseudo_rt_shift`.
    pub pseudo_rt_shift: PseudoRtShift,
    /// The directory prefix, [`FEATURE_DEBUG_PATH`] in the algorithm.
    pub path: &'a str,
}

/// `TextFile::store`: every line followed by `"\n"`, a line that already ends
/// in `"\r\n"` with that ending replaced by `"\n"`, one that ends in `"\n"`
/// unchanged.
fn text_file(lines: &[Vec<u8>]) -> Vec<u8> {
    let mut out = Vec::new();
    for line in lines {
        if line.ends_with(b"\r\n") {
            out.extend_from_slice(&line[..line.len() - 2]);
            out.push(b'\n');
        } else if line.ends_with(b"\n") {
            out.extend_from_slice(line);
        } else {
            out.extend_from_slice(line);
            out.push(b'\n');
        }
    }
    out
}

/// One `.dta` file: `toStr(shift * k + rt) + "\t" + intensity` per peak.
fn dta(traces: &MassTraces, shift: f64) -> String {
    let lines: Vec<Vec<u8>> = traces
        .iter()
        .enumerate()
        .flat_map(|(k, trace)| {
            trace.peaks.iter().map(move |peak| {
                let mut line = to_str(shift * k as f64 + peak.rt);
                line.push('\t');
                line.push_str(&to_str_f32(peak.intensity));
                line.into_bytes()
            })
        })
        .collect();
    // Every piece is ASCII.
    String::from_utf8_lossy(&text_file(&lines)).into_owned()
}

/// `writeFeatureDebugInfo_` (`FeatureFinderAlgorithmPicked.cpp:2130-2203`),
/// with the shift already read.
///
/// - `.dta`: every peak of the extended traces, trace `k` shifted by
///   `k * pseudo_rt_shift`, as `StringUtils::toStr(double)`, a tab and the
///   `float` intensity as `StringUtils::toStr(float)`.
/// - `_cropped.dta`: the same for the cropped traces, written only when they
///   hold a peak.
/// - `.plot`: one gnuplot function per extended trace, named `'f' + k` as a C++
///   `char` (so the 27th trace gets byte `0x80`), from
///   [`TraceFitter::gnuplot_formula`] with the traces' baseline and shift
///   `k * pseudo_rt_shift`; the axis labels, `set samples 1000`, the `plot`
///   command and `pause -1`. The `plot` command titles the curves `before fit
///   (RT: <centre, 2 decimals> m/z: <seed m/z, 4 decimals>)`, the cropped
///   curve `feature  - <reason>` for a rejected candidate or `feature <n>
///   (score: <3 decimals>)` for an accepted one, where `n` is the output map's
///   size plus one at that moment, and each trace `Trace <k> (m/z: <average,
///   4 decimals>)`.
///
/// # Errors
///
/// Returns [`Error::Unsupported`] for a [`PseudoRtShift::HeapAddress`] shift
/// at the first value, in the order the source writes them, whose text the
/// address can change: a peak of trace `k >= 1` in the extended or the
/// cropped traces, or the fitted centre in the formula of trace `k >= 1`. The
/// source writes a text there that depends on the heap address and differs
/// from process to process. Every other text of such a shift is the text of
/// shift `0`, which is what the source writes.
pub fn write_feature_debug_info(input: FeatureDebugInput<'_>) -> Result<FeatureDebugFiles> {
    let FeatureDebugInput {
        fitter,
        traces,
        new_traces,
        feature_ok,
        error_msg,
        final_score,
        plot_nr,
        seed_mz,
        features_len,
        pseudo_rt_shift,
        path,
    } = input;
    let pseudo_rt_shift = match pseudo_rt_shift {
        PseudoRtShift::Value(value) => value,
        PseudoRtShift::HeapAddress => {
            check_absorbed(traces, new_traces, fitter.center(), plot_nr)?;
            0.0
        }
    };
    let mut script = format!(
        "plot \"{path}{plot_nr}.dta\" title 'before fit (RT: {} m/z: {})' with points 1",
        number(fitter.center(), 2),
        number(seed_mz, 4)
    )
    .into_bytes();
    let before = dta(traces, pseudo_rt_shift);
    let mut cropped = None;
    if new_traces.peak_count() != 0 {
        cropped = Some(dta(new_traces, pseudo_rt_shift));
        script.extend_from_slice(
            format!(", \"{path}{plot_nr}_cropped.dta\" title 'feature ").as_bytes(),
        );
        if feature_ok {
            script.extend_from_slice(
                format!(
                    "{} (score: {})",
                    features_len as u64 + 1,
                    number(final_score, 3)
                )
                .as_bytes(),
            );
        } else {
            script.extend_from_slice(format!(" - {error_msg}").as_bytes());
        }
        script.extend_from_slice(b"' with points 3");
    }
    let mut lines: Vec<Vec<u8>> = Vec::with_capacity(traces.len() + 5);
    for (k, trace) in traces.iter().enumerate() {
        // `char fun = 'f'; fun += (char)k;`: arithmetic modulo 256.
        let name = (usize::from(b'f') + k) as u8;
        let mut formula = fitter
            .gnuplot_formula(trace, 'f', traces.baseline, pseudo_rt_shift * k as f64)
            .into_bytes();
        // The formula starts with the one-byte name 'f'.
        if let Some(first) = formula.first_mut() {
            *first = name;
        }
        lines.push(formula);
        script.extend_from_slice(b", ");
        script.push(name);
        script.extend_from_slice(
            format!("(x) title 'Trace {k} (m/z: {})'", number(trace.avg_mz(), 4)).as_bytes(),
        );
    }
    lines.push(b"set xlabel \"pseudo RT (mass traces side-by-side)\"".to_vec());
    lines.push(b"set ylabel \"intensity\"".to_vec());
    lines.push(b"set samples 1000".to_vec());
    lines.push(script);
    lines.push(b"pause -1".to_vec());
    Ok(FeatureDebugFiles {
        plot_nr,
        path: path.to_owned(),
        dta: before,
        cropped_dta: cropped,
        plot: text_file(&lines),
    })
}

/// The check of [`write_feature_debug_info`] for a heap-address shift, in the
/// order the source writes the values: the `.dta` peaks, the `_cropped.dta`
/// peaks, then the centre in each trace formula of the `.plot` file. Trace 0
/// is shifted by `shift * 0 = 0` and never depends on the address.
fn check_absorbed(
    traces: &MassTraces,
    new_traces: &MassTraces,
    center: f64,
    plot_nr: i64,
) -> Result<()> {
    let refuse = |file: &str, trace: usize, what: &str, value: f64| {
        Error::Unsupported(format!(
            "write_debug: debug:pseudo_rt_shift holds a string or a list, so the source shifts \
             trace k by k times the bits of a heap address; in {plot_nr}{file} trace {trace} has \
             the {what} {value:e}, where such a shift changes the written number, which then \
             differs from process to process"
        ))
    };
    for (file, set) in [(".dta", traces), ("_cropped.dta", new_traces)] {
        for (k, trace) in set.iter().enumerate().skip(1) {
            if let Some(peak) = trace.peaks.iter().find(|peak| !dta_absorbs(peak.rt, k)) {
                return Err(refuse(file, k, "retention time", peak.rt));
            }
        }
    }
    if let Some(k) = (1..traces.len()).find(|&k| !formula_absorbs(center, k)) {
        return Err(refuse(".plot", k, "fitted centre", center));
    }
    Ok(())
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Finds mass spectrometric peaks in profile mass spectra.
//!
//! Ports `OpenMS4-topp/src/PeakPickerHiRes.cpp` (topp `174b576`), the TOPP
//! wrapper of the high-resolution peak picker
//! [`crate::processing::peak_picking::PeakPickerHiRes`]. The in-memory mode is
//! ported: load the mzML input, warn once about per-peak ion mobility, refuse an
//! input without spectra and chromatograms or with unsorted records, pick every
//! selected spectrum and every chromatogram, attach one shared `peak picking`
//! processing record and store the mzML output. The tool's `algorithm`
//! subsection is the picker's source parameter tree.
//!
//! The input is read with [`PeakPickerHiRes::read_options`], the size-derived
//! library ceilings plus the source-compatibility switches a tool path takes
//! (decision D10). Instrument-sized profile runs read: the support document
//! records a 2.3 GB, 40,856-spectrum Q Exactive file picked end to end and
//! compared against the C++ tool.
//!
//! `-processOption lowmemory` is ported too: the source's `PPHiResMzMLConsumer`
//! driven by `MzMLFile::transform`, which streams the input past a picking
//! consumer that writes each record as it is produced, so no experiment is ever
//! held. That mode is not the in-memory mode with a smaller footprint - it
//! omits five of the in-memory checks, decides automatic mode on the stored
//! spectrum type alone, and can therefore write a different file from the same
//! input. [`LowMemoryPicker`] and [`run_low_memory`](fn@run_low_memory) carry
//! the list and the evidence.
//!
//! Not ported yet: progress logging.
//!
//! `docs/TOPP_PEAK_PICKER_HI_RES_SUPPORT.md` lists the source members, the
//! preserved conventions, the native differences and the evidence.

use crate::cli::{ExitCode, Tool, ToolContext, ToolError, ToolResult, ToolSpec};
use crate::format::PeakFileOptions;
use crate::format::file_handler::FileHandler;
use crate::format::file_types::FileType;
use crate::format::ms_data_writing_consumer::{
    MSDataWritingConsumer, MSDataWritingProcessor, ReferencePolicy,
};
use crate::format::mzml;
use crate::kernel::{MSChromatogram, MSExperiment, MSSpectrum, SpectrumType};
use crate::metadata::{
    DataProcessing, ImTypes, IonMobilityFormat, IonMobilityPeakType, MetaValue, MetaValueData,
    ProcessingAction, im_peak_type_to_string,
};
use crate::param::Param;
use crate::processing::peak_picking::{
    CENTROIDED_INPUT_MESSAGE, PeakPickerHiRes as Picker, PickedExperimentReport,
    PickingCompatibility,
};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::Write;

/// The `PeakPickerHiRes` TOPP tool: centroids profile spectra and chromatograms
/// with the high-resolution peak picker.
///
/// Registration follows `registerOptionsAndFlags_` (`OpenMS4-topp/src/PeakPickerHiRes.cpp:152-163`):
/// `-in` and `-out` restricted to mzML, the advanced `-processOption` with the
/// valid strings `inmemory` (the default) and `lowmemory`, and the `algorithm`
/// subsection filled from
/// [`PeakPickerHiRes::defaults`](crate::processing::peak_picking::PeakPickerHiRes::defaults),
/// as `getSubsectionDefaults_` returns `PeakPickerHiRes().getDefaults()`. An INI
/// the C++ tool writes with `-write_ini` is accepted unchanged, and so is every
/// `-algorithm:<name>` command-line override the strict parameter update
/// accepts. See [`PeakPickerHiRes::run_io`] for the run.
pub struct PeakPickerHiRes;

impl PeakPickerHiRes {
    /// The mzML load options of this tool path: [`mzml::ReadOptions::source`].
    ///
    /// Source `main_` loads through `FileHandler::loadExperiment`, whose
    /// `MzMLFile` has neither a source-compatibility switch nor a resource
    /// ceiling, so a tool path reproduces it in both respects:
    ///
    /// * **Source compatibility (decision D10).** `source_dangling_references`
    ///   drops a `softwareRef` or `dataProcessingRef` that names no definition,
    ///   which source `std::map::operator[]` does silently
    ///   (`MzMLHandler.cpp:920-952`); the library default refuses it.
    ///   `source_invalid_timestamps` keeps an unparseable `startTimeStamp`
    ///   non-fatal, which is both the source behaviour and the library default.
    ///   Using `ReadOptions::source()` rather than naming the switches here
    ///   means a switch added to that constructor reaches this tool with it.
    /// * **Resource ceilings.** The limits are the library defaults, which are
    ///   size-derived ([`mzml::InputScaling`], `src/format/mzml_scaling.rs`):
    ///   the ceiling for every cumulative quantity grows with the XML bytes already
    ///   consumed, so work and storage stay linear in the input while a
    ///   document of any realistic size fits. The ceilings this tool shipped
    ///   with were fixed, and refused instrument-sized profile data outright: a
    ///   2.3 GB, 40,856-spectrum Q Exactive run has 197,765,338 raw points
    ///   against a fixed ceiling of 10,000,000, and 2.3 GB of XML against a
    ///   fixed 512 MiB. `docs/TOPP_PEAK_PICKER_HI_RES_SUPPORT.md` measures that
    ///   input through this tool and against the C++ tool.
    ///
    /// Scientific selection stays at [`PeakFileOptions::default`], as the source
    /// tool sets no `PeakFileOptions`; that default still sorts every record by
    /// position, as source `MzMLHandler` does.
    pub fn read_options() -> mzml::ReadOptions {
        mzml::ReadOptions::source()
    }
}

/// Source warning for per-peak ion mobility (`OpenMS4-topp/src/PeakPickerHiRes.cpp:222-226`),
/// around the name of the spectrum's ion mobility peak type.
fn ion_mobility_warning(peak_type: IonMobilityPeakType) -> String {
    format!(
        "Warning: Input contains per-peak ion mobility data (IM_PEAK, {}). PeakPickerHiRes picks in m/z only and reports intensity-weighted mean ion mobility. This produces incorrect results on unbinned data. Consider IonMobilityBinning or PeakPickerIM first.",
        im_peak_type_to_string(peak_type)
    )
}

/// Source warning for an input without spectra and chromatograms
/// (`OpenMS4-topp/src/PeakPickerHiRes.cpp:233-234`). The log stream ends the line.
const EMPTY_INPUT_WARNING: &str = "The given file does not contain any conventional peak data, but might contain chromatograms. This tool currently cannot handle them, sorry.";

/// Source error for an unsorted spectrum (`OpenMS4-topp/src/PeakPickerHiRes.cpp:243`).
const UNSORTED_SPECTRA_ERROR: &str = "Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!";

/// Source error for an unsorted chromatogram (`OpenMS4-topp/src/PeakPickerHiRes.cpp:253`),
/// which also says m/z.
const UNSORTED_CHROMATOGRAMS_ERROR: &str = "Error: Not all chromatograms are sorted according to peak m/z positions. Use FileFilter to sort the input!";

/// The per-record hook of `-processOption lowmemory`: source nested class
/// `PPHiResMzMLConsumer` (`OpenMS4-topp/src/PeakPickerHiRes.cpp:107-150`).
///
/// The source class derives from `MSDataWritingConsumer` and overrides its two
/// template-method hooks; this port is the hook pair alone
/// ([`MSDataWritingProcessor`]), handed to
/// [`MSDataWritingConsumer`], which is that base class. It carries nothing but
/// the picker: the source keeps a second copy of `ms_levels` read out of
/// `pp.getParameters()`, which is by construction the picker's own
/// `ms_levels_`, so this reads [`Picker::ms_levels`] and the two cannot drift.
///
/// # How the selection differs from the in-memory mode
///
/// The source spells the automatic-mode test here as `s.getType()`, the
/// **stored** spectrum type: `MSSpectrum` re-exposes the base accessor with
/// `using SpectrumSettings::getType` (`MSSpectrum.h:655`), and the no-argument
/// overload is that base one, which returns the `type_` member and nothing
/// else. `pickExperiment`, the in-memory path, spells the same test
/// `getType(true)` (`CENTROIDING/PeakPickerHiRes.cpp:510` and `531`) - stored type, then a
/// scan of the record's data-processing history for a `PEAK_PICKING` action,
/// then `PeakTypeEstimator` over the samples.
///
/// So a spectrum whose stored type is `UNKNOWN` - which is every spectrum whose
/// mzML carries no `MS:1000127`/`MS:1000128` term, including everything written
/// by a converter that only sets `MS:1000525` - is **picked by the low-memory
/// mode and copied by the in-memory mode** whenever the slower test would have
/// called it centroided. That is not a rounding difference; it is a different
/// output file. It is reproduced here, not repaired: the source's low-memory
/// path is the specification for the low-memory path.
///
/// The manual-mode arm diverges in the same direction: `pickExperiment` refuses
/// a centroided spectrum on a selected MS level unless `-force` is given, and
/// this class runs `pp_.pick` straight away. **`-force` is inert in the
/// low-memory mode**, in the source and here, because
/// [`Picker::check_spectrum_type`] is read only by the experiment entry points
/// and this hook calls the single-record [`Picker::pick_spectrum`].
struct LowMemoryPicker {
    /// The configured picker, source member `pp_`.
    picker: Picker,
}

impl MSDataWritingProcessor for LowMemoryPicker {
    /// Source `processSpectrum_` (`OpenMS4-topp/src/PeakPickerHiRes.cpp:120-137`).
    ///
    /// Automatic mode (no `ms_levels`) leaves a spectrum whose **stored** type
    /// is centroided untouched; manual mode leaves every spectrum whose MS level
    /// is not selected untouched. Everything else is picked, with no type check.
    /// An untouched spectrum is still written, and still receives the tool's
    /// `peak picking` processing record, because the consumer appends that to
    /// every record after this hook returns.
    ///
    /// # Errors
    ///
    /// [`Picker::pick_spectrum`]'s errors. The source's hook returns `void` and
    /// can only throw; a throw there escapes through `MzMLFile::transform` and
    /// ends the run with a partially written output file, which is what this
    /// error does too - see [`run_low_memory`](fn@run_low_memory).
    fn process_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<()> {
        if self.picker.ms_levels.is_empty() {
            if spectrum.spectrum_type == SpectrumType::Centroid {
                return Ok(());
            }
        } else if !self.picker.ms_levels.contains(&spectrum.ms_level) {
            return Ok(());
        }
        *spectrum = self.picker.pick_spectrum(spectrum)?.spectrum;
        Ok(())
    }

    /// Source `processChromatogram_` (`OpenMS4-topp/src/PeakPickerHiRes.cpp:139-144`): every
    /// chromatogram is picked, unconditionally, exactly as `pickExperiment`
    /// picks every chromatogram in the in-memory mode. `ms_levels` does not
    /// apply to chromatograms in either mode.
    ///
    /// # Errors
    ///
    /// [`Picker::pick_chromatogram`]'s errors.
    fn process_chromatogram(&mut self, chromatogram: &mut MSChromatogram) -> Result<()> {
        *chromatogram = self.picker.pick_chromatogram(chromatogram)?.chromatogram;
        Ok(())
    }
}

/// Source `doLowMemAlgorithm` (`OpenMS4-topp/src/PeakPickerHiRes.cpp:170-186`).
///
/// Builds the writing consumer on `output`, gives it the `peak picking`
/// processing record, and streams `input` past it with
/// [`mzml::transform_with_options`], the port of `MzMLFile::transform`.
///
/// # What the mode does, in order
///
/// 1. **The output file is created before the input is read**, because the
///    source consumer opens its `std::ofstream` in its constructor
///    (`MSDataWritingConsumer.cpp:33`) and the constructor runs before
///    `transform`. An input that cannot be parsed therefore still leaves a file
///    behind - empty, if the failure comes before the first record. The
///    in-memory mode writes nothing until the whole run has succeeded.
/// 2. **The input is read twice.** `transform` runs `transformFirstPass_` and
///    then a second full parse (`MzMLFile.cpp:178-191`). The first pass hands
///    the consumer the record counts declared by the document and the
///    experimental settings; the second hands it the records. Both passes are
///    complete parses of the file - the mode trades I/O for memory, and a
///    low-memory run reads roughly twice the bytes an in-memory run does.
/// 3. **The header is written from the settings of pass one plus the first
///    record**, and each list tag announces the count pass one declared, not
///    the records that follow. A header field derived from the records
///    therefore describes that one record: on the 2.3 GB benchmark input the
///    `fileContent` of a low-memory output carries `MS1 spectrum` alone where
///    the in-memory output carries `MS1 spectrum` and `MSn spectrum`. **The C++
///    Release build does the same** on the same input - it is the one-record
///    dummy map of `MSDataWritingConsumer.cpp:76-83` - and apart from that one
///    line the two modes' mzML bodies are byte-identical on both sides
///    (`../oracle/p4-lowmemory`, `logs/verify2_06.log`).
///    The source notes that the list counts are not
///    enforced and that a wrong one "will lead to an inconsistent mzML".
///    [`CountPolicy::Checked`](crate::format::ms_data_writing_consumer::CountPolicy::Checked),
///    the consumer's default, is kept: a document whose declared counts and
///    actual records disagree ends this run with [`Error::InvalidValue`] after
///    the output has been closed, where the source silently writes an mzML
///    whose `count` attributes lie. Native difference, recorded in the support
///    document.
/// 4. **A record needing header entries the first record did not contribute
///    is written with the source's dangling reference.** The header is the
///    first record's, so a later record's `dataProcessing` cannot be numbered
///    against it; the source numbers the reference by the record's position in
///    the stream instead, which names nothing the header declares. This path
///    selects
///    [`ReferencePolicy::SourceDangling`](crate::format::ms_data_writing_consumer::ReferencePolicy::SourceDangling)
///    to reproduce that exactly, because the alternative - the library
///    default, which refuses the record - would stop the mode on any
///    `FileMerger` output, where every merged part carries its own
///    `dataProcessing`. Native difference, recorded in the support document
///    with the C++ issue the dangling reference itself deserves.
/// 5. **Each record is picked and written immediately**, then dropped. Peak
///    memory is one read batch
///    ([`PeakFileOptions::max_data_pool_size`], 100 records, as upstream) plus
///    the rendered text of one record, not the experiment.
///
/// # What the mode does *not* do
///
/// None of the in-memory mode's four input phases exists on this path, in the
/// source or here: no per-peak ion mobility warning, no
/// [`ExitCode::IncompatibleInputData`] for an input without spectra and
/// chromatograms, no sortedness refusal, and no per-MS-level summary on stdout.
/// An empty input produces an empty output file and exit 0.
///
/// # Threads
///
/// The mode is serial, in the source and here. The source's consumer dispatch
/// loop calls `consumeSpectrum` one record at a time
/// (`MzMLHandler.cpp:259-274`); its only OpenMP region on this path decodes
/// binary arrays, which this port's reader does not parallelise either. So
/// `-threads` reaches nothing here and the written bytes are identical at every
/// value of it - trivially, rather than by the batch-order argument the
/// in-memory mode needs. The test suite pins that at 1, 8 and 32.
///
/// # Errors
///
/// [`Error::Io`] when `output` cannot be created; the read errors of
/// [`mzml::transform_with_options`] under [`PeakPickerHiRes::read_options`];
/// the picker's errors from the hooks; and the consumer's own refusals -
/// a duplicate native identifier, a record whose header references cannot be
/// satisfied by the one header already written, and the count mismatch of
/// point 3.
///
/// Whichever of them ends the run, the document is closed before the error is
/// returned. The source's `~MSDataWritingConsumer` calls `doCleanup_` on every
/// path (`MSDataWritingConsumer.cpp:37-40`), which closes the open list and
/// writes the footer whenever writing started (`:151-173`), so a failed source
/// run leaves a closed, indexed document holding the records it got to.
/// [`MSDataWritingConsumer::finish`](crate::format::ms_data_writing_consumer::MSDataWritingConsumer::finish)
/// is this port's destructor - Rust cannot report an error from a drop - so
/// the failing path calls it too and discards its own result, which keeps the
/// failure that caused it rather than the count mismatch that a half-written
/// run raises by construction. The records already written stay where they
/// are, under the `count` attributes pass one declared: a streaming writer
/// cannot take bytes back. A failure before the first record leaves the file
/// empty, because `doCleanup_` writes nothing while `started_writing_` is
/// false.
///
/// `run_io` classifies the error as `TOPPBase::main` classifies the exceptions
/// this path lets through uncaught - `doLowMemAlgorithm` catches nothing
/// (`OpenMS4-topp/src/PeakPickerHiRes.cpp:170-186`). A reader failure is
/// `Error: Unable to read file (<reason>)` with
/// [`ExitCode::InputFileCorrupt`], as the source's `ParseError` arm
/// (`TOPPBase.cpp:460-465`) and as in the in-memory mode, whose loader errors
/// propagate the same way; [`Error::Unsupported`] is
/// `INCOMPATIBLE_INPUT_DATA`. Only [`Error::InvalidValue`],
/// [`Error::InvalidRange`] and [`Error::MissingInformation`] are reported as
/// `Error: Unexpected internal error (<reason>)` with
/// [`ExitCode::UnknownError`], because the source exceptions they stand for
/// derive straight from `BaseException` and take that arm there, while this
/// port's framework mapping reserves parameter codes for them.
///
/// One error has no source counterpart and no agreement between the two modes:
/// the transform's own administrative ceiling
/// ([`mzml::TransformOptions::max_bytes`] and `max_work`) is an
/// [`Error::InvalidValue`], so it takes the `Unexpected internal error` arm
/// here, where the in-memory mode's loader ceilings reach the framework's
/// `InvalidValue` mapping and exit `ILLEGAL_PARAMETERS`. The source has no
/// ceiling of either kind, so neither code is its answer, and nothing
/// distinguishes this error from the picker's by kind alone.
fn run_low_memory(
    ctx: &ToolContext,
    input: &str,
    output: &str,
    picker: Picker,
) -> Result<ExitCode> {
    let mut processing = ctx.processing_info(&[ProcessingAction::PeakPicking])?;
    render_list_parameters(&mut processing);
    let mut consumer = MSDataWritingConsumer::create(output, LowMemoryPicker { picker })?
        .with_reference_policy(ReferencePolicy::SourceDangling);
    consumer.add_data_processing(processing)?;
    let options = mzml::TransformOptions {
        load: mzml::LoadOptions {
            scientific: PeakFileOptions::default(),
            ..mzml::LoadOptions::default()
        },
        read: PeakPickerHiRes::read_options(),
        ..mzml::TransformOptions::default()
    };
    match mzml::transform_with_options(input, &mut consumer, &options) {
        Ok(_) => {
            consumer.finish()?;
            Ok(ExitCode::ExecutionOk)
        }
        Err(error) => {
            let _ = consumer.finish();
            Err(error)
        }
    }
}

/// The in-memory input checks of `main_` before picking
/// (`OpenMS4-topp/src/PeakPickerHiRes.cpp:216-256`), in source order.
///
/// Writes the ion mobility warning once, for the first spectrum whose format
/// is per-peak, then returns the terminal exit code, if any, after writing its
/// diagnostic: [`ExitCode::IncompatibleInputData`] for an experiment without
/// spectra and chromatograms, or with an unsorted spectrum or chromatogram.
///
/// The lines reach `err` here, where the source writes them — before picking
/// starts — and not at the end of the run, so a diagnostic survives whatever
/// the rest of the run does, including an error that propagates. That is only
/// possible because this runs on the calling thread: see
/// [`pick_experiment`](fn@pick_experiment) for what does not.
///
/// The source prints the stored ion mobility peak type
/// (`imPeakTypeToString(spec.getIMPeakType())`). The native spectrum has no
/// stored peak type; the name printed is `im_profile`, which is what the source
/// mzML reader stores for ion mobility data without the `MS:1003441` term
/// (`MzMLHandler.cpp:253-255`). An input carrying that term reads
/// `im_centroided` in the source.
///
/// Both sortedness checks are unreachable through the tool's loader, which
/// sorts every record by position, as the source `FileHandler::loadExperiment`
/// does with default `PeakFileOptions` (`MzMLHandler.cpp:218-221`, `299-302`);
/// they are kept because the source keeps them.
///
/// The two sortedness errors are `writeLogError_` lines in the source, so
/// they also reach the `-log` file ([`ToolContext::write_log_error`]); the
/// ion mobility and empty-input warnings are `OPENMS_LOG_WARN` lines, which
/// do not. Every refusal here is an exit code `main_` returns, so the
/// lifecycle's closing `PeakPickerHiRes took … .` line follows it.
///
/// # Errors
///
/// Returns [`Error::Io`] when a line cannot be written to `err`.
fn check_input(
    ctx: &ToolContext,
    experiment: &MSExperiment,
    err: &mut dyn Write,
) -> Result<Option<ExitCode>> {
    if experiment
        .spectra
        .iter()
        .any(|spectrum| ImTypes::determine_im_format(spectrum) == IonMobilityFormat::PerPeak)
    {
        writeln!(
            err,
            "{}",
            ion_mobility_warning(IonMobilityPeakType::Profile)
        )?;
    }
    if experiment.spectra.is_empty() && experiment.chromatograms.is_empty() {
        writeln!(err, "{EMPTY_INPUT_WARNING}")?;
        return Ok(Some(ExitCode::IncompatibleInputData));
    }
    if !experiment
        .spectra
        .iter()
        .all(|spectrum| spectrum.is_sorted())
    {
        ctx.write_log_error(err, UNSORTED_SPECTRA_ERROR)?;
        return Ok(Some(ExitCode::IncompatibleInputData));
    }
    if !experiment
        .chromatograms
        .iter()
        .all(|chromatogram| chromatogram.is_sorted())
    {
        ctx.write_log_error(err, UNSORTED_CHROMATOGRAMS_ERROR)?;
        return Ok(Some(ExitCode::IncompatibleInputData));
    }
    Ok(None)
}

/// Pick every selected spectrum and every chromatogram of `experiment` in
/// place, as `PeakPickerHiRes::pickExperiment`, on the `-threads` pool.
///
/// The run's policy
/// ([`ToolContext::thread_policy`](crate::cli::ToolContext::thread_policy))
/// reaches
/// [`PeakPickerHiRes::pick_experiment_in_place_with_threads`](crate::processing::peak_picking::PeakPickerHiRes::pick_experiment_in_place_with_threads):
/// the spectrum loop runs on the pool
/// [`ToolContext::in_thread_pool`](crate::cli::ToolContext::in_thread_pool)
/// sizes from that policy, and the centroids are bit-identical at every worker
/// count (the determinism contract of [`crate::concept::parallel`]). The source
/// picks serially — `PeakPickerHiRes.cpp` has no OpenMP — so there is no C++
/// parallel baseline here; what the source does parallelise for this workload
/// is the mzML reader (`MzMLHandler.cpp:206`), which this tool's loader does
/// not.
///
/// **The pool wraps this call, not the whole tool body**, which is where the
/// five earlier ported tools put it. Picking is this tool's only parallel
/// region: reading, the input checks, the summary, the processing record and
/// storing are serial in this port and in the source. Scoping the pool to the
/// region that uses it has three consequences, all wanted:
///
/// * the diagnostics and the summary are written where the source writes them,
///   to the real streams, because the phases that produce them run on the
///   calling thread — `run_io`'s streams are not [`Send`] and so cannot cross
///   onto a pool thread;
/// * the mzML read, which is about 45% of an instrument-scale run and the
///   heaviest allocator client in it, keeps the calling thread's malloc arena;
/// * the pool exists only while it is used.
///
/// **At one worker no pool is built** and picking runs on the calling thread.
/// A pool of one worker is a pure cost: glibc gives the worker a second malloc
/// arena, which measures at +0.6 s on the 2.3 GB benchmark run — a net
/// regression at `-threads 1`, which is the TOPP default. Nothing is given up,
/// because the bound a pool provides is a bound on rayon work and the only
/// rayon work here is this call, whose width is the `Threads` value passed to
/// it: at one worker the picker's batch loop maps on the calling thread and
/// builds no pool of its own. A second parallel region added to this tool
/// belongs inside the pool the same way this one is.
///
/// The **in-place** entry point is used rather than the borrowing one because
/// this tool writes the picked experiment and never reads the profile data
/// again: picking in place releases each spectrum's profile samples as its
/// centroids appear and copies no record it does not pick, where the borrowing
/// form holds a second experiment beside the first and clones every unpicked
/// record. The two produce the same experiment and the same reports; only the
/// peak memory differs. The in-place form is not atomic, which costs this tool
/// nothing: its only reaction to a picking error is to report it and exit
/// without writing an output file.
///
/// # Errors
///
/// The picker's errors, and [`Error::Io`] when the operating system refuses the
/// worker threads. `run_io` reports both the same way — as
/// `Error: Unexpected internal error (<reason>)` with
/// [`ExitCode::UnknownError`], see its `# Errors` — because this call, and with
/// it the pool, is inside the body rather than around it.
fn pick_experiment(
    ctx: &ToolContext,
    picker: &Picker,
    experiment: &mut MSExperiment,
) -> Result<PickedExperimentReport> {
    let threads = ctx.thread_policy();
    if threads.get() <= 1 {
        return picker.pick_experiment_in_place_with_threads(experiment, threads);
    }
    ctx.in_thread_pool(|| picker.pick_experiment_in_place_with_threads(experiment, threads))?
}

/// The per-level summary `pickExperiment` logs (`CENTROIDING/PeakPickerHiRes.cpp:559-563`):
/// for each MS level in ascending order, the spectra picked and the spectra
/// seen. The header is written even without spectra.
///
/// `experiment` is the picked experiment, whose spectra are those of the input
/// in input order and carry the MS level they were read with: picking replaces a
/// record's samples and leaves its metadata, so the per-level counts are the
/// same before and after.
fn pick_summary(experiment: &MSExperiment, report: &PickedExperimentReport) -> Vec<String> {
    let mut levels: BTreeMap<u32, (u64, u64)> = BTreeMap::new();
    for (spectrum, boundaries) in experiment.spectra.iter().zip(&report.spectrum_boundaries) {
        let entry = levels.entry(spectrum.ms_level).or_insert((0, 0));
        entry.0 += u64::from(boundaries.is_some());
        entry.1 += 1;
    }
    let mut summary = vec!["#Spectra that needed to and could be picked by MS-level:".to_owned()];
    for (level, (count, total)) in levels {
        summary.push(format!("  MS-level {level}: {count} / {total}"));
    }
    summary
}

/// Render the list-valued and empty parameters of a processing record as the
/// source mzML writer renders them.
///
/// `getProcessingInfo_` records every resolved parameter as a `DataValue`
/// outside `-test` (`TOPPBase.cpp:556-568`), and the source mzML writer turns a
/// list one into text: the oracle's non-test output carries
/// `<userParam name="parameter: algorithm:ms_levels" type="xsd:string" value="[]"/>`
/// (P3 oracle `notest`). This port's mzML writer refuses list and empty
/// metadata instead (`Empty/list metadata has no lossless mzML scalar
/// encoding`), so without this a run outside `-test` could not store its own
/// output: `ms_levels` is an empty integer list by default. The rendering uses
/// [`MetaValue`]'s source-style list format, which reproduces the C++ text.
///
/// This belongs in `cli::processing::processing_info`, which every tool shares,
/// or in the writer; both are outside this package, so it is an integrator
/// request and is recorded as a native difference in
/// `docs/TOPP_PEAK_PICKER_HI_RES_SUPPORT.md`.
fn render_list_parameters(processing: &mut DataProcessing) {
    for value in processing.metadata.values_mut() {
        if matches!(
            value.data(),
            MetaValueData::Empty
                | MetaValueData::StringList(_)
                | MetaValueData::IntegerList(_)
                | MetaValueData::FloatList(_)
        ) {
            *value = MetaValue::from(value.to_string());
        }
    }
}

impl Tool for PeakPickerHiRes {
    const NAME: &'static str = "PeakPickerHiRes";
    const DESCRIPTION: &'static str = "Finds mass spectrometric peaks in profile mass spectra.";

    fn register(spec: &mut ToolSpec) -> Result<()> {
        spec.register_input_file(
            "in",
            "<file>",
            "",
            "input profile data file ",
            true,
            false,
            &[],
        )?;
        spec.set_valid_formats("in", &["mzML"])?;
        spec.register_output_file("out", "<file>", "", "output peak file ", true, false)?;
        spec.set_valid_formats("out", &["mzML"])?;
        spec.register_string_option(
            "processOption",
            "<name>",
            "inmemory",
            "Whether to load all data and process them in-memory or whether to process the data on the fly (lowmemory) without loading the whole file into memory first",
            false,
            true,
        )?;
        spec.set_valid_strings("processOption", &["inmemory", "lowmemory"])?;
        spec.register_subsection("algorithm", "Algorithm parameters section")?;
        Ok(())
    }

    /// The picker's source parameter tree for every section, as
    /// `getSubsectionDefaults_` returns `PeakPickerHiRes().getDefaults()`
    /// whatever section it is asked for.
    fn subsection_defaults(_section: &str) -> Result<Option<Param>> {
        Ok(Some(Picker::defaults()?))
    }

    /// Run with the process's standard streams; see [`PeakPickerHiRes::run_io`].
    fn run(ctx: &ToolContext) -> ToolResult {
        Self::run_io(ctx, &mut std::io::stdout(), &mut std::io::stderr())
    }

    /// Source `main_` (`OpenMS4-topp/src/PeakPickerHiRes.cpp:188-273`).
    ///
    /// 1. The `algorithm` subsection configures the picker, as
    ///    `setParameters`, with every source behaviour of
    ///    [`PickingCompatibility::source`](crate::processing::peak_picking::PickingCompatibility::source)
    ///    enabled, because the source picks duplicate positions, negative
    ///    intensities and non-positive spline maxima without refusing them.
    /// 2. `-processOption lowmemory` leaves for `run_low_memory` in this
    ///    module before anything is read, as `main_` returns from
    ///    `doLowMemAlgorithm`. Steps 3 to 5 below are the in-memory mode only,
    ///    and the low-memory mode runs none of them.
    /// 3. The input is loaded as mzML with [`PeakPickerHiRes::read_options`] —
    ///    the source's dangling header references accepted (decision D10) and
    ///    the size-derived library ceilings, which admit an instrument-sized
    ///    run — then checked as described at `check_input`.
    /// 4. Spectra are picked with the spectrum type checked unless `-force` is
    ///    given (`check_spectrum_type = !force`). A centroided spectrum in manual
    ///    mode without `-force` ends the run with
    ///    `Error: Unexpected internal error (Error: Centroided data provided but profile spectra expected.)`
    ///    and [`ExitCode::UnknownError`], as the source's `IllegalArgument`
    ///    takes the `BaseException` arm of `TOPPBase::main`
    ///    ([`ToolError::unexpected`]): the exception unwinds past the closing
    ///    `PeakPickerHiRes took … .` line, so there is none, and the message
    ///    reaches the `-log` file. Nothing is written.
    /// 5. The per-level summary goes to `out`, one shared `peak picking`
    ///    processing record is attached to every spectrum and chromatogram
    ///    (`addDataProcessing_`), and the experiment is stored as mzML.
    ///
    /// Picking replaces the loaded experiment record by record rather than
    /// building a second one (see `pick_experiment` in
    /// this module). Neither that nor the worker pool changes a written byte:
    /// the output of a run is fixed by its input and its parameters, at every
    /// worker count and with or without the `parallel` feature.
    ///
    /// `-threads` reaches the picking, which is this tool's only parallel
    /// region, as `TOPPBase::main` applies the setting before `main_`
    /// (`TOPPBase.cpp:408-415`): `pick_experiment` in this module opens
    /// the [`ToolContext::in_thread_pool`] pool around that call. The five
    /// earlier ported tools wrap their whole body instead, from
    /// [`Tool::run`](crate::cli::Tool::run); a wrapper around `run` would never
    /// execute here, because this tool overrides `run_io`, which is what
    /// [`run_with`](crate::cli::run_with) calls. Every other phase runs on the
    /// calling thread, so `out` and `err`, which are not [`Send`], take each
    /// line where the source writes it: the input warnings before picking, the
    /// per-level summary before the output file is stored.
    ///
    /// # Errors
    ///
    /// Loading, storing and processing-record failures propagate and are
    /// mapped by the framework. A picker failure other than the centroided
    /// refusal is reported as `Error: Unexpected internal error (<reason>)`
    /// with [`ExitCode::UnknownError`] ([`ToolError::unexpected`], no closing
    /// line): those are native bounds (points and
    /// work per record, metadata copies) and the FWHM search that never
    /// terminates in the source, not parameter errors. A refusal by the
    /// operating system to start the `-threads` workers is reported the same
    /// way, as `Error: Unexpected internal error (cannot start <n> worker
    /// threads: <reason>)` with [`ExitCode::UnknownError`]: the pool is built
    /// inside the picking call now (`pick_experiment` in this module), so the
    /// [`Error::Io`] it raises reaches the same arm as a picker failure instead
    /// of propagating out of `run_io` as it did while the pool wrapped the whole
    /// body. No output file is written in either case. A low-memory run's
    /// failures, which `run_low_memory` in this module lists, are split
    /// differently: only the three kinds standing for source exceptions that
    /// derive straight from `BaseException` take that arm there, and a reader
    /// failure propagates to the framework's `ParseError` mapping exactly as
    /// this path's loader call does, so both modes answer a corrupt input with
    /// `Error: Unable to read file (<reason>)` and `INPUT_FILE_CORRUPT`. What
    /// the low-memory mode cannot do is take back the bytes it has already
    /// streamed; it closes the document over them instead. The exception is
    /// [`Error::Unsupported`], which propagates (`INCOMPATIBLE_INPUT_DATA`):
    /// the picker returns it for `SignalToNoise:auto_mode` 1 when noise
    /// estimation runs on a record outside the narrow input domain where that
    /// mode is defined, which every real spectrum is; the source writes out of
    /// bounds there and crashes (oracle `PPHR_auto_mode_1`, SIGBUS or SIGSEGV).
    /// With
    /// `signal_to_noise` 0 the estimator never runs in either implementation
    /// and the run succeeds.
    fn run_io(ctx: &ToolContext, out: &mut dyn Write, err: &mut dyn Write) -> ToolResult {
        let input = ctx.string("in")?;
        let output = ctx.string("out")?;
        let process_option = ctx.string("processOption")?;

        // `writeDebug_("Parameters passed to PeakPickerHiRes", pepi_param, 3)`
        // (`PeakPickerHiRes.cpp:200-201`): the -log file, from debug level 3.
        let algorithm = ctx.subsection("algorithm")?;
        ctx.write_debug_param("Parameters passed to PeakPickerHiRes", &algorithm, 3);
        let mut picker = Picker::from_param(&algorithm)?;
        picker.compatibility = PickingCompatibility::source();

        if process_option == "lowmemory" {
            return match run_low_memory(ctx, input, output, picker) {
                Ok(code) => Ok(code),
                // `doLowMemAlgorithm` catches nothing, so every failure of
                // this path is classified by `TOPPBase::main` on the exception
                // type alone, wherever it was raised. The three kinds below
                // stand for source exceptions that derive straight from
                // `BaseException` and therefore take its `UNKNOWN_ERROR` arm
                // (`TOPPBase.cpp:495-499`) where `run_failure` would map them
                // to a parameter code. Everything else - a reader failure
                // above all - propagates to `run_failure`, whose arms are that
                // same catch chain.
                Err(
                    error @ (Error::InvalidValue(_)
                    | Error::InvalidRange(_)
                    | Error::MissingInformation(_)),
                ) => Err(ToolError::unexpected(error)),
                Err(error) => Err(error.into()),
            };
        }

        let mut experiment = FileHandler::load_experiment_with_read_options(
            input,
            &[FileType::MzMl],
            &PeakFileOptions::default(),
            &Self::read_options(),
        )?;
        if let Some(code) = check_input(ctx, &experiment, err)? {
            return Ok(code);
        }
        // `!getFlag_("force")` where the source reads it, before picking; the
        // accessor writes its `Value of … option` line to the -log file.
        picker.check_spectrum_type = !ctx.flag("force")?;

        let report = match pick_experiment(ctx, &picker, &mut experiment) {
            Ok(report) => report,
            Err(error @ Error::Unsupported(_)) => return Err(error.into()),
            // The centroided refusal is reported with the source's bare
            // message, without this port's `invalid value: ` prefix.
            Err(Error::InvalidValue(reason)) if reason == CENTROIDED_INPUT_MESSAGE => {
                return Err(ToolError::unexpected(reason));
            }
            Err(error) => return Err(ToolError::unexpected(error)),
        };
        for line in pick_summary(&experiment, &report) {
            writeln!(out, "{line}")?;
        }

        let mut processing = ctx.processing_info(&[ProcessingAction::PeakPicking])?;
        render_list_parameters(&mut processing);
        ctx.add_data_processing(&mut experiment, &processing);
        FileHandler::store_experiment(output, &experiment, Some(FileType::MzMl))?;
        Ok(ExitCode::ExecutionOk)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::{ChromatogramPeak, MSChromatogram, MSSpectrum, Peak1D};

    fn spectrum(positions: &[f64]) -> MSSpectrum {
        MSSpectrum {
            peaks: positions.iter().map(|&mz| Peak1D::new(mz, 1.0)).collect(),
            ..MSSpectrum::default()
        }
    }

    fn chromatogram(positions: &[f64]) -> MSChromatogram {
        MSChromatogram {
            peaks: positions
                .iter()
                .map(|&rt| ChromatogramPeak::new(rt, 1.0))
                .collect(),
            ..MSChromatogram::default()
        }
    }

    /// The exit code, the error-stream text and the `-log` file lines (their
    /// time stamps taken off) of [`check_input`], run with a log file.
    fn checked(experiment: &MSExperiment) -> (Option<ExitCode>, String, Vec<String>) {
        let dir = crate::system::file::TempDir::new(false).expect("a temporary directory");
        let path = dir.path().join("log.txt");
        let log = std::sync::Arc::new(crate::cli::ToolLog::new());
        log.set_location("PeakPickerHiRes:1:");
        log.set_destination(path.to_str());
        let ctx = ToolContext::new(
            PeakPickerHiRes::NAME,
            "1.0.0",
            "PeakPickerHiRes:1:",
            Param::new(),
            log.clone(),
        );
        let mut err = Vec::new();
        let code =
            check_input(&ctx, experiment, &mut err).expect("writing to a vector cannot fail");
        log.finish();
        let lines = std::fs::read_to_string(&path)
            .unwrap_or_default()
            .lines()
            .map(|line| line.get(20..).unwrap_or(line).to_owned())
            .collect();
        (
            code,
            String::from_utf8(err).expect("the source messages are text"),
            lines,
        )
    }

    /// The unsorted branches cannot be reached through the tool's loader, which
    /// sorts; they are exercised here on constructed experiments, with the
    /// source's messages and `INCOMPATIBLE_INPUT_DATA`. The source writes both
    /// with `writeLogError_` (`PeakPickerHiRes.cpp:243`, `253`), so each is
    /// also a line of the `-log` file, under the INI location.
    #[test]
    fn unsorted_records_are_refused_with_the_source_messages() {
        let mut experiment = MSExperiment::default();
        experiment.spectra.push(spectrum(&[100.0, 101.0]));
        experiment.spectra.push(spectrum(&[101.0, 100.0]));
        experiment.chromatograms.push(chromatogram(&[2.0, 1.0]));
        assert_eq!(
            checked(&experiment),
            (
                Some(ExitCode::IncompatibleInputData),
                format!("{UNSORTED_SPECTRA_ERROR}\n"),
                vec![format!("PeakPickerHiRes:1:: {UNSORTED_SPECTRA_ERROR}")]
            )
        );
        experiment.spectra.pop();
        assert_eq!(
            checked(&experiment),
            (
                Some(ExitCode::IncompatibleInputData),
                format!("{UNSORTED_CHROMATOGRAMS_ERROR}\n"),
                vec![format!(
                    "PeakPickerHiRes:1:: {UNSORTED_CHROMATOGRAMS_ERROR}"
                )]
            )
        );
        experiment.chromatograms[0] = chromatogram(&[1.0, 1.0, 2.0]);
        assert_eq!(checked(&experiment), (None, String::new(), Vec::new()));
    }

    /// The empty-input refusal needs no spectra and no chromatograms; one
    /// chromatogram is enough to proceed, as in the source. The source writes
    /// it with `OPENMS_LOG_WARN` (`PeakPickerHiRes.cpp:233-234`), which does
    /// not reach the `-log` file.
    #[test]
    fn only_an_input_without_spectra_and_chromatograms_is_empty() {
        let mut experiment = MSExperiment::default();
        assert_eq!(
            checked(&experiment),
            (
                Some(ExitCode::IncompatibleInputData),
                format!("{EMPTY_INPUT_WARNING}\n"),
                Vec::new()
            )
        );
        experiment.chromatograms.push(chromatogram(&[]));
        assert_eq!(checked(&experiment), (None, String::new(), Vec::new()));
    }
}

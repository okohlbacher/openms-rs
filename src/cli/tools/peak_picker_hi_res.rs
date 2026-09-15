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
//! Not ported yet: `-processOption lowmemory`, the source's
//! `PPHiResMzMLConsumer` path through `MzMLFile::transform`, is refused with
//! `INCOMPATIBLE_INPUT_DATA` until package P4-PICKER-LOWMEM ports it; the
//! debug dump of the algorithm parameters at `-debug 3`, and progress logging.
//!
//! `docs/TOPP_PEAK_PICKER_HI_RES_SUPPORT.md` lists the source members, the
//! preserved conventions, the native differences and the evidence.

use crate::cli::{ExitCode, Tool, ToolContext, ToolSpec};
use crate::format::PeakFileOptions;
use crate::format::file_handler::FileHandler;
use crate::format::file_types::FileType;
use crate::format::mzml;
use crate::kernel::MSExperiment;
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
/// Registration follows `registerOptionsAndFlags_` (`PeakPickerHiRes.cpp:152-163`):
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

/// Source warning for per-peak ion mobility (`PeakPickerHiRes.cpp:222-226`),
/// around the name of the spectrum's ion mobility peak type.
fn ion_mobility_warning(peak_type: IonMobilityPeakType) -> String {
    format!(
        "Warning: Input contains per-peak ion mobility data (IM_PEAK, {}). PeakPickerHiRes picks in m/z only and reports intensity-weighted mean ion mobility. This produces incorrect results on unbinned data. Consider IonMobilityBinning or PeakPickerIM first.",
        im_peak_type_to_string(peak_type)
    )
}

/// Source warning for an input without spectra and chromatograms
/// (`PeakPickerHiRes.cpp:233-234`). The log stream ends the line.
const EMPTY_INPUT_WARNING: &str = "The given file does not contain any conventional peak data, but might contain chromatograms. This tool currently cannot handle them, sorry.";

/// Source error for an unsorted spectrum (`PeakPickerHiRes.cpp:243`).
const UNSORTED_SPECTRA_ERROR: &str = "Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!";

/// Source error for an unsorted chromatogram (`PeakPickerHiRes.cpp:253`),
/// which also says m/z.
const UNSORTED_CHROMATOGRAMS_ERROR: &str = "Error: Not all chromatograms are sorted according to peak m/z positions. Use FileFilter to sort the input!";

/// The refusal of `-processOption lowmemory` until package P4 ports it.
const LOW_MEMORY_UNSUPPORTED: &str = "PeakPickerHiRes -processOption lowmemory is not ported yet (package P4-PICKER-LOWMEM of the early TOPP bundle ports it); use -processOption inmemory";

/// The in-memory input checks of `main_` before picking
/// (`PeakPickerHiRes.cpp:216-256`), in source order.
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
/// # Errors
///
/// Returns [`Error::Io`] when a line cannot be written to `err`.
fn check_input(experiment: &MSExperiment, err: &mut dyn Write) -> Result<Option<ExitCode>> {
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
        writeln!(err, "{UNSORTED_SPECTRA_ERROR}")?;
        return Ok(Some(ExitCode::IncompatibleInputData));
    }
    if !experiment
        .chromatograms
        .iter()
        .all(|chromatogram| chromatogram.is_sorted())
    {
        writeln!(err, "{UNSORTED_CHROMATOGRAMS_ERROR}")?;
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

/// The per-level summary `pickExperiment` logs (`PeakPickerHiRes.cpp:559-563`):
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
    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        Self::run_io(ctx, &mut std::io::stdout(), &mut std::io::stderr())
    }

    /// Source `main_` (`PeakPickerHiRes.cpp:188-273`).
    ///
    /// 1. The `algorithm` subsection configures the picker, as
    ///    `setParameters`, with every source behaviour of
    ///    [`PickingCompatibility::source`](crate::processing::peak_picking::PickingCompatibility::source)
    ///    enabled, because the source picks duplicate positions, negative
    ///    intensities and non-positive spline maxima without refusing them.
    /// 2. `-processOption lowmemory` is refused with [`Error::Unsupported`]
    ///    (`INCOMPATIBLE_INPUT_DATA`) before anything is read.
    /// 3. The input is loaded as mzML with [`PeakPickerHiRes::read_options`] —
    ///    the source's dangling header references accepted (decision D10) and
    ///    the size-derived library ceilings, which admit an instrument-sized
    ///    run — then checked as described at `check_input`.
    /// 4. Spectra are picked with the spectrum type checked unless `-force` is
    ///    given (`check_spectrum_type = !force`). A centroided spectrum in manual
    ///    mode without `-force` ends the run with
    ///    `Error: Unexpected internal error (Error: Centroided data provided but profile spectra expected.)`
    ///    and [`ExitCode::UnknownError`], as the source's `IllegalArgument`
    ///    takes the `BaseException` arm of `TOPPBase::main`. Nothing is written.
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
    /// refusal is written as `Error: Unexpected internal error (<reason>)` and
    /// returns [`ExitCode::UnknownError`]: those are native bounds (points and
    /// work per record, metadata copies) and the FWHM search that never
    /// terminates in the source, not parameter errors. A refusal by the
    /// operating system to start the `-threads` workers is reported the same
    /// way, as `Error: Unexpected internal error (cannot start <n> worker
    /// threads: <reason>)` with [`ExitCode::UnknownError`]: the pool is built
    /// inside the picking call now (`pick_experiment` in this module), so the
    /// [`Error::Io`] it raises reaches the same arm as a picker failure instead
    /// of propagating out of `run_io` as it did while the pool wrapped the whole
    /// body. No output file is written in either case. The exception is
    /// [`Error::Unsupported`], which propagates (`INCOMPATIBLE_INPUT_DATA`):
    /// the picker returns it for `SignalToNoise:auto_mode` 1 as soon as noise
    /// estimation runs, where the source reads out of bounds and crashes
    /// (oracle `PPHR_auto_mode_1`, SIGBUS or SIGSEGV). With
    /// `signal_to_noise` 0 the estimator never runs in either implementation
    /// and the run succeeds.
    fn run_io(ctx: &ToolContext, out: &mut dyn Write, err: &mut dyn Write) -> Result<ExitCode> {
        let input = ctx.string("in")?;
        let output = ctx.string("out")?;
        let process_option = ctx.string("processOption")?;

        let mut picker = Picker::from_param(&ctx.subsection("algorithm")?)?;
        picker.compatibility = PickingCompatibility::source();
        picker.check_spectrum_type = !ctx.force();

        if process_option == "lowmemory" {
            return Err(Error::Unsupported(LOW_MEMORY_UNSUPPORTED.into()));
        }

        let mut experiment = FileHandler::load_experiment_with_read_options(
            input,
            &[FileType::MzMl],
            &PeakFileOptions::default(),
            &Self::read_options(),
        )?;
        if let Some(code) = check_input(&experiment, err)? {
            return Ok(code);
        }

        let report = match pick_experiment(ctx, &picker, &mut experiment) {
            Ok(report) => report,
            Err(error @ Error::Unsupported(_)) => return Err(error),
            // The centroided refusal is reported with the source's bare
            // message, without this port's `invalid value: ` prefix.
            Err(Error::InvalidValue(reason)) if reason == CENTROIDED_INPUT_MESSAGE => {
                writeln!(err, "Error: Unexpected internal error ({reason})")?;
                return Ok(ExitCode::UnknownError);
            }
            Err(error) => {
                writeln!(err, "Error: Unexpected internal error ({error})")?;
                return Ok(ExitCode::UnknownError);
            }
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

    /// The exit code and the error-stream text of [`check_input`], captured
    /// from the stream it writes to.
    fn checked(experiment: &MSExperiment) -> (Option<ExitCode>, String) {
        let mut err = Vec::new();
        let code = check_input(experiment, &mut err).expect("writing to a vector cannot fail");
        (
            code,
            String::from_utf8(err).expect("the source messages are text"),
        )
    }

    /// The unsorted branches cannot be reached through the tool's loader, which
    /// sorts; they are exercised here on constructed experiments, with the
    /// source's messages and `INCOMPATIBLE_INPUT_DATA`.
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
                format!("{UNSORTED_SPECTRA_ERROR}\n")
            )
        );
        experiment.spectra.pop();
        assert_eq!(
            checked(&experiment),
            (
                Some(ExitCode::IncompatibleInputData),
                format!("{UNSORTED_CHROMATOGRAMS_ERROR}\n")
            )
        );
        experiment.chromatograms[0] = chromatogram(&[1.0, 1.0, 2.0]);
        assert_eq!(checked(&experiment), (None, String::new()));
    }

    /// The empty-input refusal needs no spectra and no chromatograms; one
    /// chromatogram is enough to proceed, as in the source.
    #[test]
    fn only_an_input_without_spectra_and_chromatograms_is_empty() {
        let mut experiment = MSExperiment::default();
        assert_eq!(
            checked(&experiment),
            (
                Some(ExitCode::IncompatibleInputData),
                format!("{EMPTY_INPUT_WARNING}\n")
            )
        );
        experiment.chromatograms.push(chromatogram(&[]));
        assert_eq!(checked(&experiment), (None, String::new()));
    }
}

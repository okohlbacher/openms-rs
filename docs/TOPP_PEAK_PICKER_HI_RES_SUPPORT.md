# The `PeakPickerHiRes` TOPP tool

`cli::tools::PeakPickerHiRes` ports `OpenMS4-topp/src/PeakPickerHiRes.cpp`
(topp `174b576`, sha256 `e253bb08…`), the command-line wrapper of the
high-resolution peak picker. The picker itself, its parameter contract and its
signal-to-noise estimator are `processing::peak_picking`
([PEAK_PICKING_SUPPORT](PEAK_PICKING_SUPPORT.md)); the lifecycle around the
tool is [TOPP_CLI_SUPPORT](TOPP_CLI_SUPPORT.md). No C++ library is called or
built by the crate.

| Source | Rust |
| --- | --- |
| `topp/src/PeakPickerHiRes.cpp`: `TOPPPeakPickerHiRes` | `src/cli/tools/peak_picker_hi_res.rs` |
| its `main(argc, argv)` | `src/bin/PeakPickerHiRes.rs` |
| `PROCESSING/CENTROIDING/PeakPickerHiRes.h`, `.cpp` | `src/processing/peak_picking.rs` (read-only here) |
| `IONMOBILITY/IMTypes.h` `determineIMFormat`, `imPeakTypeToString` | `src/metadata/im_types.rs` (read-only here) |

```console
$ PeakPickerHiRes -in profile.mzML -out centroided.mzML
$ PeakPickerHiRes -ini picker.ini -in profile.mzML -out centroided.mzML -threads 16
$ PeakPickerHiRes -write_ini picker.ini
```

## API mapping

Every member of `TOPPPeakPickerHiRes`, in source order.

| Source member | Rust | Notes |
| --- | --- | --- |
| `TOPPPeakPickerHiRes()` (name, description) | `Tool::NAME`, `Tool::DESCRIPTION` | Same two strings. |
| `class PPHiResMzMLConsumer` (`processSpectrum_`, `processChromatogram_`, members `pp_`, `ms_levels_`) | `LowMemoryPicker` + `MSDataWritingConsumer` | The source class is the consumer base plus the two hooks; the port is the hooks alone, handed to `format::ms_data_writing_consumer::MSDataWritingConsumer`, which is that base. `ms_levels_` is read from the picker rather than copied out of its parameters a second time. See *The low-memory mode*. |
| `registerOptionsAndFlags_` | `Tool::register` | `-in`, `-out` (mzML only), advanced `-processOption` restricted to `inmemory,lowmemory`, subsection `algorithm`. |
| `getSubsectionDefaults_(section)` | `Tool::subsection_defaults` | Returns `PeakPickerHiRes::defaults()` for every section, as the source ignores its argument. |
| `doLowMemAlgorithm(pp)` | `run_low_memory` | Consumer on `-out`, `addDataProcessing`, `MzMLFile::transform` → `mzml::transform_with_options` with the tool's read options. |
| `main_(int, const char**)` | `Tool::run_io` | `Tool::run` forwards with the process streams. |
| members `in`, `out` | `ToolContext::string("in")`, `("out")` | Read where the source reads them. |
| `getParam_().copy("algorithm:", true)` + `pp.setParameters` | `ToolContext::subsection("algorithm")` + `PeakPickerHiRes::from_param` | |
| `writeDebug_("Parameters passed to PeakPickerHiRes", pepi_param, 3)` | not ported | No debug writer exists on `ToolContext`; no registered test compares debug output. |
| `pp.setLogType(log_type_)`, `mz_data_file.setLogType(log_type_)`, `FileHandler::loadExperiment(..., log_type_)` | not ported | The picker and the readers have no progress logging; `ToolContext::progress_log_type` exists for when they do. |
| `FileHandler().loadExperiment(in, exp, {MZML}, log_type_)` | `FileHandler::load_experiment_with_read_options(in, &[FileType::MzMl], &PeakFileOptions::default(), &PeakPickerHiRes::read_options())` | `read_options` is `mzml::ReadOptions::source()`: decision D10's source-compatibility switches, which `TOPP_PeakPickerHiRes_5` needs for its dangling header references, over the library's size-derived ceilings, which an instrument-sized run needs. Library defaults stay strict. |
| the `IMTypes::determineIMFormat` warning loop | `check_input` | Warns once, as the source `break` does. |
| `ms_exp_raw.empty() && getChromatograms().empty()` → `INCOMPATIBLE_INPUT_DATA` | `check_input` | Same message and code. |
| the two `isSorted()` loops → `INCOMPATIBLE_INPUT_DATA` | `check_input` | Same messages and code; unreachable through the loader in both implementations (see *Preserved source conventions*). |
| `pp.pickExperiment(ms_exp_raw, ms_exp_peaks, !getFlag_("force"))` | `PeakPickerHiRes::pick_experiment` with `check_spectrum_type = !ctx.force()` | The per-level `OPENMS_LOG_INFO` summary is written by the tool, because the native picker does not log. |
| `addDataProcessing_(ms_exp_peaks, getProcessingInfo_(DataProcessing::PEAK_PICKING))` | `ToolContext::processing_info(&[ProcessingAction::PeakPicking])` + `ToolContext::add_data_processing` | One shared record on every spectrum and chromatogram. |
| `FileHandler().storeExperiment(out, ms_exp_peaks, {MZML})` | `FileHandler::store_experiment(out, exp, Some(FileType::MzMl))` | |
| `return EXECUTION_OK` | `Ok(ExitCode::ExecutionOk)` | |
| — | `render_list_parameters` | Native: renders list-valued parameters of the processing record as the source mzML writer does (see *Native differences*). |

### Capabilities

| Mode | State |
| --- | --- |
| `-processOption inmemory` (the default) | ported, `TOPP_PeakPickerHiRes_1`, `_2`, `_5`, `_6` |
| `-processOption lowmemory` | ported, `TOPP_PeakPickerHiRes_3`, `_4` |
| `-force` | ported (`check_spectrum_type = !force`); **inert under `-processOption lowmemory`**, as in the source |
| `-test` | ported (processing record and unique-id seed through the framework) |
| `-write_ini`, `-ini`, `-threads`, `-debug`, `-no_progress` | through the framework; `-debug` prints nothing here, `-no_progress` has nothing to suppress |
| `-write_ctd` and the CWL/JSON writers | refused by the framework |
| `-log`, `-instance` | inert, as in the rest of the port |
| ion mobility (`IM_PEAK`) | picked with the source's warning; the mean ion mobility array is reproduced |
| input formats | mzML only, as the source registers |

## Preserved source conventions

- **Order of `main_`.** Parameters, then the low-memory branch, then loading,
  the ion mobility warning, the empty-input refusal, the two sortedness
  refusals, picking, the processing record and the output. A refusal writes
  nothing.
- **The `algorithm` subsection is the picker's `getDefaults()`**, so an INI
  written by the C++ tool loads unchanged, and `-algorithm:<name>` values reach
  the picker. `PeakPickerHiRes -write_ini` and the C++ `-write_ini` agree line
  by line under the upstream `FuzzyDiff` settings (P3 oracle `write_ini`), and
  `TOPPWRITEINI_OVERWRITE` matches the retained `WRITE_INI_OUT.ini` with the
  registered `version` whitelist.
- **Source acceptance of degenerate input.** The picker runs with
  `PickingCompatibility::source()`: duplicate positions, negative intensities
  and a non-positive spline maximum are picked rather than refused, as the C++
  picker does. The library default stays strict.
- **The centroided refusal exits 8.** `pickExperiment` throws
  `Exception::IllegalArgument` for a selected centroided spectrum when the type
  is checked; `TOPPBase::main` catches it in the `BaseException` arm
  (`TOPPBase.cpp:495-499`), so the process prints
  `Error: Unexpected internal error (Error: Centroided data provided but profile spectra expected.)`
  and exits `UNKNOWN_ERROR`, not `ILLEGAL_PARAMETERS`. Reproduced, including
  the wording (oracle `PPHR_6_noforce`).
- **Sortedness is checked after a loader that sorts.** `FileHandler` loads mzML
  with default `PeakFileOptions`, whose `sort_spectra_by_mz` and
  `sort_chromatograms_by_rt` are set, so `MzMLHandler` (`:218-221`, `:299-302`)
  and this port's reader both sort each record before the tool sees it. The two
  refusals are therefore unreachable on this path in **both** implementations;
  they are ported because the source has them, exercised by the module's unit
  test, and the P3 oracle's `unsorted_spectrum` and `unsorted_chromatogram`
  cases show the C++ tool picking such an input and producing the sorted
  workflow's output, which this port reproduces.
- **The empty-input diagnostic is a warning, not an error line**, and the tool
  still exits 11.
- **The per-level summary** (`#Spectra that needed to and could be picked by
  MS-level:` and one `  MS-level <n>: <picked> / <total>` line per level in
  ascending order) is written to standard output, where the source's
  `OPENMS_LOG_INFO` writes it. The header appears even for an input without
  spectra.
- **`-threads` changes no output byte.** `PeakPickerHiRes.cpp` has no OpenMP, so
  the source picks serially; this port picks in parallel and still writes the
  same file at every worker count (tested at 1, 2, 8, 16, 32 and 0, and on the
  2.3 GB benchmark run at 1, 8 and 32). See native difference 10.

## The low-memory mode

`-processOption lowmemory` is not the in-memory mode with a smaller footprint.
It is a second, shorter program: `doLowMemAlgorithm` builds a
`PPHiResMzMLConsumer` on the output file, gives it the `peak picking` processing
record, and calls `MzMLFile::transform`, which streams the input past it. None
of `main_`'s input handling runs. The source's low-memory path is the
specification for this port's low-memory path, including where the two modes
disagree.

### What it does, in order

1. **The output file is created before the input is read**, because the source
   consumer opens its `std::ofstream` in its constructor
   (`MSDataWritingConsumer.cpp:33`). An input that cannot be parsed still leaves
   a file behind — empty, when the failure comes before the first record. The
   in-memory mode writes nothing until the whole run has succeeded.
2. **The input is read twice.** `transform` runs `transformFirstPass_`, which
   parses the whole file for the declared record counts and the experimental
   settings and hands them to `setExpectedSize`/`setExperimentalSettings`, then
   parses it again for the records (`MzMLFile.cpp:178-191` and `212-231`). The mode trades I/O
   for memory: a low-memory run reads roughly twice the bytes an in-memory run
   does.
3. **The header comes from those settings plus the first record**, and each list
   tag announces the count the first pass declared, not the records that follow.
4. **Each record is picked and written immediately**, then dropped. Peak memory
   is one read batch (`max_data_pool_size`, 100 records, as upstream) plus the
   rendered text of one record, not the experiment.
5. **The document is closed with an index.** See *Native differences* 5 and the
   writing consumer's own support document.

### Where it differs from the in-memory mode, by the source's own construction

| | `-processOption inmemory` | `-processOption lowmemory` |
| --- | --- | --- |
| automatic-mode type test | `getType(true)`: stored type, then a `PEAK_PICKING` entry in the record's processing history, then `PeakTypeEstimator` over the samples (`PeakPickerHiRes.cpp:510`, `531`) | `s.getType()`, the `SpectrumSettings` accessor `MSSpectrum` re-exposes with `using` (`MSSpectrum.h:655`): the **stored type only** (`PeakPickerHiRes.cpp:124`) |
| centroided data on a selected MS level | `IllegalArgument` unless `-force` | picked; **there is no check at all**, so `-force` is inert |
| per-peak ion mobility | warns once | silent |
| input with neither spectra nor chromatograms | `INCOMPATIBLE_INPUT_DATA` | exit 0; nothing is written, because no record ever reaches the consumer |
| unsorted records | two `INCOMPATIBLE_INPUT_DATA` checks (unreachable through either loader, which sorts) | no checks; the reader sorts, so the two modes agree here |
| per-MS-level summary on stdout | written by `pickExperiment` | none |
| chromatograms | all picked | all picked |
| input reads | one | two |
| a record needing header entries the first record did not contribute | numbered against a header written from the whole experiment | written with a dangling reference, numbered by the record's position in the stream (`MzMLHandler.cpp:5252-5272`) |
| what a failing run leaves on disc | nothing | the batches already sent — floor(N / 100) × 100 records, closed and indexed, under the count pass one declared |
| an `-out` that names an existing directory | `CANNOT_WRITE_OUTPUT_FILE` | in the source, **exit 0 having written nothing**: the consumer's constructor never checks its `std::ofstream` (`MSDataWritingConsumer.cpp:33`) |

The first row is the one that changes numbers. A spectrum whose mzML carries
`MS:1000525` but neither `MS:1000127` nor `MS:1000128` has an unknown stored
type, which is common — it is what a converter that only records the
representation term produces. The in-memory mode falls through to the estimator
and may call it centroided and copy it; the low-memory mode sees `UNKNOWN`,
which is not `CENTROID`, and picks it. Upstream's own workflow 6 input is such a
file: the in-memory mode reports `MS-level 1: 0 / 1` and writes the 33 input
samples, and the low-memory mode writes 4 centroids. Both are reproduced
(`low_memory_automatic_mode_tests_only_the_stored_spectrum_type`,
`low_memory_never_refuses_centroided_data_and_force_is_inert`).

### Threads

The mode is serial in the source and here. The source's consumer dispatch loop
hands over one record at a time (`MzMLHandler.cpp:259-274`), and the one OpenMP
region on that path decodes binary arrays, which this port's reader does not
parallelise either. `-threads` therefore reaches nothing on this path, and the
written bytes cannot depend on it — trivially, rather than by the batch-order
argument the in-memory mode needs. Pinned at 1, 8 and 32 on an input with
spectra and on one with chromatograms
(`the_low_memory_output_is_bit_identical_at_every_thread_count`).

### What the mode is worth, measured

`ibminode06` (128 cores, 995 GB, load 0.04 per core at the start), the C++
Release install at the pins
(`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`) and the
Rust release binary, over the benchmark's 2.3 GB Q Exactive profile run
`UK222.mzML` (2,317,975,830 bytes, 40,856 spectra, one chromatogram), with the
INI the C++ tool writes with `-write_ini`, `-threads 1` and `-no_progress`.
**Peak RSS and output size are the measurement**; the wall times are context and
not a timing claim — the timing node is not used by this lane. Drivers, logs and
hashes: `../oracle/p4-lowmemory/`.

| | exit | peak RSS | wall | output |
| --- | --- | --- | --- | --- |
| C++ `-processOption lowmemory` | 0 | 165,544 KiB (161.7 MiB) | 26.6 s | 549,528,621 B |
| C++ in-memory | 0 | 3,977,240 KiB (3.79 GiB) | 25.2 s | 549,528,688 B |
| Rust `-processOption lowmemory` | 0 | **83,632 KiB (81.7 MiB)** | 66.8 s | 535,615,900 B |
| Rust in-memory | 0 | 3,451,904 KiB (3.29 GiB) | 25.0 s | 535,615,963 B |

The mode does what it exists for, in both implementations and more so here:
**41 times less resident memory than the in-memory mode on a 2.3 GB input**, and
half the C++ low-memory footprint. What it costs is time: this port's streaming
path renders each record as a one-record document and splices the record element
out of it, which is how one definition of the mzML encoding is kept, and that
shows up as 2.7 times the in-memory wall time where the C++ path, which writes
the record element directly, is level with its own in-memory run. Nothing about
the output depends on it.

**The two modes write the same data.** In each implementation the mzML body,
`<run …>` through `</mzML>`, is byte-identical between the modes: 532,278,411
bytes (`e2141e3a…`) for this port, 546,109,040 bytes (`ae37b166…`) for C++. All
four outputs carry 40,856 spectra, one chromatogram and 22,784,372 summed
`defaultArrayLength`.

The headers differ by one line, 63 bytes here and 67 in C++, and **the source
does it too**: the consumer's header comes from the experimental settings plus
the first record (`MSDataWritingConsumer.cpp:76-83`), so the `fileContent` terms
that are derived from records come from that one record. On this input the
in-memory mode writes `MS1 spectrum` and `MSn spectrum` and the low-memory mode
writes `MS1 spectrum` alone — in the C++ Release build exactly as here. A
low-memory output therefore understates the file's own content by design, and
that is reproduced, not repaired.

### How a failure ends, and what it leaves behind

`doLowMemAlgorithm` catches nothing (`PeakPickerHiRes.cpp:170-186`), so every
failure of the mode is classified by `TOPPBase::main` on the exception type
alone. A reader failure is therefore `Error: Unable to read file (…)` and
`INPUT_FILE_CORRUPT`, exactly as in the in-memory mode, whose `loadExperiment`
failures take the same arm (`TOPPBase.cpp:460-465`). This port reproduces that:
`run_low_memory` lets a reader error out to the framework's own mapping of that
catch chain, and keeps `Error: Unexpected internal error (…)` / `UNKNOWN_ERROR`
for the three error kinds whose source counterparts — `InvalidValue`,
`InvalidRange` and `MissingInformation` — derive straight from `BaseException`
and take its arm in the source while this port's framework maps them to
parameter codes. Executed on `ibminode06` against the C++ Release build over a
truncated mzML, a non-XML input and an mzML whose third record carries base64
of an impossible length: both implementations exit 3 in both modes, with
`Unable to read file`, and both leave a zero-byte output behind in the
low-memory mode and none in the in-memory mode
(`../oracle/p4-lowmemory/logs/fixdiff_06.log`, cases `trunc`, `garbage` and
`badb64`; pinned by `a_corrupt_input_exits_3_in_both_modes`).

The document is closed on the failing path as well. The source's
`~MSDataWritingConsumer` runs `doCleanup_` on every path
(`MSDataWritingConsumer.cpp:37-40`), which closes the open list and writes the
footer whenever writing started (`:151-173`), so a source run that throws after
its first record still leaves a closed, indexed document. `finish` is this
port's destructor — Rust cannot report an error from a drop — and
`run_low_memory` calls it on the failing path too, discarding its own result so
that the failure which ended the run is the one reported. The records already
written stay where they are, under the `count` the first pass declared: a
streaming writer cannot take bytes back. A failure before the first record
leaves the created file empty, as `doCleanup_` writes nothing while
`started_writing_` is false.

**That is measured, not only argued.** A C++ low-memory run can indeed end in a
closed, partial document, and it takes more than a hundred records to see it.
The source's second pass decodes and hands over in batches of
`maximal_data_pool_size_`, 100 by default (`PeakFileOptions.h:248`):
`populateSpectraWithData_` decodes a whole batch under OpenMP, throws
`ParseError` if any record of it failed (`MzMLHandler.cpp:198-244`), and only
then loops over the batch calling `consumeSpectrum` (`:259-274`). A corrupt
record therefore kills its entire batch before any of it reaches the writer, so
on a five-record file the C++ run writes nothing whichever record is corrupt —
first, third or last, all measured (`logs/fixdiff2_06.log`). On a 110-spectrum
file built on the node by the Release `FileMerger`, with the base64 of spectrum
104 corrupted, the C++ low-memory run exits 3 and leaves **875,255 bytes: the
100 records of the first batch, `</mzML>`, an index and a `fileChecksum`, under
the announced `count="110"`** (`logs/fixdiff3_06.log`, case `many_104`). That is
the document this port now leaves on its own failing paths.

Pinned by `a_failure_after_the_first_record_still_closes_the_document`, which
reaches the case through `SignalToNoise:auto_mode 1` with `ms_levels 2`, so the
MS1 spectrum is copied and written before the first MS2 spectrum reaches the
estimator. The source cannot be compared on *that* command line, because it does
not fail there in any orderly way: the same run on `ibminode06` writes all five
records and then dies of SIGSEGV (exit 139), leaving 404,915 bytes with no
index and no footer, since a signal runs no destructor
(`logs/fixdiff_06.log`, case `am1`; native difference 2). This port exits 11
with its diagnosis over a closed, indexed document holding the one record it
had written.

**Neither counting pass reads record contents, and this port leaves exactly
the same partial document.** The source's runs with `LD_RAWCOUNTS` and sets
`skip_spectrum_` (`MzMLHandler.cpp:966-974`); this port's sets `state.raw` and,
at the list's start tag, a `skip_depth` that skips to the matching end tag
(`src/format/mzml_counts.rs:859`, `:465-469`). A record that is well-formed XML
but wrong inside is therefore invisible to both first passes and is discovered
by both second passes, which deliver records to the writer in batches of
`max_data_pool_size` — the same 100 on both sides
(`src/format/mzml_consumer.rs:193`, `PeakFileOptions.h:248`). So **a low-memory
run that fails at record N leaves floor(N / 100) × 100 records behind**,
whichever implementation runs it, in a closed, indexed and reloadable document
announcing the count pass one declared. Measured on the 110-record `batches`
fixture with the base64 of one record corrupted: at index 4 both sides write
nothing, and at index 104 the C++ leaves 150,901 bytes holding 100 records and
this port leaves 137,074 bytes holding the same 100, both closed, indexed and
announcing `count="110"`. This port's partial reloads through the tool at exit
0 with 100 spectra (`logs/closediff2_06.log` section D,
`logs/closediff3_06.log` section C).

**Where the two readers disagree, that partial document is this mode's most
consequential divergence, and it is silent.** Three inputs measured on
`ibminode06` make it concrete: a record whose `defaultArrayLength` declares one
value too many, a non-numeric `scan start time`, and two records sharing a
native id. The C++ build exits 0 on all three in **both** process options —
warning once about the array length and saying nothing at all about the other
two — while this port refuses all three with `Unable to read file` in **both**
(`logs/fixdiff_06.log`, cases `badlen`, `badrt`, `dupid`;
`logs/closediff2_06.log` section D). The refusal is the mzML reader's
strictness, not the mode's: [MZML_SUPPORT](MZML_SUPPORT.md) already states that
duplicate records and incorrect array lengths are errors, and a CV value that
cannot be converted is one too. The **exit code** is therefore the same in both
modes. The **artefact** is not, and that part does belong to this mode: the
in-memory mode writes no file at all, while the low-memory mode, on an input
longer than one batch, leaves a complete, reloadable mzML holding only the
completed batches under a `spectrumList count` announcing every record — a
silently truncated output on a file the C++ tool accepts and writes in full.
Measured for each of the three kinds at index 4 and at index 104 of the
`batches` fixture: the C++ exits 0 with all 110 records in both modes at both
positions; this port leaves 0 records at index 4 and 100 at index 104.

Pinned by `a_low_memory_failure_leaves_the_batches_already_written`, which
runs the boundary from both sides on that fixture.

### The one place this port is stricter

The consumer's `CountPolicy::Checked` is kept: a document whose declared list
counts and actual records disagree ends the run with
`Error: Unexpected internal error (invalid value: mzML list counts announce …)`
and `UNKNOWN_ERROR`, after the output has been closed. The source writes such a
document silently. Executed on `ibminode06` on `PeakPickerHiRes_input.mzML`
with `<spectrumList count="5"` rewritten to `count="9"`: the C++ low-memory run
exits 0 and writes `<spectrumList count="9">` over its five records, and both
in-memory runs are unaffected (`../oracle/p4-lowmemory/logs/fixdiff_06.log`,
case `badcount`). The file this port leaves behind is the same complete,
indexed document carrying the same lying count — the refusal is a diagnosis
added to the source's output, not a change to it. Pinned by
`a_low_memory_run_reports_a_lying_list_count_over_a_closed_document`.

## Native differences

1. **List-valued parameters in the processing record.** Outside `-test`,
   `getProcessingInfo_` records every resolved parameter; the source mzML writer
   renders a list-valued one as text
   (`parameter: algorithm:ms_levels` = `[]`, P3 oracle `notest`). This port's
   mzML writer refuses list and empty metadata
   (`Empty/list metadata has no lossless mzML scalar encoding`), so the tool
   renders those values itself before attaching the record
   (`render_list_parameters`), which reproduces the C++ text. Without it no
   run outside `-test` could store its output, because `ms_levels` is an empty
   list by default. The rendering belongs in `cli::processing::processing_info`
   or in the writer, which are outside this package: integrator request.
2. **`SignalToNoise:auto_mode 1` is refused instead of crashing.** The source
   reads out of bounds in `AUTOMAXBYPERCENT` and dies from a memory-access
   signal (oracle `PPHR_auto_mode_1`: SIGBUS 138, SIGSEGV 139 on repetition).
   The picker returns `Error::Unsupported` as soon as noise estimation runs, so
   the tool exits `INCOMPATIBLE_INPUT_DATA` with an explicit message and writes
   nothing. With `signal_to_noise` 0 the estimator never runs in either
   implementation and both produce the ordinary output (P3 oracle
   `auto_mode_1_without_estimation`).
3. **The one difference the source has between its two modes' outputs does not
   appear here.** Upstream retains two output files for workflow 1 because of
   it, with the comment that the low-memory output "SHOULD be identical to
   'PeakPickerHiRes_output.mzML', but due to a missing 'Dataprocessing' entry
   (which is not known when writing the mzML header), we need an extra output
   file" (`CMakeLists.txt:2534`). The two retained files differ in exactly one
   byte: `<dataProcessingList count="3">` against `count="2">`. The C++ writer
   puts `max(1, dps.size() + <float data arrays of the whole experiment>)` in
   that attribute (`MzMLHandler.cpp:5161-5170`), and the consumer's
   `writeHeader_` is handed a dummy map holding only the first record, so it
   counts one record's arrays. This port's writer counts the processing
   histories it writes (`mzml_header/write.rs:329`), which does not depend on
   the records at all, so both modes write the same count and the low-memory
   output here is **byte-identical to the in-memory output** on both upstream
   registrations. A pre-existing native difference of the mzML writer, pinned
   from the retained pair by
   `the_two_retained_cpp_outputs_differ_only_in_the_data_processing_count`.
4. **The C++ low-memory mode cannot read an mzML holding both spectra and
   chromatograms while progress logging is on; this port can.** The C++ tool
   exits 3 with
   `Error: Unable to read file (- due to that error of type Precondition failed
   in: StopWatch.cpp@43-void OpenMS::StopWatch::start())` and writes nothing.
   `MzMLFile::transform` parses the file twice through one `ProgressLogger`
   (`MzMLFile.cpp:178-190`); in the first pass `MzMLHandler` runs with
   `LD_RAWCOUNTS` and, at `<chromatogramList>`, calls
   `logger_.startProgress("loading chromatogram list")`
   (`MzMLHandler.cpp:997`) and then immediately throws `EndParsingSoftly`
   because it now has both counts (`MzMLHandler.cpp:1001-1006`), so the
   `endProgress()` at `</chromatogramList>` (`MzMLHandler.cpp:1493-1498`, the `endProgress()` at `1497`) is
   never reached and the shared stopwatch is still running when the second pass
   calls `startProgress` again. One record kind alone means no early throw and a
   balanced pair, which is why the upstream registrations never see it: their
   inputs hold spectra only and chromatograms only. `-no_progress` avoids it;
   `-test` does not. Reproduced on a 450 KB two-kind input made on the node by
   the Release `FileMerger`, and on the 2.3 GB benchmark input; the control runs
   (spectra only, chromatograms only, and the in-memory mode on the same input)
   all exit 0. This port has no shared progress logger and completes on every
   one of those inputs. It reaches every caller of `MzMLFile::transform`, so
   `NoiseFilterGaussian`, `NoiseFilterSGolay`, `FileConverter -process_lowmemory`
   and `PeakPickerIM` are affected the same way: **integrator request for
   `OpenMS_CPP_ISSUES.md`**. Evidence:
   `../oracle/p4-lowmemory/logs/probe_06.log` and that oracle's `manifest.json`.
5. **Container differences are documented, not compared** (decision D6). The
   C++ output is an `indexedmzML` with an ISO-8859-1 declaration, the software
   alias `MS:1002135 TOPP PeakPickerHiRes`, a `dataProcessingList count`
   computed as `max(1, histories + float arrays)` (CPP-019) and the constant
   `0` for `fileChecksum` (CPP-049); this port writes an `indexedmzML` in UTF-8
   with the exact software name, one `dataProcessing` entry and the real SHA-1.
   Both modes write the same container here, as both do in the source. Both
   decode to the same content, which is what the tests compare.
6. **The ion mobility peak type in the warning is `im_profile`.** The source
   prints `imPeakTypeToString(spec.getIMPeakType())`; the native spectrum has no
   stored peak type. `im_profile` is what the source mzML reader stores for ion
   mobility data without `MS:1003441` (`MzMLHandler.cpp:253-255`), which is the
   oracle's text; an input carrying that term would read `im_centroided` in the
   source.
7. **No debug dump and no progress logging**, as listed in the API mapping.
8. **Bounded work on the input.** The C++ tool has no resource ceilings; this
   port bounds every cumulative quantity, and an input beyond a ceiling is
   refused before anything is written.

   On the **load path** the ceilings are the library defaults, which are
   size-derived (`mzml::InputScaling`, `src/format/mzml_scaling.rs`): each
   grows with the XML bytes already consumed, so work and storage stay linear
   in the input while a document of any realistic size fits. `PeakPickerHiRes::read_options` is
   `mzml::ReadOptions::source()` — those ceilings plus the source-compatibility
   switches a tool path takes (decision D10), named once so that the tests use
   the same definition and a switch added to that constructor reaches the tool.
   The ceilings the tool shipped with were fixed, 10,000,000 raw points and
   512 MiB of XML, and no instrument-sized profile run fits under either: the
   2.3 GB `UK222.mzML` has 197,765,338 raw points. They became size-derived in
   the library (`fix/mzml-reader-scale`, `MZML_READER_SCALE_SUPPORT.md`), which
   this branch merges, and the file now loads through the tool (measured
   below).

   The **picker's** ceilings are still fixed and are not this tool's to set:
   `max_points` 1,000,000 and `max_work` 10,000,000 per record, and the pooled
   acquisition-copy ledger of `pick_experiment` (`src/processing.rs`:
   50,000,000 work and 256 MiB of metadata, charged over the whole experiment
   before any record is picked). That ledger, not the loader, is what refuses an
   instrument-sized run on `integrate/wave2`; the measurements below record both
   the refusal and the complete run with the sibling lane that lifts it.
9. **A picker failure that is not the centroided refusal** is reported as
   `Error: Unexpected internal error (<reason>)` with `UNKNOWN_ERROR`, the code
   `TOPPBase` gives an unmapped exception, rather than the framework's default
   mapping of `Error::InvalidValue` to `ILLEGAL_PARAMETERS`: these are data and
   resource conditions, not parameter errors. One of them, a non-converging
   FWHM bisection, is where the source loops forever.
10. **`-threads` reaches the picking, where the source's does not.** The source
   tool applies the setting before `main_` (`TOPPBase.cpp:408-415`), but
   `PeakPickerHiRes.cpp` has no OpenMP anywhere, so in the C++ the setting
   cannot share any picking. It is not inert for the C++ tool — on a fixed set
   of eight CPUs, `-threads 8` buys it 7.1%, measured below — but whatever that
   gain is, it is outside the picker. This port runs the picker's spectrum loop
   on the rayon
   pool `ToolContext::in_thread_pool` sizes from the policy. The written file is
   fixed by the input and the parameters: byte-identical at `1`, `2`, `8`, `16`,
   `32` and `0` workers in the test suite, and byte-identical over the 2.3 GB
   benchmark run at `1`, `8` and `32`, to a build with the `parallel` feature
   off.

   **The pool is opened around the picking, not around the whole tool body**,
   which is where the five earlier ported tools put it (from `Tool::run`; a
   wrapper around `run` would never execute here, because this tool overrides
   `run_io` and that is what `run_with` calls). Picking is this tool's only
   parallel region — reading, the input checks, the summary, the processing
   record and storing are serial in this port and in the source — and scoping
   the pool to it has three consequences, all wanted:

   * every line reaches the real stream *where the source writes it*, because
     the phase that produces it runs on the calling thread. `run_io`'s `out` and
     `err` are not `Send` and cannot cross onto a pool thread, so a body that
     runs entirely on the pool has to collect its lines and write them at the
     end — and then loses them on any path that returns an error, which the
     input warnings and the per-level summary both have. Measured against the
     C++ tool on an input whose store fails: both write the ion mobility warning
     to the error stream before picking and the summary to standard output
     before the store, and both are lost by a collect-and-write-at-the-end body
     (evidence below);
   * the mzML read, which is about 45% of an instrument-scale run and the
     heaviest allocator client in it, keeps the calling thread's malloc arena;
   * the pool exists only while it is used.

   **At one worker no pool is built** and picking runs on the calling thread.
   `-threads 1` is the TOPP default, and a pool of one worker is a pure cost:
   glibc gives the worker its own malloc arena, which measured at +0.6 s on the
   2.3 GB run. Nothing is given up, because the bound a pool buys is a bound on
   rayon work, and the only rayon work here is the picking call, whose width is
   the `Threads` value it is given; at one worker the picker's batch loop maps on
   the calling thread and builds no pool of its own. This is the one place where
   this tool's thread handling differs from the five earlier ones, and
   `tests/topp_threads.rs` records it as `expected_workers`: the sampled worker
   count is exactly `-threads n` for every tool, and **none** at one worker for
   this one. What that suite had to change for this tool is not the assertion but
   the input — it samples the picker on profile records it really picks, because
   the pool now lives only as long as the picking and a centroid-like record is
   copied faster than the sampler polls.

   A second consequence is in the reporting: the `Error::Io` that
   `ToolContext::in_thread_pool` raises when the operating system refuses the
   worker threads is now raised *inside* the body, so `run_io` reports it like a
   picker failure — `Error: Unexpected internal error (cannot start <n> worker
   threads: …)` with `UNKNOWN_ERROR` — where it previously propagated out of
   `run_io` to the framework. The source has no such path at all: libgomp aborts
   the process when it cannot create a thread.
11. **Picking is in place.** The tool picks with
    `PeakPickerHiRes::pick_experiment_in_place_with_threads` rather than the
    borrowing `pick_experiment`: it writes the picked experiment and never reads
    the profile data again, so replacing each record as its centroids appear
    releases that record's profile samples and copies no record it does not
    pick (33,945 of the benchmark run's 40,856 spectra are copied by the
    borrowing form). The two produce the same experiment and the same reports;
    only the peak memory differs, by 886 MB on the benchmark run. The in-place
    form is not atomic, which costs this tool nothing: its only reaction to a
    picking error is to report it and exit without writing an output file.
12. **A record needing header entries the first record did not contribute is
    written with the source's dangling reference.** The consumer writes its
    header from the settings plus the first record, so the `sourceFileList`,
    `dataProcessingList` and `softwareList` it declares are that record's. A
    later record needing different ones cannot be numbered against that header;
    the source numbers its references by the record's own position in the
    stream instead (`MzMLHandler.cpp:5252-5272`, with `dps_` holding the one
    entry `writeHeader_` filled it with), which names nothing the header
    declares. **That is not valid mzML.** `dataProcessingRef` and
    `sourceFileRef` on a `spectrum` are `xs:IDREF` against `xs:ID` on
    `DataProcessingType` and `SourceFileType` (`mzML_1_10.xsd:851`, `:856`),
    and `dataProcessingRef` additionally carries `KEYREF_DPREF`, whose `refer`
    is `KEY_DP_ID`, the `id` of a `dataProcessingList/dataProcessing`
    (`:1064-1071`, `:983-990`). Requested as a C++ issue.

    It is reachable on ordinary data — every `FileMerger` output carries one
    `dataProcessing` per merged part — so this port reproduces it rather than
    refusing. The library type keeps both answers:
    `ReferencePolicy::Checked`, still the default, refuses the record, and
    `ReferencePolicy::SourceDangling`, which this tool's low-memory path
    selects, writes what the source writes (see
    [MS_DATA_WRITING_CONSUMER_SUPPORT](MS_DATA_WRITING_CONSUMER_SUPPORT.md)).
    The dangling identifier is the source's own spelling, `dp_sp_<s>` and
    `sf_sp_<s>`: a bare position in this writer's single zero-padded namespace
    would alias a declared entry instead of dangling.

    Measured on `ibminode06` against the C++ Release build at the pins. On the
    file the Release `FileMerger` builds from 22 copies of
    `PeakPickerHiRes_input.mzML`, the C++ low-memory run exits 0 over all 110
    records; its output declares three referenceable ids and carries 109
    references, of which **105 dangle** — `dp_sp_5` through `dp_sp_109`,
    records 1 to 4 needing none because their `dataProcessing` is the first
    record's (`logs/closediff1_06.log` sections B, C and E). The C++
    **in-memory** run of the same file declares two `dataProcessing` entries
    and dangles nothing. This port's low-memory run now writes all 110 records
    with the same references; before this round it stopped after five with
    `Error: unsupported: record needs a different mzML sourceFileList or
    dataProcessingList than the header written for the first record`.

    The five-record `refs` fixture pins every cell of the rule, because the
    110-record file exercises only one of them. Against the C++ output
    (`logs/closediff2_06.log` section A):

    | record | its `dataProcessing` | its `sourceFile` | C++ writes |
    | --- | --- | --- | --- |
    | 0 | — | — | `sourceFileRef="sf_sp_0" dataProcessingRef="dp_sp_0"`, both declared |
    | 1 | not the first's | the first's | `sf_sp_1`, `dp_sp_1` |
    | 2 | not the first's | not the first's | `sf_sp_2`, `dp_sp_2` |
    | 3 | the first's | not the first's | `sf_sp_3`, no `dataProcessingRef` |
    | 4 | the first's | the first's | `sf_sp_4`, no `dataProcessingRef` |

    Two things to read off it. `dataProcessingRef` is numbered by the record's
    *position*, not by which entry it would have matched — records 1 and 2
    share one `dataProcessing` in the input and come out as `dp_sp_1` and
    `dp_sp_2`. And `sourceFileRef` is renumbered for **every** record after the
    first that carries one, whether or not the source file is the first
    record's: record 4's is, and it still gets the dangling `sf_sp_4`, because
    `MzMLHandler.cpp:5252-5255` never consults the header. This port reproduces
    both (`the_low_memory_mode_writes_the_sources_dangling_references`).

    A chromatogram is the silent case on both sides: `writeChromatogram_`
    (`MzMLHandler.cpp:5879`) writes `id`, `index` and `defaultArrayLength` and
    no reference at all, so a chromatogram whose history the header does not
    declare is written under the list's default and its own history is lost.
    This port does the same
    (`the_source_policy_writes_a_chromatogram_without_any_reference`).

    **What the dangling `sourceFileRef` costs here, and is raised rather than
    decided.** This port's reader refuses an unregistered spectrum
    `sourceFileRef` under *either* dangling-reference policy
    (`src/format/mzml_header/read.rs:113-121`), where the source's warns once
    and carries on (`MzMLHandler.cpp:899-906`, leniency this port deliberately
    does not have). So on an input with per-record source files this port now
    writes a low-memory output it will not read back — exactly as it will not
    read the C++ output of the same run, measured both ways. The dangling
    `dataProcessingRef` has no such problem: the reader's source policy, which
    this tool selects, reads it as an empty history with one warning, so the
    `FileMerger` case round-trips.


## Checked boundaries and evidence

Evidence tier 1 (executed differential) throughout, except the load-option row
of the table, which is tier 4 (Rust-only, on a synthetic document): every other
expectation is a retained upstream output or a C++ output produced by the
product SDK (Debug, core `4fdec46`, decision D7), and the instrument-scale
measurement below is against the optimised C++ Release build. Hashes are in
`tests/data/topp_peak_picker_hi_res_provenance.json`.

| Case | Source of the expectation | What is compared |
| --- | --- | --- |
| `TOPP_PeakPickerHiRes_1`, `_2`, `_5`, `_6` | retained `PeakPickerHiRes*_output.mzML` (`CMakeLists.txt:2523-2549`) | exit 0, stdout, and decoded content: record counts, native ids, MS levels, spectrum types, float array names, every m/z, RT, intensity and array value bit for bit, history lengths 6/2/1/3, and the whole decoded walk with exact numbers |
| `TOPP_PeakPickerHiRes_3` | retained `PeakPickerHiRes_output_lowMem.mzML` (`CMakeLists.txt:2533-2536`), reproduced byte for byte by the C++ Release build (oracle `p4-lowmemory`, `w1_lowmem`) | exit 0, empty stdout and stderr, decoded content, history lengths 6, and byte equality with this port's in-memory run |
| `TOPP_PeakPickerHiRes_4` | retained `PeakPickerHiRes_2_output.mzML` (`CMakeLists.txt:2538-2540`), which the C++ Release build writes in both modes | as above, plus the five chromatogram peak counts |
| the retained C++ pair itself | `PeakPickerHiRes_output.mzML` against `PeakPickerHiRes_output_lowMem.mzML` | exactly one differing byte, `dataProcessingList count="3"` against `"2"` |
| the automatic-mode divergence | C++ Release output `oracle_lowmem_w6_auto.mzML` (oracle `p4-lowmemory`) | 4 centroids in the low-memory mode against the in-memory mode's 33 copied samples and `MS-level 1: 0 / 1`, and decoded equality with the C++ output |
| the absent centroided refusal, `-force` inert | the same C++ output, and the in-memory refusal | exit 0 and 4 centroids where the in-memory mode exits 8 with the source message; the same bytes with and without `-force` |
| the absent input checks | C++ Release outputs `oracle_lowmem_im_peak.mzML`, `oracle_lowmem_unsorted_spectrum.mzML`, and the empty-input run | no ion mobility warning, exit 0 and a zero-byte output for an input without records, and the workflow outputs for the unsorted inputs |
| the index | the retained low-memory output, and this port's own | both are `indexedmzML`; every offset in this port's index lands on the `<spectrum` element whose `id` it names |
| a corrupt input in either mode | the C++ Release build on a truncated, a non-XML and a malformed-base64 input, both process options (oracle `p4-lowmemory`, `logs/fixdiff_06.log`, `trunc`, `garbage`, `badb64`) | exit 3 and `Error: Unable to read file (…)` on both sides and in both modes; a zero-byte low-memory output and no in-memory output |
| a `spectrumList count` that overstates its records | the C++ Release build on the same derived input (`badcount`): exit 0, `count="9"` written over five records | this port writes the same count over a closed, indexed, reloadable document and then exits 8 with the count-mismatch message; its in-memory run is unaffected |
| a failure after the first record | the source's `~MSDataWritingConsumer`/`doCleanup_` contract; the C++ run on that command line dies of SIGSEGV and so cannot be compared (`am1`) | exit 11, and the output left behind is closed, indexed and reloadable, holding the one record written under the first pass's count |
| the batch boundary a failing run leaves | the C++ Release build on the 110-record `batches` fixture with four corruptions at index 4 and at index 104 (`logs/closediff2_06.log` section D, `logs/closediff3_06.log` section C) | malformed base64: both sides write nothing at index 4 and 100 records at index 104; a bad `defaultArrayLength`, a non-numeric `scan start time` and a duplicate native id: the C++ exits 0 with all 110 records in both modes, this port leaves 0 and 100 records and its in-memory run writes no file |
| the source's dangling header references | the C++ Release build on the 110-record `FileMerger` output (105 dangling `dataProcessingRef`s over 110 records) and on the five-record `refs` fixture, whose start tags give the numbering rule for every combination (`logs/closediff1_06.log`, `logs/closediff2_06.log` section A, `logs/closediff3_06.log` section A) | the same references, in the source's own spelling, over the same records; the encoded arrays of this port's two modes; and the reader refusal an unregistered `sourceFileRef` still causes here, against the source's warn-and-continue |
| an `-out` naming an existing directory | the C++ Release build in both modes, with a read-only `-out` and an `-out` under a missing directory as controls (`logs/closediff1_06.log` section F, `logs/closediff3_06.log` section D) | the C++ low-memory run exits 0 with empty standard error having written nothing, its in-memory run exits 5 `Error: Unable to write file (…could not be created. )`, and both controls exit 5 with `Cannot write output file given from parameter '-out'!` on both sides and in both modes; this port exits 8 with the operating system's refusal in both modes |
| `-threads` on the low-memory path | this port at 1, 8 and 32 on two inputs; both implementations at 1, 8 and 32 on the node (oracle `p4-lowmemory`, `threads`) | bit-identical bytes |
| instrument scale, both modes | the 2.3 GB `UK222.mzML` through both implementations and both modes on `ibminode06` | peak RSS, and the byte-identical mzML body between the modes on each side (see *The low-memory mode*) |
| `TOPP_INI_INVALIDVALUE`, `TOPP_CLI_INVALIDVALUE`, `_SECTION` (both), `TOPP_INI_INVALIDNAME`, `TOPP_CLI_INVALIDNAME` | `CMakeLists.txt:107-131` plus the C1 oracle | exit 6 and every diagnostic `ExpectToolFailure.cmake` requires |
| `TOPPWRITEINI_OVERWRITE` | retained `WRITE_INI_OUT.ini` (`CMakeLists.txt:104-106`) | line-level `FuzzyDiff` with the registered `version` whitelist, plus five of the C++ update diagnostics |
| `-write_ini` defaults | P3 oracle `write_ini` | line-level `FuzzyDiff` with the upstream settings and no whitelist |
| C++ INI fed back at `-threads 1`, `16`, `0` | P3 oracle `cpp_ini_threads_*` | exit 0, no diagnostics, decoded equality, and byte-identical Rust output at every thread count; the same run without `-test` |
| `PPHR_6_noforce` | C1 oracle | exit 8, the exact stderr line, no output |
| `PPHR_auto_mode_1`, `auto_mode_1_without_estimation` | C1 and P3 oracles | the explicit refusal (exit 11, nothing written) and the successful run with `signal_to_noise` 0 |
| `PPHR_cli_signal_to_noise_2` | C1 oracle output | 102 centroids in `scan=12663` and decoded equality |
| `empty` | P3 oracle | exit 11, the C++ stderr byte for byte, nothing on stdout, no output |
| `unsorted_spectrum`, `unsorted_chromatogram` | P3 oracle | the loader sorts; the outputs are the workflow 6 and workflow 2 outputs |
| `im_peak` | P3 oracle | the warning byte for byte and the mean ion mobility array bit for bit |
| `notest` | P3 oracle | the processing record: version, action, parameter keys in order and every value but the two paths |
| unsorted refusals | module unit test | the source messages and `INCOMPATIBLE_INPUT_DATA` on constructed experiments |
| the tool's load options | synthetic 42,776,703-byte profile document (tier 4, Rust-only) | the library default refuses its dangling `softwareRef`, the former fixed ceilings refuse its size (`parameter bytes exceed configured limit`), and `PeakPickerHiRes::read_options` reads all 20,000 spectra and picks them end to end |

**Release measurements** (`ibmi` node `dax`, 384 cores, release profile,
`-test -ini <C++ default INI>`; wall time including process start, best of
five). The C++ side is the optimised build at the same source pins
(`../oracle/release-build/manifest.json`), never the Debug product SDK.

| Input | Rust `-threads 1` | Rust `-threads 16` | C++ `-threads 1` | C++ `-threads 16` |
| --- | --- | --- | --- | --- |
| `PeakPickerHiRes_input.mzML` (432 KB, 25,818 raw points), the largest C1 input | 0.088 s | 0.086 s | 0.159 s | 0.155 s |
| 100 copies of it (42.5 MB, 2.6 M points) | 0.467 s | 0.471 s | 0.440 s | 0.422 s |
| 300 copies (127.5 MB, 7.7 M points) | 1.254 s | 1.266 s | 1.002 s | 0.974 s |
| 500 copies (212.6 MB, 12.9 M points) | refused, exit 3 | refused, exit 3 | 1.586 s | 1.515 s |

The Rust output is byte-identical at every thread count, and identical to the
`-threads 1` run; the C++ output is identical across thread counts too, and the
optimised C++ output on the smallest input equals the Debug product SDK's
(`0e63f534…`), so the two builds agree here. The refusal at 500 copies was the
mzML reader's former fixed 10,000,000-point ceiling, which
`PeakPickerHiRes::read_options` replaces with the size-derived ones: after
212.6 MB of consumed XML the point allowance stands at 10,000,000 + 8 per byte,
which is 130 times those 12.9 M points. That table was not re-measured; the
executed evidence for the load path is the instrument-scale run below, on a file
ten times larger again.

### Instrument scale, against the C++ tool

Executed on `ibminode06` (128 cores, 995 GB, a foreign load of 24 to 35 runnable
processes throughout, so the wall times are indicative), against the optimised
C++ Release build at these pins
(`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`). Input:
the benchmark's Q Exactive SILAC profile run `UK222.mzML`, 2,317,975,830 bytes
(sha256 `bd6f6e19…`), 40,856 spectra, one chromatogram, 197,765,338 raw points,
staged on node-local `/scratch`. Both tools were driven by the INI the C++ tool
wrote with `-write_ini` (sha256 `dc0f7b60…`: `threads` 1, automatic mode,
`signal_to_noise` 0, no FWHM), outside `-test`. Drivers, logs and hashes:
`../oracle/topp-peak-picker-scale/`.

| | exit | wall | peak RSS | output |
| --- | --- | --- | --- | --- |
| C++ `PeakPickerHiRes` | 0 | 26.8 s | 3,977,188 KiB | 549,528,678 B |
| Rust, this branch | 8 | 18.1 s | 3,491,772 KiB | none |
| Rust, with the picker ledger lifted | 0 | 38.7 s | 4,412,352 KiB | 535,615,957 B |

The **load path is no longer the limit**: the middle row reads the whole 2.3 GB
file — its peak RSS is within 0.01% of the reader's own measurement of the same
input (3,491,464 KiB,
`../oracle/mzml-reader-scale/hpc_scale_ibminode06_round2.log`), so the load is
the high-water mark of that run —
and then exits 8 with `Error: Unexpected internal error (invalid value: data
array description resource limit exceeded)` from the picker's pooled
acquisition-copy ledger (`src/processing.rs`), which charges all 40,856 spectra
before picking any of them. That ledger is the picker library's, not this
tool's; the third row is the same tool over the sibling lane that lifts it
(`fix/picker-scale`), measured to show what the tool does once it can run, and
is not the state of this branch.

Both implementations print the same per-MS-level summary
(`MS-level 1: 6911 / 6911`, `MS-level 2: 0 / 33945`); the C++ adds progress
logging and a timing line, which this port does not write (native difference 7).

### Instrument scale across thread counts (`perf/peak-picker`)

Same node, same input, same C++ Release build and the same C++-written INI, this
time with `-test` so the processing record carries no wall-clock time and the
written file can be compared by hash. Every configuration was run three times
and both tables are the **median**; the node carried a foreign load of 4.6 to
17.8 runnable processes, recorded per run. The first table pins each
configuration with `taskset` to as many CPUs as it was given workers (`-c 0`,
`0-7`, `0-31`), which is how the port is used; the second table pins **every**
run to the same eight CPUs, so that the worker count is not confounded with the
CPU count. Four binaries:

* **C++** — `openms4-release-bc9cc12-c19e494-174b576`.
* **Rust before** — the baseline build of `fabd4b9`.
* **Rust after** — this branch, same recipe and toolchain
  (see *Rebuilding the measured binaries* below).
* **Rust after, `parallel` off** — this branch built
  `--no-default-features --features mzml,paramxml`.

| | wall `-threads 1` | `-threads 8` | `-threads 32` | peak RSS | output |
| --- | --- | --- | --- | --- | --- |
| C++ `PeakPickerHiRes` | 28.33 s | — | — | 3,953,640 KiB | 549,526,348 B |
| Rust before (`fabd4b9`) | 37.10 s | — | — | 4,414,500 KiB | 535,613,726 B |
| Rust after | 37.32 s | **24.40 s** | **22.53 s** | 3,527,720 KiB (t1), 3,519,504 KiB (t32) | 535,613,726 B |
| Rust after, `parallel` off | 36.19 s | — | 36.48 s | 3,527,716 KiB | 535,613,726 B |

The C++ tool is measured at more than one worker in the **second** table only,
because a per-worker pinning cannot say anything about its `-threads`.

**The C++ worker count, without the CPU-count confound.** Pinning each
configuration to as many CPUs as it has workers confounds *more workers* with
*more CPUs*, and the first record of this benchmark read the C++ tool's better
time at eight as "the OpenBLAS server threads having room to spin". That
mechanism was never measured, and this benchmark cannot see it. Repinned so
that `-threads` is the only thing that changes:

| on `taskset -c 0-7` | `-threads 1` | `-threads 8` |
| --- | --- | --- |
| C++ `PeakPickerHiRes` | 27.39 s | 25.45 s |
| Rust after | 36.76 s | 23.94 s |

The C++ gain at eight workers is real — 1.94 s, 7.1%, with no extra CPU to spin
on — and it is **not** picking being shared: `PeakPickerHiRes.cpp` at `bc9cc12`
contains no `#pragma omp` at all (verified by grep; its spectrum loop at `:504`,
chromatogram loop at `:548` and on-disc loop at `:584` are plain `for`
statements). The one OpenMP region this run passes through is the mzML loader,
`MzMLHandler.cpp:205-206` — `#pragma omp parallel for` over
`populateSpectraWithData_`, and `:287-288` for chromatograms — and the C++
tool's minor page faults do rise with the setting, 1,144,276 at one worker to
1,176,911 at eight, so a thread team is created somewhere. Which region the
1.94 s comes from is not separated by this benchmark and is not claimed here.
What the same eight CPUs also show is that the port's parallel picking is worth
more than the C++ tool's whole `-threads` response on this input: 36.76 → 23.94 s
against 27.39 → 25.45 s.

**Bit-identity.** All **25** Rust runs behind the two tables — nine
configurations, three repetitions each but one: this branch at 1, 8 and 32
workers, this branch at one worker under `MALLOC_ARENA_MAX=1`, this branch at 1
and 8 on the fixed CPU set, `fabd4b9` at 1, and the `parallel`-off build at 1
and 32 — wrote the same 535,613,726 bytes, sha256
`bb13eecfe092a272b08ddc71feec3780c7bc876e8847e8a45b173cda9d2dad52`. That is the
hash the first round of this branch recorded and the hash `fabd4b9` wrote: the
worker count changes nothing, the `parallel` feature changes nothing, the arena
setting changes nothing, and rescoping the pool to the picking changed nothing.
The 22,776,198 centroids are still the C++ Release build's, bit for bit.

**Scaling, and its ceiling.** The Rust tool goes 37.32 → 24.40 → 22.53 s,
**1.66x at 32 workers**, and passes the C++ tool between one and eight workers.
The ceiling is Amdahl's and is known: the pick loop is about 47% of the run and
the serial mzML read is about 45%, so no worker count takes this below roughly
20 s while the reader is serial. The reader is also the larger half of the
remaining gap to C++ (+5.4 s of the +8.6 s measured by the profiling lane) and
is not this lane's file.

**Memory.** Peak RSS falls from 4,414,500 KiB to 3,527,720 KiB, 886 MB and 20%,
below the C++ tool's 3,953,640 KiB — the in-place entry point and the
metadata-only record construction. It is **flat in the worker count**:
3,527,720 KiB at one worker against 3,519,504 KiB at 32. The parallel path holds
one bounded batch of centroids beyond what a serial pick holds
(`PARALLEL_BATCH_RECORDS` / `PARALLEL_BATCH_POINTS`), not one per worker.

**`-threads` reaches the workers.** The running tool's `/proc/<pid>/task/*/comm`
was sampled every 50 ms over real 2.3 GB runs:

| `-threads` | 1 | 2 | 8 | 32 | 0 |
| --- | --- | --- | --- | --- | --- |
| tasks at their peak | 1 | 3 | 9 | 33 | 129 |
| of them `openms-<i>` workers | **0** | 2 | 8 | 32 | 128 |

Every count above one starts exactly that many workers named `openms-<i>`, the
names `ToolContext::in_thread_pool` gives its pool, beside the main thread and
nothing else; the pool is the only source of threads in the process, and
`-threads 0` reaches all 128 processors of the node. At `-threads 1` the process
is **one thread**: no pool is built and the picking runs on the calling thread
(native difference 10). All five runs wrote the same
`bb13eecf…` output. These counts are exact because this pick takes about
fifteen seconds. `tests/topp_threads.rs` samples the same executables the same
way on a synthetic input and holds this tool to exactly `-threads n` workers and
to none at one (`expected_workers`); its input for this tool is profile data
rather than the centroid-like records the other five tools share, because the
pool is open only while the picking runs. On that profile input the picking
takes 654 ms at one worker, 339 ms at two and 184 ms at four in a debug build on
the gate host, against a 200 µs poll — hundreds of samples — where a
centroid-like record is copied in single-digit milliseconds and can be over
between two samples, which is what once made that assertion flaky.

**The error paths write what the C++ tool writes.** Executed on the same node,
same fixture (`p3_im_peak.mzML`, per-peak ion mobility, tool defaults), same
`-out`: an existing **directory**, which the framework's writability pre-check
accepts (`file::writable` probes a directory and answers yes) and which the
store then cannot use. It is the one failure *after* the summary that is easy to
provoke from a command line: the reviewer's probe showed that an output under a
non-existent directory is stopped by the pre-check in both implementations and
never reaches the store, so this is what is left.

| `-out` = an existing directory | standard error | standard output | exit |
| --- | --- | --- | --- |
| C++ Release | the ion mobility warning, then `Error: Unable to write file (…could not be created. )` | the per-level summary | 5 |
| Rust `fabd4b9` | the warning, then `Error: Unexpected internal error (Is a directory (os error 21))` | the summary | 8 |
| Rust, this branch | byte for byte as `fabd4b9` | the summary | 8 |

Both lines survive the failure in all three, which is the point: the warning is
written where the input is checked and the summary before the file is stored.
The first version of this branch ran the whole body on the pool and so had to
collect its lines and write them at the end, which dropped them with any error
that returned from the body. Two tests pin both halves —
`the_summary_is_written_before_the_output_is_stored`, which is this same
directory case run in process, and
`a_refusal_after_the_input_checks_keeps_the_ion_mobility_warning`, which is the
`auto_mode 1` refusal, a different post-check failure. Both were executed
against that version with only the tests applied: the first found standard
output carrying nothing but the INI version warning, the second found standard
error carrying the refusal without the mobility warning ahead of it. Both pass
here. The exit code and the wording of the store failure itself
differ between the C++ tool and the port and always have (native difference 9:
an unmapped error takes `UNKNOWN_ERROR`); what this case is evidence for is
where the two lines are written, not how the store failure is reported. A
control run of the same fixture with a writable `-out` gives exit 0 and the same
warning and summary from all three builds.

**What the single-thread changes are worth, measured without the node.** Wall
time on this node cannot resolve a half-second on a 37 s run: a paired,
alternating A/B of eight reps ran while the foreign load swung between 16 and 41
and the same binary spread over 35.3 s to 42.5 s, so that A/B is reported as
inconclusive rather than as a number. Instruction and cache counts are
deterministic and do not care about the load. Callgrind (valgrind 3.27.1,
`--cache-sim=yes --branch-sim=yes`) over the 682-spectrum, 49.7 MB slice
`part_part01of60.mzML` with the same INI, at `-threads 1`, all three builds
pinned to one CPU:

| at `-threads 1` | instructions | data refs | D1 misses | branches |
| --- | --- | --- | --- | --- |
| `fabd4b9` | 5,772,056,967 | 2,269,009,920 | 25,296,992 | 1,333,319,077 |
| this branch, `parallel` off | 5,665,953,334 | 2,154,334,607 | 22,950,812 | 1,274,789,232 |
| | **-1.84%** | **-5.05%** | **-9.27%** | **-4.39%** |
| this branch, **as shipped** (`parallel` on) | 5,673,177,994 | 2,159,056,928 | 22,974,469 | 1,275,134,417 |
| | **-1.71%** | **-4.85%** | **-9.18%** | **-4.36%** |

The third row is the configuration an ordinary `PeakPickerHiRes -in … -out …`
runs in, and it is the row that matters: **the default keeps the saving.** It
did not before — with the pool around the whole body, the `parallel` build at
`-threads 1` measured 5,770,823,494 instructions, -0.02% against `fabd4b9`, so
the single-thread work was given back in full unless a worker count above one
was used. What is left between the second and third rows, 7.2 M instructions or
0.13%, is rayon being linked and the one-worker decision being taken.

All three runs wrote the same 5,841,655 bytes (sha256 `2670a1d738bf135e…`). The
data-reference and D1-miss figures are where the two removed passes live: the
whole-record clone that copied every profile sample of every picked spectrum
only to drop it, and the second `MSSpectrum::validate` over samples the
experiment-level pass had already visited. Instructions fall by less than data
references do, which is what removing copies rather than computation looks like.

**The one-worker cost, and where it went.** The first version of this branch ran
the whole tool body on the pool and built one even for `-threads 1`, and that
cost 0.68 s on this run — entirely glibc's second malloc arena, since
`MALLOC_ARENA_MAX=1` recovered it exactly. The pool is now opened around the
picking and not built at all at one worker (native difference 10), and the cost
is gone. The deterministic counter says so more clearly than wall time can:
minor page faults over the 2.3 GB run at `-threads 1` are **1,701,679** for this
branch against **1,701,678** for the same code built without `parallel`, where
`fabd4b9` needs **1,850,830**; and running this branch under
`MALLOC_ARENA_MAX=1` now changes nothing at all (1,701,681, and the same output
hash). In wall time the shipped default is 37.32 s against `fabd4b9`'s 37.10 s,
a 0.22 s difference on a 37 s run and inside `fabd4b9`'s own 0.85 s spread over
its three repetitions, where the first version of this branch was 0.66 s slower
with no overlap.

What is still on the table is crate-wide and not this lane's: raising the
allocator's mmap, trim and top-pad thresholds (`MALLOC_MMAP_THRESHOLD_`,
`MALLOC_TRIM_THRESHOLD_`, `MALLOC_TOP_PAD_`) took the earlier measurement to
35.38 s, a further 2.2 s that a `Cargo.toml`/`src/lib.rs` allocator decision
would collect for every tool.

**Evidence.** Drivers, per-run logs, the results table, the probe outputs and
the callgrind output of the round above are archived at
`/ceph/ibmi/abi/oliver/bench/openms4/results/2026-09-15-pph-parallel-minors/`
(`wbuild.sh`, `wbuild.log`, `wrun.sh`, `wrun.log`, `wresults.tsv`, `probe/`,
`callgrind/`). The round that first measured the parallel loop, before the pool
was rescoped, is beside it in
`…/2026-09-15-pph-parallel/` (`bench_remote.sh`, `arena_probe.sh`,
`ab_probe.sh`, `ir_probe.sh`, `results.tsv`, `bench_run.log`, `ab_probe.log`,
`ir_probe.log`, `ir_{base,ser,new}.log`, `bench_build.log`).

**Rebuilding the measured binaries.** `cargo build --release --locked
--offline`, rustc 1.96.0, no `RUSTFLAGS`, from this worktree — plus one node
detail the first record left out, without which the recipe does not reproduce:
on `ibminode06` `/usr/local/bin/cc` is a Ceph-quota shell script that shadows
the C compiler, so rustc's link step exits 0 and writes no binary at all.
`wbuild.sh` puts a `cc` → `/usr/bin/gcc` symlink first on `PATH`; without it
`cargo build` reports success and produces nothing. A binary sha256 is recorded
as an identifier for the artefact that was measured, not as a provenance token.
What is on record from this node, with this recipe, is three independent
rebuilds and one unexplained difference between them: a rebuild of `3b943e4`
reproduced its `parallel` binary byte for byte (6,135,424 B, sha256
`8680d951a6595840…`) and a rebuild of `fabd4b9` reproduced its
(`a37a043fe652dcc7…`), while an earlier rebuild of the same `3b943e4` tree gave
a 6,135,416 B `parallel` binary — eight bytes smaller, same behaviour, same
output hash — and reproduced the `parallel`-off binary of that round exactly.
That is one observation, not the outcome to expect, and its cause was not
identified; it is recorded so that the next rebuild can recognise it if it
recurs. The reviewer of the follow-up round supplied the mechanism class: on
this toolchain a rustdoc-only change to a single file moved this tool's release
binary by 16 bytes on the gate host — `.strtab` and the build-id note — while
`.text`, `.rodata`, `.data`, `.data.rel.ro`, `.eh_frame` and
`.gcc_except_table` stayed byte-identical. So check a small size or hash
difference section by section (`readelf -S`, `objcopy --only-section`, `cmp -l`
on the stripped images) before treating it as a behaviour difference. Either
way the provenance is the commit plus the recipe, and the hash only names the
file.

**Comparison of the two outputs.** The C++ `FuzzyDiff` from the same prefix
(`-ratio 1.001 -absdiff 1e-5`) fails at line 1, column 31 — the XML declaration,
ISO-8859-1 against UTF-8 — so it never reaches the data: the container
difference is documented, not compared (native difference 5), and D6 makes the
decoded content the contract. Decoded (`../oracle/topp-peak-picker-scale/decoded_compare_probe.rs`):

- 40,856 spectra on both sides, with equal native ids, MS levels and peak
  counts, and **22,776,198 centroids whose every m/z and every intensity are
  bit-identical**.
- Spectrum retention times differ in 10,671 of 40,856 records, by at most
  9.09e-13 s (one to two ULP of `f64`) on times of 60 to 4,400 s.
- The picked TIC **chromatogram is bit-identical**: 8,174 points on both sides,
  every retention time and every intensity equal. This run first showed it
  differing (8,173 of 8,174 retention times, at most 3.19e-3 s; 7,891
  intensities beyond 1e-6 relative, the worst 1.75e-3 at 4,393.5 s), which lane
  `fix/picked-chromatogram` then root-caused in the mzML reader rather than in
  the picker: the source narrows a unit-converted 32-bit time array back to
  `f32` (`MzMLHandlerHelper.cpp:217-222`), which this input's TIC time array —
  32-bit float in minutes — hits, and the spline apex amplifies the 9.5e-7 s
  per-point difference. `ReadOptions::source_time_array_precision`, which
  `ReadOptions::source()` sets and this tool therefore uses, reproduces the
  narrowing; the library default keeps the precision. See
  [MZML_SUPPORT](MZML_SUPPORT.md) and native difference 13 in
  [PEAK_PICKING_SUPPORT](PEAK_PICKING_SUPPORT.md).

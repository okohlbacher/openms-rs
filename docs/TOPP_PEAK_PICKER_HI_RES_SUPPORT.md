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
| `class PPHiResMzMLConsumer` (`processSpectrum_`, `processChromatogram_`, members `pp_`, `ms_levels_`) | not ported | The low-memory path; package P4-PICKER-LOWMEM ports it against `format::ms_data_writing_consumer` and `mzml::transform`. |
| `registerOptionsAndFlags_` | `Tool::register` | `-in`, `-out` (mzML only), advanced `-processOption` restricted to `inmemory,lowmemory`, subsection `algorithm`. |
| `getSubsectionDefaults_(section)` | `Tool::subsection_defaults` | Returns `PeakPickerHiRes::defaults()` for every section, as the source ignores its argument. |
| `doLowMemAlgorithm(pp)` | not ported | Refused, see below. |
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
| `-processOption lowmemory` | refused explicitly, `INCOMPATIBLE_INPUT_DATA` (package P4) |
| `-force` | ported (`check_spectrum_type = !force`) |
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
  2.3 GB benchmark run at 1, 8 and 32). See native difference 9.

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
3. **`-processOption lowmemory` is refused** with `Error::Unsupported`
   (`INCOMPATIBLE_INPUT_DATA`) and a message naming package P4, before anything
   is read. The source runs its `PPHiResMzMLConsumer` there, whose spectrum
   selection deliberately differs from the in-memory mode.
4. **Container differences are documented, not compared** (decision D6). The
   C++ output is an `indexedmzML` with an ISO-8859-1 declaration, the software
   alias `MS:1002135 TOPP PeakPickerHiRes` and a `dataProcessingList count`
   computed as `max(1, histories + float arrays)` (CPP-019); this port writes a
   plain mzML in UTF-8 with the exact software name and one `dataProcessing`
   entry. Both decode to the same content, which is what the tests compare.
5. **The ion mobility peak type in the warning is `im_profile`.** The source
   prints `imPeakTypeToString(spec.getIMPeakType())`; the native spectrum has no
   stored peak type. `im_profile` is what the source mzML reader stores for ion
   mobility data without `MS:1003441` (`MzMLHandler.cpp:253-255`), which is the
   oracle's text; an input carrying that term would read `im_centroided` in the
   source.
6. **No debug dump and no progress logging**, as listed in the API mapping.
7. **Bounded work on the input.** The C++ tool has no resource ceilings; this
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
8. **A picker failure that is not the centroided refusal** is reported as
   `Error: Unexpected internal error (<reason>)` with `UNKNOWN_ERROR`, the code
   `TOPPBase` gives an unmapped exception, rather than the framework's default
   mapping of `Error::InvalidValue` to `ILLEGAL_PARAMETERS`: these are data and
   resource conditions, not parameter errors. One of them, a non-converging
   FWHM bisection, is where the source loops forever.
9. **`-threads` reaches the picking, where the source's does not.** The source
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
   `tests/topp_threads.rs` records it as `expected_workers`: for the five the
   sampled worker count is exactly `-threads n`, for this tool it is at most
   `n`, and none at one worker. The upper bound rather than the exact count,
   because the pool now lives only as long as the picking, and on that suite's
   few-MB centroid-like input the picking can be over between two samples.
10. **Picking is in place.** The tool picks with
    `PeakPickerHiRes::pick_experiment_in_place_with_threads` rather than the
    borrowing `pick_experiment`: it writes the picked experiment and never reads
    the profile data again, so replacing each record as its centroids appear
    releases that record's profile samples and copies no record it does not
    pick (33,945 of the benchmark run's 40,856 spectra are copied by the
    borrowing form). The two produce the same experiment and the same reports;
    only the peak memory differs, by 886 MB on the benchmark run. The in-place
    form is not atomic, which costs this tool nothing: its only reaction to a
    picking error is to report it and exit without writing an output file.

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
logging and a timing line, which this port does not write (native difference 6).

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
(native difference 9). All five runs wrote the same
`bb13eecf…` output. These counts are exact because this pick takes about
fifteen seconds; `tests/topp_threads.rs` checks the same executables on a
few-MB input, where the pool is short-lived, and so holds this tool to *at
most* `-threads n` workers and none at one (`expected_workers`).

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
differ between the C++ tool and the port and always have (native difference 8:
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
picking and not built at all at one worker (native difference 9), and the cost
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
as an identifier for the artefact that was measured, not as a provenance token:
an independent rebuild of the same tree reproduced the `parallel`-off binary
byte for byte but not the `parallel` one (6,135,416 B against 6,135,424 B, same
behaviour, same output hash), so the provenance is the commit plus the recipe
and the hash only names the file.

**Comparison of the two outputs.** The C++ `FuzzyDiff` from the same prefix
(`-ratio 1.001 -absdiff 1e-5`) fails at line 1, column 31 — the XML declaration,
ISO-8859-1 against UTF-8 — so it never reaches the data: the container
difference is documented, not compared (native difference 4), and D6 makes the
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

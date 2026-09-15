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
| `FileHandler().loadExperiment(in, exp, {MZML}, log_type_)` | `FileHandler::load_experiment_with_read_options(in, &[FileType::MzMl], &PeakFileOptions::default(), &mzml::ReadOptions { source_dangling_references: true, .. })` | Decision D10: the tool path accepts the dangling header references `MzMLHandler` drops, which `TOPP_PeakPickerHiRes_5` needs; library defaults stay strict. |
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
- **Serial.** `PeakPickerHiRes.cpp` has no OpenMP and neither has the picker, so
  `-threads` changes no output byte (tested at 1, 2, 16 and 0).

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
7. **Bounded work on the input.** The mzML reader's ceilings apply
   (`mzml::ReadOptions`: 512 MiB of XML, 10,000,000 raw points, 512 MiB of
   decoded arrays), and the picker's own (`max_points` 1,000,000 per record,
   `max_work` 10,000,000, the acquisition-copy ledger). The C++ tool has no
   ceilings. An input beyond a ceiling is refused before anything is written;
   the point ceiling is reached at about 213 MB of profile mzML, which a
   benchmark on instrument-sized data will hit (see *Checked boundaries and
   evidence*).
8. **A picker failure that is not the centroided refusal** is reported as
   `Error: Unexpected internal error (<reason>)` with `UNKNOWN_ERROR`, the code
   `TOPPBase` gives an unmapped exception, rather than the framework's default
   mapping of `Error::InvalidValue` to `ILLEGAL_PARAMETERS`: these are data and
   resource conditions, not parameter errors. One of them, a non-converging
   FWHM bisection, is where the source loops forever.

## Checked boundaries and evidence

Evidence tier 1 (executed differential) throughout: every expectation is a
retained upstream output or a C++ output produced by the product SDK
(Debug, core `4fdec46`, decision D7). Hashes are in
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
(`0e63f534…`), so the two builds agree here. The refusal at 500 copies is the
mzML reader's 10,000,000-point ceiling, reported as
`Error: Unable to read file (parse error on line 0: peak count exceeds
configured limit)` with `INPUT_FILE_CORRUPT`. A benchmark over
instrument-sized profile data needs that ceiling raised for the tool path (and
a code other than 3 for a resource refusal); the ceilings are shared policy
across the bundle's tools, so this is an integrator decision, not a local one.

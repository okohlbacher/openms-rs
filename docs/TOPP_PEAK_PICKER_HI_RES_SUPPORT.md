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
- The picked TIC **chromatogram differs beyond round-off**: 8,174 points on both
  sides, but 8,173 of 8,174 retention times differ (at most 3.19e-3 s), and
  7,891 of 8,174 intensities differ by more than 1e-6 relative, 71 of them by
  more than 1e-3, the worst 1.75e-3 (252,284,752 against 252,727,072 at
  4,393.5 s). The retained workflow 2 fixture (five short chromatograms) is
  bit-exact, so this shows up only here. What has **not** been established is
  where the two part: a candidate is the input's TIC time array, which is
  32-bit float in **minutes** where every spectrum's retention time is a decimal
  `cvParam`, so the two conversions to seconds could differ before the picker
  sees the data and the apex interpolation would amplify that; that is a
  hypothesis, not a measurement. The chromatogram path is the picker library's
  (`PeakPickerHiRes::pick_chromatogram`, package P1), not this tool's, so this
  is reported and not chased further here.

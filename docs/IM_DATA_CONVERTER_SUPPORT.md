# IM data converter support

Native coverage of `IMDataConverter::splitByFAIMSCV` from
`IONMOBILITY/IMDataConverter.h` and `IONMOBILITY/IMDataConverter.cpp` at Core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Work package B8-IMSPLIT of the early
TOPP bundle; decision D5 limits it to this one member, so the header stays
**partial**.

| Artifact | Path |
| --- | --- |
| Implementation | `src/kernel/im_data_converter.rs` |
| Tests | `tests/im_data_converter.rs` |
| Fixtures | `tests/data/im_data_converter/oracle_cases.tsv`; `tests/data/im_data_converter/c2_split_records.tsv` |
| Reused fixtures (in place) | `tests/data/faims_helper/IM_FAIMS_test.mzML` (A1); `tests/data/mzml_mobility/FeatureFinderCentroided_1_input.mzML`, `FAIMS_test_data.mzML`, `FAIMS_CV-60C_V-45_Interleaved.mzML` (A3) |
| Manifest | `tests/data/im_data_converter_provenance.json` |
| Oracle driver | `../oracle/im-data-converter/` (outside this repository) |

The module lives in `kernel`, next to `faims_helper`, `spectrum_mobility` and
`experiment_mobility`. It uses `kernel -> metadata` (the drift time unit) and
`kernel -> format` (the C++ stream number formatting of
`format::file_info::text_format::ostream_g`, for one warning text); both edges
exist already, so `tools/check_module_cycles.py` records no new edge.
`kernel -> processing`, which would close a cycle, is not used.

The consumer in the pinned TOPP package is FeatureFinderCentroided
(`FeatureFinderCentroided.cpp:238`, topp `174b576`), which splits its input and
runs FeatureFinderAlgorithmPicked once per group.

## API mapping

| C++ member | Rust |
| --- | --- |
| `class IMDataConverter` | `kernel::im_data_converter::ImDataConverter`, a zero-sized unit struct: the source class holds no state and its members are static |
| `IMDataConverter()`, `~IMDataConverter()` (implicit) | derived `Default`; no drop glue |
| `static std::vector<std::pair<double, MSExperiment>> splitByFAIMSCV(PeakMap&& exp)` | `ImDataConverter::split_by_faims_cv(&mut MSExperiment) -> Result<FaimsSplit>` |
| its return value `std::vector<std::pair<double, MSExperiment>>` | `FaimsSplit::groups: Vec<FaimsGroup>` |
| `std::pair<double, MSExperiment>` | `FaimsGroup { key: FaimsGroupKey, experiment: MSExperiment }` |
| the `double` key, NaN for non-FAIMS data | `FaimsGroupKey::{NotFaims, Voltage(CompensationVoltage)}`; `FaimsGroupKey::volts()` returns the source `double` (NaN for `NotFaims`), `FaimsGroupKey::voltage()` the `Option` |
| the rvalue argument, moved-from and cleared afterwards | `&mut MSExperiment`, left as `MSExperiment::default()` on success and unchanged on error |
| its `OPENMS_LOG_INFO` and `OPENMS_LOG_WARN` records (`IMDataConverter.cpp:35`, `:60`, `:81`), and the forwarded `FAIMSHelper.cpp:52` warning | `FaimsSplit::messages: Vec<FaimsSplitMessage>`; `FaimsSplitMessage::{level, text, spectrum_index}`, `Display`; `FaimsSplitLogLevel::{Info, Warning}`; the texts as `ImDataConverter::NO_COMPENSATION_VOLTAGES_INFO`, `UNEXPECTED_COMPENSATION_VOLTAGE_WARNING` (prefix) and `SPECTRUM_WITHOUT_COMPENSATION_VOLTAGE_WARNING` |
| the spectra and chromatograms `exp.clear(true)` destroys (`IMDataConverter.cpp:84`) | `FaimsSplit::skipped_spectra`, `FaimsSplit::dropped_chromatograms` (native) |
| native only | `ImDataConverter::MAX_SPECTRA`; `FaimsSplit::has_faims()` |
| used: `FAIMSHelper::getCompensationVoltages` | `kernel::faims_helper::FaimsHelper::get_compensation_voltages` (A1, reused unchanged) |
| used: `MSExperiment::getExperimentalSettings`, `getSqlRunID` | `MSExperiment::settings`, `MSExperiment::sql_run_id` |
| `static MSExperiment reshapeIMFrameToMany(MSSpectrum im_frame)` | not ported (outside D5) |
| `static std::tuple<std::vector<MSExperiment>, Math::BinContainer> splitExperimentByIonMobility(MSExperiment&& in, UInt number_of_IM_bins, double bin_extension_abs, double mz_binning_width, MZ_UNITS mz_binning_width_unit)` | not ported: needs the unported `SpectraMerger` |
| `static MSExperiment reshapeIMFrameToSingle(const MSExperiment& in)` | not ported (outside D5) |
| `static void setIMUnit(DataArrays::FloatDataArray& fda, const DriftTimeUnit unit)` | not ported. It forwards to `IMDataArrayUtils::setIMUnit`, which is not ported either |
| `static bool getIMUnit(const DataArrays::FloatDataArray& fda, DriftTimeUnit& unit)` | not ported as a public function. It forwards to `IMDataArrayUtils::getIMUnit`, whose name-to-unit rule exists privately in `src/kernel/spectrum_mobility.rs` (`array_unit`) behind `MSSpectrum::im_data` |
| file-local `annotateAsIM`, `processDriftTimeStack` | not ported (helpers of the unported members) |

## Preserved source conventions

All verified against the executed oracle; the case names are the driver's.

- **The voltages are `FAIMSHelper::getCompensationVoltages`'s**: every FAIMS
  drift time except the sentinel `-1`, and its missing-voltage warning comes
  first (`sentinel_faims_spectrum`, `only_sentinel`).
- **No voltages, no split.** The whole input, chromatograms and settings
  included, is the single group, keyed NaN in the source and `NotFaims` here,
  with the information message (`ffc1_input`, `class_non_faims`,
  `settings_and_chromatograms_non_faims`). An empty experiment is one empty group
  (`empty`), and an experiment whose only FAIMS drift time is the sentinel is one
  unsplit group with both messages (`only_sentinel`).
- **Ascending voltage order**, whatever the run order (`interleaved_descending_cvs`,
  `infinities`: `-inf`, `-45`, `+inf`). Adjacent doubles are separate groups
  (`adjacent_doubles`).
- **Signed zero.** `-0.0` and `+0.0` are one group keyed with the sign stored
  first; each spectrum keeps its own drift time bits
  (`signed_zero_negative_first`, `signed_zero_positive_first`).
- **A FAIMS spectrum joins its own voltage at any MS level** and becomes the
  context (`ms2_with_own_cv`: an MS2 at `-60` after an MS1 at `-45` takes the
  following MS2 to `-60`).
- **The context rule** (`IMDataConverter.cpp:52-82`): a spectrum of any other
  unit joins the context group only when its MS level is above 1. MS1 and MS
  level 0 are skipped (`ms_level_zero_and_three`, `unit_none_with_value`), MS3 is
  assigned like MS2. MS2 before any FAIMS spectrum is skipped
  (`ms2_before_any_faims`). A non-FAIMS spectrum never changes the context, so
  an MS1 without voltage between FAIMS spectra does not detach the MS2 after it
  (`cv_less_ms1_keeps_context`).
- **The sentinel is an undetected voltage, not a missing context.** A FAIMS
  spectrum at `-1` among real voltages is skipped with the unexpected-voltage
  warning, it becomes the context, and the MS2 after it is skipped as well
  (`sentinel_faims_spectrum`).
- **Settings copies.** Every voltage group holds a copy of the input's
  `ExperimentalSettings` (`IMDataConverter.cpp:46`) and the sqMass run id, which
  the source stores as a meta value of the settings (`MSExperiment.cpp:768-780`)
  and the port as `MSExperiment::sql_run_id` (`settings_and_chromatograms_faims`:
  comment, date, identifiers, user parameter and run id 42 in both groups).
- **No chromatograms in voltage groups** (`faims_test_data`,
  `settings_and_chromatograms_faims`).
- **The input is consumed**: empty afterwards, chromatograms gone, settings reset
  as `clear(true)` resets them.
- **Message texts**, including the source's wording "Not FAIMS compensation
  voltages found", and the voltage of the unexpected-voltage warning formatted as
  `std::ostream << double` at precision 6 (`-1`).

## Native differences

- **Ranges are not stale.** The source builds voltage groups with `addSpectrum`
  and never calls `updateRanges`, so `spectrumRanges().byMSLevel(1)` throws
  `InvalidValue` "No ranges for this MS level" on every one of them (recorded for
  all 35 voltage groups of the oracle, and by C2 facts 1a and 1b), and
  FeatureFinderAlgorithmPicked, which reads those ranges
  (`FeatureFinderAlgorithmPicked.cpp:242`), makes FeatureFinderCentroided exit 8
  on every FAIMS input. The port keeps no range cache
  (`MSExperiment::spectrum_range_manager` computes ranges on demand), so every
  group's ranges are its own. The test asserts the recorded C++ `true` and the
  Rust MS1 ranges side by side. The crash is not emulated.
- **Skipped spectra and FAIMS chromatograms are returned**, not destroyed:
  `FaimsSplit::skipped_spectra` (one per skip message, in input order) and
  `FaimsSplit::dropped_chromatograms`. The groups are the source's.
- **Messages are returned, not logged**, one record per event. The source log
  stream prints a repeated line once and later `<line> occurred N times`, with a
  two-entry repetition cache (`LogStream.cpp:184-253`, `:299-360`). The test
  feeds the returned texts through an emulation of that cache and compares the
  printed lines with the captured C++ lines.
- **NaN voltages are refused** with `Error::InvalidValue` by
  `FaimsHelper::get_compensation_voltages`, before anything moves. The source
  puts the NaN into `std::set<double>` and `std::map<double, MSExperiment>`:
  executed, a NaN first leaves the whole input unsplit under the NaN key with a
  spurious missing-voltage warning (`nan_first`); a later NaN spectrum is found
  in the `-50` group by `std::map::find` and the MS2 after it is skipped
  (`nan_middle`, `nan_last`). The test asserts the refusal, the unchanged input
  and the recorded C++ groups.
- **The input's moved-from state** after the unsplit case is unspecified in the
  source: its settings differed from the default in
  `settings_and_chromatograms_non_faims`. The port always leaves
  `MSExperiment::default()`.
- **Atomicity and ceilings.** Everything that can fail happens before a spectrum
  moves: the voltage scan (`MAX_SPECTRA`, 100,000,000, the `FaimsHelper`
  ceiling), `ExperimentalSettings::validate` on the settings about to be copied
  (only when there are voltages; a sample tree deeper than 64 is refused where
  the source copies anything), and fallible reservations of the destination
  list, each group's spectrum list at its exact size, the group list and the
  skipped list. On error the input is unchanged.
- **The NaN key is explicit.** `FaimsGroupKey::NotFaims` replaces a NaN `double`
  that cannot be compared, ordered or hashed; `volts()` gives the NaN back.

## Performance

Serial, as the source. One voltage scan (A1's, `O(n log k)` for `k` voltages),
one classification pass that looks each FAIMS voltage up by binary search in the
sorted voltage list, and one move pass. Spectra are moved, never cloned; the
only per-group copy is the settings, as in the source. Input without voltages is
moved into its group with no copy at all, which is the FeatureFinderCentroided
path for ordinary data.

## Checked boundaries and evidence

| Evidence | Tier | What it covers |
| --- | --- | --- |
| `IMDataConverter_test.cpp:31-76`, `:201-269` literals | 3 (source review) | the constructor and destructor sections; `IM_FAIMS_test.mzML` split into -65/-55/-45 with 4/9/6 spectra, drift times equal to the key, the input emptied, group 1 dated 2019-09-07T09:40:04; MS2 assigned to the last seen voltage (3 and 2 spectra, 2 and 1 MS2); one group for a millisecond-unit spectrum |
| `oracle_cases.tsv`, file cases | 1 (executed differential) | six mzML files through `MzMLFile::load`: per spectrum the native ID, MS level, drift time bits and unit the C++ reader produced (the Rust reader agrees for all 369), then the split of each (see below). The two C2 FAIMS fixtures are rebuilt in the test from the FFC_1 input |
| `oracle_cases.tsv`, synthetic cases | 1 (executed differential) | the two class-test experiments and 18 edge cases listed under the preserved conventions and native differences |
| `oracle_cases.tsv`, per case | 1 | group count; keys bit for bit; sizes; membership and order by native ID, MS level, drift time bits and unit; chromatogram counts; settings equality with the input, date, comment and run id; the emptied input; `byMSLevel(1)` throwing; every printed log line |
| `c2_split_records.tsv` | 1 (executed differential, extracted) | C2's `faims_corrected` source sequence on `faims_one_cv.mzML` (-45, 112 spectra), `faims_two_cv.mzML` (-60 and -45, 56 each) and the FFC_1 input (NaN, 112, ranges present), their CV sets, and `faims_facts` facts 1a (-60 and -45, 2 each) and 1b (56 each), all with `byMSLevel(1)` throwing on voltage groups |
| `tests/im_data_converter.rs` native tests | 4 | group keys; message levels, texts and indices; returned skipped spectra and chromatograms; settings and run id copies; the empty experiment; refusal of settings beyond the resource limits with the input unchanged; the log-stream emulation's eviction rule; the ceiling |

The oracle is product-sdk (Debug, core `4fdec46`), the accepted development-time
oracle (D7). `git diff 4fdec46 bc9cc12` is empty for `IMDataConverter.{h,cpp}`,
`IMDataConverter_test.cpp`, `FAIMSHelper.cpp`, `MSExperiment.{h,cpp}`,
`LogStream.{h,cpp}`, `MzMLHandler.cpp`, `ExperimentalSettings.cpp`,
`SpectrumRangeManager.h` and `MSSpectrum.cpp`, and the installed
`IMDataConverter.h` and `FAIMSHelper.h` hash equal to the pin. The driver ran
twice with byte-identical output. No Debug-only precondition lies on this path.
Every comparison is exact (integers, strings and float bits), so it holds on
every platform.

## Class-test section accounting

`IMDataConverter_test.cpp` has 10 sections.

| Section | Rust test | Status |
| --- | --- | --- |
| `IMDataConverter()` | `constructor_and_destructor_sections` | passes |
| `~IMDataConverter()` | `constructor_and_destructor_sections` | passes |
| `splitByFAIMSCV(PeakMap& exp)` | `split_by_faims_cv_section` (feature `mzml`) | passes unchanged, through the Rust mzML reader |
| `splitByFAIMSCV assigns MS2 without explicit CV to last seen FAIMS CV` | `split_by_faims_cv_assigns_ms2_without_explicit_cv_to_last_seen_faims_cv_section` | passes unchanged |
| `splitByFAIMSCV returns single-element group for non-FAIMS dataset` | `split_by_faims_cv_returns_single_element_group_for_non_faims_dataset_section` | passes unchanged |
| `setIMUnit(DataArrays::FloatDataArray&, const DriftTimeUnit)` | none | not ported |
| `getIMUnit(const DataArrays::FloatDataArray&, DriftTimeUnit&)` | none | `NOT_TESTABLE` in the source (tested with `setIMUnit`); not ported |
| `reshapeIMFrameToMany(MSSpectrum im_frame)` | none | not ported |
| `splitExperimentByIonMobility(...)` | none | not ported (needs `SpectraMerger`) |
| `reshapeIMFrameToSingle(const MSExperiment& in)` | none | `NOT_TESTABLE` in the source (tested with `reshapeIMFrameToMany`); not ported |

Five sections ported, five recorded as not ported, none unaccounted. The ledger
entry for `IMDataConverter.h` stays `partial`.

## Deferrals

- **The rest of the header** (`reshapeIMFrameToMany`, `reshapeIMFrameToSingle`,
  `splitExperimentByIonMobility`, `setIMUnit`, `getIMUnit`) waits for its
  consumers (IonMobilityBinning, FeatureFinderMetabo, FeatureFinderMultiplex) and
  for `SpectraMerger` and `IMDataArrayUtils`.
- **Logging.** When kernel modules are wired to `LogStream`, the messages can be
  forwarded there; until then callers relay `FaimsSplit::messages`.
- **The FeatureFinderCentroided FAIMS closure** (B11, D5) decides whether the
  tool reproduces the source's exit 8 or runs the groups; this module gives
  groups with correct ranges either way.

## Source defects observed

C++ issue candidates, each with executed evidence in `oracle_cases.tsv`:

- **Split groups carry no ranges**, so FeatureFinderCentroided cannot process any
  FAIMS input (all 35 voltage groups throw at `byMSLevel(1)`; C2 facts 1a, 1b).
- **Chromatograms of FAIMS input are destroyed without a message**
  (`faims_test_data`, `settings_and_chromatograms_faims`).
- **NaN voltages corrupt the grouping** (`nan_first`, `nan_middle`, `nan_last`).
- **Message wording:** "Not FAIMS compensation voltages found" for "No FAIMS
  compensation voltages found".

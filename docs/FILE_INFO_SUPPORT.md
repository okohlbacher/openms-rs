# FileInfo library support

Native coverage of `FORMAT/FileInfo.h` and `FORMAT/FileInfo.cpp` at Core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`: the peak-file branch for DTA, DTA2D
and mzML and the featureXML branch, each with `-m`, `-p` and `-s`, and the text
and TSV reports. Work package A4-FILEINFO-CORE of the early TOPP bundle. The
FileInfo tool (`topp/src/FileInfo.cpp`) is package A5; `-i`, `-d` and `-c` are
A6, which has landed — see
[FILE_INFO_CHECKS_SUPPORT](FILE_INFO_CHECKS_SUPPORT.md) and
`src/format/file_info/checks.rs`; the consensusXML, identification and FASTA
branches A7; `-v`, mzXML, mzData and trafoXML A8. The ledger row for
`FileInfo.h` stays partial until those land.

| Artifact | Path |
| --- | --- |
| Options and result model | `src/format/file_info/model.rs` |
| Run, section order, shared rendering | `src/format/file_info/report.rs` |
| Peak-file branch | `src/format/file_info/peaks.rs` |
| featureXML branch | `src/format/file_info/features.rs` (feature `featurexml`) |
| Numeric text formatting (A2) | `src/format/file_info/text_format.rs`, whose module documentation is its support document; manifest `tests/data/file_info_text_format_provenance.json` |
| Tests | `tests/file_info.rs` |
| Fixtures | `tests/data/file_info/` (`inputs/`, `expected/`, `retained/`) |
| Manifest | `tests/data/file_info_provenance.json` |
| Oracle drivers | `../oracle/topp-early-bundle/` (C1) and `../oracle/file-info-core/` (this package), outside this repository |

The four Rust files cover the one header together with A2's `text_format.rs`.
They add one module edge, `format -> math`, for `SummaryStatistics`; it closes
no cycle (`tools/check_module_cycles.py`). The computation is serial, as the
source is.

## API mapping

Every public member of `FileInfo.h`, and the file-local helpers of
`FileInfo.cpp` the reports depend on.

| C++ member | Rust |
| --- | --- |
| `class FileInfo` | `format::file_info::report::FileInfo`, a zero-sized unit struct: the source class holds no state |
| `FileInfo()`, `~FileInfo()` | `FileInfo::new`, derived `Default`; no drop glue |
| `struct Range { bool present; double min, max; }` | `model::Range { min, max }`; `present == false` is `None` in the containing `RangeSet` |
| `struct RangeSet { rt, mz, mobility, intensity, has_mobility }` | `model::RangeSet` with `Option<Range>` dimensions |
| `struct Ranges { combined, spectra_overall, per_ms_level, chromatograms, is_experiment }` | `model::Ranges`; `std::map<UInt, RangeSet>` is `BTreeMap<u32, RangeSet>` |
| `struct FileMeta { file_name, file_type, file_type_name }` | `model::FileMeta`; `FileTypes::Type` is `format::FileType` |
| `struct ExperimentMeta` and its `struct Contact` | `model::ExperimentMeta`, `model::Contact`; never filled, as in the source |
| `struct ProcessingStep` | `model::ProcessingStep` |
| `struct NamedStats { title, stats }` | `model::NamedStats`; `Math::SummaryStatistics<std::vector<double>>` is `math::statistic_functions::SummaryStatistics`; never filled, as in the source |
| `struct PeakInfo` (every field) | `model::PeakInfo`; `Int` keys `i32`, `UInt64` counts `u64`, `std::map` keys `BTreeMap` |
| `PeakInfo::activationMethodsFlat() const` | `PeakInfo::activation_methods_flat` |
| `struct FeatureInfo` and its `struct MapColumn` | `model::FeatureInfo`, `model::MapColumn`; the consensus-only fields stay empty until A7 |
| `struct IdentInfo` | `model::IdentInfo`; declared, filled by A7 |
| `struct FastaInfo` | `model::FastaInfo`; `std::map<char, UInt64>` is `BTreeMap<u8, u64>`; declared, filled by A7 |
| `struct MzTabInfo` | `model::MzTabInfo`; the member `type` is `kind`; filled by the mzTab branch, which no package ports yet |
| `struct ValidationInfo` | `model::ValidationInfo`, `supported` defaulting to `true`; the four index fields are filled by A6 (`-i`) from `src/format/indexed_mzml.rs`, which departs from the source's decoder at two measured boundaries (native differences 11 and 12 of [FILE_INFO_CHECKS_SUPPORT](FILE_INFO_CHECKS_SUPPORT.md)), and the rest awaits A8 (`-v`); `schema_version` and `detail` are never written by the source run |
| `struct CorruptionInfo`, `struct DetailInfo` | `model::CorruptionInfo`, `model::DetailInfo`; never filled, as in the source (the `-c` and `-d` output goes only into the text) |
| `struct Result` | `model::FileInfoResult` (`Result` is the crate's error alias); native field `warnings` |
| `struct Options` (all eight members) | `model::Options`; `ProgressLogger::LogType` is `concept::progress_logger::ProgressLogType` |
| `Result run(const std::string&, const Options&)` | `FileInfo::run(&self, impl AsRef<Path>, &Options) -> Result<FileInfoResult>` |
| `Result run(const std::string&)` | `FileInfo::run_default` |
| `Result runAll(const std::string&)` | `FileInfo::run_all` |
| `static std::string toText(const Result&, const Options&)` | `FileInfo::to_text_with_options`, returning `&str` |
| `static std::string toText(const Result&)` | `FileInfo::to_text` |
| `static std::string toTSV(const Result&, const Options&)` | `FileInfo::to_tsv_with_options` |
| `static std::string toTSV(const Result&)` | `FileInfo::to_tsv` |
| private `report_` | private: `FileInfo::run` writes the header and dispatches; `peaks::report` and `features::report` write their branch and its `-m`, `-p` and `-s` sections |
| helper `printChargeDistribution` | `report::write_charge_distribution` (crate-private) |
| helper `operator<<(ostream&, const SummaryStatistics&)` | `report::write_summary_text` |
| helper `writeRangesHumanReadable_` (map and `MSExperiment`) | `report::write_ranges_text`, called per block by `peaks` and `features` |
| helper `writeRangesMachineReadable_` (map and `MSExperiment`) | `report::write_ranges_tsv`, with the block's key prefix |
| helper `writeSummaryStatisticsMachineReadable_` | `features::write_summary_tsv` (only the featureXML branch writes statistics TSV) |
| helpers `extractRangeSet_`, `extractRanges_`, `extractRangesExp_` | `report::range_set`; `features::feature_map_ranges`; the range part of `peaks::Summary::compute` |
| helper `struct IdData` | `identifications::IdData` (A7) |
| native only | `FileInfo::MAX_STATISTICS_VALUES`; `peaks::PEAK_TYPE_ESTIMATION_MIN_PEAKS`; `FileInfoResult::warnings` |

## Preserved source conventions

- **Report order** (`FileInfo.cpp:662-2445`): header, `-v`, `-i`, content, `-m`,
  `-p`, `-s`, then `"\n\n"`. The TSV has no section titles and no trailing
  newlines, and many text lines have no TSV twin (peak statistics, `-p`'s
  first-spectrum note, `-m` sample and instrument lists and contacts).
- **Unknown type**: forced type first, else `FileHandler::getType`; an unknown
  type returns a result with only `meta` filled (`file_type_name` `unknown`) and
  empty reports, before any flag is looked at. A directory whose name gives no
  type is unknown too: the source content check reads it as an empty file
  (`FileHandler.cpp:398-411`), where the native `FileHandler::get_type` returns
  the I/O error, so `FileInfo::run` maps that one case
  (`a_directory_without_a_recognised_name_is_an_unknown_type`; the product-SDK
  tool reports `Could not determine input file type!` for it).
- **Forced type** selects the branch and is the loader's only allowed type; the
  loader still detects the type by name and then content, so a featureXML map
  named `.tmp` loads when featureXML is forced (`FileInfo_test.cpp:120-145`).
- **Loading**: peak files with default `PeakFileOptions`; featureXML with
  convex hulls and subordinates off (`FileInfo.cpp:1080-1081`).
- **SRM spectra of mzML files** become chromatograms before anything is
  counted. The source mzML load step ends with
  `ChromatogramTools().convertSpectraToChromatograms<PeakMap>(exp, true)`
  (`FileHandler.cpp:906-911`), which the native loader leaves out, so the
  peak-file branch calls `ChromatogramTools::convert_spectra_to_chromatograms`
  with `remove_spectra` set and `force_conversion` unset on mzML input. Each
  SRM spectrum with one precursor and at least one peak adds one point per peak
  to the chromatogram of its (precursor, product) pair, appended after the
  stored chromatograms; every SRM spectrum is then removed, including those the
  conversion skips (no precursor, two precursors, no peaks). The spectrum
  counts, ranges, activation methods, charges, data-array names, `-p` (the first
  remaining spectrum) and `-s` all see the converted experiment
  (`a4_srm_spectra_become_chromatograms_all_flags`,
  `a4_srm_spectra_among_ordinary_spectra`).
- **Ranges**: two-decimal `StringUtils::number` for every bound and one decimal
  for the span in minutes, through `text_format::fixed_truncated` (A2's note:
  `fixed` would refuse what `number` cuts). An absent dimension prints
  `<none>`, never zero. Experiment blocks print an `ion mobility` line in the
  combined, overall-spectrum and per-level blocks, never in the chromatogram
  block; the TSV prints `ion-mobility` lines only when that range is present.
  Chromatogram m/z is each chromatogram's product m/z, even without points.
  Feature-map ranges extend positions and intensities first, then hull boxes,
  with the keep-first rule on equal endpoints.
- **Peak type per MS level**: the stored type of the first spectrum of the
  level (`getType(false)`, which also honours a peak-picking processing step),
  and `PeakTypeEstimator` on the first spectrum of the level with more than
  ten peaks; a level without one prints `Unknown` (the source map's
  `operator[]` default).
- **Counts**: activation methods over every precursor, keyed by `(level,
  enum)`; precursor charges from the first precursor per spectrum; data-array
  names from float, integer and string arrays together, counted per array,
  bytewise-ordered, padded by byte length; chromatogram types in enum order;
  total peaks include chromatogram points.
- **FAIMS**: `FAIMSHelper::getCompensationVoltages` on every peak file; a
  non-empty set prints `IM (FAIMS_CV): [..]` with `StringUtils::toStr` text.
- **Numbers**: the total ion current is the `float` intensities summed as
  `double` in file order; resolutions, the total ion current and statistics go
  through the default stream at the tracked precision (`ostream_g`); integers
  print in decimal.
- **Stream precision** is tracked per report: 6 at the start, set to
  `writtenDigits<float>()` before each statistics block, and kept afterwards.
  Both values are 6 on these branches, so the persistence is not observable
  yet; the consensus branch (A7) switches to 15.
- **Statistics**: `SummaryStatistics` (sort; mean; variance with `n - 1`, 0 for
  one value; all zero when empty) over `float` or `int` values promoted to
  `double`. featureXML: intensity, width (`Feature FWHM in RT dimension`),
  overall, RT and m/z quality, text and TSV. Peak files: MS1 intensities, then
  one block per data-array name over the float and integer arrays of that name
  (a string array's name gives an all-zero block); text only. As in the source,
  one name's values are collected, summarised and released before the next
  name's, so at most one block is held.
- **`-p`**: the map's processing, or the first spectrum's with the note
  `Note: The data is taken from the first spectrum!`; an empty list prints
  `No information about data processing available!`; actions in enum order;
  an unset completion time prints `0000-00-00 00:00:00`. The structured steps
  are collected only when `-p` is set.
- **`-m`**: featureXML prints `Document ID` and the TSV `meta: document ID`;
  peak files the document identifier, `DateTime::get` date, sample, instrument
  and contacts.
- **Label spacing**: the source at the pin writes one space after
  `intensity:` (`FileInfo.cpp:140`, `:189`). The retained `FileInfo_3_output.txt`
  (and `FileInfo_7`) has six; only FuzzyDiff's whitespace rule accepts it, and
  the port follows the source (`retained_file_info_3_passes_fuzzy_diff`
  asserts both spellings).
- **`-d` and `-c` on featureXML** have no effect in the source (read only in the
  peak-file branch); the port produces the default report
  (`detailed_listing_and_corrupt_check_do_not_change_a_featurexml_report`).

## Native differences

1. **Refusals instead of partial support.** `-v` for every type, the
   pepXML, mzTab, trafoXML and PQP branches, and peak files of mzXML, mzData, MGF, MS2,
   sqMass, XMass, MSP, Thermo RAW and Bruker TDF return `Error::Unsupported`
   naming the branch. The source reports them; it loads RAW and TDF when built
   with its default `WITH_THERMO_RAW` and `WITH_OPENTIMS` options (the
   product-SDK oracle is built without both). MGF and MS2 have native loaders
   but no FileInfo oracle yet. The refusal comes once the type is known and
   before the file is loaded. The type is known without file access when it is
   forced or recognised from the name (every refusal test uses such names);
   otherwise type detection reads the start of the file first, so a missing
   file with an unrecognised name gives `Error::Io` even with `-v` set
   (`refusals_follow_content_detection_of_an_unrecognised_name`).
2. **Errors of unloadable types.** The source writes the header and then throws
   from `FileHandler::loadExperiment`; the port returns no result:
   `Error::Parse` for a type no experiment loader of the source handles (source
   `ParseError`), `Error::InvalidValue` for imzML (source `InvalidFileType`)
   and for a detected type other than the forced one (source `ParseError`,
   mapped as `FileHandler` maps it).
3. **Absent ranges** are `Option<Range>`, not `present` plus zeros.
4. **Result shape.** `FileInfoResult` adds `warnings`: the source logs the
   FAIMS missing-voltage warning with `OPENMS_LOG_WARN`, twice, because it asks
   for the voltages twice (`FileInfo.cpp:1673`, `:1742`); the library computes
   them once and returns the warning once.
   `to_text` and `to_tsv` borrow the cached reports instead of copying them.
   The file name must be UTF-8 (`Error::InvalidValue`).
5. **Loader strictness is inherited.** The native readers refuse inputs the
   source loads; see *Known reader gaps*. The library runs every load with the
   strict default reader options; `Options` has no field for the
   source-compatibility load options of D10 yet (see *Deferrals*).
6. **Non-finite and signed-zero values.** The readers refuse NaN and infinite
   coordinates and intensities, and a FAIMS spectrum with a NaN voltage is
   refused; the source would print `nan` or break its set order.
   `SummaryStatistics` used to sort with `f64::total_cmp`, so `-0.0` sorted
   before `0.0` where `std::sort` calls the two equivalent and libstdc++ leaves
   them in input order, and a sample holding both zeros could print `-0` where
   the source prints `0` as minimum or maximum. That was **measured in both
   orders and pinned** in wave 8 as native difference 6 of
   [A7](FILE_INFO_A7_SUPPORT.md) (oracle cases `c_nan_one_s` and
   `c_zero_swapped_s`), and it is **closed** since the shared-math wave of
   2026-09-19: under lead decision D16 the crate's private `sort_ascending` is
   the Release build's own `std::sort`
   ([`crate::math::source_sort`](STATISTIC_FUNCTIONS_SUPPORT.md#nan-policy)), so
   the port keeps the input order the source keeps and reproduces **both**
   members of the measured pair. The test is renamed
   `consensus_a_signed_zero_sample_keeps_the_release_builds_order` and now
   asserts equality rather than a divergence. Kernel hull boxes merge with `f64::min` and `f64::max` (open B1
   follow-up); FileInfo loads no hulls. All of that is about non-finite
   *inputs*.

   A statistics block can also compute a non-finite value from finite input, and
   that is a separate matter. The consensusXML `-s` relative intensity error
   divides (`FileInfo.cpp:2310`) and inverts every ratio below 1
   (`:2312-2315`), so a sub-feature of intensity zero makes the sample
   `{1, +inf}`, whose mean is `+inf` and whose variance is a NaN; one of
   intensity `-0.0` next to one of `0.0` contributes `(-inf) + (+inf)`, which
   puts a NaN into the per-consensus-feature *sample* itself. The port
   reproduces those values bit for bit **and spells them as the reference build
   spells them**. That was native difference 5 of
   [the A7 document](FILE_INFO_A7_SUPPORT.md) — the port wrote every NaN `nan`
   where glibc writes a sign-bit NaN `-nan` — and it is **closed**. Every
   compared A7 oracle report is now byte-identical to the Release build's.

   Closing it took two steps, in this order, because the spelling is only
   honest once the value is:

   - *The value.* "Bit for bit" became true of every host, and not only of
     x86_64, in the shared-math wave of 2026-09-19: every operation of
     `src/math/statistic_functions.rs` that can *generate* a NaN is built on
     `crate::math::x86_64`, so `inf - inf` is the Release build's
     `0xfff8000000000000` on an arm64 host too. One generator was outside that,
     `src/format/file_info/consensus.rs`'s `it_aad += it_ratio` — the
     `(-inf) + (+inf)` above — and it now goes through `math::x86_64::add`. A
     sweep of `src/format/file_info/` found no second instance: the consensusXML
     reader's `map.validate()` refuses a non-finite coordinate, and the
     featureXML branch's TIC loop runs after `RangeBase::extend_value` has
     refused one, so nothing else there can reach an invalid operation.
   - *The spelling.* `text_format::nonfinite` writes `-nan` for a NaN whose
     sign bit is set. All three of its callers reach C `printf`; the
     `StringUtils::toStr` path does not come through it and still writes `NaN`
     for either sign, which is `NumericFormatting.h:29`. The measurement is
     `../oracle/a2-textfmt-linux`, which re-ran A2's driver against the Linux
     Release install: of 1018 rows exactly one differs from the macOS capture,
     and it differs in the five `printf` columns and not in the `toStr` column.

   *What is measured and what is generalised.* That corpus holds exactly three
   NaN bit patterns and one sign-bit NaN row, so it pins `number(-NaN, n)` at
   `n` in `{0, 1, 2}` and `ostream(-NaN, p)` at `p` in `{6, 15}`, and no
   sign-bit NaN `float` at all. That gap is **closed by measurement**, not by
   the argument that glibc writes the sign before `__printf_fp` dispatches on
   the class: the companion sweep `../oracle/a2-textfmt-nan-sweep` covers six
   NaN shapes by both signs by `double` and `float` over 22 digit counts and 16
   precisions — **912 sign-bearing rows, none printing a sign that disagrees
   with the argument's sign bit**, and 24 `toStr` rows all printing `NaN`. It is
   registered in `SOURCE_PROVENANCE.json` with a sha256 and a tier-1 label, and
   asserted row for row by `the_negative_nan_sweep_is_reproduced_row_for_row`. A macOS C++ build
   writes `nan` for the same bits, so a macOS comparison must not count the
   difference as a port defect — the same caveat the `%g` tie class carries.

   `SummaryStatistics::new` no longer refuses any NaN. It used to summarise only
   the two shapes in which the permutation `std::sort` leaves behind cannot be
   observed — a one-value sample and an all-NaN sample — and refuse a NaN next
   to a number. Under lead decision D16 the crate reproduces the permutation
   itself, so every shape is summarised and the only refusal left is the
   out-of-bounds one `crate::math::source_sort` raises; see section 5.2 of the
   A7 document for the measurement.
7. **Bounded work.** A statistics block holds at most
   `FileInfo::MAX_STATISTICS_VALUES` (2^27) values, checked before collection
   and allocated fallibly; the kernel range managers, the FAIMS scan, the
   stored-type query and `PeakTypeEstimator` apply their own ceilings (the
   estimator refuses a spectrum above one million peaks, where the source
   classifies it). An MS level above `i32::MAX` cannot become a `PeakInfo` key
   and is refused; the source casts. The SRM conversion of mzML input runs with
   the default `ChromatogramConversionLimits` and refuses, with
   `Error::InvalidValue`, an input above them; a file without SRM spectra costs
   it one step per spectrum and never reaches them.
8. **`log_type`** is accepted and has no effect: the native loaders on this
   path report no progress.
9. **Default-stream ties.** `ostream_g` follows the C standard and glibc for
   the class of integer-valued ties below `1e15` where Apple libc keeps trailing
   zeros (A2's platform note, `text_format.rs`). None of the 34 oracle report
   pairs here holds such a value; a macOS oracle text that differs from the
   port only there is not a port defect.

## Known reader gaps

The C++ loader accepts these inputs and the strict native mzML reader refuses
them. Each gap lies outside this package; the affected oracle tests are
`#[ignore]`d with the reason, `reader_gaps_behind_the_ignored_cases_are_still_present`
fails once a gap closes, and a derived input that avoids the gap is compared
instead where the C++ report on it is identical apart from the file name.

| Input | Refusal | C++ | Ignored test | Derived replacement |
| --- | --- | --- | --- | --- |
| `FileInfo_9_input.mzML` | `duplicate userParam name name` (repeated spectrum-level `name`) | `MetaInfo` overwrites | `c1_file_info_9_mzml_mps`, `a4_file_info_9_default_flags` | `FileInfo_9_strict_reader.mzML` |
| same | `primary-array processing has no independent native owner` (`dataProcessingRef` on m/z and intensity arrays) | kept on the array | same | same |
| same | `canonical auxiliary array binary type` (64-bit float `charge array`) | loaded | same | same (re-encoded as 32-bit integers) |
| `empty.mzML` (C1 derived) | `unresolved dataProcessingRef` (dangling `defaultDataProcessingRef`), with the strict default | ignored | `c1_empty_mzml_mps` | `empty_resolved_ref.mzML`, or `Options::source_dangling_references` (A5, P2, D10), which the FileInfo tool sets |
| `FileInfo_12_input.mzML` | `canonical auxiliary array binary type` | loaded | `a4_indexed_file_info_12_all_flags` | none |
| `MzMLFile_1.mzML` | loads, but the selected-ion drift time is not copied onto the MS2 spectrum (A3 request 5) | copied (`MzMLHandler.cpp:1871-1875`), ranges end at 8.10 | `a4_mzml_file_1_all_flags` | `MzMLFile_1_no_selected_ion_drift.mzML` |

## Checked boundaries and evidence

Tier 1, executed differential (product-SDK FileInfo, Debug, core `4fdec46`).
The installed `FileInfo.h` and `FileHandler.h` hash equal to the pin, and C1's
manifest records that on the traced FileInfo call path only
`FEATUREFINDER/MassTraceDetection.h` and `IONMOBILITY/IMTypesExperiment.cpp`
change between the oracle build and the pin, with numeric text output
identical (commit `74526a8`). No case exits non-zero, so no Debug-only
precondition is involved. The TOPP tool writes `toText` and `toTSV` of
the library result, so its `-out` and `-out_tsv` files are the library
reports. Text and TSV are compared byte for byte, only the `File name` line and
the `general: file name` TSV line normalised; 34 report pairs in all:

- C1 (`../oracle/topp-early-bundle`, run1): FileInfo_1 (`-in_type dta`),
  FileInfo_2, FileInfo_3 (`-m -s -p`) and FileInfo_9 (`-m -p -s`) with
  `-out_tsv`; the empty featureXML and empty mzML with `-m -p -s`; the retained
  FeatureFinderCentroided_1 output with `-m -p -s`.
- This package (`../oracle/file-info-core`, reproduced twice): the class-test
  featureXML with all flags, default flags and `-s` only, and forced under a
  `.tmp` name; the class-test minimal mzML; FileInfo_3 and FileInfo_9 with
  default flags; FileInfo_12; `FAIMS_test_data.mzML`,
  `FAIMS_CV-60C_V-45_Interleaved.mzML` and `IM_FAIMS_test.mzML` (FAIMS voltages
  and mobility ranges); `FeatureFinderCentroided_1_input.mzML`; 18 SRM
  chromatograms (`mzml_numpress_source_original.mzML`); `MzMLFile_1.mzML`
  (two processing steps, three activation methods, drift times, two contacts);
  `precursor_purity_input.mzML`; a profile and a centroid DTA; the three DTA2D
  header variants; the three derived inputs above; and two generated SRM-spectra
  mzML files (`srm_spectra.mzML` with all flags; `srm_spectra_mixed.mzML` with
  all and with default flags), whose generator and rules are in `oracle.py` and
  its manifest.

Tier 1, retained upstream output: TOPP_FileInfo_1, _2, _3 and _9
(`topp/CMakeLists.txt:881-883`, `:884-886`, `:887-889`, `:902-904` at test-data
`0cb15f2`) through C3's FuzzyStringComparator with the pinned `FuzzyDiff.ini`
(ratio 1.01, absdiff 0.01) and the registered whitelist `File name`. FileInfo_9
runs on the derived input.

Tier 4: every refusal and its message, before file access for forced and
recognised names and after content detection otherwise; unknown type, a
directory among them;
unloadable, imzML and forced-type mismatch errors; missing files (`Error::Io`);
truncated featureXML and mzML (`Error::Parse`; C++ exits 3); a non-UTF-8 name;
the file name printed as given; the `Options` and `Result` defaults; the
structured `PeakInfo`, `FeatureInfo`, ranges and processing steps checked
against the executed reports.

Feature lines: `tests/file_info.rs` gates mzML cases on `mzml` and featureXML
cases on `featurexml`; it runs with `--no-default-features` (22 tests), with
`mzml` (47 tests, 5 ignored), with `featurexml` (38 tests), with
`mzml,featurexml` (63 tests, 5 ignored) and with `--all-features` (64 tests,
5 ignored).

## Class-test section accounting

`src/tests/class_tests/openms/source/FileInfo_test.cpp` at `bc9cc12`, nine
sections.

| Section | Status |
| --- | --- |
| `FileInfo()` | ported: `class_test_constructor_and_destructor` |
| `~FileInfo()` | ported: same test |
| `run ... - featureXML` | ported: `class_test_run_featurexml` |
| `run ... - consensusXML` | pending (A7); `class_test_run_consensusxml_is_pending` asserts the explicit refusal meanwhile |
| `run ... - mzML peaks` | ported: `class_test_run_mzml_peaks` |
| `run ... - FASTA` | pending (A7); `class_test_run_fasta_is_pending` asserts the explicit refusal meanwhile |
| `toText` and `toTSV` | ported: `class_test_to_text_and_to_tsv` |
| `Options gating` | ported: `class_test_options_gating` |
| `forced type selects the parse branch for an unrecognized extension` | ported: `class_test_forced_type_selects_the_parse_branch` (temporary directory instead of a fixed name) |

## Deferrals

- `-v`, mzXML, mzData, trafoXML (A8); pepXML, mzTab, PQP, sqMass, XMass, MSP,
  MGF and MS2 have no package yet. `-i`, `-d` and `-c` are no longer deferred:
  A6 ported them ([FILE_INFO_CHECKS_SUPPORT](FILE_INFO_CHECKS_SUPPORT.md)), and
  its two open items are the mzML reader and kernel refusals that keep two `-c`
  lines out of reach. The consensusXML, idXML/mzIdentML and FASTA branches are
  no longer deferred either: A7 ported them
  ([FILE_INFO_A7_SUPPORT](FILE_INFO_A7_SUPPORT.md)), refusing only the three
  places the source's behaviour is an out-of-bounds `std::vector` access.
- The tool wrapper, exit codes and output routing are A5's.
- Source-compatibility load options (D10): closed for dangling header
  references, open for the rest. A5 added the native
  `Options::source_dangling_references`, which defaults to `false`, so every
  library default stays strict, and with the `mzml` feature reaches the mzML
  reader through `FileHandler::load_experiment_with_read_options` (P2's
  `mzml::ReadOptions::source_dangling_references`). The FileInfo tool sets it,
  which is what makes the C1 `empty.mzML` report reachable there
  (`tests/topp_file_info.rs`, `c1_empty_mzml_with_a_dangling_reference`); the
  library case `c1_empty_mzml_mps` stays `#[ignore]`d, because it runs with the
  strict default. The other three ignored cases need reader work no option
  covers yet (see *Known reader gaps*).

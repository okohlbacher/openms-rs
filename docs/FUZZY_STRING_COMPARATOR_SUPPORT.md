# FuzzyStringComparator, FuzzyDiff and decoded comparison

The upstream TOPP suite judges almost every tool output with
`FuzzyDiff -test -ini FuzzyDiff.ini [-whitelist ...]`, which runs
`OpenMS::FuzzyStringComparator`. This group ports that comparator, first as
shared test support and, since the `FuzzyDiff` tool was ported, as the library
module `concept::fuzzy_string_comparator`. The move changed no behaviour: the
module is the test-support code, the test support re-exports it, and every
verdict and log byte of the executed corpus below is reproduced as before. The
test support keeps an emulation of the `FuzzyDiff` tool contract for tests
built without `paramxml`, and a decoded-content comparator for XML outputs
(decision D6). The tool itself is documented in
[TOPP_FUZZY_DIFF_SUPPORT](TOPP_FUZZY_DIFF_SUPPORT.md). The source header is
outside the registered SDK (it is installed by `src/testframework`), so it has
no ledger key.

| Rust file | Covers |
|---|---|
| `src/concept/fuzzy_string_comparator.rs` | `FuzzyStringComparator.h/.cpp` (core bc9cc12), and the two input transforms of `FuzzyDiff::main_`, `sorted_lines` and `parse_matched_whitelist`; `std` only, built without any feature |
| `tests/support/fuzzy_string_comparator.rs` | re-exports the module; the `FuzzyDiff` contract emulation (topp 174b576 `src/FuzzyDiff.cpp` with the TOPPBase behaviour it relies on) and a bounded ParamXML reader for `FuzzyDiff.ini` |
| `tests/support/decoded_compare.rs` | Field-by-field comparison of `FeatureMap` and `MSExperiment` values with the comparator's number rule |
| `tests/fuzzy_string_comparator.rs` | Class-test port, executed C++ differential, `FuzzyDiff` exit-code parity, retained-vs-current pairs, decoded comparator tests |
| `tests/data/fuzzy_string_comparator/` | Pinned `FuzzyDiff.ini`, retained upstream fixtures, oracle outputs and synthetic inputs |
| `tests/data/fuzzy_string_comparator_provenance.json` | Hashes, pins and evidence tiers |

Source pins: core `bc9cc12` (`src/testframework/include/OpenMS/CONCEPT/FuzzyStringComparator.h`
`43e9d1be…`, `src/testframework/source/CONCEPT/FuzzyStringComparator.cpp` `9b7c5bd5…`,
`src/tests/class_tests/openms/source/FuzzyStringComparator_test.cpp` `737e1fd4…`), topp
`174b576` (`src/FuzzyDiff.cpp` `3809921f…`), cli `c19e494` (`source/APPLICATIONS/TOPPBase.cpp`
`326b96f8…`) and test-data `0cb15f2` (`topp/CMakeLists.txt` `0bf7d0ef…`, `topp/FuzzyDiff.ini`
`238d0598…`). The comparator source is byte-identical at the oracle's core revision 4fdec46.

## Using it

```text
#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;

// A registered text comparison, e.g. TOPP_FileInfo_3_out1 (CMakeLists.txt 887-889):
let settings = fuzzy::FuzzyDiffSettings::upstream()?.with_whitelist(&["File name"]);
settings.compare_bytes(&actual, &expected)?;          // Err carries the comparator log

// The whole tool contract on files, returning the C++ exit code:
let outcome = fuzzy::fuzzy_diff(&in1, &in2, &settings); // outcome.exit, outcome.log
```

```text
#[path = "support/fuzzy_string_comparator.rs"]
mod fuzzy;
#[path = "support/decoded_compare.rs"]
mod decoded; // uses crate::fuzzy; never load the comparator a second time

// An XML output, e.g. TOPP_FeatureFinderCentroided_1_out1 (-whitelist "id="):
let tolerance = decoded::Tolerance::from_settings(&fuzzy::FuzzyDiffSettings::upstream()?);
let options = decoded::DecodedOptions::new(tolerance).ignoring_unique_ids();
decoded::compare_feature_maps(&actual, &expected, &options)?; // Err: path + detail
```

`FuzzyDiffSettings::upstream()` loads the pinned INI: ratio 1.01, absdiff 0.01, whitelist
`<?xml-stylesheet`, verbose 1. A registration's `-whitelist` **replaces** that list, as
TOPPBase does for list parameters; executed on the oracle (`cli_whitelist_replaces_ini`).
Never widen a tolerance to make a comparison pass. A test that needs different values
sets them explicitly and names the reason next to the call.

The upstream tolerance is loose, at 1 % or 0.01 absolute. The retained
FeatureFinderCentroided_1 expectation already differs from current C++ by about 1e-10
relative (`Tolerance::exact()` rejects the pair). Tight scientific comparisons should
therefore use executed oracle outputs, not only retained expectations.

## API mapping

### `FuzzyStringComparator`

| Source member | Rust |
|---|---|
| `FuzzyStringComparator()` | `FuzzyStringComparator::new`, `Default` |
| `virtual ~FuzzyStringComparator()` | implicit `Drop` |
| copy constructor, `operator=` (declared, not implemented) | not provided; the type does not implement `Clone` |
| `getAcceptableRelative` / `setAcceptableRelative` | `acceptable_relative` / `set_acceptable_relative` (reciprocal below one; `normalise_relative`) |
| `getAcceptableAbsolute` / `setAcceptableAbsolute` | `acceptable_absolute` / `set_acceptable_absolute` (negated below zero; `normalise_absolute`) |
| `getWhitelist() const` / `getWhitelist()` / `setWhitelist` | `whitelist` / `whitelist_mut` / `set_whitelist` |
| `getMatchedWhitelist` / `setMatchedWhitelist` | `matched_whitelist` / `set_matched_whitelist` |
| `getVerboseLevel` / `setVerboseLevel` | `verbose_level` / `set_verbose_level` |
| `getTabWidth` / `setTabWidth` | `tab_width` / `set_tab_width` |
| `getFirstColumn` / `setFirstColumn` | `first_column` / `set_first_column` |
| `getLogDestination` / `setLogDestination` | `log_destination` / `set_log_destination` with `LogDestination::{Stdout, Stderr, Buffer}`; `log`, `take_log`, `log_truncated` read the buffer |
| `compareStrings` | `compare_strings` (`&str`), `compare_bytes` (`&[u8]`, the bytes of a `std::string`) |
| `compareStreams` | `compare_streams<R1: BufRead, R2: BufRead>` |
| `compareFiles` | `compare_files(&Path, &Path)` |
| `compareLines_` (protected) | `compare_lines(&[u8], &[u8])`, public so single-line behaviour can be tested |
| `reportSuccess_`, `reportFailure_`, `writeWhitelistCases_` (protected) | private `report_success`, `report_failure`, `whitelist_cases_text`; output byte-identical to C++ (see evidence) |
| `readNextLine_`, `openInputFileStream_` (protected) | private `LineSource::read_next_line`, `open_input_file` |
| `struct AbortComparison` | private `AbortComparison` returned through `Result` |
| `struct InputLine` (`setToString`, `updatePosition`, `seekGToSavedPosition`, `ok`) | private `InputLine`, emulating the stream's position, `eofbit`, `failbit` and saved position, including `tellg` setting `failbit` on a stream at EOF |
| `struct StreamElement_` (`reset`, `fillFromInputLine`) | private `StreamElement` |
| `struct PrefixInfo_` | private `PrefixInfo` |
| protected data members (`log_dest_` … `matched_whitelist_`) | private fields. `is_absdiff_small_` is not stored: it is only read on a fall-through that always continues. `use_prefix_` is kept and always false, since the source has no setter. |
| friends `Internal::ClassTest::testStringSimilar`, `isFileSimilar` | not ported; they back the C++ `TEST_STRING_SIMILAR`/`TEST_FILE_SIMILAR` macros, and Rust tests call the comparator directly |
| file-local `getLine`, `absolutePath`, `suffix`, `prefixOf`, `to_path` | `LineSource::get_line`, `absolute_display`, slicing |
| protected `input_1_name_`, `input_2_name_` (set by `compareFiles` only) | `input_names`; native `set_input_names`, with which the tool names its in-memory `-sort` texts |
| a stream exception escaping `compareFiles` (no source API) | native `InputFailure`, `input_failure` and `log_without_input_failure`: which comparison stopped on a read error or the input bound, and the log without its failure line; the log line itself is unchanged |
| file-local `tryParseNaN`, `parseFloat`, `extractDouble` | `extract_double` (public), following the `std::from_chars` contract |
| file-local `fromCharsFloat` (libc++ `strtod` fallback) | not ported; see native differences |
| number branch of `compareLines_` (556-689) | `compare_numbers` (public; also used by the decoded comparator) |
| `std::ostream << double` in both reports | `format_g` |

### `FuzzyDiff` (topp `src/FuzzyDiff.cpp`)

| Source | Rust |
|---|---|
| `registerOptionsAndFlags_`: `in1`, `in2`, `ratio` (1, min 1), `absdiff` (0, min 0), `whitelist` (`<?xml-stylesheet`), `matched_whitelist` (empty), `verbose` (2, 0-3), `tab_width` (8, min 1), `first_column` (1, min 0), `sort` | `FuzzyDiffSettings::registered_defaults` and its fields; `in1`/`in2` are the arguments of `fuzzy_diff` |
| `-ini <file>` | `FuzzyDiffSettings::load_ini` / `from_ini` (items of `FuzzyDiff:1:`; unknown items and unparsable values go to `ini_errors`) |
| `-whitelist`, `-matched_whitelist` on the command line | `with_whitelist`, `with_matched_whitelist` (replace the lists) |
| `main_`: matched-whitelist split (`IllegalArgument`), comparator setup, `-sort` temporary files, `compareFiles`, `EXECUTION_OK`/`PARSE_ERROR` | `fuzzy_diff`, `FuzzyDiffSettings::comparator`, the library's `parse_matched_whitelist` and `sorted_lines`, `FuzzyDiffOutcome` |
| TOPPBase exit codes used by the tool | `FuzzyDiffExit` (0, 1, 2, 4, 6, 7, 8, 10) |
| in-memory convenience without file checks | `FuzzyDiffSettings::compare_bytes` (native) |

Order of `fuzzy_diff`, as executed: an INI error or an out-of-range parameter, from the INI
or the command line, gives 6 during initialisation, before the file checks. Then come
`in1` and then `in2`: empty name 7, missing 1, unreadable 2, empty regular file 4. A
matched-whitelist entry that does not split into exactly two parts at `:` gives 8. The
comparison then gives 0 or 10.

## Preserved source conventions

Each of these decides verdicts or report text and is covered by an executed oracle case
(names in parentheses):

- `\n`, `\r\n` and a lone `\r` all end a line; empty and whitespace-only lines (C-locale
  `isspace`, including `\v` and `\f`) are skipped, and every read counts toward the line
  number (`cr_stream_*`, `whitespace_only_lines`).
- A whitespace run equals any whitespace run, but whitespace against the end of a line
  fails: `"a "` against `"a"` is a difference (`line_trailing_space`, `fd_trailing_space`).
- An exhausted line reads as a NUL letter, so `"ab"` against `"abc"` reports `different
  letters`. The "line … is shorter" message follows only at verbose level 3
  (`line_exhausted_v2`, `line_exhausted_v3`).
- A carriage return that starts a whitespace run inside a line is skipped when the other
  side holds a letter. It does not help against a number, where the number check comes
  first (`cr_line_*`).
- Whitelist: identical lines return before the whitelist is consulted and are never
  counted. The first entry contained in both lines is counted and skips the pair, and an
  empty entry matches every line. Matched pairs apply in either order and are not
  counted (`wl_*`, `mw_*`).
- Numbers: equal values (including `0 == -0` and equal infinities) pass. Two NaNs pass,
  and a NaN against anything else fails, as do opposite infinities and an infinity
  against a finite value. Otherwise a pair passes when `absdiff <= absdiff_allowed`, or
  when both are non-zero with the same sign and the ratio (inverted when below one) does
  not exceed `ratio_allowed` (`num_*`, `ct_nonfinite_*`).
- Number tokens: `inf`, `infinity` and `nan` are recognised anywhere a token starts, so
  `information` begins with `inf` and `Nancy` with NaN. A single leading `+` is skipped
  (`+-5` is -5, `++5` is not a number). An exponent without digits is not consumed.
  Overflow reads as letters and underflow as a number (`tok_*`).
- `compareStreams` resets only the success flag. Line numbers, maxima and whitelist counts
  carry over to the next comparison on the same instance (`reuse_*`).
- `compareFiles` refuses equal name strings ("That's cheating!") and reports both open
  failures as "Error opening first input file" (`file_*`).
- Reports: a failure report is 35 lines without whitelist hits, and all texts,
  tab-width-aware columns, raw letter bytes, `boolalpha` and `%g` numbers are reproduced
  byte for byte. The success report names the inputs as given, while the failure report
  prints absolute paths (all 137 cases).
- FuzzyDiff: `TOPP_FuzzyDiff_1` fails because both names are the same file ("cheating");
  `_2` fails on numbers; `_3` passes; `_4` exits 1 because `lorem_ipsum.featureXML` does
  not exist at the pin. The comments in `CMakeLists.txt` 142-144 have the reasons for 1
  and 2 swapped.

## Native differences

- **Number grammar.** `extract_double` follows the contract the source states for
  `std::from_chars`, general format: no hexadecimal floats, overflow rejected, underflow
  accepted. The Apple libc++ build that serves as the oracle parses through a `strtod`
  fallback. That fallback also accepts hexadecimal floats, so `"0x10"` equals `"16"` on
  macOS only. This single case (`tok_hex_vs_decimal`) is asserted as a known divergence.
  Whether libstdc++'s `std::from_chars` rejects underflow to zero, contrary to the source
  comment, is now executed: the Linux Release build rejects `1e-400` and
  accepts `1e-310` (`token_underflow`, `token_subnormal` in
  [TOPP_FUZZY_DIFF_SUPPORT](TOPP_FUZZY_DIFF_SUPPORT.md)), so the port, which
  accepts underflow as the source comment says, matches the libc++ build there
  and not the libstdc++ one.
- **Value of a non-number in the report.** The source resets it to NaN. The libc++
  fallback then overwrites it with `strtod`'s result for the letter (usually 0), while
  this port keeps NaN. The differential masks only this field, and only for elements
  whose `is_number` is false. NaN prints as `nan` whatever its sign, as on the oracle.
- **Read errors and bounds.** A read error, an input stream beyond `MAX_INPUT_BYTES`
  (1 GiB) or a larger file fails the comparison with a log line. The source treats a
  stream error as the end of input, which could let two broken inputs compare equal.
  Failure reports stop after `MAX_FAILURE_REPORTS` (10,000, verbose 3 only) and the
  buffered log after `MAX_LOG_BYTES` (256 MiB); verdicts are unaffected. INI input is
  limited to `MAX_INI_BYTES` (16 MiB), nesting to 64 and items to 100,000.
- **Log destination.** An enum replaces the `std::ostream*`. `Stdout` and `Stderr` go
  through `print!`/`eprint!`, so the test harness captures them, with invalid UTF-8
  replaced; `Buffer` keeps exact bytes.
- **Integers.** Line numbers are `i64` (C++ `int` overflow is undefined). A zero tab width
  gives a zero quotient in the column computation; in C++ the division is undefined, and
  arm64 also yields zero.
- **Maximum lines** are copied once per line instead of at every new maximum. The report
  is the same, without quadratic copying on long lines.
- **Whitelist entries** are UTF-8 `String`s; C++ uses byte strings. INI values are decoded
  as ISO-8859-1 when the file declares it, otherwise as UTF-8.
- **FuzzyDiff emulation.** The test-support `fuzzy_diff` does not reproduce the TOPPBase
  log messages (version warning, timing, "Cannot read input file"), and
  `FuzzyDiffOutcome::log` carries the comparator log or a one-line reason; the ported tool
  (`openms::cli::tools::FuzzyDiff`) writes them. `tests/topp_fuzzy_diff.rs` holds the
  emulation's exit code to the tool's and the Release build's on every oracle case it
  models. It keeps one answer of its own: an unreadable input such as a directory is a
  failed comparison (10), where the tool and the source exit 12. `-sort` compares in memory: the verdict and exit code match, and the
  log names the inputs instead of temporary files. The INI reader handles the ParamXML
  subset INI files use (`NODE`, `ITEM`, `ITEMLIST`, `LISTITEM`, quoted attributes, the
  five predefined entities and numeric character references). It is cross-checked
  against the crate's full `format::paramxml` reader on the pinned file.

## Decoded comparison

Line-level `FuzzyDiff` cannot compare Rust-written XML with C++-written XML. The
declaration encoding, attribute order, metadata order and index offsets all differ
without any difference in content (plan risk 11). Decision D6 therefore compares decoded
values:

- **Numbers** use `compare_numbers` with a fresh relative maximum and the invocation's
  tolerance, e.g. `Tolerance::from_settings(&FuzzyDiffSettings::upstream()?)`. A test
  pins this rule to the executed C++ verdicts of every single-number case in the corpus.
- **Structure is exact:** counts, order, metadata key sets, value types (`float` against
  `int` fails), hull point order, the presence of optional values, and string fields
  (names, native ids, software). Charges and MS levels are exact. Integer metadata
  values and integer data arrays follow the number rule, as their text tokens do in
  `FuzzyDiff`. This is stricter than line-level `FuzzyDiff` for strings, whose embedded
  numbers FuzzyDiff compares fuzzily, and for charges and MS levels.
- **Nested metadata without a dedicated walk** is compared through its pretty `Debug`
  rendering, line by line, with the same comparator and tolerance. This covers
  experimental settings, instrument settings, acquisition information, source files,
  products, CV terms, identifications and identification-graph references. Field names,
  variants, strings and element counts must agree, and numbers follow the rule. A
  mismatch reports the field path and the `Debug` line.
- **`ignore_unique_ids`** is the decoded counterpart of `-whitelist "id="`. It skips the
  unique ids of maps, features and subordinates, the native ids of spectra and
  chromatograms, precursor spectrum references and the SQL run id. Line-level `id=`
  skipping drops the whole line, and with it the other attributes on those lines; the
  decoded comparison still compares them.
- **Load provenance** (`loaded_file_path`, `loaded_file_type`) is not content and is never
  compared.
- **Reports** give the first mismatch in walk order, with a path such as
  `features[1].convex_hulls[1].points[3].mz` or `spectra[1].float_data_arrays[0].data[1]`.

Index offsets (`-whitelist offset indexListOffset` for mzML registrations) have no decoded
counterpart, because an index is not decoded content.

## Checked boundaries and evidence

| Evidence | Tier | What |
|---|---|---|
| `oracle/comparator_cases.tsv`, `oracle/comparator_results.tsv` | 1, executed differential | 137 cases run by `../oracle/fuzzy-string-comparator/driver.cpp` on the unmodified `FuzzyStringComparator.cpp.o` in product-sdk `libOpenMSTestFramework.a`. They cover every class-test literal, number tokens, the number rule, line structure, CR, whitelists, report formatting, instance reuse and files. Rust matches every verdict and every log byte, after replacing the working directory and masking non-number values, except the one documented hexadecimal case. |
| `oracle/tool_runs.tsv` | 1, executed differential | Exit codes of 39 `FuzzyDiff` invocations on product-sdk: TOPP_FuzzyDiff_1..4 as registered, the retained-vs-current pairs with and without their whitelists, one-digit mutations, parameter and INI errors, empty and missing files, matched-whitelist splitting and sorting. Rust `fuzzy_diff` returns the same code for each. |
| `oracle/FileInfo_3.current.txt`, `oracle/FileInfo_17.current.txt`, `oracle/FeatureFinderCentroided_1.current.featureXML` | 1, oracle-generated | Current C++ outputs of the registered invocations, byte-identical across two runs. They pass against the retained expectations under the registered whitelists (FileInfo_3: intensity-line width; FileInfo_17: CRLF; FFC_1: `id=`), in C++ `FuzzyDiff` and in Rust, and fail without the whitelists. The decoded FFC_1 comparison passes with ids ignored. |
| `fuzzydiff/*.beyond.*`, `*.within.*`, `*.id_change.*` | 1, oracle-judged | One-byte mutations of the retained expectations: beyond tolerance fails, within tolerance and an `id=` line change pass, in C++ and in Rust, line-level and decoded. |
| `retained/*`, `FuzzyDiff.ini` | retained upstream | Byte copies of the test-data 0cb15f2 fixtures. |
| class-test sections | 3, source review | All 25 `START_SECTION`s are accounted for. 15 sections carry assertions and are transcribed: constructor, destructor, the six setters, the three whitelist accessors, `compareStrings`, `compareStreams`, `compareFiles` and `[EXTRA]` non-finite numbers, including the "magic" log line counts 17, 36 and 246 and the seven failure headers. 8 are `NOT_TESTABLE`: copy and assignment (not implemented, and the Rust type is not `Clone`) and six getters "tested along with set-method". 2 are commented out (`reportFailure_`, `reportSuccess_`), and the differential checks their output byte for byte. |
| bounds, malformed INI, read errors, decoded paths | 4, Rust-only | Report and byte limits, INI rejection and entity decoding, read errors failing the comparison, `fuzzy_diff` on a directory or empty name, and the mismatch path for every decoded field family. |

The oracle environment, compiler, binary and library hashes, and the pinned input
hashes are in `../oracle/fuzzy-string-comparator/manifest.json`. No case depends on a
Debug-only precondition.

## Upstream issue candidates

Recorded for `OpenMS_CPP_ISSUES.md` (integrator-owned). Each is confirmed by an executed
oracle case against core 4fdec46, whose source is identical to bc9cc12:

1. `compareLines_` raises `ratio_max_` only when it reports a failure (665-683). A passing
   report therefore always prints `relative_max: 1`, and "Maximum relative error was
   attained at these lines" names the last line with any ratio above one, not the maximum
   (`rep_success_last_ratio_line`).
2. `compareStreams` resets only `is_status_success_` (805). On a reused instance a failure
   leaves `ratio_max_` raised, so a later comparison accepts ratios beyond the tolerance
   but below that maximum (`reuse_ratio_max_carries`: `1` against `1.2` passes after `1`
   against `1.5` failed).
3. A negative quotient that underflows to `-0.0` passes the sign test (`ratio < 0` is
   false), and its reciprocal `-inf` never exceeds the maximum. The pair `-1e-300` against
   `1e300` is accepted with absdiff 0 (`num_ratio_underflow_negative_zero`).
4. The libc++ `strtod` fallback in `fromCharsFloat` (167-211) accepts hexadecimal floats,
   which `std::from_chars` does not, so verdicts depend on the standard library
   (`tok_hex_vs_decimal`). It also writes `strtod`'s result into the reset number of a
   letter element, which changes the failure report.
5. `openInputFileStream_` reports "Error opening first input file" for the second file
   too (896; `file_missing_second`).
6. `setAcceptableRelative(NaN)` stores NaN, which then accepts every ratio
   (`num_ratio_nan_setter`).
7. `FuzzyDiff -sort` compares temporary copies, so the same-name "cheating" check never
   fires and a file compared with itself passes (`fd_sort_same_file`).

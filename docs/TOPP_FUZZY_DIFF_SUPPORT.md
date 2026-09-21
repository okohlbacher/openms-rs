# FuzzyDiff (TOPP tool)

`FuzzyDiff` compares two files and tolerates numeric differences. The upstream
test suite judges almost every tool output with it (`${DIFF}` at
`topp/CMakeLists.txt:59`, test-data `0cb15f2`). This port runs it on the TOPP
lifecycle over the library comparator. `python3 tools/core_sdk_coverage.py`
counts it as the ninth validated TOPP workflow, through
`tests/data/topp_fuzzy_diff_provenance.json`, whose tier-1 evidence is not
retained upstream output: the four upstream registrations retain no output
file, and `WILL_FAIL` on three of them is their whole expectation. The
evidence is executed Release runs: the pinned Release build ran `FuzzyDiff`
80 times (the four registrations among them), and each run's exit code and
both streams are retained and compared (*Checked boundaries and evidence*).

| Rust file | Covers |
|---|---|
| `src/cli/tools/fuzzy_diff.rs` | `OpenMS4-topp/src/FuzzyDiff.cpp` (topp `174b576`): registration and `main_` |
| `src/bin/FuzzyDiff.rs` | the executable, `openms::cli::run::<FuzzyDiff>`; `[[bin]]` requires `paramxml` |
| `src/concept/fuzzy_string_comparator.rs` | `CONCEPT/FuzzyStringComparator.h/.cpp` (core `bc9cc12`, `src/testframework`), plus the tool's two input transforms `sorted_lines` and `parse_matched_whitelist`; see [its support document](FUZZY_STRING_COMPARATOR_SUPPORT.md) |
| `tests/topp_fuzzy_diff.rs` | the executed differential against the Release build, the four upstream registrations, the deliberate divergences and the native bounds |
| `tests/data/topp_fuzzy_diff/` | synthetic inputs, the oracle case list, the Release build's streams with the oracle's placeholders, and its `-write_ini` file |
| `tests/data/topp_fuzzy_diff_provenance.json` | pins, hashes, evidence tiers |

The tool needs `paramxml` only, as the TOPP framework does. The comparator and
the two transforms are always built, so the test suite's own `${DIFF}`
reproductions (`tests/support/fuzzy_string_comparator.rs`) need no feature.

## API mapping

### Parameters (`registerOptionsAndFlags_`, `FuzzyDiff.cpp:54-86`)

Registered in source order, with the source's descriptions, defaults and
restrictions; `-write_ini` equals the Release build's file.

| Parameter | Default | Restriction | Advanced | Rust |
|---|---|---|---|---|
| `in1`, `in2` | empty | required input files | no | `ToolContext::string`; the framework's input check exits 7 (missing), 1 (not found), 2 (unreadable) or 4 (empty regular file) first |
| `ratio` | 1 | min 1 | no | `set_acceptable_relative` |
| `absdiff` | 0 | min 0 | no | `set_acceptable_absolute` |
| `whitelist` | `<?xml-stylesheet` | - | yes | `set_whitelist`; a command-line list replaces the INI list |
| `matched_whitelist` | empty (`ListUtils::create<std::string>("")` splits into nothing) | `first:second` entries | yes | `parse_matched_whitelist`, then `set_matched_whitelist` |
| `verbose` | 2 | 0-3 | no | `set_verbose_level` |
| `tab_width` | 8 | min 1 | no | `set_tab_width` |
| `first_column` | 1 | min 0 | no | `set_first_column` |
| `sort` | false | flag | no | `sorted_lines` on both texts, compared in memory |

The four empty lines the source adds, before `in1`, `ratio`, `whitelist` and
`sort`, are registered too (`add_empty_line`), so the usage text is the
source's line for line.

### `main_` (`FuzzyDiff.cpp:88-208`)

| Source | Rust |
|---|---|
| `getStringOption_`/`getDoubleOption_`/`getStringList_`/`getIntOption_`/`getFlag_` | `ToolContext` accessors; the framework validates eagerly what the source validates lazily, in the same order (`docs/TOPP_CLI_SUPPORT.md`, *Validation order*) |
| `writeDebug_` of both whitelists at debug level 1 (106-108) | not written: the framework does not port `writeDebug_` |
| split of each `matched_whitelist` entry, `IllegalArgument` on a part count other than two (112-126) | `parse_matched_whitelist`; the tool prints `Error: Unexpected internal error (<entry> does not have the format String1:String2)` and exits 8, as `TOPPBase`'s `BaseException` arm (`TOPPBase.cpp:495-499` at cli c19e494) |
| comparator setup (128-134) | the same seven setters on `FuzzyStringComparator`, logging into its buffer |
| `-sort`: `sortFile` into temporary files, then `compareFiles` on those (136-197) | `read_for_sort` and `sorted_lines`, compared with `compare_bytes` under the input names (native difference 1) |
| `compareFiles(in1, in2)` (190) | `compare_files`; the report goes to the output stream |
| `EXECUTION_OK` / `PARSE_ERROR` (199-207, "TODO think about better exit codes") | `ExitCode::ExecutionOk` / `ExitCode::ParseError` |

### Exit codes

| Code | When | Source |
|---|---|---|
| 0 | no difference | `EXECUTION_OK` |
| 1 | an input does not exist | framework input check (`FileNotFound`) |
| 2 | an input cannot be read | framework input check (`FileNotReadable`) |
| 4 | an input is an empty regular file | framework input check (`FileEmpty`) |
| 6 | a parameter outside its restriction, an unknown or malformed option, an invalid INI | framework |
| 7 | `in1` or `in2` not given or empty | framework |
| 8 | a malformed `matched_whitelist` entry | `IllegalArgument` → `UNKNOWN_ERROR` |
| 10 | a difference, or the same name twice ("That's cheating!") | `PARSE_ERROR` |
| 11 | an input beyond the comparator's 1 GiB bound, or a `-sort` input beyond 2^26 lines | native (difference 3) |
| 12 | an input that cannot be read, such as a directory | the source's `std::ios_base::failure` → `INTERNAL_ERROR` (difference 2) |

Only 0 and 10 are exit codes `main_` returns, and only they are followed by
the closing line `FuzzyDiff took … .`, as in the Release build; every other
code is a thrown exception (`ToolError`) or a refusal before `main_`, and ends
without it. The `IllegalArgument` of 8 and the `FileNotFound` of an unopenable
`-sort` input are written to the `-log` file as well; the 12 of an unreadable
input reaches only the error stream, as the source's initialisation catch
writes it (Release oracles `fd_*_log` of `../oracle/topp-exception-exits`).

## Preserved source conventions

- The report is the comparator's, byte for byte: every verbose level (0 writes
  nothing, 1 only failures, 2 adds the summary, 3 continues after a failure),
  tab-width-aware columns and the first-column offset, whitelist counts
  including an empty entry, the raw bytes of a non-UTF-8 line, and the
  absolute path the failure report gives a relative name.
- Only one of `ratio` and `absdiff` has to be satisfied. A passing report
  always prints `relative_max: 1` and names the last line with any ratio above
  one as the maximum (CPP-230).
- The same name twice is refused before anything is read, even through a
  relative path; the same file under two spellings is compared and passes.
- A command-line `-whitelist` replaces the INI list; `-matched_whitelist` pairs
  apply in either order; `a:` is a valid pair and `abc`, `a:b:c` and the empty
  entry are not.
- `-sort` keeps the first line and sorts the rest bytewise, splitting at `\n`
  only, so `-sort` on a file and itself passes (CPP-236).
- A directory given as an input passes the framework's input check, as the
  source's `inputFileReadable_` lets it (it checks emptiness only for regular
  files), and fails when it is read.
- The tool runs serially, as the source does: the comparator has no parallel
  region, so `-threads` has nothing to size and no worker pool is built.

## Native differences

1. **`-sort` compares in memory.** The source writes each sorted input to
   `<tmp>/<name>.sorted.<unique>.tmp`, compares those and deletes them, so its
   report names files that no longer exist. The port compares the same texts
   in memory and names the inputs as given. The verdict, the line and column
   numbers and every other byte of the report are the source's (executed:
   `sort_fail`, `sort_pass_verbose_2`, `sort_reordered_rows`,
   `sort_same_file`). There is no temporary-file failure (`UnableToCreateFile`,
   exit 5) to reproduce.
2. **An input that cannot be read exits 12 in both modes.** Without `-sort`
   the source's line reader lets libstdc++'s `std::ios_base::failure` escape
   `compareFiles`, and `TOPPBase` reports `Unable to initialize or run
   FuzzyDiff: basic_filebuf::underflow error reading the file: Is a
   directory` with `INTERNAL_ERROR` (executed: `directory_vs_file`,
   `file_vs_directory`, `directory_vs_directory`). The port prints the same
   line, the reports written before the failure on the output stream, and
   exits 12. With `-sort` the source's `std::getline` swallows the error, so a
   directory sorts as an empty text and is compared (`sort_directory`, exit
   10); the port refuses with 12 and `error reading the file '<path>'`, because
   it does not compare a text it could not read in full. The comparator
   reports the failure through `input_failure` (native API; its log line is
   what it always was).
3. **Bounds.** An input beyond the comparator's `MAX_INPUT_BYTES` (1 GiB) is
   refused before it is read, with the comparator's message on the error
   stream and exit 11 (`INCOMPATIBLE_INPUT_DATA`), not the 10 that would claim
   a difference; the source has no limit. A `-sort` input of more than
   `FuzzyDiff::MAX_SORT_LINES` (2^26) lines is refused the same way before it
   is sorted, because sorting holds a slice reference per line; each input is
   sorted and its unsorted text dropped before the second is read. The report is buffered up to
   `MAX_LOG_BYTES` (256 MiB) and a cut report is announced on the error
   stream; only verbose level 3 can reach it, and the verdict is unaffected.
4. **A number that rounds to zero.** The comparator accepts `1e-400` as 0, as
   the source's comment at `FuzzyStringComparator.cpp:205-208` states and its
   libc++ build does. The Linux Release build's `std::from_chars` rejects it,
   so there `1e-400` reads as letters and the comparison fails
   (`token_underflow`, `token_underflow_absdiff`: C++ 10, port 0). A subnormal
   (`1e-310`), the smallest normal and an overflow (`1e400`, letters in both)
   agree. The comparator keeps the behaviour it had as test support; the C++
   verdict depends on its standard library (CPP-233).
5. **Integer parameters** are `i64` in the parameter tree and converted to the
   source's `int`; both of the framework's readers already refuse a value
   beyond `i32`, so the conversion's refusal is defensive.

The framework's lines in this tool's streams are the Release build's: the
closing line after `main_` returns, the missing-file diagnostic, and a failed
strict update's `Parameters passed to 'FuzzyDiff' are invalid...` before its
diagnostic. The source's two `writeDebug_` lines about the whitelists reach
the `-log` file from debug level 1 (Release oracle `fd_debug1_log`). One
framework difference is left, of the environment rather than the tool: a tool
driven in process writes to explicit streams, which are not probed, so the
`stty` line the Linux Release build prints before its usage text has no
counterpart (`docs/TOPP_CLI_SUPPORT.md`, *On a console*).

## Checked boundaries and evidence

| Evidence | Tier | What |
|---|---|---|
| `tests/data/topp_fuzzy_diff/oracle/` | 1, executed differential | The pinned Release build `openms4-release-bc9cc12-c19e494-174b576` ran `FuzzyDiff` for 80 invocations on ibminode06, twice, with identical exit codes and streams after the oracle's own masking of paths, timing figures and temporary names (`../oracle/fuzzy-diff-tool`, binary sha256 `ee789a71…`). 75 run in process and match on exit code, output stream and error stream byte for byte, the Release streams as recorded apart from the Linux `stty` probe line and `-sort`'s temporary names (`tests/topp_fuzzy_diff.rs::release_streams`), the port's closing line masked as the oracle masks the Release build's; 2 run through the executable with a working directory and match the same way; 3 are the deliberate divergences above, asserted on their own; `-write_ini` matches the C++ file line by line and as a decoded tree. Until the TOPPBase completion, the comparison also dropped the closing line, rewrote the missing-file text and reordered the strict-update lines to the framework of the day; those normalisations are gone. |
| `tests/data/topp_exception_exits/fd_*_log` | 1, executed differential | Three more Release runs with `-log` (`../oracle/topp-exception-exits`): the malformed whitelist (8, its line in the log, no closing line), a directory as `-in1` (12, the initialisation catch's line on the error stream only, no log file) and a comparison at debug level 1 (the whitelist lines in the log). |
| TOPP_FuzzyDiff_1..4 | 1 | The four registrations at `topp/CMakeLists.txt:138-144` on their pinned inputs: `_1` 10 (cheating), `_2` 10 (ratio at line 22), `_3` 0, `_4` 1 (missing input), each checked against its `WILL_FAIL` expectation and then against the Release build's streams. No registration retains an output file. |
| the 35 product-SDK invocations | 1 | Re-run on the Release build with their streams; every exit code equals the product-SDK one recorded in `tests/data/fuzzy_string_comparator/oracle/tool_runs.tsv`. |
| emulation agreement | 4 | `tests/support/fuzzy_string_comparator.rs::fuzzy_diff`, which tests without `paramxml` use, gives the tool's and the Release build's exit code on every oracle case it models, except that it keeps its own answer (10) for a directory. |
| bounds | 4 | 1 GiB+1 sparse inputs with and without `-sort` exit 11 before reading; a `-sort` input of 2^26+1 lines exits 11 before sorting; an INI integer beyond `i32` never reaches the tool. |
| library unit tests | 4 | `InputFailure` recording and its unchanged log line, reports kept before a failed read, input names, `sorted_lines`, `parse_matched_whitelist`. |

## C++ issue candidates

For `OpenMS_CPP_ISSUES.md` (integrator-owned):

1. CPP-233's libstdc++ half is now executed and retained: the Release build
   rejects `1e-400` (`token_underflow`), accepts `1e-310` and the smallest
   normal (`token_subnormal`, `token_min_normal`) and rejects `0x10`
   (`token_hex`), so the comment at `FuzzyStringComparator.cpp:205-208` is
   wrong for libstdc++ exactly for values that round to zero.
2. A directory input makes `FuzzyDiff` exit `INTERNAL_ERROR` without `-sort`
   (an escaped `std::ios_base::failure`) and compare as an empty text with it
   (exit 10 against a readable file, executed; by the same reading two
   unreadable inputs would compare equal, not executed): the verdict for
   unreadable input depends on a formatting flag.
3. `FuzzyDiff -sort` reports and "Easy Access" lines name temporary files that
   `main_` deletes before it returns.

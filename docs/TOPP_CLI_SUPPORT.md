# The TOPP command-line framework

This group ports the `OpenMS4-cli` package — `TOPPBase` and the parameter
records around it — and the TOPP tools built on it. It is the layer that turns
the library into executables.

The framework needs the `paramxml` feature, because every TOPP tool supports
`-ini` and `-write_ini`. Sources are read at cli `c19e494`
(`TOPPBase.h` 576cea45…, `TOPPBase.cpp` 326b96f8…, `TOPPBase_test.cpp`
01e7de9f…).

| Source | Rust |
| --- | --- |
| `APPLICATIONS/TOPPBase.h`, `TOPPBase.cpp`: lifecycle | `src/cli.rs` |
| registration and `getDefaultParameters_` parameter part | `src/cli/spec.rs` |
| `get*_` accessors, file checks, `parseRange_`, run services | `src/cli/context.rs` |
| `getProcessingInfo_`, `addDataProcessing_` | `src/cli/processing.rs` |
| `printUsage_` | `src/cli/usage.rs` |
| `APPLICATIONS/ParameterInformation.h`, `.cpp`; `TOPPBase::ExitCodes` | `src/cli/parameter.rs` |

## API mapping

| Source | Native |
| --- | --- |
| `TOPPBase(name, description, official, citations, toolhandler_test)` | `trait Tool` with `NAME`, `DESCRIPTION`, `VERSION`; `official`, `citations` and `toolhandler_test` not ported |
| `main(argc, argv)` | `cli::run::<T>()`, or `run_with::<T>(args, out, err)` |
| `registerOptionsAndFlags_` | `Tool::register(&mut ToolSpec)` |
| `main_` | `Tool::run(&ToolContext)`, or `Tool::run_io(&ToolContext, out, err)` for tools that write to the streams |
| `getSubsectionDefaults_(section)` | `Tool::subsection_defaults` |
| `getSubsectionDefaults_()`, `getDefaultParameters_` | internal; `ToolSpec::to_param` is the parameter part |
| `getToolUserDefaults_` | not ported (reads `~/<Tool>.ini`, which would make runs depend on the account) |
| `parseCommandLine_` | internal |
| `handleWriteCommands_` | internal; `-write_ini` ported, CTD/CWL/JSON refused |
| `checkIfIniParametersAreApplicable_` | internal |
| `checkParam_` | not ported (warnings only; see native differences) |
| `TOPPBase::ExitCodes` | `ExitCode`, same 15 variants and discriminants |
| `ParameterInformation`, `ParameterTypes` | `ParameterInformation`, `ParameterType` |
| `paramEntryToParameterInformation_`, `getParamArgument_`, `paramToParameterInformation_` | `ParameterInformation::from_param_entry` |
| `register{String,Int,Double}Option_`, `registerFlag_` | `ToolSpec::register_*` |
| `registerInputFile_`, `registerOutputFile_`, `registerOutputPrefix_`, `registerOutputDir_` | `ToolSpec::register_*` |
| `register{String,Int,Double}List_`, `register{Input,Output}FileList_` | `ToolSpec::register_*_list` |
| `registerFullParam_`, `registerParamSubsectionsAsTOPPSubsections_`, `registerTOPPSubsection_` | `ToolSpec::register_full_param`, `ToolSpec::topp_subsections` |
| `setValidStrings_`, `setValidFormats_`, `setMin/MaxInt_`, `setMin/MaxFloat_` | `ToolSpec::set_*`; `setValidFormats_`'s `force_OpenMS_format` not ported |
| `registerSubsection_`, `addText_`, `addEmptyLine_` | `ToolSpec::register_subsection`, `add_text`, `add_empty_line` |
| `findEntry_`, `getParameterByName_` | `ToolSpec::find` |
| `get{String,Int,Double}Option_`, `get{String,Int,Double}List_`, `getFlag_`, `getParam_` | `ToolContext::{string,int,double,string_list,int_list,double_list,flag,param}` |
| `getOutputDirOption` | not ported; no ported tool registers an output directory |
| `getParamAs*_`, `getSubsection_` | internal |
| `getIniLocation_`, `toolName_` | `ToolContext::ini_location`, `tool_name` |
| `test_mode_`, `debug_level_`, `log_type_` | `ToolContext::test_mode`, `debug_level`, `progress_log_type` |
| `setMaxNumberOfThreads` | `ToolContext::thread_policy` |
| `UniqueIdGenerator::setSeed(19991231235959)` under `-test` | `ToolContext::unique_id_generator`, `TEST_MODE_UNIQUE_ID_SEED` |
| `getProcessingInfo_` (both overloads) | `ToolContext::processing_info` |
| `addDataProcessing_` (`ConsensusMap`, `FeatureMap`, `PeakMap`) | `ToolContext::add_data_processing`, trait `AddDataProcessing` |
| `inputFileReadable_`, `outputFileWritable_` | `cli::input_file_readable`, `cli::output_file_writable` |
| `parseRange_(text, double&, double&)` | `parse_range` |
| `parseRange_(text, Int&, Int&)` | not ported; no ported tool uses it |
| `version_`, `verboseVersion_` | `Tool::VERSION`, `TOPP_PRODUCT_VERSION`; the verbose usage line is W2.2 |
| `printUsage_` | internal, reached by `--help` / `--helphelp` and command-line errors |
| `writeLogInfo_`, `writeLogWarn_`, `writeLogError_`, `writeDebug_`, `enableLogging_` | not ported; diagnostics go to the `out`/`err` streams, and `-log` is inert |
| `getDocumentationURL`, `Citation`, `cite_openms` | not ported |

Registration is a builder rather than protected methods on the tool, so the
registered set is inspectable without running the tool, and a tool cannot mutate
its own resolved parameters mid-run.

## Lifecycle and exit codes

`run_with` follows `TOPPBase::main` phase by phase. The exit code depends on the
phase, as in the source, where errors inside the run-phase `try`
(`TOPPBase.cpp:258-425`) take the inner catch (430-499) and anything else takes
the initialisation catch (505-514).

The *Test* column names the cases that assert each row, in
`tests/topp_cli_lifecycle.rs` unless another file is named. *Evidence* says what
backs the expected exit code:

* **tier 1**: an executed C++ case on the product SDK, named as in
  `../oracle/topp-cli-lifecycle/manifest.json`, or as in
  `../oracle/topp-cli-lifecycle/ini_read_failures/manifest.json` where marked †;
* **tier 4, derived**: read from the cited source lines, with no executed case;
* **tier 4, native**: a native mapping. The Rust `Error` enum is coarser than the
  source's exceptions, so there is no C++ case to execute; the row says what the
  source does instead where the two differ.

| Phase | Condition | Exit | Source | Test | Evidence |
| --- | --- | --- | --- | --- | --- |
| registration | tool registration fails | 6 | 505-508 | `an_initialisation_failure_is_illegal_parameters` | tier 4, derived |
| parse | more than `MAX_ARGUMENTS` tokens or `MAX_ARGUMENT_BYTES` bytes | 6 | — | `an_oversized_command_line_is_refused_before_parsing` | tier 4, native bound |
| parse | a flag followed by text | 6 | 186-191, 2336-2348 | `a_flag_followed_by_text_is_refused`, `upstream_flag_with_trailing_arguments` | tier 1: `flag_with_trailing_text` |
| parse | a number that does not convert, including `-section:name` values | 6 | 186-191, 2365-2407 | `a_non_numeric_integer_is_refused`, `a_non_numeric_subsection_integer_is_refused` | tier 1: `int_option_not_numeric`, `algorithm_peakcount_not_int` |
| parse | no arguments at all | 6 | 227-232 | `a_bare_invocation_is_refused`; also `tests/topp_dta_extractor.rs`, `tests/topp_baseline_filter.rs`, `tests/topp_map_normalizer.rs` | tier 1: `no_arguments` |
| parse | `--help`, `--helphelp` | 0 | 235-239 | `usage_and_exit_codes_follow_the_source_contract` in `tests/topp_dta_extractor.rs` | tier 1: `help` (exit code; the stream differs, see *Usage text*) |
| parse | an unknown option, including an unknown `-section:name` | 6 | 242-247 | `an_unknown_option_is_refused`, `an_unknown_subsection_parameter_is_refused` | tier 1: `unknown_option`, `algorithm_bogus` |
| parse | trailing text | 6 | 250-255 | `trailing_text_is_refused` | tier 1: `trailing_text` |
| parse | trailing text after several options, listed in command-line order | 6 | 2436-2444 | `a_command_line_at_the_argument_bound_is_parsed` | tier 4, derived |
| write | `-write_ini` target not writable | 5 | 2609 | `write_ini_to_an_unwritable_path_is_refused` | tier 1: `write_ini_unwritable` |
| write | `-write_ini` | 0 | 2606-2641 | `write_ini_ignores_command_line_values`, `write_ini_with_an_invalid_ini_value_keeps_the_default` | tier 1: `write_ini_with_cli_value`, `write_ini_with_ini` |
| write | `-write_ctd`, `-write_cwl`, `-write_nested_cwl`, `-write_json`, `-write_nested_json` | 12 | 2643-2683 | `tool_description_writers_are_refused_explicitly` | tier 1: `write_cwl`, `write_nested_cwl`, `write_json`, `write_nested_json` exit 12 without TDL; `write_ctd` exits 0 there, see *Refused writers* |
| INI | `-ini` file missing | 1 | 296, 436-441 | `a_missing_ini_is_input_file_not_found` | tier 1: `ini_missing` |
| INI | `-ini` file not readable, before a run or with `-write_ini` | 2 | 296, 2630, 448-453 | `an_unreadable_ini_is_input_file_not_readable`, `a_fifo_ini_this_user_cannot_open_is_input_file_not_readable` | tier 1: `ini_unreadable`†, `write_ini_ini_unreadable`† (a regular file); `ini_fifo_denied`†, `write_ini_ini_fifo_denied`† (a FIFO with mode 000: not queried for readability before the load, whose open is refused) |
| INI | `-ini` file malformed | 3 | 296, 460-465 | `a_malformed_ini_is_input_file_corrupt` | tier 1: `ini_malformed` |
| INI | `-ini` names a directory, before a run or with `-write_ini` | 3 | 296, 2630, 460-465 | `an_ini_directory_is_input_file_corrupt` | tier 1: `ini_directory`†, `write_ini_ini_directory`† |
| INI | `-ini` names the character device `/dev/null`, before a run or with `-write_ini` | 3 | 296, 2630, 460-465 | `a_character_device_ini_is_input_file_corrupt` | tier 1: `ini_dev_null`†, `write_ini_ini_dev_null`† |
| INI | `-ini` names a FIFO this user can open | as the INI it carries once a writer opens it; blocks until then | 296, 2630 | `an_ini_fifo_is_read_once_a_writer_opens_it` (exit 0 with `-write_ini`) | tier 4, native: with no writer the C++ tool blocks opening it too, so there is no exit code to execute |
| INI | INI without a section for this tool | 0, warning | 1957-1966 | `an_ini_for_another_tool_warns_and_applies_defaults` | tier 1: `ini_foreign_section` |
| update | unknown parameter (INI item, `-instance`, top-level `common:` value) | 6 | 338-343 | `an_unknown_ini_item_is_refused`, `instance_on_the_command_line_is_refused`, `a_common_top_level_value_is_refused_as_in_the_source` | tier 1: `ini_unknown_item`, `instance_on_command_line`, `ini_common_top_level` |
| update | value outside its restrictions | 6 | 338-343 | `an_invalid_subsection_value_on_the_command_line_is_refused`, `an_invalid_subsection_value_in_an_ini_is_refused` | tier 1: `algorithm_movetype_sideways`, `ini_invalid_subsection_value` |
| update | value of another type | 6 | 338-343 | `an_ini_value_of_the_wrong_type_is_refused` | tier 1: `ini_value_type_mismatch` |
| update | INI written by another version | 0, notice on stdout | 355-366 | `an_ini_from_another_version_is_noted` | tier 1: `ini_version_mismatch` |
| validation | required value missing or empty | 7 | 1398-1406, 466-474 | `a_missing_required_output_is_missing_parameters` | tier 1: `missing_required_output` |
| validation | input missing / unreadable / empty | 1 / 2 / 4 | 1968-1995 | `a_missing_input_is_input_file_not_found`, `an_unreadable_input_is_input_file_not_readable`, `a_zero_byte_input_is_input_file_empty` | tier 1: `missing_input`, `unreadable_input`, `zero_byte_input` |
| validation | output not writable | 5 | 1997-2013 | `an_unwritable_output_is_cannot_write_output_file` | tier 1: `unwritable_output` |
| validation | output extension not registered | 6 | 1595-1608 | `a_wrong_output_extension_is_refused` | tier 1: `output_wrong_extension` |
| run | `Error::Parse` from the tool | 3 | 460-465 | `a_corrupt_input_is_input_file_corrupt`, `run_phase_errors_map_like_the_source_inner_catch` | tier 1: `corrupt_input` |
| run | `Error::Io`, not found | 1 | 436-441 | `run_phase_errors_map_like_the_source_inner_catch` | tier 4, native: as `FileNotFound` |
| run | `Error::Io`, permission denied | 5 | 430-435 | `run_phase_errors_map_like_the_source_inner_catch` | tier 4, native: an I/O error does not say whether it read or wrote; a tool's inputs are checked for readability before its body runs, so a denial there is taken as a failed write. INI files never take this row |
| run | any other `Error::Io` | 8 | 495-499 | `run_phase_errors_map_like_the_source_inner_catch` | tier 4, native |
| run | `Error::InvalidValue`, `InvalidRange` | 6 | 475-480 | `run_phase_errors_map_like_the_source_inner_catch` | tier 4, native: as `InvalidParameter`; the source's own `InvalidValue` and `InvalidRange` exceptions reach its `BaseException` arm, exit 8 |
| run | `Error::MissingInformation` | 7 | 466-474 | `run_phase_errors_map_like_the_source_inner_catch` | tier 4, native: as `RequiredParameterNotGiven`; the source's `MissingInformation` exception reaches its `BaseException` arm, exit 8 |
| run | `Error::Unsupported`, `UnsortedData` | 11 | the tools' own `INCOMPATIBLE_INPUT_DATA` returns | `run_phase_errors_map_like_the_source_inner_catch` | tier 4, native |

A malformed INI file is exit 3, not 6: the source loads it inside the run-phase
`try` (`TOPPBase.cpp:258`, `296`), so its `ParseError` is `INPUT_FILE_CORRUPT`.
The early plan classified it as an initialisation error; the oracle settles it.

An unreadable INI file is exit 2 for the same reason. `XMLFile::parse_` checks
only that the file exists (`XMLFile.cpp:124`); xerces cannot open it, and
`XMLHandler::fatalError` asks `FileHandler::getTypeByContent` for a file-type
hint before raising its `ParseError` (`XMLHandler.cpp:49-50`). For an
uncompressed file that reads through `TextFile` (`FileHandler.cpp:402`), which
throws `FileNotReadable` (`TextFile.cpp:40-43`; all at core bc9cc12). A readable
directory does reach the `ParseError`, exit 3. `-write_ini` loads its `-ini`
inside the same `try` (2630), with the same codes. This port checks existence,
then readability for a regular file or a directory, then directories, before
loading, and prints the source's `FileNotReadable` wording; for a directory it
prints the source's `ParseError` wording without the file-type hint.

Anything else is not queried for readability, because `file::readable` answers
`false` for a device or FIFO without opening it, which reported the readable
`/dev/null` as unreadable. The load opens such a file itself: a denied open is
exit 2 with the same wording, as for a FIFO with mode 000 in the oracle's
`ini_fifo_denied` and `write_ini_ini_fifo_denied`; any other open or read
failure is exit 3, and a document that does not parse exit 3. `/dev/null` reads
as an empty document, exit 3 with an `Error: Unable to read file (` line on both
paths, as the oracle's `ini_dev_null` and `write_ini_ini_dev_null` do; the
source's line also names the file and xerces's message. A FIFO this user can
open is opened like a file, so it waits for a writer; the C++ tool also blocks
on such a FIFO with no writer.

## Preserved source conventions

The sixteen common parameters are registered in source order and with source
descriptions: `ini`, `log`, `instance`, `debug`, `threads`, `write_ini`,
`write_ctd`, `write_nested_cwl`, `write_cwl`, `write_nested_json`, `write_json`,
`no_progress`, `force`, `test`, `-help`, `-helphelp`. As in source, the two help
flags are registered with a leading dash, so their command-line tokens carry two.

**Command line.** Tokens are read from last to first, so an option finds its
values already queued (`parseCommandLine_`). Every entry of a registered
subsection is an option too, as `-section:name`, typed by its default: a string
that defaults to `false` and is restricted to `true`/`false` becomes a flag, and
numbers convert while parsing with `StringUtils::toInt32` and `toDouble` rules,
whose messages are reproduced. An option given twice keeps its last value with a
warning. Unknown options are reported before trailing text, and both after
`--help`.

**Defaults.** The tree the update starts from is `getDefaultParameters_`: the
registered parameters under `<tool>:1:` without `ini`, `instance`, the help flags
and the `write_*` requests; flags as `false` restricted to `true`/`false`; the
`<tool>:version` item; the tool and instance section descriptions; and the
algorithm subsections, ordered by name as the source's `std::map` orders them.

**INI merge and strict update.** The command line comes first, then the INI
instance section, the `common:<tool>:` section and the `common:` section; each
adds only what is not yet present. The result updates the defaults with
`fail_on_invalid_values` and `fail_on_unknown_parameters` set
(`TOPPBase.cpp:339-342`). A value not found by its full key is matched by a
unique leaf name, as `Param::update` does. The diagnostics use the source's
words: `Unknown (or deprecated) Parameter '…' given in outdated parameter file!`,
`Invalid string parameter value '…' for parameter '…' given! Valid values are:
'…'. Updating failed!`, `Parameter '…' has changed value type!` and
`Parameters passed to '<tool>' are invalid. To prevent usage of wrong defaults,
please update/fix the parameters!`.

**`-write_ini`** writes the defaults, updated leniently and verbosely from `-ini`
when given, and never command-line values (`TOPPBase.cpp:2606-2641`).

**File checks** follow `inputFileReadable_` and `outputFileWritable_`: an input
must exist, be readable and, unless it is a directory, hold a byte; an output
must be writable, and an output prefix is probed as `<prefix>_0`. The
writability query never creates or deletes a file under the caller's name.

**Run services.** `-test` seeds the unique-id generator with 19991231235959 and
makes the processing record machine-independent: software version
`version_string`, completion time 1999-12-31 23:59:59 and the single metadata
entry `parameter: mode` = `test_mode`. Otherwise the record carries the product
version, the current time and `parameter: <name>` for every resolved parameter.
A feature or consensus map gains one record, every spectrum and chromatogram of
an experiment gains one shared record, and under `-test` a consensus map's column
headers keep only file base names. `-no_progress` selects `ProgressLogType::None`
over `Cmd`.

**Threads.** The core is multithreaded: `rayon` is a default dependency behind
the `parallel` feature, and `src/concept/parallel.rs` holds the contract that a
parallel result is bit-identical to the serial one. `-threads` reaches a
computation through `ToolContext::thread_policy` (`Threads::from_cli`, where 0
means every core), instead of the source's process-wide `omp_set_num_threads`.

**Source behaviours reproduced but questionable** (C++ issue candidates, all
confirmed on the oracle):

1. `-instance` cannot be used. It is registered, but `getDefaultParameters_`
   leaves it out of the defaults while the command line keeps it, so the strict
   update rejects every run that passes it (`instance_on_command_line`). An INI
   instance section other than 1 is therefore unreachable; the context always
   reads `<tool>:1:`.
2. A `common:<tool>:` value of a subsection parameter overrides the instance
   section. The `common:` copy keeps the nested `<tool>:section:name` key, which
   the update finds by leaf name after the instance value has been applied
   (`ini_instance_and_common`).
3. A `common:<tool>:` value of a top-level parameter is rejected: the nested
   copy has no leaf match outside the root (`ini_common_top_level`).

## Native differences

**Usage text.** A successful `--help` prints to the output stream, because
`tests/topp_dta_extractor.rs` asserts it there; the source prints usage to
standard error in every case, as this port already does after a command-line
error. The usage layout, the `Version:` line, the documentation URL, citations
and the subsection list differ from the source (W2.2).

**Validation order.** The source checks each option lazily when `main_` reads it;
this port checks every registered option after the update, in registration
order. Exit codes agree whenever a tool reads all its options, which every ported
tool does.

**Formats.** Input formats are checked by file extension; the source detects the
type from content with `FileHandler::getType` and only warns about an unknown
type. Output extensions must be registered; the source accepts an unknown
extension. Both are W2.2.

**`-write_ini` parity** beyond the defaults-only rule — `supported_formats`,
the ISO-8859-1 declaration and a FuzzyDiff match with the C++ file — is W2.2.

**Refused writers.** The CTD, CWL and JSON tool-description writers are not
ported. Each request exits 12 with an explicit message, which is what the oracle
build without TDL reports for the four CWL and JSON writers; its `-write_ctd`
succeeds.

**Not ported:** JSON INI files (a `.json` `-ini` is read as XML and fails with 3),
`getToolUserDefaults_`, `checkParam_` warnings (unknown subsection, wrong type,
unknown parameter), the `-type` INI item, the run timing and peak-memory line,
`-log`, debug output helpers, `ToolHandler` and the `.tools.tsv` manifest
discovery, `INIUpdater`, `SearchEngineBase`, `MapAlignerBase` and
`TOPPExternalToolBase`.

**UpdateCheck** is ported in `src/system/update_check.rs`, behind the
non-default `network` feature, and is not called by the lifecycle. The source
skips it under `-test` and when `OPENMS_DISABLE_UPDATE_CHECK` is set; the oracle
runs with it disabled.

**Accessors are typed.** `ToolContext::double` does not widen an integer, and
`ToolContext::flag` rejects strings other than `true` and `false`, as the source
throws `WrongParameterType` and `InvalidParameter`. The strict update already
guarantees that every resolved value has its registered type.

**Non-finite values in processing records.** Outside `-test`,
`processing_info` records every resolved parameter as metadata, and `MetaValue`
holds only finite floating-point values, so a resolved `inf` or `nan` double
makes it fail with `Error::InvalidValue`, exit 6 once a tool propagates it. The
source's `DataValue` records any double: the product SDK's
SpectraFilterWindowMower accepts `-algorithm:windowsize inf` and `nan` and exits
0 (verifier probe, not in a manifest), although this port's window mower refuses
both values itself. No ported tool calls `processing_info` yet; the D4 retrofit
in W2.2 must record such values in a representable form or keep this difference.

**Diagnostics.** `Param::update_with_options` decides the update and applies it;
its report is worded natively, so the CLI derives the source wording from the
same entries. A missing file's parenthetical reads `does not exist`, where the
source says `could not be found`, because an existing assertion fixes the phrase.

**Errors** are typed `Result` values mapped at the boundary. A registration
error is exit 6 (the source's initialisation catch) rather than 12.
`parse_range` leaves both bounds unchanged on error. A command line is bounded
by `MAX_ARGUMENTS` and `MAX_ARGUMENT_BYTES` before parsing, and parsed in time
linear in its length: the text left after each option is gathered in reverse and
ordered once, where the source inserts each chunk at the front of its list
(`TOPPBase.cpp:2436-2439`). Number conversion
follows `std::from_chars`, so hexadecimal floats are rejected even though the
oracle's libc++ `strtod` fallback accepts them. `unique_id_generator` returns an
independent generator per call instead of seeding a process-wide singleton.

## Checked boundaries and evidence

`tests/topp_cli_lifecycle.rs` holds 65 cases, none ignored.

* **Oracle cases (tier 1 executed differential).** `../oracle/topp-cli-lifecycle/run.sh`
  runs 38 cases of the C++ product SDK (core 4fdec46, Debug, AppleClang 21) in a
  clean environment with a per-case directory; `manifest.json` records the build,
  binary hashes, inputs, exit codes, diagnostics and output hashes. Two executions
  agreed on all 37 shared cases. The subsection override is compared against the
  C++ output `tests/data/topp_cli_lifecycle/swm_algorithm_peakcount_1.mzML`; no
  case reaches a Debug-only precondition. `ini_read_failures.sh` in the same
  directory adds ten cases on the same binaries, recorded in
  `ini_read_failures/manifest.json`: an INI written by the tool, a readable
  control run with it, that INI unreadable before a run and with `-write_ini`,
  a directory as `-ini` on both paths, the character device `/dev/null` as
  `-ini` on both paths, and a FIFO with mode 000 as `-ini` on both paths. Each
  time cases were added, first the device cases and then the FIFO cases, the
  earlier cases ran again and matched every previous run in exit code, stderr
  and outputs; only timing figures in the control run's stdout changed. A FIFO
  this user can open is not a case, because with no writer the C++ tool blocks
  opening it.
* **Upstream class test (tier 3).** Transcribed with their literals:
  `getIniLocation_` (default), `getStringOption_` (default, command line, wrong
  type, unregistered, required), `getIntOption_`, `getDoubleOption_`,
  `getIntList_`, `getDoubleList_`, `getStringList_`, `getFlag_`,
  `inputFileReadable_`, `outputFileWritable_`, `parseRange_`, data processing
  methods, `getParam_`, misc options, subsection parameters, duplicate
  parameters and flag with trailing arguments. Not transcribed: the constructor,
  destructor and `main` sections (no object lifecycle; `main` is NOT_TESTABLE
  upstream), `setMaxNumberOfThreads` (NOT_TESTABLE; the native policy is
  tested), `getIniLocation_` with `-instance 5` and the `getStringOption_` INI
  cases (the source rejects `-instance` and top-level `common:` values in its
  strict update and only reads values left behind), the `-write_ini` file
  comparison (W2.2), `-log` and `Citation::toString` (not ported).
* **Native cases (tier 4)** cover `run_io` routing, the run-phase mapping of
  every `Error` variant, the initialisation failure, the argument bound (refused
  beyond it, parsed in full at it, with trailing-text chunks in command-line
  order) and the context services.

## DTAExtractor and executed differential evidence

`src/bin/DTAExtractor.rs` is the first TOPP tool. Its three upstream tests are
reproduced in `tests/topp_dta_extractor.rs` against the retained C++ outputs:

| Upstream test | Arguments | Result |
| --- | --- | --- |
| `TOPP_DTAExtractor_1` | `-rt :61` | `DTAExtractor_RT60.0.dta` byte-identical |
| `TOPP_DTAExtractor_2` | `-level 1` | `DTAExtractor_RT60.0.dta` byte-identical |
| `TOPP_DTAExtractor_3` | `-level 2 -mz :1000` | `DTAExtractor_RT140.0_MZ5.0.dta` byte-identical |

These fixtures were produced by the C++ tool, so agreement is **tier 1 executed
differential evidence** under `docs/DIFFERENTIAL_VALIDATION.md` for the whole
chain: command line, parameter validation, mzML reading, the source number
formatter that names the output files, and DTA writing. This is the first
validated TOPP workflow in the port.

Two source behaviors had to be matched to get there, and both were real gaps:

1. **Header list counts are advisory on reading.** The port rejected an mzML
   whose `count` attribute disagreed with the number of children. The upstream
   fixture `DTAExtractor_1_input.mzML` declares `softwareList count="5"` with
   four entries and `dataProcessingList count="3"` with one, and C++ loads it —
   so the port could not read its own reference data. Reading now ignores the
   declared value; writing still emits the true count.
2. **`DTAFile::store` discards what DTA cannot represent.** The native writer
   refused to drop spectrum metadata, extra precursors, identifications or
   precursor acquisition metadata. Source writes the precursor mass and charge
   plus the peaks and silently ignores the rest, and uses the legacy proton mass
   (`(mz - 1.0) * charge + 1.0`), not the exact one. `dta::WriteOptions::source()`
   selects that behavior; the checked default is unchanged for library callers.

**A reversed retention-time range is not matched yet.** The C++ tool passes its
`-rt` bounds to `DRange<1>(rt_l, rt_u)`, whose constructor swaps reversed bounds
(`DTAExtractor.cpp:144` at topp 174b576, `DIntervalBase.h:85-90` at core
bc9cc12), so `DTAExtractor -test -rt 70:50` on the product SDK exits 0 and
writes `DTA_RT60.0.dta` (verifier probe, not in a manifest). This port's tool
compares retention times against the bounds as `parse_range` returns them and
writes nothing for that command line. The m/z bounds are compared as given in
both. The fix belongs to `src/cli/tools/dta_extractor.rs` and is recorded for
that tool's owner; no test covers it yet.

## MzMLSplitter

`src/bin/MzMLSplitter.rs` is the second tool, and the first whose output is
compared **canonically rather than byte for byte**. Both upstream
`TOPP_MzMLSplitter_*` invocations agree with the retained C++ parts on record
counts, native identifiers, MS levels, retention times and peak values. The two
mzML writers differ in serialisation detail, and the upstream test itself uses
FuzzyDiff rather than a byte comparison, so byte equality is not the contract.
Source conventions preserved: the part count derived from a file size in
KB/MB/GB base 1024, the remainder spread over the parts still to come, zero
padding to the width of the part count, and the refusal of `no_chrom` together
with `no_spec`.

## MapNormalizer and SpectraFilterWindowMower

`MapNormalizer` scales MS1 peak intensities to a percentage of the run maximum;
its upstream test reproduces the retained C++ output, and the most intense MS1
peak lands on 100. Higher MS levels are untouched and the source's commented-out
chromatogram branch is not ported.

`SpectraFilterWindowMower` is the first tool with an **algorithm subsection**,
the shape most remaining TOPP tools take. `Tool::subsection_defaults` ports
`getSubsectionDefaults_`: the algorithm's parameter tree is merged beneath the
tool's own defaults, so `-write_ini` emits it and an INI or a `-algorithm:name`
option can override it. Its upstream test matches the retained output, and three
further tests cover the subsection itself — that the C++ `WindowMower` defaults
(`windowsize` 50, `peakcount` 2, `movetype` slide) reach the INI, that changing
`peakcount` changes the result, and that a value violating a registered
restriction is **rejected with exit 6**, as the source's strict update does
(`TOPPBase.cpp:339-342`). An earlier version of this port ignored such a value in
favour of the default; the C++ product SDK exits 6 for the same INI.

One ordering detail the subsection support forced: a section description cannot
be set on a section that holds no entries, so the algorithm subsection
descriptions are applied after the subsection defaults are inserted.

`BaselineFilter` removes the baseline by morphological filtering and reproduces
its retained output with zero difference across all 132 intensities. It takes
the filter's three parameters as ordinary options rather than a subsection,
exactly as the source does, and keeps the source's two refusals: a run holding
only chromatograms, and spectra that are not sorted by m/z.

None of the five tools adds the processing record the source attaches to its
output yet; the context services for it exist, and the retrofit is a separate
wave-2 step.

The second DTA finding above is the first concrete instance of the port's
"checked boundaries" convention blocking C++ parity. The resolution pattern —
keep the guard as the library default, add an explicit source-behavior option,
and have the tool opt in — is the one to apply as further tools meet their own
guards.

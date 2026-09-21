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
| `APPLICATIONS/TOPPBase_defs.h` (`Citation`, `cite_openms`) | `src/cli/defs.rs` |
| `enableLogging_`, `writeLog*_`, `writeDebug_`, the destructor | `src/cli/logging.rs` |
| `APPLICATIONS/ToolHandler.h`, `.cpp` | `src/cli/tool_handler.rs` |
| core `FORMAT/ParamCTDFile.h`, `.cpp` | `src/cli/param_ctd.rs` |
| core `FORMAT/ToolDescriptionFile.h`, `HANDLERS/ToolDescriptionHandler.h` (read side) | `src/cli/tool_description_file.rs` |
| core `DATASTRUCTURES/ToolDescription.h`, `ToolInfo.h` | `src/data_structures/tool_description.rs`, `tool_info.rs` |

The pinned product manifest the Release build installs,
`share/openms4/tools/topp.tools.tsv`, is `resources/tools/topp.tools.tsv`
(sha256 `5a90f7c1…`, byte-identical to the Release install).

## API mapping

| Source | Native |
| --- | --- |
| `TOPPBase(name, description, official, citations, toolhandler_test)` | `trait Tool` with `NAME`, `DESCRIPTION`, `CITATIONS`; `official` and `toolhandler_test` are stored and never read by the source, so they are not ported |
| `~TOPPBase` (removes an empty log file) | `ToolLog::finish`, at the end of `run_with` |
| `main(argc, argv)` | `cli::run::<T>()`, or `run_with::<T>(args, out, err)` |
| `registerOptionsAndFlags_` | `Tool::register(&mut ToolSpec)` |
| `main_` | `Tool::run(&ToolContext)`, or `Tool::run_io(&ToolContext, out, err)` for tools that write to the streams; both return `ToolResult`: `Ok` with the exit code `main_` returns, or a `ToolError` where `main_` throws (see *Run time*) |
| `getSubsectionDefaults_(section)` | `Tool::subsection_defaults` |
| `getSubsectionDefaults_()`, `getDefaultParameters_` | internal; `ToolSpec::to_param` is the parameter part |
| `getToolUserDefaults_` | internal: `<user directory>/<Tool>.ini` updates the defaults leniently, with the source's messages; see *Per-user defaults* |
| `parseCommandLine_` | internal |
| `handleWriteCommands_` | internal; `-write_ini` with the source's ISO-8859-1 declaration (`paramxml::WriteOptions::source`), `-write_ctd` through `ParamCtdFile`, CWL and JSON refused as the Release build refuses them; see *Tool descriptions* |
| `fileParamValidityCheck_` (both overloads) | internal; input formats through `FileHandler::get_type`, output formats through `type_by_file_name` |
| `checkIfIniParametersAreApplicable_` | internal |
| `checkParam_` | internal; its three calls after the strict update, with the source's warnings |
| `TOPPBase::ExitCodes` | `ExitCode`, same 15 variants and discriminants |
| `ParameterInformation`, `ParameterTypes` | `ParameterInformation`, `ParameterType` |
| `paramEntryToParameterInformation_`, `getParamArgument_`, `paramToParameterInformation_` | `ParameterInformation::from_param_entry` |
| `register{String,Int,Double}Option_`, `registerFlag_` | `ToolSpec::register_*` |
| `registerInputFile_`, `registerOutputFile_`, `registerOutputPrefix_`, `registerOutputDir_` | `ToolSpec::register_*` |
| `register{String,Int,Double}List_`, `register{Input,Output}FileList_` | `ToolSpec::register_*_list` |
| `registerFullParam_`, `registerParamSubsectionsAsTOPPSubsections_`, `registerTOPPSubsection_` | `ToolSpec::register_full_param`, `ToolSpec::topp_subsections` |
| `setValidStrings_`, `setValidFormats_`, `setMin/MaxInt_`, `setMin/MaxFloat_` | `ToolSpec::set_*`; `setValidFormats_`'s `force_OpenMS_format` not ported; `ParameterInformation::accepted_formats` reads a file parameter's formats wherever they are kept |
| `registerSubsection_`, `addText_`, `addEmptyLine_` | `ToolSpec::register_subsection`, `add_text`, `add_empty_line` |
| `findEntry_`, `getParameterByName_` | `ToolSpec::find` |
| `get{String,Int,Double}Option_`, `get{String,Int,Double}List_`, `getFlag_`, `getParam_` | `ToolContext::{string,int,double,string_list,int_list,double_list,flag,param}` |
| `getOutputDirOption` | `ToolContext::output_dir`, which creates the directory |
| `getParamAs*_`, `getSubsection_` | internal |
| `getIniLocation_`, `toolName_` | `ToolContext::ini_location`, `tool_name` |
| `test_mode_`, `debug_level_`, `log_type_` | `ToolContext::test_mode`, `debug_level`, `progress_log_type` |
| `setMaxNumberOfThreads` | `ToolContext::thread_policy` |
| `UniqueIdGenerator::setSeed(19991231235959)` under `-test` | `ToolContext::unique_id_generator`, `TEST_MODE_UNIQUE_ID_SEED` |
| `getProcessingInfo_` (both overloads) | `ToolContext::processing_info` |
| `addDataProcessing_` (`ConsensusMap`, `FeatureMap`, `PeakMap`) | `ToolContext::add_data_processing`, trait `AddDataProcessing` |
| `inputFileReadable_`, `outputFileWritable_` | `cli::input_file_readable`, `cli::output_file_writable` |
| `parseRange_(text, double&, double&)` | `parse_range` |
| `parseRange_(text, Int&, Int&)` | `parse_range_int` |
| `version_`, `verboseVersion_` | `product_versions(registry, name)`, `cli::verbose_version`; the registry's product version, else the core version |
| `getToolPrefix`, `getIniLocation_` with `-instance` | internal; the instance number selects `<tool>:<n>:` |
| `printUsage_` | internal (`src/cli/usage.rs`), reached by `--help` / `--helphelp` and command-line errors; writes to the error stream, shaped and coloured for a console by an executable (*On a console*); `cli::print_usage` for a tool body that prints its own |
| `writeLogInfo_`, `writeLogWarn_`, `writeLogError_`, `writeDebug_` (both), `enableLogging_` | `ToolContext::write_log_info`, `write_log_warn`, `write_log_error`, `write_debug`, `write_debug_param`; `ToolLog`; see *Log file and debug levels*. A tool body on the `-threads` pool records its lines in the crate's `PoolLines`, which logs each at once and writes them to the console when the pool returns. `cli::log_error` and `cli::log_warning` write one record of the error and warning log streams, red and yellow on a terminal |
| the run-phase and initialisation catch blocks (`TOPPBase.cpp:428-513`) as a tool body reaches them | `ToolError`: `Error` (mapped by the lifecycle), `Caught { code, message }` (a catch block's own text, `ToolError::unexpected` for the `BaseException` arm, `ToolError::file_not_found`), `Escaped { what }` (a `std::exception`, the initialisation catch) |
| `getDocumentationURL` | internal; the release URL of the core version, used by the usage text and the CTD |
| `Citation`, `Citation::toString`, `cite_openms` | `Citation`, `Citation::to_source_string` (also `Display`), `CITE_OPENMS`; tool citations as `Tool::CITATIONS` |
| `TOPPBase_defs.h` exceptions `UnregisteredParameter`, `WrongParameterType`, `RequiredParameterNotGiven` | not types: the required-value check reports `RequiredParameterNotGiven`'s message and exit 7; the accessors return `Error::InvalidValue` for the other two (see *Native differences*) |
| `ToolHandler::getTOPPToolList`, `getToolVersion`, `findExecutable`, `getTypes`, `getCategory`, `getInternalToolsPath` | `ToolHandler::get_topp_tool_list`, `get_tool_version`, `find_executable`, `get_types`, `get_category`, `get_internal_tools_path`; `ToolListType` |
| `ToolHandler` private `getInternalTools_`, `getInternalToolConfigFiles_`, anonymous `prefixes`, `packageTools` | internal; `ToolHandler::internal_tool_config_files`, `package_tools`; the environment as `ToolRegistrySources` |
| `ParamCTDFile::store`, `writeCTDToStream`, private `escapeXML`, `replace` | `ParamCtdFile::store`, `write_ctd_to_stream`, `to_ctd_string`, `escape_xml`, `replace` |
| `ToolDescriptionFile::load` | `ToolDescriptionFile::load`, with the handler's non-fatal messages in `LoadedToolDescriptions::diagnostics`; `store` throws `NotImplemented` in the source and is not ported |

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
  `../oracle/topp-cli-lifecycle/manifest.json`, as in
  `../oracle/topp-cli-lifecycle/ini_read_failures/manifest.json` where marked †,
  or as in `../oracle/topp-cli-lifecycle/cli2/manifest.json` where marked ‡;
* **tier 4, derived**: read from the cited source lines, with no executed case;
* **tier 4, native**: a native mapping. The Rust `Error` enum is coarser than the
  source's exceptions, so there is no C++ case to execute; the row says what the
  source does instead where the two differ.

One row has no source phase at all and runs before all of them; *The x86_64 FMA
requirement* below says why.

| Phase | Condition | Exit | Source | Test | Evidence |
| --- | --- | --- | --- | --- | --- |
| CPU check | this build requires the x86-64 FMA3 instructions and the processor does not have them | 12 | — | `only_a_requiring_build_on_a_processor_without_fma_is_refused` and the other cases in `src/system/cpu_features.rs` | tier 4, native: the C++ build has no such requirement, so there is nothing to execute |
| registration | tool registration fails | 6 | 505-508 | `an_initialisation_failure_is_illegal_parameters` | tier 4, derived |
| parse | more than `MAX_ARGUMENTS` tokens or `MAX_ARGUMENT_BYTES` bytes | 6 | — | `an_oversized_command_line_is_refused_before_parsing` | tier 4, native bound |
| parse | a flag followed by text | 6 | 186-191, 2336-2348 | `a_flag_followed_by_text_is_refused`, `upstream_flag_with_trailing_arguments` | tier 1: `flag_with_trailing_text` |
| parse | a number that does not convert, including `-section:name` values | 6 | 186-191, 2365-2407 | `a_non_numeric_integer_is_refused`, `a_non_numeric_subsection_integer_is_refused` | tier 1: `int_option_not_numeric`, `algorithm_peakcount_not_int` |
| parse | no arguments at all | 6 | 227-232 | `a_bare_invocation_is_refused`; also `tests/topp_dta_extractor.rs`, `tests/topp_baseline_filter.rs`, `tests/topp_map_normalizer.rs` | tier 1: `no_arguments` |
| parse | `--help`, `--helphelp`: usage on the error stream | 0 | 235-239, 631-888 | `usage_text_matches_the_cpp_text_for_every_ported_tool`; `usage_and_exit_codes_follow_the_source_contract` in `tests/topp_dta_extractor.rs` | tier 1: `help`; `help_<tool>`‡ and `helphelp_<tool>`‡ for the five tools, text byte for byte apart from the revision, see *Usage text* |
| parse | an unknown option, including an unknown `-section:name` | 6 | 242-247 | `an_unknown_option_is_refused`, `an_unknown_subsection_parameter_is_refused` | tier 1: `unknown_option`, `algorithm_bogus` |
| parse | trailing text | 6 | 250-255 | `trailing_text_is_refused` | tier 1: `trailing_text` |
| parse | trailing text after several options, listed in command-line order | 6 | 2436-2444 | `a_command_line_at_the_argument_bound_is_parsed` | tier 4, derived |
| write | `-write_ini` target not writable | 5 | 2609 | `write_ini_to_an_unwritable_path_is_refused` | tier 1: `write_ini_unwritable` |
| write | `-write_ini` | 0 | 2606-2641, 2097-2256 | `write_ini_ignores_command_line_values`, `write_ini_with_an_invalid_ini_value_keeps_the_default`, `write_ini_matches_the_cpp_file_for_every_ported_tool`, `upstream_write_ini_matches_the_retained_files` | tier 1: `write_ini_with_cli_value`, `write_ini_with_ini`, `write_ini_<tool>`‡ for the five tools; the retained `TOPPBase_test_write_ini_out.ini` and `TOPPBase_test_write_ini_subsec_out.ini` |
| startup | a tool manifest the registry refuses (a row without four fields, an unsafe path, a duplicate name) or cannot read; a `.ttd` registry that collides or does not parse | 6 | 141-163, 505-508 | `a_manifest_the_registry_refuses_ends_the_run_as_in_the_release_build`, `an_unreadable_manifest_or_directory_ends_the_run_as_in_the_release_build`, `the_internal_tool_registry_behaves_as_in_the_release_build`, `a_malformed_ttd_file_ends_the_run` in `tests/topp_cli_completion.rs` | tier 1 (Release build): `reg_*`, `ttd_*` in `../oracle/toppbase-completion` |
| write | `-write_ctd <dir>` | 0 | 2643-2649, 2551-2603 | `write_ctd_matches_the_release_build_for_the_ported_tools` and the other `write_ctd_*` cases in `tests/topp_cli_completion.rs` | tier 1 (Release build): `ctd_<tool>` byte for byte, `ctd_notest_*`, `ctd_with_ini`, `ctd_with_cli`, `ctd_existing`, `ctd_trailing_slash`, `ctd_cwd_*`, `ctd_relative_*` |
| write | a tool-description target that is not writable | 5 | 2566 | `write_ctd_target_failures_exit_as_in_the_release_build` | tier 1 (Release build): `ctd_missing_dir`, `ctd_ro_dir`, `write_json_missing_dir` |
| write | a `-write_ctd` target that is a directory | 12 | 2601, 510-514 | `write_ctd_target_failures_exit_as_in_the_release_build` | tier 1 (Release build): `ctd_target_is_directory` |
| write | `-write_cwl`, `-write_nested_cwl`, `-write_json`, `-write_nested_json` | 12 | 2651-2683; `ParamCWLFile.cpp:331`, `ParamJSONFile.cpp:326` | `cwl_and_json_writers_are_refused_as_the_release_build_refuses_them` | tier 1 (Release build): `write_cwl`, `write_nested_cwl`, `write_json`, `write_nested_json`, `write_cwl_cwd`, `write_nested_cwl_and_json`; the Release build is built without TDL, see *Tool descriptions* |
| INI | a per-user `<Tool>.ini` that does not parse / is a directory | 3 / 12 | 2276-2286 | `unusable_user_defaults_end_the_run_as_in_the_release_build` | tier 1 (Release build): `ud_malformed`, `ud_directory` |
| run | an `is_executable` input found neither as given nor on `PATH` | 14 | 1534-1549, 442-447 | `an_is_executable_input_is_resolved_on_path` | tier 4, derived |
| INI | `-ini` file missing | 1 | 296, 436-441 | `a_missing_ini_is_input_file_not_found` | tier 1: `ini_missing` |
| INI | `-ini` file not readable, before a run or with `-write_ini` | 2 | 296, 2630, 448-453 | `an_unreadable_ini_is_input_file_not_readable`, `a_fifo_ini_this_user_cannot_open_is_input_file_not_readable` | tier 1: `ini_unreadable`†, `write_ini_ini_unreadable`† (a regular file); `ini_fifo_denied`†, `write_ini_ini_fifo_denied`† (a FIFO with mode 000: not queried for readability before the load, whose open is refused) |
| INI | `-ini` file malformed | 3 | 296, 460-465 | `a_malformed_ini_is_input_file_corrupt` | tier 1: `ini_malformed` |
| INI | `-ini` names a directory, before a run or with `-write_ini` | 3 | 296, 2630, 460-465 | `an_ini_directory_is_input_file_corrupt` | tier 1: `ini_directory`†, `write_ini_ini_directory`† |
| INI | `-ini` names the character device `/dev/null`, before a run or with `-write_ini` | 3 | 296, 2630, 460-465 | `a_character_device_ini_is_input_file_corrupt` | tier 1: `ini_dev_null`†, `write_ini_ini_dev_null`† |
| INI | `-ini` exists but cannot be opened for another reason than a missing file or permissions (a Unix socket, `/dev/tty` without a controlling terminal), before a run or with `-write_ini` | 8 | 296, 2630, 495-499; `TextFile.cpp:44-47` | `a_socket_ini_is_an_unexpected_internal_error`, `a_terminal_ini_without_a_controlling_terminal_is_an_unexpected_internal_error` | tier 1: `ini_socket`†, `write_ini_ini_socket`†, `ini_tty`†, `write_ini_ini_tty`† |
| INI | `-ini` names a FIFO this user can open | as the INI it carries once a writer opens it; blocks until then | 296, 2630 | `an_ini_fifo_is_read_once_a_writer_opens_it` (exit 0 with `-write_ini`) | tier 4, native, a deliberate difference: the C++ tool opens the INI twice and blocks with a single writer (observation `write_ini_ini_fifo_single_writer`†); with no writer both block |
| INI | INI without a section for this tool | 0, warning | 1957-1966 | `an_ini_for_another_tool_warns_and_applies_defaults` | tier 1: `ini_foreign_section` |
| update | unknown parameter (INI item, `-instance`, top-level `common:` value) | 6 | 338-343 | `an_unknown_ini_item_is_refused`, `instance_on_the_command_line_is_refused`, `a_common_top_level_value_is_refused_as_in_the_source` | tier 1: `ini_unknown_item`, `instance_on_command_line`, `ini_common_top_level` |
| update | value outside its restrictions | 6 | 338-343 | `an_invalid_subsection_value_on_the_command_line_is_refused`, `an_invalid_subsection_value_in_an_ini_is_refused` | tier 1: `algorithm_movetype_sideways`, `ini_invalid_subsection_value` |
| update | value of another type | 6 | 338-343 | `an_ini_value_of_the_wrong_type_is_refused` | tier 1: `ini_value_type_mismatch` |
| update | INI written by another version | 0, notice on stdout | 355-366 | `an_ini_from_another_version_is_noted` | tier 1: `ini_version_mismatch` |
| validation | required value missing or empty | 7 | 1398-1406, 466-474 | `a_missing_required_output_is_missing_parameters` | tier 1: `missing_required_output` |
| validation | input missing / unreadable / empty | 1 / 2 / 4 | 1968-1995 | `a_missing_input_is_input_file_not_found`, `an_unreadable_input_is_input_file_not_readable`, `a_zero_byte_input_is_input_file_empty` | tier 1: `missing_input`, `unreadable_input`, `zero_byte_input` |
| validation | input format detected by name, then by content, and not accepted | 6 | 1575-1593 | `the_input_format_is_detected_by_name_then_content`; `usage_and_exit_codes_follow_the_source_contract` in `tests/topp_dta_extractor.rs` | tier 1: `in_txt_extension_mzml_content`‡, `in_dta_extension`‡ |
| validation | input format detected and accepted (no extension, an unknown one, or another case) | runs | 1575-1593 | `the_input_format_is_detected_by_name_then_content` | tier 1: `in_no_extension_mzml_content`‡, `in_unknown_extension_mzml_content`‡, `in_uppercase_extension`‡ (exit 0) |
| validation | input format undetermined | warning, then runs | 1580-1583 | `an_undetermined_input_format_only_warns` | tier 1: `in_unknown_name_and_content`‡ (the warning; the later load's exit code differs, see *Formats*) |
| validation | output not writable | 5 | 1997-2013 | `an_unwritable_output_is_cannot_write_output_file` | tier 1: `unwritable_output` |
| validation | output name of a known type the parameter does not accept | 6 | 1595-1608 | `a_wrong_output_extension_is_refused`, `the_output_format_is_checked_by_name_only` | tier 1: `output_wrong_extension`, `out_txt_extension`‡ |
| validation | output name of no known type, or of an accepted type in another case | runs | 1598-1601 | `the_output_format_is_checked_by_name_only` | tier 1: `out_unknown_extension`‡, `out_no_extension`‡, `out_uppercase_extension`‡ (exit 0) |
| run | `Error::Parse` from the tool | 3 | 460-465 | `a_corrupt_input_is_input_file_corrupt`, `run_phase_errors_map_like_the_source_inner_catch` | tier 1: `corrupt_input` |
| run | `Error::Io`, not found | 1 | 436-441 | `run_phase_errors_map_like_the_source_inner_catch` | tier 4, native: as `FileNotFound` |
| run | `Error::Io`, permission denied | 5 | 430-435 | `run_phase_errors_map_like_the_source_inner_catch` | tier 4, native: an I/O error does not say whether it read or wrote; a tool's inputs are checked for readability before its body runs, so a denial there is taken as a failed write. INI files never take this row |
| run | any other `Error::Io` | 8 | 495-499 | `run_phase_errors_map_like_the_source_inner_catch` | tier 4, native |
| run | `Error::InvalidValue`, `InvalidRange` | 6 | 475-480 | `run_phase_errors_map_like_the_source_inner_catch` | tier 4, native: as `InvalidParameter`; the source's own `InvalidValue` and `InvalidRange` exceptions reach its `BaseException` arm, exit 8 |
| run | `Error::MissingInformation` | 7 | 466-474 | `run_phase_errors_map_like_the_source_inner_catch` | tier 4, native: as `RequiredParameterNotGiven`; the source's `MissingInformation` exception reaches its `BaseException` arm, exit 8 |
| run | `Error::Unsupported`, `UnsortedData` | 11 | the tools' own `INCOMPATIBLE_INPUT_DATA` returns | `run_phase_errors_map_like_the_source_inner_catch` | tier 4, native |
| run | `ToolError::Caught` from a tool body: a source exception whose catch-block text the tool knows | the tool's code | 428-499 | `a_caught_exception_prints_no_closing_line_and_reaches_the_log`; the tools' `release_*` tests | tier 1 (Release build): `fi_out_is_directory_log`, `fd_matched_whitelist_log`, `pphr_centroided_log`, `ffc_profile_noforce_log`, `ffc_ms2_only_log` in `../oracle/topp-exception-exits` |
| run | `ToolError::Escaped` from a tool body: a `std::exception` | 12 | 510-513 | `an_escaped_exception_prints_no_closing_line_and_stays_out_of_the_log` | tier 1 (Release build): `fd_directory_log` in `../oracle/topp-exception-exits` |
| run | a tool body returns an exit code after its own message | that code, closing line printed | 413-424 | the tools' `release_*` tests | tier 1 (Release build): `fi_index_on_dta_log`, `fi_notype_log`, `pphr_empty_log`, `ms_*`, `dta_bad_*`, `bf_empty_log` in `../oracle/topp-exception-exits` |

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
`ini_fifo_denied` and `write_ini_ini_fifo_denied`.

An open that fails for any other reason is exit 8 with the source's
`Error: Unexpected internal error (IO error for file '<path>')`. The file exists
and passes `File::readable`, so `TextFile` throws `IOException`
(`TextFile.cpp:44-47`), which reaches the `BaseException` arm (495-499). The
oracle records this for a Unix socket (`ini_socket`, `write_ini_ini_socket`) and
for `/dev/tty` opened without a controlling terminal (`ini_tty`,
`write_ini_ini_tty`). The port tells such an open failure from a later read
failure by opening the file once more after the load has failed; a FIFO is not
opened again, because that open would wait for a writer, so its I/O failure is
taken as a read failure. A read failure after a successful open, and a document
that does not parse, stay exit 3. `/dev/null` reads as an empty document, exit
3 with an `Error: Unable to read file (` line on both paths, as the oracle's
`ini_dev_null` and `write_ini_ini_dev_null` do; the source's line also names the
file and xerces's message.

A FIFO this user can open is opened once and read like a file, so it waits for
a writer; the C++ tool also blocks on such a FIFO with no writer. With a writer
that opens the FIFO once, the two differ, deliberately: the source opens the INI
twice, once to look for compression (`XMLFile.cpp:141-147`) and once for xerces
(`:166`), so after its first open closes no writer is left and the second open
waits. The oracle observation `write_ini_ini_fifo_single_writer` records it: the
writer sends the whole INI, and the C++ tool blocks until a 10-second alarm ends
it (exit 142) without writing anything. This port reads what the one writer
sends (`an_ini_fifo_is_read_once_a_writer_opens_it`, tier 4). Reproducing the
second open would make the port hang on an input it can read.

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
when given, and never command-line values (`TOPPBase.cpp:2606-2641`). The file
declares ISO-8859-1, as `ParamXMLFile::store` does. The tree is
`getDefaultParameters_`'s: an input file carries the `input file` tag, an output
file `output file`, an output prefix `output prefix` and an output directory
`output dir`, and their formats become `*.<format>` restrictions, which the
writer emits as `supported_formats`; flags are typed `bool`; a parameter's own
tags other than `is_executable` on an input file are not written (2097-2256).
For each of the five ported tools the file matches the product SDK's line by
line with only `version` lines skipped, and entry for entry once decoded; the
class test's `TOPPBaseTest` and `TOPPBaseCmdParseSubsectionsTest` files match
the retained C++ files at `TEST_FILE_SIMILAR`'s tolerances.

**File checks** follow `inputFileReadable_` and `outputFileWritable_`: an input
must exist, be readable and, unless it is a directory, hold a byte; an output
must be writable, and an output prefix is probed as `<prefix>_0`. The
writability query never creates or deletes a file under the caller's name.

**Formats** follow `fileParamValidityCheck_` (1498-1614). An input file's format
is what `FileHandler::get_type` detects: the name first, then the content for a
name no type claims. A detected format the parameter accepts runs; one it does
not accept is `Invalid parameter: Input file '<path>' has invalid format
'<type>'. Valid formats are: '<formats>'.`, exit 6; an undetermined format only
warns, `Warning: Could not determine format of input file '<path>'!`. An output
file's format comes from its name alone: a known type the parameter does not
accept is `Invalid parameter: Invalid output file extension for file '<path>'.
Valid file extensions are: '<formats>'.`, exit 6, and a name no type claims is
accepted. Types compare with their names without regard to case. Input file
lists are checked file by file; output file lists and output prefixes are not
checked for format, and a parameter without formats is not checked at all.

**Usage text** is `printUsage_`'s, on the error stream for `--help` and
`--helphelp` as after a command-line error: the tool line, the documentation
URL, the version line `<product version> (OpenMS core <version>, revision
<short revision>)`, the OpenMS citation, the options with a description column
at six past the longest shown name and argument, `(default: '…')`,
`(valid: …)`, `(valid formats: …)` and `(min: …)`/`(max: …)` addons, and for a
tool with subsections either their summary (`--help`) or every subsection
parameter under its section description (`--helphelp`). Descriptions start with
an upper-case letter, and a line break in a description continues at the
description column (`IndentedStream`, `ConsoleUtils::breakString_`).

**On a console** (`src/cli/console.rs`). A tool executable (`cli::run`) lays
its usage text out for the console as the source does. The width is `COLUMNS`
read as `StringUtils::toInt32` reads it, or, when `COLUMNS` is unset, the
second field of `stty size` run through `/bin/sh` with the process's standard
input, whose complaint about a standard input that is not a terminal reaches
standard error as the source's does; the width is that number less one, and
anything unreadable or below 10 turns shaping off (`readConsoleSize_`). Each
item the source inserts into its `IndentedStream` is broken as
`breakString_` breaks it: from the current column, continuation lines at the
current indentation, a line that fills the width and ends in a word shorter
than four bytes gives that word to the next line, an item of more than ten
lines keeps eight, an indented `...` and its last line, and the source's
unsigned arithmetic is kept, so an option column wider than the console
stops shortening continuation lines and a line that ends exactly at the edge
is followed by an empty one (CPP-353). On a terminal standard error the usage
text carries the source's `Colorizer` codes (the tool name and `Usage:`
inverted, headings bright, the documentation URL underlined, required options
and the mandatory note green, defaults cyan, restrictions magenta, empty
addons included), every `writeLogError_`/`OPENMS_LOG_ERROR` line of the
lifecycle is red and every `writeLogWarn_`/`Param::update` line yellow
(`LogStreamBuf::distribute_`), and the run ends with `\x1b[0m` on every
terminal among standard output and error (`InitConsole`). A tool driven in
process through `run_with` writes to explicit streams, which are not a
console: no probing, no shaping, no colours, which is the text the source
writes when standard error is not a terminal and `COLUMNS` is unset. Evidence:
100 cases (every tool's `--help` and `--helphelp` at `COLUMNS` 20 to 120, 16
`COLUMNS` values, 4 failures), byte-identical to the Release build apart from
`FeatureFinderCentroided`'s missing citations (`tests/topp_cli_console.rs`),
and 15 cases in a
pseudo-terminal (`tests/data/topp_cli_console/tty_*`), whose usage text and
log lines are compared in the unit tests of `src/cli/usage.rs` and
`src/cli/console.rs`; the executables in a pseudo-terminal match 13 of them
byte for byte (`../oracle/toppbase-completion/compare_tty.py`), and the other
two are tool bodies: `FeatureFinderCentroided` lacks its citations, and
`BaselineFilter` prints no progress lines. (`FileInfo`'s own `Error: Can only
validate indices for mzML files`, the third until the tools wrote their
messages through the log streams, is red now, as in the Release build.) Not ported: the
Windows console-buffer width (`GetConsoleScreenBufferInfo`) and virtual
terminal mode, which the standard library cannot reach; there `COLUMNS` still
applies.

**Processing records (decision D4).** BaselineFilter, MapNormalizer and
SpectraFilterWindowMower attach `getProcessingInfo_` with `BASELINE_REDUCTION`,
`NORMALIZATION` (`Intensity normalization`) and `FILTERING` (`Data filtering`)
to every spectrum and chromatogram after processing, as their sources do.
MzMLSplitter calls the same code while each part still holds no spectra and no
chromatograms, so its parts carry no new record; that is reproduced (C++ issue
candidate). DTAExtractor writes DTA files, which hold no processing record.
Against the retained outputs `BaselineFilter_output.mzML`,
`MapNormalizer_output.mzML`, `SpectraFilterWindowMower_1_output.mzML` and
`MzMLSplitter_output_part1/2.mzML`, every record holder has the same records in
the same order, with the same software, version, actions, completion time and
metadata; outside `-test`, SpectraFilterWindowMower's record matches the product
SDK's (`window_mower_notest`‡) key for key. Two comparisons follow the C++ mzML
writer rather than the tool: it stores completion times to the minute
(`MzMLHandler.cpp:3947`), and a software name as the PSI-MS term found by the
name, the name plus ` software` or `TOPP ` plus the name (3763-3787), so
`SpectraFilterWindowMower` reloads from a C++ file as
`TOPP SpectraFilterWindowMower`. This port's mzML writer keeps the seconds and
the exact name.

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
parallel result is bit-identical to the serial one. `-threads` reaches a tool
through `ToolContext::thread_policy`, and a tool opts into a scoped rayon pool
of that size by wrapping its body in `ToolContext::in_thread_pool`, where the
source calls the process-wide `omp_set_num_threads` before `main_`. A positive
count `n` gives `n` workers; **zero and every negative count give every
available processor**, which is the source's own `if (num_threads <= 0)
num_threads = omp_get_num_procs()` and was established by executed C++
(`../oracle/tool-threads`). `OMP_NUM_THREADS` and `RAYON_NUM_THREADS` do not
size the body, exactly as the source's `omp_set_num_threads` overrides
`OMP_NUM_THREADS`; the registered default stays 1, so a run without `-threads`
is serial. Without the `parallel` feature the body runs on the calling thread.
The full contract, the API mapping and the native differences of the scoped
pool are in [TOPP_THREADS_SUPPORT](TOPP_THREADS_SUPPORT.md).

Opting in is per tool, because the `out`/`err` streams of `Tool::run_io` are
not `Send` and no central wiring in `run_with` was possible. The five wave-2
tools (`BaselineFilter`, `DTAExtractor`, `MapNormalizer`, `MzMLSplitter`,
`SpectraFilterWindowMower`) wrap their bodies. `PeakPickerHiRes` opens the same
pool but **around its picking call rather than around its body**
(`src/cli/tools/peak_picker_hi_res.rs:249`), and skips it altogether at one
worker: that is measured at 0.68 s cheaper on a gigabyte-scale run than putting
the whole body on a worker, and a body wrapper in `Tool::run` would in any case
never execute for this tool, which overrides `run_io`. It is the **sixth**
tool in `tests/topp_threads.rs::TOOLS` and in
`every_executable_runs_its_body_on_the_requested_pool`, with one stated
exception: it starts no worker at `-threads 1`, by design.

`FeatureFinderCentroided` takes `ctx.thread_policy()` into its algorithm without
wrapping its body, and `FileInfo` uses neither, so `-threads` sizes no pool for
those two. Nothing is computed wrongly — both are byte-identical at every thread
count — but they do not have the contract this section describes. The wave-4
benchmark sampled all eight executables on full-size data: two use their pool
(`PeakPickerHiRes` at CPU utilisation 2.26 and `FeatureFinderCentroided` at
4.72, through its algorithm), five build it and leave it idle at utilisation
1.00, and `FileInfo` builds none. Carried forward in
[the work packages](EARLY_TOPP_WORK_PACKAGES.md#wave-4-status).

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

**The x86_64 FMA requirement.** The C++ SDK is built for the baseline
architecture with no `-march`, so a C++ tool binary runs on any processor its
operating system runs on. This port's x86_64 builds do not: `.cargo/config.toml`
sets `-C target-feature=+fma`, because the ported glibc `exp`, `log` and `powf`
are built out of `f64::mul_add` and a baseline x86-64 target has to call out of
line for every one of them. An x86_64 binary built in this checkout therefore
needs an FMA3-capable processor — Intel Haswell, AMD Piledriver or newer — and a
tool refuses to start on anything older, printing the requirement and the
command that rebuilds without the flag and exiting `INTERNAL_ERROR` (12) rather
than dying on `SIGILL` at whatever arithmetic came first.

The check is the **first statement of `cli::run`** (`src/cli.rs`), which is what
every tool executable calls, ahead of reading the command line and ahead of every
source phase; everything after it is in `run_from_environment`, which is
`#[inline(never)]`, so no instruction of the body can be hoisted above it. That
placement is the shipped one and it was arrived at by measurement: with the check
one level lower, in `run_with`, `cli::run`'s own prologue emitted `vpxor`,
`vmovdqu` and `vmovups %ymm0` while building the argument vector, all of them
ahead of the check. In the release `FileInfo` on kim, `cli::run` inlines into the
tool's `main` and the check is the first call, with no VEX or SSE instruction
before it. `cli::run_with` holds a second copy of the same check, for a caller
that drives a tool in process; no tool binary reaches the guard through that
copy. A build without the flag compiles both away. See `docs/FMA_BUILD_FLAG.md`
section 8 for what this does and does not guarantee, and
`src/system/cpu_features.rs`; the flag changes no result, only speed.

**Usage text.** Driven in process, the layout is the one the source writes
when standard error is not a terminal and `COLUMNS` is unset, as in the
oracles; the executables shape and colour it for a console (above). For the
seven ported tools that declare their citations, `--help` and `--helphelp`
match the Release
build byte for byte (`help_<tool>`, `helphelp_<tool>` in
`../oracle/toppbase-completion`; the Release build's revision is the port's,
`bc9cc12`), apart from the line its console-width probe prints first when
standard input is not a terminal (`stty: 'standard input': Inappropriate ioctl
for device` from GNU `stty`), which only an executable prints. Tool citations
follow the OpenMS citation under
`To cite <tool>:`, and a tool no manifest lists gets `Common UTIL options:`.
`FeatureFinderCentroided` does not yet
declare its two citations (`FeatureFinderCentroided.cpp:124-136`), so its usage
and its CTD lack them. An earlier version of this port printed `--help` usage
on the output stream; the oracle's `help` case writes nothing there.

**Validation order.** The source checks each option lazily when `main_` reads it;
this port checks every registered option after the update, in registration
order. Exit codes agree whenever a tool reads all its options, which every ported
tool does.

**Formats.** When `FileHandler::get_type` fails, for example on a directory with
an unknown name (for which the source reports an unknown type; an open
`src/format/file_handler.rs` follow-up), the format counts as undetermined and
only warns. After that warning the source's load of a file whose content no
type claims throws `ParseError`, exit 3 ("type: unknown is not allowed for
loading an experiment", oracle `in_unknown_name_and_content`, and with
`-log` oracle `log_undetermined_format`, whose warning and error both reach the
log in both builds); this port's `FileHandler::load_experiment` returns
`Error::InvalidValue`, exit 6. That difference belongs to the loader, not to
the check.

**Input checks on a directory.** On the Linux Release build a directory given
as an input file whose name no type claims does not warn: `FileHandler::getType`
sniffs its content, libstdc++ throws `std::ios_base::failure` reading a
directory, and the initialisation catch ends the run with `Unable to initialize
or run <tool>: basic_filebuf::underflow error reading the file: Is a
directory`, exit 12, before the tool body runs (oracle `fi_undetermined_log` of
`../oracle/topp-exception-exits`). The macOS product SDK reads the directory as
an empty file, warns, and lets the tool decide, which is what this port does
(FileInfo exits 10). Recorded as an open difference
(`a_directory_input_is_a_recorded_difference_from_the_linux_release_build`),
not closed: the Release build's answer is a libstdc++ artefact, and FuzzyDiff's
reproduction of the same message for its own read (`ToolError::escaped`) is a
tool-level choice. `setValidFormats_`'s
`force_OpenMS_format` check of registered format names is not ported.

**`-write_ini`** is compared with the C++ files by number, not by text: this
writer prints `3` where the source prints `3.0`. The ISO-8859-1 declaration is
kept consistent with the bytes (`paramxml::OutputEncoding::Latin1`): the source
copies UTF-8 bytes under that declaration, so its non-ASCII text reads back as
different characters, while this writer writes ISO-8859-1 bytes and character
references. The five ported tools' files are ASCII.

**Tool descriptions.** `-write_ctd <dir>` writes `<dir>/<tool><type>.ctd`
for each `-type` the registry lists, once without a type for every ported tool:
the defaults — with the per-user defaults, never `-ini` or command-line values,
and the same with or without `-test` — and the tool's product version, name,
documentation URL, registry category, description and citation DOIs, the
OpenMS DOI first. The writer is byte-identical to the Release build for the
seven ported tools that declare their citations and to the retained
`ParamCTDFile_test_writeCTDToStream.ctd`. It reproduces three source defects a
CTD consumer sees, each executed through the Release library by
`../oracle/toppbase-completion/ctd_driver.cpp`: the replacement loop skips the
character after each replacement, so `&&` is written `&amp;&` and a second
consecutive line break stays raw; a tab is written `&amp;#x9;`; the `<tool>`
attributes, description and citations are not escaped (C++ issue candidates).
An empty directory is the current directory.

One difference cannot be observed yet: for a tool with `-type` values, the
source recomputes `getDefaultParameters_()` for each type after setting `type`
on its command-line parameters (`handleWriteCommands_`), so a subsection whose
defaults depend on the type differs per CTD; this writer computes the defaults
once and sets only `<location>type`. No ported tool registers types, so every
CTD written today is the Release build's; the first ported tool with types
must recompute them per type.

The four CWL and JSON writers need TDL, which the Release build (and the
source's default `ENABLE_TDL=OFF`) does not have: each checks its target, then
ends with the Release build's line `Unable to initialize or run <tool>: TDL
support is not available. Rebuild with -DENABLE_TDL=ON to enable this
feature.` and exit 12. The Release build opens the target before it fails and
so leaves an empty file behind, emptying one that existed (oracle
`write_cwl_existing`); this port leaves the target untouched, a deliberate
difference, because emptying a user's file is a loss.

**Tool registry.** The product version, the `TOPP`/`UTIL` heading of the
common options, the CTD category and the `-type` values come from
`ToolHandler`, which reads the package manifests of every
`OPENMS_TOOL_PREFIX_PATH` entry and of the executable's prefix, with the
source's validation and messages. A Rust tool is not installed by the source's
CMake rules, so the pinned product manifest is compiled in and read **in
place of the executable's own prefix when that prefix has no
`share/openms4/tools`**. The source reads a manifest reached twice once, and
knows it by its canonical path (`ToolHandler.cpp:90-104`): the Release build
runs with its own prefix in `OPENMS_TOOL_PREFIX_PATH`, once, twice, or through
a symbolic link (oracles `reg_install_prefix_twice`, `reg_install_prefix_once`,
`reg_symlink_to_install`). The built-in manifest has no path, so it is known by
its name and bytes: a prefix supplying `share/openms4/tools/topp.tools.tsv`
byte for byte is that manifest reached again, and the built-in copy is not
read a second time. So a port tool runs with `OPENMS_TOOL_PREFIX_PATH` naming
an OpenMS4 installation, which before this exited 6 for every tool, even
`--help` (finding F6 of `.planning/VERIFY-phase3w1.md`). Any other manifest that
lists a product tool again — another name, other bytes such as another product
version — is a duplicate, exactly as for two C++ installations (oracle
`reg_dup_help`). One recorded difference follows from the identity: the
Release build refuses a byte-identical copy of its manifest under another
prefix (oracle `reg_copy_of_install_manifest`, exit 6, because its own
installation's file is a second path), where the port reads that copy as the
built-in manifest. Manifests are read ordered by file name, where the
source takes the file system's order; which of two duplicate rows is reported
can differ. The source also consults the prefix of the loaded `libOpenMS_CLI`
through `dladdr`; the crate is linked into the executable, so that adds
nothing. The internal-tool (`.ttd`) directory is the shared-data directory's
`TOOLS/INTERNAL` when a shared-data directory resolves; the source's
`getOpenMSDataPath` throws when none does, which would stop every tool, while
this port's tools need no shared-data tree, so a missing one only means no
internal tools. `OPENMS_TTD_INTERNAL_PATH` is read as in the source. A `.ttd`
file that is not well formed ends the run with exit 6 as in the source, with
this reader's wording of the reason where the source prints Xerces's.

A tool no manifest lists reports the core version (as `ToolManifest_test.cpp`
asserts) and `Common UTIL options:`. Its verbose version line is the core
version and `, Revision: <short revision>`; the source puts the library's
build time between them, which this crate does not record.

**Per-user defaults.** `<user directory>/<Tool>.ini` — `OPENMS_HOME_PATH`, the
`home_dir` of `OpenMS.ini`, or the home directory — updates the defaults
leniently whenever they are built (a run, `-write_ini`, `-write_ctd`; not
usage), with `Param::update`'s verbose messages. An unreadable file is
skipped; one that does not parse ends the run with 3; a directory by that name
with 12 and libstdc++'s message, as on the Linux Release build. A run therefore
depends on the account, as in the source.

**Log file and debug levels.** `-log <file>` (from the command line, or from
the INI once it is merged) receives the source's lines: every info, warning and
error line the lifecycle writes through `writeLog*_`, and from debug level 1
the source's `writeDebug_` lines, including the accessors' `Value of … option`
lines and, from level 2, the parameter dumps and file checks. A line is
`<local time> <tool>:<instance>:: <text>`. At debug level 1 or more, opening
the file prints `Writing to '<file>'` on standard output. As in the Release
build, debug text never reaches the console (`OPENMS_LOG_DEBUG` is compiled
out there), and an empty log file named by `-log` is removed at the end. Two
orders differ from the source: the file checks run before the tool body, where
the source checks each file when the tool first reads it, so at level 2 the
`Checking … file` lines come before the body's `Value of` lines; and the
`Value of` lines follow the order the ported tool reads its options, which
for BaselineFilter is not the source's. A source debug line that carries a C++
source file and line (`Error occurred in line …`) has no counterpart.

The tool bodies' own lines reach the log file as the source's do: every
`writeLogError_`, `writeLogWarn_` and `writeLogInfo_` line of a ported tool
(FileInfo's two refusals, PeakPickerHiRes's sortedness errors, MzMLSplitter's
refusals and its report, DTAExtractor's `Invalid boundary`, BaselineFilter's
warnings), the catch-block line of a `ToolError::Caught`, and the tools'
`writeDebug_` lines (FuzzyDiff's two list lines at level 1, FileInfo's
detected type at 2, DTAExtractor's ranges and levels at 1, the parameter dumps
of PeakPickerHiRes, FeatureFinderCentroided and SpectraFilterWindowMower at 3).
`OPENMS_LOG_WARN` and `OPENMS_LOG_ERROR` lines of a body (PeakPickerHiRes's and
BaselineFilter's empty-input warnings, FeatureFinderCentroided's ion-mobility
refusal) and the initialisation catch's `Unable to initialize or run` line do
not, as in the source. Evidence: 28 Release runs, `../oracle/topp-exception-exits`,
retained under `tests/data/topp_exception_exits` and replayed through
`tests/support/release_runs.rs`. A tool body on the `-threads` pool
(BaselineFilter, DTAExtractor, MzMLSplitter) logs each line at once and puts it
on the console when the pool returns, before the closing line or the report of
a failure: the pool's worker cannot write to `run_io`'s streams, and only the
calling thread knows whether they are a terminal to colour.

**Run time.** After the tool body returns, whatever its exit code, the source's
closing line `<tool> took <wall> (wall), <cpu> (CPU), <system> (system), <user>
(user); Peak Memory Usage: <n> MB.` goes to standard output. A component this
port's `StopWatch` cannot read on the platform prints as `n/a`, and the memory
part appears only where the platform reports a peak (Linux). A body that ends
where the source's `main_` throws prints no closing line, as the source's
exception unwinds past it to a catch block. The distinction is in the type a
body returns, `ToolResult`: `Ok(code)` is a `main_` that returned `code`, even
an error code after the tool's own message (FileInfo's `-i` refusal, exit 6,
closing line printed); a `ToolError` is a `main_` that threw — `Error` for this
crate's errors, which the run-phase mapping reports, `Caught` for a catch
block's text the tool knows exactly (the `BaseException` arm's `Error:
Unexpected internal error (<what>)` for PeakPickerHiRes's centroided refusal,
FileInfo's unwritable report file, FuzzyDiff's malformed whitelist and the
FeatureFinderCentroided refusals; `FileNotFound`, `FileEmpty` and
`UnableToCreateFile` texts where the source throws those), written through
`writeLogError_`, and `Escaped` for a `std::exception`, which only the
initialisation catch handles (`Unable to initialize or run <tool>: <what>`,
exit 12, on the error stream only). Before the framework printed the closing
line, the tools stood in for thrown exceptions by writing the catch text and
returning a code, which after it printed the line where the Release build
prints none (finding F0 of `.planning/VERIFY-phase3w1.md`); the five earlier
tools did the opposite and raised errors where the source returns a code.

**Update diagnostics** of the strict update, the lenient `-write_ini` update and
the per-user defaults reach standard error at the end of the run: the source
writes them to its warning log stream without flushing it, so the Release build
prints them after everything else, after `Parameters passed to '<tool>' are
invalid. …` on a failed strict update.

**Not ported:** JSON INI files (`ParamJSONFile::load`, which the Release
build reads without TDL: a `.json` `-ini` is read as XML here and fails with 3,
where the Release build applies it, oracle `json_ini_*`; it needs `serde_json`
behind the `paramxml` feature), `INIUpdater`, `SearchEngineBase`,
`MapAlignerBase` and `TOPPExternalToolBase`.

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
makes it fail with `Error::InvalidValue`, exit 6. The source's `DataValue`
records any double: the product SDK's SpectraFilterWindowMower accepts
`-algorithm:windowsize inf` and exits 0 with a `parameter:
algorithm:windowsize` record of type `xsd:double` and value `inf` (oracle
`window_mower_notest_windowsize_inf`‡). The D4 retrofit keeps this documented
difference rather than recording a non-finite value in another type: a string
would change the record's type, which a decoded comparison with the C++ output
would report anyway. Among the ported tools the difference is not reached
through the record. BaselineFilter's `-struc_elem_length inf` fails the range
check first, exit 6 in both (`baseline_filter_notest_length_inf`‡); this
port's window mower refuses a non-finite window size before the record is
built, exit 6 where the C++ tool exits 0; MapNormalizer and MzMLSplitter have no
floating-point parameter.

**Diagnostics.** `Param::update_with_options` decides the update and applies it;
its report is worded natively, so the CLI derives the source wording from the
same entries. A missing file's parenthetical reads `does not exist`, where the
source says `could not be found`, because an existing assertion fixes the phrase.

**Errors** are typed `Result` values mapped at the boundary. A registration
error is exit 6 (the source's initialisation catch) rather than 12. A tool that
reads a parameter it did not register, or with the wrong accessor, gets
`Error::InvalidValue`, which a body that lets it escape ends with exit 6 where
the source's `UnregisteredParameter` and `WrongParameterType` reach
`INTERNAL_ERROR` (12); both are programming errors, not reachable from a
command line. `Error::InvalidValue` and `InvalidRange` from a tool body exit 6
where the source's own `InvalidValue`/`InvalidRange` exceptions exit 8. The
source's explicit `return ILLEGAL_PARAMETERS` after a message is now an exit
code the tool returns (MzMLSplitter's option checks and DTAExtractor's caught
`ConversionError` were the ones that raised `InvalidValue`); the crate-wide
arm still maps `InvalidValue` to 6, because other callers, among them the
parameter readers and MapNormalizer's non-positive scale, rely on that code,
and moving it is a separate decision.
`parse_range` leaves both bounds unchanged on error. A command line is bounded
by `MAX_ARGUMENTS` and `MAX_ARGUMENT_BYTES` before parsing, and parsed in time
linear in its length: the text left after each option is gathered in reverse and
ordered once, where the source inserts each chunk at the front of its list
(`TOPPBase.cpp:2436-2439`). Number conversion
follows `std::from_chars`, so hexadecimal floats are rejected even though the
oracle's libc++ `strtod` fallback accepts them. `unique_id_generator` returns an
independent generator per call instead of seeding a process-wide singleton.

## Checked boundaries and evidence

`tests/topp_cli_lifecycle.rs` holds 75 cases, `tests/topp_cli_completion.rs`
49, `tests/param_ctd.rs` 10, `tests/tool_description.rs` 7 and
`tests/parameter_information.rs` 3, none ignored.

* **Release build (tier 1 executed differential).**
  `../oracle/toppbase-completion/cases.sh` ran 124 cases of the C++ Release
  build at the port's pins (core `bc9cc12`, cli `c19e494`, topp `174b576`) on
  ibminode06, each in a clean environment with its own `HOME`,
  `OPENMS_HOME_PATH` and working directory; 97 of them are retained under
  `tests/data/topp_cli_completion/<case>/` (argv, environment, exit code,
  streams, log file, file tree) with `fixtures.sha256.json`, and the CTDs under
  `tests/data/topp_cli_lifecycle/release/`. `ctd_driver.cpp` wrote five CTDs
  through the Release `libOpenMS` (two identical runs), retained under
  `tests/data/param_ctd/oracle/`. `../oracle/topp-exception-exits/cases.sh`
  ran 28 more on the same build, twice with identical results: how each ported
  tool ends a returned and a thrown `main_`, and what reaches `-log`
  (retained under `tests/data/topp_exception_exits/`, provenance
  `tests/data/topp_exception_exits_provenance.json`).
* **Upstream class tests (tier 3).** `ToolHandler_test.cpp`,
  `ToolManifest_test.cpp`, `ParameterInformation_test.cpp` and the `-log` and
  `Citation::toString` sections of `TOPPBase_test.cpp` (cli `c19e494`);
  `ParamCTDFile_test.cpp` (with its retained output, tier 1),
  `ToolDescriptionFile_test.cpp` and `ToolDescription_test.cpp` (core
  `bc9cc12`). The command-line part of `getIniLocation_` (`-instance 5`) is
  transcribed through the log file, whose lines carry the location. The
  `getStringOption_` INI cases (`instance1`, `instance5`, `toolcommon`,
  `common`) read what the source's non-atomic strict update left behind
  after it failed; this port's update is atomic and the run ends before any
  tool reads a value, so they have no observable counterpart. The runs' exit
  codes and messages are covered by the oracle's `instance_on_command_line`
  and `ini_common_top_level`.

* **Oracle cases (tier 1 executed differential).** `../oracle/topp-cli-lifecycle/run.sh`
  runs 38 cases of the C++ product SDK (core 4fdec46, Debug, AppleClang 21) in a
  clean environment with a per-case directory; `manifest.json` records the build,
  binary hashes, inputs, exit codes, diagnostics and output hashes. Two executions
  agreed on all 37 shared cases. The subsection override is compared against the
  C++ output `tests/data/topp_cli_lifecycle/swm_algorithm_peakcount_1.mzML`; no
  case reaches a Debug-only precondition. `ini_read_failures.sh` in the same
  directory adds fifteen cases on the same binaries, recorded in
  `ini_read_failures/manifest.json`: an INI written by the tool, a readable
  control run with it, that INI unreadable before a run and with `-write_ini`,
  a directory as `-ini` on both paths, the character device `/dev/null` as
  `-ini` on both paths, a FIFO with mode 000 as `-ini` on both paths, a Unix
  socket as `-ini` on both paths, `/dev/tty` opened from a new session (without
  a controlling terminal) on both paths, and the observation
  `write_ini_ini_fifo_single_writer`. Each time cases were added, first the
  device cases, then the FIFO cases, then the socket, terminal and single-writer
  cases, the earlier cases ran again and matched every previous run in exit
  code, stderr and outputs; only timing figures in the control run's stdout
  changed. A FIFO with no writer is not a case, because the C++ tool blocks
  opening it. `cli2_cases.sh` adds 34 cases on the five ported tools' binaries,
  recorded in `cli2/manifest.json`: `-test -write_ini` of every tool,
  `--help` and `--helphelp` of every tool, ten input and output format cases,
  the registered BaselineFilter, MapNormalizer, SpectraFilterWindowMower and
  MzMLSplitter invocations, two runs without `-test` and two with a non-finite
  floating-point parameter. A second execution (`cli2_rerun_manifest.json`)
  agreed on every exit code and stderr; the outputs differed only in the
  completion times of the three non-test mzML files, and stdout only in timing
  figures. The retained fixtures `write_ini_<tool>.ini`, `help_<tool>.txt`,
  `helphelp_<tool>.txt` (without the `stty` line) and `swm_notest_output.mzML`
  come from the first execution.
* **Retained C++ outputs (tier 1).** `TOPPBase_test_write_ini_out.ini` and
  `TOPPBase_test_write_ini_subsec_out.ini` (cli c19e494), compared with the
  ported `FuzzyStringComparator` as `TEST_FILE_SIMILAR` does. The processing
  records of the retained `BaselineFilter_output.mzML`,
  `MapNormalizer_output.mzML`, `SpectraFilterWindowMower_1_output.mzML` and
  `MzMLSplitter_output_part1/2.mzML`, compared as decoded mzML in
  `tests/topp_baseline_filter.rs`, `tests/topp_map_normalizer.rs`,
  `tests/topp_spectra_filter_window_mower.rs` and `tests/topp_mzml_splitter.rs`
  beside their existing peak comparisons, which are unchanged.
* **Upstream class test (tier 3).** Transcribed with their literals:
  `getIniLocation_` (default), `getStringOption_` (default, command line, wrong
  type, unregistered, required), `getIntOption_`, `getDoubleOption_`,
  `getIntList_`, `getDoubleList_`, `getStringList_`, `getFlag_`,
  `inputFileReadable_`, `outputFileWritable_`, `parseRange_`, data processing
  methods, `getParam_`, misc options, subsection parameters, duplicate
  parameters and flag with trailing arguments, and the `-write_ini` part of
  `getStringOption_` (483-536): `TEST_EQUAL(p1, p2)` with `Param::operator==`
  semantics, where the expected version is the port's `Tool::VERSION` instead
  of `VersionInfo::getVersion()` (the class-test tool has no installed manifest
  and so reports the core version), and both `TEST_FILE_SIMILAR` comparisons.
  Not transcribed: the constructor, destructor and `main` sections (no object
  lifecycle; `main` is NOT_TESTABLE upstream), `setMaxNumberOfThreads`
  (NOT_TESTABLE; the native policy is tested), `getIniLocation_` with
  `-instance 5` and the `getStringOption_` INI cases (the source rejects
  `-instance` and top-level `common:` values in its strict update and only
  reads values left behind), `-log` and `Citation::toString` (not ported).
* **Native cases (tier 4)** cover `run_io` routing, the run-phase mapping of
  every `Error` variant, the initialisation failure, the argument bound (refused
  beyond it, parsed in full at it, with trailing-text chunks in command-line
  order), the context services, and a FIFO read once by a single writer.
* **Environment-dependent cases.** The socket case runs where a Unix socket
  can be bound in the temporary directory and does not open for reading; the
  terminal case runs only in a process without a controlling terminal, as under
  the gate and in CI. Both print `ran:` or `skipped:` with the reason; both ran
  on the Linux gate host.

## DTAExtractor and executed differential evidence

`src/bin/DTAExtractor.rs` is the first TOPP tool. Its three upstream tests plus
a real-instrument slice are reproduced in `tests/topp_dta_extractor.rs` against
the retained C++ outputs:

| Upstream test | Arguments | Result |
| --- | --- | --- |
| `TOPP_DTAExtractor_1` | `-rt :61` | `DTAExtractor_RT60.0.dta` byte-identical |
| `TOPP_DTAExtractor_2` | `-level 1` | `DTAExtractor_RT60.0.dta` byte-identical |
| `TOPP_DTAExtractor_3` | `-level 2 -mz :1000` | `DTAExtractor_RT140.0_MZ5.0.dta` byte-identical |
| (native) real Velos slice | none | all three produced names and all three files byte-identical |

These fixtures were produced by the C++ tool, so agreement is **tier 1 executed
differential evidence** under `docs/DIFFERENTIAL_VALIDATION.md` for the whole
chain: command line, parameter validation, mzML reading, the source number
formatter that names the output files, and DTA writing. This is the first
validated TOPP workflow in the port.

Every `.dta` fixture under `tests/data/` was rewritten on 2026-09-15 from the
pinned Release build at
`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`, run on
`ibminode06`. The three upstream `TOPP_DTAExtractor_*_output` files could not be
used: they predate the C++ move from Boost.Karma to `std::to_chars`
(`NumericFormatting.h`) and spell a peak m/z `120` where the pinned tool writes
`120.0`. The upstream TOPP suite compares with FuzzyDiff, so the stale spelling
still passes there; a byte comparison needs bytes the pinned tool actually
wrote. The added `dta_extractor_velos_slice.mzML` is three spectra of the
benchmark's LTQ Orbitrap Velos centroided run, cut with the pinned `FileFilter`
(`-rt 2.0:4.0 -mz 350:450`); its numbers are irregular enough to tell the
source's two formatting rules apart, which the round values of the upstream
cases cannot.

Three source behaviors had to be matched to get there, and all three were real
gaps:

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
3. **`DTAFile::store` writes two different numeric formats, neither of them
   Rust's.** The port wrote every number as Rust's shortest round-trip text,
   which is about 8 significant digits for an `f32` intensity. The source sets
   `os.precision(writtenDigits<double>(0.0))` — 15 — and then writes the
   precursor `MH+` mass and each intensity through the stream's *default float
   field* (`%.15g`, 15 **significant** digits, an `f32` promoted to `double`
   first) but each peak m/z through `DPosition<1>`'s `operator<<`, which calls
   `precisionWrapper` and so `StringUtils::toStr(double, true)` — 15 **fraction**
   digits via `NumericFormatting::appendNumeric`. One peak line therefore carries
   both `104.115715026855469` and `260.789154052734`, and an intensity of zero
   prints `0` on a line whose m/z prints `350.133800546540385`. Both rules are
   reused from `src/format/file_info/text_format.rs`
   (`ostream_g(·, WRITTEN_DIGITS_F64)` and `to_str`), which already ports them
   with oracle evidence; the tool's file names go through the same `to_str`.

   This was found by the TOPP benchmark's adversarial review: on the 1.2 GB
   Velos centroided run the port wrote 658,901,890 bytes where C++ wrote
   836,505,793, and the benchmark's "12% faster at one thread" was comparing
   21% less text. Re-measured on `ibminode06` after the fix, at `-level 2
   -threads 1`, both write 36,443 files and **836,505,793 bytes**, every file
   name and every byte identical (aggregate SHA-256
   `d1da4273041e597b57c466548f331702a8540bdfaa4d80aa57b0b1ba4c281c87`). The port
   costs it in wall clock: median of five interleaved repetitions 36.37 s
   against C++ 29.70 s (ratio 1.22), where the port before the fix took 26.33 s
   (0.89) for the smaller output. The fix adds 10.0 s over 48.0 M formatted
   numbers, about 200 ns each, because `ostream_g` formats every value twice —
   once as `{:.14e}` to find the decimal exponent, once in the chosen field —
   and both formatters return an owned `String` per number. That is a cost of
   the shared formatter, not of the DTA writer, and is recorded for its owner.
   Peak RSS is unchanged by the fix: 1.86 GiB against the C++ tool's 1.46 GiB.
   The same run at `-threads 32` writes bytes identical to the `-threads 1`
   output and to the C++ output, so the fix is thread-invariant.

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
with `no_spec`, which, like a `-parts 1` run without `-size`, is the source's
`writeLogError_` line followed by exit 6 as a code `main_` returns.

## MapNormalizer and SpectraFilterWindowMower

`MapNormalizer` scales MS1 peak intensities to a percentage of the run maximum;
its upstream test reproduces the retained C++ output. Higher MS levels are
untouched and the source's commented-out chromatogram branch is not ported.

**The run maximum is the combined one, chromatograms included.** `main_` calls
`exp.updateRanges()` and then `exp.getMaxIntensity()`.
`MSExperiment::updateRanges()` extends `combined_ranges_` first with the
spectrum ranges over every MS level and then with the chromatogram ranges
(`MSExperiment.cpp:665-719` at core bc9cc12), and `getMaxIntensity()` returns
that combined maximum (`MSExperiment.h:1059`). `MSExperiment::ranges` in this
crate covers spectra alone — the deliberate, documented scope of that query — so
taking it as the run maximum was wrong whenever a chromatogram carried the most
intense value. A TIC chromatogram sums a whole spectrum, so it routinely does.
The tool now asks `MSExperiment::combined_ranges_with_limits`, which is that
fold in that order and already existed; nothing about the C++ semantics is
restated in the tool.

That was a real divergence on real data, not a corner case. On the benchmark's
1.2 GB LTQ Orbitrap Velos centroided run the TIC chromatogram peaks at
1,788,496,256 against a spectrum maximum of 135,038,560, so every one of the
7,302 MS1 spectra came out 1,788,496,256 / 135,038,560 = 13.2443x too high while
the 36,443 MS2 spectra, which the tool never rescales, agreed bitwise. Both
maxima are what the C++ `FileInfo` prints for that file under *Combined Ranges*
and *Spectrum Ranges*. `MapNormalizer::run_maximum_intensity` now folds the
chromatogram ranges in, and the full run agrees with the C++ Release build on
every decoded array: 43,746 holders, 88,478,237 intensity points, 88,434,492 m/z
points and the chromatogram's 43,745 time points, all bitwise equal, zero
differences, with the port's output byte-identical at 1 and 32 threads. That
comparison was re-run after the empty-range refusal below was restored, and the
port's output digest is unchanged.

Two shapes distinguish the two formulas, and
`tests/topp_map_normalizer.rs` covers both against executed C++ output
(`tests/data/map_normalizer_chromatogram_{above,below}_*.mzML`): a chromatogram
above the spectrum maximum, which sets the scale on its own, and one below it,
which leaves the scale at the spectrum maximum. In the second fixture that
maximum sits in the MS2 spectrum, so MS1 normalizes to 6 rather than to 100 —
"the most intense MS1 peak becomes 100" holds only for a run with no
chromatogram whose most intense peak is an MS1 one, which is what the upstream
fixture happens to be and why its test never saw this.

**The ceiling on that query is derived from the map.**
`SummaryLimits::max_work`, which `combined_ranges_with_limits` charges one unit
per spectrum, per peak, per chromatogram and per chromatogram point against,
defaults to 50,000,000 — a figure no real LC-MS run fits. The Velos run charges
43,745 + 88,434,492 + 1 + 43,745 = 88,521,983 and `combined_ranges()` would
refuse it. The map is fully materialized by the time the tool asks for its
ranges, and what bounds the input is the reader's own size-derived ceilings in
`src/format/mzml_scaling.rs`, so `MapNormalizer::range_limits` charges the map's
own size: still bounded work, and it can never refuse a map the reader already
admitted. `the_derived_ceiling_is_exactly_what_the_query_charges` pins it at
exactly the charge — one unit less is refused — so it is neither slack nor
capable of refusing. The fixed default remains a defect for every other caller
of `combined_ranges`, `chromatogram_ranges` and `calculate_tic_binned`, and is
reported rather than changed here; `src/kernel/experiment_summary.rs` is outside
this tool.

**One deliberate deviation.** The source divides by `getMaxIntensity() / 100`
with no guard. Executed at the pinned Release build: on a run whose intensities
are all zero it exits 0 and writes `NaN` into every MS1 peak (the MS2 spectrum,
which it never rescales, keeps its zeros); on a run whose intensities are all
negative the maximum is the least negative value, so the scale is negative and
it exits 0 having flipped the sign of every MS1 peak — `[-10, -200, -3000]` came
back as `[1000, 20000, 300000]`. This port refuses a non-positive scale with
exit 6 instead, and writes nothing.

**An empty combined intensity range is refused, as the source refuses it** —
only the exit code differs. `main_` asks for `getMaxIntensity()` unconditionally,
one line before its peak loop, and `RangeBase::getMax()`
(`RangeManager.h:139-146` at core bc9cc12) throws `Exception::InvalidRange` on an
empty range with no assertion guard, so a Release build throws as well; whether
the loop body would have run is irrelevant, the throw happens first. Executed at
the pinned Release build on the two shapes these fixtures cover — three scans
carrying a retention time but no point, and an empty `spectrumList`; the
condition is simply that the run holds no spectrum peak and no chromatogram
point — the C++ exits **8**
with *Empty or uninitialized range object. Did you forget to call
updateRanges()?* and writes **no output file**; its `FileInfo` prints
`intensity: <none> .. <none>` under *Combined Ranges* for both. The port refuses
both and likewise writes nothing, but exits **6**: `Error::InvalidValue` and
`Error::InvalidRange` map to `ILLEGAL_PARAMETERS` crate-wide, where the source's
same-named exceptions derive from `BaseException` and reach its `UNKNOWN_ERROR`
arm — the mapping already documented at `run_failure` in `src/cli.rs`, not a
choice this tool makes. `an_empty_combined_intensity_range_is_refused` covers
both fixtures end to end; the executed commands, exit codes and digests are in
`tests/data/topp_map_normalizer_provenance.json`.

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
only chromatograms, and spectra that are not sorted by m/z. Both are exit codes
`main_` returns after its message (`INCOMPATIBLE_INPUT_DATA`), so the closing
line follows them.

BaselineFilter, MapNormalizer and SpectraFilterWindowMower add the processing
record the source attaches to their outputs, and MzMLSplitter reproduces the
source's attaching it to parts that hold nothing yet (decision D4; see
*Processing records* under *Preserved source conventions*). BaselineFilter's
warning when peak type estimation finds the first spectrum centroided and the
lines MzMLSplitter writes about the file size, the totals and each part are
ported, on the console and in the log (Release oracles `bf_centroided_log`,
`ms_parts2_log`, `ms_size_kb_log`). Not ported: SpectraFilterWindowMower's
`clearMetaDataArrays` and its warning (the tool keeps the arrays aligned
instead), and the source-compatible load options (decision D10) in the five
earlier tools, whose loaders are still the strict library default: an input
with a dangling mzML header reference, which the source reads, fails to load
there (seen on the Release run `bf_empty_log`'s input).

The second DTA finding above is the first concrete instance of the port's
"checked boundaries" convention blocking C++ parity. The resolution pattern —
keep the guard as the library default, add an explicit source-behavior option,
and have the tool opt in — is the one to apply as further tools meet their own
guards.

# Lane NOTES (TOPPBase finish) -- delete before the final commit

Branch: worktree wf_673154fa-1e1-1 off main a660a0c.
Ledger: headers/packages/cli/include/OpenMS/APPLICATIONS/TOPPBase.h = partial (see docs/core-sdk-reviewed-apis.json scope).
Rust: src/cli.rs, src/cli/{context,defs,logging,processing,spec,usage,param_ctd,tool_handler,console,outcome}.rs, src/cli/tools/**
serde_json =1.0.150 float_roundtrip optional; features rna-json, proforma-json use it. paramxml = quick-xml + file-compression.

Gaps (ranked):
(a) JSON INI ParamJSONFile::load -> wire serde_json into paramxml feature; add THIRD_PARTY_CRATE_DECISIONS entry.
(b) InvalidValue exit 6->8; first per-tool ILLEGAL_PARAMETERS returns (MzMLSplitter, MapNormalizer, FileInfo, SpectraFilterWindowMower, lifecycle).
(c) FeatureFinderCentroided CITATIONS.
(d) tool bodies' error/warn lines -> cli::log_error/log_warning (file_info.rs:185, peak_picker_hi_res, feature_finder_centroided).
(e) UnregisteredParameter/WrongParameterType exit 12.
(f) Windows COLUMNS -- record only.

## Findings
- C++: core bc9cc12 src/openms/source/FORMAT/ParamJSONFile.cpp load() (lines ~38-160); header src/openms/include/OpenMS/FORMAT/ParamJSONFile.h;
  test src/tests/class_tests/openms/source/ParamJSONFile_test.cpp. NO ledger entry for ParamJSONFile.h (propose one in fragment).
- cli c19e494 source/APPLICATIONS/TOPPBase.cpp: -ini JSON branch 283-296 (getTypeByFileName==JSON -> param_inifile_=getDefaultParameters_(); load; false->ILLEGAL_PARAMETERS),
  -write_ini branch 2619-2629 same. Exceptions caught 430-500: ParseError -> 'Error: Unable to read file (<what>)' exit 4 (INPUT_FILE_CORRUPT).
- load semantics: namespace = first trace node name + ':1:' (ignores -instance!); keys child.key() with '__'->':'; null skipped;
  object w/o string 'class' -> recurse; unknown key -> ParseError 'Parameter <key> passed to '<tool>' is invalid. To prevent usage of wrong defaults, please update/fix the parameters!';
  bool (valid_strings exactly true/false) get<bool>; input file: is_executable -> object?path:string, else node["path"]; int get<int64_t> (no bool);
  double get<double> (no bool); string list (input file -> ["path"]); int list vector<int> (bool accepted, 32-bit cast); double list; EMPTY -> OPENMS_LOG_WARN; setValue(key,value).
  nlohmann what(): '[json.exception.type_error.302] type must be X, but is Y', 305 'cannot use operator[] with a string argument with Y'.
- nlohmann vendored 3.11.2 (core src/openms/extern/nlohmann_json/single_include/nlohmann/json.hpp), JSON_DIAGNOSTICS 0.
  get<int64_t>/get<double> -> get_arithmetic_value (4610-4640): bool REJECTED 'type must be number, but is boolean';
  vector<int> elements -> generic arithmetic from_json (~4919): bool ACCEPTED. Empty input msg: 'syntax error while parsing value - unexpected end of input; expected '[', '{', or a literal' (no 'empty input' text in 3.11.2).
  ParseError what() = message + ' in: ' + expression (Exception.cpp:160-161) -> '<msg> in: '.
  serde divergences to record: '-0' -> f64 -0.0 (nlohmann int 0); '1e400' refused (nlohmann inf); syntax-error text differs; BOM: strip it (nlohmann skips).
- Rust Param: src/param.rs (entry, exists, set_value, valid_strings, tags; iter() items with trace), ParamValue in src/param/value.rs (IntegerList Vec<i32>).
  file_types::type_by_file_name (src/format/file_types.rs:238). Error::Parse Display 'parse error on line N: msg' -> cli must print message itself.
  Module goes in src/format/param_json.rs under #[cfg(feature="paramxml")] (src/format/mod.rs:35).
- Upstream test ParamJSONFile_test.cpp: FileNotFound, 'not a json' ParseError, all-types doc (bool1-4, int 5, double 6.1, string, int_list 10,11,12, double_list 13.25,15.125, string_list, file_output, is_executable_v1/v2).
- Rust: src/cli.rs load_ini (line ~1329) used at ~992 (-ini) and ~1190 (-write_ini); only paramxml::load today.

- DONE (c3cc1b9, builds): src/format/param_json.rs, cli.rs load_any_ini/load_json_ini + 'Assuming INI is a cwl file' debug; Cargo paramxml += dep:serde_json.
- Oracle: remote dir /ceph/ibmi/abi/oliver/oracle-runs/json-ini (ibminode06); local ../oracle/json-ini; retain in tests/data/topp_json_ini/<case>/{argv,exit_code,stdout,stderr,out.ini}.
  Replay pattern: tests/support/release_runs.rs (run_with_registry in-process; map oracle paths).

- ORACLE run1 fetched: ../oracle/json-ini/run1/results (ibminode06, Release). run2 (adds 3 debug-log cases) -> remote results2.
  Exit codes: ParseError -> 3 (INPUT_FILE_CORRUPT=3! enum: 0 OK,1 NOT_FOUND,2 NOT_READABLE,3 CORRUPT,4 EMPTY,5 CANNOT_WRITE,6 ILLEGAL,7 MISSING,8 UNKNOWN,...,12 INTERNAL).
  missing .json -> 1 'Error: File not found (the file '<p>' could not be found)'  (XML path in Rust says 'does not exist' = latent bug, FileNotFound ctor Exception.cpp:166-167 says 'could not be found').
  directory .json -> 12 'Unable to initialize or run PeakPickerHiRes: basic_filebuf::underflow error reading the file: Is a directory' (std::exception outer catch).
  1e400 -> '[json.exception.out_of_range.406] number overflow parsing '1e400'' exit 3 (reproducible!). -0 -> 0.0 (Rust gives -0.0: divergence).
  1e30 int -> -9223372036854775808 confirmed; lenient AND strict update accept it for min 0 (int cast quirk). -debug 1 without -log prints nothing.
- FFC citations DONE (c4a693b), completion 49 + console 4 tests pass (default features).

- run2 retained in tests/data/topp_json_ini (62 cases). tests/topp_json_ini.rs replays; tests/param_json.rs unit (9 pass).
  INI load errors now via session.error (log file) for XML + JSON; XML missing wording -> 'could not be found' (needs run3 XML evidence).
- FOUND, NOT MINE (report): paramxml writer float text '4' vs Release '4.0'/'1.0e20'/'0.0:' (all -write_ini!); paramxml writer + Param::validation_error
  refuse i64 outside i32 (src/param.rs:251 to_i32) where C++ truncates -> OPEN cases w_pphr_int_huge, r_pphr_huge.
- debug2 order difference documented in docs/TOPP_CLI_SUPPORT.md ~585 (checks before body).
- Gap (e) plan: ToolContext accessors -> record programming error in context (Cell), lifecycle maps escaping one to 12 'Internal error: Request for unregistered parameter '<name>'' / 'wrong type'. what() = name (TOPPBase_defs.h:53-75). catch TOPPBase.cpp:482-493.
- scratch helpers: scratchpad/dbgall.sh (replays cases with target/debug bins).

## Next step
2. Write src/format/param_json.rs (load) behind paramxml (+dep:serde_json); wire into load_ini.
4. tests/topp_json_ini.rs replaying run2 (retain in tests/data/topp_json_ini); unit tests in tests/param_json.rs (ParamJSONFile_test.cpp).
5. Gaps (b),(d),(e); docs; THIRD_PARTY entry; ledger fragment; sweep.
- Gap (b) trial (arm flipped to 8 in run_failure, uncommitted at first): failing = lifecycle run_phase_errors (asserts native 6), file_info a_forced_type_the_file_contradicts_is_refused,
  map_normalizer an_empty_combined_intensity_range_is_refused (Release: 8 + RangeManager.h:142 text), sfwm a_non_finite_window_size_is_refused_before_the_processing_record. Splitter passes.
- DONE: (b) arm -> 8 (fe1ce53, 924aea0: MapNormalizer InvalidRange+source text, SFWM/MapNormalizer native refusals via ToolError::invalid_parameter, FileInfo forced type -> ParseError 3 (report.rs forced_type_refused + peaks.rs precheck; FileHandler check_allowed still InvalidValue = FileHandler lane!)).
  (c) FFC citations. (e) 0b6a73a context.rs refusals -> run_body escaped_refusal. (d) was already done on main by 436878b; FuzzyDiff truncation warning now log_warning.
- run3 oracle fetched: ../oracle/json-ini/run3/results (x_missing_log, x_missing_write_ini, x_malformed_log, j_missing_log, fi_forced_*). Not yet retained in tests.
- NEXT: retain run3 + test; docs (TOPP_CLI_SUPPORT JSON section, exit table rows 133-160, FILE_INFO item 3, THIRD_PARTY serde_json); provenance json; ledger fragment; gates.
- DONE (6e88010): batch3 retained+tests, docs, provenance, ledger fragment (.planning/ledger-fragments/toppbase-finish.json, untracked/ignored), store side, fmt.
  Local checks (core_sdk, doc_coverage --report, module_cycles, source_citations) rc=0.
- NEXT: sweep on dax: OPENMS_GATE_HOST=dax OPENMS_GATE_JOBS=64 ~/.local/bin/openms-ci-sweep.sh toppbase-finish-1 <logdir> (background); slices for new tests; final commit deleting NOTES.md.
- Slices measured: paramxml/param_json 10; mzml paramxml/topp_json_ini 6 (+lifecycle 76); mzml paramxml featurexml/topp_json_ini 8 (completion 49, console 4, topp_file_info 43).
- Sweep toppbase-finish-1 on dax started at 0e504dd (log scratchpad/sweep1.out). Commits after: dedupe refusals + cli module doc (re-run lifecycle locally).
- file_info_a8 failed in sweep1 (old InvalidValue assertion) -> fixed (commit above); re-run failing lines after sweep ends.
- file_info_a8 failed in sweep1 (old InvalidValue assertion) -> fixed; re-run failing lines after sweep ends.
- sweep1 (at 0e504dd): 45/54 pass; 7 fails = file_info_a8 (fixed 57409b0), clippy (3 lints + complex type, fixed), doc (links, fixed). sweep2 next.

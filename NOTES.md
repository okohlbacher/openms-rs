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

## Next step
2. Write src/format/param_json.rs (load) behind paramxml (+dep:serde_json); wire into load_ini.
3. Oracle: Release runs with JSON INIs on ibminode06 (nlohmann messages, instance quirk, bool-in-int-list).

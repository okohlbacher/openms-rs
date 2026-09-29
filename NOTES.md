# B8c lane notes (BLOCK-05: TOPPBase remainder)

Merged prerequisite worktree-wf_ee239125-5be-1 (fast-forward to 9fcb962).
No earlier incarnation notes.

## Scope (from brief, ranked) and findings
1. -write_ini numbers. C++ ParamXMLFile.cpp (core bc9cc12) writeXMLToStream:
   - ITEM int/double value: `it->value.toString()` = ParamValue::toString(true):
     INT -> std::to_string(int64 ssize_); DOUBLE -> StringUtils::toStr(d,true) =
     NumericFormatting::appendNumeric(d, 15, false). Rust port: `ParamValue::to_text(true)`
     (src/param/value.rs format_float) already implements it.
   - restrictions: StringUtils::toStr(min_float) (same 15-digit rule), ints StringUtils::toStr(int).
   - double LISTITEM: `os << list[i]` with os.precision(writtenDigits<double>=15) -> `%.15g`
     (evidence: tests/data/topp_cli_lifecycle/TOPPBase_test_write_ini_out.ini has `4`, `2.33`).
   - int LISTITEM: IntList (32-bit) as is.
   Plan: WriteOptions gets `numbers` (Native default = shortest round-trip & refuse >i32;
   Source = the above). WriteOptions::source() selects Source (used by -write_ini, cli.rs:1204).
   Native test tests/paramxml.rs precision_nan... asserts lossless round trip -> keep native default.
   Integer width: C++ ParamValue::operator int() is a modular truncation (ParamValue.cpp:445-451);
   ParamEntry::isValid `int tmp = value` (Param.cpp:122). Rust: src/param.rs:251 uses checked
   to_i32 -> change to `as i32` (minimal B6-file change, report it). cli.rs source_validity_message
   `narrow` saturates -> truncate. getIntOption_ in src/cli/context.rs:455 check too.
   C++ reader asInt_ = toInt32 -> refuses its own -9223372036854775808 (ParseError) -> the test's
   assert_ini must compare text when both reads fail.
   Oracle: tests/data/topp_json_ini/{w_pphr_int_huge,r_pphr_huge} retained; OPEN list in
   tests/topp_json_ini.rs:46; assert_ini tolerance for double lines at :206-230 -> make exact.
2. JSON syntax-error text (nlohmann parse_error.101) -> src/format/param_json.rs; SYNTAX_ERRORS
   in tests/topp_json_ini.rs; needs a nlohmann-style lexer error reproduction + more oracle cases.
3. run_failure (src/cli.rs ~2322) prints Error::Parse as `parse error on line N: msg`; source
   ParseError::what() = `<msg> in: <file>`. Error::Parse has no file; Display in src/error.rs:32.
4. run_failure maps MissingInformation -> 7; source BaseException arm -> 8 + `Error: Unexpected
   internal error (<what>)` (TOPPBase.cpp:495-499 cli c19e494). No producer stands for
   RequiredParameterNotGiven (lifecycle handles that). Test tests/topp_cli_lifecycle.rs:1746
   asserts the native 7 -> change to source 8 (source-derived, comment at :1711 already says 8).
5. ParamCWLFile.h (unmapped): Release build has no TDL -> store throws runtime_error
   TDL_UNAVAILABLE (ParamCWLFile.cpp:331). Port as ParamJsonFile store refusal. TDL writer =
   named backlog (extern/tool_description_lib convertToCWL.h + yaml).

## Progress (commits on this branch)
- DONE item 1: paramxml NumberText::{RoundTrip,Source}; WriteOptions::source() -> Source.
  Param validation truncates (src/param.rs:250). cli source_int() (context.rs) used by
  ToolContext::int, new(), cmd_int, validate int_range, source_validity_message.
  Tests: tests/topp_write_ini_text.rs (new, byte-for-byte vs 9 retained INIs), topp_json_ini
  OPEN removed + exact INI, tests/param.rs 2 tests renamed to source semantics,
  tests/trace_fitter.rs to_param_and_from_param_agree_beyond_i32 (oracle gauss-trace-fitter).
  gate1+gate2 (full default nextest on spock): all green after fixes (5536 tests).
- DONE item 4: MissingInformation -> exit 8 "Unexpected internal error"; lifecycle test updated.
- DONE item 3: parse_what(): Parse{line:0} prints the message alone; 3 tests now exact.
- DONE item 5: src/format/param_cwl.rs ParamCwlFile refusal; cli write_commands stores through
  ParamCwlFile/ParamJsonFile; tests/param_cwl.rs.
- DONE item 2 (commit d4401e9a): src/format/param_json/nlohmann.rs (nlohmann 3.11.2 lexer/parser
  error path), tests/param_json_syntax.rs, tests/param_int_width.rs; Release oracle
  ../oracle/B8c-numbers (ibminode06, Release install, run twice, differ.txt empty) retained
  in tests/data/b8c_release + tests/data/b8c_release_provenance.json. Docs updated there too.
- Restart (incarnation 3, 2026-09-29): merged -36 (own earlier run) + -1 into
  worktree-wf_cd066dab-0a7-40. d4401e9a NOT yet gated.
## Next step
1. Gate on spock: nextest default + features touched (paramxml, json?), clippy -D warnings,
   doc, MSRV +1.85.0, then CI sweep (openms-ci-sweep.sh B8c-sweep <logdir>).
2. Count tests per slice for the new test files; write report + ledger; delete NOTES.md.

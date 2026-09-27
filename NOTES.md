# Lane B0 notes — integrate and verify the five Phase-B branches

Branch: worktree-wf_ee239125-5be-1 (off main f0990c0), merged worktree-wf_bcdd783d-f3f-1 (earlier incarnation).

## Branches merged (done by earlier incarnation, commits 9e52998..89cde73)
1. -1e1-1 TOPPBase tip 2ce5bf0   2. -1e1-2 MSExperiment tip 6203012   3. -1e1-3 LogStream tip b4ab801
4. -1e1-4 FileHandler tip d37eab0   5. -1e1-5 ProteinIdentification tip 6e544f9
No conflict markers remain.

## Tooling facts (took calls to find)
- git -C only works with a LITERAL absolute path, one git per command, no shell vars:
  git -C /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages/cli show c19e494:source/APPLICATIONS/TOPPBase.cpp
  topp pin paths: src/<Tool>.cpp ; core pin checkout: .reference/openms4-core-bc9cc12 (main repo)
- Oracle dir: /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/oracle/ ; template run.sh/remote_run.sh in sne-completion, TOPP-run template topp-exception-exits/cases.sh
- Sweep 1 log: $SCRATCH/sweep1/sweep.log (scratchpad)
- Lane fragments: .claude/worktrees/wf_673154fa-1e1-{1..5}/.planning/ledger-fragments/*.json (summarised with scratchpad/frag.py)

## Findings so far
F1 run_failure now maps Error::InvalidValue/InvalidRange -> 8 crate-wide (was 6). Source: TOPPBase.cpp:475-499.
   Reachable native parameter refusal that silently moved 6->8: BaselineFilter DataPoints with negative
   struc_elem_length (src/cli/tools/baseline_filter.rs:195-199); source has no restriction (MorphologicalFilter.h:139,177).
   Stale docs: fuzzy_diff.rs int_option doc ("reported as an invalid parameter"), peak_picker_hi_res.rs:312-324
   ("in-memory mode ... exit ILLEGAL_PARAMETERS").  Fix: invalid_parameter for BaselineFilter; docs.
F2 FileHandler check_allowed returns InvalidValue where FileHandler.cpp:859-863 throws ParseError
   ("type: X is not allowed for loading an experiment. Allowed types are: , a in: <file>"), now exit 8 vs source 3.
   FileInfo lane copied the check (src/format/file_info/report.rs forced_type_refused, peaks.rs load_experiment).
   Oracle row: tests/data/file_handler_load_release.tsv:222 (dta2d_bogus_allowed). Features: InvalidFileType (1261),
   consensus/ids ParseError (1409, 1522). Fix: make FileHandler raise the source error, FileInfo call it; one copy.
   tests/file_handler_load.rs:137-148 and 789-808 assert the divergence (update to source value from oracle row).
F3 serde_json float_roundtrip is on (Cargo.toml:91); no float round-trip test exists for JSON INI -> need oracle.

F4 Oracle ../oracle/B0-attacks DONE (ibminode06, Release, gcc 14.4.0, 2 runs): 51 JSON tokens bits; unknown_* 7 tools exit 3;
   bf_dp_negative SIGSEGV 139; bf_dp_negative_fraction exit 0. Retained tests/data/b0_release (+ provenance json).
F5 Verified OK (no change): compute_hash/rewrite_source_file vs file_handler_load_release.tsv (sf rows all fields);
   MSExperiment operator<< whole-dump vs msexperiment-close stdout; PI engine-order claim (registry order, MetaInfo.cpp:162-170);
   LogStream exit flush (ProcessLogs in run_from_environment, no process::exit, panic=unwind); LogStream.h members in doc.
   MSExperiment complete: own members all in docs/EXPERIMENT_MOBILITY_SUPPORT.md; inherited MetaInfoInterface is partial (UInt/registry) -> note for lead.

## Done (commit 5486102 WIP)
- FileHandler check_allowed(kind, allowed, file, LoadTarget) -> source ParseError/InvalidFileType texts; FileInfo peaks/report/identifications/features use it
- BaselineFilter: MorphologicalFilter::new + invalid_parameter (exit 6), after input checks
- docs: TOPP_CLI_SUPPORT, FILE_HANDLER_SUPPORT (div 4 closed), FILE_INFO_A7 (diff 4 closed), A8 (diff 1 closed), TOPP_FILE_INFO (featureXML exit 8), PPHR doc 9, pphr/fuzzy rustdoc
- tests: file_handler_load (2), file_info_a7, topp_cli_lifecycle, topp_cli_completion, topp_file_info (featureXML exit 8), topp_json_ini (+fi_forced_featurexml_on_dta), param_json (+51 bits), NEW tests/topp_b0_release.rs

## HAZARD: the scratchpad dir is SHARED with other agents (gate1.log, frag.py clobbered). Use scratchpad/B0lane/ only.

## Status
- [x] merges 1..5
- [x] sweep 1 of merged tree (before fixes): 54 lines, 0 failed (scratchpad/sweep1/sweep.log)
- [x] targeted gate of fixes: all green (scratchpad/B0lane/gate_targeted2.log): file_handler_load 11, file_info 61, file_info_a7 34,
      file_info_a8 25, param_json 11, topp_b0_release 2, topp_baseline_filter 3, topp_cli_completion 49, topp_cli_lifecycle 76,
      topp_file_info 46, topp_fuzzy_diff 17, topp_json_ini 9, topp_peak_picker_hi_res 36
- [x] lead files (CHANGELOG, OpenMS_CPP_ISSUES CPP-355, PORTING_STATUS, VALIDATION) restored to main (aaeadf1); lanes' text in
      main-repo .planning/lane-reports/B0-lead-files.patch
- [x] ledger fragment written: main-repo .planning/ledger-fragments/wave4/B0.json (scratchpad/make_fragment.py -> copy in B0lane)
- [ ] sweep 2 running (slot B0-sweep, log scratchpad/B0lane/sweep2/)
- [ ] slice measurements running (scratchpad/B0lane/slices.sh -> slices/summary.txt)
- [ ] python gates: ci_matrix.py --check and core_sdk_coverage.py fail = lead's regeneration (new test files, provenance manifests)
- [ ] report -> main repo .planning/lane-reports/B0.md

## Next step
1. When sweep2 + slices done: write report B0.md (per header, slices+counts, Cargo.toml, provenance, CHANGELOG line, defect candidates, attacks).
2. Final commit deleting NOTES.md.

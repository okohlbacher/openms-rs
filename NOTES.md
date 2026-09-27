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

## Status
- [x] merges 1..5
- [ ] sweep 1 of merged tree running (dax, slot B0-sweep)
- [ ] attacks: exit codes (F1,F2), JSON float round-trip (F3), compute_hash/rewrite_source_file bytes, exit flush paths, MSExperiment operator<<, PI engine ordering
- [ ] ledger fragment consolidation -> .planning/ledger-fragments/wave4/B0.json
- [ ] report -> main repo .planning/lane-reports/B0.md

## Next step
1. Survey lane 2-5 evidence (tests vs oracle results) for the remaining attacks.
2. Build ../oracle/B0-attacks (driver: JSON doubles %a; TOPP cases: BaselineFilter negative DataPoints) on ibminode06.
3. Fix F1, F2.

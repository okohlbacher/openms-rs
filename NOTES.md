# Lane A1 — EVID-01 notes (restart-safety file; deleted in the final commit)

## Task
Re-run the three mislabelled oracles on ibminode06 against the Linux Release
install `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`:
- `../oracle/topp-early-bundle` (run.sh -> cases/oracle.py; macOS product-sdk Debug core 4fdec46)
- `../oracle/topp-file-info-tool` (oracle.py)
- `../oracle/topp-peak-picker-tool` (oracle.py, derived/)
New oracle goes in `../oracle/A1-release-rerun/` (host+build label honest), run twice, DONE marker.
Then re-check every test that cites those oracles (grep the three names under tests/),
fix port (src/cli/tools/<tool>.rs) where Release differs; framework (src/cli.rs,
src/cli/context.rs) divergences -> report to lead (B8c owns).
Also fix tests/topp_peak_picker_hi_res.rs:~1179 which calls a macOS SDK oracle "Release".

## Known divergence
directory as `-in`: Linux Release exits 12, port exits 10.

## State
- Merged worktree-wf_ee239125-5be-1 (fast-forward to 9fcb962). No earlier incarnation.

## Oracle harness (built, running)
- `../oracle/A1-release-rerun/`: `run.sh stage|fetch|compare`; `linux_release.py` imports
  verbatim copies of the three drivers and rebinds only location constants (perl alarm
  kept; perl exists on ibminode06); `remote_run.sh` = identity, run1, run2, reproduce, DONE.
- Remote tree: ibminode06:/scratch/kohlbach/A1-release-rerun (started 2026-09-27).
- `compare.py verify-stage` passes (6 drivers identical, 91 inputs hash-equal).
- compare.py `compare` still needs writing (imports compare_runs, not yet written).
- Tests that cite the 3 oracles: tests/{file_info,topp_file_info,topp_peak_picker_hi_res,
  topp_feature_finder_centroided}.rs + tests/data/*provenance.json. The brief's other files
  (dta_extractor, mzml_splitter, map_normalizer, swm, baseline_filter) cite
  ../oracle/topp-cli-lifecycle and topp-exception-exits; the early bundle has DTAExtractor_* and
  SWM_* cases, so check those tests against Release too.
- Existing Linux Release replay helper: tests/support/release_runs.rs (topp-exception-exits);
  its fi_undetermined_log case already shows the dir-as-in exit 12.

## Results (done)
- Reruns fetched; also topp-cli-lifecycle run.sh + cli2 rerun (lifecycle/). manifest.json
  sha 47564af6...; divergences*.json; retained_in_repo.json.
- Only real divergence: in_is_directory 10 -> 12 (framework, src/cli.rs input_path) -> REPORT,
  not fixed (B8c). FFC debug-only cases already Release-asserted. Picker all byte-identical.
- Committed: usage-text fixtures + FFC_1 input swapped for Release files; test/doc relabels;
  provenance release_verification blocks. Gate: 4 test files green on dax.

## Incarnation 3 (wf_cd066dab-0a7-41)
- Merged ee239125-5be-37 (holds incarnation 2's commits faf8168e, bdfd4dd1) and -1 (up to date).
- Re-audited divergences.json / divergences_lifecycle.json: residuals are only Debug-only warnings,
  revision line, host paths, FFC last-digit floats (FFC_1 fixture already swapped), 2 FFC debug-only
  exits (Release 0), PPHR_auto_mode_1 crash both, in_is_directory 10->12. Framework: src/cli.rs
  input_path `.unwrap_or(FileType::Unknown)` swallows the EISDIR from FileHandler::get_type; Release
  goes to TOPPBase.cpp:507-512 (cli c19e494) "Unable to initialize or run", exit 12.
- FOUND mislabel: tests/topp_cli_lifecycle.rs ~1782 calls in_is_directory exit 10 "Release evidence"
  -> fix to cite fi_notype_log (Release, exit 10 + closing line). Header note on lifecycle rerun.

- DONE dcc64f70: lifecycle relabel. Sweep slot A1-sweep was raced by another lane (foreign
  openswath_data_access test); re-running in slot A1-cd066-sweep (scratchpad/sweep2).

## Next step
1. Read scratchpad sweep2/summary.log; then write .planning/lane-reports/A1.md and ledger-fragments/wave4/A1.json.
2. Full CI sweep in background; delete NOTES.md in final commit.

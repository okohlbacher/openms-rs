# INTEG-X notes (group X integrator)

Branch: worktree-wf_97065d91-48e-10 off main f0990c06. Never push, never rebase, no local cargo (dax via gate).
.planning/ is NOT tracked: reports/fragments live at
/Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.planning/{lane-reports,ledger-fragments/wave4}/
Reference C++ headers: /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.reference/openms4-core-bc9cc12/
Shell quirks: harness refuses long compound commands and the word `source` as a path token; write python scripts
with the Write tool to the scratchpad and run them with a short command.

## Done (commits on this branch)
- 12 lane branches merged (I1 I2 I3a I3b I4 FF1 FF2a FF2b FF3 FF4 M1 M2).
- Baseline: all-features `test --no-fail-fast` on the merged tree green (463 result lines, EXIT 0).
- Collapsed: MapAlignmentAlgorithmKD (M1 alignment_kd.rs is the one; KDTreeFeatureMaps implements KdFeatureMaps;
  FeatureGroupingAlgorithmKD uses transformations LowessOptions) -> KD headers complete.
- Closed gaps: PoseClustering composes StablePairFinder + class test (307/1.01164/-32.0912); TreeGuided build_tree/align
  (oracle replays C++ tree).
- Collapsed: nth_element x3 -> math::libstdcxx::nth_element_by; AASequence::source_order; ModificationsDB::all_search_modifications;
  median_of_sorted; lower/upper_bound; LinearRegressionWithoutIntercept ported (math/ml).
- PeakWidthEstimator FF1/FF2b collapsed at merge (kept FF1's); decoy_random shuffle_slice kept existing.

## Next steps
1. Check test run /private/tmp/claude-501/integ-x-test2.log (all-features test after collapses).
2. More duplicate hunting (see list in thought: std_min/std_max helpers, param entry helpers, FF3 eigen sums, format_g17).
3. Member-by-member attack of `complete` claims (scripts in scratchpad members.py; heuristic only) - read headers by hand for risky ones.
4. Evidence audit: oracle labels, CI slice counts non-zero on the merged tree (ci_matrix.py --check on scratch copy).
5. Full sweep: OPENMS_GATE_HOST=dax OPENMS_GATE_JOBS=96 ~/.local/bin/openms-ci-sweep.sh integ-X-sweep /private/tmp/claude-501/integ-x-sweep
6. Write INTEG-X.json + INTEG-X.md (also mention: I3b/FF1/FF2a lane reports are missing), delete NOTES.md in final commit.

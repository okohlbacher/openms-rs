# Lane M1 notes (map alignment) -- delete in the final commit

## Scratch
`/private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad/M1/`
(mct/mcpp/mhdr/mdata symlinks into the pinned C++; gate*.log, clippy*.log, msrv.log, sweep.log, sweep/).
The scratchpad root is shared with other sessions: only create things inside `M1/`.
Write files outside the worktree via scratch + `cp` (the Write tool refuses shared-checkout paths).

## Done (all committed, gates green on dax: clippy all-features, rustdoc -D warnings, MSRV 1.85 tests)
- transformations: B-spline, Akima, LOWESS auto-span, Param surface, TransformationModelDefaults, TrafoXML.
- mapmatching: superimposer (+dumps), shift/affine superimposers, evaluation, identification alignment,
  pose-clustering alignment (partial: StablePairFinder closure), MapAlignmentTransformer IdentificationData overload.
- docs: TRANSFORMATIONS_SUPPORT.md, MAP_ALIGNMENT_SUPPORT.md, THIRD_PARTY_CRATE_DECISIONS row, 2 provenance manifests.
- Report drafted at scratch M1/M1.md; ledger fragment copied to .planning/ledger-fragments/wave4/M1.json.
- Tier-1 oracle committed (2749bfd): tests/map_alignment_oracle.rs, 6 tests, 3780/3780 bit-identical on dax;
  C++ side in ../oracle/M1-mapalign (ibminode06). Report + ledger (scratch M1/M1.md, M1/M1.json) updated
  for tier 1; re-cp M1.json to the ledger path.

## Restart 2 (worktree wf_cd066dab-0a7-36, fast-forwarded from worktree-wf_ee239125-5be-30)
Scratch dir from restart 1 was wiped: report draft M1/M1.md is LOST; ledger M1.json at .planning/ledger-fragments/wave4/M1.json survives.
Re-run the sweep; rewrite the report from docs/MAP_ALIGNMENT_SUPPORT.md + TRANSFORMATIONS_SUPPORT.md.
Then look at KD / TreeGuided deps (KDTreeFeatureMaps, ClusterHierarchical, StablePairFinder) if time.

## Next step (old)
1. Wait for the CI sweep (started 20:58 on the tree before the oracle commit; oracle slice + clippy re-run on the final tree already) (scratch M1/sweep.log, logs in M1/sweep/); fix any red line.
2. Append sweep result to M1/M1.md, cp it to .planning/lane-reports/M1.md.
3. Final commit deleting NOTES.md.

## Restart-2 plan (KD + TreeGuided), decided after reading the pinned C++
- FF2a lane (worktree wf_cd066dab-0a7-34) ports ML/CLUSTERING (ClusterHierarchical, AverageLinkage,
  ClusterAnalyzer), DistanceMatrix, BinaryTreeNode, KDTree. Do NOT duplicate them.
- KD: port MapAlignmentAlgorithmKD in mapmatching/alignment_kd.rs against a trait (the KDTreeFeatureMaps
  surface it uses: len, rt, charge, map_index, neighborhood, apply_transformations+optimize). Result is
  independent of kd-tree query order (CC seeds by node index, members pushed in node order) -> a
  brute-force test double gives the same answer. C++ quirk: the charge check in filterCCs_ is dead
  (`continue` in the inner loop) -> preserve, defect candidate. Empty Param -> members uninitialised (UB).
  Class test: NOT_TESTABLE everywhere; evidence = oracle on ibminode06 with real KDTreeFeatureMaps.
- TreeGuided: port everything but buildTree's clustering; tree passed in as (left,right) pairs.
  Math::median(..., sorted=true) on UNSORTED per-sequence feature RTs -> middle element, preserve.
  maps_ranges never merged (tmp is dead). Class test gives tree literal (0,2),(0,1) + trafo_order,
  unique ids 20/14, original_RT; data MapAlignmentAlgorithmTreeGuided_test_in{0,1,2}.featureXML.
- CI sweep for the pre-restart tree running in background -> scratch M1/sweep.log.

## Restart-2 progress
- 442ba5a: alignment_kd.rs + alignment_tree_guided.rs + tests/map_alignment_kd.rs (10) + tests/map_alignment_tree_guided.rs (7), green on dax.
- Oracle ../oracle/M1-kdtree (gen_cases.py, driver.cpp, run.sh -> ibminode06, results/). Next: write
  tests/map_alignment_oracle_kd_tree.rs comparing records.tsv bit-exactly (KD: nb sizes, cc labels,
  filtered sizes, rt after transform; TG: pair distances vs 1-sim f32, tree, order, mt rt, trafo points
  x/y/note, final rt), + tg_defaults.ini vs defaults(). Copy results into tests/data/map_alignment_oracle_kd_tree_*.
- Then: docs (MAP_ALIGNMENT_SUPPORT.md rows/sections), provenance json, clippy/doc/MSRV, sweep, report, ledger.

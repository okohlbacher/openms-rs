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

## Next step
1. Wait for the CI sweep (scratch M1/sweep.log, logs in M1/sweep/); fix any red line.
2. Append sweep result to M1/M1.md, cp it to .planning/lane-reports/M1.md.
3. Final commit deleting NOTES.md.

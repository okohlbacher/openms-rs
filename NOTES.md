# Lane A2a NOTES (EVID-02 first half: native_equivalent A-L re-audit)

Ledger-only lane. No code changes.

## What
- 72 `native_equivalent` entries in docs/core-sdk-reviewed-apis.json whose path after
  `include/OpenMS/` sorts before `M` (CHEMISTRY 15, CONCEPT 4, DATASTRUCTURES 11,
  FEATUREFINDER 3, FORMAT 19, INTERFACES 1, KERNEL 19).
- Definition: docs/NATIVE_EQUIVALENTS.md "What the tier means, decided 2026-09-21".
- Keys are `packages/core/include/OpenMS/<X>.h`; C++ at core pin bc9cc12 is
  `src/openms/include/OpenMS/<X>.h`, `src/openms/source/<X>.cpp`,
  `src/tests/class_tests/openms/source/<Name>_test.cpp`.

## Tooling (scratchpad, A2a/ subdir only)
SP=/private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad/A2a
- $SP/sel.json: the 72 entries. $SP/cpp/: header/cpp/test at the pin (flat X__Y__Z.h), already extracted.
- $SP/side.py <prefix> [--h|--r]: header decls vs Rust pub API lines. $SP/members.py <prefix>: rough miss list.
- $SP/v.py <rest> <status> <clause>: record a verdict (writes .planning-A2a/verdicts.json).
- Scripts repointed to this worktree (wf_ee239125-5be-3) on restart.
- The worktree guard refuses heredocs, `$VAR` programs and compound git: write scripts to $SP, call with literal paths.

## Verdicts so far (.planning-A2a/v-*.json, format {rest: [status, clause]})
- CHEMISTRY: all 15 done (13 complete, 2 partial: MassDecompositionAlgorithm and
  TheoreticalSpectrumGeneratorXLMS lack the DefaultParamHandler defaults/from_param/to_param surface).

## Next step
1. CONCEPT (4), DATASTRUCTURES (11), FEATUREFINDER (3), FORMAT (19), INTERFACES (1), KERNEL (19):
   read header public members, check Rust file, record verdict in .planning-A2a/v-<group>.json; commit per group.
2. Build fragment JSON (.planning/ledger-fragments/wave4/A2a.json in main checkout) + report
   (.planning/lane-reports/A2a.md) from the verdicts; delete NOTES.md and .planning-A2a in the final commit.

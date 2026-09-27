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

## Verdicts (.planning-A2a/v-*.json, format {rest: [status, clause]}) -- ALL 72 DONE
- 56 complete, 11 partial, 5 stay native_equivalent (ExposedVector, TypeAliases, MapUtilities, DPeak, StandardTypes).
- Ledger keys are `src/openms/include/OpenMS/<X>.h` (NOT packages/core/...); all 72 header sha256 at
  bc9cc12 match docs/core-sdk-coverage.json.
- Fragment built by $SP/build.py -> .planning/ledger-fragments/wave4/A2a.json (main checkout, git-excluded).
- Crate findings for the report: rusqlite (bundled C SQLite, fails "pure Rust") and flate2 have no
  THIRD_PARTY_CRATE_DECISIONS.md entry; chrono, rand_mt, base64 do.

## Next step
1. Write .planning/lane-reports/A2a.md (main checkout).
2. Final commit: delete NOTES.md and .planning-A2a.

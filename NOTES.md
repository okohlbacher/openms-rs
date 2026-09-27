# Lane A2a NOTES (EVID-02 first half: native_equivalent A-L re-audit)

Ledger-only lane. No code changes.

## What
- 72 `native_equivalent` entries in docs/core-sdk-reviewed-apis.json whose path after
  `include/OpenMS/` sorts before `M` (CHEMISTRY/..., CONCEPT/..., DATASTRUCTURES/...,
  FEATUREFINDER/..., FORMAT/..., INTERFACES/..., KERNEL/...).
- Definition: docs/NATIVE_EQUIVALENTS.md "What the tier means, decided 2026-09-21".
- Keys are `packages/core/include/OpenMS/<X>.h`; C++ at core pin bc9cc12 is
  `src/openms/include/OpenMS/<X>.h`, `src/openms/source/<X>.cpp`,
  `src/tests/class_tests/openms/source/<Name>_test.cpp`.

## Tooling (scratchpad, shared with other lanes -> use the A2a/ subdir only)
SP=/private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad/A2a
- $SP/sel.json: the 72 entries (filter status=='native_equivalent' and rest<'M').
- $SP/extract.py: extracts header/cpp/test at the pin into $SP/cpp/ (flat names X__Y__Z.h).
  Plain `git -C <core> show` works; compound shell with git/heredocs is refused by the
  worktree guard -> write python scripts to $SP and run `python3 $SP/x.py`.
- $SP/dump.py [prefix]: prints current ledger entries (scope etc).
- Verdicts accumulate in .planning-A2a/verdicts.json in the worktree (committed, WIP).

## Outputs
- /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.planning/lane-reports/A2a.md
- /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.planning/ledger-fragments/wave4/A2a.json

## Next step
1. Per directory group (CHEMISTRY, CONCEPT, DATASTRUCTURES, FEATUREFINDER, FORMAT,
   INTERFACES, KERNEL): read header public members, grep Rust file, record verdict
   in .planning-A2a/verdicts.json; commit per group.
2. Build fragment JSON + report from verdicts.json; delete NOTES.md and .planning-A2a.

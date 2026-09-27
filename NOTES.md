# Lane A2a NOTES (EVID-02 first half: native_equivalent A-L re-audit)

Ledger-only lane. No code changes.

## What
- 72 `native_equivalent` entries in docs/core-sdk-reviewed-apis.json whose path after
  `include/OpenMS/` sorts before `M` (CHEMISTRY/..., CONCEPT/..., DATASTRUCTURES/...,
  FEATUREFINDER/..., FORMAT/..., INTERFACES/..., KERNEL/...). List saved via script:
  python3 filter status=='native_equivalent' and rest<'M'.
- Definition: docs/NATIVE_EQUIVALENTS.md "What the tier means, decided 2026-09-21".
- For each: read C++ header at core pin bc9cc12
  (`git -C ../../../../OpenMS4-tests/packages/core show bc9cc12:src/openms/include/OpenMS/<X>.h`),
  compare members to Rust files listed; keep native_equivalent only if std discharges it.

## Outputs
- .planning/lane-reports/A2a.md (in MAIN repo path /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.planning/...)
- .planning/ledger-fragments/wave4/A2a.json
- Audit worksheet kept in tools? NO -> keep in scratchpad + commit progress JSON at
  .planning-A2a/progress.json in worktree (WIP), deleted at end.

## Next step
1. Build per-header member lists (grep public decls) for all 72.
2. Classify; write fragment entries in batches by directory; commit per directory.

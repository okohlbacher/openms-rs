# Lane A4 notes (DEBT-05, EVID-03, EVID-04)

Worktree wf_ee239125-5be-7 (branch worktree-wf_ee239125-5be-7), fast-forwarded
onto worktree-wf_bcdd783d-f3f-7 (1165047, this file only). Scratch scripts live in
the session scratchpad; regenerate them if lost.

## State found
- OpenMS_CPP_ISSUES.md: 6601 lines, CPP-001..354, no duplicates, no gaps, but not
  in file order (e.g. CPP-046 after "Maintaining this log"). Index table covers
  001..058 only. Header pins core `82ce5b3`; that commit EXISTS in
  ../OpenMS4-tests/packages/core, is an ancestor of bc9cc12 (8 commits back), and
  76 blob links in entries 1-58ish do use it. Later entries cite bc9cc12 (core),
  c19e494 (cli), 174b576 (topp) and "executed on product SDK at core 4fdec46".
  So the brief's "none uses it" is wrong; the header is wrong because it claims
  ALL entries were checked at 82ce5b3.
- Entry formats vary: early "Affected file/Issue and reproduction/Evidence/Rust
  handling"; later "Source revision/Status/Affected file/function/Trigger/Issue/
  Proposed C++ fix/Evidence/Rust handling"; FORMAT wave "Status and source".
- tools/check_source_citations.py --report --verbose: exit 0, 4136 checked,
  35 unreachable (names only, no locations), 3 bare ranges (not listed).
  Unreachable: Eigen (lmpar.h 7, Redux.h 3, InnerProduct.h 3,
  GeneralMatrixVector.h 3, Dot.h 2, TriangularSolverVector.h 2,
  ProductEvaluators.h 1), Boost (normal.hpp 8, erf_inv.hpp 1), libstdc++
  bits/random.h 1, TOPPBase.cpp 2 (issue-log entries naming core only?),
  driver.cpp 2 (../oracle/progress-format-readers/driver.cpp).
  Retained-source mechanism: SOURCE_PROVENANCE.json retained_sources.directories
  (lead's file) -> give lead the JSON block; retain headers under ../oracle/A4-*/.
- CI runs only test_source_citations.py (no pins in CI), and check_core_sdk.py
  without --source.

## Next step
1. Parse every CPP entry (scratch parse_cpp.py) -> classify -> index table.
2. Fix header of OpenMS_CPP_ISSUES.md; replace index; commit.
3. Checker: --verbose should list each unreachable citation and bare range with
   document:line (strengthening). Resolve: retain Eigen/Boost/libstdc++ headers
   from ibminode06 Release deps into ../oracle/A4-retained/ with sha256; check
   TOPPBase.cpp/driver.cpp cases; record unresolvable ones with reasons.
4. Measure sha256 recomputability for REPRODUCIBILITY.md (check_core_sdk.py).
5. Lane report + ledger fragment; delete NOTES.md.

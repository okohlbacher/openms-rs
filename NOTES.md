# Lane A5b NOTES (restart checkpoint)

Scope list: `/Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.planning/a5-lists.json` key `A5b` (37 headers).
Outputs: report -> main repo .planning/lane-reports/A5b.md ; ledger -> main repo .planning/ledger-fragments/wave4/A5b.json

## Tooling (worktree guard quirks)
- Scratch: `<scratchpad>/A5b/` (scratchpad is SHARED with other lanes; top-level files there are not ours except
  dump_cpp.sh, ledger_a5b.*). C++ sources for all 37 headers: `<scratchpad>/cpp/<Name>.h|.cpp|_test.cpp`.
- Print ledger scopes: `python3 <scratchpad>/A5b/scopes.py <filter...>` from the worktree root.
- Shell variables as program paths, `cd x && ...` chains and heredocs are refused by the guard: use literal paths,
  Write tool for scripts, one plain command per call.
- Ledger statuses allowed: complete | partial | native_equivalent (tools/core_sdk_coverage.py:234). No 'impossible'.

## Triage (second incarnation, from ledger scopes)
- JavaInfo/PythonInfo/RWrapper: every member mapped; 'partial' reasons are API-shape differences -> COMPLETE after
  member-by-member check.
- Gumbel: fitWeighted (h:81) and GumbelDistributionFitResult::eval (h:51) have NO definition at bc9cc12 -> complete
  with impossibility evidence (report: lead may want an 'impossible' status).
- MultipleTesting: only gap = Pi0Smoother seam. Close with a library adapter in processing implementing Pi0Smoother over
  processing::spline::smoothing::BSplineSmoothingSpline; switch tests/multiple_testing.rs to it.
- SysInfo/StopWatch/BuildInfo: platform bindings (macOS/Windows) -> leave partial, name crates.
- fancy-regex =0.19.2 IS a dependency -> SpectrumNativeIDParser/SpectrumLookup 'no regex engine' scopes are STALE.
- Cheap metadata gaps: Software (Ord + Hash), Precursor (4 activation-method string helpers + Hash),
  PeptideHit/ProteinHit (comparator functors, hash), ChromatogramSettings (hash, Display).
- Genuine ports: USI.h (1 TOPP consumer), MetaInfoRegistry.h.
- No Rust at all: Matrix.h, SIMDe.h, SimpleSVM.h, MetaInfoRegistry.h, USI.h.

## Done
(none yet)

## Next step
1. Confirm JavaInfo/PythonInfo/RWrapper member-by-member; write ledger fragment entries.
2. Gumbel entry; MultipleTesting adapter.
3. USI.h port; MetaInfoRegistry.h port.
4. Cheap metadata gaps (Software, Precursor, PeptideHit, ProteinHit).

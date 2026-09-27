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
- MultipleTesting: src/processing/spline/pi0_smoother.rs (SplinePi0Smoother, pi0_est_smoothed); tests/multiple_testing.rs
  20 pass on spock (default features). New acyclic edge processing->math (tool exit 0; docs/module-cycles.json NOT written).
- JavaInfo/PythonInfo/RWrapper confirmed member-by-member -> complete; tests/system_process.rs 34 pass on spock.
- Gumbel: verified at bc9cc12 (h:51 eval, h:81 fitWeighted undefined; test :175 uses GumbelMaxLikelihoodFitter).
- SpectrumNativeIDParser: src/metadata/spectrum_native_id_parser.rs over concept::boost_regex::BoostRegex (tokens =
  sregex_token_iterator); tests/spectrum_native_id_parser.rs; docs/SPECTRUM_NATIVE_ID_PARSER_SUPPORT.md.
- USI: src/metadata/usi.rs; tests/usi.rs (all 36 sections); docs/USI_SUPPORT.md. (gate2-4 logs in scratch A5b/)

## Next step
1. Check gate2-4 (usi, parser, lib unit, doctests); fix; commit.
2. Write ledger fragment entries so far (main repo .planning/ledger-fragments/wave4/A5b.json).
3. MetaInfoRegistry.h port.
4. Cheap metadata gaps (Software, Precursor, PeptideHit, ProteinHit); SpectrumLookup over the new parser.

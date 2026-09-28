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

- MetaInfoRegistry done (11 tests spock). Software/Precursor/PeptideHit/ProteinHit value ops (src/identification/hits.rs),
  ChromatogramSettings/SpectrumSettings Hash+Display (tests/settings_hash.rs) - all tests + clippy green on spock.
- Ledger builder: <scratchpad>/A5b/ledger.py (run from worktree root) -> writes main-repo fragment. ADD entries for:
  Software, Precursor, PeptideHit, ProteinHit, ChromatogramSettings, SpectrumSettings (all now complete by audit).

## Third incarnation (branch worktree-wf_cd066dab-0a7-4, merged worktree-wf_ee239125-5be-9 = fbf17f3)
- Old scratchpad (scopes.py, ledger.py, cpp dumps) is GONE; new scratch:
  /private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad/A5b
- Ledger fragment already exists in main repo (.planning/ledger-fragments/wave4/A5b.json) - edit it directly.
- Done since: IMDataArrayUtils + IMDataConverter frame reshaping (76b773c); MetaInfo index overloads (42783ef);
  SpectrumLookup + concept/string_conversion WIP (fbf17f3, UNGATED).

## Third incarnation progress
- Gate2 (spock) green: spectrum_lookup, usi, sni parser, meta_info_indexed, settings_hash, software,
  multiple_testing, matrix 15, imaging_types 29, on_disc_imzml 53, im_data_array_utils 7, im_data_converter 19.
- NEW: src/data_structures/matrix.rs + tests/matrix.rs + docs/MATRIX_SUPPORT.md (Matrix.h complete).
- NEW: tests/imaging_types.rs (3 imaging class tests) + Default for ImagingRegion -> IonImage/Geometry/Region complete.
- NEW: src/metadata/spectrum_meta_data_lookup.rs + tests + docs/SPECTRUM_META_DATA_LOOKUP_SUPPORT.md (gate3 pending).
- Scratch scripts: <scratch>/A5b/gate3.sh (tests+doc+clippy), scopes.py, dump.py (C++ -> <scratch>/A5b/cpp).
- Ledger status: IMDataArrayUtils, MetaInfo, MetaInfoInterface, SpectrumLookup -> complete (not yet in fragment).

## Next step
1. Finish gate3; commit SpectrumMetaDataLookup milestone.
2. MetaInfoDescription standalone + DataArrays source ordering (needs registry; name-only difference is
   unordered in source: defect candidate). Plan: metadata/meta_info_description.rs + DataArray::source_lt.
3. Ledger fragment update for all; SysInfo/StopWatch/BuildInfo stay partial (platform crates).
4. Final: CI sweep, report, delete NOTES.md.

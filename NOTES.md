# A3b NOTES (DEBT-03: decide every open item)

Deliverable: `docs/OPEN_ITEMS.md` (new). No code changes.
Inputs (main checkout only, gitignored .planning):
- /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.planning/OPEN-ITEMS.md (962 lines;
  §1.1-1.6 = OI-01..OI-116, §1.7 = ST-01..ST-12 stale, §3B stale refusals, §4 C++ defects, §5 DEC-1..DEC-18)
- .../.planning/REQUIREMENTS.md (180 lines) - M1 phase + requirement IDs for "scheduled"
- inventory was compiled against c3ac959; worktree base f0990c0.

Method: for each OI, grep current tree for the cited site; classify
resolved(commit) / scheduled(phase+REQ) / deferred(reason) / decision owed(question+options+costs).
Stale-record fixes: one-sentence edits only, never in lead-owned files (VALIDATION.md, PORTING_STATUS, README, CHANGELOG...).

Outputs also: .planning/lane-reports/A3b.md, .planning/ledger-fragments/wave4/A3b.json (likely empty entries: no headers).

## Next step
1. Read OPEN-ITEMS.md lines 219-962 and REQUIREMENTS.md, ROADMAP.md.
2. Verify OI-01..OI-116 against tree (batch greps), record verdicts in scratch file
   /private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad/verdicts.md
3. Write docs/OPEN_ITEMS.md, commit per section.

## Learned (incarnation 2)
- No Agent tool available: do it serially with batched greps.
- ROADMAP phases: A=DEBT/EVID, B=BLOCK/SHARED/CLI, C=CORE-01..05, D=CORE-06..12, E=CORE-13..18, F=STRUCT, G=UPSTREAM.
- Commits since c3ac959 that close items: 9466368 (D18 recorded in VALIDATION register -> OI-56),
  3942b51 (ST-01..04 fixed in FILE_INFO docs + VALIDATION), cf899a0 (citation item 7 closed -> OI-60),
  c242e71/6ffcb39 (cycle gate sees braced imports -> OI-69?), 1712a03 (native_equivalent def).
  eac9ac6/39fdeb6 (A8: mzXML/mzData/MGF/MS2 FileInfo -> OI-09 partly), 1a773e4 (CTD writer -> OI-37 partly),
  b0c9765 (XSD validation -> OI-30?), 2ec5beb mzXML.
- Verdicts accumulate below in "## Verdicts" (restart-safe; commit often).

## Incarnation 3 (merged f3f-6 FF at 22273d0)
Wave-4 sibling lanes (branch worktree-wf_ee239125-5be-N; their NOTES.md in .claude/worktrees/wf_ee239125-5be-N/):
1 B0 integrate Phase-B branches (FileHandler/LogStream/MSExperiment/TOPPBase/ProteinId) BLOCK-01..05;
2 B6 DefaultParamHandler trait SHARED-01; 3 A2a / 4 A2b EVID-02 native_equivalent re-audit;
5 A3a DEBT-04: range_utils(OI-07,B1), spectrum_mobility(B2), iterative(B3), faims(B5), ffp debug(B6),
  chromatogram_tools(B7), compute_rank(OI-01,B4), feature_overlap_filter stable_sort (OI-08) + oracle A3a-source-order;
7 A4 DEBT-05/EVID-03/EVID-04 (CPP index, citation tail -> OI-60?); 8 A5a / 9 A5b DEBT-01 partial re-read;
10 F1a Arrow/Parquet; 11 F2 MzTab/MzTab-M/mzQC/OMS (OI-28, OI-29); 12 F3a IdXML + modified-hit budget (OI-20);
13 F3b search-engine files; 14 F4 sqMass/CachedMzML/SwathFile/MSP/EDTA... (CORE-01).
Items worked by an in-flight lane = "scheduled" (Phase + REQ + lane name), not resolved.

## Verdicts
Header->req map: MzMLFile.h/StringUtils/MathFunctions -> B BLOCK-06; FileHandler B BLOCK-01; LogStream B BLOCK-02; MSExperiment B BLOCK-03 (branch 673154fa-1e1-2 closes OI-40); TOPPBase B BLOCK-05; FORMAT/* incl FileInfo.h, MzMLHandler.h -> C CORE-01; IMAGING -> E CORE-18.
Unmerged Phase B branches: 673154fa-1e1-{1 TOPPBase,2 MSExperiment,3 LogStream,4 FileHandler,5 ProteinIdentification}.
- OI-01 open (:1378 refuse, :1388 total_cmp) -> sched A DEBT-04 (leg a) + oracle row for leg b
- OI-02 open (UnsortedData :516/675/733/810) -> DECISION (C++ _sorted does not verify)
- OI-03/04 open (unmeasured, doc :528, :1128) -> sched A DEBT-04? no: one ibminode06 oracle row; MATH StatisticFunctions.h is native_equivalent -> EVID-02 re-audit. sched A EVID-02
- OI-05 open -> sched with OI-03 row
- OI-06 RESOLVED 6aed3c3 (VALIDATION fixed by 3942b51); WP:1901-1906 stale sentence -> fix
- OI-07 open range_utils.rs:721 total_cmp -> A DEBT-04 (B1)
- OI-08 open feature_overlap_filter.rs:1105 sort_by -> ?
- OI-09 PARTIAL: A8 first half landed (39fdeb6, eac9ac6: mzXML/mzData/MGF/MS2); left -v, pepXML/mzTab/trafoXML/PQP, sqMass/XMass/MSP/RAW/TDF -> C CORE-01
- OI-10 open -> DEC-4
- OI-11 open (ignores tests/file_info.rs:168,325,351; tripwire :1237) -> C CORE-01 MzMLHandler
- OI-12 open (ignore :494) decision taken -> C CORE-01
- OI-13 open: lib default strict is correct (Series-A D10); tool passes true (file_info.rs:205); test :210 uses default -> C CORE-01 (run case with option)
- OI-19 RESOLVED: -log e64ade6, timing fb0ae8e, writeDebug_ 287c89f; -instance reproduces source (TOPP_CLI_SUPPORT:400)

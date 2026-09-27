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
- OI-14 open (FILE_INFO_CHECKS nd 11/12 at :225,:248; IndexedMzMLDecoder.h native_equivalent) -> DEFERRED (reviewed departure; CPP-337 -> G UPSTREAM-01)
- OI-15/16/17 open (VALIDATION :1357, :1173, :933) -> DEFERRED to next a6/a7 oracle re-execution (cost stated)
- OI-18 open (TOPP_FILE_INFO :136-144) -> sched B BLOCK-01 (FileHandler, lane 673154fa-1e1-4 / B0)
- OI-20 open main; F3a lane (5bb4d01 WIP parse_metered) -> sched C CORE-01 (F3a)
- OI-21 open (MSSTATS :258-265) -> sched C CORE-01 (ConsensusXMLFile.h) ; shares cause with OI-20 (AASequence parse charge)
- OI-22 open -> DECISION? (relaxing column check) -> actually deferred: port reproduces source; upstream PERCIN-007 -> G
- OI-23 open (mzml.rs:2010 transcoding, :2270 "1.1." prefix) -> sched B BLOCK-06 MzMLFile.h / C CORE-01 ImzMLFile.h
- OI-24 open -> sched C CORE-01 ImzMLFile.h (partial)
- OI-25 IMAGING: MSImagingExperiment.h now COMPLETE, Geometry/Region/IonImage partial ("owned by later IMAGING pkg") -> E CORE-18 + DEC-8 still owed (where types live)
- OI-26 open -> sched E CORE-18 (MSImagingGeometry) + C CORE-01 ImzMLFile
- OI-27 open (QcMLFile partial) -> D CORE-10? QcMLFile is FORMAT -> C CORE-01; DECISION on defects (DEC-19 new?)
- OI-28/29 -> sched C CORE-01, lane F2 (MzTab, MzTab-M, OMS). MzTabMFile.h native_equivalent -> EVID-02
- OI-30 RESOLVED 3897138 (+4d59f17) xml_schema.rs bundles mzXML_idx_3.1.xsd
- OI-31/32 open -> sched C CORE-01 (MzDataFile.h partial)
- OI-33 open (transformation_xml.rs:175,:352) -> sched C CORE-05 (TransformationModelBSpline unmapped, Interpolated partial)
- OI-34 open; reason "no regex engine" STALE (concept::boost_regex since 854d996 2026-09-13; MASCOT_XML_SUPPORT.md:190 written 09-12) -> fix sentence; sched D CORE-08
- OI-35 open (IndexedMzMLFileLoader partial) -> sched C CORE-01 + E CORE-18 INTERFACES
- OI-36 open (RMT evidence_requires_review) -> sched A DEBT-02 / C CORE-02
- OI-37 RESOLVED 1a773e4 (CTD writer) + c36d087; CWL/JSON refused exactly as Release w/o TDL (src/cli.rs:1169-1174)
- OI-38 open (src/cli.rs:2008-2011) -> sched B BLOCK-01
- OI-39 open -> sched E CORE-18 (MSSpectrum.h partial)
- OI-40 -> sched B BLOCK-03 (branch 673154fa-1e1-2, B0 integrating)
- OI-41 open -> E CORE-18 IONMOBILITY; waits SpectraMerger etc
- OI-43 open -> DECISION (keep) -> sched B BLOCK-06 MathFunctions.h
- KERNEL/DimMapper.h unmapped & counted in CORE-18 though OI-83 says TOPPView excluded -> DECISION owed (M1 def vs exclusion)

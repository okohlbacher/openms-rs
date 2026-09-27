# Lane A2b NOTES (EVID-02 second half: native_equivalent re-audit M-Z)

Ledger source: docs/core-sdk-reviewed-apis.json (lead's file, do NOT edit). Output goes to
../../../.planning/ledger-fragments/wave4/A2b.json and .planning/lane-reports/A2b.md in MAIN repo
(/Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.planning/...).

19 native_equivalent headers at/after M (path after include/OpenMS/):
MATH/MISC/EmgGradientDescent.h, MATH/STATISTICS/{BasicStatistics,Histogram,RankData}.h,
MATH/StatisticFunctions.h, METADATA/{ContactPerson,DocumentIdentifier,ExperimentalSettings,
Gradient,HPLC,IdentifierMSRunMapper,Instrument,IonDetector,IonSource,MassAnalyzer,
PeptideEvidence,PeptideIdentificationList,Sample}.h, SYSTEM/CurlInit.h

Rule (docs/NATIVE_EQUIVALENTS.md, 2026-09-21): native_equivalent = Rust lang + std discharge
the header, no port owed, no crate. Any hand-written port => complete/partial on member coverage.
Genuine candidates: PeptideIdentificationList.h (Vec alias). CurlInit.h: no Rust type, but
ureq crate replaces curl -> a crate discharges it => not native (decide: complete w/ not-ported rationale?).

Next step:
1. For each header, check C++ header members (git -C ../OpenMS4-tests/packages/core show bc9cc12:src/openms/include/OpenMS/<path>)
   vs Rust file; decide complete/partial.
2. Write ledger fragment + report; commit per group.

## Incarnation 2 (2026-09-27)
Merged worktree-wf_bcdd783d-f3f-4 (fast-forward, NOTES only). Current ledger entries dumped from
docs/core-sdk-reviewed-apis.json 'headers' (19 NE at/after M confirmed).
Plan: per header read C++ header at bc9cc12 + Rust file + existing *_SUPPORT.md API table;
re-record complete/partial. Keep PeptideIdentificationList NE if header is a pure alias.
StatisticFunctions: its own scope says a NaN gap in computeRank/Spearman -> partial (named gap).
Next step: 3. check tools/core_sdk_coverage.py rules for complete (rust [] allowed?), then header-by-header.
- All 19 headers + cpp + class tests are byte-identical between audit pins (82ce5b3/6bfc0e4) and bc9cc12
  (git log bc9cc12 -- paths: last change ea35651/70c8251, both ancestors of the audit pins).
- Verdicts so far (member-by-member vs Rust): Sample/IonSource/IonDetector/MassAnalyzer/Instrument complete
  (experiment_values.rs fields+enums+Hash field selection match); ContactPerson/HPLC/Gradient complete
  (experiment_support.rs); DocumentIdentifier complete (setLoadedFileType delegate parity -> FileHandler.h);
  ExperimentalSettings complete (fields, eq incl. DocId-id-only, Display = 2 marker lines);
  IdentifierMSRunMapper complete (13 members); PeptideEvidence complete (derived Ord = source order, CPP-045 fix);
  PeptideIdentificationList stays NE (pure ExposedVector alias + 3 ctor/3 assign conversions).
  CurlInit -> complete (discharged by ureq crate choice, not std; THIRD_PARTY_CRATE_DECISIONS lacks a
  libcurl->ureq transport row: report it).
- StatisticFunctions -> partial (own scope names NaN gap in computeRank/Spearman).
- Git guard: use one plain `git -C <abs core> show bc9cc12:<path> > <scratch file>` per call; no loops/vars.
Next step: 4. MATH headers (Emg, BasicStatistics, Histogram 'resize', RankData helpers) then write fragment+report.
- DONE: ledger fragment written (.planning/ledger-fragments/wave4/A2b.json, 16 complete / 2 partial / 1 NE) by
  scratchpad make_fragment.py. EmgGradientDescent -> partial (DefaultParamHandler Param surface gap, B6 trait closes it).
- Gate run of all cited tests on dax started (log: scratchpad/gate_cited.log).
Next step: 5. wait for gate pass counts, write .planning/lane-reports/A2b.md, delete NOTES.md in final commit.

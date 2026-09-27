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

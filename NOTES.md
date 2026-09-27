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

## Verdicts

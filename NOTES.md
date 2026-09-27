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

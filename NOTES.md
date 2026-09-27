# Lane A4 notes (DEBT-05, EVID-03, EVID-04)

Worktree off main f0990c0, branch worktree-wf_bcdd783d-f3f-7.

## State found
- OpenMS_CPP_ISSUES.md: 6601 lines, 354 `## CPP-` headings; index table covers
  CPP-001..058 only; header pins `82ce5b3` (claimed unused by entries).
- tools/check_source_citations.py (1135 lines): exit 0, reports
  "4136 citations checked ... 35 skipped for an unreachable file and 3 bare ranges".
  (brief said 31 unreachable -> now 35).
- docs/REPRODUCIBILITY.md does not exist.

## Next step
1. Parse every CPP entry (id, title, affected revision, evidence line, Rust handling)
   with a script in scratchpad; generate index table; check which revisions are used.
2. Fix header of OpenMS_CPP_ISSUES.md; replace index; commit.
3. Checker: list the 35 unreachable + 3 bare ranges (verbose flag?), resolve each.
4. Measure sha256 recomputability for REPRODUCIBILITY.md (check_core_sdk.py).
5. Lane report + ledger fragment; delete NOTES.md.

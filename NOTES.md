# Lane B0 notes — integrate and verify the five Phase-B branches

Branch: worktree-wf_bcdd783d-f3f-1 (off main f0990c0).

## Branches to merge, in this order
1. worktree-wf_673154fa-1e1-1  TOPPBase (JSON INI, exit 6->8, citations, accessor misuse)  tip 2ce5bf0
2. worktree-wf_673154fa-1e1-2  MSExperiment (set2DData mass traces, getFirstProductSpectrum, operator<<) tip 6203012
3. worktree-wf_673154fa-1e1-3  LogStream (process-exit flush, Windows VT, CPP-355) tip b4ab801
4. worktree-wf_673154fa-1e1-4  FileHandler (log type, rewrite_source_file, compute_hash, SRM) tip d37eab0
5. worktree-wf_673154fa-1e1-5  ProteinIdentification member review tip 6e544f9

## Status
- [ ] merges 1..5
- [ ] full sweep of merged tree (dax, slot B0-sweep)
- [ ] attacks: exit codes, JSON float round-trip, compute_hash/rewrite_source_file bytes, exit flush paths, MSExperiment operator<<, PI engine ordering
- [ ] FileInfo vs FileHandler check_allowed duplicate
- [ ] ledger fragment consolidation -> .planning/ledger-fragments/wave4/B0.json
- [ ] report -> main repo .planning/lane-reports/B0.md

## Next step
1. Merge branch 1.

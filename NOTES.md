# Lane B0 notes — integrate and verify the five Phase-B branches

Branch: worktree-wf_ee239125-5be-1 (off main f0990c0), merged worktree-wf_bcdd783d-f3f-1 (earlier incarnation).

## Branches merged (done by earlier incarnation, commits 9e52998..89cde73)
1. -1e1-1 TOPPBase tip 2ce5bf0
2. -1e1-2 MSExperiment tip 6203012
3. -1e1-3 LogStream tip b4ab801
4. -1e1-4 FileHandler tip d37eab0
5. -1e1-5 ProteinIdentification tip 6e544f9
No conflict markers remain.

## Status
- [x] merges 1..5
- [ ] full sweep of merged tree (dax, slot B0-sweep)
- [ ] attacks: exit codes, JSON float round-trip, compute_hash/rewrite_source_file bytes, exit flush paths, MSExperiment operator<<, PI engine ordering
- [ ] FileInfo vs FileHandler check_allowed duplicate
- [ ] ledger fragment consolidation -> .planning/ledger-fragments/wave4/B0.json
- [ ] report -> main repo .planning/lane-reports/B0.md

## Next step
1. Launch full sweep in background (log dir in scratchpad), then read lane fragments + reports.

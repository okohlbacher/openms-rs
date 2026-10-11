# INTEG-F working notes (integrator, group F, worktree wf_97065d91-48e-4)

Base: main f0990c0. F1c is stacked on F1a (contains 8bb62562). All lanes ahead of base.

Lane branches (merge in this order; F1a then F1c is a fast chain):
F1a worktree-wf_cd066dab-0a7-5 | F1b worktree-wf_cd066dab-0a7-42 | F1c worktree-wf_cd066dab-0a7-43
F2 worktree-wf_a10fdd71-2d2-1 | F3a worktree-wf_cd066dab-0a7-7 | F3b worktree-wf_a10fdd71-2d2-2
F4 worktree-wf_cd066dab-0a7-9 | F5 worktree-wf_a10fdd71-2d2-3 | F6 worktree-wf_cd066dab-0a7-11
F7a worktree-wf_a10fdd71-2d2-4 | F7b worktree-wf_a10fdd71-2d2-5 | F7c worktree-wf_a10fdd71-2d2-6

Reports: .planning/lane-reports/<lane>.md, fragments .planning/ledger-fragments/wave4/<lane>.json
(read from the main checkout /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R, or from the lane branch).

## Progress
- [ ] merge all 12 branches (commit each merge)
- [ ] dax sweep: OPENMS_GATE_HOST=dax OPENMS_GATE_JOBS=96 ~/.local/bin/openms-ci-sweep.sh integ-F <logdir>  (background, poll)
- [ ] collapse duplicated logic (shared Parquet writer helpers F1a/F1b/F1c; gzip/zip F7b vs others; mzML writer F7a vs F4 cached)
- [ ] member-by-member check of every `complete`
- [ ] consolidated INTEG-F.json fragment + INTEG-F.md report
- [ ] delete NOTES.md in the final commit

## Next step
1. Merge F1a, F1c, F1b ... then re-read each report.

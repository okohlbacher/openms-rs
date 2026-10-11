# INTEG-X notes (group X integrator)

Branch: worktree-wf_97065d91-48e-10 off main f0990c06. Never push, never rebase, no local cargo (dax via gate).
.planning/ is NOT tracked: reports/fragments live at
/Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.planning/{lane-reports,ledger-fragments/wave4}/

## Lane branches (merge order = list order)
I1 worktree-wf_cd066dab-0a7-28 | I2 -24 | I3a -29 | I3b -31 | I4 -32 | FF1 -33 | FF2a -34
FF2b worktree-wf_a10fdd71-2d2-9 | FF3 worktree-wf_97065d91-48e-8 | FF4 worktree-wf_cd066dab-0a7-35
M1 worktree-wf_cd066dab-0a7-36 | M2 worktree-wf_97065d91-48e-2   (all prefixes worktree-wf_cd066dab-0a7- unless stated)

## Next step
1. merge each branch (git merge --no-ff), commit each merge; resolve conflicts.
2. run sweep on dax: OPENMS_GATE_HOST=dax OPENMS_GATE_JOBS=96 ~/.local/bin/openms-ci-sweep.sh integ-X <logdir> (background).
3. attack complete claims; collapse duplicated logic; fix; commit.
4. write INTEG-X.json + INTEG-X.md in main repo .planning/.  Delete this file in the final commit.

## Progress log
- (none yet)

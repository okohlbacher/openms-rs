# INTEG-B notes (integrator group B, wave 4)

Branch: worktree-wf_a10fdd71-2d2-11, base main f0990c06.

Lanes to merge (commits ahead of main):
- B0  worktree-wf_ee239125-5be-1   (112)
- B6  worktree-wf_cd066dab-0a7-1   (15)
- B7  worktree-wf_cd066dab-0a7-37  (124, contains B0 + -33)
- B8a worktree-wf_cd066dab-0a7-38  (127)
- B8b worktree-wf_cd066dab-0a7-39  (132)
- B8c worktree-wf_cd066dab-0a7-40  (124, contains B0 + -36)

Constraints: no local cargo (dax gate slot integ-B, sweep slot integ-B-sweep);
git commands must be plain single commands (worktree isolation guard).

## Next step
1. Merge B0, then B6, B7, B8a, B8b, B8c; commit each merge.

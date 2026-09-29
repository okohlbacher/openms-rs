# INTEG-F notes (integrator, group F)

Base: main f0990c06. Branch: worktree-wf_a10fdd71-2d2-14.

## Lane branches (tips at start)
- F1a  worktree-wf_cd066dab-0a7-5   8bb62562
- OA   worktree-wf_cd066dab-0a7-27  ce18035b (F1b already merged it; must come in with F1b)
- F1b  worktree-wf_cd066dab-0a7-42  c73aa407 (contains F1a + OA)
- F1c  worktree-wf_cd066dab-0a7-43  14dd8d68
- F2   worktree-wf_a10fdd71-2d2-1   75891deb
- F3a  worktree-wf_cd066dab-0a7-7   39f6fde9
- F3b  worktree-wf_a10fdd71-2d2-2   242a24f3
- F4   worktree-wf_cd066dab-0a7-9   8a60f6e3
- F5   worktree-wf_a10fdd71-2d2-3   1dc4a0ea
- F6   worktree-wf_cd066dab-0a7-11  726720ce (sweep NOT run by lane)
- F7a  worktree-wf_a10fdd71-2d2-4   7ecc49e8
- F7b  worktree-wf_a10fdd71-2d2-5   0756c4e1
- F7c  worktree-wf_a10fdd71-2d2-6   b9dd94c4

## Plan / merge order
F1a, F1b(+OA), F1c, F2, F3a, F3b, F4, F5, F6, F7a, F7b, F7c.

## Status
(none merged yet)

## Next step
1. Merge lanes in the order above, commit each merge.
2. Look for duplicates: string_utils (F3b + F4 both?), line_reader, PTMXMLFile (F3b + F6 both!), zip input (F7b vs B8a).
3. Run sweep on dax slot integ-F (background).
4. Attack complete claims.

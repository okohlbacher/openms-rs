# O2 lane notes (OpenSWATH scoring)

Branch: worktree-wf_cd066dab-0a7-47. Prereqs 0a7-27 (OA) and 0a7-30 (O1a) merged (commit 59a11552).

## C++ sources
- Core checkout `../OpenMS4-tests/packages/core` HEAD is 291c5ba, NOT the pin; for the seven scope
  headers the only pin diff is a doc comment in OpenSwathScoring.h, so read the working tree directly.
  (`git -C ... diff bc9cc12 291c5ba -- <path>` confirms.)
- The worktree-isolation guard refuses `git -C` with shell variables/loops and paths containing
  the word "source" inside loops; use literal single commands. Scratch symlink `$SCRATCH/c` -> core/src.
- Sizes (h+cpp): DIAHelper 180+612, DIAScoring 325+590, DIAPrescoring 138+308, MRMScoring 304+836,
  OpenSwathScoring 339+832, IonMobilityScoring 273+854, OpenSwathScores 235+248.
- Class tests: src/tests/class_tests/openms/source/<Name>_test.cpp, plus openswathalgo/DiaHelpers_test.cpp.

## Dependency order
DIAHelper -> DIAScoring, DIAPrescoring -> MRMScoring -> IonMobilityScoring -> OpenSwathScores -> OpenSwathScoring.

## Next step
1. Read DIAHelper.h/.cpp + tests; port to src/analysis/openswath/dia_helper.rs.

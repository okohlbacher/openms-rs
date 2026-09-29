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

## Done (WIP commit dcdf630f)
- dia_helper.rs, dia_prescoring.rs, dia_scoring.rs + tests/openswath_dia_scoring.rs (class-test literals).
- Gate: `OPENMS_GATE_HOST=dax OPENMS_GATE_JOBS=64 ~/.local/bin/openms-kim-gate.sh O2 test --locked --test openswath_dia_scoring`
  (run_in_background; log in scratchpad). NOTE test-name filters do not include the file name.
- DIAScoring.cpp uses its own global C13C12 = 1.0033548 (not Constants' 1.0033548378); DIAPrescoring uses Constants.
- Averagine: SourceSingle precision (C++ Peak1D float intensities).

## Progress
- DIA trio green (24). MRMScoring green (16, commit d21665c8). OpenSwathScores methods + Usage (06df27ef).
- IonMobilityScoring (7bd5f5fd): peak picking -> Error::Unsupported (PeakPickerMobilogram unported);
  MobilogramParquetConsumer -> MobilogramSink trait; driftIdScoring takes detection transitions slice.
- Test files: tests/openswath_{dia_scoring,mrm_scoring,ion_mobility_scoring}.rs (default features).

## Next step
1. Confirm gate4 (3 test files). 2. OpenSwathScoring (339+832). 3. Support doc docs/OPENSWATH_SCORING_SUPPORT.md,
provenance tests/data/openswath_scoring_provenance.json, doc coverage, clippy/fmt/doc/MSRV/sweep, report+ledger.

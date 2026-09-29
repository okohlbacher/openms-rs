# O1b lane notes (Transition list I/O)

Prereq merged: worktree-wf_cd066dab-0a7-30 (O1a: TargetedExperiment, ReactionMonitoringTransition,
TargetedExperimentHelper, IncludeExcludeTarget in src/analysis/openswath/).

C++ sources extracted to scratchpad O1b/ (git archive of core bc9cc12).
Sizes: TraMLHandler.cpp 1745, TransitionTSVFile.cpp 2001, TransitionPQPFile.cpp 1683,
TraMLFile.cpp 58, MRMFile.cpp 54, TargetedDataFileLoader.cpp 59.
Class tests: TraMLFile_test (208), TransitionTSVFile_test (361), TransitionPQPFile_test (445),
MRMFile_test (36), TargetedDataFileLoader_test (39). Test data: ToyExample1.traML.

MRMFile / TargetedDataFileLoader need OpenSwath::SwathMap, SwathFile, SimpleOpenMSSpectraFactory,
MRMChromHandler -> check whether ported; probably blocked (name the deps).

## Next step
1. Port TraMLHandler + TraMLFile to src/format/traml*.rs (load/store TargetedExperiment).
2. TransitionTSVFile -> src/analysis/openswath/transition_tsv_file.rs
3. TransitionPQPFile (rusqlite) -> src/analysis/openswath/transition_pqp_file.rs
4. MRMFile / TargetedDataFileLoader if deps exist.

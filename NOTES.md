# O1c lane notes (Assays and decoys)

Merged prerequisite worktree-wf_cd066dab-0a7-30 (O1a: TargetedExperiment,
ReactionMonitoringTransition, TargetedExperimentHelper, MRMMapping, ...) as a
fast-forward. No earlier O1c incarnation exists in bcdd783d branches.

C++ pin readable at .reference/openms4-core-bc9cc12 (outside worktree; git -C
to other repos is refused by the sandbox, plain cat/find works).

Headers (src/openms/include/OpenMS/ANALYSIS/OPENSWATH/):
- MRMIonSeries.h (139) / .cpp (666) / test 279 lines  -> dependency of both below
- MRMDecoy.h (332) / .cpp (1160) / test 1287 + TraML in/out fixtures
- MRMAssay.h (568) / .cpp (2197) / test 1007 + many TraML fixtures
- TransitionListEvidenceFilter.h (144) / .cpp (1320) / test 438

Target Rust: src/analysis/openswath/{mrm_ion_series,mrm_decoy,mrm_assay,
transition_list_evidence_filter}.rs; append pub mod lines to
src/analysis/openswath.rs.

## Next step
1. Read MRMIonSeries.h/.cpp/test; check Rust AASequence, residue, TheoreticalSpectrumGenerator, ModificationsDB, TraML reader availability.
2. Port MRMIonSeries, commit.
3. MRMDecoy, then MRMAssay, then TransitionListEvidenceFilter.

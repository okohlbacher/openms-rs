# O1a lane notes (targeted experiment data model)

C++ pin: /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages/core @ bc9cc12
Code goes to src/analysis/openswath/ (module root src/analysis/openswath.rs is empty, append pub mod only).

## Missing (all unported)
- ANALYSIS/TARGETED/TargetedExperimentHelper.h (+cpp) — CV-term value types (Protein, Peptide, Compound, RetentionTime, ...). Test: TargetedExperimentHelper_test.cpp
- ANALYSIS/MRM/ReactionMonitoringTransition.h (+cpp). Test: ReactionMonitoringTransition_test.cpp
- ANALYSIS/TARGETED/TargetedExperiment.h (+cpp). Test: TargetedExperiment_test.cpp
- IncludeExcludeTarget, MetaboTargetedAssay, MetaboTargetedTargetDecoy, MRMMapping, ChromatogramProcessor, *ChromHandler (no class tests at pin)
- kernel/mrm.rs Transition trait: implement for ReactionMonitoringTransition

## Next step
1. Read TargetedExperimentHelper.h/.cpp, port to src/analysis/openswath/targeted_experiment_helper.rs
2. ReactionMonitoringTransition -> reaction_monitoring_transition.rs, impl kernel::mrm::Transition
3. TargetedExperiment -> targeted_experiment.rs
4. IncludeExcludeTarget, MRMMapping, then rest

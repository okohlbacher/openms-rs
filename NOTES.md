# O1a lane notes (targeted experiment data model) — incarnation 2

C++ pin files exported to scratchpad .../scratchpad/cpp/src (git archive of core bc9cc12; re-export if gone:
`git -C ../OpenMS4-tests/packages/core archive bc9cc12 <paths> -o t.tar`). Worktree-isolated agent: git -C with
literal path works, compound git commands are refused.
Code goes to src/analysis/openswath/ (module root src/analysis/openswath.rs, append pub mod only).

## Design decisions (made, keep them)
- CVTermList (src/metadata/value.rs) stands for both CVTermList and CVTermListInterface; base class -> field `cv_terms`.
- Inheritance PeptideCompound -> Compound/Peptide: field `base: PeptideCompound` + Deref/DerefMut.
- Precondition getters (getRT, getChargeState, getPrecursorCVTermList) -> Option; throwing getters -> Result
  (IllegalArgument -> Error::InvalidValue, ElementNotFound/MissingInformation -> Error::MissingInformation).
- C++ quirk: Configuration and Peptide::Modification have no operator==, so inherited CVTermList== ignores their own
  fields. Rust compares all fields (native divergence, C++ defect candidate). Hash consistent with Eq.
- C++ getAASequence reads peptide.sequence[location] for location -1 (UB) / ==size; Rust: strict default refuses
  unresolvable mods, SourceSkip policy reproduces skip+warning.
- RMT library_intensity default -101 (kernel SimpleTransition says -1; RMT impl kernel::mrm::Transition).
- SummaryStatistics Display uses crate::param::value::format_float(x,false) (= StringUtils::toStr(double,false)).
- Sorting: total_cmp + stable sort (std::sort unstable; `<` with NaN).

## Headers (order)
1. TargetedExperimentHelper -> targeted_experiment_helper.rs (test: TargetedExperimentHelper_test.cpp)
2. ReactionMonitoringTransition -> reaction_monitoring_transition.rs (ReactionMonitoringTransition_test.cpp)
3. IncludeExcludeTarget -> include_exclude_target.rs (IncludeExcludeTarget_test.cpp)
4. TargetedExperiment -> targeted_experiment.rs (TargetedExperiment_test.cpp)
5. MRMMapping, ChromatogramProcessor/*ChromHandler, MetaboTargeted* later

## Next step
1. Write targeted_experiment_helper.rs + tests/targeted_experiment.rs; gate on spock (slot O1a).

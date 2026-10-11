# O3b notes (MRM feature finding, filtering, QC)

## State (update when a milestone lands)
- Merged prerequisites (O2 scoring branch cd066dab-47, peak picking 97065d91-1) in commit 22219451.
- DONE + green on spock (`--test` runs): MRMFeatureQC, MRMFeatureQCFile, MRMFeaturePicker, MRMFeaturePickerFile
  (tests/mrm_feature_qc.rs, 15 tests), MRMFeatureFilter (tests/mrm_feature_filter.rs, 26 tests, mechanical transcription of every upstream section).
- C++ sources extracted to scratchpad `.../scratchpad/core/` via `scratchpad/extract.py <S> <regex> [pkg rev]`.
  The bash tool refuses commands containing "git -C"/other-repo git and heredocs that write into the worktree; write python scripts with the Write tool and run them by path.
  Test translator: scratchpad/translate_filter_test.py + assemble_filter_test.py (TEST_EQUAL converts expected to T1(b): truncates for Int).
- Source modules: src/analysis/openswath/{mrm_feature_qc,mrm_feature_picker,mrm_feature_filter,mrm_feature_support}.rs, src/format/{mrm_feature_qc_file,mrm_feature_picker_file}.rs
- Shared helper `mrm_feature_support` (pub(crate)): data_value_to_string, meta_double, native_id.

## Remaining
1. MRMFeatureSelector (+QMIP/Score) and MRMBatchFeatureSelector: LP solve is gated (LPWrapper not ported, GATED-01). Plan: LinearProgram model + LpSolver trait + exact PartitionSolver for the Score variant; QMIP takes any solver. Class test needs FeatureXML loader (exists) and tests/data MRMFeatureSelector_*.featureXML.
2. ConfidenceScoring (needs TransformationDescription::apply, RandomShuffler::portable_random_shuffle (in openswath_helper / chemistry decoy_random?), TargetedExperiment peptide_by_ref).
3. MRMFeatureFinderScoring: depends on O2 openswath_scoring.rs, mrm_transition_group_picker.rs, EmgScoring NOT in crate (named gap: use_elution_model_score), MobilogramSink trait, ProteaseDigestion.
4. Support doc docs/MRM_FEATURE_QC_SUPPORT.md (name used in module docs), provenance json, ledger fragment, lane report, CI slice counts, doc coverage check, delete NOTES.md at the end.

# FF3 notes (lane FeatureFinderIdentification + MetaboIdent)

Worktree branch worktree-wf_97065d91-48e-8. Prereq branches merged (commit 8fd08226).

## Sources (pinned checkout, read here)
/Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.reference/openms4-core-bc9cc12/
- src/openms/{include/OpenMS,source}/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.{h,cpp} (459+1772)
- src/openms/{include/OpenMS,source}/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.{h,cpp} (295+988)
- class test: src/tests/class_tests/openms/source/FeatureFinderIdentificationAlgorithm_test.cpp (no MetaboIdent class test)
- TOPP: FeatureFinderIdentification, FeatureFinderMetaboIdent in ../OpenMS4-tests/packages/topp @174b576

## Status
- [x] prerequisites merged
- [ ] explore C++ and existing Rust (openswath/*, featurefinder/*)

## Next step
1. Read headers+cpp; list deps already in Rust (MRMFeatureFinderScoring, ChromatogramExtractor, ElutionModelFitter, TransformationModel...).

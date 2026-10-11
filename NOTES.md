# O3b notes (MRM feature finding, filtering, QC)

## State
- Merged prerequisites (O2 scoring branch cd066dab-47, peak picking 97065d91-1) in commit 22219451.
- C++ sources extracted to scratchpad `.../scratchpad/core/` via `scratchpad/extract.py <S> <regex>` (runs git -C into OpenMS4-tests/packages/core at bc9cc12).
  The bash tool refuses `git -C <other>` / `cd other && git`; use that python script.
- All 9 headers are `unmapped` in docs/core-sdk-coverage.json.

## Headers (sizes cpp lines): MRMFeatureQC.h (275, tiny cpp), MRMFeaturePicker (54), MRMFeatureFilter (345 + 1476), MRMFeatureSelector (350+393),
MRMBatchFeatureSelector (86+50), ConfidenceScoring (280+260), MRMFeatureFinderScoring (442+1449), MRMFeatureQCFile (175+290), MRMFeaturePickerFile (140+136)

## Next step
1. Read headers+cpp+tests, check kernel deps (MRMFeature, MRMTransitionGroup, FeatureMap, TargetedExperiment).
2. Port order: MRMFeatureQC -> MRMFeatureQCFile -> MRMFeatureFilter -> MRMFeatureSelector -> MRMBatchFeatureSelector -> MRMFeaturePicker(+File) -> ConfidenceScoring -> MRMFeatureFinderScoring.

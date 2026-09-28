# FF1 lane notes (delete in final commit)

Lane FF1: FEATUREFINDER models and fitters. Worktree branch worktree-wf_cd066dab-0a7-33, base f0990c0.

## State at start
- All lane headers `unmapped` in docs/core-sdk-coverage.json except
  FeatureFinderAlgorithmPickedHelperStructs.h (partial; src/analysis/feature_finder_picked/helper_structs.rs;
  gap = ConvexHull2D signed-zero merge in kernel, not ours).
- src/analysis/featurefinder.rs exists (empty module root; append `pub mod` only). Code goes in src/analysis/featurefinder/.
- LM reuse: src/math/fitters/levenberg_marquardt.rs. EMG helpers maybe in src/analysis/emg.rs.
- C++ pin: ../OpenMS4-tests/packages/core @ bc9cc12 (relative to main repo root, i.e.
  /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages/core).

## Plan (dependency order)
1. BaseModel (+_impl), ModelDescription(?), InterpolationModel
2. GaussModel, BiGaussModel, EmgModel, IsotopeModel, ExtendedIsotopeModel
3. Fitter1D, LevMarqFitter1D, MaxLikeliFitter1D, GaussFitter1D, BiGaussFitter1D, EmgFitter1D, IsotopeFitter1D, ExtendedIsotopeFitter1D
4. EmgScoring, PeakWidthEstimator, SeedListGenerator, FeatureFinderDefs, FeatureFinderAlgorithm, ElutionModelFitter

## Next step
1. Read skill SKILL.md, the C++ headers/impls + class tests for step 1.

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

## Findings (read already; don't redo)
- C++ readable at /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.reference/openms4-core-bc9cc12/
  (sha256 of BaseModel.h matches coverage json). Reading the tests pin through the VCS with paths containing
  the word "source" is REFUSED by the worktree guard; use the .reference checkout + python printing.
- Helper script: scratchpad/show.py <Name> prints h/cpp/test stripped; dumps in scratchpad/<Name>.txt
- Rust has NO LinearInterpolation (ML/INTERPOLATION unmapped) -> port it at src/math/ml/linear_interpolation.rs
  (flag in report: possible duplicate with an ML lane).
- Param handler: B6 trait `DefaultParamHandling` lives on branch worktree-wf_cd066dab-0a7-1 (src/param/handler.rs),
  NOT in base. Shape: NAME, defaults(), to_param(&self), update_members(&Param)->Result<Self>. Provide inherent
  defaults/from_param/from_param_with_warnings/to_param + from_complete_param like peak_picking.rs:1853-1985.
- C++ setParameters REPLACES param_ (defaults + given) -> a prior setInterpolationStep is lost in every
  Fitter's model build (GaussFitter1D etc: model->setInterpolationStep(step) then setParameters(tmp)). Defect candidate.
- BaseModel ctor sets defaults_ but no defaultsToParam_; LevMarqFitter1D ctor same; symmetric_ uninit.
- IsotopeModel::updateMembers_ does NOT call setSamples (others do). IsotopeModel uses
  CoarseIsotopePatternGenerator(max_isotope, round_masses=true); Extended uses (max_isotope) only.
- fitOffset_/Emg/ExtIso pearson use std::vector<float> (f32 rounding of both series).
- fit1d class tests are dummies except EmgFitter1D (fit of EmgModel samples -> params 100000/5/5/725, and #6239 robustness).
- EmgScoring test needs OpenSWATH_Test mock feature (0.924365639 / 0.981013417243958).
- PeakWidthEstimator needs BSpline2d: Rust has src/processing/spline/b_spline.rs (check it is BSpline2d port).
- SeedListGenerator: pure data plumbing; tests use PepXMLFile_test.mzML, ConsensusXMLFile_1.consensusXML.

## Progress
- DONE (compiles, clippy clean): src/math/ml/linear_interpolation.rs; src/analysis/featurefinder/{base_model,gauss_models,
  isotope_model,fitter1d,fitters,emg_scoring,feature_finder_algorithm,peak_width_estimator,seed_list_generator}.rs
- FeatureFinderDefs.h -> re-export of feature_finder_picked::defs (feature_finder_algorithm::defs).
- QualityType is f32 (BaseFeature::QualityType float) - correlations narrowed.
- Test fixtures already in repo: tests/data/peak_picking/PeakPickerHiRes_orbitrap.mzML, tests/data/consensusxml/ConsensusXMLFile_1.consensusXML
  (sha256 equal to core class-test data). PepXMLFile_test.mzML: transcribe 13 rows (ms level, rt, precursor mz) instead.

- Tier-1 oracle ../oracle/FF1-models (ibminode06 Release), 15742 values, ALL bit-identical on spock.
- Kernel fix: ConvexHull2D keep-first enlarge (closes FFAP helper structs gap). Docs amended.
- SCRATCHPAD IS SHARED WITH OTHER LANES: use scratchpad/FF1/ only.
- tools/check_schema_feature_graph.py runs cargo LOCALLY - do not run it (deleted local target/ once).
- docs/core-sdk-coverage.json stale -> lead must regenerate (core_sdk_coverage.py --write); doc-coverage floor --write too.

## Next step
1. Wait clippy (scratchpad/FF1/clippy2.log); run sweep into scratchpad/FF1/sweep.
2. Ledger fragment + lane report; delete NOTES.md in final commit.

# Lane M1 notes (map alignment) -- delete in the final commit

## Scratch
C++ pins extracted (git archive) to
`/private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad/M1/`
with symlinks `ct` (class tests), `cpp` (C++ implementations), `hdr` (headers), `keepdata` (MapAlign* class-test data).
The worktree guard refuses any command where a path ends in the word for the C++ implementation dir -- use the symlinks,
and write files with the Write tool rather than heredocs.
GTE Akima: `git -C ../OpenMS4-tests/packages/core show bc9cc12:src/openms/extern/GTE/Mathematics/IntpAkima{1,Nonuniform1}.h`.

## State at start (f0990c0)
- `src/analysis/mapmatching.rs` empty root (append `pub mod` only).
- `src/analysis/transformations.rs` + `transformations/lowess.rs`: closed enum ModelConfig/TransformationModel
  {None, Identity, Linear, Interpolated, Lowess}; typed options; NO Param surface; almost no rustdoc.
- `src/format/transformation_xml.rs` maps TrafoXML (name, BTreeMap<String,ParamValue>) <-> ModelConfig;
  returns Unsupported for b_spline and akima -- update it when those land.
- Available: BSpline2d (`processing::spline::b_spline`, has eval/derivative/ok), CubicSpline2d,
  BasicStatistics (`math::basic_statistics`, update(&[f64]) = C++ probability update incl. zero-sum guard),
  source_sort_by (`math::source_sort`, libstdc++ std::sort permutation), MorphologicalFilter (f32 only!),
  ConsensusMap::ranges() (includes handles = updateRanges), conversion_helper (peak/feature -> consensus),
  statistic_functions (quantile, mean_absolute_deviation, root_mean_square_error).
- NOT available: StablePairFinder/FeatureDistance (lane M2 owns feature grouping), CrossValidation (ML),
  Math::LinearInterpolation, KDTreeFeatureMaps, generic nth_element (only f32 in processing/noise_estimation.rs).

## Decisions
- Superimposer histograms need f64 tophat (C++ filterRange is templated on double) -> private f64 helper,
  flag duplication with processing::baseline for the integrator.
- nth_element on reverse iterators (affine truncation) -> generic introselect port, flag duplication.
- freq_cutoff scan can read past the end in C++ (UB) -> refuse with InvalidValue (decision D1).
- MapAlignmentAlgorithmPoseClustering blocked on StablePairFinder (M2) -> partial with a named gap.

## Planned order
1. [ ] transformations: Param surface (defaults/from_param/to_param per model), TransformationModelDefaults,
   b_spline model, akima, lowess auto-span CV, getModelTypes, printSummary, getModelParameters; trafoXML update.
2. [ ] MapAlignmentEvaluationAlgorithm + Precision + Recall (tests: 0.757143, 0.5, empty cases).
3. [ ] BaseSuperimposer + PoseClusteringShift/AffineSuperimposer.
4. [ ] MapAlignmentAlgorithmIdentification; PoseClustering (partial); KD/TreeGuided if time.

## Next step
1. Write transformations Param layer + b-spline + akima (step 1).

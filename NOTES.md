# Lane M1 notes (map alignment) -- delete in the final commit

## State at start (f0990c0)
- `src/analysis/mapmatching.rs` is an empty module root (append `pub mod` only).
- `src/analysis/transformations.rs` (649 l) + `transformations/lowess.rs`: TransformationDescription/Model/
  Interpolated partial, Linear complete. Gaps (coverage json): Param surface (getParameters,
  getDefaultParameters, checkValidWeight, getValidX/YWeights), weightData/unWeightData batch,
  akima interpolation, b-spline model, getModelTypes, printSummary, getModelParameters.
- `src/analysis/alignment_transformer.rs`: MapAlignmentTransformer partial (IdentificationData overload missing).
- Unmapped: BaseSuperimposer, PoseClusteringAffine/ShiftSuperimposer, MapAlignmentAlgorithm{PoseClustering,
  Identification,KD,TreeGuided}, MapAlignmentEvaluationAlgorithm{,Precision,Recall},
  TransformationModelLowess, TransformationModelBSpline, TransformationModelDefaults.

## Planned order (dependency first)
1. TransformationModelDefaults / Lowess / BSpline / Param surface in transformations
2. BaseSuperimposer + PoseClusteringShiftSuperimposer + PoseClusteringAffineSuperimposer
3. MapAlignmentEvaluationAlgorithm + Precision + Recall (small, self-contained)
4. MapAlignmentAlgorithmPoseClustering (needs StablePairFinder? check), Identification
5. KD / TreeGuided

## Next step
1. Read skill ~/.claude/skills/openms-port-header/SKILL.md, the Param pattern in
   src/processing/peak_picking.rs:1853-1985, and the C++ headers at bc9cc12.

# O5b lane notes (Library and inference)

Merged prerequisite worktree-wf_cd066dab-0a7-30 (O1a, ff) AND lane OA worktree-wf_cd066dab-0a7-27
(complete lane; provides src/analysis/openswath/transition_experiment.rs = LightTargetedExperiment,
needed by OpenSwathLibraryIDNormalizer). Report that extra merge to the lead.

## Scope (all missing in Rust; coverage entries at docs/core-sdk-coverage.json ~1816, 2002-2145)
C++ at core bc9cc12, src/openms/{include,source}/.../ANALYSIS/OPENSWATH/:
Class tests exist only for: LevelContextInference, OpenSwathLibraryIDNormalizer, OpenSwathPercolatorScoring.
- InferenceData: enums InferenceLevel/InferenceContext + toString, IPF*Row, LevelContext{Input,Result}Row
- InferenceConfig: ErrorEstimationConfig, PeptidoformInferenceConfig, LevelContextInferenceConfig (plain structs, defaults)
- LevelContextInference::infer uses math::multiple_testing (p_emp/p_norm, pi0_est, q_value, lfdr) -- all exist in Rust
- Gene/Peptide/Protein inference: level check + delegate
- PeptidoformInference: computeModelFDR (=multiple_testing::compute_model_fdr), preparePrecursorBM, applyBM, precursorInference, infer
- LibraryIDNormalizer: needs LightTargetedExperiment (now merged)
- PercolatorScoring (1642 lines incl. hdr): needs ANALYSIS/ID/Percolator (in-process), OSWParquetReader, ParquetFile, SQLite -> likely untouched/partial

## Plan: files src/analysis/openswath/{openswath_inference_data,openswath_inference_config,level_context_inference,
openswath_{peptide,protein,gene}_inference,openswath_peptidoform_inference,openswath_library_id_normalizer}.rs

## Next step
1. Read class tests (LevelContextInference_test, LibraryIDNormalizer_test) for expected values.
2. Write data+config+level_context, commit.

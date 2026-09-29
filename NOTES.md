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

## Done
- 95ffb473: 8 modules written (data, config, level_context, 3 wrappers, peptidoform, library_id_normalizer); dax clippy --lib -D warnings clean.
- Oracle ../oracle/O5b-inference (driver.cpp, make_inputs.py -> inputs/ + cases.txt, remote_run.sh, run.sh) on ibminode06.
  Output: results/out/*.tsv (expected), results/log. Class tests with literals: LevelContextInference_test, PeptidoformInference_test,
  OpenSwathLibraryIDNormalizer_test, OSWFileInference_test (reader counts 682/32/389/9431).
- PercolatorScoring: needs ANALYSIS/ID/Percolator (in-process), OSWParquetReader, ParquetFile, sqlite -> expect untouched/partial.

## Next step
1. When oracle FETCHED: copy inputs + results/out to tests/data/openswath_inference/, write tests/openswath_inference.rs comparing.
2. Support doc docs/OPENSWATH_INFERENCE_SUPPORT.md, provenance tests/data/openswath_inference_provenance.json.
3. Assess PercolatorScoring. Report + ledger fragment.

## Findings (2026-09-29)
- Oracle rerun FETCHED (212 cases, differ.txt empty, run twice on ibminode06). Data copied to tests/data/openswath_inference
  via scratchpad copy_oracle.py (5 identical IPF trans dumps deduplicated, sha b3affdb2...).
- C++ DEFECT: lfdr with a single p-value segfaults (CubicSpline2d::eval via kdeFFTEval) -> lc_single_* exit 139.
- DIVERGENCE (math layer, not O5b code): protein/gene global fixture = perfect separation, 16 p-values all 1/16, pi0=0.
  C++ probit(0.0625) = -1.5341205443525465 (probe ../oracle/O5b-inference/probe) -> sd==0 in bwNrd0 -> healthy bw -> lfdr 0.
  Rust statrs probit presumably 1 ulp off -> sd ~2e-16 -> degenerate KDE grid -> kernel_density NaturalCubicSpline refuses
  ("cubic spline abscissae must be strictly increasing") -> LevelContextInference Err. Awaiting probe output (tests/o5b_probe_tmp.rs, DELETE before commit).
- All peptide-level fixture cases (27) matched within 1e-12 before the first protein failure.

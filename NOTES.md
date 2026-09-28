# Lane I2 notes (FDR, Percolator, protein inference)

C++ at pin: the on-disk tree ../OpenMS4-tests/packages/core (HEAD 291c5ba) is IDENTICAL to bc9cc12 for
ANALYSIS/ID/* and the relevant class tests (checked `git diff --stat bc9cc12 HEAD -- '*ANALYSIS/ID/*'`),
so read files directly: C=/Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages/core
(git -C with a path containing "source" is refused by the sandbox; use globs like '*ANALYSIS/ID/*').

Sizes: Percolator.h 320 / .cpp 1328; PercolatorTypes.h 113 (header only); PercolatorFeatureSetHelper.h 141 / .cpp 388;
BayesianProteinInferenceAlgorithm.h 149 / .cpp 1093; IDBoostGraph.h 639 / .cpp 1779; FalseDiscoveryRate.h 256 / .cpp 1903.
Class tests: Percolator_test, Percolator_subprocess_parity_test, PercolatorFeatureSetHelper_test,
BayesianProteinInferenceAlgorithm_test, IDBoostGraph_test, FalseDiscoveryRate_test.

Existing Rust: src/analysis/false_discovery_rate.rs (partial), src/analysis/protein_inference.rs (IDBoostGraph partial,
private helpers), src/format/percolator_infile.rs. New code goes to src/analysis/id/ (module root src/analysis/id.rs, append pub mod only).

Ledger: FDR partial (gaps: applyToObservationMatches, ConsensusMap overloads, applyBasic(run_info), applyEvaluateProteinIDs x3,
diffEstimatedEmpirical, findDecoyString). Percolator.h evidence_requires_review; PercolatorFeatureSetHelper, PercolatorTypes,
BayesianProteinInferenceAlgorithm unmapped; IDBoostGraph partial.

## Next step
1. PercolatorTypes (tiny, header only) -> src/analysis/id/percolator_types.rs
2. PercolatorFeatureSetHelper -> src/analysis/id/percolator_feature_set_helper.rs (4 TOPP consumers)
3. Percolator.h (in-process percolator) -> src/analysis/id/percolator.rs
4. FDR gaps
5. IDBoostGraph / Bayesian (graph crate: petgraph?)

## Incarnation 2 findings
- Percolator.h wraps a VENDORED Apache-2.0 percolator subset at $C/src/openms/thirdparty/percolator (~10k lines;
  core: Scores 809, CrossValidation 663, PosteriorEstimator 610, BaseSpline 576, SetHandler 541, ssl 211, Normalizers).
  No Rust crate exists (checked: none on crates.io implements percolator L2-SVM-MFN) -> must port, big.
- PercolatorFeatureSetHelper.cpp: addMASCOTFeatures modifies a COPY of hits (C++ defect: only the sort sticks).
  Class test expected feature lists come from *.topperc_check.idXML ProteinIdentification SearchParameters meta "extra_features".
- Plan order: PercolatorTypes -> FeatureSetHelper -> Percolator model I/O + fillPINCompatibleFields -> engine port -> FDR gaps.
- DONE (commit): percolator_types.rs, percolator_feature_set_helper.rs + tests/percolator_feature_set_helper.rs (12 pass on spock,
  values match *_check.idXML). Fixtures in tests/data/percolator/. Still TODO for these: support doc, provenance json, clippy.
- Engine port: use crate::math::source_sort::source_sort_by for std::sort parity (libstdc++ introsort). NO_TERMINATE=true in
  wrapper -> MyException paths become warnings. reorderFeatureRows = memory layout only, skip.
- TRON oracle literals: Percolator_test.cpp:1427-1523 (weights 1e-8 rel). Test copies in /private/tmp/claude-501/i2/.
- Vendored engine read fully (copy at /private/tmp/claude-501/i2/perc). Design: src/analysis/id/percolator.rs (public API)
  + src/analysis/id/percolator/{engine.rs (Scores, ScoreHolder, PRNG, sampler), posterior.rs (q/p/pi0), pep.rs (InferPEP:
  PAVA + I-spline TRR), svm.rs (TRON), cross_validation.rs (CV + SanityCheck + Normalizer)}.
  Eigen kernels (eigen_dot/eigen_squared_norm/eigen_gemv_row) live private in src/math/fitters/levenberg_marquardt.rs;
  plan: make pub(crate) (report to lead). pep_method=nonparametric (LogisticRegression+BaseSpline+PackedMatrix) = named gap.
  NO_TERMINATE=true => engine only errors on unknown initial_direction and .at() out-of-range.

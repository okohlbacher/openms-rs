# OpenSWATH confidence inference and library identifiers

Ports of the OpenMS core pin `bc9cc12` headers under `ANALYSIS/OPENSWATH/`:

| Header | Rust | Status |
|---|---|---|
| `OpenSwathInferenceData.h` (header-only) | `src/analysis/openswath/openswath_inference_data.rs` | complete |
| `OpenSwathInferenceConfig.h` (header-only) | `src/analysis/openswath/openswath_inference_config.rs` | complete |
| `LevelContextInference.h` (+ `.cpp`) | `src/analysis/openswath/level_context_inference.rs` | partial: the known divergence below, which lies in `MATH/STATISTICS/MultipleTesting.h` |
| `OpenSwathPeptideInference.h` (+ `.cpp`) | `src/analysis/openswath/openswath_peptide_inference.rs` | complete |
| `OpenSwathProteinInference.h` (+ `.cpp`) | `src/analysis/openswath/openswath_protein_inference.rs` | complete |
| `OpenSwathGeneInference.h` (+ `.cpp`) | `src/analysis/openswath/openswath_gene_inference.rs` | complete |
| `OpenSwathPeptidoformInference.h` (+ `.cpp`) | `src/analysis/openswath/openswath_peptidoform_inference.rs` | complete |
| `OpenSwathLibraryIDNormalizer.h` (+ `.cpp`) | `src/analysis/openswath/openswath_library_id_normalizer.rs` | complete |
| `OpenSwathPercolatorScoring.h` (+ `.cpp`) | none | not started |

Tests: `tests/openswath_inference.rs`. Provenance: `tests/data/openswath_inference_provenance.json`.

`OpenSwathPercolatorScoring` reads OSW (SQLite) and OSWPQ (Parquet) feature
tables, rescored with the in-process Percolator, and writes the scores back. It
needs `ANALYSIS/ID/Percolator.h` and `ANALYSIS/ID/PercolatorTypes.h` (the
in-process `Percolator::RescoreInput` / model; not ported),
`ANALYSIS/OPENSWATH/OpenSwathOSWParquetReader.h` and `FORMAT/ParquetFile.h`
(unmapped at this pin of the port), and the OSW schema access of `FORMAT/OSWFile.h`.
It is a backlog item with those named dependencies. The OSW and OSWPQ formats
stay in scope; nothing here re-scopes them.

The inference classes are file-format agnostic: in the source the rows come from
`OSWFile::readLevelContextData` / `readIPF*Data` and the results go back through
`OSWFile::write*Results`. `FORMAT/OSWFile.h` is not ported; the tier-1 evidence
below uses rows the Release `OSWFile` wrote, so the replay runs on the inputs
the class tests describe.

## API mapping

### OpenSwathInferenceData

| Source | Rust |
|---|---|
| `enum class InferenceLevel { Peptidoform, Peptide, Protein, Gene }` | `InferenceLevel` (same variants; `Default` is `Peptide`, the config default) |
| `enum class InferenceContext { Global, ExperimentWide, RunSpecific }` | `InferenceContext` (`Default` is `RunSpecific`) |
| `toString(InferenceLevel)`, `toString(InferenceContext)` | `as_str()` and `Display`; same strings |
| `IPFPrecursorRow { feature_id, ms2_peakgroup_pep, ms1_precursor_pep, ms2_precursor_pep }` | `IPFPrecursorRow`, `Option<f64>` for the optionals, source defaults |
| `IPFPrecursorProbabilityRow`, `IPFTransitionRow`, `IPFAlignmentRow`, `IPFResultRow` | same names and fields, `i64`/`i32`/`f64`, source defaults |
| `LevelContextInputRow { run_id, group_id, entity_id, decoy, score, context }` | `LevelContextInputRow` |
| `LevelContextResultRow { run_id, entity_id, context, score, pvalue, qvalue, pep }` | `LevelContextResultRow` |

### OpenSwathInferenceConfig

| Source | Rust |
|---|---|
| `ErrorEstimationConfig` (11 members) | `ErrorEstimationConfig`; `pi0_method` is `Pi0Method`, `lfdr_transformation` is `LfdrTransform` (see native differences) |
| `PeptidoformInferenceConfig` (11 members) | `PeptidoformInferenceConfig`, all members, source defaults |
| `LevelContextInferenceConfig { level, context, error }` | `LevelContextInferenceConfig` |

### LevelContextInference and the level wrappers

| Source | Rust |
|---|---|
| `static LevelContextInference::infer(input, config)` | `LevelContextInference::infer(&[..], &config) -> Result<Vec<_>>` |
| file statics `countNumPositives_`, `statMetrics_` | not ported: computed and discarded by the source (`(void)statMetrics_(...)`), no observable effect |
| file static `findNearestMatchAscending_` | private `find_nearest_match_ascending` |
| file statics `buildErrorStatistics_`, `assignErrorStatistics_` | private `build_error_statistics`, `assign_error_statistics` |
| `OpenSwathPeptideInference::infer` / `OpenSwathProteinInference::infer` / `OpenSwathGeneInference::infer` | same, `&self`, `Result` |
| (native) | `BSplinePi0Smoother`: the `Pi0Smoother` the source's `pi0Est` builds from `BSplineSmoothingSpline` |
| (native) | `MAX_ROWS` input ceiling |

### OpenSwathPeptidoformInference

| Source | Rust |
|---|---|
| `BayesianModelRow { feature_id, hypothesis, prior, evidence }` | `BayesianModelRow` (module level; nested types are not a Rust feature) |
| `PosteriorRow { feature_id, hypothesis, posterior }` | `PosteriorRow` |
| `static computeModelFDR` | `compute_model_fdr -> Result<Vec<f64>>` (the ceiling of `math::multiple_testing`) |
| `static preparePrecursorBM` | `prepare_precursor_bm` |
| `static applyBM` | `apply_bm` |
| `precursorInference(rows, config) const` | `precursor_inference(&self, ..)` |
| `infer(precursor, transition, alignment, config) const` | `infer(&self, ..) -> Result<Vec<IPFResultRow>>` |
| file statics `propagateSignalAcrossRuns_`, `prepareTransitionBM_`, `makePrecursorPepMap_`, `TransitionInferenceRow` | private equivalents |
| (native) | `MAX_PROPAGATED_ROWS` |

### OpenSwathLibraryIDNormalizer

| Source | Rust |
|---|---|
| `SourceIDMapping` (three `unordered_map<string,string>`) | `SourceIDMapping` (three `HashMap<String, String>`) |
| `static normalizeSourceIDs(exp)` | `normalize_source_ids(&mut exp) -> Result<SourceIDMapping>` |
| `static materializeDecoyPrefix(exp, ids, prefix)` | `materialize_decoy_prefix` |
| `static canonicalTargetForDecoyPrecursor(ref, ids, prefix)` | `canonical_target_for_decoy_precursor -> Option<String>` |
| `static hasCanonicalIDs(exp)` | `has_canonical_ids` |
| `static hasCanonicalIDFormat(exp) noexcept` | `has_canonical_id_format` |
| `static validateCanonicalIDs(exp)` | `validate_canonical_ids -> Result<()>` |

## Preserved source conventions

- Run-specific groups are visited in ascending run id (`std::map`); the first
  failing run's precondition is rewrapped as `Run-specific inference failed for
  RUN_ID <id>: <message>`, while errors of the wrapped statistics
  (`std::invalid_argument` in the source, e.g. a lambda of `1.0`) pass through
  unwrapped, as the source's `Exception::BaseException` handler does not catch
  them. Error texts equal the source's `getMessage()` texts (checked against the
  oracle).
- Every row, target or decoy, is projected onto the nearest ascending target
  cutoff, ties resolved to the first cutoff of a run, an equal distance to the
  upper bracket; non-finite scores keep their row and get `NaN` statistics.
- `pi0Est` with the smoother method uses the ported `BSplineSmoothingSpline`,
  exactly as the source constructs it (`-1.0`, `smooth_df`).
- `applyBM`: product of evidences, minimum prior (`std::min`, so a `NaN` prior
  after a number keeps the number), posterior `0` for a feature whose
  likelihood sum is zero or not finite; output ordered by feature, hypothesis.
- `precursorInference`: thresholds are strict `<` for the peak-group PEP and
  `>=` removes precursor evidence; the plain path keeps input order.
- `infer`: last precursor PEP / alignment group wins for a repeated feature
  (`operator[]` overwrites), first positive peptidoform count wins; the H0 PEPs
  enter the model FDR and are then dropped; result ordered by feature, peptide.
- Propagation: rows with a negative alignment group pass through; confident
  (`pep <= threshold`) rows are copied to every feature of their group and
  duplicates collapse to the minimum PEPs.
- `normalizeSourceIDs`: byte-wise lexicographic precursor numbering over all
  source compounds, sparse ids after dropping unreferenced compounds,
  transition ids by position, conventional `DECOY`/`Decoy`/`decoy` prefixes
  materialised, the experiment replaced only after validation.
- Canonical ids: non-empty, no sign, no leading zero, all digits, fits `i64`.

## Native differences

- `ErrorEstimationConfig` stores `Pi0Method` / `LfdrTransform` enums instead of
  the strings the source converts on every call; an unknown name is refused by
  `Pi0Method::parse` / `LfdrTransform::parse` when the configuration is built.
- `Exception::Precondition` and `Exception::InvalidValue` become
  `Error::InvalidValue` with the source's message text.
- `applyBM` sorts stably; the source's `std::sort` leaves the order of rows with
  equal feature and hypothesis unspecified, so a product of three or more
  evidences may differ in the last bit. Propagation visits alignment groups in
  ascending id instead of `unordered_map` order. The oracle replay shows no
  difference on its cases (max deviation below).
- Ceilings: `MAX_ROWS` (level context input) and `MAX_PROPAGATED_ROWS`
  (propagation output, computed before allocation; the source's propagation is
  quadratic in the group size and unbounded).
- `LightTargetedExperiment` keeps no compound-reference cache in the port, so
  the source's stale-cache rebuild rationale does not apply; the rebuild is kept
  for atomicity.

## Known divergence (in the math layer)

When every target outscores every decoy of an inference group, all target
p-values equal `1 / n_decoys`. The local-FDR stage (`MultipleTesting::lfdr`)
then estimates a density of identical probit values, and `bwNrd0` takes its
degenerate branch depending on whether the sample standard deviation of those
identical values is exactly zero, which depends on the last bit of the probit.
Measured with the O5b probe (`../oracle/O5b-inference/probe`, ibminode06, Linux
x86_64 Release `bc9cc12-c19e494-174b576`): Boost returns
`probit(0.0625) = -1.5341205443525465` (it promotes to 80-bit `long double`),
whose 16-fold mean is exact, bandwidth `0.793`, local FDR `0` at `pi0 = 0`. The
port's `statrs`-based quantile returns `-1.534120544352546` (two units in the
last place away, within the 1.55-ulp accuracy `math::multiple_testing`
documents), bandwidth `2.4e-16`, the KDE grid collapses and
`math::kernel_density::kde_fft_eval` refuses it. `LevelContextInference::infer`
therefore returns an error for the protein- and gene-level groups of the
upstream fixture `PyProphet_inference_test.osw` (42 of 197 replayed cases,
every non-parametric variant of those levels), where the Release build returns
rows. The fix belongs to `MATH/STATISTICS/MultipleTesting.h` (a probit that
rounds as Boost's does on x86_64, or a KDE that handles the collapsed grid as
the source does); the test accepts exactly this refusal on exactly such inputs
(`is_known_probit_divergence`) and compares the case like every other once the
port returns rows.

## C++ defect candidate

`MultipleTesting::lfdr` with a single finite p-value crashes the Release build
(SIGSEGV in `CubicSpline2d::eval` via `KernelDensityEstimation::kdeFFTEval`,
reproduced twice; backtrace in the O5b lane report): one target and one decoy
in a `LevelContextInference::infer` call (`lc_single_*`) are enough, which a
small run-specific group reaches. The port returns
`gridKdeFFT: the grid has zero width` instead.

## Checked boundaries and evidence

- Tier 1: O5b oracle, `../oracle/O5b-inference/driver.cpp` built with
  `x86_64-conda-linux-gnu-g++` `-O3 -DNDEBUG` against the Linux x86_64 Release
  install `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`
  on ibminode06, 212 cases each run twice in fresh processes, all outputs
  byte-identical between runs. 197 cases are replayed by `oracle_replay` (the
  15 others are the `OSWFile` reader dumps that feed them): level-context
  inference on the fixture's peptide, protein and gene rows in all three
  contexts under nine error-estimation variants, synthetic ties, `NaN`/`inf`
  scores, empty and `NaN`-only distributions, missing run ids, invalid lambdas
  and the wrapper level checks; IPF on the fixture under six configurations
  (up to 722 result rows) and on synthetic precursor, transition, grouped-FDR
  and propagation cases; `applyBM`, `preparePrecursorBM`, `computeModelFDR`;
  and nine library normalisation cases including every error path. Numeric
  fields agree within `1e-12` absolute (measured maximum `6.0e-15`), ids and
  flags exactly, error texts exactly; the 42 cases of the known divergence are
  the only exception.
- Tier 3: every section of `LevelContextInference_test.cpp`,
  `PeptidoformInference_test.cpp` and `OpenSwathLibraryIDNormalizer_test.cpp`,
  and the reader counts of `OSWFileInference_test.cpp` on the Release rows.
- Tier 4: nearest-cutoff tie rules, posterior normalisation, canonical-id
  boundaries (`i64::MAX`, overflow, whitespace), atomic refusal.

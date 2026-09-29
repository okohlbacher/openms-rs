# O4b NOTES (lane: TargetedSpectraExtractor, CalibrationWorkflow, MRMRTNormalizer, RANSAC*)

Prereqs merged: 0a7-27, 0a7-30 (commit "[MERGE] O4b prerequisites").

## Next step
1. Read skill ~/.claude/skills/openms-port-header/SKILL.md, C++ headers at bc9cc12.
2. Port RANSAC family into src/math/ml/ (sink; check deps).
3. MRMRTNormalizer -> src/analysis/openswath/mrm_rt_normalizer.rs
4. TargetedSpectraExtractor, CalibrationWorkflow.

## Findings (incarnation 1)
- C++ sources readable as plain files: /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.reference/openms4-core-bc9cc12 (use literal paths; the guard refuses `$VAR` paths and `for ... git` loops).
- RANSAC needs Math::RandomShuffler (MathFunctions.h:477): boost mt19937_64 + boost uniform_int bucket mapping + descending Fisher-Yates.
  Already implemented privately in src/chemistry/decoy_random.rs (rand_mt::Mt64). Plan: src/math/random_shuffler.rs public, decoy_random delegates to its generic helpers.
  RANSAC(uint64_t seed) ctor passes seed to RandomShuffler(int) => truncation to i32 then sign-extend to u64; setSeed does not truncate.
- RANSACModelLinear needs ML/REGRESSION/LinearRegression (unmapped; GTE ApprHeightLine2 fit: mean = sum*(1/n), slope=cov01/cov00, intercept=-slope*mx+my;
  goodness uses Math::mean/variance/covariance (Rust: src/math/statistic_functions.rs) + boost students_t quantile (statrs 0.18 already a dep)).
  RANSACModelQuadratic needs QuadraticRegression (unmapped). Nobody else is assigned ML/REGRESSION -> port into src/math/ml/.
- MRMRTNormalizer::removeOutliersRANSAC default-seeds RANSAC with time(nullptr) (nondeterministic). computeBinnedCoverage: negative bin => UB (defect candidate).
- CalibrationWorkflow needs MRMFeatureFinderScoring, OpenSwathWorkflow, SwathMapMassCorrection, TransitionPQP/TSV (all unmapped) -> probably untouched/partial.
- TSE needs AccurateMassSearchEngine (unmapped), PeakPickerHiRes/Gauss/SGolay (partial), Deisotoper, S/N median, BinnedSpectralContrastAngle (complete).

## Done (incarnation 1)
- src/math/random_shuffler.rs (API identical to lane B8b's concept::math_functions::random_shuffler on branch
  worktree-wf_cd066dab-0a7-39 -> report: lead collapses them; math must own it (sink)). decoy_random.rs UNTOUCHED (B8b moves it).
- src/math/ml/{linear_regression,linear_regression_without_intercept,quadratic_regression,ransac}.rs
- src/analysis/openswath/mrm_rt_normalizer.rs; tests/ransac.rs (13), tests/mrm_rt_normalizer.rs (9) pass on spock.
- Oracle ../oracle/O4b-ransac (gen_data.py -> datasets.tsv, cases.tsv; driver.cpp interposes time()); run via `bash run.sh` (background).

## Next step
1. Fetch oracle results; write tests/ransac_oracle.rs comparing bits (t-quantile: tolerance), fixtures in tests/data/ransac_oracle_*.tsv.
2. docs/RANSAC_SUPPORT.md, tests/data/ransac_provenance.json, THIRD_PARTY entry? (statrs already a dep; add row for students_t).
3. TSE partial (deps: GaussFilter, SavitzkyGolay, PeakPickerHiRes, S/N median, BinnedSpectrum, Deisotoper, TargetedExperiment from O1a).
4. CalibrationWorkflow: likely untouched (needs MRMFeatureFinderScoring, OpenSwathWorkflow, SwathMapMassCorrection...).
5. Report + ledger fragment; delete NOTES.md.

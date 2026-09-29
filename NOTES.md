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

## Done later (incarnation 1)
- RANSAC group complete: tests/ransac_oracle.rs (tier 1, 128 cases), docs/RANSAC_SUPPORT.md, tests/data/ransac_provenance.json.
  All gates on spock green (clippy all-features lib+tests, +1.85 no-default tests, doc -D warnings, doctest).
- TSE: src/analysis/openswath/targeted_spectra_extractor.rs compiles clippy-clean (untested). Gaps: searchSpectrum (AMSE), storeSpectraMSP (MSP writer).

## Next step
1. tests/targeted_spectra_extractor.rs: class-test sections (pickSpectrum inline data; params; mzML-based ones need feature mzml + a test-side
   TargetedExperiment built from TargetedSpectraExtractor_13CFlux_TraML.csv [transition_group_id=peptide ref, RetentionTime minutes, PrecursorMz]).
   Test data at .reference/openms4-core-bc9cc12/src/tests/class_tests/openms/data/TargetedSpectraExtractor_*.
2. docs/TARGETED_SPECTRA_EXTRACTOR_SUPPORT.md + provenance; CalibrationWorkflow assessment (likely untouched).
3. Report + ledger fragment; delete NOTES.md.

## Incarnation 3 (restart after lead stopped wave)
- Merged 0a7-50 (own earlier work incl. TSE class tests 12 pass, executed differential 2682 records bit-identical, support doc). 27/30 already in.
- Ledger fragment wave4/O4b.json already exists (05:43) — review, amend.
- Lane report O4b.md does NOT exist yet.
## Next step (incarnation 3)
1. Verify TSE gaps (searchSpectrum/AMSE, storeSpectraMSP) still open; CalibrationWorkflow assessment.
2. Re-run gates on spock for the merged tree (clippy, tests, +1.85, doc) and record pass counts per slice.
3. Write report, check ledger, delete NOTES.md, final commit.

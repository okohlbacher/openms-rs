# O4a lane notes (Extraction and SWATH windows)

Merged prerequisite worktree-wf_cd066dab-0a7-27 (fast-forward, OA lane: openswath data access, SwathMap, spectrum access, helpers).

Headers (all under src/openms/include/OpenMS/ANALYSIS/OPENSWATH/ at core bc9cc12, tests in
src/tests/class_tests/openms/source/<Name>_test.cpp):
- ChromatogramExtractorAlgorithm.h  (base: others depend on it)
- ChromatogramExtractor.h
- SwathWindowLoader.h
- PeakMapExtractor.h
- SwathQC.h
- SwathMapMassCorrection.h
TOPP: OpenSwathChromatogramExtractor? (check topp 174b576 / cli c19e494)

Rust target: src/analysis/openswath/<snake>.rs, pub mod appended to src/analysis/openswath.rs.
Never run cargo locally; gate: OPENMS_GATE_HOST=dax OPENMS_GATE_JOBS=64 ~/.local/bin/openms-kim-gate.sh O4a ... (background).

## Findings (read before coding)
- Rust deps present: SpectrumAccess trait (spectrum_access.rs), light Spectrum/Chromatogram (data_structures.rs),
  SwathMap (swath_map.rs), LightTargetedExperiment (transition_experiment.rs), openswath_helper::compute_precursor_id /
  compute_transition_group_id, data_access_helper::convert_to_openms_chromatogram, processing::deisotoping::Deisotoper,
  processing::peak_picking::PeakPickerHiRes, math::statistic_functions::adaptive_quantile, kernel::mrm::MRMTransitionGroup,
  analysis::transformations::TransformationDescription, constants C13C12_MASSDIFF_U.
- ChromatogramExtractor TOPP consumer: FeatureFinderMetaboIdent only (not ours to port).
- C++ quirks found: tophat left-walk never adds mz_start after >1 step; when mz_it==end the last peak is added twice
  (current + first left step) -> defect candidate, preserve & document; prepare_coordinates MS2 derefs tr2pep[missing]=null (UB)
  -> Rust errors. SwathWindowLoader: trailing blank line keeps previous lower/upper (sentry fails, values untouched) -> duplicate row
  (verify on ibminode06 oracle). annotate non-atomic in C++.
- SwathMapMassCorrection (819 lines .cpp) needs MRMFeatureFinderScoring::MRMTransitionGroupType -> likely partial.

## Progress (committed)
- chromatogram_extractor_algorithm.rs, swath_window_loader.rs, peak_map_extractor.rs, chromatogram_extractor.rs
  + tests/openswath_extraction.rs (21 tests pass on dax: gate `test --locked --test openswath_extraction`).
- Oracle ../oracle/O4a-extraction (driver.cpp, node.sh, run.sh push|pull) started on ibminode06, results/DONE marker.
  Covers tophat, xic file/IM, peakmap, prepare_coordinates/return_chromatogram, SwathWindowLoader edge files,
  SwathQC subsample grid + charge distribution + JSON, SMC estimateWindow/params/findSwathMapsPasef.

## Next step
1. Pull oracle (run.sh pull), compare vs Rust; add tier-1 tests from records.tsv (tests/data/openswath_extraction_oracle.tsv).
2. SwathQC -> swath_qc.rs (Deisotoper + PeakPickerHiRes exist).
3. SwathMapMassCorrection -> partial: params, getters, estimate_window, find_swath_maps_pasef; correctMZ/IM gaps
   (need OpenSwathScoring::fetchSpectrumSwath, DIAHelpers, Linear/QuadraticRegression).
4. Support doc docs/OPENSWATH_EXTRACTION_SUPPORT.md, provenance json, doc coverage, full CI sweep, report + ledger.

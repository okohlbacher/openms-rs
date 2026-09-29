# O3a NOTES (chromatogram peak picking)

Branch: worktree-wf_cd066dab-0a7-48. Prereqs 0a7-27 (OA: OPENSWATHALGO, data access) and 0a7-30 (TargetedExperiment, MRMMapping) merged at 5adedcb (conflict in src/analysis/openswath.rs resolved keeping both).

## State of headers (from docs/core-sdk-coverage.json)
- PeakPickerChromatogram.h: partial, src/processing/chromatogram.rs. Gaps: updateMembers_, DefaultParamHandler surface (11 params), write_sn_log_messages_, copy-assign.
- PeakIntegrator.h: partial, src/analysis/peak_integrator.rs. Gaps: 6 ConstIterator overloads; getDefaultParameters/Param surface.
- PeakPickerMobilogram.h: unmapped.
- MRMTransitionGroupPicker.h: unmapped (1726 lines). Depends on PeakPickerChromatogram, PeakIntegrator, MRMFeature/MRMTransitionGroup kernel types?
- MasstraceCorrelator.h: unmapped; TOPP ClusterMassTraces, ClusterMassTracesByPrecursor.

## Param pattern: src/processing/peak_picking.rs:1853-1985 (defaults/from_param/from_param_with_warnings/to_param)

## Next step
1. Close PeakPickerChromatogram Param surface.
2. PeakIntegrator Param surface + iterator overloads (slice/range entry points).
3. MRMTransitionGroupPicker.
4. PeakPickerMobilogram, MasstraceCorrelator.

## Findings (session 1)
- C++ pin readable at /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.reference/openms4-core-bc9cc12 (sha256 verified vs coverage json). Use Read tool on it: Bash with `$VAR` paths or `git -C` compound commands is refused by the worktree guard.
- MRMTransitionGroupPicker is header-mostly (.h 1579 lines, .cpp 147). Needs: MRMTransitionGroup/MRMFeature (src/kernel/mrm.rs, complete), PeakPickerChromatogram (+Param), PeakIntegrator (+Param), LinearResamplerAlign raster(iter) = free `resample_to_grid`, Scoring (normalizedCrossCorrelation, xcorrArrayGetMaxPeak, computeAndAppendRank, rankedMutualInformation) in src/analysis/openswath/scoring.rs, StatsHelpers mean_and_stddev.
- MRMTransitionGroupPicker params: stop_after_feature -1, stop_after_intensity_ratio 1e-4, min_peak_width 0.001, peak_integration original|smoothed, background_subtraction none|original|exact, recalculate_peaks false, use_precursors false, use_consensus true, recalculate_peaks_max_z 1.0, minimal_quality -10000, resample_boundary 15, compute_peak_quality false, compute_peak_shape_metrics false, compute_total_mi false, boundary_selection_method largest|widest, subsections PeakPickerChromatogram:, PeakIntegrator:.
- C++ defect candidates: (a) computeQuality_ MRMTransitionGroupPicker.h:1322 `max_int = intensity > max_int;` stores a bool (1.0) -> after first in-range peak, any later peak with intensity > 1 replaces the borders (last, not most intense). (b) pickFragmentChromatograms smoothed_chroms[k] indexed by transition k while smoothed_chroms holds only detecting chromatograms (in chromatogram order) -> misaligned when non-detecting transitions exist or transition order != chromatogram order. (c) PeakPickerChromatogram updateMembers_ never re-seeds pp_ (seed S/N fixed at construction).
- Helpers float_entry etc. are pub(super) in src/processing/peak_picking/noise.rs; need widening to pub(crate).

## Progress
- DONE (795edde + next): PeakPickerChromatogram Param surface (defaults/to_param/from_param[_with_warnings]), IDX_* consts, noise_log in result, inactive_smoothing field; default snt write_log_messages=false (source). PeakIntegrator Param surface, 7 source string consts, six *_between index overloads. tests/openswath_peak_picking.rs 10 pass on spock (gate: `--test openswath_peak_picking`).
- Helpers float_entry etc. now pub(crate), re-exported from processing::peak_picking.
- cargo fmt reordered src/analysis/openswath.rs mod lines (merged file was unformatted) -> mention in report.
- NEXT: MRMTransitionGroupPicker in src/analysis/openswath/mrm_transition_group_picker.rs. Needs: local source-exact raster (dist_right/(dl+dr), f32 accumulate), local stable sortByIntensity with arrays, UniqueIdGenerator param for ensureUniqueId, MetaInfo on Feature.base.metadata.
- WIP: src/analysis/openswath/mrm_transition_group_picker.rs written (all public members incl. pick_apex/pick_fragment/pick_precursor/remove_overlapping/find_*; params; meta consts). Tests appended to tests/openswath_peak_picking.rs (class-test transcription). Scratch files in scratchpad/O3a/ (shared scratchpad has other lanes' g*.log!).
- More defect candidates: (d) computeQuality_ outlier index counts only traces with signal but indexes picked_chroms; (e) total_mi divides the already-averaged MI by detecting count again; (f) stop_after_feature collects N+1 features (cnt > N); (g) mean_and_stddev functor reused: coel mean averages shapes+coel; (h) pickPrecursorChromatograms prec_idx = all fragment chroms + k (misaligned when non-detecting fragments/!use_precursors); (i) pickApex integrates intensity-sorted picked chroms with PosBegin/PosEnd binary searches (reproduced via libstdc++ lower/upper_bound emulation).
- DONE (f0e8543): MRM picker bit-exact vs ibminode06 Release oracle (../oracle/O3a-mrm-picker, 32 cases, 3358 records), after switching scoring.rs standardize_data + calculate_cross_correlation to Eigen SSE2 redux order (eigen_sum). Provenance tests/data/openswath_peak_picking_provenance.json; docs MRM_TRANSITION_GROUP_PICKER_SUPPORT.md, CHROMATOGRAM_PICKING_SUPPORT.md, PEAK_INTEGRATION_SUPPORT.md updated. clippy/doc/MSRV clean at ee923e0 (rerun before final).
- NEXT: (1) python gates + CI sweep (log dir outside repo); (2) ledger fragment .planning/ledger-fragments/wave4/O3a.json + report .planning/lane-reports/O3a.md; (3) PeakPickerMobilogram (needs mobilogram smoothing/HiRes via chromatogram conversion; source remove_overlapping_ never initialised!); (4) MasstraceCorrelator.
- DONE through b3be24f + test fix: PeakPickerMobilogram (complete), MasstraceCorrelator (complete; map_xml.rs dataProcessing-without-software leniency; test substitutes rt="4the_second_seed202.0"->0.0 as source asDouble_ does). Ledger fragment written to .planning/ledger-fragments/wave4/O3a.json (MRM partial: pickApex+fit_EMG Unsupported, pools not ported).
- Final CI sweep running: slot O3a-sweep, logs scratchpad/O3a/sweep2. Slice counts: scratchpad/O3a/slice{1,2,3}.log.
- NEXT: when sweep passes -> write report .planning/lane-reports/O3a.md, delete NOTES.md, final commit.

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

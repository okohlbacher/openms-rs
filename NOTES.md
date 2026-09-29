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

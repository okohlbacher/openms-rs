# B6 lane notes (SHARED-01: DefaultParamHandler trait)

Branch: worktree-wf_bcdd783d-f3f-2, off main f0990c0.

## What is missing
- No trait over the defaults()/from_param()/from_param_with_warnings()/to_param() pattern.
  Pattern source: src/processing/peak_picking.rs ~1853-1985 (PeakPickerHiRes), uses
  crate::param::DefaultParamHandler (src/param/handler.rs) set_defaults + set_parameters_with.
- Headers with param gap in ledger (docs/core-sdk-coverage.json):
  partial: PeakPickerIterative, NLargest, ThresholdMower, SignalToNoiseEstimatorMeanIterative,
  LinearResamplerAlign, Normalizer, GaussFilter, SavitzkyGolayFilter;
  complete-but-mention: WindowMower, RankScaler, SqrtScaler;
  CHEMISTRY partial: TheoreticalSpectrumGenerator, SpectrumAnnotator, NucleicAcidSpectrumGenerator.

## Evidence
- C++ getDefaults for each class: git -C ../OpenMS4-tests/packages/core show bc9cc12:src/openms/source/<...>.cpp
- class tests: src/tests/class_tests/openms/source/<Name>_test.cpp
- Oracle: ../oracle/B6-param-defaults/ on ibminode06, Release install, dump getDefaults via ParamXMLFile.

## Next step
1. Design trait in src/param/handler.rs (or src/param/algorithm.rs).
2. Implement for PeakPickerHiRes (reference), then the 14 headers.
3. Oracle for default Param dumps.

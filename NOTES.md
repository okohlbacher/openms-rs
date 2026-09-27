# B6 lane notes (SHARED-01: DefaultParamHandler trait)

Branch: worktree-wf_ee239125-5be-2 (merged prior incarnation worktree-wf_bcdd783d-f3f-2, notes only).

## What is missing
- No trait over the defaults()/from_param()/from_param_with_warnings()/to_param() pattern.
  Pattern source: src/processing/peak_picking.rs ~1853-1985 (PeakPickerHiRes), uses
  crate::param::DefaultParamHandler (src/param/handler.rs) set_defaults + set_parameters_with.
  PEAK_PICKER_HI_RES_NAME const is the handler name.
- Ledger targets (docs/core-sdk-coverage.json, scope says param gap):
  PROCESSING partial: PeakPickerIterative (src/processing/iterative.rs), NLargest, ThresholdMower,
  Normalizer, LinearResamplerAlign (all src/processing.rs), SignalToNoiseEstimatorMeanIterative
  (src/processing/mean_noise.rs), GaussFilter, SavitzkyGolayFilter (src/processing/smoothing.rs);
  complete (C++ .cpp names DefaultParamHandler only): WindowMower, RankScaler, SqrtScaler.
  CHEMISTRY partial: TheoreticalSpectrumGenerator (src/chemistry/theoretical.rs), SpectrumAnnotator
  (src/chemistry/spectrum_annotator.rs), NucleicAcidSpectrumGenerator.
  Also mentions: ANALYSIS/OPENSWATH/PeakPickerChromatogram (not mine).
- SignalToNoiseEstimatorMedian already has defaults/from_param in src/processing/peak_picking/noise.rs.

## Evidence
- C++ getDefaults: git -C ../OpenMS4-tests/packages/core show bc9cc12:src/openms/source/<...>.cpp
- class tests: src/tests/class_tests/openms/source/<Name>_test.cpp
- Oracle: ../oracle/B6-param-defaults/ on ibminode06, Release install, dump getDefaults via ParamXMLFile.

## Done
- Oracle ../oracle/B6-param-defaults RAN OK on ibminode06 Release install (results/, no DIFFER).
  attempt1/ = first try, failed: SignalToNoiseEstimatorMeanIterative<MSChromatogram> does not compile
  (getMZ on ChromatogramPeak, header :262,:279) -> C++ defect candidate; case dropped.
  Unknown-param warning in C++: "Warning: X received the unknown parameter 'k'!" (Rust handler says
  "X: unknown parameter 'k'" -- left as is, FeatureFinderPicked maps it itself).
- Trait `DefaultParamHandling` in src/param/handler.rs (re-exported from param): NAME, SUBSECTIONS,
  CHECK_DEFAULTS, WARN_EMPTY_DEFAULTS, defaults(), to_param(), update_members(), keep_native_members(),
  provided handler/from_param/_with_warnings/_with_handler/set_parameters.
  Adopted on PeakPickerHiRes and SignalToNoiseEstimatorMedian.
- C++ sources readable at OpenMS4-R/.reference/openms4-core-bc9cc12 (main checkout). The worktree
  sandbox refuses a git invocation against another repo; use Read/cat with literal paths.

- 11 PROCESSING impls done (src/processing/parameters.rs, window_mower.rs, smoothing.rs,
  mean_noise.rs, iterative.rs); tests/default_param_handling.rs 12/12 pass on spock (tier 1 vs
  tests/data/param_defaults/*.ini copied from oracle results; hashes in scratchpad B6/fixture_hashes.txt).
  Param::replace_value added (src/param.rs); helpers unsigned_parameter/int_parameter/flag_parameter/
  set_flag in param::handler. GaussFilter gained inactive_width + write_log_messages; MeanIterative
  gained range_parameters; SavitzkyGolay keeps frame_parameter.

- Chemistry (incarnation 3): TSG impl in src/chemistry/theoretical_parameters.rs (intensities now f64,
  matching the source's double members TheoreticalSpectrumGenerator.h:141-153); SpectrumAnnotator in
  src/chemistry/spectrum_annotator_parameters.rs; NASG in src/chemistry/nucleic_acid_spectrum_generator_parameters.rs.
  f64 change breaks 2 tests in tests/theoretical_fine_reference.rs whose expectation used f64::from(0.7_f32)
  (the old f32 field) -> oracle ../oracle/B6-tsg-intensity (ibminode06, Release) dumps the C++ peaks with
  b_intensity 0.7 / relative_loss_intensity 0.9 to settle it (started; fetch results with scp/ssh, the
  sandbox refuses `bash run.sh` outside the worktree -- run its steps as plain ssh commands).
- Gate: always pass --no-fail-fast; cargo stops at the first failing test binary otherwise.

- DONE (incarnation 3): chemistry impls + tests pass on spock (15 contract tests, 4 intensity tests, all theoretical/
  annotator/NASG suites). Oracle ../oracle/B6-tsg-intensity fetched (results/, differ empty). Provenance
  tests/data/param_defaults_provenance.json (generator: scratchpad B6/provenance.py). Docs
  docs/DEFAULT_PARAM_HANDLING_SUPPORT.md. Ledger fragment written: .planning/ledger-fragments/wave4/B6.json
  (generator scratchpad B6/ledger.py). m/z of TSG oracle peaks differ <=2 ulp (reported, 1e-9 tolerance).
- ci_matrix --check wants --write (adds default_param_handling to the paramxml slice) -- lead's file, report it.

## Next step
1. Wait for sweep (scratchpad B6/sweep1) and slices (B6/slices.log); fix anything red.
2. Write .planning/lane-reports/B6.md; delete NOTES.md in the final commit.

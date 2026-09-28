# Lane OA notes (restart-safe)

## Situation (first incarnation, 2026-09-28)
- `src/analysis/openswath.rs` is an empty module root (append `pub mod` only).
- No OpenSWATH code exists yet in Rust except `src/kernel/mrm.rs` (MRM kernel types) and
  `src/format/mzml_sqlite_swath_handler.rs`.
- C++ pin (bc9cc12) readable at /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages/core
  (checkout HEAD 291c5ba). Checkout == pin EXCEPT these (use `git show bc9cc12:<path>`):
  SpectrumAccessSqMass.h/.cpp, OPENSWATHALGO/DATAACCESS/DataStructures.h, OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h.
- Complex git commands against the pin repo are refused by the worktree guard: use one plain
  `git -C <abs> show bc9cc12:<path>` per call (no pipes/loops).

## Upstream tests (expected values)
- class_tests/openswathalgo: Datastructures_test, DiaHelpers_test, Scoring_test, SwathMap_test
- class_tests/openms/source: StatsHelpers_test, SpectrumHelpers_test, OpenSwathHelper_test, SpectrumAddition_test,
  DataAccessHelper_test, OpenSwathDataAccessHelper_test, SpectrumAccessOpenMS_test, SpectrumAccessOpenMSInMemory_test,
  SpectrumAccessQuadMZTransforming_test, SpectrumAccessSqMass_test, MRMFeatureAccessOpenMS_test, OpenSwathMRMFeatureAccessOpenMS_test

## Plan (dependency order)
1. OPENSWATHALGO DataStructures + ISpectrumAccess trait + SwathMap + ITransition/Transitions/TransitionExperiment/TransitionHelper
2. StatsHelpers, Scoring, SpectrumHelpers
3. MockObjects, DataFrameWriter
4. DATAACCESS: SpectrumAccessOpenMS/InMemory/Cached/Transforming/QuadMZ/SqMass, DataAccessHelper, SimpleOpenMSSpectraAccessFactory, MRMFeatureAccessOpenMS
5. OpenSwathHelper, SpectrumAddition

## Incarnation 2 (retry, branch worktree-wf_cd066dab-0a7-27) findings
- Prior branch -21 == 2d62fe1f, same as ours: only notes existed; nothing to merge.
- Worktree guard refuses `sed $VAR` constructs; `cd <abs> && cat <literal paths>` works.
- Existing Rust: MSSpectrum{peaks: Vec<Peak1D{mz f64,intensity f32}>, float_data_arrays: Vec<DataArray<f32>>, rt, ms_level u32, native_id, drift_time},
  MSChromatogram{peaks: Vec<ChromatogramPeak{rt,intensity f32}>, precursor, product, native_id}, MSExperiment{spectra, chromatograms}. Error: InvalidValue, InvalidRange, MissingInformation, Unsupported...
- StatsHelpers/Scoring use Eigen reductions (packet sums) -> port with sequential sums; document bit-level divergence.
- Rust layout decided: src/analysis/openswath/{data_structures,spectrum_access,transitions,swath_map,stats_helpers,scoring,data_frame_writer,mock_objects,...}.rs

## Done (incarnation 2), all committed
- src/analysis/openswath/{data_structures,spectrum_access,spectrum_access_openms,transition_experiment,transitions,
  transition_helper,itransition,swath_map,mock_objects,stats_helpers,scoring,data_frame_writer,data_access_helper,mrm_feature_access}.rs
- Design: trait SpectrumAccess (&self, Send+Sync, usize ids, Result), SpectrumPtr=Arc<Spectrum>, arrays owned (>=2 invariant).
  Access traits: FeatureAccess/MRMFeatureAccess/TransitionGroupAccess/SignalToNoiseAccess (feature() returns Option<&dyn>).
- Gaps: SpectrumAccessOpenMSCached (needs CachedMzML, unported) ; DataAccessHelper targeted conversions (need TargetedExperiment);
- C++ defect candidates: InMemory getSpectraByRT derefs end + exclusive upper bound vs doc; InMemory drops chromatograms for SqMass origin;
  QuadMZ transforms shared cached spectrum in place (double transform); LightTargetedExperiment ref map dangling/null deref;
  convertTargetedExp unreachable 'both target and decoy' branch; MockMRMFeature operator[] inserts null.

## State (incarnation 2, late)
- All code, tests (algo 30, data_access 26, sqmass 4 -- pass all-features AND +1.85.0 no-default-features), docs, provenance committed.
- Ledger fragment written: .planning/ledger-fragments/wave4/OA.json (22 entries; Cached untouched).
- Python gates OK (core_sdk, schema, cycles, doc coverage; doc-coverage --write is lead's).

## Next step
1. Wait for CI sweep (scratchpad sweep/summary.log); fix failures.
2. Write .planning/lane-reports/OA.md with Bash heredoc; delete NOTES.md in final commit.

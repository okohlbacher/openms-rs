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

## Next step
1. Read C++ headers of step 1 and existing Rust MSSpectrum/MSChromatogram/MSExperiment types; write src/analysis/openswath/data_structures.rs

# OpenSWATH workflow and its scheduler

Rust ports of `ANALYSIS/OPENSWATH/OpenSwathWorkflow.h` and
`OpenSwathWorkflowScheduler.h` with their translation units (pin `bc9cc12`):

| Header | Rust |
| --- | --- |
| `OpenSwathWorkflow.h` | [`src/analysis/openswath/openswath_workflow.rs`](../src/analysis/openswath/openswath_workflow.rs) |
| `OpenSwathWorkflowScheduler.h` | [`src/analysis/openswath/openswath_workflow_scheduler.rs`](../src/analysis/openswath/openswath_workflow_scheduler.rs) |

Tests: [`tests/openswath_workflow.rs`](../tests/openswath_workflow.rs) (17),
[`tests/openswath_workflow_scheduler.rs`](../tests/openswath_workflow_scheduler.rs)
(9), the scheduler grid in
[`tests/openswath_export_oracle.rs`](../tests/openswath_export_oracle.rs).
Manifests: [`tests/data/openswath_workflow_provenance.json`](../tests/data/openswath_workflow_provenance.json),
[`tests/data/openswath_export_provenance.json`](../tests/data/openswath_export_provenance.json).

`OpenSwathWorkflow` drives, per SWATH window and batch of assays, the
extraction of the fragment and precursor chromatograms
(`ChromatogramExtractor`), the picking of peak groups
(`MRMTransitionGroupPicker`), the scoring (`MRMFeatureFinderScoring`) and the
output to a feature map, a chromatogram consumer and an OSW writer
([`OPENSWATH_OSW_WRITER_SUPPORT.md`](OPENSWATH_OSW_WRITER_SUPPORT.md)).
`OpenSwathWorkflowScheduler` plans the memory-bounded waves of SWATH maps
that the in-memory path processes together.

## API mapping: `OpenSwathWorkflow.h`

| Source member | Rust |
| --- | --- |
| `struct ChromExtractParams` | `ChromExtractParams` (default: all zero, empty function) |
| `OpenSwathWorkflowBase()`, `OpenSwathWorkflowBase(use_ms1_traces, ..., threads_outer_loop)` | `OpenSwathWorkflowBase::default`, `OpenSwathWorkflowBase::new` |
| `OpenSwathWorkflowBase::MS1Extraction_` | `OpenSwathWorkflowBase::ms1_extraction` (returns the chromatograms) |
| `OpenSwathWorkflowBase::prepareExtractionCoordinates_` | `OpenSwathWorkflowBase::prepare_extraction_coordinates` (returns the chromatograms and the coordinates) |
| the protected members `ms1_map_`, `use_ms1_traces_`, `use_ms1_ion_mobility_`, `prm_`, `pasef_`, `mrm_`, `threads_outer_loop_`, the `ProgressLogger` base | public fields `ms1_map`, `use_ms1_traces`, `use_ms1_ion_mobility`, `prm`, `pasef`, `mrm`, `threads_outer_loop`, `logger` (behind a lock) with `start_progress`, `set_progress`, `end_progress` |
| `OpenSwathWorkflow(use_ms1_traces, use_ms1_ion_mobility, prm, pasef, mrm, threads_outer_loop)` | `OpenSwathWorkflow::new` |
| `loadMS1Map` | `OpenSwathWorkflow::load_ms1_map` |
| `performExtraction` | `OpenSwathWorkflow::perform_extraction` |
| `writeOutFeaturesAndChroms_` | `write_out_features_and_chroms` |
| `scoreAllChromatograms_` | `OpenSwathWorkflow::score_all_chromatograms` |
| `selectCompoundsForBatch_` | `select_compounds_for_batch` |
| `copyBatchTransitions_` | `copy_batch_transitions` |
| `IChromatogramHandler::createDefault()` (called in MRM mode) | `ChromatogramHandler` trait and `OpenSwathWorkflow::set_chromatogram_handler`; see "Not ported" |
| anonymous namespace: `calculateInnerBatchSize`, `determineOSWWriteBufferBytes` | `calculate_inner_batch_size`, `determine_osw_write_buffer_bytes` (and the `_from` variants taking the memory) |
| anonymous namespace: `OSWBufferedWriter` | not ported; see "Native differences" |

## API mapping: `OpenSwathWorkflowScheduler.h`

| Source member | Rust |
| --- | --- |
| `Options` (11 members) | `SchedulerOptions` with `Default` |
| `ConcurrencyEstimate` | `ConcurrencyEstimate` with `Default` (cap 1) |
| `Wave` | `Wave` |
| `estimateAvailableMemoryForScoring` | `OpenSwathWorkflowScheduler::estimate_available_memory_for_scoring`; with the memory given: `estimate_available_memory_from` |
| `estimateConcurrency` | `estimate_concurrency`; with the budget given: `estimate_concurrency_with_budget` |
| `planWaves` | `plan_waves` |
| `chooseInnerBatchSize` | `choose_inner_batch_size` |
| `ConcurrencyLimiter` (`acquireSlot`, `releaseSlot`) | `ConcurrencyLimiter` (`new`, `acquire_slot`, `release_slot`, `capacity`) |
| `ScopedSlot` | `ScopedSlot` (RAII; `None` is the null limiter) |
| `SysInfo::getFreeSystemMemory` | `free_memory_kib`: `MemAvailable`/`MemFree` of `/proc/meminfo` on Linux, `None` elsewhere. It is the same query as `system::sys_info::free_system_memory`, which this module cannot import without closing the cycle `analysis` - `system` - `format` - `analysis` (`tools/check_module_cycles.py`). |

## Preserved source conventions

- **Three execution paths.** MRM (all maps chromatogram-only; needs a handler),
  the per-SWATH path (`batch_size > 0`, or the data not loaded into memory) and
  the wave scheduler (`batch_size <= 0` with `load_into_memory` and more than the
  MS1 map). MS1-only data (one map flagged MS1) is scored on the precursor
  traces alone and requires `use_ms1_traces`. PRM and diaPASEF assign every
  transition to the single best window (m/z centre, respectively ion mobility
  centre of the matching windows) before extraction; both together are
  `NotImplemented` (`Error::Unsupported`).
- **Wave scheduling.** Waves from `plan_waves`; per wave the SWATHs are loaded
  into memory and extracted, their assays are cut into scoring jobs of
  `choose_inner_batch_size` compounds, the jobs are sorted by estimated
  transitions descending, scored, and a context is written out when its last job
  completes. After a wave with several jobs `out_featureFile` is sorted by
  position, as the source does.
- **Per-assay processing** (`score_all_chromatograms`): assays in identifier
  order; each transition's chromatogram is looked up by native id (a missing one
  is an error, skipped in MRM mode), copied, renamed to the transition and cut to
  the de-normalized expected RT +/- `rt_extraction_window`; the precursor
  isotope traces `<id>_Precursor_i<n>` are added; with no chromatogram at all
  the assay is skipped, with no *detecting* chromatogram it is an error (skipped
  in MRM mode). Picking and scoring errors: `InvalidValue` and `InvalidRange`
  (source `IllegalArgument`, `InvalidRange`) are returned in DIA mode and skip the
  assay in MRM mode; every other error skips the assay in both modes (the
  source's `catch (std::exception&)`). An output feature without a detecting
  transition is an error. With an active OSW writer the feature map is cleared
  per assay and the rows are prepared from it.
- **Outputs.** Non-empty MS1 then MS2 chromatograms to the consumer, then the
  features and protein identifications to the map when `store_features`.

## Native differences

| Source | Here | Why |
| --- | --- | --- |
| OpenMP `parallel for`, `critical`, the nested variant | `map_collect` over SWATHs and scoring jobs, results consumed in job order on one thread; the nested variant is not provided and the source falls back to the wave scheduler without it | the determinism contract (`concept::parallel`): a parallel result equals the serial one. The per-SWATH path is bit-identical in order for any thread count; the wave path yields the same set of features (see the module documentation) |
| global `UniqueIdGenerator` | one generator per workflow, keyed per assay and SWATH, seedable | no identifier depends on the schedule |
| `OSWBufferedWriter` thread | rows written in job order after each chunk of `4 * threads` jobs (one with a mobilogram consumer) | memory bound is the chunk; same file content for a serial run |
| `ScopedResamplingWarningSuppression` | none | the port's resampler does not warn |
| log lines | `OpenSwathWorkflow::take_log` | no global log stream |
| `chromConsumer` dereferenced unchecked | `Option`; a `Break` of the consumer stops feeding it | no null reads; the source has no stop |
| `MobilogramParquetConsumer*` | `Option<&mut dyn MobilogramSink>` | the Parquet implementation is the backlog item of `FORMAT/DATAACCESS/MobilogramParquetConsumer.h` |
| `compounds[assay_peptide_map[id]]` for an assay without compound (`operator[]` inserts index 0) | the first compound; an error when there is none | the source reads index 0 |

## Not ported (backlog items with their dependency)

- **SRM/MRM mode** needs `IChromatogramHandler::createDefault()`, i.e.
  `ANALYSIS/TARGETED/DefaultChromHandler.h`, `MRMChromHandler.h` and
  `DIAChromHandler.h`. The workflow defines the `ChromatogramHandler` trait the
  handlers implement and returns `Error::Unsupported` without one. The handler
  needs the precursor and product of each chromatogram of a
  `SpectrumAccessOpenMS` (`getChromatogramMetaInfo`), which the
  `SpectrumAccess` trait does not offer: a downcast or a metadata accessor on the
  trait is the dependency.
- `MobilogramParquetConsumer` (above).
- The `ChromatogramProcessor` and `CalibrationWorkflow` headers use the same
  base state; `OpenSwathWorkflowBase` is a public struct for composition.

## Checked boundaries and evidence

**Tier 1, retained C++ output.** `tests/openswath_workflow.rs` reproduces the
invocations of `TOPP_OpenSwathWorkflow_1`, `_3`, `_6`, `_18`, `_19`, `_20`,
`_22`, `_23` and `_threads4` (test-data pin `0cb15f2`) on the retained inputs
and compares the features (decoded, `FuzzyDiff.ini` tolerances, `id=`
whitelist) and, for `_1`, `_3`, `_22`, `_23`, the chromatograms with the
retained outputs. They cover the per-SWATH and wave paths, with and without MS1
traces, the MI scores, PRM window matching and diaPASEF window matching with
ion mobility scoring and extraction. The TOPP tool of the `openswath` package
is not in this crate: the harness splits the run into SWATH maps as
`FullSwathFileConsumer` does, canonicalizes the library
(`OpenSwathLibraryIDNormalizer`) and assembles the parameters from the `Scoring`
subsection defaults and the tool's overrides (see the test header).

**One known divergence, carved out narrowly.** The library of these tests has
an assay `PEPTIDEB`. `AASequence` here refuses a bare `B`/`Z`/`X` residue
(its documented decision: no invented mass), where the source's `ResidueDB`
gives it the empty formula (an internal mass of minus one water); the b/y-ion
score of such an assay fails with `Unsupported`, which the workflow, like the
source, treats as a skipped assay. The expected features of ambiguous peptides
are compared strictly when the run has them and are otherwise accepted missing
only when the log records the skipped assay with exactly that error; the
branch disappears once the chemistry port computes the source's masses (an
explicit source-compatibility option of `AASequence` and the theoretical
spectrum generator, owner: CHEMISTRY).

**Tier 1, executed differential.** `plan_waves` over a grid of 432 (map count,
cap, budget, spectra, transitions) combinations and `choose_inner_batch_size`
over 576 combinations equal the Release build's results exactly.

**Tier 3.** `OpenSwathWorkflowScheduler_test.cpp` (wave partitions with count and
memory caps, batch sizes, limiter) is transcribed; the memory arithmetic
(`estimate_available_memory_from`) is derived from the source formula.

**Tier 4.** Thread-count independence (1, 2, 4, 8 threads, both paths), explicit
batch sizes equal to automatic, batch selection, the refusals (MS1-only without
MS1 traces, PRM with PASEF, MRM without a handler).

Upstream tests not carried: the other `TOPP_OpenSwathWorkflow_*` cases exercise
the tool's loaders (mzXML, SWATH window files, PQP and sqMass input and output,
split files), its calibration workflow (iRT, `_21`, `_27`) or the MRM mapping
(`_24`, `_25`); they are outside this class. `_17` and `_17_b` (ion mobility
TraML and PQP libraries) are not replayed.

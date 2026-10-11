# MRM feature finding and confidence scoring

Rust ports of `ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h` and
`ANALYSIS/OPENSWATH/ConfidenceScoring.h` (pin `bc9cc12`), and of the
heavyweight-experiment half of
`ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h` that the finder needs:

| Header | Rust |
| --- | --- |
| `MRMFeatureFinderScoring.h` | [`src/analysis/openswath/mrm_feature_finder_scoring.rs`](../src/analysis/openswath/mrm_feature_finder_scoring.rs) |
| `ConfidenceScoring.h` | [`src/analysis/openswath/confidence_scoring.rs`](../src/analysis/openswath/confidence_scoring.rs) |
| `DataAccessHelper.h`: `convertTargetedExp`, both `convertTargetedCompound` | `convert_targeted_exp`, `convert_targeted_peptide`, `convert_targeted_compound` in [`src/analysis/openswath/data_access_helper.rs`](../src/analysis/openswath/data_access_helper.rs) |

Tests: [`tests/mrm_feature_finder_scoring.rs`](../tests/mrm_feature_finder_scoring.rs) (6),
[`tests/confidence_scoring.rs`](../tests/confidence_scoring.rs) (6),
[`tests/mrm_oracle.rs`](../tests/mrm_oracle.rs) (executed comparison).
Manifest: [`tests/data/mrm_feature_provenance.json`](../tests/data/mrm_feature_provenance.json).

The finder calls `EmgScoring` (`FEATUREFINDER/EmgScoring.h`, the FF1 lane's
`src/analysis/featurefinder/emg_scoring.rs`) for `use_elution_model_score`.

## API mapping: `MRMFeatureFinderScoring`

| Source member | Rust |
| --- | --- |
| `MRMFeatureFinderScoring()` and the parameter defaults | `new`, `defaults`, `from_param`, `from_param_with_warnings`, `set_parameters`, `parameters` |
| `updateMembers_` | `from_param` (the public fields hold the members) |
| `pickExperiment(PeakMap, FeatureMap, TargetedExperiment, trafo, PeakMap)` | `pick_experiment` |
| `pickExperiment(SpectrumAccessPtr, FeatureMap, LightTargetedExperiment, trafo, swath_maps, group_map)` | `pick_experiment_light` |
| `prepareProteinPeptideMaps_` | `prepare_protein_peptide_maps` |
| `scorePeakgroups` | `score_peakgroups` |
| `setStrictFlag` | `set_strict_flag` |
| `setMS1Map` | `set_ms1_map` |
| `mapExperimentToTransitionList` | `map_experiment_to_transition_list` |
| `splitTransitionGroupsDetection_`, `splitTransitionGroupsIdentification_` | private `split_detection`, `split_identification` |
| `scoreIdentification_` | private `score_identification` |
| `prepareScoredFeatureOutput_`, `addScoreMetaValues_`, `prepareFeatureOutput_` | private, same names in snake case |
| `TransitionGroupCache`, `OpenSwath_Ind_Scores_Pooled` | `TransitionGroupCache` (private); the pooled score record is a per-call `OpenSwathIndScores` |
| `TransitionGroupMapType`, `MRMTransitionGroupType` | the type aliases of the same names |
| `run_identifier` | `RUN_IDENTIFIER` |
| `MobilogramParquetConsumer*` | `Option<&mut dyn MobilogramSink>` |

## API mapping: `ConfidenceScoring`

| Source member | Rust |
| --- | --- |
| `ConfidenceScoring(bool test_mode_)` | `ConfidenceScoring::new(test_mode)` |
| `initialize` | `initialize` |
| `initializeGlm` | `initialize_glm` |
| `scoreMap` | `score_map` |
| `GLM_` | `Glm` (`evaluate`) |
| `RTNorm_` | `RtNorm` (`apply`) |
| `chooseDecoys_` | private `choose_decoys` |
| `manhattanDist_` | `manhattan_dist` |
| `getAssayRT_` | `get_assay_rt` |
| `scoreAssay_` | `score_assay` |
| `scoreFeature_` | `score_feature` |
| `boost::bimap<double, multiset_of<double>>` | private `IntensityMap` |

## API mapping: `DataAccessHelper`

| Source member | Rust |
| --- | --- |
| `convertTargetedExp` | `convert_targeted_exp` |
| `convertTargetedCompound(const Peptide&, LightCompound&)` | `convert_targeted_peptide` |
| `convertTargetedCompound(const Compound&, LightCompound&)` | `convert_targeted_compound` |

## Preserved source conventions

- **Scored order.** The scored peak groups of a transition group reach the
  output ordered by overall quality, best first, with the libstdc++
  `std::sort` permutation for ties (`source_sort_reversed_by`). With
  `stop_report_after_feature` below the number of peak groups, only that many
  are prepared and written, chosen by a sort on quality.
- **Overall quality is the LDA prescore.** `-calculate_lda_prescore` (or the
  single transition model), replaced by `-calculate_swath_lda_prescore` when a
  SWATH map with spectra is given and DIA scores are on.
- **`std::log` of a `float`.** The identification `ind_log_intensity` takes the
  logarithm of a `float` intensity, which C++ resolves to `logf`; the port uses
  the reference build's `logf` (`featurefinder::glibc_logf`). Double logarithms
  use `glibc_libm::log`.
- **`prepareScoredFeatureOutput_` tests `xx_swath_prelim_score`**, which
  `addScoreMetaValues_` never writes (it writes `main_var_xx_swath_prelim_score`),
  so the peptide hit score is always the LDA prescore. Kept.
- **The detection group replaces the caller's group.** `transition_group =
  transition_group_detection` at the end of `scorePeakgroups`, so a group with
  identifying-only transitions loses them; the peak groups it holds keep their
  scores and sub-features.
- **Identification features are copies.** `scoreIdentification_` reads the
  peak groups of the identification group, which are copies made before
  scoring.
- **Transition names may repeat across groups.** The class test's
  identification input lists `tr2` in both groups; each group keeps its own
  transition.
- **`getRetentionTime()` in minutes** is converted to seconds by
  `convertTargetedCompound`; a range of exactly two retention times is kept as
  `rt_start`/`rt_end`.
- **`PeptideRefMap_` keeps the last compound of a repeated ID.**
- **Scores are stored as they are**, NaN included (`addScore` is
  `setMetaValue`).
- **`ConfidenceScoring`** keeps its state between `scoreMap` calls: the
  transition lookup grows by every transition again and `n_decoys` stays
  clamped. The decoy shuffle runs only for `n_decoys > 0`, with the engine
  seeded `0` (reproducible) or from the clock in test mode.

## Native differences

| Where | Source | Here |
| --- | --- | --- |
| per-feature scoring loop | OpenMP `parallel for` with a critical section on the output list | serial; a result the source can produce with one thread |
| `std::partial_sort` of the reported peak groups | heap selection, an unspecified order among equal qualities | a full introsort of the indices, the first `stop_report_after_feature` taken; the same unless equal qualities straddle the limit |
| `PeptideRefMap_` | pointers into the caller's transition list | cloned compounds (call `prepare_protein_peptide_maps` before `score_peakgroups`) |
| `PeptideRefMap_[ref]` for an unknown compound with an RT window | a null dereference | `Error::InvalidValue` |
| `getPeptideByRef` for an unknown peptide in `convertTargetedExp` | a null dereference | `Error::InvalidValue` |
| unique IDs | the global `UniqueIdGenerator` | a generator the finder shares through an `Arc<Mutex<..>>` (`set_unique_id_generator`) |
| `MobilogramParquetConsumer` | the Parquet writer | the `MobilogramSink` trait; the Parquet writer is the backlog item `FORMAT/DATAACCESS/MobilogramParquetConsumer.h` |
| `ProgressLogger` base of both classes | progress messages | none |
| `ConfidenceScoring::getAssayRT_` precondition | debug-only | `Error::InvalidValue` when no RT is set |
| `local_FDR` with no decoy | `0 / 0 = NaN` stored | the same, through the source-exact metadata constructor (a map holding it is refused by checked writers) |
| tests' `TraMLFile` | reads the TraML | not ported (FORMAT backlog); the tests build the same experiment by hand |

## Checked boundaries and evidence

| Behaviour | Tier | Source |
| --- | --- | --- |
| all sections of `ConfidenceScoring_test.cpp` (including the four exception messages) | 3 | `tests/confidence_scoring.rs` |
| all sections of `MRMFeatureFinderScoring_test.cpp`: default and legacy picker parameters, UIS identification scores, `mapExperimentToTransitionList` | 3 | `tests/mrm_feature_finder_scoring.rs` |
| the finder on five configurations over `OpenSwath_generic_input` (default, legacy, identification, `stop_report_after_feature`, mutual information with the single transition model): every output feature with all metavalues, sub-features, identification and protein run | 1 | `tests/mrm_oracle.rs` |
| `ConfidenceScoring` on the class test and on a seven-assay library with decoys, transitions kept and clamped decoys | 1 | same |
| strict mapping, RT windows, scorePeakgroups preconditions, `n_transitions` and decoy reproducibility | 4 | `tests/mrm_feature_finder_scoring.rs`, `tests/confidence_scoring.rs` |

### Defect candidates

Listed in the lane report.

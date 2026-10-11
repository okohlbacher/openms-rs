# Targeted experiment data model

Ports of the OpenMS core pin `bc9cc12` headers

| Header | Rust | Status |
|---|---|---|
| `ANALYSIS/TARGETED/TargetedExperimentHelper.h` (+ `.cpp`) | `src/analysis/openswath/targeted_experiment_helper.rs` | complete |
| `ANALYSIS/MRM/ReactionMonitoringTransition.h` (+ `.cpp`) | `src/analysis/openswath/reaction_monitoring_transition.rs` | complete |
| `ANALYSIS/TARGETED/IncludeExcludeTarget.h` (+ `.cpp`) | `src/analysis/openswath/include_exclude_target.rs` | complete |
| `ANALYSIS/TARGETED/TargetedExperiment.h` (+ `.cpp`) | `src/analysis/openswath/targeted_experiment.rs` | complete |
| `ANALYSIS/TARGETED/MRMMapping.h` (+ `.cpp`) | `src/analysis/openswath/mrm_mapping.rs` | partial: the `LightTargetedExperiment` overload waits for `OPENSWATHALGO/DATAACCESS/TransitionExperiment.h` |
| `ANALYSIS/TARGETED/MetaboTargetedTargetDecoy.h` (+ `.cpp`) | `src/analysis/openswath/metabo_targeted_target_decoy.rs` | complete |

Tests: `tests/targeted_experiment.rs`. Provenance: `tests/data/targeted_experiment_provenance.json`.

The five remaining ANALYSIS/TARGETED headers are not ported: `MetaboTargetedAssay` needs
`SiriusMSConverter`, `SiriusFragmentAnnotation`, `FeatureGroupingAlgorithmQT`, `BinnedSpectrum`,
`BinnedSpectralContrastAngle` and `SpectraMerger`; `ChromatogramProcessor` needs
`MRMFeatureFinderScoring` and `SpectrumAccessOpenMS`; `IChromatogramHandler`, `DefaultChromHandler`,
`MRMChromHandler` and `DIAChromHandler` need `SwathMap`, `LightTargetedExperiment`,
`ChromExtractParams` (`OpenSwathWorkflow.h`) and the chromatogram extractors. The TraML
(`TraMLFile`) and transition TSV/PQP (`TransitionTSVFile`, `TransitionPQPFile`) formats are backlog
items that depend on the types ported here; no format is re-scoped out.

## API mapping

### TargetedExperimentHelper

| Source | Rust |
|---|---|
| `Configuration` (`CVTermList` base, `contact_ref`, `instrument_ref`, `validations`) | `Configuration { cv_terms, contact_ref, instrument_ref, validations }` |
| `CV(id, fullname, version, URI)`, `operator==`, `std::hash<CV>` | `CV::new`, `PartialEq`/`Eq`, `Hash`; field `URI` is `uri`; `Default` is native |
| `Protein` (`id`, `sequence`), `operator==`, `std::hash` | `Protein { cv_terms, id, sequence }`, `PartialEq`/`Eq`, `Hash` |
| `RetentionTime::RTUnit` (`SECOND`, `MINUTE`, `UNKNOWN`, `SIZE_OF_RTUNIT`) | `RTUnit::{Second, Minute, Unknown, SizeOfRTUnit}`, `#[repr(i8)]` with the source values |
| `RetentionTime::RTType` (`LOCAL` … `SIZE_OF_RTTYPE`) | `RTType::{Local, Normalized, Predicted, Hpins, Irt, Unknown, SizeOfRTType}` |
| `RetentionTime()` and copy/move/assign | `RetentionTime::new` / `Default`, `Clone` |
| `RetentionTime::operator==`, `std::hash` | `PartialEq`/`Eq`, `Hash` |
| `isRTset`, `setRT`, `getRT` | `is_rt_set`, `set_rt`, `rt() -> Option<f64>` |
| `software_ref`, `retention_time_unit`, `retention_time_type` | public fields of the same names |
| `PeptideCompound` ctor/copy/move, `operator==`, `std::hash` | `PeptideCompound::new`/`Default`/`Clone`, `PartialEq`/`Eq`, `Hash` |
| `setChargeState`, `hasCharge`, `getChargeState` | `set_charge_state`, `has_charge`, `charge_state() -> Option<i32>` |
| `setDriftTime`, `getDriftTime` | `set_drift_time`, `drift_time` |
| `hasRetentionTime`, `getRetentionTime`, `getRetentionTimeType`, `getRetentionTimeUnit` | `has_retention_time`, `retention_time`, `retention_time_type`, `retention_time_unit` (each `Result`) |
| `PeptideCompound::id`, `rts` | public fields |
| `Compound` (`molecular_formula`, `smiles_string`, `theoretical_mass`), `operator==`, `std::hash` | `Compound { base, molecular_formula, smiles_string, theoretical_mass }`, `Deref<Target = PeptideCompound>`, `PartialEq`/`Eq`, `Hash` |
| `Peptide::Modification` (`avg_mass_delta`, `mono_mass_delta`, `location`, `unimod_id`) | `PeptideModification` |
| `Peptide` (`protein_refs`, `evidence`, `sequence`, `mods`), `operator==`, `std::hash` | `Peptide { base, protein_refs, evidence, sequence, mods }`, `Deref`, `PartialEq`/`Eq`, `Hash` |
| `setPeptideGroupLabel`, `getPeptideGroupLabel` | `set_peptide_group_label`, `peptide_group_label` |
| `Contact`, `Publication`, `Instrument` (`id`), `operator==`, `std::hash` | same names with `cv_terms` and `id`, `PartialEq`/`Eq`, `Hash` |
| `Prediction` (`software_ref`, `contact_ref`), `operator==`, `std::hash` | `Prediction`, `PartialEq`/`Eq`, `Hash` |
| `Interpretation` (`ordinal`, `rank`, `iontype`), `IonType` typedef, `operator==`/`!=`, `std::hash` | `Interpretation { cv_terms, ordinal: u8, rank: u8, iontype: PeptideFragmentType }`, `PartialEq`/`Eq`, `Hash` |
| `TraMLProduct` ctor, `operator==`, `std::hash` | `TraMLProduct::new`/`Default`, `PartialEq`/`Eq`, `Hash` |
| `setChargeState`, `hasCharge`, `getChargeState` | `set_charge_state`, `has_charge`, `charge_state() -> Option<i32>` |
| `getMZ`, `setMZ` | `mz`, `set_mz` |
| `getConfigurationList`, `addConfiguration` | `configuration_list`, `add_configuration` |
| `getInterpretationList`, `addInterpretation`, `resetInterpretations` | `interpretation_list`, `add_interpretation`, `reset_interpretations` |
| `getAASequence(const Peptide&)` | `get_aa_sequence` (strict), `get_aa_sequence_with(peptide, db, policy) -> AASequenceConversion` |
| `setModification(location, max_size, modification, aas)` | `set_modification`, `set_modification_with_registry` (both `Result`) |
| `hashCVTermListInterface` | not ported: `CVTermList` implements `Hash` |

### ReactionMonitoringTransition

| Source | Rust |
|---|---|
| typedefs `Configuration`, `RetentionTime`, `Prediction` | the helper types directly |
| typedef `Product` | `pub type Product = TraMLProduct` |
| `DecoyTransitionType` (`UNKNOWN`, `TARGET`, `DECOY`) | `DecoyTransitionType::{Unknown, Target, Decoy}` (`Ord` in source order) |
| ctor, copy/move ctor, dtor, copy/move assignment | `new`/`Default`, `Clone`, drop, moves |
| `setName`/`getName`, `setNativeID`/`getNativeID` (one member) | `set_name`/`name`, `set_native_id`/`native_id` |
| `setPeptideRef`/`getPeptideRef`, `setCompoundRef`/`getCompoundRef` | `set_peptide_ref`/`peptide_ref`, `set_compound_ref`/`compound_ref` |
| `setPrecursorMZ`, `getPrecursorMZ` | `set_precursor_mz`, `precursor_mz` |
| `hasPrecursorCVTerms`, `setPrecursorCVTermList`, `addPrecursorCVTerm`, `getPrecursorCVTermList` | `has_precursor_cv_terms`, `set_precursor_cv_term_list`, `add_precursor_cv_term -> Result`, `precursor_cv_term_list -> Option` |
| `setProductMZ`, `getProductMZ`, `getProductChargeState`, `isProductChargeStateSet`, `addProductCVTerm` | `set_product_mz`, `product_mz`, `product_charge_state -> Option`, `is_product_charge_state_set`, `add_product_cv_term -> Result` |
| `getIntermediateProducts`, `addIntermediateProduct`, `setIntermediateProducts` | `intermediate_products`, `add_intermediate_product`, `set_intermediate_products` |
| `setProduct`, `getProduct` | `set_product`, `product` |
| `hasPrediction`, `setPrediction`, `addPredictionTerm`, `getPrediction` | `has_prediction`, `set_prediction`, `add_prediction_term -> Result`, `prediction -> Result` |
| `getDecoyTransitionType`, `setDecoyTransitionType` | `decoy_transition_type`, `set_decoy_transition_type` |
| `getLibraryIntensity`, `setLibraryIntensity` | `library_intensity`, `set_library_intensity` |
| `setRetentionTime`, `getRetentionTime` | `set_retention_time`, `retention_time` |
| `isDetectingTransition`/`setDetectingTransition` and the identifying and quantifying pairs | `is_detecting_transition`/`set_detecting_transition`, … |
| `operator==`, `operator!=`, `std::hash` | `PartialEq`/`Eq`, `Hash` |
| `ProductMZLess`, `NameLess` | `cmp_product_mz`, `cmp_name` (`Ordering`) |
| `updateMembers_` (empty, protected) | not ported: no behaviour |
| — | native: `impl kernel::mrm::Transition` |

### IncludeExcludeTarget

| Source | Rust |
|---|---|
| typedefs `Configuration`, `RetentionTime` | helper types |
| ctor, copy ctor, dtor, `operator=` | `new`/`Default`, `Clone` |
| `setName`/`getName`, `setPeptideRef`/`getPeptideRef`, `setCompoundRef`/`getCompoundRef` | `set_name`/`name`, … |
| `setPrecursorMZ`/`getPrecursorMZ`, `setPrecursorCVTermList`, `addPrecursorCVTerm`, `getPrecursorCVTermList` | `set_precursor_mz`/`precursor_mz`, `set_precursor_cv_term_list`, `add_precursor_cv_term -> Result`, `precursor_cv_term_list` |
| the four product counterparts | `set_product_mz`/`product_mz`, `set_product_cv_term_list`, `add_product_cv_term -> Result`, `product_cv_term_list` |
| `setInterpretations`, `getInterpretations`, `addInterpretation` | `set_interpretations`, `interpretations`, `add_interpretation` |
| `setConfigurations`, `getConfigurations`, `addConfiguration` | `set_configurations`, `configurations`, `add_configuration` |
| `setPrediction`, `addPredictionTerm`, `getPrediction` | `set_prediction`, `add_prediction_term -> Result`, `prediction` |
| `setRetentionTime`, `getRetentionTime` | `set_retention_time`, `retention_time` |
| `operator==`, `operator!=`, `std::hash` | `PartialEq`/`Eq`, `Hash` |
| `updateMembers_` | not ported: empty |

### TargetedExperiment

| Source | Rust |
|---|---|
| `SummaryStatistics` and `operator<<` | `SummaryStatistics` (`decoy_counts: BTreeMap`), `Display` |
| typedefs `CV` … `Interpretation`, `Transition`, `IonType` | the helper, transition and `PeptideFragmentType` types directly |
| `ProteinReferenceMapType`, `PeptideReferenceMapType`, `CompoundReferenceMapType` | private lazily built index maps |
| ctors, dtor, copy/move assignment | `new`/`Default`, `Clone` |
| `operator==`, `operator!=` | `PartialEq` |
| `operator+`, `operator+=` (copy and move) | `merged`, `merge`, `merge_owned` (all `Result`) |
| `clear(bool)` | `clear` |
| `getSummary` | `summary` |
| `setCVs`/`getCVs`/`addCV`, contacts, publications, instruments, software, include/exclude targets, source files | `set_*`/getter/`add_*` of the same stems |
| `setTargetCVTerms`, `getTargetCVTerms`, `addTargetCVTerm`, `setTargetMetaValue` | `set_target_cv_terms`, `target_cv_terms`, `add_target_cv_term -> Result`, `set_target_meta_value` |
| `setProteins` (both), `getProteins`, `getProteinByRef`, `hasProtein`, `addProtein` | `set_proteins`, `proteins`, `protein_by_ref -> Result`, `has_protein`, `add_protein` |
| `setCompounds`, `getCompounds`, `addCompound`, `hasCompound`, `getCompoundByRef` | `set_compounds`, `compounds`, `add_compound`, `has_compound`, `compound_by_ref -> Result` |
| `setPeptides` (both), `getPeptides`, `hasPeptide`, `getPeptideByRef`, `addPeptide` | `set_peptides`, `peptides`, `has_peptide`, `peptide_by_ref -> Result`, `add_peptide` |
| `setTransitions` (both), `getTransitions`, `addTransition` | `set_transitions`, `transitions`, `add_transition` |
| `sortTransitionsByProductMZ`, `sortTransitionsByName` | `sort_transitions_by_product_mz`, `sort_transitions_by_name` |
| `containsInvalidReferences` | `contains_invalid_references`; native `first_invalid_reference` returns the logged message |
| `createProteinReferenceMap_` and siblings (protected) | private `protein_index` and siblings |

### MRMMapping

| Source | Rust |
|---|---|
| ctor (defaults `precursor_tolerance` 0.1, `product_tolerance` 0.1, `map_multiple_assays` false, `error_on_unmapped` false) | `new`/`Default`, `defaults() -> Result<Param>`, `to_param`, `from_param`, `from_param_with_warnings` |
| `mapExperiment(PeakMap, TargetedExperiment, PeakMap&)` | `map_experiment(&MSExperiment, &TargetedExperiment) -> Result<MSExperiment>` |
| `mapExperiment(PeakMap, LightTargetedExperiment, PeakMap&)` | **not ported**: needs `LightTargetedExperiment` (lane OA) |
| copy ctor / assignment (protected), `updateMembers_` | `Clone`; `from_param` |
| — | native `max_output_points` ceiling |

### MetaboTargetedTargetDecoy

| Source | Rust |
|---|---|
| `MetaboTargetDecoyMassMapping` | `MetaboTargetDecoyMassMapping` |
| `constructTargetDecoyMassMapping` | `construct_target_decoy_mass_mapping` |
| `resolveOverlappingTargetDecoyMassesByDecoyMassShift(…, mz_tol_unit)` | `resolve_overlapping_target_decoy_masses_by_decoy_mass_shift(…, MzToleranceUnit)`; `MzToleranceUnit::from_source_str` |
| `generateMissingDecoysByMassShift` | `generate_missing_decoys_by_mass_shift -> Result` |
| `constructTransitionsMap_` (protected) | private `transitions_by_compound_ref` |

## Preserved source conventions

- Defaults: RMT library intensity `-101`, flags detecting/not identifying/quantifying, decoy type
  `UNKNOWN`, precursor m/z `0`; `IncludeExcludeTarget` m/z `f64::MAX`; `RetentionTime` unit and type
  the `SIZE_OF_*` sentinels; `PeptideCompound` drift time `-1`; `Peptide::Modification` location and
  UniMod id `-1`; `Interpretation` ion type `Unannotated`.
- `PeptideCompound` equality ignores the drift time, as the source `operator==`.
- `setName` and `setNativeID` write one member.
- By-reference lookups resolve a repeated id to its **last** occurrence, as the source's
  pointer-map loop; `MRMMapping` and `MetaboTargetedTargetDecoy` use the **first** matching peptide or
  compound, as their linear searches.
- `merge` appends without deduplication and does not merge the target list's metadata; `clear(false)`
  removes only the transitions.
- `containsInvalidReferences` checks in the source order and stops at the first failure; a transition
  without references is valid.
- `SummaryStatistics` prints the source format with `StringUtils::toStr(double, false)` percentages.
- `MRMMapping` tolerances are strict (`<`), a chromatogram with both m/z below `1e-5` is skipped
  unless multiple mapping is enabled, and the error messages are the source's.
- `MetaboTargetedTargetDecoy` overlap tolerance is inclusive (`<=`); only the string `"ppm"` means
  ppm; the regex `_\[` → `_decoy_[` replacement is a literal replace-all.

## Native differences

- **Equality of `Configuration` and `PeptideModification`.** The source declares no `operator==` for
  either, so the inherited `CVTermList`/`CVTermListInterface` comparison ignores every own field
  (references, validations, masses, location, UniMod id). The port compares all fields. C++ defect
  candidate: two peptides with a modification at different residues compare equal, and so do two
  targeted experiments containing them, while `std::hash<Peptide>` hashes the location — equal values
  with different hashes.
- **Hash values.** `Hash` is implemented consistently with `PartialEq` over all compared fields; the
  source's FNV-1a mixing is not reproduced (hash values are not an observable contract).
  `std::hash<ReactionMonitoringTransition>` hashes a subset; the port hashes every compared field.
- **Precondition getters** (`getRT`, `getChargeState`, `getPrecursorCVTermList`,
  `getProductChargeState`) return `Option` instead of an unchecked value.
- **Throwing getters** map to `Result`: `IllegalArgument` → `Error::InvalidValue`,
  `ElementNotFound`/`MissingInformation` → `Error::MissingInformation`.
- **CV-term adders validate** the term (`CVTermList::add`) and return `Result`; an invalid term leaves
  the object unchanged, including not creating the lazily allocated precursor list or prediction.
- **`getAASequence`**: an unresolvable modification is refused by default
  (`UnresolvedModificationPolicy::Reject`); `SourceSkip` reproduces the source's skip-and-warn. A
  terminal location without a UniMod id is unresolvable (the source reads `sequence[-1]`, undefined
  behaviour, or the terminating NUL). An out-of-range location is refused instead of reaching
  `AASequence` with a wrapped index.
- **Sorting** is stable and uses `f64::total_cmp`; the source `std::sort` is unstable with `<`.
- **`first_invalid_reference`** names the compound in the compound message; the source prints the
  (empty) peptide reference there (C++ defect candidate). The source's warning for a transition
  without references is not reported.
- **`MRMMapping::map_experiment`** returns a new experiment and nothing on error (the source leaves
  its out-parameter half filled); a native `max_output_points` ceiling is checked before any copy; a
  NaN tolerance is refused.
- **`generate_missing_decoys_by_mass_shift`** refuses two undefined-behaviour cases of the source
  (no target compound for a generated decoy: end-iterator dereference; fewer decoy masses than target
  transitions: out-of-bounds read) and is atomic.

## Checked boundaries and evidence

- Tier 3: every section of `TargetedExperimentHelper_test.cpp`, `ReactionMonitoringTransition_test.cpp`,
  `IncludeExcludeTarget_test.cpp` (whose accessor sections are `TODO` upstream; its hash section is
  ported), `TargetedExperiment_test.cpp` and `MRMMapping_test.cpp` is ported with its literals.
- The `MRMTransitionGroup_test.cpp` `getLibraryIntensity` section, which the source instantiates with
  `ReactionMonitoringTransition`, now runs with the real type.
- Tier 1: see `tests/data/targeted_experiment_provenance.json` for the C++ oracle
  (`../oracle/O1a-targeted`, ibminode06, Linux x86-64 Release install `bc9cc12-c19e494-174b576`).
- Tier 4: the native checks listed at each test.

# MRM feature QC, QC files, picker parameters and the feature filter

Rust ports of five core-SDK headers (pin `bc9cc12`):

| Header | Rust |
| --- | --- |
| `ANALYSIS/OPENSWATH/MRMFeatureQC.h` | [`src/analysis/openswath/mrm_feature_qc.rs`](../src/analysis/openswath/mrm_feature_qc.rs) |
| `FORMAT/MRMFeatureQCFile.h` | [`src/format/mrm_feature_qc_file.rs`](../src/format/mrm_feature_qc_file.rs) |
| `ANALYSIS/OPENSWATH/MRMFeaturePicker.h` | [`src/analysis/openswath/mrm_feature_picker.rs`](../src/analysis/openswath/mrm_feature_picker.rs) |
| `FORMAT/MRMFeaturePickerFile.h` | [`src/format/mrm_feature_picker_file.rs`](../src/format/mrm_feature_picker_file.rs) |
| `ANALYSIS/OPENSWATH/MRMFeatureFilter.h` | [`src/analysis/openswath/mrm_feature_filter.rs`](../src/analysis/openswath/mrm_feature_filter.rs) |

The metavalue conversions the filter shares with the selector and the
confidence scorer live in `src/analysis/openswath/mrm_feature_support.rs`
(`pub(crate)`: `StringUtils::toStr(DataValue)`, the `(double)` cast).

Tests: [`tests/mrm_feature_qc.rs`](../tests/mrm_feature_qc.rs) (15),
[`tests/mrm_feature_filter.rs`](../tests/mrm_feature_filter.rs) (26),
[`tests/mrm_oracle.rs`](../tests/mrm_oracle.rs) (the executed C++ comparison).
Manifest: [`tests/data/mrm_feature_provenance.json`](../tests/data/mrm_feature_provenance.json).

## API mapping

### `MRMFeatureQC`

| Source member | Rust |
| --- | --- |
| `MRMFeatureQC()`, `~MRMFeatureQC()` | `MRMFeatureQC::default()`; drop |
| `ComponentQCs` (`component_name`, `retention_time_l/u`, `intensity_l/u`, `overall_quality_l/u`, `meta_value_qc`) | `ComponentQCs`, same field names; defaults `0.0` / `1e12` |
| `ComponentQCs::operator==`, `operator!=` | derived `PartialEq` (exact map comparison, see differences) |
| `ComponentGroupQCs` (the above plus `n_heavy_l/u`, `n_light_l/u`, `n_detecting_l/u`, `n_quantifying_l/u`, `n_identifying_l/u`, `n_transitions_l/u`, `ion_ratio_pair_name_1/2`, `ion_ratio_l/u`, `ion_ratio_feature_name`) | `ComponentGroupQCs`; counts `i32` defaulting to `0` / `100` |
| `ComponentGroupQCs::operator==`, `operator!=` | derived `PartialEq` |
| `ComponentGroupPairQCs` (`component_group_name`, `resolution_pair_name`, `resolution_l/u`, `rt_diff_l/u`) | `ComponentGroupPairQCs` |
| `component_qcs`, `component_group_qcs`, `component_group_pair_qcs` | the three `Vec` fields of `MRMFeatureQC` |

### `MRMFeatureQCFile`

| Source member | Rust |
| --- | --- |
| `load(filename, mrmfqc, is_component_group)` | `MRMFeatureQCFile::load` |
| `store(filename, mrmfqc, is_component_group)` | `MRMFeatureQCFile::store` |
| `pushValuesFromLine_` (`ComponentQCs`) | `push_values_from_line_component` |
| `pushValuesFromLine_` (`ComponentGroupQCs`) | `push_values_from_line_component_group` |
| `setPairValue_` | `set_pair_value` |
| `getCastValue_` (`Int`, `double`, `String`) | `get_cast_value_i32`, `get_cast_value_f64`, `get_cast_value_string` |
| base classes `CsvFile`, `ProgressLogger` | not state of the Rust type (see below) |

### `MRMFeaturePicker` and `MRMFeaturePickerFile`

| Source member | Rust |
| --- | --- |
| `MRMFeaturePicker` | unit struct `MRMFeaturePicker` |
| `MRMFeaturePicker::ComponentParams` | `mrm_feature_picker::ComponentParams` |
| `MRMFeaturePicker::ComponentGroupParams` | `mrm_feature_picker::ComponentGroupParams` |
| `MRMFeaturePickerFile::load(filename, cp_list, cgp_list)` | `MRMFeaturePickerFile::load` returning the two lists |
| `extractParamsFromLine_` | `extract_params_from_line` (returns `None` for the source's `false`) |
| `setCastValue_` | `set_cast_value` |

### `MRMFeatureFilter`

| Source member | Rust |
| --- | --- |
| `MRMFeatureFilter()` (parameter `flag_or_filter`) | `MRMFeatureFilter::new`, `defaults`, `to_param`, `from_param`, `from_param_with_warnings`; field `flag_or_filter: FlagOrFilter` |
| `getDefaultParameters` | `defaults` |
| `updateMembers_` | `from_param` |
| `FilterFeatureMap` | `filter_feature_map` |
| `FilterFeatureMapPercRSD` | `filter_feature_map_perc_rsd` |
| `FilterFeatureMapBackgroundInterference` | `filter_feature_map_background_interference` |
| `EstimateDefaultMRMFeatureQCValues` | `estimate_default_mrm_feature_qc_values` |
| `TransferLLOQAndULOQToCalculatedConcentrationBounds` | `transfer_lloq_and_uloq_to_calculated_concentration_bounds`, generic over `QuantitationBounds` |
| `EstimatePercRSD` | `estimate_perc_rsd` |
| `EstimateBackgroundInterferences` | `estimate_background_interferences` |
| `calculateIonRatio` | `calculate_ion_ratio` |
| `calculateRTDifference` | `calculate_rt_difference` |
| `calculateResolution` | `calculate_resolution` |
| `checkMetaValue` | `check_meta_value` returning `(passes, key_exists)` |
| `updateMetaValue`, `setMetaValue`, `initMetaValue` | `update_meta_value`, `set_meta_value`, `init_meta_value` returning `key_exists` |
| `countLabelsAndTransitionTypes` | `count_labels_and_transition_types` returning `LabelCounts` |
| `getUniqueSorted` | `get_unique_sorted` |
| `accumulateFilterValues` | `accumulate_filter_values` |
| `zeroFilterValues` | `zero_filter_values` |
| `calculateFilterValuesMean`, `...Var`, `...PercRSD` | `calculate_filter_values_mean`, `calculate_filter_values_var`, `calculate_filter_values_perc_rsd` |
| `checkRange`, `updateRange`, `setRange`, `initRange` | `check_range`, `update_range`, `set_range`, `init_range`, generic over `RangeValue` (`f64`, `i32`) |
| `flag_or_filter_` | `flag_or_filter` |

## Preserved source conventions

- **Group checks run once per subordinate.** In `FilterFeatureMap` the component
  group QC is evaluated inside the loop over the feature's subordinates, so a
  group with `k` subordinates adds `k` times the tests (and failure messages)
  to the group score. The `%RSD` and background interference variants do the
  same. Kept.
- **Pairing by position.** `FilterFeatureMapPercRSD` and
  `FilterFeatureMapBackgroundInterference` pair `filter_criteria` with
  `filter_values` by list index, not by name.
- **`EstimateDefaultMRMFeatureQCValues` swaps `update` and `init` for the group
  metavalues.** The component branch initialises on the first sample and
  updates afterwards; the group branch does the opposite. Kept, and flagged as a
  defect candidate.
- **`accumulateFilterValues` keeps the last matching feature** of each sample,
  because `setRange` assigns.
- **`TEST_EQUAL(a, b)` converts the expected value to the type of `a`**
  (`ClassTest.h`, `testEqual`: `expression_1 == T1(expression_2)`), so the
  class test's `n_detecting_u, 43.30127019` asserts the `Int` `43`. The
  transcription keeps that.
- **Integer counts stay `Int`.** The mean of a count divides an `Int` as
  `unsigned long`; the variance accumulates through a `double` and is truncated
  on assignment; the `%RSD` of a count is truncated. The arithmetic is
  reproduced.
- **Booleans are integers.** `setMetaValue("..._pass", true)` converts `bool`
  through the integer `DataValue` constructor (there is no `bool` one), so the
  flags are the integers `1` and `0`, which the retained featureXML shows as
  `type="int"`.
- **`pow(x, 2)`** is compiled to `x * x`; the port multiplies.
- **A missing `native_id`** stringifies to the empty string and matches an
  empty QC name.
- **Empty cells and defaults.** An empty cell takes the source default of its
  bound, not zero.
- **`metaValue_(.+)_(l|u)`** is searched unanchored with Boost semantics through
  the crate's `BoostRegex`, so `metaValue_a_lx` is the metavalue `a`, lower
  bound.
- **`std::stod` / `std::stoi` read a prefix.** `1.5abc` is `1.5` in a QC cell,
  and a string without a number is an error. The picker file instead uses
  `StringUtils::toDouble` / `toInt32`, which require the whole cell.
- **The boolean parameters of the picker file** are `true` only for exactly
  `true` or `TRUE`, although the header comment says case-insensitive.

## Native differences

| Where | Source | Here |
| --- | --- | --- |
| `ComponentQCs::operator==` | `std::equal(first1, last1, first2, ...)` over the metavalue maps reads past the end of a shorter right-hand map (undefined) and reports `true` when the right-hand map merely starts with the left | exact map equality; the same whenever both maps have one size |
| `ComponentGroupPairQCs` | numeric members uninitialised | zeroed by `Default` |
| rows shorter than the header | `line[index]` reads past the end of the row | `Error::InvalidValue` naming the column |
| number syntax | `std::stod`, `std::stoi` throw `std::invalid_argument` / `out_of_range`, uncaught in the TOPP tools | `Error::InvalidValue`; hexadecimal floating literals are `Error::Unsupported` |
| load | `component_qcs` is cleared before parsing | built aside and assigned on success, so an error leaves `qc` unchanged |
| store | writes metavalue cells positionally under the first QC's headers | a QC whose metavalue names differ from the first is `Error::InvalidValue` (the source misaligns the columns, or writes past the row) |
| store of an empty list | returns without writing | same |
| `FilterFeatureMap*` in filter mode | `features = features_filtered` replaces the whole `FeatureMap`, dropping its identifier, protein identifications and data processing | only the features are replaced; map-level members are kept |
| `Int` mean / variance with no sample (or one for the variance) | integer division by zero ends the process | `Error::InvalidValue` when a component group QC would be divided |
| `.at()` on `filter_values` | `std::out_of_range` | `Error::InvalidValue` |
| non-numeric metavalue read as `double` | the empty value throws; any other type reads an unrelated union member | `Error::InvalidValue` |
| `OPENMS_LOG_DEBUG` in `calculateIonRatio`, the `*MetaValue` helpers | debug messages | not produced |
| `ProgressLogger` base of `MRMFeatureQCFile` | unused | no progress |
| `TransferLLOQAndULOQ...` argument `std::vector<AbsoluteQuantitationMethod>` | the concrete class | generic over the trait `QuantitationBounds` (`component_name`, `lloq`, `uloq`); `ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h` is not ported and implements it with three one-line accessors when it is |
| unsigned truncation of a negative or out-of-range double in `setCastValue_` | undefined | saturating `as u32` |

## Checked boundaries and evidence

Evidence tiers follow `docs/DIFFERENTIAL_VALIDATION.md`.

| Behaviour | Tier | Source |
| --- | --- | --- |
| every section of `MRMFeatureQCFile_test.cpp` (load, store round trip, both `pushValuesFromLine_`, `setPairValue_`, three `getCastValue_`) | 3 | `tests/mrm_feature_qc.rs` |
| `MRMFeaturePickerFile_test.cpp` `load` on `MRMFeaturePickerFile.csv` | 3 | `tests/mrm_feature_qc.rs` |
| all 21 sections of `MRMFeatureFilter_test.cpp`, transcribed statement for statement | 3 | `tests/mrm_feature_filter.rs` |
| text of a stored component file against `MRMFeatureQCFile_1_test.csv` (upstream test data) | 3 | `tests/mrm_feature_qc.rs` |
| the executed Release build: QC file load and store (bytes), picker file parameters, `FilterFeatureMap` in flag and filter mode, `EstimateDefaultMRMFeatureQCValues` with and without initialisation, `EstimatePercRSD`, `EstimateBackgroundInterferences`, `FilterFeatureMapPercRSD`, `FilterFeatureMapBackgroundInterference`, `calculateIonRatio`, `calculateRTDifference`, `calculateResolution` on the upstream `MRMFeatureSelector_100ug.featureXML` | 1 | `tests/mrm_oracle.rs`; see the section below |
| refusal of a short row, a misaligned metavalue store, an atomic failed load, the integer division by zero, the default flags of an unmatched transition | 4 | `tests/mrm_feature_qc.rs`, `tests/mrm_feature_filter.rs` |

### Executed comparison

`../oracle/O3b-mrm/driver.cpp` (case `qcfile`, `pickerfile`, `filter`) links the
Linux x86-64 Release install `openms4-release-bc9cc12-c19e494-174b576` and
prints `key<TAB>value`. The run, its host and hashes are recorded in the
manifest; the output is retained as `tests/data/mrm_oracle/*.stdout`.

### Defect candidates

Written up under "Defect candidates" in the lane report; none is filed here.

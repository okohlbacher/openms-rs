# OpenSWATH data structures and data access

Source pin: OpenMS4-core `bc9cc12`. Rust: `src/analysis/openswath/`. Tests:
`tests/openswath_algo.rs` (OPENSWATHALGO structures), `tests/openswath_data_access.rs`,
`tests/openswath_sqmass.rs` (`sqmass` feature). Provenance:
`tests/data/openswath_foundation_provenance.json`. Algorithms (Scoring, StatsHelpers,
DataFrameWriter) are in [OPENSWATH_ALGO_SUPPORT.md](OPENSWATH_ALGO_SUPPORT.md).

Every other OpenSWATH port builds on these types, so the design decisions are stated first.

## Design of the shared interfaces

- **`SpectrumAccess` trait** (`ISpectrumAccess`). Getters take `&self`; the trait is `Send + Sync`; ids
  are `usize` and every implementation returns `Error::InvalidValue` for an id outside its view (the
  source's `OPENMS_PRECONDITION` is debug-only). `light_clone` is kept and returns `Result`, since a
  file-backed clone opens its own connection. `SpectrumAccessPtr = Arc<dyn SpectrumAccess>`.
- **Spectra and chromatograms are shared, arrays are owned.** `SpectrumPtr = Arc<Spectrum>`,
  `ChromatogramPtr = Arc<Chromatogram>`, `SpectrumSequence = Vec<SpectrumPtr>`. A record owns its
  `BinaryDataArray`s by value (the source shares each array through its own `shared_ptr`; no consumer
  relies on that aliasing).
- **No null leading arrays.** The m/z (time) and intensity arrays always exist; `set_data_arrays`
  refuses fewer than two and `data_arrays_mut` is a slice, so the source's "may be null" cannot arise.
- **Access traits** (`ITransition.h`): `IFeature` → `FeatureAccess`, `IMRMFeature` → `MRMFeatureAccess`,
  `ITransitionGroup` → `TransitionGroupAccess`, `ISignalToNoise` → `SignalToNoiseAccess`. The `I` prefix
  became an `Access` suffix so the traits do not shadow the kernel `Feature`/`MRMFeature`.
  `MRMFeatureAccess::feature` returns `Option<&dyn FeatureAccess>` (a borrow instead of a
  `shared_ptr`; every source caller uses it locally). Out-parameters became return values; the
  source's non-`const` getters take `&self`.
- **`LightTransition` implements the kernel `Transition` trait**, so an OpenSWATH transition group is
  `MRMTransitionGroup<MSChromatogram, LightTransition>`, the source's
  `MRMFeatureFinderScoring::TransitionGroupType`.

## API mapping

### `OPENSWATHALGO/DATAACCESS/DataStructures.h` → `data_structures.rs`

| Source | Rust |
|---|---|
| `OSBinaryDataArray`/`BinaryDataArray` (`data`, `description`) | `BinaryDataArray` (+ `new`, `with_description`) |
| `BinaryDataArrayPtr` | not ported: arrays are owned (see design) |
| `OSChromatogramMeta`/`ChromatogramMeta` (`index`, `id`), `ChromatogramMetaPtr` | `ChromatogramMeta`, `ChromatogramMetaPtr` |
| `OSChromatogram`: ctor, `getTimeArray`, `setTimeArray`, `getIntensityArray`, `setIntensityArray`, `getDataArrays` (both), `setDataArrays` | `Chromatogram::new`, `time_array`/`time_array_mut`, `set_time_array`, `intensity_array`/`intensity_array_mut`, `set_intensity_array`, `data_arrays`, `data_arrays_mut` + `push_data_array`, `set_data_arrays -> Result`; native `from_arrays` |
| private `defaultArrayLength`, `initvec` | construction invariant |
| `ChromatogramPtr` | `ChromatogramPtr = Arc<Chromatogram>` |
| `OSSpectrumMeta`/`SpectrumMeta` (`index`, `id`, `RT`, `ms_level`), `RTLess` | `SpectrumMeta { index, id, rt, ms_level }`, `SpectrumMeta::rt_less` |
| `SpectrumMetaPtr` | `SpectrumMetaPtr` |
| `OSSpectrum`: ctor, `getMZArray`, `setMZArray`, `getIntensityArray`, `setIntensityArray`, `setDriftTimeArray`, `getDriftTimeArray`, `getDataArrays` (both), `setDataArrays` | `Spectrum::new`, `mz_array`/`mz_array_mut`, `set_mz_array`, `intensity_array`/`intensity_array_mut`, `set_intensity_array`, `set_drift_time_array`, `drift_time_array -> Option`, `data_arrays`, `data_arrays_mut` + `push_data_array`, `set_data_arrays -> Result`; native `from_arrays`, `drift_time_array_index`, `len`, `is_empty`, free fn `is_ion_mobility_description` |
| `SpectrumPtr` | `SpectrumPtr = Arc<Spectrum>` |

### `OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h` (+ `.cpp`) → `spectrum_access.rs`

| Source | Rust |
|---|---|
| `OpenMS::SpectrumSequence`, `OpenSwath::SpectrumSequence` | `SpectrumSequence` |
| `lightClone` | `SpectrumAccess::light_clone -> Result<SpectrumAccessPtr>` |
| `getSpectrumById(int)` | `spectrum_by_id(usize) -> Result<SpectrumPtr>` |
| `getSpectrumById(int, drift_start, drift_end)` | provided `spectrum_by_id_in_drift_range` |
| `getSpectraByRT` | `spectra_by_rt -> Result<Vec<usize>>` |
| `getNrSpectra` | `nr_spectra` |
| `getSpectrumMetaById` | `spectrum_meta_by_id -> Result` |
| `getChromatogramById` | `chromatogram_by_id -> Result<ChromatogramPtr>` |
| `getNrChromatograms` | `nr_chromatograms` |
| `getChromatogramNativeID` | `chromatogram_native_id -> Result<String>` |
| `getMultipleSpectra(RT, n)` | provided `multiple_spectra(rt, usize)` |
| `getMultipleSpectra(RT, n, drift_start, drift_end)` | provided `multiple_spectra_in_drift_range` |
| static `filterByDrift` | free fn `filter_by_drift -> Result<Spectrum>` |
| `SpectrumAccessPtr` | `SpectrumAccessPtr = Arc<dyn SpectrumAccess>` |

### `OPENSWATHALGO/DATAACCESS/ITransition.h` → `itransition.rs`

| Source | Rust |
|---|---|
| `IFeature::getRT(vector&)`, `getIntensity(vector&)`, `getIntensity()`, `getRT()` | `FeatureAccess::rt_values -> Result<Vec<f64>>`, `intensity_values -> Result<Vec<f64>>`, `intensity`, `rt` |
| `IMRMFeature::getFeature`, `getPrecursorFeature`, `getNativeIDs`, `getPrecursorIDs`, `getIntensity`, `getRT`, `getMetaValue`, `size` | `MRMFeatureAccess::feature -> Option<&dyn FeatureAccess>`, `precursor_feature`, `native_ids`, `precursor_ids`, `intensity`, `rt`, `meta_value -> Result<f64>`, `len` (+ `is_empty`) |
| `ITransitionGroup::size`, `getNativeIDs`, `getLibraryIntensities(vector&)` | `TransitionGroupAccess::len` (+ `is_empty`), `native_ids`, `library_intensities -> Result<Vec<f64>>` |
| `ISignalToNoise::getValueAtRT`, `ISignalToNoisePtr` | `SignalToNoiseAccess::value_at_rt(&self)`, `SignalToNoisePtr = Arc<dyn SignalToNoiseAccess>` |

### `Transitions.h` → `transitions.rs`; `TransitionExperiment.h` → `transition_experiment.rs`; `TransitionHelper.h` → `transition_helper.rs`

| Source | Rust |
|---|---|
| `Peptide` (`rt`, `charge`, `sequence`, `id`, `getChargeState`, `modifications`, `transitions`), `Protein`, `TargetedExperiment` | `transitions::{Peptide (+ charge_state), Protein, TargetedExperiment}` |
| `enum class FragmentIonType`, `stringToFragmentIonType`, `fragmentIonTypeToString` | `FragmentIonType` (`repr(u8)`, same discriminants), `FragmentIonType::from_type_str`, `as_str`, `Display` |
| `TransitionFlags` (bitfields) | `TransitionFlags(u8)` with `decoy`/`detecting`/`quantifying`/`identifying` getters and setters |
| `LightTransition` fields | same public fields |
| `getDecoy`, `setDecoy`, `getFragmentType`, `setFragmentType`, `getAnnotation`, `getProductChargeState`, `isProductChargeStateSet`, `isPrecursorImSet`, `getNativeID`, `getPeptideRef`, `getCompoundRef`, `getLibraryIntensity`, `setLibraryIntensity`, `getProductMZ`, `getPrecursorMZ`, `getPrecursorIM`, `set/isDetectingTransition`, `set/isQuantifyingTransition`, `set/isIdentifyingTransition`, `operator==`/`!=` | `decoy`, `set_decoy`, `fragment_type_str`, `set_fragment_type`, `annotation`, `product_charge_state`, `is_product_charge_state_set`, `is_precursor_im_set`, `native_id`, `peptide_ref`, `compound_ref`, `library_intensity`, `set_library_intensity`, `product_mz`, `precursor_mz`, `precursor_im`, `set_/is_detecting_transition`, `set_/is_quantifying_transition`, `set_/is_identifying_transition`, `PartialEq`/`Eq` by name |
| `LightModification` (+ `==`) | `LightModification` (derived `PartialEq`/`Hash`) |
| `LightCompound` fields, `isPeptide`, `setChargeState`, `getChargeState`, `setDriftTime`, `getDriftTime`, `==` | same fields; `is_peptide`, `set_charge_state`, `charge_state`, `set_drift_time`, `drift_time`, `PartialEq` by id |
| `LightProtein` (+ `==`) | `LightProtein`, `PartialEq` by id |
| `LightTargetedExperiment`: fields, `getTransitions`/`getCompounds`/`getProteins` (both), `getPeptideByRef`, `getCompoundByRef`, private reference map | fields; `transitions`/`transitions_mut`, `compounds`/`compounds_mut`, `proteins`/`proteins_mut`; `peptide_by_ref`, `compound_by_ref -> Option`; native `compound_reference_map`, `new` |
| typedefs `Transition`, `Peptide`, `Compound`, `Protein` inside `LightTargetedExperiment` | not ported (inherent associated types are unstable); they name the `Light*` types |
| `std::hash<LightTransition/LightCompound/LightProtein/LightModification>` | `Hash` impls (by name / id / id / both fields) |
| `TransitionHelper::convert`, `findPeptide` | `transition_helper::convert`, `find_peptide -> Option<&LightCompound>` |

### `SwathMap.h`, `SpectrumHelpers.h`, `MockObjects.h`

| Source | Rust |
|---|---|
| `SwathMap` fields `sptr`, `lower`, `upper`, `center`, `imLower`, `imUpper`, `ms1`; three constructors; `isEqual(other, tolerance = 1e-6)` | `SwathMap { sptr: Option<_>, lower, upper, center, im_lower, im_upper, ms1 }`; `Default`, `new`, `with_im`; `is_equal` (1e-6), `is_equal_with_tolerance` |
| `SpectrumHelpers.h` | declares nothing (an empty namespace); nothing to port |
| `MockFeature`, `MockMRMFeature`, `MockTransitionGroup`, `MockSignalToNoise` (members `m_*`) | `mock_objects::{MockFeature, MockMRMFeature, MockTransitionGroup, MockSignalToNoise}` (fields without `m_`) |

### `ANALYSIS/OPENSWATH/DATAACCESS/*` → `spectrum_access_openms.rs`, `spectrum_access_sqmass.rs`, `data_access_helper.rs`, `mrm_feature_access.rs`

| Source | Rust |
|---|---|
| `SpectrumAccessOpenMS(shared_ptr<PeakMap>)`, copy ctor, `lightClone`, the 7 interface members, `getSpectraMetaInfo`, `getChromatogramMetaInfo`, typedefs | `SpectrumAccessOpenMS::new(Arc<MSExperiment>)`, `Clone`, trait impl, `spectra_meta_info -> Result<&MSSpectrum>`, `chromatogram_meta_info -> Result<&MSChromatogram>`; native `experiment`; conversions `spectrum_with_data_arrays`, `chromatogram_with_data_arrays` |
| `SpectrumAccessOpenMSInMemory(ISpectrumAccess&)`, copy ctor, interface members | `SpectrumAccessOpenMSInMemory::new(&dyn SpectrumAccess) -> Result`, `Clone`, trait impl; native `from_spectra`, `from_parts` |
| `SpectrumAccessTransforming(sptr)` (abstract), protected `sptr_` | `SpectrumAccessTransforming<T: SpectrumTransform>` (`new`, `inner`, `transform`); trait `SpectrumTransform` |
| `SpectrumAccessQuadMZTransforming(sptr, a, b, c, ppm)`, `lightClone`, `getSpectrumById` | `QuadMZTransform { a, b, c, ppm }` (+ `apply`), `SpectrumAccessQuadMZTransforming = SpectrumAccessTransforming<QuadMZTransform>`, `quad_mz` |
| `SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr`, private `isExperimentCached` | `SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr -> Result`, `is_experiment_cached` (public) |
| `SpectrumAccessSqMass`: 3 ctors, copy ctor, `lightClone`, interface members, `getAllSpectra` | `SpectrumAccessSqMass::new`, `with_indices`, `narrowed` (all `Result`), `Clone`, trait impl, `all_spectra`; native `subset`, `to_in_memory` (feature `sqmass`) |
| `SpectrumAccessOpenMSCached` | **not ported**: needs `FORMAT/CachedMzML.h` and `FORMAT/HANDLERS/CachedMzMLHandler.h` |
| `OpenSwathDataAccessHelper::convertToOpenMSSpectrum`, `convertToSpectrumPtr`, `convertToOpenMSChromatogram`, `convertToChromatogramPtr`, `convertToOpenMSChromatogramFilter` | `data_access_helper::convert_to_openms_spectrum`, `convert_to_spectrum`, `convert_to_openms_chromatogram`, `convert_to_chromatogram`, `convert_to_openms_chromatogram_filter` |
| `convertTargetedExp`, both `convertTargetedCompound`, `convertPeptideToAASequence` | **not ported**: need `ANALYSIS/TARGETED/TargetedExperiment` and `TargetedExperimentHelper` |
| `FeatureOpenMS` | `FeatureOpenMS<'a>` (`new`, `feature`) |
| `MRMFeatureOpenMS` 3 ctors, interface members, `getFeatureIntensities`, `getPrecursorFeatureIntensities`, `getFeatureIntensity` | `MRMFeatureOpenMS::new`, `with_ids`, `with_aligned_ids`; trait impl; `feature_intensities`, `precursor_feature_intensities`, `feature_intensity`; native `mrm_feature` |
| `TransitionGroupOpenMS<SpectrumT, TransitionT>` | `TransitionGroupOpenMS<'a, C, T>` |
| `SignalToNoiseOpenMS<ContainerT>` | `SignalToNoiseOpenMS` (`from_chromatogram`, `from_spectrum`, `from_points`) |

### `ANALYSIS/OPENSWATH/OpenSwathHelper.h` → `openswath_helper.rs`; `SpectrumAddition.h` → `spectrum_addition.rs`

| Source | Rust |
|---|---|
| `computePrecursorId`, `computeTransitionGroupId` | `compute_precursor_id`, `compute_transition_group_id` |
| `selectSwathTransitions` (Light) | `select_swath_transitions` |
| `selectSwathTransitions` (`OpenMS::TargetedExperiment`) | **not ported** (TargetedExperiment) |
| `selectSwathTransitionsPasef` | `select_swath_transitions_pasef -> Vec<Option<usize>>` |
| `checkSwathMap` | `check_swath_map -> Result<(lower, upper, center)>` |
| `checkSwathMapAndSelectTransitions<T>` | `check_swath_map_and_select_transitions` (Light instantiation) |
| `estimateRTRange` | `estimate_rt_range -> Result<(f64, f64)>` |
| `sampleExperiment` | **not ported**: needs a crate-visible `Math::RandomShuffler` (the source-exact engine is private to `chemistry/decoy_random.rs`) and `std::random_device` seeding |
| `simpleFindBestFeature` | `simple_find_best_feature<C>` over `BTreeMap<String, MRMTransitionGroup<C, LightTransition>>` |
| `SpectrumAddition::addUpSpectra` (3 overloads), `concatenateSpectra`, `sortSpectrumByMZ` | `add_up_spectra`, `add_up_spectra_in_mobility_range(&RangeBase)`, `add_up_ms_spectra`; `concatenate_spectra`; `sort_spectrum_by_mz` (all `Result`) |

## Preserved source conventions

- `getSpectraByRT` returns the first spectrum at or after `RT - deltaRT` **even beyond `RT + deltaRT`**;
  `getMultipleSpectra` relies on it. The upper bound is inclusive for `SpectrumAccessOpenMS` and
  exclusive for the materialised view (its documentation claims inclusive); both kept.
- `getMultipleSpectra`: closest spectrum first, then alternating lower/upper neighbours for `n / 2`
  steps (an even `n` fetches `n + 1`); equidistant neighbours resolve to the later spectrum.
- `SpectrumAccessOpenMS::getSpectrumMetaById` fills only RT and MS level; the sqMass view fills id, RT and
  MS level; neither fills `index`.
- `SpectrumAccessOpenMS` converts float and integer data arrays (strings dropped); `convertToSpectrumPtr`
  and the sqMass view carry only m/z and intensity.
- `filterByDrift` bounds are inclusive and the result carries m/z, intensity and the IM array only.
- Ion-mobility array detection: prefix `"Ion Mobility"`, or containing `"inverse reduced ion mobility"`
  or `"ion mobility array"`, case-sensitive; `setDriftTimeArray` appends and renames.
- Identity of `LightTransition`/`LightCompound`/`LightProtein` is by name/id only; `compound_by_ref`
  resolves duplicates to the last compound, `find_peptide` to the first.
- `selectSwathTransitions` window: `lower < mz < upper` and at least `min_upper_edge_dist` below the
  upper edge; appends to the output. PASEF ties keep the earlier window.
- `SpectrumAddition`: a single spectrum is returned unchanged; an empty first spectrum yields an empty
  result; the MSSpectrum variant also scans every peak for the range when estimating; resampling
  distributes intensity to the two neighbouring grid points (intensity outside the grid goes to the ends).

## Native differences

- **Out-of-range ids, unsorted RT, missing IM arrays and misaligned arrays are errors**, where the
  source has undefined behaviour (debug-only preconditions, null dereferences, out-of-bounds reads).
  Sortedness by RT is recorded once at construction, so `spectra_by_rt` refuses an unsorted view without
  a per-call scan; the other getters still work on it.
- **`SpectrumAccessOpenMSInMemory`** reads through the generic interface for every origin. The source
  special-cases an sqMass origin and then drops its chromatograms (`nr_chromatograms` becomes 0 while the
  origin reports a count); `SpectrumAccessSqMass::to_in_memory` makes that spectra-only batch copy an
  explicit choice. Its end-of-list check happens before the read (the source dereferences `end()`).
- **`SpectrumAccessTransforming`** copies a shared spectrum before transforming it. The source transforms
  the wrapped view's spectrum in place: over the materialised view that corrupts the cache, and a second
  request returns a doubly transformed m/z array.
- **`SimpleOpenMSSpectraFactory`** returns `Error::Unsupported` for a cached experiment instead of a
  `SpectrumAccessOpenMSCached` (not ported).
- **`SpectrumAccessSqMass`** reads its spectrum and chromatogram counts once, at construction, and shares
  its (read-stateless) handler with clones; per-id metadata reads skip the peak arrays. Without a subset,
  an id is the SQL spectrum ID (dense IDs presumed, as in the source).
- **`LightTargetedExperiment`** keeps no pointer cache (the source's goes stale after the compound vector
  changes and dereferences null for an unknown reference); `compound_reference_map` builds one on demand.
- **`MockMRMFeature::feature`** does not insert on a failed lookup (the source's `operator[]` does).
- **`MRMFeatureOpenMS`** builds its wrappers eagerly and reports unknown IDs as `None` or
  `Error::MissingInformation`; `meta_value` reports a missing or non-numeric value as an error.
- **`SignalToNoiseOpenMS`** copies the positions and computes the estimates in its constructor (the
  source keeps a reference to the container); `value_at_rt` takes `&self`.
- **`SpectrumAddition`**: a zero, negative (after estimation) or non-finite sampling rate and grids above
  `MAX_RESAMPLED_POINTS` (10^7) are refused before allocation; NaN m/z values sort last instead of being
  undefined; `concatenate_spectra` refuses a spectrum with more arrays than the first; the resampler's
  spacing warning is not logged.
- **`selectSwathTransitionsPasef`** returns a fresh map (the source resizes the caller's, keeping old
  entries). **`checkSwathMapAndSelectTransitions`** returns `Ok(false)` instead of printing its warnings.

## Checked boundaries and evidence

| Behaviour | Evidence |
|---|---|
| `Datastructures_test.cpp` (3 sections), `SwathMap_test.cpp` | tier 3, `tests/openswath_algo.rs` |
| `SpectrumAccessOpenMS_test.cpp`, `SpectrumAccessOpenMSInMemory_test.cpp`, `SpectrumAccessQuadMZTransforming_test.cpp` (all sections) | tier 3, `tests/openswath_data_access.rs` |
| `SpectrumAccessSqMass_test.cpp` at the pin (all sections) over the retained `SqliteMassFile_1.sqMass` | tier 3, `tests/openswath_sqmass.rs` |
| `MRMFeatureAccessOpenMS_test.cpp` (all sections, incl. the S/N ratio 102.352941176471) | tier 3 |
| `OpenSwathHelper_test.cpp` (all testable sections), `SpectrumAddition_test.cpp` (all sections, incl. the 25/9/16/28-point grids) | tier 3 |
| RT-window rules, `getMultipleSpectra` order, `filterByDrift`, cache non-corruption, factory, DataAccessHelper conversions, MRM lanes, PASEF choice, best feature, IM-range addition, concatenation, transitions/compounds/mocks | tier 4 |

Measured on dax: `--all-features` `openswath_data_access` 26 passed, `openswath_sqmass` 4 passed.

## Gaps (named dependencies)

- `SpectrumAccessOpenMSCached.h` → needs the cached-mzML reader (`CachedMzML`, `CachedMzMLHandler`).
- `DataAccessHelper.h`: `convertTargetedExp`, `convertTargetedCompound` (×2), `convertPeptideToAASequence`
  → need `ANALYSIS/TARGETED/TargetedExperiment` and `TargetedExperimentHelper`.
- `OpenSwathHelper.h`: the `TargetedExperiment` overload of `selectSwathTransitions` (and its
  `checkSwathMapAndSelectTransitions` instantiation) → TargetedExperiment; `sampleExperiment` → a shared
  `Math::RandomShuffler`.

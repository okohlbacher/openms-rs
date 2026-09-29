# MRM transition group picking

`src/analysis/openswath/mrm_transition_group_picker.rs` ports
`ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h` and `.cpp` at OpenMS4-core
`bc9cc12`. The class finds peak groups ("features") in the chromatograms of one
precursor: every detecting chromatogram is smoothed and picked with
`PeakPickerChromatogram` (`src/processing/chromatogram.rs`), the largest (or
widest) remaining picked peak defines the boundaries of a group, every
chromatogram is resampled onto the sampling of the seed's raw chromatogram and
integrated with `PeakIntegrator` (`src/analysis/peak_integrator.rs`) between
them, and the loop repeats until no picked peak is left.

```rust
use openms::analysis::openswath::mrm_transition_group_picker::MRMTransitionGroupPicker;
use openms::concept::unique_id::UniqueIdGenerator;
use openms::kernel::mrm::{MRMTransitionGroup, SimpleTransition};
use openms::kernel::{ChromatogramPeak, MSChromatogram};

let mut group = MRMTransitionGroup::<MSChromatogram, SimpleTransition>::new();
group.set_transition_group_id("PEPTIDE/2");
let peaks: Vec<_> = (0..30)
    .map(|i| {
        let x = f64::from(i);
        ChromatogramPeak::new(x, (1000.0 * (-0.5 * ((x - 15.0) / 2.0).powi(2)).exp() + 3.0) as f32)
    })
    .collect();
let chromatogram = MSChromatogram { native_id: "t1".into(), peaks, ..Default::default() };
group.add_chromatogram(chromatogram, "t1")?;
group.add_transition(SimpleTransition::new("t1", -1.0), "t1")?;

let mut param = MRMTransitionGroupPicker::defaults()?;
param.set_value("PeakPickerChromatogram:gauss_width", openms::param::ParamValue::Float(5.0), "", &[])?;
let picker = MRMTransitionGroupPicker::from_param(&param)?;
picker.pick_transition_group(&mut group, &mut UniqueIdGenerator::from_seed(1))?;
assert_eq!(group.features().len(), 1);
# Ok::<(), openms::Error>(())
```

## API mapping

| Source member | Rust |
| --- | --- |
| `MRMTransitionGroupPicker()`, `~MRMTransitionGroupPicker()` | `MRMTransitionGroupPicker::default()`, `Drop` |
| `DefaultParamHandler` base: `getDefaults`, `getParameters`, `setParameters` | `defaults()`, `to_param()`, `from_param()` / `from_param_with_warnings()` |
| `updateMembers_()` (protected) | the private `from_complete_param`, run by `from_param` |
| `operator=` (protected; copies nothing) | `Clone` (copies every field) |
| `pickTransitionGroup(MRMTransitionGroup<SpectrumT, TransitionT>&)` | `pick_transition_group(&mut MRMTransitionGroup<MSChromatogram, T: Transition>, &mut UniqueIdGenerator)` |
| `createMRMFeature(group, picked, smoothed, chr_idx, peak_idx)` | `create_mrm_feature(...)` |
| `createMRMFeature(group, picked, smoothed, xics, detecting_xic, transition_chroms, picked_inputs, chr_idx, peak_idx)` | `create_mrm_feature_with_inputs(..., &TransitionGroupInputs, ...)`; `TransitionGroupInputs::new` builds the precomputed arguments |
| `pickApex(...)` | `pick_apex(...)` |
| `pickFragmentChromatograms(...)` | `pick_fragment_chromatograms(...)`; the four `double&` totals are `&mut PeakGroupTotals` |
| `pickPrecursorChromatograms(...)` | `pick_precursor_chromatograms(...)` |
| `remove_overlapping_features(picked, best_left, best_right)` | `remove_overlapping_features(...)` (associated function) |
| `findLargestPeak(picked, chr_idx&, peak_idx&)` | `find_largest_peak(picked) -> Option<(usize, usize)>` |
| `findLargestPeakCached_(picked, cache, chr_idx&, peak_idx&)` | `find_largest_peak_cached(picked, &mut [Option<usize>])` |
| `findWidestPeakIndices(picked, chrom_idx&, point_idx&)` | `find_widest_peak_indices(picked) -> Result<Option<(usize, usize)>>` |
| `selectChromHelper_` (protected) | private `select_chrom_helper` |
| `computeQuality_` (protected) | private `compute_quality` |
| `recalculatePeakBorders_` (protected) | private `recalculate_peak_borders` |
| `prepareMasterContainer_`, `resampleChromatogram_` (protected) | private `prepare_master_container`, `resample_chromatogram` (with the source-exact `raster`) |
| members `peak_integration_`, `background_subtraction_`, `boundary_selection_method_` | typed fields `peak_integration: PeakIntegrationSignal`, `background_subtraction: BackgroundSubtraction`, `boundary_selection_method: BoundarySelection` |
| members `recalculate_peaks_`, `use_precursors_`, `use_consensus_`, `compute_peak_quality_`, `compute_peak_shape_metrics_`, `compute_total_mi_`, `min_qual_`, `stop_after_feature_`, `stop_after_intensity_ratio_`, `min_peak_width_`, `recalculate_peaks_max_z_`, `resample_boundary_` | public fields of the same names (`min_qual_` is `minimal_quality`, the parameter's name) |
| members `picker_`, `pi_` | `picker: PeakPickerChromatogram`, `integrator: PeakIntegrator` |
| namespace `MRMTransitionGroupPickerMeta` (29 registered names) | module `meta` (29 string constants) |
| `PickerPoolShared`, `g_picker_pool_shared()` | none: thread-local allocation pools that only avoid reallocation |
| `lin_interpolate` (defined in the `.cpp`, declared nowhere, unused) | none |

Every parameter (`stop_after_feature` -1, `stop_after_intensity_ratio` 0.0001,
`min_peak_width` 0.001, `peak_integration` original, `background_subtraction`
none, `recalculate_peaks` false, `use_precursors` false, `use_consensus` true,
`recalculate_peaks_max_z` 1.0, `minimal_quality` -10000, `resample_boundary`
15, `compute_peak_quality` false, `compute_peak_shape_metrics` false,
`compute_total_mi` false, `boundary_selection_method` largest) has its source
type, description, valid strings and `advanced` tag, and the
`PeakPickerChromatogram:` and `PeakIntegrator:` subsections are those classes'
own parameter sets.

## Preserved source conventions

- The picker's arithmetic is the Linux x86-64 Release build's, bit for bit on
  the 32 oracle cases (below). Resampling is the source's
  `LinearResamplerAlign::raster(iterators)`: samples before the grid are added
  to its first point in `f32`, each interior sample is split
  `dist_right / (dist_left + dist_right)` and `dist_left / (...)` in `f64` and
  stored as `f32`, and the rest goes to the last point.
  `crate::processing::resample_to_grid` computes the split as `1 - fraction`,
  which rounds differently, so the picker keeps its own copy.
- Picked chromatograms are sorted by ascending intensity (stable, permuting
  the data arrays) before the loop, as `sortByIntensity` does; the `largest`
  seed search walks a per-chromatogram cursor down from the end.
- Sub-features are keyed by chromatogram native ID with the source's
  duplicate-key behaviour (`DuplicateKeyPolicy::SourceOverwrite`).
- `width_at_50` is stored even without `compute_peak_shape_metrics`, and
  identifying (non-detecting) transitions are integrated but do not add to
  `total_xic`.
- Peak groups whose `[leftWidth, rightWidth]` lies inside an earlier group's
  are dropped, as are zero-intensity groups.

Source defects that change results are reproduced, because every OpenSWATH
result downstream depends on them; each is documented at its item:

1. `stop_after_feature` compares `cnt > stop_after_feature`, so one more peak
   group is collected than the number given (oracle case `g3_stop1`: two).
2. `pickFragmentChromatograms` indexes the smoothed chromatograms by
   transition, although they hold only the detecting chromatograms in
   chromatogram order.
3. `pickPrecursorChromatograms` looks the precursor's smoothed chromatogram and
   edges up at "number of fragment chromatograms + k", so smoothed integration
   of a group with a precursor chromatogram throws unless every fragment is
   detecting and `use_precursors` is set (oracle cases `g1_smoothed`,
   `g3_smoothed`).
4. With `compute_total_mi`, a detecting transition's mutual information,
   already averaged over the detecting transitions, is divided by their number
   again before it is added to the peak group's `total_mi`.
5. `computeQuality_` averages the co-elution scores in the accumulator that
   already holds the shape scores, and its outlier index counts only
   chromatograms with signal but indexes all picked chromatograms.
6. `pickApex` bounds its integration of the intensity-sorted picked
   chromatograms with `PosBegin`/`PosEnd` binary searches; the port reproduces
   libstdc++'s `lower_bound`/`upper_bound` on that unsorted data.

A seventh, `computeQuality_` storing `max_int = intensity > max_int` (a
comparison) where it means an intensity, only affects border lists the source
never reads; the port computes the peak counts alone.

## Native differences

- Errors instead of undefined behaviour: an index past a picked chromatogram's
  arrays, a seed whose raw chromatogram is empty, an `original` background with
  no raw sample between the boundaries, a smoothed chromatogram that is not
  there, and a crossed binary-search range in `pick_apex` all return
  `Error::InvalidValue` where the source reads out of bounds.
- A failed pick leaves the group unchanged; the source may already have added
  nothing but has computed everything in place.
- Unique IDs come from a caller-owned `UniqueIdGenerator` instead of the
  process-wide one; `from_seed` gives a reproducible stream.
- A non-finite metadata value (a NaN peak quality when no chromatogram has
  signal around the peak) is an error: native `MetaValue` holds finite floats
  only.
- `pick_apex` with EMG fitting returns `Error::Unsupported`: the source would
  fit the EMG model to the centroided picked peaks, which this port does not
  reproduce.
- `max_chromatograms` (100,000) bounds the group before any work; every other
  bound is the picker's and the integrator's (`max_points`, `max_work`).
- The integrator refuses a shape-metric apex with no sample at or above it
  inside the bounds, where the source returns metrics from its sampled search.
- `PeakPickerChromatogram` refuses an apex without two neighbouring samples,
  where the source's boundary walk reads before the first sample.

## Checked boundaries and evidence

- **Tier 1**: `group_picker_matches_the_release_build_bit_for_bit` replays 32
  cases over five groups (the two class-test groups and three synthetic ones:
  a five-transition group with an identifying and a non-quantifying
  transition, a precursor chromatogram and three peaks; a group with two
  precursor isotopes; a precursor-only group) against the records the Release
  build wrote on ibminode06 (`../oracle/O3a-mrm-picker`, run twice, identical).
  The cases cover legacy and corrected picking, Savitzky-Golay smoothing,
  overlap removal, consensus and non-consensus boundaries, `largest` and
  `widest` seeds, recalculated borders, peak quality with and without a
  minimal quality, shape metrics, mutual information, `original` and `exact`
  background with all three integration types and a vertical-division
  baseline, smoothed integration, precursors (used and not), a minimal width
  that drops everything, a stop count and a stop ratio, and three cases where
  the source throws. Every one of the 3,358 records (feature and sub-feature
  positions, intensities, qualities, every metadata key and value, every hull
  point) agrees bit for bit, and the Rust pick fails exactly where the source
  throws.
- **Tier 3**: every section of `MRMTransitionGroupPicker_test.cpp` is
  transcribed with its literals and `TEST_REAL_SIMILAR` semantics:
  `pickTransitionGroup` on both groups (consensus, non-consensus, only
  quantifying, only detecting), `createMRMFeature` with and without
  background subtraction and Simpson integration, `findLargestPeak`,
  `findWidestPeakIndices` and `remove_overlapping_features`. The class test
  spells four keys `PeakPickerChromatogram::method`, which match no parameter;
  the transcription uses the defaults they leave in place.
- **Tier 4**: unit tests of the source-exact raster, the libstdc++ bound
  emulation and the padded range.

Source files at `bc9cc12`:

| Source file | SHA-256 |
| --- | --- |
| `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h` | `a4ad5a43d42eaa665bf0f3e69d75a036750b266d9bda31e3fe20b231a3ff86ba` |
| `src/openms/source/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.cpp` | see `tests/data/openswath_peak_picking_provenance.json` |
| `src/tests/class_tests/openms/source/MRMTransitionGroupPicker_test.cpp` | see `tests/data/openswath_peak_picking_provenance.json` |

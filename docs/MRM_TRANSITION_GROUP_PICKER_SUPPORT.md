# MRM transition group picking, mobilogram picking and mass trace correlation

This document covers three OpenSWATH headers at OpenMS4-core `bc9cc12`:
`MRMTransitionGroupPicker.h` (below), `PeakPickerMobilogram.h` (section
"Mobilogram peak picking") and `MasstraceCorrelator.h` (section "Mass trace
correlation").

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

## Mobilogram peak picking

`src/analysis/openswath/peak_picker_mobilogram.rs` ports
`ANALYSIS/OPENSWATH/PeakPickerMobilogram.h` and `.cpp`. The source class is
`PeakPickerChromatogram` on ion mobility: its smoothing
(`GaussFilter::filter(Mobilogram&)`, `SavitzkyGolayFilter::filter(Mobilogram&)`),
seed picking (`PeakPickerHiRes::pick(const Mobilogram&, ...)`), boundary walk,
overlap split and raw sum are the same template code as for a chromatogram, with
the S/N gate left out. The port runs `PeakPickerChromatogram` on the
mobilogram's `(mobility, intensity)` samples with `signal_to_noise` 0 and
`report_sn` off, and converts the result back.

| Source member | Rust |
| --- | --- |
| `PeakPickerMobilogram()`, `~PeakPickerMobilogram()` | `Default`, `Drop` |
| `getDefaults` / `getParameters` / `setParameters` + `updateMembers_` | `defaults()`, `to_param()`, `from_param()` / `from_param_with_warnings()` |
| enum `FLOATINDICES` (six) | `PeakPickerMobilogram::IDX_*`, `SIZE_OF_FLOATINDICES` |
| struct `PeakPositions` | `PeakPositions` |
| `pickMobilogram(mobilogram, picked)`, `pickMobilogram(mobilogram, picked, smoothed)` | `pick_mobilogram(&Mobilogram) -> MobilogramPickingResult` |
| public members `integrated_intensities_`, `left_width_`, `right_width_` | `MobilogramPickingResult::{integrated_intensities, left_widths, right_widths}` |
| `filterTopPeak(picked, vector&, peak_pos&)`, `filterTopPeak(picked, mobilogram&, peak_pos&)` | `filter_top_peak(&picked, &mut [Mobilogram])`, `filter_top_peak_single(&picked, &mut Mobilogram)` |
| protected `pickMobilogram_`, `integratePeaks_`, `findClosestPeak_`, `removeOverlappingPeaks_` | `PeakPickerChromatogram`'s private equivalents |
| protected `findHighestPeak_`, `filterPeakIntensities_` (2), `extractFloatValues_`, `extractIntValues_`, `filterTopPeak_` (2) | private `highest_peak`, `filter_peak_intensities`, `float_to_size` |
| members `sgolay_*`, `gauss_width_`, `use_gauss_`, `peak_width_`, `signal_to_noise_`, `sn_win_len_`, `sn_bin_count_`, `write_sn_log_messages_`, `method_`, `remove_overlapping_` | the public fields of the same meaning |
| members `pp_`, `sgolay_`, `gauss_` | built per call |

Parameters: `sgolay_frame_length` 9, `sgolay_polynomial_order` 3, `gauss_width`
0.002, `use_gauss` `false`, `peak_width` -1, `signal_to_noise` 1 (at least 0),
`sn_win_len` 1 (an integer parameter), `sn_bin_count` 4,
`write_sn_log_messages` `false`, `remove_overlapping_peaks` `false`, `method`
`corrected` (valid `legacy`, `corrected`, `crawdad`). `signal_to_noise`,
`sn_win_len`, `sn_bin_count` and `write_sn_log_messages` are stored and, as in
the source, change nothing: there is no boundary S/N gate, and the seed picker
took the constructor's S/N of one.

Source defects, and what the port does:

- `updateMembers_` never assigns `remove_overlapping_`, so `pickMobilogram`
  reads an uninitialised `bool`. The port assigns it from
  `remove_overlapping_peaks`, as that parameter's description promises.
- `method` `crawdad` is accepted and then handled in neither branch, so the
  per-peak vectors of an earlier call (or none) are read. The port refuses it
  in `from_param`.
- `findHighestPeak_` returns the picked peak's index as `apex`, although it
  documents all three positions as mobilogram indices; the port keeps the
  source value and says so on `PeakPositions`.

Evidence: **tier 1**, `mobilogram_and_hull_scores_match_the_release_build_bit_for_bit`
replays eight cases (the class-test mobilogram and two synthetic ones;
corrected, legacy, a forced width, Gaussian and Savitzky-Golay smoothing)
against the Release build on ibminode06 (`../oracle/O3a-mobilogram-masstrace`,
run twice, identical): picked and smoothed mobilograms, all six arrays, the
public per-peak vectors and both `filterTopPeak` overloads, bit for bit. The
driver constructs the picker in zeroed storage, so the never-assigned
`remove_overlapping_` reads `false`, which is the port's default. Tier 3:
`PeakPickerMobilogram_test.cpp:67-84` (apex position and intensity,
integrated intensity and both border mobilities of the 60-point class-test
mobilogram, through `getParameters`/`setParameters`); tier 4 for the
defaults, the refusals and the empty-input branches of the top-peak filters.

## Mass trace correlation

`src/analysis/openswath/masstrace_correlator.rs` ports
`ANALYSIS/OPENSWATH/MasstraceCorrelator.h` and `.cpp`, the algorithm behind
`ClusterMassTraces` and `ClusterMassTracesByPrecursor`.

| Source member | Rust |
| --- | --- |
| `MasstraceCorrelator()`, `~MasstraceCorrelator()` | `Default`, `Drop` |
| `DefaultParamHandler` base (`getDefaults`, `getParameters`, `setParameters`) | `defaults()`, `to_param()`, `from_param()` / `from_param_with_warnings()` |
| `ProgressLogger` base | the `_with_progress` methods (`CACHE_PROGRESS_LABEL`, `CORRELATION_PROGRESS_LABEL`) |
| typedef `MasstracePointsType` | `MasstracePoints` |
| `createPseudoSpectra(map, pseudo_spectra&, min_peak_nr, min_correlation, max_lag, max_rt_apex_difference)` | `create_pseudo_spectra(...) -> Vec<PseudoSpectrum>` and `create_pseudo_spectra_with_progress` |
| `scoreHullpoints(h1, h2, lag&, lag_intensity&, pearson_score&, min_corr, max_lag, mindiff = 0.1)` | `score_hullpoints(h1, h2, min_corr, max_lag, mindiff) -> HullpointScores` |
| `createConsensusMapCache(map, feature_points&, max_intensities&, rt_cache&)` | `create_consensus_map_cache(map) -> ConsensusMapCache` and `create_consensus_map_cache_with_progress` |
| protected `matchMassTraces_(h1, h2, vec1&, vec2&, mindiff, padEnds = true)` | private `match_mass_traces` (both `pad_ends` modes, unit-tested) |

The three parameters (`sgolay_frame_length` 15, `sgolay_polynomial_order` 3,
`gauss_width` the integer 50) are declared under the handler name
`MRMFeatureFinderScoring`, a copy-paste slip the port keeps because warnings
carry it, and nothing reads them (the class has no `updateMembers_`).

Native differences:

- The source stores each joined trace's `RT_apex`, `RT_diff`, `lag`,
  `pearson_score` and `lag_intensity` in float arrays on the pseudo-spectrum,
  one entry per joined trace, so the arrays are one shorter than the peaks. A
  native spectrum's arrays must match its peaks; the port returns them as
  `PseudoSpectrum::correlations`, in `f32` as stored.
- When the Pearson score is not above `min_corr`, the source leaves its `lag`
  and `lag_intensity` out-parameters as they were; the port returns `None`.
- A non-finite RT is an error, where the source's matching loop never advances;
  a seed feature without handles is an error, where the source reads the first
  point of an empty trace; `max_work` bounds the all-vs-all loop.
- Equal RTs within a trace keep their handle order (a stable sort); the
  source's `std::sort` may order them differently for traces longer than 16
  points.

Evidence: **tier 1**, `pseudo_spectra_match_the_release_build_bit_for_bit`
(five `createPseudoSpectra` cases on the class-test map and a synthetic
14-feature map, varying `min_peak_nr`, the correlation, lag and RT limits)
and six `scoreHullpoints` pairs in
`mobilogram_and_hull_scores_match_the_release_build_bit_for_bit`: every
pseudo-spectrum peak, every score array and every score bit-identical to the
Release build (`../oracle/O3a-mobilogram-masstrace`). Its stderr carries the
source's non-fatal `Double conversion error of "4the_second_seed202.0"`,
which confirms the 0.0 the test substitutes for that RT. Tier 3: every
section of `MasstraceCorrelator_test.cpp`: the matching cases (`38-209`, on
the private function), the cross-correlation and the three `scoreHullpoints`
calls (`211-283`), and `createPseudoSpectra` with `min_peak_nr` 0 and 1 on
the retained `Masstraces_Testdata.consensusXML` (`285-315`).

Source files at `bc9cc12`:

| Source file | SHA-256 |
| --- | --- |
| `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h` | `a4ad5a43d42eaa665bf0f3e69d75a036750b266d9bda31e3fe20b231a3ff86ba` |
| `src/openms/source/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.cpp` | `8f88e5dd65b9916c49ba8a7bb4977503441c4a065d401a389f263d6f02303c1b` |
| `src/tests/class_tests/openms/source/MRMTransitionGroupPicker_test.cpp` | `bb65f0664bfa589d4cd3a7f33e23d3b475434d3d991ea1a09259365d9d0782f2` |
| `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/PeakPickerMobilogram.h` | `e9116b9c00eddf76efd53a75ab7bf01c99729e3992928e7567cf8f89e945a5f7` |
| `src/openms/source/ANALYSIS/OPENSWATH/PeakPickerMobilogram.cpp` | `e4525f6f0e2a9bfdbb60c824529afac9f5948319edd5439311eb266786318d0d` |
| `src/tests/class_tests/openms/source/PeakPickerMobilogram_test.cpp` | `794cae20b1350166a5d9c9a873791863f64b9d034f16e9f32a3587485de849f8` |
| `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MasstraceCorrelator.h` | `5014897ad693b08963b556fcb93971277ce5d8026d19a09f5fad38c7f325a67e` |
| `src/openms/source/ANALYSIS/OPENSWATH/MasstraceCorrelator.cpp` | `3cec77ff3750f3b4151dc042e4e61996264f9309e15dca6f29d463a6fc047a7c` |
| `src/tests/class_tests/openms/source/MasstraceCorrelator_test.cpp` | `3a17dafa044a17d27a6ae702292130f5b32fa1fc047052b3097c8bdf9ffe3b5e` |

Every other source the port reproduces, the fixtures and the oracle artifacts
are hashed in `tests/data/openswath_peak_picking_provenance.json`.

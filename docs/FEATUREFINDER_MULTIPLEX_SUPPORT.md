# FeatureFinderMultiplex support (`FEATUREFINDER/Multiplex*`, `FeatureFinderMultiplexAlgorithm`)

Source revision: OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`;
TOPP test data `0cb15f2` of the test-data package.

| Source header | Rust | Status |
|---|---|---|
| `FEATUREFINDER/MultiplexDeltaMasses.h` | `src/analysis/featurefinder/multiplex/delta_masses.rs` | complete |
| `FEATUREFINDER/MultiplexDeltaMassesGenerator.h` | `src/analysis/featurefinder/multiplex/delta_masses_generator.rs` | complete |
| `FEATUREFINDER/MultiplexIsotopicPeakPattern.h` | `src/analysis/featurefinder/multiplex/isotopic_peak_pattern.rs` | complete |
| `FEATUREFINDER/MultiplexSatelliteCentroided.h` | `src/analysis/featurefinder/multiplex/satellite.rs` | complete |
| `FEATUREFINDER/MultiplexSatelliteProfile.h` | `src/analysis/featurefinder/multiplex/satellite.rs` | complete |
| `FEATUREFINDER/MultiplexFilteredPeak.h` | `src/analysis/featurefinder/multiplex/filtered_peak.rs` | complete |
| `FEATUREFINDER/MultiplexFilteredMSExperiment.h` | `src/analysis/featurefinder/multiplex/filtered_ms_experiment.rs` | complete |
| `FEATUREFINDER/MultiplexFiltering.h` | `src/analysis/featurefinder/multiplex/filtering.rs` | complete |
| `FEATUREFINDER/MultiplexFilteringCentroided.h` | `src/analysis/featurefinder/multiplex/filtering_centroided.rs` | complete |
| `FEATUREFINDER/MultiplexFilteringProfile.h` | `src/analysis/featurefinder/multiplex/filtering_profile.rs` (private filters in `filtering.rs`) | complete |
| `FEATUREFINDER/MultiplexClustering.h` | `src/analysis/featurefinder/multiplex/clustering.rs` | complete |
| `FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h` | `src/analysis/featurefinder/multiplex/algorithm.rs` | complete |
| `FEATUREFINDER/PeakWidthEstimator.h` (dependency) | `src/analysis/featurefinder/peak_width_estimator.rs` | complete |
| `PROCESSING/MISC/SplinePackage.h` (dependency) | `src/processing/spline/spline_package.rs` | complete |
| `PROCESSING/MISC/SplineInterpolatedPeaks.h` (dependency) | `src/processing/spline/spline_interpolated_peaks.rs` | complete |

`ML/REGRESSION/LinearRegressionWithoutIntercept.h` is used by the algorithm
and is **not** ported here: its one computation (`sum(x*y) / sum(x*x)`, NaN for
fewer than two points) is a private helper, `slope_without_intercept`, in
`algorithm.rs`. When the regression header is ported, that helper should be
replaced by it.

## Placement

Everything is `FEATUREFINDER` and needs `kernel`, `chemistry`, `processing`
(peak picking, splines), `math` (statistics, `source_sort`, grid clustering)
and `format` (the debug consensusXML and the ostream number formatting), so
it lives in `src/analysis/featurefinder/multiplex/`, one file per header except
that the two satellite headers share `satellite.rs`. `PeakWidthEstimator` is a
`FEATUREFINDER` header and sits beside it in `src/analysis/featurefinder/`. The
two `PROCESSING/MISC` spline headers depend only on `CubicSpline2d` and the
kernel, and join the other spline ports in `src/processing/spline/`.

## API mapping

### MultiplexDeltaMasses

| Source | Rust |
|---|---|
| `LabelSet` (`std::multiset<std::string>`) | `LabelSet`: `new`, `single`, `insert`, `iter`, `as_slice`, `len`, `is_empty`, `count`, `FromIterator` |
| `DeltaMass` struct, `delta_mass`, `label_set` | `DeltaMass { delta_mass, label_set }` |
| `DeltaMass(double, LabelSet)` | `DeltaMass::new` |
| `DeltaMass(double, const std::string&)` | `DeltaMass::with_label` |
| `MultiplexDeltaMasses()` | `new()`, `Default` |
| `MultiplexDeltaMasses(const std::vector<DeltaMass>&)` | `from_delta_masses` |
| `getDeltaMasses()` (both overloads) | public field `delta_masses`, `delta_masses()` |
| `labelSetToString` (static) | `label_set_to_string` (associated) |
| free `operator<` | `source_less`, `source_cmp` (not `PartialOrd`; see the item) |

### MultiplexDeltaMassesGenerator

| Source | Rust |
|---|---|
| `Label` struct and constructor | `Label`, `Label::new` |
| `MultiplexDeltaMassesGenerator()` (the `labels` param handler) | `LabelMassShifts` with `defaults`, `to_param`, `from_param`, `from_param_with_warnings` |
| `MultiplexDeltaMassesGenerator(labels, max_nr_labelled_aas, label_mass_shift)` | `new(&str, i32, BTreeMap<String, f64>) -> Result` |
| `generateKnockoutDeltaMasses` | `generate_knockout_delta_masses -> Result` |
| `printSamplesLabelsList(std::ostream&)` | `print_samples_labels_list(&mut impl Write)`, `samples_labels_text` |
| `printDeltaMassesList(std::ostream&)` | `print_delta_masses_list(&mut impl Write)`, `delta_masses_text` |
| `getDeltaMassesList` (both overloads) | `delta_masses_list() -> &[MultiplexDeltaMasses]` |
| `getSamplesLabelsList` (both overloads) | `samples_labels_list() -> &[Vec<String>]` |
| `getLabelShort`, `getLabelLong` | `label_short`, `label_long -> Option<&str>` |
| `extractLabelSet(const AASequence&)` | `extract_label_set`, and `extract_label_set_from_str` |
| private `labels_`, `labels_list_`, `max_nr_labelled_aas_`, `label_delta_mass_` | `labels()`, `labels_list()`, `max_nr_labelled_aas()`, `label_delta_mass()` |
| private `fillLabelMasterList_`, `label_master_list_` | `label_master_list()` |
| private `label_short_long_`, `label_long_short_` | private maps behind `label_long`/`label_short` |
| (native) | `Labelling` and `labelling()`, `KNOWN_LABELS`, `MAX_DELTA_MASS_PATTERNS` |

### MultiplexIsotopicPeakPattern

| Source | Rust |
|---|---|
| `MultiplexIsotopicPeakPattern(int c, int ppp, MultiplexDeltaMasses ms, int msi)` | `new(...) -> Result` |
| `getCharge`, `getPeaksPerPeptide`, `getMassShifts`, `getMassShiftIndex` | `charge`, `peaks_per_peptide`, `mass_shifts`, `mass_shift_index` |
| `getMassShiftCount`, `getMZShiftCount` | `mass_shift_count`, `mz_shift_count` |
| `getMassShiftAt`, `getMZShiftAt` | `mass_shift_at`, `mz_shift_at -> Option<f64>` |
| (native) | `mz_shifts()`, `MAX_MZ_SHIFTS` |

### MultiplexSatelliteCentroided, MultiplexSatelliteProfile

| Source | Rust |
|---|---|
| `MultiplexSatelliteCentroided(size_t rt_idx, size_t mz_idx)`, `getRTidx`, `getMZidx` | `new`, public fields `rt_idx`, `mz_idx` |
| `MultiplexSatelliteProfile(float rt, double mz, float intensity)`, `getRT`, `getMZ`, `getIntensity` | `new`, public fields `rt`, `mz`, `intensity` |

### MultiplexFilteredPeak

| Source | Rust |
|---|---|
| constructor | `new(mz, rt, mz_idx, rt_idx)` |
| `getMZ`, `getRT`, `getMZidx`, `getRTidx` | `mz`, `rt`, `mz_idx`, `rt_idx` |
| `addSatellite(size_t, size_t, size_t)`, `addSatellite(const MultiplexSatelliteCentroided&, size_t)` | `add_satellite`, `add_satellite_from` |
| `addSatelliteProfile(float, double, float, size_t)`, `addSatelliteProfile(const MultiplexSatelliteProfile&, size_t)` | `add_satellite_profile`, `add_satellite_profile_from` |
| `checkSatellite` | `check_satellite` |
| `getSatellites`, `getSatellitesProfile` (`std::multimap`) | `satellites`, `satellites_profile` (`SatelliteMultimap<T> = BTreeMap<usize, Vec<T>>`), and `satellite_entries`, `satellite_profile_entries` |
| `size`, `sizeProfile` | `satellite_count`, `satellite_profile_count` |

### MultiplexFilteredMSExperiment

| Source | Rust |
|---|---|
| constructor | `new`, `Default` |
| `addPeak`, `getPeak` | `add_peak`, `peak -> Option<&_>`; `peaks()` |
| `getMZ(size_t)`, `getMZ()`, `getRT(size_t)`, `getRT()` | `mz_at`, `mz`, `rt_at`, `rt` |
| `size` | `len`, `is_empty` |
| `writeDebugOutput(exp_picked, debug_out)` | `write_debug_output(exp_picked, path, &mut UniqueIdGenerator)` (feature `consensusxml`), built on `debug_consensus_map` |

### MultiplexFiltering

| Source | Rust |
|---|---|
| `White2Original` | `White2Original = Vec<Vec<usize>>` |
| constructor | `new(...) -> Result` (`usize` isotope counts) |
| `getCentroidedExperiment` | `centroided_experiment` |
| `getBlacklist` | `blacklist_experiment` |
| protected `updateWhiteMSExperiment_` | `update_white_ms_experiment` |
| protected `checkForSignificantPeak_` | `check_for_significant_peak -> Option<usize>` (white spectrum by index) |
| protected `filterPeakPositions_` | `filter_peak_positions(mz, band: Range<usize>, pattern, &mut peak) -> Result<bool>` |
| protected `blacklistPeak_(peak, pattern_idx)` | `blacklist_peak -> Result` |
| protected `blacklistPeak_(peak)` | not ported: declared, never defined in the source |
| protected `filterAveragineModel_`, `filterPeptideCorrelation_` | `filter_averagine_model`, `filter_peptide_correlation -> Result<bool>` |
| protected members `exp_centroided_`, `blacklist_`, `exp_centroided_white_`, `exp_centroided_mapping_`, `patterns_`, the parameters | `pub(crate)` fields; readers `blacklist`, `white_experiment`, `white_mapping`, `patterns`, `intensity_cutoff` |
| (native) | `AveragineType`, `MAX_FILTERED_PEAKS` |

### MultiplexFilteringCentroided, MultiplexFilteringProfile

| Source | Rust |
|---|---|
| constructors | `new(...) -> Result`; the base is reachable through `Deref` |
| `filter` | `filter`, `filter_with_progress` |
| `MultiplexFilteringProfile::getPeakBoundaries` | `peak_boundaries` |
| private `filterAveragineModel_`, `filterPeptideCorrelation_` (profile) | `MultiplexFiltering::filter_averagine_model_profile`, `filter_peptide_correlation_profile` |
| private `exp_spline_profile_`, `boundaries_` | `spline_profile`, `peak_boundaries` |

### MultiplexClustering

| Source | Rust |
|---|---|
| `Point` typedef | `grid_based_cluster::Point` |
| constructor (profile, picked, boundaries, rt_typical) | `from_profile -> Result` |
| constructor (exp, mz_tolerance, unit, rt_typical) | `from_centroided -> Result` |
| `cluster` | `cluster`, `cluster_with_progress` |
| `MultiplexDistance(double)`, `MultiplexDistance()`, `operator()` | `MultiplexDistance::new`, `Default`, `distance` |
| private grid, scaling and typical RT | `grid_spacing_mz`, `grid_spacing_rt`, `rt_scaling`, `rt_typical` |
| (native) | `MAX_GRID_LINES` |

### FeatureFinderMultiplexAlgorithm

| Source | Rust |
|---|---|
| constructor, `DefaultParamHandler` | `new`, `Default`, `with_parameters`, `defaults`, `to_param`, `from_param`, `from_param_with_warnings`; typed `FeatureFinderMultiplexParameters` |
| `run(MSExperiment&, bool)` | `run(&mut MSExperiment, bool)`, `run_with_progress` |
| `getFeatureMap`, `getConsensusMap`, `getBlacklist` | `feature_map`, `consensus_map`, `blacklist`, and `feature_map_mut`, `consensus_map_mut`, `blacklist_mut` for the non-const references the tool edits in place |
| protected `generatePeakPatterns_` | `generate_peak_patterns` (method and free function) |
| protected `correctPeptideIntensities_` | `correct_peptide_intensities` |
| protected `determinePeptideIntensitiesCentroided_`, `...Profile_` | `determine_peptide_intensities_centroided`, `determine_peptide_intensities_profile` |
| protected `generateMapsCentroided_`, `generateMapsProfile_` | `generate_maps_centroided`, `generate_maps_profile` |
| protected members | `profile_experiment`, `centroided_experiment`, `is_centroided`, `reports_progress`, `charge_range`, `isotopes_per_peptide_range`, `label_mass_shift` |
| `ProgressLogger` base, `prog_log_` | `ProgressReporter` arguments |
| process-wide `UniqueIdGenerator` | public field `unique_ids` |

### PeakWidthEstimator, SplinePackage, SplineInterpolatedPeaks

| Source | Rust |
|---|---|
| `PeakWidthEstimator(exp_picked, boundaries)`, `getPeakWidth` | `new -> Result`, `peak_width -> Result<f64>` |
| `SplinePackage(pos, intensity)`, `getPosMin`, `getPosMax`, `getPosStepWidth`, `isInPackage`, `eval` | `new -> Result`, `pos_min`, `pos_max`, `pos_step_width`, `is_in_package`, `eval -> Result<f64>` |
| `SplineInterpolatedPeaks(pos, intensity)`, `(const MSSpectrum&)`, `(const MSChromatogram&)` | `new`, `from_spectrum`, `from_chromatogram -> Result` |
| `getPosMin`, `getPosMax`, `size`, `getNavigator(scaling = 0.7)` | `pos_min`, `pos_max`, `size`, `navigator(scaling) -> Result<Navigator<'_>>`, `DEFAULT_NAVIGATOR_SCALING` |
| `Navigator::eval`, `Navigator::getNextPos` | `Navigator::eval`, `Navigator::next_pos -> Result<f64>` |
| `Navigator()` (pyOpenMS) | not ported: a null package pointer that every call dereferences |

## Preserved source conventions

- **Sort orders are the Release build's.** Every `std::sort` whose ties are
  observable goes through `crate::math::source_sort::source_sort_by` (the
  libstdc++ introsort): the mass shift list, the peak patterns
  (`lessPattern`), the filtering experiment's spectra, the feature map, and the
  median m/z of the clustering. Stable sorts (peaks by m/z, chromatograms by RT,
  the consensus map) use Rust's stable sort.
- **The filters run serially.** The source's `omp parallel for` over the peaks
  reads the blacklist while a `critical` section writes it, so its result
  depends on scheduling; the one-thread result is the reference, and the
  upstream TOPP tests run with `-threads 1`.
- **Averagine intensities are binary32**, as the source's `IsotopeDistribution`
  stores them (`ProbabilityPrecision::SourceSingle`).
- **White versus original indices.** The zeroth-peak test of
  `filterPeakPositions_` compares a *white* peak index with satellites that
  store *original* indices; reproduced.
- **Parsing.** Labels split exactly as `boost::split` does, keeping empty
  tokens; the known-label test is a substring test (`Arg` passes); `std::stod`
  reads a numeric prefix (`3.5abc` is 3.5). `Label:13C(6)` maps back to `Arg6`,
  never `Lys6`.
- **Map building.** A cluster's centroided satellites stop merging at the
  first duplicate of a peak (`break`, not `continue`); the RT span test of a
  peptide uses the source's minimum/maximum bookkeeping that resets on every
  satellite of the lightest peptide's first isotope; column sizes are counted
  before a multiplet can be dropped; the maps are not cleared between runs;
  `run` swaps the input into the algorithm's buffer.
- **Knock-outs** keep the source's missing quadruplet triplet `(0,1,3)`.

## Native differences

Each is documented at its item.

- Undefined behaviour of the source becomes an error: reading past a pattern's
  m/z shifts or the averagine distribution, a Leu labelling with one sample,
  a knock-out pattern shorter than the first, satellite indices outside the
  data, empty inputs to `front()`, a peptide without a mono-isotopic trace.
- Infinite loops of the source become errors: a non-advancing profile scan,
  a zero grid step in the clustering, `getNextPos(NaN)`.
- Wrapping conversions become errors: negative `max_nr_labelled_aas`, negative
  isotope counts; a zero charge is refused where the source divides by zero.
- `SplinePackage` (through `CubicSpline2d`) refuses non-increasing positions,
  where the source builds a NaN spline. In the algorithm's centroided
  quantification such a chromatogram gets no spline: the source's NaN spline
  evaluates to zero everywhere and contributes no ratio pair, so the result is
  the same.
- `SplineInterpolatedPeaks` with fewer than two points left after the removal
  of surplus zeros builds no package (the source reads past its work vector).
- `std::stod` hexadecimal floats are not recognised in numeric labels.
- `getLabelShort`/`getLabelLong` return `None` for an unknown label where the
  source inserts and returns an empty string.
- Resource ceilings: `MAX_MZ_SHIFTS`, `MAX_DELTA_MASS_PATTERNS`,
  `MAX_FILTERED_PEAKS`, `MAX_GRID_LINES`, `MAX_RATIO_SCAN_STEPS` and the
  spline `MAX_POINTS`, each checked before the work it bounds.
- Unique IDs come from a caller-owned generator.

## Checked boundaries and evidence

- **Tier 1** (retained C++ output): `tests/featurefinder_multiplex_topp.rs`
  runs the algorithm on the inputs and `.ini` files of the eleven upstream
  `TOPP_FeatureFinderMultiplex_*` tests (`topp/CMakeLists.txt:463-526`) and
  compares feature and consensus maps with the retained outputs under the
  upstream FuzzyDiff contract (ratio 1.01 or absolute 0.01; counts, charges,
  map indices and labels exact). The cases cover centroided and profile input,
  SILAC, dimethyl doublets and triplets with knock-outs, unlabelled RNA and DNA
  averagines, numeric labels and a label-free run.
- **Tier 3** (class tests): `tests/featurefinder_multiplex.rs` (data
  structures, delta masses and their generator), `tests/featurefinder_multiplex_filtering.rs`
  (base, centroided and profile filtering, clustering, `PeakWidthEstimator`),
  `tests/featurefinder_multiplex_algorithm.rs` (the algorithm class test) and
  `tests/spline_interpolated_peaks.rs` (both spline headers): every section of
  the fourteen class tests.
- **Tier 4**: the pattern formula, the multiset order, the documented SILAC
  mass shift list, the parameter tree, error paths.

No oracle driver was run for this group: the retained TOPP outputs cover the
algorithm end to end, and every class-test literal is asserted. An exact
(bitwise) differential against the Release build remains open.

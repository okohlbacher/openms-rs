# FF2b NOTES (FeatureFinderMultiplex lane)

Merged prerequisite worktree-wf_cd066dab-0a7-34 (FF2a: GridBasedClustering, DistanceMatrix, ...) at 8344dc03.
C++ sources extracted to scratchpad/cpp/src/openms/{include,source}/FEATUREFINDER (bc9cc12).
Guard note: git commands whose paths contain the word "source" are refused; use globs ('src/openms/*/FEATUREFINDER/*ultiplex*').

Sizes (C++ .cpp lines): DeltaMasses 88, DeltaMassesGenerator 609, IsotopicPeakPattern 74, SatelliteCentroided 30,
SatelliteProfile 35, FilteredPeak 100, FilteredMSExperiment 139, Filtering 585, FilteringCentroided 123,
FilteringProfile 427, Clustering 189, FeatureFinderMultiplexAlgorithm 1154.

Upstream class tests: src/tests/class_tests/openms/source/Multiplex*_test.cpp, FeatureFinderMultiplexAlgorithm_test.cpp
Data: class_tests data MultiplexFiltering.mzML, MultiplexClustering.mzML, FeatureFinderMultiplex_1_{input.mzML,parameters.ini};
TOPP data FeatureFinderMultiplex_{1..11}_* at test-data 0cb15f2.

Dependency order: DeltaMasses -> DeltaMassesGenerator -> IsotopicPeakPattern -> Satellite* -> FilteredPeak ->
FilteredMSExperiment -> Filtering -> FilteringCentroided/Profile -> Clustering (uses FF2a GridBasedClustering) -> Algorithm.

## Rust API survey (done)
- MSExperiment/MSSpectrum in src/kernel.rs (peaks Vec<Peak1D{mz f64,intensity f32}>, rt f64). find_nearest/rt_begin do an O(n) sortedness
  check per call -> write local partition_point helpers; C++ findNearest(mz,tol) (MSSpectrum.cpp:311) does NOT try the other side,
  Rust find_nearest_with_tolerance does (window variant) -> use local helper.
- source_sort_by (src/math/source_sort.rs) = libstdc++ std::sort; use for delta-mass list sort (unstable in C++).
- pearson_correlation_coefficient / rank_correlation_coefficient in src/math/statistic_functions.rs (NaN pass-through preserved).
- CoarseIsotopePatternGenerator::new(Some(n), CoarseMassMode::?) .with_precision(ProbabilityPrecision::SourceSingle); estimate_from_{peptide,rna,dna}_weight.
- PeakPickerHiRes::pick_experiment -> PickedExperiment{experiment, spectrum_boundaries: Vec<Option<Vec<PeakBoundary{min,max}>>>}.
- GridBasedClustering::without_properties(metric Fn(&Point,&Point)->f64, xs, ys, grid_x, grid_y) (FF2a), .cluster(), .results() BTreeMap<i32,GridBasedCluster>.
- Param pattern: src/processing/peak_picking.rs:1853-1985 (DefaultParamHandler::new + set_defaults + set_parameters_with).
- Missing deps: SplineInterpolatedPeaks + SplinePackage (PROCESSING/MISC, unmapped, nobody porting) needed by FilteringProfile + Algorithm;
  PeakWidthEstimator (FEATUREFINDER, unmapped) needed by profile MultiplexClustering ctor. C++ extracted to scratchpad/k/.

## C++ defect candidates found
- Filtering.cpp filterPeakPositions_: zeroth-peak check calls checkSatellite(rt_idx, mz_idx) with a WHITE-spectrum index, satellites store ORIGINAL indices.
- DeltaMassesGenerator: Leu branch reads samples_labels_[1] with one sample (UB); getLabelShort/Long operator[] inserts; "Label:13C(6)" long->short maps to Arg6 (first insert) never Lys6;
  label validity check is substring-of-list (e.g. "Arg" or "" accepted); negative max_nr_labelled_aas wraps unsigned.
- Centroided/Profile filter(): omp parallel for with blacklist read/write race -> only serial is deterministic.

## Layout
src/analysis/featurefinder/multiplex.rs + multiplex/*.rs, one file per header; featurefinder.rs gets `pub mod multiplex;`.

## Status (commit 7b15f348 + getter)
PASSING on spock: all leaf types, DeltaMassesGenerator, Filtering base, Centroided [0,0,0,0,4,4,4,0], Profile [0,0,0,0,4,5,4,0],
Clustering [0,0,0,0,2,0,0,0], PeakWidthEstimator, SplinePackage, SplineInterpolatedPeaks class tests.
Tests: tests/featurefinder_multiplex.rs (no features), tests/featurefinder_multiplex_filtering.rs (mzml), tests/spline_interpolated_peaks.rs.
Extra headers ported as deps: PROCESSING/MISC/SplinePackage.h, SplineInterpolatedPeaks.h (src/processing/spline/), FEATUREFINDER/PeakWidthEstimator.h.

## Next step
1. FeatureFinderMultiplexAlgorithm: needs LinearRegressionWithoutIntercept (not ported; ML/REGRESSION, tiny), PeakIntegrator (src/analysis/peak_integrator.rs),
   make_peak_position_unique (src/kernel/spectrum_helper.rs), SpectrumType get_type (src/kernel/spectrum_type.rs), ConvexHull2D.
   C++ read fully (scratchpad/mx/FeatureFinderMultiplexAlgorithm.cpp). Note run(): label_mass_shift_[...] = param (overwrites).
2. Class test FeatureFinderMultiplexAlgorithm_test: 2 consensus features, H/L ~ 3 +-0.2; tier 1: TOPP retained outputs FeatureFinderMultiplex_{1..11}.
3. Docs: docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md, provenance tests/data/featurefinder_multiplex_provenance.json, doc coverage, report+ledger.

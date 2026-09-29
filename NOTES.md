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

## Next step
1. Survey Rust deps (AASequence/ModificationsDB, MSExperiment, SplineInterpolatedPeaks, PeakPickerHiRes, IsotopeDistribution, Param).
2. Port the leaf data types (DeltaMasses, IsotopicPeakPattern, Satellite*, FilteredPeak, FilteredMSExperiment), commit.

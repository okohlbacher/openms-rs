# FeatureHypothesis support

`analysis::feature_hypothesis` implements the complete public `FeatureHypothesis` operation group in OpenMS4-core `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. The same header's `CmpMassTraceByMZ` and `CmpHypothesesByScore` become `mass_trace_mz_less` and `hypothesis_score_greater`; its two-field `Range` value becomes `MetaboIsotopeMassWindow`. This does not complete `FeatureFindingMetabo`: candidate assembly, model prediction, scoring and feature-map publication remain the next group.

`FeatureHypothesis<'a>` owns an ordered vector of immutable borrowed `MassTrace` references. `add_mass_trace` retains duplicate references. Rust lifetimes prevent dropping or changing a trace while the hypothesis uses it; there are no raw pointers or hidden trace copies. Clone and assignment preserve reference identity, score, charge and limits. `checked_clone` provides a bounded copy of the reference vector. Ordinary Rust Clone, Debug, borrowing, destruction and direct scalar access keep their standard costs.

The public operation mapping is:

| Source | Native |
|---|---|
| default/copy/assignment/destructor | Default/new, Clone/checked_clone, assignment, ownership/drop |
| addMassTrace/getSize | add_mass_trace, len/is_empty/traces |
| getLabel/getLabels | label/labels |
| getScore/setScore, getCharge/setCharge | score/set_score, charge/set_charge |
| getAllIntensities | all_intensities(smoothed) |
| getAllCentroidMZ/RT/IM | all_centroid_mz/rt/im |
| getIsotopeDistances | isotope_distances |
| getCentroidMZ/RT, getFWHM | centroid_mz/rt, fwhm |
| getMonoisotopicFeatureIntensity | monoisotopic_feature_intensity(smoothed) |
| getSummedFeatureIntensity/getMaxIntensity | summed_feature_intensity/max_intensity(smoothed) |
| getNumFeatPoints | number_of_feature_points |
| getConvexHulls/getChromatograms | convex_hulls/chromatograms(feature_id) |

Rust calls supply the source default `false` explicitly for optional smoothing arguments. First-trace queries use cached values exactly as stored, without re-sorting traces or recomputing centroids. Isotope distances are consecutive signed m/z differences without division by charge. The point total counts duplicate references repeatedly; it is cached safely because shared borrowing fixes trace lengths. IM vectors include cached zero from records without an IM flag. Label joining retains empty components and underscores. Scores retain all `f64` values, including infinities and NaNs; the source predicates use ordinary `<`/`>` and return false for unordered comparisons. No Eq/Ord, total floating order or source object-layout claim is made.

Intensity queries reuse each trace's quantification choice and established MassTrace behavior. Area uses cached FWHM borders, Median uses raw peak intensities even when `smoothed=true`, and MaxHeight selects raw or smoothed heights. Hypothesis maximum intensity instead compares each trace's apex and starts from zero. Cached FWHM and summed/max intensity are zero for an empty hypothesis; monoisotopic intensity and centroid queries return errors. The existing MassTrace finite-data, unset-smoothing and resource boundaries continue to apply.

Hull export produces one scan-envelope `ConvexHull2D` per reference using every raw RT/m/z point, ignoring intensity and cached centroids. Chromatogram export uses every raw RT and `f32` intensity, ignoring raw m/z, trace labels, smoothing and quantification. It sorts each chromatogram by RT, uses BasePeak type, sets both name and native ID to decimal `feature_id` plus underscore plus zero-based trace index, and gives every precursor the **first** trace's cached m/z. Precursor metadata `peptide_sequence` is the decimal feature ID. The common hypothesis charge is applied to every precursor. Product, acquisition settings, source file, processing history and arrays retain their new-record defaults; no unrelated metadata is copied.

The following native boundaries are deliberate:

- Equal-RT chromatogram points retain encounter order; source `std::sort` does not specify ties.
- Empty-hypothesis chromatogram export returns a checked error instead of source unchecked indexing. A nonempty hypothesis may contain empty traces and returns empty chromatograms for them.
- Charge uses platform-independent `i64` storage. Chromatogram export checks its conversion to the narrower precursor `i32`; other operations do not reject unused charge values.
- Consumed nonfinite coordinates/intensities and overflowing distance/sum arithmetic return errors. Unused raw fields remain uninspected: for example, a hull can use peaks with NaN intensity, and a chromatogram can use peaks with NaN raw m/z. Finite negative coordinates/intensities and signed zero remain represented.
- Checked mutations publish only after fallible work. Queries return owned results, so a failure cannot publish a partial vector or change any borrowed trace.

`FeatureHypothesisLimits` defaults to one million references, ten million total member peaks, 50 million weighted work visits/comparisons and 256 MiB new logical payload per checked call. Adding a reference checks aggregate membership including duplicates before reserving; geometric vector growth bounds repeated additions without cloning all existing references. Multitrace queries recheck current membership limits and share one cumulative work/byte allowance. Monoisotopic intensity only charges and checks the first trace it consumes; cached scalar/point-count access is O(1).

Allocating queries charge vector slots, copied strings, sorting scratch and sparse precursor metadata-map storage before construction. Existing MassTrace median/hull calls receive conservative cumulative precharges before invocation, so each member cannot reset the enclosing hypothesis allowance. Hull precharges include small-vector growth as well as input points, output scans and sorting scratch. The individual MassTrace limits also remain effective. Budgets describe checked work and conservative new logical storage, not allocator internals or a hard resident-memory reservation.

`tests/feature_hypothesis.rs` exercises every operation, duplicate borrowing/copy identity, all quantification modes, source seven-peak literals, independent raw-envelope/chromatogram expectations, exact first-trace/cache semantics, signed zero and full score storage, resource failures and unconsumed invalid fields. The compile-fail lifetime example prevents a hypothesis from outliving a trace. `tests/data/feature_hypothesis_source.tsv` is an exact projection of the MassTrace class-test input, including its original double intensity before the source `f32` narrowing. The FeatureFindingMetabo class test contains no separate active FeatureHypothesis assertions: helper expectations are source-expression-derived, with the numeric MassTrace reference anchors identified in the provenance manifest. No expected value is generated from this Rust implementation. This subgroup adds no dependency, model resource, registry or global state.

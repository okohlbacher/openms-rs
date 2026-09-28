# Feature finder models and fitters

Ports, at core SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4`, the
one-dimensional models and fitters of `FEATUREFINDER/` and the classes built on
them:

| Header | Rust |
| --- | --- |
| `ML/INTERPOLATION/LinearInterpolation.h` (header only) | [`src/math/ml/linear_interpolation.rs`](../src/math/ml/linear_interpolation.rs) |
| `FEATUREFINDER/BaseModel.h`, `BaseModel_impl.h`, `ModelDescription.h`, `InterpolationModel.h` | [`src/analysis/featurefinder/base_model.rs`](../src/analysis/featurefinder/base_model.rs) |
| `GaussModel.h`, `BiGaussModel.h`, `EmgModel.h` | [`src/analysis/featurefinder/gauss_models.rs`](../src/analysis/featurefinder/gauss_models.rs) |
| `IsotopeModel.h`, `ExtendedIsotopeModel.h` | [`src/analysis/featurefinder/isotope_model.rs`](../src/analysis/featurefinder/isotope_model.rs) |
| `Fitter1D.h`, `LevMarqFitter1D.h`, `MaxLikeliFitter1D.h` | [`src/analysis/featurefinder/fitter1d.rs`](../src/analysis/featurefinder/fitter1d.rs) |
| `GaussFitter1D.h`, `BiGaussFitter1D.h`, `EmgFitter1D.h`, `IsotopeFitter1D.h`, `ExtendedIsotopeFitter1D.h` | [`src/analysis/featurefinder/fitters.rs`](../src/analysis/featurefinder/fitters.rs) |
| `EmgScoring.h` | [`src/analysis/featurefinder/emg_scoring.rs`](../src/analysis/featurefinder/emg_scoring.rs) |
| `PeakWidthEstimator.h` | [`src/analysis/featurefinder/peak_width_estimator.rs`](../src/analysis/featurefinder/peak_width_estimator.rs) |
| `SeedListGenerator.h` | [`src/analysis/featurefinder/seed_list_generator.rs`](../src/analysis/featurefinder/seed_list_generator.rs) |
| `FeatureFinderAlgorithm.h`, `FeatureFinderDefs.h` (header only) | [`src/analysis/featurefinder/feature_finder_algorithm.rs`](../src/analysis/featurefinder/feature_finder_algorithm.rs) |

Tests: [`tests/featurefinder_models.rs`](../tests/featurefinder_models.rs) (class
tests, no feature), [`tests/featurefinder_io.rs`](../tests/featurefinder_io.rs)
(class tests reading fixtures; `mzml` + `consensusxml`) and
[`tests/featurefinder_oracle.rs`](../tests/featurefinder_oracle.rs) (executed
Release build, no feature). Fixture:
[`tests/data/featurefinder_models/oracle.tsv`](../tests/data/featurefinder_models/oracle.tsv).
Manifest: [`tests/data/featurefinder_models_provenance.json`](../tests/data/featurefinder_models_provenance.json).

Everything is serial, as the source is (none of these files has an OpenMP
pragma). No crate was added: the Lorentzian peak shape uses `statrs`' Cauchy
density (already a dependency; its `pdf` is Boost's `1 / (pi * scale * (1 +
xs * xs))` term for term), the Levenberg-Marquardt path is
`math::fitters::levenberg_marquardt` (Eigen's `maxfev` emulated), and the
B-spline is `processing::spline::b_spline::BSpline2d`.

## API mapping

### LinearInterpolation

`template <typename Key = double, typename Value = Key> class LinearInterpolation`
becomes the non-generic `LinearInterpolation` over `f64`: the library
instantiates only `<double, double>`.

| Source | Rust |
| --- | --- |
| `LinearInterpolation(KeyType scale = 1., KeyType offset = 0.)` | `new(scale, offset)`, `Default` (1, 0) |
| copy constructor, `operator=`, destructor | `Clone`, `PartialEq`, drop |
| `value(pos)`, `addValue(pos, value)`, `derivative(pos)` | `value`, `add_value`, `derivative` |
| `getData()` const / non-const, `setData`, `empty` | `data`, `data_mut`, `set_data`, `is_empty` |
| `key2index`, `index2key` | `key2index`, `index2key` |
| `getScale`, `setScale`, `getOffset`, `setOffset` | `scale`, `set_scale`, `offset`, `set_offset` |
| `setMapping(scale, inside, outside)` | `set_mapping` |
| `setMapping(inside_low, outside_low, inside_high, outside_high)` | `set_mapping_by_points` |
| `getInsideReferencePoint`, `getOutsideReferencePoint` | `inside_reference_point`, `outside_reference_point` |
| `supportMin`, `supportMax` | `support_min`, `support_max` |
| typedefs `value_type`, `key_type`, `container_type` and aliases | `f64`, `Vec<f64>` |

### BaseModel, BaseModel_impl, ModelDescription, InterpolationModel

| Source | Rust |
| --- | --- |
| `class BaseModel : public DefaultParamHandler` (abstract) | `trait BaseModel` |
| typedefs `IntensityType`, `CoordinateType` (double), `PositionType` (`DPosition<1>`), `PeakType` (`Peak1D`), `SamplesType` | `f64`, `f64`, `f64`, `Peak1D`, `Vec<Peak1D>` |
| `BaseModel()` (default `cutoff` 0), copy, `operator=`, destructor | implementors' `new`/`defaults`, `Clone` |
| `getIntensity(const PositionType&)` (pure) | `intensity(f64)` |
| `isContained(pos)` | `is_contained` (provided) |
| `fillIntensity(PeakType&)`, `fillIntensities(begin, end)` (templates) | `fill_intensity(&mut Peak1D)`, `fill_intensities(&mut [Peak1D])` |
| `getCutOff`, `setCutOff` | `cut_off`, `set_cut_off` |
| `getSamples(SamplesType&)` (pure), `getSamples(std::ostream&)` | `samples() -> Vec<Peak1D>`, `write_samples(&mut dyn Write)` |
| `getName`, `getParameters`, `setParameters`, `getDefaults` (inherited) | `name`, `to_param`, `set_parameters`; each model's associated `defaults()` |
| protected `cut_off_`, `updateMembers_` | `InterpolationCore::cut_off`, the private `set_parameters_impl` |
| `BaseModel_impl.h` (includes the five models) | `enum Model1D { Gauss, BiGauss, Emg, Isotope, ExtendedIsotope }`, implementing both traits |
| `template <UInt D> class ModelDescription` | `struct ModelDescription { name, param }` (the unused `D` dropped) |
| `ModelDescription()`, copy, `operator=`, `ModelDescription(const BaseModel*)` | `Default`, `Clone`, `from_model(&M)` |
| `getName`/`setName` (const and mutable), `getParam`/`setParam` (const and mutable) | public fields `name`, `param` |
| `operator==`, `operator!=` | `source_equal` (name and source `Param` equality) |
| (class documentation: "allows reconstruction of the model") | native `create_model() -> Result<Model1D>` |
| `class InterpolationModel : public BaseModel` (abstract) | `trait InterpolationModel: BaseModel` + `struct InterpolationCore` |
| `InterpolationModel()` (`interpolation_step` 0.1, `intensity_scaling` 1), copy, `operator=` | `InterpolationCore::default`, `Clone` |
| `getIntensity(pos)`, `getIntensity(CoordinateType)` | `BaseModel::intensity` (`interpolation().value(pos)`) |
| `getInterpolation` | `interpolation` |
| `getScalingFactor`, `setScalingFactor` | `scaling_factor`, `set_scaling_factor` |
| `setOffset` (virtual) | `set_offset` |
| `getSamples(SamplesType&)` | `InterpolationCore::samples` (via `BaseModel::samples`) |
| `getCenter` (virtual, throws `NotImplemented`) | `center() -> Result<f64>` |
| `setSamples` (virtual, throws `NotImplemented`) | `set_samples() -> Result<()>` |
| `setInterpolationStep` | `set_interpolation_step`; native `interpolation_step` getter |
| protected `interpolation_`, `interpolation_step_`, `scaling_` | `InterpolationCore` fields; `core`/`core_mut` |
| typedefs `KeyType`, `LinearInterpolation` | `f64`, `LinearInterpolation` |

### The five models

Every model has `NAME`, `new() -> Result<Self>` (the source default
constructor, which samples the defaults - except `IsotopeModel`, whose
`updateMembers_` does not sample), `defaults()`, `from_param`,
`from_param_with_warnings`, and implements `BaseModel` and
`InterpolationModel`. Copy constructor and `operator=` are `Clone`.

| Source | Rust |
| --- | --- |
| `GaussModel`: `setOffset`, `setSamples`, `getCenter` (the mean); protected `min_`, `max_`, `statistics_` | `GaussModel`: trait methods; `bounding_box()`, `statistics()` |
| `BiGaussModel`: `setOffset`, `setSamples`, `getCenter` (`statistics2_` mean); protected `min_`, `max_`, `statistics1_`, `statistics2_` | `BiGaussModel`: trait methods; `bounding_box()`, `statistics()` |
| `EmgModel`: `setOffset`, `setSamples`, `getCenter` (the `statistics:mean`); protected `min_`, `max_`, `statistics_`, `height_`, `width_`, `symmetry_`, `retention_`; typedef `ContainerType` | `EmgModel`: trait methods; `bounding_box()`, `emg_parameters()` |
| `IsotopeModel`: `enum Averagines {C, H, N, O, S, AVERAGINE_NUM}`, `getCharge`, `setOffset`, `getOffset`, `getFormula`, `setSamples(const EmpiricalFormula&)`, `using setSamples` (throws), `getCenter` (monoisotopic m/z), `getIsotopeDistribution`; protected `isotope_stdev_`, `isotope_lorentz_fwhm_`, `charge_`, `mean_`, `monoisotopic_mz_`, `averagine_[]`, `max_isotope_`, `trim_right_cutoff_`, `isotope_distance_`, `isotope_distribution_` | `IsotopeModel`: `charge`, `set_offset`, `offset`, `formula`, `set_samples_from_formula`, `set_samples` (`Unsupported`), `center`, `isotope_distribution`; `peak_shape()` and `IsotopePeakShape` for `isotope:mode:mode`; the averagine order is the private `AVERAGINE_KEYS` |
| `ExtendedIsotopeModel`: `enum Averagines`, `getCharge`, `setOffset`, `getOffset`, `setSamples`, `getCenter`; protected `isotope_stdev_`, `charge_`, `monoisotopic_mz_`, `averagine_[]`, `max_isotope_`, `trim_right_cutoff_`, `isotope_distance_` | `ExtendedIsotopeModel`: `charge`, `set_offset`, `offset`, `set_samples`, `center` |

### Fitters

| Source | Rust |
| --- | --- |
| `class Fitter1D : public DefaultParamHandler` (instantiable) | `trait Fitter1D` (`name`, `to_param`, `set_parameters`, `fit1d`); `struct BaseFitter1D` for the class itself |
| typedefs `IndexSet`, `ChargedIndexSet` | `feature_finder_algorithm::defs::{IndexSet, ChargedIndexSet}` |
| typedefs `CoordinateType` (double), `QualityType` (`float`), `PeakType`, `RawDataArrayType`, `PeakIterator` | `f64`, `fitter1d::QualityType = f32`, `Peak1D`, `&[Peak1D]`, slices |
| `Fitter1D()`, copy, `operator=`, destructor | `Default`, `Clone` |
| `fit1d(range, std::unique_ptr<InterpolationModel>&)` (throws `NotImplemented` in the base) | `fit1d(&self, &[Peak1D]) -> Result<(QualityType, Model1D)>`; `BaseFitter1D` returns `Unsupported` |
| protected `tolerance_stdev_box_`, `statistics_`, `interpolation_step_`, `updateMembers_` | `Fitter1DCore` fields and `read_param` |
| `LevMarqFitter1D`: `ContainerType`, `GenericFunctor` (`inputs`, `values`, `operator()`, `df`, `m_inputs`, `m_values`), constructor (`max_iteration` 500), copy, `operator=`; protected `symmetric_`, `max_iteration_`, `optimize_`, `updateMembers_`; file-local `GenericFunctorConstEigenAdapter` | `lev_marq_optimize(x, values, residual, jacobian, max_iteration)` with closures for the functor; `MAX_ITERATION_DEFAULT`; `EmgFitter1D::max_iteration`; `symmetric_` is internal to the EMG fit |
| `MaxLikeliFitter1D`: constructor, copy, `operator=`; protected `fitOffset_`, `updateMembers_` | `max_likeli_fit_offset(model, set, stdev1, stdev2, offset_step)` |
| `GaussFitter1D`, `BiGaussFitter1D` (protected `statistics1_`, `statistics2_`), `IsotopeFitter1D` (protected `charge_`, `isotope_stdev_`, `max_isotope_`), `ExtendedIsotopeFitter1D` (protected `charge_`, `isotope_stdev_`, `monoisotopic_mz_`, `max_isotope_`): constructor, copy, `operator=`, `fit1d`, `updateMembers_` | structs of the same names with public fields, `NAME`, `defaults`, `from_param`, `Default`, `Clone`, and `Fitter1D` |
| `EmgFitter1D`: `fit1d`; protected `Data {n, set}`, `EgmFitterFunctor` (with the static constants `c`, `sqrt2pi`, `emg_const`, `sqrt_2`), `setInitialParameters_`, `setInitialParametersMOM_`, `height_`, `width_`, `symmetry_`, `retention_`, `updateMembers_` | `EmgFitter1D` (`max_iteration`, `init_mom`); the functor is the private `emg_residuals`/`emg_jacobian`; the start values the private `initial_parameters(_mom)` |

### EmgScoring, PeakWidthEstimator, SeedListGenerator, FeatureFinderAlgorithm, FeatureFinderDefs

| Source | Rust |
| --- | --- |
| `EmgScoring()`, destructor | `EmgScoring::new`, `Default` |
| `setFitterParam`, `getDefaults` | `set_fitter_param`, `defaults` |
| `template calcElutionFitScore(MRMFeature&, MRMTransitionGroup<SpectrumType, TransitionT>&)` | `calc_elution_fit_score<C: NativeIdentified, T>(&MRMFeature, &MRMTransitionGroup<C, T>)` (`MSChromatogram` or `MSSpectrum`) |
| `elutionModelFit(points, smooth_data)` | `elution_model_fit(&[Point2D], bool)` |
| protected `fitRT_`, `prepareFit_`, `fitter_emg1D_params_` | private `fit_rt`, `prepare_fit`, `fitter_emg1d_params` |
| `PeakWidthEstimator(exp_picked, boundaries)` (hidden default constructor), destructor | `PeakWidthEstimator::new(&MSExperiment, &[Vec<PeakBoundary>])` |
| `getPeakWidth(mz)`; private `bspline_`, `mz_min_`, `mz_max_` | `peak_width(mz)`; private fields |
| `SeedListGenerator()`, typedef `SeedList` | `SeedListGenerator::new`, `type SeedList = Vec<DPosition2>` |
| `generateSeedList(const PeakMap&, SeedList&)` | `generate_seed_list_from_experiment` |
| `generateSeedList(PeptideIdentificationList&, SeedList&, bool)` | `generate_seed_list_from_peptides(&mut [PeptideIdentification], bool)` |
| `generateSeedLists(const ConsensusMap&, std::map<UInt64, SeedList>&)` | `generate_seed_lists -> BTreeMap<u64, SeedList>` |
| `convertSeedList(const SeedList&, FeatureMap&)`, `convertSeedList(const FeatureMap&, SeedList&)` | `seeds_to_features`, `features_to_seeds` |
| `struct Summary` (`exception`, `no_exceptions`, `mz_model`, `mz_stdev`, `charge`, `corr_mean`, `corr_max`, `corr_min`, constructor) | `Summary` (same fields, `mz_stdev: Vec<(f32, u32)>` kept sorted; `count_mz_stdev`), `Default` |
| `FeatureFinderAlgorithm` (`MapType`, `CoordinateType`, `IntensityType`, constructor, `run`, `getDefaultParameters`, `setData`, `setSeeds`, protected `map_`, `features_`, private copy) | `trait FeatureFinderAlgorithm` (`default_parameters`, `set_seeds` with the source default, `run(map, features)`); `SEEDS_NOT_SUPPORTED` |
| `FeatureFinderDefs` (`IndexPair`, `ChargedIndexSet`, `IndexSet`, `enum Flag`, `NoSuccessor`) | `feature_finder_algorithm::defs`, re-exporting the existing port in `feature_finder_picked::defs` |

## Preserved source conventions

- **`setParameters` replaces, it does not merge.** Every `set_parameters`
  completes the given tree from the defaults and forgets earlier setter calls.
  In particular the Gaussian, bi-Gaussian, EMG and isotope fitters call the new
  model's `setInterpolationStep(interpolation_step_)` and then
  `setParameters(tmp)`, so every fitted model is sampled at the model default
  step 0.1 whatever the fitter's `interpolation_step`; the fitter's step still
  sets the offset scan. The isotope fitters' `isotope_model:` subsection is
  likewise passed to the model and then replaced; it can only make a fit fail
  validation. Both are reproduced and listed as defect candidates below.
- **Sampling loops.** The Gaussian-family models sample at `min + i * step`
  for `i = 0, 1, ...` and include the first position at or past `max`; the
  isotope peak shapes accumulate `coord += step`. `min == max` clears the data
  and keeps the previous scale and offset.
- **Normalisation.** Gauss and bi-Gauss scale by `scaling / step / sum`; the
  isotope model by `scaling / (step * sum)`; the extended isotope model by
  `scaling / step / sum`; the EMG is not normalised. Sums are left to right.
- **`float` where the source has `float`.** Fit qualities are
  `Feature::QualityType` (`float`); correlations are computed on
  `std::vector<float>` copies of intensities and predictions; the EMG median
  start value accumulates in `float`; the isotope mean offset multiplies
  `float * Int` in `float`; isotope distributions use `Peak1D` binary32
  arithmetic (`ProbabilityPrecision::SourceSingle`).
- **NaN behaviour.** Bounding boxes start at the first position and ignore
  later NaNs; `fitOffset_` never lets a NaN correlation win; a NaN quality
  becomes -1 in every fitter; the EMG fitter's issue-#6239 guards are
  transcribed with their NaN-sensitive negated comparisons.
- **`EmgModel::setOffset` does not move `emg:retention`**, and
  `EmgModel::getCenter` is the `statistics:mean` parameter, not the EMG's
  maximum.
- **`IsotopeModel`** samples only on `set_samples_from_formula`; its
  `set_samples()` is the base class's `NotImplemented`. Its formula is built
  from `Int(0.5 + mass * atoms_per_dalton)` counts, zero counts omitted.
- **`SeedListGenerator::generateSeedLists`** pops the last seed of a
  sub-feature's map even when that seed belongs to an earlier consensus
  feature (two sub-features from one map); creating a list for a map without
  a column header leaves it empty.
- **`PeakWidthEstimator`** pairs spectra with boundary lists by position and
  uses the first and last *visited* peak's m/z as its range.

## Native differences

Each is stated at the item as well.

- **Errors instead of undefined behaviour.** An empty data set (the source
  reads `set[0]`), a zero charge (division by zero), a negative
  `isotope:maximum` (converted to a huge `Size`), an interpolation step larger
  than the isotope distance (the source reads past its shrunken pattern), a
  non-positive or NaN step with a non-empty range (the source loops forever),
  a zero-width Lorentzian (Boost throws `std::domain_error`) and an averagine
  count beyond 32 bits all return `Error::InvalidValue`/`InvalidRange`. Models
  and fitters are left unchanged by every error; the source's
  `setParameters` has already replaced `param_` when `checkDefaults` throws.
- **Resource ceilings.** `base_model::MAX_SAMPLES` (10,000,000 interpolation
  samples), `fitter1d::MAX_OFFSET_WORK` (offset scan points times data points),
  the Levenberg-Marquardt point ceiling and `MAX_RESIDUAL_WORK`, and
  `seed_list_generator::MAX_SEEDS` are checked before allocation.
- **Unknown parameters are not kept** in `to_param()` (the source keeps them in
  `param_` after warning), except the isotope fitters' `isotope_model:`
  subsection, which is kept because the source reads it.
- **`LinearInterpolation`** returns 0 (or NaN for a NaN position) where the
  source reads an empty container or converts a non-finite index; the
  derivative treats samples out of range as 0.
- **`EmgModel` and the EMG fitter** use `feature_finder_picked::glibc_libm::exp`,
  the reference build's `exp`, so their values are bit-identical to the Linux
  Release build on every host; the source's `pow(x, 2)` is computed as `x * x`
  (glibc's `pow` is within 0.52 ulp, so the two can differ in the last bit near
  a rounding midpoint). The method-of-moments cube root uses the host's `pow`.
  The Gaussian models use `BasicStatistics`, which calls the host's `exp`.
- **`ModelDescription::create_model`** is native.
- **`EmgScoring::calc_elution_fit_score`** refuses a feature without a convex
  hull (the source's precondition is debug-only), and propagates `MRMFeature::feature`'s error for
  a native ID without a feature, where the source's `feature_map_[key]`
  silently returns the first feature.
- **`SeedListGenerator`** maps an unset peptide RT or m/z (`None`) to NaN,
  the source's unset value, and does not pop an empty seed list (undefined in
  the source).
- **`Summary::mz_stdev`** is a sorted vector (a `float`-keyed `std::map` has no
  `Ord` key in Rust); NaN is refused.
- **`FeatureFinderAlgorithm`** passes the map and features to `run` instead of
  storing pointers with `setData`; the default `set_seeds` returns
  `Error::InvalidValue` for the source's `IllegalArgument`.

## Checked boundaries and evidence

- **Tier 1, executed differential** (`tests/featurefinder_oracle.rs`): 15,660
  values from `../oracle/FF1-models/driver.cpp` on **ibminode06** against the
  Linux x86_64 **Release** install `openms4-release-bc9cc12-c19e494-174b576`
  (conda gcc 14.4.0, libOpenMS core flags `-O3 -DNDEBUG ... -ffp-contract=off`,
  glibc 2.39), run twice with identical output: LinearInterpolation value,
  derivative and resampling; Gauss (with scaling and a shifted offset),
  bi-Gauss and two EMG models (all samples, raw data, offset, scale); the
  isotope model for the class-test formula with both peak shapes and for an
  averagine; the extended isotope model; Gauss, bi-Gauss, EMG (median and
  method-of-moments starts, and a skewed peak with a 37-evaluation budget),
  isotope (charge 2 and the charge-0 Gaussian path) and extended isotope fits
  (quality, center, every numeric model parameter, all samples); and
  EmgScoring's raw and smoothed fits of the three mock elution profiles.
  Values are compared within a relative 1e-9; the EMG models and the
  interpolation are asserted bit-identical.
- **Tier 3, class-test literals** (`tests/featurefinder_models.rs`,
  `tests/featurefinder_io.rs`): every section of the 20 class tests at
  `bc9cc12` with its persisting ClassTest tolerances, including the 173- and
  160-value isotope patterns, the EMG fit to exact samples, the four
  issue-#6239 robustness cases, EmgScoring's `0.924365639` and
  `0.981013417243958`, PeakWidthEstimator's three widths on the picked
  `PeakPickerHiRes_orbitrap.mzML`, and SeedListGenerator's seeds (the 13
  spectra of `PepXMLFile_test.mzML` are transcribed, as `tests/pepxml.rs`
  does, rather than committing the 1.4 MB file).
- **Not transcribed:** the second half of the EMG fitter class test perturbs
  the samples with `boost::mt19937` and `boost::normal_distribution`, which are
  not ported; EmgModel's three `ABORT_IF(isinf(...))` lines skip the rest of a
  section without asserting anything.
- **Tier 4:** the native refusals and ceilings above, `ModelDescription`
  round trips, `write_samples`, `Summary` and the default `set_seeds`.

## C++ defect candidates (unnumbered)

1. The model-building fitters' `setInterpolationStep` before `setParameters`
   is discarded (`GaussFitter1D.cpp:72-79`, `BiGaussFitter1D.cpp:73-80`,
   `IsotopeFitter1D.cpp:79-86` and `94-103`, `ExtendedIsotopeFitter1D.cpp:81-88`
   and `96-104`, `EmgFitter1D.cpp:159-169`): every fitted model is sampled at
   0.1.
2. `IsotopeFitter1D` and `ExtendedIsotopeFitter1D` pass the `isotope_model:`
   parameters to the model (`IsotopeFitter1D.cpp:92-94`) and then replace them
   with the second `setParameters` (`:103`), so they never take effect.
3. `BaseModel`'s constructor sets `defaults_` without `defaultsToParam_`, and
   `LevMarqFitter1D`'s adds `max_iteration` without it, leaving `cut_off_`,
   `max_iteration_` and `symmetric_` uninitialised until a derived constructor
   runs (`BaseModel.h:32-35`, `LevMarqFitter1D.h:63-66`).
4. `IsotopeModel::setSamples` and `ExtendedIsotopeModel::setSamples` read past
   the resized pattern when the interpolation step exceeds the isotope
   distance, and divide by a zero charge (`IsotopeModel.cpp:127-143`).
5. `EmgScoring::calcElutionFitScore` scores the first sub-feature for a native
   ID the MRM feature does not hold (`EmgScoring.h:67` through
   `MRMFeature::getFeature`).
6. `SeedListGenerator::generateSeedLists` pops from an empty seed list, and
   pops another consensus feature's seed, when two sub-features come from one
   map (`SeedListGenerator.cpp:83-88`).
7. `PeakWidthEstimator` pairs picked spectra with boundary lists by position,
   although `PeakPickerHiRes::pickExperiment` reports boundaries for picked
   spectra only (`PeakWidthEstimator.cpp:20-33`).

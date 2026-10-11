# OpenSWATH scoring

Ports of the OpenSWATH scoring headers at the core pin `bc9cc12`:

| C++ header | Rust | Status |
|---|---|---|
| `ANALYSIS/OPENSWATH/DIAHelper.h` | `src/analysis/openswath/dia_helper.rs` | complete |
| `ANALYSIS/OPENSWATH/DIAPrescoring.h` | `src/analysis/openswath/dia_prescoring.rs` | complete |
| `ANALYSIS/OPENSWATH/DIAScoring.h` | `src/analysis/openswath/dia_scoring.rs` | complete |
| `ANALYSIS/OPENSWATH/MRMScoring.h` | `src/analysis/openswath/mrm_scoring.rs` | complete |
| `ANALYSIS/OPENSWATH/OpenSwathScores.h` | `src/kernel/mrm.rs` (records), `src/analysis/openswath/openswath_scores.rs` (methods, usage) | complete |
| `ANALYSIS/OPENSWATH/IonMobilityScoring.h` | `src/analysis/openswath/ion_mobility_scoring.rs` | partial: ion-mobilogram peak picking (needs `PeakPickerMobilogram.h`) |
| `ANALYSIS/OPENSWATH/OpenSwathScoring.h` | `src/analysis/openswath/openswath_scoring.rs` | partial: inherits the ion-mobility peak-picking gap |

`convertPeptideToAASequence` of `ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h`
was added to `src/analysis/openswath/data_access_helper.rs`, because the DIA
and precursor scores need it.

## API mapping

### DIAHelper (`namespace DIAHelpers`)

| C++ | Rust |
|---|---|
| `integrateWindow(SpectrumPtr, mz, im, intensity, RangeMZ, RangeMobility, centroided)` | `integrate_window(&Spectrum, &RangeBase, &RangeBase, bool) -> Result<Option<WindowIntegral>>` |
| `integrateWindow(SpectrumSequence, ...)` | `integrate_window_sequence` |
| `integrateWindows(SpectrumPtr, centers, width, int, mz, im, RangeMobility, remove_zero)` | `integrate_windows(...) -> Result<IntegratedWindows>` |
| `integrateWindows(SpectrumSequence, ...)` | `integrate_windows_sequence` |
| `integrateWindow_` (public helper) | `integrate_window_sums(..., &mut WindowSums)` |
| `adjustExtractionWindow(right, left, window, ppm)` | `adjust_extraction_window(left, right, window, ppm) -> Result<(left, right)>` |
| `getBYSeries(AASequence, b, y, generator, charge)` | `by_series(&AASequence, &TheoreticalSpectrumGenerator, charge) -> Result<(b, y)>` |
| `getTheorMasses` | `theoretical_masses` |
| `getAveragineIsotopeDistribution(mz, out, charge, nr_isotopes, mannmass)` | `averagine_isotope_distribution(mz, charge, nr_isotopes, mannmass)` |
| `simulateSpectrumFromAASequence` | `simulate_spectrum_from_sequence -> Result<SimulatedSpectrum>` |
| `addPreisotopeWeights(vector, ...)` | `add_preisotope_weights` |
| `addPreisotopeWeights(double, ...)` | `add_preisotope_weights_single` |
| `addIsotopes2Spec` | `add_isotopes_to_spec` |
| `addSinglePeakIsotopes2Spec` | `add_single_peak_isotopes_to_spec` |
| `sortByFirst` / `extractFirst` / `extractSecond` | `sort_by_first` / `extract_first` / `extract_second` |
| `createMZRangePPM` | `create_mz_range_ppm` |
| `DEFAULT_MANN_MASS` (file-local) | `DEFAULT_MANN_MASS`; the header default `1.000482` of `addPreisotopeWeights` is `DEFAULT_PREISOTOPE_MANN_MASS` |

### DIAPrescoring (`class DiaPrescore`)

| C++ | Rust |
|---|---|
| `DiaPrescore()` / `defineDefaults` | `DiaPrescore::new` / `defaults` |
| `DiaPrescore(window, nr_isotopes, nr_charges)` | `DiaPrescore::with_window` |
| `updateMembers_` / `getParameters` | `from_param`, `from_param_with_warnings` / `to_param` |
| `TransitionGroupTheoreticalSpectrumCache` (+ `clear`) | `TransitionGroupTheoreticalSpectrumCache` (+ `clear`) |
| `score(spec, lt, im_range, dotprod, manhattan)` | `score(...) -> Result<PrescoreResult>` |
| `buildTheoreticalSpectrum(lt, cache)` | `build_theoretical_spectrum(lt) -> Result<cache>` |
| static `scorePrepared` | associated `score_prepared` |
| `operator()(swath_ptr, exp, im_range, IDataFrameWriter*)` | `score_all(&dyn SpectrumAccess, &LightTargetedExperiment, &RangeBase, &mut dyn DataFrameWriter)` |
| file-local `getNormalizedLibraryIntensities`, `getMZIntensityFromTransition` | not ported: defined in `DIAPrescoring.cpp`, declared in no header and called nowhere |

### DIAScoring (`class DIAScoring`)

| C++ | Rust |
|---|---|
| `DIAScoring()`, parameters, `updateMembers_` | `DiaScoring::new`, `defaults`, `from_param(_with_warnings)`, `to_param` |
| `dia_isotope_scores(transitions, spectrum, IMRMFeature*, im_range, corr, overlap)` | `dia_isotope_scores(..., &dyn MRMFeatureAccess, ...) -> Result<IsotopeScores>` |
| `dia_massdiff_score(..., ppm_score, ppm_score_weighted, diff_ppm)` | `dia_massdiff_score -> Result<MassDiffScores>` |
| `dia_ms1_massdiff_score(...) -> bool` | `dia_ms1_massdiff_score -> Result<Ms1MassDiffScore>` (`signal_found` is the return value) |
| `dia_ms1_isotope_scores_averagine` | `dia_ms1_isotope_scores_averagine -> Result<IsotopeScores>` |
| `dia_ms1_isotope_scores(..., EmpiricalFormula)` | `dia_ms1_isotope_scores(..., &EmpiricalFormula)` |
| `dia_by_ion_score` | `dia_by_ion_score -> Result<ByIonScores>` |
| `score_with_isotopes` | `score_with_isotopes -> Result<PrescoreResult>` |
| `getDIAExtractionWindow` / `isDIAExtractionPPM` | `dia_extraction_window` / `is_dia_extraction_ppm` |
| private `diaIsotopeScoresSub_`, `largePeaksBeforeFirstIsotope_`, `scoreIsotopePattern_` (3), `getIsotopeIntysFromExpSpec_` | private helpers of the same behaviour |
| private copy constructor / assignment (deleted semantics) | `Clone` (the scorer has no shared state) |

### MRMScoring (`class OpenSwath::MRMScoring`)

| C++ | Rust |
|---|---|
| `XCorrMatrixType` (`OpenMS::Matrix<XCorrArrayType>`) | `XCorrMatrixType = ScoreMatrix<XCorrArrayType>` |
| `getXCorrMatrix`, `getXCorrContrastMatrix`, `getXCorrPrecursorContrastMatrix`, `getXCorrPrecursorCombinedMatrix` | `xcorr_matrix`, `xcorr_contrast_matrix`, `xcorr_precursor_contrast_matrix`, `xcorr_precursor_combined_matrix` (+ native `xcorr_precursor_matrix`) |
| `getMIMatrix`, `getMIContrastMatrix`, `getMIPrecursorContrastMatrix`, `getMIPrecursorCombinedMatrix` | `mi_matrix`, `mi_contrast_matrix`, `mi_precursor_contrast_matrix`, `mi_precursor_combined_matrix` (+ native `mi_precursor_matrix`) |
| `initializeXCorrMatrix(data)` / `(IMRMFeature*, ids)` | `initialize_xcorr_matrix_from_data` / `initialize_xcorr_matrix` |
| `initializeXCorrContrastMatrix`, `initializeXCorrPrecursorMatrix`, `initializeXCorrPrecursorContrastMatrix` (feature and data overloads), `initializeXCorrPrecursorCombinedMatrix` | same names in snake case; the data overload is `..._from_data` |
| `initializeMIMatrix`, `initializeMIContrastMatrix`, `initializeMIPrecursorMatrix`, `initializeMIPrecursorContrastMatrix`, `initializeMIPrecursorCombinedMatrix` | same names in snake case |
| every `calcXcorr*Score` (13), `calcSeparateXcorrContrast*` (2), `calcMI*Score` (5), `calcSeparateMIContrastScore` | same names in snake case |
| static `calcLibraryScore(..., six out-parameters)` | `calc_library_score -> Result<LibraryScores>` |
| static `calcRTScore`, `calcSNScore`, `calcSeparateSNScore` | `calc_rt_score`, `calc_sn_score`, `calc_separate_sn_score` |
| typedefs `SpectrumType`, `TransitionType`, `PeptideType`, `ProteinType`, `FeatureType` | the Rust types are used directly |

### OpenSwathScores

| C++ | Rust |
|---|---|
| `OpenSwath_Scores` data members | `kernel::mrm::OpenSwathScores` (earlier port) |
| `get_quick_lda_score(...) const` | associated `OpenSwathScores::get_quick_lda_score` |
| `calculate_lda_prescore(scores) const`, `calculate_lda_single_transition`, `calculate_swath_lda_prescore` | `&self` methods; the source ignores `this` and every caller passes the object itself |
| `OpenSwath_Scores_Usage` (20 switches) | `OpenSwathScoresUsage` (trailing underscores dropped) |
| `OpenSwath_Ind_Scores` | `kernel::mrm::OpenSwathIndScores` (earlier port) |

### IonMobilityScoring

| C++ | Rust |
|---|---|
| `driftScoring(spectra, transitions, scores, drift_target, im_range, window, ppm, drift_extra, peak_picking, consumer, feature_rt, feature_id)` | `drift_scoring(spectra, transitions, &mut scores, drift_target, im_range, &DriftExtraction, peak_picking, Option<MobilogramExport>)` |
| `driftScoringMS1`, `driftScoringMS1Contrast` | `drift_scoring_ms1`, `drift_scoring_ms1_contrast` |
| `driftIdScoring(..., MRMTransitionGroupType& detection, ...)` | `drift_id_scoring(..., detection_transitions: &[LightTransition], ...)` |
| `computeIonMobilogram(spectra, mz, im, im&, intensity&, Mobilogram&, eps)` | `compute_ion_mobilogram -> Result<IonMobilogram>` |
| private `computeGrid_`, `alignToGrid_`, `extractIntensities` (2) | public `compute_grid`, `align_to_grid`, `extract_intensities`, `extract_intensities_single` |
| free `sumAlignedMobilograms` | `sum_aligned_mobilograms` |
| `MobilogramParquetConsumer*` | `MobilogramSink` trait (the Parquet implementation is a backlog item of `FORMAT/DATAACCESS/MobilogramParquetConsumer.h`) |
| `PeakPickerMobilogram` use for `apply_im_peak_picking` | **not ported**: `Error::Unsupported` naming `ANALYSIS/OPENSWATH/PeakPickerMobilogram.h` |

### OpenSwathScoring

| C++ | Rust |
|---|---|
| private enums `SpectrumAdditionMethod`, `SpectrumMergeMethodType` | public `SpectrumAdditionMethod`, `SpectrumMergeMethodType` |
| private members | public fields |
| `initialize(...)` | `initialize(...) -> Result<()>` |
| `calculateChromatographicScores`, `calculateChromatographicIdScores`, `calculateLibraryScores` | same, snake case, `-> Result<()>` |
| `calculateDIAScores(..., masserror_ppm&, drift_target, range_im, consumer, feature_id)` | `calculate_dia_scores(..., Option<&mut dyn MobilogramSink>, feature_id) -> Result<Vec<f64>>` (the mass errors) |
| `calculatePrecursorDIAScores`, `calculateDIAIdScores` | same, snake case (`calculate_dia_id_scores` takes the detection transitions) |
| `fetchSpectrumSwath(SpectrumAccessPtr, ...)` returning / out-parameter | `fetch_spectrum_swath` / `fetch_spectrum_swath_into` |
| `fetchSpectrumSwath(vector<SwathMap>, ...)` returning / out-parameter | `fetch_spectrum_swath_maps` / `fetch_spectrum_swath_maps_into` |

## Preserved source conventions

- Window integration takes peaks with `min <= mz < max` (binary-search lower
  bound, exclusive upper bound) and the closed ion-mobility range; the
  intensity-weighted means and the `-1` ion-mobility sentinel without an IM
  range are kept, as the scores consume them.
- `DIAScoring.cpp` uses its own global `C13C12_MASSDIFF_U = 1.0033548`, which
  unqualified lookup inside `namespace OpenMS` resolves to; `DIAPrescoring.cpp`
  uses `Constants::C13C12_MASSDIFF_U = 1.0033548378`. Both are kept where the
  source uses them (`DIA_SCORING_C13C12_MASSDIFF_U`, `chemistry::C13C12_MASSDIFF_U`).
- Averagine distributions use the source's binary32 `Peak1D` intensities
  (`ProbabilityPrecision::SourceSingle`); the oracle pins those bits.
- The weight given to the averagine model is `product_mz * |charge|` (the source
  notes this leaves the charge carriers in).
- `score_with_isotopes` passes the extraction window as Th even when the unit
  is ppm, as the source does.
- The contrast shape score (`calcXcorrContrastShapeScore`) is a **sum**, not a
  mean; MI triangle means divide by `n*n/2 + (n+1)/2` in integer arithmetic.
- `sumAlignedMobilograms` adds `f32` intensities; ion-mobilogram bins sum in
  `f64` and are stored as the mobilogram's `f32`.
- `alignToGrid_`'s maximum index is taken from the next profile point after
  each step, as in the source.
- The two `fetchSpectrumSwath` forms differ as in the source: the returning
  single-map form with resampling returns one merged spectrum even when
  nothing was fetched (an empty spectrum; the source's `addUpSpectra` indexes
  the empty vector there, which is undefined behaviour); the out-parameter form
  returns nothing.

## Native differences

- Out-parameters became return values; failures are `Result` errors, and every
  scoring call leaves its scores untouched on error.
- Checked where the source reads out of bounds or checks only in debug builds:
  array lengths, `normalized_library_intensity` length, charge 0 and
  non-positive fragment charges, empty transition or SWATH-map lists, missing
  features (`Error::MissingInformation` instead of a null dereference), IM bin
  indices beyond `i32`.
- `addPreisotopeWeights(double, ...)` refuses charge 0 (the source divides by
  zero); `DiaPrescore` refuses negative `nr_isotopes` (the source converts it to
  a huge unsigned count).
- Charge-0 isotope windows (non-finite positions) find no signal, reproducing
  the source's NaN/infinite windows (pinned by the oracle's `glucose0` case).
- Thread-local pools and caches (averagine cache of 100000 entries, b/y series,
  prescore spectra, formulas, AASequence, a persistent `MRMScoring`) are not
  ported; results are recomputed and identical.
- The source's debug-only `OPENSWATH_PRECONDITION`s on matrix shapes are not
  enforced; the release build computes on and so does the port.
- `ScoreMatrix` is a column-major stand-in for `OpenMS::Matrix`
  (`DATASTRUCTURES/Matrix.h`, not ported); every initialisation starts zeroed.
- `sortByFirst` is stable (`std::sort` is not); NaN sorts last.
- `DiaPrescore::score_all` does not print the per-spectrum progress line.
- Log output (`OPENMS_LOG_WARN`/`DEBUG`) is not reproduced.

## Checked boundaries and evidence

- **Tier 1 (executed differential):** `tests/openswath_scoring_oracle.rs`
  against `tests/data/openswath_scoring/oracle_all.stdout`, produced by
  `../oracle/O2-scoring/driver.cpp` built with x86_64-conda-linux-gnu-g++
  14.4.0 `-O3 -DNDEBUG` against the Linux x86-64 Release install
  `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`
  (libOpenMS.so sha256 `abd4fc99...66d6c`) on **ibminode06**, run twice with
  byte-identical output. 1014 records cover every DIAHelper function (including
  the averagine values of the `#if 0` class-test section), all DIAScoring and
  DiaPrescore scores (including the untested sum-formula isotope scores), every
  MRMScoring score (including the untested contrast-shape, precursor and MI
  precursor scores), the IonMobilityScoring routines and `computeIonMobilogram`,
  and OpenSwathScoring's DIA, precursor, chromatographic, ID and library scores
  in both addition modes.
- **Tolerance:** records compare **exactly** (bit for bit), NaN to NaN, except one declared
  class compared with relative tolerance 1e-12 (observed 1-3 ulp): fragment m/z from the chemistry
  port's `TheoreticalSpectrumGenerator` (`by2.*`, `theor1.*`): the source starts from
  `PROTON_MASS_U * charge` and adds residue masses (`TheoreticalSpectrumGenerator.cpp:724-794`), the
  port adds the protons last. This is outside the OpenSWATH port and is reported as a chemistry-port
  divergence; the b/y *counts* built from those m/z values match exactly. Two further classes the
  lane declared (cross-correlation apex values of `Scoring.cpp`, and the `StatsHelpers.cpp`
  `dotProd`/`norm`/`manhattanDist`/sums behind the prescore, `dotprod_score_dia`,
  `manhatt_score_dia`, `lib.dot` and `lib.man`) were accumulation-order differences against Eigen's
  blocked reductions; integration (group O) moved both modules onto the measured Eigen 5.0.1
  two-lane kernels and they are compared exactly.
- **Tier 3 (source review):** every section of `DIAHelper_test.cpp`,
  `DIAPrescoring_test.cpp`, `DIAScoring_test.cpp`, `MRMScoring_test.cpp`,
  `IonMobilityScoring_test.cpp`, `OpenSwathScores_test.cpp` and
  `OpenSwathScoring_test.cpp` (the four `calculate*` sections are NOT_TESTABLE
  upstream) in `tests/openswath_dia_scoring.rs`, `tests/openswath_mrm_scoring.rs`,
  `tests/openswath_ion_mobility_scoring.rs` and `tests/openswath_scoring.rs`.
- **Tier 4:** native guards and composition checks marked "derived".

## Source defect candidates

- `OpenSwathScoring.h` forward-declares `struct RangeMZ; struct RangeMobility;`
  in the **global** namespace (before `namespace OpenMS`), so any translation
  unit with `using namespace OpenMS` that includes it cannot name `RangeMZ` or
  `RangeMobility` unqualified (ambiguous reference; hit by the oracle driver).
- `DIAScoring.cpp` shadows `Constants::C13C12_MASSDIFF_U` (1.0033548378) with a
  file-global 1.0033548, so DIAScoring and DiaPrescore space isotopes
  differently.
- `DIAPrescoring.cpp` comments "WARNING: This was spectrumWIso ... Was this a
  bug?" about the signed spectrum; `nr_charges_` is stored but never used.
- `IonMobilityScoring.h` declares `sumAlignedMobilograms` without
  `OPENMS_DLLAPI`, so the free function is not exported from libOpenMS.
- `DiaPrescore(double, int, int)` does not call `defineDefaults`, so its
  `getParameters()` is empty.

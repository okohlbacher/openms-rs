# TargetedSpectraExtractor: native header equivalent (partial)

[`src/analysis/openswath/targeted_spectra_extractor.rs`](../src/analysis/openswath/targeted_spectra_extractor.rs)
ports core SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4`
`ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h` with its `.cpp`, except two
members that need components the crate does not have yet (see
[Gaps](#gaps)).

Tests: [`tests/targeted_spectra_extractor.rs`](../tests/targeted_spectra_extractor.rs)
(8 tests without features; 12 with `mzml`, which loads the class test's mzML).
Manifest: [`tests/data/targeted_spectra_extractor_provenance.json`](../tests/data/targeted_spectra_extractor_provenance.json).

## API mapping

| Source member | Native representation |
| --- | --- |
| `TargetedSpectraExtractor()` | `TargetedSpectraExtractor::new() -> Result<Self>` |
| `DefaultParamHandler` base: `getDefaults`, `getParameters`, `setParameters`, `updateMembers_` | `defaults()`, `to_param()`, `from_param()` / `from_param_with_warnings()` (the crate's parameter pattern; the fields are the source's cached members) |
| `getDefaultParameters(Param&)` | `default_parameters() -> Result<Param>` |
| `struct Match` (`spectrum`, `score`, both constructors) | `Match` (`spectrum`, `score`), `Default`, `Match::new` |
| `class Comparator` (`generateScores`, `init`, `getLibrary`, protected `library_`) | trait `Comparator` (`generate_scores`, `init`, `library`) |
| `class BinnedSpectrumComparator` (`generateScores`, `init`, private `cmp_bs_`, `bs_library_`, `bin_size_`, `peak_spread_`, `bin_offset_`) | `BinnedSpectrumComparator` (`new`, `Comparator` impl, accessors `bin_size`, `peak_spread`, `bin_offset`) |
| `annotateSpectra(spectra, targeted_exp, annotated, features, compute_features = true)` | `annotate_spectra(spectra, targeted_exp) -> Result<(Vec<MSSpectrum>, FeatureMap)>` |
| `annotateSpectra(spectra, targeted_exp, annotated)` and `compute_features = false` | `annotate_spectra_only(..) -> Result<Vec<MSSpectrum>>` |
| `annotateSpectra(spectra, ms1_features, ms2_features, annotated)` | `annotate_spectra_with_ms1_features(spectra, ms1_features, &mut UniqueIdGenerator)` |
| `searchSpectrum(feat_map, feat_map_output, add_unidentified_features)` | **not ported** (gap 1) |
| `pickSpectrum(spectrum, picked)` | `pick_spectrum(&MSSpectrum) -> Result<MSSpectrum>` |
| `scoreSpectra(annotated, picked, features, scored, compute_features = true)` | `score_spectra(annotated, picked, &mut FeatureMap) -> Result<Vec<MSSpectrum>>` |
| `scoreSpectra(annotated, picked, scored)` and `compute_features = false` | `score_spectra_only(annotated, picked)` |
| `selectSpectra(scored, features, selected, selected_features, compute_features = true)` | `select_spectra(scored, features) -> Result<(Vec<MSSpectrum>, FeatureMap)>` |
| `selectSpectra(scored, selected)` and `compute_features = false` | `select_spectra_only(scored)` |
| `extractSpectra(experiment, targeted_exp, spectra, features, compute_features = true)` | `extract_spectra(experiment, targeted_exp)` |
| `extractSpectra(experiment, targeted_exp, spectra)` | `extract_spectra_only(experiment, targeted_exp)` |
| `extractSpectra(experiment, ms1_features, spectra, features)` and the private five-argument overload | `extract_spectra_with_ms1_features(experiment, ms1_features, ids)` |
| `extractSpectra(experiment, ms1_features, spectra)` | `extract_spectra_with_ms1_features_only(..)` |
| `matchSpectrum(input, cmp, matches)` | `match_spectrum(&MSSpectrum, &dyn Comparator) -> Result<Vec<Match>>` |
| `targetedMatching(spectra, cmp, features)` | `targeted_matching(spectra, cmp, &mut FeatureMap) -> Result<Vec<usize>>` (the unmatched indices the source logs) |
| `untargetedMatching(spectra, cmp, features)` | `untargeted_matching(spectra, cmp) -> Result<(FeatureMap, Vec<usize>)>` |
| `constructTransitionsList(ms1, ms2, t_exp)` | `construct_transitions_list(ms1, ms2, &mut TargetedExperiment)` |
| `storeSpectraMSP(filename, experiment)` | **not ported** (gap 2); its two preparation steps are |
| `mergeFeatures(fmap_input, fmap_output)` | `merge_features(input, &mut output, &mut UniqueIdGenerator)` |
| protected `deisotopeMS2Spectra_(experiment)` | public `deisotope_ms2_spectra(&mut MSExperiment)` |
| protected `removeMS2SpectraPeaks_(experiment)` | public `remove_ms2_spectra_peaks(&mut MSExperiment)` |
| protected `organizeMapWithSameIdentifier(fmap_input, fmapmap)` | public `organize_map_with_same_identifier(&FeatureMap) -> Result<BTreeMap<..>>` |
| private members `mz_unit_is_Da_` ... `max_precursor_mass_threashold_` | public fields `mz_unit_is_da` ... `max_precursor_mass_threashold`, `deisotoping: DeisotopingParameters`, `subsection_parameters` |

## Preserved source conventions

- **Parameters.** The tree is the source's: `getDefaultParameters` entries in
  their order with their descriptions, then the constructor's subsections
  `SavitzkyGolayFilter:` (frame length 15, order 3), `GaussFilter:`
  (width 0.2), `PeakPickerHiRes:` (signal-to-noise 1.0), and
  `AccurateMassSearchEngine:` (that class's constructor defaults, transcribed,
  with `keep_unidentified_masses` set to `false`). Only the first three are
  registered subsections, so their keys are not checked; `pick_spectrum`
  applies them over the filter's and the picker's own defaults on every call.
- **Annotation.** The targeted overload matches the peptide retention time
  (converted from minutes) within `rt_window / 2` and the precursor m/z within
  `±mz_tolerance`; the MS1-feature overload uses `mz_tolerance / 2` in m/z
  units regardless of `mz_unit_is_Da`, skips MS1 spectra and features whose
  `PeptideRef` is `"null"`, and uses subordinates instead of a feature when it
  has any. A spectrum without precursors has m/z 0, which disables the m/z
  check; the targeted overload's lower limit is then
  `numeric_limits<double>::min()`, the smallest positive double.
- **Picking.** Gaussian (or Savitzky-Golay) smoothing, then `PeakPickerHiRes`
  with spacing constraints disabled and absolute FWHM; peaks filtered by
  height and by the FWHM in the first float data array; no surviving peak
  gives a fully cleared spectrum.
- **Scoring.** `log10(TIC) * tic_weight + 1/avgFWHM * fwhm_weight + avgSNR *
  snr_weight` with the median S/N estimator at a 40-Th window, noise 2 for
  sparse windows and at least 10 points; the scored copy's float arrays are
  resized to five and arrays 1-4 carry `score`, `log10_total_tic`,
  `inverse_avgFWHM`, `avgSNR` as one-element arrays.
- **Selection** keeps, per name, the first spectrum unless a later one scores
  strictly higher, skipping scores below `min_select_score`, in name order.
- **Matching** sorts with the source's `std::sort` (`math::source_sort`), so
  tied scores come out in libstdc++ order.
- **Transitions list**: first-wins maps for proteins and transitions; a
  repeated peptide collects retention times; native ids stream the retention
  time as `std::ostream` does (`110`, `1.23457e+06`).
- **Merging**: `std::to_string` (six decimals) in group keys; weights from
  `peak_apex_int` or intensity; the `PeptideRef` before the first `_`.

## Native differences

- **Unique ids** come from a caller-owned `UniqueIdGenerator` instead of the
  source's process-wide one, so they are reproducible.
- **Return values instead of out-parameters.** The MS1-feature annotation
  returns its MS2 features; the source appends them to the caller's map.
- **Atomic failures.** `score_spectra`, `targeted_matching`,
  `construct_transitions_list`, `merge_features`, `deisotope_ms2_spectra` and
  `remove_ms2_spectra_peaks` leave their in/out argument unchanged on error.
  The source replaces the target experiment's contents before validating
  them, and `mergeFeatures` swallows every exception (logging it) after
  appending the groups processed so far; here the error is returned.
- **Checked instead of undefined**: mismatched picked/annotated counts, a
  scored spectrum without a score array, a picked spectrum without an FWHM
  array, an MS2 spectrum without a precursor in `remove_ms2_spectra_peaks`, a
  comparator reporting an index outside its library, and non-string
  `PeptideRef`/`native_id`/`Comments`/`identifier` values (the source's
  `ConversionError`) are errors.
- **Bounded work**: an annotation that would copy more than
  `MAX_ANNOTATIONS` (1,000,000) spectra is refused by a counting pass before
  anything is copied.
- **Logging**: the source's debug, info and warning lines are not written;
  `targeted_matching` returns the unmatched indices it would have logged.
- **Unsorted input**: `pick_spectrum`'s refusal of an unsorted spectrum is
  `Error::InvalidValue` with the source's message
  ("Spectrum must be sorted by position").
- Non-finite scores (an empty picked or annotated spectrum) are stored as the
  source stores them, in meta values the crate otherwise refuses; its
  validation reports them downstream.

## Gaps

1. **`searchSpectrum`** runs `ANALYSIS/ID/AccurateMassSearchEngine` (not in the
   crate; header status `unmapped`). The parameter subsection is carried.
2. **`storeSpectraMSP`** writes through `FileHandler` to MSP with
   `MSPGenericFile::store`. Lane F4 has ported `MSPGenericFile`
   (`src/format/msp_generic.rs`, branch `worktree-wf_cd066dab-0a7-9`); once it
   is integrated, `storeSpectraMSP` is `deisotope_ms2_spectra` (when
   `deisotoping:use_deisotoper`), `remove_ms2_spectra_peaks` and that store.
3. The class test's matching sections (`matchSpectrum`, `targetedMatching`,
   `untargetedMatching`) and the MS1-feature `extractSpectra` section load MSP
   libraries with `MSPGenericFile`, and the targeted sections load the target
   list with `TransitionTSVFile` (unmapped). The targeted sections run here on a
   target list built in the test from the same CSV (see the test's
   `targeted_experiment`); the matching logic is tested with a fixed-score
   comparator and a self-match of `BinnedSpectrumComparator`.

## Checked boundaries and evidence

| Tier | Source | What it covers |
| --- | --- | --- |
| 3, source review | `TargetedSpectraExtractor_test.cpp` literals, with its `TEST_REAL_SIMILAR` tolerance | parameters; `pickSpectrum` on the inline C-aconitate spectrum (three filter settings, the unsorted refusal); `annotateSpectra` (30 spectra, names, sizes, feature RT/m/z); `scoreSpectra` (all 24 asserted values and the size refusal, both overloads); `selectSpectra` and `extractSpectra` (both overloads); `mergeFeatures`; MS1-feature `annotateSpectra`; `constructTransitionsList`; `storeSpectraMSP`'s in-place peak removal |
| 4, independently derived | fixed-score comparator; contrast angle of a spectrum with itself | `match_spectrum` ordering, cut-off and top-N; `targeted_matching` meta values and unmatched indices; option parsing of `BinnedSpectrumComparator` |

The class test's mzML (`TargetedSpectraExtractor_13C1_spectra0to100.mzML`) and
target list (`TargetedSpectraExtractor_13CFlux_TraML.csv`) are copied from the
core checkout at the pin, with their sha256 in the manifest.

## Candidate C++ issues

1. `annotateSpectra` (targeted) turns a ppm `mz_tolerance` into
   `mz_tolerance / 1e6` without multiplying by the m/z, so a ppm tolerance acts
   as a tiny absolute one; `pickSpectrum` does the same with `fwhm_threshold`
   against an absolute FWHM.
2. `annotateSpectra` (targeted) uses `numeric_limits<double>::min()` (the
   smallest positive double) as the "no lower limit" value, which excludes a
   target m/z of 0.
3. `pickSpectrum` reads the FWHM from float data array 0, which is the
   ion-mobility array when the input carries one.
4. `scoreSpectra` resizes the float data arrays to five, truncating any the
   annotated spectrum had beyond that and appending to arrays 1-4 if they
   already held data.
5. `targetedMatching` and `untargetedMatching` temporarily mutate
   `top_matches_to_report_`, which makes them non-`const` and not re-entrant.
6. `removeMS2SpectraPeaks_` reads `getPrecursors()[0]` without a check.
7. `mergeFeatures` catches every exception, logs it and returns normally.

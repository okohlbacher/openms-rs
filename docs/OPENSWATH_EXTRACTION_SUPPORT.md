# OpenSWATH extraction and SWATH windows

Source pin: OpenMS4-core `bc9cc12`. Rust: `src/analysis/openswath/`. Tests:
`tests/openswath_extraction.rs`. Provenance: `tests/data/openswath_extraction_provenance.json`.
Oracle: `../oracle/O4a-extraction/` (outside the crate; C++ driver, scripts, results).

The data types these ports build on (`SpectrumAccess`, light `Spectrum`/`Chromatogram`,
`SwathMap`, `LightTargetedExperiment`) are described in
[OPENSWATH_DATAACCESS_SUPPORT.md](OPENSWATH_DATAACCESS_SUPPORT.md).

| Header | Rust | Status |
|---|---|---|
| `ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h` | `chromatogram_extractor_algorithm.rs` | complete |
| `ANALYSIS/OPENSWATH/ChromatogramExtractor.h` | `chromatogram_extractor.rs` | partial: `return_chromatogram<TargetedExperiment>` |
| `ANALYSIS/OPENSWATH/PeakMapExtractor.h` | `peak_map_extractor.rs` | complete |
| `ANALYSIS/OPENSWATH/SwathWindowLoader.h` | `swath_window_loader.rs` | complete |
| `ANALYSIS/OPENSWATH/SwathQC.h` | `swath_qc.rs` | complete |
| `ANALYSIS/OPENSWATH/SwathMapMassCorrection.h` | `swath_map_mass_correction.rs` | complete (`correctMZ`, `correctIM` ported at group O integration) |

## API mapping

### `ChromatogramExtractorAlgorithm.h` → `chromatogram_extractor_algorithm.rs`

| Source | Rust |
|---|---|
| class (a `ProgressLogger`), default ctor/dtor | `ChromatogramExtractorAlgorithm` (unit struct), `new`, `Default` |
| `ExtractionCoordinates { mz, ion_mobility, mz_precursor, rt_start, rt_end, id }` | `ExtractionCoordinates` (same fields, all default 0 / empty) |
| `ExtractionCoordinates::SortExtractionCoordinatesByMZ` | `ExtractionCoordinates::sort_by_mz` |
| `ExtractionCoordinates::SortExtractionCoordinatesReverseByMZ` | `ExtractionCoordinates::sort_reverse_by_mz` |
| — | native `ExtractionCoordinates::cmp_by_mz` (an `Ordering` for `sort_by`) |
| `extractChromatograms(input, output, coords, mz_window, ppm, im_window, filter)` | `extract_chromatograms(&dyn SpectrumAccess, &mut [Chromatogram], ..) -> Result<()>`; `extract_chromatograms_with_progress(.., &mut ProgressLogger)` reports as the `ProgressLogger` base does |
| `extract_value_tophat(mz_start, mz_it, mz_end, int_it, mz, integrated, window, ppm)` | `extract_value_tophat(&[f64] mz, &[f64] int, &mut usize pos, mz, window, ppm) -> Result<f64>` |
| `extract_value_tophat(.., im_it, mz, im, integrated, mz_window, im_window, ppm)` | `extract_value_tophat_im(mz, int, im, &mut pos, mz, im, mz_window, im_window, ppm) -> Result<f64>` |
| private `getFilterNr_` | `ExtractionFilter::from_name` (`Tophat`, `Bartlett`) |
| — | `MAX_EXTRACTED_POINTS`, `EXTRACTION_PROGRESS_LABEL` |

### `ChromatogramExtractor.h` → `chromatogram_extractor.rs`

| Source | Rust |
|---|---|
| class (a `ProgressLogger`), ctor | `ChromatogramExtractor` (unit struct), `new`, `Default` |
| typedef `ExtractionCoordinates` | re-export `ExtractionCoordinates` |
| `extractChromatograms` (6 args; IM window -1) | `ChromatogramExtractor::extract_chromatograms` |
| `extractChromatograms` (7 args) | `ChromatogramExtractor::extract_chromatograms_im` |
| static `prepare_coordinates(output, coords, LightTargetedExperiment, rt_window, ms1 = false, ms1_isotopes = 0)` | free fn `prepare_coordinates(&mut Vec<Chromatogram>, &mut Vec<ExtractionCoordinates>, &LightTargetedExperiment, f64, bool, i32) -> Result<()>` |
| static template `return_chromatogram<LightTargetedExperiment>` (also `const`) | free fn `return_chromatogram(.., settings: &MSSpectrum, &mut Vec<MSChromatogram>, ms1, im_extraction_width) -> Result<()>` |
| static template `return_chromatogram<TargetedExperiment>` | **not ported**: the heavyweight `ANALYSIS/TARGETED/TargetedExperiment` is not ported |
| private `extract_id_<LightTargetedExperiment>` (+ `const`) | private `extract_id` |
| private `extract_id_<TargetedExperiment>` | not ported (heavyweight `TargetedExperiment`) |
| private `prepareSpectra_`, `outsideExtractionWindow_`, `getFilterNr_`, `populatePeptideRTMap_`, `PeptideRTMap_` | not ported: private, and no public member of the pinned source calls them (the source marks three "TODO deprecate") |
| file-local `populateMS1Transition`, `populateMS2Transition` | inlined in `prepare_coordinates` |
| — | `MAX_COORDINATES` |

### `PeakMapExtractor.h` → `peak_map_extractor.rs`

| Source | Rust |
|---|---|
| class (a `ProgressLogger`) | `PeakMapExtractor` (unit struct), `new`, `Default` |
| using `ExtractionCoordinates` | `chromatogram_extractor_algorithm::ExtractionCoordinates` |
| `ExtractedPeakMap { native_id, target_mz, target_rt = NaN, target_ion_mobility = -1, rt_start = 0, rt_end = -1, mz, rt, ion_mobility, intensity }` | `ExtractedPeakMap` (same fields and defaults); native `for_coordinate`, `len`, `is_empty` |
| `extractPeakMaps(input, output, coords, mz_window, ppm, im_window, filter = "tophat")` | `extract_peak_maps(..) -> Result<Vec<ExtractedPeakMap>>`; `extract_peak_maps_with_progress` |
| private `getFilterNr_` | inline check |
| — | `MAX_PEAK_MAP_POINTS`, `PEAK_MAP_PROGRESS_LABEL` |

### `SwathWindowLoader.h` → `swath_window_loader.rs`

| Source | Rust |
|---|---|
| static `annotateSwathMapsFromFile(filename, swath_maps, do_sort, force)` | `annotate_swath_maps_from_file(path, &mut Vec<SwathMap>, do_sort, force) -> Result<Vec<String>>` (the forced-override messages the source logs) |
| static `readSwathWindows(filename, lower, upper)` | `read_swath_windows(path) -> Result<(Vec<f64>, Vec<f64>)>` |
| — | native `parse_swath_windows(&str)`, `annotate_swath_maps(&[f64], &[f64], ..)`, `MAX_SWATH_WINDOW_FILE_BYTES` |

### `SwathQC.h` → `swath_qc.rs` (namespace `OpenSwath`)

| Source | Rust |
|---|---|
| typedef `ChargeDistribution = std::map<int,int>` | `ChargeDistribution = BTreeMap<i32, i32>` |
| deleted default ctor | none |
| `SwathQC(cd_spectra, decon_ms1_mz_tol)` | `SwathQC::new` |
| `getExpSettingsFunc()` | `exp_settings_func()` (closure over `apply_experimental_settings`) |
| `getSpectraProcessingFunc()` | `spectra_processing_func()` (closure over `process_spectrum`) |
| static `getChargeDistribution(swath_maps, nr_samples, mz_tol)` | `SwathQC::charge_distribution_of_maps` (distinct name: Rust has no overloading) |
| `storeJSON(filename)` | `store_json(path) -> Result<()>`; native `to_json() -> String` |
| `getChargeDistribution() const` | `charge_distribution()` |
| `setNrMS1Spectra(nr)` | `set_nr_ms1_spectra` |
| protected static `isSubsampledSpectrum_` | `is_subsampled_spectrum` (public; a pure function the class test reaches through a subclass) |
| members `cd_`, `nr_ms1_spectra_`, `cd_spectra_`, `decon_ms1_mz_tol_`, `ms1_spectra_seen_` | private fields |

### `SwathMapMassCorrection.h` → `swath_map_mass_correction.rs`

| Source | Rust |
|---|---|
| ctor (a `DefaultParamHandler`), dtor | `new`, `Default`; `defaults()`, `from_param`, `from_param_with_warnings`, `to_param` (the `peak_picking.rs` shape) |
| `updateMembers_` | `from_param` (private `from_complete_param`) |
| parameters `mz_extraction_window` … `debug_mz_file` | public fields of the same names |
| `correctMZ(transition_group_map, targeted_exp, swath_maps, pasef)` | `correct_mz(&mut self, &BTreeMap<String, &MRMTransitionGroupType>, &LightTargetedExperiment, &mut [SwathMap], pasef) -> Result<()>`; the source's `std::map` of group pointers is a `BTreeMap` of references, the maps are all wrapped or none |
| `correctIM(transition_group_map, targeted_exp, swath_maps, pasef, im_trafo)` | `correct_im(&mut self, groups, targeted_exp, &[SwathMap], pasef, &mut TransformationDescription) -> Result<()>`; serial where the source runs `#pragma omp parallel for` |
| `findSwathMapsPasef(transition_group, swath_maps)` | `find_swath_maps_pasef(&self, &MRMTransitionGroupType, &[SwathMap]) -> Result<Vec<SwathMap>>` |
| static `estimateWindow(residuals, quantile = 0.99, full_width = true, padding_factor = 1.0)` | `estimate_window(&[f64], f64, bool, f64) -> Result<f64>` |
| `get/setFragmentMzWindow`, `get/setFragmentImWindow`, `get/setPrecursorMzWindow`, `get/setPrecursorImWindow` | `fragment_mz_window`/`set_fragment_mz_window`, … (default -1) |
| file-local `findBestFeature`, `findSwathMaps` | private `find_best_feature`, `find_swath_maps` of `swath_map_mass_correction.rs` |
| — | `MRMTransitionGroupType`, `MZ_CORRECTION_FUNCTIONS`, `IM_CORRECTION_FUNCTIONS`, `SWATH_MAP_MASS_CORRECTION_NAME` |

## Preserved source conventions

- **Top-hat walk, both quirks** (`ChromatogramExtractorAlgorithm.cpp:63-119`). A target above every
  peak adds the last peak twice (current peak plus first leftward step); a leftward walk of two or
  more steps never adds the first peak of the spectrum. The class test's `400.05 -> 8400` encodes
  the second. Both are reproduced and oracle-checked (`quirk_end` = 114, `quirk_first` = 110,
  `single#1` = 6 on a one-peak spectrum).
- Open windows in the top-hat (`>`/`<`), **closed** windows in `PeakMapExtractor` (`>=`/`<=`).
- The walk position is shared by all coordinates of one spectrum; coordinates must be sorted by m/z.
  Descending targets silently give wrong sums (the class test documents this).
- `extractChromatograms` returns before any argument check when the input has no spectra; the
  Bartlett filter fails only when the first point would be extracted with it.
- `prepare_coordinates` appends and then stably sorts the **whole** coordinates vector; a compound
  without transitions gets one coordinate (m/z 0, ion mobility 0, whole RT range) and no isotopes;
  a non-positive `ms1_isotopes` adds none; NaN RT window copies the compound's `rt_start`/`rt_end`.
- `return_chromatogram`: base-peak type for MS1, SRM type for MS2, `peptide_sequence` from the
  compound's sequence or else its name, isolation offsets from the first precursor of the settings,
  drift time only for a non-negative coordinate IM and positive width.
- SWATH windows: header iff the first token is not a number; only two columns read; CRLF accepted;
  maps matched by position after an optional sort by `upper`; `-force` applies a widening row.
- `isSubsampledSpectrum_` arithmetic (`floor`/`ceil` steps, `lround`) reproduced; the whole
  13 × 13 × 14 grid matches the oracle.
- `SwathQC` deisotopes with the source's argument list (0.04 Th, charges 1–10, 3–10 peaks,
  keep only deisotoped, annotate charge) **and shares isotopes between envelopes**
  (`allow_shared_isotopes = true`), which the Rust `Deisotoper` makes opt-in; the charge
  distributions match the oracle for five sample counts and two tolerances.

## Native differences

- **Atomic results.** Every function either succeeds or leaves its outputs unchanged. The source
  appends partially before throwing (`extractChromatograms`, `prepare_coordinates`,
  `readSwathWindows`, `annotateSwathMapsFromFile`, `extractPeakMaps`).
- **Undefined behaviour becomes an error:** a transition referencing an unknown compound in
  `prepare_coordinates` (null dereference) and an unknown transition or compound id in
  `return_chromatogram` → `MissingInformation`; array length mismatches in the top-hat and extraction
  → `InvalidValue`; fewer coordinates than chromatograms → `InvalidValue`; an MS1 map without
  spectrum access in `SwathQC` → `MissingInformation`; an empty transition group in
  `findSwathMapsPasef` → `InvalidValue`.
- **`findSwathMapsPasef` precondition checked**: a first transition with ion mobility -1 is an
  error; the source's `OPENMS_PRECONDITION` is debug-only.
- **SWATH-window rows are parsed strictly.** The source's two unchecked `operator>>` leave both values
  unchanged on a blank row, so a blank line **repeats the previous window** (oracle:
  `400 425\n425 450\n\n` reads three windows, the last duplicated) and an empty or header-only file
  reads uninitialised values; a one-column row keeps the previous upper bound. Here these are
  `Error::Parse`. A non-numeric column is `Parse` here where the source reports "illegal ranges".
- **Non-finite SWATH bounds** are refused (`InvalidValue`); NaN passes the source's `!(lower < upper)`
  test only as a failure, infinity is accepted there.
- **Stable sort** of SWATH maps by `upper` (the source's `std::sort` is unstable); NaN `upper` with
  `do_sort` is refused (undefined in the source's sort).
- `nr_ms1_spectra` meta value: negative or out of `i32` → `InvalidValue` (source truncates to `int`
  and wraps to `size_t`).
- `return_chromatogram` marks **copies** of the settings' data-processing entries with
  `performed_on_spectra`; the source mutates the caller's shared objects.
- `return_chromatogram` takes the settings as an `MSSpectrum` (the source's `SpectrumSettings`,
  which `MSSpectrum` derives from); `SpectrumAccessOpenMS::spectra_meta_info` returns one.
- Logging: the source's `OPENMS_LOG_*` output is not reproduced; forced SWATH overrides are
  returned as strings. The unreachable "Drift time of ion is negative" warning of
  `extractChromatograms` (the branch requires a non-negative drift time) has no counterpart.
- Ceilings: `MAX_EXTRACTED_POINTS` (2^30 spectrum×coordinate pairs), `MAX_PEAK_MAP_POINTS` (2^30),
  `MAX_COORDINATES` (2^28), `MAX_SWATH_WINDOW_FILE_BYTES` (16 MiB), all preflighted or checked before
  results are committed.
- The extraction is serial, as the source's is (no OpenMP in these files except `correctMZ`/`correctIM`).

## C++ defect candidates (unnumbered)

1. `ChromatogramExtractorAlgorithm::extract_value_tophat` (both overloads) counts the last peak twice
   when the target m/z lies above every peak and the last peak is inside the window
   (`ChromatogramExtractorAlgorithm.cpp:63-103`, `:179-222`; the IM overload repeats the walk). Oracle: `[100, 200, 200.02]` at
   200.05 ± 0.1 gives 114 (7 counted twice); a one-peak spectrum gives twice its intensity.
2. Same walk: the leftward loop (`:98`, `:217`) stops *on* `mz_start` without adding it, so the
   first peak is dropped whenever the walk reaches it in two or more steps. Oracle:
   `[100, 100.05, 100.1]` at 100.1 ± 0.2 gives 110, not 111.
3. `SwathWindowLoader::readSwathWindows` never checks its stream extractions: a blank line repeats
   the previous window, a one-column row reuses the previous upper bound, an empty or header-only
   file reads uninitialised `lower`/`upper` (`SwathWindowLoader.cpp:118-131`).
4. `ChromatogramExtractor::prepare_coordinates` dereferences `tr2pep[ref]` (a default-inserted null)
   for a transition whose peptide reference is not a compound (`ChromatogramExtractor.cpp:123`);
   `return_chromatogram` does the same for `trans_map[coord.id]` (`ChromatogramExtractor.h:212`) and
   `getCompoundByRef`.
5. `return_chromatogram` mutates the caller's `DataProcessing` objects through the copied
   `SpectrumSettings` (shared pointers), adding `performed_on_spectra` to the input's history.

## Checked boundaries and evidence

- **Tier 1** (executed C++ Release, ibminode06, x86_64 Linux, gcc 14.4.0 `-O2 -DNDEBUG`, install
  `openms4-release-bc9cc12-c19e494-174b576`; two runs byte-identical): `oracle_tophat` (34 calls),
  `oracle_extract_chromatograms` (3 × 59 points from ChromatogramExtractor_input.mzML, 3 IM windows),
  `oracle_peak_maps` (4 window combinations), `oracle_prepare_coordinates` (7 cases, and the
  TraML-converted light experiment), `oracle_return_chromatogram` (MS2 and MS1),
  `oracle_swath_windows` (9 files, 5 annotation runs), `oracle_swath_qc_subsample_grid`,
  `oracle_swath_qc_charge_distribution` (10 runs + JSON), `oracle_estimate_window` (96 values),
  `oracle_swath_map_mass_correction_defaults`. All compared bit for bit.
- **Tier 3**: every section of the six class tests at the pin, and of
  `SwathMapMassCorrection_test.cpp` (`correctMZ`: none, unweighted, quadratic, weighted quadratic,
  both delta-ppm functions and the PASEF flag; `correctIM`: MS2 with and without the MS1 map, a
  single peptide, PASEF, MS1 without MS2 (`UnableToFit`), MS1 calibration), in
  `swath_map_mass_correction_correct_mz_class_test` and
  `swath_map_mass_correction_correct_im_class_test` of `tests/openswath_extraction.rs`. Two native
  guards: a group with features but no transition is `InvalidValue`, and the debug files are
  created before the loop (a file that cannot be created is an error; the source's `ofstream`
  silently writes nothing).
- **Tier 4**: atomicity, argument order of checks, ceilings, guards.

## Gaps (named dependencies)

- None for `SwathMapMassCorrection.h`: `correctMZ` / `correctIM` were ported at group O integration
  once `OpenSwathScoring::fetchSpectrumSwath`, `DIAHelpers`, `LinearRegression` and
  `QuadraticRegression` were in the tree.
- `ChromatogramExtractor::return_chromatogram<TargetedExperiment>` needs
  `ANALYSIS/TARGETED/TargetedExperiment.h`.
- The class tests load `ChromatogramExtractor_input.TraML` through `TraMLFile` and
  `OpenSwathDataAccessHelper::convertTargetedExp`; neither is ported, so the tests transcribe the
  converted experiment (checked against the oracle's print of it).

# mzML spectrum mobility and FileHandler loading options

Work package A3-FORMAT-IO of the early TOPP bundle
(`docs/EARLY_TOPP_BUILD_PLAN.md`). It closes four narrow routes that
FileInfo and FeatureFinderCentroided depend on: spectrum- and scan-level ion
mobility including the FAIMS compensation voltage, the `MS:1000525`
representation reset, units on ion-mobility data arrays, and FileHandler type
detection plus option-taking loaders. It does **not** close `MzMLFile.h`,
`MzMLHandler.h` or `FileHandler.h`.

Pins: core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`, TOPP `174b576`,
test data `0cb15f2`. Oracle: product-sdk (Debug, core `4fdec46`), whose traced
sources are unchanged against the core pin.

| File | Role |
|---|---|
| `src/format/mzml.rs` | Reader routes `Record::spectrum_mobility`, the `MS:1000525` arm and `Binary::mobility_array_unit`; writer preflight `validate_spectrum_mobility` |
| `src/format/mzml_precursor.rs` | `mobility_term` and `mobility_cv`, the accession/unit tables of the existing selected-ion route, extracted unchanged and now shared |
| `src/format/mzml_settings.rs` | `write_scan` emits the spectrum mobility term |
| `src/format/file_handler.rs` | `FileHandler::get_type`, `load_experiment_with_options`, `load_feature_map_with_options` |
| `tests/mzml_mobility.rs`, `tests/file_handler_type_detection.rs` | Integration tests |
| `tests/data/mzml_mobility/` | Upstream and synthetic inputs, oracle output |
| `../oracle/a3-format-io/` | C++ driver, case generator, run script, hashed manifest |

## API mapping

### MzMLHandler reader routes

| Source route (`MzMLHandler.cpp`) | Rust | Status |
|---|---|---|
| spectrum `MS:1000127`, `MS:1000128` (1634-1641) | `MSSpectrum::spectrum_type` | existing |
| spectrum `MS:1000525` resets the type to `UNKNOWN` (1642-1645) | `SpectrumType::Unknown` | ported |
| spectrum `MS:1003441` sets `IM_CENTROIDED` (1646-1649) | none | not ported: `MSSpectrum` has no ion-mobility peak type |
| spectrum `MS:1001581`, FAIMS compensation voltage (1731-1738) | `drift_time`, `DriftTimeUnit::FaimsCompensationVoltage` | ported |
| spectrum `MS:1002476`, `MS:1002815`, `MS:1002954` | none | not a source route: warned and ignored, as here |
| scan `MS:1002476`, `MS:1002815`, `MS:1001581`, `MS:1002954` (2279-2311) | `drift_time` and `drift_time_unit` | ported |
| selectedIon mobility on the precursor (1873, 1876) | `Precursor::drift_time`, `drift_time_unit` | existing |
| selectedIon mobility copied to the spectrum (1874-1875) | none | not ported, see native differences |
| binaryDataArray unit stored as `unit_accession` (`MzMLHandlerHelper.cpp:291-295`) | `DataArray::metadata["unit_accession"]` | ported for ion-mobility arrays only |

### MzMLHandler writer routes

| Source route | Rust | Status |
|---|---|---|
| scan mobility term after the start time (5412-5440) | `mzml_settings::write_scan` | ported, with the preflight below |
| array `unit_accession` written on the array cvParam (5781-5786) | written as a `unit_accession` userParam | not ported; decoded content is identical |
| type `UNKNOWN` written as `MS:1000525` (5285-5288) | nothing written | existing; reads back as `Unknown` either way |
| precursor mobility (4607-4628) | `mzml_precursor::write_end` | existing |

### FileHandler.h

| Source member | Rust | Status |
|---|---|---|
| `getType` | `FileHandler::get_type` | ported |
| `getTypeByFileName` | `file_types::type_by_file_name` | existing |
| `hasValidExtension` | `file_types::has_valid_extension` | existing |
| `stripExtension`, `swapExtension` | `file_types::strip_extension`, `swap_extension` | existing |
| `getConsistentOutputfileType` | `file_types::consistent_output_type` | existing |
| `getTypeByContent(filename)` | `file_handler::type_by_content(&[u8])`; the path form runs inside `get_type` | existing |
| `isSupported` | none; `can_read_experiment` and siblings report native availability | not ported |
| `getOptions`, `setOptions` | the `&PeakFileOptions` argument of `load_experiment_with_options` | ported as an argument |
| `getFeatOptions`, `setFeatOptions` | the `&FeatureFileOptions` argument of `load_feature_map_with_options` | ported as an argument |
| `loadExperiment` | `load_experiment`, `load_experiment_with_options` | partial: DTA, DTA2D, MGF, MS2, mzML; no log type, `rewrite_source_file` or `compute_hash` |
| `storeExperiment` | `store_experiment` | existing |
| `loadSpectrum`, `storeSpectrum` | none | not ported |
| `loadFeatures` | `load_feature_map`, `load_feature_map_with_options` | partial: featureXML only |
| `storeFeatures` | `store_feature_map` | existing |
| `loadConsensusFeatures`, `storeConsensusFeatures` | `load_consensus_map`, `store_consensus_map` | existing |
| `loadIdentifications`, `storeIdentifications` | `load_identifications`, `store_identifications` | existing, idXML only |
| `loadTransitions`, `storeTransitions`, `loadTransformations`, `storeTransformations`, `storeQC` | none | not ported |
| `computeFileHash` | none | not ported |
| `applyPostLoadOptions_` (protected) | none | not needed: only the TDF branch uses it |

`PeakFileOptions.h` is the existing `format::peak_options::PeakFileOptions`,
unchanged. `load_experiment_with_options` forwards it the way
`FileHandler::loadExperiment` does:

| Detected type | Options the source reader consumes | Rust |
|---|---|---|
| mzML | all | `mzml::read_with_load_options`, default XML and binary limits |
| DTA2D | RT, m/z and intensity ranges (`DTA2DFile.h:208-243`) | `dta2d::ReadOptions` ranges |
| DTA, MGF, MS2 | none | options ignored, as in the source |
| anything else | | `Error::Unsupported` |

## Preserved source conventions

- All four mobility accessions map to the same units on every route, through
  one table. Only `MS:1001581` is read directly below `<spectrum>`; the other
  three are ignored there, as the source warns and ignores them.
- `MS:1000525` after `MS:1000127` or `MS:1000128` leaves the type unknown;
  before them it is overridden. Oracle rows `spectrum_representation`.
- A FAIMS voltage of `-1` equals the unset sentinel. It is stored as `-1`,
  `determine_im_format` reports `None`, and the writer still writes it because
  the unit is FAIMS, exactly as the source checks the unit first.
- An ion-mobility array keeps its unit as `unit_accession` string metadata,
  the source representation. Whether an array is ion mobility follows
  `MSSpectrum::contains_im_data` (`IMDataArrayUtils::getIMUnit`), so CV
  children of `MS:1002893` and the `Ion Mobility` user-parameter names qualify;
  the parent term `ion mobility array` does not.
- `get_type` strips trailing `/` and `\`, trusts a known extension without
  touching the file, maps every Bruker TDF name (`.d`, `.d/`, `.d.zip`) to
  unknown as a build without OpenTIMS does, and otherwise recognises content.
- Ranges are half-open (`DRange::encloses`): the minimum is kept and the
  maximum dropped. DTA2D drops a spectrum the ranges leave empty.
- FeatureFinderCentroided's intensity range is **executed** as `[0, DBL_MAX)`:
  `std::numeric_limits<DPosition<1>>::min()` has no specialisation and is
  `DPosition()`. The oracle keeps `0.0`, a subnormal and `DBL_MIN` and drops
  `-1`. The source comment intends `[DBL_MIN, DBL_MAX)`; both give 112 spectra
  and 3084 peaks on `FeatureFinderCentroided_1_input.mzML`, which has no zero or
  negative MS1 peak. A faithful tool port must pass `0.0`.

## Native differences

1. **Selected-ion mobility is not copied to the spectrum.**
   `MzMLFile_test.cpp:418-421` expects spectrum 1 of `MzMLFile_1.mzML` to carry
   its precursor's drift time 8.1 ms, and the oracle shows it. The native reader
   leaves the spectrum unset, because the existing
   `tests/precursor_workflow.rs:157,188` asserts that a precursor-only mobility
   round-trips without a spectrum value. Adopting the source behaviour flips that
   assertion and needs an integrator decision; the new test pins the current
   state.
2. **Conflicting repeats are refused.** The source keeps the last of several
   mobility terms on one spectrum (oracle: 8.5 ms after 7.5 ms; -60 V after a
   spectrum-level -50 V). The native reader returns `Error::Unsupported`
   for a differing value or unit and accepts an identical repeat, matching its
   refusal of repeated scan start times.
3. **Mobility unit attributes are checked.** The source ignores them (oracle:
   7.5 read as milliseconds despite `UO:0000010`). The native reader requires the
   quantity's unit, as its selected-ion route already did, and additionally
   reads `UO:000218` as volts for the FAIMS voltage: the upstream
   `FAIMS_CV-60C_V-45_Interleaved.mzML:324` and `FAIMS_test_data.mzML:171` carry
   that spelling from an earlier OpenMS writer, while the pinned writer emits
   `UO:0000218`.
4. **Units stay unsupported on other auxiliary arrays.** The source keeps a
   `unit_accession` on every non-default array (oracle:
   `im_arrays_non_mobility_unit`); the native reader still refuses them, the
   contract of `tests/mzml_auxiliary_review.rs:324-336`.
5. **Array units are written as a userParam.** The source writes them on the
   array cvParam. Both forms decode to the same `unit_accession` metadata in
   either implementation; moving the unit onto the cvParam needs the header
   array-parameter writer (`src/format/mzml_header/write.rs`) to skip the key.
6. **Lossy spectrum mobility is refused on write.** The source writes a drift
   time without a unit as milliseconds with a warning, and drops a non-FAIMS unit
   without a drift time. The native writer rejects both, and a non-finite drift
   time, before emitting any byte.
7. **Loader errors differ in kind.** A detected type outside the allowed list is
   `Error::InvalidValue` (source `ParseError`); a missing unknown-extension file
   in `get_type` is the I/O error (source `FileNotFound`).
   `load_feature_map_with_options` has no allowed-type list (featureXML is the
   only native reader, where the source throws `InvalidFileType` for a
   disallowed type). `get_type` reads at most 64 KiB, and returns
   `Error::Unsupported` for ZIP content and for gzip/bzip2 content without the
   `file-compression` feature, where the source looks inside.
8. **No SRM conversion after mzML loading.** `loadExperiment` moves SRM spectra
   into chromatograms (`ChromatogramTools::convertSpectraToChromatograms`);
   neither native loader does.
9. **NaN never reaches a range filter.** The source keeps a NaN intensity, since
   both `encloses` comparisons fail; the native readers reject non-finite values
   first.

## Checked boundaries and evidence

Tier 1, oracle-generated executed differential
(`tests/data/mzml_mobility/a3_format_io_oracle.tsv`, sha256 `1a7009ed…772a`,
reproduced byte for byte in two runs; `../oracle/a3-format-io/manifest.json`):

| Oracle records | Rust test |
|---|---|
| `FAIMS_CV-60C_V-45_Interleaved`: 12 spectra, -45/-60 V, one chromatogram | `upstream_faims_interleaved_file_reads_spectrum_level_compensation_voltages` |
| `FAIMS_test_data`: -65 V below `<scan>` on MS1 and MS2 | `upstream_faims_test_data_reads_the_scan_level_compensation_voltage` |
| `MzMLFile_1`: 7.1 ms scan value; precursor 8.1 ms and unset | `class_test_fixture_scan_drift_time_and_the_precursor_propagation_gap` |
| `scan_mobility`: 13 routes, units, sentinel, repeat, second scan, legacy unit | `scan_and_spectrum_level_mobility_terms_match_the_oracle` |
| `spectrum_representation`: five orderings | `spectrum_representation_reset_matches_the_oracle` |
| `im_arrays`: five IM arrays with and without units, all `im_peak` | `mobility_arrays_with_and_without_units_are_per_peak` |
| `FeatureFinderCentroided_1_input`: 112 spectra, unknown type | `feature_finder_input_matches_the_oracle` |
| `scan_mobility_conflict`, `scan_mobility_unit_mismatch`, `im_arrays_non_mobility_unit` | `source_lenient_inputs_stay_explicit_native_errors` |
| `type`: 13 names, directories, content, gzip, missing file | `get_type_matches_the_oracle_for_names_directories_and_content` |
| `filtered`: FeatureFinderCentroided options on the FFC_1 input, per spectrum | `feature_finder_options_load_the_ffc_1_input` |
| `spectrum_peak`, `filtered_detail_peak`: intensity edges | `intensity_range_edges_match_the_oracle` |
| `features`: FFC_1 output 8/30/0 and 8/0/0; `featurexml_source_1` 2/1/2 and 2/0/0 | `feature_file_options_match_the_oracle` |

Every compared value is bitwise: drift times and peaks through the driver's
`%a` output, counts and names exactly.

Tier 3, source review: `MzMLFile_test.cpp:368-369, 451-452, 464-465` literals;
DTA2D range routing (`DTA2DFile.h:208-243`) and DTA/MGF option neglect
(`FileHandler.cpp:869-936`) in `dta2d_receives_only_the_three_ranges_the_source_reads`
and `dta_and_mgf_ignore_peak_file_options_like_the_source`.

Tier 4, native invariants: writer round trips of all four units, the FAIMS
sentinel and a spectrum without retention time, plain and zlib; round trips of
`scan_mobility`, `im_arrays`, `spectrum_representation`, `FAIMS_test_data` and
`MzMLFile_1`; refusal of lossy writes with the output untouched; unit, value and
repeat validation on the reader.

Regression: the 69 existing test targets that load or write mzML or featureXML,
use FileHandler or ion-mobility types passed before the change, and together
with the two new targets (71 binaries, 1177 tests) after it, on Rust stable and
1.85.0 (`--all-features`, kim). The new targets also pass with
`--no-default-features --features mzml` (19 tests) and
`--features mzml,featurexml` (file-handler target, 8 tests) on both toolchains
where run; clippy `-D warnings`, `cargo doc -D warnings` and a
`--no-default-features` check of the library and tests are clean.

## C++ issue candidates

- FeatureFinderCentroided's intensity filter keeps zero intensities
  (`FeatureFinderCentroided.cpp:182-184`): `numeric_limits<DPosition<1>>::min()`
  is zero; `DPosition<1>::minPositive()` was meant. Executed.
- A FAIMS voltage of -1 V is indistinguishable from an unset drift time
  (`IMTypes::DRIFTTIME_NOT_SET`): reported as no ion mobility. Executed.
- The precursor writer has no FAIMS case (`MzMLHandler.cpp:4607-4628`): a
  precursor FAIMS voltage falls into `default`, is written as milliseconds with
  a warning and reads back with the wrong unit. Source review only.
- Upstream FAIMS fixtures spell the volt unit `UO:000218`. Fixture defect.

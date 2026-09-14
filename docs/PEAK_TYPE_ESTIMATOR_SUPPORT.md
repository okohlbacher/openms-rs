# Peak type estimator support

Native coverage of `FORMAT/PeakTypeEstimator.h` at Core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The header is header-only: one
static template, `estimateType(begin, end)`, plus a commented-out alternative
kept "for reference".

| Artifact | Path |
| --- | --- |
| Implementation | `src/format/peak_type_estimator.rs` (public entry point); `src/kernel/spectrum_type.rs` (`estimate`, the shared shoulder transcription, unchanged) |
| Tests | `tests/peak_type_estimator.rs` |
| Fixtures | `tests/data/peak_type_estimator/oracle_estimates.tsv`; the three upstream DTA spectra are read from `tests/data/spectrum_type/` (identical sha256 to the pinned class-test data) |
| Manifest | `tests/data/peak_type_estimator_provenance.json` |
| Oracle driver | `../oracle/pte-faims-helper/` (outside this repository) |

The estimator algorithm was already transcribed line by line in
`src/kernel/spectrum_type.rs`, where `MSSpectrum::get_type(true)` uses it, but
only as a crate-private function. The only public classifier was
`processing::peak_picking::estimate_spectrum_type`, which validates the whole
spectrum and refuses negative intensities and unsorted or duplicate m/z values
where the source classifies them. FileInfo (`FileInfo.cpp:1599`) and
PeakPickerHiRes call the source estimator on raw spectrum peaks, so they need a
public entry point with source semantics. This work package adds that entry
point and delegates to the existing transcription; it does not add a second
copy of the algorithm.

## API mapping

| C++ member | Rust |
| --- | --- |
| `class PeakTypeEstimator` | `format::peak_type_estimator::PeakTypeEstimator`, a zero-sized unit struct: the source class holds no state |
| `PeakTypeEstimator()`, `~PeakTypeEstimator()` (implicit; exercised by the two `[EXTRA]` sections) | derived `Default`; no drop glue |
| `template<typename PeakConstIterator> static SpectrumSettings::SpectrumType estimateType(const PeakConstIterator& begin, const PeakConstIterator& end)` | `PeakTypeEstimator::estimate_type(&[Peak1D]) -> Result<SpectrumType>` with the default `SpectrumTypeQueryLimits`, and `PeakTypeEstimator::estimate_type_with_limits(&[Peak1D], SpectrumTypeQueryLimits)`. An iterator sub-range is a sub-slice |
| the literal `5` of the `@note` and of `end - begin < 5` | `PeakTypeEstimator::MIN_PEAKS` |
| `SpectrumSettings::SpectrumType::{UNKNOWN, CENTROID, PROFILE}` | `kernel::SpectrumType::{Unknown, Centroid, Profile}` |
| commented-out `estimateType` (quartile-of-distances variant, `PeakTypeEstimator.h:158-242`) | not ported: it is inside a comment and never compiled; the header says it "does not work reliably" |

The template accepts any peak type with `getMZ`/`getIntensity`/`setIntensity`.
Every caller in the pinned core (`MSSpectrum.cpp:162`, `FileInfo.cpp:1599`)
passes `MSSpectrum` iterators, that is `Peak1D`, so the Rust entry point takes a
`Peak1D` slice. No other instantiation exists to port.

## Preserved source conventions

- **Fewer than five peaks are `Unknown`**, before any value is looked at, so a
  short range containing NaN is still `Unknown` (`PeakTypeEstimator.h:40`, `:47`).
- **Up to five maxima, stopping early once more than 50% of the total intensity
  is explained.** The explained intensity is accumulated as the source does:
  the left shoulder adds the maximum once, it is subtracted again, and the right
  shoulder adds it back.
- **Strict shoulder conditions.** A point joins a shoulder while it is at most
  the previous point's intensity, positive, *strictly* more than 10% of the
  previous point (`10/100` does not qualify, oracle `ratio_exact_tenth`), and
  *strictly* less than 1 Th from the maximum (`mz + 1 > max_mz`, which fails at
  m/z 1e17 where `mz + 1 == mz`, oracle `huge_mz_profile_shape`).
- **The begin element is never scanned** by the left shoulder loop
  (`it != data.begin()` ends the loop there), and the sink-restoring rule
  `(it+1)->setIntensity(int_last)` runs on whatever point the loop stopped at.
- **Fewer than two points on either shoulder is centroid evidence**; at least
  two on both is profile evidence.
- **The evidence ratio is `float`**: `profile / float(profile + centroid)`, and
  `> 0.75` is Profile. With no evidence the ratio is `0/0 = NaN`, not greater
  than 0.75, so a spectrum without any positive intensity is `Centroid`
  (oracle `all_zero_five`).
- **No sortedness, distinctness or sign requirement.** Descending m/z,
  duplicate m/z and negative intensities are classified (oracle
  `descending_mz_profile`, `duplicate_mz_profile`, `negative_edges_profile`;
  all `Profile`).
- **The input is not modified.** The source copies the range; the Rust slice is
  borrowed immutably and the scratch copy is internal.
- **Precision.** The source keeps `float` intensities and promotes them to
  `double` in every comparison, division and sum. The Rust transcription holds
  the same values in `f64`. Every value written back into the scratch copy is
  zero or an intensity previously read, so the scratch values remain exactly
  the input `f32` values and every decision is the source's.

## Native differences

- **Non-finite values are refused.** Once at least five peaks are present, a NaN
  or infinite m/z or intensity returns `Error::InvalidValue`. The source
  classifies such input as `Centroid` (oracle `nan_intensity`: a NaN is never a
  maximum and makes the total NaN; `inf_intensity`: every shoulder ratio against
  an infinite maximum is NaN). A class derived from meaningless comparisons is
  not reported. `MSSpectrum::get_type(true)` refuses the same input, so the two
  public paths agree.
- **Resource ceilings.** The ceilings and their order are those
  `MSSpectrum::get_type_with_limits` applies to the peaks when it reaches
  estimation: `max_points`, then 32 work units per peak against `max_work`,
  then 16 bytes per peak against `max_bytes`, then the finiteness check, then a
  fallible allocation. Short input returns `Unknown` before any ceiling is
  consulted. For a spectrum with no stored type and no data-processing records,
  queried with `query_data` set, the test asserts both paths accept and refuse
  at the same boundaries. The source has no ceiling.
- **Data-processing records spend work the slice does not.** Before the peaks
  are charged, `get_type_with_limits` searches the spectrum's data-processing
  records for a peak-picking step. It charges each record `1 + 12 * h` work
  units against the same `max_work`, where `h` is the bit length of the
  record's action count, so a record without actions costs one unit
  (`get_type_with_budget` in `src/kernel/spectrum_type.rs`). A spectrum with
  records can therefore fail under limits its peak slice fits. The 7-peak
  `ascending_profile` shape at `max_work` 224 (7 × 32) is `Profile` as a slice
  and a resource error as a spectrum carrying one default record; at 225 both
  are `Profile` (`a_data_processing_record_spends_work_the_peak_slice_does_not`).
  A record with a peak-picking action instead makes the spectrum `Centroid`
  without estimating. `estimate_type_with_limits` receives no records and
  charges none.
- **`Result` return type.** The source cannot fail. The Rust function returns
  `Result<SpectrumType>` for the two refusals above; the early-TOPP-bundle plan
  sketched `-> SpectrumType`, which cannot express them.

## Checked boundaries and evidence

| Evidence | Tier | What it covers |
| --- | --- | --- |
| `PeakTypeEstimator_test.cpp:33-57` literals | 3 (source review) | the two `[EXTRA]` sections; `PeakTypeEstimator_raw.dta` and `_rawTOF.dta` are `PROFILE`, `_peak.dta` is `CENTROID`, and each is `UNKNOWN` after `resize(4)` |
| `oracle_estimates.tsv`, `pte_dta` / `pte_dta_first4` | 1 (executed differential) | the same three spectra through the unmodified C++ `DTAFile` and `estimateType`: sizes 66, 99386 and 121, and the same classes |
| `oracle_estimates.tsv`, `pte_mzml` | 1 (executed differential) | `estimateType` on every spectrum (65) of eight mzML fixtures, with the native ID, MS level and peak count; the Rust estimate on the Rust-loaded peaks matches for all 65, and `MSSpectrum::get_type(true)` with nothing stored equals the estimate for all 65 |
| `oracle_estimates.tsv`, `pte_synthetic` | 1 (executed differential) | eleven shapes: short input, flat, all-zero, profile in ascending and descending m/z order, negative edges, duplicate m/z, the exact 10% ratio, m/z 1e17, NaN and infinite intensity |
| `tests/peak_type_estimator.rs` native tests | 4 | finiteness after the short-input gate; ceiling parity with `get_type_with_limits` for a spectrum without stored type or data-processing records; the one-unit charge of a default data-processing record that makes the spectrum fail at `max_work` 224 where its 7-peak slice succeeds; the stricter picker refusing what the estimator classifies |
| `tests/spectrum_type.rs` (unchanged, earlier package) | 3 | the transcription's own boundary tests through `MSSpectrum::get_type` |

The stored and queried types (`getType(false)`, `getType(true)`) recorded in the
same oracle rows are reader evidence rather than estimator evidence.
`stored_and_queried_types_match_the_executed_oracle` asserts them too, and the
Rust reader matches C++ on all 65 spectra. None of the eight fixtures places
`MS:1000525` inside a spectrum after `MS:1000128`, so the C++ reset to `UNKNOWN`
(`MzMLHandler.cpp:1642-1645`), which A3-FORMAT-IO ports, is not exercised here.

The oracle is product-sdk (Debug, core `4fdec46`), accepted as a
development-time oracle. `git diff 4fdec46 bc9cc12` is empty for
`PeakTypeEstimator.h`, `PeakTypeEstimator_test.cpp` and the DTA fixtures, and the
installed `PeakTypeEstimator.h` hashes equal to the pin.

## Class-test section accounting

| Section | Rust test |
| --- | --- |
| `[EXTRA]PeakTypeEstimator()` | `extra_constructor_and_destructor_sections` |
| `[EXTRA] ~PeakTypeEstimator()` | `extra_constructor_and_destructor_sections` |
| `estimateType(begin, end)` | `estimate_type_section_classifies_the_three_upstream_spectra`, plus the oracle comparison `the_upstream_spectra_match_the_executed_oracle_record_for_record` |

All three sections pass unchanged.

## Source defects observed

None in the estimator. The header comment above `if (evidence_ratio > 0.75)`
says "80% are profile" while the threshold is 75%; the code is authoritative
and the Rust documentation states 0.75.

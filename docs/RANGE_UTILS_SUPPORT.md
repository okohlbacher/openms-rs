# RangeUtils

The native `kernel::range_utils` module ports `OpenMS/KERNEL/RangeUtils.h` at
Core SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The header is header-only:
fifteen functor templates and no `.cpp`, so the inline bodies are the
implementation. See the [implementation](../src/kernel/range_utils.rs),
[tests](../tests/range_utils.rs) and
[hashed provenance](../tests/data/range_utils_provenance.json). No C++ was
built or executed; every transcribed expectation is tier 3 (source review) and
every native check is tier 4.

The source functors are consumed by `FileFilter`, `FeatureFinderCentroided`,
`FeatureFinderMultiplex`, `AssayGeneratorMetabo` and
`AssayGeneratorMetaboSirius` (`docs/core-sdk-coverage.json`), always through
`std::remove_if`/`std::erase_if`.

## API mapping

Every source functor exposes a constructor and `operator()`; the protected
fields (`min_`, `max_`, `reverse_`, …) have no accessors and are private here.
`operator()` is `SpectrumPredicate::keep` or `PeakPredicate::keep`. The
`reverse` constructor flag is preserved on every predicate.

| Source | Native | Notes |
| --- | --- | --- |
| `@defgroup RangeUtils` examples (`erase_if` with `InRTRange`, `InIntensityRange`) | Two module-level doctests | Rewritten to `retain_spectra` / `retain_peaks_where`; removal is `reverse = true` or a negated closure |
| `HasMetaValue<MetaContainer>(metavalue, reverse)`, `operator()` | `HasMetaValue::new(key, reverse)`, `keep`, `evaluate(&MetaInfo)` | Source is generic over any `MetaInfoInterface`; `keep` serves `MSSpectrum`, `evaluate` any `MetaInfo` owner |
| `InRTRange<SpectrumType>(min, max, reverse)`, `operator()` | **Mapped**: `MSExperiment::spectra_in_rt_range(min, max, 0)` for the closed-range query (sorted RT required, borrowing); removal is `retain_spectra(&\|s\| !(min <= s.rt && s.rt <= max))` | No dedicated struct, to avoid a second RT-range API; `spectra_in_rt_range` also rejects `min > max` with `InvalidValue` |
| `InMSLevelRange<SpectrumType>(levels, reverse)`, `operator()` | **Mapped**: `MSExperiment::ms_levels()` / `contains_scan_of_level(level)` for membership; removal is `retain_spectra(&\|s\| !levels.contains(&s.ms_level))` | No dedicated struct; `ms_level` is `u32` so the source `IntList` negative levels cannot be expressed |
| `HasScanMode<SpectrumType>(Int mode, reverse)`, `operator()` | `HasScanMode::new(ScanMode, reverse)`, `keep` | Typed enum instead of the source cast integer |
| `HasScanPolarity<SpectrumType>(IonSource::Polarity, reverse)`, `operator()` | `HasScanPolarity::new(Polarity, reverse)`, `keep` | `POLNULL` is `Polarity::Unknown` |
| `IsEmptySpectrum<SpectrumType>(reverse)`, `operator()` | `IsEmptySpectrum::new(reverse)`, `keep` | Peak list only, as the source `empty()` |
| `IsZoomSpectrum<SpectrumType>(reverse)`, `operator()` | `IsZoomSpectrum::new(reverse)`, `keep` | Reads `instrument_settings.zoom_scan` |
| `HasActivationMethod<SpectrumType>(StringList methods, reverse)`, `operator()` | `HasActivationMethod::new(methods, reverse)`, `from_names(names, reverse)`, `keep` | Typed `ActivationMethod` set; `from_names` accepts the source long names and the short names, and rejects unknown names (source: silently never match) |
| `InPrecursorMZRange<SpectrumType>(mz_left, mz_right, reverse)`, `operator()` | `InPrecursorMZRange::new(mz_left, mz_right, reverse) -> Result`, `keep` | All precursors must be inside; none is vacuously inside |
| `HasPrecursorCharge<SpectrumType>(IntList charges, reverse)`, `operator()` | `HasPrecursorCharge::new(charges, reverse)`, `keep` | Any precursor with a listed charge matches |
| `InMzRange<PeakType>(min, max, reverse)`, `operator()` | `InMzRange::new(min, max, reverse) -> Result`, `keep(&Peak1D)`, `contains(mz)` | `contains` serves other peak types |
| `InIntensityRange<PeakType>(min, max, reverse)`, `operator()` | `InIntensityRange::new(min, max, reverse) -> Result`, `keep(&Peak1D)`, `contains(f32)` | `f32` widened to `f64` before comparing, as the source |
| `IsInCollisionEnergyRange<SpectrumType>(min, max, reverse)`, `operator()` | `IsInCollisionEnergyRange::new(min, max, reverse) -> Result`, `keep`; helper `collision_energy(&Precursor)` | Reads `COLLISION_ENERGY_KEY` in `precursor.cv_terms.metadata`, then the `MS:1000045` term |
| `IsInIsolationWindowSizeRange<SpectrumType>(min_size, max_size, reverse)`, `operator()` | `IsInIsolationWindowSizeRange::new(min_size, max_size, reverse) -> Result`, `keep` | Width = lower + upper offset |
| `IsInIsolationWindow<SpectrumType>(vector<double> vec_mz, reverse)`, `operator()` | `IsInIsolationWindow::new(mz_values, reverse) -> Result`, `mz_values()`, `keep` | Sorted on construction, as the source |
| `std::remove_if` / `std::erase_if` over spectra | `MSExperiment::retain_spectra(&P) -> Result<usize>`, `retain_spectra_with_limits` | Keeps matches; returns the removed count |
| `std::remove_if` / `std::erase_if` over peaks | `MSSpectrum::retain_peaks_where(&P) -> Result<usize>`, `retain_peaks_where_with_limits` | Drives the existing `retain_peaks`, so aligned data arrays move with the peaks |
| — | `SpectrumPredicate`, `PeakPredicate` traits; blanket impls for `Fn(&MSSpectrum) -> bool` / `Fn(&Peak1D) -> bool` | Native: closures are predicates |
| — | `RangeFilterLimits { max_items }` (default 50,000,000) | Native preflight ceiling |
| — | `COLLISION_ENERGY_KEY`, `COLLISION_ENERGY_ACCESSION` | Native constants naming the source storage |

Doxygen carried: every `@param` (all constructor parameters, including each
`reverse` sentence) is on the constructor; the `@note` on `InMzRange` (m/z is
dimension 0) is neutralised by `Peak1D::mz`; the three `@note`s claiming the
MSn predicates "return true" for MS1 spectra, and the collision-energy note
about spectra without a collision energy, are carried with the correction that
the code returns `false` (see below). The `@ingroup` tags are dropped; the
`@see`-free header has no unresolved link targets. No `@exception` exists in
the source.

## Preserved source conventions

- All ranges are closed intervals: `min <= x && x <= max` for RT, m/z,
  intensity and precursor m/z (`RangeUtils.h:379`); `!(x > max || x < min)`
  for collision energy (`:552`) and isolation width (`:601`). The class-test
  boundaries 5.0/10.0 (in) and 4.9/10.1 (out) are reproduced for m/z, RT and
  intensity, the latter via `f32` widening (`4.9f32 < 5.0`, `10.1f32 > 10.0`).
- `reverse` is an XOR on the result for every predicate except where the source
  returns early: `IsInCollisionEnergyRange`, `IsInIsolationWindowSizeRange` and
  `IsInIsolationWindow` return `false` for `ms_level == 1` **regardless of
  `reverse`** (`:542`, `:595`, `:640`), and `IsInCollisionEnergyRange` returns
  `false` regardless of `reverse` when no precursor carries a collision energy
  (`:557`). Only level 1 is special; level 0 falls through, as in the source.
- The source `@note`s say those predicates "return true" for MS1 spectra. The
  code returns `false`; under the source's `remove_if` usage `false` is what
  keeps the spectrum, which is the note's intent. The port preserves the code
  and documents the mismatch at each item.
- `HasActivationMethod`: any precursor with any listed method (`:331`);
  `HasPrecursorCharge`: any precursor with a listed charge (`:427`);
  `InPrecursorMZRange`: every precursor inside, vacuously true without
  precursors (`:379`); `IsInIsolationWindowSizeRange` without precursors is
  not in range, so the answer is `reverse`.
- `IsInIsolationWindow` sorts its m/z list on construction (`:634`) and uses
  `lower_bound` semantics (`:651`): the first listed value not below
  `mz - lower_offset` must be at or below `mz + upper_offset` (`:655`). The
  class-test triple (200.3 ± 0.5 hits 200.0; 201.1 misses; a second precursor
  at 299.9 hits 300.0) is reproduced.
- `IsInCollisionEnergyRange` reads the metavalue `"collision energy"`, the key
  the source mzML handler fills from `MS:1000045`
  (`MzMLHandler.cpp:1900-1903`); integer-typed values convert to `f64` as the
  source `DataValue` conversion does.

## Native differences

- Range constructors validate: nonfinite bounds → `Error::InvalidValue`;
  `min > max` → `Error::InvalidRange`. The source stores any pair; an inverted
  pair then never matches (or, reversed, always matches), which is a silent
  no-op filter. `IsInIsolationWindow::new` rejects nonfinite m/z values, which
  the source would sort into an unspecified position.
- `HasActivationMethod::from_names` rejects unknown names with
  `Error::InvalidValue`; the source `StringList` keeps them and they never
  match. It also accepts the short names (`"CID"`), which the source does not.
- `InRTRange` and `InMSLevelRange` are not separate structs (see the mapping);
  closures over `rt`/`ms_level` implement `SpectrumPredicate` directly.
- `HasScanMode` compares `ScanMode` values, so the source's out-of-range
  integer mode (never matching) cannot be expressed.
- `IsInIsolationWindow` does not emit the source `OPENMS_LOG_WARN` for a
  precursor with a zero lower or upper offset (`:647`); it evaluates the same
  way. The module has no log sink; the condition is documented at the item.
- A string-typed collision-energy metavalue counts as absent here. The source
  `DataValue::operator double()` returns the union's `dou_` field for a
  `STRING_VALUE` (`DataValue.cpp:466-478`, union at `DataValue.h:409-413`) —
  the bit pattern of a `std::string*` — instead of throwing. Unconfirmed by
  execution; recorded as a C++ issue candidate for the shared log.
- `retain_spectra` / `retain_peaks_where` keep matches (like `Vec::retain`)
  whereas the source pairs the functors with `remove_if`; the module docs show
  the `reverse = true` / negated-closure translation. Both return the number of
  removed elements, which the source callers compute separately.
- The port is serial, as the source header is (no `#pragma omp` here).

## Checked boundaries and evidence

`retain_spectra_with_limits` and `retain_peaks_where_with_limits` check the
element count against `RangeFilterLimits::max_items` before any mutation, so a
rejected call leaves the container unchanged. `retain_peaks_where` also inherits
the aligned-array consistency check of `retain_peaks`, which rejects a
misaligned nonempty data array before selecting. Predicates are infallible
`bool` functions, so removal cannot stop halfway and no half-filtered container
can result. The predicates allocate nothing per call; `IsInIsolationWindow`
does one `partition_point` per precursor over its sorted list.

### Class-test sections (tier 3)

All 33 `START_SECTION`s of `RangeUtils_test.cpp` are ported to
`tests/range_utils.rs`; the destructor sections drop a constructed predicate.

| # | Section (`RangeUtils_test.cpp` line) | Test |
| --- | --- | --- |
| 1 | `InRTRange(min, max, reverse)` (32) | `in_rt_range_constructor` |
| 2 | `~InRTRange()` (37) | `in_rt_range_destructor` |
| 3 | `InRTRange::operator()` (41) | `in_rt_range_operator` — 4.9/5.0/7.5/10.0/10.1 truth table via closure and `spectra_in_rt_range` |
| 4 | `MSLevelRange(levels, reverse)` (67) | `in_ms_level_range_constructor` |
| 5 | `~InMSLevelRange()` (73) | `in_ms_level_range_destructor` |
| 6 | `InMSLevelRange::operator()` (77) | `in_ms_level_range_operator` — levels 1..5 against {2,3,4} via closure, `ms_levels`, `contains_scan_of_level` |
| 7 | `HasScanMode(mode, reverse)` (106) | `has_scan_mode_constructor` |
| 8 | `~HasScanMode()` (111) | `has_scan_mode_destructor` |
| 9 | `HasScanMode::operator()` (115) | `has_scan_mode_operator` — SIM/MASSSPECTRUM |
| 10 | `InMzRange(min, max, reverse)` (131) | `in_mz_range_constructor` |
| 11 | `~InMzRange()` (136) | `in_mz_range_destructor` |
| 12 | `InMzRange::operator()` (140) | `in_mz_range_operator` |
| 13 | `IntensityRange(min, max, reverse)` (165) | `in_intensity_range_constructor` |
| 14 | `~InIntensityRange()` (170) | `in_intensity_range_destructor` |
| 15 | `InIntensityRange::operator()` (174) | `in_intensity_range_operator` — `f32` literals |
| 16 | `IsEmptySpectrum(reverse)` (200) | `is_empty_spectrum_constructor` |
| 17 | `~IsEmptySpectrum()` (205) | `is_empty_spectrum_destructor` |
| 18 | `IsEmptySpectrum::operator()` (209) | `is_empty_spectrum_operator` — empty, then `resize(5)` |
| 19 | `IsZoomSpectrum(reverse)` (224) | `is_zoom_spectrum_constructor` |
| 20 | `~IsZoomSpectrum()` (229) | `is_zoom_spectrum_destructor` |
| 21 | `IsZoomSpectrum::operator()` (233) | `is_zoom_spectrum_operator` |
| 22 | `HasActivationMethod(methods, reverse)` (248) | `has_activation_method_constructor` — the source's `""` list is rejected here (documented difference) |
| 23 | `~HasActivationMethod()` (253) | `has_activation_method_destructor` |
| 24 | `HasActivationMethod::operator()` (257) | `has_activation_method_operator` — PSD/BIRD, BIRD, +LCID, +PD |
| 25 | `InPrecursorMZRange(mz_left, mz_right, reverse)` (313) | `in_precursor_mz_range_constructor` |
| 26 | `~InPrecursorMZRange()` (318) | `in_precursor_mz_range_destructor` |
| 27 | `InPrecursorMZRange::operator()` (322) | `in_precursor_mz_range_operator` — 150, 444, 444+150 |
| 28 | `IsInIsolationWindow(...)` (362) | `is_in_isolation_window_constructor` |
| 29 | `~IsInIsolationWindow()` (367) | `is_in_isolation_window_destructor` |
| 30 | `IsInIsolationWindow::operator()` (371) | `is_in_isolation_window_operator` — unsorted 300/100/200/400; 200.3, 201.1, +299.9 |
| 31 | `HasScanPolarity(polarity, reverse)` (412) | `has_scan_polarity_constructor` |
| 32 | `~HasScanPolarity()` (417) | `has_scan_polarity_destructor` |
| 33 | `HasScanPolarity::operator()` (421) | `has_scan_polarity_operator` |

### Native tests (tier 4)

`HasPrecursorCharge`, `IsInCollisionEnergyRange`,
`IsInIsolationWindowSizeRange` and `HasMetaValue` have no class-test section;
their tests are derived from the source bodies (any-precursor semantics, MS1
and missing-energy early returns regardless of `reverse`, closed boundaries,
metadata-key and `MS:1000045` lookups, integer conversion, string-as-absent).
Further tests cover inverted/nonfinite range rejection, the removed-count
helpers, aligned data-array movement, limit and misalignment errors leaving
the container unchanged, and MS level 0 falling through.

43 tests pass with all features and without default features; the two module
doctests pass; strict Clippy and rustfmt pass; the module compiles on Rust
1.85.0; rustdoc coverage of the module is 100% (41/41) and `src/kernel.rs`
stays at 100%.

33 ported, 0 mapped-with-evidence, 0 mapped-without-evidence, 0 unaccounted

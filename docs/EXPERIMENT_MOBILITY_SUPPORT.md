# Experiment ion mobility, rasterization and the closing audit of `MSExperiment.h` / `AreaIterator.h`

Source: Core SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4`,
`src/openms/include/OpenMS/KERNEL/MSExperiment.h`,
`src/openms/source/KERNEL/MSExperiment.cpp` and
`src/openms/include/OpenMS/KERNEL/AreaIterator.h`.

This work package closes the two headers' recorded residuals and then audits
both of them member by member. The new code is
[`src/kernel/experiment_mobility.rs`](../src/kernel/experiment_mobility.rs)
(`IMBegin`, `IMEnd`, `isIMFrame`, `get2DPeakDataIM`,
`get2DPeakDataIMPerSpectrum`, `rasterizeRTMZ`) and the mobility extension of
[`src/kernel/area_iteration.rs`](../src/kernel/area_iteration.rs)
(`Param::lowIM`/`highIM`, `getDriftTime`, and both `RangeManager` overloads of
`areaBegin`/`areaBeginConst`). Everything else in the two headers was ported by
earlier work packages and is **mapped**, not re-ported; the tables below name the
Rust counterpart for every public member so an auditor can check the claim
without reading the Rust.

Tests: [`tests/experiment_mobility.rs`](../tests/experiment_mobility.rs).
Provenance: [`tests/data/experiment_mobility_provenance.json`](../tests/data/experiment_mobility_provenance.json).
The pre-existing scalar area traversal has its own document,
[`AREA_ITERATION_SUPPORT.md`](AREA_ITERATION_SUPPORT.md), whose API table this
package extended rather than duplicated.

## Two mobility representations, one owner each

The source keeps ion mobility in two places and this port keeps that split:

* a **scalar** drift time for a whole scan, with `-1`
  (`IMTypes::DRIFTTIME_NOT_SET`) meaning unset. `IMBegin`, `IMEnd`, `isIMFrame`
  and the area iterator's mobility filter all read this one.
* a **per-peak float data array** whose *name* marks the spectrum as an ion
  mobility frame. `get2DPeakDataIM` and `get2DPeakDataIMPerSpectrum` read this
  one.

Neither rule is re-derived here. The scalar sentinel and the array-name lookup
live in [`src/kernel/spectrum_mobility.rs`](../src/kernel/spectrum_mobility.rs)
(`MSSpectrum::drift_time_if_set`, `contains_im_data`, `im_data`,
`maybe_im_data`), and the mobility *range* of a run — which uses the array when
present and the scalar otherwise — lives in
[`src/kernel/ranges.rs`](../src/kernel/ranges.rs). This module calls both.

## API mapping: `KERNEL/MSExperiment.h`

Ported here (the residual "scan-mobility/range overloads" plus the previously
unmapped rasterizer):

| Source member | Rust counterpart |
| --- | --- |
| `ConstIterator IMBegin(CoordinateType im) const` | `MSExperiment::im_begin(im) -> Result<usize>` |
| `ConstIterator IMEnd(CoordinateType im) const` | `MSExperiment::im_end(im) -> Result<usize>` |
| `bool isIMFrame() const` | `MSExperiment::is_im_frame() -> Result<bool>` |
| `void get2DPeakDataIM(min_rt, max_rt, min_mz, max_mz, ms_level, rt&, mz&, intensity&, ion_mobility&) const` | `MSExperiment::get_2d_peak_data_im(bounds, ms_level) -> Result<FlatPeakDataIm>`; `append_2d_peak_data_im(..., &mut FlatPeakDataIm)` for the source's append semantics; `*_with_limits` for explicit ceilings |
| `void get2DPeakDataIMPerSpectrum(..., rt&, mz&, intensity&, ion_mobility&) const` | `MSExperiment::get_2d_peak_data_im_per_spectrum(bounds, ms_level) -> Result<SpectrumPeakDataIm>`; `append_2d_peak_data_im_per_spectrum(...)`; `*_with_limits` |
| `enum class RasterAggregation { SUM, MAX }` | `kernel::spectrum_mobility::RasterAggregation::{Sum, Max}` — one enum serves both rasterizers, as the source declares one per class with identical members |
| `void rasterizeRTMZ(float* output, rt_bins, mz_bins, min_rt, max_rt, min_mz, max_mz, ms_level, aggregation) const` | `MSExperiment::rasterize_rt_mz(&RtMzRaster) -> Result<Vec<f32>>`, with `RtMzRaster::{new, with_aggregation, pixels}` holding the scalar parameters and the owned `Vec` replacing the caller-supplied `float*` |
| `AreaIterator areaBegin(const RangeManagerType& range, UInt ms_level = 1)` | `MSExperiment::area_begin_mut_from_ranges(&RangeManager, ms_level)` |
| `ConstAreaIterator areaBeginConst(const RangeManagerType& range, UInt ms_level = 1) const` | `MSExperiment::area_begin_from_ranges(&RangeManager, ms_level)`; `AreaOptions::from_range_manager` exposes the dimension extraction on its own |

Mapped to code from earlier work packages:

| Source member | Rust counterpart |
| --- | --- |
| `MSExperiment()` | `MSExperiment::new()` / `Default` (`src/kernel.rs`) |
| `MSExperiment(const MSExperiment&)`, `operator=(const MSExperiment&)` | `Clone` |
| `MSExperiment(MSExperiment&&)`, `operator=(MSExperiment&&) &` | native move; Rust moves by default and needs no member |
| `~MSExperiment()` | `Drop` glue; no member |
| `operator=(const ExperimentalSettings& source)` | assign `MSExperiment::settings` (`src/metadata/experimental_settings.rs`) |
| `operator==`, `operator!=` | `PartialEq` on `MSExperiment` |
| `Size size() const noexcept`, `bool empty() const noexcept` | `MSExperiment::len`, `is_empty` |
| `void resize(Size n)`, `void reserve(Size n)` | `experiment.spectra.resize_with(n, Default::default)` / `Vec::reserve` on the public `spectra` field |
| `SpectrumType& operator[](Size)` (both) | `Index`/`IndexMut` on the public `spectra` field |
| `begin`, `cbegin`, `end`, `cend` (all overloads) | `experiment.spectra.iter()` / `iter_mut()` |
| `void reserveSpaceSpectra(Size)`, `void reserveSpaceChromatograms(Size)` | `Vec::reserve` on `spectra` / `chromatograms` |
| `template<class Container> void get2DData(Container&) const` | `MSExperiment::get_2d_data() -> Result<Vec<Peak2D>>`, `append_2d_data(&mut Vec<Peak2D>)` (`src/kernel/experiment_2d.rs`). The C++ template is duck-typed on `push_back`/`setRT`/`setMZ`/`setIntensity`; the Rust form fixes the element type to `Peak2D`, which every source instantiation in the SDK uses |
| `template<class Container> void set2DData(const Container&)` | `MSExperiment::set_2d_data(&[Peak2D]) -> Result<MSExperiment>` |
| `template<class Container> void set2DData(const Container&, const StringList&)` | `MSExperiment::set_2d_data_rich(&[RichPeak2D], &[String])` |
| `template<bool add_mass_traces, class Container> void set2DData(const Container&)` | **not ported**: the `add_mass_traces = true` specialization, which expands `num_of_masstraces` / `masstrace_intensity_<i>` meta values into peaks at ¹³C spacing, has no Rust counterpart. See *Deferrals* below |
| `struct SumIntensityReduction` | the default reduction of `MSExperiment::aggregate` / `extract_xics` (`src/kernel/experiment_aggregation.rs`); the source needs a named functor only for template deduction |
| `template<class MzReductionFunctionType> aggregate(ranges, ms_level, func)` | `MSExperiment::aggregate_with(regions, ms_level, closure)`, `aggregate_with_limits` |
| `aggregate(ranges, ms_level)` | `MSExperiment::aggregate(regions, ms_level)` |
| `template<class MzReductionFunctionType> extractXICs(ranges, ms_level, func)` | `MSExperiment::extract_xics_with(...)`, `extract_xics_with_limits` |
| `extractXICs(ranges, ms_level)` | `MSExperiment::extract_xics(regions, ms_level)` |
| `aggregateFromMatrix(ranges, ms_level, mz_agg)` | `MSExperiment::aggregate_from_matrix(rows, ms_level, MzAggregation)` |
| `extractXICsFromMatrix(ranges, ms_level, mz_agg)` | `MSExperiment::extract_xics_from_matrix(rows, ms_level, MzAggregation)` |
| `AreaIterator areaBegin(min_rt, max_rt, min_mz, max_mz, UInt ms_level)` | `MSExperiment::area_begin_mut(min_rt, max_rt, min_mz, max_mz, ms_level)` |
| `ConstAreaIterator areaBeginConst(min_rt, max_rt, min_mz, max_mz, UInt ms_level) const` | `MSExperiment::area_begin(...)` |
| `AreaIterator areaEnd()`, `ConstAreaIterator areaEndConst() const` | `AreaIter::default()` / an exhausted iterator; Rust iteration has no end sentinel to hand out |
| `void get2DPeakDataPerSpectrum(...)` | `MSExperiment::get_2d_peak_data_per_spectrum(bounds, ms_level)` (`src/kernel/peak_data.rs`) |
| `void get2DPeakData(...)` | `MSExperiment::get_2d_peak_data(bounds, ms_level)` |
| `ConstIterator RTBegin(CoordinateType) const`, `Iterator RTBegin(CoordinateType)` | `MSExperiment::rt_begin(rt) -> Result<usize>` (one index serves both) |
| `ConstIterator RTEnd(CoordinateType) const`, `Iterator RTEnd(CoordinateType)` | `MSExperiment::rt_end(rt) -> Result<usize>` |
| `void clearRanges()` | no counterpart by design: ranges are computed on demand, so there is no cache to clear. `src/kernel/ranges.rs` records the decision |
| `getMinRT`, `getMaxRT`, `getMinMZ`, `getMaxMZ`, `getMinIntensity`, `getMaxIntensity`, `getMinMobility`, `getMaxMobility` | `MSExperiment::combined_range_manager()` then `min_rt()`/`max_rt()`/`min_mz()`/… on the returned `RangeManager`; the legacy triple `combined_ranges()` gives RT/m/z/intensity only |
| `void updateRanges()` | no counterpart: `spectrum_range_manager()`, `chromatogram_range_manager()` and `combined_range_manager()` each compute their manager on demand |
| `UInt64 getSize() const` | `MSExperiment::total_peak_count() -> Result<u64>` (`src/kernel/experiment_summary.rs`) |
| `std::vector<UInt> getMSLevels() const` | `MSExperiment::ms_levels() -> Vec<u32>` |
| `UInt64 getSqlRunID() const`, `void setSqlRunID(UInt64)` | the public `MSExperiment::sql_run_id` field. The source stores it as the `"sqMassRunID"` meta value and answers `0` when absent; a typed field keeps that default without a string key |
| `void sortSpectra(bool sort_mz = true)` | `MSExperiment::sort_spectra(sort_mz) -> Result<()>` |
| `void sortChromatograms(bool sort_rt = true)` | `MSExperiment::sort_chromatograms(sort_rt) -> Result<()>` |
| `bool isSorted(bool check_mz = true) const` | `MSExperiment::is_sorted(check_mz) -> bool` |
| `void reset()` | `MSExperiment::clear(true)`; the source's `reset()` and `clear(true)` differ only in that `reset()` also clears the range cache this port does not keep |
| `bool clearMetaDataArrays()` | `MSExperiment::clear_meta_data_arrays() -> Result<bool>` |
| `const ExperimentalSettings& getExperimentalSettings() const`, non-const | the public `MSExperiment::settings` field |
| `void getPrimaryMSRunPath(StringList& toFill) const` | `ExperimentalSettings::primary_ms_run_path` handling in `src/metadata/experimental_settings.rs`, reached through `experiment.settings` |
| `ConstIterator getPrecursorSpectrum(ConstIterator) const`, `int getPrecursorSpectrum(int) const` | `MSExperiment::precursor_spectrum_index(index) -> Result<Option<usize>>`; `None` replaces both the past-the-end iterator and the `-1` |
| `ConstIterator getFirstProductSpectrum(ConstIterator) const`, `int getFirstProductSpectrum(int) const` | **not ported**. See *Deferrals* below |
| `void swap(MSExperiment& from)` | `std::mem::swap(&mut a, &mut b)` |
| `void setSpectra(const std::vector<MSSpectrum>&)`, `setSpectra(std::vector<MSSpectrum>&&)` | assign the public `spectra` field (by clone or by move) |
| `void addSpectrum(const MSSpectrum&)`, `addSpectrum(MSSpectrum&&)` | `experiment.spectra.push(...)` |
| `const std::vector<MSSpectrum>& getSpectra() const`, non-const | the public `spectra` field |
| `ConstIterator getClosestSpectrumInRT(double) const` and its three siblings | `MSExperiment::closest_spectrum_in_rt(rt, ms_level)`, with `ms_level == 0` meaning "any level" as in the source's single-argument overload |
| `void setChromatograms(...)` (both), `addChromatogram(...)` (both), `getChromatograms()` (both) | the public `chromatograms` field and `Vec::push` |
| `MSChromatogram& getChromatogram(Size)` / const, `MSSpectrum& getSpectrum(Size)` / const | `experiment.chromatograms[id]` / `experiment.spectra[id]`, or `get(id)` for a checked lookup |
| `Size getNrSpectra() const`, `Size getNrChromatograms() const` | `experiment.spectra.len()`, `experiment.chromatograms.len()` |
| `const MSChromatogram calculateTIC(float rt_bin_size = 0, UInt ms_level = 1) const` | `MSExperiment::calculate_tic(ms_level)` for the unbinned form and `calculate_tic_binned(rt_bin_size, ms_level)` for the resampling form |
| `void clear(bool clear_meta_data)` | `MSExperiment::clear(clear_metadata)` |
| `bool containsScanOfLevel(size_t) const` | `MSExperiment::contains_scan_of_level(ms_level) -> Result<bool>` |
| `bool hasZeroIntensities(size_t) const` | `MSExperiment::has_zero_intensities(ms_level) -> Result<bool>` |
| `const SpectrumRangeManagerType& spectrumRanges() const` | `MSExperiment::spectrum_range_manager() -> Result<SpectrumRangeManager>` |
| `const ChromatogramRangeManagerType& chromatogramRanges() const` | `MSExperiment::chromatogram_range_manager() -> Result<RangeManager>` |
| `const RangeManagerType& combinedRanges() const` | `MSExperiment::combined_range_manager() -> Result<RangeManager>` |
| `std::ostream& operator<<(std::ostream&, const MSExperiment&)` (free function) | **not ported**: there is no `Display` for `MSExperiment`. See *Deferrals* below |
| `using MSRun = MSExperiment` (namespace alias) | not ported; one name, `MSExperiment` |
| `PeakMap` (the `StandardTypes.h` alias this header pulls in) | not ported; one name, `MSExperiment` |

Type aliases, nested helpers and protected state:

| Source member | Rust counterpart |
| --- | --- |
| `PeakT`, `PeakType` → `Peak1D`; `ChromatogramPeakT`, `ChromatogramPeakType` → `ChromatogramPeak` | `kernel::Peak1D`, `kernel::ChromatogramPeak`; the aliases carry no behaviour |
| `CoordinateType` (`double`), `IntensityType` (`float`) | `f64` and `f32` throughout |
| `RangeManagerType` = `RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>` | `RangeManager::experiment()` (`src/kernel/ranges.rs`), whose dimension set is a run-time value because Rust has no variadic generics |
| `SpectrumRangeManagerType`, `ChromatogramRangeManagerType` | `ranges::SpectrumRangeManager`, `RangeManager::chromatogram_manager()` |
| `SpectrumType`, `ChromatogramType`, `Base`, `value_type` | `MSSpectrum`, `MSChromatogram`, `Vec<MSSpectrum>` |
| `Iterator`, `ConstIterator`, `iterator`, `const_iterator` | `slice::IterMut<MSSpectrum>` / `slice::Iter<MSSpectrum>`, obtained from the public field |
| `AreaIterator`, `ConstAreaIterator` | `kernel::AreaIterMut`, `kernel::AreaIter` |
| protected `spectra_`, `chromatograms_` | the public `spectra`, `chromatograms` fields |
| protected `spectrum_ranges_`, `chromatogram_ranges_`, `combined_ranges_` | no counterpart: these three caches are what on-demand computation replaces |
| private `ContainerAdd_<T, false>::addData_` (both overloads) | the private `PointSource` trait of `src/kernel/experiment_2d.rs` |
| private `ContainerAdd_<T, true>::addData_` | not ported; the mass-trace deferral above |
| private `createSpec_(rt)`, `createSpec_(rt, metadata_names)` | private row construction inside `src/kernel/experiment_2d.rs` |

## API mapping: `KERNEL/AreaIterator.h`

`Internal::AreaIterator` is a five-parameter template used at exactly two
instantiations, which become `AreaIter` (shared) and `AreaIterMut`
(exclusive) in [`src/kernel/area_iteration.rs`](../src/kernel/area_iteration.rs).
The iterator pair collapses into one Rust iterator, so the `Param` builder
becomes the owned `AreaOptions` value and the end iterator becomes exhaustion.

| Source member | Rust counterpart |
| --- | --- |
| `Param::Param(first, begin, end, uint8_t ms_level)` | `AreaOptions::new(bounds, ms_level)`; the three iterators are replaced by the experiment plus the RT bounds, since the Rust iterator borrows the run |
| `static Param Param::end()` | `AreaIter::default()` |
| `Param& Param::operator=(const Param&)` | `Copy` on `AreaOptions` |
| `Param& Param::lowMZ(CoordinateType)`, `highMZ` | `AreaBounds::new(min_rt, max_rt, min_mz, max_mz)` / the `AreaBounds::mz` field |
| `Param& Param::lowIM(CoordinateType)`, `highIM` | **ported here**: `AreaOptions::with_mobility(min_im, max_im)` and the `AreaOptions::mobility` field |
| `Param& Param::msLevel(int8_t)` | the `AreaOptions::ms_level` field; `AreaOptions::source_compatible` reproduces the byte narrowing |
| protected `first_`, `current_scan_`, `end_scan_`, `current_peak_`, `end_peak_`, `low_mz_`, `high_mz_`, `low_im_`, `high_im_`, `ms_level_`, `is_end_` | the private interval plan and cursor of `AreaIter`/`AreaIterMut` |
| `iterator_category` (`forward_iterator_tag`) | `Iterator` + `FusedIterator`; the plan also makes `ExactSizeIterator` sound |
| `value_type`, `reference`, `pointer` | `AreaPeak` / `AreaPeakMut`, which bundle the peak with its indices |
| `difference_type` (`unsigned int`) | `usize` in `size_hint`/`len` |
| `explicit AreaIterator(const Param&)` | `MSExperiment::area_iter(options)` / `area_iter_mut(options)` |
| `AreaIterator()` | `AreaIter::default()` |
| `~AreaIterator()` | `Drop` glue; no member |
| `AreaIterator(const AreaIterator&)` | `Clone` on `AreaIter`; deliberately absent on `AreaIterMut`, which would alias |
| `AreaIterator& operator=(const AreaIterator&)` | assignment/move; the source's hand-written operator exists only to avoid copying the iterators of an end iterator |
| `bool operator==`, `bool operator!=` | `PartialEq`/`Eq` on `AreaIter`, comparing the current peak's address, with all exhausted iterators equal |
| `AreaIterator& operator++()`, `AreaIterator operator++(int)` | `Iterator::next`, which yields the current item and then advances |
| `reference operator*() const`, `pointer operator->() const` | the yielded `AreaPeak::peak` / `AreaPeakMut::peak`; `AreaIter::peek()` is the non-consuming form |
| `CoordinateType getRT() const` | `AreaPeak::spectrum.rt`, `AreaPeakMut::rt` |
| `CoordinateType getDriftTime() const` | **ported here**: `AreaPeak::drift_time()` and the `AreaPeakMut::drift_time` field |
| `const SpectrumT& getSpectrum() const` | `AreaPeak::spectrum`. `AreaPeakMut` deliberately exposes RT, MS level and drift time instead of a whole-spectrum reference, which would alias its exclusive peak borrow |
| `PeakIndex getPeakIndex() const` | `AreaPeak::spectrum_index` / `peak_index` (and the same two fields on `AreaPeakMut`), which are exactly `PeakIndex`'s components; `kernel::PeakIndex::new(spectrum, peak)` builds the value type |
| private `nextScan_()` | the interval planner of `area_iteration.rs`, which resolves every scan and peak window before the first item is yielded |
| private `p_` | the iterator's plan and cursor |

## Preserved source conventions

* **Scalar drift time only, in the searches and the area filter.**
  `MSSpectrum::IMLess` compares `getDriftTime()`, and `AreaIterator::nextScan_`
  tests `containsMobility(current_scan_->getDriftTime())`. A spectrum that
  carries a per-peak ion mobility array but no scalar drift time therefore
  presents the sentinel `-1`, and any mobility window above `-1` excludes it.
  The port does the same and says so at the item.
* **An unrestricted mobility window selects every scan.** The source defaults
  `low_im_`/`high_im_` to `lowest()`/`max()` and the scalar `areaBegin`
  overloads pass `RangeMobility{}.getNonEmptyRange()`, a full range. `mobility:
  None` is that case.
* **`getNonEmptyRange()` semantics for the `RangeManager` overload.** An empty
  dimension does not restrict the area. A dimension the Rust manager does not
  carry behaves identically, because the source's manager is a fixed template
  and cannot express absence. The intensity dimension is ignored, as the source
  iterator has no intensity filter.
* **`isIMFrame` compares against the previous scan only.** Drift times
  `1, 2, 1` report a frame; `1, 1, 2` do not. The reference retention time is
  the *first* spectrum's, compared with `!=`.
* **`-1` as the missing-mobility sentinel** in both bulk exports, which is the
  same value a genuine mobility of `-1` would produce. The source cannot
  distinguish them and neither does the port.
* **Row grouping by narrowed retention time** in `get2DPeakDataIMPerSpectrum`:
  the `f64` retention time is compared against the `f32` cursor, so scans whose
  retention times narrow to the same `f32` merge into one row. The ion mobility
  array is fetched only when a row *starts*, so a merged second spectrum's peaks
  are looked up in the first spectrum's array — preserved, and pinned by
  `merged_rows_use_the_first_spectrum_ion_mobility_array`.
* **The array is indexed by the peak's position inside its own spectrum**
  (`getPeakIndex().peak`), not by its position in the selection, so an m/z
  filter does not shift the lookup.
* **Both exports append.** The source never clears the vectors it is handed;
  `append_*` is the faithful form and `get_*` is the convenience wrapper.
* **Exact MS-level matching with source byte narrowing.** `Size`/`UInt` levels
  reach an `uint8_t` parameter stored in an `int8_t` field, so `256` selects
  level `0` and `255` selects `u32::MAX`. Reproduced by
  `AreaOptions::source_compatible` and by every wrapper that takes a source
  level, including the two `RangeManager` ones.
* **Raster geometry.** Row-major with m/z as the slow axis; a peak exactly at
  `max_rt` or `max_mz` is clamped into the last bin instead of overflowing; the
  RT and m/z windows are half-inclusive at the top through `RTEnd`/`MZEnd`, so
  peaks outside them are excluded rather than clamped in.

## Native differences

* **No cached ranges.** The source's `updateRanges()`/`clearRanges()` cycle has
  no counterpart; the backward-compatible `getMinMobility()` family is answered
  by `combined_range_manager()`. This was settled crate-wide before this package
  and is only restated here.
* **Indices instead of iterators.** `IMBegin`/`IMEnd`/`RTBegin`/`RTEnd` return
  `usize`, with `spectra.len()` for the past-the-end iterator. The source
  declares only const `IMBegin`/`IMEnd` even though its own class test names a
  mutable overload; one index serves both.
* **Checked preconditions.** The source states "make sure the spectra are
  sorted" as a `@note` and does not check it; `im_begin`/`im_end` return
  `Error::UnsortedData`, and `rasterize_rt_mz` checks both the RT order of the
  run and the m/z order of every contributing spectrum. `OPENMS_PRECONDITION`
  compiles to nothing in a release build (`Macros.h:91`), so the source's
  `areaBegin` swapped-bound and sortedness guards do not exist in shipped
  builds; the port's are always on.
* **Non-finite values are refused, not silently dropped.** `RangeBase::contains`
  answers `false` for NaN, so a source scan with a NaN drift time is skipped
  without trace; with a mobility window set, the port returns
  `Error::InvalidValue` instead. Without a window the drift time is never read,
  which is why a NaN drift time elsewhere in the run does not affect an
  unfiltered traversal. `rasterize_rt_mz` rejects non-finite bounds,
  coordinates and intensities, where the source would reach
  `static_cast<Int64>(NaN)` — undefined behaviour in C++.
* **Bounds-checked ion mobility arrays.** The source indexes the array with the
  peak index and never compares lengths, so an array shorter than the peak list
  reads out of bounds. The port returns `Error::InvalidValue` and writes
  nothing.
* **The flat export's retention-time `-1` defect is not reproduced.**
  `get2DPeakDataIM` declares its row cursor `float t = -1.0;` *inside* the
  per-peak loop (`MSExperiment.cpp:240-242`), so a spectrum whose retention time
  is exactly `-1` never has its mobility array fetched and reports `-1` for
  every peak. `-1` is OpenMS's unset retention time, so the case is reachable,
  and `-1` is also the "no mobility" sentinel, so the loss is undetectable
  downstream. The port always fetches. No source-compatibility option is
  offered, because the source path only destroys information — there is nothing
  a caller could want it for. The row-grouping export keeps the source's cursor
  semantics, because there the cursor is a real grouping key.
* **An unopenable first row is an error, not a write past the end.** When the
  first selected retention time is exactly `-1` and the output holds no row,
  the source appends through `mz.back()` into an empty vector. The port returns
  `Error::InvalidValue`; when the output *does* hold a row, it continues that
  row exactly as the source does, mobility sentinel included.
* **Owned raster buffer.** `rasterize_rt_mz` returns a `Vec<f32>` instead of
  filling a caller-supplied `float*`, which removes the source's
  `Exception::NullPointer` path and its `numpy.empty` advice; that `@note` and
  the Python `@code` block describe the pyOpenMS buffer and have no counterpart.
  The source multiplies `rt_bins * mz_bins` unchecked; `RtMzRaster::pixels`
  returns an error on overflow.
* **Serial rasterizer.** The source runs `rasterizeRTMZ` under
  `#pragma omp parallel for` with per-thread accumulation buffers and a
  thread-count heuristic (`MSExperiment.cpp:333-515`). The port reproduces the
  source's own single-threaded branch, which computes the same image; it
  introduces no threads, so the performance gap is real and stated. The same
  applies to `aggregate`/`extractXICs`, which the aggregation package already
  records.
* **One `RasterAggregation`.** The source declares a separate nested enum in
  `MSSpectrum` and in `MSExperiment` with the same two members; the port has
  one, in `spectrum_mobility`.
* **A typed `sql_run_id` field** instead of the `"sqMassRunID"` meta value, with
  the same `0`-when-unset default.

## Checked boundaries and evidence

`ExperimentMobilityLimits::default()` allows one million spectra, ten million
visited peaks, ten million written points, one million written rows, 50 million
work units and 256 MiB of newly allocated payload. Every ceiling counts what is
already in the caller's output, is checked before anything is allocated, and is
configurable per call. The area selection shares one budget with the export that
consumes it, so an empty selection cannot bypass the whole-input validation the
area planner performs. Both exports build a complete replacement and assign it
at the end, so a late failure leaves the caller's vectors byte-identical — pinned
by `short_or_nonfinite_ion_mobility_arrays_are_refused_atomically` and
`mobility_exports_append_and_continue_the_callers_last_row`.

`MSExperiment::MAX_MOBILITY_ITEMS` (100,000,000) caps the spectra a scan-mobility
search or `is_im_frame` will walk and the peaks `rasterize_rt_mz` will visit;
`MSExperiment::MAX_RASTER_PIXELS` (16,777,216) caps one image at 64 MiB. Both
match the `MSSpectrum` constants of the same name, so a spectrum whose IM frame
can be rasterized can also be rasterized as part of a run.

Errors use only existing variants: `Error::InvalidValue` for a non-finite or
out-of-range value, an exhausted ceiling, a short ion mobility array and a
misaligned output; `Error::UnsortedData` for a violated ordering precondition;
`Error::InvalidRange` for the two reversed-range throws of `rasterizeRTMZ`.

**Evidence tier 3 (source review)** for every literal: they are transcribed from
`MSExperiment_test.cpp` and `AreaIterator_test.cpp` at the pinned revision. No
C++ was built or executed and no retained C++ output exists for these headers, so
tier 1 is not available. The checks the source does not perform — the ordering
and finiteness guards, the array-length guard, the ceilings and the atomicity —
are tier 4 (Rust-only invariants) derived from the source lines anchored in the
provenance manifest.

### Class-test section accounting

`AreaIterator_test.cpp` has 15 sections; all 15 are accounted for. Six are
ported in `tests/experiment_mobility.rs` or `tests/area_iteration.rs` because
they carry more than five assertion macros:

| Section | Where |
| --- | --- |
| `getDriftTime()` (7 macros) | ported: `source_area_iterator_drift_time_and_scan_mobility_window` |
| `getRT()` (9) | ported: `source_area_iterator_rt_and_ms_level_cases` |
| `[EXTRA] Overall test` (53) | ported across `tests/area_iteration.rs::source_area_fixture_quadrants_indices_rt_and_empty_scan` (every quadrant and empty window) and, for the ion-mobility and MS-level cases, `source_area_iterator_drift_time_and_scan_mobility_window` and `source_area_iterator_rt_and_ms_level_cases` |
| `operator==` (6), `operator!=` (6) | ported: `tests/area_iteration.rs::const_clones_have_independent_cursors_address_equality_and_fused_end` |
| `getPeakIndex()` (17) | ported: `tests/area_iteration.rs::source_area_fixture_quadrants_indices_rt_and_empty_scan`, which asserts the same eight `(spectrum, peak)` pairs `(0,0) (0,1) (1,0) (1,1) (3,0) (3,1) (4,0) (4,1)` |

The remaining nine (`AreaIterator()`, the `Param` constructor, `~AreaIterator()`,
the copy constructor, `operator=`, `operator*`, `operator->`, and both `++`
forms; 0–4 macros each) are mapped to
`tests/area_iteration.rs::const_clones_have_independent_cursors_address_equality_and_fused_end`
and `source_area_fixture_quadrants_indices_rt_and_empty_scan`; the latter
reproduces the `operator*`/`operator->`/`++` value `510.0` as the first peak of
the `(0, 7, 505, 520)` window and then `506.0`.

`MSExperiment_test.cpp` has 66 sections. 62 are accounted for: 8 are ported in
`tests/experiment_mobility.rs` (`areaBeginConst(range)`, `areaBegin(range)`,
`IMBegin`, `IMEnd`, `isIMFrame`, `get2DPeakDataIM`,
`get2DPeakDataIMPerSpectrum`, `Backward compatibility tests`), and 54 are
ported or mapped by the earlier packages named in the tables above —
`tests/kernel.rs`, `tests/area_iteration.rs`, `tests/peak_data.rs`,
`tests/experiment_2d.rs`, `tests/experiment_summary.rs`,
`tests/experiment_aggregation.rs`, `tests/experimental_settings.rs`,
`tests/precursor_workflow.rs` and `tests/ranges.rs`. Four sections are **not**
accounted for, covering the three deferrals below:

| Section | Macros | Why |
| --- | --- | --- |
| `template<bool add_mass_traces, class Container> void set2DData(...)` | 11 | the mass-trace specialization is not ported |
| `ConstIterator getFirstProductSpectrum(ConstIterator) const` | 11 | not ported |
| `int getFirstProductSpectrum(int) const` | 10 | not ported |
| `std::ostream& operator<<(std::ostream&, const MSExperiment&)` | 5 | no `Display` for `MSExperiment` |

Four sections, three deferrals: the two `getFirstProductSpectrum` overloads
share one. Across both class tests: 81 sections, 14 ported as new tests here,
63 mapped to existing Rust tests with a cited value, 4 unaccounted.

## Deferrals

These three are outside this work package's file ownership. Each needs a file
this package does not own, so the gap is recorded rather than half-closed:

1. **`set2DData<add_mass_traces = true>`** — expanding a feature's
   `num_of_masstraces` and `masstrace_intensity_<i>` meta values into peaks
   spaced by `Constants::C13C12_MASSDIFF_U / charge` (charge `0` treated as `1`),
   and throwing `Exception::Precondition` when an expected
   `masstrace_intensity_<i>` is missing. Belongs in
   `src/kernel/experiment_2d.rs`, next to the other 2D import. The constant
   already exists in `src/concept/constants.rs` and the meta-value names in
   `src/analysis/feature_finding_metabo.rs`.
2. **`getFirstProductSpectrum`** (both overloads) — the forward search for the
   first spectrum of the next higher MS level whose first precursor names the
   current scan, bounded by the next scan of a lower level. Belongs in
   `src/kernel.rs` beside `precursor_spectrum_index`.
3. **`operator<<(std::ostream&, const MSExperiment&)`** — the
   `MSEXPERIMENT BEGIN` / `MSSPECTRUM BEGIN` / `MSCHROMATOGRAM BEGIN` diagnostic
   stream. `MSSpectrum` and `MSChromatogram` already have their `Display`
   counterparts; the experiment-level one belongs with them in `src/kernel.rs`.

## Source defect candidates

Recorded here because `OpenMS_CPP_ISSUES.md` is owned by the integrating agent.

1. **`get2DPeakDataIM` never reads the ion mobility array of a spectrum at
   retention time `-1`** (`MSExperiment.cpp:240-242`). `DriftTimeUnit unit`,
   `std::vector<float> im` and `float t = -1.0` are all declared inside the
   per-peak loop, so the `it.getRT() != t` guard is dead for every retention
   time except exactly `-1`, where it turns the guard off and makes the function
   report `-1` mobility for peaks that have a real value. `-1` is OpenMS's own
   unset retention time. Proposed fix: hoist the three declarations out of the
   loop, as the sibling `get2DPeakDataIMPerSpectrum` already does. The port
   always reads the array and documents the difference.
2. **`isIMFrame()` reports `true` for a single spectrum with no ion mobility at
   all** (`MSExperiment.cpp:1322-1333`). `last_drift` starts at
   `numeric_limits<double>::lowest()`, so the unset sentinel `-1` differs from
   it and the loop ends with "RT stable, IM changing". Proposed fix: require at
   least two spectra, or reject the `DRIFTTIME_NOT_SET` sentinel. The port
   reproduces the source answer and documents it.
3. **Both mobility exports index the ion mobility array without a bounds
   check** (`MSExperiment.cpp:196-197`, `MSExperiment.cpp:252-253`). The index is
   the peak's position inside its spectrum, and nothing requires the array to be
   as long as the peak list — `containsIMData` inspects only the array's *name*.
   An array shorter than the peak list reads out of bounds. Proposed fix: compare
   the lengths and skip or throw. The port returns `Error::InvalidValue`.
4. **A merged retention-time row reads the wrong spectrum's array**
   (`MSExperiment.cpp:186-190`). Two spectra whose retention times narrow to the
   same `f32` share one row, and the array is fetched only when the row opens,
   so the second spectrum's peaks are looked up in the first's. Proposed fix:
   fetch per spectrum rather than per row. The port preserves the behaviour and
   pins it with a test.

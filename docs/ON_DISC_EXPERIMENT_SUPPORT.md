# On-disc MS experiment

Rust: `src/kernel/on_disc_experiment.rs` (feature `mzml`), registered in
`src/kernel.rs`
Tests: `tests/on_disc_experiment.rs`
Manifest: `tests/data/on_disc_experiment_provenance.json`

Source, pinned at `bc9cc12514c768385ce121d6ca4bb710fe1983c4`:

- `src/openms/include/OpenMS/KERNEL/OnDiscMSExperiment.h` (274 lines) — ported here
- `src/openms/source/KERNEL/OnDiscMSExperiment.cpp` (233 lines) — ported here
- `src/openms/include/OpenMS/FORMAT/IndexedMzMLFileLoader.h` and its `.cpp` — read
  for the consumer contract; only the materialisation loop of `store` is ported,
  as `load_experiment`

The engine is `src/format/indexed_mzml_handler.rs`
(`docs/INDEXED_MZML_HANDLER_SUPPORT.md`), which already owns index parsing,
counts, fetch by index and by native id, the raw record XML, the checked byte
ranges and the `PeakFileOptions` execution. The source class body is almost pure
delegation to `Internal::IndexedMzMLHandler`, so this module is a thin facade: no
byte-range arithmetic and no XML handling of its own.

## What the facade adds

Two things, both the source's own:

1. **A metadata experiment.** `loadMetaData_` loads the whole file through
   `FileHandler` with `setFillData(false)` and **no** other filter, deliberately,
   so that metadata index *i* is file record *i*
   (`OnDiscMSExperiment.cpp:57-63`). Every field except peak data is therefore
   available without touching the record: retention time, MS level, precursors,
   instrument settings, data arrays descriptors, identifiers.
2. **`PeakFileOptions` filtering around each fetch**, in two groups:
   - RT range, MS level and precursor m/z are tested against the metadata
     **before** any peak I/O, and an excluded record is returned as its metadata
     spectrum with no peaks. That preserves every index — the source's `@note`
     contrasts it with in-memory loading, where a filtered spectrum disappears
     from the container.
   - m/z range and intensity range select peaks **inside** a record, so they
     require the read.

## API mapping — `OnDiscMSExperiment.h`

Every public, private and protected member of the header appears here.

| C++ member | Rust | Notes |
|---|---|---|
| `OnDiscMSExperiment() = default` | `OnDiscMSExperiment::new()`, `Default` | Holds no file; the source's own constructor documentation says to call `openFile`. |
| implicit `~OnDiscMSExperiment()` | `Drop` on the owned `File` inside the handler | Not declared in the header, but the class test has a section for it. |
| `bool openFile(const std::string&, bool skipMetaData = false)` | `open_file(path, skip_metadata) -> Result<bool>` | Same `bool`: whether the index parsed. An unreadable file or a broken metadata load is `Err`, where the source lets the `FileHandler` exception escape. Built into temporaries and committed at the end. |
| — | `OnDiscMSExperiment::open(path) -> Result<Self>` | Native. The commented-out `OnDiscMSExperiment(const std::string&)` section of the class test, with the source's `false` return turned into `Error::Parse`. |
| — | `OnDiscMSExperiment::with_limits(OnDiscLimits)` | Native. Explicit ceilings. |
| `OnDiscMSExperiment(const OnDiscMSExperiment&)` | `try_clone() -> Result<Self>` | The source copy reopens the file ("critical for parallel access"); reopening can fail, so this is not `Clone`. It also rebuilds the two native-id caches, which the source copy loses (defect 1 below). |
| `bool operator==(const OnDiscMSExperiment&) const` | `PartialEq::eq` | File path and metadata only, with the source's null-pointer fallback preserved. |
| `bool operator!=(const OnDiscMSExperiment&) const` | `PartialEq::ne` | Derived from `eq`, as the source derives it from `operator==`. |
| `bool isSortedByRT() const` | `is_sorted_by_rt()` | `MSExperiment::is_sorted(false)` on the metadata; `false` without metadata, as in the source. |
| `Size size() const` | `len()` | |
| `bool empty() const` | `is_empty()` | Spectra only; a chromatogram-only file is empty, as in the source. |
| `Size getNrSpectra() const` | `spectrum_count()` | From the index, so zero when the index did not parse even if the metadata lists spectra. |
| `Size getNrChromatograms() const` | `chromatogram_count()` | |
| `std::shared_ptr<const ExperimentalSettings> getExperimentalSettings() const` | `experimental_settings() -> Option<&ExperimentalSettings>` | `None` is the source's null pointer after `openFile(..., true)`; its own class test asserts that null. |
| `std::shared_ptr<PeakMap> getMetaData() const` | `metadata() -> Option<&MSExperiment>` | Shared reference, not a mutable handle: see "Native differences". |
| `MSSpectrum operator[](Size n)` | `spectrum(n)` | `Index` cannot express it — the value is produced on demand and the fetch needs `&mut self`. The source's `operator[]` is a one-line forward to `getSpectrum`. |
| `MSSpectrum getSpectrum(Size id)` | `spectrum(index) -> Result<MSSpectrum>` | Full filtering, in the source's order. |
| `OpenMS::Interfaces::SpectrumPtr getSpectrumById(Size id)` | not ported | The crate has no OpenSWATH `Interfaces` layer (`INTERFACES/DataStructures.h` is `unmapped` in `docs/core-sdk-coverage.json`). Its `MZArray` and `IntensityArray` are the one `peaks` vector of `spectrum(index)`. Its second property — bypassing both the metadata and `options_` — has no separate Rust entry point. |
| `MSChromatogram getChromatogram(Size id)` | `chromatogram(index) -> Result<MSChromatogram>` | RT and intensity ranges select points; no before-I/O check, as in the source. |
| `OpenMS::Interfaces::ChromatogramPtr getChromatogramById(Size id)` | not ported | As the spectrum pointer overload; `TimeArray`/`IntensityArray` are `chromatogram(index)`'s `peaks`. |
| `MSChromatogram getChromatogramByNativeId(const std::string&)` | `chromatogram_by_native_id(&str) -> Result<MSChromatogram>` | Unfiltered, as in the source. Unknown id is `Error::InvalidValue` with the source's message text. |
| `MSSpectrum getSpectrumByNativeId(const std::string&)` | `spectrum_by_native_id(&str) -> Result<MSSpectrum>` | Unfiltered, as in the source. |
| `void setSkipXMLChecks(bool)` | `set_skip_xml_checks(bool)` | **Ported.** It sets the handler's `options_mut().skip_xml_checks`, which reaches the Base64 whitespace strip at `src/format/mzml.rs:2386`. That is the whole effect in the source too: `skip_xml_checks_` is forwarded to `MzMLHandlerHelper::decodeBase64Arrays` and never to XML syntax checking. Remembered across `open_file`, because in the source it is a handler member that `openFile` does not reset. |
| — | `skip_xml_checks() -> bool` | Native getter; the source has none. |
| `PeakFileOptions& getOptions()` | `options_mut() -> &mut PeakFileOptions` | |
| `const PeakFileOptions& getOptions() const` | `options() -> &PeakFileOptions` | |
| `void setOptions(const PeakFileOptions&)` | `set_options(PeakFileOptions)` | Replaces the value and, as in the source, does **not** forward `skip_xml_checks` to the decoder. |
| `OnDiscMSExperiment& operator=(const OnDiscMSExperiment&)` (private) | not ported | Declared private and never defined, because the handler's file streams cannot be copied. `try_clone` is the supported copy. |
| `void loadMetaData_(const std::string&)` (private) | inside `open_file` | `fill_data = false`, no other filter. |
| `MSSpectrum getMetaSpectrumById_(const std::string&)` (private) | inside `spectrum_by_native_id`, over an eagerly built `BTreeMap` | |
| `MSChromatogram getMetaChromatogramById_(const std::string&)` (private) | inside `chromatogram_by_native_id` | |
| `typedef ChromatogramPeak ChromatogramPeakT` (private) | `crate::kernel::ChromatogramPeak` | Private alias, unused by the class body. |
| `typedef Peak1D PeakT` (private) | `crate::kernel::Peak1D` | Private alias, unused by the class body. |
| `std::string filename_` (protected) | `path()` | |
| `Internal::IndexedMzMLHandler indexed_mzml_file_` (protected) | owned `Option<IndexedMzMLHandler>`; `is_indexed()` | `None` reproduces the handler whose `parsing_success_` is false. |
| `std::shared_ptr<PeakMap> meta_ms_experiment_` (protected) | owned `Option<MSExperiment>`; `metadata()` | |
| `std::unordered_map<std::string, Size> chromatograms_native_ids_` (protected) | private `BTreeMap<String, usize>` | Built eagerly with the metadata; first entry wins, as `emplace` does. |
| `std::unordered_map<std::string, Size> spectra_native_ids_` (protected) | private `BTreeMap<String, usize>` | |
| `PeakFileOptions options_` (protected) | private `PeakFileOptions`; `options()` / `options_mut()` / `set_options()` | |
| `typedef OpenMS::OnDiscMSExperiment OnDiscPeakMap` | `pub type OnDiscPeakMap = OnDiscMSExperiment` | The name the class test uses throughout. |

Native additions with no source counterpart: `OnDiscLimits` (and its five
fields), `with_limits`, `open`, `limits()`, `path()`, `is_indexed()`,
`skip_xml_checks()`, `try_clone()` and `load_experiment()`.

`load_experiment` is the loop of `IndexedMzMLFileLoader::store`
(`IndexedMzMLFileLoader.cpp:50-59`) — `getNrSpectra()` through `getSpectrum` and
`getNrChromatograms()` through `getChromatogram` — collected into an
`MSExperiment` instead of fed to a writing consumer.

## Preserved source conventions

- **Filter order and effect.** RT, then MS level, then precursor m/z, each
  returning the metadata spectrum before any peak read; then the record; then the
  m/z and intensity selection (`OnDiscMSExperiment.cpp:136-192`).
- **Half-open ranges.** `DRange<1>::encloses` rejects `value < min` and
  `value >= max` (`DRange.h:152-159`), so a value equal to the maximum is
  outside. `encloses` in this module is that comparison, matching
  `src/format/mzml_load.rs` for a whole-file load.
- **Only the first precursor is tested,** and a spectrum without precursors is
  not tested against the precursor range at all (`:156`).
- **Chromatograms are never skipped before I/O**, and the RT range selects their
  points rather than dropping the record (`:197-227`).
- **The by-native-id fetches consult no option** (`:84`, `:116`). A caller with
  an m/z range configured gets filtered peaks from `spectrum` and unfiltered
  peaks from `spectrum_by_native_id`. Documented below as an inconsistency rather
  than repaired.
- **Unknown native identifier is an error** — `Exception::IllegalArgument` with
  the message `Could not find spectrum with id '<id>'.`, reproduced verbatim as
  `Error::InvalidValue`.
- **The metadata load is never filtered**, so metadata indices stay aligned with
  the index sections; the facade's own options do not reach it.
- **`skipMetaData` consequences**, none of which the header documents:
  `experimental_settings()` and `metadata()` are `None`, `is_sorted_by_rt()` is
  `false`, and the RT, MS-level and precursor filters are not applied at all,
  because the values they test are only known from the metadata (`:169-173`).
- **`empty()` and `size()` ask only about spectra.**
- **Equality ignores the index** and falls back to pointer comparison when either
  metadata is absent (`OnDiscMSExperiment.h:112-124`).
- **`setSkipXMLChecks` is independent of `options_`** in both directions.
- **Duplicate native identifiers:** the first entry wins, as
  `unordered_map::emplace` does.

## Native differences

- **A failed open changes nothing.** The source assigns `filename_` first and
  replaces `meta_ms_experiment_` inside `loadMetaData_`, so a throwing metadata
  load leaves the object naming the new file with the old or a half-built
  metadata. Here the handler, metadata and both caches are built into temporaries
  and committed together.
- **Filtering preserves the record.** The source's m/z or intensity path builds a
  fresh `MSSpectrum` and assigns only the `SpectrumSettings` base
  (`OnDiscMSExperiment.cpp:178-180`), which does not carry retention time, MS
  level, drift time, name or the three data arrays — all of them are lost
  (defect 2). This port selects peaks in place with `MSSpectrum::retain_peaks`,
  so the record keeps its metadata *and* its aligned annotation arrays. The
  chromatogram path is the same, losing name and the data arrays.
- **The identifier caches are built eagerly and replaced with the metadata.** The
  source builds them lazily and only when empty, and `openFile` never clears
  them, so a second `openFile` leaves stale identifiers in place (defect 3).
- **A copy resolves identifiers.** The source's copy constructor omits the two
  caches, and the handler's copy constructor omits its own two, so every
  by-native-id lookup on a copy throws (defect 1). `try_clone` copies them.
- **Out-of-range access is an error.** The source indexes `meta_ms_experiment_`
  with `std::vector::operator[]` and `getChromatogram(id)` with no bound check;
  only the handler's own check stops an out-of-range id, and only for the offset
  vector.
- **No mutable metadata handle.** `getMetaData` hands out a non-const
  `shared_ptr<PeakMap>`, so a caller can mutate the metadata the facade is using
  and silently desynchronise the lazily built identifier maps. `metadata()`
  returns a shared reference; clone it for an owned copy.
- **No unopened-and-unusable state.** `is_indexed()` reports what the source
  keeps only in the handler's `parsing_success_`, and a fetch without an index is
  a clear `Error::InvalidValue` rather than an exception from inside the handler.
- **Exclusive access is enforced.** Every fetch takes `&mut self`, so the
  source's "not thread-safe" `@note` becomes a compile-time property.
- **No threads.** The source parallelises nothing in this class, but its class
  documentation recommends `#pragma omp parallel for firstprivate(ondisc_map)`.
  This port is serial; concurrency is a caller-side decision and costs one extra
  `open` or `try_clone` per reader. The performance gap is stated, not closed.
- **`load_experiment` never dereferences absent metadata.** The source's own
  materialisation loop does (defect 4); here the settings stay default.
- **Ordered maps.** `BTreeMap` instead of `unordered_map`: same first-wins rule,
  reproducible iteration.

## Checked boundaries and evidence

`OnDiscLimits`, all checked before anything is allocated or committed:

| Field | Default | Guards |
|---|---|---|
| `record` | `RecordReadLimits::default()` | handed to the handler: one record's byte range, the cached header, the index |
| `read` | `mzml::ReadOptions::default()` | XML and binary-array ceilings for the metadata load and for each record |
| `max_records` | 1 000 000 | metadata spectra plus chromatograms at open; index spectra plus chromatograms in `load_experiment` |
| `max_native_id_bytes` | 64 MiB | stored identifiers, counted once for the key and once for the record it names |
| `max_materialized_points` | 10 000 000 | peaks plus chromatogram points `load_experiment` holds at once |

The record ceilings are the handler's, so a hostile index still cannot cause an
unbounded read through this facade; `tests/on_disc_experiment.rs` asserts that a
1 KiB `max_record_bytes` rejects the fixture's first spectrum and leaves the
index and metadata intact. `load_experiment` reserves both vectors with
`try_reserve_exact` after the record-count preflight and checks the running point
total before each push, so a rejection leaves this object and any destination
untouched. A refused `open_file` leaves the previous file open and fetchable.

**Evidence tier 3 (source review).** Literals transcribed from
`OnDiscMSExperiment_test.cpp` and its unmodified fixtures: 2 spectra, 1
chromatogram, 19914 and 19800 peaks, 48 chromatogram points, instrument
`LTQ FT` with 1 mass analyzer, native ids
`controllerType=0 controllerNumber=1 scan=1`, `…scan=2` and `TIC`, the m/z window
400–600, the intensity window 1000–1000000, and from `MzMLFile_4_indexed.mzML`
the MS2 at index 1 with precursor m/z 5.55. Scan start times 0.2961 s and
0.4738 s and the peak counts 15 and 10 are read from the same fixtures through
this crate's own metadata load.

**Evidence tier 4 (independently derived).** The resource ceilings, the
atomicity of a failed `open_file`, the half-open endpoints exercised at a
spectrum's own retention time, `load_experiment`, `try_clone` resolving native
identifiers where a source copy cannot, and the assertion that the by-native-id
fetches ignore the options.

No C++ was built or executed and no retained C++ output was used, so nothing here
is tier 1 or tier 2. An oracle driver over `OnDiscMSExperiment` in
`../oracle/drivers/` would be the way to reach tier 1, and is recorded as a
deferral.

## Class-test section coverage

`OnDiscMSExperiment_test.cpp` has 29 `START_SECTION` occurrences: 28 live and one
commented out upstream. All 29 are ported; none are merely mapped.

| Section | Rust test |
|---|---|
| `OnDiscMSExperiment()` | `default_construction_holds_no_file` |
| `~OnDiscMSExperiment()` | `dropping_an_experiment_releases_the_file` |
| `OnDiscMSExperiment(const OnDiscMSExperiment&)` | `a_reopened_clone_sees_the_same_file` |
| `OnDiscMSExperiment(const std::string&)` (commented out at `:53`) | `open_is_the_constructor_from_a_filename` |
| `bool operator==` | `equality_compares_the_file_and_the_metadata` |
| `bool operator!=` | `inequality_is_the_negation_of_equality` |
| `bool openFile(const std::string&, bool)` | `open_file_reports_whether_the_index_parsed` |
| `bool isSortedByRT() const` | `sorted_by_rt_reads_the_metadata_only` |
| `Size size() const` | `size_is_the_indexed_spectrum_count` |
| `bool empty() const` | `empty_asks_only_about_spectra` |
| `Size getNrSpectra() const` | `spectrum_count_matches_the_index` |
| `Size getNrChromatograms() const` | `chromatogram_count_matches_the_index` |
| `getExperimentalSettings() const` | `experimental_settings_come_from_the_metadata` |
| `MSSpectrum operator[](Size)` | `indexing_is_an_alias_for_spectrum` |
| `MSSpectrum getSpectrum(Size)` | `spectrum_by_index_merges_peaks_into_the_metadata` |
| `Interfaces::SpectrumPtr getSpectrumById(Size)` | `spectrum_arrays_are_the_peak_vector` |
| `MSChromatogram getChromatogram(Size)` | `chromatogram_by_index_merges_points_into_the_metadata` |
| `Interfaces::ChromatogramPtr getChromatogramById(Size)` | `chromatogram_arrays_are_the_point_vector` |
| `getChromatogramByNativeId(const std::string&)` | `chromatogram_by_native_id_resolves_through_the_metadata` |
| `getSpectrumByNativeId(const std::string&)` | `spectrum_by_native_id_resolves_through_the_metadata` |
| `PeakFileOptions& getOptions()` | `mutable_options_persist` |
| `const PeakFileOptions& getOptions() const` | `shared_options_read_back_the_same_values` |
| `void setOptions(const PeakFileOptions&)` | `set_options_replaces_the_whole_value` |
| `[EXTRA] m/z range filtering on getSpectrum` | `mz_range_selects_peaks_after_loading` |
| `[EXTRA] intensity range filtering on getSpectrum` | `intensity_range_selects_peaks_after_loading` |
| `[EXTRA] copy constructor copies options` | `a_reopened_clone_copies_the_options` |
| `[EXTRA] RT range filter skips loading peak data` | `rt_range_filter_skips_loading_peak_data` |
| `[EXTRA] MS level filter skips loading peak data` | `ms_level_filter_skips_loading_peak_data` |
| `[EXTRA] precursor m/z range filter skips loading peak data for MS2` | `precursor_mz_range_filter_skips_loading_peak_data` |

Native tests beyond the sections: `a_failed_open_leaves_the_previous_file_in_place`,
`by_native_id_ignores_the_options`, `load_experiment_materialises_every_record`,
`materialisation_ceilings_are_checked`,
`metadata_ceilings_are_checked_before_the_caches_are_built`,
`record_ceilings_are_carried_to_the_handler`,
`skip_xml_checks_reaches_the_decoder_and_survives_reopening`,
`the_metadata_is_loaded_without_peak_data`,
`half_open_range_endpoints_follow_drange_encloses`,
`a_duplicate_native_identifier_keeps_the_first_record`.

Fixture substitution: the failure sections open the upstream non-indexed
`MzMLFile_1.mzML`, which this crate's mzML reader rejects with "binary array
count mismatch" for an unrelated reason. The already-committed upstream
`MzMLFile_2_minimal.mzML` (`tests/data/mzml_upstream_minimal.mzML`) stands in;
like `MzMLFile_1.mzML` it is plain mzML with no `indexListOffset` footer, which
is the only property those sections rely on. Its metadata is empty where
`MzMLFile_1.mzML`'s is not, so the failed-open experiment's metadata *content* is
not exercised.

## Consumers

`docs/core-sdk-coverage.json` records exactly one direct TOPP consumer of this
header: **TICCalculator**. Inside the SDK the header is also included by
`IndexedMzMLFileLoader`, `PeakPickerHiRes`, `SpectrumMetaDataLookup`,
`GNPSMGFFile` and `OnDiscImzMLExperiment`, none of which is ported.

## Source defects found while porting

Reported to the integrator for `OpenMS_CPP_ISSUES.md`; numbering is assigned
there.

1. **A copy cannot resolve native identifiers.**
   `OnDiscMSExperiment.h:97-103` copies `filename_`, `indexed_mzml_file_`,
   `meta_ms_experiment_` and `options_`, but not `spectra_native_ids_` or
   `chromatograms_native_ids_`. Those are lazily rebuilt when empty, so the
   omission is only wasted work — except that copying `indexed_mzml_file_` by
   value invokes `IndexedMzMLHandler`'s own copy constructor
   (`IndexedMzMLHandler.cpp:76-88`), which omits *its* two native-id maps from
   its initialiser list. `getSpectrumByNativeId` and
   `getChromatogramByNativeId` both end in a handler call keyed by identifier,
   so on any copy they throw `Exception::IllegalArgument` for every identifier,
   including ones the file certainly has. The class exists to be copied: the
   header's own `@note` at line 61 recommends
   `#pragma omp parallel for firstprivate(ondisc_map)`, which produces exactly
   these broken per-thread copies.

2. **Filtering by m/z or intensity discards retention time and MS level.**
   `OnDiscMSExperiment.cpp:178-180` builds `MSSpectrum filtered;` and then
   `filtered.SpectrumSettings::operator=(spectrum)`. `SpectrumSettings`
   (`SpectrumSettings.h:172-179`) holds the native id, comment, instrument
   settings, source file, acquisition info, precursors, products and data
   processing — while `retention_time_`, `drift_time_`, `drift_time_unit_`,
   `ms_level_`, `name_` and the float/integer/string data arrays are members of
   `MSSpectrum` itself (`MSSpectrum.h:743-764`). A spectrum fetched with an m/z
   or intensity range therefore comes back with `rt == -1`, `ms_level == 1`, no
   drift time, no name and no annotation arrays, whatever the file said. The
   chromatogram path (`:213-214`) loses `name_` and the three data arrays the
   same way. The class test's own m/z and intensity sections check only the peak
   count and the peak values, and its RT-preservation assertion is on the
   early-return path, so nothing catches this. Proposed fix: copy the whole
   spectrum and erase the rejected peaks, or assign the remaining members
   explicitly.

3. **`openFile` leaves stale native-identifier maps.**
   `OnDiscMSExperiment.cpp:16-25` replaces `filename_` and, via `loadMetaData_`,
   `meta_ms_experiment_`, but never clears `spectra_native_ids_` or
   `chromatograms_native_ids_`; `getMetaSpectrumById_` (`:100`) and
   `getMetaChromatogramById_` (`:68`) rebuild them only when they are empty.
   Opening a second file on an object that has already resolved one identifier
   therefore keeps the first file's identifier-to-index mapping, so
   `getSpectrumByNativeId` silently returns the wrong spectrum's metadata, or
   throws for an identifier the new file does have. The underlying handler
   compounds this by appending rather than replacing its offset vectors
   (`IndexedMzMLHandler.cpp:37-46`), so counts are the sum of both files.

4. **`IndexedMzMLFileLoader::store` dereferences a null metadata pointer.**
   `IndexedMzMLFileLoader.cpp:47` calls
   `consumer.setExperimentalSettings(*exp.getExperimentalSettings().get())` with
   no null check. `getExperimentalSettings` returns a null `shared_ptr` whenever
   the experiment was opened with `skipMetaData == true`
   (`OnDiscMSExperiment.h:170-173`), a state the header documents nowhere, so
   storing such an experiment dereferences null. The loader's `load` always
   passes `skipMetaData == false`, but the two calls are independent public API.

5. **The by-native-id fetches ignore the configured `PeakFileOptions`.**
   `getSpectrum` (`:130`) and `getChromatogram` (`:197`) both consult `options_`;
   `getSpectrumByNativeId` (`:116`) and `getChromatogramByNativeId` (`:84`) never
   do. Two fetches of the same record therefore disagree, and the unfiltered one
   is the one whose name suggests a targeted lookup. Reported as an
   inconsistency, not a crash: it may be deliberate, but neither the header nor
   the implementation says so. This port preserves the behaviour and documents
   it at the two methods.

## Deferrals

- Tier 1 evidence would need an oracle driver over `OnDiscMSExperiment` in
  `../oracle/drivers/`, printing full-precision peak arrays for each record of
  `IndexedmzMLFile_1.mzML` plus the five filtered variants. Not built here.
- `IndexedMzMLFileLoader` is not ported; only the materialisation loop of its
  `store`, as `load_experiment`. Its own `store` overloads, its options and its
  `load` are out of scope.
- The `Interfaces::Spectrum` / `Interfaces::Chromatogram` pointer overloads are
  not ported; the crate has no OpenSWATH `Interfaces` layer.
- `OnDiscImzMLExperiment` (the imzML sibling) is not ported.
- `docs/core-sdk-reviewed-apis.json` records
  `FORMAT/HANDLERS/IndexedMzMLHandler.h` with "Not ported: … setSkipXMLChecks",
  which is wrong — that member maps to `options_mut().skip_xml_checks` and
  reaches the decoder. Correcting that scope string belongs to the integrator,
  who owns the ledger.

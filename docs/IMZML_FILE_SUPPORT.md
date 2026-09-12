# imzML file adapter — `FORMAT/ImzMLFile.h`

Stage 3 of the imzML family: the user-facing file API that ties the reader, the
writer and the on-disc facade together, plus the family's entire class-test
suite.

- C++ source revision: `bc9cc12514c768385ce121d6ca4bb710fe1983c4`
- Header / implementation: `src/openms/include/OpenMS/FORMAT/ImzMLFile.h` (283
  lines), `src/openms/source/FORMAT/ImzMLFile.cpp` (645 lines)
- Class tests: `ImzMLFile_test.cpp` (1,769 lines, 41 sections),
  `ImzMLFile_all_modes_test.cpp` (258 lines, 9 sections)
- Rust: `src/format/imzml_file.rs`
- Tests: `tests/imzml_file.rs` (42 section tests for 41 sections + 15 native + 13 for the float oracle, 70 in all), `tests/imzml_all_modes.rs` (9 sections)
- Provenance: `tests/data/imzml_file_provenance.json`
- Fixtures: reused unmodified from stage 1 —
  `tests/data/ImzMLFile_1_Example_Continuous.{imzML,ibd}`,
  `tests/data/ImzMLFile_2_Example_Processed.{imzML,ibd}`

The stages this one builds on: `docs/IMZML_HANDLER_SUPPORT.md` (index and
`.ibd` reads), `docs/IMZML_WRITER_SUPPORT.md` (the writer),
`docs/ON_DISC_IMZML_SUPPORT.md` (the grid, the regions, the ion image and the
on-disc facade).

## What imzML is

An imzML 1.1.0 dataset is **two files**. The `.imzML` is mzML XML carrying the
metadata and, per spectrum, the IMS CV params giving the byte offset
(`IMS:1000102`), element count (`IMS:1000103`) and stored byte length
(`IMS:1000104`) of that spectrum's m/z and intensity arrays inside the
companion `.ibd`, together with the pixel's image coordinates (`IMS:1000050`
position x, `IMS:1000051` position y, `IMS:1000052` position z). The `.ibd`
begins with a 16-byte UUID that must equal the `IMS:1000080`
universally-unique-identifier param in the XML. Two storage modes exist:
*continuous* (`IMS:1000030`, one shared m/z array stored once) and *processed*
(`IMS:1000031`, a private m/z array per spectrum). Both are exercised by the
upstream fixtures and by every round trip here.

---

## API mapping

Every public member declared in `ImzMLFile.h`, in header order. The private
members and the file-static helpers of `ImzMLFile.cpp` follow, because they
carry behaviour a caller can observe.

### Declared public members

| C++ | Rust | Notes |
|---|---|---|
| `ImzMLFile()` | `ImzMLFile::new`, `impl Default for ImzMLFile` | The source constructor's only work is registering `mzML_1_10.xsd` / `"1.1.0"` with its `XMLFile` base (`ImzMLFile.cpp:208`); the Rust equivalent is which schema `is_valid` selects. |
| `~ImzMLFile() override = default` | not ported: nothing to release | The source destructor is defaulted. `ImzMLFile` owns no handle; the `.ibd` handle lives in `ImzMLHandler` and is dropped when a load returns. |
| `PeakFileOptions& getOptions()` | `ImzMLFile::options_mut` | |
| `const PeakFileOptions& getOptions() const` | `ImzMLFile::options` | |
| `void setOptions(const PeakFileOptions&)` | `ImzMLFile::set_options` | Takes the value; the source copies. |
| `void load(const std::string&, MSImagingExperiment&)` | `ImzMLFile::load`, `ImzMLFile::load_with_ibd` | Out-parameter becomes the return: `Result<(ImagingExperiment, LoadReport)>`. The `_with_ibd` variant exposes the `ibd_path_override` that the source reaches only through `loadSpectraIndex`. See **the mzML metadata gap** below for the one behavioural divergence. |
| `static void buildImagingGeometry(const MSExperiment&, MSImagingGeometry&)` | `build_imaging_geometry_from_experiment` | Returns `Result<(ImagingGeometry, MetaGeometryReport)>`. The source has no report and logs instead. |
| `static void buildImagingGeometry(const std::vector<ImzMLSpectrumIndex>&, const ImzMLMeta&, MSImagingGeometry&)` | `build_imaging_geometry_from_index`, a re-export of `crate::kernel::on_disc_imzml_experiment::build_imaging_geometry` | Ported by the on-disc stage, because the on-disc reader needs it too; re-exported here under the name this header declares. One function, two paths. |
| `void load(const std::string&, Interfaces::IMSDataConsumer&)` | `ImzMLFile::load_into_consumer`, `ImzMLFile::load_into_consumer_with_ibd` | Consumer trait is `crate::interfaces::MSDataConsumer`. Returns a `LoadReport`; the source returns `void`. |
| `void loadSpectraIndex(const std::string&, ImzMLMeta&, std::vector<ImzMLSpectrumIndex>&, const std::string& ibd_path = "")` | `ImzMLFile::load_spectra_index`, `ImzMLFile::load_spectra_index_with_ibd`, `ImzMLFile::load_spectra_index_checked` | The two out-parameters are the two fields of the returned `ImzMLIndex`. The defaulted fourth parameter becomes the `_with_ibd` overload, per the crate's rule against defaulted arguments. `_checked` additionally returns the UUID verdict that the source logs. |
| `bool isValid(const std::string&, std::ostream&)` | `ImzMLFile::is_valid`, gated on the existing `mzml-schema` feature | Returns `Result<SchemaValidationReport>`: `report.is_valid()` is the source's `bool`, `report.diagnostics` is what the source wrote to `os`. Feature-gated because the libxml2 validator is optional in this crate and Xerces is not in OpenMS. |
| `void store(const std::string&, const MSExperiment&) const` | `ImzMLFile::store` | Forwards to `imzml_writer::store_with_options` with this adapter's options, as the source forwards to `Internal::ImzMLWriter::store`. Returns the writer's `StoreReport`. |
| `void store(const std::string&, const MSImagingExperiment&) const` | `ImzMLFile::store_imaging` | Synthesises `imzml:x` / `imzml:y` / `imzml:z` and the grid metadata from the geometry onto a clone and reuses the writer, exactly as the source does (`ImzMLFile.cpp:500-532`). |

### Private members and file-static helpers

| C++ | Rust | Notes |
|---|---|---|
| `void loadImpl_(filename, consumer, meta_exp, out_meta, out_index, index_only, ibd_path_override)` | private `ImzMLFile::load_experiment_with_ibd` + `decode_all` / `stream` / `decode_one` | The source's one function with six switches is three private methods here, because the three public entry points want three different products. The switch semantics are preserved: `index_only` sets `FillData(false)` and `AlwaysAppendData(false)`, and a null consumer sets `AlwaysAppendData(true)` (`ImzMLFile.cpp:550-559`). |
| `static std::string inferIbdPath_(const std::string&)` | `imzml_handler::infer_ibd_path`, ported by stage 1 | Case-insensitive `.imzml` suffix replaced with `.ibd`, any other name gains `.ibd` (`ImzMLFile.cpp:639`). Checked in `the_ibd_sibling_is_inferred_case_insensitively`. |
| `PeakFileOptions options_` | private `ImzMLFile::options` field | |
| file-static `uuidStringToBytes_` | `imzml_handler::uuid_bytes` | Strips `-`, `{` and `}` and requires exactly 32 hex digits (`ImzMLFile.cpp:45`). |
| file-static `bytesToHex_` | private, inside `ImzMLHandler::uuid_status` | Lower-case hex; surfaces as the two strings of `UuidStatus::Mismatch`. |
| file-static `verifyIbdUuid_` | `ImzMLHandler::uuid_status`, surfaced as `LoadReport::uuid` | Source logs a warning for every non-matching outcome and loads anyway (`ImzMLFile.cpp:88-123`); this returns the verdict. |
| file-static `attachImzMLMeta_` | `attach_dataset_meta`, public | Public here because it is the documented contract between a load and a later store, and the round trip is testable only if a caller can invoke it (`ImzMLFile.cpp:125-184`). |
| file-static `struct XercesPlatformGuard` | not ported: no counterpart | Xerces needs a process-global `Initialize()` and must deliberately never be `Terminate()`d from a destructor, or MSVC fastfails (`ImzMLFile.cpp:186-205`). `quick_xml` has no global state, so the whole hazard is absent. The source's `@note` about `Terminate()` is therefore neutralised rather than carried. |

### Surface inherited through public base classes

`class ImzMLFile : public Internal::XMLFile, public ProgressLogger` makes both
bases' public members public members of `ImzMLFile`. Rust has no
implementation inheritance, so they are mapped as members or recorded as absent.

| C++ (inherited) | Rust | Notes |
|---|---|---|
| `XMLFile::isValid(filename, os)` | `ImzMLFile::is_valid` | The derived class re-declares and forwards to it; mapped once, above. |
| `XMLFile::getVersion() const` | not ported: constant `"1.1.0"` | The value is fixed by the constructor and no caller of `ImzMLFile` reads it. Recorded rather than exposed as a getter that can only return one string. |
| `XMLFile::~XMLFile()` | not ported | |
| `ProgressLogger::setLogType` / `getLogType` | `ImzMLFile::set_log_type` / `ImzMLFile::log_type` | `loadImpl_` copies the log type onto the logger it hands the handler; `store` onto the writer's. Checked in `the_progress_log_type_round_trips_and_a_logged_load_still_works`. |
| `ProgressLogger::startProgress` / `setProgress` / `endProgress` / `nextProgress` | not ported on `ImzMLFile`: driven internally | The ported `crate::concept::progress_logger::ProgressLogger` has all four. Exposing them on the file adapter would let a caller desynchronise the range this adapter owns; the source's `mutable` log type and `static int recursion_depth_` are exactly how it leaks that state. |
| `ProgressLogger::setLogger` | not ported on `ImzMLFile` | As above; `ProgressLogger::set_logger` exists on the logger type. |
| `ProgressLogger::LogType` | `crate::concept::progress_logger::ProgressLogType` | `CMD` / `GUI` / `NONE` become `Cmd` / `Gui` / `None`. |
| `ProgressLogger(const ProgressLogger&)`, `operator=` | `#[derive(Clone)]` on `ImzMLFile` | |

### `MSImagingExperiment`, which is a different header

`ImzMLFile::load` and one `store` overload are typed on
`IMAGING/MSImagingExperiment.h`. That header has its own class test and is
**not** this package's to port; leaving those two members unported was the
alternative, and it was worse. `ImagingExperiment` in this module reproduces
the surface those two members need, and nothing else is claimed for the header.
For the record, of `MSImagingExperiment`'s public members:

| C++ | Covered here |
|---|---|
| `explicit MSImagingExperiment(MSExperiment)` | `impl From<MSExperiment> for ImagingExperiment` |
| `MSImagingExperiment& operator=(MSExperiment)` | `ImagingExperiment::set_ms_experiment_and_clear` — the overload that clears the geometry |
| `getMSExperiment()` / `const` | `ms_experiment_mut` / `ms_experiment` |
| `setMSExperiment(MSExperiment)` | `set_ms_experiment` — keeps the geometry |
| `getGeometry()` / `const` | `geometry_mut` / `geometry` |
| `setGeometry(MSImagingGeometry)` | `set_geometry` |
| `getNumberOfPixels()` | `number_of_pixels` |
| `getNumberOfSpectra()` | `number_of_spectra` |
| `hasPixel(UInt, UInt)` | `has_pixel` |
| `getSpectrum(UInt, UInt)` / `const` | `spectrum_mut` / `spectrum` |
| `extractIonImage(double, double)` | `extract_ion_image` |
| `extractIonImage(double, double, Size)` | `extract_ion_image_in_region` |
| `validate()` | `validate` |
| `getRegionSpectrumIndices(Size)` | `region_spectrum_indices` |

That is the whole declared surface, and a later IMAGING package owns the header
itself along with `MSImagingGeometry.h`, `MSImagingRegion.h`, `IonImage.h` and
`IonImageExtraction.h` (`ImagingGeometry`, `ImagingRegion` and `IonImage` were
reproduced to the same standard by the on-disc stage; see
`docs/ON_DISC_IMZML_SUPPORT.md`). No claim is made on any of them.

---

## Preserved source conventions

1. **Two files, one adapter.** Every read path opens the companion `.ibd` — so
   a missing one is an error, not a metadata-only load — and records its path
   in `ImzMLMeta::ibd_file_path` (`ImzMLFile.cpp:568-570`).
2. **An explicit `.ibd` beats the inferred sibling**, and it reaches the index
   load *and* the UUID check, not just the reads (`ImzMLFile.cpp:565-569`).
3. **The `.ibd` UUID check is advisory.** A mismatch, a missing or unparsable
   `IMS:1000080`, or an `.ibd` shorter than 16 bytes is reported and the dataset
   still loads, deliberately, so legacy and non-conformant writers stay readable
   (`ImzMLFile.cpp:82-87`).
4. **The grid comes from the parsed index, not from meta values.**
   `load(MSImagingExperiment&)` uses the index builder, the same one the on-disc
   reader uses, so the in-memory and on-disc geometries cannot diverge
   (`ImzMLFile.cpp:229-233`). `tests/imzml_all_modes.rs` asserts the two agree
   pixel for pixel, and `the_geometry_builders_agree_on_the_upstream_fixtures`
   asserts the meta-value builder agrees with both.
5. **Coordinates are 1-based in the file and 0-based in the grid**, in both
   directions: a load subtracts one (`ImzMLFile.cpp:297`), the geometry-driven
   store adds one (`ImzMLFile.cpp:518`).
6. **Only the `z == 1` plane is placed.** The grid is two-dimensional; other
   planes stay reachable by index (`ImzMLFile.cpp:287-295`, `:397-400`).
7. **Four ways a spectrum can fail to reach the grid, all non-fatal**: no
   coordinates at all (skipped silently), a coordinate below 1, a pixel outside
   the declared `IMS:1000042` x `IMS:1000043` grid, and a pixel an earlier
   spectrum already claimed. Only the first spectrum per pixel is mapped; every
   spectrum stays reachable by index (`ImzMLFile.cpp:270-322`).
8. **An undeclared grid is derived from the observed maxima**, and a declared
   one larger than any pixel survives (`ImzMLFile.cpp:349-361`).
9. **The two geometry builders disagree on two details, and both are
   reproduced.** The meta-value builder tests the coordinates before `z`
   (`ImzMLFile.cpp:277` then `:287`), so a spectrum at `(0, 0, 2)` counts as a
   non-positive coordinate; the index builder tests `z` first (`:397` then
   `:401`), so the same spectrum counts as another plane. And the meta-value
   builder copies the pixel size when both keys merely *exist* (`:363`) where the
   index builder copies it only when both are strictly positive (`:472`), so the
   two can disagree about the pixel size of one dataset; neither C++
   `setPixelSize` nor the ported `set_pixel_size` rejects a zero, and only the
   latter rejects a non-finite one. The ordering difference is asserted in
   `the_meta_value_geometry_builder_reports_every_kind_of_unplaceable_spectrum`.
10. **A duplicate pixel is storable again.** The writer writes both spectra with
    their coordinates and only the reader's grid drops the second, so a dataset
    that loads can be stored (`ImzMLFile.h`'s `store` documentation, and
    `duplicate_pixel_coordinates_survive_a_store_and_a_reload`).
11. **Streaming retains nothing** and delivers only after the `spectrumList`
    section is parsed, never per `startElement` (`ImzMLFile.cpp:478-482` with
    the header's note).
12. **`FillData(false)` is what the index-only path is built on**
    (`ImzMLFile.cpp:551-554`); honoured for any load here, so a caller can have
    the coordinates without the arrays.
13. **The dataset metadata is mirrored onto the experiment** with the exact
    presence rules of `attachImzMLMeta_`: grid counts and `.ibd` path always,
    pixel sizes and extents only when positive, identifier, checksums, array
    data types and acquisition-geometry terms only when the file declared them
    (`ImzMLFile.cpp:125-184`). A key present therefore means the file said so.
14. **`isValid` validates against the unmodified mzML 1.1.0 schema.** imzML adds
    CV terms, never elements, which is why the source registers no imzML schema
    of its own (`ImzMLFile.cpp:209`, `:490`).

---

## Native differences

### The mzML metadata gap — the one behavioural divergence

Source `ImzMLFile` parses an `.imzML` with `Internal::ImzMLHandler`, which
**derives from `MzMLHandler`** and intercepts only the IMS terms its base does
not know. Everything else — the instrument configuration, the software list,
the data processing, and per spectrum `MS:1000016` scan start time and
`MS:1000511` ms level — is the base class's work.

This port cannot reach that base for an imzML file, and the reason is in
`src/format/mzml.rs`, which this package does not own. Measured in this
worktree with a probe against the merged tree:

| Input | `mzml::read` |
|---|---|
| `ImzMLFile_1_Example_Continuous.imzML` | `parse error: unresolved softwareRef` — the fixture puts `<softwareList>` *after* `<dataProcessingList>`, so the single-pass reference resolution sees a forward reference |
| `ImzMLFile_2_Example_Processed.imzML` | `unsupported: only mzML 1.1 is supported` — the fixture declares `version="1.1"` where the reader requires a `1.1.` prefix; behind that, its `ISO-8859-1` bytes are refused as "requires transcoding" |
| any imzML this port writes | `unsupported: binary array CV IMS:1000101` — the external-data param is the one thing an mzML reader cannot be asked to interpret |

The third row is the decisive one: it is not a fixture defect, it is what makes
an imzML an imzML. `read_metadata` clears the third obstacle but not the first
two, and it maps `MS:1000031` "instrument model" carried as a param value to
nothing, which is the same as C++ (OpenMS writes a model as a child term).

**Consequence.** After `ImzMLFile::load`, an experiment carries its peaks,
pixel coordinates, auxiliary float arrays, native identifiers and the whole
`imzml:*` dataset metadata — and carries the OpenMS unset defaults for
retention time (`-1.0`) and MS level (`1`), an empty instrument and no data
processing. `store` writes RT and MS level correctly from whatever the
experiment holds, so a store is unaffected; only a load cannot recover them. The
upstream section `void store metadata round-trip` asserts
`reloaded[0].getRT() == 12.34`; that half is checked here against the written
document instead, and the test says so in place.

**Why not work around it.** Two options were rejected. Re-parsing the `.imzML`
in this module for RT and MS level would be a second XML parser next to
`mzml.rs` and `imzml_handler.rs`. Attempting `mzml::read` and silently falling
back would degrade unpredictably and differ per file. The gap is therefore
recorded, proven by a probe, and left for whoever relaxes the three
strictnesses in `src/format/mzml.rs`.

### Other differences

1. **Reports instead of a log.** `LoadReport`, `MetaGeometryReport` and the
   writer's `StoreReport` carry what the source sends to `OPENMS_LOG_WARN`: the
   `.ibd` UUID verdict, every spectrum that missed the grid, every skipped
   auxiliary array, the filtered-out count. The policy is unchanged — none of
   them stops a load — but a caller that wants strict conformance can now
   enforce it. `LoadReport::is_clean` and `MetaGeometryReport::is_clean`
   summarise.
2. **Filters whose input is not parsed are refused, not applied.** A
   `PeakFileOptions` carrying a retention-time range, an MS-level selection or a
   precursor m/z range makes `load` and `load_into_consumer` return
   `Error::Unsupported`. The source applies all three inside `MzMLHandler`;
   here every spectrum would carry the unset defaults and such a filter would
   discard the whole dataset. Refusing is the crate's rule for a silently lossy
   operation. `sort_spectra_by_mz`, the m/z and intensity peak filters,
   `metadata_only` and `fill_data` *are* honoured, through the writer's
   `apply_store_options`, which is the same `applyStoreOptions_` port and the
   same `DRange<1>::encloses` endpoints.
3. **Sortedness is checked, not assumed.** `MSImagingExperiment::extractIonImage`
   documents sorted peaks as a precondition — "Phase 1 callers must ensure it
   manually" — and would return a wrong sum otherwise.
   `ImagingExperiment::extract_ion_image` returns `Error::UnsortedData`, because
   the check is one pass and the alternative is an ion image that is quietly
   wrong. Loads sort by default, so this only fires for a hand-built
   experiment.
4. **Bounded loads.** `ImzMLLoadLimits` sums the parsed index's declared element
   counts and refuses the dataset before the first array is read. The source has
   no ceiling on a load at all.
5. **`AlwaysAppendData` has no counterpart.** Every load here returns a fresh
   experiment, so there is nothing to append to; the source's flag exists
   because `loadImpl_` writes into a caller-owned `MSExperiment&`.
6. **A defaulted argument becomes an overload.** `loadSpectraIndex`'s
   `ibd_path = ""` is `load_spectra_index` / `load_spectra_index_with_ibd`.
7. **`Exception::InvalidParameter` has no variant.** A declared continuous mode
   the spectra cannot support is `Error::InvalidValue`.
8. **A consumer may stop the load.** `MSDataConsumer::consume_spectrum` returns
   `ControlFlow`; `Interfaces::IMSDataConsumer::consumeSpectrum` returns `void`
   and cannot. Recorded in `LoadReport::stopped_early`.
9. **Serial.** The source parallelises the array decode inside `ImzMLHandler`
   with `#pragma omp parallel for`. Nothing here starts a thread, so a large
   image decodes on one core. The performance gap is stated rather than
   discovered.
10. **No Xerces lifecycle.** See `XercesPlatformGuard` above.

---

## Checked boundaries and evidence

| Boundary | Rust outcome |
|---|---|
| missing `.imzML` or `.ibd`, on any load path | `Error::Io` (source `Exception::FileNotFound`) |
| malformed XML, unparsable IMS value, spectrum with no pixel coordinate | `Error::Parse` (source `Exception::ParseError`) |
| compressed external m/z, intensity or auxiliary array | `Error::Unsupported` (source `Exception::ParseError`) |
| declared peaks across the dataset > `ImzMLLoadLimits::max_loaded_peaks` | `Error::InvalidValue`, before the first array is read |
| declared auxiliary values > `max_loaded_float_values` | `Error::InvalidValue`, same preflight |
| every per-array and per-file ceiling of `ImzMLReadLimits` | `Error::InvalidValue`, per `docs/IMZML_HANDLER_SUPPORT.md` |
| declared or derived grid above `MAX_IMAGE_PIXELS` | `Error::InvalidValue` |
| `imzml:x` / `imzml:y` / `imzml:z` / `imzml:max_count_*` of the wrong type or outside `u32` | `Error::InvalidValue` (source `Exception::ConversionError`, or a silent `static_cast` wrap) |
| a geometry pixel referencing a spectrum index at or above the spectrum count | `Error::InvalidValue`, from `validate`, `spectrum`, `store_imaging` and both extractions (source `Exception::InvalidValue`) |
| no pixel at a requested coordinate | `Error::InvalidValue` (source `Exception::ElementNotFound`) |
| unknown region id | `Error::InvalidValue` (source `Exception::ElementNotFound`) |
| negative or non-finite `mz` / `tolerance_ppm` | `Error::InvalidValue` (source `Exception::InvalidValue`) |
| a referenced spectrum not sorted by m/z, on the in-memory extraction | `Error::UnsortedData` (source: unchecked precondition) |
| RT, MS-level or precursor filter in the load options | `Error::Unsupported` (native) |
| a store with no spectra, or a spectrum without `imzml:x` / `imzml:y` | `Error::MissingInformation` (source `Exception::MissingInformation`) |
| a declared continuous mode the spectra cannot support | `Error::InvalidValue` (source `Exception::InvalidParameter`) |
| a 1-based coordinate derived from a grid pixel overflowing `u32` | `Error::InvalidValue` (native) |

**Evidence tier 3, source review**, with two independent cross-checks.

The transcribed literals of both suites are listed in the two test files'
module documentation. No C++ was built or executed and no C++ output was
retained, so **this is not a tier 1 differential**.

The two checks that do not depend on the suite are stronger than a
transcription:

- **Every round trip closes the loop through the `.ibd` offsets.** A dataset is
  written and read back through the stage-1 reader, and the decoded m/z,
  intensities, auxiliary arrays and pixel coordinates are compared against what
  went in — for every spectrum, not only the first. In processed mode the
  per-spectrum offsets are additionally asserted to be distinct and increasing,
  which is what distinguishes the two modes in the binary file.
- **Every ion image is compared against a sum computed directly from the
  decoded peaks**, not through either extraction path. That is the same
  independent cross-check the C++ sections make between their in-memory and
  on-disc paths, and it is made here against both.
- **The three access paths are compared bit for bit.** `tests/imzml_all_modes.rs`
  fetches pixel (2,3) through the full load, the on-disc reader and the RAM
  lookup and compares the `f64`/`f32` bit patterns, where the upstream section
  compares only the peak counts.

The resource ceilings, the refusal of an unparseable filter, the `fill_data`
load, the UUID-mismatch path and the error-variant choices are independently
derived (tier 4): no upstream fixture reaches them.

Numeric comparisons use a faithful port of `TEST_REAL_SIMILAR` —
`ClassTest::isRealSimilar`, similar when the absolute difference is within
`1e-5` **or** the ratio within `1 + 1e-5` — because the suite's decimal
literals are roundings of float32-stored values. Comparisons between two
computations inside this crate use `1e-12` relative instead, since they must
agree to f64 accumulation error and not to a reporting tolerance.

---

## Class-test section accounting

**All 50 sections of the family's two suites are ported.** None is mapped and
none is unaccounted for. A section may exercise a header other than
`ImzMLFile.h`; the *Covers* column names the header it actually tests, which is
what the ledger entries credit.

### `ImzMLFile_test.cpp` — 41 sections

| # | Line | Upstream section | Covers | Rust test in `tests/imzml_file.rs` |
|---|---|---|---|---|
| 1 | 370 | `void load(filename, MSExperiment&)` | ImzMLFile | `the_continuous_fixture_loads_with_its_first_pixel_at_one_one` |
| 2 | 389 | `OnDiscImzMLExperiment random access` | OnDiscImzMLExperiment | `the_on_disc_reader_serves_the_same_pixel_by_index_and_by_coordinate` |
| 3 | 409 | `const MSImagingGeometry& getGeometry() const` | OnDiscImzMLExperiment | `the_on_disc_grid_matches_the_dataset_and_maps_every_spectrum` |
| 4 | 429 | `IonImage extractIonImage(mz, tolerance_ppm)` | ImzMLFile + OnDiscImzMLExperiment | `the_in_memory_and_on_disc_whole_image_extractions_agree` |
| 5 | 476 | `[EXTRA] open honours an explicit .ibd path override` | OnDiscImzMLExperiment + ImzMLFile | `an_explicit_ibd_override_reaches_the_index_load_and_the_uuid_check` |
| 6 | 502 | `void load(filename, IMSDataConsumer&)` | ImzMLFile | `a_streaming_load_delivers_every_spectrum_and_retains_none` |
| 7 | 520 | `void load Example_Processed imzML` | ImzMLFile | `the_processed_fixture_loads` |
| 8 | 532 | `dataset metadata mirrored on MSExperiment after load` | ImzMLFile + ImzMLHandler | `the_dataset_metadata_of_the_continuous_fixture_is_mirrored_on_the_experiment`, `..._of_the_processed_fixture_...` (the section's two blocks) |
| 9 | 634 | `void load(filename, MSImagingExperiment&)` | ImzMLFile | `the_imaging_load_gives_pixel_random_access_over_the_whole_grid` |
| 10 | 657 | `void buildImagingGeometry(exp, geom)` | ImzMLFile | `the_meta_value_geometry_builder_places_every_fixture_spectrum` |
| 11 | 672 | `static void buildImagingGeometry(index, meta, geom)` | ImzMLFile + OnDiscImzMLExperiment | `the_index_geometry_builder_needs_no_experiment_and_no_meta_values` |
| 12 | 703 | `void store round-trip continuous imzML` | ImzMLWriter + ImzMLFile | `a_continuous_dataset_round_trips_through_store_and_load` |
| 13 | 729 | `void store(filename, const MSImagingExperiment&)` | ImzMLFile | `the_imaging_store_takes_its_coordinates_from_the_geometry` |
| 14 | 774 | `void store round-trip processed imzML` | ImzMLWriter + ImzMLFile | `a_processed_dataset_round_trips_through_store_and_load` |
| 15 | 796 | `void store round-trip FloatDataArray ion mobility and non-standard` | ImzMLWriter + ImzMLHandler + OnDiscImzMLExperiment | `an_ion_mobility_and_a_free_text_float_array_round_trip_on_both_paths` |
| 16 | 921 | `bool isValid(filename, os)` | ImzMLFile | `a_stored_imzml_passes_the_mzml_schema` |
| 17 | 936 | `void store rejects missing pixel coordinates` | ImzMLWriter | `a_store_of_a_spectrum_without_pixel_coordinates_is_refused` |
| 18 | 951 | `void store tolerates duplicate pixel coordinates by default` | ImzMLWriter + ImzMLFile | `duplicate_pixel_coordinates_survive_a_store_and_a_reload` |
| 19 | 1006 | `void store rejects incompatible continuous mode` | ImzMLWriter | `a_declared_continuous_mode_the_spectra_cannot_support_is_refused` |
| 20 | 1021 | `void buildImagingGeometry tolerates duplicate pixels by default` | ImzMLFile | `the_meta_value_geometry_builder_keeps_the_first_spectrum_of_a_shared_pixel` |
| 21 | 1038 | `void load and OnDisc reject zlib-compressed external m/z and intensity` | ImzMLHandler + OnDiscImzMLExperiment | `zlib_compressed_external_peak_arrays_are_refused_on_both_paths` |
| 22 | 1064 | `void store writes a single allowed unit for PSI-MS aux arrays` | ImzMLWriter | `a_psi_term_with_one_allowed_unit_gets_it_and_a_term_with_two_gets_none` |
| 23 | 1107 | `OnDiscImzMLExperiment rejects compressed zero-length aux arrays` | OnDiscImzMLExperiment + ImzMLHandler | `a_compressed_zero_length_auxiliary_array_is_still_refused` |
| 24 | 1135 | `void load rejects spectrum missing pixel coordinates` | ImzMLHandler + ImzMLFile | `a_spectrum_declaring_no_pixel_coordinate_fails_the_load` |
| 25 | 1157 | `void load tolerates duplicate pixel coordinates by default` | ImzMLFile | `a_duplicated_pixel_coordinate_loads_with_only_the_first_in_the_grid` |
| 26 | 1185 | `OnDiscImzMLExperiment tolerates duplicate pixel coordinates by default` | OnDiscImzMLExperiment | `the_on_disc_reader_tolerates_a_duplicated_pixel_coordinate` |
| 27 | 1210 | `void store applies PeakFileOptions m/z range filter` | ImzMLWriter + ImzMLFile | `the_store_mz_range_filter_shrinks_every_spectrum` |
| 28 | 1240 | `void store honors PeakFileOptions binary precision` | ImzMLWriter + ImzMLFile | `the_stored_binary_precision_follows_the_peak_file_options` |
| 29 | 1264 | `void store metadata round-trip` | ImzMLWriter + ImzMLFile | `the_acquisition_metadata_round_trips_through_store_and_load` — the RT and MS-level assertions are checked at the written document; see **the mzML metadata gap** |
| 30 | 1304 | `IonImage extractIonImage(mz, tolerance_ppm, region_id)` | ImzMLFile + OnDiscImzMLExperiment | `the_in_memory_and_on_disc_region_extractions_agree` |
| 31 | 1354 | `void load sorts external peaks when getSortSpectraByMZ is true` | ImzMLFile + OnDiscImzMLExperiment | `a_load_sorts_external_peaks_and_moves_the_annotations_with_them` |
| 32 | 1407 | `void load and OnDisc reject numpress-compressed external m/z` | ImzMLHandler + OnDiscImzMLExperiment | `numpress_compressed_external_peak_arrays_are_refused_on_both_paths` |
| 33 | 1430 | `void load skips one bad aux length and still returns all spectra` | ImzMLHandler + OnDiscImzMLExperiment | `one_auxiliary_array_of_the_wrong_length_does_not_abort_the_load` |
| 34 | 1469 | `void load and OnDisc drop zero-length aux without a ghost IM array` | ImzMLHandler + OnDiscImzMLExperiment | `a_zero_length_auxiliary_array_leaves_no_ghost_ion_mobility_array` |
| 35 | 1503 | `void store skips unnamed FloatDataArray` | ImzMLWriter | `an_unnamed_float_array_is_not_written_at_all` |
| 36 | 1533 | `void load drops phantom IntegerDataArray for integer-typed aux` | ImzMLHandler + OnDiscImzMLExperiment | `an_integer_typed_auxiliary_array_becomes_a_float_array_and_no_integer_array` |
| 37 | 1565 | `void load and OnDisc drop inline aux when peaks are external` | ImzMLHandler + OnDiscImzMLExperiment | `an_inline_auxiliary_array_is_dropped_when_the_peaks_are_external` |
| 38 | 1597 | `void store writes MS:1000576 on mz and intensity referenceableParamGroups` | ImzMLWriter | `both_peak_reference_groups_declare_no_compression` |
| 39 | 1626 | `void store skips FloatDataArrays named after the peak arrays` | ImzMLWriter + ImzMLHandler | `float_arrays_named_after_the_peak_arrays_are_skipped` |
| 40 | 1680 | `void store drops integer and string data arrays` | ImzMLWriter | `integer_and_string_data_arrays_are_dropped_and_the_float_array_survives` |
| 41 | 1728 | `void load and OnDisc skip aux array without a supported binary data type` | ImzMLHandler + OnDiscImzMLExperiment | `an_auxiliary_array_without_a_decodable_type_is_skipped_not_fatal` |

Section 8 is one `START_SECTION` containing two independent blocks, one per
fixture; it is ported as two tests so a failure names the fixture.

### `ImzMLFile_all_modes_test.cpp` — 9 sections

| # | Line | Upstream section | Covers | Rust test in `tests/imzml_all_modes.rs` |
|---|---|---|---|---|
| 1 | 65 | `mode 1: full load into MSExperiment` | ImzMLFile | `mode_1_a_full_load_gives_every_spectrum_and_the_dataset_metadata` |
| 2 | 84 | `mode 2: FileHandler and FileTypes` | FileHandler / FileTypes | `mode_2_the_file_handler_recognises_imzml_and_refuses_to_load_it` |
| 3 | 96 | `mode 3: streaming IMSDataConsumer` | ImzMLFile | `mode_3_a_streaming_load_reaches_the_consumer_once_per_spectrum` |
| 4 | 122 | `mode 4: loadSpectraIndex for on-disc access` | ImzMLFile + ImzMLHandlerHelper | `mode_4_the_index_only_load_gives_the_offsets_and_the_dataset_metadata` |
| 5 | 142 | `mode 5: OnDiscImzMLExperiment random access` | OnDiscImzMLExperiment | `mode_5_the_on_disc_reader_serves_one_pixel_at_a_time` |
| 6 | 167 | `mode 6: MSImagingExperiment in-memory pixel lookup` | ImzMLFile | `mode_6_the_imaging_load_maps_every_spectrum_to_a_pixel` |
| 7 | 190 | `mode 7: buildImagingGeometry from MSExperiment` | ImzMLFile | `mode_7_the_meta_value_geometry_builder_agrees_with_document_order` |
| 8 | 206 | `cross-mode consistency: full load vs on-disc vs RAM lookup` | ImzMLFile + OnDiscImzMLExperiment | `the_three_access_paths_return_the_same_pixel` |
| 9 | 237 | `processed imzML encoding` | ImzMLFile + OnDiscImzMLExperiment | `a_processed_dataset_works_on_the_in_memory_and_the_on_disc_path` |

### Beyond the suites

Fifteen further tests in `tests/imzml_file.rs` cover what no upstream section
reaches: the options and log-type accessors, a missing file on every load path,
the refusal of an unparseable filter, the `fill_data` and `metadata_only` loads,
the load peak ceiling, a mismatched `.ibd` UUID, a dangling geometry pixel, the
whole `ImagingExperiment` accessor surface, the unsorted-extraction refusal, all
five unplaceable-spectrum reports with the two builders' ordering difference, a
wrongly typed meta value, the two builders' agreement on both fixtures, and the
`.ibd` path rule.

Thirteen further tests cover the suite's own float oracle; see below.

---

## The float oracle — `ClassTest::isRealSimilar`

`TEST_REAL_SIMILAR` is what 38 assertions in `tests/imzml_file.rs` rest on, so
the helper that stands in for it is load-bearing and is tested directly.

The first version of that helper claimed to port
`src/testframework/source/CONCEPT/ClassTest.cpp:364-489` exactly and did not:
it omitted the opposite-sign branch at `:439-451`. Because the reciprocal of a
negative quotient is still negative, `ratio <= 1 + 1e-5` was then trivially
true, so the helper accepted **any sign error whose magnitudes matched** —
`close(-1.0, 1.0)`, `close(-100.0, 100.0)` and `close(-1e9, 3.0)` all passed
where C++ returns `false`. Every one of the 38 assertions was weaker than it
read.

`is_real_similar` in `tests/imzml_file.rs` is now the whole decision tree,
branch for branch, and `close` asserts it. **Re-running the file with the
faithful oracle changed no result: all 57 pre-existing tests still pass.** The
weak oracle was hiding no sign error — which is a real finding, not an excuse
for it, because nothing but re-running it could have established that.

Porting it faithfully surfaced two defects in upstream's own oracle. Both were
confirmed by compiling and executing upstream's control flow, not by reading it:

1. **Any two infinities are "similar".** The quotient is NaN and so is the
   absolute difference, so neither the sign test nor the ratio test fires and
   control reaches "ratio of numbers is small" — including `+inf` against
   `-inf`.
2. **`isRealSimilar` is not symmetric.** It decides "opposite signs" from the
   sign of the quotient, and that quotient underflows to `-0.0` once the
   magnitudes are far enough apart. `-0.0 < 0.` is false, so the sign branch is
   skipped; `-0.0 < 1.` then takes the reciprocal to `-inf`, and `-inf > 1 +
   1e-5` is false. So `isRealSimilar(1e-300, -1e300)` is `true` while
   `isRealSimilar(-1e300, 1e-300)` is `false`, and the same-sign pair
   `isRealSimilar(1e-300, 1e300)` is correctly `false` — which is what makes
   this a defect and not a deliberate tolerance.

Both are reproduced by `is_real_similar`, whose contract is to *be* the C++
oracle, and both are pinned by tests so they cannot be mistaken for port bugs.
`close`, the assertion the call sites use, then refuses both: it requires each
side to be finite, and decides "opposite signs" from the operands rather than
from the quotient. Those two guards can only reject pairs upstream would have
accepted, so **no assertion in the file is weaker than its C++ original**, and
neither guard can fire on a correct run, since every call site compares one
finite measured quantity against another.

The file's other hand-rolled comparison, `identical_holds`, was audited for the
same defect and does not have it: it is a plain absolute-plus-relative band
around its second argument with no quotient, so a sign flip lands a full
`2 * |right|` away and fails. A test pins that. No assertion in either file
compares computed floats with `==` where the oracle is the right tool: the one
exact float comparison, `reloaded.spectra[0].rt == -1.0` in section 29, asserts
that RT is still exactly its untouched default and so is deliberately exact,
and the float comparisons in `tests/imzml_all_modes.rs` use `to_bits()`, which
is stronger than either.

---

## Two defects in the upstream test fixtures

Both were found by running the ported `isValid` and the crate's mzML reader
against the fixtures. Neither contradicts the C++ suite, which validates only a
file it has just written, and neither affects any assertion here — the imzML
reader in this crate reads both fixtures correctly. Candidates for
`OpenMS_CPP_ISSUES.md`, which this package does not own.

1. **`ImzMLFile_1_Example_Continuous.imzML` is not schema-valid mzML 1.1.0.**
   Eleven errors from libxml2 against `mzML_1_10.xsd`: ten `cvParam` elements
   omit the required `cvRef` attribute, and `<scanSettingsList>` (line 35)
   precedes `<softwareList>` (line 57) while the schema's sequence is
   `softwareList`, `scanSettingsList`, `instrumentConfigurationList`,
   `dataProcessingList`. The same misordering is what makes this crate's mzML
   reader report `unresolved softwareRef`. `ImzMLFile::isValid` would return
   `false` for this fixture in C++ as well; no upstream section asks it to.
2. **`ImzMLFile_2_Example_Processed.imzML` declares `version="1.1"`, not
   `1.1.0`, and is `ISO-8859-1` with three Latin-1 bytes** (`ü` in a contact
   affiliation, `ß` twice in a contact address). The mzML 1.1.0 schema requires
   the three-component version; `MzMLFile` in C++ warns on an unexpected version
   rather than refusing, which is why the fixture has survived.

---

## Deferred

- **The mzML metadata layer**, above. Closing it needs three changes in
  `src/format/mzml.rs`, which this package does not own: tolerate a forward
  `softwareRef`, accept `version="1.1"`, and either transcode `ISO-8859-1` or
  ignore `IMS:1000101` on a binary array. Until then `ImzMLFile::load` produces
  no instrument, software or data-processing metadata and no per-spectrum RT or
  MS level.
- **`IMAGING/MSImagingExperiment.h` and its four siblings** are separate headers
  with their own class tests (`MSImagingExperiment_test.cpp`, `IonImage_test.cpp`,
  `MSImagingGeometry_test.cpp`, `MSImagingRegion_test.cpp`). `ImagingExperiment`
  reproduces the surface `ImzMLFile`'s two `MSImagingExperiment`-typed members
  need and nothing more; those suites were not ported section by section and no
  claim is made on the headers. A later IMAGING package owns them and should
  fold `ImagingExperiment`, `ImagingGeometry`, `ImagingRegion` and `IonImage`
  into its own module.
- **`FileHandler::loadExperiment` still refuses imzML**, as in C++, which is
  correct: the generic loader cannot know about the `.ibd`. `FileType::ImzMl`
  is recognised by extension and by content already. Routing `FileHandler`
  through `ImzMLFile` would change `src/format/file_handler.rs`, outside this
  package.
- **`docs/doc-coverage.json` was not rewritten.** `src/format/imzml_file.rs`
  measures 100.0% (48/48); recording the new floor with
  `check_doc_coverage.py --write` is the integrator's step.
- **CI wiring** is the integrator's step: append `--test imzml_file --test
  imzml_all_modes` to a `--features mzml` line of the `minimum-rust` job in
  `.github/workflows/rust.yml`. Both tests pass under `--no-default-features
  --features mzml` and under `cargo +1.85.0`; `is_valid` and its one test are
  gated on the existing `mzml-schema` feature, so they are simply absent there.
- **`OpenMS_CPP_ISSUES.md`** was not edited; the two fixture defects above are
  reported to the integrating agent, which owns that log.

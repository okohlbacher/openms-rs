# imzML writer support

Port of `src/openms/include/OpenMS/FORMAT/HANDLERS/ImzMLWriter.h` (87 lines)
with `src/openms/source/FORMAT/HANDLERS/ImzMLWriter.cpp` (1,523 lines), at core
SDK revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

- Rust: `src/format/imzml_writer.rs`
- Tests: `tests/imzml_writer.rs` (47 cases)
- Fixtures: the stage-1 `tests/data/ImzMLFile_1_Example_Continuous.{imzML,ibd}`
  and `tests/data/ImzMLFile_2_Example_Processed.{imzML,ibd}`, reused unchanged.
  This package adds none.
- Provenance: `tests/data/imzml_writer_provenance.json`
- Depends on: `docs/IMZML_HANDLER_SUPPORT.md` (stage 1), whose `.ibd` array
  writers, UUID codec, `.ibd` path rule and `ImzMLMeta` this package reuses and
  whose reader every round-trip test reads the output back with.

## What the writer produces

An imzML dataset is two files, and the writer emits both.

The `.ibd` opens with a 16-byte UUID header and then carries every binary array
back to back: uncompressed, little-endian, no padding, no per-array header. In
**continuous** mode one shared m/z array sits immediately after the header and
every spectrum's m/z offset points at it; each spectrum then contributes its
intensity array and its auxiliary arrays. In **processed** mode each spectrum
contributes its own m/z array first.

The `.imzML` is mzML 1.1.0 XML. Its `binaryDataArray` elements are *external*:
`<binary/>` is empty, and three IMS params say where the payload is —
`IMS:1000102` byte offset, `IMS:1000103` element count, `IMS:1000104` stored
byte length. Spectra carry 1-based pixel coordinates (`IMS:1000050` x,
`IMS:1000051` y, `IMS:1000052` z), and the file declares the `.ibd` identifier
(`IMS:1000080`) plus both its checksums (`IMS:1000091` SHA-1, `IMS:1000090`
MD5).

## Why this is not the mzML writer

The source writes the XML by hand, with raw `std::ostream` inserters in a
file-static `writeImzMLXml_`, rather than delegating to `MzMLHandler`. It has
to: an mzML serialiser writes base64 inside `<binary>`, and every array here is
external with an empty `<binary/>` and a byte offset into a second file. This
port keeps that division for the same reason.

What it does reuse:

| Reused from | What |
|---|---|
| `crate::format::imzml_handler` | `write_float32_array`, `write_float64_array`, `write_mz_as_float32`, `write_mz_as_float64` — the whole `.ibd` payload encoder; `uuid_bytes` for decoding a declared identifier; `infer_ibd_path`; `IBD_UUID_BYTES`; `ImzMLMeta`, `ImagingMode`, `ImzMLDataType`, `ImzMLReadLimits` |
| `crate::format::controlled_vocabulary` | `ControlledVocabulary::psi_ms` and `first_child_with_name`, which is the budgeted equivalent of the source's `iterateAllChildren("MS:1000513", …)` |
| `crate::format::peak_options` | `PeakFileOptions`, unchanged |
| `crate::concept::progress_logger` | `ProgressLogger`, threaded through as the source's fourth parameter |
| `quick_xml::escape` | the same escape `crate::format::mzml` uses |
| `sha1` | `IMS:1000091`, and the identifier derivation |

Nothing here encodes base64, parses XML, or re-implements an array codec.

## API mapping — `ImzMLWriter.h`

The header declares one class with one public member. Everything else in the
translation unit is in an anonymous namespace and is therefore not public API;
the parts of it that carry documented behaviour are mapped below the line
because a reader of the Rust needs them, not because the header exports them.

### `class ImzMLWriter`

| C++ member | Rust | Notes |
|---|---|---|
| `static void store(const std::string& imzml_path, const MSExperiment& exp, const PeakFileOptions& options, ProgressLogger& logger)` | `store_with_options(imzml_path, exp, options, write, logger) -> Result<StoreReport>` | The four source parameters in the same order, plus `write` for the ceilings and the shared-m/z tolerance. Returns what the source logs instead of `void`. `store(imzml_path, exp, options)` is the convenience form with default ceilings and a silent logger. |
| *(implicit ctor / dtor / copy)* | none | The source class is a namespace for one static function; it holds no state and the Rust equivalent is a free function. |
| `OPENMS_DLLAPI` on the class | none | Export visibility is a build-system concern; `pub` carries it. |

Every documented behaviour of that one `@brief` block is accounted for:

| Header statement | Rust |
|---|---|
| "Supports continuous (shared m/z array) and processed (per-spectrum m/z) modes" | `storage_mode`, `ImagingMode`, `StoreReport::mode` |
| "Pixel coordinates are read from `imzml:x/y/z` MetaValues on each spectrum (required for export)" | private `pixel_coord`; `Error::MissingInformation` when absent |
| "dataset metadata is taken from experiment MetaValues when present" | `dataset_meta` |
| "external binary arrays with a 16-byte UUID header in the `.ibd`" | `IBD_UUID_BYTES`, `derive_uuid`, `StoreReport::meta.uuid` |
| "Binary precision … follows PeakFileOptions (`getMz32Bit`, `getIntensity32Bit`)" | `PeakFileOptions::mz_32_bit` / `intensity_32_bit`; `StoreReport::meta.mz_data_type` / `int_data_type` |
| "`FloatDataArray` values are exported as additional external binary arrays after m/z and intensity" | the `aux` plan; `StoreReport::aux_arrays_written` |
| "`FloatDataArray` only (not integer/string data arrays)" | `StoreReport::dropped_data_array_count` / `dropped_data_array_names` |
| "unnamed arrays are skipped with a warning" | `FloatArraySkipReason::Unnamed` |
| "arrays named after the peak CV terms (`MS:1000514`, `MS:1000515`) are skipped with a warning: such an array would be read back as the spectrum's peak metadata" | `FloatArraySkipReason::ReservedPeakArrayName` |
| "arrays must have the same length as the spectrum peaks (others are skipped)" | `FloatArraySkipReason::LengthMismatch`, plus `FloatArraySkipReason::Empty` for the source's silent skip of an empty array |
| "always 32-bit float, uncompressed (`MS:1000576` on the array and on the m/z and intensity `referenceableParamGroup` entries)" | `external_aux_array` and the two group writers |
| "PSI-MS accession resolved via the ontology (children of `MS:1000513`)… exactly one allowed unit is written… several allowed units get no unit attributes… Unknown names become `MS:1000786`" | `resolve_float_array_cv`, `ResolvedArrayCv` |
| "Viewers can rely on `MSSpectrum::containsIMData()` after load for IM arrays" | Not ported: `containsIMData` belongs to `KERNEL/MSSpectrum.h` and is not in this package. The writer's half of the contract — the `MS:1003006` accession and its single `MS:1002814` unit — is written and asserted in `tests/imzml_writer.rs`. |
| "PeakFileOptions spectrum/peak filters (MS level, RT, precursor m/z, m/z and intensity ranges, metadata-only, sort-by-m/z) are applied to a temporary copy before export" | `apply_store_options`, called by `store_with_options` on a clone |
| "Spectra sharing a pixel coordinate are written out as-is with a warning" | `StoreReport::duplicate_pixel_count` / `duplicate_pixels`, `DuplicatePixel` |
| `@throws Exception::MissingInformation` (no spectra, or no `imzml:x`/`imzml:y`) | `Error::MissingInformation` |
| `@throws Exception::InvalidValue` (invalid pixel coordinates) | `Error::InvalidValue` |
| `@throws Exception::InvalidParameter` (incompatible continuous export) | `Error::InvalidValue` — this crate has no `InvalidParameter` variant and may not add one. Stated at `storage_mode`. |
| `@throws Exception::UnableToCreateFile` | `Error::Io` |
| `@throws Exception::ParseError` (binary serialisation failure) | `Error::Io` for a failed write, `Error::InvalidValue` for a ceiling. The source's `ParseError` here always wraps an I/O failure or the array-element ceiling. |
| `@ingroup FileIO` | dropped; module structure carries it |
| `@param[in] imzml_path / exp / options / logger` | the four arguments, documented at `store_with_options` |

### Anonymous-namespace helpers with public Rust counterparts

These are not header API. They are listed because a documented source behaviour
lives in each one and a reader needs to know where it went.

| C++ file-static | Rust | Notes |
|---|---|---|
| `resolveFloatArrayCv_` + `struct ResolvedArrayCv` + `ArrayCvCache` | `resolve_float_array_cv`, `ResolvedArrayCv` | The cache is private to one `store`, as the source's is |
| `extractMeta_` | `dataset_meta` | |
| `spectraShareMz_` | `spectra_share_mz` | the source's default `tolerance = 1e-5` is `SOURCE_SHARED_MZ_TOLERANCE` |
| `isContinuousMode_` | `storage_mode` | |
| `applyStoreOptions_` | `apply_store_options` | |
| `md5Hex_`, `md5ProcessBlock_`, `md5DigestToHex_` | `md5_hex` and a private streaming `Md5` | |
| `ensureUuidBytes_`, `uuidBytesToString_` | `derive_uuid`, `IBD_UUID_NAMESPACE`, private `uuid_to_string` | derivation replaces the source's RNG; see **Native differences** |
| `uuidStringToBytes_` | stage-1 `imzml_handler::uuid_bytes` | not re-implemented |
| `inferIbdPath_` | stage-1 `imzml_handler::infer_ibd_path` | not re-implemented |
| `elementByteSize_` | private `element_bytes` | |
| `validatePixelMetadataForStore_`, `readPixelCoord_`, `struct PixelCoord` | private `pixel_coord` + the duplicate pass in `Plan::build` | merged: the source validates as `Int` and re-reads as `uint32_t` |
| `warnOnDroppedDataArrays_` | the same pass, into `StoreReport` | |
| `updateGridFromPixels_` | the grid block at the end of `Plan::build` | |
| `instrumentModelForExport_` | private `instrument_model` | |
| `appendAndWriteFloatDataArrays_`, `struct AuxArrayWritePlan` | private `Planner::aux`, `AuxPlan` | |
| `struct SpectrumWritePlan` | private `SpectrumPlan` | `share_mz` is `Plan::shared_mz` instead of a per-spectrum flag |
| `writeMzArray_`, `writeIntArray_` | private `write_mz` and the intensity branch of `Plan::write_ibd` | |
| `writeIbdUuidHeader_` | private `write_ibd_uuid_header` | |
| `writeCvParam_` | private `cv_param` | |
| `writeImsGeometryCvParams_` | private `geometry_terms` | |
| `writeExternalBinaryArray_`, `writeExternalAuxBinaryArray_` | private `external_array`, `external_aux_array`, `external_extents` | |
| `writeImzMLXml_` | private `Plan::write_xml`, `Plan::write_spectrum` | |
| `sha1Hex_`, `sha1LeftRotate_`, `sha1DigestToHex_` | the `sha1` crate, in private `ibd_checksums` | the crate already depends on it for indexed mzML |
| `struct UniqueFile_` | `std::fs::File` | RAII is the language's |
| `floatDataArrayAsFloat_`, `intensitiesAsFloat_`, `mzAsDouble_` | inline `iter().map().collect()` | |
| `IBD_UUID_HEADER_BYTES` | stage-1 `imzml_handler::IBD_UUID_BYTES` | |

### Native API with no source counterpart

| Rust | Why |
|---|---|
| `StoreReport`, `DuplicatePixel`, `SkippedFloatArray`, `FloatArraySkipReason` | the source's warnings, returned instead of logged |
| `ImzMLWriteLimits` | resource ceilings; the source has none on the write side |
| `ImzMLWriteOptions`, `ImzMLWriteOptions::source`, `SOURCE_SHARED_MZ_TOLERANCE` | the lossy/lossless choice, in the established `dta::WriteOptions::source()` shape |
| `derive_uuid`, `IBD_UUID_NAMESPACE` | a deterministic identifier where the source uses an RNG |
| `store` (three-argument) | convenience over `store_with_options` |

## Preserved source conventions

- **The `.imzML` is written last.** A dataset that has an `.imzML` has a
  complete `.ibd`, because the offsets in the XML are only meaningful once the
  payload they name exists.
- **Mode selection.** An explicit `imzml:imaging_mode` of `"continuous"` is
  honoured only if the spectra really share an m/z axis and raises an error
  otherwise; `"processed"` is honoured unconditionally; anything else, including
  an unrecognised string, falls through to auto-detection. A single non-empty
  spectrum trivially shares its own axis, so a one-pixel dataset is written
  continuous unless told otherwise.
- **The grid is only ever raised.** `IMS:1000042` / `IMS:1000043` start from the
  experiment's declared `imzml:max_count_x` / `_y` and rise to the largest pixel
  actually written, so a declared grid larger than the data survives.
  `max_count_z` has no CV term and is the largest observed `z`, floored at 1.
- **The physical extents are recomputed**, not trusted: `max_dim_x = pixel_size_x
  × max_count_x` whenever the pixel size is positive, overwriting any
  `imzml:max_dim_x` the experiment carried. The upstream suite asserts the
  resulting 300 µm for the 3×3, 100 µm continuous fixture.
- **Duplicated pixels are written, not rejected.** The readers accept them, so a
  dataset that loads must be storable again. Readers map only the first spectrum
  per pixel into the imaging geometry.
- **Auxiliary arrays are always 32-bit float and uncompressed**, whatever the
  peak precision, and `MS:1000576` is written on them and on both
  `referenceableParamGroup`s. External arrays must be uncompressed for the
  readers to decode them at all.
- **A unit is written only when the ontology term allows exactly one.**
  `MS:1000821` "pressure array" gets `UO:0000110` pascal; `MS:1003007` "raw ion
  mobility array" allows milliseconds and seconds and gets none.
- **The array-identity name comes from the CV term, not the caller.** A resolved
  term writes its own `name`; only `MS:1000786` carries the caller's string, as
  the param's `value`.
- **`MS:1000016` "scan start time" is written for every spectrum**, including
  the OpenMS unset default of `-1`, and in seconds.
- **A missing native identifier becomes `spectrum=<index+1>`**, 1-based.
- **The instrument model falls back** to the instrument name and then to the
  literal `"OpenMS export"`, because `MS:1000031` must carry a value.
- **Skip order inside a spectrum's float arrays**: empty, then unnamed, then
  length mismatch, then reserved peak-array name. The order is observable
  because only the first matching reason is reported.
- **`DRange<1>::encloses` semantics** for every `PeakFileOptions` range: closed
  at the minimum, open at the maximum.
- **The XML layout is byte-for-byte the source's** in element order, attribute
  order and tab indentation, including the `encodedLength="0"` attribute on
  every external `binaryDataArray` and the empty `<binary/>` child.

## Native differences

Each is a deliberate divergence, stated at the item in the rustdoc as well.

1. **The identifier is derived, not random.** `ensureUuidBytes_` draws 16 bytes
   from `std::random_device` and stamps them as an RFC 4122 version 4 UUID.
   This crate has no random-number dependency and adding one is a `Cargo.toml`
   change outside this package, so `derive_uuid` computes a version **5**
   identifier instead: `SHA-1(IBD_UUID_NAMESPACE || payload)` over the `.ibd`
   bytes after the header, stamped. Consequences: a dataset written twice from
   the same arrays carries the same identifier, which makes the writer
   reproducible and testable; and two datasets with byte-identical payloads
   share an identifier. The identifier's job in imzML is to bind one `.imzML` to
   one `.ibd`, which it still does exactly. A caller that wants a specific
   identifier sets `imzml:uuid`, which both implementations honour verbatim.
2. **Continuous mode is not lossy by default.** `spectraShareMz_` compares m/z
   with an absolute tolerance of `1e-5` and then stores one axis, silently
   replacing every other spectrum's m/z within that window.
   `ImzMLWriteOptions::default()` uses `0.0`, so a dataset whose axes are not
   bit-identical is written processed instead, and an explicit `"continuous"`
   over such a dataset is refused. `ImzMLWriteOptions::source()` selects the
   source's tolerance and its loss. This follows `dta::WriteOptions::source()`.
3. **Nothing is written until everything is checked.** The source opens the
   `.ibd`, writes the UUID header and streams arrays, discovering a bad pixel
   coordinate or an incompatible mode part-way through and leaving a truncated
   `.ibd` with no `.imzML`. Here the whole plan — coordinates, offsets, lengths,
   resolved ontology terms, every ceiling — is built before either file is
   created, so a rejected experiment leaves nothing behind. `tests/imzml_writer.rs`
   asserts the absence of both files on every refusal path.
4. **Every offset is checked arithmetic against a ceiling.** The source
   accumulates `.ibd` offsets in a plain `uint64_t` and compares them against
   nothing. `ImzMLWriteLimits` bounds the spectrum count, the per-spectrum and
   total peak counts, the auxiliary-array counts, the `.ibd` size, the
   experiment-derived XML text, one string's length, the vocabulary lookups and
   the checksummed bytes; `Planner::reserve` is the single place an offset is
   produced, so the overflow check cannot be bypassed.
5. **Warnings are returned.** Every `OPENMS_LOG_WARN` becomes a field of
   `StoreReport`. The listings are capped at
   `ImzMLWriteLimits::max_reported_items`, defaulting to the source's own
   `max_listed` of 20, while the counts are exact.
6. **An empty `FloatDataArray` is reported.** The source skips it before any
   other check and without a warning, so a named but valueless array vanishes
   silently. `FloatArraySkipReason::Empty` records it.
7. **Wrong-typed and out-of-range metadata are errors.** `static_cast<UInt>` on
   a negative or above-`2^32` pixel count wraps silently in the source;
   `dataset_meta` and `pixel_coord` reject both. A meta value of the wrong type
   throws `Exception::ConversionError` in the source and is `Error::InvalidValue`
   here.
8. **Annotation arrays are kept aligned.** The source's sort and peak filters
   reorder and subset `spectrum` without touching its data arrays;
   `MSSpectrum::sort_by_position` and `MSSpectrum::select` move them together and
   refuse an array whose length is neither zero nor the peak count, leaving the
   spectrum unchanged on refusal.
9. **Strings are validated as XML 1.0.** `XMLHandler::writeXMLEscape` escapes
   the five entity characters and nothing else, so a native identifier holding,
   say, `U+0001` produces a document no parser can read back. Such a string is
   `Error::InvalidValue` here.
10. **Floating-point values round-trip.** The source formats doubles through
    `std::ostringstream` at `writtenDigits<double>()` = 15 significant digits,
    which is not always exact. Rust's `{}` writes the shortest representation
    that parses back to the same `f64`, so the port loses nothing the source
    would.
11. **Two file passes instead of three.** The source re-reads the `.ibd` once for
    MD5 and once for SHA-1. Here both digests are computed in one re-read, and
    the identifier's seed digest is accumulated in the pass that writes the
    payload.
12. **Progress nesting is balanced on the error paths.** `ProgressLogger`'s
    nesting depth is decremented whatever `store` returns. The source's
    `startProgress` / `endProgress` pair leaks its recursion counter when `store`
    throws between them — see **Source findings**.
13. **The port is serial.** The source's writer is serial too, so there is no
    OpenMP gap here; the parallel loop in this family is on the *read* side and
    is recorded in `docs/IMZML_HANDLER_SUPPORT.md`.
14. **`exp.updateRanges()` has no counterpart.** This crate computes ranges on
    demand in `MSExperiment::ranges`, so there is no cache to refresh. Nothing in
    the source's `store` reads the refreshed ranges either.
15. **MD5 is implemented here.** The source implements RFC 1321 inline for the
    same reason: OpenMS has no MD5 elsewhere and imzML 1.1.0 mandates
    `IMS:1000090`. `md5_hex` is public so the implementation can be checked
    against the published RFC 1321 appendix A.5 vectors, which
    `tests/imzml_writer.rs` does. MD5 is used only as a file-integrity
    fingerprint and is documented as unsuitable for authentication.

## Source findings

**`metadata_only` cannot store a continuous dataset.** `applyStoreOptions_`
clears every peak for `PeakFileOptions::getMetadataOnly()`, and
`isContinuousMode_` then runs `spectraShareMz_` over the emptied spectra, which
reports no shared axis because it cannot find a non-empty reference. An
experiment declaring `imzml:imaging_mode = "continuous"` — which is what loading
any continuous imzML produces — therefore fails its own metadata-only store with
`Exception::InvalidParameter`. The port reproduces this rather than quietly
downgrading the declared mode, and
`metadata_only_over_a_declared_continuous_dataset_is_refused_as_in_the_source`
pins it. A C++ fix would skip the shared-axis check when every spectrum is
empty, or honour the declared mode for a metadata-only store. Anchor:
`ImzMLWriter.cpp:485` (`clear(false)`) with `ImzMLWriter.cpp:738`
(`isContinuousMode_`).

**`ProgressLogger` nesting leaks on a throw.** `store` calls
`logger.startProgress(...)` at `ImzMLWriter.cpp:1437` and `logger.endProgress()`
only on the success path at `ImzMLWriter.cpp:1519`. Every throw in between —
`UnableToCreateFile` on the `.ibd`, `ParseError` from the array writers or the
flush, `UnableToCreateFile` on the `.imzML` — leaves `ProgressLogger`'s recursion
depth incremented, so subsequent progress output in the same process is indented
wrongly and, at the depth cap, suppressed. A C++ fix would use a scope guard.
The port decrements unconditionally.

**A partially written `.ibd` outlives a failed store.** Every throw after
`fopen` at `ImzMLWriter.cpp:1440` leaves a truncated `.ibd` next to no `.imzML`,
which a later run of the same pipeline can mistake for a stale companion. The
port's preflight removes every non-I/O cause of this.

Neither finding was reproduced against running C++ — see **Evidence** — so both
are labelled candidates in `tests/data/imzml_writer_provenance.json` and neither
is claimed as an applied upstream fix.

## Checked boundaries and evidence

### Ceilings

| Field | Default | Guards |
|---|---|---|
| `max_spectra` | 5,000,000 | the spectrum count and therefore the plan vector |
| `max_peaks_per_spectrum` | 100,000,000 | one array; also becomes `ImzMLReadLimits::max_array_elements`, the source's `MAX_IBD_ARRAY_ELEMENTS` |
| `max_total_peaks` | 4,000,000,000 | the dataset |
| `max_aux_arrays_per_spectrum` | 256 | one spectrum's auxiliary plan |
| `max_total_aux_arrays` | 10,000,000 | the dataset's |
| `max_ibd_bytes` | 1 TiB | every reserved offset, checked with `checked_add` |
| `max_text_bytes` | 1 GiB | cumulative experiment-derived XML text |
| `max_name_bytes` | 64 KiB | one identifier, array name or meta string |
| `max_cv_lookups` | 100,000 | distinct array names resolved against the ontology |
| `max_checksum_bytes` | 16 GiB | the `.ibd` re-read for its digests, checked in the preflight |
| `max_reported_items` | 20 | each `StoreReport` listing; the source's own `max_listed` |

`every_ceiling_refuses_before_a_file_is_created` drives all eleven and asserts
that neither output file exists afterwards.

### Evidence tier

**Tier 3 (source review) for the behavioural literals, with a Rust-only
round-trip closing the loop.** No C++ was built or executed and no C++ output
was retained, so this is not a tier 1 differential.

`ImzMLWriter.h` has no class test upstream. `ImzMLFile_test.cpp` exercises the
writer through `ImzMLFile::store`, and these literals are transcribed from its
store sections: the round trip of both fixtures with the mode preserved and a
non-empty `imzml:ibd_md5`; `MS:1003006` with `MS:1002814` and `MS:1000786` with
the free-text name, at `binaryDataArrayList count="4"`; `MS:1000821` "pressure
array" with `unitAccession="UO:0000110"`, `unitName="pascal"`,
`unitCvRef="UO"`; `MS:1003007` written with no `unitAccession`;
`binaryDataArrayList count="2"` and no `MS:1000786` for an unnamed array;
`count="2"` and an empty `aux` for arrays named after the peak arrays;
`MS:1000576` inside both `referenceableParamGroup`s; `Exception::MissingInformation`
for a spectrum with no coordinates; two spectra written and one geometry pixel
for a duplicated coordinate; `Exception::InvalidParameter` for an incompatible
continuous mode; `FLOAT64` for both arrays after `setMz32Bit(false)` /
`setIntensity32Bit(false)`; a shrunken first spectrum inside an m/z range; and
the geometry, instrument-model and 300 µm extent round trip.

Two strands do not depend on that suite:

- **The round trip through the stage-1 reader.** Every store test reads its own
  output back with `crate::format::imzml_handler` and compares decoded m/z,
  intensities, auxiliary arrays and pixel coordinates against what went in, for
  both upstream fixtures and for every synthetic case. That is what checks the
  `.ibd` offsets, which no transcribed literal can. The continuous round trip
  additionally asserts the exact `.ibd` size from first principles
  (`16 + peaks × 8 + 9 × peaks × 4`) and that the file on disk is that long, that
  every spectrum's `IMS:1000102` equals 16, and that the recomputed
  `IMS:1000091` SHA-1 verifies.
- **RFC 1321 appendix A.5.** The seven published MD5 vectors check the ported
  digest against a specification rather than against OpenMS.

The ceilings, the atomicity of a rejected store, the identifier derivation, the
exact-versus-source tolerance, the XML 1.0 validation and the error-variant
choices are independently derived (tier 4): the source has no analogue for any
of them.

### Comparison policy

Decoded m/z and intensities are compared **exactly**, not within a tolerance.
With the default `PeakFileOptions` the m/z path is `float32` in the fixture →
`f64` in memory → `float64` on disk → `f64`, and the intensity path is `float32`
throughout; both are exact widenings, so any difference would be a defect rather
than accumulation. Where a written value is re-read as text — pixel sizes,
extents, retention time — the comparison uses `1e-9`, which is looser than the
round-trip formatting guarantees and is there only to keep the assertion
readable.

## Deferred

- **Ledger, coverage and CI wiring are not written**, because
  `docs/core-sdk-reviewed-apis.json`, `docs/core-sdk-coverage.json`,
  `docs/doc-coverage.json`, `docs/VALIDATION.md`, `README.md`, `CHANGELOG.md`,
  `OpenMS_CPP_ISSUES.md` and `.github/workflows/rust.yml` are outside this
  package. The integrating agent should add `--test imzml_writer` to a
  `--features mzml` line of the `minimum-rust` job — the test passes under
  `--no-default-features --features mzml` and under `cargo +1.85.0` — run
  `python3 tools/core_sdk_coverage.py --write` and
  `python3 tools/check_doc_coverage.py --write`, and record both source findings
  above in `OpenMS_CPP_ISSUES.md`. `docs/core-sdk-coverage.json` was already
  stale at this branch point: stage 1's `ImzMLHandler.h` is still recorded as
  `unmapped` with no reference manifests.
- **No tier 1 differential.** An oracle driver under `../oracle/drivers/` linking
  the prebuilt `libOpenMS` could call `ImzMLWriter::store` on a fixed experiment
  and retain the `.imzML` and `.ibd` for byte comparison of the XML and exact
  comparison of the decoded arrays. Two fields would have to be excluded: the
  `<software version=...>` and the `IMS:1000080` identifier, which is random in
  the source and derived here. That is the single highest-value upgrade for this
  package.
- **`ImzMLFile::store`'s own layer is not ported.** The `MSImagingExperiment`
  overload that derives pixel coordinates from an `MSImagingGeometry` when the
  spectra carry no `imzml:x`/`imzml:y` lives in `FORMAT/ImzMLFile.h` and
  `IMAGING/MSImagingGeometry.h`, neither of which this package owns. The
  upstream `store(const MSImagingExperiment&)` section is therefore unreachable
  here, and an experiment whose coordinates exist only in a geometry is refused
  rather than exported.
- **`ProgressLogger` output is not compared.** The progress range and the
  balanced nesting are exercised; the rendered lines are
  `crate::concept::progress_logger`'s contract, not this package's.
- **The written `.imzML` is not read back by this crate's whole-file mzML
  reader.** The stage-1 reader validates the imzML half. Reading the written
  document as a complete `MSExperiment` needs the `FORMAT/ImzMLFile.h` stage,
  and `docs/IMZML_HANDLER_SUPPORT.md` lists the three mzML-reader strictnesses
  that stage must resolve first. The document this writer emits declares
  `<softwareList>` before `<dataProcessingList>`, so it does not reproduce the
  forward-`softwareRef` problem the upstream continuous fixture has.

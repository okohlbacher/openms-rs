# imzML handler support

Port of `src/openms/include/OpenMS/FORMAT/HANDLERS/ImzMLHandlerHelper.h`
(251 lines) with `src/openms/source/FORMAT/HANDLERS/ImzMLHandlerHelper.cpp`
(381 lines), and `src/openms/include/OpenMS/FORMAT/HANDLERS/ImzMLHandler.h`
(212 lines) with `src/openms/source/FORMAT/HANDLERS/ImzMLHandler.cpp`
(731 lines), at core SDK revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

- Rust: `src/format/imzml_handler.rs`
- Tests: `tests/imzml_handler.rs` (56 cases)
- Fixtures: `tests/data/ImzMLFile_1_Example_Continuous.{imzML,ibd}`,
  `tests/data/ImzMLFile_2_Example_Processed.{imzML,ibd}`
- Provenance: `tests/data/imzml_handler_provenance.json`

## What imzML is

An imzML dataset is two files. The `.imzML` is mzML 1.1.0 XML carrying metadata
and, per spectrum, IMS ontology params that give the byte offset
(`IMS:1000102`), element count (`IMS:1000103`) and stored byte length
(`IMS:1000104`) of that spectrum's arrays inside a companion `.ibd` binary file.
The `.ibd` opens with a 16-byte UUID that must equal the XML's `IMS:1000080`.
Spectra carry 1-based image coordinates (`IMS:1000050` x, `IMS:1000051` y,
`IMS:1000052` z). Two storage modes exist: **continuous** (`IMS:1000030`, one
shared m/z array stored once) and **processed** (`IMS:1000031`, a private m/z
array per spectrum). Both upstream fixtures are exercised here.

## Division of labour, and what this package is not

C++ `ImzMLHandler` derives from `MzMLHandler` and, by its own documentation,
intercepts only "the ~15 IMS:* CV terms that MzMLHandler does not know".
Everything else — instrument, data processing, retention time, MS level,
inline binary arrays — is the base class's work.

This port keeps that split. `crate::format::mzml` is the base class, and
`src/format/imzml_handler.rs` is the interception layer: it scans an `.imzML`
for the IMS vocabulary and nothing else, and reads `.ibd` ranges. It does not
parse mzML: no binary payload, no header metadata, no namespace or schema
validation, and no `PeakFileOptions` filtering. Joining the two halves into one
`MSExperiment` is `FORMAT/ImzMLFile.h`'s job and belongs to the stage that ports
it; see **Deferred** below for the two mzML-reader strictnesses it must deal
with first.

`IMAGING/MSImagingGeometry.h` and `IMAGING/IonImage.h` are separate headers in a
separate directory and are **not** part of this package. The geometry-shaped
type in the two headers this package owns is `ImzMLMeta`, which carries the
image dimensions, pixel size, scan pattern, scan direction and line-scan
direction; it is ported in full.

## API mapping — `ImzMLHandlerHelper.h`

Every public member of the header, in declaration order.

### `struct ImzMLMeta` → `ImzMLMeta`

| C++ member | Rust | Notes |
|---|---|---|
| `uint32_t max_count_x` | `ImzMLMeta::max_count_x` | `IMS:1000042`, raised to the largest observed `x` |
| `uint32_t max_count_y` | `ImzMLMeta::max_count_y` | `IMS:1000043`, raised to the largest observed `y` |
| `uint32_t max_count_z` | `ImzMLMeta::max_count_z` | no CV term; largest observed `z`. Source default is 1; `Default` here is 0 and `read_index` raises it to at least 1 for any dataset with a spectrum |
| `double pixel_size_x` | `ImzMLMeta::pixel_size_x` | `IMS:1000046`, µm |
| `double pixel_size_y` | `ImzMLMeta::pixel_size_y` | `IMS:1000047`, µm |
| `double max_dim_x` | `ImzMLMeta::max_dim_x` | `IMS:1000044`, µm |
| `double max_dim_y` | `ImzMLMeta::max_dim_y` | `IMS:1000045`, µm |
| `std::string imaging_mode` | `ImzMLMeta::imaging_mode: Option<ImagingMode>` | `ImagingMode::{Continuous, Processed}`; `None` is the source's empty string. `ImagingMode::as_str` reproduces the source spelling |
| `std::string ibd_file_path` | `ImzMLMeta::ibd_file_path: PathBuf` | the source's loader assigns it after construction; `ImzMLHandler::open*` fills in the file it opened |
| `std::string ibd_sha1` | `ImzMLMeta::ibd_sha1` | `IMS:1000091`, verifiable with `ImzMLHandler::verify_ibd_sha1` |
| `std::string ibd_md5` | `ImzMLMeta::ibd_md5` | `IMS:1000090`, parsed only — see **Checksums** |
| `std::string uuid` | `ImzMLMeta::uuid` | `IMS:1000080` as written, dashes and braces included |
| `std::string mz_data_type` | `ImzMLMeta::mz_data_type: ImzMLDataType` | `ImzMLDataType::name` gives the source's string, `Unknown` its empty state; **when** it is filled differs, see native differences |
| `std::string int_data_type` | `ImzMLMeta::int_data_type: ImzMLDataType` | as above |
| `std::string scan_pattern` | `ImzMLMeta::scan_pattern` | `"top down"` / `"bottom up"` |
| `std::string scan_direction` | `ImzMLMeta::scan_direction` | `"flyback"` / `"meander"` / `"horizontal"` / `"vertical"` |
| `std::string line_scan_direction` | `ImzMLMeta::line_scan_direction` | `"left-right"` / `"right-left"` |
| `std::string polarity` | `ImzMLMeta::polarity` | `"positive"` / `"negative"` |

### `struct ImzMLSpectrumIndex` → `ImzMLSpectrumIndex`

| C++ member | Rust | Notes |
|---|---|---|
| `enum class DataType : uint8_t { FLOAT32, FLOAT64, INT32, INT64, UNKNOWN }` | `ImzMLDataType::{Float32, Float64, Int32, Int64, Unknown}` | plus `width()`, `name()` (the source's `dtStr_` strings) and `from_accession()` (the source's four accession comparisons) |
| `struct AuxArray` | `ImzMLAuxArray` | all eight fields below |
| `AuxArray::name` | `ImzMLAuxArray::name` | CV term name for an `MS:1000513` child, `value` for `MS:1000786` |
| `AuxArray::accession` | `ImzMLAuxArray::accession` | |
| `AuxArray::unit_accession` | `ImzMLAuxArray::unit_accession` | |
| `AuxArray::offset` | `ImzMLAuxArray::offset` | `IMS:1000102` |
| `AuxArray::length` | `ImzMLAuxArray::length` | `IMS:1000103` |
| `AuxArray::encoded_bytes` | `ImzMLAuxArray::encoded_bytes` | `IMS:1000104` |
| `AuxArray::type` | `ImzMLAuxArray::data_type` | renamed: `type` is a Rust keyword |
| `AuxArray::compressed` | `ImzMLAuxArray::compressed` | child of `MS:1000572` other than `MS:1000576` |
| `int32_t index` | `ImzMLSpectrumIndex::index: u32` | 0-based document order; a count cannot be negative |
| `uint32_t x` | `ImzMLSpectrumIndex::x` | `IMS:1000050`, 1-based |
| `uint32_t y` | `ImzMLSpectrumIndex::y` | `IMS:1000051`, 1-based |
| `uint32_t z` | `ImzMLSpectrumIndex::z` | `IMS:1000052`, 1-based, default 1 |
| `uint64_t mz_offset` | `ImzMLSpectrumIndex::mz_offset` | |
| `uint64_t mz_length` | `ImzMLSpectrumIndex::mz_length` | |
| `DataType mz_type` | `ImzMLSpectrumIndex::mz_type` | |
| `bool mz_compressed` | `ImzMLSpectrumIndex::mz_compressed` | |
| `uint64_t int_offset` | `ImzMLSpectrumIndex::int_offset` | |
| `uint64_t int_length` | `ImzMLSpectrumIndex::int_length` | |
| `DataType int_type` | `ImzMLSpectrumIndex::int_type` | |
| `bool int_compressed` | `ImzMLSpectrumIndex::int_compressed` | |
| `std::vector<AuxArray> aux` | `ImzMLSpectrumIndex::aux` | named arrays only, as the source's builder |
| — | `native_id`, `mz_encoded_bytes`, `int_encoded_bytes`, `mz_external`, `int_external`, `unnamed_aux`, `inline_aux_names` | native additions, each justified under **Native differences** |

### `class ImzMLBinaryIO` → `ImzMLBinaryIO` and four free functions

| C++ member | Rust | Notes |
|---|---|---|
| `static void readMzArray(FILE*, offset, count, dt, std::vector<double>& out, ibd_path)` | `ImzMLBinaryIO::read_mz_array(offset, count, data_type) -> Result<Vec<f64>>` | the `FILE*` and `ibd_path` become the receiver; the out-parameter becomes the return |
| `static void readIntArray(FILE*, offset, count, dt, std::vector<float>& out, ibd_path)` | `ImzMLBinaryIO::read_intensity_array(offset, count, data_type) -> Result<Vec<f32>>` | |
| `static void readAuxArray(FILE*, offset, count, dt, std::vector<float>& out, ibd_path, array_name)` | `ImzMLBinaryIO::read_aux_array(offset, count, data_type, name) -> Result<Vec<f32>>` | an empty `name` becomes `"auxiliary array"` in messages, as the source does |
| `static void writeFloat32Array(FILE*, const float*, count, ibd_path)` | `write_float32_array(out: impl Write, data: &[f32], limits) -> Result<()>` | free function: writing targets a stream the writer owns, not this reader. Pointer + count become a slice, so the source's null-`FILE*` error cannot arise |
| `static void writeMzAsFloat32(FILE*, const std::vector<double>&, ibd_path)` | `write_mz_as_float32(out, mz: &[f64], limits)` | narrows with the same `as f32` the source's `static_cast` performs |
| `static void writeFloat64Array(FILE*, const double*, count, ibd_path)` | `write_float64_array(out, data: &[f64], limits)` | |
| `static void writeMzAsFloat64(FILE*, const std::vector<double>&, ibd_path)` | `write_mz_as_float64(out, mz: &[f64], limits)` | forwards to `write_float64_array`, as the source does |

File-static helpers of `ImzMLHandlerHelper.cpp` (not header members, ported
because they carry the behaviour): `MAX_IBD_ARRAY_ELEMENTS` →
`ImzMLReadLimits::max_array_elements`; `validateCount_` → the count half of
`ImzMLBinaryIO::preflight`; `seekIbd_` → `std::io::Seek`, whose 64-bit offsets
make the source's platform `off_t` ceiling check unnecessary; `throwReadError_`
→ `Error::Parse`; `swapU32_` / `swapU64_` / `decodeLittleEndian_` →
`from_le_bytes` / `to_le_bytes`, which need no host test.

## API mapping — `ImzMLHandler.h`

| C++ member | Rust | Notes |
|---|---|---|
| `class ImzMLInterceptConsumer;` (forward declaration, defined in the `.cpp`) | not ported | It exists to re-associate batched `MzMLHandler` deliveries with per-spectrum IMS state through a running counter. This port builds the index in one pass and returns it by value, so there is no delivery order to re-synchronise and no consumer to bridge. Its decode body is `ImzMLHandler::spectrum` |
| `ImzMLHandler(PeakMap& exp, const std::string& filename, const ProgressLogger& logger)` | `ImzMLHandler::open`, `open_with_ibd`, `open_with_limits` | The source binds a `PeakMap` that its base class fills; nothing here writes into a caller's container. `ProgressLogger` is not ported anywhere in this module — see **Deferred** |
| `~ImzMLHandler()` (closes the `.ibd` `FILE*`) | `Drop` of the owned `File` | |
| `ImzMLHandler(const ImzMLHandler&) = delete` | no `Clone` impl | |
| `ImzMLHandler& operator=(const ImzMLHandler&) = delete` | not applicable | |
| `void openIBD(const std::string& ibd_path)` | `ImzMLBinaryIO::open` / `open_with_limits`, called by `ImzMLHandler::open*` | The source can be called twice and closes the previous handle; a handler here owns exactly one `.ibd` for its lifetime, so a second dataset is a second handler |
| `const ImzMLMeta& getImzMLMeta() const noexcept` | `ImzMLHandler::meta` | |
| `ImzMLMeta& getImzMLMeta() noexcept` | not ported as a mutable accessor | Its one documented use is for the loader to set `ibd_file_path` before the parse; `open_with_limits` sets it from the file it opened, so no caller needs to reach in |
| `const std::vector<ImzMLSpectrumIndex>& getIndex() const noexcept` | `ImzMLHandler::index` (slice), `parsed` (index and metadata together), `entry(i)` (one entry, checked) | `entry` matches `OnDiscImzMLExperiment::getIndex(i)`, which raises `Exception::IndexOverflow` |
| `void connectDecodeConsumer(IMSDataConsumer*, bool append_spectra_to_map, bool decode_ibd = true)` | not ported | It wires the intercept consumer, flips `PeakFileOptions::setAlwaysAppendData` and records `decode_ibd`. The port has no consumer chain: index-only is `read_index`, decoding is `spectrum`, and both are explicit calls rather than a mode flag |
| `void onStartElement(const char16_t*, const XMLAttributes&) override` | `read_index_with_limits`'s `Event::Start` arm | SAX override becomes a pull-parser arm |
| `void onEndElement(const char16_t*) override` | `read_index_with_limits`'s `Event::End` arm | including the missing-coordinate check on `</spectrum>` |
| `friend class ImzMLInterceptConsumer` | not applicable | |
| inherited `MzMLHandler` public surface (`setOptions`, `getOptions`, `setMSDataConsumer`, `parse_`, …) | belongs to `MzMLHandler.h` | Not this package's headers. The Rust equivalent of that base is `crate::format::mzml`, already ported |

Private members, for completeness, because they carry behaviour that is now
public: `handleIMSCvParam_` → `Parser::handle` (the same three context-ordered
blocks); `applyRefGroup_` → `Parser::apply_ref_group`; `struct ArrayMeta` →
private `ArrayMeta`, whose `is_ext`, `encoded_bytes` and name fields are
promoted into the public index; `struct SpecIMS` → the per-spectrum fields of
`Parser` plus `inline_aux_names` on the public entry; `struct CvEntry` →
private `CvEntry`; `ref_groups_`, `cur_*`, `in_*`, `spec_ims_`, `index_`,
`decode_bridge_` → `Parser` state.

## Preserved source conventions

- **Three-block CV dispatch, in order.** Coordinates only inside a `<scan>`
  inside a `<spectrum>`, and only the three coordinate accessions return from
  that block, so an `MS:1000129` polarity term seen inside a scan still reaches
  the dataset block. Then the `binaryDataArray` block, which returns
  unconditionally, so a dataset-level term inside an array is ignored and
  cannot overwrite the array's name. Then dataset-level metadata.
  (`ImzMLHandler.cpp:637-720`.)
- **Compression is "any child of `MS:1000572` that is not `MS:1000576`"**, a
  vocabulary question rather than a list, so numpress and zlib+numpress are
  caught without being enumerated. `MS:1000576` clears the flag.
- **Array identity is `MS:1000786` or a child of `MS:1000513`**, and a
  `MS:1000513` child takes the CV term's own name, not the XML `name`
  attribute, so a decoded array is named the way `MzMLHandler` names it.
- **Referenceable parameter groups are captured and replayed at the reference**,
  in the context prevailing there; an undeclared reference is ignored, as the
  source's failed map lookup is.
- **A spectrum with no `IMS:1000050` or no `IMS:1000051` is a parse error**
  (`ImzMLHandler.cpp:561`), which is the upstream `load rejects spectrum missing
  pixel coordinates` section.
- **`z` defaults to 1** per spectrum, and `max_count_z` has no CV term at all.
- **The declared bounding box is a floor.** `IMS:1000042`/`43` are raised to the
  largest observed coordinate (`ImzMLHandler.cpp:193-195`).
- **Duplicate pixel coordinates are tolerated**, and the first spectrum in
  document order owns the pixel, as the source's geometry builder and the
  upstream `tolerates duplicate pixel coordinates by default` sections.
- **Compressed external arrays are rejected at decode time, never inflated**,
  with the source's own advice in the message ("need uncompressed MS:1000576.
  Re-export without compression."). A compressed auxiliary array's message
  names its `IMS:1000104` length, as the source's does.
- **Auxiliary arrays are warned about and skipped, not fatal**, for an unnamed
  array, a zero length, a length that is not the peak count, and a missing
  binary data type. The upstream `load skips one bad aux length and still
  returns all spectra` section depends on this.
- **An unnamed external auxiliary array is dropped from the index entry**, so
  `aux.len()` counts what the source counts.
- **A zero-length array never touches the file.** `readMzArray` returns before
  its seek, so a garbage offset on a zero-length array is not an error.
- **Numeric conversions.** m/z widen to `f64` and intensities narrow to `f32`
  with the same casts the source applies, including the precision a 64-bit
  integer m/z above 2^53 loses in both.
- **`.ibd` arrays are little-endian** (imzML 1.1.0). The source byte-swaps under
  `OPENMS_IS_BIG_ENDIAN`; `from_le_bytes` does it on every host with no
  conditional compilation, so there is no untested big-endian path.
- **`.ibd` path inference.** A case-insensitive `.imzML` suffix becomes `.ibd`;
  any other name gains `.ibd` (`ImzMLFile::inferIbdPath_`). An explicit `.ibd`
  override is used for both the index and the UUID check, as
  `OnDiscImzMLExperiment::open(imzml, ibd)` threads it through.
- **UUID string parsing** strips dashes and braces and requires exactly 32 hex
  digits (`ImzMLFile.cpp:45`), so reader and writer agree on the header bytes.
- **An empty `IMS:1000080` keeps the previous value** (`ImzMLHandler.cpp:707`).

## Native differences

Each one is also stated at the item in rustdoc.

1. **Bounds are preflighted, not discovered.** Source `readMzArray` checks the
   element count against `MAX_IBD_ARRAY_ELEMENTS`, then `resize`s the output,
   then lets a short `fread` report the truncation; `readFloatVector_` does the
   same. Nothing compares the range against the `.ibd`'s actual length, so a
   hostile `IMS:1000103` at that ceiling commits 800 MB for a float64 array, or
   1.2 GB for a float32 one, which `readMzArray` stages through a second
   `std::vector<float>` before widening. Here
   `ImzMLBinaryIO::preflight` rejects an unknown data type, a count above the
   ceiling, a count above `usize`, a stored size above `max_array_bytes`, an
   `offset + bytes` that overflows `u64`, and a range that leaves the measured
   file — all before a single byte is allocated, and the allocation itself uses
   `try_reserve_exact`.
2. **The index scan has ceilings of its own.** `ImzMLReadLimits` bounds spectra,
   auxiliary arrays per spectrum and in total, XML bytes, referenceable groups
   and their parameters, stored string bytes, vocabulary lookups and checksum
   bytes. The source bounds none of these.
3. **A cvParam no imzML rule can act on is skipped without its value being
   decoded.** The source reaches the same conclusion one step later — its
   dispatch ignores every unlisted accession, explicitly so inside a
   `binaryDataArray` "so they cannot overwrite the array name"
   (`ImzMLHandler.cpp:694`) — but Xerces has already transcoded the value by
   then. Not decoding it is why the upstream processed fixture, which declares
   `ISO-8859-1` and carries a non-UTF-8 byte in an `MS:1000590` contact
   affiliation, still indexes. A non-UTF-8 byte in a value this parser *does*
   read is `Error::Parse`, not a silent replacement: the source relies on Xerces
   to transcode from the declared encoding and this port implements no
   transcoder.
4. **Data types are recorded at parse time, not decode time.** The source sets
   `mz_data_type`/`int_data_type` inside the decode branch of
   `consumeSpectrum`, so `loadSpectraIndex` — its index-only path — leaves both
   empty. Taking the first declared type instead keeps a metadata-only parse
   informative, and agrees with the source on any file whose first spectrum
   declares a type, which includes both fixtures.
5. **Omissions are returned, not logged.** The source writes five distinct
   warnings to `OPENMS_LOG_WARN`: an unnamed auxiliary array, a zero-length one,
   a length mismatch, a missing data type, and an inline auxiliary array
   alongside external peaks. `DecodedSpectrum::skipped_aux` carries all five as
   `AuxSkipReason` values and `inline_peaks` reports peak arrays this handler
   could not supply, so a caller can notice what a log would have buried. This
   is also why `unnamed_aux` and `inline_aux_names` are public.
6. **The UUID and checksum verdicts are returned, not warned.**
   `verifyIbdUuid_` logs for a mismatch, a missing or unparsable XML UUID, an
   unopenable `.ibd` and a too-short `.ibd`, then loads the dataset anyway —
   deliberately, so a non-conformant file still opens. `uuid_status` returns
   `UuidStatus`, and `open*` rejects nothing, so the source's tolerance is kept
   and the policy is the caller's.
7. **`Option<ImagingMode>` replaces a two-valued string**, `ImzMLDataType`
   replaces the data-type strings, and both keep the source spelling behind
   `as_str` / `name`.
8. **Non-finite `IMS:1000044`–`47` values are rejected.** `std::stod` accepts
   `inf` and `nan`; a NaN pixel size cannot be used by any consumer and would
   propagate silently.
9. **IMS integers require digits only.** The source's `std::stoull` tolerates
   leading whitespace and a `+` sign and then rejects trailing characters;
   `parseImsUInt64_` separately rejects any value containing `-` because
   `stoull` would wrap it. This trims ASCII whitespace and then requires digits,
   which keeps every source rejection and additionally rejects `+`. No imzML
   writer emits one and no upstream fixture carries one.
10. **Errors are `Result`, and one `ParseError` becomes `Unsupported`.** The
    source raises `Exception::ParseError` for a compressed external array and
    for an unsupported data type; both are `Error::Unsupported` here, because
    the input is valid imzML asking for a feature this port does not implement,
    which is what that variant means. Every other `ParseError` maps to
    `Error::Parse`, an out-of-range index to `Error::InvalidValue` (matching
    `Exception::IndexOverflow`), a missing pixel to `Error::InvalidValue` on the
    coordinate lookup (matching `Exception::ElementNotFound`), and a failed
    `fopen` to `Error::Io` (matching `Exception::FileNotFound`).
11. **`Error::Parse` carries line 0.** quick-xml reports a byte position rather
    than a line, and the crate's other XML readers do the same.
12. **Exclusive access is enforced.** Every decode takes `&mut self`. The source
    shares one `FILE*` and only documents that this is not thread-safe.
13. **Coordinate lookup is a linear scan.** `getSpectrumAtCoord` is served by an
    `MSImagingGeometry` grid built during `open()`, O(1) but 2-D, so its own
    documentation says only the `z == 1` plane is addressable.
    `index_at_coord` scans the index, so every `z` in the file is addressable;
    the grid type belongs to a header this package does not own.
14. **Peaks are not sorted.** The source re-sorts external peaks when
    `PeakFileOptions::getSortSpectraByMZ()` is set, because its base class's
    sort ran before the external arrays overwrote the peaks. This handler holds
    no `PeakFileOptions`; a caller that needs the guarantee calls
    `MSSpectrum::sort_by_position`.
15. **`native_id` is recorded.** The source leaves native identifiers to its
    base class. A stand-alone index cannot be matched against mzML records
    without one, and stage 3 needs exactly that.
16. **No `default_array_length_` zeroing and no empty-array pruning.** Two of
    the source's most intricate passages exist only to stop `MzMLHandler` from
    logging a size mismatch for every spectrum from inside an OpenMP loop — a
    heap corruption on Windows — and to delete the empty data arrays that base
    class had already created (`ImzMLHandler.cpp:542-553`, `595-611`,
    `133-146`). Neither applies: this port never hands the empty inline
    `<binary/>` placeholders to the mzML reader, and it creates a data array
    only when it has values for it.
17. **Serial.** The source parallelises `populateSpectraWithData_` with
    `#pragma omp parallel for`. Nothing here starts a thread, so a large image
    decodes on one core.

## Checksums

The `.ibd` UUID header is read and compared on request
(`ImzMLHandler::uuid_status`), which is what the source does, one layer up, in
`ImzMLFile.cpp`.

`IMS:1000091` **SHA-1 is parsed and verifiable**:
`ImzMLHandler::verify_ibd_sha1` re-hashes the whole `.ibd` in 64 KiB chunks and
compares case-insensitively, bounded by `max_checksum_bytes`. No OpenMS read
path recomputes it; the upstream processed fixture's declared digest is
reproduced in `tests/imzml_handler.rs`. Verification is never automatic,
because hashing a multi-gigabyte `.ibd` is a full extra pass.

`IMS:1000090` **MD5 is parsed only, never verified.** The reason is the
dependency set: this crate has a SHA-1 implementation because indexed mzML
needs one, and no MD5 implementation. Adding a dependency is a `Cargo.toml`
change, which is outside this package's scope. The declared string is preserved
verbatim in `ImzMLMeta::ibd_md5` so a later stage can verify it without
re-reading the XML.

## Checked boundaries and evidence

| Boundary | Rust outcome |
|---|---|
| element count > `max_array_elements` (source's 100,000,000) | `Error::InvalidValue` before any allocation |
| element count > `usize` | `Error::InvalidValue` |
| count × element width > `max_array_bytes`, or overflowing | `Error::InvalidValue` |
| `offset + bytes` overflowing `u64` | `Error::InvalidValue` |
| `offset + bytes` past the measured `.ibd` length | `Error::Parse` naming the range and the file length |
| short read inside a checked range (file shrank) | `Error::Parse` |
| `ImzMLDataType::Unknown` on a decode | `Error::Unsupported` |
| compressed m/z, intensity or auxiliary array | `Error::Unsupported` with the source's message |
| decoded m/z and intensity lengths differ | `Error::Parse` naming the pixel, as the source's message |
| spectrum index out of range | `Error::InvalidValue` |
| auxiliary index out of range | `Error::InvalidValue` |
| no spectrum at a coordinate | `Error::InvalidValue` |
| spectra > `max_spectra`, aux > `max_aux_arrays` / `max_total_aux_arrays` | `Error::InvalidValue` |
| XML bytes > `max_xml_bytes`; groups > `max_param_groups`; group params > `max_group_params`; stored text > `max_text_bytes`; vocabulary lookups > `max_cv_lookups` | `Error::InvalidValue` |
| `.ibd` longer than `max_checksum_bytes` on a SHA-1 request | `Error::InvalidValue` |
| spectrum count > `u32` | `Error::InvalidValue` |
| malformed XML, empty/negative/non-numeric/out-of-range IMS value, non-finite dimension, missing pixel coordinate, non-UTF-8 in a read value | `Error::Parse` |
| missing `.imzML` or `.ibd` | `Error::Io` |

**Evidence tier 3, source review**, with two independent cross-checks. Neither
of these headers has a class test upstream; the imzML family is tested through
`ImzMLFile_test.cpp` and `ImzMLFile_all_modes_test.cpp`, whose literals are
transcribed here: nine spectra on a 3x3 grid, 100 µm pixels, 300 µm extents,
`continuous` / `processed`, the two UUIDs, the MD5
`4b5dd9fa84fafc955cfdd301f9ed55d7` and SHA-1
`7e8fdb93053915d3edb51b70aa0619ac209964df`, `float32` for both arrays,
`negative` polarity with `top down` / `horizontal` / `left-right` geometry, and
the first m/z of pixel (1,1) — 100.0 continuous, 100.083336 processed.

The two independent checks are stronger than a transcription. The decoded
intensities of processed pixel (1,1) sum to 121.85039039868471, which is exactly
the `MS:1000285` total ion current that same fixture's XML declares — a number
written by the instrument's exporter, not by a test. And the fixture's declared
`IMS:1000091` SHA-1 is recomputed from the `.ibd` bytes and matches, which
confirms that the file this port reads is the file the fixture describes.

No C++ was built or executed and no C++ output was retained, so this is not a
tier 1 differential. The synthetic documents, the resource ceilings and the
hostile offsets are independently derived (tier 4): no upstream fixture reaches
them.

## Deferred

- **`FORMAT/ImzMLFile.h`, `KERNEL/OnDiscImzMLExperiment.h` and
  `FORMAT/HANDLERS/ImzMLWriter.h`** are separate headers and separate stages.
  This package is their foundation: `ImzMLIndex`, `ImzMLMeta`,
  `ImzMLSpectrumIndex`, `ImzMLHandler::{open, open_with_ibd, open_with_limits,
  spectrum, spectrum_at_coord, mz_array, intensity_array, aux_array}` and the
  four `write_*` functions are what they build on.
- **`IMAGING/MSImagingGeometry.h`, `IMAGING/IonImage.h`,
  `IMAGING/MSImagingExperiment.h`, `IMAGING/MSImagingRegion.h`,
  `IMAGING/IonImageExtraction.h`** are a separate directory with their own class
  tests (`IonImage_test.cpp`, `MSImagingExperiment_test.cpp`) and are unported.
  `index_at_coord` is the minimum this package needed and is not a substitute
  for that grid.
- **Two mzML-reader strictnesses block a combined load**, and both are in
  `src/format/mzml.rs`, which this package does not own. `mzml::read` rejects
  `ImzMLFile_1_Example_Continuous.imzML` with `unresolved softwareRef` — the
  fixture puts `<softwareList>` *after* `<dataProcessingList>`, so the
  single-pass reference resolution sees a forward reference that C++ tolerates —
  and rejects `ImzMLFile_2_Example_Processed.imzML` twice over, first because it
  declares `version="1.1"` where the reader requires a `1.1.` prefix, and then
  because the reader refuses non-ASCII `ISO-8859-1` ("requires transcoding").
  The stage that ports `ImzMLFile` must resolve these three before an imzML
  dataset can become an `MSExperiment` with its mzML metadata attached.
- **`ProgressLogger`** is a constructor parameter of source `ImzMLHandler` and
  is not threaded through this module; the crate has the type but this package
  reports no progress.
- **`PeakFileOptions` filtering** (RT, MS level, m/z range, `SortSpectraByMZ`,
  `FillData`) is applied by `ImzMLFile` / `OnDiscImzMLExperiment` around the
  handler, not by it, and is left to those stages.
- **`docs/doc-coverage.json` was not rewritten.** `src/format/imzml_handler.rs`
  measures 100.0% (63/63) and `check_doc_coverage.py` reports "coverage improved
  or the surface changed"; recording the new floor with `--write` is the
  integrator's step, because that file is outside this package's scope.
- **CI wiring** is the integrator's step for the same reason: append
  `--test imzml_handler` to a `--features mzml` line of the `minimum-rust` job
  in `.github/workflows/rust.yml`. The test passes under
  `--no-default-features --features mzml` and under `cargo +1.85.0`.

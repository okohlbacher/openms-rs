# mzXML support in the Rust port

Port of `src/openms/include/OpenMS/FORMAT/MzXMLFile.h` (111 lines) with
`src/openms/source/FORMAT/MzXMLFile.cpp` (112 lines), and
`src/openms/include/OpenMS/FORMAT/HANDLERS/MzXMLHandler.h` (186 lines) with
`src/openms/source/FORMAT/HANDLERS/MzXMLHandler.cpp` (1296 lines), at core SDK
revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

- Rust: `src/format/mzxml.rs` (gated on the existing `mzml` feature, for
  quick-xml, base64 and flate2)
- Tests: `tests/mzxml.rs` (59 cases)
- Fixtures: `tests/data/MzXMLFile_1.mzXML`,
  `tests/data/MzXMLFile_1_compressed.mzXML`,
  `tests/data/MzXMLFile_2_minimal.mzXML`, `tests/data/MzXMLFile_3_64bit.mzXML`
  (all unmodified upstream) and `tests/data/MzXMLFile_5_nested.mzXML`
  (hand-written for this package)
- Provenance: `tests/data/mzxml_provenance.json`

## What mzXML is, and where it differs from mzML

mzXML 3.1 is mzML's predecessor. Five differences matter for the port:

| | mzXML | mzML |
|---|---|---|
| Peak arrays | one base64 payload, m/z and intensity interleaved as **pairs** | two independent `binaryDataArray` elements |
| Element width | `precision="32"` / `"64"` attribute | `MS:1000521` / `MS:1000523` CV term |
| Byte order | big-endian, `byteOrder="network"` | little-endian, fixed |
| Compression | `compressionType="none"`/`"zlib"` attribute | `MS:1000574` CV term |
| MS2 placement | `<scan>` **nested inside** its MS1 parent | sibling `<spectrum>` elements |

Metadata is attribute-driven rather than CV-driven: `<msIonisation value="ESI">`
instead of `MS:1000073`. Five small string tables in `MzXMLHandler::init_`
(`MzXMLHandler.cpp:1271-1294`) map those strings to OpenMS enums by **position**,
so the empty leading entry is each enum's `Unknown`. Those tables are
transcribed verbatim into `src/format/mzxml.rs` and a unit test checks that each
term's index equals the corresponding `crate::metadata` enum's index.

### Nesting carries no data

An MS2 `<scan>` inside an MS1 `<scan>` still becomes an ordinary spectrum in
document order; upstream's `nesting_level_` only decides when the pending batch
is flushed. This port reproduces that: spectra are appended in start order and a
batch is committed only when the scan stack is empty, so a whole nested tree is
committed together. Reading `MzXMLFile_1.mzXML` therefore yields `scan=10`,
`scan=11`, `scan=12`, `scan=13` in that order, with `scan=13` (the MS2 child of
`scan=12`) last — which is what `MzXMLFile_test.cpp:125-137` asserts.

The writer reproduces the nesting from MS levels alone: it opens one `<scan>` per
spectrum and closes `ms_level - next_ms_level + 1` of them when the next
spectrum's level is not higher (`MzXMLHandler.cpp:1082-1096`). All six MS-level
patterns of `MzXMLFile_test.cpp:500-568` are exercised and each one produces a
balanced document that reloads to five spectra with the same levels.

## API mapping — `FORMAT/MzXMLFile.h`

Every declared member, plus the surface it inherits.

| C++ member | Rust | Notes |
|---|---|---|
| `MzXMLFile()` | `MzXMLFile::new` / `Default` | registers `SCHEMA` and `SCHEMA_VERSION` |
| `~MzXMLFile() override` | derived `Drop` | nothing to release |
| `PeakFileOptions& getOptions()` | `MzXMLFile::options.peaks` (public field) | mutable access is field access |
| `const PeakFileOptions& getOptions() const` | `MzXMLFile::options.peaks` | same field, shared borrow |
| `void setOptions(const PeakFileOptions&)` | assign `MzXMLFile::options.peaks` | |
| `void load(const std::string&, MapType&)` | `MzXMLFile::load`, `load`, `load_with_options`, `load_into`, `load_into_with_options`, `read`, `read_with_options`, `read_with_report`, `read_metadata`, `read_scan_count` | value-returning; the `_into` forms replace a destination only after a complete parse |
| `void store(const std::string&, const MapType&) const` | `MzXMLFile::store`, `store`, `store_with_options`, `write`, `write_with_options` | atomic publication through the shared temp-file writer |
| `void transform(filename, consumer, skip_full_count)` | `MzXMLFile::transform`, `transform_with_options` | |
| `void transform(filename, consumer, map, skip_full_count)` | `MzXMLFile::transform_into`, `transform_into_with_options` | forces `always_append_data`, as the source does |
| `protected void transformFirstPass_(...)` | first pass inside `run_transform` | private; its effect is visible through `TransformReport::expected_spectra` |
| `private PeakFileOptions options_` | `MzXMLFile::options` | `ReadOptions` wraps it together with `ReadLimits` |
| `private typedef PeakMap MapType` | `crate::kernel::MSExperiment` | |
| inherited `XMLFile::getVersion()` | `MzXMLFile::version`, `SCHEMA_VERSION` | |
| inherited `XMLFile::isValid(filename, ostream)` | **not ported** | needs `mzXML_idx_3.1.xsd` and a Xerces grammar pool; this crate embeds no mzXML schema. `crate::format::mzml_schema` covers mzML only |
| inherited `XMLFile::parse_` / `parseBuffer_` / `save_` | the module's reader and writer | compression sniffing comes from `crate::format::path_io` |
| inherited `XMLFile::schema_location_` / `schema_version_` | `SCHEMA`, `SCHEMA_VERSION` | |
| inherited `XMLFile::enforceEncoding_` | reader's declaration handling | UTF-8, US-ASCII and ASCII-only ISO-8859-1 are accepted |
| inherited `ProgressLogger` (`setLogType`, `startProgress`, `setProgress`, `endProgress`, `nextProgress`) | **not ported** | the source constructor takes a logger and reports per scan; this port has no log stream. Non-fatal observations come back in `ReadReport` instead |

## API mapping — `FORMAT/HANDLERS/MzXMLHandler.h`

The handler is `Internal` and its own documentation says "Do not use this class.
It is only needed in MzXMLFile." There is no Rust type corresponding to it; its
behaviour is the module's reader (`Run`) and writer (`write_engine`). Every
declared member is accounted for.

| C++ member | Rust | Notes |
|---|---|---|
| `MzXMLHandler(MapType&, filename, version, ProgressLogger&)` | `Run::new` (private) | read-only handler |
| `MzXMLHandler(const MapType&, filename, version, const ProgressLogger&)` | `write_engine` (private) | write-only handler; no state object is needed |
| `~MzXMLHandler() override` | derived `Drop` | |
| `LOADDETAIL getLoadDetail() const` | `Detail` (private) | |
| `void setLoadDetail(LOADDETAIL)` | `read_scan_count` selects `Detail::RawCounts` | only `LD_ALL` and `LD_RAWCOUNTS` are reachable for this handler |
| `void onStartElement(qname, attributes)` | `Run::start` (private) | |
| `void onEndElement(qname)` | `Run::end` (private) | |
| `void onCharacters(chars, length)` | `Run::characters` (private) | |
| `void writeTo(std::ostream&)` | `write`, `write_with_options` | |
| `void setOptions(const PeakFileOptions&)` | `ReadOptions::peaks`, `WriteOptions::from_peak_options` | |
| `UInt getScanCount() const` | `ReadReport::scan_count`, `read_scan_count` | |
| `void setMSDataConsumer(IMSDataConsumer*)` | the `consumer` argument of the `transform*` functions | |
| `protected typedef MapType::PeakType PeakType` | `crate::kernel::Peak1D` | |
| `protected typedef MSSpectrum SpectrumType` | `crate::kernel::MSSpectrum` | |
| `protected MapType* exp_` | `Run::experiment` | owned, not borrowed |
| `protected const MapType* cexp_` | the `experiment` argument of `write_engine` | |
| `protected PeakFileOptions options_` | `ReadOptions::peaks` / `WriteOptions` | |
| `protected Int nesting_level_` | `Run::scans` depth | a stack of frames, not a counter — see **Native differences** |
| `protected struct SpectrumData` | `Pending` + `PeaksBlock` (private) | `peak_count_` → `Pending::declared_peaks`, `precision_`/`compressionType_`/`char_rest_` → `PeaksBlock`, `spectrum` → `Pending::spectrum` |
| `SpectrumData::skip_data` | **not ported** | dead upstream: written nowhere, read nowhere |
| `protected std::vector<SpectrumData> spectrum_data_` | `Run::pending` | |
| `protected bool skip_spectrum_` | `Frame::slot == None` | per-scan rather than global |
| `protected UInt spec_write_counter_` | **not ported** | upstream initialises it to 1 and resets it to 1 at the end of `writeTo`; nothing reads it. `write_scans`'s local `written` counter is the renumbering counter that is actually used |
| `protected IMSDataConsumer* consumer_` | `Sink::consumer` (private) | |
| `protected UInt scan_count_` | `ReadReport::scan_count` | |
| `protected const ProgressLogger& logger_` | **not ported** | see the `ProgressLogger` row above |
| `protected writeAttributeIfExists_(os, meta, metakey, attname)` | `write_scan_statistics` (private) | |
| `protected writeUserParam_(os, meta, indent, tag)` | `write_user_param` (private) | |
| `protected doPopulateSpectraWithData_(SpectrumData&)` | `Run::decode` (private) | |
| `protected populateSpectraWithData_()` | `Run::flush` (private) | serial; see **Native differences** |
| `protected std::vector<shared_ptr<DataProcessing>> data_processing_` | `Run::data_processing` (`Vec<Arc<DataProcessing>>`) | |
| `private MzXMLHandler()` | none | upstream declares it and does not implement it |
| `private void init_()` | the five `*_TERMS` constants | |
| file-static `struct IndexPos` | the `(id, offset)` pairs `write_scans` returns | |
| file-static `writeKeyValue(os, key, value)` | inlined in `write_scan_statistics` | |
| namespace `typedef PeakMap MapType` / `typedef MSSpectrum SpectrumType` | `MSExperiment` / `MSSpectrum` | |

### Public Rust surface with no C++ counterpart

`ReadLimits`, `ReadOptions`, `ReadReport`, `WriteLimits`, `WriteOptions`,
`PeakPrecision`, `TransformOptions`, `TransformReport`, `SCHEMA`,
`SCHEMA_VERSION`, `NAMESPACE`, `COMMENT_KEY`, `PHONE_KEY`,
`PROCESSING_TYPE_KEY`, `INTENSITY_CUTOFF_KEY`, `FILTER_STRING_KEY` and
`ATTRIBUTE_METADATA_KEYS`. The first eight exist because this port bounds its
work and returns diagnostics instead of logging them; the key constants name the
`#`-prefixed metadata conventions the source spells as string literals inside
the handler.

## Preserved source conventions

- **`peaksCount` is the pair count.** The decode loop runs to `2 * peak_count_`
  and reads `data[n]`, `data[n+1]` (`MzXMLHandler.cpp:1177-1187`), so a scan
  with `peaksCount="5"` needs ten values. The stale comment at `:300` and the
  `reserve(peak_count_ / 2 + 1)` at `:302` say otherwise and are wrong; only the
  loop is authoritative.
- **`xs:duration` parsing.** The text after the **last** `T` is scanned for `H`,
  `M` and `S` in that order and the components are summed
  (`MzXMLHandler.cpp:254-279`). A day component is therefore dropped, and a
  duration with no `T` and no `H`/`M`/`S` contributes nothing. A missing
  `retentionTime` gives 0.0, not the `MSSpectrum` sentinel `-1`.
- **`windowWideness` is a full width.** It is stored in the lower offset at the
  start tag and both offsets are set to half of it when the m/z text arrives
  (`:219-222`, `:574-579`). An empty `<precursorMz/>` therefore leaves the full
  width in the lower offset and no m/z at all — preserved.
- **Non-fatal attribute defaults.** `precision` defaults to `"32"`, `byteOrder`
  to `"network"`, `contentType` to `"m/z-int"` and `compressionType` to
  `"none"`. `MzXMLFile_1.mzXML` writes `pairOrder=` rather than `contentType=`,
  which upstream never reads, so the default carries it.
- **`scanType` mapping**, including `Q1`, `Q3` and the three non-standard ABI
  Sashimi values `EMS`, `EPI` and `ER`; `EPI` also rewrites the MS level to 2
  (`:370-383`). `Full` selects `MSnSpectrum` above MS level 1.
- **`msLevel="0"`** warns and is read as level 1 (`:248-252`).
- **`centroided`, `chargeDeconvoluted`, `deisotoped`, `collisionEnergy`,
  `basePeakMz`, `lowMz`, `highMz`, `basePeakIntensity`, `totIonCurrent`** are
  written but never read back, as upstream states at `:305`.
- **`<operator>` is singular.** The source `resize(1)` plus `back()` means a
  second `<operator>` overwrites the first; preserved.
- **`#`-prefixed metadata is internal** and is not emitted as `<nameValue>`
  (`:1141`). `#type`, `#intensity_cutoff`, `#phone` and `#comment` keep their
  source spellings.
- **`fileType` is an enumeration in mzXML**, so the writer searches the OpenMS
  free-text type for `raw` and emits `RAWData` or `processedData` (`:662-671`),
  and replaces a checksum that is not 40 characters of SHA-1 with forty zeros
  (`:672-681`).
- **`scanCount` of an empty experiment is one**, because an empty mzXML is not
  schema-valid (`:634`).
- **`startTime`/`endTime` are the first and last spectrum's retention times**,
  not the minimum and maximum (`:637-642`), so a descending file writes an
  `endTime` before its `startTime`.
- **`<index>` spelling.** `<index name = "scan" >` and
  `<offset id = "N" >`, with the spaces around `=` the source emits
  (`:1108-1114`), and `<indexOffset>` giving the byte position of `<index`.
- **MaxQuant compatibility** reproduces every branch of
  `PeakFileOptions::force_mq_compatibility`: empty spectra are skipped, a Thermo
  `<msManufacturer>` is inferred from an Xcalibur acquisition software name, an
  unset scan type becomes `Full`, `lowMz`/`highMz`/`basePeakIntensity`/
  `totIonCurrent` are added, the `<peaks>` tag is broken across lines, a missing
  activation method falls back to `CID`, and an index is written even when it
  was not requested.
- **Base64 leniency.** Upstream's decoder keeps only whole bytes and never
  inspects the surplus bits of the final symbol. `MzXMLFile_3_64bit.mzXML`
  depends on that: its first payload ends `AA1=`, whose two surplus bits are
  set, and `MzXMLFile_test.cpp:313` asserts the resulting intensity
  100.0000991821289. The module therefore decodes `<peaks>` with a
  trailing-bits-permissive, padding-indifferent engine. Whitespace inside the
  payload is stripped as at `:1158-1160`.

## Native differences

Each of these changes observable behaviour and is documented at the item in
`src/format/mzxml.rs`.

1. **Metadata attribution under nesting is per-scan.** Upstream attaches
   `<nameValue>`, `<comment>`, `<precursorMz>` and `<peaks>` to
   `spectrum_data_.back()` — the most recently *started* scan — and nothing pops
   that on `</scan>`. A parent's `<nameValue>` written after a nested child
   therefore lands on the child. This port keeps a stack of open scans, so each
   level keeps its own metadata. `MzXMLFile_5_nested.mzXML` and
   `nested_scan_metadata_stays_on_its_own_scan` pin the corrected behaviour.
   Recorded as a C++ issue candidate.
2. **`skip_spectrum_` is per-scan.** Upstream's single flag stays set after a
   filtered child's `</scan>`, so the parent's remaining children are dropped
   too. A frame's `slot == None` replaces it.
3. **`peaksCount` must match the payload.** Upstream guards the decode loop with
   an `assert` that is compiled out in release builds and then iterates to
   `2 * peak_count_` regardless, reading past the decoded buffer for a short
   payload. This port compares the two and returns `Error::Parse`. An *empty*
   payload still follows the source's early return and yields no peaks with a
   diagnostic, because upstream reaches that return before the count is used.
4. **The four undefined `<peaks>` attribute values are refused.** Upstream logs
   a non-fatal error for an unknown `precision`, a non-`network` `byteOrder`, a
   non-`m/z-int` `contentType` or an unknown `compressionType`, and then decodes
   anyway — reading little-endian data as big-endian, or intensity-first pairs as
   m/z-first. This port returns `Error::Unsupported` and records the source's own
   message in `ReadReport`.
5. **`<peaks>`, `<precursorMz>`, `<software>`, `<processingOperation>` and
   `<msResolution>` outside their expected parent are errors or diagnostics, not
   undefined behaviour.** Upstream indexes `spectrum_data_.back()`,
   `data_processing_.back()` or `getMassAnalyzers()[0]` unconditionally, and
   `*(open_tags_.end() - 2)` even at depth one. Recorded as a C++ issue
   candidate.
6. **A negative `xs:duration` keeps its sign.** Upstream takes the text after
   the last `T` before looking at anything, so `-PT1S` — which its own writer
   emits for the `MSSpectrum` default retention time `-1` — reads back as `+1`.
   `MzXMLFile_4_long.mzXML` is exactly such a file. This port negates, which
   makes a store/load cycle lossless. Recorded as a C++ issue candidate.
7. **Writer defaults are lossless.** `WriteOptions::default` writes
   `precision="64"` and Rust's shortest round-tripping numbers, and keeps a
   fractional `precursorIntensity`. `WriteOptions::source()` selects the source
   behaviour — 32-bit peaks, six significant digits and
   `(int)precursor.getIntensity()` — which narrows every m/z to `f32` and rounds
   retention times and precursor m/z to six digits. This follows the
   `dta::WriteOptions::source()` precedent: the native default refuses to
   discard information and an explicit option selects the source behaviour.
8. **`startTime`/`endTime` are signed.** Upstream writes `PT-1S` there, which
   `xs:duration` does not permit; this port writes `-PT1S`, the same spelling it
   uses for a scan's `retentionTime`.
9. **The XML declaration follows the bytes.** Upstream hard-codes
   `encoding="ISO-8859-1"` and streams `std::string` bytes unchanged, so a UTF-8
   experiment produces a document whose declaration contradicts its content —
   which this crate's own readers then refuse. This port keeps
   `ISO-8859-1` whenever every written string is ASCII and switches to `UTF-8`
   when it is not.
10. **Attribute-backed metadata is not duplicated.** Upstream writes
    `filterLine`, `lowMz`, `highMz`, `basePeakMz`, `basePeakIntensity` and
    `totIonCurrent` as `<scan>` attributes *and* again as `<nameValue>` children,
    because `writeUserParam_` does not skip them. This port skips the six keys
    named in `ATTRIBUTE_METADATA_KEYS`.
11. **MaxQuant mode looks at the next *written* spectrum.** Upstream reads
    `cexp_[s + 1]`'s MS level even when that spectrum is about to be skipped for
    being empty, which can nest an MS1 scan inside another MS1 scan. This port
    skips over the empty ones. Recorded as a C++ issue candidate.
12. **MaxQuant mode refuses an unsorted spectrum.** Upstream logs a non-fatal
    error and then writes `begin()->getMZ()` and `rbegin()->getMZ()` as
    `lowMz`/`highMz` anyway, producing wrong attributes. This port returns
    `Error::UnsortedData` before writing anything.
13. **A malformed `fileSha1` is stored with `ChecksumType::Unknown`.** Upstream
    tags any string as SHA-1; this crate's `SourceFile::validate` requires 40
    hexadecimal characters for that tag, so the text is kept verbatim and the
    algorithm is downgraded. The writer then emits the source's forty-zero
    placeholder, as it does for any non-SHA-1 checksum.
14. **Split text nodes are concatenated before use.** Upstream's `onCharacters`
    applies `setComment` and `setMZ` per chunk, so a `<comment>` split by the
    parser keeps only its last chunk and a split `<precursorMz>` halves the
    isolation window twice. This port accumulates the whole node first.
15. **A duplicated XML attribute, a second `<peaks>` in one scan, a
    `startMz` above its `endMz`, a DTD and a general entity reference are
    refused.** None appears in a conformant document; the last two are refused
    for the same reason `crate::format::mzml` refuses them.
16. **A scan comment lives in metadata.** `MSSpectrum` has no comment field in
    this port, so `<comment>` is stored under `COMMENT_KEY` (`#comment`), the
    same internal key upstream already uses for the instrument comment. Round
    trips are exact because `#`-prefixed keys are not written as `<nameValue>`.
17. **`<nameValue>` order is alphabetical** on write, because metadata is a
    `BTreeMap`; upstream emits `MetaInfoRegistry` index order. Content is
    identical.
18. **Peaks are decoded at `</peaks>`, not at the batch flush.** The result is
    the same — the same filters and the same optional sort — but the base64 text
    of at most one array is retained instead of one per pooled spectrum.
19. **`collisionEnergy` comes from `Precursor::activation_energy`.** This port's
    `Precursor` has no `MetaInfo`, so the source's
    `precursor.metaValueExists("collision energy")` maps to the typed field,
    which is also where `crate::format::mzml` puts `MS:1000509`.
20. **A peak value that is not finite, or an intensity beyond `f32`, is
    refused.** A 64-bit payload can hold an intensity the kernel's `f32` field
    cannot; upstream narrows it with an implicit conversion and stores the
    resulting infinity, which then fails `MSSpectrum::validate` downstream.
21. **Serial.** `populateSpectraWithData_` decodes the batch under
    `#pragma omp parallel for` (`MzXMLHandler.cpp:1224`) and re-throws a single
    `ParseError` for any failure. This port is serial, per
    `docs/REPOSITORY_ANALYSIS.md`, so a large file decodes on one core and the
    error that is reported is the first one encountered rather than an arbitrary
    one.

## Checked boundaries and evidence

### Resource ceilings

`ReadLimits` is checked before any allocation derived from a file-declared
length, following `src/identification/run_mapping.rs` and
`src/format/indexed_mzml_handler.rs`.

| Field | Default | Bounds |
|---|---|---|
| `max_xml_bytes` | 512 MiB | decompressed XML consumed; the input is capped at the ceiling plus one byte so a single oversized event cannot be buffered whole |
| `max_scans` | 1,000,000 | `<scan>` elements, and `msRun`'s declared `scanCount` |
| `max_peaks_per_scan` | 10,000,000 | one scan's `peaksCount`, checked before `try_reserve_exact` |
| `max_total_peaks` | 20,000,000 | declared peaks summed over the file, filtered scans included |
| `max_encoded_bytes` | 128 MiB | retained base64 characters of one `<peaks>` |
| `max_decoded_bytes` | 64 MiB | decoded bytes of one `<peaks>`, before and after inflation |
| `max_depth` | 64 | open elements, which bounds nested `<scan>` recursion |
| `max_source_files` | 100,000 | `<parentFile>` |
| `max_data_processing` | 100,000 | `<dataProcessing>` |
| `max_precursors_per_scan` | 10,000 | `<precursorMz>` |
| `max_metadata_entries` | 1,000,000 | metadata entries across the experiment |
| `max_text_bytes` | 1 MiB | one non-`<peaks>` text node |
| `max_diagnostics` | 1,000 | retained `ReadReport` messages |

`WriteLimits` bounds the writer the same way: `max_spectra`, `max_total_peaks`,
`max_precursors_per_spectrum`, `max_ms_level` (which bounds the indentation the
source derives from the MS level, where `std::string(ms_level + 1, '\t')` would
otherwise allocate gigabytes for a hostile level) and `max_metadata_entries`.

### Atomicity and panic freedom

- The reader returns a value; `load_into` and `transform_into` replace their
  destination only after a complete parse.
- The writer runs a full preflight — counts, MS levels, finiteness of every
  retention time, peak, precursor and scan window, XML validity of every string,
  and the sortedness MaxQuant mode needs — before emitting a byte, and
  `path_io::write` publishes through a temporary file, so a rejected experiment
  leaves both a fresh stream and an existing destination untouched.
- No string that came from a file is byte-sliced. `rsplit_once`/`split`/`rsplit`
  do the duration parsing, and a `completionTime` is clipped with
  `chars().take(19)` rather than the source's `substr(0, 19)`.
  `non_ascii_input_is_handled_without_slicing_bytes` reads a document whose
  metadata, source-file name, filter string, comment and `completionTime` are
  all multi-byte, from a path named `日本語.mzXML`, and writes it back.
- Every index into file-derived data goes through `get`, `chunks_exact` or a
  checked length comparison; arithmetic on declared lengths uses `checked_*`.

### Evidence

**Tier 3, source review.** Expectations come from `MzXMLFile_test.cpp` literals
and from the four unmodified upstream fixtures. No C++ was built or executed and
no C++ output was retained, so this is **not** a tier 1 differential.

All 15 `START_SECTION`s of `MzXMLFile_test.cpp` are ported; none is merely
mapped. `tests/mzxml.rs` carries 59 cases. The one substitution is
`[EXTRA] static bool isValid(...)` (1 assertion): Xerces XSD validation is
unavailable here, so `the_stored_document_is_structurally_valid` asserts instead
that the writer's output carries the 3.1 namespace and schema location, has
balanced `<scan>` elements, and reloads to the same experiment.

`MzXMLFile_4_long.mzXML` (10,641,595 bytes) is **not** copied into this
repository; its sha256 is recorded in the manifest. An equivalent document with
the same 997530 peaks is generated in
`load_reads_a_very_long_spectrum_split_across_text_events`, with the base64
wrapped at 76 characters so the payload arrives as thousands of text events. The
peak count is the upstream literal; the coordinates are independently derived.

Everything else — the resource ceilings, the attribute rejections, the writer's
lossless defaults, the nesting attribution, the retention-time sign, the
non-ASCII inputs and the malformed documents — is tier 4, independently derived,
because no upstream fixture reaches it.

## PeakFileOptions fields this format does not use

`force_tpp_compatibility`, `write_supplemental_data`, `mz_32_bit`,
`intensity_32_bit`, `skip_chromatograms`, `sort_chromatograms_by_rt` and the
three Numpress configurations are mzML- or mzData-only; mzXML has no
chromatograms, one shared `precision` attribute for both halves of a pair, and
no Numpress. `skip_xml_checks` bypasses only upstream's base64 whitespace
removal; this port always strips whitespace, because stripping is what makes
`MzXMLFile_2_minimal.mzXML` readable and the cost is one pass. `zlib_compression`
has no effect in the source mzXML writer and is honoured here as the native
extension described above.

## Deferred

- `XMLFile::isValid` against `mzXML_idx_3.1.xsd`: no mzXML schema is embedded and
  `crate::format::mzml_schema` is mzML-only.
- `ProgressLogger`: the source constructor threads a logger through and reports
  progress per scan and per stored spectrum. `ReadReport` carries the non-fatal
  observations; there is no progress callback.
- `FileHandler` dispatch: `crate::format::file_types::FileType::MzXml` already
  exists and content detection already recognises `<mzXML`, but
  `src/format/file_handler.rs` does not route to this module. That file is
  outside this package.
- The mzML metadata mzXML cannot express — sample description, HPLC gradient,
  instrument configurations, multiple ion sources or analyzers, chromatograms,
  auxiliary data arrays, ion mobility, peptide identifications — is lost by the
  format, not by the port. A `store` of an mzML-derived experiment drops it, as
  upstream does.
- Parallel batch decoding, deliberately: see native difference 21.

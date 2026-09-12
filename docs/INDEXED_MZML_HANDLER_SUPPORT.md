# Indexed mzML record handler

Rust: `src/format/indexed_mzml_handler.rs` (feature `mzml`)
Tests: `tests/indexed_mzml_handler.rs`
Manifest: `tests/data/indexed_mzml_handler_provenance.json`

Source headers, pinned at `bc9cc12514c768385ce121d6ca4bb710fe1983c4`:

- `src/openms/include/OpenMS/FORMAT/HANDLERS/IndexedMzMLHandler.h` — ported here
- `src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h` — behaviourally
  absorbed; see the second table
- `src/openms/include/OpenMS/KERNEL/OnDiscMSExperiment.h` — read for the
  consumer contract only; the facade is a separate work package

The footer index itself is `src/format/indexed_mzml.rs`
(`docs/INDEXED_MZML_SUPPORT.md`); this module consumes it.

## What the module does

`IndexedMzMLHandler::open` parses the `indexListOffset` footer and the
`<indexList>` sections with `IndexedMzMLDecoder`, then caches three things:

1. the two ordered offset vectors and an ordered native-id lookup for each,
2. the raw bytes of the document from byte zero up to the first record list's
   opening tag, and
3. the `defaultDataProcessingRef` of each record list's opening tag, located by a
   bounded backwards scan from that list's first record.

Fetching record *i* then reads only that record's byte range, trims it at the
record's own closing tag, and assembles

```text
<cached header bytes><spectrumList count="1" defaultDataProcessingRef="…">
  <the one record>
</spectrumList></run></mzML></indexedmzML>
```

which goes to `mzml::read_with_load_options`. The closing suffix is derived from
the elements still open at the end of the cached header, so it matches whatever
wrapper and namespace prefixes the file actually uses.

There is no second XML parser. The consequence is that a record fetched here
carries everything the whole-file reader produces — retention time, MS level,
polarity, precursors, products, scan settings, float/integer/string data arrays
and metadata — where the source decodes only the `binaryDataArray` payloads and
the `id` attribute.

## API mapping — `IndexedMzMLHandler.h`

Every public member of the header appears here.

| C++ member | Rust | Notes |
|---|---|---|
| `IndexedMzMLHandler()` | not ported | The default constructor leaves `parsing_success_` false, a state the header itself calls invalid to use. There is no unopened handler here. |
| `explicit IndexedMzMLHandler(const std::string&)` | `IndexedMzMLHandler::open(path)` | Returns `Result<Self>`; a handler exists only when its index parsed. |
| — | `IndexedMzMLHandler::open_with_limits(path, RecordReadLimits, mzml::ReadOptions)` | Native. Explicit ceilings and mzML decoding limits. |
| `IndexedMzMLHandler(const IndexedMzMLHandler&)` | not ported | The source copy reopens the file "critical for parallel access". A handler owns a seek position, so an independent reader is an independent `open`. The source copy also silently loses both native-id maps (issue below). |
| `~IndexedMzMLHandler()` | `Drop` on the owned `File` | The source destructor is `= default`. |
| `void openFile(const std::string&)` | `IndexedMzMLHandler::open` | No reopen on an existing handler: `parseFooter_` never clears its containers, so a second `openFile` on one object appends (issue below). |
| `bool getParsingSuccess() const` | not ported as a predicate | It is the `Ok`/`Err` of `open`. |
| `size_t getNrSpectra() const` | `spectrum_count()`, also `count(RecordKind::Spectrum)` | |
| `size_t getNrChromatograms() const` | `chromatogram_count()`, also `count(RecordKind::Chromatogram)` | |
| `Interfaces::SpectrumPtr getSpectrumById(int)` | not ported | The crate has no OpenSWATH `Interfaces` layer. `spectrum(i)` returns the same peak data in `MSSpectrum`. |
| `const MSSpectrum getMSSpectrumById(int)` | `spectrum(index) -> Result<Option<MSSpectrum>>` | `None` means the configured `PeakFileOptions` excluded the whole record. |
| `void getMSSpectrumById(int, MSSpectrum&)` | `spectrum(index)` | The out-parameter overload returns the value instead. The source overload *merges* into the caller's spectrum, keeping metadata the caller already had; `OnDiscMSExperiment` relies on that, and the facade work package reproduces it by overlaying this return value on its own metadata record. |
| `void getMSSpectrumByNativeId(const std::string&, MSSpectrum&)` | `spectrum_by_native_id(&str) -> Result<Option<MSSpectrum>>` | Unknown id is `Error::InvalidValue`, the source's `Exception::IllegalArgument`. |
| `Interfaces::ChromatogramPtr getChromatogramById(int)` | not ported | As the spectrum pointer overload. |
| `const MSChromatogram getMSChromatogramById(int)` | `chromatogram(index) -> Result<Option<MSChromatogram>>` | |
| `void getMSChromatogramById(int, MSChromatogram&)` | `chromatogram(index)` | As the spectrum out-parameter overload. |
| `void getMSChromatogramByNativeId(const std::string&, MSChromatogram&)` | `chromatogram_by_native_id(&str)` | |
| `void setSkipXMLChecks(bool)` | `options_mut().skip_xml_checks` | Reaches the decoder: `PeakFileOptions` is passed into `mzml::read_with_load_options`, which suppresses the four-character Base64 whitespace strip (`src/format/mzml.rs:1335`, `:2386`). That is the whole effect in the source too — `skip_xml_checks_` is forwarded to `MzMLHandlerHelper::decodeBase64Arrays`, never to XML syntax checking. |
| `filename_` (private) | `path()` | |
| `spectra_offsets_`, `chromatograms_offsets_` (private) | `offset(kind, index)` | |
| `spectra_native_ids_`, `chromatograms_native_ids_` (private) | `index_of(kind, native_id)`, `native_id(kind, index)` | `BTreeMap`, so lookups are ordered and reproducible; the source uses `unordered_map`. |
| `index_offset_` (private) | `index_list_offset()` | |
| `spectra_before_chroms_` (private) | `spectra_before_chromatograms()` | |
| `filestream_` (private) | owned `File`; every fetch takes `&mut self` | |
| `parsing_success_` (private) | not represented | See `getParsingSuccess`. |
| `skip_xml_checks_` (private) | `PeakFileOptions::skip_xml_checks` | See `setSkipXMLChecks`. |
| `parseFooter_()` (private) | inside `open_with_limits` | |
| `getSpectrumById_helper_`, `getChromatogramById_helper_` (private) | `record_xml(kind, index) -> Result<Vec<u8>>` | Public here, and trimmed at the record's closing tag. |

Native additions with no source counterpart: `RecordKind`, `RecordReadLimits`,
`limits()`, `read_options()`, `options()`, `options_mut()`, `set_options()`,
`is_empty()`.

`options()` / `set_options()` carry the names of `OnDiscMSExperiment::getOptions`
and `setOptions`; in the source the options live one layer up, in the facade.

## API mapping — `MzMLSpectrumDecoder.h`

This header has no ledger entry of its own here. Its job — turn one
`<spectrum>`/`<chromatogram>` text into a record — is done by
`src/format/mzml.rs` through the envelope described above.

| C++ member | Rust | Notes |
|---|---|---|
| `explicit MzMLSpectrumDecoder(bool skip_xml_checks = false)` | not ported as a type | The decoder is `mzml::read_with_load_options` over an assembled document; the constructor flag is `PeakFileOptions::skip_xml_checks`. |
| `void domParseSpectrum(const std::string&, Interfaces::SpectrumPtr&)` | not ported | No `Interfaces` layer. |
| `void domParseSpectrum(const std::string&, MSSpectrum&)` | covered by `IndexedMzMLHandler::spectrum` | Strictly more is decoded: the source reads only `id`, `defaultArrayLength` and `binaryDataArray`. |
| `void domParseChromatogram(const std::string&, MSChromatogram&)` | covered by `IndexedMzMLHandler::chromatogram` | |
| `void domParseChromatogram(const std::string&, Interfaces::ChromatogramPtr&)` | not ported | No `Interfaces` layer. |
| `void setSkipXMLChecks(bool)` | `options_mut().skip_xml_checks` | As above; the decoder's copy of the flag is the same option value. |
| `domParseString_`, `handleBinaryDataArray_`, `decodeBinaryData*_` (protected) | `src/format/mzml.rs` | The whole-file reader's own binary-array handling. |

Two source checks in `MzMLSpectrumDecoder.cpp:25-50` (`checkData_`) deserve
naming, because the replacement must not lose them and does not. Integer-encoded
m/z, RT or intensity arrays are rejected there; `src/format/mzml.rs:462` rejects
any primary peak array that is not floating point. Primary arrays of unequal
length are rejected there; `src/format/mzml.rs:592` is stricter still, requiring
every decoded array to match the record's `defaultArrayLength`, so two primary
arrays cannot disagree in the first place.

`domParseString_` carries a `@pre` that the input must have `<spectrum>` or
`<chromatogram>` as its root element, but in the source that precondition is
only an `OPENMS_PRECONDITION`, which `CONCEPT/Macros.h:91` expands to nothing
outside debug builds. Here `record_xml` checks the element name unconditionally,
before anything is decoded.

## Preserved source conventions

- **Record byte ranges.** Record *i* spans `[offset[i], offset[i+1])`. The last
  record of a kind ends at the first record of the other kind when that kind
  follows it in the file, and at `indexListOffset` otherwise. This reproduces
  the branch at `IndexedMzMLHandler.cpp:142` and `:202` exactly, including the
  `spectra_before_chroms_` rule below.
- **`spectra_before_chroms_`.** Defaults to true, and is computed only when both
  index sections are non-empty, by comparing the two first offsets
  (`IndexedMzMLHandler.cpp:48`).
- **Duplicate native identifiers.** `parseFooter_` uses `unordered_map::emplace`,
  which does not overwrite, so the *first* entry wins. `index_of` does the same.
- **Index order.** Offsets stay in the order the index lists them; no sorting.
- **Out-of-range access is an error**, not a clamp or an empty record.
- **Unknown native identifier is an error** (`Exception::IllegalArgument` →
  `Error::InvalidValue`).
- **`PeakFileOptions` semantics** follow `OnDiscMSExperiment::getSpectrum`: RT
  range, MS level and precursor m/z decide the record; m/z and intensity ranges
  select peaks within it. `DRange::encloses` half-open endpoints are preserved
  by `src/format/mzml_load.rs`.

## Native differences

- **Failure is a `Result`, not a flag.** `open` returns `Err` where the source
  records `parsing_success_ == false`, so no handler can be used in the state
  the header calls invalid.
- **Filtered-out records are `None`.** `OnDiscMSExperiment::getSpectrum` returns
  the metadata-only spectrum for an excluded record, explicitly to preserve the
  caller's index mapping. Here the index mapping is the caller's own argument,
  so absence is reported as absence rather than as an empty record.
- **The record range is trimmed at the record's closing tag.** The source hands
  the untrimmed range — including `</spectrumList>` and the following
  `<chromatogramList …>`, or `</run></mzML>` — to a `XercesDOMParser` created
  without an error handler, so the resulting well-formedness errors are
  discarded. Trimming means no parser is ever fed trailing content, and a
  genuinely unterminated record is an error instead of a silent partial read.
- **Full record metadata.** As described above.
- **The index is verified against the record.** After decoding, the record's `id`
  attribute must equal the native identifier the index recorded for that offset.
  The source never compares them.
- **The offset must point at the element.** `record_xml` requires the range to
  begin with `<spectrum`/`<chromatogram` followed by whitespace, `>` or `/`. The
  source seeks wherever the index says.
- **Ordered maps.** `BTreeMap` instead of `unordered_map`: same first-wins rule,
  reproducible iteration.
- **Encoding.** The assembled document keeps the original `<?xml …?>`
  declaration from the cached header, so an `ISO-8859-1` or `US-ASCII` document
  is decoded as the file declares. The source parses the extracted record text
  through `MemBufInputSource` with no declaration at all, which Xerces treats as
  UTF-8 regardless of what the file said.
- **No threads.** The source parallelises nothing in this class, but its own
  class documentation and `OnDiscMSExperiment`'s recommend copying the object per
  OpenMP thread. This port is serial; concurrency is a caller-side decision and
  costs one extra `open` per reader. The performance gap is stated rather than
  closed.
- **Per-fetch header cost.** Each fetch re-parses the cached header bytes so that
  `dataProcessingRef`, `sourceFileRef`, `instrumentConfigurationRef` and
  `referenceableParamGroupRef` on the record resolve. The source skips the header
  entirely and therefore also skips every reference. The header is read once and
  bounded by `max_header_bytes`; the extra work per record is proportional to the
  header, not to the file.

## Checked boundaries and evidence

`RecordReadLimits` (all checked before any allocation or read):

| Field | Default | Guards |
|---|---|---|
| `max_record_bytes` | 256 MiB | the span of one record's byte range |
| `max_header_bytes` | 16 MiB | the cached document prefix |
| `max_list_tag_bytes` | 64 KiB | the backwards scan for a record list's opening tag |
| `max_records` | 1 000 000 | index entries per kind |
| `max_native_id_bytes` | 64 MiB | stored identifiers, counted for both the vector and the map |
| `max_depth` | 64 | open elements in the cached header, bounding the closing suffix |
| `index` | `IndexReadLimits::default()` | handed to the footer decoder |

Additional unconditional checks: `indexListOffset` must not exceed the file
length; every index entry must lie within the file; a record range must not
decrease; the record element must be present at the offset and closed inside the
range; the reads use `try_reserve_exact` and verify the byte count actually read.
Every failure returns before mutating the handler, so a rejected fetch leaves the
index and the cached header untouched — the tests assert this for the
`max_record_bytes` case.

**Evidence tier 3 (source review).** Literals transcribed from
`IndexedMzMLFile_test.cpp` and its unmodified fixture: 2 spectra, 1 chromatogram,
`indexListOffset` 667742, record offsets 24146 / 345745 / 665563, native ids
`controllerType=0 controllerNumber=1 scan=1`, `…scan=2` and `TIC`, 19914 / 19800 /
48 points, MS level 1, scan start times 0.2961 s and 0.4738 s. Where the class
test compares the handler against `MzMLFile().load` of the same file, the port
compares against this crate's own whole-file reader on the same file.

**Evidence tier 4 (independently derived).** The chromatograms-before-spectra
ordering branch, the hostile-index rejections, the resource ceilings and the
`PeakFileOptions` behaviour use synthetic documents built in the test, because
the upstream fixture cannot reach those paths. One test loads each synthetic
document with the whole-file reader first, so the handler is never the only
reader that accepts it.

No C++ was built or executed and no retained C++ output was used, so nothing here
is tier 1 or tier 2. An oracle driver over `IndexedMzMLHandler` would be the way
to reach tier 1, and is recorded as a deferral.

## Class-test section coverage

All 16 `START_SECTION` blocks of `IndexedMzMLFile_test.cpp` are ported; none are
merely mapped.

| Section | Rust test |
|---|---|
| `IndexedMzMLHandler(std::string filename)` | `constructor_from_filename_parses_the_index` |
| `~IndexedMzMLHandler()` | `dropping_a_handler_releases_the_file` |
| `IndexedMzMLHandler()` | `there_is_no_unopened_handler` |
| `IndexedMzMLHandler(const IndexedMzMLHandler&)` | `a_second_handler_on_one_file_reads_the_same_data` |
| `bool getParsingSuccess() const` | `parsing_success_is_the_result_of_open` |
| `void openFile(std::string)` | `opening_replaces_rather_than_accumulates` |
| `size_t getNrSpectra() const` | `spectrum_count_matches_the_index` |
| `size_t getNrChromatograms() const` | `chromatogram_count_matches_the_index` |
| `Interfaces::SpectrumPtr getSpectrumById(int)` | `spectrum_arrays_match_a_whole_file_load` |
| `MSSpectrum getMSSpectrumById(int)` | `spectrum_by_index_carries_full_record_metadata` |
| `void getMSSpectrumByNativeId(std::string, MSSpectrum&)` | `spectrum_by_native_id_resolves_through_the_index` |
| `Interfaces::ChromatogramPtr getChromatogramById(int)` | `chromatogram_arrays_match_a_whole_file_load` |
| `MSChromatogram getMSChromatogramById(int)` | `chromatogram_by_index_carries_its_native_id` |
| `void getMSChromatogramByNativeId(std::string, MSChromatogram&)` | `chromatogram_by_native_id_resolves_through_the_index` |
| `[EXTRA] load broken file` (2^64) | `footer_offset_beyond_sixty_three_bits_is_rejected` |
| `[EXTRA] load broken file` (2^63-1) | `footer_offset_past_end_of_file_is_rejected` |

## Source defects found while porting

Reported to the integrator for `OpenMS_CPP_ISSUES.md`; numbering is assigned
there.

1. **The copy constructor loses both native-id maps.**
   `IndexedMzMLHandler.cpp:76-88` initialises `filename_`, both offset vectors,
   `index_offset_`, `spectra_before_chroms_`, a fresh `filestream_`,
   `parsing_success_` and `skip_xml_checks_`, but not `spectra_native_ids_` or
   `chromatograms_native_ids_`, which are therefore default-constructed empty.
   Every `getMSSpectrumByNativeId` / `getMSChromatogramByNativeId` on a copy
   throws `Exception::IllegalArgument` for every identifier. The class exists to
   be copied — its own comment at line 82 says reopening rather than copying the
   stream "is critical for parallel access to the same file", and
   `OnDiscMSExperiment`'s copy constructor copies the handler by value, so the
   `#pragma omp parallel for firstprivate(ondisc_map)` pattern that
   `OnDiscMSExperiment.h:66` recommends produces per-thread handlers whose
   `getSpectrumByNativeId` cannot work.

2. **`openFile` accumulates instead of replacing.** `parseFooter_`
   (`IndexedMzMLHandler.cpp:37-46`) pushes into `spectra_offsets_`,
   `chromatograms_offsets_` and both maps without clearing them, and `openFile`
   (`:92`) does not clear them either. Opening two valid indexed files on one
   object leaves `getNrSpectra()` reporting the sum, with offsets from the first
   file read through the second file's stream. The class test does call
   `openFile` repeatedly (`IndexedMzMLFile_test.cpp:100-106`) but its first two
   calls fail before any append, so the defect is not observed there.

3. **The record read length is unchecked in both directions.**
   `IndexedMzMLHandler.cpp:162-167` and `:222-227` compute
   `std::streampos readl = endidx - startidx` and pass it to
   `new char[readl + 1]`. Nothing checks that `endidx >= startidx`, that either
   lies inside the file, or that the difference is a sane size — all three come
   straight from the file's own footer index. A decreasing pair makes the length
   negative, so the `new[]` throws `std::bad_array_new_length`, a non-OpenMS
   exception that escapes past every `Exception::` handler a caller installed; an
   oversized pair is an unbounded allocation. The subsequent
   `filestream_.read` result is never inspected, so when the range runs past the
   end of the file the buffer keeps uninitialised bytes and only `buffer[readl]`
   is set to `'\0'`, after which `std::string text(buffer)` reads that
   uninitialised memory.

4. **The chromatogram bound error reports the spectrum count.**
   `IndexedMzMLHandler.cpp:134-136`: `getChromatogramById_helper_` checks against
   `getNrChromatograms()` but its message reads "id needs to be smaller than the
   number of spectra" and interpolates `getNrSpectra()`. On the class test's own
   fixture (2 spectra, 1 chromatogram) chromatogram index 1 is rejected with
   "maximal allowed is 2". Both helpers additionally say "maximal allowed is N"
   where the largest accepted index is N-1.

5. **Record XML errors are discarded.**
   `MzMLSpectrumDecoder.cpp:538-543` constructs a `XercesDOMParser` and calls
   `parse()` without `setErrorHandler`, so errors go to the default handler and
   are dropped. This is load-bearing rather than incidental: the handler's byte
   range for the last record of each kind deliberately includes
   `</spectrumList>` and the next list's opening tag (or `</run></mzML>`), so
   that record is *always* parsed from input that is not well-formed. The same
   silence turns a truncated or corrupt record into an empty spectrum rather than
   a reported error.

6. **`getMSChromatogramById(int)` updates ranges twice.**
   `IndexedMzMLHandler.cpp:278-284` calls the two-argument overload, which
   already ends with `c.updateRanges()` at `:301`, and then calls
   `c.updateRanges()` again — a second full pass over the chromatogram. The
   spectrum counterpart at `:246-251` does not.

## Deferrals

- `OnDiscMSExperiment` itself (the container facade, its metadata load through
  `FileHandler` with `fillData=false`, `isSortedByRT`, `operator==`, the
  `OnDiscPeakMap` alias) is a separate work package. This module deliberately
  exposes what that facade needs: counts, per-index and per-native-id fetch,
  offsets, identifiers and the options.
- `IndexedMzMLFileLoader` is not ported.
- The `Interfaces::Spectrum` / `Interfaces::Chromatogram` pointer overloads are
  not ported; the crate has no OpenSWATH `Interfaces` layer.
- Tier 1 evidence would need an oracle driver over `IndexedMzMLHandler` in
  `../oracle/drivers/`, printing full-precision arrays for each record of the
  fixture. Not built here.

# mzData 1.05 support

Rust port of `FORMAT/MzDataFile.h` (87 lines) and its
`FORMAT/HANDLERS/MzDataHandler.h` (191 lines, with a 1510-line `.cpp`) at
source revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

| | |
|---|---|
| Rust | `src/format/mzdata.rs` |
| Tests | `tests/mzdata.rs` |
| Provenance | `tests/data/mzdata_provenance.json` |
| Feature | `mzml` — the module reuses the `base64` crate that feature already pulls in; it does **not** use the mzML reader, writer or schema code |
| Evidence | Tier 3 (transcribed class-test literals) plus one independent cross-check; no C++ was built or executed |

mzData predates mzXML and mzML. A document is a `<description>` block of
experiment metadata followed by a `<spectrumList>` of `<spectrum>` elements,
each with a `<spectrumDesc>` and one `<mzArrayBinary>` / `<intenArrayBinary>`
pair plus any number of `<supDataArrayBinary>` annotation arrays. Metadata is
`<cvParam accession="PSI:1000xxx" value="…"/>` against a vocabulary the handler
hard-codes in `init_` rather than loading from an OBO file, plus free
`<userParam name value/>` pairs.

The property that distinguishes it from every later PSI format: **each binary
array declares its own `precision`, `endian` and `length` as XML attributes**,
and the source honours all three. `endian="big"` is read for the m/z array, the
intensity array and every annotation array independently, so one spectrum may
mix both byte orders. `tests/mzdata.rs::big_endian_arrays_are_honoured` and
`byte_order_changes_the_decoded_values` hold that behaviour down, and
`big_endian_input_survives_a_round_trip` stores a big-endian document back out
— little endian, as `writeBinary_` writes — and asserts every value returns.

---

## API mapping

### `FORMAT/MzDataFile.h`

| C++ member | Rust counterpart |
|---|---|
| `typedef PeakMap MapType` (private) | `crate::kernel::MSExperiment` throughout; no alias is introduced |
| `MzDataFile()` | `MzDataFile::new` / `MzDataFile::default`. The source constructor also pins `XMLFile("/SCHEMAS/mzData_1_05.xsd", "1.05")`; those are the constants `SCHEMA_FILE` and `SCHEMA_VERSION` |
| `~MzDataFile() override` | not ported: `Drop` is derived. The source destructor is `= default` |
| `PeakFileOptions& getOptions()` | `MzDataFile::options_mut` |
| `const PeakFileOptions& getOptions() const` | `MzDataFile::options` |
| `void setOptions(const PeakFileOptions&)` | `MzDataFile::set_options` |
| `void load(const std::string&, MapType&)` | `MzDataFile::load` (returns the experiment), `MzDataFile::load_into` (replaces a destination), `MzDataFile::load_report` (also returns `LoadReport`). Free forms: `load`, `load_with_options`, `load_into`, `read`, `read_with_options` |
| `void store(const std::string&, const MapType&) const` | `MzDataFile::store`, `MzDataFile::store_report`. Free forms: `store`, `store_with_options`, `store_report_with_options`, `write`, `write_with_options` |
| `bool isSemanticallyValid(const std::string&, StringList& errors, StringList& warnings)` | **not ported**: `MzDataFile::is_semantically_valid` returns `Error::Unsupported`. It needs `/MAPPING/mzdata-mapping.xml` and `/CV/psi-mzdata.obo`, neither of which ships with this crate, and `MzDataFile_test.cpp:849-854` marks its own section `NOT_TESTABLE` because "the mapping file was hand-crafted by Marc Sturm". The generic machinery is `crate::format::semantic_validator` behind the `semantic-validation` feature |
| `PeakFileOptions options_` (private) | private `MzDataFile::options` field |
| inherited `Internal::XMLFile::isValid(filename, os)` | **not ported**: XSD validation. `crate::format::mzml_schema` is mzML-specific and `mzData_1_05.xsd` is not shipped. `tests/mzdata.rs::stored_documents_are_wellformed` reproduces what can be reproduced — both documents are well-formed and reload |
| inherited `Internal::XMLFile::getVersion()` | `MzDataFile::version`, and the `SCHEMA_VERSION` constant |
| inherited `XMLFile::parse_` / `save_` / `schema_location_` / `schema_version_` (protected) | replaced by the free `read_*` / `write_*` functions and the `SCHEMA_*` constants. `store_with_options` publishes through `crate::format::path_io::write`, which is this crate's equivalent of `save_`: a sibling temporary renamed onto the destination, with `.gz`/`.bz2` suffix compression |
| inherited `XMLFile::parseBuffer_` (protected) | `read` / `read_with_options`, which take any `BufRead` and so cover both the file and the in-memory case; `parse_` is the only one `MzDataFile` itself calls |
| inherited `XMLFile::enforceEncoding_` / `enforced_encoding_` (protected) | **not ported**: the override exists for X!Tandem output whose declaration the Xerces parser stumbles on, and no `MzData*` code path sets it. The declared encoding is honoured instead — ISO-8859-1, UTF-8 and US-ASCII are read and anything else is refused rather than misread |
| inherited `ProgressLogger` (`setLogType`, `getLogType`, `startProgress`, `setProgress`, `endProgress`, `nextProgress`) | **not ported here**: `crate::concept::progress_logger` exists but is not threaded through this module. The counters the handler spends on progress are returned instead, in `LoadReport` and `StoreReport` |

### `FORMAT/HANDLERS/MzDataHandler.h`

The handler is a SAX subclass whose class documentation says "Do not use this
class. It is only needed in MzDataFile." Its observable behaviour is the whole
of mzData reading and writing, so all of it is ported — as a one-pass,
value-returning parser and a writer, not as a subclass. The three
`typedef`s the header declares at namespace scope, above the class, are rows of
the table like any other member.

| C++ member | Rust counterpart |
|---|---|
| `typedef PeakMap MapType` (namespace scope) | `crate::kernel::MSExperiment`; no Rust alias is introduced |
| `typedef MSSpectrum SpectrumType` (namespace scope) | `crate::kernel::MSSpectrum`; no Rust alias is introduced |
| `typedef MSChromatogram ChromatogramType` (namespace scope) | **not ported**: the alias is declared and never used — mzData 1.05 has no chromatogram element, `MzDataHandler` never mentions a `ChromatogramType` value, and neither reading nor writing has anything to map it onto. The writer refuses `MSExperiment::chromatograms` (see **What the writer refuses by default**) rather than silently dropping them, which is where a caller meets the gap |
| `MzDataHandler(MapType& exp, filename, version, ProgressLogger&)` | the reading path: `read_with_options` and the private `Parser`. The header comments on the two constructors are **swapped**: this one assigns `exp_` and is used by `load`, but is documented "Constructor for a write-only handler" |
| `MzDataHandler(const MapType& exp, filename, version, const ProgressLogger&)` | the writing path: `write_with_options`. Documented "Constructor for a read-only handler", and it is the store constructor |
| `~MzDataHandler() override` | not ported: empty body, `Drop` is derived |
| `void onStartElement(qname, attributes) override` | private `Parser::start` |
| `void onEndElement(qname) override` | private `Parser::end` |
| `void onCharacters(chars, length) override` | private `Parser::push_text` plus the per-element dispatch in `Parser::end`. The source applies character data at every chunk, overwriting the previous chunk for every tag except `data`; this port accumulates per element and applies once at the closing tag, so a value split across chunks is not truncated |
| `void writeTo(std::ostream&) override` | private `write_document` / `write_description` / `write_spectra`, reached through `write_with_options` |
| `void setOptions(const PeakFileOptions&)` | `MzDataFile::set_options`, and the `options` argument of `read_with_options` / `load_with_options` |
| `void init_()` (private) | the module constants `SAMPLE_STATE`, `IONIZATION_MODE`, `RESOLUTION_METHOD`, `RESOLUTION_TYPE`, `SCAN_DIRECTION`, `SCAN_LAW`, `REFLECTRON_STATE`, `ACQUISITION_MODE`, `IONIZATION_TYPE`, `INLET_TYPE`, `DETECTOR_TYPE`, `ANALYZER_TYPE`, `ACTIVATION_METHOD` — the 13 non-empty of the 19 `cv_terms_` tables. Tables 4 (`ScanFunction`), 12 (`TandemScanningMethod`), 15 (`EnergyUnits`), 16 (`ScanMode`) and 17 (`Polarity`) are commented "no longer used" and left empty upstream, so they have no Rust counterpart |
| `typedef MapType::PeakType PeakType` (protected) | `crate::kernel::Peak1D` |
| `typedef MSSpectrum SpectrumType` (protected) | `crate::kernel::MSSpectrum` |
| `MapType* exp_` (protected) | the owned `Parser::experiment`, returned by value |
| `const MapType* cexp_` (protected) | the `&MSExperiment` argument of `write_with_options` |
| `PeakFileOptions options_` (protected) | `Parser::options` (a borrow) and `WriteOptions` |
| `UInt peak_count_` (protected) | the declared `length` of the m/z array, kept per `<data>` element in `Encoded::declared`. The source keeps it as a handler member that is never reset per spectrum |
| `SpectrumType spec_` (protected) | `Parser::spectrum` |
| `std::vector<std::pair<std::string, MetaInfoDescription>> meta_id_descs_` (protected) | `Parser::descriptions`, a `Vec<Description { reference, metadata }>` |
| `std::vector<std::string> data_to_decode_` (protected) | `Encoded::payload`, one per `<data>` element in `Parser::arrays` |
| `std::vector<float> data_to_encode_` (protected) | a local `Vec<f64>` per array in `write_spectra`; `Sink::binary` narrows to the declared precision at the point of encoding |
| `std::vector<std::vector<float>> decoded_list_` (protected) | not ported as state: `Parser::decode` returns one `Vec<f64>` per array and the 32/64-bit split never needs two parallel vectors |
| `std::vector<std::vector<double>> decoded_double_list_` (protected) | as above |
| `std::vector<std::string> precisions_` (protected) | `Encoded::precision`, typed as `Precision` |
| `std::vector<std::string> endians_` (protected) | `Encoded::endian`, typed as `Endian` |
| `bool skip_spectrum_` (protected) | `Parser::skip` |
| `const ProgressLogger& logger_` (protected) | not ported; see the `ProgressLogger` row above |
| `void fillData_()` (protected) | private `Parser::fill_data` |
| `writeCVS_(os, double value, acc, name, indent=4) const` | private `Sink::cv_number`. Nothing is written when the value is exactly zero, as upstream |
| `writeCVS_(os, const std::string& value, acc, name, indent=4) const` | private `Sink::cv_string`. Nothing is written for an empty value, as upstream |
| `writeCVS_(os, UInt value, UInt map, acc, name, indent=4)` | private `Sink::cv_enum` plus `cv_term`. A map or value index outside the table is skipped; upstream warns, and here the table is a compile-time constant so the index cannot be out of range |
| `writeUserParam_(os, const MetaInfoInterface& meta, indent=4)` | private `Sink::user_params`. Keys whose first character is `#` are skipped, as upstream |
| `cvParam_(const std::string& accession, const std::string& value)` | private `Parser::cv_param`, split per section into `cv_spectrum_instrument`, `cv_ion_selection`, `cv_activation`, `cv_detector`, `cv_source`, `cv_sample`, `cv_analyzer`, `cv_additional` and `cv_processing`. The declared parameter name is `name` in the header and `accession` in the `.cpp`; it is an accession |
| `writeBinary_(os, size, tag, name="", id=-1)` | private `Sink::binary` |
| `std::shared_ptr<DataProcessing> data_processing_` (protected) | `Parser::processing`, cloned into an `Arc` per spectrum |
| header `@improvement` "Add implementation and tests of 'supDataArray' to store IntegerDataArray and StringDataArray" | still open. `MSSpectrum::integer_data_arrays` and `string_data_arrays` have no mzData element, and the writer refuses them unless `WriteOptions::discard_unrepresentable` is set |

### Native additions with no source counterpart

| Rust item | Why |
|---|---|
| `SCHEMA_VERSION`, `SCHEMA_FILE`, `SCHEMA_LOCATION` | the three literals the source constructor and writer embed |
| `Endian`, `Precision` | the two per-array transport attributes, typed instead of compared as strings. `parse` reproduces the source's exact tolerance and reports whether the spelling was one of the two it knows |
| `ReadLimits` | explicit ceilings. The source has none; see **Checked boundaries** |
| `LoadReport`, `Loaded` | the source reports everything through `XMLHandler::warning`, which `#ifdef`s down to a debug log in a release build, so a caller cannot see it |
| `WriteOptions`, `WriteOptions::source`, `StoreReport` | the refuse-or-discard policy for everything mzData cannot represent |
| `MzDataFile::limits` / `set_limits` / `discards_unrepresentable` / `set_discard_unrepresentable` / `write_options` | plumbing for the two above |

---

## Preserved source conventions

* **Per-array byte order and precision.** `precision="32"` selects 32-bit and
  *every other spelling* selects 64-bit; `endian="big"` selects big endian and
  *every other spelling* selects little (`MzDataHandler.cpp:503-527`). Both
  fallbacks are reproduced, each with a warning the source does not raise.
* **Base64 tolerance.** Whitespace is stripped from a payload before decoding,
  because "line breaks inside the base64 data are unfortunately no exception"
  (`MzDataHandler.cpp:492-494`). A payload shorter than four characters decodes
  to nothing, a length that is not a multiple of four is an error, and a
  trailing partial element is dropped (`Base64.h:315-344`). The upstream
  fixtures `MzDataFile_3_minimal.mzData` and `MzDataFile_4_64bit.mzData` depend
  on that last tolerance: their 36-character m/z payload is 25 bytes where the
  declared three 64-bit values need 24.
* **Declared lengths are advisory.** `length` is spent as a resource ceiling and
  then the payload wins, with a warning, exactly as `fillData_` does
  (`MzDataHandler.cpp:542-547`). `<spectrumList count>` is likewise not compared
  with the number of children — the source spends it only on `reserve` and a
  progress range.
* **Scan-window zero sentinel.** A window is kept only when at least one bound
  is nonzero, so a spectrum with `mzRangeStart="110"` and no `mzRangeStop`
  yields the window `(110, 0)` (`MzDataHandler.cpp:392-396`). That is what
  `MzDataFile_1.mzData` contains, and `ScanWindow::validate` rejects it, so the
  reader deliberately does not validate what it produces. `MSExperiment::validate`
  on the loaded fixture returns `Err`, asserted in
  `load_instrument_settings_and_acquisition`.
* **Filters.** `metadata_only` abandons the parse at `<spectrumList>`; an
  MS-level outside the selected set, a retention time outside the RT range and a
  precursor m/z outside the precursor range each drop the whole spectrum; the
  m/z and intensity ranges drop individual peaks together with the aligned
  annotation values. Range membership is `DRange<1>::encloses`, inclusive at the
  minimum and exclusive at the maximum (`DRange.h:152-161`).
* **Retention-time units.** `PSI:1000038` is minutes and is multiplied by 60;
  `PSI:1000039` is seconds (`MzDataHandler.cpp:1143`, `:1152`).
* **Polarity spellings.** `Positive`/`positive`/`+` and
  `Negative`/`negative`/`-` are all accepted, "be flexible here, actually only
  the first one is correct"; anything else warns and leaves the polarity unknown.
* **Vocabulary fallback.** `cvStringToEnum_` resolves an unknown term to index 0
  (`XMLHandler.cpp:131-145`). That is the enum's unknown value for every table
  except `cv_terms_[18]`, which has no leading empty entry, so an unrecognised
  activation method becomes `CID`. Reproduced, with the fallback named at each
  call site.
* **Scan-mode aliases.** `Zoom` and `EnhancedResolutionScan` both set the zoom
  flag with the `MassSpectrum` mode, and `ProductIonScan` overrides the
  `msLevel` attribute with 2.
* **Multiple precursor charges.** A second `PSI:1000041` resets the charge to
  zero and warns, rather than overwriting it.
* **XML strictness.** A duplicate attribute, an unbalanced tree, a mismatched
  end tag, a DTD and an undeclared entity are all refused, as Xerces refuses
  them for the source. CDATA sections and the five predefined entities and
  numeric character references resolve, because Xerces hands the source handler
  their content through `characters()` like any other text. quick-xml reports
  each `&…;` in character data as its own `Event::GeneralRef`, so the event
  loop resolves that variant explicitly and refuses `DocType` and an
  unexpanded `Empty` — there is no catch-all arm that could swallow a
  reference, the defect that turned `1&#46;5` into `15` in the Mascot XML
  reader (`src/format/mzml.rs:2465` is the established shape; commit
  `c218077` applied it to `src/format/imzml_handler.rs`).
  `tests/mzdata.rs::entity_references_in_character_data_never_vanish` splices a
  reference into a base64 payload — where a dropped one would change numbers
  with no error at all — and into a metadata value, and asserts the undeclared
  entity is refused.
* **Deliberately ignored terms.** `PSI:1000017` (`ScanFunction`),
  `PSI:1000020` (`TandemScanningMethod`), `PSI:1000035` (`PeakProcessing`),
  `PSI:1000043` (intensity unit) and `PSI:1000046` (energy unit, "we assume
  electronvolt") are read and dropped, and `<supSourceFile>`'s three children
  are read and dropped.
* **Writer layout.** Tab indentation, `\n` line endings, the
  `encoding="ISO-8859-1"` declaration, the element order, the placeholder
  `<contact>` and `<analyzerList count="1">` a document with none needs, the
  single document-wide `<dataProcessing>` taken from `spectra[0]`, and the
  zero-length placeholder spectrum an empty experiment produces. The whole
  empty-experiment document is asserted byte for byte in
  `stored_documents_are_wellformed`.
* **Two round-trip losses that are the source's own**, reproduced rather than
  repaired, and pinned by
  `tests/mzdata.rs::the_two_round_trip_losses_that_are_the_sources_own`:
  `writeCVS_` writes nothing for a numeric value that is exactly zero
  (`MzDataHandler.cpp:1442-1448`), so a retention time of exactly 0 s emits no
  `TimeInSeconds` and reads back as the `MSSpectrum` default of −1 (every other
  zero-valued field has zero as its default and so is unaffected); and the
  placeholder spectrum above means an empty experiment reloads with one empty
  spectrum in it. Both are also what `MzDataFile_test.cpp` gets away with: its
  three spectra sit at 60, 120 and 180 s and its stored experiments are never
  empty.
* **Native-ID renumbering.** `<spectrum id>` is an integer, so the writer
  reproduces the three-flag analysis of `MzDataHandler.cpp:764-800`: ids are
  taken from a `spectrum=`-prefixed number, else from a bare number, else the
  spectra are renumbered from 1 — silently when every native ID is empty, with
  a warning otherwise.
* **`<precursor spectrumRef>`** is written from the last id seen at each MS
  level, and `-1` when no parent level has been written. The reader has no
  handler for `msLevel` or `spectrumRef` on `<precursor>` and ignores both, so
  the value does not survive a round trip either way.

---

## Native differences

Each is documented at the Rust item as well.

1. **No panic on untrusted input.** Nine places in `MzDataHandler.cpp` call
   `.back()` on a container a hostile document can leave empty —
   `exp_->getContacts()` (`:117`, `:121`, `:125`), `exp_->getSourceFiles()`
   (`:150`, `:158`, `:166`), `meta_id_descs_` (`:265`, `:292`),
   `spec_.getPrecursors()` (`:284`, `:1186`…), the three instrument component
   vectors, and `exp_->getSpectra()` (`:1136`) — and dereference
   `data_processing_` (`:113`, `:129`, `:133`, `:246`, `:316`) without a null
   check. Each of those is an `Error::Parse` here.
2. **`fillData_`'s out-of-bounds reads are refused.** The source reads
   `precisions_[0]` and `precisions_[1]` *before* its
   `data_to_decode_.size() < 2` guard (`:522-533`), and after reporting an
   m/z-versus-intensity length disagreement through the non-fatal
   `error(LOAD, …)` (`:539`) it indexes the intensity array and every
   annotation array at every position below the m/z length (`:562-572`). This
   port returns `Ok` with no peaks for a spectrum with no binary array,
   `Error::Parse` for a mismatched pair, and `Error::Parse` for an annotation
   array whose length differs from the spectrum's.
3. **`peak_count_` does not leak between spectra.** It is a handler member
   assigned only from `<data>` inside `<mzArrayBinary>`, so upstream a spectrum
   without an m/z array inherits the previous one's count.
4. **m/z is written at 64-bit precision by default.** `writeBinary_` hardcodes
   `precision="32"` and stages every coordinate through a
   `std::vector<float>`, losing about nine significant digits of every mass.
   mzData declares the precision per array and both readers honour it, so the
   default here is `precision="64"`; `WriteOptions::mz_32_bit` (set by
   `WriteOptions::source`) selects the source behaviour.
5. **Text is XML-escaped and the declared encoding is kept truthful.**
   `writeTo` streams every `std::string` raw, so an `&` or a `<` in any
   metadata value produces a document that is not well-formed, and UTF-8
   metadata produces bytes `encoding="ISO-8859-1"` cannot describe. Here the
   five XML delimiters are escaped, characters in `U+0080..U+00FF` are written
   as single ISO-8859-1 bytes, and anything above `U+00FF` as a numeric
   character reference. Reading transcodes ISO-8859-1 exactly, accepts UTF-8
   and US-ASCII, and refuses any other declared encoding rather than
   misreading it.
6. **Unparsable numbers are refused.** `asDouble_` and `asInt_` log a
   non-fatal error and substitute `0` (`XMLHandler.h:264-317`). Substituting
   zero for a mass, an intensity or a charge is worse than refusing the
   document, so both are `Error::Parse` here. Nonfinite values are refused for
   the same reason, in an attribute and in a *decoded* array alike: a NaN or
   infinite m/z makes every later sort and binary search on that spectrum
   undefined, and an intensity outside the f32 range cannot be stored at all.
   `asDateTime_`'s failure remains non-fatal: an unparsable `completionTime`
   leaves `None` and raises a warning, as the source leaves its
   invalid-`DateTime` sentinel.
7. **Character data is accumulated per element.** The source applies each
   chunk as it arrives, so a value split across `characters()` calls keeps only
   the last chunk for every tag except `data`. Accumulating and applying at the
   closing tag also makes CDATA sections and character references work, which
   Xerces gives the source handler for free.
8. **`<supDataArrayBinary>` cannot shift its neighbours.** The source pushes
   the decode slot when `<arrayName>` opens (`:426-429`), so an annotation array
   without that child shifts every following payload by one slot and the
   trailing array reads past the end of `precisions_`. Here the `<data>`
   element creates its own slot, `<arrayName>` only names it, and each slot
   carries the index of the float data array it belongs to, so a
   `<supDataArrayBinary>` with no `<data>` child or with two is an error rather
   than a silent mispairing.
9. **The m/z and intensity arrays are identified by their parent element.**
   `fillData_` identifies them by *position* among the `<data>` elements —
   `precisions_[0]` and `precisions_[1]` (`:522-527`) — so a document that
   writes `<intenArrayBinary>` before `<mzArrayBinary>` has its two arrays
   silently swapped. This port classifies by parent tag and refuses a second
   array of either kind.
10. **A document with no `<dataProcessing>` leaves the spectra without one.**
   The source pushes `data_processing_` unconditionally at `<spectrum>`
   (`:345`), so such a document gives every spectrum a null
   `DataProcessingPtr`.
11. **The unknown-scan-mode fallback stays on its own spectrum.**
    `cvParam_` writes `MSNSPECTRUM` onto `exp_->getSpectra().back()` — the
    *previous* spectrum — and warns only on the MS1 branch (`:1136-1142`).
12. **The spectrum comment lives in metadata.** `spec_.setComment` (`:137`)
    has no `MSSpectrum` counterpart in this crate; the comment is on
    `SpectrumSettings`, which `MSSpectrum` does not embed, so the value is kept
    as a `comment` metadata entry. The source writer never emits the element,
    so it is read-only either way.
13. **No progress logging and no OpenMP.** The source handler takes a
    `ProgressLogger&` and keeps its scan counter in a **function-local
    `static UInt`** (`:436`), shared by every handler instance and thread, so
    concurrent or successive loads interleave the count. Neither the logger nor
    any threading is ported; the source parallelises nothing in this file, so
    there is no OpenMP gap to record beyond that.
14. **Nothing mzData cannot represent is discarded silently.** The default
    `WriteOptions` refuses; `WriteOptions::source` discards, warning where the
    source warns. The refusals are enumerated below.
15. **Every document this writer produces, this reader loads.**
    `<supDataArrayBinary>` is aligned with the peak array element for element,
    and mzData has no way to say otherwise. `writeTo` nevertheless writes an
    annotation array whose length differs from the spectrum's, reporting it
    through the non-fatal `error(LOAD, …)` (`MzDataHandler.cpp:1032-1037`) —
    and `fillData_` then reads `decoded_list_[2 + i][n]` past the end of that
    array for every peak it does not have (`:562-572`), which this reader
    refuses. A store followed by a load therefore failed on any experiment
    holding a misaligned array, including the empty placeholder array
    `DataArray` documents and `MSExperiment::validate` accepts. The array
    length is now checked in the preflight, before anything is written, and
    `WriteOptions::discard_unrepresentable` drops such an array — with the
    source's warning and the `<supDesc>` / `<supDataArrayBinary>` numbering of
    the surviving arrays kept in step — instead of writing a document that
    cannot be read back. A nonfinite m/z, intensity or annotation value is
    refused in *both* modes for the same reason: it is not a representability
    gap that could be discarded, and the reader refuses a decoded array that
    holds one (difference 6). `tests/mzdata.rs::annotation_array_length_must_match_the_peak_count`,
    `::misaligned_annotation_arrays_are_dropped_rather_than_written` and
    `::nonfinite_values_are_refused_by_the_writer` hold all three down.

### What the writer refuses by default

Experiment level: chromatograms; an SQL run id; more than one source file; a
source-file checksum, size, native-ID type or CV term; a contact email, URL,
address or metadata, and a first or last name containing a space (`<name>` is
one field the reader splits); a sample organism, comment or subsample; more than
one ion source or ion detector; an instrument ion-optics type or instrument
software; a nonzero component `order`; HPLC settings; an experiment date,
comment or fraction identifier; named instrument configurations;
experiment-level metadata; float data arrays while `write_supplemental_data` is
false; native IDs that are neither numbers nor `spectrum=` followed by a number.

Spectrum level: a spectrum name, per-spectrum source file, product, peptide
identification, spectrum-level metadata or drift time; integer or string data
arrays; more than one scan window; a scan mode outside
`Unknown`, `MassSpectrum`, `SelectedIonMonitoring`,
`SelectedReactionMonitoring`, `ConsecutiveReactionMonitoring`,
`ConstantNeutralGain`, `ConstantNeutralLoss` and `Precursor`; a zoom flag
outside `MassSpectrum`; an unknown spectrum type together with acquisitions, or
a spectrum type without them; acquisition-info metadata; an acquisition
identifier that is not a 32-bit integer; more than one activation method; a
precursor isolation window, charge-state list, CV term, drift time or spectrum
reference; an annotation array whose length differs from the spectrum's peak
count, the empty placeholder array included; data-array processing history;
more than one data-processing record per spectrum, spectra whose records differ,
software metadata, and a processing action outside `Deisotoping`,
`ChargeDeconvolution` and `PeakPicking`.

Refused in **both** modes, `WriteOptions::source` included: a nonfinite m/z,
intensity or annotation value. That is not something mzData cannot represent —
it is a value the reader refuses on the way back in, so writing it would
produce a document this crate cannot load (difference 15).

`ScanMode::Ms1Spectrum` and `MsnSpectrum` are refused because `writeTo` maps
both to the same `MassScan` value as `MassSpectrum`. `Absorption`,
`EnhancedMultiplyCharged` and `TimeDelayedFragmentation` are refused because the
source writer emits `PhotodiodeArrayDetector`, `EnhancedMultiplyChargedScan` and
`TimeDelayedFragmentationScan` — three values its own `cvParam_` does not
recognise, so each reads back as `MassSpectrum`.

---

## Checked boundaries and evidence

### Ceilings

The source has none: `<spectrumList count>` goes straight into
`MSExperiment::reserve` (`MzDataHandler.cpp:352-353`) after
`attributeAsInt_` returns a signed `Int` that is assigned to a `UInt`, so a
huge or negative count asks for roughly a terabyte before any data is read.
Every ceiling below is checked before the allocation it bounds, so a refused
document never commits the memory it asked for, and
`tests/mzdata.rs::read_ceilings_are_enforced` exercises each one.

| `ReadLimits` field | Default | Bounds |
|---|---|---|
| `max_xml_bytes` | 512 MiB | input bytes, counted before any transcoding |
| `max_xml_depth` | 64 | open elements |
| `max_spectra` | 1 000 000 | retained spectra, and the `<spectrumList count>` attribute |
| `max_arrays_per_spectrum` | 4096 | `<data>`, `<supDataArrayBinary>`, `<supDesc>` and `<acquisition>` elements per spectrum |
| `max_array_bytes` | 64 MiB | decoded bytes of one array, checked against the encoded character count before the decoder allocates |
| `max_array_elements` | 20 000 000 | elements of one array, and the declared `length` attribute |
| `max_total_peaks` | 50 000 000 | retained peaks across the document |
| `max_text_bytes` | 512 MiB | cumulative character data, charged whether or not the spectrum is kept |
| `max_metadata_entries` | 1 000 000 | metadata entries created across all records |
| `max_warnings` | 256 | retained warnings; `LoadReport::warning_count` keeps counting past it |

Atomicity: the parsed experiment is committed only when the whole document has
been read, each spectrum is assembled in a temporary and pushed at
`</spectrum>`, peak selection runs before anything is appended, and
`store_report_with_options` serialises into memory and refuses before any file
is created — `a_refused_store_leaves_the_destination_alone` and
`store_refuses_to_drop_the_software_comment` both assert that.

No `unsafe`. No index, slice or arithmetic on file-derived data is unchecked:
lengths go through `checked_sub`/`checked_add`, allocations through
`try_reserve_exact`, byte groups through `chunks_exact`, and the one prefix
strip is `str::strip_prefix`. No string this module did not construct is
byte-sliced; the `completionTime` truncation to 19 characters uses
`chars().take(19)`, not `&text[..19]`, which is exactly the shape that aborted
a process in this project's earlier audit.

### Class-test sections

All fifteen `START_SECTION`s of `MzDataFile_test.cpp` (435 assertion macros)
are accounted for. Thirteen are ported; two — both at or below the
five-macro threshold — are mapped with the Rust function and the concrete
asserted value they reproduce.

| # | Section (line range, macros) | Status | Rust | Asserted value reproduced |
|---|---|---|---|---|
| 1 | `MzDataFile()` (36-41, 1) | ported | `default_construction_and_drop` | the adapter exists and pins version `1.05` |
| 2 | `~MzDataFile()` (43-47, 0) | ported | `default_construction_and_drop` | the adapter is dropped; the C++ section only `delete`s |
| 3 | `const PeakFileOptions& getOptions() const` (49-57, 2) | ported | `option_accessors` | `hasMSLevels() == false` on a fresh adapter and on a const copy |
| 4 | `setOptions(const PeakFileOptions&)` (59-71, 2) | ported | `option_accessors` | `hasMSLevels()` false, then true after `setOptions` |
| 5 | `PeakFileOptions& getOptions()` (73-79, 1) | ported | `option_accessors` | `hasMSLevels() == true` after `addMSLevel(1)` in place |
| 6 | `load(...)` (81-497, 299) | ported | `load_document_identity_and_scan_axis`, `load_annotation_array_descriptions`, `load_precursors`, `load_instrument_settings_and_acquisition`, `load_peaks_and_annotation_values`, `load_experiment_metadata`, `load_special_cases` | `e.size() == 3`; RT 60/120/180; `spectrum=10`/`11`/`12`; `lsid`; `MS-Sample`/`0-815`/gas; resolutions 22.33 and 12.3; 997530 peaks; the 64-bit and big-endian arrays |
| 7 | `[EXTRA] load with metadata - only flag` (499-524, 10) | ported | `load_metadata_only` | `e.size() == 0` with `MzDataFile_test_1.raw` and `MS-Instrument` still present |
| 8 | `[EXTRA] load with selected MS levels` (526-555, 14) | ported | `load_selected_ms_levels` | `e.size() == 2`, native IDs `spectrum=10` and `spectrum=12`, 1 and 5 peaks |
| 9 | `[EXTRA] load with RT range` (557-577, 5) | ported | `load_rt_range` | `e.size() == 2`, RT 120 and 180 |
| 10 | `[EXTRA] load with MZ range` (579-614, 14) | ported | `load_mz_range` | peak counts 1/2/2 with coordinates 120; 120,130; 120,130 |
| 11 | `[EXTRA] load with intensity range` (616-648, 12) | ported | `load_intensity_range` | peak counts 0/1/3 with intensities 200; 200,300,200 |
| 12 | `store(...)` (650-667, 3) | ported | `store_round_trip_equals_the_loaded_experiment` | `e2.getIdentifier() == "lsid"` and `e1 == e2` after restoring the software comment |
| 13 | `[EXTRA] storing / loading of meta data arrays` (669-828, 69) | ported | `store_and_load_annotation_arrays` | array counts 1/0/2, names `MDA1`/`MDA2`, values 1.1…1.5 and −2.1…−2.5, and 1.3/1.4/1.5 after the [2.5, 7.0) filter |
| 14 | `[EXTRA] static bool isValid(...)` (830-847, 2) | **mapped** | `stored_documents_are_wellformed` | `isValid(...) == true` for both stored documents, reproduced as "both are well-formed, both reload, and the empty-experiment document matches the exact byte layout". The XSD itself is not shipped; see the gap below |
| 15 | `bool isSemanticallyValid(...)` (849-854, 1) | **mapped** | `semantic_validation_is_unsupported` | `NOT_TESTABLE` — the section asserts nothing, and the method reports that the mapping file and OBO are unavailable |

**Sections ported: 13. Mapped with a cited asserted value: 2. Unaccounted: 0.**

### Evidence tier

Tier 3, source review with transcribed literals, plus one independent
cross-check. No C++ was built or executed, no C++ output was retained, and no
retained mzData output exists anywhere in the pinned tree (`MzDataFile_1.mzData`
is hand-written: it carries `ScanFunction` and `TandemScanningMethod` cvParams
the writer never emits), so this is **not** a tier 1 differential.

The independent check is
`tests/mzdata.rs::upstream_payload_holds_a_partial_trailing_element`: the m/z
payload of `MzDataFile_3_minimal.mzData` is decoded outside the reader and shown
to hold 25 bytes where the declared three 64-bit values need 24. That the
upstream fixtures depend on `Base64::decodeUncompressed_`'s truncation is
established by the bytes, not by anything the class test asserts.

Tier 4, independently derived: the big-endian, ISO-8859-1 and UTF-8 fixtures;
every resource ceiling; the hostile `length` and `count` attributes; the writer
refusals; the round-trip closure of difference 15
(`annotation_array_length_must_match_the_peak_count`,
`misaligned_annotation_arrays_are_dropped_rather_than_written`,
`nonfinite_values_are_refused_by_the_writer`,
`big_endian_input_survives_a_round_trip`); the character-reference regression
`entity_references_in_character_data_never_vanish`; and the byte-exact
empty-experiment document, which is transcribed from
`MzDataHandler.cpp:583-1072` rather than captured from a run.

---

## Gaps

* **XSD validation** (`XMLFile::isValid`, class-test section 14) is not ported:
  `mzData_1_05.xsd` does not ship with this crate and
  `crate::format::mzml_schema` is mzML-specific. A caller with the schema can
  validate with any XML tool; the port guarantees well-formedness only.
* **Semantic validation** (`MzDataFile::isSemanticallyValid`) is not ported:
  `mzdata-mapping.xml` and `psi-mzdata.obo` do not ship either, and the upstream
  section is `NOT_TESTABLE` by its own admission.
* **`ProgressLogger`** is a base class of `MzDataFile` and a constructor
  parameter of the handler; it is not threaded through this module.
  `LoadReport` and `StoreReport` carry the counters instead.
* **`MzDataFile_2_long.mzData`** (10.6 MB, 997530 peaks) is not copied into
  `tests/data`: it would nearly double the 15 MB of fixture data. Its sha256 is
  recorded in the manifest and `load_special_cases` writes and reads an
  equivalent 997530-peak document instead. What upstream calls "CDATA
  splitting" is character data split across parser callbacks, which
  `split_base64_payload_is_concatenated` reaches directly with an XML comment
  inside a payload.
* **`FileHandler` is not wired to mzData.** `FileType::MzData` already exists
  and `type_by_content` already recognises `<mzData`
  (`src/format/file_handler.rs:355`), but `FileHandler`'s load and store
  dispatch does not call this module. `src/format/file_handler.rs` is outside
  this package's scope.
* **`FORMAT/VALIDATORS/MzDataValidator.h`** is a separate header with its own
  class test and is unported.
* **The spectrum comment** has no `MSSpectrum` field in this crate and is kept
  as metadata; see native difference 11.
* **`docs/doc-coverage.json` was not rewritten.** `src/format/mzdata.rs`
  measures 100.0% (44/44) and `src/format/mod.rs` improves from 4/28 to 5/29,
  so `check_doc_coverage.py` passes, but recording the new floor with `--write`
  touches a file outside this package's scope.
* **`docs/core-sdk-coverage.json` and `docs/CORE_SDK_COMPLETION.md` are stale**
  until the integrator runs `python3 tools/core_sdk_coverage.py --write`: the
  generator now detects `src/format/mzdata.rs` as a candidate for
  `MzDataFile.h`, which moves that header from `unmapped` to
  `evidence_requires_review` before the ledger entries are added. Both files are
  outside this package's scope.
* **CI wiring** was not added for the same reason: append `--test mzdata` to
  `.github/workflows/rust.yml:87`, the second `--no-default-features --features
  mzml` line of the `minimum-rust` job (the one that already carries the imzML
  tests). Verified: all 63 tests pass under
  `cargo nextest run --locked --no-default-features --features mzml --test mzdata`
  and `cargo +1.85.0 check --locked --all-features --all-targets` is clean.

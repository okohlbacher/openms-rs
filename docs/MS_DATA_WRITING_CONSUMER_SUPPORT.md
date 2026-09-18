# Streaming mzML writing consumer

Source pin: `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. This increment ports
`FORMAT/DATAACCESS/MSDataWritingConsumer.h` (251 lines) and its `.cpp`
(175 lines) to `src/format/ms_data_writing_consumer.rs`, gated on the existing
`mzml` feature. Ten TOPP tools depend on the header — `CometAdapter`,
`FileConverter`, `NoiseFilterGaussian`, `NoiseFilterSGolay`,
`OpenSwathMzMLFileCacher`, `OpenSwathWorkflow`, `PeakPickerHiRes`,
`PeakPickerIM`, `SageAdapter` and `TICCalculator` — and every one of them needs
the same thing: a sink that turns a stream of records into an mzML file without
holding the experiment in memory, with a hook to transform each record on the
way past.

The source class is an abstract `Internal::MzMLHandler` subclass that also
implements `Interfaces::IMSDataConsumer`. It writes the mzML header when the
first record arrives, appends each record immediately, and closes the document
in its destructor. Its two private pure virtuals are the template-method hooks.

## API mapping

Every public and protected member of the three classes in the header.

### `MSDataWritingConsumer`

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `typedef PeakMap MapType` | `crate::kernel::MSExperiment` | Used internally to build the one-record document. |
| `typedef MapType::SpectrumType SpectrumType` | `crate::kernel::MSSpectrum` | |
| `typedef MapType::ChromatogramType ChromatogramType` | `crate::kernel::MSChromatogram` | |
| `explicit MSDataWritingConsumer(const std::string& filename)` | `MSDataWritingConsumer::create(path, processor)` | Also `PlainMSDataWritingConsumer::create_plain(path)`. `Error::Io` when the file cannot be created; the source does not check `is_open`. |
| *(no source counterpart)* | `MSDataWritingConsumer::new(writer, processor)`, `PlainMSDataWritingConsumer::plain(writer)` | Any `W: Write`. The source is bound to an owned `std::ofstream`. |
| `~MSDataWritingConsumer()` (and the `virtual ~` the test names) | `MSDataWritingConsumer::finish() -> Result<W>` | The destructor calls `doCleanup_`. A Rust drop cannot report an I/O failure, so this is explicit. |
| `void setExperimentalSettings(const ExperimentalSettings& exp) override` | `MSDataWritingConsumer::set_experimental_settings(&ExperimentalSettings)` | `Result`: refused once the header is written, where the source's late value is silently unused. Also reachable through `MSDataConsumer`. |
| `void setExpectedSize(Size, Size) override` | `MSDataWritingConsumer::set_expected_size(usize, usize)` | `Result`: refused past a written list tag, or over `WritingLimits`. Also reachable through `MSDataConsumer`. |
| `void consumeSpectrum(SpectrumType& s) override` | `MSDataWritingConsumer::consume_spectrum(&mut MSSpectrum)` | `MSDataConsumer::consume_spectrum` wraps it as `Ok(ControlFlow::Continue(()))`. |
| `void consumeChromatogram(ChromatogramType& c) override` | `MSDataWritingConsumer::consume_chromatogram(&mut MSChromatogram)` | Likewise. |
| `virtual void addDataProcessing(DataProcessing d)` | `MSDataWritingConsumer::add_data_processing(DataProcessing)` | `Result`: refused once the header is written, because the header declares the `dataProcessingList`. |
| `virtual Size getNrSpectraWritten()` | `MSDataWritingConsumer::spectra_written()` | |
| `virtual Size getNrChromatogramsWritten()` | `MSDataWritingConsumer::chromatograms_written()` | |
| `private virtual void processSpectrum_(SpectrumType&) = 0` | `MSDataWritingProcessor::process_spectrum` | Returns `Result<()>`; the source's `void` hook can only throw. |
| `private virtual void processChromatogram_(ChromatogramType&) = 0` | `MSDataWritingProcessor::process_chromatogram` | Likewise. |
| `private virtual void doCleanup_()` | the body of `finish()` | Not separately overridable; `NoopMSDataWritingConsumer` is the only source class that overrides it, and that class is a distinct type here. |
| `protected std::ofstream ofs_` | the `W` type parameter | |
| `protected bool started_writing_` | `MSDataWritingConsumer::started_writing()` | Read-only accessor; protected upstream. |
| `protected bool writing_spectra_` | `MSDataWritingConsumer::writing_spectra()` | Read-only. |
| `protected bool writing_chromatograms_` | `MSDataWritingConsumer::writing_chromatograms()` | Read-only. |
| `protected Size spectra_written_` | `MSDataWritingConsumer::spectra_written()` | The public getter and the member are one thing here. |
| `protected Size chromatograms_written_` | `MSDataWritingConsumer::chromatograms_written()` | |
| `protected Size spectra_expected_` | `MSDataWritingConsumer::expected_size().0` | |
| `protected Size chromatograms_expected_` | `MSDataWritingConsumer::expected_size().1` | |
| `protected bool add_dataprocessing_` | `MSDataWritingConsumer::additional_data_processing().is_some()` | The flag and the pointer are one `Option` here. |
| `protected DataProcessingPtr additional_dataprocessing_` | `MSDataWritingConsumer::additional_data_processing() -> Option<&Arc<DataProcessing>>` | |
| `protected Internal::MzMLValidator* validator_` | not ported: the crate's semantic validation is a separate operation behind the `mzml-validation` feature and is not run per record | |
| `protected ExperimentalSettings settings_` | `MSDataWritingConsumer::settings()` | Read-only. |
| `protected std::vector<std::vector<ConstDataProcessingPtr>> dps_` | not ported as a field: the equivalent state is the rendered `dataProcessingList` text and the set of identifiers it declares, both captured from the first record's header | See *Native differences*. |
| inherited `MzMLHandler::writeHeader_` / `writeSpectrum_` / `writeChromatogram_` | `crate::format::mzml::write_with_options`, driven once per record | The port composes rather than inherits. |
| inherited `MzMLHandlerHelper::writeFooter_` | the document-closing tags and the index in `finish()` | Indexed, as the inherited `write_index_` is. See *Indexed output*. |
| inherited `ProgressLogger` base | not ported | |

### `PlainMSDataWritingConsumer`

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `class PlainMSDataWritingConsumer : public MSDataWritingConsumer` | `type PlainMSDataWritingConsumer<W> = MSDataWritingConsumer<W, PlainProcessor>` | |
| `explicit PlainMSDataWritingConsumer(std::string filename)` | `PlainMSDataWritingConsumer::create_plain(path)` | |
| `void processSpectrum_(...) override {}` | `PlainProcessor::process_spectrum` | |
| `void processChromatogram_(...) override {}` | `PlainProcessor::process_chromatogram` | |

### `NoopMSDataWritingConsumer`

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `class NoopMSDataWritingConsumer : public MSDataWritingConsumer` | `NoopMSDataWritingConsumer` | A separate type, not a specialisation: it shares no state with the writing consumer. |
| `explicit NoopMSDataWritingConsumer(std::string filename)` | `NoopMSDataWritingConsumer::new()` | **Takes no path.** The source still runs the base constructor, which creates and truncates the named file. |
| `void setExperimentalSettings(...) override {}` | `MSDataConsumer::set_experimental_settings`, `Ok(())` | |
| `void consumeSpectrum(...) override {}` | `MSDataConsumer::consume_spectrum`, counts and continues | |
| `void consumeChromatogram(...) override {}` | `MSDataConsumer::consume_chromatogram`, counts and continues | |
| `private void doCleanup_() override {}` | nothing to clean up | |
| `private void processSpectrum_(...) override {}` | not applicable: no hook | |
| `private void processChromatogram_(...) override {}` | not applicable | |
| *(inherited `getNrSpectraWritten`, always 0)* | `NoopMSDataWritingConsumer::spectra_written()`, counts what arrived | Native: the source's overrides discard the records without counting them. |
| *(inherited `getNrChromatogramsWritten`, always 0)* | `NoopMSDataWritingConsumer::chromatograms_written()` | |
| *(inherited `setExpectedSize`, not overridden)* | `MSDataConsumer::set_expected_size`, `Ok(())` | |

Native additions with no source counterpart: `WritingLimits` and its three
constants, `CountPolicy`, `with_limits`, `with_write_options`,
`with_count_policy`, `limits`, `write_options`, `count_policy`, `processor`.

## Preserved source conventions

- **Nothing is written until the first record arrives**, and that record's
  arrival writes the header (`MSDataWritingConsumer.h:49`). A consumer that
  receives no record produces an empty output, because `doCleanup_` writes the
  footer only under `started_writing_` (`.cpp:167`).
- **The header describes the settings plus the first record only.** The source
  builds a dummy one-record map for `writeHeader_` (`.cpp:76`); this port
  renders a one-record document the same way.
- **The list tags are written by hand, once, carrying the expected counts**
  (`.cpp:89` and `.cpp:134`), with `defaultDataProcessingRef` pointing at the
  first data-processing entry.
- **Spectra may not follow chromatograms.** `Error::InvalidValue` carries the
  source's own message, "cannot write spectra after writing chromatograms"
  (`.cpp:57`), and the reason is that two `spectrumList` elements cannot appear
  in one mzML file.
- **A chromatogram closes an open `spectrumList` first** (`.cpp:102`), so the
  lists never interleave and `doCleanup_` only ever has one open.
- **The record is copied before processing** (`.cpp:62`), so neither the hook
  nor the added data-processing entry is visible to the caller. A test asserts
  the caller's record is unchanged after both `consume_spectrum` and
  `consume_chromatogram`.
- **`addDataProcessing` replaces rather than appends** (`.cpp:143`): two calls
  leave one entry.
- **The running record number is the `index` attribute** and, in this port, the
  fallback identifier — `writeSpectrum_` takes `spectra_written_++`
  (`.cpp:95`).
- **`doCleanup_`'s order**: close the open list, then the footer only if writing
  started, then close the stream (`.cpp:151-173`).
- **The expected size is not enforced** (`MSDataWritingConsumer.h:56`). It is
  still not *enforced* here — the announced value is what goes into the file —
  but the mismatch is reported; see below.

## Native differences

- **The port composes instead of inheriting.** The source *is* an
  `MzMLHandler` and calls its protected `writeHeader_`, `writeSpectrum_`,
  `writeChromatogram_` and `writeFooter_`. The equivalent internals of
  `src/format/mzml*` are private to that module, and this package does not own
  it, so the consumer drives the public whole-document writer once per record
  and splices the record element out from between the list tags. Three tag
  spellings and the `index="0"` marker are held as module constants and every
  one is checked: a change to the writer's layout makes the consumer return
  `Error::Unsupported` rather than emit wrong XML. Two tests pin the result —
  a streamed three-spectrum document is byte-identical to
  `mzml::write` over the same experiment, and every streamed document is read
  back with `mzml::read` and compared record by record.
- **Cost.** Each record re-renders the header text and discards it, so the work
  per record is the header plus the record, where the source pays the record
  alone. Memory stays at one record, which is the property the class exists
  for. Rendering later records without the settings is *not* an option: a
  record's `sourceFileRef` indexes a `sourceFileList` that begins with the
  settings' own source files, so dropping them would renumber the reference and
  point it at the run's source file instead of the record's. A test with two
  run-level source files pins that the record reference stays `sf_…2`.
- **Dangling references are refused, not written.** Only the first record
  contributes to the header, so a later record needing a different
  `sourceFileList` or `dataProcessingList` — its own, or one on an auxiliary
  array — has nothing correct to point at. This port compares the rendered
  declaration lists and the set of declared identifiers and returns
  `Error::Unsupported`. The source emits the dangling reference:
  `MzMLHandler.cpp:5254` builds `sourceFileRef="sf_sp_<n>"` from the running
  index, and `MzMLHandler.cpp:5258-5272` falls back to
  `dataProcessingRef="dp_sp_<n>"` when no entry in its one-element `dps` matches
  — the gap the source's own `// TODO ... assert this here` at `.cpp:93` marks.
  The consequence for a caller is real: a stream with per-record source files or
  histories cannot be written by this port at all, and must go through the
  whole-document writer.
- **An empty `native_id` is filled in** with `index=N` or `chromatogram=N`,
  which is what the whole-document writer would have produced at that position.
  The source writes the empty identifier through, giving a file whose records
  all carry `id=""`.
- **Duplicate identifiers are refused.** The whole-document writer rejects them;
  the streaming path must not be the weaker one. Retained identifier text is
  charged against `WritingLimits::max_native_id_bytes`.
- **`CountPolicy`.** The source's own `@note` says a wrong expected size leaves
  an inconsistent mzML. `CountPolicy::Checked`, the default, closes the document
  and *then* returns `Error::InvalidValue` naming both pairs, so the caller
  learns what it shipped; `CountPolicy::SourceInconsistent` accepts it silently,
  as upstream. This follows the crate's established pattern of refusing a lossy
  source behaviour by default and offering it explicitly.
- **`finish()` instead of a destructor.** Dropping a consumer without calling it
  leaves the `run` and `mzML` elements unclosed. That is deliberate: a
  best-effort `Drop` would have to swallow the I/O failure on a file the caller
  believes is complete.
- **`NoopMSDataWritingConsumer` takes no path**, so asking for a consumer that
  does nothing cannot truncate an existing output — the source's does, through
  its base constructor.
- **`MzMLValidator` and `ProgressLogger` are not threaded through.**
- **Serial.** `MzMLHandler` carries `#pragma omp` in its binary encoding; this
  port introduces no threads, so a large record encodes on one core.

## Indexed output

The source consumer is an `MzMLHandler`, so it carries that handler's
`PeakFileOptions`, whose `write_index_` defaults to true, and no TOPP consumer
reaches those options to change it. `doCleanup_` therefore hands over to
`MzMLHandlerHelper::writeFooter_`, which emits `indexList`, `indexListOffset`,
`fileChecksum` and the closing `indexedmzML` from the offsets
`writeSpectrum_`/`writeChromatogram_` recorded; `writeHeader_` opened the
`indexedmzML` element to match. The retained upstream output
`PeakPickerHiRes_output_lowMem.mzML`, which `TOPP_PeakPickerHiRes_3` compares
against, is an `indexedmzML` accordingly.

This port does the same, and reuses the whole-document writer's own output
adapter (`mzml::IndexedOutput`, `src/format/mzml_write_options.rs`) rather than
restating the layout: the adapter writes the `indexedmzML` opening after the XML
declaration, counts bytes, and hashes them as they go. This consumer keeps its
own identifier-and-offset table, because it learns its records one at a time and
cannot reserve the adapter's table up front, and passes it to
`Output::footer_ids` at the end.

Two consequences are asserted by the tests:

* **A streamed document is byte-identical to the document the whole-document
  indexed writer would have produced** from the same records, when the announced
  list counts are the real ones. The record blocks already came from that writer;
  the index entries, `indexListOffset` and the SHA-1 `fileChecksum` fall out of
  the same byte positions. The file a caller gets from streaming is the file it
  would have got from holding the experiment.
* **A consumer that received no record still writes nothing at all**, index
  included, because `doCleanup_` writes a footer only when `started_writing_` is
  set. The source's `writeFooter_` would emit a dummy `-1` index entry, but it is
  never reached on that path.

`fileChecksum` is the real SHA-1 over every byte through the opening
`<fileChecksum>` tag, as the indexed mzML schema specifies and as this crate's
whole-document writer already did; the source writes the constant `0` there
(CPP-049).

## Checked boundaries and evidence

`WritingLimits` bounds the spectra and chromatograms one consumer may write
(one million each), the intermediate document rendered for one record (256 MiB,
enforced inside a bounded `Write` sink so the ceiling stops the render rather
than measuring it afterwards, and surfacing as `Error::Io`) and the retained
identifier text (64 MiB). Every check runs before a byte of the record reaches
the writer, so a rejected record leaves the file exactly as it was; bytes
already written for earlier records are not rolled back, which no streaming
writer can do. Counters and flags are not advanced on a rejected record — the
tests assert that for the ordering rule, both count ceilings, the byte ceiling,
the duplicate identifier, the undeclared reference and a refusing processor.

The splice never byte-slices a string this module did not construct: the
document comes from the crate's own writer, and every offset used is either a
`str::find` result on an ASCII marker or such a result plus that marker's
length, so no slice can fall inside a multi-byte character in a caller-supplied
identifier. `get` is used throughout rather than indexing, so a layout
mismatch is an error and not a panic. A test streams
`scan=日本語` and reads it back.

Evidence is tier 3 for the structure and tier 4 for everything else. The
upstream class test contributes no expected value at all: ten of its eleven
sections contain only `// TODO`, the eleventh calls
`new MSDataWritingConsumer()` on an abstract class that has no default
constructor, and `executables.cmake:239` comments the whole test out so it is
never built. All eleven sections are still accounted for in
`tests/ms_data_writing_consumer.rs`, with expectations read from the header and
`.cpp` and each pinned to a line-level anchor in
[the provenance record](../tests/data/ms_data_writing_consumer_provenance.json).
The strongest checks depend on no C++ at all: byte-identity with
`mzml::write` and a full read-back comparison. No full format implementation was built or executed and
no format output was retained, so this is not a tier-1 differential.

Two candidate defects are specific to this class: dangling `sf_sp_`/`dp_sp_`
references in streamed files and a disabled class test whose constructor call
no longer matches the class. Both are source-review findings; no full consumer
C++ execution is claimed. The separate `SVOutStream` manipulator-state finding
is recorded in [its support document](SV_OUT_STREAM_SUPPORT.md).

The earlier additional claim that `std::string + Size` narrows numeric IDs to
characters is withdrawn. The pinned `StringUtils.h` supplies numeric string
operators, so IDs such as `sf_sp_0` and `dp_sp_0` contain decimal numbers. An
isolated overload-resolution probe is recorded in the manifest; it does not
validate streaming-reference correctness or execute the full SDK.

# MS data consumers of `FORMAT/DATAACCESS`

Source pin: core `bc9cc12`. This increment ports seven of the eight consumer
headers of `FORMAT/DATAACCESS` onto the crate's consumer trait
[`MSDataConsumer`](../src/interfaces.rs) (the source's
`Interfaces::IMSDataConsumer`):

| Header | Rust | Feature | Status |
|---|---|---|---|
| `NoopMSDataConsumer.h` | `src/format/noop_ms_data_consumer.rs` | none | complete |
| `MSDataStoringConsumer.h` | `src/format/ms_data_storing_consumer.rs` | none | complete |
| `MSDataChainingConsumer.h` | `src/format/ms_data_chaining_consumer.rs` | none | complete |
| `MSDataTransformingConsumer.h` | `src/format/ms_data_transforming_consumer.rs` | none | complete |
| `MSDataAggregatingConsumer.h` | `src/format/ms_data_aggregating_consumer.rs` | none | partial: the merge is caller-supplied until `ANALYSIS/OPENSWATH/SpectrumAddition.h` is ported |
| `MSDataSqlConsumer.h` | `src/format/ms_data_sql_consumer.rs` | `sqmass` | complete |
| `MSDataCachedConsumer.h` | `src/format/ms_data_cached_consumer.rs` | none | complete (wave-4 integration, on lane F4's `CachedMzMLHandler`) |
| `MSDataWritingConsumer.h` | `src/format/ms_data_writing_consumer.rs` | `mzml` | partial, see [its own document](MS_DATA_WRITING_CONSUMER_SUPPORT.md) |

## The trait

`IMSDataConsumer` has four `void` members. `MSDataConsumer` keeps the four and
makes each fallible (`Result`), and `consume_*` additionally answers
`ControlFlow<()>`: `Break` asks the producer to stop before the record is
retained, a soft stop the source interface does not have. Every consumer here
answers `Continue` unless it forwards another consumer's `Break`.

## API mapping

### `NoopMSDataConsumer`

| C++ | Rust | Notes |
|---|---|---|
| `NoopMSDataConsumer()` | `NoopMSDataConsumer::new()`, `Default` | Zero-sized. |
| `setExperimentalSettings`, `setExpectedSize`, `consumeSpectrum`, `consumeChromatogram` (all `{}`) | the four trait methods, each `Ok` / `Ok(Continue)` | |
| file-scope `default_noop_consumer` in the `.cpp` | none | Declared in no header, so no caller can name it. |

### `MSDataStoringConsumer`

| C++ | Rust | Notes |
|---|---|---|
| `MSDataStoringConsumer()` | `new()`, `Default` | |
| `setExperimentalSettings(settings)` | trait method | `exp_ = settings` assigns only the `ExperimentalSettings` base: records are kept. |
| `setExpectedSize(s, c)` | trait method | A reservation hint capped at `MAX_RESERVED_RECORDS` (65,536) per list; a failed reservation is ignored. The source reserves whatever the producer passes. |
| `consumeSpectrum` / `consumeChromatogram` | trait methods | Store a copy; the caller's record is unchanged, as `addSpectrum(s)` copies. |
| `const PeakMap& getData() const` | `data() -> &MSExperiment` | Plus `into_data()` to take it without a copy. |
| private `PeakMap exp_` | private `experiment` | |

### `MSDataChainingConsumer`

| C++ | Rust | Notes |
|---|---|---|
| `MSDataChainingConsumer()` | `new()`, `Default` | |
| `MSDataChainingConsumer(std::vector<IMSDataConsumer*>)` | `from_consumers(Vec<&'a mut dyn MSDataConsumer>)` | Borrowed, not owned, as the source "does not transfer ownership"; the lifetime is the source's "keep every consumer alive". |
| `~MSDataChainingConsumer()` | drop | Deletes nothing, in both. |
| `appendConsumer(IMSDataConsumer*)` | `append_consumer(&'a mut dyn MSDataConsumer)` | |
| the four overrides | trait methods | In chain order. The first error stops the loop (an exception leaves the source's loop). A `Break` stops the chain for that record and is returned: native, see above. |
| private `consumers_` | private; `len()`, `is_empty()` | |

### `MSDataTransformingConsumer`

| C++ | Rust | Notes |
|---|---|---|
| `MSDataTransformingConsumer()` | `new()`, `Default` | The three functions start unset (`nullptr`). |
| `setSpectraProcessingFunc(std::function<void(SpectrumType&)>)` | `set_spectra_processing_func(impl FnMut(&mut MSSpectrum) -> Result<()> + 'a)` | Fallible, so a function can fail the consumption instead of throwing. Passing `nullptr` is `clear_spectra_processing_func()`. |
| `setChromatogramProcessingFunc(...)` | `set_chromatogram_processing_func`, `clear_chromatogram_processing_func` | |
| `setExperimentalSettingsFunc(...)` | `set_experimental_settings_func`, `clear_experimental_settings_func` | |
| `setExpectedSize` | trait method | Ignored. |
| `setExperimentalSettings` | trait method | Calls the settings function when set. The header's comment says "ignored"; the implementation calls `lambda_exp_settings_`, and this follows the implementation. |
| `consumeSpectrum` / `consumeChromatogram` | trait methods | Apply the function in place when set. |
| protected `lambda_spec_`, `lambda_chrom_`, `lambda_exp_settings_` | private; `has_*_func()` | The lifetime `'a` of the stored closures is the source's "captured variables must still be in scope". |

### `MSDataAggregatingConsumer`

| C++ | Rust | Notes |
|---|---|---|
| `MSDataAggregatingConsumer(IMSDataConsumer* next_consumer)` | `new(&'a mut dyn MSDataConsumer, add_up: F)` | `F: FnMut(&[MSSpectrum]) -> Result<MSSpectrum>` is the merge the source hard-codes as `SpectrumAddition::addUpSpectra(s_list, -1, true)` (`MSDataAggregatingConsumer.cpp:23`, `:44`). **This is the header's one open gap**: `SpectrumAddition` is an `ANALYSIS/OPENSWATH` header, not ported, and `format` may not depend on `analysis` (`tools/check_module_cycles.py`). When it lands, a constructor that passes it is one line in the analysis domain. |
| `~MSDataAggregatingConsumer()` | `finish() -> Result<ControlFlow<()>>`, and the same flush on drop | A drop cannot report an error, so it logs a warning; `finish` reports it. |
| `setExpectedSize`, `setExperimentalSettings` (`{}`) | trait methods | Not forwarded to the next consumer, as in the source. |
| `consumeSpectrum` | trait method | Groups consecutive spectra with `|RT - previous_rt| < 1e-5` (`RT_TOLERANCE`), to the *previous* spectrum, so a drift chains a group; merges each group and copies the first spectrum's metadata (`copySpectrumMeta(s_list[0], tmps, false)`). A failed merge or next consumer leaves the group, the incoming spectrum and the RT state unchanged. |
| `consumeChromatogram` | trait method | Passed on unchanged and immediately. |
| private `next_consumer_`, `previous_rt_`, `rt_initialized_`, `s_list` | private; `pending()` | `s_tmp` is an unused scratch member in the source. |

### `MSDataSqlConsumer`

| C++ | Rust | Notes |
|---|---|---|
| `MSDataSqlConsumer(sql_filename, run_id, buffer_size = 500, full_meta = true, lossy_compression = false, linear_mass_acc = 1e-4)` | `new(path, run_id)` with the defaults, `with_config(path, run_id, buffer_size, full_meta, lossy, accuracy)` | Creates the tables, as the source does. `Error::InvalidValue` for a zero or over-limit buffer or a non-finite accuracy; the source passes them unchecked. |
| `~MSDataSqlConsumer()` | `finish() -> Result<()>`, and the same on drop (warning on failure) | Flushes, then writes the run-level information unless a run was added. |
| `flush()` | `flush() -> Result<()>` | Spectra, then chromatograms, each one atomic handler write; a failed write keeps its buffer. |
| `addRun(filename, run_id)` | `add_run(&str, u64)` | |
| `setRunId(run_id)` | `set_run_id(u64)` | Bit 63 cleared by the handler. |
| `consumeSpectrum` / `consumeChromatogram` | trait methods | Buffer a copy, `clear(false)` the caller's record, keep the cleared record for the snapshot under full metadata, flush at `buffer_size`. |
| `setExpectedSize`, `setExperimentalSettings` | trait methods | Ignored, as in the source: the snapshot carries default settings. |
| protected `filename_`, `handler_`, `flush_after_`, `full_meta_`, `spectra_`, `chromatograms_`, `peak_meta_`, `wrote_any_run_` | `filename()`, `handler()`, `flush_after()`, `full_meta()`, `buffered_spectra()`, `buffered_chromatograms()`, `peak_meta()`, `wrote_any_run()` | Read-only. |
| *(none)* | `set_loss_policy`, `loss_policy` | The handler's `LossPolicy`. `Reject`, the default, refuses metadata sqMass would drop; `Source` is the source's silent drop. Under full metadata the handler cannot check a batch against a snapshot that is written only at the end, so the consumer writes batches under `Source` and applies the handler's own full-metadata refusals (auxiliary arrays under `Reject`, peptide identifications, empty native IDs) to each record as it arrives, before the record is taken. |

### `MSDataCachedConsumer`

| C++ | Rust | Notes |
|---|---|---|
| `MSDataCachedConsumer(filename, clearData = true)` | `new(path)` with the default, `with_clear_data(path, clear_data)` | Creates the file and writes the identifier, as the source does. `Error::Io` when the file cannot be created; the source's `ofstream` fails silently and loses every later write. |
| `~MSDataCachedConsumer()` | `finish() -> Result<()>`, and the same on drop (warning on failure) | Writes the spectrum and chromatogram counts and closes the file; until then the file is not a valid cache. |
| `consumeSpectrum` | trait method | Writes the record with the handler's `writeSpectrum_` (`CachedMzMLHandler::write_spectrum`, crate-private), counts it, and with `clearData` calls `clear(false)` on the caller's spectrum. After a chromatogram: `Error::InvalidValue` with the source's message "Cannot write spectra after writing chromatograms." (`Exception::IllegalArgument`), nothing written. |
| `consumeChromatogram` | trait method | As `consumeSpectrum`, with `writeChromatogram_`; no ordering check, as in the source. |
| `setExpectedSize`, `setExperimentalSettings` (both `{;}`) | trait methods | Ignored, as in the source. |
| base `Internal::CachedMzMLHandler` | `handler()`, `handler_mut()` | Composition in place of the source's public inheritance: the handler's reading and indexing members. |
| protected `ofs_`, `clearData_`, `spectra_written_`, `chromatograms_written_` | `filename()`, `clear_data()`, `spectra_written()`, `chromatograms_written()` | Read-only. |
| typedefs `SpectrumType`, `ChromatogramType` (private) | `MSSpectrum`, `MSChromatogram` | |
| `OPENMS_POSTCONDITION`s after clearing | none | The kernel's `clear(false)` empties the peaks and all three data-array families, which is what they assert. |

## Evidence

- `tests/ms_data_consumers.rs` (feature `mzml`, 29 tests): every section of
  `NoopMSDataConsumer_test.cpp`, `MSDataStoringConsumer_test.cpp`,
  `MSDataTransformingConsumer_test.cpp`, `MSDataChainingConsumer_test.cpp` and
  `MSDataAggregatingConsumer_test.cpp`, transcribed with line numbers (tier 3),
  on the unmodified upstream `MzMLFile_1.mzML`; plus `native_*` tests for the
  port's own contract (tier 4). The aggregating sections run with a merge
  double that adds intensities on the shared m/z grid; the class test's
  17/120/230 are reached by any merge that adds on that grid, so they check the
  grouping and metadata plumbing, not `SpectrumAddition`. The class tests'
  `setComment` is a metadata entry, since `MSSpectrum` has no comment field.
- `tests/ms_data_sql_consumer.rs` (feature `sqmass`, 11 tests): no upstream
  class test exists. Tier 4, derived from `MSDataSqlConsumer.cpp`: streamed
  records read back identically to a `write_experiment` of the same experiment
  for buffers of 1, 2 and 500, with and without full metadata; buffer and
  flush behaviour; the caller's record cleared; `addRun`/`setRunId`; the run
  written once at the end, named after the sqMass file; drop equivalent to
  `finish`; the refusals.
- `tests/ms_data_cached_consumer.rs` (feature `mzml`, 13 tests): every
  section of `MSDataCachedConsumer_test.cpp` (tier 3) on the unmodified
  `MzMLFile_1.mzML`. `release_build_bytes` streams that file through the
  consumer and requires the bytes the Linux x86_64 Release build's
  `CachedMzMLHandler::writeMemdump` wrote for it
  (`tests/data/f4_progress/MzMLFile_1.cached` and `empty.cached`, from
  `../oracle/F4-progress` on ibminode06): `writeMemdump`
  (`CachedMzMLHandler.cpp:33-58`) and the consumer's constructor,
  `consumeSpectrum`, `consumeChromatogram` and destructor
  (`MSDataCachedConsumer.cpp:18-75`) write the same identifier, records and
  counts (tier 1 bytes, tier 4 equivalence). The `release_consumer_*` tests
  compare with the Release build's own `MSDataCachedConsumer` directly
  (`../oracle/integ-F-cachedconsumer`, ibminode06, run twice, identical): all
  records with `clearData`, spectra 0 and 1 without it, nothing, and a
  spectrum after a chromatogram give the same four files byte for byte
  (`tests/data/ms_data_cached_consumer/`), the same cleared or unchanged
  records, and the source's `IllegalArgument` message (tier 1). Native tests
  cover the cached and cleared data arrays and an uncreatable file.

## C++ observations

- `MSDataTransformingConsumer.h` documents `setExperimentalSettings` as
  "ignored", but the `.cpp` calls `lambda_exp_settings_` when set.
- `MSDataAggregatingConsumer_test.cpp` constructs one with a null next consumer
  ("don't do that ..."); the first flush dereferences it. A reference cannot be
  null here.
- `MSDataSqlConsumer`'s `buffer_size` is an `int` stored into `size_t
  flush_after_`: a negative value becomes a huge flush threshold and a huge
  `reserve`.
- `MSDataCachedConsumer.h` describes a consumer that "is able to transform a
  spectrum on the fly ... using a function pointer"; the class has no such
  member (that is `MSDataTransformingConsumer`).
- `MSDataCachedConsumer` opens its `ofstream` unchecked and never tests the
  stream: an uncreatable output file, or a full disk, loses the data silently.

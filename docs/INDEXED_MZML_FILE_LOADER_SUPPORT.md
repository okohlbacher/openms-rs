# `IndexedMzMLFileLoader`

Source pin: core `bc9cc12`. `FORMAT/IndexedMzMLFileLoader.h` and its `.cpp`
are ported as `src/format/indexed_mzml_file_loader.rs`, feature `mzml`. The
header's one TOPP consumer is `TICCalculator`.

Before this port each operation existed separately (`OnDiscMSExperiment::open_file`,
`OnDiscMSExperiment::load_experiment`, the streaming writer, `mzml::store`); the
type now joins them behind the loader's own `PeakFileOptions`, as the source
does.

## API mapping

| C++ | Rust | Notes |
|---|---|---|
| `typedef OnDiscMSExperiment OnDiscPeakMap` | `crate::kernel::on_disc_experiment::OnDiscPeakMap` | Already present. |
| `IndexedMzMLFileLoader()` / `~IndexedMzMLFileLoader()` | `new()`, `Default`, drop | |
| `PeakFileOptions& getOptions()` | `options_mut()` | |
| `const PeakFileOptions& getOptions() const` | `options()` | |
| `void setOptions(const PeakFileOptions&)` | `set_options(PeakFileOptions)` | |
| `bool load(filename, OnDiscPeakMap&)` | `load(path, &mut OnDiscPeakMap) -> Result<bool>` | `exp.openFile(filename)`: `Ok(true)` when an index parsed, `Ok(false)` for a readable file without a usable index. An unreadable file is `Err`, not an escaping exception. The loader's options are not passed on, as in the source. The source's `@note` that a corrupt `indexListOffset` throws `ConversionError` does not apply: the on-disc index parser treats every unusable index as "not indexed", so that is `Ok(false)`. |
| `void store(filename, OnDiscPeakMap&)` | `store_on_disc(path, &mut OnDiscPeakMap)` | Streams record by record through `PlainMSDataWritingConsumer`, so one record is in memory at a time, exactly the source loop. Sets `options.write_index = true` persistently first, then hands the options to the consumer (`MSDataWritingConsumer::set_peak_file_options`, the source's `consumer.setOptions(options_)`), so precision, zlib, Numpress and TPP compatibility apply to every record. The output is byte-identical to `store` of the same records with the same options. See the native differences. |
| *(none)* | `store_on_disc_with_policy(path, exp, ReferencePolicy)` | The streaming writer's reference policy; `SourceDangling` reproduces the source's dangling references. |
| `void store(filename, PeakMap&)` | `store(path, &MSExperiment)` | `MzMLFile::store` with the loader's options, after `write_index = true` (persistent). |
| private `PeakFileOptions options_` | private `options` | |

## Native differences

1. **A store of an on-disc map without metadata is refused**
   (`Error::MissingInformation`) before a byte is written. The source
   dereferences the null `getExperimentalSettings()` pointer
   (`IndexedMzMLFileLoader.cpp:47`).
2. **Records whose history differs from the first record's are refused by
   default** (`ReferencePolicy::Checked`), where the source writes a dangling
   reference; `store_on_disc_with_policy(.., SourceDangling)` is the source's
   output. See `docs/MS_DATA_WRITING_CONSUMER_SUPPORT.md`.
3. **The on-disc store is published atomically**, plain even for a `.gz`/`.bz2`
   name (the source writes through an `std::ofstream`), so a failure leaves no
   partial file. The in-memory store compresses by suffix, as `XMLFile::save_`
   does.
4. **An in-memory experiment with no records is refused** by the mzML writer
   (`Error::InvalidValue`): an indexed mzML needs a record to index, and the
   source writes a dummy `-1` offset (CPP-050).

## Evidence

`tests/indexed_mzml_file_loader.rs` (feature `mzml`, 13 tests): every section of
`IndexedMzMLFileLoader_test.cpp`, with line numbers (tier 3), on the unmodified
upstream `IndexedmzMLFile_1.mzML` and `MzMLFile_1.mzML`, including full record
and settings equality between the on-disc map and a whole-file load after both
stores. The upstream `VALIDATE_TMP_FILES` step (XSD validation of the written
files) runs under the `mzml-schema` feature. The `native_*` tests (tier 4) pin
the differences above, and one invariant the source implies without testing:
its two stores reach the same `MzMLHandler` writer with the same `options_`
(`IndexedMzMLFileLoader.cpp:49`, `:66`), so the streamed on-disc store and the
whole-document store of the same records are compared byte for byte under nine
writer configurations (defaults, zlib, 32-bit m/z, 64-bit intensity, TPP
compatibility and four Numpress combinations).

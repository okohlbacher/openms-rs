# Cached mzML support

Ports of `FORMAT/HANDLERS/CachedMzMLHandler.h` (`OpenMS::Internal`) and
`FORMAT/CachedMzML.h` (core `bc9cc12`): the binary on-disk cache of spectra
and chromatograms and the random-access reader over a cached mzML pair
(metadata `.mzML` plus binary `.mzML.cached`). Rust:
`src/format/cached_mzml_handler.rs`, `src/format/cached_mzml.rs` (feature
`mzml`; the handler's binary half builds without it). Tests:
`tests/cached_mzml.rs` (14, file-gated on `mzml`) and the progress half of
`tests/f4_progress.rs` (3, file-gated on `mzml` and `consensusxml`).
Provenance: `tests/data/cached_mzml_provenance.json`,
`tests/data/f4_progress_provenance.json`.

## API mapping

| Source member | Rust |
|---|---|
| `CachedMzMLHandler()`, `~CachedMzMLHandler()`, `operator=` | `CachedMzMLHandler::new()` / `Default`, drop, `Clone` |
| typedefs `MapType`, `SpectrumType`, `ChromatogramType`, `DatumSingleton`, `Datavector` | `MSExperiment`, `MSSpectrum`, `MSChromatogram`, `f64`, `Vec<f64>` |
| `writeMemdump(exp, out) const` | `write_memdump(&exp, out)`, `write_memdump_to(&exp, writer)` |
| `writeMetadata(exp, out_meta, addCacheMetaValue)` | `write_metadata(&exp, out_meta, add_cache_meta_value)` (mzML) |
| `writeMetadata_x(exp, out_meta, addCacheMetaValue)` | `write_metadata_x` (mzML) |
| `readMemdump(exp_reading, filename) const` | `read_memdump(filename) -> MSExperiment` |
| `createMemdumpIndex(filename)` | `create_memdump_index(filename)` |
| `getSpectraIndex()`, `getChromatogramIndex()` | `spectra_index() -> &[u64]`, `chromatogram_index() -> &[u64]` |
| `readSpectrumFast` (both overloads) | `read_spectrum_fast(reader) -> (Vec<BinaryDataArray>, ms_level, rt)` |
| `readChromatogramFast` (both overloads) | `read_chromatogram_fast(reader) -> Vec<BinaryDataArray>` |
| `readSpectrum(spectrum, ifs)`, `readChromatogram(chromatogram, ifs)` | `read_spectrum(&mut spectrum, reader)`, `read_chromatogram(&mut chromatogram, reader)` |
| protected `writeSpectrum_`, `writeChromatogram_`, `readDataFast_`, `spectra_index_`, `chrom_index_` | private functions and fields |
| `OpenSwath::BinaryDataArrayPtr` | `BinaryDataArray` (`data`, `description`); the OpenSWATH type itself is not ported |
| file identifier `8094` | `CACHED_MZML_FILE_IDENTIFIER` |
| `ProgressLogger` base: the calls of `writeMemdump`, `readMemdump`, `createMemdumpIndex` | `write_memdump_with_progress`, `read_memdump_with_progress`, `create_memdump_index_with_progress`, taking the `&mut ProgressLogger` the source object is |
| `CachedmzML()`, `CachedmzML(filename)`, copy constructor, destructor | `CachedmzML::new()`, `CachedmzML::open(filename)`, `Clone` (reopens the side-car), drop |
| `getSpectrum(id)`, `getChromatogram(id)` | `spectrum(id)`, `chromatogram(id)` |
| `getNrSpectra()`, `getNrChromatograms()`, `getMetaData()` | `nr_spectra()`, `nr_chromatograms()`, `meta_data()` |
| `static store(filename, map)`, `static load(filename, CachedmzML&)` | `CachedmzML::store`, `CachedmzML::load` |
| protected `load_`, `meta_ms_experiment_`, `ifs_`, `filename_`, `filename_cached_`, `spectra_index_`, `chrom_index_` | private; `filename()`, `cached_filename()` accessors |

## Preserved source conventions

- The binary layout field for field (identifier, per-record counts, MS level
  and RT for spectra, m/z-or-RT then intensity doubles, named arrays, the two
  trailing counts); little-endian with 64-bit sizes, which is every platform
  the source supports.
- Float arrays before integer arrays, all as doubles; everything reads back
  as float arrays; string arrays and metadata other than MS level and RT live
  in the metadata mzML.
- `writeMetadata` with the cache flag adds one shared data-processing entry
  ("File format conversion", `cached_data` = `true`).
- The side-car name is the mzML path with `.cached` appended.
- Progress: `startProgress(0, spectra + chromatograms, label)` with the
  source's three labels, `setProgress(i)` per record with the index
  restarting at 0 for chromatograms, `endProgress()` after the last record
  (after the close when writing); a refused header makes no call; the file's
  counts reach the logger with the source's unsigned-to-signed wrap-around.

## Native differences

- The source skips the arrays of a record without peaks while counting them
  in the header, so its reader then misreads the next record; this port
  writes them, which both readers read correctly.
- Write failures are `Error::Io`; the source writes to an unchecked stream.
  With progress, an output that cannot be created fails before the start,
  where the source still makes every call.
- An array name longer than `MAX_ARRAY_NAME_BYTES` (1023) is refused; the
  source skips it and reuses the previous name.
- A record extending past the end of the file is refused in
  `create_memdump_index`; the source records the offset. The previous indices
  are kept on error; `CachedmzML::load` leaves the target unchanged on error.
- `spectrum(id)`/`chromatogram(id)` check `id` (a debug-only assertion in the
  source).

## Checked boundaries and evidence

- Tier 1: every progress call (arguments, nesting depth), the outcome and
  the depth afterwards of the three members on `MzMLFile_1.mzML`, an empty
  experiment and a file with the wrong magic number, and the written dumps
  byte for byte, against the Linux x86_64 Release build run twice on
  ibminode06 (`../oracle/F4-progress`, `tests/f4_progress.rs`).
- Tier 3: every `START_SECTION` of `CachedMzMLHandler_test.cpp` and
  `CachedMzML_test.cpp` over the class-test file `MzMLFile_1.mzML` (held
  byte-identical as `tests/data/mzml_validator/MzMLFile_1.mzML`).
- Tier 4: truncated and forged dumps are refused; records without peaks keep
  their arrays and chromatogram arrays keep their values. The array-name
  ceiling is not exercised by a test.

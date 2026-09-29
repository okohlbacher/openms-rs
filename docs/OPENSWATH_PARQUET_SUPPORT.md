# OpenSWATH Parquet files and the MSExperiment Arrow export

Lane F1b (wave 4). Everything here is behind the optional `arrow` feature,
which now also enables `numpress` (the `.xic`/`.xim`/`.xipm` arrays are
MS-Numpress + zlib). C++ pin: core `bc9cc12`.

| Header | Status | Rust |
|---|---|---|
| `FORMAT/XICParquetFile.h` | complete | `src/format/xic_parquet_file.rs` |
| `FORMAT/XIMParquetFile.h` | complete | `src/format/xim_parquet_file.rs` |
| `FORMAT/XIPMParquetFile.h` | complete | `src/format/xipm_parquet_file.rs` |
| `FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h` | complete | `src/format/ms_chromatogram_parquet_consumer.rs` |
| `FORMAT/DATAACCESS/MobilogramParquetConsumer.h` | complete | `src/format/mobilogram_parquet_consumer.rs` |
| `FORMAT/DATAACCESS/XIPMParquetConsumer.h` | complete | `src/format/xipm_parquet_consumer.rs` |
| `FORMAT/MSExperimentArrowExport.h` | partial | `src/format/msexperiment_arrow_export.rs` |
| `ANALYSIS/OPENSWATH/PeakMapExtractor.h` | partial (result type only) | `src/analysis/openswath/peak_map_extractor.rs` |
| `FORMAT/ConsensusMapArrowExport.h` | untouched | needs `QPXIdentity.h`, `QPXValueValidation.h`, `IDScoreSwitcherAlgorithm.h`, full `SpectrumLookup.h`/`SpectrumNativeIDParser.h` |
| `FORMAT/ProteinGroupArrowExport.h` (+ `_impl.h`) | untouched | needs `ConsensusMapArrowExport.h`, `QPXIdentity.h` |

The private `FORMAT/DATAACCESS/ParquetReaderHelpers.{h,cpp}` is ported as the
crate-private `src/format/parquet_reader_helpers.rs`; the three writers'
repeated anonymous-namespace code (transition maps, id assignment, Numpress
encoding, the ZSTD-11 `FileWriter`) is one crate-private
`src/format/openswath_parquet_writer.rs`.

## Readers (`XICParquetFile`, `XIMParquetFile`, `XIPMParquetFile`)

| C++ | Rust |
|---|---|
| `XICParquetFile(filename)` / `(filenames)` | `XicParquetFile::new` / `from_files` (`Error::Io` NotFound for `FileNotFound`; empty list `InvalidValue`, XIPM: NotFound "No .xipm files given") |
| `getFilename()`, `getFilenames()` | `filename()`, `filenames()` |
| `load(out)` | `load()` |
| `getChromatograms(out, precursor_id, …, run_id, filter)` | `chromatograms(&XicQuery)`; `-1`/`""` defaults are `None`/`""`, a negative value is ignored as in the source |
| `getChromatograms(out, const ParquetFilter&)` / `(…, const ParquetFilterBuilder&)` | `chromatograms_with_filter` / `chromatograms_with_builder` |
| `getRuns(out)` | `runs()` |
| `getAnalytes(out, columns, nest_transitions)` | `analytes(&columns, nest)` (the source's defaults: `&[]`, `true`) |
| `getColumns(out)` | `columns()` |
| `XICChromatogram`, `XICRunInfo`, `XICAnalyte` (and XIM/XIPM twins) | `XicChromatogram`, `XicRunInfo`, `XicAnalyte`, …; each `has_x`/`x` pair is one `Option` |
| XIM: `getMobilograms(…, mobilogram_type, feature_id, feature_rt, filter)` | `XimParquetFile::mobilograms(&XimQuery)`; `feature_rt` is matched within 1e-9, row by row |
| XIPM: `getPeakMaps(…, peakmap_type)` | `XipmParquetFile::peak_maps(&XipmQuery)` |

`XICRunInfo`/`XIMRunInfo`/`XIPMRunInfo` and `XICAnalyte`/`XIMAnalyte` are
member-for-member identical in the source; the port has one `RunInfo` and one
`Analyte`, re-exported under the three names.

**Filter evaluation.** The source tries Arrow Dataset push-down, then Arrow
compute, and falls back to its own row evaluator when Arrow cannot bind or run
the expression. Every expression the parser builds is a left-to-right chain of
comparisons joined by AND/OR, without negation, so Kleene logic with null rows
dropped selects the same rows as the row evaluator with a null comparison
counted false. The port evaluates row by row with the operand rules the two
paths give together: an integer condition compares against any integer column
(unsigned values above `i64::MAX` are an error) or a boolean one; a string
condition compares byte-wise against `string`, `large_string`, `binary` and
`large_binary`; anything else is an error once a non-null value must be
compared. Native difference: an integer condition on a floating-point column,
which Arrow compute would evaluate after casting the literal, is refused (no
OpenSWATH schema has a filterable floating-point column). The Arrow Dataset
push-down is not ported; it is an optimisation with the same result.

**Resource ceilings (native).** `ReaderLimits`: each file's footer is checked
(rows, uncompressed bytes) before decoding, the combined table is capped at
2^32 rows, each binary cell decodes within the Numpress coder limits (which
bound zlib inflation), and one call decodes at most 2^31 values.

**Other native differences.** A file without row groups reads as an empty
table (the source's `getColumn` refuses its chunk-less columns). Raw
compression-0 arrays are read little-endian (the source `memcpy`s host order).

## Writers (the three consumers)

| C++ | Rust |
|---|---|
| `MSChromatogramParquetConsumer(filename, run_id, source_file, transition_exp)` | `MSChromatogramParquetConsumer::new` (infallible, as the source's constructor) |
| `consumeSpectrum`, `consumeChromatogram`, `setExpectedSize`, `setExperimentalSettings` | the `MSDataConsumer` trait; `consume` is `consumeChromatogram` |
| `finalize()`, destructor | `finalize() -> Result<()>`; `Drop` finalises and discards the error, as the destructor logs it |
| `MobilogramParquetConsumer(…)`, `consumeMobilogram(m, type, ms_level, transition_id, native_id, feature_rt, feature_id)`, `setExpectedSize`, `finalize` | `new`, `consume_mobilogram(&m, &MobilogramAnnotation)`, `set_expected_size`, `finalize` |
| `XIPMParquetConsumer(filename, transition_exp)`, `consumePeakMap(peak_map, run_id, source_file, ms_level)`, `setExpectedSize`, `finalize` | `new`, `consume_peak_map`, `set_expected_size`, `finalize` |
| `PeakMapExtractor::ExtractedPeakMap` | `analysis::openswath::peak_map_extractor::ExtractedPeakMap` (source defaults: target RT NaN, IM -1, window 0..-1) |

The transition maps follow `buildTransitionMaps_` exactly: compounds are
numbered first (an integer id that `StringUtils::toInt64` accepts and is free
keeps its value, anything else takes the next free id), then transitions;
duplicated ids keep their first entry; a precursor is a decoy only if it has
decoy transitions and no target transition. The mobilogram writer rejects
negative parsed ids and clears the run id's top bit, as its source does.

**Native differences.**

- *Arrays that the lossy encoding cannot represent are stored losslessly*
  under identifier 1 (zlib-compressed doubles), which every source reader
  decodes: a Numpress-linear fixed point that is not positive (the source's
  mass-accuracy estimator returns 0 for fewer than three values, so the
  source stores a one- or two-point chromatogram with fixed point 0), first
  values that do not quantise into the unsigned 32-bit range Numpress linear
  stores them in (negative retention times, for example), negative or
  non-finite values for slof, and encoder rejections.
- *Staged output.* The file is written next to the destination and renamed
  onto it by `finalize`; the source opens and truncates the destination at its
  first flush, so a later failure leaves a truncated file.
- *Atomic consume.* A chromatogram or peak map whose native id the library
  does not know is refused before anything is buffered; the source has
  already appended `RUN_ID`, `SOURCE_FILE` and `MS_LEVEL` (and `PEAKMAP_TYPE`)
  when it throws, leaving its builders ragged.
- *No writes after `finalize`.* The chromatogram writer refuses them; the
  source's buffers them and would reopen and truncate the finished file at its
  next flush.

The Parquet bytes differ from the source's (arrow-rs vs Parquet C++ footers;
`miniz_oxide` vs zlib streams at the same level 6); decoded values are
identical (see Evidence).

## `MSExperimentArrowExport`

| C++ | Rust |
|---|---|
| `ArrowExportFormat`, `ArrowSpectraExportConfig`, `ArrowChromatogramExportConfig` | same names, public fields, source defaults |
| `ParquetWriteConfig` | `arrow_io_helpers::ParquetWriteConfig` (lane F1a), re-exported |
| `getSpectraArrowColumnNames`, `getChromatogramArrowColumnNames` | `spectra_arrow_column_names` (`Result`: an RT bound needs sorted spectra), `chromatogram_arrow_column_names` |
| `exportSpectraToArrowCDataInterface`, `exportChromatogramsToArrowCDataInterface` | `spectra_to_arrow`, `chromatograms_to_arrow` return the `RecordBatch` the source exports. **Gap:** the C ABI export itself (`ArrowSchema`/`ArrowArray` structs) needs arrow-rs's `ffi` feature (`arrow-array/ffi`, which adds the `libc` edge to `Cargo.lock`) |
| `exportSpectraToParquet`, `exportChromatogramsToParquet` | `export_spectra_to_parquet`, `export_chromatograms_to_parquet` (`Result<()>` for `bool` + log), through F1a's `write_table_to_parquet` |

Native differences: an RT or m/z bound on unsorted spectra/peaks is
`Error::UnsortedData` (the source's binary search assumes sorted input); an
inverted window selects nothing (the source's iterator loop runs past the
end); an ion-mobility array shorter than its peak list is an error (the source
reads past it); a column beyond Arrow's 32-bit offsets is an error before
anything is built (arrow-rs builders would panic; the source's fail).

## Evidence

- Tier 3: every section of `XICParquetFile_test.cpp`, `XIMParquetFile_test.cpp`,
  `XIPMParquetFile_test.cpp`, `MSChromatogramParquetConsumer_test.cpp` (the
  TraML-based sections with the library written out by hand, since
  `TraMLFile` is unported), `MobilogramParquetConsumer_test.cpp`,
  `XIPMParquetConsumer_test.cpp` and `MSExperimentArrowExport_test.cpp`.
- Differential against retained C++-written files: the 18 chromatograms of
  `MSChromatogramParquetConsumer_1_output.chrom.mzML`, written here, decode to
  bit-identical RT and intensity arrays and identical metadata as
  `XICParquetFile_1_input.xic`, which the C++ consumer wrote; the two peak maps
  of the XIPM class test, written here, equal `XIPMParquetFile_1_input.xipm`
  field by field.
- Tier 1: the executed Release oracle `../oracle/F1b-openswath` (ibminode06,
  Linux x86_64 Release `openms4-release-bc9cc12-c19e494-174b576`, g++ 14.4
  `-O2 -DNDEBUG`, Arrow C++ 23.0.1; run twice, records and files
  byte-identical), `tests/openswath_parquet_oracle.rs`: 28 textual filters, 7
  argument/typed filters, 12 analyte listings, columns, runs and a decoded row
  of the `.xic` reader; 6 filters, 2 argument queries, 4 analyte listings and
  two decoded rows of the `.xim` reader; the three writers on one library
  (numeric and non-numeric ids, a negative id, a decoy-only precursor, an
  orphan reference, empty and one-point arrays), where the C++-written files
  read here equal the C++ reader's rows bit for bit and the Rust-written files
  do too (except the arrays below); and seven `MSExperimentArrowExport`
  tables (long/semi-wide, RT/m/z/level windows, column selection), whose
  C++-written Parquet files equal the Rust tables column by column. See
  `tests/data/openswath_parquet_provenance.json`.

Measured on the Release build and recorded in the oracle:

- A one- or two-point chromatogram written by `MSChromatogramParquetConsumer`
  reads back with RT `-nan` (fixed point 0), and negative RTs come back as
  about 4.29e8 (wrapped unsigned integers). The Rust writer stores both
  losslessly (identifier 1); the Rust reader decodes such source-written
  payloads as the source does (NaN), so these files stay readable.
- One chromatogram with an unknown native id makes the source's `finalize()`
  fail ("Column 3 named PRECURSOR_ID expected length 2 but got length 1") and
  the file holds no chromatogram at all.
- `MSChromatogramParquetConsumer` accepts a chromatogram after `finalize()`
  and silently drops it; the mobilogram and peak-map writers throw
  `IllegalArgument`.
- Arrow compute is not initialised in this build, so every filter falls back
  to the row evaluator, logging two warnings per call.

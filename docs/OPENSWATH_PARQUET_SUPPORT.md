# OpenSWATH Parquet files, the MSExperiment Arrow export and the `.oswpq` containers

Part 1 is group F (lane F1b), part 2 group O (lanes O5a to O5c).

## Part 1: the XIC/XIM/XIPM files and the MSExperiment Arrow export

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
| `ANALYSIS/OPENSWATH/PeakMapExtractor.h` | complete (the result type here, the extractor by lane O4a; see `docs/OPENSWATH_EXTRACTION_SUPPORT.md`) | `src/analysis/openswath/peak_map_extractor.rs` |
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

## Part 2: the `.oswpq` containers and the OpenSWATH exports (group O)

Ports the OpenSWATH Parquet I/O of `ANALYSIS/OPENSWATH` at OpenMS4-core
`bc9cc12`: the `.oswpq` assay library (`TransitionParquetFile.h`), the OSW
results container written by OpenSwathWorkflow and read by PyProphet-style
consumers (`OpenSwathOSWParquetWriter.h`, `OpenSwathOSWParquetReader.h`), and
the flat scored-table exports (`OpenSwathParquetExporter.h` with the row and
table types of `OpenSwathExportData.h`). The container code sits behind the
`arrow` feature (Arrow/Parquet, lane F1a) together with `file-compression`
(the ZIP archive layer, `ZipArchiveFile`, lane F7b), which the default feature
set already enables; the exporter needs `arrow` only, and the export data
types need no feature.

| Header | Rust | Status |
|---|---|---|
| `TransitionParquetFile.h` | `src/analysis/openswath/transition_parquet_file.rs` | complete |
| `OpenSwathOSWParquetWriter.h` | `src/analysis/openswath/openswath_osw_parquet_writer.rs` | complete |
| `OpenSwathOSWParquetReader.h` | `src/analysis/openswath/openswath_osw_parquet_reader.rs` | complete |
| `OpenSwathExportData.h` | `src/analysis/openswath/openswath_export_data.rs` | complete |
| `OpenSwathParquetExporter.h` | `src/analysis/openswath/openswath_parquet_exporter.rs` | complete |

### The `.oswpq` container

A directory, or a store-only ZIP archive of one:

```text
<name>.oswpq
├── metadata.json              container summary (results containers only)
├── library/
│   ├── precursors.parquet     OSWPrecursorSchema
│   ├── transitions.parquet    OSWTransitionSchema
│   └── metadata.json          QC counts
└── runs/                      results containers only
    ├── runs.parquet
    └── run_id=<id>/
        ├── features.parquet
        ├── feature_precursor.parquet
        └── feature_transition.parquet
```

A zipped container written by `TransitionParquetFile` also carries a JSON size
index as an entry named after the **staging** archive, `<name>.oswpq.tmp.idx.json`
(the source indexes `<path>.tmp` before renaming it into place; measured on the
Release build, and reproduced).

### TransitionParquetFile: API mapping

| C++ member | Rust | Notes |
|---|---|---|
| `TransitionParquetFile()` / `~TransitionParquetFile()` | `TransitionParquetFile::new()`, `Default` | stateless |
| `convertParquetToTargetedExperiment(oswpq_dir, targeted_exp)` | `convert_parquet_to_targeted_experiment(path, &mut exp)` | |
| `convertParquetToTargetedExperiment(oswpq_dir, targeted_exp, source_ids*)` | `convert_parquet_to_targeted_experiment_with_source_ids(path, &mut exp, Option<&mut SourceIDMapping>)` | `nullptr` is `None` |
| `convertLightTargetedExperimentToParquet(oswpq_path, targeted_exp)` | `convert_light_targeted_experiment_to_parquet(path, &exp)` | the compatibility entry point that normalises source-style ids |
| `convertLightTargetedExperimentToParquet(oswpq_path, targeted_exp, source_ids*)` | `convert_light_targeted_experiment_to_parquet_with_source_ids(path, &exp, Option<&SourceIDMapping>)` | |
| file-static `OpenMSLibraryStats`, `writeLibraryMetadata_`, `jsonMapByClass_`, `jsonEscape_` | private `LibraryStats`, `library_metadata`, `json_map_by_class`; `ParquetFile::json_escape` | |
| file-static `PrecursorInfo`, `addModifications_`, `joinProteinAccessions_` | private `PrecursorInfo`, `add_modifications`, `join(";")` | |
| entry names `library/precursors.parquet`, `library/transitions.parquet`, `library/metadata.json` | `PRECURSORS_ENTRY`, `TRANSITIONS_ENTRY`, `METADATA_ENTRY` | native constants |

### Preserved source conventions

- Reader: tables validated against the registry schemas in **subset** mode:
  a column the schema does not know is ignored, a known column of the wrong
  type is refused, and a missing column is not reported by validation (the
  optional columns `traml_id`, `annotation`, `library_drift_time`, … may be
  absent; a missing *required* column then fails in `getColumn`, "Missing
  required column 'type'").
- `protein_accessions` is read with `getStringList`: `"P2;;P9;"` gives
  `P2`, `""`, `P9`, `""`, and the empty accession becomes a protein of its
  own, exactly as the Release build does. A `list<string>` column, which the
  header documents as accepted, is refused by the schema validation first
  (measured; see the defect list).
- Compound `sequence` is `modified_sequence`, or `unmodified_sequence` when
  that is empty; modifications are N-terminal (`-1`), C-terminal (length),
  then residues, with UniMod record id `-1` for a mass tag.
- Fragment type is `type`, or `annotation` when `type` is empty (an
  annotation such as `b4^2` therefore becomes the unknown type).
- `charge`, fragment charge and ordinal are narrowed with the source's
  wrapping casts (`300` becomes fragment charge `44`, `70000` ordinal `4464`;
  measured).
- Defaults for null or absent optional values: drift time -1, library RT and
  intensity 0, decoy false, detecting and quantifying true, identifying false,
  ordinal -1.
- `traml_id` provenance: first non-empty value wins for transitions; an
  ambiguous precursor mapping is an error only when provenance is requested.
- Writer: precursor m/z from the first transition of the compound, precursor
  decoy when any of its transitions is, `traml_id` from the provenance or the
  canonical id, `protein_accessions` joined by `;`, `unmodified_sequence`
  from `AASequence`; `metadata.json` byte for byte, charges in ascending
  order, fragment classes `b`/`y`/other on the lower-cased type string
  (`unknown` counts as other), proteins starting with `DECOY_` as decoys.
- Writer to a path that is not a directory: staged in a temporary
  directory, entries added in creation order (precursors, transitions,
  metadata), sidecar index, then renamed over an existing file.

### Native differences

- **Compound and protein order on read.** The source materialises precursor
  rows in a `std::unordered_map<int64_t, …>` and emits compounds (and the
  proteins they introduce) in its iteration order: on the Release build's
  libstdc++ that is the reverse of insertion for the fixtures measured, and
  it depends on the standard library in general. The port emits precursor
  row order. Every field and the transition order are identical.
- **Unparseable sequences.** The source catches only
  `Exception::InvalidValue` from `AASequence::fromString` (an unknown
  modification); a `ParseError` (a missing `)`, an unknown residue letter)
  escapes and aborts the read or the write (measured: `PEPT(UniMod:35`). The
  port's `AASequence` reports both as `Error::InvalidValue`, so every
  unparseable sequence is handled as the source handles an unknown
  modification: no modifications, empty `unmodified_sequence`.
- **`build_time`** is written as `""`: the crate records no build time, and
  two writes of the same library give identical bytes.
- **Archive entries** are read by extracting them to a temporary file, the
  source's fallback path; its first attempt (`ZipRandomAccessFile`, a reader
  over a stored entry) reads the same bytes and is not ported.
- **Failure leaves the outputs empty.** The source resets the experiment
  first and leaves it partially filled when it throws; the port leaves the
  experiment and the provenance mapping empty on every error.
- **Error types.** `Exception::FileNotFound` is `Error::Io` of kind
  `NotFound`; `InvalidValue` and `MissingInformation` map to the variants of
  the same name, their message and offending value in one string.

### OSW results container: API mapping

| C++ member | Rust | Notes |
|---|---|---|
| `OpenSwathOSWParquetWriter()` | `OpenSwathOSWParquetWriter::new()`, `Default` | |
| `write(output_path, assay_library, feature_map, run_id, input_filename, enable_uis_scoring)` | `write(path, &library, &features, run_id, input_filename, uis) -> Result<Vec<String>>` | the returned strings are the source's `OPENMS_LOG_WARN` lines |
| `write(..., source_ids*)` | `write_with_source_ids(..., Option<&SourceIDMapping>)` | `nullptr` is `None` |
| `setPreserveExisting(bool)` | `set_preserve_existing(bool)`, getter `preserve_existing()` | |
| `OpenSwathOSWParquetReader()`, `(oswpq_dir)` | `new()`, `from_path(path) -> Result<Self>` | |
| `load`, `oswpqPath`, `rows` | `load`, `oswpq_path`, `rows` | |
| `fetchPeakGroupFeatures(dir, level, main_score)` | `fetch_peak_group_features` | the main-score warning is returned in `warnings` |
| `fetchTransitionFeatures`, `fetchUnscoredData` | `fetch_transition_features`, `fetch_unscored_data` | `UnscoredResult` members in snake case |
| `Row`, `PeakGroupFeatureScoresResult`, `TransitionFeaturesResult`, `UnscoredResult` | same names | |

### Scored-table exports: API mapping

| C++ member | Rust | Notes |
|---|---|---|
| `OpenSwathParquetExporter::writeFeatureScores(filename, table)` | `write_feature_scores(path, &table) -> Result<()>` | |
| `OpenSwathParquetExporter::writeTransitionScores(filename, table)` | `write_transition_scores(path, &table) -> Result<()>` | |
| enums `OpenSwathExportFileFormat`, `OpenSwathIPFExportMode`, `OpenSwathMatrixLevel`, `OpenSwathMatrixNormalization` and their `toString` | same enums (`TSV` is `Tsv`), `as_str()` and `Display` | `Default` is the first enumerator (a value-initialised `enum class`) |
| `OpenSwathExportRow`, `OpenSwathFeatureScoreRow`, `OpenSwathFeatureScoreTable`, `OpenSwathTransitionScoreRow`, `OpenSwathTransitionScoreTable`, `OpenSwathQuantMatrix` | same structs, `Default` = member initialisers | `std::optional` is `Option` |

### OSW results container: native differences

- **Scores read from meta values (C++ defect candidate).** The source writes
  every score it reads from a meta value as
  `appendOptionalFloat_(b, extractMetaDouble_(x, key, value), value, ..)`
  with `value` (one `double value = 0.0;` per feature) passed by value. C++
  leaves the evaluation order of function arguments unspecified, and the
  executed Linux Release build (GCC 14.4, `-O2`) copies `value` before the
  extraction runs: every such cell of `features.parquet`,
  `feature_transition.parquet` and `feature_precursor.parquet` holds the
  previous successful extraction on the same feature (`0.0` for the first),
  while its null flag is right. Measured on a feature with 31 successful
  extractions, all 31 cells shifted by one (for example the first
  transition's `total_area_intensity` is 7.25, the feature's
  `im_log_intensity`, where its `total_xic` is 20). The port writes each
  cell's own value, which is what a compiler evaluating left to right
  produces. `tests/openswath_osw_parquet.rs` carries the chain (`STALE`) and
  requires the Release and port files to differ in exactly those cells.
- **Preserved archives.** The Release reader cannot read the archive its own
  writer extended with `preserve_existing` (libzip rewrites the replaced
  entries deflated, and `ZipRandomAccessFile` then fails to seek:
  "Failed to create parquet reader from RandomAccessFile", with no fallback);
  the port reads entries by extraction and reads that archive.
- Archive entries are written in sorted path order (the source follows the
  file system's `recursive_directory_iterator` order); the three run tables
  are written one after another (the source uses `std::async`); `runs.parquet`
  is read through the type-checking accessors (the source's
  `static_pointer_cast` is undefined for another type); `build_time` is `""`.
- Reproduced: an error raised while the features are converted leaves the
  created run directory behind, and the same `run_id` is then refused.
- Reader: both sorts replay the Release `std::sort` comparison by
  comparison, so ties keep the Release order; a failed `load` leaves `rows`
  empty (the source keeps the rows read before the error).

### Scored-table exports: native differences

- A failure to write is `Error::Io` carrying the reason; the source maps the
  `false` of `ArrowIOHelpers::writeTableToParquet` to
  `Exception::FileNotWritable` and loses it.

### Checked boundaries and evidence

Tier 1, executed differential: `../oracle/O5c-parquet` runs a driver against
the Release install `openms4-release-bc9cc12-c19e494-174b576` on ibminode06
(g++ 14.4.0, `-O2 -DNDEBUG`, Arrow C++ 23.0.1), twice, with identical record
streams and files. `tests/transition_parquet_file.rs` replays it:

- 24 reader cases on containers the Release build wrote or read: the class
  test's TraML library written through the 2-argument and the 3-argument
  writer, directory and ZIP; a synthetic canonical library with sparse ids
  (`0`, `3`, `5`, `8`, `11`, `17`), decoys, a metabolite, terminal, residue
  and nested-label modifications, every counted fragment class, negative and
  zero charges; an empty library; and hand-built tables with extra, missing
  and mistyped columns, duplicate and negative ids, unknown references, a
  null required value, an ambiguous `traml_id`, wrapping casts, a missing
  table and a missing container. Every compound, transition, protein and
  provenance entry is compared field by field, every error by kind and
  message.
- Writer: every table value by value (schema, nulls and bits) and
  `metadata.json` byte by byte except `build_time`, for the 2-argument and
  3-argument paths, with and without provenance, empty, and zipped (entry
  list and index keys); the four writer errors.
- Both sections of `TransitionParquetFile_test.cpp` transcribed.

Provenance: `tests/data/transition_parquet_file_provenance.json`.

Tier 1, executed differential for the results container:
`../oracle/O5c-oswpq` (same host, install and flags, twice, identical records
and files). `tests/openswath_osw_parquet.rs` replays it:

- Writer: the class-test inputs and a rich feature map (tricky numbers,
  NaN/inf texts, lists, subordinates by name and by number, MS1 isotopes,
  UIS target and decoy identification transitions, 64-bit unique and run
  ids, source-id mapping), written to a directory, appended as a second
  run, rejected as a duplicate run, to an empty feature map, to an archive,
  extended with `preserve_existing`, and replaced; every table compared
  cell by cell (schema, nulls, bits) with the Release files, except the
  documented stale cells, which must differ exactly as modelled;
  `metadata.json` byte for byte except `build_time`; archive entry sets and
  index keys; 10 writer error cases by kind and message, and which run
  directories stay behind.
- Reader: `load`, four peak-group variants (MS2, MS1+MS2, a main score, a
  missing main score), transition features and unscored data on the
  upstream fixture `OpenSwathWorkflow_tworuns_1_17.output.oswpq` (verbatim
  from core `bc9cc12` class-test data) and on every Release-written
  container, as directories and as archives, including sort ties and the
  unknown-reference errors; every returned column compared value by value.
- `OpenSwathOSWParquetWriter_test.cpp`, `OpenSwathOSWParquetReader_test.cpp`
  and `OpenSwathOSWParquetRoundTrip_test.cpp` transcribed.

Tier 1 for the exports: `../oracle/O5c-export` (same host, install and
flags, twice, identical). `tests/openswath_parquet_exporter.rs` writes the
class-test tables, empty tables, a table using every optional member with
dynamic columns and value vectors shorter and longer than the column lists
(NaN, infinity, -0.0, 1e-300, wrapping-sized ordinals), and a column-only
table, and compares every column (name, type, nullability, values, nulls)
with the Release file; the unwritable destination; both class-test sections
transcribed. `OpenSwathExportData.h` defaults and `toString` are unit tests
in `openswath_export_data.rs`.

Provenance: `tests/data/openswath_osw_parquet_provenance.json`,
`tests/data/openswath_parquet_exporter_provenance.json`.

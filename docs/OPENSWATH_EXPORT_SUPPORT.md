# OpenSWATH export data, configuration and writers

Rust ports of the user-facing OpenSWATH export classes of
`ANALYSIS/OPENSWATH` (pin `bc9cc12`):

| Header | Rust |
| --- | --- |
| `OpenSwathExportData.h` | [`src/analysis/openswath/openswath_export_data.rs`](../src/analysis/openswath/openswath_export_data.rs) |
| `OpenSwathExportConfig.h` | [`src/analysis/openswath/openswath_export_config.rs`](../src/analysis/openswath/openswath_export_config.rs) |
| `OpenSwathResultsExporter.h` | [`src/analysis/openswath/openswath_results_exporter.rs`](../src/analysis/openswath/openswath_results_exporter.rs) |
| `OpenSwathMatrixExporter.h` | [`src/analysis/openswath/openswath_matrix_exporter.rs`](../src/analysis/openswath/openswath_matrix_exporter.rs) |

Tests: [`tests/openswath_export.rs`](../tests/openswath_export.rs) (22 with `sqlite`, 23 without; the
upstream class tests and hand-derived cases),
[`tests/openswath_export_oracle.rs`](../tests/openswath_export_oracle.rs) (4,
executed comparison). Manifest:
[`tests/data/openswath_export_provenance.json`](../tests/data/openswath_export_provenance.json).

The exporters take *rows*; reading them out of an OSW or Parquet file is
`FORMAT/OSWFile.h` (`readOpenSwathExportRows`), which is not part of these
headers. The Parquet arm of both exporters is `ExportTable::write_parquet` over
`table()` behind the `arrow` feature (group O integration); see "Parquet" below.

## API mapping: `OpenSwathExportData.h`

| Source member | Rust |
| --- | --- |
| `enum class OpenSwathExportFileFormat { TSV, Parquet }` | `OpenSwathExportFileFormat::{Tsv, Parquet}` |
| `enum class OpenSwathIPFExportMode { Peptidoform, Augmented, Disable }` | `OpenSwathIPFExportMode` |
| `enum class OpenSwathMatrixLevel { Precursor, Peptide, Protein, Gene }` | `OpenSwathMatrixLevel` |
| `enum class OpenSwathMatrixNormalization { None, Median, MedianMedian }` | `OpenSwathMatrixNormalization` |
| `toString(OpenSwathExportFileFormat)`, `toString(OpenSwathIPFExportMode)`, `toString(OpenSwathMatrixLevel)`, `toString(OpenSwathMatrixNormalization)` | `as_str()` and `Display` of each enum |
| `struct OpenSwathExportRow` (59 members) | `OpenSwathExportRow`; `std::optional<T>` is `Option<T>`; `Int64` is `i64`, `Int32` `i32` |
| `struct OpenSwathFeatureScoreRow` (88 members) | `OpenSwathFeatureScoreRow`; `std::vector<double>` is `Vec<f64>` |
| `struct OpenSwathFeatureScoreTable` | `OpenSwathFeatureScoreTable` |
| `struct OpenSwathTransitionScoreRow` (20 members) | `OpenSwathTransitionScoreRow` |
| `struct OpenSwathTransitionScoreTable` | `OpenSwathTransitionScoreTable` |
| `struct OpenSwathQuantMatrix` | `OpenSwathQuantMatrix` |

Member names are the source's. A script compares the member list and order of
every struct with the header (`ALL OK`, 59 + 88 + 3 + 20 + 2 + 4 members).
Defaults are the source's: the id members are -1, everything else zero, `false`,
empty or absent.

## API mapping: `OpenSwathExportConfig.h`

| Source member | Rust |
| --- | --- |
| `OpenSwathExportFilterConfig` (14 members, PyProphet defaults) | `OpenSwathExportFilterConfig` with `Default` |
| `OpenSwathResultsExportConfig { filters, format }` | `OpenSwathResultsExportConfig` |
| `OpenSwathMatrixExportConfig { filters, level, normalization, top_n, consistent_top, format }` | `OpenSwathMatrixExportConfig`; `top_n` is `usize` (source `Size`) |
| `OpenSwathParquetExportConfig { filters, include_transition_data }` | `OpenSwathParquetExportConfig` |

## API mapping: exporters

| Source member | Rust |
| --- | --- |
| `OpenSwathResultsExporter::write(filename, rows, config)` | `OpenSwathResultsExporter::write`; the TSV text alone: `tsv_string`, `write_tsv`; the Parquet columns: `table` |
| `OpenSwathMatrixExporter::buildMatrix(rows, config)` | `OpenSwathMatrixExporter::build_matrix` |
| `OpenSwathMatrixExporter::writeMatrix(filename, matrix, config)` | `OpenSwathMatrixExporter::write_matrix`; `write_tsv`; `table` |
| anonymous-namespace helpers (`median_`, `buildSampleNames_`, `selectBestPeakGroups_`, `summarizePeptides_`, `validateMatrixDimensions_`, `applyNormalization_`) | private functions of the same names in snake case |

## Preserved source conventions

- **Number formatting.** TSV doubles are written with the default stream state
  (`os << double`: `%g` at precision 6, so `450.7`, `1.23457e+06`, `1e-05`);
  optional doubles through `StringUtils::toStr(double)` at full precision
  (`321.5`, `1.0e-05`, `3.6739e04`); integers in decimal; flags as `1` and `0`;
  an absent optional as the empty field. Both formatters are the crate's
  existing reproductions (`format::file_info::text_format`).
- **Matrix construction.** Sample columns are the runs in increasing run id;
  the label is the run name, `RUN_ID_<id>` when empty, and `_RUN<id>` appended
  when the label is already used by a run met earlier in row order. The best
  peak group of a (run, transition group) has the lowest `m_score`, then the
  highest `intensity`, then the lowest `feature_id`. Peptide level averages the
  `top_n` precursors with the highest median intensity over the runs
  (`consistent_top`) or the `top_n` most intense peak groups of each run; protein
  and gene level average the peptide summaries per sample, over the `top_n`
  peptides of highest median when `consistent_top`. Rows with an empty protein
  or gene name are skipped at that level. Ties between precursors of equal
  median intensity are resolved in the order of the Release build's
  `std::sort` (`math::source_sort`), as the source's selection depends on it.
- **Checks.** No rows, `top_n` of zero above precursor level, no row surviving
  the peptide summary, and no protein or gene annotation are
  `Exception::Precondition` in the source (`Error::InvalidValue` here); so are
  an empty matrix and a matrix whose row widths differ from its column names at
  write time, which are refused before anything is written.
- **Normalization.** Per sample column: the median of the present values;
  `Median` divides each value by it (a column with a zero or non-finite median
  is left as it is), `MedianMedian` then multiplies by the median of the finite
  column medians.

## Native differences

| Source | Here | Why |
| --- | --- | --- |
| Parquet output built with Arrow | `ExportTable::to_record_batch` and `write_parquet` (feature `arrow`; `ArrowIOHelpers::writeTableToParquet` is `format::arrow_io_helpers::write_table_to_parquet` with its default configuration); without the feature `write` / `write_matrix` return `Error::Unsupported` and create nothing | the exporters build the typed columns themselves (`table()`) and hand the batch to the shared writer; every field is nullable as `arrow::field(name, type)` makes it |
| `FileNotWritable` | `Error::Io` | the crate's error set |
| Parquet-first column helper `ParquetFile::appendOrThrow` | arrays are built from the typed columns (`Vec<Option<T>>`), which cannot fail | no Arrow builders are involved |
| partial file on a failed write | same | not atomic in the source either |

## Checked boundaries and evidence

**Tier 1, executed differential.** `tests/data/openswath_export/` holds the
output of the C++ Release build for 34 matrix scenarios (every level,
normalization and selection mode, ties, label collisions, sparse and
unannotated data), one results export of 40 rows over all 59 columns, the
row preparation of the OSW writer and a grid of scheduler plans. They were
produced by `../oracle/O5a-export/driver.cpp` on ibminode06 (Linux x86_64,
`openms4-release-bc9cc12-c19e494-174b576`, compiler and flags in the manifest),
run twice with identical output. The replay compares text byte for byte:
34 of 34 matrix scenarios and the results table are equal.

The Parquet arm has its own tier-1 run, `../oracle/integ-O-parquet/driver.cpp`
(the same case file; ibminode06, the same Release install, run twice, the 35
Parquet files byte-identical between runs): `OpenSwathMatrixExporter::writeMatrix`
and `OpenSwathResultsExporter::write` with format Parquet for every scenario,
retained as `tests/data/openswath_export/parquet/*.parquet`. The replay
(`matrix_exporter_parquet_arm_matches_the_release_build`,
`results_exporter_parquet_arm_matches_the_release_build`, feature `arrow`)
writes each scenario through this crate, reads both files back with
`ParquetFile::read_table` and compares the tables: schema (names, Arrow types,
nullability) and every value. 34 of 34 matrix files and the 59-column results
file are equal.

**Tier 3, source review.** `OpenSwathResultsExporter_test.cpp`,
`OpenSwathMatrixExporter_test.cpp` and the matrix and shape sections of
`OpenSwathExport_test.cpp` are transcribed; the Parquet reading sections of the
upstream tests are transcribed as `results_exporter_writes_the_parquet_of_the_upstream_test`
(feature `arrow`): the two-row export read back column by column with
`ParquetFile`, nulls and 64-bit values included, and the empty export's 59
columns.

**Tier 4.** Hand-derived peptide, protein and gene summaries and the
normalization; the label rules; the number formats.

Upstream tests not carried: `OpenSwathExport_test.cpp` sections that read
`PyProphet_inference_test.osw` through `OSWFile` and write Parquet through
Arrow (`OSWFile`, `OpenSwathParquetExporter`, `ParquetFile`).

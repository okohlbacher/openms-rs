# OpenSWATH OSW writer

Rust port of `ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h` and its translation
unit (pin `bc9cc12`):
[`src/analysis/openswath/openswath_osw_writer.rs`](../src/analysis/openswath/openswath_osw_writer.rs).

The writer turns a `FeatureMap` into buffered rows of typed cells and inserts
them into the tables `RUN`, `FEATURE`, `FEATURE_MS1`, `FEATURE_MS2`,
`FEATURE_PRECURSOR` and `FEATURE_TRANSITION` of an SQLite file, the input of
PyProphet. The row preparation is pure Rust and always present; the members
that touch the database need the `sqlite` feature and return
`Error::Unsupported` without it (so the workflow compiles either way).

Tests: [`tests/openswath_osw_writer.rs`](../tests/openswath_osw_writer.rs)
(7, `sqlite`), the row preparation in
[`tests/openswath_export.rs`](../tests/openswath_export.rs) and
[`tests/openswath_export_oracle.rs`](../tests/openswath_export_oracle.rs).
Manifest: [`tests/data/openswath_export_provenance.json`](../tests/data/openswath_export_provenance.json).

## API mapping

| Source member | Rust |
| --- | --- |
| `OpenSwathOSWWriter(output_filename, uis_scores)` | `OpenSwathOSWWriter::new` |
| `isActive` | `is_active` |
| `writeHeader` | `write_header` (`sqlite`) |
| `addRun(run_id, input_filename)` | `add_run` (`sqlite`) |
| `setRunId` | `set_run_id` (the sign bit is cleared) |
| `getScore(feature, name)` | `get_score` |
| `getSeparateScore(feature, name)` | `get_separate_score` |
| `prepareRows(pep, transition, output, id)` | `prepare_rows` |
| `prepareRowsInto(rows, pep, transition, output, id)` | `prepare_rows_into` |
| `prepareLine(pep, transition, output, id)` | `prepare_line` |
| `writeRows(osw_output)` | `write_rows` (`sqlite`) |
| `writeLines(to_osw_output)` | `write_lines` (`sqlite`) |
| `OSWValue` (`Type`, `Storage`, the constructors, `null`, `text`, `isNull`, `asInt`, `asDouble`, `asText`) | `OswValue` (enum `Null`, `Int64`, `Double`, `Text`), `OswValue::{null, text, real, parse, is_null, as_int, as_double, as_text}`, `From<i64/u64/i32/f64>`; the string constructor is `parse` |
| `FEATURE_COLUMN_COUNT` .. `FEATURE_TRANSITION_COLUMN_COUNT` | the constants of the same names |
| `FeatureRow` .. `FeatureTransitionRow` | the type aliases of the same names |
| `OSWData` (`reserve` x2, `append`, `empty`, `rowCount`, `estimateMemoryUsage`) | `OswData` (`reserve`, `reserve_rows`, `append`, `is_empty`, `row_count`, `estimate_memory_usage`) |
| the column lists, the `CREATE TABLE` statements | `FEATURE_COLUMNS` .. `FEATURE_TRANSITION_COLUMNS`, `CREATE_SQL`, `insert_statement` |
| `SqliteHelper::clearSignBit` | `clear_sign_bit` |
| `SqliteConnector` persistent connection, `conn_mutex_` | a `Mutex<Option<SqliteConnector>>` opened on the first write |

## Preserved source conventions

- **Cell classification.** A text is classified as the source's string
  constructor does: empty, `NULL`, `nan`, `-nan`, `NaN`, `-NaN` are NULL; text
  `strtoll` consumes completely (base 10, sign, surrounding white space) is an
  integer; text `strtod` consumes completely a real (hexadecimal, `inf`,
  `infinity`, `nan(...)` included, a NaN being NULL); anything else, and any
  number that over- or underflows, is text. So a numeric transition group id
  is an integer cell and `PRECURSOR_ID` of `NIAGGLEK_2` a text cell.
- **Metavalues.** An integer metavalue is an integer cell, a real a real cell
  (NaN as NULL), a string is classified as above, and a list is classified as
  its `toString`. Score lists (`masserror_ppm`, the `id_target_*` and
  `id_decoy_*` lists) yield one cell per element; an index outside the list,
  or an undefined score, is NULL. The `masserror_ppm` of an MS2 transition row
  is indexed by the position of the subordinate among the MS2 subordinates.
- **Transition rows.** Without UIS scoring one row per MS2 subordinate. With
  UIS scoring the rows are those of the target then the decoy lists, as many as
  `id_*_num_transitions` says but no more than there are names; the MS2
  subordinate rows are used only when no UIS row exists at all. The peak shape
  columns of UIS rows are filled only when the *target* `start_position_at_5`
  list is non-empty and its first element does not print as `0` (a real zero
  prints `0.0`, so it enables them, an integer zero disables them), and that
  decision also applies to the decoy rows. The `VAR_MASSDEV_SCORE` list is
  written to both `MASSERROR_PPM` (index 7) and `VAR_MASSDEV_SCORE` (index 15)
  of a UIS row, as the source does.
- **Precursor rows.** A subordinate with `FeatureLevel` `MS1` and a positive
  intensity; `ISOTOPE` is the text after `Precursor_i` of its `native_id`,
  NULL when absent.
- **`RUN.ID` is an INTEGER** (bound as int64, sign bit cleared), so the
  `FEATURE.RUN_ID = RUN.ID` join works.
- **Write mode.** The connection is opened on the first write with
  `journal_mode=WAL` and `synchronous=OFF` (a failure of either is ignored);
  each `write_rows` / `write_lines` is one transaction, rolled back on error.
  `prepare_line` produces the literal `INSERT` statements in the order
  `FEATURE`, `FEATURE_MS1`, `FEATURE_PRECURSOR`, `FEATURE_MS2`,
  `FEATURE_TRANSITION`; a text cell is quoted with doubled single quotes and a
  real written through `StringUtils::toStr(double)`.

## Native differences

| Source | Here | Why |
| --- | --- | --- |
| `asInt`/`asDouble`/`asText` of another type (undefined) | `Option` | no undefined reads |
| `OSWValue` union with a manually managed string | an enum | memory safety |
| `estimateMemoryUsage` row size (`sizeof(OSWValue)` 40) | the size of this port's cell (32) plus the source's 16 bytes of headroom per cell | the budget is a heuristic; only its units matter |
| a failed conversion of `norm_RT`, `delta_rt` or `id_*_num_transitions` (`ConversionError`) | `Error::InvalidValue` | the crate's error set |
| the caller guards `writeRows` with an OpenMP critical section | the writer holds a lock | no caller discipline needed |
| hexadecimal significands of more than 16 digits | rounded twice | `strtod` rounds once; no identifier column holds such text |
| `executeStatement` failure type (`IllegalArgument`) | `Error::Io` | the crate's error set |

## Checked boundaries and evidence

**Tier 1.** The five OSW scenarios of the oracle (`osw.out`: the upstream
masserror case with `prepareLine`, UIS lists with peak shape metrics and decoys,
UIS with a zero start position and the MS2 fallback, MS1 scores with precursor
rows and several native id forms, and a table of scores and 39 text
classification cases) are replayed bit for bit: every cell type and value, as
the hex of each double, and the SQL text. See
`tests/data/openswath_export_provenance.json`.

**Tier 3.** The three database sections of `OpenSwathOSWWriter_test.cpp`
(`RUN.ID` INTEGER, masserror indexing, UIS name count) are replayed against a
real SQLite file, reading back the values the upstream test queries.

**Tier 4.** The schema (`PRAGMA table_info` against the column lists), the
transaction rollback on a duplicate key, the WAL mode, the NULL and text cell
binding.

Not carried: `OpenSwathOSWParquetWriter`, `OpenSwathOSWParquetReader` and
`OSWFile` are separate headers.

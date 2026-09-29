# TransitionPQPFile support

Ports `ANALYSIS/OPENSWATH/TransitionPQPFile.h` (with its `.cpp`) at
OpenMS4-core `bc9cc12`, TraML-structure (heavy) path. Rust:
`src/analysis/openswath/transition_pqp_file.rs`, behind the existing `sqlite`
feature (`rusqlite =0.40.2`, already pinned; no new dependency). SQLite is
reached through `format::sqlite_connector::SqliteConnector`.

## API mapping

| C++ member | Rust | Notes |
|---|---|---|
| `TransitionPQPFile()` / dtor, base `TransitionTSVFile` | `TransitionPQPFile { tsv: TransitionTSVFile }`, `new()`, `Default` | the base becomes a field |
| inherited parameters | `defaults()`, `to_param()`, `from_param()`, `from_param_with_warnings()` | those of `TransitionTSVFile` |
| `PQPSqlQueryInfo` | `PqpSqlQueryInfo` | all four flags and the query |
| `buildPQPSelectQuery_(conn, legacy, stable_order)` | `TransitionPQPFile::build_pqp_select_query(&SqliteConnector, legacy, stable_order)` | identical SQL text |
| `readPQPInput_` | `read_pqp_input(path, legacy) -> Vec<TsvTransition>` | |
| `writePQPOutput_` | private `write_pqp_output` | |
| `convertTargetedExperimentToPQP` | `convert_targeted_experiment_to_pqp(path, &exp)` | |
| `convertPQPToTargetedExperiment` (TraML structure) | `convert_pqp_to_targeted_experiment(path, legacy) -> (TargetedExperiment, warnings)` | |
| `getPQPIDToTraMLIDMap` | `TransitionPQPFile::pqp_id_to_traml_id_map(path, table)` | `HashMap`, as the source's `unordered_map` |
| `getPQPCurrentIDToTraMLIDMap` | `TransitionPQPFile::pqp_current_id_to_traml_id_map(path, table)` | |
| inherited `validateTargetedExperiment` | `validate_targeted_experiment` | |
| `convertPQPToTargetedExperiment` (`LightTargetedExperiment`), `streamPQPToLightTargetedExperiment_`, both `convertLightTargetedExperimentToPQP` | **not ported** | need `OPENSWATHALGO/DATAACCESS/TransitionExperiment.h` and `ANALYSIS/OPENSWATH/OpenSwathLibraryIDNormalizer.h` (both unmapped) |
| anonymous-namespace `PreparedInsert`, `tableHasRows_`, `validateTraMLIDTableName_`, `configurePQPReadConnection_` | rusqlite prepared statements, private `has_rows`, `validate_traml_id_table_name`; pragmas not set | the cache/mmap/`query_only`/`threads` pragmas only tune performance |

## Preserved source conventions

- Schema (13 tables, `VERSION` 3), the select statement text including the
  `GENE` and peptidoform aggregates, the stable `ORDER BY` of the heavy
  reader, reading until a `NULL` precursor m/z, `NA` gene names read as empty.
- `sqlite3_column_*` conversions: an integer id reads as its decimal text;
  charges are only taken from integer columns (`extractValueIntStr`).
- Writer: ids are the ranks of the sorted, deduplicated peptide (UniMod
  string, plus IPF peptidoforms), compound, protein and group id lists; gene
  ids in first-seen order; a precursor's m/z and decoy state from its first
  (detecting) transition; decoy flags propagated by the four `UPDATE`s;
  compound names from `CompoundName` or the id; an existing file is removed
  first (`TransitionPQPFile.cpp:653`).
- An unset peptide charge is written as `0`, the value the source's
  release-mode `getChargeState` returns (`TransitionPQPFile.cpp:1079`); a
  compound without charge writes `NULL`.

## Native differences

- Everything that can fail on the data (invalid references, missing retention
  times, unbuildable sequences) is checked before the file is removed or
  created, so a refused experiment leaves an existing file untouched. A SQLite
  error midway still leaves a partial file, as upstream.
- The id-map functions open the file read-only; the source opens it
  read-write-or-create (`TransitionPQPFile.cpp:1118`, `:1151`), creating an
  empty database for a missing path.
- The mapping tables are written with prepared statements instead of batched
  SQL text; the rows and their order are the same.
- `ReadLimits::max_transitions` of the TSV parameters bounds the row count.

## Checked boundaries and evidence

`tests/transition_pqp_file.rs` (7 tests; `--features
semantic-validation,sqlite --test transition_pqp_file` on dax: 7 passed):

- Tier 1 against retained C++ outputs at test-data `0cb15f2`
  (`topp/CMakeLists.txt:1762-1851`), via the `TargetedFileConverter` heavy path
  (openswath `e29e7ee`, `TargetedFileConverter.cpp:229-276`):
  TargetedFileConverter 3 (TraML to PQP to TraML, numeric ids) and 4 (legacy
  ids) reproduce `_3_output`/`_4_output`; 11 and 12 (proteomics, 12 with
  `override_group_label_check`) reproduce their inputs with `id=` lines
  skipped as the upstream whitelist does; 10 reads the **C++-written**
  `TargetedFileConverter_10_input.pqp` (legacy ids) and reproduces
  `_10_output.tsv`.
- Tier 3: `TransitionPQPFile_test.cpp:71-137` (empty GENE table, issue #8687)
  and `:163-206` (two genes do not duplicate a transition), heavy halves.
- Tier 4: id maps over both tables, invalid table names refused (the source
  checks the same at `TransitionPQPFile_test.cpp:313-316`), invalid references
  refused with no file created.

## C++ defect candidates

1. `TransitionPQPFile.cpp:1079` binds an unset peptide charge as `0` (a
   release build skips the precondition) instead of `NULL`, unlike compounds.
2. `TransitionPQPFile.cpp:1118`/`:1151`: read-only queries open the database
   in create mode.

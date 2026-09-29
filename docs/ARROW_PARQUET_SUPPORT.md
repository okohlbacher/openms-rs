# Arrow and Parquet support

Lane F1a (wave 4) ports the core of OpenMS's Apache Arrow / Parquet family:
the table schemas, the shared read/write helpers, the filter expressions and
the protein-identification tables. Everything Arrow-typed is behind the
optional `arrow` feature (see the `arrow-rs` row of
[THIRD_PARTY_CRATE_DECISIONS](THIRD_PARTY_CRATE_DECISIONS.md));
`ParquetFilter` needs no Arrow type and is in the default build.

| Source header (core `bc9cc12`) | Rust | Status |
|---|---|---|
| `FORMAT/ArrowSchemaRegistry.h` | `src/format/arrow_schema_registry.rs` | complete |
| `FORMAT/ParquetFile.h` | `src/format/parquet_file.rs` | complete |
| `FORMAT/ParquetFilter.h` | `src/format/parquet_filter.rs` | complete |
| `FORMAT/ArrowIOHelpers.h` (+ `ParquetWriteConfig` of `MSExperimentArrowExport.h`) | `src/format/arrow_io_helpers.rs` | complete (QPX value validation by lane F1c, see [QPX_SUPPORT](QPX_SUPPORT.md)) |
| `FORMAT/ProteinIdentificationArrowIO.h` | `src/format/protein_identification_arrow_io.rs` | complete |
| `FORMAT/PSMArrowIO.h` | — | not ported (`QPXFile`, which it builds on, is ported: [QPX_SUPPORT](QPX_SUPPORT.md)) |
| `FORMAT/FeatureMapArrowIO.h` | `src/format/feature_map_arrow_io.rs` | partial: features table; PSM table and bundle need `QPXFile` |
| `FORMAT/ConsensusMapArrowIO.h` | — | not ported: needs `QPXFile`, ProForma parsing of peptidoforms |

## Representation

* A source `std::shared_ptr<arrow::Table>` is an `arrow_array::RecordBatch`
  holding every row: one array per column. The source's readers combine
  chunks anyway; a caller with several batches concatenates them with
  `arrow_select::concat::concat_batches`. A multi-chunk table therefore
  cannot reach `ParquetFile::getColumn`, whose source returns only
  `chunk(0)` (measured: 2 of 5 rows) — that defect cannot occur here.
* A possibly-null `std::shared_ptr<arrow::Array>` is an `Option<&dyn Array>`,
  a row index a `usize`. The source casts arrays to the type a getter expects
  without checking (`std::static_pointer_cast`, undefined behaviour for any
  other type) and reads rows without bounds checks; every getter here checks
  both.
* `arrow::KeyValueMetadata` is `arrow_io_helpers::KeyValueMetadata`, ordered
  and allowing repeated keys, as the source's type. An arrow-rs `Schema`
  holds its metadata in a map; converting keeps the first value of a
  repeated key.
* Functions that return `bool`/`nullptr` and log in the source return a
  `Result`; where the source logs a warning and continues, the port returns
  the warning lines.

## API mapping

### `ArrowSchemaRegistry.h`

| Source | Rust |
|---|---|
| `ArrowSchemaValidation::Mode` | `validation::Mode` |
| `ArrowSchemaValidation::ValidationResult` (+ `toString`) | `validation::ValidationResult` (+ `Display`) |
| `ArrowSchemaValidation::validate(table, expected, mode)` | `validation::validate`; `validation::validate_schema` for a schema alone (native) |
| the 22 schema structs (`ProteinSchema` ... `XIPMSchema`): column-name constants, `*Type()` helpers, `schema()` | the same-named structs: associated `&'static str` constants, `*_type()` and `schema()` |
| `arrow::DataType::Equals` as validation uses it | `validation::types_equal` (native helper) |
| `DataType::ToString` / `Field::ToString` in messages | `validation::cpp_type_string` / `cpp_field_string` (native helpers) |

### `ParquetFile.h`

| Source | Rust |
|---|---|
| `appendOrThrow(status, column)` | `ParquetFile::append_or_error` |
| `finishArray<BuilderT>(builder, name)` | `ParquetFile::finish_array` |
| `writeTable(table, filename, row_group_size = 262144)` | `ParquetFile::write_table`, `write_table_with_row_group_size` |
| `readTable(filename)` | `ParquetFile::read_table` (+ `read_table_with_limits`) |
| `readTable(RandomAccessFile)` | `ParquetFile::read_table_from` (+ `_with_limits`) |
| `getColumn(table, name)` / `getOptionalColumn` | `ParquetFile::get_column` / `get_optional_column` |
| `getInt64`, `getDouble`, `getBool`, `getString`, `getStringList` | `ParquetFile::get_int64`, `get_double`, `get_bool`, `get_string`, `get_string_list` |
| `jsonEscape(input)` | `ParquetFile::json_escape` |
| `rowCount(filename)` | `ParquetFile::row_count` |
| `throw_finish_error_` (private) | folded into `finish_array` |

### `ParquetFilter.h`

| Source | Rust |
|---|---|
| `ColumnType`, `Condition`, `FilterExpression` | same-named types |
| `ParquetFilter` constructors, `andNext`/`orNext`, `eq`/`ne`/`lt`/`le`/`gt`/`ge`/`in` for `Int64` and strings | `ParquetFilter::new`, `and_next`/`or_next`, `eq_int` ... `in_int`, `eq_str` ... `in_str` |
| `getExpression`, `empty`, `merge` | `expression`, `is_empty`, `merge` |
| `ParquetFilterBuilder` and its forwarding members, `build` | `ParquetFilterBuilder` (same names), `filter` |

### `ArrowIOHelpers.h` and `ParquetWriteConfig`

| Source | Rust |
|---|---|
| `ParquetWriteConfig` (+ `Compression`) | `ParquetWriteConfig`, `ParquetCompression` |
| `generateUuidV4()` | `generate_uuid_v4` |
| `tableRowCount(table)` | `table_row_count` |
| `writeTableToParquet(table, filename, config)` | `write_table_to_parquet` |
| `concatenateAndWriteToParquet(tables, filename, config, metadata)` | `concatenate_and_write_to_parquet` |
| `qpxFileMetadata(file_type, config, extra)` | `qpx_file_metadata` → `KeyValueMetadata` |
| `qpxScanFormat(native_id)` / `qpxScanFormat(native_ids)` | `qpx_scan_format` / `qpx_scan_format_of_all` (returns the warning) |
| `qpxRunFileName(path)` | `qpx_run_file_name` |
| `qpxScanComponents(reference)` | `qpx_scan_components` |
| `QPXRunFileNameKeys` | `QpxRunFileNameKeys` |
| `qpxPsmRunFileName(hit, identification, resolved, keys)` | `qpx_psm_run_file_name` |
| `qpxWarnOnRunNameCollisions(context, paths)` | `qpx_warn_on_run_name_collisions` (+ `qpx_run_name_collisions`, the warning lines) |
| `qpxIntensityLabel(label, channel)` / `qpxIntensityLabels(cmap)` | `qpx_intensity_label` / `qpx_intensity_labels` |
| `qpxIsCanonicalIntensityLabel(label)` | `qpx_is_canonical_intensity_label` |
| `getColumn(table, name, required)` | `get_column` |
| `getStringValue`, `getDoubleValue`, `getFloatValue`, `getInt32Value`, `getInt64Value`, `getBoolValue`, `isNull` | `get_string_value`, `get_double_value`, `get_float_value`, `get_int32_value`, `get_int64_value`, `get_bool_value`, `is_null` |
| `readMetaValues(array, row, target, excluded_keys)` | `read_meta_values` |
| `validateQPXValuesForWrite` (anonymous) | QPX view tables refused (see gaps) |

### `ProteinIdentificationArrowIO.h`

| Source | Rust (`ProteinIdentificationArrowIO::`) |
|---|---|
| `exportProteinsToArrow` / `exportProteinsToParquet` | `export_proteins_to_arrow` / `export_proteins_to_parquet` |
| `exportProteinGroupsToArrow` / `...ToParquet` | `export_protein_groups_to_arrow` / `..._to_parquet` |
| `exportSearchParamsToArrow` / `...ToParquet` (+ `definitions_by_run`) | `export_search_params_to_arrow` / `..._to_parquet` |
| `importFromParquet` | `import_from_parquet` (+ `_with_registry`) |
| `importSearchParamsFromArrow` / `...FromParquet` | `import_search_params_from_arrow` (+ `_with_registry`) / `..._from_parquet` |
| `importProteinsFromArrow` / `...FromParquet` | `import_proteins_from_arrow` (+ `_with_registry`) / `..._from_parquet` |
| `importProteinGroupsFromArrow` / `...FromParquet` | `import_protein_groups_from_arrow` / `..._from_parquet` |
| `synthesizeRunIdentifiers` | `synthesize_run_identifiers` (also returns the warnings) |
| `applyRunIdentifierRename` | `apply_run_identifier_rename` |
| `checkUniqueIdentifiers` | `check_unique_identifiers` |
| `writeArrowTableToParquet_` (anonymous) | `write_table` |

### `FeatureMapArrowIO.h`

| Source | Rust (`FeatureMapArrowIO::`) |
|---|---|
| `exportFeaturesToArrow(feature_map)` | `export_features_to_arrow` |
| `importFeaturesFromArrow(table, feature_map)` | `import_features_from_arrow` (returns the warnings) |
| `exportPSMsToArrow`, `importPSMsFromArrow` | not ported: `QPXFile::exportToArrow` / `importFromArrow` |
| `exportToParquet`, `importFromParquet` | not ported: the PSM table above, and the `DataProcessing`/metadata JSON of `features.parquet`'s file-level metadata |

## Preserved source conventions

* **Parquet writing.** `ParquetFile::writeTable` stores no Arrow schema;
  `writeTableToParquet` and the protein exports store it. Both write lists in
  the Parquet-compliant layout (element `element`), a `timestamp[s]` column
  as milliseconds (Parquet has no seconds; both source readers return
  `timestamp[ms]`, measured), and a table with no rows as one empty row group
  whose column chunk has no statistics (measured). The row group size is a row
  count, capped at 1,048,576 rows (measured: 1,100,000 rows at the default
  become 1,048,576 + 51,424), whatever `ParquetWriteConfig` documents about
  bytes. LZ4 is Parquet's `LZ4_RAW`, as Arrow C++ writes it. GZIP levels -1
  and 0..=9 succeed and 10 fails; ZSTD accepts every level (measured).
* **Parquet reading** combines chunks, returns `timestamp[ms]` for a stored
  `timestamp[s]`, restores `large_string` from the stored schema, and keeps
  the schema metadata (`ARROW:schema` removed), all as the source's reader
  does (measured on the C++-written files).
* **`readMetaValues`** converts `int` with `std::stoll` and `double`/`float`
  with glibc `strtod` semantics: prefix parsing, hexadecimal floats,
  `inf`/`infinity`/`nan(...)`, and `ERANGE` for overflow and for an inexact
  result below the smallest normal double (82 cases measured, including
  `1e-310` and `0x1.8p-1074` falling back to strings while `0x1p-1074` and
  `0e-400` convert). Lists strip one surrounding `[` `]` and use
  `ListUtils::create`; string items are trimmed of space, tab, newline and
  carriage return but not vertical tab.
* **QPX naming.** `qpxScanComponents` follows `SpectrumNativeIDParser`
  (the port reuses `percolator_infile::ScanNumberPattern`), including the
  last-match rule and the `Int32` range. The isobaric channel vocabulary is
  transcribed from the executed `IsobaricQuantitationMethod` registry (TMT 6
  to 35-plex, iTRAQ 4/8-plex).
* **Features table.** Features are flattened depth first, each followed by
  its subordinates, with the parent's id and the depth; `width` 0 is null, the
  bounding box is the union of the hulls' boxes, `FWHM` is left out of
  `metavalues` and restored from `width` by the import's `setWidth`
  (measured: a feature exported with `FWHM` 2.5 and width 3.5 reads back with
  `FWHM` 3.5). The import attaches subordinates deepest first and turns a row
  whose parent is missing into a top-level feature with the source's warning
  (measured on a sliced table).
* **Protein tables.** Every column's null rule is the source's (coverage
  unknown, empty sequence/description/version/database fields, unknown
  enzyme or specificity, pre-1970 dates → null). Ranks and missed cleavages
  pass through `int32` with the source's wrapping casts. Group indices are
  0-based per run and group type and are not read back. The import creates a
  run for a hit or group whose run identifier has no search parameters, and
  the hit form warns, as the source does.

## Native differences

* **Atomic writes.** The source opens and truncates the destination before
  encoding and removes the fragment on failure, which also destroys a file
  that existed there; the port encodes in memory, writes a sibling temporary
  file and renames it, so a failure leaves the destination unchanged.
* **Deterministic bytes.** Schema metadata is written as footer key-value
  pairs in key order and left out of the stored Arrow schema; both readers
  merge the pairs back. arrow-rs keeps schema metadata in a randomly seeded
  hash map, and encoding it inside the stored schema made the bytes vary
  between runs. `tests/arrow_parquet.rs` requires byte-identical output for
  equal tables.
* **Metadata entry order.** The source lists `metavalues` entries in
  `MetaInfoRegistry` order (the order the process first registered each
  name); the port lists them by key. Import is order-independent.
* **ZSTD levels below 1** have no arrow-rs counterpart and are written at
  level 1 (0 means zstd's default 3); the file decodes to the same table.
* **`QPXRunFileNameKeys`.** The source resolves registry indices and silently
  falls back to the resolved run path when the names were not registered
  before the keys were built (measured); the port has no registry and always
  consults the names.
* **Run paths and inference engine.** The source keeps the primary MS run
  paths as `spectra_data` metadata; the port keeps them as
  `ProteinIdentification::primary_ms_run_paths` and writes them into
  `metavalues` under `spectra_data` when non-empty. A path list the source
  sets explicitly empty (present, zero entries) cannot be told from an unset
  one here.
* **Null inputs.** The source's null tables and arrays are `Option`s where a
  function can meaningfully receive none; a `RecordBatch` reference cannot be
  null, so `writeTableToParquet(nullptr)` has no counterpart.
* **Protein modifications** that are mass tags (not named
  `ResidueModification`s) are refused on export (`Error::Unsupported`); the
  source's `ProteinHit` cannot hold them.
* **Features table.** Export refuses subordinates nested deeper than
  `Feature::MAX_SUBORDINATE_DEPTH` (128); the source recurses without a bound.
  A hull with no points is left out of the bounding box (the source's union
  with an empty `DBoundingBox` was not measured), and a hull point that is
  not finite is refused on import.
* **Registries.** Modification definitions read from `sp_metavalues` are
  registered into a caller-owned `ModificationsDB` (`*_with_registry`), not a
  process-wide database.

## Gaps (named dependencies)

* (Closed by lane F1c.) `writeTableToParquet` validates a QPX `psm_file`,
  `feature_file` or `pg_file` table with `QPXValueValidation` and refuses one
  that violates the QPX value contract (`Error::InvalidValue`), as the source
  does.
* `PSMArrowIO`, the PSM and bundle half of `FeatureMapArrowIO`, and
  `ConsensusMapArrowIO` build and read
  their PSM tables through `QPXFile::exportToArrow`/`importFromArrow`
  (`FORMAT/QPXFile.h`, now `QpxFile::export_to_arrow`/`import_from_arrow`),
  reconstruct peptidoforms with ProForma, and store `DataProcessing` as JSON;
  their feature and consensus-feature tables need the kernel metadata JSON
  helpers of those files. They are backlog items; `QPXFile` no longer blocks
  them.

## Evidence

Tier 1, executed differential against the Linux x86_64 Release install
`openms4-release-bc9cc12-c19e494-174b576` (g++ 14.4 `-O2 -DNDEBUG`, Arrow C++
23.0.1) on ibminode06, each driver run twice with identical records. Drivers
live outside the repository (`../oracle/F1a-arrow`, `F1a-arrow2`,
`F1a-protein`, `F1a-features`); their sha256s and the record sha256s are in
`tests/data/arrow_parquet_provenance.json`.

* `tests/arrow_parquet.rs` (23 tests): every column of the 22 schemas and the
  32 nested-type helpers, type equality, 20 validation cases, the
  `ParquetFile` getters and I/O, `writeTableToParquet` codecs, levels, row
  groups, QPX gate and metadata, concatenation, QPX metadata, the 82
  `readMetaValues` cases, the QPX naming and label helpers, nested columns,
  and byte determinism. Where a C++-written file exists, the port reads it
  and its own output with the same checks.
* `tests/protein_identification_arrow_io.rs` (6 tests): the protein tables of
  six retained C++-written bundles (test-data `0cb15f2`) import field by field
  as the executed source imports them; every table of those bundles validates
  strictly against the registry (a feature or consensus bundle's PSM table
  with its owner-id column in front of `PSMSchema`); the port's exports of a constructed set of
  runs equal the C++-written tables cell by cell; both round trips
  reconstruct the source's runs; identifier synthesis, renaming and the
  uniqueness check.
* `tests/feature_map_arrow_io.rs` (4 tests): the retained C++-written
  features table of `FileConverter_9_output.featureparquet` imports feature
  by feature as the executed source imports it; the port's features table of
  a constructed map with a three-level hierarchy, hulls, `FWHM` and a
  `UInt64` id above `i64::MAX` equals the C++-written one; both read back as
  the source reads them; orphaned subordinates and an empty table.
* `tests/parquet_filter.rs` (19 tests, default build): every section of the
  upstream class test `ParquetFilter_test.cpp` (tier 3).

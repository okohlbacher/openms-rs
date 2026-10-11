# QPX support

Lane F1c (wave 4) ports the QPX (Quantitative Proteomics eXchange, quantms)
layer of OpenMS's Arrow/Parquet family: the opaque row identities, the value
contract every QPX view is checked against before it is written, the PSM
exporter and importer, and the scope guard of a whole collection. Everything is
behind the optional `arrow` feature, like the rest of the family (see
[ARROW_PARQUET_SUPPORT](ARROW_PARQUET_SUPPORT.md)). The source marks all four
classes experimental.

| Source header (core `bc9cc12`) | Rust | Status |
|---|---|---|
| `FORMAT/QPXIdentity.h` | `src/format/qpx_identity.rs` | complete |
| `FORMAT/QPXValueValidation.h` | `src/format/qpx_value_validation.rs` | complete |
| `FORMAT/QPXFile.h` | `src/format/qpx_file.rs` | complete |
| `FORMAT/QPXCollectionExport.h` | `src/format/qpx_collection_export.rs` | partial: `requireExportable` |

Porting `QPXValueValidation` also closes the last gap of
`FORMAT/ArrowIOHelpers.h`: `write_table_to_parquet` now validates a QPX
`psm_file`, `feature_file` or `pg_file` table instead of refusing it.

## API mapping

### `QPXIdentity` (namespace)

| Source | Rust |
|---|---|
| `FeatureLinks` (`unordered_map<Int64, Int64>`) | `qpx_identity::FeatureLinks` (`HashMap<i64, i64>`) |
| `Null`, `Float32`, `Value` (`std::variant`) | `qpx_identity::Value` (`Null`, `String`, `Int64`, `Float32`, `Int32List`, `StringList`) |
| `formatFloat(float)` | `format_float(f32)` |
| `canonical(values, unordered_list_indices)` | `canonical(&[Value], &[usize])` |
| `deriveId(values, unordered_list_indices)` | `derive_id(&[Value], &[usize])` |
| `FEATURE_COMPOSITE`, `PSM_COMPOSITE`, `PG_COMPOSITE` | the same constants |
| `featureId(run, peptidoform, charge, rt, scan, observed_mz)` | `feature_id(&str, &str, i64, Option<f32>, &[i32], f32)` |
| `psmId(run, scan, peptidoform, charge)` | `psm_id(&str, &[i32], &str, i64)` |
| `pgId(pg_accessions, grouped_runs, label)` | `pg_id(&[A], &[R], Option<&str>)` |

### `QPXValueValidation`

| Source | Rust |
|---|---|
| `View::{PSM, FEATURE, PROTEIN_GROUP}` | `View::{Psm, Feature, ProteinGroup}` |
| `Result` (`valid`, `errors`) | `ValidationResult` (same fields, default valid) |
| `Result::toString()` | `Display` (errors joined by `"; "`) |
| `QPXValueValidation(View)`, destructor, move constructor/assignment | `QpxValueValidation::new(View)`; drop and move are Rust's |
| deleted copy constructor/assignment | not `Clone` |
| `validate(const shared_ptr<Table>&)` | `validate(&mut self, &RecordBatch)` |
| `reset()` | `reset(&mut self)` |
| — | `view()` (native accessor) |

### `QPXFile`

| Source | Rust |
|---|---|
| `exportToArrow(proteins, peptides, export_all_psms)` | `QpxFile::export_to_arrow(..) -> Result<QpxOutput<RecordBatch>>` |
| `exportPSMsToQPXArrow(proteins, peptides, export_all_psms, feature_links)` | `QpxFile::export_psms_to_qpx_arrow(.., Option<&FeatureLinks>)` |
| `exportToParquet(proteins, peptides, filename, export_all_psms, config, feature_links)` | `QpxFile::export_to_parquet(..) -> Result<Vec<String>>` (the warnings) |
| `exportToParquet(table, filename, config, scan_format)` | `QpxFile::export_table_to_parquet(&RecordBatch, path, &ParquetWriteConfig, &str)` |
| `exportToParquetStreaming(proteins, ptrs, filename, export_all_psms, batch_size, config, n_threads, feature_links)` | `QpxFile::export_to_parquet_streaming(proteins, &[&PeptideIdentification], path, &StreamingOptions, Option<&FeatureLinks>)`; `StreamingOptions` holds the four defaulted knobs with the source's defaults |
| `requireResolvableMergeIndices(proteins, PeptideIdentificationList)` | `QpxFile::require_resolvable_merge_indices(&[ProteinIdentification], &[PeptideIdentification])` |
| `requireResolvableMergeIndices(proteins, vector<const PeptideIdentification*>)` | `QpxFile::require_resolvable_merge_indices_of(.., &[&PeptideIdentification])` |
| `importFromArrow(table, proteins, peptides)` | `QpxFile::import_from_arrow(..) -> Result<Vec<String>>`; `import_from_arrow_with_registry` resolves against a caller-owned `ModificationsDB` |

### `QPXCollectionExport`

| Source | Rust |
|---|---|
| `requireExportable(cmap, design)` | **not ported** (see Gaps) |
| `Transaction(directory)` | `Transaction::new(path)` |
| `Transaction::commit()` | `Transaction::commit(self)` |
| `~Transaction()` | `Drop` (logs to the global log) |
| — | `Transaction::roll_back(self) -> Vec<RollbackMessage>`, `directory()`, `preexisting()` (native) |
| the three file names | `QPX_COLLECTION_FILES` |

## Preserved source conventions

* **Identities** are BLAKE2b with an 8-byte digest (the digest length is in the
  parameter block: BLAKE2b-64, not truncated BLAKE2b-512) over compact,
  ASCII-escaped JSON, read big-endian and signed. Floats are widened to double
  and written as Python's `repr` writes them; set-valued lists
  (`pg_accessions`, `grouped_runs`) are deduplicated and ordered by their JSON
  encoding; `scan` keeps its order.
* **Validation** reports every diagnostic of the source, word for word and in
  the source's order; keys accumulate across batches only when a whole batch is
  valid.
* **Exports**: one row per exported hit (the first hit, not the best-scoring
  one, unless `export_all_psms`); identifications without hits add no row but
  count in `peptide_identification_index`. A PEP column from the main score or
  the first PEP-named metadata value of the first hit; additional scores are the
  known numeric score names (`_score` suffix stripped) with the orientation of
  the category the name parses as (null for engine names such as `expect`).
  `charge` narrows to `int16` wrapping, `calculated_mz` is 0 for a charge <= 0
  in the QPX view and null in the internal one; `mass_error_ppm` is computed in
  `float` then multiplied by the double `1e6`, as the source's expression is.
  Modifications are grouped by name in name order; a mass the sequence parser
  could not resolve (`M[+123.456]`) is written with an empty name and a null
  accession (measured). `run_file_name` prefers the hit's
  `reference_file_name` or `run_file_name`, then the identification's
  `reference_file_name`, then the merged-run origin file, stemmed.
* **Streaming** writes each batch as row groups of at most `row_group_size`
  (capped at 1,048,576) rows and a batch without rows as one empty row group,
  as Parquet C++'s `WriteTable` does (measured layouts `3:2,1,2` and `3:1,0,1`).
* **Import** groups rows by `peptide_identification_index` value, converts
  peptidoforms best-effort with the source's warnings, restores the dedicated
  columns as the source's metadata keys, and appends shell runs for unknown
  identifiers with the first non-empty `reference_file_name` as their path.
* The warnings, refusals and messages are the source's text; the refusals of
  `IdentifierMSRunMapper` that QPXFile surfaces were aligned with the source as
  part of this lane (`Error::MissingInformation` with the source's three
  messages; the duplicate-path message).

## Native differences

* **Tables** are `RecordBatch`es: the source's chunk handling, null tables and
  unchecked casts have no counterpart; a column of an unexpected type is a
  reported error, not undefined behaviour.
* **Metadata order.** The source lists `additional_scores`, `psm_metavalues`
  and `spectrum_metavalues` entries in `MetaInfoRegistry` order (the order the
  process first registered each name); this crate lists them in key order. The
  content is the same; the differential compares them as sorted lists.
* **Rank.** The source keeps a hit's rank as the `rank` metadata value (only
  when non-zero); `PeptideHit::rank` is a field here, written as that entry
  and restored from it.
* **Non-numeric metadata** under a numeric column (`predicted_RT`,
  `ion_mobility`, a PEP score, ...) is refused. The source's
  `DataValue::operator double` reads a string's storage as a double.
* **Errors instead of logs.** `false`/`nullptr` with a logged error is an
  `Err`; logged warnings are returned. The source's `ConversionError` for a
  non-integer `missed_cleavages` is `Error::InvalidValue` with the same text.
* **Streaming output** is written to a temporary sibling and renamed into
  place, so a refused batch leaves no file and never destroys an existing one
  (the source truncates first and removes its partial output). With more than
  one thread the source writes each partition as its own row groups, so its
  row-group boundaries depend on the thread count; the port writes the batch
  whole, so the file is the same for every count and equals the source's
  serial output. A non-positive `n_threads` means every core, as in the
  source.
* **Registries.** The import resolves peptidoform modifications against a
  caller-owned `ModificationsDB` (a copy of the global one by default). The
  source resolves against the process-wide database, which
  `AASequence::fromString` extends with every unknown mass it parses: a C++
  process that parsed `PEPTM[+123.456]IDEK` imports `PEPTM[+123.4560]IDEK`
  with the mass, a fresh process drops it with a warning. The port always
  behaves as the fresh process (measured both ways, see Evidence).
* **Transaction** joins paths instead of concatenating strings, so an empty
  directory is the current directory, not the root; its explicit
  `roll_back` returns the lines the destructor would log.
* **Arithmetic.** `calculated_mz` is the sequence's `mono_mass_for(Full,
  charge) / charge`, summed in the source's order. The sequence port sums a
  modified residue as its internal formula's mass plus the modification's,
  where the source subtracts water from the full residue's formula mass, so the
  double can differ by up to 2 ULP (measured on `PEM(Oxidation)TIDER`, charge
  2); the QPX view's `float` does not differ.

## Gaps (named dependencies)

* `QPXCollectionExport::requireExportable` needs
  `ConsensusMapArrowExport::requireUnambiguousIdentities` and
  `requireResolvableIdRuns`, the quantification-unit checks of
  `ProteinGroupArrowExport` (`FORMAT/ProteinGroupArrowExport_impl.h`) and
  `METADATA/ExperimentalDesign`. The class test's `requireExportable` sections
  are carried over when those land.
* The class test sections of `QPXFile_test.cpp` that exercise
  `ProteinGroupArrowExport` and `ConsensusMapArrowExport` are carried over for
  the psm view only.

## Evidence and tolerances

Tier 1, executed differential. The oracle driver `F1c-qpx`
(`../oracle/F1c-qpx/driver.cpp`) links the Linux x86_64 Release install
`openms4-release-bc9cc12-c19e494-174b576` (g++ 14.4 `-O2 -DNDEBUG`, Arrow C++
23.0.1) on ibminode06 and ran twice with identical records; a second process
imports the C++-written `internal_all.parquet` fresh. Its outputs are
`tests/data/qpx/` (hashes in `tests/data/qpx_provenance.json`):

* `tests/qpx_file.rs` compares every cell of both export layouts (best hit and
  all hits, 29 and 26 columns), the feature links, the warnings, the import
  (in-process and fresh), the re-export, the refusals, and five Parquet files
  (row-group layout, footer pairs, every cell; C++-written files are read back
  by the port). Tolerances: exact, except `calculated_mz` of the internal
  layout (relative 1e-15, see Native differences) and the ProForma parser's own
  message inside the "failed to parse peptidoform" warning, whose text is the
  parser port's (the source's reads `Expected ']' in: <input>`); QPXFile's
  framing is compared exactly.
* `tests/qpx_value_validation.rs` checks all 28 validation cases (psm, feature,
  pg, a second batch, a table of another view) against the executed
  diagnostics.

Tier 3: every section of `QPXIdentity_test.cpp` (its values come from the qpx
Python reference implementation, independent of both), of
`QPXValueValidation_test.cpp`, of the `QPXFile_test.cpp` sections that need
only `QPXFile`, and the `Transaction` section of `QPXCollectionExport_test.cpp`.

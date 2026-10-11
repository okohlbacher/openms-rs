# SqMassFile support

Port of `FORMAT/SqMassFile.h` / `SqMassFile.cpp` (core `bc9cc12`), the
adapter over the sqMass (SQLite) storage of `HANDLERS/MzMLSqliteHandler`
(`docs/MZML_SQLITE_HANDLER_SUPPORT.md`). Rust: `src/format/sqmass_file.rs`
(feature `sqmass`: `rusqlite` plus `mzml`). Tests: `tests/sqmass_file.rs` (7,
file-gated on `sqmass`). Provenance: `tests/data/sqmass_file_provenance.json`.

**Status: partial** — `convertToXICParquet` is not ported (see below).

## API mapping

| Source member | Rust |
|---|---|
| `SqMassConfig { write_full_meta, use_lossy_numpress, linear_fp_mass_acc }` | `SqMassConfig` with the same fields and defaults, plus the native `loss_policy` |
| `typedef MSExperiment MapType` | `MSExperiment` |
| `SqMassFile()`, `~SqMassFile()` | `SqMassFile::new()` / `Default`, drop |
| `load(filename, MapType&) const` | `load(filename) -> MSExperiment` |
| `store(filename, const MapType&) const` | `store(filename, &MSExperiment)` |
| `transform(filename_in, IMSDataConsumer*, skip_full_count, skip_first_pass) const` | `transform(filename_in, &mut impl MSDataConsumer, skip_full_count, skip_first_pass)` |
| `convertToXICParquet(filename_in, xic_filename, run_id, source_file, transition_exp) const` | **not ported**: needs an Arrow/Parquet writer, `MSChromatogramParquetConsumer` and `OpenSwath::LightTargetedExperiment` |
| `setConfig(const SqMassConfig&)` | `set_config(SqMassConfig)`; `config()` reads it back |
| protected `config_` | private field |

## Preserved source conventions

- One configuration applied to the handler for every call; the `RUN` ID is the
  source's `getSqlRunID` (from `sql_run_id`, else a non-negative integer
  `sqMassRunID` meta value, else 0).
- `transform` announces counts and settings, then streams spectra and
  chromatograms in input order in batches of 500 (the source's
  `batch_size`); `skip_full_count` and `skip_first_pass` are ignored, as in the
  source.
- `peptide_sequence` of the first precursor round-trips with the data.

## Native differences

- `LossPolicy::Reject` (the default) refuses to store metadata that neither
  the SQL tables nor the snapshot keep; `LossPolicy::Source` discards it
  silently as the source does.
- The source issues an empty last batch when a count is a positive multiple of
  500, which makes its handler select every record and throw
  `IllegalArgument`; empty batches are skipped here.
- A consumer may end the stream early with `ControlFlow::Break`.
- The lossless Numpress path releases its transient zlib state inside the
  codec budget (`src/format/numpress_coder.rs`), a fix made for this port.

## Checked boundaries and evidence

- Tier 3: every `START_SECTION` of `SqMassFile_test.cpp` with its tolerances,
  over the byte-identical class-test fixtures `SqliteMassFile_1.sqMass` and
  `MzMLSqliteHandler_1.mzML`.
- Tier 4: `transform` against `load`, and a record count that is a
  multiple of the batch size of 500. The `LossPolicy::Reject` refusal is
  exercised inside the store section. Early stop is not tested.

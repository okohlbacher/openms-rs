# OpenSWATH Parquet support

Ports the OpenSWATH Parquet I/O of `ANALYSIS/OPENSWATH` at OpenMS4-core
`bc9cc12`: the `.oswpq` assay library (`TransitionParquetFile.h`). All of it
sits behind the `arrow` feature (Arrow/Parquet, lane F1a) together with
`file-compression` (the ZIP archive layer, `ZipArchiveFile`, lane F7b), which
the default feature set already enables.

| Header | Rust | Status |
|---|---|---|
| `TransitionParquetFile.h` | `src/analysis/openswath/transition_parquet_file.rs` | complete |

## The `.oswpq` container

A directory, or a store-only ZIP archive of one:

```text
<name>.oswpq
└── library/
    ├── precursors.parquet     OSWPrecursorSchema
    ├── transitions.parquet    OSWTransitionSchema
    └── metadata.json          QC counts
```

A zipped container written by `TransitionParquetFile` also carries a JSON size
index as an entry named after the **staging** archive, `<name>.oswpq.tmp.idx.json`
(the source indexes `<path>.tmp` before renaming it into place; measured on the
Release build, and reproduced).

## TransitionParquetFile: API mapping

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

## Preserved source conventions

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

## Native differences

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

## Checked boundaries and evidence

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

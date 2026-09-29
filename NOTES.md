# O5c lane notes (OpenSWATH Parquet I/O)

## Merges done (all committed)
- worktree-wf_cd066dab-0a7-5 (F1a arrow core), worktree-wf_cd066dab-0a7-44 (O1b TraML/TSV/PQP) -- the brief's prereqs.
- EXTRA, needed and complete lanes (report to lead!):
  - worktree-wf_cd066dab-0a7-27 (OA: LightTargetedExperiment, data_access_helper convertTargetedExp)
  - worktree-wf_cd066dab-0a7-46 (O5b: OpenSwathLibraryIDNormalizer + SourceIDMapping; contains OA)
  - worktree-wf_a10fdd71-2d2-5 (F7b: ZipArchiveFile, zip =6.0.0 under file-compression)
  Cargo.lock resolved by union (arbitrary + arrow + half blocks).

## C++ sources
Extracted to scratchpad/cpp via `git -C <core> archive -o X.tar bc9cc12 '*Parquet*' ...` (the word
"source" in a path trips the worktree guard; use pathspec globs).

## Design facts learned
- TransitionParquetFile reader: C++ iterates unordered_map<int64,PrecursorInfo> to emit compounds ->
  order is stdlib-hash order. Rust: precursor row order (document divergence, C++ defect candidate).
- Missing entry in a directory layout throws FileNotFound (via ZipArchiveFile::extractEntryToTempFile),
  header docs say MissingInformation (doc defect candidate).
- Directory output: writes <dir>/library/{precursors,transitions}.parquet + metadata.json.
  Non-dir output: stage in temp dir, addOrReplaceFromFile each file into <abs>.tmp, writeSidecarIndex, rename.
- metadata.json uses VersionInfo getVersion/getTime.
- Class test TransitionParquetFile_test uses MRMAssay_detectingTransistionCompound_input.TraML.
- OpenSwathParquetExporter needs OpenSwathExportData.h (NOT ported anywhere) -> likely gap.

## Status
- TransitionParquetFile DONE (commit "[FEATURE,TEST,DOC] O5c: TransitionParquetFile"): module, tests (8 pass on spock,
  `--features arrow --test transition_parquet_file`), docs/OPENSWATH_PARQUET_SUPPORT.md, provenance json
  (generator: scratchpad/o5c/prov.py + tpf_spec.json).
- Oracle ../oracle/O5c-parquet done (ibminode06, identical runs).

## Next step
1. Writer+Reader code committed (WIP; clippy clean). Oracle ../oracle/O5c-oswpq launched (run.sh pull; fixture tworuns
   copied into fixtures/). Records F/Z/X/OK/E/L/PGC/PG/TFC/TF/UC/U (see driver header); files flattened case__rel with '/'->'#'.
2. Copy results into tests/data/openswath_osw_parquet/ + provenance; write tests/openswath_osw_parquet.rs:
   Rust writer on same inputs (mirror richFeatures()/secondRunFeatures()/classFeatures() exactly) -> compare every table
   value by value + metadata.json modulo build_time; reader outputs vs L/PG/TF/U dumps; X facts (run dir left behind).
3. Class tests: Writer_test, Reader_test, RoundTrip_test (needs traml_input records from transition_parquet_file data).
4. Exporter: OpenSwathExportData.h unported -> decide.

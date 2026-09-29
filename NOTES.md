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

## Next step (incarnation 3, after a restart)
Committed so far: TPF complete; Writer+Reader code (WIP); ExportData+ParquetExporter code (WIP);
tests/openswath_osw_parquet.rs draft (1341 lines, NOT yet run on spock) + tests/data/openswath_osw_parquet/
(oracle O5c-oswpq results, ibminode06 Release, runs identical).
Oracle O5c-export ran too (ibminode06, identical, results pulled to ../oracle/O5c-export/results): 7 files
fs_{class,empty,rich,empty_cols}, ts_{class,empty,rich} + err_unwritable FileNotWritable.
1. Run the osw parquet test on spock (--features arrow --test openswath_osw_parquet); fix.
2. Copy O5c-export results into tests/data/openswath_parquet_exporter/, write tests/openswath_parquet_exporter.rs.
3. Provenance json, support doc section, ledger, report, delete NOTES.

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

## Next step
1. Write src/analysis/openswath/transition_parquet_file.rs (cfg arrow + file-compression), tests/transition_parquet_file.rs
   transcribing the class test.
2. Then OSWParquetWriter/Reader (read headers + tests), exporter last.

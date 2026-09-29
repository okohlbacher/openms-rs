# O5c lane notes (OpenSWATH Parquet I/O)

Prereqs merged: worktree-wf_cd066dab-0a7-5 (F1a arrow core), worktree-wf_cd066dab-0a7-44 (O1b TraML/TSV/PQP).

C++ sources extracted (git -C <core> archive bc9cc12 '*Parquet*' '*OSWFile*' ...) into
scratchpad/cpp/ (source-path words trip the worktree guard; use pathspec globs).

Headers (bc9cc12): ANALYSIS/OPENSWATH/{TransitionParquetFile 189+759, OpenSwathOSWParquetWriter 111+1309,
OpenSwathOSWParquetReader 260+850, OpenSwathParquetExporter 34+373}.
Class tests: TransitionParquetFile_test(444), OpenSwathOSWParquetWriter_test(163), ...Reader_test(311),
...RoundTrip_test(86), OpenSwathParquetExporter_test(116).
TOPP data: OpenSwathAssayGenerator_parquet_output.tsv, OpenSwathDecoyGenerator_parquet_output.tsv,
TargetedFileConverter_parquet_roundtrip.tsv.

## Next step
1. Read the four headers + tests; order: TransitionParquetFile first (TOPP consumers), then writer/reader, exporter.

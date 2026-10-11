# O5a notes (OpenSWATH workflow + writers)

Status at start: prerequisites merged (0a7-44, 0a7-47, 48e-6, 0a7-49; no earlier bcdd783d incarnation existed for this lane).

Headers in scope (all `unmapped` in docs/core-sdk-coverage.json at start), C++ at bc9cc12 under
`src/openms/include/OpenMS/ANALYSIS/OPENSWATH/` and `src/openms/source/ANALYSIS/OPENSWATH/`:
OpenSwathExportConfig.h (83 lines), OpenSwathExportData.h (352), OpenSwathResultsExporter (.h 36/.cpp 354),
OpenSwathMatrixExporter (.h 40/.cpp 532), OpenSwathOSWWriter (.h 339/.cpp 1515), OpenSwathWorkflowScheduler (.h 213/.cpp 258),
OpenSwathWorkflow (.h 461/.cpp 1627).
Out of scope but adjacent: OpenSwathOSWParquetReader/Writer, OpenSwathParquetExporter (Parquet needs arrow/parquet crates, not in Cargo.toml).

Upstream tests (class_tests/openms/source/): OpenSwathExport_test, OpenSwathMatrixExporter_test, OpenSwathOSWWriter_test,
OpenSwathResultsExporter_test, OpenSwathWorkflowScheduler_test. No OpenSwathWorkflow_test (TOPP-tested only).

Tooling hint: guard refuses `git -C ../x`; use /private/tmp/claude-501/.../scratchpad/dump.py <regex> to extract pinned C++ into scratchpad/cpp/.

Rust ground: src/analysis/openswath/* (transition_tsv_file, targeted_experiment, openswath_scoring, mrm_feature_finder_scoring,
chromatogram_extractor, swath_map_mass_correction, spectrum_access_*), src/format/sqlite_connector.rs, `sqlite` feature = rusqlite.

## Next step
1. Port ExportConfig + ExportData (enums, rows) -> src/analysis/openswath/openswath_export.rs (or split).
2. ResultsExporter (TSV; Parquet = typed error naming the missing ParquetFile dependency).
3. MatrixExporter. 4. OSWWriter (rusqlite). 5. WorkflowScheduler. 6. Workflow.

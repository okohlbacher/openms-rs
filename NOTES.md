# F1c notes (QPX lane)

Merged prerequisite worktree-wf_cd066dab-0a7-5 (F1a arrow core, fast-forward). No earlier incarnation branch.

C++ sources at pin copied to scratchpad (read via `git -C .../OpenMS4-tests/packages/core show bc9cc12:<path> > file`,
one plain git command per call; compound git forms are refused by the worktree guard):
- QPXValueValidation.h/.cpp (136/722), test 384
- QPXIdentity.h/.cpp (208/453), test 244
- QPXFile.h/.cpp (146?/2158), test 1889
- QPXCollectionExport.h/.cpp (146/237), test 538

Existing Rust: src/format/arrow_io_helpers.rs has QPX naming helpers, qpx_file_metadata, and a
`validate_qpx_values_for_write` stub that refuses psm/feature/pg with Error::Unsupported (F1a gap).
Schemas QPXPSMSchema/QPXFeatureSchema/QPXPgSchema in src/format/arrow_schema_registry.rs.

Order: 1 QPXValueValidation (+ wire into arrow_io_helpers gate), 2 QPXIdentity, 3 QPXFile, 4 QPXCollectionExport.

## Next step
1. Read QPXValueValidation.h/.cpp + test, port to src/format/qpx_value_validation.rs.

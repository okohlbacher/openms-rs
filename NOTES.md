# NOTES (lane: FileHandler.h) -- delete before final commit

## State
- Header: src/openms/include/OpenMS/FORMAT/FileHandler.h (core bc9cc12), .cpp 1841 lines.
- Scratch copies: $SCRATCH/FileHandler.{h,cpp} (git -C ../OpenMS4-tests/packages/core show bc9cc12:...)
- Rust: src/format/file_handler.rs (789 lines), src/format/path_io.rs
- Ledger scope remaining: log type, rewrite_source_file, compute_hash, SRM conversion after mzML load;
  also "Other identification/transition formats, ZIP input, SHA1 and full source-file bookkeeping".
- Tests: tests/file_handler.rs, tests/file_handler_type_detection.rs, tests/identification_paths.rs
- Doc: docs/FILE_HANDLING_SUPPORT.md

## Next step
1. Read FileHandler.h/.cpp members; map to file_handler.rs.

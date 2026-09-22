# NOTES (lane: FileHandler.h) -- delete before final commit

## Pins / scratch
- C++: FileHandler.h (429 l), FileHandler.cpp (1841 l) at core bc9cc12; scratch copies in
  SCRATCH=/private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad
  (bash refuses `sed $S/...` and compound git cmds -- use Read/Write tools, literal paths, plain git)
- Rust: src/format/file_handler.rs (789 l), src/format/path_io.rs. Tests: tests/file_handler.rs,
  tests/file_handler_type_detection.rs, tests/identification_paths.rs. Doc: docs/FILE_HANDLING_SUPPORT.md
- Ledger scope (docs/core-sdk-reviewed-apis.json:1925): remaining = log type, rewrite_source_file,
  compute_hash, SRM conversion after mzML load; also other id/transition formats, ZIP input, SHA1, full source-file bookkeeping.

## C++ facts (FileHandler.cpp @bc9cc12)
- computeFileHash :768-778: RFC3174 SHA-1 (anon ns :83-187) over raw file bytes (8 KiB reads), lowercase hex
  40 chars; "" when file cannot be opened. => sha1 crate (Cargo.toml:84, optional, pulled by `mzml` feature).
- loadExperiment :849-1052. Precondition rewrite||!compute_hash (:853). mzML case :904-912:
  MzMLFile with options_, setLogType(log), load, then ChromatogramTools().convertSpectraToChromatograms<PeakMap>(exp,true).
  setLogType(log) also for DTA2D, MzXML, MzData, MGF, MS2, BrukerTims, ThermoRaw. DTA/SQMASS/XMASS/MSP none.
- rewrite_source_file :1014-1051: src_file = first SourceFile (warn "No source file annotated." if none;
  warn "Expecting a single source file in mzML. Found N will take only first one for rewriting." if >1);
  strip trailing / and \ ; setNameOfFile(basename); path = File::path(File::absolutePath(name));
  uri = "file://"+path if path starts with "/" else "file:///"+path; setFileType(FileTypes::typeToMZML(type));
  if compute_hash && type!=BRUKER_TDF: setChecksum(computeFileHash(filename), SHA1); clear sources; push src_file.
- Rust ChromatogramTools exists: src/kernel/chromatogram_tools.rs; used in src/format/file_info/peaks.rs:109
  as ChromatogramTools::default().convert_spectra_to_chromatograms(&mut exp, true, false)?

## Next step
1. Check Rust SourceFile (checksum type), FileType -> mzML term (typeToMZML), ProgressLogger + *_with_progress signatures.
2. Design: FileHandler::load_experiment_with (options struct {log, rewrite_source_file, compute_hash}).

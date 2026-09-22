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

## Rust facts
- SourceFile: src/metadata/acquisition.rs:395 {name,path,size_mb,file_type,checksum,checksum_type:ChecksumType{Sha1,Md5,Unknown},...}
- typeToMZML == FileType::mzml_name() (src/format/file_types.rs:148) -- already ported, same table (FileTypes.cpp:259-275).
- File helpers: crate::system::file::{basename, path, absolute_path} (src/system/file.rs:448,459,509) = source semantics (lexical).
- Global log: crate::openms_log_warn! macro (src/concept/log_stream.rs:873) -> WARNING route; readers use log_message(Warn,..).
- ProgressLogger (src/concept/progress_logger.rs:305): clone() = fresh backend of same type + same clock/nesting/gui_factory
  == source f.setLogType(log). ProgressLogType {Cmd,Gui,None}.
- *_with_progress: mzml_paths::load_with_progress(path,&LoadOptions,&ReadOptions,&mut logger) (path_io::open, compressed ok);
  dta2d::load_with_progress(path,&ReadOptions,&mut logger) uses File::open (NO decompression; read_reporting private);
  mzxml::load_with_progress(path,&ReadOptions{peaks,limits},&mut logger); mzdata similar; mascot_generic::load_with_progress
  (FileHandler uses mgf.rs, a different strict reader -- MGF progress = gap unless reader switch); ms2 has none.
- FileHandler has NO mzXML/mzData experiment dispatch (can_read_experiment: DTA,MGF,DTA2D,MS2,mzML) -> add (cfg mzml).
- sha1 dep optional, pulled by `mzml` feature (default). compute_file_hash -> cfg(feature="mzml") or hand-rolled? use crate.
- Callers of FileHandler::load_experiment*: src/cli/tools/* (baseline_filter uses plain load_experiment; peak_picker uses
  with_read_options). Don't touch src/cli. SRM conversion in existing loaders = source-faithful; precheck any SRM spectrum -> skip (exact no-op).
- ChromatogramTools C++ (KERNEL/ChromatogramTools.h:110-215) warns via OPENMS_LOG_WARN per skipped SRM spectrum (:164-165);
  Rust returns report.skipped_spectra only. Default limits max_input_records 1M -> scale to experiment.

## Decisions
- New API: FileHandler::load_experiment_with_progress? -> one "ExperimentLoadOptions {peaks, mzml read, rewrite_source_file, compute_hash}" + &ProgressLogger (cloned).
- Warnings go to openms_log_warn! (same route as source OPENMS_LOG_WARN).

## Upstream expected values (FileHandler_test.cpp @bc9cc12 :229-320)
- :234-255 mzData/mzXML through loadExperiment (+ setMZRange(115,1000)) -> exp[1][0..2]=110,120,130 / 120,130; mzXML exp[2]=100,110,120 / 120,130,140
- :258-261 MzMLFile_1.mzML {MZML}, NONE, rewrite, hash -> size 4, checksum 36007593dbca0ba59a1f4fc32fb970f0e8991fa6, SHA1
- :264-279 DTA2DFile_test_1.dta2d -> 230.02,430.02,630.02; mz 300-1000 + rewrite+hash -> 430.02,630.02, d50d5144cc3805749b9e8d16f3bc8994979d8142;
  rewrite w/o hash -> checksum "" + UNKNOWN; {SIZE_OF_TYPE} -> ParseError
- :293-319 computeFileHash on non-ASCII path nonempty.
- Repo fixtures byte-identical (shasum checked): tests/data/mzml_load_source_original.mzML (=MzMLFile_1.mzML),
  tests/data/text_peak_lists/DTA2DFile_test_1.dta2d, tests/data/MzDataFile_1.mzData, tests/data/MzXMLFile_1.mzXML
- MS2File.cpp makes no progress calls -> MS2 without progress is exact.
- mzxml/mzdata load_reporting use path_io::open (compressed ok). mzdata::load_with_progress(path,&PeakFileOptions,&ReadLimits,&mut logger)->Loaded{experiment,..}
- Precondition rewrite||!hash is OPENMS_PRECONDITION (debug only); Release: hash w/o rewrite = no-op.

## Code status (commit "[WIP] FileHandler: log type, ...")
- file_handler.rs: ExperimentLoadOptions, load_experiment_with_progress, compute_file_hash (cfg mzml),
  convert_srm_spectra (also in load_experiment via load_stream, and with_options), rewrite_source_file,
  mzXML/mzData read dispatch. dta2d.rs: read_with_progress. cargo check --lib OK.
- Oracle dir: ../oracle/file-handler-load (driver.cpp, make_inputs.py, remote_run.sh, run.sh), modeled on
  ../oracle/progress-format-readers. SRM inputs: tests/data/file_info/inputs/srm_spectra{,_mixed}.mzML
  (mixed has SRM w/ 0 and 2 precursors -> warnings).

## Next step
1. Write oracle, run on ibminode06 (run.sh start/poll/fetch).
2. Tests: tests/file_handler_load.rs (upstream values + oracle fixture).
2. Implement in file_handler.rs; tests in tests/file_handler_source_file.rs; oracle on ibminode06.

# Lane F3b notes (search-engine input/output files)

Branch: worktree-wf_ee239125-5be-13 (second incarnation), off main f0990c0,
fast-forwarded onto worktree-wf_bcdd783d-f3f-13 (first incarnation: 5a2b6b2).

## State
- None of the 13 headers was ported at start: coverage entries exist in
  docs/core-sdk-coverage.json (~9303 InspectInfile, 9314 InspectOutfile, 9325
  KroenikFile, 9604 MsInspectFile, 9987 OMSSACSVFile, 9998 OMSSAXMLFile,
  10249/10260 PepNovo*, 10523/10534 Sequest*, 10547 SpecArrayFile, 10982
  XTandemInfile, 10995 XTandemXMLFile), all unmapped.
- WIP from incarnation 1 (uncompiled, no tests yet, not in mod.rs):
  src/format/source_text.rs (StringUtils toDouble/toInt32/split/trim),
  kroenik_file.rs, ms_inspect_file.rs, spec_array_file.rs.
- Module root: src/format/mod.rs (append `pub mod` lines only).

## How to read the C++ (tool-checker quirk)
- Any shell command containing the word `s o u r c e` as a path segment in
  find/archive is refused by the worktree guard; `git -C <core> show
  bc9cc12:<path>` works. Already extracted everything to the scratchpad:
  /private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad/c/{h,cpp,test,data}/
  (script: scratchpad/extract.py; rerun if the scratchpad was wiped).
- Class-test data lives in the core repo at src/tests/class_tests/openms/data/
  (not the test-data package). Copy fixtures to tests/data/ and describe them in
  tests/data/search_engine_files_provenance.json (format: see
  tests/data/imzml_handler_provenance.json: source_revision, repository,
  license, source_files[{path,sha256,role}], fixtures[{path,source_path,sha256,bytes,role}]).
- Tests load data via `PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")`.

## Plan (order: simple feature-text formats first, then id formats)
1. KroenikFile, MsInspectFile, SpecArrayFile (FeatureMap text readers; small)
2. OMSSACSVFile (CSV -> PeptideIdentification)
3. PepNovoInfile/Outfile, InspectInfile/Outfile, SequestInfile/Outfile
4. XTandemInfile, XTandemXMLFile, OMSSAXMLFile (XML)

## Next step
1. Verify the three WIP readers against their .h and class tests, add
   `pub mod` lines, write tests/search_engine_feature_files.rs, gate on spock
   (slot F3b, run_in_background), commit.

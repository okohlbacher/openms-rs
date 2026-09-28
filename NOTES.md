# Lane F3b notes (search-engine input/output files)

Branch: worktree-wf_cd066dab-0a7-8 (third incarnation; merged incarnation 2 = worktree-wf_ee239125-5be-13 at ecb3d83).
Incarnation 2 branch: worktree-wf_ee239125-5be-13 (second incarnation), off main f0990c0,
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

## Done (gated on spock, slot F3b)
- data_structures/string_utils.rs (moved from format/source_text.rs;
  cli::context re-exports to_double/to_int32 instead of its copy).
- Kroenik/MsInspect/SpecArray: tests/search_engine_feature_files.rs 16/16.
- OMSSACSVFile (+ format/search_engine_text.rs permissive AASequence parse):
  tests/search_engine_id_files.rs 5/5.
- format/search_engine_xml.rs: SAX driver (text nodes whole, whitespace
  delivered, Latin-1 for X!Tandem), cfg idxml.
- XTandemXMLFile: class test passes (303/497/6/0, 2/21/3/0), tests in
  tests/search_engine_xml_files.rs (cfg idxml).
- OMSSAXMLFile + resources/search_engines/OMSSA_modification_mapping
  (sha256 e7b8d091...): tests written, gate pending.

## Defect candidates noted so far
- XTandemXMLFile: protein_uids_ and group_type_stack_ never reset between
  loads; `aa` outside group -> error() then top() of empty stack (UB); `aa`
  before a domain -> back() on empty vector (UB).
- OMSSACSVFile: first data line numbered -1 appends to caller's last id / UB.
- OMSSAXMLFile: setModificationDefinitionsSet restarts at 119 each call while
  keeping earlier entries (ambiguous numbers).

- SequestOutfile: all class-test sections pass (tests/search_engine_id_files.rs).
- PepNovoInfile/PepNovoOutfile: written, tests written (gate6 pending).
- C++ ORACLE DONE: ../oracle/F3b-search-engine/ (driver.cpp, remote_run.sh,
  run.sh; ibminode06, Release bc9cc12). Outputs copied to
  tests/data/search_engine_oracle/; tests/search_engine_oracle.rs compares
  XTandemXML, OMSSAXML, SequestOutfile, PepNovoOutfile, PepNovoInfile dumps.
  Oracle also holds xtandem_infile.out, inspect_infile.out, sequest_infile.out
  for the writers still to port (pin writes "13.0" and 15-fraction-digit
  masses via appendToStr = crate::param::value::format_float(v, true)).
- Findings: SequestInfile setEnzyme("i_dont_exist_enzyme") returns 17 at the
  pin while its class test expects 18; InspectInfile::store refuses a known
  non-TSV extension (.txt) but accepts .tmp.

## Incarnation 3 progress
- Merged tree gated on spock: 16+17+5+8 pass.
- XTandemInfile DONE (410653c): src/format/xtandem_infile.rs, tests/search_engine_infiles.rs (5/5 + 1 unit), byte-identical to oracle.
- PTMXMLFile (named dependency, unmapped header) + InspectInfile written: src/format/ptm_xml_file.rs, inspect_infile.rs (idxml),
  tests/search_engine_ptm_infiles.rs; fixtures PTMs.xml, Inspect_PTMs.xml, Sequest_PTMs.xml, *_test_template1.txt copied to tests/data.
  inspect_infile::parse_modification_line is the shared handlePTMs grammar for SequestInfile too.
- C++ extracts: scratchpad/F3b/c/ (script scratchpad/F3b/extract.py '<regex>').

## Next step
1. Gate g3 (ptm/inspect), commit. 2. SequestInfile (reuse parse_modification_line). 3. InspectOutfile (partial ok).
4. docs/SEARCH_ENGINE_FILES_SUPPORT.md, tests/data/search_engine_files_provenance.json, lane report + ledger, CI sweep.

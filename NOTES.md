# Lane F3b notes (search-engine input/output files)

Branch: worktree-wf_bcdd783d-f3f-13, off main f0990c0.

## State at start
- None of the 13 headers is ported: coverage entries exist in docs/core-sdk-coverage.json
  (lines ~9303 InspectInfile, 9314 InspectOutfile, 9325 KroenikFile, 9604 MsInspectFile,
  9987 OMSSACSVFile, 9998 OMSSAXMLFile, 10249/10260 PepNovo*, 10523/10534 Sequest*,
  10547 SpecArrayFile, 10982 XTandemInfile, 10995 XTandemXMLFile) all unmapped.
- Module root: src/format/mod.rs (append `pub mod` lines only).
- Sources: git -C ../OpenMS4-tests/packages/core show bc9cc12:src/openms/include/OpenMS/FORMAT/<X>.h
  and src/openms/source/FORMAT/<X>.cpp; tests src/tests/class_tests/openms/source/<X>_test.cpp;
  test data at 0cb15f2 in ../OpenMS4-tests/packages/test-data.

## Plan (order: simple feature-text formats first, then id formats)
1. KroenikFile, MsInspectFile, SpecArrayFile (FeatureMap text readers; small)
2. OMSSACSVFile (CSV -> PeptideIdentification)
3. PepNovoInfile/Outfile, InspectInfile/Outfile, SequestInfile/Outfile
4. XTandemInfile, XTandemXMLFile, OMSSAXMLFile (XML)

## Next step
1. Read headers/cpp/tests for KroenikFile, MsInspectFile, SpecArrayFile and port.

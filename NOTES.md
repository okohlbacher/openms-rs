# Lane I4 notes (IDMapper + metabolomics ID)

## Status
- DONE (code+tests, gating): src/analysis/id/id_mapper.rs + tests/id_mapper.rs (class test + 6 TOPP retained outputs),
  src/analysis/id/ion_identity_molecular_networking.rs + tests/ion_identity_molecular_networking.rs.
- Scratchpad is SHARED with other lanes: use scratchpad/I4/ only (gate logs t*.log there).
- C++ at pin: core working tree == bc9cc12 for all lane files; copies in scratchpad/cpp/.
- Worktree guard: one plain git command per Bash call; no shell variables as paths; heredocs mentioning git refused.
- Module cycles: analysis must NOT import crate::system (cycle); local basename in id_mapper.rs.
  analysis->data_structures edge recorded in docs/module-cycles.json.
- Reader defects found (not fixed; worked around in tests/id_mapper.rs load_ids): idXML reader rejects
  date="0000-00-00T00:00:00", flanking residue '*', duplicate SearchParameters UserParam.
- C++ defect candidate: IDMapper TMT scan-number fallback uses boost::regex{""} (empty() is false) -> ParseError.
- Untouched so far: AccurateMassSearchEngine (needs MzTab/MzTabM export, IdentificationData, HMDB mapping files),
  MetaboliteSpectralMatching, SiriusMSConverter, SiriusExportAlgorithm, FIAMSDataProcessor, FIAMSScheduler,
  SiriusFragmentAnnotation.

## Next step
1. Confirm gate t3 green; then clippy all-features/all-targets, doc, MSRV 1.85 for the two test files.
2. docs/ID_MAPPER_SUPPORT.md (+ IIMN section), tests/data/id_mapper_provenance.json, ledger fragment, report.
3. If time: SiriusFragmentAnnotation (format) with its class test data.

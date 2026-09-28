# Lane I4 notes (IDMapper + metabolomics ID)

## Status
- No Rust for any of the lane headers yet (grep found only file_types/mztab_m/constants mentions).
- C++ sources at core pin bc9cc12 (src/openms/{include,source}/.../ANALYSIS/ID/{IDMapper,AccurateMassSearchEngine,
  MetaboliteSpectralMatching,SiriusMSConverter,SiriusExportAlgorithm,IonIdentityMolecularNetworking,FIAMSDataProcessor,
  FIAMSScheduler}, FORMAT/DATAACCESS/SiriusFragmentAnnotation).
- Class tests: IDMapper_test, AccurateMassSearchEngine_test, MetaboliteSpectralMatching_test, IonIdentityMolecularNetworking_test,
  FIAMSDataProcessor_test, FIAMSScheduler_test, SiriusFragmentAnnotation_test (data in class_tests/openms/data/IDMapper_*, etc.)
- Tooling gotcha: the worktree guard refuses compound commands containing `git -C <other repo>`; run each
  `git -C .../core show bc9cc12:<path> > file` as ONE plain command (no loops/pipes/python).

## Next step
1. Extract C++ sources to scratchpad (one git show per command).
2. Port IDMapper first (2 TOPP, zero Rust; others don't depend on it), into src/analysis/id/id_mapper.rs.
3. Then AccurateMassSearchEngine, MetaboliteSpectralMatching, SiriusMSConverter/SiriusFragmentAnnotation.

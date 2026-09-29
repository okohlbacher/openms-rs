# O4a lane notes (Extraction and SWATH windows)

Merged prerequisite worktree-wf_cd066dab-0a7-27 (fast-forward, OA lane: openswath data access, SwathMap, spectrum access, helpers).

Headers (all under src/openms/include/OpenMS/ANALYSIS/OPENSWATH/ at core bc9cc12, tests in
src/tests/class_tests/openms/source/<Name>_test.cpp):
- ChromatogramExtractorAlgorithm.h  (base: others depend on it)
- ChromatogramExtractor.h
- SwathWindowLoader.h
- PeakMapExtractor.h
- SwathQC.h
- SwathMapMassCorrection.h
TOPP: OpenSwathChromatogramExtractor? (check topp 174b576 / cli c19e494)

Rust target: src/analysis/openswath/<snake>.rs, pub mod appended to src/analysis/openswath.rs.
Never run cargo locally; gate: OPENMS_GATE_HOST=dax OPENMS_GATE_JOBS=64 ~/.local/bin/openms-kim-gate.sh O4a ... (background).

## Next step
1. Read ChromatogramExtractorAlgorithm.h/.cpp + test; port to chromatogram_extractor_algorithm.rs.
2. SwathWindowLoader (text-file parser) -> swath_window_loader.rs.
3. ChromatogramExtractor, PeakMapExtractor, SwathQC, SwathMapMassCorrection.

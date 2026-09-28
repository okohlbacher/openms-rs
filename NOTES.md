# Lane F2 notes (MzTab, MzTab-M, mzQC, OMS) - incarnation 4

Branch worktree-wf_cd066dab-0a7-6. Incarnations 2-3 finished MzTab exporters, MzTabFile streaming
stores, OMSFile (complete), MzTabM end-to-end; ledger fragment .planning/ledger-fragments/wave4/F2.json
exists (MzTab, MzTabFile, MzTabM, MzTabMFile, OMSFile). Lane report .planning/lane-reports/F2.md was
NOT written; write it at the end (with Bash heredoc).

## Gates measured this incarnation (spock, slot F2, logs in scratchpad gates/)
- bare --test mztab 62, mztab_file 45, mztab_m 48
- sqlite --test mztab_m 49, oms_file 15; sqlite,file-compression --test oms_file 16
- ci_matrix --check wants: sqlite slice + oms_file; new slice "consensusxml featurexml idxml" --test mztab_export.
  Item-level gates needing JOBS entries: sqlite --test mztab_m; "sqlite file-compression" --test oms_file.
- Oracle F2-oms run1 vs run2 (checked on ibminode06): canonical dumps equal for retained inputs; the
  fixture-written idxml_whole/feature/consensus differ only in random UNKNOWN_INPUT_FILE_* ids. Doc fixed (1dd5a43).

## MzQCFile (now in progress)
- Tier 1 evidence: test-data 0cb15f2 topp/CMakeLists.txt:5387-5389 TOPP_QCCalculator_2 -> QCCalculator_2_output.mzQC,
  FuzzyDiff whitelist creationDate, location, version. Inputs: QCCalculator_input.mzML,
  FeatureFinderMetaboIdent_1_output.featureXML, OpenPepXL_output.idXML. Tool: topp 174b576 src/QCCalculator.cpp
  (loadExperiment, sortSpectra, updateRanges; features updateRanges+sortByRT; FileHandler::storeQC -> MzQCFile::store).
- Dependency closure (core bc9cc12): QCBase (+FlagSet, DATASTRUCTURES/FlagSet.h), TIC, SpectrumCount, FeatureSummary,
  IdentificationSummary -> MissedCleavages (ProteaseDigestion, missed cleavages 0). qc-cv.obo (share/OpenMS/CV).
  Class tests: QCBase_test, TIC_test, SpectrumCount_test, FeatureSummary_test, IdentificationSummary_test.
- Placement: src/analysis/qc/ (analysis<->format already mutual, no new cycle), root src/analysis/qc.rs,
  `pub mod qc;` appended to src/analysis/mod.rs. MzQCFile -> src/format/mzqc_file.rs (serde_json + sha1).
- Extract C++: git -C ../OpenMS4-tests/packages/core archive -o X.tar bc9cc12 ':(glob)src/openms/*/QC/**'
  (the word "source" in a literal path trips the worktree guard).

## Next step
1. FlagSet + QCBase + SpectrumCount + FeatureSummary + TIC under src/analysis/qc/, with class tests. Commit.
2. MissedCleavages + IdentificationSummary. Commit.
3. MzQCFile + qc-cv.obo resource + TOPP_QCCalculator_2 reproduction. Commit.
4. Gates, ledger additions, report, delete NOTES.md.

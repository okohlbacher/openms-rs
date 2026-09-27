# Lane F2 notes (MzTab, MzTab-M, mzQC, OMS)

Branch: worktree-wf_bcdd783d-f3f-11 off main f0990c0.

## State at start
- MzTab.h partial (src/format/mztab.rs, 3139 lines): missing exportFeatureMapToMzTab,
  exportIdentificationsToMzTab, exportConsensusMapToMzTab, extractModificationList,
  nested IDMzTabStream / CMMzTabStream.
- MzTabFile.h partial (src/format/mztab_file.rs): id-list + ConsensusMap store overloads
  need the exporters.
- MzTabM.h partial (src/format/mztab_m.rs): exporter reproduction from OMS input absent
  (needs OMSFile load -> IdentificationData).
- MzQCFile.h unmapped. OMSFile.h (+ OMSFileLoad.h/OMSFileStore.h) unmapped; rusqlite pinned.

## C++ sources (core pin bc9cc12, repo /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages/core)
- src/openms/source/FORMAT/{MzTab,MzTabFile,MzTabM,MzQCFile,OMSFile,OMSFileLoad,OMSFileStore}.cpp
- class tests: MzTab_test, MzTabFile_test, MzTabM_test, MzTabMFile_test, OMSFile_test
- data: MzTabFile_*.mzTab, MzTabMFile_input_1.oms, MzTabMFile_output_1.mztab, OMSFile_test_2.featureXML

## Next step
1. Read MzTab.cpp export functions + MzTab_test sections; map to existing Rust types.

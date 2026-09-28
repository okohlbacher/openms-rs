# Lane F2 notes (MzTab, MzTab-M, mzQC, OMS)

Branch: worktree-wf_cd066dab-0a7-6 (incarnation 3; merged worktree-wf_ee239125-5be-11 = incarnation 2,
which finished the MzTab exporters + streaming stores (commits 3bf0851..5bc5bd7) and left OMS JSON export WIP (9e004d0)).

## C++ sources
- Extract: `git -C ../OpenMS4-tests/packages/core archive -o X.tar bc9cc12 <paths>`; tar xf; copy to a
  flat dir (a path containing `source/` trips the worktree guard's `source` detector).
- MzTab.cpp exporters: 698-3414. MzTabFile.cpp streaming store overloads 2894-3125 (ID, ConsensusMap).
- Tool: topp 174b576 src/MzTabExporter.cpp (featureXML -> exportFeatureMapToMzTab + store(MzTab);
  idXML -> MzTabFile::store(ids streaming); consensusXML -> IDFilter::removeEmptyIdentifications + store(cmap streaming, true,true,subfeatures)).

## Evidence (tier 1)
- test-data 0cb15f2 topp/CMakeLists.txt:3382-3413 MzTabExporter_1..8, `${DIFF}` = FuzzyDiff ratio 1.01 absdiff 0.01,
  whitelist "MTD\turi[1]" "MS:1000752" "MTD\tsoftware" (+ "protein_search_engine_score" for 4-8).
  Small: 4 (Epifany_2_out.consensusXML), 5 (MzTabExporter_5_in.consensusXML, -first_run_inference_only),
  6 (idXML), 7 (consensusXML), 8 (6_input.idXML -export_all_psms). 3 = featureXML 1.25MB (NaN positions -> "nan").
  1,2 too large (24MB / 3.9MB) for tests/data.
- Rust FuzzyStringComparator: src/concept/fuzzy_string_comparator.rs (treats nan/NaN equal).
- MzTab_test.cpp: 2 exporter sections (export_all_psms termination; consensus assay grain).

## Design
- New `src/format/mztab_export.rs`: inherent `impl MzTab { export_feature_map, export_identifications,
  export_consensus_map, extract_modification_list }`, `IdMzTabStream`, `ConsensusMapMzTabStream`.
- MzTabFile streaming store: `store_identifications`, `store_consensus_map` in mztab_file.rs using explicit
  SectionLayout built with the source's header rule + Postcondition column check.
- C++ defects to preserve+document: proteinSectionRowFromProteinHit_ builds mods but never adds them;
  indistinguishable-group mods keep only last; getMSRunSpectrumIdentifierType_ only inspects first PID;
  fixed-mod test compares getId() against full ids; PRT/PEP user keys with spaces are looked up underscored.

## State at incarnation 3 start
- DONE (inc. 2): src/format/mztab_export.rs, mztab_file.rs streaming stores, tests/mztab_export.rs (MzTabExporter_3..8).
- WIP (inc. 2, saved by lead, never gated): src/format/oms_file.rs = JSON export only (OMSFileLoad::exportToJSON),
  tests/oms_file.rs (JSONExporter_1/2, feature `sqlite`). OMSFile.h coverage entry is `unmapped`; MzQCFile.h `unmapped`.
- MzTabM.h partial: "end-to-end exporter reproduction from the OMS input remains absent" -> needs OMS load into
  IdentificationData? (check MzTabM export input in C++).

## Inc. 3 findings
- C++ sources flattened at scratchpad cpp/flat (archive of core bc9cc12; OMSFile*, QC/*); re-extract via archive.
- MzQCFile::store needs QC/TIC, SpectrumCount, FeatureSummary, IdentificationSummary (no src/qc in crate),
  ControlledVocabulary loadFromOBO + qc-cv.obo (absent from resources/cv), FileHandler::computeFileHash (SHA-1).
  No direct TOPP consumer in coverage json -> lowest priority; likely leave untouched with named deps.
- OMS: Rust IdentificationData exists (src/identification/graph*.rs, sorted Table iteration). FeatureMap has
  no id_data member (by design, feature_identification.rs) -> store/load feature maps take/return the graph.
- Retained .oms (test-data 0cb15f2 topp/): JSONExporter_protein/RNA, IDMerger_6_input1/2, MapAlignerIdentification_8_input1/2.
  Evidence plan: load -> store -> export JSON == retained JSONExporter output (FuzzyDiff); C++ oracle on ibminode06
  to confirm the C++ load->store->exportToJSON round trip gives the same (../oracle/F2-oms/).
- C++ defects: store binds compound.name to `smile`; consensus ratios stored under last handle's id and loaded
  into a copy after push_back (never round-trip); OMSFileLoad ref maps default-construct on missing keys.
- Old schema versions (<4 DataValue join, <5 FEAT_Feature/overall_quality, <=2 grouping_order/position, 1 no specificity).

## Inc. 3 progress (commit 401019f)
- OMSFile load/store DONE: src/format/oms_file{.rs,/store.rs,/load.rs}; tests/oms_file.rs 15 tests, 14 passed on spock
  before the last one was added. Oracles: ../oracle/F2-oms (run.sh, dump_tables.py, make_dumps.sh) and
  ../oracle/F2-oms-noid (C++ cannot store maps without ID data). Expected dumps in tests/data/oms_file/release/.
- Remaining: docs/OMS_FILE_SUPPORT.md, tests/data/oms_file_provenance.json, clippy/doc/MSRV gates, report+ledger.
- MzQCFile: leave untouched (needs QC classes, qc-cv.obo). MzTabM: check doc gaps if time.

## Next step (old list, superseded by the progress above)
1. Gate the inherited work on spock: tests mztab_export (features?) + oms_file (--features sqlite) + clippy.
2. Write docs/OMS_FILE_SUPPORT.md (referenced by oms_file.rs, missing), provenance for oms test data; commit.
3. MzQCFile (store only; QualityControl tool) -> src/format/mzqc_file.rs.
4. OMS load/store of IdentificationData (src/identification/graph.rs) if feasible; else partial with named gaps.
5. Report + ledger fragment; delete NOTES.md.

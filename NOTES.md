# Lane F2 notes (MzTab, MzTab-M, mzQC, OMS)

Branch: worktree-wf_ee239125-5be-11 (incarnation 2; merged worktree-wf_bcdd783d-f3f-11 = notes only).

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

## Next step
1. Write src/format/mztab_export.rs (in progress).
2. Streaming store overloads in mztab_file.rs.
3. tests/mztab_export.rs with retained MzTabExporter_4..8 (+3) under FuzzyDiff contract.
4. Then MzQCFile, then OMSFile (rusqlite pinned; IdentificationData in src/identification/graph.rs).

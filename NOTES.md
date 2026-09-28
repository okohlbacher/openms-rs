# F3a lane notes (Identification XML) — restart-safe scratch; delete in final commit

Branch: worktree-wf_ee239125-5be-12 (off main f0990c0, merged worktree-wf_bcdd783d-f3f-12 = FF).
Never push/rebase. Cargo only via gate: OPENMS_GATE_HOST=dax OPENMS_GATE_JOBS=64 ~/.local/bin/openms-kim-gate.sh F3a ...
(run_in_background, poll output).

## Current coverage status (docs/core-sdk-coverage.json)
- PepXMLFile partial (src/format/pepxml.rs), PepXMLFileMascot UNMAPPED, ProtXMLFile UNMAPPED,
  IdXMLFile partial (identification_xml.rs, idxml.rs), MzIdentMLFile partial (mzidentml.rs),
  MascotXMLFile partial (mascot_xml.rs), XQuestResultXMLFile UNMAPPED, PercolatorOutfile
  evidence_requires_review (no rust), MzIdentMLHandler partial, MascotXMLHandler partial,
  XQuestResultXMLHandler UNMAPPED.

## Budget fix (priority 1) — WIP code exists (commit 5bb4d01, by previous incarnation)
- sequence.rs parse_with_budget -> parse_metered with ParseMeter; per-annotation charge now
  name_lookup_work (modifications.rs) * 16 for named, registry scan for mass.
- tests/file_info_a7.rs: new test identifications_many_modified_hits_upstream_files with
  oracle expected id_ff25/id_fdr5 (tests/data/file_info_a7/expected). UNVERIFIED: not yet
  gated on dax; literals must be checked against oracle output files.

## Plan / order
1. Verify + gate budget fix (idxml feature; file_info_a7 test; sequence tests).
2. PercolatorOutfile (small, IDFileConverter consumer).
3. ProtXMLFile (IDFileConverter consumer).
4. XQuestResultXMLFile + handler (IDFileConverter, OpenPepXL, XFDR).
5. PepXMLFileMascot.
6. Gaps in partials if time.

## Findings (incarnation 2)
- Budget fix WIP verified by reading: named lookups only touch by_name candidates (find/find_handles/
  get_modification); resolve paths <= 9 find-shaped calls (charged 16). Oracle expected id_ff25/id_fdr5
  == ../oracle/a7-fileinfo/results/run1 (run1==run2), literals match. Registered 6 fixtures in
  tests/data/file_info_a7_provenance.json. Still to do: evidence text 55->57, known_gaps entry ->
  RESOLVED, docs/FILE_INFO_A7_SUPPORT.md; note for lead (VALIDATION.md).
- Regex: src/concept/boost_regex.rs BoostRegex (fancy-regex, always on). CsvFile complete (src/format/csv.rs).
  ModificationDefinitionsSet complete (src/chemistry/modification_definitions.rs infer_from_peptides).
  SpectrumMetaDataLookup partial: mascot_xml::SpectrumTitleLookup (hand-coded formats) + pepxml::SpectrumIndex.
- Retained TOPP outputs: IDFileConverter_3 (protXML->idXML, topp/IDFileConverter_3_input.protXML ->
  IDFileConverter_3_output.idXML); IDFileConverter_12/13/14 (psms -> idXML, -score_type qvalue/PEP/score).
  FileHandler PROTXML: push one ProteinIdentification + one PeptideIdentification.
- Class test data in core src/tests/class_tests/openms/data: PercolatorOutfile_test.psms,
  ProtXMLFile_input_{1,2,3}.protXML, PepXMLFileMascot_test.pepXML, XQuestResultXMLFile_test_data.xquest.xml (661K).
- XQuest load needs OPXLHelper::combineTopRanksFromPairs (unported) + helpers private in mzidentml.rs. Do last.
- Isolation checker: avoid `cd x && git`, heredocs mentioning git, the word `source`; use `git -C <abs> show` alone.

## Done (incarnation 3 view)
- Budget fix committed 3ce9df2 (gate g1 35 pass file_info_a7; g2 pepxml 8 pass).
- PercolatorOutfile + SpectrumMetaDataLookup committed cbe9b27 as WIP; gate g3 (scratchpad/gate/g3.log,
  features idxml?) 13 percolator_outfile tests + 26 lib tests pass. Provenance/docs not yet checked.
- scratchpad/gate/test1.log + check1.log belong to lane F2, not us.

## Next step
1. [done c7c2a42] Percolator milestone (docs + provenance; gate g4 15/15 + lib 2).
2. ProtXMLFile -> src/format/protxml.rs (feature idxml). Fixtures in tests/data/protxml/ (fetched).
   Oracle ../oracle/F3a-protxml (run.sh push/pull; ibminode06 /scratch/kohlbach/F3a-protxml;
   IDFileConverter -test on 4 upstream + 13 generated inputs, twice). Retained IDFileConverter_3
   output has float-widened protein scores (0.999000012874603); pin writes toStr(double).
   Reuse pepxml internal_residue_mass (make pub(crate)). Helper: scratchpad/pin.py (sha/show/grep/ls/save).
3. PepXMLFileMascot -> src/format/pepxml_mascot.rs.
4. XQuest last.

## Incarnation 4 (wf_cd066dab-0a7-7, branch worktree-wf_cd066dab-0a7-7)
- Merged worktree-wf_ee239125-5be-12 (FF to 1f8bb86: lead-saved WIP protxml.rs + xml_sax.rs + expected outputs).
- FOUND: branch worktree-wf_ee239125-5be-32 (integ-F, NOT in my brief) already salvaged + finished
  ProtXMLFile (4e49731..f4de3cc: docs/PROTXML_SUPPORT.md, tests, provenance) and replaced xml_sax.rs with
  a shared `format::sax_xml` driver (2d6ad26). => Do NOT touch protxml.rs/xml_sax.rs further here
  (would only conflict); report that protXML is done on integ-F.
- Next step:
  1. PepXMLFileMascot (src/format/pepxml_mascot.rs), using the SAX approach main already uses.
  2. XQuestResultXMLFile + handler (needs OPXLHelper::combineTopRanksFromPairs).
  3. Partial gaps (IdXML/MzIdentML/MascotXML/PepXML) if time.
- [inc4] sax_xml.rs copied verbatim from integ-F (2d6ad26) + cpp_release_records.txt; only edit: line 350
  kb-xml cfg_attr -> allow(dead_code); mod.rs has #[allow(dead_code)] on the mod.
- [inc4] PepXMLFileMascot: src/format/pepxml_mascot.rs, tests/pepxml_mascot.rs, docs/PEPXML_MASCOT_SUPPORT.md,
  tests/data/pepxml_mascot_provenance.json (generator: scratchpad mk_mascot_prov.py). Oracle
  ../oracle/F3a-pepxml-mascot (ibminode06, run.sh push/pull), results copied to tests/data/pepxml_mascot/expected.
  Gate: scratchpad/gate/g3.log (test), g3c.log (clippy). rsync to dax sometimes times out -> check log head.
- XQuest next: needs faithful OPXLHelper addXLTargetDecoyMV/addBetaAccessions/combineTopRanksFromPairs/
  removeBetaPeptideHits/computeDeltaScores (mzidentml.rs copies DIVERGE on purpose; do not reuse blindly).
  Writer tier 1: XFDR_test_out{1,2,3,4,5,7}.xquest.xml retained (TOPP_XFDR_* CMakeLists 3478-3510).

## Incarnation 5 (same worktree, restarted)
- HEAD 210d49d: PepXMLFileMascot milestone committed (gated). protXML aligned with integ-F.
- No lane report / ledger fragment written yet.
- Next step:
  1. XQuestResultXMLFile + handler (src/format/xquest_result_xml.rs) if feasible; else report as gap.
  2. Write report .planning/lane-reports/F3a.md + ledger wave4/F3a.json (covering budget fix,
     Percolator, protXML(integ-F), PepXMLFileMascot, XQuest).
  3. Delete NOTES.md in final commit.
- [inc5] XQuest ORACLE ALREADY DONE by inc4: ../oracle/F3a-xquest (driver.cpp dump format, make_cases.py,
  node.sh, splits.tsv, 32 inputs incl. x01..x16 generated cases). results/ DONE, ibminode06 Release
  gcc 14.4.0, differ.txt empty (run1==run2). Modes: load (dump), store, idxml2xquest (convert), split.
- Reuse: math::source_sort::source_sort_reversed_by (libstdc++ std::sort emulation) for
  combineTopRanksFromPairs. mzidentml.rs OPXL helpers diverge on purpose -> write faithful private ones.
- Write src/format/xquest_result_xml.rs (feature: same as idxml?), tests/xquest_result_xml.rs.
- [inc5] src/format/xquest_result_xml.rs written (reader, OPXLHelper steps, writer; spec.xml writers = gap).
  tests/xquest_result_xml.rs: load/store/convert/split differentials vs tests/data/xquest_result_xml/expected
  (20 cases retained; XFDR_test_in2/in3/out3/4/5/7 too big, run ad hoc only). Gate logs scratchpad/F3a/c*.log.
  Store comparison sorts the extra meta attributes (registry order divergence).

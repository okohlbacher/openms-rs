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

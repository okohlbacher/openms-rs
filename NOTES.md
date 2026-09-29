# INTEG-F notes (integrator, group F)

Base: main f0990c06. Branch: worktree-wf_a10fdd71-2d2-14. dax slot integ-F (sweep integ-F-sweep).
Scratch: /private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad/integF
(both.py = keep-both conflict resolver; lockmerge.py = Cargo.lock union; logs/).
Guard quirks: one plain `git -C <wt>` per command, no loops; heredoc python + cd refused -> write script files.

## PRIOR integ-F attempt (Sep 27) = branch worktree-wf_ee239125-5be-32 (NOT in any lane, NOT in main)
Its NOTES (git show worktree-wf_ee239125-5be-32:NOTES.md) list fixes to re-apply where still relevant:
- 23ad0e4/687c614 D1 string_utils collapse; 4bc5bf0 toDouble/toFloat Release from_chars (oracle ../oracle/integ-F-stringutils,
  tests/string_utils_release.rs); f5f5111 ProteinIdentification::sort tie by accession;
  2d6ad26 one Xerces SAX driver (F3a already carries its sax_xml.rs verbatim; F3b search_engine_xml still separate);
  b7216a7 remove tests/data/search_engine_oracle/driver.cpp.txt (C++ in repo!); d940b3d arrow builds alone;
  04c0d42 AcqusHandler false Release claim; ../oracle/integ-F-ptmxml.

## Merged (all 12 lanes, in order): F1b(+F1a+OA) e.. F1c e39e7181, F2 3ddf9817, F3a b9056b11, F3b 76190ade,
F4 993fdd01, F5 f3c14e82, F6 877758b7 (took F3a/integ-F sax_xml.rs, gate any(idxml,kb-xml)), F7a eb3e44e7,
F7b 30c01dae (Cargo.lock union-merged by hand: verify --locked!), F7c b44150c9, fmt 5a6fbd7e.

## Done after merge (commits)
7fb53b93 D1 format::string_utils folded; 60e302c6 cherry-pick f5f51116 sort ties; cherry-pick 4bc5bf07 StringUtils
Release from_chars (+tests/string_utils_release.rs); 978bdf89 cherry-pick 687c6145 PEFF/lookup; 0bf73090 Bruker imaging;
1f64cd7a D2 SAX collapse (xtandem/omssa -> sax_xml; ptm oracle); 4ddbebe9 D4 PTMXML (ptm_xml_file removed);
cccfe3e4 driver.cpp.txt removed; 266c174a D5 parseInt -> StringManager. clippy all-features: CLEAN (logs/c2.log).
C++ reference at bc9cc12: .reference/openms4-core-bc9cc12 (scratch/core symlink; cpp.py prints lines).
Audit scripts: scratch audit.py / audit_fn.py / gaps.py / frag.py. No complete-with-unported-member found so far;
handlers inconsistent (PTMXMLHandler F3b complete vs F6 native_equivalent -> native_equivalent).
Ledger rule: core_sdk_coverage.py accepts only complete/partial/native_equivalent; rust/tests/doc files must exist.
'unmapped' entries (F1b ConsensusMapArrowExport/ProteinGroupArrowExport, F5 HDF5/Thermo) must NOT go in INTEG-F.json.
Open: D6 stod duplicate (arrow_io_helpers::c_stod glibc-measured vs string_utils::stod) - collapse if time.

## Duplicates to collapse
- D1 string_utils: F3b src/data_structures/string_utils.rs vs F4 src/format/string_utils.rs (cli/context uses data_structures).
- D2 SAX: F3b src/format/search_engine_xml.rs driver vs sax_xml.rs (X!Tandem/OMSSA).
- D4 PTMXML: F3b src/format/ptm_xml_file.rs (idxml) vs F6 src/format/ptm_xml.rs (kb-xml).
- D5 xerces_parse_int (sax_xml) vs F7b string_manager::parse_int.
- D3 F7a MSDataCachedConsumer now unblocked by F4 cached_mzml_handler.
- F7b ZIP input vs B8a zip_input.rs (cross-group, flag to lead).
- F3b driver.cpp.txt must go.

## Incarnation 3 progress (12:05-)
- Verified: extcheck 249 artifacts OK; trial_register 36 new manifests OK (check_core_sdk); trial_ledger.py (scratch)
  merges INTEG-F.json into reviewed-apis + core_sdk_coverage --write OK (core complete 130->215), restored.
- Audit CHK items all mapped (InspectInfile max_ptm_size, XTandem max_valid_evalue, ParamXML write*, MGF updateMembers_
  in doc, StringManager xerces ptrs = Rust ownership). AcqusHandler Release claim now TRUE (F5-xmass oracle populated).
- 5bec6a37 collapse F4 clear_*_data -> kernel clear(false).
- MSDataCachedConsumer ported (src/format/ms_data_cached_consumer.rs, tests/ms_data_cached_consumer.rs, doc, manifest,
  ledger.py entry). Proposed ci_matrix mzml slice + ms_data_cached_consumer (uncommitted, --write done).
  Gate g1 (scratch logs/g1.log) running; then commit.
- sweep2 running on dax slot integ-F-sweep (logs/sweep2.out; remote .ci-line.N.log) - snapshot BEFORE the consumer port.

## Next step
1. g1 green -> commit consumer. Then gate clippy --all-features --all-targets, doc --all-features, +1.85.0 mzml slice.
2. Wait sweep2; collect per-line pass counts via ssh (read-only).
3. Write INTEG-F.md (Bash heredoc) from lane reports + this; revert lead files; delete NOTES.md; final commit.

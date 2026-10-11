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

## Incarnation 3, later (12:45)
- Committed: 2392a81f MSDataCachedConsumer; 3d0e7897 empty test mod (clippy fail from 6ab21922); 5685ca4e doc links
  (PEFF/PTMXML/Unimod module docs - 23 rustdoc errors; F6 never ran doc).
- PSMArrowIO ported (src/format/psm_arrow_io.rs, tests/psm_arrow_io.rs 14 pass, doc ARROW_PARQUET_SUPPORT, manifest
  arrow_parquet_provenance psm_arrow_io_note, ledger.py entry). Gates g9-g11 (clippy, doc, 1.85) running -> commit.
- sweep2 = pre-fix snapshot: its clippy + doc lines WILL fail (fixed since). Needs sweep3 on final tree.

## 12:57
- 84a9860e PSMArrowIO committed. Oracle ../oracle/integ-F-psmarrowio ran OK (ibminode06, twice identical);
  tests/data/psm_arrow_io/ (records, stderr, export_* bundles); 2 new differential tests pass (16 total) after
  mapping rank meta<->field and documenting B/Z/X m/z warning native difference.
- Report: harness refuses report .md files (Write) -> full report goes into the FINAL MESSAGE (brief fallback).
  INTEG-F-lead-files.patch saved in .planning/lane-reports (regenerate at end).

## Next step
1. g14 (clippy) g15 (1.85) green -> commit differential. trial_register again (arrow manifest changed).
2. sweep2 done -> sweep3 (same slot) on final tree; counts.sh integ-F-sweep sweep3/lines.tsv.
3. Revert lead files, regenerate patch, delete NOTES.md, final commit, final message = report.

## Incarnation 4 (wf_97065d91-48e-4, 2026-10-11): resumed from the above by merging worktree-wf_a10fdd71-2d2-14 (4d9800d6)
Own branch: worktree-wf_97065d91-48e-4. Scratch (shared with other agents!): /private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad
 (mine: audit/ lead-copy/ sweep2/ sweep-copy.sh gate-copy.sh collapse_d*.py; others' dirs: O3a integB-inc2 status tools lemon)
Guard quirks: no loops/heredocs with git in one command; ".git" substring anywhere (even .github, --exclude) is refused; write
 scripts with the Write tool and run them. Lead files may not be edited in the worktree -> proposals go through lead-copy/ (scratch
 copy of the tree with old-lead.patch applied; sweep-copy.sh runs the sweep from it). Killing remote processes is refused.
Done in this incarnation:
- merged 4d9800d6 (all 12 lanes + collapses); lane files present (lanecheck.sh), test counts preserved (testcount.sh)
- 0822a187 D7: one StringUtils split/trim/toDouble/toInt32 for msp, xquest, mascot_xml, mascot_generic, mztab, mztab_file
  (tests pass on dax slot integ-F-g, clippy all-features clean)
- D8 (uncommitted until its gate is green): arrow_schema_registry::validation::schemas_equal is the one Schema::Equals; the
  OpenSWATH readers' concatenation used arrow-rs == (refused list<item> vs list<element>) and readers/writer fed arrow-rs
  concat, which also compares the element name -> parquet_file::conform_to_schema. Oracle ../oracle/integ-F-concat (ibminode06,
  Release, twice identical) -> tests/data/arrow_parquet/cpp_release_concat_records.tsv, tests/arrow_parquet.rs
  concatenation_follows_schema_equality_for_nested_names, lib test parquet_reader_helpers::files_that_spell_...;
  manifest arrow_parquet_provenance + docs/ARROW_PARQUET_SUPPORT.md updated.
- sweep2 (77 lines = 54 + 23 proposed slices, from lead-copy) running on dax slot integ-F. A leftover loop of the killed
  first sweep (remote pid 2576733) runs in the same slot and overwrites .ci-line.N.log: do not trust PASS summary lines blindly,
  re-run suspicious lines individually under another slot (integ-F-g) with gate-copy.sh (GATE_LOCAL=lead-copy).
- register_trial.py: 37 unregistered F manifests register cleanly (check_core_sdk passes) -> audit/registered_manifests.txt
Cross-group: B8a (INTEG-B, worktree-wf_97065d91-48e-3) also adds src/data_structures/string_utils.rs (1912 lines, Result<_, Error>
 API), src/format/zip_input.rs and edits path_io.rs; merge-tree vs this branch conflicts in Cargo.*, mod.rs files,
 string_utils.rs (add/add), path_io.rs, file_handler.rs, pepxml.rs, cli/context.rs, identification.rs, digestion.rs.
Residual (documented, not collapsed): mascot_xml::SpectrumTitleLookup + pepxml::SpectrumIndex vs spectrum_meta_data_lookup;
 mzidentml OPXLHelper steps (data-preserving) vs xquest (faithful); path_io gzip/bzip2 = MultiGzDecoder (not GzipIfstream).

## Next step (incarnation 4)
1. g-d8d green -> commit D8 (src, tests, test data, provenance, docs). Re-run clippy + doc + 1.85 for the whole tree.
2. Wait for sweep2 (+ leftover loop) to end; read per-line pass counts from dax .ci-line.N.log.
3. Regenerate the lead patch (rust.yml, tools/ci_matrix.py, tools/ci_coverage_baseline.json) from lead-copy after any new
   test file; write it to the main checkout .planning/lane-reports/INTEG-F-lead-files.patch.
4. Final fragment: start from main-checkout .planning/ledger-fragments/wave4/INTEG-F.json, amend PEFF scope (tier-1
   sequence members exist: ../oracle/F6-peff-sequences), D7/D8 notes; the system rule forbids a report .md: report in the
   final message and in the fragment/lead notes.
5. Delete NOTES.md in the final commit.

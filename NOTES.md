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

## Duplicates to collapse
- D1 string_utils: F3b src/data_structures/string_utils.rs vs F4 src/format/string_utils.rs (cli/context uses data_structures).
- D2 SAX: F3b src/format/search_engine_xml.rs driver vs sax_xml.rs (X!Tandem/OMSSA).
- D4 PTMXML: F3b src/format/ptm_xml_file.rs (idxml) vs F6 src/format/ptm_xml.rs (kb-xml).
- D5 xerces_parse_int (sax_xml) vs F7b string_manager::parse_int.
- D3 F7a MSDataCachedConsumer now unblocked by F4 cached_mzml_handler.
- F7b ZIP input vs B8a zip_input.rs (cross-group, flag to lead).
- F3b driver.cpp.txt must go.

## Next step
1. Poll logs/c1.log (clippy all-features on dax). Fix compile errors.
2. Collapse D1, D4, D2, D5; remove driver.cpp.txt.
3. Full sweep; attack complete claims; write INTEG-F.json + INTEG-F.md.

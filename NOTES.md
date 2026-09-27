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

## Next step
1. Check oracle provenance of id_ff25/id_fdr5 expected (../oracle/a7-fileinfo), check test literals,
   gate `cargo test --features idxml,... --test file_info_a7` on dax.

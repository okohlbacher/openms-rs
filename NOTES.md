# F7c NOTES (lane working notes; delete in final commit)

Branch: worktree-wf_cd066dab-0a7-19 (retry; merged -14 notes).

## Ledger gaps as recorded (docs/core-sdk-coverage.json)
- DTA2DFile (src/format/dta2d.rs): MISSING DocumentIdentifier loaded-path/type state; C++ stream numeric text (Rust shortest float).
- DTAFile (src/format/dta.rs): no whole-header completion claim; default_ms_level_ protected.
- MS2File (src/format/ms2.rs): MapType template param fixed to MSExperiment (arguably native).
- MzDataFile (src/format/mzdata.rs): isSemanticallyValid unsupported; FileHandler dispatch (B0, not ours).
- MzXMLFile (src/format/mzxml.rs): transform overloads report no progress; truncated-after-</scan> accepted where C++ ParseError.
- MascotGenericFile (src/format/mascot_generic.rs, mgf.rs): "generic base/template/API behavior should be reviewed".
- ConsensusXMLFile (src/format/consensusxml.rs): source ZIP input unsupported; full passive PeakFileOptions separate.
- TransformationXMLFile (src/format/transformation_xml.rs): b_spline/polynomial/akima => Error::Unsupported.

## Next step
1. Read each C++ header at pin (git -C ../OpenMS4-tests/packages/core show bc9cc12:src/openms/include/OpenMS/FORMAT/X.h) member by member vs Rust.
2. Pick tractable gaps: MzXML truncated doc ParseError; TrafoXML model resolution (depends on transformations models in src/analysis/transformations/); DTA2D numeric text; DocumentIdentifier.

## Incarnation 3 (wf_a10fdd71-2d2-6) findings
- Merged worktree-wf_cd066dab-0a7-19 (fast-forward, notes only).
- TrafoXML b_spline/akima: ALREADY CLOSED on lane M1 branch worktree-wf_cd066dab-0a7-36
  (commits 48b83c7d, c060e480: src/format/transformation_xml.rs + tests/transformation_xml.rs;
  M1 says `polynomial` is not a source interpolation type). DO NOT edit transformation_xml.rs
  (would conflict with M1). Only do the member-by-member check vs M1's version and report.
- Order of work: MzXML truncated-doc ParseError (FileConverter consumer) -> DTA2D numeric text
  -> DocumentIdentifier state -> MzData isSemanticallyValid -> MGF review -> MS2 template -> ConsensusXML zip.
- MzXML truncated-after-scan: ALREADY FIXED (docs/FILE_INFO_A8_SUPPORT.md s.6, repair F2). Remaining MzXML gap:
  transform overloads report no progress (needs oracle: extend ../oracle/progress-format-readers driver pattern
  into ../oracle/F7c-transform-progress/). FileHandler dispatch is B0's, not an MzXMLFile member.
- DTA2D tier-1 bytes ALREADY RETAINED: ../oracle/progress-format-readers/results/out_1/dta2d_store_cmd.dta2d and
  dta2d_store_tic_cmd.dta2d (ibminode06 Release, 3 reps identical). C++ store = precisionWrapper -> StringUtils::toStr:
  RT/mz via text_format::to_str, intensity via to_str_f32. tests/progress_format_readers.rs:634 exempts their sizes
  (remove exemption once bytes match). Class test DTA2DFile_test.cpp:72-73 + MS2File_test.cpp:49-50 give
  DocumentIdentifier expectations (loaded path, type name dta2d/ms2).
- Plan DTA2D: WriteOptions -> struct {limits, discard_unrepresentable}; strict default refuses text that does not
  read back exactly; source() = C++ behaviour. MS2: record DocumentIdentifier, doc ms2.rs to 100%.

## Done (incarnation 3)
- 9b55dd03 DTA2D source number text (tier 1 vs retained Release bytes) + DocumentIdentifier for DTA2D/MS2.
- 629a82cf DTAFile path load/store, name=basename, source line split, DEFAULT_MS_LEVEL; class test ported.
- Oracle ../oracle/F7c-transform-progress (driver.cpp, run.sh start|poll|fetch) launched on ibminode06
  for MzXMLFile::transform progress (8 cases + filtered). Next: fetch, add transform_with_progress to mzxml.rs.
- Gate cmd that works: OPENMS_GATE_HOST=spock OPENMS_GATE_JOBS=64 ~/.local/bin/openms-kim-gate.sh F7c test --locked --all-features --test <t> ...
  (no shell variables in the command; the sandbox refuses them). Edit scripts go in scratchpad as .py files.
- Do NOT run check_doc_coverage.py --write (docs/doc-coverage.json is the lead's).

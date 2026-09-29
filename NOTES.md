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

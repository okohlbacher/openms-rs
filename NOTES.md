# F7c NOTES (lane working notes; delete in final commit)

Branch: worktree-wf_cd066dab-0a7-14, base main f0990c0.

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

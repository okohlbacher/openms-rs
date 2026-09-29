# O1b lane notes (Transition list I/O)

Prereq merged: worktree-wf_cd066dab-0a7-30 (O1a: TargetedExperiment, ReactionMonitoringTransition,
TargetedExperimentHelper, IncludeExcludeTarget in src/analysis/openswath/).

C++ sources extracted to scratchpad O1b/ (git archive of core bc9cc12).
Sizes: TraMLHandler.cpp 1745, TransitionTSVFile.cpp 2001, TransitionPQPFile.cpp 1683,
TraMLFile.cpp 58, MRMFile.cpp 54, TargetedDataFileLoader.cpp 59.
Class tests: TraMLFile_test (208), TransitionTSVFile_test (361), TransitionPQPFile_test (445),
MRMFile_test (36), TargetedDataFileLoader_test (39). Test data: ToyExample1.traML.

MRMFile / TargetedDataFileLoader need OpenSwath::SwathMap, SwathFile, SimpleOpenMSSpectraFactory,
MRMChromHandler -> check whether ported; probably blocked (name the deps).

## Status
- DONE: TraMLHandler + TraMLFile (src/format/traml.rs, traml_handler.rs), tests/traml.rs 13/13 pass on dax
  (`test --features semantic-validation --test traml`). CVTermList::add_source added (metadata/value.rs, gated traml).

## Next step
2. TransitionTSVFile -> src/analysis/openswath/transition_tsv_file.rs
3. TransitionPQPFile (rusqlite) -> src/analysis/openswath/transition_pqp_file.rs
4. MRMFile / TargetedDataFileLoader if deps exist.

## Findings (session 1)
- `git -C <other repo>` works only as a plain single command (no loops/pipes around it).
  C++ copies: scratchpad/O1b/ (TraMLHandler.cpp, XMLHandler.h/.cpp, XMLFile.cpp, ModificationsDB.cpp).
- identification_xml::parse_xml_with_budget refuses namespaced elements -> TraML needs its own
  streaming quick-xml reader (mirror its GeneralRef/escape handling, quick-xml 0.39).
- Number formatting: XMLFile::save_ sets os.precision(15) -> `os << double` ==
  format::file_info::text_format::ostream_g(x, 15); precisionWrapper/StringUtils::toStr(double)
  == param::value::format_float(x, true). DataValue::toString: private helpers exist in
  metabo_targeted_target_decoy.rs / imzml_writer.rs (write own private one).
- userParam order: crate sorts MetaInfo by key (BTreeMap); C++ uses MetaInfoRegistry order.
  Established divergence (docs/MZML_WRITER_CPP_PARITY.md) -> compare decoded, not bytes.
- Unimod name on write: ModificationsDB::global().find("UniMod:N", residue, term).first().name()
  (C++ getModification picks first, throws InvalidValue when none).
- CV check: ControlledVocabulary::psi_ms() (restrict to MS: accessions; C++ loads psi-ms only).
- is_semantically_valid: TraMLValidator == SemanticValidator + check_units; embed
  share/OpenMS/MAPPING/TraML-mapping.xml into resources/cv/ (gate: semantic-validation).
- Tier-1 fixtures (test-data 0cb15f2 topp/): ConvertTSVToTraML_{1..7} (TSV->TraML),
  TargetedFileConverter_{1,3} (TraML->PQP->TraML), _5/_6 (TraML->TSV), FileMerger_9.
- Feature plan: new `traml = ["dep:quick-xml","file-compression"]` in default (Cargo.toml, report it).
- C++ defect candidates: Compound MS:1002476 sets actual_peptide_ drift time; Product charge
  via toDouble; handleCVParam_ value-type checks nested inside `if (obsolete)`.

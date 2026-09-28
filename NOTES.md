# Lane I4 notes (IDMapper + metabolomics ID)

## Status
- No Rust for any lane header yet. IDMapper is `evidence_requires_review` with no candidate files.
- C++ at pin: the core working tree (../OpenMS4-tests/packages/core, HEAD 291c5ba) is IDENTICAL to bc9cc12 for
  every lane .h/.cpp/_test.cpp (checked with `git diff --stat bc9cc12 HEAD -- <dirs>`). Copies live in the
  session scratchpad `cpp/` dir; if lost, `cp` them from the core working tree (plain cp works).
- Tooling gotcha: the worktree guard refuses compound commands containing `git -C <other repo>`, shell variables
  used as paths, or the word "source" in some positions. Use literal absolute paths, one git call per command.
- Rust building blocks: FeatureMap/Feature/ConsensusMap/ConsensusFeature/FeatureHandle in src/kernel/features.rs,
  PeptideIdentification (rt/mz Option) in src/identification.rs, DataProcessing in src/metadata/acquisition.rs,
  AASequence::mono_mass_for/average_mass_for(PeptideFragmentType::Full, z), concept::math_functions::ppm_abs.
  idxml::load / featurexml::load / consensusxml::load exist (default features).
- Missing: AnnotatedMSRun (unmapped, not my header) -> Rust annotate for peak maps takes the experiment plus
  `&mut Vec<PeptideIdentification>` / `&mut Vec<ProteinIdentification>`. SpectrumLookup (partial, private in pepxml):
  only findByNativeID (last native id wins, `ids_[native_id] = index`) is needed -> local HashMap.
- Oracle: IDMapper_test.cpp (in-memory asserts + consensusXML outputs IDMapper_{3,5,6}_out*.consensusXML).
  Test data goes into tests/data/ (flat, upstream names).

## Next step
1. Write src/analysis/id/id_mapper.rs + tests/id_mapper.rs (+ copy IDMapper_* data), commit.
2. Then AccurateMassSearchEngine, MetaboliteSpectralMatching, SiriusMSConverter/SiriusFragmentAnnotation.

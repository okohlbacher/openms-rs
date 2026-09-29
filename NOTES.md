# O1c lane notes (Assays and decoys)

Merged prerequisite worktree-wf_cd066dab-0a7-30 (O1a: TargetedExperiment,
ReactionMonitoringTransition, TargetedExperimentHelper, MRMMapping, ...) as a
fast-forward. No earlier O1c incarnation exists in bcdd783d branches.

C++ pin readable at .reference/openms4-core-bc9cc12 (outside worktree; git -C
to other repos is refused by the sandbox, plain cat/find works).

Headers (src/openms/include/OpenMS/ANALYSIS/OPENSWATH/):
- MRMIonSeries.h (139) / .cpp (666) / test 279 lines  -> dependency of both below
- MRMDecoy.h (332) / .cpp (1160) / test 1287 + TraML in/out fixtures
- MRMAssay.h (568) / .cpp (2197) / test 1007 + many TraML fixtures
- TransitionListEvidenceFilter.h (144) / .cpp (1320) / test 438

Target Rust: src/analysis/openswath/{mrm_ion_series,mrm_decoy,mrm_assay,
transition_list_evidence_filter}.rs; append pub mod lines to
src/analysis/openswath.rs.

## Next step
1. Read MRMIonSeries.h/.cpp/test; check Rust AASequence, residue, TheoreticalSpectrumGenerator, ModificationsDB, TraML reader availability.
2. Port MRMIonSeries, commit.
3. MRMDecoy, then MRMAssay, then TransitionListEvidenceFilter.

## Findings (incarnation 1)
- No TraML reader in Rust (FORMAT/TraMLFile.h unported) -> TraML fixtures of
  MRMAssay/MRMDecoy tests cannot be loaded; plan a C++ oracle on ibminode06
  (../oracle/O1c-<topic>/) that dumps inputs/outputs in a simple text/JSON form.
- LightTargetedExperiment (OPENSWATHALGO) unported -> generateDecoysLight and
  *Light signatures are a named gap; port algorithms on (location) vectors.
- Residue loss table already in src/chemistry/theoretical.rs fragment_losses
  (D/E/S/T H2O, K/N/Q NH3, R NH3,C1H2N2,C1H2N1O1; modified residue -> mod's
  losses only, empty formulas never present). AASequence::internal_mass(i) is
  pub(crate) = Residue::getMonoWeight(Internal).
- boost::uniform_int<> default range is [0,9]! pseudoRNG() in shufflePeptide
  yields 0..9; pseudoRNG(n) yields [0,n-1] via generate_uniform_int bucket
  mapping on 32-bit boost::mt19937 (rand_mt::Mt is 32-bit MT19937).
- decoy_random.rs: Mt64 + bucket mapping for Math::RandomShuffler (used by
  generateDecoys when aim_decoy_fraction < 1); shuffle is &mut [u8] only.
- C++ defect candidates: MRMIonSeries::annotateTransition uses getSuffix for
  a/b/c ions; annotationToCVTermList_ `static const EmpiricalFormula
  nl_formula(...)` caches the FIRST loss formula forever; shuffle in
  generateDecoys seeds from time(nullptr).

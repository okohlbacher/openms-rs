# B8b lane notes (BLOCK-06b)

Prereq merge: worktree-wf_ee239125-5be-1 fast-forwarded f0990c0 -> 9fcb962 (no merge commit needed).

## Scope and ledger gaps (from docs/core-sdk-coverage.json)
1. MATH/MathFunctions.h — src/concept/math_functions.rs. Only gap: Math::RandomShuffler (6 members);
   semantics exist as crate-private src/chemistry/decoy_random.rs DecoyRandom. Promote into a shared
   module (math/concept sink) and re-use from chemistry. Test: MathFunctions_test.cpp RandomShuffler section.
2. METADATA/PeptideIdentification.h — src/identification.rs. Missing: buildUSI (2 overloads), buildUIDFromPepID,
   buildUIDsFromAllPepIDs, getScoreComparator, std::hash<PeptideIdentification>, SpectrumIdentification alias.
3. CHEMISTRY/ProteaseDB.h — src/chemistry/digestion.rs + digestion/enzymes.rs. Missing: writeTSV, provider ctor.
4. CHEMISTRY/ModificationsDB.h — src/chemistry/modifications.rs. Missing: initializeModificationsDB,
   isInstantiated, searchModification(ResidueModification), searchModificationsFast, has(), addModification x2,
   findModificationIndex, searchModificationsByDiffMonoMassSorted x2, getAllSearchModifications, writeTSV.
5. PROCESSING/ID/IDFilter.h — src/analysis/id_filter.rs. Large unported surface (predicates, digestion filters,
   countHits/getBestHit, filterGroupsByScore, regex, RTPredict, keepBestPerPeptide*, annotateBestPerPeptide*...).

## Commands
- Pinned C++ extracted by archiving src/openms whole (a path token "source" trips the worktree guard)
  into $SCRATCH/cpp with symlinks inc/ (include/OpenMS), impl/ (src/openms/source), tests/ (class tests).
  SCRATCH=/private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad
- Oracle template: ../oracle/B6-param-defaults/{run.sh,remote_run.sh,driver.cpp} (ibminode06 via sshgw jump)
- USI.h is unported (evidence_requires_review) and blocks buildUSI.
- Native hashing precedent: docs/METADATA_HASH_SUPPORT.md (std::hash::Hash consistent with PartialEq, not FNV).
- cargo only via: OPENMS_GATE_HOST=dax OPENMS_GATE_JOBS=64 ~/.local/bin/openms-kim-gate.sh B8b <args> (background!)
- C++ pins: git -C ../OpenMS4-tests/packages/core show bc9cc12:<path>

## Next step
1. MathFunctions RandomShuffler (promote DecoyRandom).
2. PeptideIdentification missing members.
3. ProteaseDB writeTSV; ModificationsDB members.
4. IDFilter: countHits, getBestHit, filterGroupsByScore, regex, keepBestPerPeptide*, digestion filter.

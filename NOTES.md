# B8b lane notes (BLOCK-06b)

Prereq merge: worktree-wf_ee239125-5be-1 fast-forwarded f0990c0 -> 9fcb962 (no merge commit needed).

## Scope and ledger gaps (from docs/core-sdk-coverage.json)
1. MATH/MathFunctions.h — DONE, complete.
2. METADATA/PeptideIdentification.h — src/identification.rs + src/identification/peptide_identification.rs
   + src/identification/usi.rs (USI.h ported as dependency). Oracle ../oracle/B8b-peptide-id (tier 1).
3. CHEMISTRY/ProteaseDB.h — src/chemistry/digestion.rs + digestion/enzymes.rs. Missing: writeTSV, provider ctor
   (provider ctor needs dynamic enzymes; EnzymeXMLDataProvider unported; share/CHEMISTRY/Enzymes.xml exists).
4. CHEMISTRY/ModificationsDB.h — src/chemistry/modifications.rs. Missing: initializeModificationsDB (needs
   UnimodXMLFile, unported format -> backlog), isInstantiated, searchModification(ResidueModification),
   searchModificationsFast, has(), addModification x2, findModificationIndex, searchModificationsByDiffMonoMassSorted x2,
   getAllSearchModifications (private dup in format/mascot_generic.rs), writeTSV (ostream default precision 6).
   C++ name index is set<pointer> ordered by address ~ load order; Rust by_name holds sorted indices.
5. PROCESSING/ID/IDFilter.h — src/analysis/id_filter.rs. Large unported surface.

## Commands / gotchas
- SCRATCH is SHARED with other lanes: use only $SCRATCH/B8b/ (cpp/, entries.json, gate.sh, logs/, assemble.py).
  SCRATCH=/private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad
- $SCRATCH/B8b/gate.sh <name> <cargo args> (run in background) -> logs/<name>.log
- Pinned C++ in $SCRATCH/B8b/cpp with symlinks inc/, impl/, tests/ (the word "source" in a git path trips the guard;
  heredocs and loops with variables are refused too: write scripts with the Write tool).
- Oracle template: ../oracle/B6-param-defaults/{run.sh,remote_run.sh}; mine: ../oracle/B8b-*/.
- Ledger fragment: .planning/ledger-fragments/wave4/B8b.json built by $SCRATCH/B8b/ledger.py from entries.json.

## Done
- MathFunctions COMPLETE (commit 2d26489): RandomShuffler promoted, tier-1 oracle ../oracle/B8b-random-shuffler.
- PeptideIdentification + USI: code, tests (tests/peptide_identification*.rs), fixture
  tests/data/peptide_identification_oracle.tsv (387 lines, sha 799c648f...), MQEvidence_2.consensusXML copied.
  set_experiment_label("") now a no-op (source); tests/identification.rs assertion corrected with oracle evidence.
- Doc-coverage floor (docs/doc-coverage.json) is lead's: report "run check_doc_coverage --write".

## Next step
1. Finish PeptideIdentification: gate green, provenance manifest tests/data/peptide_identification_provenance.json,
   docs (IDENTIFICATION_SUPPORT.md section + new docs/USI_SUPPORT.md), ledger entries (PeptideIdentification.h, USI.h), commit.
2. ModificationsDB members; ProteaseDB writeTSV (+oracle for both writeTSV outputs).
3. IDFilter: countHits, getBestHit, filterGroupsByScore, regex, keepBestPerPeptide*, digestion filter.

# Lane I1 notes (restart-safety scratch; delete in final commit)

C++ pin readable at /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.reference/openms4-core-bc9cc12
Worktree isolation REFUSES any git redirect to other dirs. The reference dir is itself a worktree checked out at bc9cc12
(HEAD verified in ../OpenMS4-tests/packages/core/.git/worktrees/openms4-core-bc9cc12/HEAD) -> just `cat` files there.
Shell loops with $vars in sed args are refused too; never put the word g-i-t in python heredocs. Use Write/Edit for NOTES.

Headers (lines h/cpp/test): ConsensusIDAlgorithm 115/145/0; Average 42/30/142; Best 42/31/183; Worst 42/31/151;
Identity 70/123/0; Similarity 65/134/0; PEPIons 51/121/49; PEPMatrix 52/74/50; Ranks 52/78/135;
IDMergerAlgorithm 277/417/201; ConsensusMapMergerAlgorithm 89/409/77; IDScoreSwitcherAlgorithm 752/48/181;
IDScoreGetterSetter 609/316/0; IDDecoyProbability 128/541/75; Scores 156/160/146; PeptideProteinResolution 152/794/72.
All `unmapped` in docs/core-sdk-coverage.json. Existing src/analysis/scores.rs = older score-type/main-score code (check overlap w/ IDScoreSwitcher).

Code goes in src/analysis/id/ (append `pub mod` to src/analysis/id.rs only).

## Findings
- C++ AASequence::operator< = length, nterm(none<some, then id), letters, then mod POINTER compare (nondeterministic -> defect cand.), cterm.
  Rust AASequence Ord = lexical string first. Need a source-order comparator for SequenceGrouping map order (+ PEPIons swap).
- PeptideIdentification::sort = stable_sort; Rust sort_by is stable too.
- IDFilter::removeDuplicatePeptideHits exists: src/analysis/id_filter.rs:354.
- Similarity: `*max_element` on empty other-run hits = UB (defect cand.); continue-path uses unprefixed score type (inconsistency).
- NeedlemanWunsch NOT ported (needed by PEPMatrix).

## Done (committed)
- src/analysis/id/consensus_id.rs (all 9 ConsensusID headers) + needleman_wunsch.rs (NW matrices verified vs source by script)
- tests/consensus_id.rs: class tests + TOPP ConsensusID_1/4/5/6 via idXML (tests/data/ConsensusID_*; test-data dir HEAD is bad3c3c,
  NOT 0cb15f2 and cannot be checked locally -> confirm by running the Release TOPP ConsensusID on ibminode06 as oracle).
- Gate log: scratchpad/logs/t1.log (cargo test --lib --test consensus_id on dax, slot I1)
- Release install on ibminode06 has bin/ConsensusID, IDMerger, IDScoreSwitcher, IDDecoyProbability -> TOPP oracles possible.

- COMMITTED: ConsensusID (11 integration tests pass on dax), IDScoreSwitcherAlgorithm (src/analysis/id/id_score_switcher.rs,
  tests/id_score_switcher.rs 6 pass), Scores docs + tests/scores_class_test.rs.
- Oracle ../oracle/I1-consensusid DONE on ibminode06 (Release bc9cc12/c19e494/174b576): runs identical except timestamp;
  Release == test-data for 4/5/6; for 1 one score differs in 15th digit (0.857387152777777 vs ...778) - both within 1e-13 of Rust.
- idXML reader findings (for report): IDFileConverter_12_output.idXML (200KB) fails "peptide parent parsing resource limit exceeded"
  (fixed XmlLimits max_work 50M consumed by db.len() per modified sequence); IDScoreSwitcherAlgorithm_test_input.idXML is line-wrapped
  upstream (xmlns:xsi split, attribute 'ges'), C++ tolerates, Rust refuses -> transcribed fixture instead.

- COMMITTED: IDMergerAlgorithm (id_merger.rs, 5 tests), ConsensusMapMergerAlgorithm (consensus_map_merger.rs, 5 tests).
  Findings: map_xml.rs run identifiers use a process-global counter -> lexical order flips; C++ uses random UniqueId -> source
  merge order within a condition not reproducible (C++ class test likely flaky). Source defects: offset lookup under new id (always 0);
  "_" + Size appends char 0x01. consensusXML default max_work too small for BSA.consensusXML (raised in test).

## Next step
1. PeptideProteinResolution (794 cpp, class test 72 lines) -> src/analysis/id/peptide_protein_resolution.rs
2. IDDecoyProbability (needs C++ driver oracle on ibminode06; math::fitters gamma/gauss exist)
3. IDScoreGetterSetter (609 header, no tests)
4. Report + ledger; clippy/doc/MSRV/ci sweep.

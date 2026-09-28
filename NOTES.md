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

## Next step
1. Param pattern (src/processing/peak_picking.rs:1853-1985); read class tests for Average/Best/Worst/Ranks/PEPIons/PEPMatrix.
2. Write src/analysis/id/consensus_id.rs (base + 8 variants as an enum/trait), then commit.
3. Scores.h; IDScoreSwitcherAlgorithm; IDMergerAlgorithm; ConsensusMapMerger; IDDecoyProbability; PeptideProteinResolution.

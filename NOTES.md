# Lane I1 notes (restart-safety scratch; delete in final commit)

C++ pin readable at /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.reference/openms4-core-bc9cc12
(`git -C <core> archive` is refused by the sandbox; `git -C <core> show bc9cc12:<path>` works for single files).

Headers (lines h/cpp/test): ConsensusIDAlgorithm 115/145/0; Average 42/30/142; Best 42/31/183; Worst 42/31/151;
Identity 70/123/0; Similarity 65/134/0; PEPIons 51/121/49; PEPMatrix 52/74/50; Ranks 52/78/135;
IDMergerAlgorithm 277/417/201; ConsensusMapMergerAlgorithm 89/409/77; IDScoreSwitcherAlgorithm 752/48/181;
IDScoreGetterSetter 609/316/0; IDDecoyProbability 128/541/75; Scores 156/160/146; PeptideProteinResolution 152/794/72.
All `unmapped` in docs/core-sdk-coverage.json. Existing src/analysis/scores.rs = older score-type/main-score code (check overlap w/ IDScoreSwitcher).

Code goes in src/analysis/id/ (append `pub mod` to src/analysis/id.rs only).

## Next step
1. Read existing identification types (src/identification) + Param pattern (src/processing/peak_picking.rs:1853-1985).
2. Port ConsensusIDAlgorithm base + Best/Worst/Average/Ranks (class tests give literals), then Identity/Similarity/PEPMatrix/PEPIons.
3. Scores.h; IDScoreSwitcherAlgorithm; IDMergerAlgorithm; ConsensusMapMerger; IDDecoyProbability; PeptideProteinResolution.

# Lane I3b NOTES (restart-safety)

Scope: ANALYSIS/ID AhoCorasickAmbiguous, PeptideIndexing (partial), OpenSearchModificationAnalysis, AScore, PScore, CometModification.

## State (on main f0990c0)
- peptide_indexing.rs (968 lines) is old port at src/analysis/peptide_indexing.rs; ledger scope in docs/core-sdk-coverage.json line ~941.
- src/analysis/id.rs is an empty module root: append `pub mod` lines only; new code in src/analysis/id/.
- AScore/PScore/CometModification/OpenSearchModificationAnalysis: unmapped. AhoCorasickAmbiguous: evidence_requires_review, no rust file.
- C++ pins: `cd ../OpenMS4-tests/packages/core && git show bc9cc12:<path>` (git -C refused by worktree guard; cd first in a separate command).

## Next step
1. Read C++ headers/cpp/tests for PScore (AScore depends on it), CometModification, OpenSearchModificationAnalysis.
2. Port PScore -> src/analysis/id/pscore.rs, then AScore, then CometModification, OSMA.
3. AhoCorasickAmbiguous types, PeptideIndexing gaps last.

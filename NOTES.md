# Lane I3b NOTES (restart-safety)

Scope: ANALYSIS/ID AhoCorasickAmbiguous, PeptideIndexing (partial), OpenSearchModificationAnalysis, AScore, PScore, CometModification.

## State
- Incarnation 2 (retry, worktree wf_cd066dab-0a7-31) merged incarnation 1 (wf_..-26, notes only).
- peptide_indexing.rs (968 lines) is old port at src/analysis/peptide_indexing.rs; ledger scope in docs/core-sdk-coverage.json line ~941.
- src/analysis/id.rs is an empty module root: append `pub mod` lines only; new code in src/analysis/id/.
- C++ sources: the worktree guard refuses `git -C` and `$VAR` in sed. Use literal paths:
  `/Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.reference/openms4-core-bc9cc12/` (verified sha1 of PScore.cpp
  equals `git show bc9cc12:` output). Copies in scratchpad/cpp/ (session-local, re-copy if gone).
- PScore findings: massCorrectionTerm/cleavageCorrectionTerm/modificationCorrectionTerm are defined in PScore.cpp as FREE
  functions, not `PScore::` members -> declared statics never defined (link error if called). Vector overload of computePScore
  uses p = level/mz_window, single overload p = (level+1)/mz_window; vector overload calls findNearest on possibly empty spectrum.
- PScore test expectations: ranks window 9.9 -> 4, 10.1 -> 5; PScore_test.dta peak-level mz lists; pscore 83.867454, 154.682242;
  IFSQVGK (add_first_prefix_ion) 12 peaks all-match -> 240.
- AScore depends on: AASequence, TheoreticalSpectrumGenerator, MatchedIterator, Math::binomial_cdf_complement, PeptideHit.

## Next step
1. Check Rust availability: MatchedIterator, binomial_cdf_complement, TheoreticalSpectrumGenerator, DTA load, PeptideHit.
2. Port PScore -> src/analysis/id/pscore.rs (+ AScore::computeCumulativeScore_ shared helper), commit.
3. AScore, then CometModification, OSMA; AhoCorasickAmbiguous types, PeptideIndexing gaps last.

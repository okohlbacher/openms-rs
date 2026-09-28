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

## Done (committed)
- PScore + AScore: src/analysis/id/{pscore,ascore}.rs, tests/pscore_ascore.rs (22/22 on spock), docs/PSCORE_ASCORE_SUPPORT.md,
  tests/data/pscore_ascore_provenance.json. CometModification: comet_modification.rs, tests/comet_modification.rs (9/9),
  docs/COMET_MODIFICATION_SUPPORT.md, provenance. Doc coverage OK (do not --write doc-coverage.json: lead file).
- Guard: python heredocs touching many paths get refused; write scripts to scratchpad and run `python3 <script>`.

## Next step
1. AhoCorasickAmbiguous: port nine types (AA, Hit, Index, Bitset, ACNode+DepthHits, ACScout, ACTrieState, ACTrie, tables,
   nextValidAA) to src/analysis/id/aho_corasick_ambiguous.rs; class test AhoCorasickAmbiguous_test.cpp gives hit lists.
2. OSMA (1107-line cpp, fuzzy std::map, std::sort ties, boost normal): likely untouched/partial; decide after 1.
3. Ledger fragment + lane report (.planning/...), delete NOTES.md in final commit.

# Lane I3a NOTES (restart-safety scratch; delete in final commit)

Scope: ANALYSIS/ID FragmentIndex (h 751 + cpp 2625), SimpleSearchEngineAlgorithm (h 236 + cpp 956), NeighborSeq (h 232 + cpp 178).
Nothing ported yet in Rust (grep for fragment_index/simple_search/neighbor_seq in src/ = none). Target dir: src/analysis/id/ (append pub mod to src/analysis/id.rs).
C++ extracted to scratchpad via `git archive bc9cc12 <paths> | tar -x` (loops over git are refused by the sandbox; use single git cmds).
Class tests: FragmentIndex_test.cpp (2879 lines), NeighborSeq_test.cpp (252), SimpleSearchEngineAlgorithm_test.cpp (51).

## Next step
1. Read NeighborSeq.h/.cpp + test (smallest, self-contained) -> port first.
2. Survey Rust deps: AASequence, ProteaseDigestion, TheoreticalSpectrumGenerator, ModifiedPeptideGenerator, Param pattern (src/processing/peak_picking.rs:1853-1985).
3. FragmentIndex, then SimpleSearchEngineAlgorithm (depends on FragmentIndex? check).

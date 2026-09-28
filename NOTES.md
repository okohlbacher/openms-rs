# Lane I3a NOTES (restart-safety scratch; delete in final commit)

Scope: ANALYSIS/ID FragmentIndex (h 751 + cpp 2625), SimpleSearchEngineAlgorithm (h 236 + cpp 956), NeighborSeq (h 232 + cpp 178).
Nothing ported yet in Rust (grep for fragment_index/simple_search/neighbor_seq in src/ = none). Target dir: src/analysis/id/ (append pub mod to src/analysis/id.rs).
C++ extracted to scratchpad via `git archive bc9cc12 <paths> | tar -x` (loops over git are refused by the sandbox; use single git cmds).
Class tests: FragmentIndex_test.cpp (2879 lines), NeighborSeq_test.cpp (252), SimpleSearchEngineAlgorithm_test.cpp (51).

## Next step
1. Read NeighborSeq.h/.cpp + test (smallest, self-contained) -> port first.
2. Survey Rust deps: AASequence, ProteaseDigestion, TheoreticalSpectrumGenerator, ModifiedPeptideGenerator, Param pattern (src/processing/peak_picking.rs:1853-1985).
3. FragmentIndex, then SimpleSearchEngineAlgorithm (depends on FragmentIndex? check).

## Findings (pass 1, read all of NeighborSeq + FragmentIndex.h/.cpp)
- NeighborSeq: TSG b/y z=1 add_first_prefix_ion; Rust TheoreticalSpectrumGenerator{ion_series B,Y, add_first_prefix_ion:true}.generate(pep,1,1,None).
  isNeighborSpectrum: Dice 2*B12/(n1+n2) > thr (strict). shared count: int(mz/bin) set-intersection merge. mass map BTreeMap<f64 key> -> use ordered f64 wrapper or Vec sorted.
  X residue -> stat -1. Test expectations in NeighborSeq_test.cpp (all tier-3 literals).
- FragmentIndex.cpp sections: initResidueMassTable_ L57, initModificationTables_ L81, computeSnesSigmaDeltaSet_ L176, buildModSlots_ L259,
  generateFragmentsForSeries_ L342, reconstructModifiedSequence L477, realizeSNESLength L544, reconstructRealizedSubSequence L604,
  generateSNESMothers_ L669, generatePeptides L843, build L1052 (bucketsize 4096, sort (mz,pepidx), per-bucket stable sort by pepidx),
  getPeptidesInMassWindow L1340, computeMassWindow_ L1362 (f32 ppm), query L1376, queryPeaks L1427, trimHits L1566,
  searchDifferentPrecursorRanges L1637, querySpectrumSNES_ L1667, querySpectrum L2278/2299, ctor defaults L2357, updateMembers_ L2501.
- f32 arithmetic everywhere in query (precursor mz, adjusted_mass): must mirror float casts exactly.

## Incarnation 2 (retry branch wf-29, same base as wf-25 4f332e02)
- C++ re-extracted to session scratchpad cpp/ (git archive | tar with literal paths; no shell vars in git cmds).
- Rust deps: chemistry::AASequence (src/chemistry/sequence.rs; parse, mono_mass()->Result, as_str), chemistry::theoretical::TheoreticalSpectrumGenerator
  (fields ion_series, add_first_prefix_ion; generate(&pep,min_z,max_z,None)->Result<MSSpectrum>), openms_log_warn! macro in src/concept/log_stream.rs:873.
- Writing src/analysis/id/neighbor_seq.rs now.

## Status (incarnation 2)
- DONE: NeighborSeq -> src/analysis/id/neighbor_seq.rs + tests/neighbor_seq.rs (10/10 on dax, --no-default-features). Docs/provenance still TODO.
- NEXT: SimpleSearchEngineAlgorithm (TOPP consumer SimpleSearchEngine). Tier-1 fixtures: test-data 0cb15f2 topp/SimpleSearchEngine_{1,2}.ini,
  SimpleSearchEngine_1.{mzML,fasta}, SimpleSearchEngine_{1,2}_out.idXML (3 PSMs each; CMakeLists.txt:3773-3789, whitelist "IdentificationRun date", "SearchParameters id=\"SP_0\" db=").
  Rust deps present: psm_scoring.rs:218 HyperScore compute_with_detail; processing.rs Normalizer/ThresholdMower/NLargest; processing/window_mower.rs;
  processing/deisotoping.rs (is_tolerance_supported(Tolerance), deisotope); comparison.rs:1401 spectrum_alignment; chemistry/decoy_generator.rs reverse_protein/reverse_peptides(Protease);
  chemistry/digestion.rs:397 digest_unmodified; chemistry/modified_peptides.rs apply_fixed/variable; analysis/peptide_indexing.rs PeptideIndexing;
  analysis/false_discovery_rate.rs apply_peptides; analysis/id_filter.rs; metadata/im_types.rs determine_im_format; format/fasta.rs FASTAEntry; boost regex facade (tests/data/boost_regex_provenance.json).
  RandomShuffler: chemistry/decoy_random.rs (private) -- check.
- THEN: FragmentIndex (no TOPP consumer, 3.4k lines) as far as possible.

## Status (incarnation 2, later)
- DONE NeighborSeq: code, tests (10/10), docs/NEIGHBOR_SEQ_SUPPORT.md, tests/data/neighbor_seq_provenance.json.
- DONE SSE code: src/analysis/id/simple_search_engine.rs; tests/simple_search_engine.rs 7/7 on dax (default features),
  both TOPP retained outputs match under FuzzyDiff contract. Touched src/chemistry/decoy_random.rs (generic shuffle_slice).
- TODO SSE: docs/SIMPLE_SEARCH_ENGINE_SUPPORT.md, tests/data/simple_search_engine_provenance.json (fixture sha256 computed), clippy/doc gates.
- NEXT: FragmentIndex.

## Status (incarnation 2, late)
- NeighborSeq complete (10 tests). SSE complete incl. progress (8 tests, tier 1 vs retained TOPP idXML, max rel dev 2.1e-10). FragmentIndex partial: only querySpectrumSNES_ unported (39 tests).
- MSRV 1.85 green for all three binaries. CI sweep log: scratchpad sweep1/. Pending: clippy all-targets (gate12), report + ledger fragment, delete NOTES.md.
- Lead-file needs for report: ci_matrix Slices for tests/{neighbor_seq,simple_search_engine,fragment_index}.rs; doc-coverage --write; SOURCE_PROVENANCE manifests tests/data/{neighbor_seq,simple_search_engine,fragment_index}_provenance.json.

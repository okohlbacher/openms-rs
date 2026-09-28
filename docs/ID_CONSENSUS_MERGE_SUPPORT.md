# Consensus identification, run merging and score handling

Ports of eleven `ANALYSIS/ID` headers and `ANALYSIS/SEQUENCE/NeedlemanWunsch.h`
at core `bc9cc12`. Provenance: `tests/data/id_consensus_merge_provenance.json`.

| Header | Rust | Status |
|---|---|---|
| `ConsensusIDAlgorithm.h` (+ `Identity`, `Similarity`, `Average`, `Best`, `Worst`, `Ranks`, `PEPMatrix`, `PEPIons`) | `src/analysis/id/consensus_id.rs` | complete |
| `ANALYSIS/SEQUENCE/NeedlemanWunsch.h` | `src/analysis/id/needleman_wunsch.rs` | complete |
| `Scores.h` | `src/analysis/scores.rs` (`ScoreType`, free functions) | complete |
| `IDScoreSwitcherAlgorithm.h` | `src/analysis/id/id_score_switcher.rs` | complete |
| `IDMergerAlgorithm.h` | `src/analysis/id/id_merger.rs` | complete |
| `ConsensusMapMergerAlgorithm.h` | `src/analysis/id/consensus_map_merger.rs` | complete |
| `PeptideProteinResolution.h` | `src/analysis/id/peptide_protein_resolution.rs` | complete |
| `IDDecoyProbability.h` | `src/analysis/id/id_decoy_probability.rs` | partial: `generateDistributionImage_` not ported |
| `IDScoreGetterSetter.h` | not ported | untouched |

## API mapping

### ConsensusIDAlgorithm family

| C++ | Rust |
|---|---|
| `ConsensusIDAlgorithm::apply(ids, se_info, number_of_runs)` | `ConsensusIDAlgorithm::apply_with_search_engines` (trait) |
| `ConsensusIDAlgorithm::apply(ids, number_of_runs)` | `ConsensusIDAlgorithm::apply` |
| members `considered_hits_`, `min_support_`, `count_empty_`, `keep_old_scores_` | `ConsensusIDFilter` fields |
| `number_of_runs_` | argument of `apply*` (0 means `ids.len()`) |
| `HitInfo`, `SequenceGrouping` (protected) | private `HitInfo`, `Grouping` ordered by source `AASequence::operator<` |
| `apply_`, `preprocess_`, `getAggregateScore_`, `getSimilarity_` (protected virtual) | private functions `identity_group`, `identity_preprocess`, `ranks_preprocess`, `similarity_group`, `pep_matrix_similarity`, `pep_ions_similarity` |
| `compareChargeStates_` | private `compare_charge_states` |
| `updateMembers_` / `DefaultParamHandler` | `defaults`, `from_param`, `from_param_with_warnings`, `to_param` on each concrete type |
| `ConsensusIDAlgorithmSimilarity::similarities_` | a per-call cache inside `pep_ions_similarity` (PEPMatrix never used it) |
| `ConsensusIDAlgorithmRanks::current_number_of_runs_`, `current_considered_hits_` | locals of the Ranks preprocessing |
| `ConsensusIDAlgorithmPEPMatrix::alignment_`, `matrix`, `penalty` | `ConsensusIDAlgorithmPEPMatrix::{matrix, penalty}` |
| `ConsensusIDAlgorithmPEPIons::mass_tolerance_`, `min_shared_` | `ConsensusIDAlgorithmPEPIons::{mass_tolerance, min_shared}` |
| deleted copy constructors | `Clone` (the types hold only parameters) |

### NeedlemanWunsch

| C++ | Rust |
|---|---|
| `ScoringMatrix::{identity, PAM30MS}` | `ScoringMatrix::{Identity, PAM30MS}` |
| `NamesOfScoringMatrices` | `ScoringMatrix::NAMES` |
| constructors | `NeedlemanWunsch::new`, `Default` (PAM30MS, 5) |
| `align` | `align` (returns `i64`) |
| `setMatrix(ScoringMatrix)`, `getMatrix` | public field `matrix` |
| `setMatrix(string)` | `set_matrix_by_name`, `ScoringMatrix::from_name` |
| `setPenalty`, `getPenalty` | public field `penalty` |
| `first_row_`, `second_row_` | local rolling rows |

### Scores

| C++ | Rust |
|---|---|
| `IDType` | `ScoreType` |
| `isScoreType` | `ScoreType::matches` |
| `parseIDType` | `ScoreType::parse` |
| `isHigherBetter` | `ScoreType::higher_is_better` |
| `getAllIDScoreNames` | `all_score_names` |
| `getIDNamesForType` | `ScoreType::names` |
| `findIDTypeByName` | `ScoreType::from_name` |
| `normalizeScoreName` | `normalize_score_name` |
| `isKnownScoreType` | `is_known_score_type` |
| `Maps_`, `getMaps_` | the static tables inside `ScoreType::names` |

### IDScoreSwitcherAlgorithm

| C++ | Rust |
|---|---|
| constructor, `updateMembers_`, parameters | `IDScoreSwitcherAlgorithm` fields, `defaults`, `from_param*`, `to_param`, `higher_better`, `effective_new_score_type` |
| `ScoreType` typedef | `crate::analysis::scores::ScoreType` |
| `isScoreType`, `toScoreTypeEnum`, `isScoreTypeHigherBetter`, `getScoreNames` | associated functions of the same names |
| `ScoreSearchResult`, `findScoreType<T>` | `ScoreSearchResult`, `find_score_type<I: SwitchableIdentification>` |
| `switchScores<T>(id, counter)` | `switch_scores` (returns the count) |
| `switchToGeneralScoreType` (vector / list / ConsensusMap) | `switch_to_general_score_type`, `switch_consensus_map_to_general_score_type` |
| `determineScoreNameOrientationAndType` (list / ConsensusMap) | `determine_score_name_orientation_and_type`, `determine_consensus_map_score_name_orientation_and_type` |
| `switchScores(ConsensusMap&, …)`, `switchScores(PeptideIdentificationList&, …)` | `switch_consensus_map_scores`, `switch_peptide_scores` |
| `IDSwitchResult` | `IDSwitchResult` |
| `switchToScoreType` (ConsensusMap / list) | `switch_consensus_map_to_score_type`, `switch_to_score_type` |
| `switchBackScoreType` (ConsensusMap / list) | `switch_back_consensus_map_score_type`, `switch_back_score_type` |
| `tolerance_` | `SCORE_TOLERANCE` |

### IDMergerAlgorithm and ConsensusMapMergerAlgorithm

| C++ | Rust |
|---|---|
| `IDMergerAlgorithm(runIdentifier, addTimeStampToID)` | `IDMergerAlgorithm::new`, `Default` |
| parameters `annotate_origin`, `allow_disagreeing_settings` | `IDMergerParameters` |
| `insertRuns` (move and copy) | `insert_runs` (by value; clone to keep the input) |
| `returnResultsAndClear` | `return_results_and_clear` (returns the pair) |
| `getNewIdentifier_`, `copySearchParams_`, `checkOldRunConsistency_` (both), `insertProteinIDs_`, `updateAndMovePepIDs_`, `movePepIDsAndRefProteinsToResultFaster_` | private helpers |
| `ConsensusMapMergerAlgorithm::mergeProteinsAcrossFractionsAndReplicates`, `mergeAllIDRuns`, `mergeProteinIDRuns` | same names in snake case |
| parameter `annotate_origin` (declared, never read by the source) | field `annotate_origin` |
| inherited `ProgressLogger` (never started by either source class) | not ported |

### PeptideProteinResolution

| C++ | Rust |
|---|---|
| `ConnectedComponent`, `operator<<` | `ConnectedComponent`, `Display` |
| constructor `(statistics)` | `new` |
| `resolve`, `run` (static) | `resolve`, `run` |
| `buildGraph`, `resolveGraph`, `findConnectedComponent`, `resolveConnectedComponent` | same names in snake case |

### IDDecoyProbability

| C++ | Rust |
|---|---|
| constructor, copy, assignment, parameters | `IDDecoyProbability` fields, `Default`, `Clone`, `defaults`, `from_param*`, `to_param` |
| `apply(ids)` | `apply` |
| `apply(prob_ids, fwd_ids, rev_ids)` | `apply_forward_reverse` (returns the probabilities) |
| `Transformation_`, `normalizeBins_`, `getProbability_`, `apply_` | private |
| `generateDistributionImage_` (two overloads; gnuplot files, called only in `IDDECOYPROBABILITY_DEBUG` builds) | **not ported** |
| debug parameters `rev_filename`, `fwd_filename` | not ported (absent from release builds) |

## Preserved source conventions

- Consensus output keeps only score type and orientation; hits carry
  `consensus_support`, `target_decoy` and, with `keep_old_scores`,
  `<type>_score`; ties in the final score keep the source `SequenceGrouping`
  order (length first, as `AASequence::operator<`).
- PEPIons' ion ladder adds full (water-including) residue masses at every
  step, its similarity is a single-precision ratio, and the sequence pair is put
  in source order before the (asymmetric) comparison.
- NeedlemanWunsch's `J`/`O`/`U` entries are `INT16_MAX`, which the maximum
  recurrence rewards although the source documentation says they forbid.
- `IDScoreSwitcherAlgorithm::switch_scores` overwrites a `~` backup, and treats
  a `0/0` relative difference as equal, exactly as the source.
- `PeptideProteinResolution::resolve` keeps the source's lowest-probability-wins
  choice within target and within decoy groups.
- `IDDecoyProbability` starts its score range from `f64::MIN_POSITIVE`, uses
  `number_of_bins - 1` for the forward and `number_of_bins` for the reverse
  histogram, and writes both `<type>_Score` and `<type>_score`.

## Native differences

- Every operation is transactional: on error the input is unchanged (the
  source leaves partial results).
- Warnings the source logs are returned as strings.
- ConsensusID: a residue modification pair is ordered by name and full ID
  (the source compares pointers); an identification without hits contributes
  no best match in the similarity algorithms (the source dereferences an empty
  `max_element`); Ranks refuses more supporting identifications than runs (the
  source's unsigned subtraction wraps); a zero self-alignment score is refused
  (the source divides by zero); ceilings `MAX_CONSENSUS_HITS` and
  `MAX_SIMILARITY_EVALUATIONS`.
- NeedlemanWunsch refuses non-`A`-`Z` input (source out-of-bounds reads) and
  bounds lengths (`MAX_ALIGNMENT_LENGTH`, `MAX_ALIGNMENT_CELLS`).
- IDMergerAlgorithm and ConsensusMapMergerAlgorithm return protein hits in
  first-insertion order (the source's `unordered_set` hash order is
  implementation-defined).
- ConsensusMapMergerAlgorithm visits runs in first-matched order instead of
  identifier order (the identifiers are random, so the source order is not
  reproducible), shifts `id_merge_index` by the intended offset unless
  `source_merge_index_offsets` asks for the source's always-zero offset,
  and refuses out-of-range new run indices and negative merge indices.
- PeptideProteinResolution keeps component groups equal in probability and
  peptide count in index order (`std::sort` leaves it unspecified).
- IDDecoyProbability puts a reverse score at the overall maximum into the last
  bin (source out-of-bounds write), falls back from any Gauss fit error (source:
  only `UnableToFit`), and refuses empty or constant score sets and non-finite
  probabilities (the source computes NaN).

## Checked boundaries and evidence

| Header | Evidence (tier) |
|---|---|
| ConsensusID Average/Best/Worst/Ranks | every class-test literal (3); TOPP_ConsensusID_6 (1, retained output) |
| ConsensusID PEPMatrix, PEPIons | TOPP_ConsensusID_1, _4, _5 retained outputs to 1e-13 (1); the Release TOPP tool on ibminode06 reproduces them (only a 15th-digit difference in one score of _1) |
| NeedlemanWunsch | every class-test literal (3); matrices checked cell by cell against the source by script |
| Scores | every class-test assertion (3) |
| IDScoreSwitcherAlgorithm | every class-test assertion (3, fixture transcribed); TOPP_IDScoreSwitcher_2 retained output, exact (1); TOPP_IDScoreSwitcher_1 blocked (the native idXML reader exceeds its fixed work ceiling on `IDFileConverter_12_output.idXML`) |
| IDMergerAlgorithm | every class-test section (3) |
| ConsensusMapMergerAlgorithm | every class-test section (3); offsets derived (4) |
| PeptideProteinResolution | both class-test input/output file pairs (3); `resolve` derived (4) |
| IDDecoyProbability | executed Release oracle, three cases, max difference 7.4e-15 / 4.4e-16 / 9.1e-15 (1) |

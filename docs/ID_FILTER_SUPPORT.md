# Identification filtering

`openms::analysis::id_filter` ports the pinned `OpenMS::IDFilter` algorithms,
all but the members listed under [Coverage limits](#coverage-limits), to the existing native identification and feature
records. It uses free functions and standard-library collections. No additional
dependency, C++ build, search engine or identification graph is involved.

The source reference is OpenMS4-core revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4` (the header, implementation and
class test are unchanged since `7c029e8`, the first increment's pin):

- `src/openms/include/OpenMS/PROCESSING/ID/IDFilter.h`
- `src/openms/source/PROCESSING/ID/IDFilter.cpp`
- `src/tests/class_tests/openms/source/IDFilter_test.cpp`
- `src/tests/class_tests/openms/data/IDFilter_test.idXML`

## Implemented operations

Functions take mutable slices when filtering hits and mutable vectors when
removing whole identifications. Pass `std::slice::from_mut(&mut id)` to filter a
single record. Hit filters retain identification-level settings and metadata,
including when all hits are removed. They do not update the stored `hit.rank`.

| OpenMS operation | Native API and behavior |
| --- | --- |
| `filterHitsByScore` | `filter_peptides_by_score`, `filter_proteins_by_score`: inclusive cutoff, using each identification's score direction. No sorting. |
| `keepNBestHits` | `keep_n_best_peptide_hits`, `keep_n_best_protein_hits`: stable score sort, then at most N hits. N=0 clears hits. Ties can be split at the cutoff. |
| `filterHitsByRank` | `filter_peptides_by_rank`, `filter_proteins_by_rank`: inclusive **dense** ranks computed from sorted scores, e.g. `1,1,2,3`. Existing rank fields are unchanged, as in the pinned implementation. |
| `keepBestPeptideHits` | `keep_best_peptide_hits`: keep all best-score ties; `strict=true` removes every hit when more than one winner exists. |
| `keepNBestSpectra` | `keep_n_best_spectra`: sort each spectrum's hits, then retain N identifications by their best hit. Keep all hits in selected identifications. Empty identifications sort last. |
| `removeDecoyHits` | `remove_decoy_peptide_hits`, `remove_decoy_protein_hits`: remove exact unitless string annotations `target_decoy="decoy"` **or** `isDecoy="true"`. The legacy annotation also removes a hit explicitly marked target. Missing annotations and `target+decoy` alone are retained. |
| `removeEmptyIdentifications` | `remove_empty_peptide_identifications`, `remove_empty_protein_identifications`: discard records whose hit vector is empty. |
| `filterPeptidesByLength` | `filter_peptides_by_length`: inclusive residue-count bounds. Modification labels do not add residues. |
| `filterPeptidesByCharge` | `filter_peptides_by_charge`: inclusive signed charge bounds. Unknown zero and negative charges can be selected. |
| `keepHitsMatchingProteins`, `removeHitsMatchingProteins` | `filter_peptides_by_accessions`, `filter_proteins_by_accessions` with `MatchAction::Keep` or `Remove`. A peptide matches if **any** nonempty evidence accession is in the set; all its evidence remains attached. |
| `keepPeptidesWithMatchingModifications`, `removePeptidesWithMatchingModifications` | `filter_peptides_by_modifications` with a `MatchAction`. Match source **full IDs**, including terminal modifications, e.g. `Oxidation (M)` and `Acetyl (N-term)`. An empty set means any modification. |
| `extractPeptideSequences` | `extract_peptide_sequences`: return a set of canonical modified strings or unmodified sequences. |
| `keepPeptidesWithMatchingSequences`, `removePeptidesWithMatchingSequences` | `filter_peptides_by_sequences` with the set of canonical strings and a `MatchAction`; optionally ignore modifications. Charge is ignored. |
| `keepUniquePeptidesPerProtein` | `keep_unique_peptides_per_protein`: keep the exact unitless annotation `protein_references="unique"`. Does not infer uniqueness from evidence count. |
| `removeDuplicatePeptideHits` | `remove_duplicate_peptide_hits`: preserve the first full-record occurrence with `DuplicatePolicy::Exact`, or the first complete modified-sequence value with `Sequence`. Same-text records with different custom chemistry remain distinct. No sorting and no preference for a better score. |
| `filterPeptidesByRT`, `filterPeptidesByMZ` | `filter_peptides_by_rt`, `filter_peptides_by_mz`: inclusive coordinate range over whole identifications. Missing coordinates do not match. |
| `filterPeptidesByMZError` | `filter_peptides_by_mz_error`: inclusive precursor tolerance; `Tolerance::Ppm` is relative to observed precursor m/z. Charge zero assumes +1. Uses modification-aware `AASequence::mz`. |
| `removeUnreferencedProteins` | `remove_unreferenced_proteins`: retain proteins referenced by peptide evidence **within the same run identifier**. |
| `removeDanglingProteinReferences` | `remove_dangling_protein_references`: remove evidence with no matching surviving protein in the same run. Optionally remove hits with no remaining evidence. Repeated protein identification records with one run ID contribute their union of accessions. |
| `updateProteinGroups` | `update_protein_groups`: remove absent members and empty groups, retaining group probabilities and all sample arrays. Returns false only when a surviving group lost some members. Removing an entire group alone returns true. |
| `removeUngroupedProteins` | `remove_ungrouped_proteins`: retain hits appearing in any supplied group. |

### Feature and consensus maps

The map entry points visit **top-level features**, plus unassigned peptide
identifications. They keep feature coordinates, quantitative data, hulls,
subordinate features, and other metadata attached. Subordinate peptide IDs are
not visited. FeatureMap adapters extend the source's generic map conveniences;
ConsensusMap cleanup follows the source's run-aware algorithms.

| Purpose | FeatureMap | ConsensusMap |
| --- | --- | --- |
| Top N per peptide ID | `keep_n_best_hits_in_feature_map` | `keep_n_best_hits_in_consensus_map` |
| Remove unreferenced proteins | `remove_unreferenced_feature_proteins` | `remove_unreferenced_consensus_proteins` |
| Remove dangling evidence | `remove_dangling_feature_references` | `remove_dangling_consensus_references` |
| Remove empty peptide IDs | `remove_empty_feature_identifications` | `remove_empty_consensus_identifications` |

The unreferenced-protein functions take `include_unassigned`. Dangling-evidence
and empty-ID cleanup always visit unassigned IDs as well as top-level features.
Empty-ID cleanup does not remove features or protein identification records.

## Native validation and deliberate differences

- All functions returning `Result` validate affected identification records
  before mutation. An invalid later record cannot leave earlier records
  partially filtered. Map functions validate the complete map first.
- Numeric cutoffs must be finite. Length, charge and rank filters use
  `Option` for an absent upper bound, avoiding sentinel arithmetic. Reversed
  bounds are errors; the C++ length/charge/rank filters silently ignore some
  reversed upper bounds. Rank zero is an error in Rust.
- Sorting preserves input order for equal scores, including signed-zero ties.
  `keep_n_best_spectra` requires one score type **and direction** across all
  records, including empty ones. The pinned C++ implementation inconsistently
  checks empty score types, does not check directions, and has a comparator
  that returns true for two empty spectra. Rust gives empty spectra a valid,
  deterministic order.
- Exact duplicate removal uses the native hit's equality, preserving any
  differences in metadata floats, units and analysis results. The pinned C++
  `DataValue` compares scalar floats within `1e-6`, and `PeptideHit::operator==`
  omits analysis results. Rust does not discard those distinct native values.
  Sequence-only removal keeps modifications distinct and ignores other fields.
- Full-record duplicate comparison follows the source's quadratic scan with a
  preflight worst-case budget of `MAX_EXACT_DUPLICATE_COMPARISONS` (10 million
  comparisons per call). Exceeding it returns an error before mutation.
  Sequence-only comparison uses an ordered set and has no quadratic scan.
- Precursor-error filtering requires a positive observed m/z on every supplied
  identification. Missing coordinates, invalid calculated ion m/z and overflow
  are errors. All decisions are computed before mutation. Negative charges
  follow the native `AASequence::mz` ion convention; no absolute-charge rewrite
  is performed by the filter.
- Metadata string matching is exact, including the absence of units, as in
  source `DataValue` equality. It does not reinterpret integers as booleans or
  infer decoys from accession prefixes.

Protein filtering does not automatically update protein groups, peptide
evidence, coverage, or previously calculated protein modifications. A typical
explicit sequence is: filter hits; remove dangling evidence; remove empty
peptide IDs; remove unreferenced proteins; update each protein group vector;
recompute coverage if needed. Supply only accepted peptide hits to coverage
calculation; candidate alternatives also contribute if they remain present.

## Coverage limits

This is not yet complete `IDFilter` class parity; the member table below is the
authority. Three source members are not ported, each with a named dependency:

- `filterHitsByScore(ids, threshold, IDScoreSwitcherAlgorithm::ScoreType)`:
  needs `ANALYSIS/ID/IDScoreSwitcherAlgorithm.h` (`isScoreType`,
  `findScoreType`, `isScoreTypeHigherBetter`), which is unported.
- The three `AnnotatedMSRun` overloads (`filterHitsByScore`, `keepNBestHits`,
  `keepHitsMatchingProteins`): need `KERNEL/AnnotatedMSRun.h`, which is
  unported.
- `annotateBestPerPeptideWithData`, `annotateBestPerPeptidePerRunWithData` and
  the typedefs `ChargeToPepHitP`, `SequenceToChargeToPepHitP`,
  `RunToSequenceToChargeToPepHitP`: their state is a map of raw pointers into
  caller-owned hits, carried between calls. A borrowed-pointer index has no safe
  Rust counterpart; the index is internal, and the map-level functions (the
  reason the source exposes them) share one index across all of a map's lists.

## Member mapping (`PROCESSING/ID/IDFilter.h`, core `bc9cc12`)

The header is unchanged since the `7c029e8` pin the first increment used.
Items in `openms::analysis::id_filter`. "list" means the source's
`PeptideIdentificationList` overload, which the Rust slice functions cover.

| Source member | Rust |
| --- | --- |
| `IDFilter()`, `~IDFilter()` | none: free functions |
| typedefs `ChargeToPepHitP`, `SequenceToChargeToPepHitP`, `RunToSequenceToChargeToPepHitP` | not ported (see above) |
| `HasGoodScore` | `HasGoodScore` |
| `HasMetaValue` | `HasMetaValue` (source `DataValue` equality: same type and unit, doubles within 1e-6) |
| `HasMaxMetaValue`, `HasMinMetaValue` | `HasMaxMetaValue`, `HasMinMetaValue` (`Result`: a non-numeric value is an error) |
| `HasDecoyAnnotation` | `HasDecoyAnnotation` |
| `HasMatchingAccessionImpl`, `HasMatchingAccession`, `HasMatchingAccessionUnordered` | `HasMatchingAccession` over `AccessionSet` (`BTreeSet`/`HashSet`) |
| `GetMatchingItems<PeptideEvidence, FASTAEntry>` (ctor, `getKey`, `exists`, `getHitKey`, `getValue`) | `GetMatchingItems` (`new`, `key`, `exists`, `hit_key`, `value`) |
| `HasMinPeptideLength`, `HasMinCharge`, `HasLowMZError`, `HasMatchingModification`, `HasMatchingSequence`, `HasNoEvidence` | same names |
| `PeptideDigestionFilter` (ctor, `disabledValue`, `operator()`, `filterPeptideSequences`) | `PeptideDigestionFilter` (`new` with `Option` bounds for `-1`, `matches`, `filter_peptide_sequences`) |
| `DigestionFilter` (ctor, `operator()`, `filterPeptideEvidences`) | `DigestionFilter` (`new`, `matches`, `filter_peptide_evidences`) |
| `HasNoHits` | `HasNoHits` |
| `HasRTInRange`, `HasMZInRange` | same names |
| `removeMatchingItems`, `keepMatchingItems`, `moveMatchingItems` | `remove_matching_items`, `keep_matching_items` (+ `try_` forms), `move_matching_items` |
| `removeMatchingItemsUnroll`, `keepMatchingItemsUnroll` | `remove_matching_items_unroll`, `keep_matching_items_unroll` over `Identification` |
| `keepMatchingPeptideHits`, `removeMatchingPeptideHits`, `removeMatchingPeptideIdentifications` (map; list) | same names over `PeptideIdentificationMap`; the list form is `remove_matching_items` |
| `countHits` (template; list) | `count_hits` |
| `getBestHit` (template; list) | `best_hit` |
| `extractPeptideSequences` | `extract_peptide_sequences` |
| `extractUnassignedProteins` | `extract_unassigned_proteins` |
| `FilterPeptideEvidences`, `filterPeptideEvidences` (declared, never defined) | `filter_peptide_evidences`, `try_filter_peptide_evidences` |
| `removeUnreferencedProteins(ConsensusMap&, bool)` | `remove_unreferenced_consensus_proteins` |
| `removeUnreferencedProteins(vector<ProteinIdentification>&, peptides)` (and the list overload) | `remove_unreferenced_proteins` |
| `removeUnreferencedProteins(ProteinIdentification&, peptides)` | `remove_unreferenced_proteins(std::slice::from_mut(&mut run), peptides)` |
| `removeDanglingProteinReferences(peptides, proteins, bool)` | `remove_dangling_protein_references` |
| `removeDanglingProteinReferences(ConsensusMap&, bool)` | `remove_dangling_consensus_references` |
| `removeDanglingProteinReferences(ConsensusMap&, const ProteinIdentification&, bool)` | `remove_dangling_consensus_references_to_run` |
| `updateProteinGroups`, `removeUngroupedProteins` | `update_protein_groups`, `remove_ungrouped_proteins` |
| `removeEmptyIdentifications` (vector; list; map) | `remove_empty_peptide_identifications`, `remove_empty_protein_identifications`; `remove_empty_feature_identifications`, `remove_empty_consensus_identifications` |
| `filterHitsByScore(vector, threshold)` (and list; single ID) | `filter_peptides_by_score`, `filter_proteins_by_score` (a single ID via `std::slice::from_mut`) |
| `filterHitsByScore(vector, threshold, ScoreType)` | **not ported** (`IDScoreSwitcherAlgorithm`) |
| `filterGroupsByScore` | `filter_groups_by_score` |
| `keepNBestHits` (vector; list) | `keep_n_best_peptide_hits`, `keep_n_best_protein_hits` |
| `filterHitsByRank` (vector; list) | `filter_peptides_by_rank`, `filter_proteins_by_rank` |
| `removeDecoyHits` (vector; list) | `remove_decoy_peptide_hits`, `remove_decoy_protein_hits` |
| `removeHitsMatchingProteins`, `keepHitsMatchingProteins` (single ID; vector; list) | `filter_peptides_by_accessions`, `filter_proteins_by_accessions` with `MatchAction` |
| `keepBestPeptideHits` | `keep_best_peptide_hits` |
| `filterPeptidesByLength`, `filterPeptidesByCharge`, `filterPeptidesByRT`, `filterPeptidesByMZ`, `filterPeptidesByMZError` | `filter_peptides_by_length`, `_by_charge`, `_by_rt`, `_by_mz`, `_by_mz_error` |
| `filterPeptidesByRTPredictPValue` | `filter_peptides_by_rt_predict_p_value` (returns the hits removed for a missing value) |
| `removePeptidesWithMatchingModifications`, `keepPeptidesWithMatchingModifications` | `filter_peptides_by_modifications` with `MatchAction` |
| `removePeptidesWithMatchingRegEx` | `remove_peptides_with_matching_regex` |
| `removePeptidesWithMatchingSequences`, `keepPeptidesWithMatchingSequences` | `extract_peptide_sequences` + `filter_peptides_by_sequences` with `MatchAction` |
| `keepUniquePeptidesPerProtein`, `removeDuplicatePeptideHits` | `keep_unique_peptides_per_protein`, `remove_duplicate_peptide_hits` |
| `filterHitsByScore(AnnotatedMSRun&, ...)`, `keepNBestHits(AnnotatedMSRun&, n)`, `keepHitsMatchingProteins(AnnotatedMSRun&, fasta)` | **not ported** (`AnnotatedMSRun`) |
| `keepNBestSpectra` | `keep_n_best_spectra` |
| `keepNBestPeptideHits(map, n)` | `keep_n_best_hits_in_feature_map`, `keep_n_best_hits_in_consensus_map` |
| `keepBestPerPeptide`, `annotateBestPerPeptide` | `keep_best_per_peptide`, `annotate_best_per_peptide` |
| `keepBestPerPeptidePerRun(proteins, peptides, ...)`, `annotateBestPerPeptidePerRun(proteins, peptides, ...)` | `keep_best_per_peptide_per_run`, `annotate_best_per_peptide_per_run` (the protein runs do not affect the result) |
| `keepBestPerPeptidePerRun(map, ...)`, `annotateBestPerPeptidePerRun(map, ...)` | `keep_best_per_peptide_per_run_in_map`, `annotate_best_per_peptide_per_run_in_map` |
| `annotateBestPerPeptidePerRunWithData`, `annotateBestPerPeptideWithData` | not ported (see above) |
| `keepBestMatchPerObservation`, `filterObservationMatchesByScore`, `removeDecoys` (`IdentificationData`) | `keep_best_match_per_observation`, `filter_observation_matches_by_score`, `remove_decoys` |

### Source conventions kept in the wave-4 members

- `getBestHit` takes the orientation of the first identification with hits,
  replaces the best only on a strictly better score, and throws for two score
  types; `best_hit` does the same (`Error::InvalidValue`).
- `filterGroupsByScore` keeps a group only when its probability is strictly
  above the threshold (higher better) or at most the threshold (lower better):
  a group at the threshold is removed in one orientation and kept in the other.
- `filterPeptidesByRTPredictPValue` keeps hits whose value is at most
  `1 - threshold`.
- The best-per-peptide annotation sorts each identification, visits its first
  `nr_best_spectrum` hits, and leaves 1 on the hit that ends as its key's best,
  0 on every other visited hit and nothing on unvisited ones; a later hit wins
  only when strictly better in its own identification's orientation.
- `moveMatchingItems`, and so `extractUnassignedProteins`, reorder both parts
  exactly as `std::partition`'s two-pointer exchange does.
- `PeptideDigestionFilter` counts missed cleavages as the source's Boost
  tokenizer: a lookahead-only enzyme (Asp-N, Lys-N) at position 0 counts one
  more than `ProteaseDigestion::count_internal_cleavage_sites`.
- `DigestionFilter` passes `end - start` as the product length although the
  evidence end is inclusive, so a correct product loses its last residue and is
  usually rejected (`WVTFISLLLLFSSAYSR` at 2..=18 of the test protein); it
  applies the methionine allowance before the emptiness check. Both are kept.

### Native differences in the wave-4 members

- A non-numeric meta value is an error for `HasMaxMetaValue`,
  `HasMinMetaValue` and `filter_peptides_by_rt_predict_p_value`. The source's
  `DataValue::operator double()` returns the union's `double` member for a
  string or list, i.e. pointer bits, and compares that garbage.
- `remove_peptides_with_matching_regex` compiles the crate's Boost-compatible
  expression where the source uses `std::regex` (ECMAScript); constructs only
  one grammar has can differ. An invalid expression is an error in both.
- `filter_groups_by_score` requires a finite threshold.
- The `IdentificationData` filters iterate the native graph in semantic key
  order; the source's containers iterate in heap-address order of the
  referenced records (the oracle lists spectrum_3 before spectrum_1). The
  surviving records are the same.

## Verification

`tests/id_filter.rs` embeds the small scores/sequences/accessions from the pinned
IDFilter fixture independently of the Rust idXML reader. It reproduces source
cutoff, top-N, dense-rank, tied-best, modification, duplicate, coordinate,
accession, decoy and group-cleanup expectations. Additional tests cover signed
charges, zero-N behavior, run isolation and run union, metadata and sample-array
retention, map traversal, precursor tolerances, atomic errors, native equality
and the exact-duplicate work limit (tier 3 and 4).

`tests/id_filter_oracle.rs` replays the executed C++ oracle `../oracle/B8b-idfilter`
(Linux x86_64 Release install `bc9cc12-c19e494-174b576`, run twice on ibminode06,
byte-identical; tier 1) for every wave-4 member: `countHits`, `getBestHit` (16
variants and the score-type error), `filterGroupsByScore` (8), the RT p-value
filter (9), 13 regular expressions and an invalid one, all four best-per-peptide
functions and both map forms under 7 option combinations (1,267 annotated
rows), the meta-value and decoy predicates, `moveMatchingItems`,
`extractUnassignedProteins`, `PeptideDigestionFilter` over 5 enzymes and 7 bound
pairs, `DigestionFilter` over 24 evidences and 8 settings, and the three
`IdentificationData` filters. Hashes are in `tests/data/id_filter_provenance.json`.

Run `cargo test --features idxml --test id_filter_oracle` and
`cargo test --no-default-features --test id_filter`.

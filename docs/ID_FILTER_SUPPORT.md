# Identification filtering

`openms::analysis::id_filter` ports a substantive subset of the pinned
`OpenMS::IDFilter` algorithms to the existing native identification and feature
records. It uses free functions and standard-library collections. No additional
dependency, C++ build, search engine or identification graph is involved.

The source reference is OpenMS4-core revision
`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`:

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

This is not complete `IDFilter` class parity. Arbitrary predicate/template
wrappers, regex filtering, metadata-score-type discovery, RT prediction
p-value filtering, group score filtering, digestion-evidence predicates,
best-per-peptide/run annotation, unassigned-protein extraction and the separate
`IdentificationData` observation graph algorithms are not ported here.
The existing digestion and score-switching APIs can be used explicitly in a
larger workflow; this module does not silently invoke them.

## Verification

`tests/id_filter.rs` embeds the small scores/sequences/accessions from the pinned
IDFilter fixture independently of the Rust idXML reader. It reproduces source
cutoff, top-N, dense-rank, tied-best, modification, duplicate, coordinate,
accession, decoy and group-cleanup expectations. Additional tests cover signed
charges, zero-N behavior, run isolation and run union, metadata and sample-array
retention, map traversal, precursor tolerances, atomic errors, native equality
and the exact-duplicate work limit.

Run `cargo test --offline --no-default-features --test id_filter`.

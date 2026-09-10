# Basic protein inference

`openms::analysis::protein_inference::BasicProteinInference` ports the scientific
aggregation and grouping behavior of `OpenMS::BasicProteinInferenceAlgorithm`.
There is no separate `SimpleProteinInferenceAlgorithm` in the pinned source.
The implementation operates on native identification records and standard-library
maps and sets. It does not require a graph package or optional file-format feature.

Reference revision: `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`.

- `src/openms/include/OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h`
- `src/openms/source/ANALYSIS/ID/BasicProteinInferenceAlgorithm.cpp`
- `src/openms/source/ANALYSIS/ID/IDBoostGraph.cpp`
- `src/tests/class_tests/openms/source/BasicProteinInferenceAlgorithm_test.cpp`
- `src/tests/class_tests/openms/data/newMergerTest_out.idXML`

## Entry points and source defaults

All entry points validate and work on temporary records. Any returned error leaves
the caller's records unchanged. The result reports processed runs, removed peptide
and protein hits, resulting groups, and PSMs with greedily reduced evidence.

| Native method | Behavior |
| --- | --- |
| `run(peptides, proteins)` | Source vector overload: stable best-original-score selection to one candidate per spectrum, then optional score switching, inference per run, and original-score restoration. Pass a one-element protein slice for this convention on one run. |
| `run_single(peptides, protein)` | Source single-object selection convention: switch scores before sorting; aggregate only the best candidate, but retain other candidates unless dangling-evidence cleanup removes them. Peptide IDs belonging to other runs remain untouched. |
| `run_consensus_map(map, protein, include_unassigned)` | Use all top-level consensus-feature IDs as a union, ignoring their run identifiers, and optionally include unassigned IDs. Preserve original IDs on output and sort protein hits by inferred score. The supplied protein record represents the union; map-stored protein records are not replaced implicitly. |

Spectrum-level records, their coordinates, metadata, and surviving hits' stored
ranks remain attached. Neither vector nor single-object inference sorts protein
hits. No method recomputes peptide ranks or FDR. Inference scores and protein-group
`probability` fields are algorithm outputs, not newly calibrated confidence estimates.

| Option | Default |
| --- | --- |
| `aggregation` | `AggregationMethod::Best` |
| `min_peptides_per_protein` | 1 |
| `treat_charge_variants_separately` | true |
| `treat_modification_variants_separately` | true |
| `use_shared_peptides` | true |
| `skip_count_annotation` | false |
| `annotate_indistinguishable_groups` | true |
| `greedy_group_resolution` | false |
| `score_type` | None: current main score |

`score_type` accepts `ScoreType::Raw`, `PosteriorErrorProbability`, or `QValue`,
matching the source's explicit categories. Other main score names can be used
without switching. Category lookup and backup collision handling reuse the native
score-switching module. A switched score is restored exactly, including when
the source's relative-tolerance rule accepts a nearly equal pre-existing backup.
Score metadata created during switching remains attached.

## Aggregation

Each selected PSM enters a sequence/charge lookup. Sequence keys include named
modifications unless modification variants are collapsed; charge becomes zero
when charge variants are collapsed. The best score chooses one representative per
key. Exact score ties retain the first representative. Repeated spectra remain
in the peptide records; only protein aggregation uses these representatives.

Every representative contributes to all its protein accessions. Missing protein
accessions are skipped during scoring, as in C++; minimum filtering subsequently
removes dangling evidence. Distinct charge variants must have consistent accession
sets, making the source's lowest-charge-accession assumption explicit.

With shared peptides disabled, a missing `protein_references` annotation or the
value `non-unique` excludes a PSM from scoring. `unique` and `unmatched` are accepted;
the latter normally has no evidence and contributes to no protein. This setting
does not itself remove shared PSMs from peptide records or the grouping graph.

PEP scores are complemented to `1 - PEP` before aggregation. Protein score type
becomes `Posterior Probability`, with higher scores better. Otherwise the selected
score type and direction are retained.

| Method | Source behavior |
| --- | --- |
| `Best` | Maximum for higher-better scores, minimum for lower-better scores. Both source names `best` and `maximum` map here. |
| `Product` | Initialize at 1 and multiply strictly positive contributions. Zero and negative contributions are skipped but still counted. |
| `Mean` | Arithmetic mean of all representative contributions. This is the actual source method named `sum`, including its final division by count. |

`nr_found_peptides` counts representatives, not raw PSMs or necessarily distinct
unmodified peptide sequences. Charge and modification options change this count.
Protein minimum filtering uses computed counts even when count annotation is
suppressed. Scores and counts are not recomputed after greedy evidence resolution.

Search-parameter metadata records the four source `TOPPProteinInference:*`
settings and `InferenceEngine="TOPPProteinInference"`. `InferenceEngineVersion`
identifies this Rust package rather than reporting an unexecuted C++ version.

## Grouping and greedy resolution

Default grouping partitions referenced proteins by identical sets of selected
PSM neighbors. It does not collapse graph nodes by sequence, modification or charge.
Groups include referenced singleton proteins. Their score is `max(-1, member
scores)`, preserving the source's initialization even for negative scores or
lower-better score types. Group ordering uses the native source-compatible
`ProteinGroup::sort`: descending score, fewer members, then lexical accessions.

Greedy resolution first collapses these protein groups and groups PSMs sharing the
same parent protein/group set. Each shared cluster chooses its best parent by:

1. Highest protein/group score.
2. Target-containing group ahead of a decoy-only group at equal scores.
3. Most currently connected PSM nodes, including duplicate sequences.
4. Lexically first accession list for any remaining tie.

The source traversal for step 3 requests peptide-level nodes but accepts PSM nodes
in its default graph. The native implementation preserves that actual behavior.
Clusters are processed in their first-PSM input order, defining an order where the
C++ implementation uses hash-map iteration. Losing parent edges are removed from
the affected PSM evidence. Original groups survive when their members remain
referenced; retained proteins absent from a group receive singleton groups.

Unreferenced-protein removal considers every remaining candidate, which matters
for `run_single`. Enabling greedy resolution with group annotation disabled still
resolves evidence and removes unreferenced proteins, but emits no groups.
Existing `protein_references` and target/decoy annotations are not recalculated
after resolution, matching the source. Consumers needing resolved associations
should use the surviving evidence rather than these original annotations.

## Checked differences and current boundaries

- Run identifiers must be unique. Accessions must be nonempty and unique within
  each run. `run` requires every peptide ID to reference a supplied run. Empty
  strings can identify one exact run; they are not wildcards. `run_single` ignores
  unrelated runs instead of applying the source's potentially destructive global
  cleanup to them.
- Score names and directions must agree among nonempty IDs within each run.
  Different runs may use different score types. Empty IDs do not determine a
  run's score type. This avoids the C++ implementation's global-first-ID coupling.
- PEP and posterior probabilities require the matching direction and values in
  0..=1. Greedy resolution requires higher-better scores or PEP; the C++ graph
  explicitly lacks lower-better support and otherwise maximizes those scores too.
- Greedy resolution requires known target/decoy annotations on referenced proteins.
  Invalid sharedness metadata is rejected when that annotation is used.
- Unlike C++, native records never receive infinity or NaN. Minimum zero is valid
  when resulting scores are finite. Unreferenced `Best`/`Mean` proteins that survive
  resolution cause an error; unreferenced `Product` proteins have the source score
  1. Undefined `Best` scores compare below finite scores during higher-better greedy
  resolution, so a negative finite score still beats an unscored protein. An
  unscored `Mean` protein cannot participate in greedy comparison because the
  source's NaN ordering is undefined. Unreferenced records outside the graph can
  still be removed by greedy cleanup.
- Floating aggregation overflow is an error. Product underflow follows finite
  floating arithmetic and may yield zero. Deterministic sequence iteration can
  change final floating rounding relative to C++ unordered iteration.
- Suppressing count annotation does not disable minimum filtering. Existing count
  metadata is retained when suppressed, but is never used instead of computed
  counts. This corrects the source's missing/stale-count dependency.
- Consensus-map unassigned IDs excluded by `include_unassigned=false` remain
  untouched. Some C++ top-N and dangling-reference helpers also modify excluded
  unassigned IDs. Map-stored protein runs likewise remain caller-managed.
- Every entry point replaces the supplied run's old groups with inferred groups.
  This matches the source vector/single-object paths and makes the consensus path
  repeatable; the C++ consensus overload does not explicitly clear previous groups.
- The full `IDBoostGraph` public API, experimental-design graph hierarchy,
  probabilistic inference engines and protein quantification are separate work.
  This module implements the graph operations needed by BasicProteinInference.

Input defaults allow at most 1,000,000 total peptide/protein hits and 5,000,000
peptide evidence records. `max_input_hits` and `max_input_evidences` are positive,
adjustable limits checked before cloning selected records. Work depends on input
sorting, representative keys and evidence edges; there is no protein-pair matrix.
Temporary records preserve transactional errors at the cost of copying selected
records. The caller's pre-existing metadata and sequence storage are not counted
as bytes by these limits.

## Validation

`tests/data/protein_inference_merger.tsv` preserves all scientifically relevant
input fields used by the source BasicProteinInference tests. Its provenance JSON
records the source revision, hashes, extraction and omitted unrelated fields.
The XML source's invalid all-zero date is not routed through the checked native
idXML parser. Tests run with optional format features disabled too.

Source expectations include protein scores/counts, four default groups, three
resolved groups, shared-peptide exclusion and RAW-score restoration. The minimum-0
source's unreferenced negative infinity is tested as an atomic native error;
minimum-1 tests preserve the remaining finite scores. The original minimum-0 greedy
fixture produces the source's complete finite result unchanged.

Independent cases cover mean/product edge arithmetic, representative and candidate
ties, charge/modification collapsing, PSM-based grouping and greedy tie counts,
negative scores against undefined scores, switched-score backup preservation,
single-run and consensus-map conventions, independent runs, empty input, invalid
metadata, resource limits and mutation atomicity.

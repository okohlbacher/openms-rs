# Identification conflict resolution

`analysis::id_conflict_resolver` ports the represented feature, consensus and peptide-list operations of `IDConflictResolverAlgorithm` from OpenMS4-core revision `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. It uses standard-library collections and existing native records. It resolves supplied annotations; it does not combine search-engine scores or infer new peptide sequences.

## Within-feature resolution

`resolve_feature_map` and `resolve_consensus_map` take a `ResolutionMethod`. They annotate each top-level feature with its unique ID under `feature_id`. Previously unassigned peptide identifications receive the string `feature_id="not mapped"`. Rejected feature annotations append to that unassigned list. Subordinates, feature geometry, intensities, protein records and other metadata remain attached.

| Method | Selection and retained information |
| --- | --- |
| `BestScore` (default) | Sort hits by their declared score direction, reduce every ID to its best hit, select the best ID and move the others to unassigned. All affected IDs receive their feature's ID as string metadata. |
| `KeepMatching` | Select the best top hit after sorting, swap its ID to the first position, and keep other IDs containing that same modified sequence at any charge. Each matching secondary ID keeps its best occurrence of that sequence. The winning ID retains **all** its sorted hits, matching the actual C++ overload. Rejected IDs also retain all sorted hits and receive `feature_id`; kept IDs do not receive a new annotation. |
| `RankAggregation` | Rank each unique modified sequence by its first zero-based position in every sorted ID. Missing sequences receive a penalty equal to the largest hit count. Choose the highest `1 − total_rank / (max_hits × number_of_IDs)`, then retain the best original-scoring occurrence of that sequence. The aggregate is used for selection; the retained hit keeps its original score. Rejected IDs keep their top hit. All affected IDs receive `feature_id`. |

`resolve_identifications` applies the same selection to a peptide-ID vector and a supplied feature ID, appending rejected records to a separate vector. This direct adapter preserves previously unassigned records verbatim, since it has no enclosing map to annotate.

Score sorting and selection ties retain input order. Rank aggregation counts each sequence once per ID, ignores charge when grouping that sequence, and includes empty IDs in the missing-sequence penalty. If every ID is empty, the first survives and the rest are unassigned. If populated IDs exist, empty IDs cannot win best-score/keep-matching selection. This avoids the source's lower-better empty-vector selection and unchecked hit dereference.

Populated competing IDs must agree on score type and direction. A rank winner's original scores cannot safely be compared across arbitrary scales. Empty IDs do not impose a score convention. These checks are stronger than the source's assumption that the first record describes every other record.

## Between-feature resolution

`resolve_between_features` and `resolve_between_consensus_features` require at most one peptide identification per top-level feature. They sort its hits, group by **feature charge** and modified top-hit sequence, and retain each group on the highest-intensity feature. The first feature wins equal intensities. The complete sorted identification from a losing feature moves to unassigned, including any alternative hits. Features themselves are retained; empty IDs and features without IDs are untouched. This operation does not change `feature_id` metadata.

Run within-feature resolution first when features have multiple IDs. The `KeepMatching` mode can intentionally leave several IDs on a feature and therefore does not necessarily satisfy this precondition.

## Repeated spectrum annotations

`reduce_to_one_per_spectrum` operates on one run's peptide identifications. It keys each nonempty ID by its spectrum reference, **stored first hit** modified sequence, and hit charge. It does not sort hits. It keeps the best main score per key, retains the first exact-score tie, and preserves surviving input order.

Different peptidoforms or charges on the same spectrum remain separate. Missing references and empty IDs remain untouched. `UnresolvedIdentifications` reports removed IDs, spectra with multiple remaining peptidoforms, nonempty IDs without references, and groups with inconsistent score directions. An inconsistent group is counted as one peptidoform for the chimeric-spectrum statistic, as in the source. The native report also counts inconsistent score types; those groups remain intact. Mixed run identifiers are rejected before mutation to prevent identical local spectrum references from different files being conflated.

The first reduced group is reported as `<reference> / <sequence> / charge <n>`. Key order follows source sequence-length ordering before terminal/residue comparison. Native modification names and immutable registry IDs replace the C++ residue-modification pointer ordering, making modified-sequence ties reproducible across processes. Distinct native modification records are never silently merged because their allocation addresses or names happen to compare alike.

## Validation and evidence

Every public mutation validates inputs and computes its complete result before committing. A malformed later feature, incompatible score scale, invalid subordinate structure, or count overflow leaves all supplied data unchanged. Operations allocate owned temporary records and collections proportional to the input; they have no external I/O or C++ dependency.

`tests/id_conflict_resolver.rs` reproduces the source spectrum-reduction cases, both score directions, report counts, original-order preservation, the rank-aggregation example, and intensity-based feature resolution. The source rank example uses `SEQB`; the Rust test substitutes `SEQC` because B is outside the current mass-capable peptide model, preserving the independent rank calculation (`5/6` versus `4/6`). Tests extend the original intensity example to actually include its distinct-charge and modified-sequence features. Additional cases check missing penalties, duplicate sequences, sequence-length ties, keep-matching alternative preservation, empty-ID handling, metadata/subordinates, consensus adapters and atomic late failures.

Pinned sources: `src/openms/{include/OpenMS,source}/ANALYSIS/ID/IDConflictResolverAlgorithm.{h,cpp}`, `src/tests/class_tests/openms/source/IDConflictResolverAlgorithm_test.cpp`, and `src/openms/source/CHEMISTRY/AASequence.cpp`. These are recorded in the [source inventory](source-inventory.json). The C++ reference was inspected, not built or executed. Source and derived code retain BSD-3-Clause attribution. Broader consensus-identification score combination and the separate `IdentificationData` graph remain distinct port work.

Caller-owned chemistry is included in conflict keys: matching names, full IDs or
accessions do not equate records with different formulas or masses. The previous
length/name/ID tie order is retained, with complete annotation value ordering as
a final tie-break. Cloning these keys shares known records instead of copying
their chemical data. Anonymous tags retain owned strings.

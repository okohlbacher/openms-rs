# False discovery rates and q-values

`analysis::false_discovery_rate` ports target/decoy calculations and score application from OpenMS4-core commit `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. It supports the represented peptide/protein identification records, annotated unmodified-peptide estimates, indistinguishable protein groups, explicit-affix picked protein FDR, posterior probability estimates, and target/decoy ROC area. It uses only the Rust standard library. The wider core port remains in progress.

```rust
use openms::analysis::false_discovery_rate::FalseDiscoveryRate;
use openms::identification::PeptideIdentification;

fn calculate_psm_q_values(ids: &mut [PeptideIdentification]) -> openms::Result<()> {
    // Each selected hit must carry target_decoy metadata.
    FalseDiscoveryRate::default().apply_peptides(ids, false)
}
```

## The two source calculations differ

`FalseDiscoveryRate` exposes the actual legacy and Basic source calculations separately. Its default output is a q-value; `output: FdrOutput::Fdr` requests a raw FDR. In this table, T and D are cumulative target and decoy contributions at a score threshold.

| Calculation | Implemented formula and behavior |
| --- | --- |
| `calculate_legacy(targets, decoys, higher)` | D/T at inclusive target-score thresholds, without pseudocounts. Q-values are the cumulative minimum from worse to better thresholds, capped at one. |
| `calculate_basic(observations, higher)`, `conservative: true` | **(D+1)/(T+1)**, using the pinned implementation's denominator. |
| `calculate_basic(observations, higher)`, `conservative: false` | **(D+1)/(T+D+1)**. |
| `calculate_estimated(probabilities, higher)` | Running mean posterior error probability, or the corresponding mean error probability from posterior probabilities. Exact score ties share the value at their final cumulative endpoint. |

The Basic formulas above deliberately describe the code, whose `j + 1` denominator includes an additional observation beyond the cumulative count. Several source comments instead say `(D+1)/T` or `(D+1)/(T+D)`, and the class header reverses the conservative flag. The native implementation does not silently replace the actual computation with those descriptions. Raw FDR can exceed one; q-values use a cumulative minimum initialized to one.

Legacy thresholds use exact score equality. Basic groups scores whose absolute difference from the group's first score is at most `1e-12`. This is not transitive neighbor chaining. `ScoreLabel` permits fractional target contributions in `[0,1]` for Basic; D accumulates `1 - target_fraction`. Record-level methods use binary labels: peptide `target+decoy` contributes as a target, matching the source.

`ScoreCurve` exposes ascending `(original_score, value)` entries and a checked lookup. Legacy curves include supplied target and decoy scores and require exact lookup. Basic lookup retains the source cutoff rule: the next greater-or-equal stored score for higher-better input, or the next smaller-or-equal score for lower-better input, clamped at the smallest score. Missing upper coverage is an error rather than an unchecked iterator dereference.

Legacy decoy values follow the pinned source's target-score assignment. For q-values this is a nearest target-score lookup, with an equal-distance tie going to the worse target. Raw FDR retains a different source target order, so its decoy assignment generally selects the best target's value, or the worst target's value when the decoy is at least as good as every target. A shared raw target/decoy score can therefore overwrite the target's D/T entry. This surprising behavior is tested explicitly. Use the Basic API when that source algorithm is intended; the APIs are not interchangeable.

No-decoy legacy targets receive zero. A legacy curve with only decoys assigns one to their scores. Empty score arrays produce empty curves. Basic retains its pseudocount behavior even without decoys.

## Applying scores to identifications

| Method | Counting and output behavior |
| --- | --- |
| `apply_peptides(ids, annotate_peptide_fdr)` | Legacy concatenated search. Sorts each record by its original direction and retains its best hit unless `use_all_hits` is true. Optionally pools by run identifier and/or observed charge. Finally sorts retained hits by their new lower-better values. |
| `apply_basic_peptides(ids)` | Basic concatenated search. Counts each first hit unless `use_all_hits`, then scores all remaining hits against their pool's curve. First-hit mode requires records to be sorted already. |
| `apply_separate_peptides(targets, decoys)` | Legacy separate searches. Counts all hits, independent of `use_all_hits`, without requiring target/decoy labels. Reverse records are changed only when `add_decoy_peptides` is true. |
| `apply_proteins(ids)` | Legacy concatenated protein search over all supplied runs. Updates protein hits, retaining input order. |
| `apply_separate_proteins(targets, decoys)` | Legacy separate protein searches. Counts all hits and changes only forward records. Reverse records are read-only. |
| `apply_basic_protein(id, groups_too)` | Basic protein FDR and, when requested, a separate calculation for indistinguishable groups. |
| `apply_picked_protein(id, affix, groups_too)` | Basic calculation after target/decoy competition by paired accession, plus the supported source group calculation. |
| `apply_estimated_protein(id)` | Posterior-probability estimates for one protein identification record. |
| `apply_basic_peptide_level(ids)` | Basic best-representative calculation over unmodified peptide sequences. |

Concatenated target/decoy methods require recognized labels. By default, decoy peptide/protein hits are removed. The corresponding `add_decoy_peptides` or `add_decoy_proteins` flag retains them. The legacy peptide special case with no targets or no decoys removes decoys even when retention is requested, and assigns zero to any targets. It does not add the optional peptide-level annotation in that special case.

`split_charge_variants` and `treat_runs_separately` apply to the concatenated peptide methods. The native implementation builds pools from observed charges/identifiers, rather than iterating declared search-charge ranges. Each pool must contain one score type and direction; different independent pools may have different types/directions. Basic selection and assignment use the same pool boundaries. A pool containing secondary hits but no selected first-hit score is an error. This avoids source paths that reuse already-converted scores, annotate unselected runs, or dereference an empty lookup. The protein and separate-search methods retain their source pooling scope and do not use these peptide-only flags.

Before replacing each retained hit's score, ordinary methods store the original value as typed float metadata under `<old_score_type>_score`. Score type becomes `q-value` or `FDR`, and direction becomes lower-better. Existing metadata at that original-score key is overwritten, following the source. Ranks are retained; they are not reassigned. The record's `significance_threshold` also remains unchanged and is not converted into an FDR threshold. Other fields, acquisition references, protein evidence, analysis results, annotations, search parameters, paths, and metadata retain their values unless a hit is removed. Protein group cleanup after hit removal is a separate operation.

## Peptide-level estimates

The legacy `annotate_peptide_fdr: true` option computes an additional curve from the best score per unmodified sequence, separately for target and decoy sequences. It stores `peptide q-value` or `peptide FDR` on retained PSMs while keeping the PSM-level FDR as the primary score. Modified forms and charge states of the same unmodified sequence share that peptide representative within their selected pool.

`apply_basic_peptide_level` selects one best first-hit representative per unmodified sequence across all input records. A better duplicate replaces both the score and label; targets win exact score ties. Every retained record for that sequence receives its representative's result. This implements the pinned regression fix for input-order dependence. This API requires at most one hit per record and does not split runs or charges. The source only replaces the first hit but changes the whole record's score type, which otherwise leaves mixed score units in multi-hit records. The native API rejects that ambiguous case and correctly sets lower-better direction, which the source setter omits.

For this method alone, the source original-score key is the old score type **without** `_score`. A backup into `target_decoy` is rejected before committing changes, because a numeric backup would destroy the required label. Other existing original-score metadata is overwritten as in the source. The primary score type becomes `peptide q-value` or `peptide FDR`. A removed decoy leaves an empty identification record; its other record fields remain intact.

## Proteins, groups, and picked competition

Basic group calculation applies to `indistinguishable_groups`. `protein_groups` remains unchanged, matching the source method. A group is target if any listed accession corresponds to a target protein; all groups remain present even when decoy hits are removed. Group `probability` becomes the FDR/q-value; group quantity arrays and accession lists remain unchanged. Groups have no independent score-type field. When `groups_too` is false, group scores remain in their original units.

Picked competition requires `DecoyAffix::Prefix("decoy_")` or an explicit suffix. Labels determine whether a protein hit is decoy; each decoy accession must match the affix, and removing it must leave a nonempty accession. The best score per target/decoy accession pair contributes to the curve. Exact ties prefer a target. Both winning and losing target hits can remain in output and receive scores from the picked curve, as in the source; decoy retention is controlled separately.

Picked groups preserve a source detail: when a group's accession order encounters a picked decoy before a picked target, that group can contribute both a decoy and a target observation. Encountering the target first contributes only a target. This order-sensitive source behavior is covered by a regression test and is not described as a general statistical group model. Missing protein references, duplicate accessions, and unavailable group-score lookup coverage are errors. Automatic decoy-affix inference, with the source helper's permissive prefix/suffix fallback, is `apply_picked_protein_auto` (see "Source-exact apply members").

## Probability estimates and ROC

`apply_estimated_protein` accepts `Posterior Error Probability` with lower-better direction or `Posterior Probability` with higher-better direction, all in `[0,1]`. It processes one explicit record rather than silently ignoring all but the first record of a vector. It writes `Estimated Q-Values`, lower-better, and keeps the original scores under the usual suffixed key. Target/decoy labels are unnecessary when all hits are retained; removing decoys requires labels.

The pinned estimated-q implementation reserves a vector without sizing it before indexed writes, so its behavior is undefined. The native routine implements the intended running probability mean with safe allocation and gives equal scores the cumulative endpoint value. Its numerical tests use independent probability examples; no equivalence to execution of that invalid C++ memory access is claimed. This method always computes estimated q-values and does not use `output` or the conservative flag.

`roc_n` takes score/label observations and sorts them in the supplied direction. Labels must be binary. It integrates the target/decoy count curve by trapezoids, including an entire exact-score tie batch when it reaches the optional false-positive cutoff. `None` selects the full area. An empty input returns zero and target-only input returns one. A missing target class, or reaching the cutoff before any target, returns an error instead of division by zero. Partial target labels are supported in Basic FDR but not this ROC API.

## Evaluation members (source-exact, `bc9cc12`)

These members follow the source arithmetic exactly and are checked bit for bit
against the Linux x86_64 Release build (`tests/false_discovery_rate_evaluation.rs`,
`tests/data/fdr_evaluation_oracle.tsv`, `../oracle/I2-fdr`; 180 results).

| C++ | Rust | Notes |
|---|---|---|
| `double diffEstimatedEmpirical(const ScoreToTgtDecLabelPairs&, double pepCutoff = 1.0) const` | `diff_estimated_empirical(&[(f64, f64)], f64)` | Pairs sorted best first. Empty input returns 1. The last entry always uses `decoys / (targets + decoys)` whatever `conservative` says, as in the source. |
| `double rocN(const ScoreToTgtDecLabelPairs&, Size fpCutoff = 50) const` | `roc_n_sorted(&[(f64, f64)], usize)` | No sorting or validation; `decoys * targets` is a 32-bit product; a cutoff before any target divides by zero as the source. The older checked `roc_n` stays. |
| `double rocN(const PeptideIdentificationList&, Size) const`, and with `identifier` | `roc_n_peptides(ids, fp_cutoff, Option<&str>)` | First hits unless `use_all_hits`; direction of the first record; `fp_cutoff = 0` is the full curve. |
| `double rocN(const ConsensusMap&, Size, bool) const`, and with `identifier` | `roc_n_consensus(map, fp_cutoff, Option<&str>, include_unassigned)` | Feature identifications, then unassigned ones. The identifier overload reads the first feature's first identification unchecked in the source; absent, it is an error here. |
| `double applyEvaluateProteinIDs(const ProteinIdentification&, double, UInt, double) const` | `apply_evaluate_protein_ids` | Requires `Posterior Probability`. `(1 - diff)(1 - w) + rocN w`. |
| `double applyEvaluateProteinIDs(const std::vector<ProteinIdentification>&, ...) const` | `apply_evaluate_protein_ids_first` | First run only, and `diff w + rocN (1 - w)` (defect candidate below). |
| `double applyEvaluateProteinIDs(ScoreToTgtDecLabelPairs&, ...) const` | `apply_evaluate_score_labels(&mut [(f64, f64)], ...)` | Sorts the pairs in place, as the source. |

Defect candidate: the vector overload of `applyEvaluateProteinIDs` combines
the two measures as `diff * diffWeight + rocN * (1 - diffWeight)`, a score
that grows with the estimated-vs-empirical difference, while the other two
overloads compute `(1 - diff) * (1 - diffWeight) + rocN * diffWeight`. The port
reproduces both.

## Source-exact apply members (`bc9cc12`)

The native `apply_*` methods above check their input and refuse ambiguous
records. The members below instead compute what the source computes,
including its known quirks, and are checked bit for bit against the Linux
x86_64 Release build: `tests/false_discovery_rate_basic_oracle.rs` (350
records from `../oracle/I2-fdr/driver2.cpp`, including the cases where the
source crashes or throws, run in forked children) and
`tests/false_discovery_rate_observation_oracle.rs` (48 runs, 3,872 match
values from `driver3.cpp`). Provenance is in
[fdr_apply_provenance.json](../tests/data/fdr_apply_provenance.json).

| C++ | Rust |
|---|---|
| `applyBasic(PeptideIdentificationList&, bool, int charge, std::string identifier, bool only_best_per_pep)` | `apply_basic_selected` |
| `applyBasic(const std::vector<ProteinIdentification>& run_info, PeptideIdentificationList&)` | `apply_basic_runs` |
| `applyBasic(ConsensusMap&, bool)` | `apply_basic_consensus` |
| `applyBasicPeptideLevel(PeptideIdentificationList&)` | `apply_basic_peptide_level_source` (the checked `apply_basic_peptide_level` stays) |
| `applyBasicPeptideLevel(ConsensusMap&, bool)` | `apply_basic_peptide_level_consensus` |
| `applyPickedProteinFDR(id, "", prefix, groups_too)` | `apply_picked_protein_auto` |
| `DecoyStringHelper::findDecoyString(const ProteinIdentification&)` | `find_decoy_string`, `DecoyStringResult` |
| `applyToObservationMatches(IdentificationData&, ScoreTypeRef)` | `apply_to_observation_matches` |
| `DefaultParamHandler` parameters | `defaults`, `to_param`, `from_param`, `from_param_with_warnings` (the seven source flags; `max_*` limits are native) |

Where the source dereferences `map::end` or reads the first hit of an empty
record, the port returns an error; the oracle records those cases as
`crash`/`error` and the tests require an error there. Everything else,
including the log lines the source writes, matches.

Defect candidates in the source, reproduced by the port:

- `applyBasic(run_info, ids)` without `treat_runs_separately` iterates the
  charges from the smallest range start to the largest range *start*: it
  reads `.first` for both ends of the search charge range.
- `applyBasic(ConsensusMap&)` does not clear the collected scores between runs
  or charges, so each later curve also counts the earlier ones, and it marks
  rescored records with the data's direction instead of lower-better.
- `applyBasicPeptideLevel(ConsensusMap&)` takes the direction of the *last*
  feature with an identification (its loop never breaks).
- `applyBasic(ids, higher, charge, identifier, ...)` rescores every hit of
  every record, whatever its run or charge; `applyBasic(run_info, ids)` calls
  it per run and charge, so the later calls read the earlier q-values as
  scores.
- `applyToObservationMatches` counts the best match per observation through
  `getBestMatchPerObservation`, which breaks equal scores by
  `IdentifiedMolecule::operator<`, a comparison of heap addresses. A copy of
  the same `IdentificationData` can therefore count a different match. A
  first oracle fixture with such ties gave a different winner in the first of
  eight copies of the same graph, in two of six graphs (26 and 5 changed
  values). The port takes the first match in registration order, and the
  committed fixture keeps scores distinct within an observation (ties across
  observations stay).

## Validation, limits, and remaining scope

All source scores and numeric record fields must be finite; posterior probabilities and target fractions have explicit domain checks. Concatenated score pools require compatible score types/directions. Every record mutation is prepared on owned copies and committed only after the complete operation succeeds, including separate-search pairs and group computations. Reported errors cannot leave partially changed scores, metadata, or filtered hit lists.

Defaults permit one million hits and one million identification records per operation, plus one hundred thousand peptide pools or indistinguishable groups. Separate-search hit counts and record counts are checked across both inputs. Zero limits are invalid. Curves and grouping use standard ordered maps and sorting, with storage linear in the supplied records/hits. Limits bound these counts rather than annotation byte sizes or the cost of nested record validation. Allocation failure follows the Rust allocator's normal behavior. The checked code avoids the legacy decoy-assignment nested scan without changing its score-order rule.

Every public member of `FalseDiscoveryRate.h` now has a Rust counterpart (see the two previous sections). `applyEstimated` is the one member with a deliberate difference: the source writes through a reserved but unsized vector, which is undefined behaviour, so `apply_estimated_protein` implements the intended running mean and is not claimed equal to the source's execution. Score naming uses the existing source strings rather than a new ontology mapping.

## Verification and provenance

The tests extract numeric and target/decoy fields from pinned upstream idXML files without building C++ or depending on the in-progress Rust idXML reader. [fdr_provenance.json](../tests/data/fdr_provenance.json) records source and fixture hashes and each extraction. Peptide sequences are not needed for the extracted PSM/protein score tests; dedicated sequence tests cover modified/unmodified peptide representatives.

Tests check OMSSA's 1,534 PSMs including q-values `0.0730478589420655` and `0.409926470588235`, the target+decoy label, removed decoys, X!Tandem's `0.08` PSM boundary and `0.897384` protein value, and the six picked-protein values `0.25, 0.25, 0.25, 0.4, 0.4, 0.5`. The X!Tandem protein inputs use the upstream `withProtScores` fixtures: the unscored files used by some other source test sections have all-zero protein scores and cannot exercise their advertised nonzero score thresholds. The native tests assert that the chosen threshold conditions actually occur.

Additional tests cover the source peptide duplicate/tie regression, actual Basic formulas, raw legacy shared-score behavior, score direction, fractional labels, grouping tolerance, first-hit counting, run/charge pools, optional peptide annotations, separate searches, group metadata and picked-group ordering, probability estimates, ROC cutoffs, missing labels, limits, and atomic failure. Extracted fixtures and derived Rust code retain BSD-3-Clause provenance.

Pinned sources:

- [FalseDiscoveryRate implementation](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/ANALYSIS/ID/FalseDiscoveryRate.cpp), [API](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/include/OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h), and [tests](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/FalseDiscoveryRate_test.cpp).
- [IDScoreGetterSetter implementation](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/ANALYSIS/ID/IDScoreGetterSetter.cpp) and [template score extraction/replacement](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/include/OpenMS/ANALYSIS/ID/IDScoreGetterSetter.h).

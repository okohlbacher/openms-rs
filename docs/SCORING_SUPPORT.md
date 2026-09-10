# Identification score handling and fragment scores

`analysis::scores` and `analysis::psm_scoring` port score handling and peptide-spectrum scoring from OpenMS4-core revision `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. They use the existing native spectra, peptide/protein records and typed metadata. They are separate from false-discovery estimation and database search.

## Categories and score switching

`ScoreType` represents the source's six categories: raw search-engine score, raw E-value, posterior probability, posterior error probability, false-discovery rate and q-value. The registry contains all 29 pinned names, including PSI accessions, in the source's category and lexical set order. Higher values are better for raw scores and posterior probabilities; the other categories default to lower values being better. This is the source category convention, not an assertion about every possible external search-engine score.

Exact name lookup is case-sensitive, following `Scores.cpp`; the C++ header incorrectly describes this operation as case-insensitive. `normalize_score_name` removes one literal `_score` suffix. `ScoreType::matches` and `from_normalized_name` use that normalization; `from_name` performs exact lookup without stripping. Category `parse` follows the source's separate rules: remove the literal suffix, lowercase and remove spaces, hyphens and underscores. Category strings such as `q-value`, `Raw E-Value` and `Posterior_Error_Probability` are accepted. A search-engine name such as `hyperscore` is a score name, not a category string.

`find_peptide_score` and `find_protein_score` check the main score category first. Otherwise they inspect only the first hit's metadata, trying each registered name and then its `_score` form in lexical source order. A found value is not assumed numeric until switching. Absence returns `None`.

`ScoreSwitcher` reads the requested numeric metadata value for every hit, stores the previous main score and updates the record's main score name and direction. It does not sort candidates or rewrite ranks. The old-score metadata key defaults to the previous main score type; the new score type defaults to the new metadata key. Both names can be supplied explicitly. Integer and float metadata are accepted; numeric-looking strings are rejected by the native typed-value conversion. Metadata units are not converted or interpreted as score scaling factors.

If the old-score key already contains an equivalent numeric value, it remains intact. If it contains a different value, the previous score is stored under that key plus `~`. Equivalence uses the source's symmetric relative tolerance of `1e-6`, evaluated with scaling to avoid intermediate overflow. The native implementation rejects a conflicting occupied backup instead of overwriting it, and retains the exact type/unit of an already-compatible backup. Missing/empty old values are replaced by an owned float. The reserved `target_decoy` label cannot be used as a numeric backup key. The requested new score is read before any backup operation.

Single-category helpers examine every identification independently. This corrects the C++ list shortcut that assumes every record is already correct when the first record has the requested category. Records already in that category remain unchanged, including their declared direction. Other records must provide a matching score on every hit; the selected metadata key becomes the new score type after removing `_score`. An empty record with a different main category cannot supply a score and returns an error. Plain `ScoreSwitcher` can update the type/direction of an empty record without processing hits.

Switching peptide/protein slices, feature maps and consensus maps is transactional: changes are made to owned temporary results, then committed. A missing score or backup conflict in a later record leaves all inputs unchanged. This requires temporary memory proportional to the affected input. Map adapters visit top-level features and optionally unassigned peptide IDs. Subordinate peptide IDs and protein scores are not implicitly changed. To restore an earlier main score, switch explicitly to its saved metadata key and original direction.

## Morpheus scoring

`MorpheusScore::compute` and `compute_with_charges` preserve the source's two ordered traversals. One counts theoretical peaks, allowing an experimental peak to support multiple theoretical peaks. The other sums each matching experimental intensity once, allowing a theoretical peak to account for multiple experimental peaks. The score is `matched_theoretical_peaks + matched_experimental_intensity / total_experimental_intensity`.

Mass matches include the tolerance boundary. Ppm tolerance is relative to theoretical m/z, with the source's float64 arithmetic. The charge-aware variant requires exact charge equality, including zero; a charge mismatch still advances the pointer chosen by the source mass-matching branch. It does not search further for a same-charge alternative. Consequently the two passes can report different match counts.

`MorpheusResult` retains source float32 result precision after float64 accumulation. It reports theoretical peak count, match count, matched/total ion current and mean absolute errors in Da and ppm. The source divides both accumulated experimental error sums by the **theoretical match count**, not by the number of experimental matches. This convention is preserved. Nonempty inputs with no counted matches report error sentinel `1e10`; if either spectrum is empty, the entire result is zero, as in source.

Inputs must be finite, sorted and nonnegative in m/z and intensity. A nonempty experimental spectrum with zero total intensity returns an error instead of the source's NaN score. A matched theoretical m/z of zero cannot produce a finite ppm error and is rejected, including in absolute-tolerance mode because Morpheus computes both errors. Overflow, incompatible charge lengths and invalid tolerances are errors.

## HyperScore variants

The source has scientifically distinct overloads. The native API keeps them distinct:

| Method | Matching and ion convention |
| --- | --- |
| `compute` | Source `MatchedIterator` nearest matching with inclusive tolerance; b/y counts, including names containing `$b`/`$y`; every matched intensity contributes to the dot product |
| `compute_with_detail` | Same matching, but a/b/c prefix and x/y/z suffix counts; includes first-dollar cross-link names; returns counts and source mean mass error |
| `compute_with_charges` | Float64 nearest experimental peak, strict tolerance boundary and exact charge match; only b/y ions contribute to counts and intensity |
| `compute_with_intensity_sum` | Charge-aware matching, ordinary leading b/y names with parsed positive ordinals; duplicate ordinals combine experimental intensity and contribute once per series when their sum is positive |

The score is `ln(1 + dot_product) + ln(prefix_count!) + ln(suffix_count!)`. Integer log-factorials are summed with standard-library logarithms, avoiding a new numerical dependency. This adds linear work bounded by the matched peak count. The pinned implementation uses `lgamma`; small floating-point rounding differences are checked against explicit source tolerances.

Basic/detailed matching shares the existing comparison module's source `MatchedIterator` traversal. Tolerance, distances and ppm reference coordinates use float32. Equal-distance ties keep the lower/current target, and duplicate target positions preserve the source's forward-traversal stopping behavior. A target can be reused. Extremely large coordinates/tolerances that cannot fit the source float32 representation return errors. Charge-aware matching instead uses the kernel's float64 nearest lookup, including its lower-m/z tie rule.

For basic/detailed scoring, source float32 intensities are multiplied before float64 summation. Charge-aware overloads first promote intensities to float64. That distinction is preserved: very large finite intensities can overflow the basic product while remaining valid in the charge-aware calculation. Nonfinite results return errors.

Theoretical annotations come from the **first string data array**, following source behavior, and must be aligned and nonempty. The array's name is not used to silently select a different array. In the basic b/y classifier, a `$y` match takes precedence over b, matching the source conditions. The detailed classifier first uses the leading ion type, then the character after the first `$`. Its error numerator includes every matched peak, including unclassified annotations, but its denominator counts only classified prefix/suffix ions; this source peculiarity is preserved. Error units follow the configured tolerance.

The intensity-support overload accepts an existing array of 1–100,000 finite values. A b ordinal `i` adds at `i - 1`; a y ordinal `i` adds at `peptide_length - i`. Existing values are incremented. Out-of-range ordinals and arithmetic errors leave the entire array unchanged, correcting unchecked source indices. Cross-link names are not interpreted by this ordinal overload. Empty spectra produce zero score without changing the existing array. The charge-aware methods do not use unknown charge as a wildcard.

## Evidence and remaining scope

`tests/scores.rs` covers the complete registry, source lookup priority, switching/restoration, backup collisions, zero/opposite-sign/large-value comparisons, heterogeneous lists, typed metadata and map atomicity. `tests/psm_scoring.rs` reproduces the source PEPTIDE 11/33-peak HyperScore and Morpheus results and the squared-m/z ppm/Da examples. It independently checks two-pass reuse, charge mismatch behavior, all terminal-series detail, cross-link names, strict/inclusive boundaries, ordinal collapse, error conventions and invalid inputs. Existing comparison tests also run after the shared matching change.

Pinned source paths are `ANALYSIS/ID/{Scores,HyperScore,MorpheusScore}.cpp`, `ANALYSIS/ID/IDScoreSwitcherAlgorithm.h` and `DATASTRUCTURES/MatchedIterator.h`, with corresponding class tests. All are recorded in [source-inventory.json](source-inventory.json). Morpheus's source credits inspiration from C. Wenger's MIT-licensed C# implementation; this Rust translation is derived from the OpenMS BSD-3-Clause implementation. No third-party code or C++ build was modified or executed.

These modules do not perform database searching, infer score confidence, or implement other PSM scorers such as PScore, AScore, cross-link scoring families or consensus-ID algorithms. The full core-port goal remains ongoing; callers should use the [capability mapping](PORTING_STATUS.md) to distinguish implemented operations from remaining algorithms.

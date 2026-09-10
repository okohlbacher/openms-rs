# Averagine isotope cluster removal

`processing::deisotoping::AveragineDeisotoper` ports the Poisson/KL method from [Deisotoper.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/PROCESSING/DEISOTOPING/Deisotoper.cpp), pinned at `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. It uses the existing native `CoarseIsotopePatternGenerator::approximate_intensities` helper, not a separate isotope approximation. The [simple decreasing-intensity method](DEISOTOPING_SUPPORT.md) remains available as `Deisotoper`.

## API and defaults

`deisotope(&MSSpectrum) -> Result<DeisotopingResult>` returns a spectrum and accepted `IsotopeCluster` memberships. `SpectrumFilter::filter_spectrum` replaces its input only after success; its inherited experiment method also commits atomically. Metadata, precursors, scan settings and aligned auxiliary arrays are retained through every selection and ordering change.

| Option | Default |
| --- | --- |
| `tolerance` | 10 ppm; finite, nonnegative, at most 100 ppm or 0.1 Da |
| `min_charge`, `max_charge` | 1, 3; positive `u8` values |
| `min_isotope_peaks`, `max_isotope_peaks` | 2, 10 |
| `top_n` | `Some(5000)`; strongest peaks retained **before** clustering |
| `keep_only_deisotoped` | false; retain unassigned peaks after preprocessing |
| `make_single_charged` | true |
| `add_up_intensity` | false |
| `allow_shared_isotopes` | true, matching the source |
| `annotate_charge`, `annotate_isotope_peak_count`, `annotate_features` | false |
| `max_work` | 10,000,000 charged work units |

There is no decreasing-intensity switch or adjustable KL threshold in this API. The source's size-dependent probability checks determine accepted extensions.

## Preprocessing and original indices

Input peaks must have sorted finite nonnegative m/z values and nonnegative finite intensities. Invalid numerical values and misaligned arrays are errors before mutation. No MS-level restriction is imposed; the upstream reference spectrum is MS1.

First, retain intensities greater than or equal to **0.05**, comparing each stored f32 value as f64. This is the source `ThresholdMower` default and removes small positive peaks as well as zeros. It is applied even when `top_n` is disabled. Next, retain the `top_n` highest-intensity peaks when necessary and restore increasing m/z order. `None` disables top-N selection; `Some(0)` is invalid rather than acting as an alias for disabled. Equal-intensity boundary ties retain input order in the native implementation.

Top-N selection is preprocessing, not a cap on the final number of monoisotopic peaks. It can remove a weak isotope before cluster discovery. Existing integer, float and string arrays follow both selections.

Returned `IsotopeCluster::peak_indices` always refer to the **original input spectrum**, before thresholding or top-N selection. Clusters are listed by ascending seed m/z before optional charge conversion, independently of the final spectrum order. A valid initially empty spectrum is returned unchanged without new arrays; a nonempty spectrum emptied by thresholding produces valid empty requested annotation arrays.

## Detection and probability checks

Process unassigned seeds in ascending m/z order. Try every requested charge in ascending order. The expected isotope position is

`seed_mz + isotope_index * C13C12_MASSDIFF_U / charge`.

The nearest observed peak must be inside an inclusive tolerance window. Ppm tolerance is computed once from the seed m/z and used for every extension. Nearest-distance ties select the lower m/z. Members must advance strictly in input order; a ladder cannot reuse its seed or another member even when a tolerance window is wider than the isotope spacing.

When exactly one precursor has a known charge and positive neutral mass, the source precursor constraint applies. Its operation order is retained: `precursor_mz * precursor_charge - PROTON_MASS_U * precursor_charge`, and the candidate comparison uses `seed_mz * charge - PROTON_MASS_U * charge`. A candidate above precursor mass plus seed tolerance is skipped. Zero precursor charge, nonpositive precursor mass, no precursor or multiple precursors disable this restriction.

After finding a first extension, generate a truncated Poisson model with neutral mass `charge * (seed_mz - PROTON_MASS_U)` and lambda `mass / 1800`. This expression deliberately differs from the precursor constraint's floating-point operation order. Nonpositive model masses leave the seed unassigned; low-m/z noise does not abort detection of later valid clusters. Nonfinite mass arithmetic is an error.

The shared helper follows the source recurrence for ordinary finite inputs: start at weight 1, multiply each next weight by `lambda / isotope_index`, accumulate the sum and normalize the complete requested vector. It retains a logarithmic fallback when direct weights or their sum overflow. See [isotope support](ISOTOPE_SUPPORT.md).

For each proposed extension, normalize the observed and model prefixes separately, including that proposed peak. Compute natural-log Kullback–Leibler divergence `sum(P * ln(P / Q))`. Observed/model values and each term are f64, while the accumulator is f32: each update adds the f64 term to the promoted current accumulator and then rounds back to f32. Casting each term to f32 first changes real acceptance boundaries and is not equivalent.

| Prospective number of peaks | f32 KL threshold |
| --- | --- |
| 2 | 0.05 |
| 3 | 0.1 |
| 4 | 0.2 |
| 5 | 0.4 |
| 6 or more | 0.6 |

Reject only when finite KL is greater than its threshold, retaining source equality behavior. Missing peaks or failed model checks stop extension; the previously accepted prefix can still qualify if it reaches the minimum length. Positive observed signal against zero or invalid model-prefix probability rejects the candidate, so NaN cannot silently pass a comparison.

Choose the longest qualifying charge hypothesis for each seed; higher charge breaks equal-length ties. All requested charges are considered, unlike the simple method's first-success rule. Counts and summed intensities come from the **selected** ladder, rather than a later nonwinning trial.

## Shared isotopes and output

The source permits multiple accepted ladders to share heavy-isotope extensions, while skipping seeds already assigned to any accepted ladder. `allow_shared_isotopes: true` preserves this behavior by default. A shared original index appears in each relevant returned cluster.

With `add_up_intensity`, original observed member intensities are summed in f64 and checked before conversion to f32. A shared isotope contributes to **each** accepted cluster's sum, so summed output intensity can exceed the input's total. This is source behavior, not an intensity-conservation claim. Set `allow_shared_isotopes: false` for disjoint accepted memberships and sums. The existing simple `Deisotoper` also exposes this flag, with its earlier native default of false preserved.

Keep accepted monoisotopic peaks and, unless restricted, unassigned preprocessed peaks. Optional single-charge conversion uses `mz * charge - (charge - 1) * PROTON_MASS_U`, followed by sorting and the same permutation of all annotation arrays.

Optional integer arrays are `charge`, `iso_peak_count` and the native convenience `feature_number`. Accepted monoisotopic peaks receive the selected charge, ladder length and discovery ID; retained unassigned peaks receive 0, 1 and -1. Existing names in any array type cause a checked error instead of creating ambiguous duplicate annotations. `feature_number` is a native convenience corresponding to the simple method's array, not an extra source KL option.

## Checked differences and limits

The implementation preserves source sharing, scoring, candidate ranking and preprocessing. Deliberate checked differences are:

- Input and experiment mutations are transactional on all returned errors.
- Peak indices must advance within a ladder; a wide window cannot repeatedly select one observed peak.
- Isotope counts describe the winning cluster. C++ updates a seed's count during every trial and can retain a nonwinning or rejected trial's count.
- Nonpositive model masses remain unassigned, and nonfinite/zero model-prefix probabilities cannot produce NaN acceptance.
- Array-name collisions, invalid shapes, negative peak values, unsupported tolerances, overflow and exceeded resource limits are errors.
- Equal-intensity top-N ties use stable input ordering.

The work budget counts original input visits, preprocessed seed visits, charge hypotheses, isotope lookups, generated Poisson vector entries, and **every individual KL term**, including rejected extensions. Each loop accumulation is checked. Standard-library sorting and binary-search comparisons are not separately counted; input and discovery sizes remain bounded by the charged operations. At most 1,000,000 isotope entries are permitted by the shared helper. All isotope limits must be at least two, maximum must be at least minimum, and the work limit must be positive.

Discovery retains only the best completed candidate and the current trial, rather than allocating a padded matrix for every charge. Storage is linear in input size, the isotope vector and the total returned cluster memberships; preprocessing temporarily clones the input and retains selected arrays. Shared memberships can appear in multiple clusters and are bounded by the work budget. Repeated prefix KL evaluation can be quadratic in a trial's isotope count; the charged term budget prevents an oversized setting from running without a bound.

The method scores a C13-spaced Poisson envelope. It does not infer a molecular formula or peptide identity. Negative-ion processing and specialized nucleotide models remain outside these APIs.

## Validation and provenance

The [independent real-spectrum fixtures](../tests/data/averagine_deisotoping_provenance.json) retain all 5,407 upstream input peaks and all 104 expected output peaks, including exact binary floating values extracted from the source mzML fixtures. The default native output matches **all 104 peaks bit for bit**, with independently inferred original seeds and charge assignments. Metadata and aligned arrays are also checked.

The reference includes a genuine shared extension: source seed 4550 at charge 1 and seed 4564 at charge 3 both use original isotope peak 4578. The returned default clusters preserve both memberships. Disjoint mode produces the independently justified 103-peak subset, removing only the second overlapping seed. Source golden data remains unchanged.

`tests/averagine_deisotoping_review.rs` independently tests adjacent-f32 KL acceptance/rejection pairs for cluster sizes 2–7, correct mixed-precision accumulation, longest-ladder/highest-charge ranking, selected counts and sums, inclusive tolerance endpoints, nearest ties, precursor arithmetic and shared extensions. `tests/averagine_deisotoping.rs` covers preprocessing, original index mapping, aligned arrays, both sharing policies on both algorithms, low-m/z noise, invalid data, overflow and transactional spectrum/experiment failures. No C++ project was built or executed.

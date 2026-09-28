# PScore and AScore support

Source headers, core pin `bc9cc12`:

- `src/openms/include/OpenMS/ANALYSIS/ID/PScore.h` (`PScore.cpp`)
- `src/openms/include/OpenMS/ANALYSIS/ID/AScore.h` (`AScore.cpp`)

Rust: `src/analysis/id/pscore.rs`, `src/analysis/id/ascore.rs`.
Tests: `tests/pscore_ascore.rs`. Provenance: `tests/data/pscore_ascore_provenance.json`.

The two headers are one group because they depend on each other: `AScore` includes
`PScore.h`, and `PScore` calls `AScore::computeCumulativeScore_` through a
`friend` declaration. The TOPP consumer is `PhosphoScoring` (`AScore`); `PScore` has
none.

## API mapping

### `PScore` (a struct of static members)

| Source | Rust |
|---|---|
| `PScore()` / `~PScore()` (implicit) | `PScore` unit struct: `Default`, `Clone`, `Copy` |
| `static std::vector<Size> calculateIntensityRankInMZWindow(mz, intensities, mz_window)` | `PScore::calculate_intensity_rank_in_mz_window(&[f64], &[f64], f64) -> Result<Vec<usize>>` |
| `static std::vector<std::vector<Size>> calculateRankMap(const PeakMap&, double mz_window = 100)` | `PScore::calculate_rank_map(&MSExperiment, f64)`; default `DEFAULT_MZ_WINDOW` |
| `static std::map<Size, PeakSpectrum> calculatePeakLevelSpectra(spec, ranks, min_level = 1, max_level = 9)` | `PScore::calculate_peak_level_spectra(&MSSpectrum, &[usize], usize, usize) -> Result<BTreeMap<usize, MSSpectrum>>`; defaults `DEFAULT_MIN_LEVEL`, `DEFAULT_MAX_LEVEL` |
| `static double computePScore(tol, ppm, peak_level_spectra, const std::vector<PeakSpectrum>&, mz_window = 100)` | `PScore::compute_pscore_best_of(..)` |
| `static double computePScore(tol, ppm, peak_level_spectra, const PeakSpectrum&, mz_window = 100)` | `PScore::compute_pscore(..)` |
| `static double massCorrectionTerm(double)` | `PScore::mass_correction_term(f64) -> f64` |
| `static double cleavageCorrectionTerm(Size, bool)` | `PScore::cleavage_correction_term(usize, bool) -> f64` |
| `static double modificationCorrectionTerm(Size)` | `PScore::modification_correction_term(usize) -> f64` |

### `ProbablePhosphoSites`

| Source | Rust |
|---|---|
| `Size first, second, seq_1, seq_2, peak_depth, AScore` | `ProbablePhosphoSites { first, second, seq_1, seq_2, peak_depth, a_score }` (`usize`, `Default` zeroes, as value-initialisation) |

### `AScore` (a `DefaultParamHandler` named `AScore`)

| Source | Rust |
|---|---|
| `AScore()`, `defaults_` | `AScore::default()`, `AScore::defaults() -> Result<Param>`, `ASCORE_NAME` |
| `setParameters` + `updateMembers_` | `AScore::from_param`, `AScore::from_param_with_warnings` |
| `getParameters` | `AScore::to_param` |
| `~AScore()` | drop |
| `PeptideHit compute(const PeptideHit&, PeakSpectrum&)` | `AScore::compute(&self, &PeptideHit, &mut MSSpectrum) -> Result<PeptideHit>` |
| protected `int compareMZ_(double, double) const` | `AScore::compare_mz -> Ordering` |
| protected template `getSpectrumDifference_(first1, last1, first2, last2, out)` | `AScore::spectrum_difference(&[Peak1D], &[Peak1D]) -> Vec<Peak1D>` |
| protected `computeSiteDeterminingIons_(th_spectra, candidates, out)` | `AScore::compute_site_determining_ions -> Result<[MSSpectrum; 2]>` |
| protected `getSites_(const std::string&)` | `AScore::get_sites(&str) -> Vec<usize>` |
| protected `computePermutations_(sites, Int)` | `AScore::compute_permutations(&[usize], usize) -> Result<Vec<Vec<usize>>>` |
| protected `numberOfMatchedIons_(th, window, depth)` | `AScore::number_of_matched_ions -> Result<usize>` |
| protected `peptideScore_(scores)` | `AScore::peptide_score(&[f64]) -> Result<f64>` |
| protected `determineHighestScoringPermutations_(scores, sites&, permutations, ranking&)` | `AScore::determine_highest_scoring_permutations(scores, permutations, ranking) -> Result<Vec<ProbablePhosphoSites>>` |
| protected `computeBaseProbability_(double)` | `AScore::compute_base_probability(&self, f64) -> f64` |
| protected `computeCumulativeScore_(N, n, p)` | `AScore::compute_cumulative_score(usize, usize, f64) -> Result<f64>` (associated: it reads no member) |
| protected `numberOfPhosphoEvents_(const std::string&)` | `AScore::number_of_phospho_events(&str) -> usize` |
| protected `removePhosphositesFromSequence_(const std::string&)` | `AScore::remove_phosphosites_from_sequence(&str) -> Result<AASequence>` |
| protected `createTheoreticalSpectra_(permutations, seq)` | `AScore::create_theoretical_spectra -> Result<Vec<MSSpectrum>>` |
| protected `peakPickingPerWindowsInSpectrum_(PeakSpectrum&)` | `AScore::peak_picking_per_windows_in_spectrum(&MSSpectrum) -> Result<Vec<MSSpectrum>>` |
| protected `calculatePermutationPeptideScores_(th_spectra&, windows)` | `AScore::calculate_permutation_peptide_scores(&self, &[MSSpectrum], &[MSSpectrum], f64)` |
| protected `rankWeightedPermutationPeptideScores_(scores)` -> `std::multimap<double, Size>` | `AScore::rank_weighted_permutation_peptide_scores -> Result<Vec<(f64, usize)>>` (ascending, stable) |
| protected `updateMembers_()` | `from_param` (private `from_complete_param`) |
| member `fragment_mass_tolerance_` | field `fragment_mass_tolerance` (`fragment_mass_tolerance`, default 0.05, min 0) |
| member `fragment_tolerance_ppm_` | field `fragment_tolerance_ppm` (`fragment_mass_unit`, `Da`/`ppm`) |
| member `max_peptide_length_` | field `max_peptide_length` (`max_peptide_length`, default 40, advanced) |
| member `max_permutations_` | field `max_permutations` (`max_num_perm`, default 16384, advanced) |
| member `unambiguous_score_` | field `unambiguous_score: i32` (`unambiguous_score`, integer default 1000, advanced) |
| member `base_match_probability_` (mutable state set in `compute`) | no field: an explicit argument of `calculate_permutation_peptide_scores`, computed per call by `compute_base_probability` |
| `friend struct PScore` | `compute_cumulative_score` is public |

The protected members are public here because they are the surface the upstream
class test drives (through an `AScoreTest` subclass) and each is one documented
step of `compute`.

## Preserved source conventions

- Ranks are zero based and count strictly more intense neighbours within
  `+/- mz_window / 2`; the scan stops at the first neighbour outside the window.
- A peak level exists in the map only once a peak was added to it.
- The single-spectrum PScore uses the `MatchedIterator` walk (one match per
  theoretical peak, inclusive tolerance in `f32`, ppm from the theoretical m/z)
  and `p = (level + 1) / mz_window`; the vector overload uses `findNearest`, a
  strict `<` in `f64`, and `p = level / mz_window` (so a level-0 match scores
  `+inf`). Both are kept, each documented at its function.
- `AScore::compute` resets the score to `-1` on every aborted path, returns
  `unambiguous_score` when every S/T/Y is phosphorylated, writes
  `search_engine_sequence` (lowercase, as the code writes; the header's `@note`
  says `Search_engine_sequence`), `AScore_pep_score` and `AScore_<n>`, and scores
  the hit with the lowest `AScore_<n>`.
- Equal best and runner-up weighted scores give `AScore_<n> = 0` for every site.
- Ranking ties keep insertion order and the best is the last entry, as
  `std::multimap::rbegin()`.
- Placements are enumerated in lexicographic order of site positions, which is
  the order the source recursion produces.
- Windows are 100 Th wide from `floor(first/100)*100`, include their upper bound,
  keep at most 10 peaks by stable descending intensity.
- The peptide score weights and the `/ 7` are the source's, summed in its order.
- The spectrum difference skips every further ion within tolerance on both sides
  of a match, exactly as the source loop does.
- `compute` sorts the input spectrum in place when it is unsorted, as the source
  does through its non-const reference.
- Other modifications than phospho are kept in the returned sequence (the code's
  behaviour; the header's brief claims they are dropped).

## Native differences

Every one of these replaces undefined behaviour or an unchecked precondition in
the source; none changes a result the source defines.

- **Refusals instead of undefined reads:** mismatched rank input lengths; more
  ranks than peaks; an empty level spectrum in the vector PScore (the source's
  `findNearest` throws `Exception::Precondition`); out-of-range placement indices;
  a ranking or permutation list that is empty or has no qualifying runner-up (the
  source walks past `rend()`); a score row that is not 10 long (debug-only
  precondition upstream); negative m/z in the window split (the source casts the
  bound to `Size`).
- **Non-finite input** (m/z, intensity, window, tolerance) is refused with
  `Error::InvalidValue`; the source's comparisons misbehave silently on NaN.
- **Unsorted spectra** passed to the `MatchedIterator`-based functions are
  refused with `Error::UnsortedData`; the source documents sorted input and
  does not check.
- **More phosphorylations than S/T/Y** (a phospho on another residue): the
  source recurses into `sites[0]` of an empty vector; `compute` returns the reset
  hit (score `-1`) and `compute_permutations` returns no placement.
- **An empty placement** in `create_theoretical_spectra` adds no modification
  (the source reads `permutations[i][0]` out of bounds).
- **Resource ceilings**, checked before allocation: `MAX_RANK_PEAKS`
  (rank input), `MAX_ASCORE_PERMUTATIONS` (placement count, computed as a
  binomial coefficient before enumerating, so `max_num_perm = 0` cannot ask for
  billions of placements), `MAX_ASCORE_WINDOWS` (m/z span).
- **Binomial errors** (`p > 1`, or `n > N` when an ion is matched in two adjacent
  windows) return `Error::InvalidValue`; the source's `Math::binomial_cdf_complement`
  throws `std::invalid_argument` for the same inputs.
- **Non-finite scores** cannot be stored as metadata: a site score of `inf - inf`
  (both cumulative probabilities underflowing to 0) returns an error where the
  source stores NaN.
- `AScore::compare_mz` returns `Ordering` instead of `-1/0/1`.
- `unambiguous_score` is an `i32` because the source parameter is an integer;
  the returned score is its `f64` value, as the source's `double` member.

## Source defect candidates

1. `PScore.cpp` defines `massCorrectionTerm`, `cleavageCorrectionTerm` and
   `modificationCorrectionTerm` as free functions in `namespace OpenMS`, not as
   `PScore::` members; the declared static members have no definition, so a call
   to `PScore::massCorrectionTerm` fails to link. The class test's sections for
   them are empty ("Not tested"), which is why it was never noticed.
2. The two `computePScore` overloads use different match probabilities
   (`level / mz_window` vs `(level + 1) / mz_window`); with the level-0 spectrum
   the vector overload scores any match as `+inf`.
3. `AScore.h`'s `@note` on `compute` names the metadata key
   `Search_engine_sequence`; the code writes `search_engine_sequence`. The brief
   says other modifications are dropped; the code keeps them.
4. `ProbablePhosphoSites::AScore` is never assigned.
5. `computePermutations_` with more events than sites, and `compute` for a hit
   with a phospho on a non-S/T/Y residue, index an empty vector.

## Checked boundaries and evidence

Tier 3 (source review): every literal of `PScore_test.cpp` and `AScore_test.cpp`
at `bc9cc12`, including the six `AScore::compute` results on the four upstream
DTA inputs (`Ascore_test_input1..4.dta`, copied unchanged and hashed in the
provenance manifest) in Da and ppm, the unambiguous-score case, and the PScore
values `83.867454`, `154.682242` and `240` on `PScore_test.dta` and `IFSQVGK`.
All pass at the class test's `TEST_REAL_SIMILAR` tolerance (1e-5).

Tier 4 (native): the refusals above, `calculate_rank_map` against the
per-spectrum ranks, the vector PScore against the closed form `20 N` for a
self-match at level 1, the correction tables, parameter defaults and round trip,
the spectrum-difference skip rule, and the aborted-compute paths.

No tier 1 oracle was run for this group: the class-test literals cover every
public function with a defined result except `calculateRankMap`, the vector
`computePScore` and the three correction terms, which the upstream test leaves
empty.

The source is serial (no OpenMP); so is the port.

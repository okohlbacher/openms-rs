# PosteriorErrorProbabilityModel: native header equivalent

[`src/math/posterior_error_probability.rs`](../src/math/posterior_error_probability.rs)
provides a native equivalent for the **numerical** surface of core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`
`MATH/STATISTICS/PosteriorErrorProbabilityModel.h` and
`source/MATH/STATISTICS/PosteriorErrorProbabilityModel.cpp`. **Five public
members** of the header are not ported here — `extractAndTransformScores`
(`:72`), `updateScores` (`:95`), `initPlots` (`:216`),
`plotTargetDecoyEstimation` (`:228`) and `tryGnuplot` (`:237`), the last public
member before `private:` at `:239` — together with the two private helpers of
the first pair, `transformScore_` (`:247`) and `getScore_` (`:253`). Each is
named below with its reason.

Tests: [`tests/posterior_error_probability.rs`](../tests/posterior_error_probability.rs).
Manifest: [`tests/data/math_kde_provenance.json`](../tests/data/math_kde_provenance.json).
The component fitters are [`DISTRIBUTION_FITTERS_SUPPORT.md`](DISTRIBUTION_FITTERS_SUPPORT.md).

## API mapping

Every public member of the header appears here.

| Source member | Native representation |
| --- | --- |
| `PosteriorErrorProbabilityModel()` | `PosteriorErrorProbabilityModel::new`, `Default` |
| `~PosteriorErrorProbabilityModel()` | `Drop` is derived; nothing to release |
| `DefaultParamHandler` base, `defaults_` for `out_plot`, `number_of_bins`, `incorrectly_assigned`, `max_nr_iterations`, `neg_log_delta`, `outlier_handling` | `PepParameters` (`number_of_bins`, `incorrectly_assigned`, `max_nr_iterations`, `neg_log_delta`), `IncorrectComponent`, `OutlierHandling`, `PosteriorErrorProbabilityModel::with_parameters`, `parameters`, `set_parameters`. `out_plot` is **not ported**: nothing here writes files |
| `bool fit(std::vector<double>&, const std::string& outlier_handling)` | `fit(&mut [f64], OutlierHandling) -> Result<bool>` |
| `bool fit(std::vector<double>&, std::vector<double>& probabilities, const std::string&)` | `fit_with_probabilities(&mut [f64], OutlierHandling) -> Result<Option<Vec<f64>>>`; the out-parameter becomes the return value and the `false` return becomes `None` |
| `bool fitGumbelGauss(std::vector<double>&, const std::string&)` | `fit_gumbel_gauss(&mut [f64], OutlierHandling) -> Result<bool>` |
| `void fillDensities(const std::vector<double>&, std::vector<double>&, std::vector<double>&)` | `fill_densities(&[f64]) -> Result<(Vec<f64>, Vec<f64>)>`, returning `(incorrect, correct)` |
| `void fillLogDensities(...)` | `fill_log_densities(&[f64]) -> Result<(Vec<f64>, Vec<f64>)>` |
| `void fillLogDensitiesGumbel(...)` | `fill_log_densities_gumbel(&[f64]) -> Result<(Vec<f64>, Vec<f64>)>` |
| `double computeLogLikelihood(const std::vector<double>&, const std::vector<double>&) const` | `compute_log_likelihood(&[f64], &[f64]) -> f64` |
| `double computeLLAndIncorrectPosteriorsFromLogDensities(const std::vector<double>&, const std::vector<double>&, std::vector<double>&) const` | `compute_ll_and_incorrect_posteriors_from_log_densities(&[f64], &[f64]) -> (f64, Vec<f64>)` |
| `std::pair<double,double> pos_neg_mean_weighted_posteriors(const std::vector<double>&, const std::vector<double>&)` | associated `pos_neg_mean_weighted_posteriors(&[f64], &[f64]) -> (f64, f64)`; it reads no member, so it is not a method |
| `std::pair<double,double> pos_neg_sigma_weighted_posteriors(const std::vector<double>&, const std::vector<double>&, const std::pair<double,double>&)` | associated `pos_neg_sigma_weighted_posteriors(&[f64], &[f64], (f64, f64)) -> (f64, f64)` |
| `GaussFitter::GaussFitResult getCorrectlyAssignedFitResult() const` | `correctly_assigned_fit_result()` |
| `GaussFitter::GaussFitResult getIncorrectlyAssignedFitResult() const` | `incorrectly_assigned_fit_result()` |
| `GumbelMaxLikelihoodFitter::GumbelDistributionFitResult getIncorrectlyAssignedGumbelFitResult() const` | `incorrectly_assigned_gumbel_fit_result()` |
| `double getNegativePrior() const` | `negative_prior()` |
| `static double getGumbel_(double, const GaussFitter::GaussFitResult&)` | `PosteriorErrorProbabilityModel::gumbel_density(f64, GaussFitResult) -> Result<f64>`; public in the source despite the trailing underscore |
| `double computeProbability(double) const` | `compute_probability(f64) -> Result<f64>` |
| `double getSmallestScore() const` | `smallest_score()` |
| `const std::string getGumbelGnuplotFormula(const GaussFitter::GaussFitResult&) const` | associated `gumbel_gnuplot_formula(GaussFitResult) -> String` |
| `const std::string getGaussGnuplotFormula(const GaussFitter::GaussFitResult&) const` | associated `gauss_gnuplot_formula(GaussFitResult) -> String` |
| `const std::string getBothGnuplotFormula(const GaussFitter::GaussFitResult&, const GaussFitter::GaussFitResult&) const` | `both_gnuplot_formula(GaussFitResult, GaussFitResult) -> String`; reads the negative-formula selector, so it stays a method |
| `getNegativeGnuplotFormula_` / `getPositiveGnuplotFormula_` member function pointers | `NegativeFormula` plus `negative_formula()`; the positive pointer is always the Gaussian and needs no field |
| private `void processOutliers_(std::vector<double>&, const std::string&) const` | private `process_outliers`, driven by `OutlierHandling` |
| `static std::map<std::string, std::vector<std::vector<double>>> extractAndTransformScores(...)` | **not ported**: needs `ProteinIdentification`, `PeptideIdentificationList`, `PeptideHit` from the crate's `metadata` and `identification` modules |
| `static void updateScores(...)` | **not ported**: same dependencies |
| private `static double transformScore_(...)`, `static double getScore_(...)` | **not ported**: helpers of the two above |
| `TextFile initPlots(std::vector<double>&)` | **not ported**: returns `FORMAT/TextFile` and writes a `.txt` beside it |
| `void plotTargetDecoyEstimation(std::vector<double>&, std::vector<double>&)` | **not ported**: writes two files and shells out |
| `void tryGnuplot(const std::string&)` | **not ported**: calls `system("gnuplot ...")` |
| Not in source but added | `MAX_ITEMS`, `SCORE_SHIFT`, `PepParameters`, `IncorrectComponent`, `OutlierHandling`, `NegativeFormula` |

The five unported public members, and the two private helpers that serve them,
are all identification-, file- or process-facing.
`tools/check_module_cycles.py` freezes the crate's cross-module edge set and
`math` reaches no other top-level module today; `extractAndTransformScores` and
`updateScores` would need `metadata` and `identification`, `initPlots` and
`plotTargetDecoyEstimation` would need `format`, and `tryGnuplot` would need a
process launcher. They belong in a layer above `math` — the score-transformation
table is search-engine knowledge, not numerics — and are recorded as a deferral
rather than forced through the graph. The gnuplot *formula* builders are pure
string formatting and are ported, so the layer above has nothing numeric left to
reimplement.

## Preserved source conventions

- **The score transform.** Both fits sort the caller's scores ascending in
  place, record the smallest, and work on `score + |smallest| + 0.001`. The
  header's `@note` that "the vector is sorted from smallest to biggest value" is
  therefore load-bearing, and `computeProbability` re-applies the same shift, so
  it must be called with an untransformed score.
- **The initial parameters, literally.** For `fit`: the incorrect component's
  location is `mean(x[0 .. ceil(n/2)]) + x[0]` — the mean of the lower half
  *plus the minimum* — its width is `sd(x, that location)` about that location
  rather than about the sample mean, and its amplitude is
  `1 / sqrt(2 pi sigma^2)`. The correct component starts at
  `mean(x[s ..]) + x[s]` with `s = min(n-1, ceil(0.7 n))`, shares the incorrect
  component's width, and takes the same amplitude form. `negative_prior_` starts
  at `0.7`. `fitGumbelGauss` uses the same expressions for the Gumbel's `a` and
  `b`. These are not textbook starting values; a mixture fit lands on a
  different local optimum if they change, which is why they are transcribed
  rather than rationalised.
- **The convergence test and the cap.** The loop is a `do`/`while`, so its body
  always runs at least once. It stops when
  `new_loglikelihood - loglikelihood < 10^(-neg_log_delta)` or when
  `itns >= max_nr_iterations`, both tested *before* `itns` is incremented, so
  the body runs at most `max_nr_iterations + 1` times. A decrease in the
  likelihood stops the loop with `good_fit = false`; a `NaN` step returns
  `false` immediately.
- **The "impossible standard deviations" branch reports success.** Both fits
  `break` out of the loop there with `good_fit` still `true`, so a fit that
  aborted mid-iteration returns `true`. Reproduced, and stated at the item.
- **Two logarithm bases.** `computeLogLikelihood` accumulates base-ten
  logarithms; the EM loop uses
  `computeLLAndIncorrectPosteriorsFromLogDensities`, which accumulates natural
  ones, and it is the natural-log increase that the `10^-delta` threshold is
  compared against. Nothing in the source calls `computeLogLikelihood`.
- **The log-sum-exp responsibilities.** The larger of the two log
  responsibilities is subtracted before exponentiating, the posteriors are the
  normalised incorrect share, and the likelihood accumulates
  `max_log_resp + ln(sum)`.
- **`fit` fits two Gaussians whatever `incorrectly_assigned` says.** The
  parameter selects only which formula the plots use and which density
  `computeProbability` evaluates. The source's own `TODO: incorrect is currently
  filled with gauss as fitting gumble is not supported` says so at
  `fillDensities` and `fillLogDensities`.
- **`computeProbability` evaluates a Gumbel in every branch**, reading
  `incorrectly_assigned_fit_param_`'s `x0` and `sigma` as a Gumbel location and
  scale even when the Gaussian was selected. A second `TODO` calls this
  "confusing at best". Below the incorrect location the incorrect density is
  replaced by its peak and above the correct location the correct density is, so
  the probability stays monotone.
- **`processOutliers_`'s quartiles are the median-of-halves
  `quantile1st`/`quantile3rd` with `sorted = true`**, not the interpolating
  `Math::quantile`. The IQR rule drops outside `Q1 - 3 IQR .. Q3 + 3 IQR`; the
  clamp rule rewrites them to the nearest inside value; the percentile rule
  drops everything `<= x[floor(n / 100) + 1]` or `>= x[floor(n * 99.9 / 100)]`
  — note `99.9`, not the `99` the parameter description promises, and inclusive
  comparisons, so equal values at either end go too.
- **The gnuplot expressions**, character for character, including the space in
  `exp(( ` and the number format — see below.

## Native differences

- **`DefaultParamHandler` becomes a plain struct.** `PepParameters` carries the
  four numeric parameters with the source's defaults; `IncorrectComponent` and
  `OutlierHandling` carry the two enumerated ones, with `parse`/`as_str` for the
  source's strings. The source's `outlier_handling` dispatch tests `!= "none"`
  and then two named alternatives, so any unrecognised string silently selects
  `ignore_extreme_percentiles`; `OutlierHandling::parse` rejects it instead,
  which is what `Param::setValidStrings` was there to do.
- **Numbers in the gnuplot expressions are formatted as `printf("%g")` with six
  significant digits**, which is what an unconfigured `std::ostream` writes.
  Rust's `{}` prints the shortest round-tripping form, which would be a
  different byte sequence for the same fit, so `format_g` implements the `%g`
  rule — including the `%e` switch outside `[-4, 6)` and the trailing-zero trim
  — and is pinned by its own test.
- **`compute_probability` refuses on an unfitted model**
  (`Error::MissingInformation`) and after `fit_gumbel_gauss`. See the defect
  below: the source computes `max_incorrectly_` there from a Gaussian parameter
  set `fitGumbelGauss` never writes, so a fresh model's "peak density" is
  `-0.368`. The port declines to report a probability derived from it rather
  than returning a wrong number. The upstream class test's own probability loop
  is commented out in that section for the same reason.
- **Degenerate inputs are refused** rather than indexed: an outlier rule that
  removes every score, a `set_iqr_to_closest_valid` where nothing is inside the
  fence (the source decrements an `upper_bound` unconditionally and dereferences
  a `lower_bound`), and `ignore_extreme_percentiles` on a **single** score,
  which is the only length for which the source's unchecked
  `x_scores[floor(n / 100) + 1]` reads past the end — at `n = 1` that index is
  `1` on a one-element vector, and for every larger `n` both indices are in
  range. A sample of fewer than a hundred scores is *not* refused for being
  small: it is accepted whenever the rule leaves something behind, and refused
  by the "removed every score" guard above when it does not. A non-finite score
  is refused before the sort.
- **`gumbel_density` refuses a non-positive or non-finite scale**, where the
  source divides by it twice.
- **A single score is refused by `fit`**, because `sd` with one degree of
  freedom removed is undefined there; the source divides by zero and carries
  `NaN` into the EM loop.
- **`fit` and `fit_gumbel_gauss` return `Result<bool>`.** The `bool` is the
  source's return value; the `Err` is a refusal the source does not make.
- **Bounded work.** `MAX_ITEMS = 50,000,000` scores.
- **Serial**, as the source is: no `#pragma omp` in the translation unit.

## Source defects found

1. **`fitGumbelGauss` computes its peak from an unfitted parameter set.** Its
   last numeric statement is
   `max_incorrectly_ = getGumbel_(incorrectly_assigned_fit_param_.x0, incorrectly_assigned_fit_param_)`,
   reading the *Gaussian* parameters, which that function never writes. On a
   freshly constructed model they are still `(-1, -1, -1)` and the expression
   evaluates to `-0.3679` — a negative "peak density" that `computeProbability`
   then divides by. The upstream class test only reaches this function on an
   object a previous `fit` had already populated, which is why it has never
   surfaced. The port computes the peak only when a previous `fit` left a usable
   Gaussian set, and marks the model as not ready for `compute_probability`
   either way.
2. **The "Aborting fit" branch returns success.** Both EM loops `break` on an
   impossible standard deviation with `good_fit` still `true`, so a caller that
   checks the return value cannot tell an aborted fit from a converged one.
   Reproduced.
3. **`ignore_extreme_percentiles` uses `99.9 / 100` where its own parameter
   description says "99th and 1st percentile"**, and indexes
   `x_scores[x_scores.size() / 100 + 1]` without a bounds check. Working the
   two index expressions out, `floor(n / 100) + 1` and `floor(n * 99.9 / 100)`,
   the read is out of range for exactly one length, `n = 1`; the missing bounds
   check is still a defect, but a narrower one than "small samples".

All three are reported in the work package's C++ issue list.

## Checked boundaries and evidence

| Boundary | Where |
| --- | --- |
| `MAX_ITEMS = 50,000,000` scores | `prepare_scores` |
| every score finite | `prepare_scores`, before the sort |
| outlier handling must leave at least one score | `prepare_scores` |
| `set_iqr_to_closest_valid` needs a value inside the fence | `process_outliers` |
| `ignore_extreme_percentiles` needs both indices in range | `process_outliers` |
| Gumbel scale finite and positive, evaluation point finite | `gumbel_density` |
| mixture density non-zero, model fitted | `compute_probability` |
| every Gaussian evaluation validated | `GaussFitResult::eval`, `log_eval_no_normalize` |

Evidence:

- **Tier 3** for the class test's literals. The upstream tolerances are
  deliberately loose — `TOLERANCE_ABSOLUTE(0.5)` on the component parameters —
  because they are what an EM fit lands on rather than a closed form. What makes
  them worth reproducing is that they pin the *local optimum*: a port with a
  different initialisation, convergence test or iteration cap converges
  elsewhere, and at this tolerance that shows.
  - `GaussMix_2_1D.csv`, 2,000 draws from a mixture of `N(1.5, 0.5)` and
    `N(3.5, 1.0)`, reproduces `x0 = 3.5`, `sigma = 1.0`, `x0 = 1.5`,
    `sigma = 0.5` and a prior of `0.5`.
  - The 20 hand-written scores reproduce `4.62`, `0.87`, `1.06`, `0.77` and a
    prior of `0.546`, and `getSmallestScore` of `-0.39`.
  - `GumbelGaussMix_2_1D.csv`, 2,000 draws, reproduces the Gaussian at
    `8 - smallest` with `sigma = 3.5`, the Gumbel at `2 - smallest` with
    `b = 0.6`, and a prior of `0.6`.
  - The gnuplot substring assertions `(1/0.90` and `exp(( 1.47` are the
    sharpest of these: they pin the fitted scale and location of the small-sample
    fit to three significant digits *and* the number format at the same time.
- **Tier 4** for the values derived in the test and marked there: the
  weighted-moment helpers against hand-computed sums; the Gumbel density at its
  location being exactly `exp(-1) / sigma`; `fill_densities` and
  `fill_log_densities` shown equal to the component evaluations bit for bit;
  `compute_log_likelihood` shown equal to its own base-ten formula;
  probabilities lying in `[0, 1]` and never rising along the sorted scores; the
  ten `%g` formatting cases; and the outlier rules changing the fit.

No retained C++ output and no oracle driver exists for this header, so no tier 1
or tier 2 claim is made.

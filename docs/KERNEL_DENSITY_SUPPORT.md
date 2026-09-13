# KernelDensityEstimation: native header equivalent

[`src/math/kernel_density.rs`](../src/math/kernel_density.rs) provides a native
equivalent for the complete public surface of core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`
`MATH/STATISTICS/KernelDensityEstimation.h` and
`source/MATH/STATISTICS/KernelDensityEstimation.cpp`.

Tests: [`tests/kernel_density.rs`](../tests/kernel_density.rs). Manifest:
[`tests/data/math_kde_provenance.json`](../tests/data/math_kde_provenance.json).
The transforms it calls are [`FFT_SUPPORT.md`](FFT_SUPPORT.md).

The estimator follows Silverman (1982), Algorithm AS 176, as reimplemented in
`statsmodels.nonparametric` and, through that, in PyProphet — whose conventions
the source adopts where they differ from the paper.

## API mapping

Every public member of the header appears here.

| Source member | Native representation |
| --- | --- |
| `struct KernelDensityEstimation` | the module; the type holds no state and is not reproduced |
| `static double bwNrd0(const std::vector<double>& x)` | `bw_nrd0(&[f64]) -> Result<f64>` |
| `static std::vector<double> linBin(const std::vector<double>&, double xmin, double xmax, std::size_t nbins, const std::vector<double>* weights)` | `lin_bin(&[f64], f64, f64, usize, Option<&[f64]>) -> Result<Vec<f64>>` |
| `static std::vector<double> linBin(const std::vector<double>&, double, double, std::size_t)` (uniform-weight overload) | the same function with `None`; a null pointer and an absent slice are the same thing |
| `static std::vector<double> forRt(const std::vector<double>& X, std::size_t M = 0)` | `for_rt(&[f64], usize) -> Result<Vec<f64>>`; `0` still means "use the input length" |
| `static std::vector<double> revRt(const std::vector<double>& Xp, std::size_t M = 0)` | `rev_rt(&[f64], usize) -> Result<Vec<f64>>` |
| `static std::vector<double> silvermanKernelFFT(double bw, std::size_t M, double RANGE)` | `silverman_kernel_fft(f64, usize, f64) -> Result<Vec<f64>>` |
| `static std::pair<std::vector<double>, std::vector<double>> gridKdeFFT(const std::vector<double>&, double bw, std::size_t gridsize = 512, double cut = 3.0)` | `grid_kde_fft(&[f64], f64, usize, f64) -> Result<GridKde>`; the pair becomes the named `GridKde { density, grid }`, because `first`/`second` do not say which is which. `first` is `density`, `second` is `grid` |
| `static std::vector<double> kdeFFTEval(const std::vector<double>&, double bw, std::size_t gridsize = 512, double cut = 3.0)` | `kde_fft_eval(&[f64], f64, usize, f64) -> Result<Vec<f64>>` |
| the default arguments `gridsize = 512`, `cut = 3.0` | `DEFAULT_GRIDSIZE`, `DEFAULT_CUT` |
| file-static `fast_linbin(const std::vector<double>&, double a, double b, std::size_t M)` | private `fast_lin_bin`; not public in the source and not public here |
| `OpenMS::CubicSpline2d`, used by `kdeFFTEval` | private `NaturalCubicSpline`, a term-by-term copy of `CubicSpline2d::init_`/`eval` — see differences |
| Not in source but added | `MAX_ITEMS`, `MAX_GRID`, `MIN_GRID`, `NRD0_IQR_DIVISOR`, `GridKde` |

## Preserved source conventions

- **`bwNrd0`'s exact formula and accumulation order.** Non-finite values are
  dropped first; fewer than two survivors return `0.0`. The mean is a running
  sum in input order divided by `n`, the variance a second pass divided by
  `n - 1`. The quartiles are `Math::quantile` at `0.25` and `0.75` — Hyndman-Fan
  type 7, the `numpy.percentile` convention — **not** the median-of-halves
  `quantile1st`/`quantile3rd` that live next to it in `StatisticFunctions.h`.
  The result is `0.9 * min(sd, IQR / 1.34) * n^(-1/5)`.
- **PyProphet's fallback chain** when `min(sd, IQR/1.34)` is not positive: the
  standard deviation, then `|x[0]|` — and since the sample has already been
  sorted, PyProphet's `x[0]` is the *minimum*, not the first input — then the
  literal `1.0`. Ten zeros therefore give `0.9 * 10^(-1/5)`, which the class
  test pins.
- **`linBin` is a histogram.** Each value lands whole in
  `floor((v - xmin) / width)` with `width = (xmax - xmin) / nbins`, and `xmax`
  folds into the last bin. See the defects below; the class tests of both
  `KernelDensityEstimation` and `MultipleTesting` pin the histogram behaviour,
  so it is what callers depend on.
- **A weight slice of the wrong length is ignored, not rejected**, reproducing
  `weights != nullptr && weights->size() == x.size()`.
- **`fast_linbin` really is linear binning**, splitting each value between its
  two neighbouring grid points of `linspace(a, b, M)` — `M - 1` intervals, not
  `M` — and putting anything at or past the last interval whole into the last
  point.
- **`forRt` is unscaled and `revRt` is the plain inverse.** `statsmodels`'
  `forrt` divides by `M` and its `revrt` multiplies by `M`; the source does
  neither. Both pairs compose to the identity and the kernel product between
  them is unaffected, which is why the estimates agree.
- **The Munro packing**: `M` reals holding
  `[Re Y_0 .. Re Y_{M/2}, Im Y_1 .. Im Y_{M/2-1}]`.
- **Silverman's analytic kernel**, `exp(-2 (pi bw j / range)^2)` over the binning
  correction `1 - (j pi / M)^2 / 3`, with the correction replaced by
  `numeric_limits<double>::min()` when it is not positive — a substitution that
  drives that bin's response to the largest finite value instead of flipping its
  sign. The packed second half mirrors bins `1 .. M/2 - 1`.
- **`gridKdeFFT`'s grid.** `M = bit_ceil(max(gridsize, n, 512))`, so a caller's
  `gridsize` is a lower bound and never the answer; `a = min(x) - cut*bw`,
  `b = max(x) + cut*bw`, and an empty sample centres on zero. The binned counts
  are divided by `spacing * n`, the transform multiplied by the kernel, and the
  result **renormalised so that `sum(density) * spacing == 1`** — the source's
  final rescale, which pulls the estimate slightly away from the unrenormalised
  convolution.
- **`kdeFFTEval` interpolates with a natural cubic spline and clamps negatives to
  zero**, in that order.

## Native differences

- **`CubicSpline2d` is reimplemented privately.** The crate ports that class as
  `CubicSpline2d` in `src/processing/spline/cubic.rs`, but
  `tools/check_module_cycles.py` freezes the cross-module edge set and `math`
  reaches no other top-level module. Adding `math -> processing` would widen the
  graph, which the gate forbids. `NaturalCubicSpline` is therefore a term-by-term
  copy of `CubicSpline2d::init_` and `CubicSpline2d::eval`, private to the
  module. `tests/kernel_density.rs` asserts the two produce **bit-identical**
  values on the actual KDE grid, for a three-point sample and for a 64-point one,
  so the duplication is checked rather than assumed. Removing the duplication
  needs a restructuring of the module graph and is recorded as a deferral.
- **Errors instead of silent degeneracies.** `silverman_kernel_fft` refuses a
  zero or non-finite `range` and a non-finite `bw`, where the source divides
  unchecked and returns an all-`NaN` kernel for `bw == range == 0`.
  `grid_kde_fft` refuses a non-finite sample value — the source hands it to
  `std::minmax_element` — and a grid of zero width. `for_rt`/`rev_rt` refuse a
  non-power-of-two length, where the source silently transforms a different
  number of points (see [`FFT_SUPPORT.md`](FFT_SUPPORT.md)).
- **`kde_fft_eval` reports a query outside the grid** as `Error::InvalidValue`,
  which is the source's `Exception::IllegalArgument` from `CubicSpline2d::eval`.
- **Bounded work.** `MAX_ITEMS = 50,000,000` sample values and
  `MAX_GRID = 2^22` grid points, both checked before anything is allocated. The
  source has no ceiling and, above `2^16`, computes the wrong transform.
- **Serial**, as the source is: `KernelDensityEstimation.cpp` carries no
  `#pragma omp`, so there is no parallel behaviour to match and no performance
  gap to state.

## Source defects found

1. **`linBin` does not do linear binning.** The header says "Each data point is
   allocated to its two nearest grid points proportionally based on distance"
   and cites Scott (1992); the implementation is a plain histogram. The file
   contains real linear binning in the file-static `fast_linbin`, which is what
   `gridKdeFFT` uses, so the public function and the internal one disagree
   despite the shared name. The class test comment at
   `KernelDensityEstimation_test.cpp:94` ("first bin gets count from 0.0 and
   part of 0.5") describes the documented behaviour rather than the tested one.
2. **`revRt`'s documented scaling is wrong.** The header says the output is
   "scaled by multiplying by M (applied by the `evergreen::real_ifft`
   function)". `real_ifft` is the exact inverse: `DIF::real_ifft1d_packed`
   scales by `1 / (N/2)`, which is the inverse normalisation of the half-length
   complex transform. If the header were right, the class test's own
   `forRt`/`revRt` round-trip section would fail by a factor of `M`. The
   sentence describes `statsmodels`' `revrt`, not this code.
3. **`forRt`'s documented rounding does not happen.** "`M` Length of FFT (if 0,
   uses `X.size()` and rounds up to next power of 2)" — the code uses
   `X.size()` and rounds nothing; see [`FFT_SUPPORT.md`](FFT_SUPPORT.md).

All three are reported in the work package's C++ issue list.

## Checked boundaries and evidence

| Boundary | Where |
| --- | --- |
| `MAX_ITEMS = 50,000,000` sample values | `check_items`, at every entry point |
| `MAX_GRID = 2^22` bins or grid points | `lin_bin`, `silverman_kernel_fft`, `grid_kde_fft` |
| transform length a power of two | `for_rt`, `rev_rt`, via `math::fft::check_length` |
| `nbins > 0`, `xmax > xmin`, finite bounds | `lin_bin` |
| finite `bw`, finite non-zero `range` | `silverman_kernel_fft` |
| finite sample, finite `bw` and `cut`, positive grid width and spacing | `grid_kde_fft` |
| spline abscissae strictly increasing and finite, coefficients finite, query inside the grid | `NaturalCubicSpline` |

Evidence:

- **Tier 3** for the class test's literals, section by section, with the source
  line above each Rust test.
- **Tier 3** for the two fixtures the upstream test ships and this port copies
  into `tests/data/`: `kde_reference_data.csv`, produced by `statsmodels`, and
  the `scipy.stats.gaussian_kde` vector transcribed in the class test's last
  section. Both are third-party reference output retained by the upstream test
  rather than output of the C++ under test, so they are independent of the
  implementation but still reach this port through a transcription step. They
  are the strongest evidence in this group: `tests/kernel_density.rs` reproduces
  every `statsmodels` density in the fixture to the class test's 5% relative or
  0.005 absolute tolerance, and every scipy density to 10% relative or 0.01
  absolute.
- **Tier 4** for the values derived here and marked in the test: `bw_nrd0` of
  the symmetric nine-point ladder recomputed in closed form as
  `0.9 sqrt(15/8) 9^(-1/5)`; the five-point ladder as `0.9 (2/1.34) 5^(-1/5)`;
  the exact packed transform of `[1,2,3,4]`; the Silverman kernel evaluated from
  its closed form at every bin; `bw_nrd0` with `NaN`/`inf` shown equal to
  `bw_nrd0` of the surviving values; the estimate's integral being exactly one;
  and the bit-identical agreement with the crate's ported `CubicSpline2d`.

No retained C++ output and no oracle driver exists for this header, so no tier 1
or tier 2 claim is made.

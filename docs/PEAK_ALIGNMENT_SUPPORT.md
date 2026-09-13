# PeakAlignment support

Port of `src/openms/include/OpenMS/COMPARISON/PeakAlignment.h` and
`src/openms/source/COMPARISON/PeakAlignment.cpp` at Core SDK revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(header sha256 `88b7056ccddfc491436078fd0a0123ab2b736e8389dd3c3a6661944fe877696e`,
`.cpp` sha256 `a621dee886979eba7771ddbe52cc8c450d108e872971baef395542e04979a053`).

Rust: [`comparison::PeakAlignment`](../src/comparison.rs).
Tests: [`tests/comparison_scorers.rs`](../tests/comparison_scorers.rs).
Provenance: [`tests/data/comparison_scorers_provenance.json`](../tests/data/comparison_scorers_provenance.json).

## API mapping

Every public and private member of the header appears here.

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| `class PeakAlignment : public PeakSpectrumCompareFunctor` | `pub struct PeakAlignment` implementing [`PeakSpectrumCompareFunctor`](PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md) | composition of a `DefaultParamHandler` instead of inheritance |
| `PeakAlignment()` | `PeakAlignment::new() -> Result<Self>` | reproduces the four `defaults_.setValue` calls and `defaultsToParam_()` (`PeakAlignment.cpp:24-28`) **and the absent `setName`**: the handler keeps the base's `"PeakSpectrumCompareFunctor"` |
| `PeakAlignment(const PeakAlignment& source)` | `Clone` | the C++ copy constructor is `= default` |
| `~PeakAlignment() override` | drop glue | C++ destructor is `= default` |
| `PeakAlignment& operator=(const PeakAlignment& source)` | `Clone::clone_from`, or assignment of a clone | the C++ body is the self-assignment guard plus the base assignment |
| `double operator()(const PeakSpectrum& spec1, const PeakSpectrum& spec2) const override` | `PeakSpectrumCompareFunctor::score` | returns `Result<f64>`; the two shortcuts, the matrix, the border scan and the normalisation are all reproduced |
| `double operator()(const PeakSpectrum& spec) const override` | `PeakSpectrumCompareFunctor::self_score`, the trait default | the C++ override is `return operator()(spec, spec)` (`PeakAlignment.cpp:44-47`) |
| `std::vector<std::pair<Size, Size>> getAlignmentTraceback(const PeakSpectrum&, const PeakSpectrum&) const` | `alignment_traceback(&self, &MSSpectrum, &MSSpectrum) -> Result<Vec<(usize, usize)>>` | returns `Result`; the matrix fill is shared with `score` rather than duplicated |
| `private: double peakPairScore_(double&, double&, double&, double&, const double&) const` | private `peak_pair_score(f64, f64, f64, f64, f64) -> Result<f64>` | the source's non-const reference parameters are taken by value, which they always were in effect; returns `Result` because of the `sqrt` |

The header's `@param[in] spec1 First spectrum given in a binned representation`
is wrong: `operator()` takes `PeakSpectrum`, not `BinnedSpectrum`, and nothing in
this class bins anything. The Rust documentation says what the arguments are.

`@htmlinclude OpenMS_PeakAlignment.parameters` documents four entries,
registered with the source's own defaults and description strings:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `epsilon` | `0.2` (float) | absolute mass error, **and** the gap cost |
| `normalized` | `1` (integer) | **never read**; see below |
| `heuristic_level` | `0` (integer) | number of strongest peaks the shortcut considers; `0` disables it |
| `precursor_mass_tolerance` | `3.0` (float) | precursor distance beyond which the score is `0` |

## Preserved source conventions

- **Two shortcuts, in this order** (`PeakAlignment.cpp:54-96`): precursors more
  than `precursor_mass_tolerance` apart score `0` - with a missing precursor
  reading as m/z `0`, as everywhere in this hierarchy - and then, when
  `heuristic_level` is nonzero, two spectra whose `heuristic_level` most intense
  peaks share no m/z within `epsilon` also score `0`.
- **The gap cost is `epsilon`.** `PeakAlignment.cpp:99` reads the same parameter
  a second time under a `//TODO gapcost dependence on distance ?`.
- **The matrix.** `(n + 1) * (m + 1)` cells, first row and column pre-charged
  with `-gap * i`. A cell whose two peaks are more than `epsilon` apart can only
  come from a gap. Ties in the direction matrix keep the zero
  `Matrix<Size>(rows, cols)` value-initialised them to, which the traceback reads
  as "from the left"; that is the case the source marks
  `// TODO the cases where all or two values are equal`, and it is reproduced
  because the reported alignment is observable.
- **The reported score is not the corner cell.** It is the largest value of the
  last row or the last column, so a suffix of either spectrum may go unaligned -
  the recurrence is global but the report is semi-global.
- **`numeric_limits<double>::min()` seeds the border scan**
  (`PeakAlignment.cpp:186`). That is the smallest *positive* normal double, not
  the most negative one, so a matrix whose entire last row and column are
  negative reports `2.2e-308` rather than its real maximum. The port reproduces
  the seed, because every published score from this class carries it.
- **The position term is not the Gaussian it looks like.**
  `exp(-(fabs(pos1 - pos2)) / 2 * sigma * sigma)` (`PeakAlignment.cpp:401`)
  parses as `exp(((-|Δ|) / 2) * sigma * sigma)`: the distance enters linearly and
  sigma **multiplies** the exponent instead of dividing it, so a wider peak
  distribution makes distant peaks score *less*. The intended `exp(-Δ²/(2σ²))`
  would need different parentheses. The expression is transcribed exactly.
- **The sigma passes.** `mid` and `var` both run over every pair in row-major
  order and both divide by the exact integer pair count widened to `double`.
- `normalized` is registered and never read; `PeakAlignment.cpp:217` normalises
  unconditionally. The parameter is registered here so the surface matches, and
  is read no more than the source reads it. A test asserts that clearing it
  changes nothing.

## Native differences

- **One matrix fill serves both entry points.** The source duplicates the loop
  in `operator()` and `getAlignmentTraceback`, and the two copies already differ:
  `operator()` guards a zero variance, `getAlignmentTraceback` does not. Sharing
  the fill removes the risk that a future edit makes the score and the reported
  alignment disagree; the sigma difference is preserved, because it is the
  source's behaviour, and is what each entry point passes in.
- **The zero-variance guard is unusable upstream, and refusing is the only
  honest option.** When every pairwise distance is equal - a single peak on each
  side, say - `operator()` substitutes `numeric_limits<double>::min()` for sigma.
  The position term is then `1 / (DBL_MIN * sqrt(2 pi))`, about `1.8e307`, so the
  **product** of the two self-alignment scores overflows to infinity for every
  nonzero `f32` intensity, down to the smallest subnormal, and the quotient
  silently becomes `0` - complete dissimilarity for a spectrum compared with
  itself. Here the overflow is `Err(Error::InvalidValue)`.
- `getAlignmentTraceback` has no variance guard at all, so a zero variance
  divides by a zero sigma and fills the matrix with infinities. That is refused.
- An empty spectrum that reaches the alignment divides by a zero pair count and
  carries NaN through the whole matrix, ending at `+inf`. That is refused. The
  class test's "empty spectra should return zero" is unaffected: its empty
  spectrum has no precursor, so the precursor shortcut fires first and the port
  returns the same `0`.
- Unsorted peaks are `Err(Error::UnsortedData)`; the source checks nothing.
- Negative intensities are refused, because `peakPairScore_` takes the square
  root of the intensity product.
- A negative `heuristic_level` is refused. The source casts it to `UInt`, turning
  `-1` into 4294967295 and quietly making the heuristic scan the whole spectrum.
- `MAX_ALIGNMENT_MATRIX_CELLS = 4_000_000` bounds the score matrix, checked
  before either buffer is allocated. The source allocates that matrix plus, in
  the traceback, an `n * m` direction matrix, with no ceiling: two 5000-peak
  spectra ask it for 200 MB.
- The heuristic sorts stably by intensity. `std::sort` leaves the order among
  equal intensities unspecified, so which peaks land in a tied top-`level` set is
  undefined upstream; the stable sort makes the selection deterministic without
  changing it whenever the intensities at the cut are distinct.

## Checked boundaries and evidence

| Boundary | Behaviour |
| --- | --- |
| precursors beyond `precursor_mass_tolerance` | `Ok(0.0)`, before anything is allocated |
| `heuristic_level > 0` with no shared strong peak | `Ok(0.0)` |
| empty spectrum with the precursor shortcut armed | `Ok(0.0)`, as the class test observes |
| empty spectrum reaching the alignment | `Err(Error::InvalidValue)`; source returns `inf` |
| zero peak-distance variance, `score` | `Err(Error::InvalidValue)`; source overflows to `0` |
| zero peak-distance variance, `alignment_traceback` | `Err(Error::InvalidValue)`; source divides by zero |
| zero self-alignment product | `Err(Error::InvalidValue)` |
| unsorted peaks | `Err(Error::UnsortedData)` |
| negative intensity | `Err(Error::InvalidValue)` |
| negative `heuristic_level` | `Err(Error::InvalidValue)`; source casts to a huge `UInt` |
| matrix cells above `MAX_ALIGNMENT_MATRIX_CELLS` | `Err(Error::InvalidValue)` before allocation |
| all-negative border | `2.2e-308`, the source's `DBL_MIN` seed, reproduced |

OpenMP: no `#pragma omp` in this header or its `.cpp`; the port is serial, as the
source is. The two `O(n*m)` sigma passes and the matrix anti-diagonals would
parallelise, but the source does neither and a parallel float sum would have to
reproduce this accumulation order exactly, so there is nothing to take.

### Class-test sections

`src/tests/class_tests/openms/source/PeakAlignment_test.cpp`
(sha256 `ec9f506ac9fb62a54ae0121cda0fa057e07e132ac09ab2b02a5775fedfe2743e`),
**seven sections, all ported**. The fixture is
`PILISSequenceDB_DFPIANGER_1.dta` (127 peaks), retained byte-identical as
[`tests/data/comparison_dfpianger.dta`](../tests/data/comparison_dfpianger.dta).

| Section | Rust test | Asserted value | Tier |
| --- | --- | --- | --- |
| `PeakAlignment()` | `peak_alignment_construction_copy_and_assignment` | `name() == "PeakSpectrumCompareFunctor"` - the absent `setName` - and the four defaults with their source values and types | 4 |
| `~PeakAlignment()` | same | construction and drop | 4 |
| `(PeakAlignment(const PeakAlignment&))` | same | equal name and parameters, the upstream pair | 3 |
| `(PeakAlignment& operator=(const PeakAlignment&))` | same | equal name and parameters | 3 |
| `(double operator()(const PeakSpectrum&, const PeakSpectrum&) const)` | `peak_alignment_reproduces_the_upstream_golden_score` | `0.997477` to 1e-6 and `0.9974770186204426` to 1e-12 against the fixture with its last peak dropped; both empty-spectrum cases return `0.0` | 3 |
| `(double operator()(const PeakSpectrum&) const)` | same | `self_score(s1) == 1.0` exactly | 3, plus 4 |
| `(vector<pair<Size,Size>> getAlignmentTraceback(const PeakSpectrum&, const PeakSpectrum&) const)` | `peak_alignment_traceback_is_the_diagonal_for_a_self_alignment` | 127 pairs, each of them `(i, i)`, plus strict ascent on both indices | 3, plus 4 |

The transcribed literal `0.997477` is tier 3; it was reproduced to full precision
as `0.9974770186204426` by an independent Python model of the matrix fill, the
`DBL_MIN`-seeded border scan and the mis-parenthesised position term, written
before any Rust, which is what makes the 1e-12 tolerance defensible. The self
score `1` and the 127 traceback pairs are derived rather than transcribed: the
diagonal of a self-alignment is exactly the self-alignment sum, so the quotient
is exactly one - asserted with `assert_eq`, not a tolerance - and every diagonal
step of a self-alignment strictly beats both gaps, so the traceback is the
identity. The upstream section only checks the first index of each pair; this
test checks both.

Two further tests cover behaviour the class test never exercises:
`peak_alignment_shortcuts_and_unread_normalized_flag` pins both shortcuts and the
dead `normalized` flag, and `peak_alignment_bounds_its_matrix_and_rejects_unsorted_peaks`
pins the cell ceiling and the sortedness check.

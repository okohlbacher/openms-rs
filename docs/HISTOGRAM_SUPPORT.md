# Histogram: native header equivalent

[`src/math/histogram.rs`](../src/math/histogram.rs) provides a native equivalent
for the complete public surface of core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4` `MATH/STATISTICS/Histogram.h`. The
accompanying `src/openms/source/MATH/STATISTICS/Histogram.cpp` is an empty
translation unit; the header is a header-only template.

Tests: [`tests/histogram.rs`](../tests/histogram.rs). Manifest:
[`tests/data/math_statistics_provenance.json`](../tests/data/math_statistics_provenance.json).

## The binning convention

Everything the 23 class-test sections are about follows from four rules:

1. The value range `[min, max]` is **closed at both ends**.
2. Bins are **half-open**, `[left, right)`, laid out from `min` with width
   `bin_size`.
3. The bin count is `ceil((max - min) / bin_size)`, except that `max == min`
   yields exactly one bin. The last bin therefore usually extends past `max`.
4. `max` itself is routed into the **last bin** by an explicit special case,
   whatever rule 2 would otherwise say.

Rule 4 is why `right_border_of_bin` reports the next representable value above
`max` for the last bin instead of `min + (index + 1) * bin_size`: the last bin
really does contain its upper endpoint, so its open right border has to be one
step beyond it.

## API mapping

Every public and protected member of the header appears here.

| Source member | Native representation |
| --- | --- |
| `template <typename ValueType = UInt, typename BinSizeType = double> class Histogram` | non-generic `Histogram`; both types are `f64` (see differences) |
| `typedef ConstIterator` | `std::slice::Iter<'_, f64>`, returned by `iter()` |
| `Histogram()` | `Histogram::new()`, `Default` |
| `Histogram(const Histogram&)` | `Clone` |
| `Histogram(min, max, bin_size)` | `Histogram::with_bounds(f64, f64, f64) -> Result<Self>` |
| `template <DataIterator> Histogram(begin, end, min, max, bin_size)` | `Histogram::from_values(&[f64], f64, f64, f64) -> Result<Self>` |
| `virtual ~Histogram()` | ordinary drop; the port is not polymorphic |
| `BinSizeType minBound() const` | `min_bound(&self) -> f64` |
| `BinSizeType maxBound() const` | `max_bound(&self) -> f64` |
| `ValueType maxValue() const` | `max_value(&self) -> Option<f64>` |
| `ValueType minValue() const` | `min_value(&self) -> Option<f64>` |
| `BinSizeType binSize() const` | `bin_size(&self) -> f64` |
| `Size size() const` | `len(&self) -> usize`, plus `is_empty()` |
| `ValueType operator[](Size) const` | `bin(&self, usize) -> Result<f64>` |
| `BinSizeType centerOfBin(Size) const` | `center_of_bin(&self, usize) -> Result<f64>` |
| `BinSizeType rightBorderOfBin(Size) const` | `right_border_of_bin(&self, usize) -> Result<f64>` |
| `BinSizeType leftBorderOfBin(Size) const` | `left_border_of_bin(&self, usize) -> Result<f64>` |
| `ValueType binValue(BinSizeType) const` | `bin_value(&self, f64) -> Result<f64>` |
| `Size inc(BinSizeType val, ValueType increment = 1)` | split: `inc(&mut self, f64)` for the default increment, `inc_by(&mut self, f64, f64)` otherwise; both return `Result<usize>` |
| `Size incUntil(BinSizeType, bool inclusive, ValueType increment = 1)` | `inc_until(&mut self, f64, bool, f64) -> Result<usize>` |
| `Size incFrom(BinSizeType, bool inclusive, ValueType increment = 1)` | `inc_from(&mut self, f64, bool, f64) -> Result<usize>` |
| `static void getCumulativeHistogram(begin, end, complement, inclusive, Histogram&)` | `add_cumulative(&mut self, &[f64], bool, bool) -> Result<()>`, a method on the histogram it was mutating anyway |
| `void reset(min, max, bin_size)` | `reset(&mut self, f64, f64, f64) -> Result<()>` |
| `bool operator==(const Histogram&) const` | `PartialEq` (derived: bounds, bin width and every bin) |
| `bool operator!=(const Histogram&) const` | `PartialEq` |
| `Histogram& operator=(const Histogram&)` | ordinary assignment / `Clone` |
| `ConstIterator begin() const` | `iter()`, and `IntoIterator for &Histogram` |
| `ConstIterator end() const` | the same iterator's end |
| `void applyLogTransformation(BinSizeType multiplier)` | `apply_log_transformation(&mut self, f64) -> Result<()>` |
| `Size valueToBin(BinSizeType) const` | `value_to_bin(&self, f64) -> Result<usize>` |
| protected `min_`, `max_`, `bin_size_`, `bins_` | private fields behind the accessors above, plus `bins() -> &[f64]` |
| protected `initBins_()` | the private `bin_count` helper used by `with_bounds` and `reset` |
| free `operator<<(ostream&, const Histogram&)` | `impl Display`, one line per bin: centre, tab, count |
| Not in source but added | `MAX_BINS` resource ceiling; `bins()` slice accessor |

## Preserved source conventions

- The four binning rules above, including the one-bin `max == min` case and the
  `nextafter` right border of the last bin.
- `centerOfBin(i) = min + (i + 0.5) * bin_size` and
  `leftBorderOfBin(i) = min + i * bin_size`, computed in that form rather than
  by accumulating.
- `inc`'s increment is neither required to be positive nor integral. The class
  test increments by `45.0` and `250.3`, and `applyLogTransformation` replaces
  counts with logarithms, so "count" is a loose word on both sides.
- `incUntil` raises `[0, index)` and, when inclusive, `index`; `incFrom` raises
  `(index, len)` and, when inclusive, `index`. Both return the index of the
  value's own bin.
- `getCumulativeHistogram` **adds to** whatever the histogram already holds and
  has no rollback: a value outside the range aborts the loop with the earlier
  values already applied. That is reproduced, and said at the item.
- `applyLogTransformation` computes `multiplier * ln(x + 1)`.

## Native differences

| Difference | Reason |
| --- | --- |
| Not generic over `ValueType`/`BinSizeType`; both are `f64` | Every instantiation the library uses stays expressible, including the class test's `Histogram<float, float>` with fractional bin contents. It also removes a real defect: under the default `ValueType = UInt`, `applyLogTransformation` casts each transformed value back to an unsigned integer and truncates it. |
| `Result` instead of exceptions | `Exception::OutOfRange` on a non-positive bin width maps to `Error::InvalidValue`; `Exception::OutOfRange` on a value outside `[min, max]` and `Exception::IndexOverflow` on a bin index both map to `Error::InvalidRange`. |
| `bin(index)` instead of `Index` | The source's `operator[]` throws; a Rust `Index` implementation would have to panic, which this crate does not do on untrusted input. |
| `min_value`/`max_value` return `Option` | The source returns `*std::max_element(bins_.begin(), bins_.end())`, which dereferences the end iterator of a default-constructed histogram's empty bin vector. |
| `with_bounds` and `reset` reject an inverted range, a non-finite bound or width, and a bin count above `MAX_BINS` | The source checks only the bin width. An inverted range makes it take the ceiling of a negative quotient and convert it to an unsigned `Size`; a narrow width lets the bin count grow until the allocation fails. |
| A bin count that underflows to zero is clamped to one | The source would produce a histogram with a non-empty range and no bins, into which nothing can be counted. |
| `reset` builds the new bins before replacing anything | The source clears `bins_` and assigns the new bounds *before* throwing on a non-positive width, leaving an object whose bounds changed and whose bins are gone. |
| `value_to_bin` rejects NaN and an empty bin vector | Both of the source's bound comparisons are false for a NaN, which then reaches `floor(NaN)` and an unsigned conversion — undefined behaviour in C++. |
| `value_to_bin` clamps the computed index to the last bin | Only reachable when the division rounds up onto the boundary; it replaces an out-of-bounds index with the neighbouring bin the value belongs to. |
| `inc`, `inc_by`, `inc_until`, `inc_from` reject a non-finite increment | The source would store it, after which `min_value`/`max_value` are poisoned and `apply_log_transformation` is undefined. |
| `apply_log_transformation` rejects a bin at or below `-1` and commits atomically | The source evaluates `log` of a non-positive argument and stores the NaN in place, bin by bin. |
| `~Histogram()` is virtual in the source | Nothing derives from it in the SDK; the port is a plain value type. |
| `getCumulativeHistogram` is a method, not a static | It took the histogram by mutable reference already. |
| Serial only | The header carries no `#pragma omp`. |

## Checked boundaries and evidence

Resource boundaries: `MAX_BINS = 16_777_216` bins (128 MiB of `f64` counts),
checked from the caller-controlled `ceil((max - min) / bin_size)` before the
vector is allocated.

Numeric boundaries: `min == max`; `max < min`; non-finite bounds, widths,
values and increments; a value exactly on `min`, exactly on `max`, exactly on an
interior bin edge, and just outside either end; every bin index at and beyond
`len()`; the last bin's `nextafter` right border; a default-constructed
histogram, on which every lookup errors.

Evidence, per `docs/DIFFERENTIAL_VALIDATION.md`:

- **Tier 3 (source review, transcribed class-test literals)** for all 23
  `START_SECTION`s of `Histogram_test.cpp`, in `tests/histogram.rs`. The class
  test mutates one shared `Histogram<float, float>` down the file; each test
  here rebuilds the state its section sees, and the shared `filled()` helper
  reproduces the increments of the `inc` section that four later sections read.
- **Tier 4 (independently derived)**: the bin counts `10`, `5` and `6` from
  `ceil((max - min) / bin_size)`; the log-transformed last bin as exactly
  `ln(10001)` rather than the class test's rounded `9.21044`; the `incUntil`
  and `incFrom` bin patterns, which have no class-test section at all; and the
  `nextafter` right border asserted as `f64::from_bits(bits + 1)` with the
  accompanying `> 5.0` check.
- **Tier 4 (Rust-only)** for every guard in the differences table.

The class test's `rightBorderOfBin` expectation is
`std::nextafter(5.0f, 6.0f)` under `BinSizeType = float`; this port's bin
coordinates are `f64`, so the step is the `f64` ulp. The *behaviour* — "one
representable step above `maxBound()`" — is identical, and the literal cannot
be, which is why the test asserts the step rather than the number.

No tier 1 or tier 2 evidence exists for this group.

## Candidate C++ issues

Reported to the integrating agent rather than written into
`OpenMS_CPP_ISSUES.md`, which that agent owns:

1. `Histogram::minValue`/`maxValue` (lines 112-123) dereference
   `std::min_element`/`std::max_element` over `bins_` without checking that it
   is non-empty; the default constructor leaves it empty.
2. `Histogram::reset` (lines 274-296) clears `bins_` and assigns the new bounds
   before throwing `Exception::OutOfRange` on a non-positive `bin_size`, so a
   failed reset leaves a half-updated object.
3. `Histogram::initBins_` and `reset` compute
   `Size(ceil((max_ - min_) / bin_size_))` with no check that `max_ >= min_`;
   converting a negative `double` to an unsigned `Size` is undefined behaviour,
   and there is no ceiling on the resulting allocation.
4. `Histogram::valueToBin` (lines 371-387) compares `val < min_ || val > max_`,
   both of which are false for a NaN, and then converts `floor(NaN)` to `Size`.
5. `Histogram::applyLogTransformation` (lines 358-365) casts the result back to
   `ValueType`, so the header's default `UInt` instantiation truncates the
   transformed value to an integer, and it evaluates `log` of `*it + 1` without
   checking the argument's sign.

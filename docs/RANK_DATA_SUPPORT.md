# RankData: native header equivalent

[`src/math/rank_data.rs`](../src/math/rank_data.rs) provides a native equivalent
for the complete public surface of core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4` `MATH/STATISTICS/RankData.h`. The
header is header-only — a struct of static templates plus three free functions,
with no accompanying translation unit.

Tests: [`tests/rank_data.rs`](../tests/rank_data.rs). Manifest:
[`tests/data/math_statistics_provenance.json`](../tests/data/math_statistics_provenance.json).

The source replicates `scipy.stats.rankdata(a, method=..., nan_policy=...)` with
one-dimensional semantics, and its class test pins whole vectors against SciPy
1.17.1. Ranks are one-based.

**This is not the crate's only ranking.**
[`statistic_functions::compute_rank`](STATISTIC_FUNCTIONS_SUPPORT.md), which
Spearman correlation is built on, detects ties with a *relative tolerance* of
`1e-7` and only ever averages them. `rankdata` compares by exact equality and
offers five tie rules. The two disagree on near-ties, both have callers, and
neither replaces the other; `tests/rank_data.rs` asserts the disagreement so it
cannot be "unified" by accident.

## API mapping

Every public member of the header appears here.

| Source member | Native representation |
| --- | --- |
| `struct RankData` | the module; there is no state to hold, so the namespace struct is not reproduced as a type |
| `enum class RankData::Method { Average, Min, Max, Dense, Ordinal }` | `RankMethod` with the same five variants; `Default` is `Average`, the source's default argument |
| `enum class RankData::NaNPolicy { Propagate, Omit, Raise }` | `NanPolicy` with the same three variants; `Default` is `Propagate` |
| `template <class T> static std::vector<double> rankdata(const std::vector<T>&, Method, NaNPolicy)` | `rankdata<T: Rankable>(&[T], RankMethod, NanPolicy) -> Result<Vec<f64>>` |
| `static std::vector<double> rankdata_double(const std::vector<double>&, Method, NaNPolicy)` | `rankdata_f64(&[f64], RankMethod, NanPolicy) -> Result<Vec<f64>>` |
| `static std::vector<double> rankdata_float(const std::vector<float>&, Method, NaNPolicy)` | `rankdata_f32(&[f32], RankMethod, NanPolicy) -> Result<Vec<f64>>` |
| `static std::vector<double> rankdata_int(const std::vector<int>&, Method, NaNPolicy)` | `rankdata_i32(&[i32], RankMethod, NanPolicy) -> Result<Vec<f64>>` |
| free `rankdata_double(std::vector<double> a, ...)` (by value) | the same `rankdata_f64`; see differences |
| free `rankdata_float(std::vector<float> a, ...)` (by value) | the same `rankdata_f32` |
| free `rankdata_int(std::vector<int> a, ...)` (by value) | the same `rankdata_i32` |
| the template's `if constexpr (std::is_floating_point<T>::value)` NaN branch | the `Rankable` trait's `is_nan_value`, which is `false` for `i32` |
| the template's `static_cast<D>` promotion to `double` | the `Rankable` trait's `as_f64` |
| Not in source but added | `MAX_ITEMS` resource ceiling; the public `Rankable` trait that carries the two operations the template's `if constexpr` performed |

The source's `rankdata_double`/`_float`/`_int` exist twice — once as static
members for Python bindings and once as by-value free functions — and both
spellings forward to the same template. The port has one function per element
type. Taking a `&[T]` makes the by-value copy unnecessary: the source's
by-value overloads copy because the static members take a const reference and
the free functions were written for call sites that already owned a temporary;
neither semantics is observable, because `rankdata` never modifies its input.

## Preserved source conventions

- **One-based ranks**, matching SciPy.
- **A stable sort by value**, which is what makes `Ordinal` break ties by
  original position and what makes every other method see a tied run as one
  contiguous block.
- **Ties by exact equality** of the values widened to `f64`. An `f32` input is
  compared after widening, so an `f32` tie and an `f64` tie agree.
- **The five tie rules**, each on the block spanning sorted positions
  `[lo, hi)` whose one-based ranks are `lo+1 ..= hi`: `Min` takes `lo + 1`,
  `Max` takes `hi`, `Average` takes `0.5 * (min + max)`, `Dense` takes the
  block's ordinal among distinct values, `Ordinal` takes the position itself.
- **`Propagate` returns an all-NaN vector of the input's length** as soon as any
  NaN is present, for every method, matching SciPy 1.17.1 and the source's
  early return.
- **`Omit` ranks only the non-NaN values and leaves the NaN positions NaN**,
  and an all-NaN input under `Omit` returns all NaN rather than an error.
- **An empty input returns an empty vector**, under every policy including
  `Raise`.
- **An integer input never triggers a NaN policy**, because the source's
  `if constexpr` compiles the detection away.

## Native differences

| Difference | Reason |
| --- | --- |
| `Result` instead of `std::invalid_argument` | `NanPolicy::Raise` on an input containing a NaN returns `Error::InvalidValue`. |
| One function per element type instead of a template plus three forwarders and three by-value overloads | `&[T]` removes the reason the by-value overloads existed; `Rankable` carries what `if constexpr` decided. |
| `Rankable` is public | A caller cannot add an instantiation the way it could instantiate the C++ template, but it can see exactly which types are supported and why. |
| The sort comparator is `partial_cmp`, not `f64::total_cmp` | The source compares with `<`, under which `-0.0` and `0.0` are equal. The IEEE-754 total order would put `-0.0` first and change `Ordinal`'s answer for a vector containing both zeros. No NaN can reach the comparator: `Propagate` and `Raise` have already returned and `Omit` has filtered them out. |
| `MAX_ITEMS` ceiling | The source allocates whatever it is handed. |
| Serial only | The header carries no `#pragma omp`. |

## Checked boundaries and evidence

Resource boundaries: `MAX_ITEMS = 50_000_000` values, checked before either of
the two owned buffers is allocated.

Numeric boundaries: empty input; all-NaN input under each policy; a NaN under
each of the three policies and each of the five methods; signed zeros, which tie
under every method and keep input order under `Ordinal`; infinities, which sort
to the ends and are not NaN, so no policy applies to them; near-ties at `1e-12`,
which are distinct here and tied in `compute_rank`.

Evidence, per `docs/DIFFERENTIAL_VALIDATION.md`:

- **Tier 3 (source review, transcribed class-test literals)** for all 11
  `START_SECTION`s of `RankData_test.cpp`, in `tests/rank_data.rs`. The
  cross-validation section's vectors were generated by the source's authors with
  SciPy 1.17.1 and are transcribed; they are not an execution performed here,
  and that is the whole reason they are tier 3 rather than tier 1.
- **Tier 4 (independently derived)**: that `Average` is exactly the midpoint of
  `Min` and `Max` at every position of the ten-element ties vector; that
  `Ordinal` is a permutation of `1..=n`; and the block arithmetic behind every
  hand-computed small vector.
- **Tier 4 (Rust-only)** for the ceiling, the signed-zero ordering, the
  infinity handling, and the documented disagreement with `compute_rank`.

No tier 1 or tier 2 evidence exists for this group.

## Candidate C++ issues

None found in this header. Its one surprise — that `Propagate` discards a whole
vector because of a single NaN — is SciPy's documented behaviour and is
reproduced deliberately.

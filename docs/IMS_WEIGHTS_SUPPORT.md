# IMS integer weights

`chemistry::IMSWeights` is the owned Rust equivalent of source
`OpenMS::ims::Weights` at SDK commit
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. It pairs original f64 masses with
integer weights computed at a chosen precision. It does not enumerate mass
compositions and does not require the strictly positive alphabet used by a solver.
No dependency or C++ runtime is involved.

```rust
use openms::chemistry::IMSWeights;

let mut weights = IMSWeights::from_masses(&[3.0, 5.0, 8.0], 0.1)?;
assert_eq!(weights.weights(), &[30, 50, 80]);
assert!(weights.divide_by_gcd()?);
assert_eq!(weights.weights(), &[3, 5, 8]);
assert_eq!(weights.precision(), Some(1.0));
assert_eq!(weights.parent_mass(&[2, 1, 0])?, 11.0);
# Ok::<(), openms::Error>(())
```

## Complete public operation mapping

| Source operation | Native equivalent |
| --- | --- |
| Empty constructor | `IMSWeights::new()` / `Default` |
| Masses/precision constructor | `from_masses(&[f64], precision) -> Result<Self>` |
| Copy constructor, assignment | Standard owned `Clone` and Rust assignment |
| `size` | `len`; also `is_empty` |
| `getWeight`, `operator[]` | `weight(index) -> Result<u64>`; borrowed `weights()` slice |
| `getAlphabetMass` | `alphabet_mass(index) -> Result<f64>`; borrowed `masses()` slice |
| `getPrecision` | `precision() -> Option<f64>` |
| `setPrecision` | `set_precision(precision) -> Result<()>` |
| `back` | `back() -> Result<u64>` |
| `getParentMass` | `parent_mass(&[u32]) -> Result<f64>` |
| `swap` | `swap(first, second) -> Result<()>` |
| `divideByGCD` | `divide_by_gcd() -> Result<bool>` |
| `getMinRoundingError`, `getMaxRoundingError` | `min_rounding_error`, `max_rounding_error`, both checked |
| Stream insertion | `to_text() -> Result<String>`, one integer and newline per row |

The source empty constructor leaves precision uninitialized. Native empty values
use `None`; no invented numeric precision is reported. A successful
`from_masses` or `set_precision` stores `Some(precision)`, including for empty
arrays. Empty parent-mass and rounding-error operations return zero without
requiring a precision, matching the source loops. Empty text is an empty string;
`back` or an invalid index returns an error rather than accessing invalid memory.
The native text method does not reproduce `std::endl`'s per-line stream flush.

All returned slices are immutable. A clone owns separate mass/weight vectors.
Swapping moves both members of the pair, so subsequent precision changes and
parent-mass evaluation use the new order. The utility never sorts masses.

## Numerical conventions and finite quirks

Integer scaling evaluates the literal expression
`floor(original_mass / precision + 0.5)`. It does not use a ties-to-even conversion
or round away from zero as a replacement. The original mass vector is retained
without normalization. Calling `set_precision`, including after GCD reduction,
recomputes weights from those original masses. A failure on any mass leaves the
old precision, vectors and their allocations unchanged.

Native integer storage is fixed `u64`; the source uses platform-dependent
`unsigned long`: usually 64-bit on LP64 Linux/macOS and 32-bit on Windows. This
is an explicit portable-width native policy. Every representable result from zero through the largest
representable f64 integer below `2^64` is accepted. The range check uses the exact
exclusive `2^64` boundary; it does not compare with `u64::MAX as f64`, which rounds
up. Values above `2^53` are supported. This is deliberately broader than the
private mass-decomposition solver's exact-integer search bound.

Finite signed masses and precision are accepted when this rounding expression
has a valid unsigned result. For example, mass `-0.5` at precision `1` rounds to
zero; masses `[-1, -2]` at precision `-0.5` give `[2, 4]`. Negative rounded integers,
nonfinite inputs and nonfinite/out-of-range scaled values are checked errors.
Zero precision is valid only for an empty alphabet, where no division occurs.
Zero original masses and weights are valid here even though a decomposition
solver cannot use a zero integer denominator.

GCD reduction computes the GCD of the stored integers, multiplies precision by it,
and divides integers by it **without rerounding original masses**. This preserves
rounding history. Fewer than two entries return `false`. A source control-flow
quirk is retained: exactly two coprime entries return `true` with unchanged values;
three or more coprime entries return `false`. For two or more all-zero weights the
source divides by zero; native code returns an error before mutation. A nonfinite
rescaled precision is also rejected atomically. All masses remain unchanged.

Parent mass is the ordered f64 accumulation
`sum(original_mass[i] * u32_count[i])`. It is independent of precision and GCD
scaling, accepts finite signed results, and preserves floating-point order rather
than sorting, reassociating or compensated summing. Length mismatch is an error,
including the source expected/actual size message. A nonfinite final result is a
checked native error.

The two rounding-error getters evaluate the source expression
`(precision * integer_weight - original_mass) / original_mass` in stored order,
starting their result at zero and selecting only negative/minimum or
positive/maximum values. Consequently, zero-mass `0/0` terms are ignored by the
source comparisons, and the native code preserves that behavior. An overflowing
term in the opposite direction does not invalidate a finite requested bound;
only a nonfinite **returned** bound is an error. Coarse precision can produce zero
weights for positive masses and a minimum relative error of exactly `-1`.

## Resource boundaries and reuse

The private container invariant limits each vector to one million elements,
checked before construction allocates. Two vectors contain at most 16 million
payload bytes; rebuilding precision temporarily adds at most eight million.
Standard cloning is bounded by that invariant. These are element-payload limits,
not guarantees about allocator metadata, capacity or physical process memory.
Fallible reservations are used by constructors, setters and text generation;
ordinary Rust `Clone` has ordinary allocation behavior.

Checked operations have a 50-million-unit work limit. GCD iterations are charged
cumulatively before each step and its final mutation is preflighted. Formatting
checks its complete length and work before allocation and limits generated text
to eight MiB, including all trailing newlines. Oversized input/output, invalid
indices, undefined conversions, exhausted work and selected nonfinite results
return errors rather than truncating or partially modifying values.

This is the complete standalone `Weights` operation surface. `IMSAlphabet`,
`IMSElement`, isotope-distribution utilities and standalone decomposition engines
are separate source APIs and are not claimed here.

The existing [mass-decomposition solver](MASS_DECOMPOSITION_ALGORITHM_SUPPORT.md)
already has a small private implementation of weight scaling/GCD/error bounds.
It was frozen independently and is not changed by this increment. A later narrow
reuse can move its shared-budget numeric setup through this utility while
retaining solver-only positive weights, its `2^53 - 1` bound, alphabet order and
cumulative work ledger. Reusing the utility must not silently narrow the utility's
signed/zero/full-u64 behavior or reset the solver's budget.

## Evidence

[The source manifest](../tests/data/ims_weights_provenance.json) records pinned
header, implementation, class test and GCD helper hashes. Native tests exercise
every distinct active source literal: the five-mass fixture at all four tested
precisions, the copy/access/swap/parent-mass values, GCD examples and four rounding
error values. Independent tests cover 2,800 rational rounding cases, ordered-sum
rounding, the full u64 boundary, signed/zero behavior, the coprime pair quirk,
resource guards, and rollback on later failures.

Validation is native Rust testing plus source inspection. No C++ execution,
ABI compatibility or differential-runtime parity is claimed.

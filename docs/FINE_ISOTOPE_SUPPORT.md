# Fine isotope patterns

`chemistry::FineIsotopePatternGenerator` returns individual isotope
configurations with accurate masses. It preserves the materialized pattern
conventions of pinned
[OpenMS4-core FineIsotopePatternGenerator](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.cpp),
using a native bounded enumerator. No optional Cargo feature or external runtime
is required. For nominal isotope envelopes and averagine estimates, see
[coarse isotope support](ISOTOPE_SUPPORT.md).

## Construction and output

The generator has one public setting, `stop: FineIsotopeStop`.
`run(&EmpiricalFormula) -> Result<IsotopeDistribution>` borrows the formula and
returns owned peaks sorted by increasing mass. Configurations remain distinct;
equal or nearby masses are not merged. Errors leave the formula unchanged and
never return an unnoticed partial pattern.

```rust
use openms::chemistry::{
    EmpiricalFormula, FineIsotopePatternGenerator, FineIsotopeStop,
};

fn main() -> openms::Result<()> {
    let formula = EmpiricalFormula::parse("C6H12O6")?;
    let generator = FineIsotopePatternGenerator {
        stop: FineIsotopeStop::UnexplainedProbability(0.01),
    };
    let pattern = generator.run(&formula)?;
    for peak in pattern.peaks() {
        println!("{:.6}\t{:.8}", peak.mass, peak.probability);
    }
    Ok(())
}
```

This selects three glucose configurations. `FineIsotopePatternGenerator::default()`
uses the same setting.

`run_with_isotopes(atom_counts, isotope_masses, isotope_probabilities)` accepts
custom populations using the same checked input format as the iterator below.
It applies the same stop conditions, final binary32 probability storage and mass
sorting as `run`, while preserving the custom input weights as binary64.

| Stop setting | Meaning |
| --- | --- |
| `UnexplainedProbability(0.01)` | Select configurations until their stored probabilities cover at least 0.99 |
| `AbsoluteThreshold(1e-5)` | Retain configurations whose absolute probability reaches the cutoff |
| `RelativeThreshold(1e-5)` | Retain configurations reaching this fraction of the globally most probable configuration |

Unexplained probability must be finite and in `[0,1]`; one requests an empty
distribution. Thresholds must be finite and nonnegative. Equality is retained
as an explicit native policy. Zero threshold requests the entire finite
support, subject to limits, including configurations whose output probabilities
underflow to zero. The most probable configuration can contain heavy isotopes:
it is not necessarily the monoisotopic configuration.

## Probability precision and coverage

Fine enumeration uses two source rounding steps. Natural abundance literals
are narrowed to binary32 and widened to binary64 before multinomial calculation.
Each output configuration probability is narrowed to binary32 again, then
widened into the native `IsotopePeak::probability: f64` field. Internal log
probabilities and isotope masses use binary64. This differs from the native
coarse generator's documented binary64 probability arithmetic.

Neither the input abundances nor the returned pattern are renormalized. Their
sum may be slightly above or below one. Threshold comparisons use internal
log probabilities, allowing very small cutoffs without prematurely discarding
underflowed output values. For example, the complete C100 support retains 101
configurations even though its smallest stored probabilities are zero.

Coverage sums the **stored, rounded probabilities** in binary64, selecting a
minimal prefix in descending probability order. This follows the source
materialized wrapper's trimming precision. For `C520H817N139O147S8` at coverage
0.99999, the source requires 19,615 configurations; accumulating unrounded
binary64 probabilities instead gives a different stopping boundary. Coverage
does not rescale the target to the actual abundance sum. If the complete support
is exhausted below the target, the full distribution is returned, consistent
with source exhaustion. Reaching a resource limit first is an error.

## Isotope labels, charge and formulas

Natural isotope entries with zero abundance are skipped. Explicit labels such
as `(13)C`, `D` or `(3)H` instead specify a fixed isotope with probability one,
including naturally absent tritium. Natural and labeled atoms of the same
element remain separate categories. A neutral empty formula has the identity
configuration `(mass=0, probability=1)` when the stop condition admits it.
Negative atom counts and negative formula charge are errors.

Positive charge follows the pinned source's deprecated convention: add that
many **natural hydrogen atoms**, clear the charge, then enumerate. Thus charged
glucose with charge two has the pattern of neutral `C6H14O6`. Masses are not
divided by charge, and no electron correction is applied. Use an explicitly
neutral formula for a neutral pattern; do not interpret these direct outputs as
charge-adjusted m/z. Theoretical-spectrum generation applies its separate ion
and charge conventions, described in [theoretical spectra](THEORETICAL_SPECTRA.md).

The generator needs an empirical composition. A peptide with known mass but
unknown formula cannot acquire a fine pattern through an invented composition.
The existing element table, including its documented native corrections,
supplies all isotope masses and abundances.

## Owning configuration iterator and custom abundances

`FineIsotopeIterator::from_formula(&formula)` returns an owning, fallible iterator
in descending configuration probability order. Each successful item is a
`FineIsotopeConfiguration` with binary64 `mass`, `probability` and
`log_probability`. Callers can take a prefix, stop at their own coverage target,
and resume the same iterator. Returned values remain usable after advancement
or after the iterator is dropped. No complete pattern is allocated first.

The raw formula adapter follows `IsoSpecWrapper`: it **ignores formula charge**,
including negative charge, and enumerates exactly the supplied atom counts.
This differs from the high-level pattern generator's hydrogen-adduct convention
above. Natural abundances still originate in the source's binary32 tables;
explicit labels still have unit probability. Negative atom counts are errors.

```rust
use openms::chemistry::FineIsotopeIterator;

fn main() -> openms::Result<()> {
    let stream = FineIsotopeIterator::from_isotopes(
        &[2],
        &[vec![12.0, 13.00335483507]],
        &[vec![0.25, 0.75]],
    )?.with_absolute_threshold(0.1)?;
    for item in stream {
        let configuration = item?;
        println!("{} {}", configuration.mass, configuration.probability);
    }
    Ok(())
}
```

Each custom row describes an independent population: an atom count, matching
nonempty isotope-mass and isotope-weight vectors. The native API derives the
category count from the vector length. Masses must be finite and nonnegative;
weights must be finite and strictly positive. Weights retain binary64 precision
and are never silently normalized or narrowed before enumeration. Rows need not
sum to one. Zero-count rows are validated; empty outer vectors give the empty
configuration `(0, 1, 0)`. Equal masses and repeated populations do not merge
distinct configurations. Input vectors are copied into owned search state.

The example returns the two most probable enriched-carbon configurations,
with weights 0.5625 and 0.375. `with_relative_threshold(value)` instead measures
the cutoff against the original global mode. Both consuming methods accept
finite nonnegative cutoffs and retain equality. They replace the cutoff on the
remaining iterator; lowering it cannot restore already consumed configurations
or revive an exhausted iterator. A zero cutoff admits the full support subject
to lifetime resource limits.

Iterator absolute cutoffs compare the returned binary64 probability directly.
Relative cutoffs compare the returned probability to the original mode's raw
probability when both are finite and positive; logarithmic comparison handles
underflow or an unrepresentable mode. A relative cutoff greater than one admits
no configuration. This preserves equality when callers reuse a returned raw
value, without an arbitrary numerical tolerance. The materialized generator
retains its established logarithmic threshold comparison, so the two interfaces
can differ at a floating-point rounding boundary.

Raw probability can underflow to zero while log probability stays finite.
`configuration.to_peak()` performs a checked conversion to `Peak1D`, narrowing
only its intensity to binary32. The mass is copied into `mz` without adding
protons or dividing by charge, as in the source wrapper; callers must supply
any intended ion convention. Raw binary64 values that cannot fit a finite
binary32 intensity produce a conversion error.

An iterator error is yielded once, after which the iterator is permanently
exhausted. A caller may therefore receive a valid prefix followed by a resource
or numerical error; using `collect::<openms::Result<Vec<_>>>()` preserves the
error. Materialized generator methods remain all-or-error. No current-value
accessor can be called in an uninitialized or stale state.

The [streaming example](../examples/stream_isotopes.rs) prints the first five
natural glucose configurations and two enriched-carbon configurations. A
[workflow test](../tests/fine_isotope_stream_workflow.rs) connects custom streams
to annotated spectra, checked selection and mzML interchange, preserving distinct
configurations at identical masses.

## Enumeration and checked limits

Each element's multinomial mode is found by allocating atoms to the strongest
marginal probability gains and checking improving transfers. A maximum heap
starts at the product of these modes. Every valid one-atom transfer within an
element produces a neighboring configuration; a visited set retains complete
configuration identity, including equal-probability neighbors. Multinomial
configurations have nonincreasing-probability paths from a mode, so this traversal
can visit global top configurations without enumerating the full support first.

Probabilities and masses are recomputed from counts rather than accumulated
along a discovery path. Log-factorial arithmetic avoids direct factorial
overflow. The implementation uses standard-library containers and bounded
tables; it does not call or reproduce IsoSpec's layered backend.

| Bound | Value |
| --- | ---: |
| Total atoms, including hydrogen added for charge | 1,000,000 |
| Materialized configurations per pattern | 100,000 |
| Custom input categories, including zero-count rows | 1,000,000 |
| Unique visited configurations | 250,000 |
| Conservative retained-payload allowance | 128 MiB |
| Estimated work units per independent call | 100,000,000 |

Memory accounting includes configuration vectors, heap/set capacity overhead,
factorial storage and output, rather than counting only emitted peaks. Work
charges preparation, state copying/hashing, probability evaluation, heap
operations and final sorting, accounting for isotope dimension. A zero-threshold
request whose full combinatorial support exceeds the output cap is rejected
before expensive enumeration. Positive-cutoff searches are allowed despite a
large full support, but frontier growth remains bounded. Nonfinite arithmetic
and exceeded limits return errors; binary32 underflow to zero is permitted.

`TheoreticalIsotopeModel::Fine { unexplained_probability }` reuses this generator.
All fine envelopes within one theoretical-spectrum operation share a cumulative
work allowance, including fragment and loss envelopes; each enumeration also
retains its own allocation limits. Existing theoretical output limits and
transactional append behavior still apply.

The iterator shares the atom, visited-state, memory and work bounds over its
whole lifetime. It does not impose the materializer's 100,000-output cap or
precompute full-support size, so a small prefix can be requested from a much
larger support. It retains its frontier and visited set; streaming is not a
constant-memory algorithm. The last emitted state's neighbors are expanded only
when another configuration is requested, so stopping or dropping the iterator
does not perform unused enumeration work. Custom table cells, including those
in zero-count rows, count toward input work and storage checks before copying.

## Compatibility boundary and evidence

The generator produces materialized, mass-sorted patterns; the iterator supplies
ordered raw configurations and thresholds for natural or custom populations.
IsoSpec's layered traversal, performance-hint behavior and untrimmed layer-based
selection are not reproduced. The source total-probability stream's hint is not
a stop condition: its complete support can be obtained through the native ordered
iterator, while callers may stop explicitly at their own coverage target.
The backend source was unavailable in
the pinned snapshot, so exact threshold equality and final tied-configuration
membership are documented native policies. Equal probabilities use deterministic
configuration ordering; IsoSpec layer and pivot order are not claimed.

The [independent reference review](FINE_ISOTOPE_REFERENCE_REVIEW.md) and
[provenance manifest](../tests/data/fine_isotope_provenance.json) distinguish
upstream literals from derived checks and record twelve source hashes. Tests
cover 44 exact source count cases, fourteen fructose mass/probability references,
six bromine configurations, the insulin coverage boundary and carbon tail
underflow. Source-count examples include 2,548 complete glucose configurations,
20,503 complete C100H202 configurations and the 19,615-configuration insulin
coverage result.

[Implementation tests](../tests/fine_isotopes.rs) independently enumerate all
54 C2H2O2 configurations, check labels/charge, errors and real frontier-limit
exhaustion. Private tests check equal-probability identity and shared work.
[Theoretical integration tests](../tests/theoretical_fine_reference.rs) use a
separate direct chemical oracle for ion masses, intensities and losses. Focused
checks pass on current Rust and Rust 1.85. No C++ program was built or executed.

The [streaming reference review](FINE_ISOTOPE_STREAM_REFERENCE_REVIEW.md) records
the raw-wrapper charge and binary64 conventions, upstream traversal/count cases,
custom-input oracles and the limits of backend-order comparisons.

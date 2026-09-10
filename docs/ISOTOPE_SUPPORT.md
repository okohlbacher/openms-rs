# Isotope distribution support

`openms::chemistry::isotopes` provides native coarse isotope patterns and mass
estimation without external dependencies. It uses the ported element tables and
preserves the principal numerical conventions of [OpenMS4-core revision
7c029e8cdba6abab503708ecdd56f6ab55e38ce4](https://github.com/okohlbacher/OpenMS4-core/tree/7c029e8cdba6abab503708ecdd56f6ab55e38ce4).

## API coverage

| API | Behavior |
| --- | --- |
| `IsotopePeak`, `IsotopeDistribution` | Validated mass/weight container; identity and empty constructors; insert, resize, clear, min/max/most-abundant queries, sorting, trimming, renormalization and weighted average mass |
| `CoarseIsotopePatternGenerator::run` | Nominal isotope convolution, gap filling, exponentiation by squaring, isotope labels, positive-charge behavior and normalization of retained peaks |
| `convolve`, `convolve_power` | Raw nominal-bin convolution and powers with unnormalized weights |
| `set_isotope_override` | Generator-local isotope enrichment; shared element tables remain immutable |
| `CoarseMassMode` | Approximate carbon-13 mass spacing or rounding of those corrected masses |
| `AveragineComposition`, `FormulaEstimate` | Peptide, RNA, DNA and custom average compositions; average/monoisotopic mass estimation; explicit hydrogen-adjustment status; optional fixed sulfur count |
| `estimate_from_*` | Pattern generation from peptide/RNA/DNA average weights, peptide/RNA monoisotopic weights, custom compositions and peptide weights with fixed sulfur |
| `calc_fragment_isotope_dist` | Joint fragment/isolation probabilities given independent fragment and complementary-fragment distributions and selected precursor isotope indices |
| `estimate_fragment_from_weights` | The same fragment calculation using independently estimated fragment/complement formulas with any supported average composition |
| `approximate_intensities`, `approximate_from_peptide_weight` | Normalized truncated Poisson approximation with lambda = mass / 1800, source recurrence with log-space overflow fallback |

## Mass and probability conventions

- Isotope probabilities are grouped by **nominal neutron count**. Contributions
  such as carbon-13 and nitrogen-15 enter the same coarse bin. This does not
  generate fine structure or probability-weighted exact masses within each bin.
- `Approximate` labels peak `i` as the formula's **lightest isotope mass** plus
  `i * C13C12_MASSDIFF_U`. This may differ from its most-abundant-isotope mass;
  selenium is a tested example. Explicit isotope labels use the specified
  isotope's tabulated mass and have probability one unless locally overridden.
  For a modified sequence, pass `sequence.formula()?`: isotope labels and atom
  deltas are retained. Formula-based mass labels can differ slightly from
  `sequence.mono_mass()` because OpenMS uses separately tabulated terminal
  modification mass deltas in the latter calculation. Unresolved residues and
  mass-only tags return a formula error; they cannot define an isotope envelope.
- `Nominal` rounds each corrected approximate mass, matching the source. It does
  not simply add integer indices to a rounded starting mass. At sufficiently
  high indices, this can produce a two-dalton label increment.
- `run` retains the lowest requested bins and renormalizes them to sum to one.
  These are conditional probabilities within the retained range; a truncated
  result does not report the original discarded probability. Raw convolution
  and fragment-isolation calculations do **not** renormalize.
- Positive formula charge preserves the pinned source's deprecated behavior:
  extra natural hydrogen atoms contribute isotope probabilities, while the mass
  origin gains `charge * PROTON_MASS_U`. Returned coordinates are masses, not
  m/z: there is no division by charge. Negative charge and negative atom counts
  are rejected. Explicit adduct compositions should be represented by atom
  counts with zero charge. This implementation does not silently adopt the
  planned future OpenMS charge convention.
- Natural isotope support includes declared zero-abundance tail bins, such as
  tritium in hydrogen. Trimming is explicit. Raw convolution rounds its input
  masses and requires strictly increasing distinct nominal bins; gaps are filled
  with zeros. Sorting distributions by probability before convolution therefore
  requires sorting them back by mass.
- Fragment calculations use **positions in the provided arrays** as successive
  isotope indices starting at M0; their mass values are ignored. Supply complete
  index sequences, including zero bins. Do not pass trimmed or sparse vectors
  whose first position is not M0. Duplicate precursor indices are ignored.
  Renormalize the returned joint weights for true conditional probabilities.
- The Poisson layout starts at exactly its supplied mass and spaces peaks by
  `NEUTRON_MASS_U / charge`, matching the source. The starting mass is neither
  divided by charge nor protonated.

## Validation and intentional improvements

Probabilities use `f64`, whereas C++ stores isotope weights in `Peak1D`'s `float`
intensity field. Golden comparisons allow the source's rounding precision.
Inputs and container mutations reject negative/nonfinite masses or weights.
Weights above one are allowed before normalization. Empty normalization is a
no-op; a nonempty zero-sum or overflowed-sum distribution returns an error and
stays unchanged. An empty distribution has no min/max/most-abundant value and no
weighted average mass.

`trim_left` removes every peak when all weights fall below the cutoff, fixing an
upstream edge case that leaves the distribution unchanged. Trimming retains
weights equal to the cutoff. `resize` pads with `(0, 0)`, as in C++, and can
therefore change mass ordering. Fragment truncation bounds the complete loop,
preventing the source's potential out-of-bounds access when the configured result
length is shorter than its input fragment distribution.

Local isotope overrides must use declared isotope mass numbers, retain the
lightest declared isotope entry (zero weight is allowed), and have a finite,
positive total weight. This keeps the mass anchor and bin origin consistent;
for example, carbon enrichment can be `[(12, 0), (13, 1)]`. An override containing
only the heavy isotope is rejected rather than silently shifting its peak to the
lightest-isotope mass. Explicit labeled isotopes are distinct override keys.

Averagine uses rounded heavy-element counts followed by hydrogen adjustment to
the requested mass. When negative hydrogen would be needed, `FormulaEstimate`
returns the heavy-element formula with `hydrogen_adjustment_succeeded = false`.
Pattern convenience methods follow upstream behavior and still use that formula.
Call the composition estimator directly to inspect this status. Fixed-sulfur
estimation inserts sulfur safely and rejects a sulfur mass greater than the
target mass. Zero target mass yields an empty formula and an identity pattern.

## Resource limits

- Any isotope vector, including temporary gap-filled vectors, is limited to
  1,000,000 peaks (`MAX_ISOTOPE_PEAKS`).
- A coarse `run`, raw convolution or convolution-power call is limited to
  50,000,000 probability products (`MAX_CONVOLUTION_PRODUCTS`). Powers share the
  budget with all convolutions inside a run. Fragment weight estimation performs
  separate bounded runs for the fragment and its complement.
- Nominal convolution masses must fit the exact integer range of `f64`.
- Overflowed probabilities, invalid count estimates and a retained distribution
  that underflows entirely to zero return errors.

Request fewer peaks when an unconstrained calculation exceeds these limits.
`None` means all peaks within the limits; `Some(0)` is invalid. Limits are checked
before large result allocations or convolution loops. Poisson estimation uses the source recurrence and running normalization sum
for finite weights, falling back to log space when those values overflow.
This preserves source rounding for ordinary model-based deisotoping while
keeping extreme finite inputs numerically defined.

## Deferred capabilities

[FineIsotopePatternGenerator](FINE_ISOTOPE_SUPPORT.md) is implemented separately
with bounded native configuration enumeration and source precision conventions.
Fine support also includes ordered streaming and custom isotope populations,
without emulating IsoSpec layered traversal. Arbitrary-resolution distribution merging,
fixed-sulfur fragment-estimator convenience wrappers, and the legacy formula
isotope-count estimator remain outside this module. Fine states must not be inferred from the
coarse mass labels.

## Provenance and verification

Reviewed source SHA-256 values at the pinned revision:

| Source | SHA-256 |
| --- | --- |
| `src/openms/source/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.cpp` | `987934145e3cf0d5f5eff97fc580ef2c680e8b3a932d36c8c4855d4ef063bf38` |
| `src/openms/source/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.cpp` | `0974fabc36a93f741c210d89e940215fa773306b99b5b45309a35cb4b43a2f8c` |
| `src/openms/source/CHEMISTRY/EmpiricalFormula.cpp` | `7ebd85a51d761174faf90068997b2348d6291641d231db0a7f2131cdc3d614ad` |

`tests/isotopes.rs` ports golden glucose, charged glucose, large-formula,
bromine-gap, averagine, fixed-sulfur and conditional-fragment cases from
`CoarseIsotopeDistribution_test.cpp`. Independent tests cover probability sums,
raw versus normalized truncation, isotope enrichment, labeled formulas,
normalization failures, empty/single distributions, Poisson probabilities,
conditional-probability arithmetic, invalid inputs and resource limits.

```text
cargo test --offline --test isotopes
cargo clippy --offline --test isotopes -- -D warnings
```

No C++ build or live C++/Rust differential run was performed.

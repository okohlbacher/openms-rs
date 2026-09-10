# Precursor purity and SPS matching

`analysis::precursor_purity::PrecursorPurity` ports the scalar, fuzzy-scan,
interpolation, batch and SPS methods from pinned
[PrecursorPurity.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/ANALYSIS/ID/PrecursorPurity.cpp).
These native Rust APIs have no feature dependency. They borrow inputs and return
owned results; an error never changes the experiment, spectrum or precursor.

## API and acquisition fields

All five methods are associated functions of `PrecursorPurity` and return
`Result`:

| Method | Inputs | Successful result |
| --- | --- | --- |
| `compute` | `&MSSpectrum`, `&Precursor`, `Tolerance` | `PurityScores` |
| `count_sps_matches` | `&[Precursor]`, `&AASequence`, `Tolerance`, maximum fragment charge `u8` | Matching precursor count `usize` |
| `compute_single_scan` | `&MSExperiment`, selected spectrum index, parent index, maximum deviation in ppm | `Vec<f64>`, one value per selected-spectrum precursor |
| `compute_interpolated` | The same inputs, with `Option<usize>` next-parent index before the ppm argument | `Vec<f64>` |
| `compute_all` | `&MSExperiment`, `Tolerance`, `ignore_missing: bool` | `BTreeMap<String, PurityScores>`, keyed by MS2 native ID |

`Tolerance` is `comparison::Tolerance::{Absolute(f64), Ppm(f64)}`. There is no
implicit default tolerance. The scalar and SPS methods deliberately use
different source tolerance arithmetic, described below.

`kernel::Precursor` directly holds selected-ion m/z, intensity, charge, isolation
offsets, optional distinct isolation target and optional `spectrum_reference`,
alongside activation and mobility metadata. `PrecursorInfo` is a compatibility
wrapper over this same record. Both scalar and fuzzy purity center their windows
on **selected-ion `mz`**, using the lower and upper offsets; a distinct
`isolation_target_mz` is retained for interchange but does not recenter these
source algorithms. See [mzML support](MZML_SUPPORT.md) for acquisition-field
preservation and serialization limits.

## Scalar isolation purity

`compute` takes peaks in the inclusive window
`[mz - lower_offset, mz + upper_offset]`. `PurityScores` contains binary64
`total_intensity`, `target_intensity` and `signal_proportion`, `usize`
`target_peak_count` and `interfering_peak_count`, and an owned
`interfering_peaks: MSSpectrum`. The residual spectrum has unmatched peaks in
their original m/z order and default metadata, matching the source's newly
constructed isolated spectrum. It does not copy unrelated input annotations.

Charge uses its unsigned magnitude, with zero interpreted as one. Isotope
spacing is `1.0033548378 / charge`, the source C13/C12 mass difference. The
initial negative isotope index is `-trunc(lower_offset * charge)`; if its
position lies below the isolation window, the index is incremented once.
Positions remain anchored at `mz + isotope_index * spacing` with the source
multiply/divide operation order retained.

Each expected isotope selects the nearest remaining isolated observation.
Equal-distance ties favor the lower m/z. Matching endpoints are inclusive, and
the tolerance is doubled: `2 * value` for absolute tolerance, or
`mz * value * 2 * 1e-6` for ppm, in that order. A matched observation is removed
before the next search, preventing reuse within the call. This is a greedy
assignment, without an averagine model or a requirement that the monoisotopic
observation be present. The latter behavior differs from the source header's
description but matches its implementation.

Observed binary32 intensities accumulate in binary64, starting from positive
zero. The ratio is target/total when target intensity is positive, otherwise
zero. An empty isolation window returns default zero scores. Finite results
are not clipped; the implementation preserves source arithmetic rather than
adding probability normalization.

## SPS fragment matching

`count_sps_matches` uses the native compact theoretical helper's default b/y
series, including both the first fragment and full peptide length, for charges
one through the supplied maximum. Maximum charge zero means one. Peptide
terminal and residue modifications follow the helper's mass conventions;
known mass-only annotations work without invented formulas, while unresolved
required masses return an error. Empty precursor lists or empty peptides return
zero before requesting peptide chemistry.

Each precursor is tested independently; duplicate windows count separately.
Precursor charge and isolation offsets do not select the theoretical charges or
define the match window. Absolute tolerance is used once, without doubling.
Ppm width is `abs((value / 1e6) * mz)`, preserving division-first arithmetic.
Both binary64 endpoints are narrowed to binary32 before inclusive comparison
with the sorted binary32 theoretical values. Consequently, zero tolerance can
accept a distinct binary64 m/z that rounds to the same binary32 value.
Nonfinite endpoints, including overflow during narrowing, return errors.

## Fuzzy scan purity and interpolation

`compute_single_scan` preserves a separate source algorithm. Its selected and
parent indices must exist, but the method does not require their MS levels to
form a parent/child pair. Each initial seed is the globally nearest parent peak,
without a tolerance or isolation-window requirement. Charge zero means one;
negative charge is rejected because the source's signed stepping can fail to
terminate. Spacing uses the neutron mass `1.00866491566 / charge`, not C13/C12.

Strict isolation bounds are expanded to fuzzy bounds by multiplying the lower
bound by `1 - ppm/1e6` and the upper bound by `1 + ppm/1e6`. Isotope matching
requires a deviation strictly less than the ppm tolerance. Each successful
match anchors the next step at that observed peak. At exact strict boundaries
and beyond, contributions receive half intensity; exact outer fuzzy boundaries
are excluded from the neighboring denominator scans. The initial seed is
always included at full intensity. Numerator and denominator accumulate in
binary32, and their binary32 ratio is widened to binary64. Half contributions
retain the source's promotion and rounding order.

The source isotope lookup compares `lower_bound(expected)` with its physical
successor, favoring the successor on equal distance; it does not compare the
predecessor. This can ignore a closer lower observation. Native code preserves
physical candidates beyond the logical search interval when they exist. Such
a candidate may contribute to target intensity but not to the denominator,
producing a finite purity above one. At physical end, a sole remaining candidate
is used, or no match is reported, instead of dereferencing invalid iterators.

An empty parent yields a vector of ones after numerical validation. A zero-width
isolation window stops processing: preceding results remain, while that
precursor and all subsequent entries stay one. Nonpositive totals, invalid
boundaries, overflow and non-advancing isotope traversal return checked errors.

`compute_interpolated` calculates the earlier scan first. `None`, an out-of-range
next index, a next scan whose level is not MS1, or a nonpositive/nonfinite absolute
parent RT difference returns the earlier result without using the next scan's
peaks. Otherwise it also calculates the later purity and returns
`abs(selected_rt - earlier_rt) * ((late - early) / abs(later_rt - earlier_rt)) + early`.
Absolute RT distances preserve source extrapolation; values may be below zero
or above one and are not clipped. Other numerical errors propagate.

## Batch parent selection

`compute_all` walks scans in acquisition order, processing only MS2 spectra and
only their first precursor. It first looks for an earlier MS1 with the native ID
named by `spectrum_reference`, then falls back to the most recent preceding MS1.
Later references do not select future scans. Repeated parent IDs resolve to the
most recent preceding matching scan. Retention-time sorting is not performed.
The result uses the selected earlier parent only; the source header's claim of
combining preceding and following scans does not describe this method.

An empty experiment returns an empty map. Every scan must have a positive MS
level, and MS2 native IDs must be nonempty, control-free and unique. When
`ignore_missing` is false, the first scan must be MS1 and every MS2 must have a
parent. When true, an MS2 with no parent receives default zero scores. An MS2
with a usable parent but no precursor still returns an error. Suitability errors
replace source warning/empty-map and unchecked indexing paths; no partial map
is returned. Non-MS2 records do not create result entries.

## Validation, work and storage bounds

Validation is scoped to numerical data used by purity: checked spectra have
finite RT, positive MS level, and sorted finite nonnegative m/z and intensity.
Duplicate m/z values and signed zero are allowed. Checked precursors have finite
nonnegative selected m/z, intensity, isolation offsets and optional isolation
target. Tolerances are finite and nonnegative; there is no unrelated binary32 ppm
limit. Fuzzy calculations also require nonnegative computed isolation/fuzzy
boundaries. Charge policies differ as specified above.

Unused auxiliary arrays, peptide identifications, CV collections and other
acquisition metadata are **not traversed or validated**. This prevents repeated
parent lookup from rescanning arbitrary annotation payloads outside the work
budget. Call the kernel's full validation separately when needed. Batch checks
only the parent numerical spectra and selected first precursors used in its
calculations; fuzzy methods check the selected and parent spectra and all
selected-spectrum precursor scalars. Interpolation's documented fallback may
skip an unusable next scan.

| Bound | Value and scope |
| --- | --- |
| `MAX_PURITY_PEAKS` | 1,000,000 peaks per numerically checked spectrum |
| `MAX_PURITY_PRECURSORS` | 100,000 precursors per checked spectrum or SPS list; also 100,000 batch scans |
| `MAX_PURITY_ISOTOPES` | 1,000,000 estimated positions per scalar call; a combined conservative traversal allowance across fuzzy windows in one scan calculation |
| `MAX_PURITY_WORK` | 50,000,000 work units per public call; shared by all batch calculations or both interpolated scans |

Scalar work charges peak validation, isolation copies/sums, isotope positions,
binary searches and physical vector removals. Fuzzy planning charges conservative
search/traversal and total-intensity scan bounds before allocating its output or
walking ladders. Batch also charges scan-ID/reference processing and lookup work.
The scalar signed-int isotope initialization is checked against the source's
integer range before conversion. All paths check finite arithmetic and progress;
undefined source conversions, iterator accesses and NaN results become errors.

SPS additionally inherits the compact helper's independent limits: 4,096
residues, 100,000 theoretical values and 10,000,000 estimated generation/sort work
units. Its purity budget covers precursor processing and matching. Output storage
is proportional to isolated residual peaks, precursor count or batch results;
batch results can retain separate residual copies from the same parent, bounded
by their charged work. No whole annotated input spectrum is cloned.

## Verification and provenance

The [independent review](PRECURSOR_PURITY_REFERENCE_REVIEW.md) distinguishes
literal upstream goldens from derived numerical cases. The
[provenance manifest](../tests/data/precursor_purity_provenance.json) records ten
pinned source hashes, extraction details and fixture hashes. Its 72 observations
retain exact IEEE bits from both upstream MS1 scans, alongside all seven scan
records, seven scalar/map cases and eight literal SPS cases. Independent decoding
verified the compact scalar windows against all 3,872 original MS1 observations.
No C++ build or C++ runtime comparison was performed.

[Scalar/SPS tests](../tests/precursor_purity.rs) cover exact totals, inclusive
boundaries, missing monoisotopes, greedy ties/removal, signed zero, charge policies,
binary32 SPS windows, error bounds and input preservation. The annotation-work
regression reuses 10,000 parent placeholders across 10,000 MS2 scans.
[Fuzzy tests](../tests/precursor_purity_fuzzy.rs) cover source successor selection,
physical-end safety, half weights, binary32 arithmetic, unclipped ratios,
interpolation fallbacks and malformed scalars. The
[reference suite](../tests/precursor_purity_reference.rs) independently checks
source-derived numerical and parent-selection cases. The
[workflow suite](../tests/precursor_workflow.rs) connects the shared precursor
model, batch purity, mzML acquisition round-trips and legacy-format loss guards.

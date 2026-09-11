# MassTrace

The native `kernel::MassTrace` implements the complete public MassTrace value,
container and computational surface in Core SDK
`82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. The source header, implementation and
class test are byte-identical to revision 54a232f. See the
[implementation](../src/kernel/mass_trace.rs), [tests](../tests/mass_trace.rs) and
[hashed provenance](../tests/data/mass_trace_provenance.json). This increment
builds on the existing Peak2D and scan-envelope ConvexHull2D representations; it
adds no dependency and does not build or execute C++.

## Public mapping

| Source operations | Native operations |
| --- | --- |
| Default, vector/list construction, copy/assignment | `new`/`Default`, owned `from_peaks`, checked `from_slice`, standard collection into Vec and `Clone`/assignment |
| `operator[]`, mutable/const forward/reverse iterators, `getSize` | Index/IndexMut, checked `get`/`get_mut`, `peaks`/`peaks_mut` slices, `iter`/`iter_mut` and their double-ended iterators, borrowed IntoIterator, `len`/`is_empty` |
| Quantification enum, names, string lookup, setter/getter | `MassTraceQuantMethod::{Area,Median,MaxHeight}`, `NAMES`, `ALL`, exact `from_name`, `set_quant_method`/`quant_method` |
| Label, centroid MZ/RT/SD/IM, IM presence, FWHM and borders | Corresponding snake_case getters; checked label/SD/IM setters; `contains_im_data` |
| Public average m/z and IM FWHM fields | `fwhm_mz_avg`, `fwhm_im_avg` |
| Smoothed intensity getter/setter | Borrowed `smoothed_intensities`, checked copying `set_smoothed_intensities` |
| Trace length, average scan time | `trace_length`, `average_ms1_cycle_time` |
| Raw/smoothed peak area, intensity sum, maximum index | `compute_peak_area`, `compute_smoothed_peak_area`, `compute_intensity_sum`, `find_max_by_int_peak` |
| FWHM estimate, raw/smoothed FWHM area | `estimate_fwhm`, `compute_fwhm_area`, `compute_fwhm_area_smooth` |
| Quantified intensity, maximum intensity, hull | `intensity`, `max_intensity`, `convex_hull` |
| All eight centroid-update operations | `update_weighted_mean_rt`, `update_smoothed_weighted_mean_rt`, `update_smoothed_max_rt`, `update_median_rt`, `update_median_mz`, `update_mean_mz`, `update_weighted_mean_mz`, `update_weighted_mz_sd` |

Unknown quantification names return `None`, replacing the invalid source enum
sentinel; a typed setter cannot accept that sentinel. Peaks retain f64 coordinates
and f32 intensities. Smoothed data uses f64. IM is an explicitly set scalar
centroid and presence flag, just as in source; no per-peak mobility unit or
additional scan data is invented.

Peak mutation is exposed as a fixed-length slice, matching the source's mutable
iterator/index surface. Construction does not calculate centroids. Later peak
edits leave centroids, FWHM, borders and supplied smoothing untouched until the
caller requests an update. Clone retains every stored field, including the label,
IM flag, smoothing, cached borders, average widths and native limits. Native
PartialEq compares complete stored values; no Eq claim is made for float fields.

## Numerical and state contracts

Raw area is the ordered trapezoidal integral, including the first source
zero-width term. It is not an intensity sum. Finite signed intensities and
descending/unsorted RTs remain accepted where source accepts them: signed RT
differences can produce negative areas and average scan times, while trace length
uses the absolute endpoint difference. No operation silently sorts the peaks.

`compute_smoothed_peak_area` preserves an unusual source expression: the first
previous intensity is smoothed, subsequent previous/current intensities are raw,
and each interval is included only when its current smoothed intensity is
positive. It is not a conventional smoothed trapezoidal integral. By contrast,
`compute_fwhm_area_smooth` uses smoothed values throughout.

Weighted RT uses right-hand rectangles starting at the second peak, including
the preceding RT gap. The first peak's intensity is unused. Smoothed weighted RT
uses positive smoothed intensities without RT gaps. Weighted m/z uses all raw
intensities. Its total weight, smoothed RT's positive total weight, and the SD
weight must meet the exact source `f64::EPSILON` threshold. Raw weighted RT has
no epsilon threshold, so any finite nonzero denominator that produces a finite
result is allowed. Weighted SD preserves `exp(2 * ln(abs(mz - cached_mz)))` and
the separate square-root division; `ln(0)` followed by `exp(-infinity)` is valid
zero variance. Callers must update the desired m/z center first.

FWHM requires nondecreasing finite RT, chooses the first maximum, walks outward
over values at least half maximum, and interpolates each bracket with its paired
coordinates. It uses the endpoint when that flank never falls below half height.
Duplicate RTs and equal bracket intensities use the source fallback coordinate.
FWHM areas integrate the **whole inclusive index bracket**, not truncated
interpolated half-height endpoints.

An apex at either endpoint returns zero and resets borders to `(0, 0)` while
leaving the previous cached FWHM value unchanged. Unset/cleared borders give area
zero even if smoothing has never been supplied. Empty intensity sum and raw area
are zero; empty maximum intensity is zero, while maximum-index lookup rejects an
empty trace. Maximum intensity starts at zero, so a wholly negative trace returns
zero there while maximum-index lookup still finds its largest negative peak.
Median quantification always uses raw intensities, ignoring its smoothing flag.

Singleton centroid updates preserve the source early return and signed-zero
coordinate without reading unused intensity or m/z fields. Smoothed singleton
updates still require smoothing but do not require positive intensity. Ordinary
read-only getters do not validate unrelated fields or recompute cached values.

## Checked boundaries and work limits

Computations reject nonfinite consumed inputs/results, invalid interpolation
brackets, and division/square-root outcomes outside the finite native contract.
Previously undefined empty smoothed-area/median accesses and absent smoothing for
nonzero FWHM borders return errors. Smoothed setters require exactly one finite
value per peak; SD and IM setters require finite scalars. Peak value construction
and borrowed mutation retain the underlying Peak2D storage policy, with checks
when a numerical operation consumes those fields. All fallible updates calculate
their new state before publishing it, including both FWHM cache fields.

`MassTraceLimits` is public and configurable on the trace and on constructors.
Default per-call ceilings are one million peaks, 50 million weighted work visits,
and 256 MiB cumulative logical new payload. Point counts are checked before
loops, and scan/copy costs are precharged. Median and hull sorting use a
conservative `32 * n * bit_length(n)` comparison allowance. Hull preflight also
covers its temporary Point2D input, geometry Scan output and stable-sort scratch;
the existing geometry implementation then constructs the same scan envelope.
This is a conservative logical bound, not exact allocator or resident-memory
measurement. A large sort may hit the work ceiling below the peak-count ceiling.

Slice constructors and label/smoothing replacements reserve fallibly before
copying. Owned-vector construction needs no new peak allocation and retains the
input bits and order. `from_peaks` does not meter allocation already performed by
its caller. Ordinary standard-library Clone, indexing, iterator callbacks,
equality and caller destruction are outside these checked-operation budgets;
checked `get`/`get_mut` provide a non-panicking alternative to indexing. No source
input controls recursion and no numerical operation resizes the peak container.

## Evidence

Twelve focused tests cover source area/FWHM/centroid and quantification fixtures,
the three explicit symmetric/asymmetric/open-flank examples, missing-scan area
invariance, hull containment, singleton/empty/negative/zero rules, cached state,
borrowed ownership and late atomic/resource failures. One test independently
enumerates 100 deterministic small irregular grids, comparing trapezoidal area,
normalized-weight centroids and a directly squared variance expression.

Some scalar expectations in the source class test are approximate despite many
printed decimal places. For its seven exact input peaks, f32 intensities sum to
`69831325` (the source expected line says `69831326`); sorted m/z median is
`230.10223` (source line says `230.10198`); arithmetic mean is
`230.10205142857143` (source line says `230.101918`). The native suite verifies
these independently derived values and the source's precise area/FWHM literals,
with explicit per-assertion tolerances. It does not treat broad source comparison
macros as an exact numerical oracle or alter the source expressions to fit stale
expected constants.

The 12 direct tests plus 19 existing feature/geometry tests and six Peak2D tests
pass on Rust 1.98 with all features and Rust 1.85 without default features.
Scoped strict library/test Clippy passes on both compilers, and changed-file
rustfmt checks pass. Root reviewed the complete implementation and geometry cost
path; an independent chemistry reviewer checked the source numerical and cached
state expressions. Both reviews closed without an outstanding finding.

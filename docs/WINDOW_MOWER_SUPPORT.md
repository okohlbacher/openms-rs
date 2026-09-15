# WindowMower

`processing::window_mower::WindowMower` retains the most intense peaks in m/z
windows. Defaults are source `window_size = 50.0` Th, `peak_count = 2`, and
`WindowMowerMethod::Sliding`. `FromStr` accepts source movement names `slide` and
`jump`. The width must be finite and positive; a zero peak quota is valid.
Finite signed intensities and coordinates are accepted, and input need not be
sorted. All represented acquisition metadata and aligned float, integer and
string arrays are preserved through original-index selection, including empty
array placeholders.

`retained_indices` returns original input indices in result order.
`filtered_spectrum` returns an owned result. The `SpectrumFilter` methods mutate a
spectrum or every spectrum of an experiment only after successful computation;
`filter_sliding` and `filter_jumping` explicitly select a movement method without
changing the configuration. Experiment filtering retains chromatograms and
experiment metadata, and applies the point/work limits to each spectrum
separately (see [Checks and computational limits](#checks-and-computational-limits)).

## Source semantics

Both methods first arrange an index view by ascending position. Every window
starts at an observed position and excludes points whose distance from its start
is equal to the configured width. There is no intensity averaging or m/z change.

**Sliding** advances its start by one sorted raw peak. It stops immediately after
the first window containing the last raw peak; it does not process progressively
shorter trailing windows. The union of chosen **m/z positions** determines the
result: every peak at a selected position survives, even if its own intensity was
not selected. Output remains in original input order.

**Jumping** starts the next window at the next observed peak after the current
window, including across large gaps. Completed windows use the full peak quota.
The final window uses
`round((last_position - window_start) / window_size * peak_count)`, with positive
half-integers rounded upward, matching `std::round`. A singleton final window
therefore contributes no peak; an entirely singleton spectrum becomes empty.
The selected **(m/z, intensity) pairs** determine membership. All identical
original peaks survive when one is selected, potentially exceeding the quota;
coincident peaks with different intensities remain distinct. Output is sorted by
position. IEEE positive and negative zeros compare equal in membership, just as
source peak equality does.

The native implementation resolves equal-intensity ties deterministically using
ascending m/z and then original input index. The pinned C++ intensity sort does
not specify its tie order. Jumping output also retains original relative order
among equal-position peaks. Otherwise the source's finite window and membership
behavior is retained. In particular, its last-window quota uses the observed
span, not the distance to an ideal bin edge. The source triangle test comment
claims a single final peak for `round(0.9 * 2)`; the implementation and its total
20-peak assertion both require two.

## Checks and computational limits

Source `WindowMower` has no resource ceilings at all. The port bounds its work,
but **every ceiling is per spectrum**, and the work ceiling is derived from that
spectrum's point count:

| Field | Default | Meaning |
| --- | --- | --- |
| `max_points` | 1,000,000 | Input peaks in one spectrum. |
| `max_work` | 50,000,000 | Work units for a spectrum before its points are credited. |
| `work_per_point` | 32,768 | Work units credited per input peak, so the ceiling for `n` points is `max_work + work_per_point * n`. |

`max_points` and the work ceiling must be positive. Each spectrum's point count
is checked before processing, and every selection is planned before cloning
spectra or committing any change. Empty records still consume one work unit.
Invalid native array lengths, nonfinite peaks, overflowing coordinate spans used
during window processing, or exhausted limits return errors without partially
changing the input. Invalid configuration is rejected even for empty input. Zero
quotas select nothing without sorting or processing coordinate spans.

### Why the ceilings are per record and size-derived

Until the benchmark of 2026-09-15 the two ceilings were charged **across a whole
experiment**: `filter_experiment` summed `spectrum.len()` over the map, compared
that against `max_points`, and ran one `max_work` ledger through every spectrum.
A run-wide ledger shrinks as the run grows, so it refuses real data that each of
its records passes comfortably. Measured on the benchmark inputs:

| Input | Spectra | Total peaks | Largest spectrum | Largest spectrum's work | Peak work per point |
| --- | ---: | ---: | ---: | ---: | ---: |
| `sub_centroid_uk222_picked_first600` | 600 | 186,536 | 551 | 135,039 | 270.6 |
| `sub_profile_uk222_first600` | 600 | 2,585,718 | 7,641 | 25,233,450 | 3,492.5 |
| `sub_centroid_uk222_picked_first5000` | 5,000 | 2,139,727 | 3,247 | 6,304,927 | 1,975.9 |
| `50amol_R1` (1.2 GB LTQ Orbitrap Velos) | 43,745 | 88,434,492 | 16,766 | 45,766,084 | 3,311.2 |

The Velos run exceeds a run-wide million peaks 88-fold, so `SpectraFilterWindowMower`
refused it and the 5,000-spectrum slice outright, while the C++ tool processed the
same file in 710 s at one thread. Its largest single spectrum is 60 times inside
the same million. The metadata-copy ledger that `filter_experiment` also ran once
over the whole run has the same shape and is now metered per spectrum as well,
as [`MorphologicalFilter`](MORPHOLOGICAL_FILTER_SUPPORT.md) already does. Both
moves were needed, and raising the point ceiling alone would not have been
enough: a build differing from this one *only* in running the copy ledger once
over the run still fails with `data array description resource limit exceeded`
(exit 6) on the Velos run, after 442.7 s, and on the 547 MB, 40,856-spectrum
`centroid_lcms_qe_silac_uk222_picked/UK222_picked.mzML`, after 67.0 s. Every one
of those runs' spectra is far inside the same allowance on its own.

The work ceiling additionally grows with the spectrum, following the mzML
reader's size-derived allowances (`src/format/mzml_scaling.rs`,
[MZML_READER_SCALE_SUPPORT.md](MZML_READER_SCALE_SUPPORT.md)) with the point
count in place of consumed input bytes. A fixed per-spectrum ceiling would not
do: the Velos run's largest spectrum needs 45,766,084 units against the former
50,000,000, a 9% margin. The default rate is eight times the largest ratio in the
table above, rounded up to a power of two — the same margin the mzML reader's
allowances use.

Sliding cost is `Θ(n · w)` for mean window occupancy `w`, and `w` is a property
of the data, not of `n`: peaks spread over about twice the window width put `n/2`
points into each of `n/2` windows, which is quadratic. (Peaks packed into a
*single* window are not: sliding stops at the first window that reaches the last
peak.) A ceiling linear in `n` therefore does two things. Any spectrum at a
realistic occupancy is accepted however many points it has — the measured ratios
above are 271 to 3,493 units per point against a credited 32,768. And a spectrum
that is quadratic in its own size is stopped after work linear in that size
instead of running to completion, so no input can amplify: the total this module
will ever spend on one spectrum is `max_work + work_per_point * n`, and
`max_points` bounds `n` itself.

Work units are an explicit estimate, not exact CPU instructions or a time limit:
preflight costs `n + 1`; each window-end probe costs one; a window with a nonzero
quota costs `4 * window_length` before copying and linear selection; sorting a
list of length `n >= 2` costs `4 * n * ceil(log2(n))` before sorting. Sliding
membership costs `2 * n`. Jumping additionally charges sorting the retained keys
and binary membership searches using `n * (floor(log2(key_count)) + 2)`, or `n`
for no keys. All arithmetic in the accounting is checked. Native container
validation also examines attached metadata and identification records, whose
cost depends on their own size.

The implementation uses linear order-statistic selection within windows rather
than sorting each entire window. Sliding windows can still visit quadratically
many candidates; the work budget covers those visits before allocation/selection.
Jumping uses sorted key lookup instead of the source's quadratic scan for each
original peak. Working storage is linear in input peaks; owned-return and
experiment operations additionally clone original records and annotations
before native selection. No dependencies or custom tree implementation were added.

## Validation and provenance

The port follows pinned OpenMS4-core revision
[`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`](https://github.com/okohlbacher/OpenMS4-core/tree/7c029e8cdba6abab503708ecdd56f6ab55e38ce4),
specifically `PROCESSING/FILTERING/WindowMower.h`, its constructor/dispatch `.cpp`,
and `WindowMower_test.cpp`. Source hashes and methodology are recorded in
[window_mower_provenance.json](../tests/data/window_mower_provenance.json).

Tests reuse the unchanged source `Transformers_tests.dta`: its 121 input peaks
produce the literal source counts 56 sliding and 30 jumping. The source triangle
assertions check exact retained indices and annotations for width 50 and width
10. Further tests cover strict boundaries, gaps, the first-end sliding stop,
final singleton/rounding behavior, duplicates, signed zeros/intensities,
deterministic ties, empty/huge quotas, invalid data and transactional per-record
experiment limits. A deliberately simple independent interpreter uses full stable
sorts and linear source equality membership on 480 small cases to check the
optimized selection implementation.

The per-record, size-derived ceilings above were established against the C++
`SpectraFilterWindowMower` of the pinned Release build
`openms4-release-bc9cc12-c19e494-174b576` on the 1.2 GB LTQ Orbitrap Velos run
`centroid_lcms_velos_pxd001819_50amol_r1/50amol_R1.mzML` and on the 547 MB
Q Exactive run `centroid_lcms_qe_silac_uk222_picked/UK222_picked.mzML`, at one
thread, with the benchmark's own `SpectraFilterWindowMower.ini`; the decoded
outputs of the two implementations were compared spectrum by spectrum, and agree
bitwise on all 87,492 and 81,714 binary arrays respectively. The per-spectrum work and
point figures in the table are from a direct simulation of this module's own
accounting over each benchmark input. The rest of this document's source
semantics were established by reading the pinned C++ only.

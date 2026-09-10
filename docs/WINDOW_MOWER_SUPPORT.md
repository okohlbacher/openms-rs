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
experiment metadata, and shares the point/work limits across processed spectra.

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

`max_points` defaults to 1,000,000 processed input peaks and `max_work` to
50,000,000 work units. Both must be positive. For an experiment, their limits
apply across its spectra, not independently to each spectrum. The total peak
count is checked before processing, and every selection is planned before
cloning spectra or committing any change. Empty records still consume one work
unit. Invalid native array lengths, nonfinite peaks, overflowing coordinate
spans used during window processing, or exhausted limits return errors without
partially changing the input. Invalid configuration is rejected even for empty
input. Zero quotas select nothing without sorting or processing coordinate spans.

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
deterministic ties, empty/huge quotas, invalid data and transactional experiment
limits. A deliberately simple independent interpreter uses full stable sorts and
linear source equality membership on 480 small cases to check the optimized
selection implementation. No C++ program was built or executed.

# Iterative mean noise estimation

`processing::mean_noise::SignalToNoiseEstimatorMeanIterative` estimates noise
and signal-to-noise for ordered spectra, chromatograms or position/intensity
slices. It follows OpenMS4-core revision
`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`, runs three histogram clipping passes
per window, and returns owned estimates without modifying its input.

## API and defaults

The methods are `estimate(positions, intensities)`, `estimate_spectrum(input)`
and `estimate_chromatogram(input)`. `MeanNoiseEstimates` contains per-point noise
and signal-to-noise, the histogram maximum and the percentage of sparse windows.
Sparse-window warnings are represented by diagnostics rather than global logs.

Defaults retain the source: full window length 200, 30 bins, clipping deviation
multiplier 3, at least 10 included samples, and sparse noise 1e20. The default
histogram range is global mean plus three population standard deviations.
`MeanNoiseHistogramRange` also supports a positive manual maximum and the
checked `LegacyPercentile` behavior described below.

Coordinates retain supplied units: Th for spectra and seconds for chromatograms.
Finite signed intensities and coincident positions are accepted. Bin assignment
clamps negative intensities to zero; global statistics and signal-to-noise
numerators retain their original signed values. A negative automatic maximum
returns an error where the source warns and returns zero-filled estimates.

## Source numerical conventions

The sliding interval is left-inclusive and right-exclusive. Bin width is the
larger of 1 and maximum/bin_count. Intensities whose bin index is at least
bin_count are omitted, rather than clamped into the final bin as in the median
estimator. The width floor means the effective upper bound can exceed a small
configured maximum. Only included samples count towards the sparse minimum.

Each non-sparse window starts with all bins. Each of three passes calculates
the weighted midpoint mean and population variance, then sets the exclusive
upper bin to the integer truncation of
`(mean + stddev * stdev_multiplier - 1) / bin_width + 1`, capped at bin_count.

The mean and variance always divide by the **original included window count**.
That denominator does not decrease when high bins are excluded. Each pass can
reconsider previously excluded bins if the threshold grows. Final noise is the
third pass's mean with a floor of 1; it is not recomputed after the final
threshold. A single occupied high bin can disappear after threshold truncation
even with zero variance. This differs from ordinary renormalized sigma clipping.

Division precedes multiplication, preserving source accumulation order.
Nonfinite inputs/results, invalid integer conversions and histogram underflows
are errors. A window whose upper endpoint cannot advance beyond its center at
very large coordinates also returns an error.

## Historical percentile mode

`LegacyPercentile { percentile }` preserves the source's defined finite behavior;
it is not a conventional percentile. The C++ comparator selects the minimum
intensity where its comment claims a maximum. The 100-bin histogram uses source
binary32 division by 100 and subtraction of 1 before bin assignment. Its scan
also stops after at most as many bins as input samples. These behaviors apply
to the native f64 slice API too.

Ten intensities of 100 produce a maximum of 9.5 at percentile 95 because the
scan stops before reaching occupied bin 99. One hundred such values produce
99.5. Heterogeneous data often produces invalid source indices. Rust rejects
out-of-range histogram indices, invalid spacing, empty input and a negative
resulting maximum before unsafe access or arithmetic. Percentile 0 therefore
returns an error where the source can warn and return zero-filled estimates.
This mode does not silently substitute a corrected statistical percentile.

## Limits and validation

Point and bin limits default to 1,000,000, with 50,000,000 work units. Units
cover output points, global intensity scan visits, histogram slots/updates and
clipping mean/variance bin visits. Percentile histogram work is also charged.
Limits precede the corresponding allocations/scans, including convenience
method coordinate/intensity copies. Container validation also checks metadata
and attached arrays. Worst-case numerical work is O(n*b) and storage O(n+b).
Empty default/manual calls return empty estimates with zero sparse percentage;
historical percentile mode requires samples. No C++ or new dependency is used.

[Fixture provenance](../tests/data/mean_noise_provenance.json) records all 2,526
historical source samples and expected outputs. The source test uses window
40.1, sparse noise 2, minimum 10 and absolute output tolerance 0.5.
[Native tests](../tests/mean_noise.rs) additionally check calculated clipping
cases, strict boundaries, histogram exclusion, signed/zero/empty signals,
percentile quirks, container equivalence, input preservation and exact work
limits. No C++ code was built or executed.

An [independent review suite](../tests/mean_noise_review.rs) compares incremental
windows with fresh histogram rescans over 192 configurations, including signed
global statistics and historical-percentile rounding boundaries.

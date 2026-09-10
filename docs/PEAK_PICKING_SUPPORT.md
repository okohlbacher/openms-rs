# High-resolution peak picking

The native `processing::peak_picking` module ports `PeakPickerHiRes`, its natural cubic spline interpolation, the supported histogram-median noise modes, and the spectrum-type heuristic used for automatic experiment selection. The reference is OpenMS4-core commit `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. No C++ library is called or built.

```rust
use openms::processing::peak_picking::{FwhmUnit, PeakPickerHiRes};
use openms::MSSpectrum;

fn centroid(profile: &MSSpectrum) -> openms::Result<MSSpectrum> {
    let picker = PeakPickerHiRes {
        signal_to_noise: 1.0,
        report_fwhm: Some(FwhmUnit::Absolute),
        ..Default::default()
    };
    let result = picker.pick_spectrum(profile)?;
    // result.boundaries follows the centroid order.
    // result.omitted_arrays identifies profile annotations without an aggregation rule.
    Ok(result.spectrum)
}
```

`pick_spectrum`, `pick_chromatogram`, and `pick_experiment` return newly owned results with boundaries and the names of omitted profile arrays. `SpectrumFilter::filter_spectrum` and `filter_experiment`, plus `filter_chromatogram`, replace the input only after the complete operation succeeds. These convenience filters discard the omission report; use the returned-result APIs when annotations matter.

## Algorithm behavior

The default signal-to-noise threshold is **zero**, which disables noise estimation. A core requires a strict intensity maximum with nonzero immediate neighbors. The first and last two input samples cannot be peak cores; fewer than five samples produce an empty centroid record. The source's oscillating-satellite rejection and jump over already considered samples are retained.

Spectrum spacing checks default to factors 1.5 for a missing sample and 4.0 for a hard gap, relative to the smaller immediate apex spacing. Comparisons are strict. Zero disables either constraint; disabling both bypasses spacing checks. Extension proceeds independently to the left and right along non-increasing intensities, stops after zero intensity or a hard gap, and permits one missing sample per side by default. `allow_missing_flank` retains the source's recently added one-sided-core behavior. Chromatograms disable spacing checks by default. Explicit `*_with_spacing` methods expose the source overloads.

Centroids are maxima of a **natural cubic spline**, using the source tridiagonal recurrence and derivative bisection with a coordinate tolerance of `1e-6`. The intensity is evaluated at that spline maximum; it is not an integrated area or an intensity-weighted m/z. `CubicSpline2d` exposes interpolation, first through third derivatives, and the source peak-bracket bisection as reusable fallible helpers. Its `peak_maximum` assumes a peak bracket with a positive derivative to the left; it does not find a global maximum of an arbitrary spline.

Optional FWHM uses spline half-height crossings with a one-percent-of-half-height intensity tolerance. If support ends above half height, the support endpoint becomes the crossing, as in the source; the resulting width can underestimate the full physical peak width. `FwhmUnit::Absolute` writes the `FWHM` float array in Th or seconds. `Ppm` writes `FWHM_ppm` as width divided by centroid position times one million and requires a positive centroid position. Absolute units are normally appropriate for chromatograms.

Boundaries retain the source's **extension boundaries**, including the final missing sample even when that sample was rejected from the spline support. Thus a boundary is not necessarily an interpolation knot or a half-height crossing.

## Noise estimation

`SignalToNoiseEstimatorMedian::estimate` accepts parallel coordinate/intensity slices and returns per-sample noise and signal-to-noise values, the histogram upper range, and sparse/rightmost-bin percentages. It recomputes results for each call; there is no stale parameter cache.

The default histogram range is mean plus three population standard deviations. `NoiseHistogramRange::Manual` supplies a positive range explicitly. Defaults otherwise match the source: window length 200 coordinate units, 30 histogram bins, 10 required samples, and noise `1e20` for sparse windows. The moving window includes both coordinate endpoints. Histogram bins have width at least one, and intensities above the upper range enter the final bin. The lower median rank is `(count+1)/2`; noise is interpolated within that median bin and floored at one. This is the actual implementation at the pinned revision, rather than an exact sorted median or the older bin-center approximation.

The source percentile auto-range branch is not exposed. At this revision it uses a reversed maximum comparator and unchecked histogram indices, so reproducing it would include invalid accesses for ordinary data. Manual and standard-deviation range modes are implemented and tested. Logging is replaced by returned diagnostic percentages. Empty inputs return empty estimates and zero percentages.

## Selection, metadata, and annotations

Single-record picking marks spectra as centroided. It copies all acquisition metadata represented by the Rust kernel: identifiers, name, retention time, MS level, precursor information, and metadata map. Profile float, integer, and string arrays are removed because they have no general mapping to newly interpolated peaks; omitted names are reported.

For spectra, the names `Ion Mobility` and `raw inverse reduced ion mobility array` are recognized automatically. Alternatively, `ion_mobility_array` selects an exact custom float-array name. Its values are intensity-weighted over precisely the included spline support and retain the array name. Duplicate matching arrays, missing explicitly selected arrays, nonfinite values, or misaligned values are errors. This covers the source test cases without pretending to implement its full controlled-vocabulary ion-mobility classification. Mobility reduction uses float64 products; the C++ expression multiplies float32 values before accumulation, so minor rounding differences are expected. Chromatogram profile arrays are omitted.

With an empty `ms_levels` selection, experiments copy spectra already marked centroided and pick the others. Unknown types use the source `PeakTypeEstimator` heuristic: inspect up to five high-intensity maxima, stop after explaining over half the total intensity, and classify shoulder evidence with the source 0.75 threshold. Explicit MS-level selection copies other levels and rejects selected centroid spectra when `check_spectrum_type` is true. All chromatograms are picked. Boundary results retain an entry for every input spectrum, with `None` for a copied spectrum; this removes the source API's ambiguity when it only appends boundaries for picked spectra.

Data-processing history, which can also indicate centroiding in C++, is not represented by the current kernel and cannot influence automatic selection. Mobilograms, on-disc experiments, other peak-picker families, and full two-dimensional ion-mobility processing are outside this increment.

## Validation and limits

Input coordinates must be finite, strictly increasing, and distinct; intensities must be finite and nonnegative. Duplicate coordinates and negative intensities are explicit Rust errors, avoiding undefined or nonconvergent interpolation cases in the source. Baseline-corrected data with negative intensities must be clipped or otherwise handled deliberately before picking.

Picking and spline construction default to at most one million points per record. Peak picking defaults to ten million apex/extension visits per record; noise estimation defaults to one million bins and fifty million histogram updates/median-bin visits. Limits are configurable and checked. Coordinate, spline, histogram, FWHM, and output-float arithmetic failures return errors. Spline maximization and each FWHM crossing are bounded by 128 iterations. Maximization stops if floating-point coordinates cannot move; FWHM returns an error if it cannot reach the required intensity tolerance. Limits are per record, and the experiment APIs retain the complete owned experiment in memory.

These checks intentionally define failure for malformed or numerically unrepresentable input instead of promising identical C++ failure behavior. This is algorithm-level source parity on validated inputs, not full API or binary compatibility.

## Verification and provenance

The fixture manifest `tests/data/peak_picking_provenance.json` records upstream and derived-file SHA-256 values. Small test fixtures contain the first spectrum's binary arrays from the original Orbitrap and FTMS mzML inputs and S/N 1 and 4 reference outputs. They were decoded with Python's standard XML, base64, zlib, and struct libraries, without filtering or numerical transformation. The signal-to-noise DTA input/output files are copied unchanged.

Tests compare exact peak counts and centroid positions/intensities for Orbitrap (45 and 18 peaks) and FTMS (14 and 4 peaks), all reference noise values, published cubic-spline values/derivatives, source boundary literals, and source mobility/missing-flank examples. Additional tests cover histogram interpolation, sparse windows, zero input, FWHM, selection, chromatograms, mutation atomicity, malformed inputs, and resource limits. Only first-spectrum comparisons are claimed for those instrument fixtures; the larger simulation and other input spectra have not been exhaustively validated.

Source links at the pinned revision:

- [PeakPickerHiRes implementation](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/PROCESSING/CENTROIDING/PeakPickerHiRes.cpp) and [tests](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/PeakPickerHiRes_test.cpp).
- [SignalToNoiseEstimatorMedian](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h) and [tests](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/SignalToNoiseEstimatorMedian_test.cpp).
- [CubicSpline2d implementation](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/MATH/MISC/CubicSpline2d.cpp), [derivative bisection](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/include/OpenMS/MATH/MISC/SplineBisection.h), and [spline tests](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/CubicSpline2d_test.cpp).
- [PeakTypeEstimator heuristic](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/include/OpenMS/FORMAT/PeakTypeEstimator.h).

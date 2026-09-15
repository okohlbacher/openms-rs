# Gaussian, Savitzky–Golay, and morphological processing

The native Rust port implements three additional OpenMS signal-processing families in `processing::smoothing` and `processing::baseline`. Their numerical operations are derived from the pinned C++ source, with source-derived golden tests. They operate on Rust-owned spectra and chromatograms without a C++ runtime or additional numerical dependency.

```rust
use openms::processing::SpectrumFilter;
use openms::processing::smoothing::{GaussFilter, GaussianWidth, SavitzkyGolayFilter};
use openms::processing::baseline::{MorphologicalFilter, MorphologicalMethod, StructuringElement};
use openms::{MSSpectrum, Result};

fn process(spectrum: &mut MSSpectrum) -> Result<()> {
    GaussFilter::new(GaussianWidth::Absolute(0.2))?.filter_spectrum(spectrum)?;
    SavitzkyGolayFilter::new(11, 4)?.filter_spectrum(spectrum)?;
    MorphologicalFilter::new(
        MorphologicalMethod::TopHat,
        StructuringElement::Thomson(3.0),
    )?.filter_spectrum(spectrum)?;
    Ok(())
}
```

The example shows API composition; choose smoothing and baseline parameters for the acquisition resolution and peak width of the actual data. Running both smoothers is optional.

## Gaussian convolution

`GaussFilterAlgorithm::filter` accepts parallel f64 coordinate/intensity slices and returns `GaussianOutput { intensities, found_signal }`. The low-level default width is 0.8 and lookup-table spacing is 0.01, matching `GaussFilterAlgorithm`. `GaussFilter` wraps it with the C++ wrapper's default width of 0.2. `GaussFilter::new` accepts `GaussianWidth::Absolute(width)` or `GaussianWidth::Ppm(ppm)`; its public `algorithm` field allows the lookup spacing and coefficient limit to be configured.

The numerical implementation preserves these details from the C++ code:

- Width is **eight standard deviations**: `sigma = width / 8`. The legacy wrapper's parameter description mentions FWHM, but the actual implementation uses the eight-sigma definition.
- The one-sided Gaussian table contains `ceil(4 * sigma / kernel_spacing) + 1` coefficients. Values between entries use linear interpolation; the final table entry extends through its last bin.
- Each output integrates neighboring signal intervals with the trapezoidal rule, weights them with the interpolated Gaussian, and divides by the integrated kernel area.
- Integration uses strict boundary comparisons and clips its bounds to the first/last coordinate. Intervals whose outside endpoint exactly equals the clipped boundary are excluded. This unusual endpoint behavior is retained and tested.
- Nonpositive integrated numerators become zero. A one- or two-point input consequently produces zeros in the low-level implementation.
- Ppm width is recalculated at every m/z as `(ppm / 1e6) * mz`. Positive m/z values are required. Chromatograms reject ppm because their coordinates are time.

`GaussFilter` also retains the wrapper's important behavior: if every computed intensity is zero and the record contains at least three points, original intensities are preserved. For shorter records the zeros are written. Successful spectrum filtering marks the spectrum as `Profile`, including empty spectra and preserved all-zero-result cases. `filter_chromatogram` changes only intensities. No automatic logging is emitted; the low-level `found_signal` flag is available to callers that need diagnostics.

The algorithm supports nonuniform sample positions because it integrates coordinate intervals rather than averaging a fixed count of samples. Positions must be finite and sorted in nondecreasing order. Duplicate coordinates are accepted; degenerate integrations with no positive area return zero. Width and kernel spacing must be finite and strictly positive, and Gaussian variance/amplitude must be representable in f64. The default limit is one million one-sided coefficients; exceeding the configured limit returns an error before changing data.

## Savitzky–Golay filtering

`SavitzkyGolayFilter::new(frame_length, polynomial_order)` computes and stores reusable coefficients. Its default is frame length 11 and degree 4. Even lengths increase by one, including a supplied length of zero becoming one as in the C++ implementation. Degree must be smaller than the resulting odd frame length.

The filter fits the same polynomial least-squares problem as OpenMS. Interior points use a centered full window. The first half-window uses the first complete window evaluated at successive positions; the final half-window uses the reflected edge coefficients from the last complete window. Negative fitted values are clamped to zero. If the entire signal is shorter than the frame, its intensities remain unchanged, including signed values.

The Rust coefficient calculation uses twice-orthogonalized QR on a scaled Vandermonde basis and constructs the least-squares projection. The pinned C++ implementation uses Eigen's SVD. These solve the same full-rank least-squares problem; the Rust code does not use moving averages or normal-equation inversion as a substitute. Floating-point coefficients are not promised to be bit-identical to Eigen. Source-derived numeric tests and polynomial-reproduction tests cover both interior and asymmetric edge behavior.

Coefficient generation currently accepts frames up to 1023 and degrees up to 32. Numerically rank-deficient fits are rejected with a request to reduce degree. These explicit bounds limit work and avoid returning an unreliable high-degree fit. They are Rust implementation limits, not limits declared by OpenMS C++.

The fit uses **sample indices**, not physical coordinates. Like the source algorithm, it expects uniformly sampled profile data. Finite, sorted nonuniform or repeated coordinates are accepted and preserved, but their intensities still receive an index-space polynomial fit; the filter neither resamples data nor performs a polynomial regression against the actual m/z/time values. Spectrum representation flags remain unchanged.

## Morphological filtering

`MorphologicalFilter` supports the following operations with a flat, centered structuring element and windows clipped to the signal boundaries (with the one end-of-signal exception below):

| Rust method | C++ parameter | Operation |
| --- | --- | --- |
| `Identity` | `identity` | Preserve intensities |
| `Erosion` / `ErosionSimple` | `erosion` / `erosion_simple` | Local minimum |
| `Dilation` / `DilationSimple` | `dilation` / `dilation_simple` | Local maximum |
| `Opening` | `opening` | Dilation after erosion |
| `Closing` | `closing` | Erosion after dilation |
| `Gradient` | `gradient` | Dilation minus erosion |
| `TopHat` | `tophat` | Input minus opening |
| `BottomHat` | `bothat` | **Input minus closing**, retaining the source's signed convention |

Min/max operations use a monotonic deque with O(N) work, replacing the C++ van Herk helper's prefix and suffix blocks while preserving the same mathematical min/max operation and boundary windows. The tests compare all operations against the upstream golden table and min/max results against the source test's direct-window reference across increasing, decreasing, and irregular signals and many window lengths.

One end-of-signal case is not a clipped window, and it is reproduced deliberately: when the structuring element is **one sample** wide and the signal has **more than five** samples, the source's van Herk erosion and dilation never write the last output sample, which therefore keeps the zero its output buffer was allocated with. `erosion`, `dilation`, `opening` and `closing` zero the last sample; `tophat` and `bothat` keep it and zero everything else; `gradient` subtracts a value the source's process-wide scratch buffer carries over from an earlier spectrum, which `filter_experiment` reproduces by carrying one buffer through the spectra in order. The `Simple` variants and `identity` are unaffected. `BaselineFilter` reaches this whenever the element in Thomson is narrower than the peak spacing, which is the normal case for centroided MS2 spectra. Executed C++ evidence, the exhaustive sweep behind "and nowhere else", and the remaining native differences are in [MORPHOLOGICAL_FILTER_SUPPORT.md](MORPHOLOGICAL_FILTER_SUPPORT.md).

The default is top-hat with `StructuringElement::Thomson(3.0)`. For spectra:

- `DataPoints(n)` chooses an integer sample count and rounds even counts upward to odd.
- `Thomson(width)` computes `ceil(width * (N - 1) / (last_mz - first_mz))`, then rounds upward to odd.

Coordinate-width conversion uses **global average spacing**, including for nonuniformly spaced input, exactly as the source wrapper does. It does not choose a different physical-width window at every point. Widths must be finite and positive; sample count zero is rejected. Coordinate-based conversion requires a positive finite total span for records of two or more points. Very wide sample windows cover the entire available signal without allocating storage proportional to the requested window length.

`filter_range` operates directly on f32 intensities and requires `DataPoints` units. `filter_spectrum` preserves the source wrapper's short-record rule: empty/singleton spectra keep their intensities and become `Profile`. Consequently, a singleton top-hat range yields zero while a singleton top-hat spectrum remains unchanged. This distinction is tested.

`filter_chromatogram` is a Rust convenience extension; the C++ morphological class exposes spectrum and spectrum-experiment methods. Chromatograms accept `DataPoints` or `Seconds(width)`, with the same average-spacing conversion and short-record convention. `Thomson` is rejected for chromatograms and `Seconds` is rejected for spectra.

## Mutation and experiment behavior

All methods validate finite input, coordinate ordering, and existing data-array lengths before mutating a record. Output intensities must fit finite f32 values for kernel containers; subtraction/convolution overflow is an error. Coordinates, record names, precursor fields, metadata, and aligned auxiliary arrays remain intact. Gaussian and morphology update the spectrum representation flag to `Profile`; Savitzky–Golay preserves it.

The `SpectrumFilter::filter_experiment` implementations are atomic: failure leaves the original experiment unchanged. Gaussian and Savitzky–Golay process **both spectra and chromatograms**, matching their C++ wrappers. Morphology processes **spectra only**, matching its source wrapper, and leaves stored chromatograms unchanged. Invoke its chromatogram method explicitly when needed.

The Rust APIs replace cached string parameters, iterator output arguments, mutable global scratch buffers, progress loggers, and C++ exception behavior with typed configuration, returned values, and `Result`. They do not implement mobilogram-specific wrappers because this initial kernel does not yet expose that type.

## Verification and sources

The targeted suites currently contain 10 smoothing and 8 baseline tests. Gaussian outputs are checked against the nine fixed assertions in `GaussFilterAlgorithm_test.cpp` and the constant-signal assertion. Savitzky–Golay outputs are checked against the source's even-frame quadratic example, with independent polynomial-reproduction and boundary cases. All ten morphology operations are checked against the complete 40-row `MorphologicalFilter_test_1.txt` fixture.

Additional checks cover nonuniform Gaussian intervals, ppm equivalence at each evaluated position, all-zero preservation, short records, metadata/array preservation, failed-input atomicity, parameter/overflow errors, global-spacing conversion, signed bottom-hat behavior, and full experiment handling. Gaussian and Savitzky–Golay golden tolerances account for rounded C++ assertions and f32 storage. The smoothing tests are source-derived reference tests: no C++ build was run for them. Morphology additionally carries executed differential evidence — the OpenMS4 Release build's own `MorphologicalFilter` and `BaselineFilter` over edge shapes, an exhaustive element-length sweep and the benchmark's real runs, in `tests/baseline_filter_edges.rs` and `tests/topp_baseline_filter_edges.rs`, recorded in [baseline_filter_edges_provenance.json](../tests/data/baseline_filter_edges_provenance.json) and [MORPHOLOGICAL_FILTER_SUPPORT.md](MORPHOLOGICAL_FILTER_SUPPORT.md).

Fixture provenance, derivation, and SHA-256 hashes are recorded in [smoothing_provenance.json](../tests/data/smoothing_provenance.json) and [baseline_provenance.json](../tests/data/baseline_provenance.json). The morphology table is copied byte-for-byte; smoothing tables transcribe the pinned source's inputs and expected numbers.

All source links refer to OpenMS4-core revision `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`:

- [GaussFilterAlgorithm header](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/include/OpenMS/PROCESSING/SMOOTHING/GaussFilterAlgorithm.h) and [kernel initialization](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/PROCESSING/SMOOTHING/GaussFilterAlgorithm.cpp).
- [GaussFilter wrapper](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/PROCESSING/SMOOTHING/GaussFilter.cpp).
- [SavitzkyGolayFilter application](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/include/OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h) and [coefficient generation](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/PROCESSING/SMOOTHING/SavitzkyGolayFilter.cpp).
- [MorphologicalFilter operations and unit conversion](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/include/OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h).

# EMG reference review

This review covers `EmgGradientDescent` and its use by `PeakIntegrator` at OpenMS4-core revision `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. It uses the immutable source checkout and literal source test assertions. No C++ library was built or executed. The broader native port remains in progress.

The reference files are [the C++ implementation](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/MATH/MISC/EmgGradientDescent.cpp), [its header](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/include/OpenMS/MATH/MISC/EmgGradientDescent.h), and [its tests](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/EmgGradientDescent_test.cpp). Their SHA-256 hashes, the constants and integration sources, and fixture hashes are recorded in [emg_provenance.json](../tests/data/emg_provenance.json).

## Literal scientific references

[emg_traces.tsv](../tests/data/emg_traces.tsv) preserves 429 samples in seven traces: glutamate (107), saturated peaks in minutes and seconds (83 each), saturated cutoff peaks in minutes and seconds (66 each), and cutoff peaks in minutes and seconds (12 each). Seconds are retained as the original decimal literals, including their rounding; they are not reconstructed by multiplying the minute coordinates. Each row includes the original decimal positions and intensities, their binary64 bits, and the binary32 intensity bits used in source spectra and chromatograms.

[emg_fits.tsv](../tests/data/emg_fits.tsv) preserves the source output counts of 107, 87, 71 and 28 and the seven parameter assertions. Source fit parameters were transported through a float data array, so tests round the Rust parameters to binary32 before comparing them. The native result retains parameters as binary64 in a typed estimate. Parameter comparisons use relative tolerance `1e-5`, including widths below one; there is no unit-sized absolute tolerance floor. Both spectrum and chromatogram overloads are exercised.

The four source loss assertions need a separate calculation: they evaluate the full original binary64 intensity trace with fitted parameters rounded to binary32. They are not the optimizer's selected-training-point loss, nor the loss on the stored binary32 intensities. The independent test reconstructs that exact convention.

The source's fixed-parameter saturated-cutoff example uses `h=15515900`, `mu=14.3453`, `sigma=0.0344277`, and `tau=0.188507`. It produces 66 samples without extension and 71 with extension. The extended first point is approximately `(14.2717555076923, 108845.941990663)`. Scalar source assertions exercise positive, negative and asymptotic `z` branches in both coordinate units.

## Selection and optimizer review

Training selection retains its collection order: ascending points on the left, the final point followed by descending points on the right, then additional interior points selected from either side using derivatives. Sorting that selected set would change the floating-point accumulation order. [emg_training.tsv](../tests/data/emg_training.tsv) records ordered indices, initial means and initial losses from an independent translation of the source rules. The tests compare the selected loss in that exact order. The source test directly asserts training counts of 107, 77, 61 and 12 and four initial means; the additional ordered-index values are derived references, not upstream assertions.

The review also checks the source's `i == 1` behavior: if the second sample already exceeds the 80% intensity threshold, the derivative loop is skipped. A plateau case confirms that only the two endpoints enter training. Both endpoints are always selected, so the native sorted-input span used for the mean constraint equals the source training-set span.

Initialization depends on absolute coordinates: `sigma = initial_mu * 0.01`, `tau = sigma * 2`. No shift, normalization or guaranteed minute/second invariance is introduced. After updates, amplitude is bounded below by the observed maximum, mean by its initial value plus or minus 35% of the training span, sigma by `[1e-4, 20]`, and tau by `[sigma, 15*sigma]`. The best earlier loss determines returned parameters. The source samples loss every 50 iterations into a ten-value history initialized to zero and stops when its population standard deviation is below one. Its zero-gradient iRprop update still subtracts the current learning rate; that behavior is retained and covered by the implementation's private tests.

## Independent gradient and special-function checks

All twelve analytic gradient branch expressions were compared as arithmetic syntax trees against the C++ source, normalizing only function qualification and integer versus floating literal spellings. Their expression structure agrees. Four centered finite-difference checks per dataset independently compare the amplitude, mean, sigma and tau derivatives with an independently transcribed scalar model. Positive-`z`, negative-`z`, mixed-`z` and asymptotic datasets are recorded in [emg_gradient_review.json](../tests/data/emg_gradient_review.json), including parameters and difference steps.

The normal-branch checks differ by at most approximately `3.4e-9` relative. The asymptotic test has very small derivatives; its tau derivative is approximately `-3.98e-11`, with a cancellation-limited relative difference of `8.44e-6`. The other asymptotic derivatives agree within approximately `1.16e-8`. No gradient transcription defect was found. These are independent numerical review results, not claims of a C++ execution comparison or exhaustive conditioning guarantees.

Additional literal scalar oracles use Python's standard-library `math.erfc`, independently of the Rust `libm::erfc` implementation. They cover the `z=0` branch boundary, a small left-tail value near `6e-149`, a right-tail value near `4e-304`, subnormal output near `1.7e-321`, and ordinary underflow to zero. Normal finite cases use relative tolerances; the subnormal case permits two representable binary64 steps.

The source's positive-`z` expression multiplies an exponential by `erfc` directly. It can overflow well before the stated asymptotic threshold `6.71e7`; this port preserves the expression and reports its nonfinite result instead of substituting a scaled complementary error function. Tests distinguish that failure from the zero returned by the source's asymptotic branch. Using `libm` does not remove the overflow in surrounding exponential expressions.

## Checked boundaries and integration

Application preserves the source's average adjacent sampling difference, sampled apex, endpoint-based choice of extension side, `1e-3` stopping threshold and threefold positional bound. Added positions advance recursively by that sampling difference. The left-side implementation accumulates and reverses extra points, preserving values without repeated insertion at the front of a vector.

The independent tests exercise inclusive point and evaluation limits at exact success/failure boundaries, including the combined estimate-plus-application budget. The implementation validates finite inputs and parameters, strictly increasing positions, positive widths, finite numerical output and representable binary32 container intensities. Invalid calls return errors without mutating their input containers.

One deliberate difference is strict failure after numerical breakdown. The C++ optimizer can retain an earlier finite best when a later loss or parameter becomes nonfinite. The Rust call returns an error rather than returning that partial fit. For example, positions `[100,101,102,10000]` with intensities `[1,10,1,1]` have a finite initial model at `(h,mu,sigma,tau)=(10,101,1.01,2.02)`, but the far-right amplitude-gradient expression overflows. The independent regression confirms that initial evaluation succeeds and fitting fails. Nonpositive coordinate-derived initial widths and zero iteration limits are also checked errors.

The source attaches a four-element `emg_parameters` float array regardless of peak count. The native API returns a typed estimate and reports omitted input arrays, preserving aligned-array invariants. Optional `Some(0.0)` bounds mean a literal zero boundary; the C++ wrapper uses zero as an omitted-bound sentinel.

The separate [integration tests](../tests/emg_integration.rs) exercise explicit fit-then-integrate against `PeakIntegrator`'s optional EMG path, including both left and right extension, both container types, integration methods, background methods and shape metrics. Those tests use the full fitted span and preserve supplied apex and height arguments. They are workflow equivalence checks, not invented upstream integrated-area goldens.

## Validation

The independent [emg_reference.rs](../tests/emg_reference.rs) suite contains eight tests. It is checked with the current Rust toolchain and the declared Rust 1.85 minimum, with and without optional format features as covered by the repository validation matrix. The fixture extraction, gradient review and Rust tests do not execute the C++ implementation. Overall supported behavior and remaining limits are documented in [EMG_SUPPORT.md](EMG_SUPPORT.md).

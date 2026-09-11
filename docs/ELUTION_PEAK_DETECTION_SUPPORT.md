# ElutionPeakDetection

`analysis::elution_peak_detection` implements the complete class-specific public
scientific operation group from SDK `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`.
It reuses native `MassTrace`, `SavitzkyGolayFilter` and `ProgressLogger`, with no
new dependency. C++ inheritance, optional OpenMP scheduling and binary layout
are not native API contracts.

| Source operation | Native operation |
| --- | --- |
| Construction and six parameters | `ElutionPeakDetection::new/default`, public `options`, `limits`, `logger` |
| `detectPeaks(MassTrace&, ...)` | `detect_peaks`, `detect_peaks_into` |
| `detectPeaks(vector<MassTrace>&, ...)` | `detect_peaks_many`, `detect_peaks_many_into` |
| `filterByPeakWidth` | `filter_by_peak_width`, `filter_by_peak_width_into` |
| `findLocalExtrema` | `find_local_extrema`, returning owned `ElutionExtrema { maxima, minima }` |
| `smoothData` | `smooth_data(trace, signed_window)` |
| `computeMassTraceNoise` | `compute_mass_trace_noise` |
| `computeMassTraceSNR` | `compute_mass_trace_snr` |
| `computeApexSNR` | `compute_apex_snr` |

The six defaults are chrom_fwhm 5, chrom_peak_snr 3, width_filtering Fixed,
min_fwhm 1, max_fwhm 60, and masstrace_snr_filtering false. The width enum also
supports Off and Auto. A typed options value replaces stringly typed parameter
mutation; no implicit generic `Param` synchronization is claimed. Finite negative
settings retain defined source branches, and unrelated settings are not validated
by helpers that do not use them.

## Scientific behavior

Smoothing uses a sample-index quadratic Savitzky–Golay fit; RT must be
nondecreasing on that path, but spacing need not be uniform. Windows are at least
three and even values are incremented. Fewer than three peaks copy raw
intensities and ignore the window. A frame longer than the trace also leaves the
raw values unchanged; only fitted values are clamped to zero. The existing
smoother's odd-frame ceiling is **1023**, including frames longer than the trace.
The native QR projection replaces source SVD with the existing tested smoother;
bitwise coefficient equivalence is not claimed. Each result is narrowed to f32,
then promoted to f64 smoothing, matching the source's temporary Peak1D storage.

Detection computes ceil(chrom_fwhm / average MS1 cycle time) before smoothing.
It rejects nonfinite/negative/out-of-i32-range conversions; the source's
floating-to-Size or subsequent signed narrowing is undefined or platform
sensitive there. Empty, singleton and zero-span traces therefore need not be
valid detection inputs although direct smoothing accepts short traces. A finite
zero window uses smoothing frame three and extrema neighborhood zero.

Extrema preserve ascending-intensity seed order, index ties, positive seed
selection, used-window suppression, and the source neighborhood's **exclusive
right endpoint**. Short traces return their first maximum even when nonpositive.
The source bisection is retained rather than replaced with a global valley
minimum. Tied final valleys choose the right boundary; valley intensity is
clamped to at least one. Splits require both height ratios at least two and both
RT distances at least min_fwhm/2. This uses min_fwhm even when width filtering is
Off or Auto, despite the upstream parameter-description wording. Rejecting a
valley does not remove entries from the returned maxima vector.

Successful detection publishes source input state changes: smoothing always,
and single-maximum fixed-filter FWHM changes even on rejection. An accepted
single trace updates smoothed-apex RT and FWHM while preserving its prior m/z
centroid/SD and label. Multiple-maxima splitting creates nonoverlapping segments;
accepted segments update apex RT, weighted m/z and SD, inherit quantification,
IM centroid when present, and average m/z/IM widths. Their numeric label suffix
uses the original segment ordinal, including gaps for rejected segments.
Multiple maxima without a split valley still takes the `.1` segment branch.

The split SNR test deliberately uses the **whole original trace**, as the source
does. Width and SNR tests are both evaluated when enabled, even if width already
rejects. Auto does not invoke quantile filtering during detection: the caller
must explicitly call `filter_by_peak_width`.

Quantile filtering recalculates every input's smoothed FWHM, including rejected
traces, and returns stable width order. Inclusive ranks floor(n*0.05) through
floor(n*0.95) retain 19 of 20 traces. Empty input succeeds. This operation is
independent of the configured width-filtering mode.

Noise is source-order RMSE of raw f32-promoted minus smoothed f64 intensity.
Empty smoothing gives zero. Mass-trace SNR divides raw trapezoidal area by noise
and absolute trace length; an empty trace gives zero, but nonfinite division
results on nonempty input are checked errors. Apex SNR uses the smoothed maximum
only for positive noise, and otherwise returns zero. Existing MassTrace finite
arithmetic and FWHM checks remain effective; unused cached/scalar fields are not
normalized.

## Ownership, limits and progress

Every detection/width-filter call stages the entire input and output before
publication. Errors preserve all caller input traces and output, including late
batch failures. Successful rejection still publishes the source mutations.
`*_into` methods return the old output's ownership without scanning, cloning,
validating or destroying it. New segments inherit the input's native MassTrace
limits. Ordinary Clone/equality/caller destruction have normal Rust costs.

The default operation limits are one million input/generated traces, ten million
input/generated peaks, 50 million weighted work units and 256 MiB cumulative
allocated payload. Input and generated counts are each capped; staging and
copies share bytes/work. Work covers descriptor scans, complete trace copies,
SG coefficient construction and convolution, vector initialization, sorting,
actual neighborhood visits/marking/bisection, all segments including discarded
ones, labels, and conservatively charged MassTrace summary calls. These charges
precede calls to reused helpers; their individual allowances never reset the
operation's allowance. SG coefficient work is precharged conservatively by frame
squared. Count limits are ceilings, not promises that every count below them fits
work/byte allowances. Checked arithmetic and fallible reservations guard owned
buffers; the existing bounded SG and ordinary Clone paths retain their normal
Rust allocator behavior.

Default logging is CMD. Single detection emits no progress events. Batch
execution emits start, set before processing each trace, and end, using the
existing clock suppression and caller-replaceable backend. An error after start
attempts end to balance nesting; logger state/output cannot be rolled back.
Scientific publication also waits for successful end. Native serial trace/segment
order replaces source optional OpenMP interleaving.

## Evidence

`tests/elution_peak_detection.rs` uses an exact scalar projection of the source
333-scan fixture, generated by `tools/generate_elution_peak_detection_fixture.py`.
Default MassTraceDetection first yields T1; detection then verifies three split
traces, source labels, two maxima/minima count pairs, noise 573.8585 and all six
published SNR values. The source explicitly uses relative tolerance 1.01; native
rounded-literal comparisons allow 1% relative error. The third label T1.3 is a
control-flow-derived assertion, distinguished from active upstream assertions.

Independent tests include a closed-form RMSE/area/apex oracle, all 4^5 small
signals across six neighborhoods, equal-height/valley and threshold cases,
f32 smoothing storage, inclusive quantiles, successful rejection state,
original-trace split SNR, zero-window label gaps, old output identity, late
rollback, and experimentally located byte/work boundaries proving batch
allowances do not reset per trace. Current/minimum regression runs also include
existing MassTrace, MassTraceDetection and smoothing suites.

Source and fixture hashes are in
`tests/data/elution_peak_detection_provenance.json`. No C++ build or execution
is used to generate numerical expectations.

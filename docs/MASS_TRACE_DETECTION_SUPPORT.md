# MassTraceDetection

`analysis::mass_trace_detection` implements the complete computation surface of
`FEATUREFINDER/MassTraceDetection` from OpenMS revision
`82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. It uses the existing owned `MassTrace`
and borrowed experiment/area models. No C++ build, runtime, new dependency, or
mobility-field approximation is involved.

## Native API

| Source operation | Native operation |
|---|---|
| Construction and parameter defaults | `MassTraceDetection::new()` / `Default`, public `options` |
| Both termination enum values | `TraceTerminationCriterion::{Outlier, SampleRate}` |
| `run(PeakMap, output, max_traces)` | `run(input, max_traces)` or `run_into(input, output, max_traces)` |
| `run(ConstAreaIterator&, ConstAreaIterator&, output)` | `run_area(begin, end, output)` |
| `hasFwhmMz`, `hasFwhmIm`, `hasCentroidIm` | `has_fwhm_mz`, `has_fwhm_im`, `has_centroid_im` |
| Inherited progress selection/backend | Public `logger: ProgressLogger`; default mode is `Cmd` |
| CCS-tolerance warning | `ccs_tolerance_warning()` on the last successful run |

The eleven parameters are represented directly, with no string conversion or
silent normalization:

| Field | Default |
|---|---:|
| `mass_error_ppm` | 20 |
| `noise_threshold_int` | 10 |
| `chrom_peak_snr` | 3 |
| `ion_mobility_tolerance` | 0.01 |
| `reestimate_mt_sd` | true |
| `quant_method` | `MassTraceQuantMethod::Area` |
| `trace_termination_criterion` | `Outlier` |
| `trace_termination_outliers` | 5 |
| `min_sample_rate` | 0.5 |
| `min_trace_length` | 5 |
| `max_trace_length` | -1 |

A negative maximum length disables that check; other finite negative scalar
settings retain their arithmetic meaning. The outlier count is a native `usize`;
source signed-to-unsigned parameter conversion can be represented explicitly by
an unsigned value, without an implicit narrowing/wrapping parser. Enum names,
parameter-description trees and `DefaultParamHandler` inheritance are represented
by the typed options, rather than a second configuration state.

`run` returns a fresh vector. `run_into` atomically replaces a caller's output and
returns the previous vector by ownership. `run_area` returns `None` for equal
endpoints, or `Some(previous_output)` after successful replacement. None of these
operations traverses, validates, clones, or destroys the previous output.

## Source computation

Only MS1 spectra enter numerical detection. All MS1 scans count toward the
minimum of three scans, including empty scans and scans emptied by noise
filtering. Retained peaks satisfy `intensity > noise_threshold_int`; apices also
satisfy the strict `intensity > chrom_peak_snr * noise_threshold_int` comparison.
Apices are visited by descending intensity, reversing encounter order on equal
intensities. Zero `max_traces` means unlimited, subject to resource limits.

Each apex starts a trace. Incremental m/z and IM means count that apex twice,
exactly as the source initializes its counters. Each extension iteration tries
the preceding scan before the following scan, so both directions share the
updated mean and standard deviation. Dynamic deviation follows the source's
logarithm/exponential expression and machine-epsilon replacement condition.

Without IM, the nearest m/z wins, with the predecessor winning a midpoint tie.
With IM, inclusive m/z and mobility windows are applied; the nearest **m/z**
inside that window wins, despite the source comment referring to nearest IM.
Equal IM-window m/z distances retain the first candidate. A visited nearest peak
is rejected without searching for a second candidate. Failed traces do not mark
any peak visited. Accepted traces claim all their peaks before subsequent apices.

In `Outlier` mode the counter must be **greater than** the configured limit.
Empty scans move the scan index and increase the scan count, but do not increase
or reset consecutive misses. In `SampleRate` mode a direction is eligible to stop
only after more than five scanned positions; its decision uses the combined
hits and scans of both directions at that exact point in the down/up iteration.
Final quality divides trace length by all scanned positions minus only the
trailing nonempty-scan misses. Empty scans remain in that denominator.

Accepted traces use input scan order, the source right-rectangle weighted RT,
ordinary intensity-weighted m/z, weighted m/z deviation, optional median FWHM
values, and the incremental IM centroid. Labels are `T1`, `T2`, and so on. No
implicit smoothing, FWHM estimation, chromatogram conversion or RT sorting occurs.

### Mobility and width arrays

The first MS1 scan with a nonempty float-array vector determines the three exact
names: `FWHM_ppm`, `Ion Mobility`, and `IM Peak FWHM`. The discovered indices and
names must match across every MS1 scan. A first vector containing only unrelated
arrays prevents discovery of these names in later scans. Nonempty arrays of all
three scalar types must align with their original peaks, matching source
`MSSpectrum::select`; empty arrays remain empty during filtering.

The algorithm reads only retained indices from recognized arrays, preserving
f32-to-f64 conversion. All other array values, descriptions, acquisition fields,
processing handles, experiment metadata and chromatograms are borrowed and
unmodified. No irrelevant metadata graph is cloned.

The CCS warning is separate from scientific IM selection. Its scan includes all
MS levels and stops at the first array recognized by the pinned ontology/vendor
name rules. The nine proper descendants of `MS:1002893` have millisecond or VSSC
units; none has CCS units. CCS can be recognized by the source vendor fallback,
such as `Ion Mobility MS:1002954`. The generic parent `ion mobility array` is not
its own descendant. Prefix and accession precedence follow `IMDataArrayUtils`.
The native diagnostic getter replaces a global warning-stream side effect; it
never rescales or changes any scientific values.

### Area overload

Equal endpoint iterators return immediately, without inspecting options,
limits, prior output, or flags. Otherwise the source overload reconstructs plain
MS1 scans from peak/RT pairs, discarding source annotations and array metadata.
Exactly equal consecutive f64 RT values are grouped. A preceding RT `-1` group
is discarded when the RT changes, while a final `-1` group is retained, preserving
the source temporary-spectrum sentinel behavior. The ordinary three-scan rule
then applies. This API uses the existing scalar RT/m/z area iterator; it does not
add source scan-mobility filtering or claim that unrepresented overload.

## Checked boundaries and ownership

* Retained m/z coordinates must be sorted within each MS1 scan for source binary
  searching. The detector does not require or impose RT ordering; finite input
  order is preserved. Equal-RT traces can later fail the checked zero-area mean.
* All scalar options, consumed RT/m/z coordinates, MS1 intensities, recognized
  consumed array values, and resulting arithmetic must be finite. Logarithms of
  zero used through `exp(-inf) = 0` remain valid. A nonfinite incremental ratio,
  variance result, zero-area RT mean, or invalid total weight returns an error.
  Thus a finite negative setting can be accepted in one branch and fail later
  if a consumed arithmetic expression becomes nonfinite.
* Non-MS1 peaks and unrelated payload values are not validated. Their float-array
  names may still be inspected for the independent CCS diagnostic.
* Availability flags reset on every **successful** run. Source flags can otherwise
  survive a previous IM run and index missing arrays on reuse; this native
  correction makes mixed-array reuse safe. On failure, prior flags and the
  previous CCS diagnostic are preserved.
* Failure preserves output and, for `run_area`, the caller's iterator. A live
  endpoint that cannot be reached is an error, replacing source end dereference.
* Progress uses the existing logger's start, throttled set, and end operations.
  Successfully started progress is ended even after a computation error to
  balance nesting. Progress/backend side effects cannot be rolled back; an end
  error prevents publication of scientific output and state.

Default configurable `MassTraceDetectionLimits` allow one million input spectra,
ten million MS1 input peaks (or total visited area peaks), one million output
traces, fifty million charged operations, and 256 MiB of cumulative new logical
allocation. Counts are checked before traversal/allocation. The same allowance
covers area reconstruction, metadata descriptor/name visits, filtering, binary
searches, candidate visits and rejected traces, sorting, vector initialization
and growth copies, output creation, and conservative precharges for existing
`MassTrace` summary methods. Those methods retain their own checked per-call
limits; the shared precharge prevents a new allowance per accepted trace.
Allocated bytes are cumulative and are not refunded when temporary buffers are
released. Caller-owned input payload and unrelated shared handles are not copied
or charged as new allocation. Allocator bookkeeping and logger/caller callbacks
are outside the logical scientific allowance. Ordinary Rust `Clone`, indexing,
configuration assignment, and caller destruction retain ordinary Rust costs.

## Evidence

`tests/mass_trace_detection.rs` uses a lossless scalar projection of the upstream
133-scan, 1539-row mzML fixture, generated by
`tools/generate_mass_trace_detection_fixture.py` with standard-library XML,
base64, zlib and explicit little-endian decoding. Input intensities are explicitly
narrowed to source f32 before projection. It tests the source default two traces,
three traces after reducing the minimum length, exact lengths 86/31/16, the
source RT/m/z literals, all three precise area literals, and ignored MS2 scans.
Explicit tolerances are 0.0005 s for the rounded source RT literals, 0.0001 Da for
the approximate source m/z literals, and 1e-8 for the precise area literals. These
are declared native test tolerances, not an inference about C++ test macros.

Further cases test both termination modes, empty scans, inclusive windows,
visited-candidate rejection, strict thresholds, reverse ties, FWHM medians,
double-apex IM means, metadata coherence, repeated mixed-array calls, negative
settings, maximum counts/lengths, output ownership, area no-op/sentinel/failure
behavior, and progress errors. One hundred independently generated small grids
are checked against direct per-mass grouping and weighted RT arithmetic.
Private tests pin source iterative-mean literals, exact vector-initialization
work, and shared work exhaustion across accepted traces. The accompanying JSON
manifest hashes all source inputs and records the projection digest.

# Experiment aggregation and XIC extraction

The native `MSExperiment` implements the complete public `aggregate`,
`extractXICs`, `aggregateFromMatrix` and `extractXICsFromMatrix` operation group
from Core SDK revision `54a232fe2cae9c590d5c997fa49d20e7769860fb`. The relevant
header and class test are byte-identical to the previous `6bfc0e4` reference.
This increment does not implement the separate area/export/rasterization or
ion-mobility APIs.

## API mapping

| Source | Native |
| --- | --- |
| `pair<RangeMZ, RangeRT>` | `kernel::MzRtRegion`, with `mz` and `rt` `NumericRange` fields and checked `new(min_mz, max_mz, min_rt, max_rt)` |
| `aggregate(ranges, level)` | `MSExperiment::aggregate` |
| `aggregate(ranges, level, reducer)` | `aggregate_with`, accepting `FnMut(&[Peak1D]) -> Result<f64>` |
| `extractXICs(ranges, level)` | `extract_xics` |
| `extractXICs(ranges, level, reducer)` | `extract_xics_with` |
| `aggregateFromMatrix` | `aggregate_from_matrix(&[[f64; 4]], level, MzAggregation)` |
| `extractXICsFromMatrix` | `extract_xics_from_matrix(&[[f64; 4]], level, MzAggregation)` |
| `"sum"`, `"max"`, `"min"`, `"mean"` | `MzAggregation::{Sum, Max, Min, Mean}`; case-sensitive `FromStr` |
| caller-selected native limits | `aggregate_with_limits` and `extract_xics_with_limits`, with a reducer and `&AggregationLimits` |

Matrix rows are `[min_mz, max_mz, min_rt, max_rt]`. The Rust array type enforces
four columns. Custom reducers receive a borrowed contiguous slice, including an
empty slice when a selected scan has no peaks inside the m/z window.

```rust
use openms::kernel::{MzAggregation, MzRtRegion};
use openms::MSExperiment;

let experiment = MSExperiment::new();
let window = MzRtRegion::new(499.5, 500.5, 30.0, 60.0)?;
let sums = experiment.aggregate(&[window], 1)?;
let chromatograms = experiment.extract_xics_from_matrix(
    &[[499.5, 500.5, 30.0, 60.0]],
    1,
    MzAggregation::Sum,
)?;
# Ok::<(), openms::Error>(())
```

## Source behavior

Both m/z and RT intervals include their endpoints. Input region order, matching
scan order, duplicate m/z values and duplicate RT scans are preserved. MS level
is an exact filter: zero does not select every level. Empty regions return an
empty outer vector immediately. No spectra at the requested level also returns
an empty outer vector. When that level exists, every region produces a row;
regions without matching RT scans produce empty rows.

The default reducer accumulates intensities in f32 and then widens to f64. The
matrix sum and mean reducers accumulate in f64. This distinction is intentional:
intensities `[16777216, 1, -16777216]` give default sum zero and matrix sum one.
Min/max retain the first equally extreme value and return zero for an empty
slice. Mean and sum also return zero for empty slices. Finite negative m/z, RT
and intensity values retain their arithmetic behavior.

XICs retain the full f64 scan RT and narrow their reducer result once to f32
intensity. Every output chromatogram, including an empty row, gets product m/z
`(min_mz + max_mz) / 2.0` in that operation order. Its other fields have native
defaults. `MSChromatogram::product` uses the existing `metadata::Product` record
with target m/z, isolation offsets and CV metadata. Clone, peak selection,
sorting and `clear(false)` preserve it; `clear(true)` resets it.

## Checked boundaries and resource accounting

Selected scans must have finite, nondecreasing RT. Scans of other MS levels are
not validated. For selected scans covered by an RT window, all peak m/z values
must be finite and nondecreasing, because window lookup requires sorted data.
Only intensities inside selected m/z windows are consumed and checked. Other
arrays, precursors, chromatograms, metadata and identification graphs are not
traversed or cloned. The source assumes correctly partitioned binary-search
inputs; native unsorted input returns `UnsortedData` rather than an unreliable
selection. Publicly edited malformed region fields are checked when the
operation needs regions; matrix row construction validates before the source
no-matching-level early return.

Nonfinite inputs/results, f32 sum or XIC narrowing overflow, and midpoint
overflow return checked errors. The source can store nonfinite values in these
cases. There is no midpoint rearrangement that would silently avoid the source
addition overflow. Region inversion is rejected; zero-width regions remain
valid and include exact-coordinate peaks.

The default `AggregationLimits` are 50,000,000 units of library work,
256 MiB cumulative vector-allocation bytes and 10,000,000 output points. Limits
are caller-configurable through the two limited custom-reduction operations;
there is no separate lifetime counter. Work covers scan/region passes,
conservative binary-search comparisons, one-time coordinate validation per
visited scan, output visits and the intensity-validation/reduction traversal
for every window. An overlapping window consumes additional work. Byte
accounting includes vector elements plus allocation overhead, and allocations
use `try_reserve_exact`. All rows contribute to the same point bound, checked
before reducers run. The implementation stores scan references and interval
indices, rather than a dense scan-by-region mapping.

Outputs are owned and returned only on success. The experiment and its existing
chromatograms remain unchanged on every error. Native callbacks execute
sequentially in region order and then scan order; the source OpenMP API does
not promise callback order. Callback side effects are not rolled back, and
arbitrary caller code cannot be preempted or included in the library's work
bound. The library meters the provided slices; a custom reducer is responsible
for its own extra work/allocation. `MzAggregation::reduce` is independently
bounded for standalone use.

## Evidence

[`tests/experiment_aggregation.rs`](../tests/experiment_aggregation.rs) includes
the pinned four-scan source dataset and its aggregate, all four matrix-reducer,
and XIC literals, including product targets. Independent tests cover precision,
inclusive edges, duplicates, signed data, empty slices/rows, f64 RT retention,
finite-value failures, source early returns, callback errors, shared resource
limits and Product retention. No C++ execution or native-generated numerical
golden is used.

Source hashes and exact test anchors are recorded in
[`experiment_aggregation_provenance.json`](../tests/data/experiment_aggregation_provenance.json).

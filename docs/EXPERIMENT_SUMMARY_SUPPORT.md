# Experiment TIC, chromatogram organization and summaries

This native group follows Core SDK
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. It adds checked methods to
`MSExperiment`, with implementation in
[`experiment_summary.rs`](../src/kernel/experiment_summary.rs). It does not
certify the complete MSExperiment, range-manager or area-iterator headers.

| Native operation | Source operation and behavior |
| --- | --- |
| `calculate_tic_binned(rt_bin_size: f32, ms_level: u32)` | Complete finite calculateTIC bin branches: one f32 sum per selected spectrum, then absolute RT redistribution when the bin size is positive |
| `chromatogram_ranges()` | Chromatogram RT/intensity ranges and Product m/z extension from updateRanges |
| `combined_ranges()` | Combined RT/m/z/intensity extension from all spectra and chromatograms |
| `sort_chromatograms(sort_rt: bool)` | Product-m/z ordering, optionally sorting each chromatogram's points and all aligned arrays by RT |
| `total_peak_count()` | getSize: spectra peaks plus chromatogram points, returned as checked u64 |
| `contains_scan_of_level(level: usize)` | Exact level membership, including empty scans |
| `has_zero_intensities(level: usize)` | Any +0 or -0 intensity in a scan at exactly that level; empty scans do not qualify |
| `clear_meta_data_arrays()` | Remove float/integer/string arrays from spectra only and report whether any arrays existed |

Every method returns `Result`; each has a corresponding `_with_limits` method
whose final argument is `kernel::SummaryLimits`. Default wrappers use 50 million
work units, 256 MiB cumulative temporary-vector bytes and 10 million generated
points. These limits are independent of format-reader and older kernel limits.

## TIC arithmetic

Level zero includes all spectra for TIC. For level membership and zero-intensity
predicates it is an exact level, not a wildcard. Stored chromatograms are never
used to calculate a TIC; new chromatogram metadata, Product and precursor fields
have their usual defaults.

Zero and negative finite bin sizes preserve storage order, duplicate RTs and
empty scans as zero-intensity points. Positive bins require nondecreasing RTs
among selected scans. Empty input remains empty, including with all limits zero;
a singleton or multiple scans at exactly one RT produce one resampled point.
Finite signed coordinates and intensities are accepted. Missing matching MS
levels produce an empty chromatogram.

The grid has `ceil((end-start)/spacing + 1)` points and positions
`start + i*spacing`. Spacing is the caller's f32 value widened to f64, matching the
source argument and parameter conversion. Intensities are accumulated in f32
per spectrum. Redistribution computes both absolute distances, multiplies the
intensity by the opposite distance, divides by their sum, adds the previous bin
value, then narrows each updated bin to f32. The last grid point can be beyond
the input end, and any remaining tail is accumulated at the right boundary.

A small private loop preserves that expression order. The existing general
native resampler instead computes a fraction and `1-fraction`, and constructs
its count using `ceil(q)+1`. Those expressions need not round identically. For
example, with RTs zero and the next f64 above one, and bin size one, source
`ceil(q+1)` emits two points while `ceil(q)+1` emits three. No existing resampler
or its API is changed here. The source spacing-warning log is not emitted by
this pure checked query.

The older `calculate_tic(ms_level)` remains unchanged: it is an unbounded,
unbinned convenience method returning a chromatogram directly. The new method
checks finite consumed RTs/intensities, every f32 sum/output, source signed-int
grid representability, finite/advancing grid coordinates and allocation/work
bounds. NaN/infinite bin sizes, f32 overflow, reversed selected RT order for
positive bins, unrepresentable grids and precision-stalled grid steps are
errors. Unrelated peak m/z, arrays, identifications and acquisition metadata are
not inspected or cloned.

## Bounds and ordering

The range queries compute current values on demand, independently of input
sorting. Empty dimensions are `None`; no stale mutable cache is maintained.
Empty spectra contribute their RT, and every stored chromatogram contributes
its Product m/z even with no points. The string metadata key `product_mz` is not
a substitute for the Product field. RT/intensity from chromatograms come only
from points. Spectrum ranges are combined first, then chromatogram ranges, so
an equal zero endpoint retains the source's first sign of zero.

`ExperimentRanges` contains RT, m/z and intensity only. Ion-mobility ranges,
MS-level range-manager objects, area iterators and full RangeManager operations
remain separate gaps. The existing spectrum-only `ranges(ms_level)` query is
unchanged. The new combined query includes every level and all chromatograms.
Only consumed scalar values are validated; isolation offsets, CV payloads and
unrelated annotations are not traversed by these queries.

Chromatogram sorting moves whole records and never clones their metadata.
Product-m/z ties use stable input order, a deterministic native extension to the
source outer `std::sort`; RT ties preserve source stable order. With RT sorting
enabled, all nonempty float/integer/string arrays must have the peak count;
empty arrays stay empty. The native check also rejects malformed arrays in an
already-sorted chromatogram, where the source's early return can skip that
check. With RT sorting disabled, RTs and arrays are unconsumed and remain as-is.
Intensity values are not sort keys and are not independently validated.

All permutations, shape checks and numeric checks finish before any mutation.
The final swaps move String buffers, CV records and identification payloads
without copying their contents. `clear_meta_data_arrays` similarly preflights
all spectrum array descriptors and String destruction counts before removing
anything. It accepts malformed array lengths because clearing does not index
their contents. Chromatograms and their arrays are retained, and any array
counts as present even if its data vector is empty. Replacing vectors with empty
vectors releases the source containers' former storage.

## Resource and error boundaries

Work is shared across an entire operation, including all scans/chromatograms,
TIC raw and raster passes, permutation construction, array checks, swaps and
String destruction counts. Sorting uses a bounded iterative merge of indices;
it has no recursive input-dependent traversal or unmetered sort comparator.
Vector capacity and conservative allocation overhead are charged before
`try_reserve_exact`. Byte limits cover new temporary storage, not borrowed
payloads that are moved unchanged. Counts and membership queries allocate no
output collections. Generated-point limits include both intermediate raw TIC
points and final raster points cumulatively.

Queries return an owned complete result or an error. Mutating methods leave all
records and arrays unchanged on any checked failure. Ordinary allocator/process
failures outside Rust's fallible reservation guarantees are not converted into
application errors. This increment adds no dependencies, parallel runtime,
C++ execution, hidden metadata clone or public tuning state on MSExperiment.

## Validation and integration

[`experiment_summary.rs` tests](../tests/experiment_summary.rs) include upstream
TIC and chromatogram-order literals; source range examples; independent grid
rounding/f32 accumulation/signed-zero cases; empty/duplicate/signed inputs;
three-cycle aligned stable permutations; unchanged String-buffer ownership;
shared resource failures and atomic mutation failures. The 14 new tests plus
13 existing aggregation tests and 19 kernel tests pass under Rust 1.98 with all
features and Rust 1.85 without default features. Focused strict Clippy and
Rust 2024 formatting checks pass on both supported compilers.

[`experiment_summary_provenance.json`](../tests/data/experiment_summary_provenance.json)
records the exact source files, hashes and literal anchors. No C++ reference
program was built or run.

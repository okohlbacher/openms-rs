# Filtered bulk peak export

The non-mobility `MSExperiment::get2DPeakData` and
`get2DPeakDataPerSpectrum` paths follow Core SDK
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. The implementation is
[`peak_data.rs`](../src/kernel/peak_data.rs), reexported through `openms::kernel`.
This group does not include the mobility overloads or the separate unfiltered
`get2DData`/`set2DData` template group (including metadata arrays and mass-trace
expansion). [Concrete plain/rich 2D conversion](EXPERIMENT_2D_SUPPORT.md) is
implemented separately. Neither group completes the MSExperiment header.

## API

`MSExperiment::get_2d_peak_data(bounds: AreaBounds, ms_level: usize)` returns
`FlatPeakData`, with public parallel `rt`, `mz` and `intensity: Vec<f32>` fields.
`get_2d_peak_data_per_spectrum(bounds, ms_level)` returns `SpectrumPeakData`,
with `rt: Vec<f32>` and parallel nested `mz`/`intensity: Vec<Vec<f32>>` fields.
Both methods return `Result` and start from empty output.

`append_2d_peak_data(bounds, ms_level, &mut output)` and
`append_2d_peak_data_per_spectrum(...)` preserve the source append semantics.
Every method has a `_with_limits` variant whose final argument is
`PeakDataLimits`. Public records have ordinary `Default`, `Clone`, `Debug` and
`PartialEq`; no floating `Eq` is added.

Use `AreaBounds::new(min_rt, max_rt, min_mz, max_mz)` for the source scalar
argument order. Optional `rt`/`mz` dimensions set to `None` are unrestricted
over finite coordinates. Source `Size` MS levels pass through `UInt`, `uint8_t`
and `int8_t` before unsigned comparison. The native export explicitly applies
`((ms_level as u32) as u8 as i8) as u32` on every platform. Thus 256 selects
level zero, 255 selects `u32::MAX`, and zero is never a wildcard. For an exact
u32 level without source narrowing, the separate borrowed area API is available;
the export methods themselves preserve source compatibility.

## Selection, precision and row grouping

Selection uses the same scalar area traversal as the source: closed RT/m/z
bounds, original spectrum-then-peak order, all boundary duplicates and no empty
rows for scans without selected peaks. Full-precision f64 coordinates are
filtered **before** conversion to f32. Flat output repeats the converted scan
RT for every peak, narrows m/z once, and copies each stored f32 intensity.
Chromatograms, annotations, arrays and identifications are not exported.

The source name `PerSpectrum` does not guarantee one row per spectrum. Before
every selected **peak**, it compares the raw f64 scan RT to the previous RT
stored in f32. A mismatch creates a new row and updates the remembered f32 RT.
Consequently:

- Equal RTs exactly representable in f32 merge into one row across adjacent
  selected spectra, even across intervening empty or other-level scans.
- RT `0.1` is not exactly representable in f32, so two peaks at raw RT `0.1`
  create two rows, each with stored RT `0.1f32`.
- After such a row, a later scan whose raw RT equals the widened `0.1f32`
  joins the previous row. Grouping is neither spectrum identity nor equality
  of rounded RT keys.
- The remembered value resets to `-1f32` on each call. If the first selected
  raw RT is exactly `-1`, source appends into the caller's existing final m/z
  and intensity row without adding or changing an RT value. The native append
  method preserves this finite behavior even if the old final RT is unrelated.
  With no existing row, source `back()` is undefined; native construction
  returns a checked error. Flat output accepts RT `-1` normally. If an earlier
  selected scan has another RT, later RT `-1` creates a normal row.

These behaviors are retained intentionally and tested with exact oracles. No
ordinary scan-grouping or f64-preserving export is claimed by this API.

## Checked native boundaries

The shared area preflight validates finite nondecreasing RTs across the entire
experiment and finite nondecreasing m/z values in every scan, including scans
outside the selected region or at other MS levels. Invalid ordering returns
`Error::UnsortedData`; invalid boundaries, values or limits return
`Error::InvalidValue`. Like the scalar area mapping, this assumes represented
default scan mobility. The native spectrum model cannot yet represent or
filter source nonfinite drift times. See
[`AREA_ITERATION_SUPPORT.md`](AREA_ITERATION_SUPPORT.md) for that explicit gap.

Selected RT/m/z values must narrow to finite f32, and selected intensities must
be finite. Underflow and signed-zero narrowing follow Rust's f64-to-f32 cast;
finite signed values are accepted. Unselected intensities and unselected f32
conversions are not inspected. Source nonfinite output is a checked native
error. Existing output numeric values are retained as payload, without changing
or revalidating their bit patterns.

When at least one peak is selected, native append requires aligned existing
parallel vectors; nested m/z/intensity rows must also align. This is an explicit
native invariant beyond the independent source output-vector parameters.
Empty selection leaves the entire output unchanged, even if its unused prefix
is malformed or over an output limit. Input ordering/raw-size checks still run.

Append stages a complete bounded numeric result, then commits once. A late
conversion, missing first row, allocation, work or shape error leaves the
caller and experiment unchanged. Existing prefix storage may be replaced on a
successful append; pointer identity/capacity is not preserved. No unrelated
metadata graph is cloned. The two exports share one remaining-work and byte
budget with area preflight rather than starting a fresh budget after selection.

Default limits allow one million raw spectra, ten million raw peaks, ten
million final output points, one million final nested rows, 50 million work
units and 256 MiB cumulative newly allocated storage. Existing and new output
count together. The nested-row cap is unused for flat output. Work includes
global ordering checks, binary comparisons, area traversal, row-layout pass,
numeric copying/writing, row descriptors and nested-vector destruction. Bytes
include the area plan, growing row-size scratch, final top-level vectors and
all old/new nested numeric buffers with conservative allocation overhead.
Capacity/count arithmetic is checked before fallible reservations. Limits are
per call and configurable; ordinary caller `Clone`, comparisons and consumers
are outside these checked operations. A completely empty experiment permits
all-zero limits. No recursion or dependency is introduced.

## Evidence

[`tests/peak_data.rs`](../tests/peak_data.rs) pins all five source row-layout and
four source flat-layout literal queries. It adds append/history tests, the
raw-f64/f32 grouping corner cases, RT `-1`, signed values, exact boundary
selection before narrowing, level conversions, empty scans, atomic failure and
shared limits. Six hundred forty deterministic queries compare both layouts
against direct linear enumeration. The unchanged ten area tests are rerun for
the shared-budget extraction.
All 20 tests and strict focused Clippy pass on Rust 1.98 with all features and
Rust 1.85 without default features. Rust 2024 formatting checks pass. Independent
read-only source and budget review found no actionable discrepancy.

[`peak_data_provenance.json`](../tests/data/peak_data_provenance.json) records
the source pin, hashes and implementation/test anchors. No C++ reference
program was compiled or executed.

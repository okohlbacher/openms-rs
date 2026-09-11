# Concrete unfiltered 2D point conversion

[`experiment_2d.rs`](../src/kernel/experiment_2d.rs) ports the concrete plain
and rich-point branches of MSExperiment `get2DData`/`set2DData` at Core SDK
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. It uses the native
[Peak2D/RichPeak2D values](PEAK2D_SUPPORT.md). The Feature mass-trace template
specialization, arbitrary user container instantiations and all mobility
conversion remain separate gaps. The filtered f32 bulk exports are separate
operations with different ordering/grouping rules.

## Public operations

| Method | Native result / source mapping |
| --- | --- |
| `get_2d_data()` | `Result<Vec<Peak2D>>`: unfiltered MS1 points with f64 coordinates and f32 intensity |
| `append_2d_data(&mut Vec<Peak2D>)` | Append the same points while retaining the existing prefix |
| `set_2d_data(&[Peak2D])` | Replace current data with grouped MS1 spectra; return `Result<MSExperiment>` containing all previous ownership |
| `set_2d_data_rich(&[RichPeak2D], &[String])` | Same replacement, with requested numeric metadata copied into named f32 arrays; return the previous experiment |

Every method has a `_with_limits` variant with a final `Data2DLimits` argument.
The successful setter result is intentionally useful ownership, not a count:

```rust
use openms::kernel::{MSExperiment, Peak2D};

fn replace_points(experiment: &mut MSExperiment) -> openms::Result<MSExperiment> {
    let points = [Peak2D::new(2.0, 3.0, 1.0), Peak2D::new(5.0, 6.0, 4.0)];
    let previous = experiment.set_2d_data(&points)?;
    // The caller decides when to release or reuse all previous data.
    Ok(previous)
}
```

## Source selection and replacement semantics

Export visits only spectra at exactly MS level one, in storage order, and then
all their peaks in storage order. No RT/m/z sorting, filtering, isotope handling
or f64-to-f32 coordinate conversion occurs. Empty scans produce no points.
Chromatograms and all metadata are ignored. The output contains only RT, m/z and
intensity. Unlike filtered area extraction, unsorted RTs and m/z values are
valid. Appending repeats this complete sequence after the existing prefix.

Import requires nondecreasing RTs. Adjacent points at exactly equal f64 RT share
one new spectrum; there is no f32 RT comparison. m/z order is preserved even if
unsorted, and each new scan is MS1. Negative finite RTs/m/z/intensities and exact
duplicate positions are accepted. Equal positive/negative zero RTs share a scan
whose RT retains the first point's sign. Plain roundtrip therefore preserves
numeric equality, while later equal-zero RT bit patterns may become that first
sign. Empty input creates an empty experiment without inspecting unused names.

Successful source `clear(true)` removes spectra, chromatograms and experiment
metadata. Native success gives the current map those same cleared/default fields
plus the new scans. All previous spectra, chromatograms, metadata and nested
identifications are returned unchanged in the old map; their buffers are moved,
not copied. Existing public `clear` behavior is unchanged. Source clears the map
before checking all later inputs; native import instead prepares the complete
result and atomically replaces the current map only after success.

The source RT ordering exception is active in both release and debug builds.
Native import returns `Error::UnsortedData` on decreasing RT, and checks finite
consumed coordinates/intensity with `Error::InvalidValue`. Raw peak value types
themselves remain permissive. Export checks only emitted points; wrong-level
scans and an empty MS1 scan's unused RT are not numerically validated. Existing
append prefix values are retained as payload. No-selection export leaves that
prefix unchanged without applying an output limit to it.

## Rich metadata arrays

Requested names create one float array each, in the given order, in every new
scan. Duplicate names remain duplicate arrays, including an empty name if
requested. Each array stays aligned with the new scan's peaks. Missing metadata
keys produce quiet NaN, exactly as source; the importer deliberately does not
run a generic finite-array validator afterward. This valid source missing-value
sentinel can be rejected by other native algorithms/formats that separately
require finite auxiliary arrays; those contracts are not changed here.

Present integer values cast directly from stored i64 to f32, preserving the
source float conversion instead of routing through f64. For example
`2^62 + 2^38 + 1` rounds differently through f64; a direct-bit oracle covers this
boundary. Present f64 values narrow directly to f32. Units on numeric values
do not alter scalar conversion. Present Empty, string or list values are errors:
the source rejects Empty and its inactive-union read for other nonnumeric types
has no defined numeric behavior. Present values that cannot produce finite f32
are checked errors. Missing keys remain NaN rather than errors.

No other input metadata, unique IDs, units, annotations or descriptions are
copied into the new spectra. New arrays have default empty descriptions and
processing records. Unrequested metadata payloads are not traversed or validated.
Export does not reconstruct these arrays as rich metadata, matching the source
loss of annotations when converting back to plain points.

## Resource and ownership boundaries

Default limits allow one million scanned/generated spectra, ten million points,
one million generated metadata arrays, one MiB of aggregate requested-name bytes,
50 million work units and 256 MiB of cumulative newly allocated storage. Point
caps include old plus new output for append; import caps count all input points.
Export skips wrong-level peak payloads and counts only selected MS1 points.
Limits are configurable per operation.

Construction, contiguous RT group planning, numeric copies, name replication,
array construction, potential new-output disposal and metadata lookup charges
share one budget. Each BTree metadata query is charged conservatively using the
entire map's entry count times the requested-key length before lookup; this may
reject a very large metadata map even though its actual tree query is faster.
It avoids scanning unrelated values/units/lists or assuming unchecked key work.
Allocation and capacity arithmetic is checked before fallible reservations.
Peak/row/name/array buffers and growing group-size scratch count cumulatively.

All checked failures leave current experiment/output contents unchanged. Prior
experiment and rich-value payloads are returned by ownership to avoid a hidden
deep graph clone or unbounded destructor inside the checked replacement.
Ordinary caller drop/Clone/format/hash operations and consuming a returned old
map are outside these per-call resource guarantees. No source input controls
recursion. Empty import/export supports all-zero limits where its consumed
input is empty. No dependency or parallel runtime is added.

Eight [conversion tests](../tests/experiment_2d.rs) preserve source MS1 sequence,
plain roundtrip and `111.1`/`333.3` metadata fixtures; add signed-zero grouping,
duplicate names, NaN sentinels, direct integer narrowing, old-map buffer identity,
late atomic errors and cumulative limits; and check 100 deterministic plain
roundtrips plus independent unsorted MS1 export enumeration. Source hashes are
recorded in [peak2d_provenance.json](../tests/data/peak2d_provenance.json). No C++
build or reference execution is used.

## Focused validation

All 14 tests (six value and eight conversion tests) pass with Rust 1.98.0 /
all features and Rust 1.85.0 / no default features, locked and offline with
two build jobs. Strict focused Clippy passes in both configurations. An
independent source review of value semantics, conversion arithmetic, ownership
and cumulative preflight found no actionable issue.

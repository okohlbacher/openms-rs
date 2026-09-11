# Borrowed RT/m/z peak-area iteration

The native scalar area traversal follows Core SDK
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. It is implemented in
[`area_iteration.rs`](../src/kernel/area_iteration.rs) and reexported from
`openms::kernel`. This group covers the scalar RT/m/z `areaBegin` and
`areaBeginConst` paths. It does **not** complete `AreaIterator.h` or
`MSExperiment.h`: the scan-mobility `RangeManager` overload, `lowIM`/`highIM`,
and `getDriftTime` remain unmapped pending the spectrum mobility model and its
transport. The standalone Mobilogram type is not that spectrum model.
Even the source scalar overload internally applies a full finite mobility
range; the native scalar mapping assumes the currently represented default
scan mobility. It does not reproduce exclusion of source scans carrying a
nonfinite drift time, a value the current native spectrum model cannot store.

## API mapping

| Native API | Source operation / native contract |
| --- | --- |
| `MSExperiment::area_begin(min_rt, max_rt, min_mz, max_mz, ms_level)` | Scalar `areaBeginConst`, with source MS-level narrowing |
| `MSExperiment::area_begin_mut(...)` | Scalar mutable `areaBegin`, with the same narrowing |
| `area_iter(AreaOptions)` / `area_iter_mut(AreaOptions)` | Shared checked traversal using owned RT/m/z options and an exact native u32 MS level |
| `area_iter_with_limits(options, limits)` / `area_iter_mut_with_limits(...)` | The same operations with caller-specified resource limits |
| `AreaBounds::new(min_rt, max_rt, min_mz, max_mz)` | Checked closed scalar boundaries, RT first |
| `AreaBounds { rt, mz }` | Optional closed dimensions; `None` includes all finite coordinates in that dimension |
| `AreaOptions::new(bounds, level)` | Exact u32 level; default level is one |
| `AreaOptions::source_compatible(bounds, level)` | Explicitly reproduces the source unsigned-byte/signed-byte MS-level conversion |
| `AreaIter::peek()` / iterator items | Checked dereference, scan reference, RT and original spectrum/peak indices |
| `AreaIter::default()` / exhausted `next()` | Source end sentinels map to an empty iterator and `None` |
| `Clone`, `PartialEq`, `Eq` on `AreaIter` | Independent cursors sharing an immutable interval plan; current-peak address equality, with all exhausted iterators equal |
| `Iterator`, `ExactSizeIterator`, `FusedIterator` | Forward traversal, exact remaining peak count and permanent exhaustion |

`AreaPeak` holds `spectrum_index`, `peak_index`, `&MSSpectrum` and `&Peak1D`.
`AreaPeakMut` holds the original indices, the scan's RT/MS-level snapshots and
`&mut Peak1D`. The mutable item intentionally does not expose a whole-spectrum
reference overlapping its exclusive peak borrow. Multiple yielded mutable
references may coexist; no peak is yielded twice. Mutable iterator cloning,
raw pointer conversions and aliased mutable source iterator copies are not
exposed. Ordinary Rust assignment/moves and immutable `Clone` provide the
safe ownership equivalents.

```rust
use openms::kernel::{AreaBounds, AreaOptions, MSExperiment};

fn sum_area(experiment: &MSExperiment) -> openms::Result<f64> {
    let bounds = AreaBounds::new(30.0, 60.0, 400.0, 405.0)?;
    let points = experiment.area_iter(AreaOptions::new(bounds, 1))?;
    Ok(points.map(|point| f64::from(point.peak.intensity)).sum())
}
```

## Selection and source details

Both dimensions are inclusive: RT uses lower/upper bounds, and each qualifying
scan uses lower/upper m/z bounds. Equal RTs, duplicate m/z values and zero-width
ranges retain every matching peak. Traversal follows the original spectrum
order, then peak order. Empty scans and scans with no matching peaks are
skipped. No chromatogram is visited. Zero is an exact MS level, not a wildcard.
Finite negative coordinates, signed zero and the finite f64 extrema are valid.

`MzRtRegion` already represents closed RT/m/z boundaries; its explicit
conversion to `AreaBounds` keeps those fields. Its constructor takes m/z first,
whereas `AreaBounds::new` and source `areaBegin` take RT first. This conversion
does not equate the operation's MS-level rules, validation or budgets with the
aggregation API.

Source `MSExperiment` accepts an unsigned integer level but passes it through
`AreaIterator::Param(uint8_t)` into an `int8_t` field, then compares it as the
spectrum's unsigned level. The compatibility constructor and scalar wrappers
therefore implement `(requested as u8 as i8) as u32`: 256 selects level zero;
255 selects `u32::MAX`; 128 selects `u32::MAX - 127`. The ordinary native options
constructor never truncates a requested level. The conversion is not inferred
from typical MS1/MS2 data; dedicated boundary tests pin it.

Source mutable iteration fixes the current scan's peak endpoints before that
scan is traversed. Editing a yielded peak's m/z does not recalculate those
endpoints. Native construction fixes all selected scan intervals in advance;
safe borrowing prevents modification of unvisited scans, so the reachable
behavior agrees. Coordinates and intensities can be edited through the peak
reference, and the remaining preselected peaks are still yielded. A later new
iterator construction can reject an experiment made unsorted by those edits.
Arrays and metadata remain attached to their original peaks/scan positions and
are never copied, reordered or rewritten by area iteration.

Immutable iterator equality follows current peak identity rather than region
identity. Two different regions positioned at the same original peak compare
equal; equal-valued peaks in different experiments do not. End equality is
independent of experiment ownership. `peek()` returns `None` at end rather than
invoking source undefined dereference behavior; `next()` returns the current
item and then advances according to the usual Rust iterator convention.

## Validation, resources and errors

Construction always checks finite ordered bounds, nondecreasing RTs across the
whole experiment, and finite nondecreasing m/z coordinates in **every** scan.
This makes source `isSorted(true)` preconditions checked errors even when a
malformed scan is outside the region or has another MS level. NaN/infinite RT
or m/z values are rejected rather than entering binary search with an invalid
ordering. Intensities, precursor records, identifications, chromatograms and
array contents/alignment are unconsumed and are not validated or traversed.
In particular, a NaN intensity may be borrowed unchanged.

`AreaLimits::default()` allows one million spectra, ten million total raw
peaks, 50 million work units and 256 MiB of newly allocated plan storage. All
limits are configurable per call. Raw peak and scan caps include unselected
input; an empty region cannot bypass validation. A completely empty experiment
can use all-zero limits. Work is shared across global validation, actual
binary-search comparisons, interval construction and reserved future scan and
peak traversal. There is no per-scan reset. The byte budget conservatively
charges at most one three-index interval per RT-candidate scan plus a fixed
shared-plan header allowance, before fallible vector reservation. Neither peaks
nor nested metadata are cloned. The iterator's internal plan does not grow
during traversal, and its destructor traverses no borrowed metadata graph.

All checked failures happen during construction, before any mutable reference
escapes. Iteration itself returns borrowed items without a fallible midstream
state or additional allocations. Ordinary caller-driven repeated `peek`,
immutable iterator clones/replays, standard iterator consumers, and cloning a
borrowed spectrum are outside this per-construction budget. Rust allocator or
process failures outside fallible vector reservation remain ordinary platform
failures. Invalid bounds/resources use `Error::InvalidValue`; invalid ordering
uses `Error::UnsortedData`.

## Evidence

[`tests/area_iteration.rs`](../tests/area_iteration.rs) reproduces the source
five-scan fixture, all its scalar quadrants/empty windows, original indices and
RTs, and the source MSExperiment mutable `4711` edit. Independent tests cover
closed duplicates/extrema, exact versus narrowed MS levels, cursor identity,
coexisting mutable borrows, unrelated malformed data, whole-input ordering and
late shared-budget rollback. Eight hundred deterministic queries over small
generated grids compare both const and mutable traversal with independent
linear enumeration; unselected peaks are checked unchanged.
All ten tests and focused strict Clippy pass on Rust 1.98 with all features and
Rust 1.85 without default features; Rust 2024 formatting checks pass.

[`area_iteration_provenance.json`](../tests/data/area_iteration_provenance.json)
records source hashes and exact test/implementation anchors. No C++ build or
reference execution, new dependency, spectrum field or generic RangeManager
implementation is included.

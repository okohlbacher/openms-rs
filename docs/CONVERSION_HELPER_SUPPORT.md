# Map conversion support

[`src/kernel/conversion_helper.rs`](../src/kernel/conversion_helper.rs) ports
`KERNEL/ConversionHelper.h` at Core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The header was the last unmapped
`KERNEL` header besides the two documented deferrals; it declares one class,
`MapConversion`, with three static overloads that move data between the peak,
feature and consensus containers.

The containers themselves live in
[`src/kernel/features.rs`](../src/kernel/features.rs) and
[`src/kernel.rs`](../src/kernel.rs); their own members are documented in
[MAP_OPERATIONS_SUPPORT.md](MAP_OPERATIONS_SUPPORT.md).

## API mapping

Every public member of the header appears below.

| Source member | Rust counterpart |
| --- | --- |
| `class MapConversion` | `MapConversion`, a unit struct: the source class holds no state and exists only to scope three static functions |
| `static void convert(UInt64 input_map_index, PeakMap& input_map, ConsensusMap& output_map, Size n = -1)` | `MapConversion::peak_map_to_consensus(input_map_index, &input_map, n, &mut generator) -> Result<ConsensusMap>` |
| `static void convert(ConsensusMap const& input_map, const bool keep_uids, FeatureMap& output_map)` | `MapConversion::consensus_to_feature_map(&input_map, keep_uids, &mut generator) -> Result<FeatureMap>` |
| `static void convert(UInt64 input_map_index, FeatureMap const& input_map, ConsensusMap& output_map, Size n = -1)` | `MapConversion::feature_map_to_consensus(input_map_index, &input_map, n) -> Result<ConsensusMap>` |
| the `Size n = -1` default | `n: Option<usize>`; `None` is the source's `Size(-1)` sentinel and keeps everything |
| the `output_map` out-parameter | the returned container. The source's "previous contents are dropped" contract is automatic, and a failure cannot leave the caller with a half-written map |
| the source's `output_map.setUniqueId()` draws | the `generator: &mut UniqueIdGenerator` parameter, replacing the process-wide singleton. `feature_map_to_consensus` needs none, because it takes the container ID from its input |
| `output_map.updateRanges()` at the end of each overload | not needed: `FeatureMap::ranges` and `ConsensusMap::ranges` compute on demand, so the documented "range queries reflect the new contents" guarantee holds with no cache to refresh |
| — | `MapConversion::MAX_ITEMS`, the element ceiling, with no source counterpart |

The overload set becomes three named functions because the port's naming rule
forbids overloading; the names say which direction each one runs.

## Preserved source conventions

* **The header `size` fields disagree on purpose.** The `PeakMap` overload sets
  the column header `size` for `input_map_index` to the number of peaks it
  actually wrote (`ConversionHelper.cpp:47`); the `FeatureMap` overload sets it
  to `input_map.size()`, the *full* input size, even when `n` truncated the copy
  (`ConversionHelper.cpp:104`). The header's own `@note` records this, and the
  class test asserts it: after `convert(33, fm, cm, 2)` on a three-feature map,
  `cm.size()` is 2 and `cm.getColumnHeaders()[33].size` is 3.
* **The feature-map conversion does not sort.** The first `n` features are taken
  in input order; `n` is only useful after pre-sorting. The header says so and
  the port repeats it.
* **The feature-map conversion keeps the input's container unique ID**
  (`ConversionHelper.cpp:98`), which the source itself calls "an arguable design
  decision". Callers that need a fresh identity overwrite it afterwards.
* **The feature-map conversion stamps `map_index`.** It builds each element
  through `ConsensusFeature(UInt64, const BaseFeature&)`, which copies the
  feature with `BaseFeature(element, map_index)` and therefore writes a
  `map_index` meta value onto every copied peptide identification
  (`ConsensusFeature.cpp:35`). The port calls
  `BaseFeature::clone_with_map_index` first, for exactly that reason.
* **The consensus-to-feature conversion copies only the `BaseFeature` slice.**
  The source resizes the output with default-constructed features
  (`ConversionHelper.cpp:56`) and then assigns `f.BaseFeature::operator=(c)`
  (`:74`), so position,
  intensity, quality, charge, width, meta values, peptide identifications and
  identification-graph references survive, while the consensus handles and
  ratios are dropped and the produced features have no convex hulls, no
  per-dimension qualities and no subordinates.
* **`keep_uids` covers both levels.** `true` keeps the container unique ID and
  every element's; `false` replaces all of them.
* **The document identifier and the protein / unassigned peptide
  identifications are carried** by the consensus-to-feature and
  feature-to-consensus conversions; the peak-map conversion has nothing to
  carry.
* **Only MS level 1 peaks are converted.** The peak-map overload reads through
  `MSExperiment::get2DData`, which skips every other level
  (`MSExperiment.h:175`).
* **The peaks are written in descending intensity order**, and each handle's
  unique ID is its position in that order.

## Native differences

* **The peak count is clamped against the collected points, not
  `getSize()`.** The source clamps `n` against `MSExperiment::getSize()`, which
  sums the peaks of *every* spectrum plus every chromatogram point
  (`MSExperiment.cpp:744`), and then uses that count as the middle iterator of a
  `std::partial_sort` over the MS1-only vector `get2DData` produced
  (`ConversionHelper.cpp:24-44`). With any MS2 spectrum or chromatogram present
  — and the default `Size(-1)` reaches this path every time — the middle
  iterator is past the end and the following loop indexes past the end of the
  vector. This port clamps against the number of points actually collected,
  which is the evident intent, and `peak_map_to_consensus_clamps_against_collected_points`
  asserts it.
* **The descending sort is stable.** The source uses `std::partial_sort`, which
  is not; with ties at the cut-off, which peaks survive is unspecified. Sorting
  stably keeps spectrum-then-peak order among equal intensities, so the result
  is reproducible.
* **The input peak map is borrowed immutably.** The source takes `PeakMap&`
  only to call `updateRanges()` on it, which this port does not need.
* **Elements are validated.** `ConsensusFeature::from_feature` rejects a
  non-finite coordinate or a negative width; the source copies both unchecked.
* **No parallelism is lost.** `ConversionHelper.cpp` carries no OpenMP pragma.

## Checked boundaries and evidence

`MapConversion::MAX_ITEMS` is 10,000,000 elements per conversion, checked before
anything is allocated; the peak-map conversion additionally inherits the
ceilings of `MSExperiment::get_2d_data` (`Data2DLimits`), which bound the
spectra and point counts and the bytes staged. Each conversion builds its result
independently of its input, so a rejected conversion allocates nothing and
changes nothing.

Evidence is **tier 3 (source review)**: every literal in
[`tests/conversion_helper.rs`](../tests/conversion_helper.rs) is transcribed
from the pinned class test, hashed in
[`tests/data/conversion_helper_provenance.json`](../tests/data/conversion_helper_provenance.json).
No C++ was built or executed, and no retained C++ output exists for this header.
The ceilings and the element validation are tier 4 (Rust-only invariants).

### Class-test section coverage

`ConversionHelper_test.cpp`: 3 sections, all 3 ported.

| Section | Test |
| --- | --- |
| `convert(UInt64, FeatureMap const&, ConsensusMap&, Size n = -1)` (line 25) | `feature_map_to_consensus` |
| `convert(UInt64, PeakMap&, ConsensusMap&, Size n = -1)` (line 80) | `peak_map_to_consensus` |
| `convert(ConsensusMap const&, bool, FeatureMap&)` (line 104) | `consensus_to_feature_map` |

Native tests without a source section:
`feature_map_to_consensus_carries_records_and_stamps_map_index`,
`peak_map_to_consensus_clamps_against_collected_points`,
`consensus_to_feature_map_keeps_records_and_drops_handles`,
`conversion_ceilings`.

One transcription note: the third section asserts `TEST_EQUAL(cm[i], out_fm[i])`,
comparing a `ConsensusFeature` with a `Feature`. Only the `BaseFeature` part of
a consensus feature survives the conversion, and that part is what the port
compares.

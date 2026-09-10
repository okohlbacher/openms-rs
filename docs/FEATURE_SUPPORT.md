# Feature and geometry support

The native `kernel::features` and `kernel::geometry` modules port a defined subset of OpenMS4-core revision `7c029e8`. Their public types are also reexported from `kernel`. There is no C++ dependency, feature finder, feature grouping engine, `IdentificationData` reference graph, or FeatureXML/ConsensusXML reader in these modules.

## API mapping

| OpenMS source | Native Rust coverage |
| --- | --- |
| `BaseFeature` / measured `RichPeak2D` fields | `BaseFeature`: RT, m/z, intensity, quality, charge, width, unique ID, typed metadata, attached peptide identifications |
| `Feature` | `Feature`: composed `base`, RT/m/z qualities, mass-trace hulls, subordinate features, overall hull and containment |
| `FeatureHandle` | `FeatureHandle`: owned measured values and `(map_index, unique_id)` identity |
| `ConsensusFeature` | Sorted, unique handles; insertion, replacement, set union; mean, monoisotopic and decharge consensus; handle ranges |
| `FeatureMap` | Owned features, identifier/ID/metadata, protein runs and unassigned peptide IDs, stable sorting, checked selection, ID lookup and ranges |
| `ConsensusMap` | Owned consensus features, column headers, protein runs and unassigned peptide IDs, stable sorting, checked selection, ID lookup, ranges and column consistency |
| `ConvexHull2D` | Scan-interval hull construction, outline preservation, containment, compression, bounds and bounding-box expansion |
| Selected `DPosition<2>` / `DBoundingBox<2>` operations | `Point2D` and validated, inclusive `BoundingBox2D` with rectangular union |

`Feature` and `ConsensusFeature` compose a public `BaseFeature` and implement `Deref`/`DerefMut`: `feature.rt`, `feature.charge` and `feature.metadata` remain convenient field access. Set overall feature quality through `feature.quality`; dimension qualities are named `quality_rt` and `quality_mz`, avoiding unchecked numeric dimension indices. Width is an explicit field; `set_width` updates it and the source `FWHM` metadata alias together. Direct field assignment still requires the caller to keep the alias consistent for featureXML output.

Feature, map and column metadata uses owned typed `MetaInfo` values. Peptide/protein records are attached through the [identification module](IDENTIFICATION_SUPPORT.md) and validated recursively. Identification references, annotation-state inference, ratios, document provenance graphs, automatic unique-ID generation, map append/split and label interpretation are not implemented here. Maps additionally preserve `data_processing`, `loaded_file_path` and `loaded_file_type`; clearing metadata clears these fields too. Column headers preserve filename, label, feature count, unique ID and metadata. Consensus experiment type defaults to `label-free`; validation permits the source values `label-free`, `labeled_MS1` and `labeled_MS2`.

## Hull semantics

OpenMS `ConvexHull2D` is a scan envelope and can be non-convex. `from_points` and `add_points` group equal RT values into minimum/maximum m/z intervals; intervals are joined linearly between scans. The lower outline is emitted in ascending RT and the upper outline in descending RT, omitting duplicate endpoint vertices. It is not a mathematical convex-hull algorithm.

Containment includes boundaries. The source first checks the exact RT interval; if that fails, it also interpolates between the strict neighboring scans. Thus a point outside a narrow interior scan can still be accepted at that exact RT. This unusual source behavior is preserved and tested explicitly. Interpolation arithmetic overflow returns an error.

`set_hull_points` retains an ordered outline without inferring scan intervals, matching the C++ distinction. Containment on a nonempty outline-only hull returns `Unsupported`. Adding scan points to such a hull also returns an error, deliberately avoiding the C++ method's silent deletion of the existing outline. Explicitly `clear` it, or `expand_to_bounding_box` to create scan intervals, before adding points. Batch point addition and outline replacement validate all points before mutation.

`compress` removes an interior scan only when its m/z interval equals those of both immediate neighbors. Outlines are computed on demand: reading them never changes equality, and compression cannot leave a stale cached outline. Expanding an empty hull keeps it empty instead of constructing a rectangle from empty-range sentinel values.

A feature with one mass trace returns that hull unchanged. With multiple traces, `convex_hull` returns their rectangular bounding union, following the current C++ algorithm; empty traces contribute no bounds. `Feature::encloses` tests the individual mass traces, so gaps inside that overall rectangle remain excluded. `BoundingBox2D::union` likewise means the smallest enclosing rectangle, not an exact polygon/set union.

## Consensus calculations

Handles are stored privately in ascending `(map_index, unique_id)` order and exposed through an immutable slice. Individual insertion and replacement reject duplicate keys. `merge` follows source set-union behavior: an existing key retains the receiving consensus's measured values. Insertion never implicitly recomputes summary values. `clear` removes handles while retaining the summary.

`compute_consensus` takes arithmetic means of RT, m/z and intensity. `compute_monoisotopic_consensus` instead takes minimum m/z and means of RT and intensity. Both use the source's running charge vote: the most frequent charge wins; equal counts prefer smaller absolute charge. When both magnitude and count tie, the running winner is retained in handle order. This is not equivalent to sorting a final count table; for example, charges `[2, -2, -2, 2]` yield `-2`. The native absolute-value calculation handles `i32::MIN` safely.

`compute_decharge_consensus` sums intensity and averages RT and neutral mass, with optional intensity weighting. A handle's neutral mass is `mz * abs(charge) - adduct_mass`. The adduct defaults to `charge * PROTON_MASS_U`; a finite numeric `dc_charge_adduct_mass` value on the matching source feature overrides it. Source-map lookup uses the feature's unique ID, ignoring its map index, as in the C++ algorithm. The result stores neutral mass in `mz` and sets charge to zero. This applies to negative charges too: their default proton adduct is signed.

All calculations use source double-precision accumulation with final float32 intensity conversion. Empty consensus calculations return errors instead of dividing by zero. Decharging rejects unknown charge and absent/ambiguous source IDs. Weighted mode requires nonnegative intensities and positive total intensity; unweighted calculations accept finite signed intensities. Non-finite arithmetic or float32 intensity overflow returns an error. Summary assignments occur only after successful checks, so errors cannot partially update RT/mass/intensity/charge. Quality, width, ID and metadata are preserved.

## Container and validation policies

RT and m/z are float64; intensity, quality and width are float32; charge is int32; IDs and map indices are uint64. Coordinates, intensity and quality may be signed but must be finite. Checked operations require finite, nonnegative width. Bounding boxes have private, validated bounds; empty dimensions use `Option`, not numeric sentinels.

Zero means unassigned unique ID and may repeat. Assigned top-level feature IDs must be unique within each map; ID lookup rejects zero and ambiguous maps and returns `None` for a missing nonzero ID. These IDs are identities, not vector indices. Subordinate IDs are outside the top-level index. Selection indices are unique, zero-based vector positions and can reorder selected features; invalid or duplicate indices leave the map unchanged. Child features and metadata stay attached to the selected parent.

All map sorts are stable, including descending intensity/quality sorts. This preserves the C++ consensus-map guarantee and strengthens C++ feature-map sorting, which uses `std::sort`. Position order compares RT first, then m/z. Consensus size sort is descending; map sort compares full sequences of handle identities lexicographically, including feature IDs. Map sorting, selection and range computation validate before operating. Direct public edits can temporarily violate invariants; call `validate` afterward.

Validation traverses subordinate features iteratively. Checked operations reject a subordinate depth greater than 128 before recursively cloning a selected feature tree. This bounds checked-operation traversal/cloning, not arbitrary Rust `Clone`/`Drop` calls on a deliberately invalid tree assembled through public fields. Collections otherwise grow with owned input size; batch hull insertion uses sorting rather than repeated vector insertion. Neither geometry nor map operations allocate according to coordinate magnitude or ID value.

Feature-map RT/m/z ranges include top-level centroids and hull bounds; intensity ranges include top-level intensities. Subordinates do not contribute, following the source. Consensus-map ranges include both summaries and every handle in all three dimensions. Ranges are computed on demand, so public edits cannot leave a stale range cache.

`ConsensusMap::validate_consistency` additionally requires each handle's map index to have a column header and each `(filename, label)` pair to be unique. It does not treat header `size` as an upper bound on unique IDs. Ordinary numeric validation is separate, allowing a consensus map to be assembled before its column descriptions are complete.

## Validation and provenance

`tests/features.rs` contains source-derived containment, hull outline/compression, feature enclosure, consensus means and decharge/adduct examples, plus independent boundary, duplicate-key, running-charge-tie, signed-intensity, signed-charge, ID, atomicity, overflow, range and stable-selection checks. The decharge golden uses the proton constant explicitly: the historical C++ test called the neutral hydrogen atom mass `proton_mass`, while the implementation subtracts the actual proton constant. This avoids encoding that test's loose numerical tolerance as a scientific identity.

Validation commands:

```sh
cargo test --offline --no-default-features --test features
cargo clippy --offline --no-default-features --lib --test features -- -D warnings
```

No C++ build was performed. Numerical tests were derived from pinned source and class tests, not from executing the C++ library. Source implementations are under `src/openms/source/` in the pinned reference:

| Source path | SHA-256 |
| --- | --- |
| `DATASTRUCTURES/ConvexHull2D.cpp` | `75181e622beebcde6020cafa816363342cd3e9ed9e18d3e5f8924cac3b1572a1` |
| `KERNEL/BaseFeature.cpp` | `b0166d47181ce74ced4937a112da94788cbc51eb57f88a426a1e62578458151d` |
| `KERNEL/Feature.cpp` | `5b6ea1b5ecb99d656213c02ebfe4e9a89f5216bd6c7b535ab8d24a0bd7579f72` |
| `KERNEL/FeatureHandle.cpp` | `56a66eea0fa16b56e8c88bd9818fccfd94e2872866537c21f7f186e55579c283` |
| `KERNEL/FeatureMap.cpp` | `210ebf2fbe17259af5034b78855d0f57e16becf525aae370217a044f0c926a27` |
| `KERNEL/ConsensusFeature.cpp` | `d604cb44c4964c3f6d8526ebe95bbee10d1d9fb73115833d8a97c4aab7542c61` |
| `KERNEL/ConsensusMap.cpp` | `bc5fa8b64feb192eda86f1dd6606aba6108e97fc5cdcc7b600574189bc166449` |

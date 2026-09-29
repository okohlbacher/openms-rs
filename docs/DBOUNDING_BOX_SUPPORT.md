# DBoundingBox: the source template

`DATASTRUCTURES/DBoundingBox.h` at core `bc9cc12` (header-only) is ported as
`DBoundingBox<const D: usize>` in `src/data_structures/dbounding_box.rs`, with
the aliases `DBoundingBox1` and `DBoundingBox2`. It wraps the already-ported
`DIntervalBase<D>` (`docs/DPOSITION_SUPPORT.md`) and reaches every base
accessor and mutator through `Deref`/`DerefMut`, as the derived C++ class does.

`kernel::geometry::BoundingBox2D` remains the validated two-dimensional box the
feature hulls use (it rejects inverted and non-finite corners and represents
emptiness as `Option`); `ConvexHull2D::map_points` yields `DBoundingBox1`
values, and `ConvexHull2D::source_hash64` hashes them.

## API mapping

| Source | Rust |
|---|---|
| `template <UInt D> class DBoundingBox` | `DBoundingBox<D>`, `DBoundingBox1`, `DBoundingBox2` |
| `enum {DIMENSION = D}` | `DBoundingBox::DIMENSION` |
| `Base`, `PositionType`, `CoordinateType` | `DIntervalBase<D>`, `DPosition<D>`, `f64` |
| `DBoundingBox()` | `Default` (the empty sentinel, min `f64::MAX`, max `f64::MIN`) |
| copy constructor, `operator=`, destructor | `Copy`, `Clone`, assignment, drop |
| `operator=(const Base&)` | `From<DIntervalBase<D>>` (and back) |
| `DBoundingBox(minimum, maximum)` | `new` (normalised per dimension, as the base) |
| `enlarge(const PositionType&)`, `enlarge(x, y)` | `enlarge`, `enlarge_xy` |
| `operator==(const DBoundingBox&)`, `operator==(const Base&)` | `PartialEq`, `PartialEq<DIntervalBase<D>>` |
| `encloses(const PositionType&)`, `encloses(x, y)` | `encloses`, `encloses_xy` |
| `intersects(const DBoundingBox&)` | `intersects` |
| `isEmpty()` | `is_empty` |
| `operator<<` | `Display` |
| `std::hash<DBoundingBox<D>>` | `Hash`; `source_hash64` |
| protected `min_`, `max_` of the base | the base accessors through `Deref` |

## Preserved source conventions

- Margins are inclusive in `encloses` and `intersects`; a NaN coordinate never
  rejects, so a position with a NaN is enclosed.
- `is_empty` is true when any dimension has `max <= min`, so a zero-width box
  is empty here but not for the base's `is_empty`, which the source hides the
  same way.
- `enlarge` of the empty sentinel by one point gives the degenerate box there;
  a NaN coordinate changes nothing.
- `Display` writes the four lines `--DBOUNDINGBOX BEGIN--`, `MIN --> ...`,
  `MAX --> ...`, `--DBOUNDINGBOX END--`, coordinates through
  `precisionWrapper` (`2` prints as `2.0`).
- `source_hash64` folds `hash_float` of every minimum and then every maximum
  coordinate into a zero seed, so `-0.0` and `0.0` hash alike.

## Native differences

- The coordinate type is `f64` only, as for `DPosition`.
- `Hash` feeds `source_hash64` to the hasher; there is no `Eq` (floats).

## Checked boundaries and evidence

`tests/dbounding_box.rs` ports `DBoundingBox_test.cpp` at `bc9cc12` section by
section (tier 3) and checks the source hash against values recomputed in Python
from the header's recurrence over little-endian bytes (tier 4).
`tests/convex_hull_2d.rs` covers `getMapPoints` and the hull hash. Every
operation is constant work over `D` coordinates.

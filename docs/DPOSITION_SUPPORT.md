# DPosition, DIntervalBase and DRange: native header equivalents

`DPosition<D>`, `DIntervalBase<D>` and `DRange<D>` in
[`src/data_structures/dposition.rs`](../src/data_structures/dposition.rs),
[`dinterval.rs`](../src/data_structures/dinterval.rs) and
[`drange.rs`](../src/data_structures/drange.rs) provide native equivalents for
the complete public surface of core SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4`
`DATASTRUCTURES/DPosition.h`, `DIntervalBase.h` and `DRange.h`. All three are
header-only templates; their `.cpp` files are empty or instantiate default
objects. This describes value-type behavior, not C++ ABI, stream state or
hash-digest compatibility.

The dimension is a const generic (`DPosition<const D: usize>`), with aliases
`DPosition1`/`DPosition2`, `DIntervalBase1`/`DIntervalBase2` and
`DRange1`/`DRange2` for the two instantiations the SDK exports. The coordinate
type is fixed to `f64`; the source's `TCoordinateType` parameter and its
`Int`/`Int64`/`char` test instantiations are not ported (see differences). Every
operation runs over one or two fixed-size arrays, so no bounded-work ceiling
(`MAX_ITEMS`/preflight) applies to this group.

## API mapping

### `DPosition.h`

| Source member | Native representation |
| --- | --- |
| `CoordinateType`, `DataType`, `DIMENSION`, STL typedefs | `f64`, `[f64; D]`, `DPosition::DIMENSION`; iterator typedefs are `slice::Iter`/`IterMut`/`array::IntoIter` |
| `DPosition()` | `DPosition::new()`, `Default` (all zero) |
| `DPosition(CoordinateType x)` | `DPosition::filled(x)` |
| `DPosition(x, y)` (D == 2 `static_assert`) | `DPosition::<2>::xy(x, y)`, restricted by the `impl DPosition<2>` block |
| `DPosition(x, y, z)` (D == 3) | `DPosition::<3>::xyz(x, y, z)` |
| Copy/move constructors, `operator=`, destructor | `Copy`, `Clone`, assignment, ordinary drop |
| `swap` | `std::mem::swap` |
| `abs()` (in place, returns `*this`) | `abs(self) -> Self` on the `Copy` value |
| `operator[]` const and mutable | `Index<usize>`/`IndexMut<usize>` (panic out of range) plus checked `get`/`get_mut` (`Option`); the array is also public as `coordinates` |
| `getX/getY/setX/setY` (D == 2 precondition) | `x()/y()/set_x()/set_y()` on `impl DPosition<2>` |
| `operator==`/`!=` | `PartialEq` (derived, ordinary float equality) |
| `operator<`/`<=`/`>`/`>=` (lexicographic) | `PartialOrd` derived over `[f64; D]`, which is lexicographic from dimension 0 |
| `spatiallyLessEqual`/`spatiallyGreaterEqual` | `spatially_less_equal`/`spatially_greater_equal` |
| `operator+`, `+=`, `operator-` (binary), `-=` | `Add`, `AddAssign`, `Sub`, `SubAssign` |
| Unary `operator-` | `Neg` |
| Member `operator*(const DPosition&)` (inner product) | `Mul<DPosition<D>>` with `Output = f64`, and `dot(&other)` |
| `operator*=`/`operator/=` scalar | `MulAssign<f64>`/`DivAssign<f64>` |
| Free `operator*(pos, s)`, `operator*(s, pos)`, `operator/(pos, s)` | `Mul<f64> for DPosition`, `Mul<DPosition> for f64`, `Div<f64>` |
| Static `size()` | `const fn size() -> usize` |
| `clear()` | `clear()` |
| Static `zero()`, `minPositive()`, `minNegative()`, `maxPositive()` | `const fn zero()`, `min_positive()` (`f64::MIN_POSITIVE`), `min_negative()` (`f64::MIN`), `max_positive()` (`f64::MAX`) |
| `begin()/end()` const and mutable | `iter()`, `iter_mut()`, `as_slice()`, `IntoIterator` for `&`, `&mut` and by value |
| Free `operator<<` (space-separated `precisionWrapper`) | `Display`, space-separated, precision option forwarded |
| `std::hash<DPosition>` (`hash_float`, signed-zero normalization) | `Hash`, zero-normalized, no `Eq` |
| Not in source but added | `From<[f64; D]>`, `From<DPosition<D>> for [f64; D]` |

### `DIntervalBase.h` (`Internal::DIntervalBase<D>`)

| Source member | Native representation |
| --- | --- |
| `DIMENSION`, `PositionType`, `CoordinateType` | `DIntervalBase::DIMENSION`, `DPosition<D>`, `f64` |
| `DIntervalBase()` (empty; comment says "corners at infinity") | `Default`, `empty()`: min `f64::MAX`, max `f64::MIN` (finite extrema, not infinities) |
| Copy/move constructors, `operator=`, destructor | `Copy`, `Clone`, assignment, ordinary drop |
| `DIntervalBase(minimum, maximum)` (normalizes) | `new(minimum, maximum)` |
| Protected pair constructor (no normalization, for statics) | Private field construction inside `empty()`/`zero()` |
| `minPosition()`/`maxPosition()` | `min_position()`/`max_position()` returning `&DPosition<D>` |
| `setMin` / `setMax` (adjust the other corner; `@note` carried) | `set_min` / `set_max` |
| `setMinMax` (normalizes) | `set_min_max` |
| `assign<D2>` (copies `min(D, D2)` dimensions) | `assign<const D2>(&DIntervalBase<D2>)` |
| `operator==`/`!=` | `PartialEq` |
| `operator+`, `+=`, `-`, `-=` with a position | `Add/AddAssign/Sub/SubAssign<DPosition<D>>` |
| `clear()` | `clear()` |
| `isEmpty()` (`min == max` is not empty) | `is_empty()` |
| `isEmpty(UInt dim)` | `is_empty_dim(dim) -> Result<bool>` |
| `setDimMinMax(dim, DIntervalBase<1>)` | `set_dim_min_max(dim, &DIntervalBase<1>) -> Result<()>` |
| `center()`, `diagonal()` | `center()`, `diagonal()` |
| Static `empty`, `zero` | `const fn empty()`, `const fn zero()` |
| `minX/minY/maxX/maxY`, `setMinX/setMinY/setMaxX/setMaxY`, `width`, `height` | Same names in snake case; instantiating them with too few dimensions is a compile-time error |
| Protected `min_`/`max_`, `normalize_()` | Private fields, private `normalize()`; `DRange` reaches the corners through `pub(super) corners_mut` |
| `friend` declaration for other dimensions | Not needed; `assign` reads the private fields of any `DIntervalBase<D2>` within the module |
| Free `operator<<` | `Display`: `--DIntervalBase BEGIN--`, `MIN --> …`, `MAX --> …`, `--DIntervalBase END--`, each newline-terminated |
| No `std::hash` in source | `Hash` added (min then max, zero-normalized) so `DRange` can delegate |

### `DRange.h`

| Source member | Native representation |
| --- | --- |
| `DIMENSION`, `Base`, `PositionType`, `CoordinateType` | `DRange::DIMENSION`, `DIntervalBase<D>`, `DPosition<D>`, `f64` |
| `enum DRangeIntersection { Disjoint, Intersects, Inside }` | `DRangeIntersection` with the same three variants, in that order |
| `DRange()` (comment "all coordinates zero", body calls empty base ctor) | `Default`, `empty()`: the empty sentinel, as the source actually does |
| `DRange(lower, upper)` | `new(lower, upper)` (normalizes) |
| Copy/move constructors, `operator=(const DRange&)`, destructor | `Copy`, `Clone`, assignment, ordinary drop |
| `DRange(const Base&)`, `operator=(const Base&)` | `From<DIntervalBase<D>>`; reverse `From<DRange<D>> for DIntervalBase<D>`; `base()`/`base_mut()` views |
| `DRange(minx, miny, maxx, maxy)` (D == 2) | `DRange::<2>::xy(min_x, min_y, max_x, max_y)` |
| Inherited base members (`minPosition` … `height`, `empty`, `zero`, `+`/`-`) | Explicitly delegated methods of the same names; `Add/Sub` keep the `DRange` type where the source returns the base |
| `operator==(const DRange&)`, `operator==(const Base&)` | `PartialEq`, `PartialEq<DIntervalBase<D>>` (and the symmetric impl) |
| `encloses(const PositionType&)` | `encloses(&position)` |
| `encloses(x, y)` (2D) | `DRange::<2>::encloses_xy(x, y)` |
| `united` | `united(&other)` |
| `intersects` | `intersects(&range) -> DRangeIntersection` |
| `isIntersected` | `is_intersected(&range)` |
| `extend(double factor)` (throws `InvalidParameter` below zero) | `extend_by_factor(factor) -> Result<&mut Self>` with `Error::InvalidValue` |
| `extend(PositionType addition)` | `extend_by(addition) -> &mut Self` |
| `ensureMinSpan` | `ensure_min_span(min_span) -> &mut Self` |
| `swapDimensions` (D == 2) | `DRange::<2>::swap_dimensions() -> &mut Self` |
| `pullIn(DPosition<D>& point)` (in/out parameter) | `pull_in(point) -> DPosition<D>` |
| Free `operator<<` | `Display`: `--DRANGE BEGIN--`, `MIN --> …`, `MAX --> …`, `--DRANGE END--`, newline-terminated |
| `std::hash<DRange>` | `Hash` (all minima then all maxima, zero-normalized), no `Eq` |

Members named in the work-package brief that do not exist in the pinned
headers and were therefore not ported: `squaredDistance`/`distance` on
`DPosition`; `contains`/`intersect`/`union` on `DIntervalBase`;
`minPositive`/`maxNegative`, `nonNegative` and `intersection` on `DRange`. The
pinned `DRange` provides `united`, `intersects`, `isIntersected`, `encloses`
and `pullIn` instead.

## Preserved source conventions

- **Empty sentinel.** `DIntervalBase::empty` is built through the protected
  pair constructor with `min = maxPositive()` (`+DBL_MAX`) and `max =
  minNegative()` (`-DBL_MAX`), i.e. finite extrema rather than infinities
  and an inverted pair that violates the class invariant on purpose. `clear()`
  restores it, `isEmpty()` is equality with it, and `isEmpty(dim)` checks the
  per-dimension pair `(+DBL_MAX, -DBL_MAX)`. `min == max` is never empty.
- **Normalization.** The corner constructor, `setMinMax` and the four-coordinate
  `DRange` constructor swap `min[i]`/`max[i]` where inverted; `setMin`/`setMax`
  and `setMinX`… move the *other* corner instead; `assign` and `setDimMinMax`
  copy without normalizing. All of this is transcribed.
- **Lexicographic ordering** of `DPosition` from dimension 0, and the source's
  early-exit `spatially*` loops (a NaN pair does not fail them).
- **Half-open enclosure**: `>= min` and `< max` per dimension; NaN coordinates
  fail neither test and are enclosed.
- **`intersects` decision order** and **`united`** are transcribed line for
  line, including the source quirk that uniting two empty ranges yields the
  universal range `[-DBL_MAX, +DBL_MAX]` (the sentinel's inverted corners are
  normalized away by `setMinMax`). Uniting a non-empty range with the sentinel
  returns that range.
- **`extend`** arithmetic: `(max - min) / 2 * (factor - 1)` per side, and the
  additive form halves `addition`, translates, then collapses any inverted
  dimension to its center. The source's `@param` remark that invalid results
  "are not fixed automatically" is contradicted by that collapse; the port
  follows the implementation and says so at the item.
- **`pullIn`** transcribes `std::max(min, std::min(point, max))`, so a NaN
  coordinate becomes `min[i]` (documented; `f64::clamp` would differ).
- **Hashing** normalizes `-0.0` to `+0.0` before hashing bits, as
  `hash_float` does, in the same coordinate order.
- **`DRange()`** yields the empty sentinel: the source comment promises zeros,
  the body does not (C++ issue candidate below).

## Native differences

- `TCoordinateType` is fixed to `f64`. The `[EXTRA]` int/char test sections are
  ported with the same integral values as doubles; the `Int64` half of the
  `abs()` section (a value not representable in a double) has no counterpart.
- `DPosition::abs` returns a new value instead of mutating through a
  reference, because the type is `Copy`.
- `DPosition` `operator[]` has an `OPENMS_PRECONDITION` that is checked only in
  debug builds; `Index` always panics out of range, and `get`/`get_mut`
  offer non-panicking access.
- The 2D convenience accessors (`min_x` … `height`) are defined for every `D`
  as in the source, but instantiating them with too few dimensions is a
  compile-time error (inline `const` assertion) rather than out-of-bounds
  array access.
- `isEmpty(dim)` and `setDimMinMax(dim, …)` index unchecked in the source;
  here `dim >= D` returns `Error::InvalidValue` and leaves the value unchanged.
- `extend_by_factor` rejects NaN and infinite factors in addition to the
  source's negative check, before any mutation; the source would propagate
  NaN/infinity into the corners.
- `DRange` composes `DIntervalBase` instead of deriving from it. All inherited
  members are delegated explicitly; conversions go through `From` in both
  directions, and `PartialEq<DIntervalBase<D>>` replaces
  `operator==(const Base&)`. `DRange + DPosition` returns a `DRange` rather
  than the source's base-class result.
- `pullIn`'s in/out parameter is a return value; the source `extend` overload
  set is `extend_by_factor`/`extend_by`.
- No `Eq`, so none of the types is a `HashMap`/`HashSet` key: NaN makes float
  equality non-reflexive and `Eq` would assert a false contract. The source's
  `unordered_set`/`unordered_map` test sections are represented by
  equality/hash agreement, as for `Peak1D`/`Peak2D`
  ([`KERNEL_VALUE_TRAITS_SUPPORT.md`](KERNEL_VALUE_TRAITS_SUPPORT.md)). No
  numeric digest compatibility with the source FNV-1a combination is promised.
- `Display` uses Rust float formatting (round-trippable by default, or the
  caller's `{:.N}`) instead of `precisionWrapper`'s 15 significant digits and
  `std::endl`; the labels, separators and newline placement are the source's.
- `DIntervalBase` gains a `Hash` impl and `From<[f64; D]>` conversions that the
  source lacks; both are marked as native at the item.

Doxygen accounting: every `@note` (`setMin`/`setMax` corner adjustment,
`isEmpty` "min==max is NOT empty") and the `@invariant` are carried; the single
`@exception`-equivalent (`extend`'s `InvalidParameter`) maps to
`Error::InvalidValue`; `@param` constraints (`factor` in `[0, inf)`, negative
`addition` allowed) are carried; the `DRange` class description is carried
verbatim in substance. There are no `@see` targets in the three headers. The
related in-crate 2D point type used by kernel geometry is
`crate::kernel::geometry::Point2D`; the peak types keep their own
`[f64; 2]` positions and are not rebased on `DPosition` in this work package.

## Checked boundaries and evidence

Evidence tier 3 (source review): the class-test literals of all 48 + 35 + 21
`START_SECTION`s are transcribed into
[`tests/dposition.rs`](../tests/dposition.rs) (51 tests, including two
`should_panic` splits of the precondition sections),
[`tests/dinterval.rs`](../tests/dinterval.rs) (36) and
[`tests/drange.rs`](../tests/drange.rs) (22), each test citing its source line.
Native additions check the sentinel values, dimension errors, NaN routes,
stream output and hash normalization.
[`dposition_provenance.json`](../tests/data/dposition_provenance.json) pins the
nine inspected files with SHA-256 and thirteen source anchors. No C++ build or
execution was used; the `.cpp` files contain no logic to execute.

Self-audit:

- `DPosition.h`: 56 public members read (8 type definitions, 9 constructors/
  assignments/destructor, `swap`, `abs`, 6 accessors, 8 comparisons, 8
  arithmetic operators, `size`, `clear`, 4 statics, 4 iterator accessors, 4
  free operators, `std::hash`); 56 mapped, 0 deferred; 48/48 sections ported.
- `DIntervalBase.h`: 40 public members read (3 type definitions, 6
  constructors/assignment/destructor, 6 corner accessors/mutators, 6 operators,
  4 emptiness/dimension members, `center`, `diagonal`, 2 statics, 10 2D
  convenience members, free `operator<<`); 40 mapped, 0 deferred; 35/35
  sections ported.
- `DRange.h`: 28 public members read (5 type definitions incl. the enum, 9
  constructors/assignments/destructor, 12 predicates and mutators, free
  `operator<<`, `std::hash`) plus the inherited surface; 28 mapped, 0 deferred;
  21/21 sections ported.

C++ issue candidates (unconfirmed, source review only): `DRange()` documents
"all coordinates zero" but constructs the empty sentinel; `DRange::united`
of two empty ranges returns the universal range; `DRange::extend(addition)`'s
`@param` says inverted results are not fixed while the body collapses them.

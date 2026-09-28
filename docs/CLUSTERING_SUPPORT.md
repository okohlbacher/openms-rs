# Clustering support (`ML/CLUSTERING`, `DistanceMatrix`, `BinaryTreeNode`, `GridFeature`, `QTCluster`)

Source revision: OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

| Source header | Rust | Status |
|---|---|---|
| `DATASTRUCTURES/DistanceMatrix.h` | `src/data_structures/distance_matrix.rs` | complete |
| `DATASTRUCTURES/BinaryTreeNode.h` | `src/data_structures/binary_tree_node.rs` | complete |
| `ML/CLUSTERING/ClusterFunctor.h` | `src/math/ml/clustering/cluster_functor.rs` | complete |
| `ML/CLUSTERING/SingleLinkage.h` | `src/math/ml/clustering/single_linkage.rs` | complete |
| `ML/CLUSTERING/CompleteLinkage.h` | `src/math/ml/clustering/complete_linkage.rs` (+ `linkage.rs`) | complete |
| `ML/CLUSTERING/AverageLinkage.h` | `src/math/ml/clustering/average_linkage.rs` (+ `linkage.rs`) | complete |
| `ML/CLUSTERING/ClusterAnalyzer.h` | `src/math/ml/clustering/cluster_analyzer.rs` | complete |
| `ML/CLUSTERING/ClusterHierarchical.h` | `src/math/ml/clustering/cluster_hierarchical.rs`, `src/analysis/spectra_clustering.rs` | complete |
| `ML/CLUSTERING/EuclideanSimilarity.h` | `src/math/ml/clustering/euclidean_similarity.rs` | complete |
| `ML/CLUSTERING/ClusteringGrid.h` | `src/math/ml/clustering/clustering_grid.rs` | complete |
| `ML/CLUSTERING/GridBasedCluster.h` | `src/math/ml/clustering/grid_based_cluster.rs` | complete |
| `ML/CLUSTERING/GridBasedClustering.h` | `src/math/ml/clustering/grid_based_clustering.rs` | complete |
| `ML/CLUSTERING/HashGrid.h` | `src/math/ml/clustering/hash_grid.rs` | partial (iteration order, float `hash_value`) |
| `DATASTRUCTURES/GridFeature.h` | `src/analysis/mapmatching/grid_feature.rs` | complete |
| `DATASTRUCTURES/QTCluster.h` | `src/analysis/mapmatching/qt_cluster.rs` | partial (hash-map order: neighbour order, summation order) |
| `DATASTRUCTURES/KDTree.h` | none | not ported (see Backlog) |

## Placement

`math` is a sink in the module graph. The hierarchical and grid clustering
code needs only `math`, `data_structures` (`DistanceMatrix`, `DPosition`,
`DIntervalBase`) and `concept` (`ProgressReporter`), none of which reaches
`math`, so it lives in `src/math/ml/clustering/`. Three pieces need types above
`math` and live beside their consumers:

- the binned-spectrum overload of `ClusterHierarchical::cluster`
  (`BinnedSpectrum`, `BinnedSpectrumCompareFunctor`, `MSSpectrum`) is
  `src/analysis/spectra_clustering.rs`;
- `GridFeature` and `QTCluster` (`BaseFeature`, `AASequence`,
  `PeptideIdentification`) are in `src/analysis/mapmatching/`, beside the future
  `QTClusterFinder` port. They are `DATASTRUCTURES` headers, but
  `data_structures -> kernel` would close a cycle (`kernel -> metadata ->
  data_structures`).

`DistanceMatrix` and `BinaryTreeNode` are `DATASTRUCTURES` headers with no
upward dependency and live in `src/data_structures/`.

## API mapping

### DistanceMatrix<Value> (`DistanceMatrix<V>`)

| Source | Rust |
|---|---|
| `value_type`, `SizeType`, `ValueType` typedefs | `V`, `usize` |
| `DistanceMatrix()` | `new()`, `Default` |
| `DistanceMatrix(dimensionsize, value)` | `with_value(dimension, value) -> Result` |
| copy constructor | `Clone` |
| `~DistanceMatrix()` | drop |
| `operator()(i,j)` const and non-const, `getValue(i,j)` const and non-const | `get(i, j) -> Result<V>` |
| `setValue` | `set(i, j, v) -> Result<()>` |
| `setValueQuick` | `set_quick(i, j, v) -> Result<()>` |
| `clear` | `clear()` |
| `resize(dimensionsize, value)` | `resize(dimension, value) -> Result<()>` |
| `reduce(j)` | `reduce(j) -> Result<()>` |
| `dimensionsize` | `dimension_size()` |
| `updateMinElement` | `update_min_element() -> Result<()>` |
| `operator==` | `PartialEq` |
| `getMinElementCoordinates` | `min_element_coordinates() -> Result<(usize, usize)>` |
| private `operator=` ("unsafe") | not ported (private) |
| `operator<<` | `Display` |
| protected `matrix_`, `init_size_`, `dimensionsize_`, `min_element_` | private fields (`init_size_` is not needed: rows are owned vectors) |

### BinaryTreeNode

| Source | Rust |
|---|---|
| `BinaryTreeNode(i, j, x)` | `new(i, j, x)` |
| copy constructor, `operator=`, destructor | `Clone`, `Copy`, drop |
| `left_child`, `right_child`, `distance` | public fields |
| private default constructor | not ported (private) |

### ClusterFunctor

| Source | Rust |
|---|---|
| `InsufficientInput(file, line, function, message)` and destructor | `insufficient_input(message) -> Error` (`Error::InvalidValue`, message prefixed `ClusterFunctor::InsufficentInput`, the source's spelling) |
| default/copy constructor, destructor, `operator=` | trait, no state |
| `operator()(original_distance, cluster_tree, threshold = 1)` | `ClusterFunctor::cluster(&mut matrix, threshold) -> Result<Vec<BinaryTreeNode>>`, `cluster_with_progress(..., &mut ProgressReporter)`; default threshold `DEFAULT_CLUSTER_THRESHOLD` |

### SingleLinkage, CompleteLinkage, AverageLinkage

| Source | Rust |
|---|---|
| default constructor, copy constructor, destructor, `operator=` | `new()`, `Default`, `Clone`, `Copy` |
| `operator()(original_distance, cluster_tree, threshold)` | `ClusterFunctor::cluster` / `cluster_with_progress` |
| `ProgressLogger` base | the `ProgressReporter` passed to `cluster_with_progress` |

### ClusterAnalyzer

| Source | Rust |
|---|---|
| default/copy constructor, destructor | `new()`, `Default`, `Clone`, `Copy` |
| `averageSilhouetteWidth` | `average_silhouette_width` |
| `dunnIndices(tree, original, tree_from_singlelinkage = false)` | `dunn_indices(tree, original, bool)` |
| `cohesion` | `cohesion` |
| `averagePopulationAberration` | `average_population_aberration` |
| `cut(k, tree, clusters)` | `cut(k, tree) -> Result<Vec<Vec<usize>>>` |
| `cut(k, tree, subtrees)` | `cut_subtrees(k, tree) -> Result<Vec<Vec<BinaryTreeNode>>>` |
| `newickTree(tree, include_distance = false)` | `newick_tree(tree, Option<&dyn Fn(f32) -> String>)` (see native differences) |
| private `operator=` | not ported (private) |
| free `compareBinaryTreeNode` | `compare_binary_tree_node` |

### ClusterHierarchical, UnnormalizedComparator

| Source | Rust |
|---|---|
| default/copy constructor, destructor | `new()`, `Default`, `Clone`, `Copy` |
| template `cluster<Data, SimilarityComparator>` | `cluster<D, S: SimilarityValue>(data, comparator, clusterer, &mut matrix)` |
| `cluster(PeakSpectrum..., BinnedSpectrumCompareFunctor, sz, sp, offset, ...)` | `analysis::spectra_clustering::cluster_binned_spectra` |
| `getThreshold` / `setThreshold` | `threshold()` / `set_threshold()` |
| `UnnormalizedComparator(...)` and destructor | `unnormalized_comparator(message) -> Error` (declared, never raised by the source either) |

### EuclideanSimilarity

| Source | Rust |
|---|---|
| default/copy constructor, destructor, `operator=` | `new()`, `Default`, `Clone`, `Copy` |
| `operator()(a, b)` | `similarity(a, b) -> Result<f32>` |
| `operator()(c)` | `self_similarity(c) -> Result<f32>` |
| `setScale` | `set_scale` |

### ClusteringGrid

| Source | Rust |
|---|---|
| `CellIndex`, `Point` typedefs | `CellIndex = (i32, i32)`, `Point = DPosition<2>` |
| constructor | `new(x_lines, y_lines) -> Result` |
| `getGridSpacingX` / `getGridSpacingY` | `grid_spacing_x()` / `grid_spacing_y()` |
| `addCluster`, `removeCluster`, `removeAllClusters` | `add_cluster`, `remove_cluster`, `remove_all_clusters` |
| `getClusters` | `clusters(cell) -> Vec<i32>` |
| `getIndex` | `index(&point) -> Result<CellIndex>` |
| `isNonEmptyCell` | `is_non_empty_cell` |
| `getCellCount` | `cell_count` |

### GridBasedCluster

| Source | Rust |
|---|---|
| `Point`, `Rectangle` typedefs | `Point = DPosition<2>`, `Rectangle = DIntervalBase<2>` (+ `enlarge`) |
| five-argument constructor | `new(centre, box, points, a, b)` |
| three-argument constructor | `without_properties(centre, box, points)` |
| `getCentre`, `getBoundingBox`, `getPoints`, `getPropertyA`, `getPropertiesB` | `centre`, `bounding_box`, `points`, `property_a`, `properties_b` |
| `operator<`, `operator>`, `operator==` | `PartialOrd`, `PartialEq` (centre y only) |

### GridBasedClustering<Metric>, MinimumDistance

| Source | Rust |
|---|---|
| `MinimumDistance(cluster, neighbour, distance)` | `MinimumDistance::new` |
| `getClusterIndex`, `getNearestNeighbourIndex` | `cluster_index`, `nearest_neighbour_index` |
| `MinimumDistance::operator<`, `>`, `==` | `PartialOrd`, `PartialEq` (distance only) |
| private default constructor | not ported (private) |
| `Point`, `Rectangle`, `CellIndex`, `MultisetIterator`, `NNIterator` typedefs | `grid_based_cluster::{Point, Rectangle}`, `clustering_grid::CellIndex`; the iterator typedefs are internal |
| seven-argument constructor | `new(metric, x, y, a, b, x_lines, y_lines) -> Result` |
| five-argument constructor | `without_properties(metric, x, y, x_lines, y_lines) -> Result` |
| `cluster` | `cluster()`, `cluster_with_progress()` |
| `extendClustersY` | `extend_clusters_y() -> Result` |
| `removeSmallClustersY` | `remove_small_clusters_y` |
| `getResults` | `results() -> &BTreeMap<i32, GridBasedCluster>` |
| private `init_`, `mergeVeto_`, `findNearestNeighbour_`, `eraseMinDistance_` | private `init`, `merge_veto`, `find_nearest_neighbour`, `erase_min_distance` |

### HashGrid<Cluster>

| Source | Rust |
|---|---|
| `ClusterCenter`, `CellIndex`, `CellContent`, `Grid`, `key_type`, `mapped_type`, `value_type` | `ClusterCenter = DPosition<2>`, `CellIndex = [i64; 2]`, `CellContent<C> = Vec<(ClusterCenter, C)>`; the grid map is private |
| `Iterator`, `ConstIterator` (with `index()`) | `iter()`, `iter_mut()` yielding `(CellIndex, pair)` |
| `cell_dimension`, `grid_dimension` public members | `cell_dimension()`, `grid_dimension()` |
| constructor | `new(cell_dimension)` |
| `insert` | `insert(key, value) -> Result<(CellIndex, usize)>` |
| `erase(iterator)` | `erase_at(cell, position) -> Result` |
| `erase(key)` | `erase(&key) -> Result<usize>` |
| `clear`, `empty`, `size` | `clear`, `is_empty`, `len` |
| `begin`/`end` (const and not) | `iter` / `iter_mut` |
| `grid_begin`/`grid_end` (const and not) | `grid_iter` / `grid_iter_mut` |
| `grid_at` (const and not) | `grid_at` / `grid_at_mut` |
| `grid_find` (const and not) | `grid_find` / `grid_find_mut` |
| `cellIndexAtClusterCenter` | `cell_index_at_cluster_center` |
| free `hash_value(DPosition<N, T>)` | `hash_value(&[u64])` for integral `T`; the `double` instantiation is not ported |

### GridFeature

| Source | Rust |
|---|---|
| constructor, destructor | `GridFeature::new(&feature, map_index, feature_index)`, drop |
| `getFeature`, `getMapIndex`, `getFeatureIndex`, `getID`, `getAnnotations`, `getRT`, `getMZ` | `feature`, `map_index`, `feature_index`, `id`, `annotations`, `rt`, `mz` |
| `std::set<AASequence>` | `BTreeSet<SourceOrderedSequence>` |

### QTCluster

| Source | Rust |
|---|---|
| `NeighborList`, `NeighborMapMulti`, `Neighbor`, `NeighborMap`, `Element`, `Elements` | `NeighborList<'a>`, `Neighbor<'a>`, `Element<'a>`, `Elements<'a>`; the maps are private `BTreeMap`s |
| `BulkData(center, num_maps, max_distance, x, y, id)` | `BulkData::new` |
| `QTCluster(BulkData*, bool)` | `QTCluster::new(BulkData, bool)` (owns the body) |
| deleted default constructor; copy/move constructors and assignments; destructor | none; `Clone`; drop |
| `getCenterPoint`, `getId`, `getCenterRT`, `getCenterMZ`, `getXCoord`, `getYCoord`, `size` | `center_point`, `id`, `center_rt`, `center_mz`, `x_coord`, `y_coord`, `size` |
| `operator<` | `PartialOrd` (current quality) |
| `add`, `getElements`, `update`, `getQuality`, `getCurrentQuality`, `getAnnotations` | `add`, `elements`, `update`, `quality`, `current_quality`, `annotations` |
| `setInvalid`, `isInvalid`, `initializeCluster`, `finalizeCluster`, `getAllNeighbors` | `set_invalid`, `is_invalid`, `initialize_cluster`, `finalize_cluster`, `all_neighbors` |
| private `computeQuality_`, `optimizeAnnotations_`, `makeSeqTable_`, `recomputeNeighbors_` | private `compute_quality`, `optimize_annotations`, `make_seq_table`, `recompute_neighbors` |
| private `finalized_` | readable through `is_finalized()` (native) |

## Preserved source conventions

- **`DistanceMatrix` layout and minimum tracking.** Row `i` stores the `i` cells
  below the diagonal; `reduce` drops row and column `j` exactly as the source's
  copy loop does; `update_min_element` scans from `(1, 0)` with strict `<`, so
  ties keep the earliest row and column; `set` follows the pinned source after
  the #9488 fix.
- **Linkage.** `SingleLinkage` is SLINK with the pre-tree sorted by the Release
  build's `std::sort` (`math::source_sort::source_sort_by`), so equal distances
  come out in the source's order, then renumbered. `CompleteLinkage` and
  `AverageLinkage` share one loop (`linkage.rs`) with the source's
  Lance-Williams updates in `f32` and placeholder nodes `(first, lead, -1)`.
  `std::min` is emulated where a NaN could reach it.
- **ClusterAnalyzer's `std::map::operator[]`.** Reads insert default entries and
  loops run to the map's current size; the port does the same, so every loop
  visits the same keys. Merges append then clear, as the source does.
  `dunnIndices` keeps the source's reset of the merged-away cluster's minimum to
  `-1` before comparing it.
- **GridBasedClustering's multiset.** `std::multiset<MinimumDistance>` (distance
  order, a new element after its equals) is an ordered map keyed by distance and
  an insertion number; `lower_bound(0)` is a range from `(0, 0)`. Its two hash
  maps are only searched, never iterated observably. `extendClustersY` keeps the
  source's stale-list pair walk, its `std::map<GridBasedCluster, int>` keyed on
  centre y (equal y collapses onto the first index), and scans only y cell 1.
- **ClusteringGrid.** `getIndex` is `std::lower_bound` with `std::less_equal`
  (`math::libstdcxx::lower_bound`); a NaN coordinate passes the range check and
  lands in cell 0, as in the source.
- **QTCluster.** `std::multimap` insertion order (after equal keys), `std::min`
  argument order in `makeSeqTable_` and the unspecific fold, `best_pos`
  starting at the first sequence (so with no total below `num_maps *
  max_distance` the first sequence is chosen), and the post-finalize quirk: a
  cluster changed by `update` recomputes its quality with no collected
  neighbours and loses them (quality `(max - num_maps*max/num_other)/max`); the
  oracle's `Q1/updated` records it.
- **GridFeature annotation order.** `SourceOrderedSequence` transcribes
  `AASequence::operator<` (`AASequence.cpp:328-380` at the pin): residue count,
  N-terminal modification (none first, then by id), residue by residue the
  one-letter code and then the modification (none first), C-terminal
  modification. Oracle case `Q1` shows the source choosing `CC` over `AAA` on a
  tie, which the crate's own `AASequence` order (residue string first) would
  not.

## Native differences

- Every source exception is an `Error`: `OutOfRange` → `InvalidRange`;
  `InvalidParameter`, `IllegalArgument`, `InvalidValue`, `DivisionByZero`,
  `InsufficientInput`, `OutOfMemory` → `InvalidValue`; `NotImplemented` →
  `Unsupported`.
- Undefined behaviour in the source is an error here: `DistanceMatrix::set`
  after `reduce` removed the tracked minimum's row (recomputed first);
  out-of-range tree indices in `SingleLinkage`'s conversion and in
  `ClusterAnalyzer`; `newickTree` of an empty tree; `extendClustersY` merging two
  clusters with the same centre y (the source erases `map::end()`); a NaN
  nearest-neighbour distance (breaks the multiset order); `HashGrid` cell indices
  of `2^63` or NaN (the source converts undefinedly); `ClusteringGrid` with no
  grid lines; `GridBasedClustering` data vectors of different lengths.
- `DistanceMatrix` allocation is bounded by `MAX_DISTANCE_MATRIX_CELLS` (2^28
  cells) and refused before anything changes; the source allocates unbounded and
  leaves an empty matrix after a failed `resize`.
- Output parameters are return values. `ClusterAnalyzer::cut` returns a fresh
  vector; the source appends to the caller's vector and then sorts, truncates or
  pads over its old content too. `CompleteLinkage`/`AverageLinkage` stop the
  placeholder fill at `n - 1` nodes; the source stops at the caller's vector
  capacity, which can be larger.
- `newick_tree` takes the distance formatter: `math` may not depend on `format`,
  where `StringUtils::toStr(float)` lives (`to_str_f32`). Passing
  `Some(&to_str_f32)` is the source's `include_distance = true`; the oracle
  compares that text byte for byte.
- `HashGrid` iterates cells in ascending index and pairs in insertion order; the
  source iterates in `boost::unordered` hash order, which is not reproduced.
- `QTCluster` keeps neighbours by ascending map index; the source's
  `std::unordered_map` order is not reproduced. `getElements` /
  `getAllNeighbors` order differs, and the quality sums neighbour distances in a
  different order, which can change its last bit with three or more neighbours
  (oracle case `Q2` happens to agree). A cluster owns its `BulkData` instead of
  sharing a pointer. The debug-only `OPENMS_PRECONDITION`s are not checked, as in
  the Release build.
- `SourceOrderedSequence` falls back to the crate's `AASequence` order where the
  source compares two different residue modifications by their registry
  addresses (not reproducible by any other program), and where the source finds
  two unequal sequences equivalent (different terminal modifications sharing an
  id), so that the key's equality stays the sequence's.
- `DistanceMatrix`'s `Display` never colours the diagonal (the source colours
  only a terminal `std::cout`/`std::cerr`).
- `EuclideanSimilarity::self_similarity` returns 1, as the source computes; the
  source documents 0.
- `ClusteringGrid::is_non_empty_cell` never fails; the source documents two
  exceptions its implementation never throws.
- `GridBasedClustering::merge_veto` returns `true` for a veto, as the source's
  code does; its documentation states the opposite.

## Checked boundaries and evidence

- **Tier 1, executed differential** (`tests/ml_clustering_oracle.rs`, fixture
  `tests/data/ml_clustering_oracle.tsv`): the FF2a oracle
  (`../oracle/FF2a-clustering/driver.cpp`, sha256 `995baf50…`) linked against
  the Linux x86_64 **Release** install `openms4-release-bc9cc12-c19e494-174b576`
  (`libOpenMS.so` sha256 `abd4fc99…`) on **ibminode06**, compiled with
  libOpenMS's core flags, run twice with byte-identical output. It covers the
  three linkage functors on nine matrices (2 to 100 elements, with and without
  ties, thresholds 1, 0.5 and 0.25, including `std::sort`'s introsort path),
  every `ClusterAnalyzer` function on single, complete (0.5) and average trees up
  to 40 elements, `GridBasedClustering` on four point sets (continuous,
  quantized with ties, with properties A/B, non-uniform grid, `extendClustersY`
  and `removeSmallClustersY`), two `QTCluster` scenarios (annotation tie,
  post-update quirk, twelve neighbours without IDs) and fifty `HashGrid` cell
  indices. Floating-point values compare by bit pattern.
- **Tier 3, class tests** (`tests/ml_clustering.rs`, `tests/ml_grid_clustering.rs`,
  `tests/mapmatching_qt_cluster.rs`): every section of `DistanceMatrix_test`,
  `BinaryTreeNode_test`, `ClusterFunctor_test`, `SingleLinkage_test`,
  `CompleteLinkage_test`, `AverageLinkage_test`, `ClusterAnalyzer_test`,
  `ClusterHierarchical_test` (both overloads), `EuclideanSimilarity_test`,
  `ClusteringGrid_test`, `GridBasedCluster_test`, `GridBasedClustering_test`,
  `HashGrid_test`, `GridFeature_test` and `QTCluster_test`, with the source's
  tolerances.
- **Tier 4**: error paths listed above, the allocation ceiling, invariants
  (every point in exactly one grid cluster, no repeated property B, `n - 1`
  linkage nodes), `SourceOrderedSequence`'s length-first order.

## Backlog

- **`DATASTRUCTURES/KDTree.h`** (2,208 lines of libkdtree++, no class test).
  Consumer: `ANALYSIS/QUANTITATION/KDTreeFeatureMaps` (`insert`, `size`,
  `clear`, `optimize`, `find_within_range`), reached by
  `AssayGeneratorMetaboSirius`. Its result order and `optimize`'s median split
  (`std::nth_element`) are observable, so a crate (`kiddo`, `kdtree`) must first
  be measured against an oracle of `find_within_range`, or the header ported with
  `nth_element` emulated. Named dependency: that crate evaluation.
- `HashGrid` iteration order and `QTCluster` neighbour order follow Boost's and
  libstdc++'s hash tables; reproducing them needs an emulation of those tables'
  bucket and rehash policies. `QTClusterFinder` is the consumer to measure it
  against.

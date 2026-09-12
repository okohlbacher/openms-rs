# Feature and consensus map support

[`src/kernel/map_operations.rs`](../src/kernel/map_operations.rs) completes
`KERNEL/FeatureMap.h` and `KERNEL/ConsensusMap.h` at Core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The two containers themselves are
declared in [`src/kernel/features.rs`](../src/kernel/features.rs) and their
per-element identification surface in
[`src/kernel/feature_identification.rs`](../src/kernel/feature_identification.rs);
this document supersedes the parts of [FEATURE_SUPPORT.md](FEATURE_SUPPORT.md)
and [FEATURE_IDENTIFICATION_SUPPORT.md](FEATURE_IDENTIFICATION_SUPPORT.md) that
describe the *containers*.

Both headers were ledgered `evidence_requires_review` with no entry in
`core-sdk-reviewed-apis.json`, so this work package is first of all an audit:
the tables below enumerate every public member of both headers by hand. 38 of
the rows are marked **new** — members that had no Rust counterpart before this
work package; the rest already existed and are cited.

The module holds operations only, no container types. It sits beside
`features.rs` for the same reason as `feature_identification.rs` and
`gap_closures.rs`: that file is frozen for this work package, and Rust allows
inherent and trait implementations for a crate-local type from any module of the
same crate. Callers see one type with one API.

## Ownership: identifications attach by reference

The source containers own an `IdentificationData id_data_` member
(`FeatureMap.h:294`, `ConsensusMap.h:376`) and their elements hold iterators
into it. This port deliberately does not embed the graph:

* `IdentificationData` is not `Clone` by design — IDs are owner-tagged, so a
  copy needs `try_clone_with_translation`, not a derived clone.
* Embedding it would therefore strip `Clone` and `PartialEq` from `FeatureMap`
  and `ConsensusMap`, which 106 map-equality assertions in the existing test
  suite depend on, and would cascade into `consensusxml::ReadReport` and
  `FeatureFindingMetaboOutput`.

Every operation that must *resolve* a reference therefore takes the graph, or
the `ReferenceTranslator`, as a parameter — exactly the shape
`feature_identification.rs` established for `BaseFeature` and `Feature`:

| Source | Rust |
| --- | --- |
| `getIdentificationData()` (both overloads, both containers) | not ported: the graph is the caller's |
| `getUnassignedIDMatches()` | `unassigned_id_matches(&graph)` |
| the reference repointing inside the copy constructor, `operator=`, `operator+=`, `appendRows` and `appendColumns` | `update_id_references(&translator)`, called explicitly after `IdentificationData::merge_from` |

The observable consequence is that a stale or foreign reference is an error at
the call site instead of a dangling iterator.

## API mapping — `KERNEL/FeatureMap.h`

Every public member of the header appears below. Members this work package adds
are marked **new**.

### `struct AnnotationStatistics` and the free stream operators

| Source member | Rust counterpart |
| --- | --- |
| `std::vector<Size> states` | **new** `AnnotationStatistics::states() -> &[usize; 4]`. A fixed-length array, so it can never have the wrong length; the source's `SIZE_OF_ANNOTATIONSTATE` sentinel that sized the vector has no counterpart (see `feature_identification.rs`) |
| `AnnotationStatistics()` | **new** `AnnotationStatistics::new()` / `Default` |
| `AnnotationStatistics(const AnnotationStatistics&)` | `Clone` (derived) |
| `AnnotationStatistics& operator=(const AnnotationStatistics&)` | assignment; the source's self-assignment guard is unnecessary for a `Copy` type |
| `bool operator==(const AnnotationStatistics&) const` | `PartialEq` (derived), plus `Eq` and `Hash`, which the source does not provide |
| `AnnotationStatistics& operator+=(BaseFeature::AnnotationState)` | **new** `AnnotationStatistics::add(state) -> Result<()>` (checked) and `impl AddAssign<AnnotationState>` (saturating). The source increments an unchecked `Size` |
| `operator<<(std::ostream&, const AnnotationStatistics&)` | **new** `impl Display for AnnotationStatistics` |
| `operator<<(std::ostream&, const FeatureMap&)` | **new** `impl Display for FeatureMap` |
| — | **new** `AnnotationStatistics::count(state)`, a named accessor with no source counterpart |

### Base classes and the exposed vector

| Source member | Rust counterpart |
| --- | --- |
| `EXPOSED_VECTOR_INTERFACE(Feature)` — `size`, `empty`, `begin`/`end`/`cbegin`/`cend`, `rbegin`/`rend`/`crbegin`/`crend`, `operator[]`, `at`, `front`, `back`, `push_back`, `emplace_back`, `pop_back`, `insert` (5 overloads), `emplace`, `erase` (2), `resize`, `reserve`, `capacity`, `max_size`, `shrink_to_fit`, `assign` (3), `clear`, `swap`, `getData` (2), `operator==`/`!=`/`<`/`<=`/`>`/`>=` | the public `FeatureMap::features: Vec<Feature>` field and the standard `Vec` API |
| `ExposedVector` constructors `(n)`, `(n, val)`, `(begin, end)` | `FeatureMap::from_features(vec![Feature::default(); n])` and friends; covered by the `[EXTRA] ExposedVector Ctor` section port |
| typedefs `ExpVec`, `value_type`, `iterator`, `const_iterator`, `reverse_iterator`, `const_reverse_iterator`, `size_type`, `pointer`, `reference`, `const_reference`, `difference_type` | not ported: `Vec<Feature>`'s own associated types |
| typedefs `RangeManagerContainerType`, `RangeManagerType`, `Iterator`, `ConstIterator`, `ReverseIterator`, `ConstReverseIterator` | not ported: aliases for the above and for `FeatureRanges` |
| `MetaInfoInterface` base | the public `metadata: MetaInfo` field |
| `DocumentIdentifier` base | the public `identifier`, `loaded_file_path` and `loaded_file_type` fields |
| `RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity>` base | `FeatureMap::ranges() -> Result<FeatureRanges>`, computed on demand |
| `UniqueIdInterface` base | the public `unique_id: u64` field plus `HasUniqueId` (`src/concept/unique_id.rs`) |
| `UniqueIdIndexer<FeatureMap>` base — `uniqueIdToIndex` | `FeatureMap::unique_id_to_index` |
| `UniqueIdIndexer<FeatureMap>` — `updateUniqueIdToIndex`, `resolveUniqueIdConflicts`, `swap` | not ported as public API: `CONCEPT/UniqueIdIndexer.h` is not this work package's header. The conflict resolution the map merges depend on is reproduced privately in `map_operations.rs`; the cached hash map itself has no counterpart, because `unique_id_to_index` rebuilds the index per call rather than keeping a mutable cache |
| `MapUtilities<FeatureMap>` base — `applyFunctionOnPeptideHits` / `applyFunctionOnPeptideIDs`, const and non-const | not ported: `DATASTRUCTURES/Utils/MapUtilities.h` is a different header, and plain Rust iteration over `features` and `unassigned_peptide_identifications` is the counterpart (`map.features.iter_mut().flat_map(\|f\| &mut f.peptide_identifications)`) |

### Construction, comparison and arithmetic

| Source member | Rust counterpart |
| --- | --- |
| `FeatureMap()` | `FeatureMap::new()` / `Default` |
| `FeatureMap(const FeatureMap&)` | `Clone`. The source additionally merges `id_data_` and repoints the features; the graph is the caller's here, so a clone carries the owner-tagged IDs unchanged |
| `FeatureMap(FeatureMap&&)` | Rust move |
| `FeatureMap& operator=(const FeatureMap&)` | assignment from a `Clone` |
| `FeatureMap& operator=(FeatureMap&&)` | Rust move |
| `~FeatureMap()` | `Drop` (derived) |
| `bool operator==(const FeatureMap&) const` | `PartialEq` (derived). Two differences, both observable: the source compares the cached `RangeManager`, which this port does not store, and its `DocumentIdentifier::operator==` compares only the identifier, while the derived `PartialEq` also compares `loaded_file_path` and `loaded_file_type` |
| `bool operator!=(const FeatureMap&) const` | `PartialEq` |
| `FeatureMap operator+(const FeatureMap&) const` | **new** `FeatureMap::merged(&rhs, &mut generator)` |
| `FeatureMap& operator+=(const FeatureMap&)` | **new** `FeatureMap::append(&rhs, &mut generator) -> Result<usize>`, returning the number of redrawn unique IDs |

### Sorting, ranges and swapping

| Source member | Rust counterpart |
| --- | --- |
| `void sortByIntensity(bool reverse)` | `FeatureMap::sort_by_intensity(reverse)` |
| `void sortByPosition()` | `FeatureMap::sort_by_position()` |
| `void sortByRT()` | `FeatureMap::sort_by_rt()` |
| `void sortByMZ()` | `FeatureMap::sort_by_mz()` |
| `void sortByOverallQuality(bool reverse)` | `FeatureMap::sort_by_quality(reverse)` |
| `void updateRanges()` | `FeatureMap::ranges()`, recomputed on demand |
| `void swapFeaturesOnly(FeatureMap&)` | **new** `FeatureMap::swap_features_only(&mut from)` |
| `void swap(FeatureMap&)` | **new** `FeatureMap::swap(&mut from)` |

### Records

| Source member | Rust counterpart |
| --- | --- |
| `getProteinIdentifications()` (const and non-const), `setProteinIdentifications` | the public `protein_identifications: Vec<ProteinIdentification>` field |
| `findProteinIdentification(const std::string&) const` | **new** `FeatureMap::find_protein_identification(identifier) -> Option<&_>` |
| `findProteinIdentification(const std::string&)` | **new** `FeatureMap::find_protein_identification_mut(identifier)` |
| `getUnassignedPeptideIdentifications()` (const and non-const), `setUnassignedPeptideIdentifications` | the public `unassigned_peptide_identifications` field |
| `getDataProcessing()` (const and non-const), `setDataProcessing` | the public `data_processing: Vec<DataProcessing>` field |
| `void setPrimaryMSRunPath(const StringList&)` | **new** `FeatureMap::set_primary_ms_run_path(&paths)` |
| `void setPrimaryMSRunPath(const StringList&, MSExperiment&)` | **new** `FeatureMap::set_primary_ms_run_path_from_experiment(&paths, &experiment)` |
| `void getPrimaryMSRunPath(StringList&) const` | **new** `FeatureMap::primary_ms_run_path() -> Result<Vec<String>>` |
| `void clear(bool clear_meta_data)` | `FeatureMap::clear(clear_metadata)` |

### Traversal, statistics and identifications

| Source member | Rust counterpart |
| --- | --- |
| `template<Type> Size applyMemberFunction(Size (Type::*)())` | **new** `FeatureMap::for_each_unique_id(visit)` |
| `template<Type> Size applyMemberFunction(Size (Type::*)() const) const` | **new** `FeatureMap::count_unique_ids(visit)` |
| `AnnotationStatistics getAnnotationStatistics() const` | **new** `FeatureMap::annotation_statistics(graph)` |
| `std::set<ObservationMatchRef> getUnassignedIDMatches() const` | **new** `FeatureMap::unassigned_id_matches(&graph)` |
| `const IdentificationData& getIdentificationData() const` | not ported: the graph is the caller's (see above) |
| `IdentificationData& getIdentificationData()` | not ported: same reason |
| protected `id_data_` | not ported: same reason; `update_id_references` replaces the implicit repointing |

## API mapping — `KERNEL/ConsensusMap.h`

### `enum class SplitMeta` and `struct ColumnHeader`

| Source member | Rust counterpart |
| --- | --- |
| `SplitMeta::DISCARD` / `COPY_ALL` / `COPY_FIRST` | **new** `SplitMeta::Discard` / `CopyAll` / `CopyFirst` |
| `ColumnHeader()` / copy constructor / copy assignment | `ColumnHeader::default()` / `Clone` / assignment |
| `std::string filename`, `std::string label`, `Size size`, `UInt64 unique_id` | the public `filename`, `label`, `size`, `unique_id` fields |
| `MetaInfoInterface` base of `ColumnHeader` | the public `metadata: MetaInfo` field |
| `unsigned getLabelAsUInt(const std::string&) const` | `ColumnHeader::label_as_uint(experiment_type)` |
| typedef `ColumnHeaders` (`std::map<UInt64, ColumnHeader>`) | **new** `map_operations::ColumnHeaders`, an alias for `BTreeMap<u64, ColumnHeader>`, the type of the public `ConsensusMap::column_headers` field |
| typedef `FeatureType` | not ported: `ConsensusFeature` is named directly |
| typedefs `RangeManagerContainerType`, `RangeManagerType`, `Iterator`, `ConstIterator`, `ReverseIterator`, `ConstReverseIterator` | not ported, as for `FeatureMap` |
| `EXPOSED_VECTOR_INTERFACE(ConsensusFeature)` | the public `features: Vec<ConsensusFeature>` field |
| `MetaInfoInterface`, `DocumentIdentifier`, `RangeManagerContainer`, `UniqueIdInterface`, `UniqueIdIndexer<ConsensusMap>`, `MapUtilities<ConsensusMap>` bases | as for `FeatureMap` above |

### Construction, comparison and merging

| Source member | Rust counterpart |
| --- | --- |
| `ConsensusMap()` | `ConsensusMap::new()` / `Default` |
| `ConsensusMap(const ConsensusMap&)` | `Clone` |
| `ConsensusMap(ConsensusMap&&)` | Rust move |
| `~ConsensusMap()` | `Drop` (derived) |
| `explicit ConsensusMap(size_type n)` | **new** `ConsensusMap::with_size(n) -> Result<Self>` (checked against `MAX_ITEMS`) |
| `ConsensusMap& operator=(const ConsensusMap&)` | assignment from a `Clone` |
| `ConsensusMap& operator=(ConsensusMap&&)` | Rust move |
| `bool operator==(const ConsensusMap&) const` / `operator!=` | `PartialEq` (derived), with the same two differences as `FeatureMap`'s |
| `ConsensusMap& appendRows(const ConsensusMap&)` | **new** `ConsensusMap::append_rows(&rhs, &mut generator) -> Result<usize>` |
| `ConsensusMap& appendColumns(const ConsensusMap&)` | **new** `ConsensusMap::append_columns(&rhs, &mut generator) -> Result<usize>` |
| `void clear(bool clear_meta_data)` | `ConsensusMap::clear(clear_metadata)` |

### Columns, sorting and swapping

| Source member | Rust counterpart |
| --- | --- |
| `getColumnHeaders()` (const and non-const), `setColumnHeaders` | the public `column_headers: BTreeMap<u64, ColumnHeader>` field |
| `const std::string& getExperimentType() const` | the public `experiment_type: String` field |
| `void setExperimentType(const std::string&)` | **new** `ConsensusMap::set_experiment_type(&str) -> Result<()>`, the checked setter; `ConsensusMap::validate` also rejects an invalid value written through the field |
| `void sortByIntensity(bool reverse)` | `ConsensusMap::sort_by_intensity(reverse)` |
| `void sortByRT()` | `ConsensusMap::sort_by_rt()` |
| `void sortByMZ()` | `ConsensusMap::sort_by_mz()` |
| `void sortByPosition()` | `ConsensusMap::sort_by_position()` |
| `void sortByQuality(bool reverse)` | `ConsensusMap::sort_by_quality(reverse)` |
| `void sortBySize()` | `ConsensusMap::sort_by_size()` |
| `void sortByMaps()` | `ConsensusMap::sort_by_maps()` |
| `void sortPeptideIdentificationsByMapIndex()` | **new** `ConsensusMap::sort_peptide_identifications_by_map_index()` |
| `void updateRanges()` | `ConsensusMap::ranges()`, recomputed on demand |
| `void swap(ConsensusMap&)` | **new** `ConsensusMap::swap(&mut from)` |

### Records

| Source member | Rust counterpart |
| --- | --- |
| `getProteinIdentifications()` (const and non-const), `setProteinIdentifications(const&)` and the rvalue overload | the public `protein_identifications` field; the rvalue overload is an ordinary move assignment in Rust |
| `findProteinIdentification(const std::string&) const` | **new** `ConsensusMap::find_protein_identification(identifier)` |
| `findProteinIdentification(const std::string&)` | **new** `ConsensusMap::find_protein_identification_mut(identifier)` |
| `getUnassignedPeptideIdentifications()` (const and non-const), `setUnassignedPeptideIdentifications` | the public `unassigned_peptide_identifications` field |
| `getDataProcessing()` (const and non-const), `setDataProcessing` | the public `data_processing` field |
| `void setPrimaryMSRunPath(const StringList&)` | **new** `ConsensusMap::set_primary_ms_run_path(&paths)` |
| `void setPrimaryMSRunPath(const StringList&, MSExperiment&)` | **new** `ConsensusMap::set_primary_ms_run_path_from_experiment(&paths, &experiment)` |
| `void getPrimaryMSRunPath(StringList&) const` | **new** `ConsensusMap::primary_ms_run_path() -> Vec<String>` |

### Traversal, consistency, splitting and identifications

| Source member | Rust counterpart |
| --- | --- |
| `template<Type> Size applyMemberFunction(Size (Type::*)())` | **new** `ConsensusMap::for_each_unique_id(visit)` |
| `template<Type> Size applyMemberFunction(Size (Type::*)() const) const` | **new** `ConsensusMap::count_unique_ids(visit)` |
| `bool isMapConsistent(Logger::LogStream* stream) const` | `ConsensusMap::validate_consistency() -> Result<()>`. The source returns `false` and, when given a stream, writes a report naming every offending map index; this returns the first failure as an error and does not log, because no module of this crate logs |
| `std::vector<FeatureMap> split(SplitMeta mode) const` | **new** `ConsensusMap::split(mode) -> Result<Vec<FeatureMap>>` |
| `std::set<ObservationMatchRef> getUnassignedIDMatches() const` | **new** `ConsensusMap::unassigned_id_matches(&graph)` |
| `const IdentificationData& getIdentificationData() const` / non-const | not ported: the graph is the caller's |
| protected `id_data_` | not ported; `ConsensusMap::update_id_references` replaces the implicit repointing |
| `operator<<(std::ostream&, const ConsensusMap&)` | **new** `impl Display for ConsensusMap` |

## Preserved source conventions

These are behaviours a caller can observe, kept exactly as the source has them.

* **`swap` does not swap the meta values.** `FeatureMap::swap` and
  `ConsensusMap::swap` swap the elements, ranges, document identifier, unique ID
  and index, the records and the identification data — but never the inherited
  `MetaInfoInterface` (`FeatureMap.cpp:323`, `ConsensusMap.cpp:393`). The Rust
  `swap` leaves `metadata` alone for the same reason and says so at the item.
  `std::mem::swap` is the complete exchange.
* **`appendRows` pairs column headers positionally.** After inserting `rhs`'s
  headers, the source walks the merged map and `rhs`'s map in lockstep and, for
  the length of the shorter one, renames the merged entry to
  `mergedConsensusXMLFile` and sets its size to the sum of the two entries at
  that *position* (`ConsensusMap.cpp:85`). Which sizes are added therefore
  depends on iteration order, not on which columns describe the same run.
* **`appendColumns` shifts by the header count.** `rhs`'s column index `k`
  becomes `k + n` where `n` is the left map's header *count*
  (`ConsensusMap.cpp:168`), so sparse keys can still collide; a collision keeps
  the existing header, matching `std::map::insert`.
* **Both merges deduplicate modifications everywhere.** The fixed and variable
  modification lists of *every* protein identification in the result are sorted
  and deduplicated, including the ones the left map already held
  (`ConsensusMap.cpp:100`).
* **`ConsensusMap::setPrimaryMSRunPath` writes by position.** Paths go to the
  headers keyed `0..n-1` through `std::map::operator[]`, which default-inserts
  (`ConsensusMap.cpp:508`); a map whose headers are keyed otherwise gains
  columns instead of being renamed, and the count check applies only when
  headers already exist (`ConsensusMap.cpp:518`).
* **`FeatureMap::getPrimaryMSRunPath` falls back to `UNKNOWN`** when nothing is
  annotated; `ConsensusMap::getPrimaryMSRunPath` has no such placeholder.
* **`split` routes by `map_index` and default-inserts.** An identification whose
  `map_index` names a column the consensus feature has no handle for still
  creates a feature there (`ConsensusMap.cpp:751`), and the isobaric branch
  copies every protein identification into the first feature map, which the
  source comment itself marks as wrong (`ConsensusMap.cpp:789`).
* **`split`'s `COPY_FIRST` tests the handle map index.** The source captures the
  smallest map index over the *handles*, before routing the identifications, and
  `COPY_FIRST` throws when it is not zero — so a feature that exists at index 0
  only because an identification was routed there does not satisfy the check,
  even though the meta values are then written to whichever key is smallest
  after routing. Both halves are reproduced.
* **`split` drops the handle unique ID.** The source builds each feature through
  `BaseFeature(const FeatureHandle&)`, which slices the handle to its `Peak2D`
  base (`BaseFeature.cpp:40`), so position, intensity, charge and width survive
  and the ID does not.
* **Unassigned (zero) unique IDs never conflict.** The source's
  `updateUniqueIdToIndex` postcondition compares the number of *distinct* valid
  IDs against the number of valid IDs (`UniqueIdIndexer.h:86`), so a container
  of unassigned elements is accepted and no redraw happens.
* **`applyMemberFunction` counts the container itself.** Both maps add their own
  unique ID to the accumulated total before visiting their elements, and the
  `FeatureMap` overload recurses into subordinate features while the
  `ConsensusMap` one does not.
* **The truncating conversions keep their asymmetry**, documented in
  [CONVERSION_HELPER_SUPPORT.md](CONVERSION_HELPER_SUPPORT.md).

## Native differences

* **Ranges are computed on demand**, so `updateRanges()` has no cache to
  refresh, `swapFeaturesOnly` has no range to swap, and — the one place where a
  transcribed class-test expectation had to change — a map cleared with
  `clear(false)` compares *equal* to a default map, while the source's
  `operator==` sees the range cache `clear(false)` deliberately leaves
  populated. Both `operator==`/`operator!=` section ports say so at the
  assertion.
* **`operator==` compares more than the source does.** The derived `PartialEq`
  includes `loaded_file_path` and `loaded_file_type`, which
  `DocumentIdentifier::operator==` ignores.
* **The unique-ID generator is caller-owned.** `append`, `merged`, `append_rows`
  and `append_columns` take `&mut UniqueIdGenerator` in place of the source's
  process-wide singleton, and return the number of redrawn IDs rather than
  logging it. A redraw loop is bounded at 64 attempts per element; the source
  loops until the generator happens to produce a free value.
* **Merges are atomic.** Each merge builds the result into a temporary and
  commits only on success, so a rejected merge leaves the destination exactly as
  it was. The source mutates in place and can leave a half-merged map behind.
* **`update_id_references` verifies before it writes.** Every reference in the
  whole map (subordinate features included) is checked in a read-only pass
  before any feature is changed, so a missing translation anywhere leaves the
  map unchanged; the source updates element by element.
* **Out-of-range map indices are errors.** `split` sizes its result by the
  column count and the source then indexes it with the map index, reading out of
  bounds whenever the headers are not keyed `0..n-1` (`ConsensusMap.cpp:702`);
  this port returns `Error::InvalidValue`. The same applies to an isobaric
  consensus feature with no handle, which the source dereferences unguarded
  (`ConsensusMap.cpp:741`).
* **`sort_peptide_identifications_by_map_index` requires an integer index.** The
  source compares the `DataValue`s themselves, so a non-integer `map_index`
  yields a type-dependent order; this rejects it, and reads every key before
  reordering anything.
* **No logging.** The source writes `OPENMS_LOG_INFO`/`OPENMS_LOG_WARN` when
  document identifiers are lost in a merge, when an MS run path list is empty,
  when a path is not an mzML, when a feature map has no annotated run, when
  unique IDs are replaced, and inside `isMapConsistent`. No module of this crate
  logs; every one of those messages is either carried into the rustdoc of the
  member, returned (the redraw count), or has no behavioural effect.
* **The source is serial here.** Neither `FeatureMap.cpp` nor `ConsensusMap.cpp`
  carries an OpenMP pragma, so this port loses no parallelism. It is serial
  throughout regardless, as `docs/REPOSITORY_ANALYSIS.md` requires.

## Checked boundaries and evidence

`FeatureMap::MAX_ITEMS` and `ConsensusMap::MAX_ITEMS` are 10,000,000 elements,
records or processing entries per checked operation, and
`ConsensusMap::MAX_COLUMNS` is 1,000,000 headers. Every operation whose cost
scales with the input checks its ceiling before allocating: the two merges check
the *combined* counts, `split` checks elements and columns, the traversals and
the identification queries check the element count, and
`BaseFeature::MAX_ID_MATCHES` bounds the match sets. A single element gets at
most `MAX_UNIQUE_ID_REDRAWS` (64) generator draws while conflicts are resolved.
Shifted column and map indices, header size sums and accumulated counts all use
checked arithmetic.

Evidence is **tier 3 (source review)**: every literal in
[`tests/map_operations.rs`](../tests/map_operations.rs) is transcribed from the
two pinned class tests, hashed in
[`tests/data/map_operations_provenance.json`](../tests/data/map_operations_provenance.json).
No C++ was built or executed, and no retained C++ output exists for these
headers. The checks the port adds — the ceilings, the bounded redraw, the
rejected out-of-range map index, the atomic translation — are tier 4
(Rust-only invariants) derived from the source lines anchored in the manifest.

### Class-test section coverage

`FeatureMap_test.cpp`: 32 sections, all 32 ported.

| Section | Test |
| --- | --- |
| `FeatureMap()`, `~FeatureMap()` | `fm_default_constructor_and_destructor` |
| `getProteinIdentifications()` const / non-const, `setProteinIdentifications` | `fm_protein_identification_accessors` |
| `getUnassignedPeptideIdentifications()` const / non-const, setter | `fm_unassigned_peptide_identification_accessors` |
| `getDataProcessing()` const / non-const, `setDataProcessing` | `fm_data_processing_accessors` |
| `updateRanges()` | `fm_update_ranges` |
| copy constructor, `operator=`, move `operator=` | `fm_copy_assign_and_move` |
| `operator==`, `operator!=` | `fm_equality_and_inequality` |
| `operator+` | `fm_operator_plus` |
| `operator+=` | `fm_operator_plus_assign` |
| `sortByIntensity` | `fm_sort_by_intensity` |
| `sortByPosition` | `fm_sort_by_position` |
| `sortByMZ` | `fm_sort_by_mz` |
| `sortByRT` | `fm_sort_by_rt` |
| `swap` | `fm_swap` |
| `swapFeaturesOnly` | `fm_swap_features_only` |
| `sortByOverallQuality` | `fm_sort_by_overall_quality` |
| `clear` | `fm_clear` |
| `[EXTRA] uniqueIdToIndex` | `fm_unique_id_to_index` |
| `applyMemberFunction` (both overloads) | `fm_apply_member_function` |
| `getAnnotationStatistics` | `fm_annotation_statistics` |
| `[EXTRA] ExposedVector Ctor` | `fm_exposed_vector_constructors` |

`ConsensusMap_test.cpp`: 39 sections, all 39 ported. The eight `NOT_TESTABLE`
sorting sections are ported as real tests rather than mapped away.

| Section | Test |
| --- | --- |
| `ConsensusMap()`, `~ConsensusMap()` | `cm_default_constructor_and_destructor` |
| `getProteinIdentifications()` const / non-const, `setProteinIdentifications` | `cm_protein_identification_accessors` |
| `getUnassignedPeptideIdentifications()` const / non-const, setter | `cm_unassigned_peptide_identification_accessors` |
| `getDataProcessing()` const / non-const, `setDataProcessing` | `cm_data_processing_accessors` |
| `updateRanges()` | `cm_update_ranges` |
| `appendRows` | `cm_append_rows` |
| `appendColumns` | `cm_append_columns` |
| `operator=` const, copy constructor, move `operator=`, move constructor | `cm_copy_assign_and_move` |
| `ConsensusMap(size_type)` | `cm_size_constructor` |
| `[ColumnHeader] ColumnHeader()` | `cm_column_header_default_constructor` |
| `getColumnHeaders()` const / non-const | `cm_column_header_accessors` |
| `getExperimentType()`, `setExperimentType` | `cm_experiment_type` |
| `swap` | `cm_swap` |
| `operator==`, `operator!=` | `cm_equality_and_inequality` |
| `sortByIntensity`, `sortByRT`, `sortByMZ`, `sortByPosition`, `sortByQuality`, `sortBySize`, `sortByMaps` | `cm_sorting` |
| `sortPeptideIdentificationsByMapIndex` | `cm_sort_peptide_identifications_by_map_index` |
| `clear` | `cm_clear` |
| `applyMemberFunction` (both overloads) | `cm_apply_member_function` |
| `split` | `cm_split` |

Native tests without a source section: `fm_append_resolves_unique_id_conflicts_atomically`,
`fm_annotation_statistics_display`, `fm_display`, `fm_primary_ms_run_path`,
`find_protein_identification_returns_the_first_match`,
`primary_ms_run_path_from_experiment_falls_back`,
`cm_append_columns_shifts_map_index_annotations`,
`cm_split_copy_first_uses_the_handle_map_index`,
`cm_split_carries_processing_and_handle_fields`, `cm_split_error_paths`,
`cm_primary_ms_run_path`, `cm_display`,
`unassigned_id_matches_and_reference_translation`,
`ceilings_are_declared_and_checked`, `document_identity_moves_and_resets`.

# Feature identification support

[`src/kernel/feature_identification.rs`](../src/kernel/feature_identification.rs)
completes the identification surface of `KERNEL/BaseFeature.h`,
`KERNEL/Feature.h` and `KERNEL/ConsensusFeature.h` at Core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The containers themselves are in
[`src/kernel/features.rs`](../src/kernel/features.rs) and are described by
[FEATURE_SUPPORT.md](FEATURE_SUPPORT.md); that document's sentence "Identification
references, annotation-state inference, ratios ... are not implemented here" is
superseded by this work package for the members listed below.

The module holds operations only, no types. It exists beside `features.rs` for
the same reason as [`kernel/gap_closures.rs`](../src/kernel/gap_closures.rs):
`features.rs` was frozen for this work package, and Rust allows inherent and
trait implementations for a crate-local type from any module of the same crate.
Callers see one type with one API.

## Ownership: identifications attach by reference

The source `FeatureMap` owns an `IdentificationData id_data_` member
(`FeatureMap.h:294`) and its features hold iterators into it. This port
deliberately does not embed the graph:

* `IdentificationData` is not `Clone` by design. IDs are owner-tagged, so a copy
  needs `try_clone_with_translation`, not a derived clone.
* Embedding it would therefore strip `Clone` and `PartialEq` from `FeatureMap`
  and `ConsensusMap`, which 106 map-equality assertions in the existing test
  suite depend on.

A feature therefore stores owner-tagged IDs
(`BaseFeature::primary_id`, `BaseFeature::id_matches`) and every operation that
must *resolve* one takes the graph, or the `ReferenceTranslator`, as a
parameter: `annotation_state(Some(&graph))`, `set_primary_id_checked(&graph, …)`,
`add_id_match_checked(&graph, …)`, `update_id_references(&translator)`,
`update_all_id_references(&translator)`. The consequence a caller sees is that a
stale or foreign reference is an error at the call site instead of a dangling
iterator: the source cannot check this at all, because its reference is a
container iterator.

## API mapping

Every public member of the three headers appears below. Members already covered
by earlier work packages name their existing counterpart; members this work
package adds are marked **new**.

### `KERNEL/BaseFeature.h`

| Source member | Rust counterpart |
| --- | --- |
| `typedef float QualityType` | `f32` (`BaseFeature::quality`) |
| `typedef Int ChargeType` | `i32` (`BaseFeature::charge`) |
| `typedef float WidthType` | `f32` (`BaseFeature::width`) |
| `enum class AnnotationState` | **new** `AnnotationState` (`None`, `Single`, `MultipleSame`, `MultipleDivergent`) |
| `AnnotationState::SIZE_OF_ANNOTATIONSTATE` | not ported: a C++ array-sizing sentinel; `AnnotationState::NAMES` is a fixed-length array |
| `static const std::string NamesOfAnnotationState[]` | **new** `AnnotationState::NAMES`, `AnnotationState::name()`, `Display` |
| `BaseFeature()` | `Default` / `BaseFeature::new(rt, mz, intensity)` |
| `BaseFeature(const BaseFeature&)` | `Clone` |
| `BaseFeature(BaseFeature&&) noexcept` | Rust move; cannot fail, so `noexcept` has no counterpart |
| `BaseFeature(const BaseFeature& rhs, UInt64 map_index)` | **new** `BaseFeature::clone_with_map_index(map_index)` |
| `explicit BaseFeature(const Peak2D&)` | not ported as a conversion: write `BaseFeature { rt: point.rt(), mz: point.mz(), intensity: point.intensity, ..BaseFeature::default() }`, which is what the source constructor produces (quality, charge, width and identifications all zero/empty) |
| `explicit BaseFeature(const RichPeak2D&)` | not ported as a conversion: as above, additionally carrying `metadata` and `unique_id` from `point` |
| `explicit BaseFeature(const FeatureHandle&)` | not ported as a conversion: the source slices the handle to its `Peak2D` base, so write the struct with `rt`, `mz`, `intensity`, `charge` and `width` from the handle and **not** its `unique_id`. The inverse direction, `FeatureHandle::new(map_index, &BaseFeature)`, exists |
| `~BaseFeature()` | `Drop` (derived) |
| `getQuality()` / `setQuality(QualityType)` | `quality` field |
| `struct QualityLess` (4 `operator()` overloads) | `quality` field ordering; `FeatureMap::sort_by_quality(reverse)` is the container-level user |
| `getWidth()` | `width` field |
| `setWidth(WidthType)` | `BaseFeature::set_width` (also writes the `FWHM` metadata alias; rejects negative and non-finite, which the source stores) |
| `getCharge()` / `setCharge(const ChargeType&)` | `charge` field |
| `operator=` / move `operator=` | assignment |
| `operator==` / `operator!=` | derived `PartialEq` over all fields, including `primary_id` and `id_matches`, as the source's `operator==` |
| `getPeptideIdentifications() const` / non-const / `setPeptideIdentifications` | `peptide_identifications` field |
| `sortPeptideIdentifications()` | **new** `BaseFeature::sort_peptide_identifications()` |
| `getAnnotationState() const` | **new** `BaseFeature::annotation_state(graph)` |
| `hasPrimaryID() const` | **new** `BaseFeature::has_primary_id()` |
| `getPrimaryID() const` | **new** `BaseFeature::primary_id()` → `Result<IdentifiedMolecule>` |
| `setPrimaryID(const IdentifiedMolecule&)` | **new** `BaseFeature::set_primary_id(id)`, plus `set_primary_id_checked(&graph, id)` with no source counterpart |
| `clearPrimaryID()` | **new** `BaseFeature::clear_primary_id()` |
| `getIDMatches() const` / non-const | `id_matches` field (`BTreeSet<ObservationMatchId>`) |
| `addIDMatch(ObservationMatchRef)` | **new** `BaseFeature::add_id_match(id)`, plus `add_id_match_checked(&graph, id)` |
| `updateIDReferences(const RefTranslator&)` | **new** `BaseFeature::update_id_references(&translator)` |
| protected `quality_`, `charge_`, `width_`, `peptides_`, `primary_id_`, `id_matches_` | public fields of the same names (`peptide_identifications` for `peptides_`) |
| inherited `RichPeak2D` surface | `rt`, `mz`, `intensity`, `metadata`, `unique_id` fields; see [FEATURE_SUPPORT.md](FEATURE_SUPPORT.md) |

### `KERNEL/Feature.h`

| Source member | Rust counterpart |
| --- | --- |
| `Feature()` | `Default` / `Feature::new(rt, mz, intensity)` |
| `explicit Feature(const BaseFeature&)` | `From<BaseFeature>` |
| `Feature(const Feature&)`, `Feature(Feature&&)`, `~Feature()` | `Clone`, move, `Drop` |
| `getOverallQuality()` / `setOverallQuality(QualityType)` | `quality` field (through `Deref` to `BaseFeature`) |
| `getQuality(Size index)` / `setQuality(Size, QualityType)` | `quality_rt` and `quality_mz` fields; the source's `OPENMS_PRECONDITION(index < 2)` has no counterpart because an out-of-range dimension cannot be named |
| `typedef QualityLess OverallQualityLess` | as `QualityLess` above |
| `getConvexHulls()` const / non-const / `setConvexHulls` | `convex_hulls` field |
| `getConvexHull() const` | `Feature::convex_hull()` (computed on demand, not cached) |
| `encloses(double rt, double mz)` | `Feature::encloses(rt, mz)` → `Result<bool>` |
| `operator=`, move `operator=`, `operator==` | assignment and derived `PartialEq` |
| `getSubordinates()` const / non-const / `setSubordinates` | `subordinates` field |
| `applyMemberFunction(Size (Type::*)())` | **new** `Feature::for_each_unique_id(visit)` |
| `applyMemberFunction(Size (Type::*)() const) const` | **new** `Feature::count_unique_ids(visit)` |
| `updateAllIDReferences(const RefTranslator&)` | **new** `Feature::update_all_id_references(&translator)` |
| protected `qualities_[2]` | `quality_rt`, `quality_mz` |
| protected `convex_hulls_`, `subordinates_` | public fields |
| protected `convex_hulls_modified_`, `convex_hull_` | not ported: the overall hull is recomputed on demand, so there is no invalidation flag and no mutable cache to keep consistent |

### `KERNEL/ConsensusFeature.h`

| Source member | Rust counterpart |
| --- | --- |
| `typedef std::set<FeatureHandle, IndexLess> HandleSetType` | `Vec<FeatureHandle>` kept sorted and unique by `(map_index, unique_id)`; read through `ConsensusFeature::handles()` |
| `const_iterator`, `iterator`, `const_reverse_iterator`, `reverse_iterator` | `handles().iter()` and `.iter().rev()`; there is no mutable iterator, because editing a handle's identity in place would break the set invariant. Copy out, edit, `set_handles` |
| `struct SizeLess` (4 `operator()` overloads) | `ConsensusFeature::len()` ordering; `ConsensusMap::sort_by_size()` is the container-level user (descending) |
| `struct MapsLess` | lexicographic ordering of `handles().iter().map(FeatureHandle::key)`; `ConsensusMap::sort_by_maps()` |
| `struct Ratio` | `Ratio { ratio_value, denominator_ref, numerator_ref, description }` |
| `Ratio::Ratio()`, copy ctor, `operator=`, `~Ratio()` | `Default` (a defined `0.0`, where the source leaves `ratio_value_` uninitialised), `Clone`, assignment, `Drop`. The source's `virtual` destructor on a value type stored in a `std::vector` has no counterpart |
| `Ratio::ratio_value_`, `denominator_ref_`, `numerator_ref_` | `ratio_value`, `denominator_ref`, `numerator_ref` |
| `Ratio::description_` (`std::vector<std::string>`) | the public `Ratio::description` field (`src/kernel/features.rs`), with `Ratio::description()`, `add_description`, `set_description` and `validate_description` as the checked path over it (`src/kernel/consensus_display.rs`). This package left the field out and recorded the gap as a deferral; the residual-closure package added it. See [CONSENSUS_DISPLAY_SUPPORT.md](CONSENSUS_DISPLAY_SUPPORT.md) |
| `ConsensusFeature()`, copy, move, `~ConsensusFeature()` | `new` / `Default`, `Clone`, move, `Drop` |
| `explicit ConsensusFeature(const BaseFeature&)` | `From<BaseFeature>` |
| `ConsensusFeature(UInt64 map_index, const Peak2D&, UInt64 element_index)` | `From<BaseFeature>` plus `insert(FeatureHandle::from_peak(map_index, point, element_index))` |
| `ConsensusFeature(UInt64 map_index, const BaseFeature&)` | `ConsensusFeature::from_feature(map_index, &feature)`; the source stamps `map_index` on the peptide identifications through `BaseFeature(rhs, map_index)` and `from_feature` does not, so pass `feature.clone_with_map_index(map_index)` when that is wanted |
| `insert(const ConsensusFeature&)` / `insert(ConsensusFeature&&)` | `ConsensusFeature::merge(&other)` |
| `insert(const FeatureHandle&)` / `insert(FeatureHandle&&)` | `ConsensusFeature::insert(handle)` |
| `insert(const HandleSetType&)` / `insert(HandleSetType&&)` | `ConsensusFeature::set_handles(handles)` |
| `insert(UInt64, const Peak2D&, UInt64)` | `insert(FeatureHandle::from_peak(…))` |
| `insert(UInt64, const BaseFeature&)` | `insert(FeatureHandle::new(map_index, &feature))` |
| `getFeatures()` | `ConsensusFeature::handles()` |
| `getFeatureList()` | `handles().to_vec()` |
| `setFeatures(HandleSetType)` | `ConsensusFeature::set_handles` |
| `getPositionRange()` / `getIntensityRange()` | `ConsensusFeature::handle_ranges()` (`rt`, `mz`, `intensity` as `Option<NumericRange>`) |
| `computeConsensus()` | `ConsensusFeature::compute_consensus()` |
| `computeMonoisotopicConsensus()` | `ConsensusFeature::compute_monoisotopic_consensus()` |
| `computeDechargeConsensus(const FeatureMap&, bool)` | `ConsensusFeature::compute_decharge_consensus(&map, weighted)` |
| `addRatio(const Ratio&)` | **new** `ConsensusFeature::add_ratio(ratio)` |
| `setRatios(std::vector<Ratio>&)` | **new** `ConsensusFeature::set_ratios(ratios)`, taking the vector by value |
| `getRatios() const` | **new** `ConsensusFeature::ratios()` |
| `getRatios()` (non-const) | the public `ratios` field, which bypasses the checks exactly as the source accessor does |
| `size()`, `empty()`, `clear()` | `len()`, `is_empty()`, `clear()` |
| `begin()`, `end()`, `rbegin()`, `rend()` (const and non-const) | `handles()` slice as above |
| private `handles_`, `ratios_` | private `handles`, public `ratios` |
| `operator<<(std::ostream&, const ConsensusFeature&)` | `impl Display for ConsensusFeature` (`src/kernel/consensus_display.rs`). This package had no `Display` and recorded the gap as a deferral; the residual-closure package added it. `FeatureHandle`'s stream operator is ported in `gap_closures.rs`. See [CONSENSUS_DISPLAY_SUPPORT.md](CONSENSUS_DISPLAY_SUPPORT.md) |

## Preserved source conventions

* **Annotation state dispatch.** A non-empty match set wins outright: the legacy
  peptide identifications are then ignored, exactly as `getAnnotationState`
  branches. That branch never reports `None`.
* **Annotation state from one sequence.** Two identifications of which only one
  has hits collect a single sequence and therefore report `MultipleSame`, not
  `Single`. The source's `size() == 1` shortcut only applies to a single
  identification, so this is its behaviour, and it is preserved.
* **Best hit only.** Where an identification carries several hits, only its best
  hit is considered, and the comparison is by sequence string.
* **Sort order.** `sort_peptide_identifications` puts the best identification
  first and identifications without hits last, which is what the source's
  reverse-iterator `std::sort` produces.
* **`map_index` metadata key.** `clone_with_map_index` writes the source key
  `map_index` on every attached identification.
* **Ratios are experimental.** The source `@note` that the consensus feature
  handler ignores ratios still holds; nothing in this port reads them either.
* **Charge, position and decharge arithmetic** are unchanged from
  [FEATURE_SUPPORT.md](FEATURE_SUPPORT.md); this module adds no numerics.

## Native differences

* **The graph is a parameter, not a member.** See the ownership section above.
  `annotation_state` takes `Option<&IdentificationData>`; passing `None` with a
  non-empty match set is `Error::MissingInformation` rather than a silent
  answer.
* **Mixed score directions are refused.** The source comparator reads
  `isHigherScoreBetter()` from its left operand only, so identifications that
  disagree give `std::sort` an asymmetric comparator and therefore undefined
  behaviour. `sort_peptide_identifications` returns `Error::InvalidValue`
  instead.
* **Every identification's hits are sorted.** The source sorts hits inside the
  comparator, so an identification that never takes part in a comparison keeps
  unsorted hits — with a single attached identification, none are sorted at all.
* **Empty identifications are ordered, not compared.** The source comparator
  answers "less" for two empty identifications, which is not a strict weak
  ordering; here they simply sort last and keep their relative order.
* **Stable order.** Identifications with equal best scores keep their relative
  order. The source's `std::sort` may permute them, and its unstable hit sort
  means a tie between two differently-sequenced best hits can change the
  annotation state from run to run; the port always takes the first best hit in
  stored order.
* **Atomic reference updates.** `BaseFeature::update_id_references` translates
  into temporaries and commits only when every translation succeeded; the source
  swaps `id_matches_` out first and leaves a partially translated set behind
  when `translate` throws. `Feature::update_all_id_references` extends that to
  the whole subordinate tree with a read-only verification pass before any
  feature changes, where the source updates as it descends.
* **No `allow_missing` translation.** The source `RefTranslator` can be told to
  keep an untranslatable reference; `ReferenceTranslator` has no such mode,
  because a reference the new graph does not own can never be resolved.
* **Iterative traversal.** `for_each_unique_id`, `count_unique_ids` and
  `update_all_id_references` walk the subordinate tree with an explicit stack in
  the source's pre-order, so a deep tree does not recurse on the machine stack.
* **Checked ratios.** `add_ratio` and `set_ratios` reject a non-finite ratio
  value; the source appends and assigns unchecked, and its `Ratio()` leaves
  `ratio_value_` indeterminate.
* **`set_ratios` takes ownership.** The source signature
  `setRatios(std::vector<Ratio>&)` cannot be called with a temporary or a
  `const` vector at all.
* **Serial.** 36 source files carry `#pragma omp`; none of the members ported
  here is one of them, so there is no parallel-execution gap for this module
  specifically. The port is serial by policy regardless.

## Checked boundaries and evidence

Ceilings, all checked before anything is allocated or mutated, following
`src/identification/run_mapping.rs`:

| Constant | Value | Guards |
| --- | --- | --- |
| `BaseFeature::MAX_PEPTIDE_IDENTIFICATIONS` | 1 000 000 | `sort_peptide_identifications`, `annotation_state`, `clone_with_map_index` |
| `BaseFeature::MAX_PEPTIDE_HITS` | 4 000 000 | `sort_peptide_identifications`, `annotation_state` (cumulative across identifications) |
| `BaseFeature::MAX_ID_MATCHES` | 1 000 000 | `add_id_match`, `update_id_references`, `annotation_state` |
| `ConsensusFeature::MAX_RATIOS` | 1 000 000 | `add_ratio`, `set_ratios` |
| `Feature::MAX_TRAVERSED_FEATURES` | 1 000 000 | `for_each_unique_id`, `count_unique_ids`, `update_all_id_references` |
| `Feature::MAX_SUBORDINATE_DEPTH` | 128 (existing) | the same three traversals, matching `Feature::validate` |

Atomicity: `sort_peptide_identifications` performs every fallible check first
and then only infallible work, so a failure leaves the feature untouched;
`clone_with_map_index` never mutates the receiver; `update_id_references`,
`update_all_id_references`, `add_ratio` and `set_ratios` build into temporaries
and commit. `for_each_unique_id` is the exception and says so at the item: the
caller's closure may already have run on the features visited before a limit was
reached, because the limits are checked during the walk.

Accumulator overflow in `for_each_unique_id` / `count_unique_ids` is an error,
not a wrap: the closure's return value is caller-supplied.

**Evidence: tier 3 (source review).** Every expectation in
[tests/feature_identification.rs](../tests/feature_identification.rs) is
transcribed from the pinned class tests or derived from the pinned
implementation. No C++ was built or executed; no retained C++ output exists for
these headers, so tier 1 is not available for them. Hashes and source anchors
are in
[tests/data/feature_identification_provenance.json](../tests/data/feature_identification_provenance.json).

### Class-test section accounting

All 92 sections of the three class tests are ported, one Rust test each.

| Source section | Rust test |
| --- | --- |
| `BaseFeature()` | `bf_default_constructor` |
| `~BaseFeature()` | `bf_destructor` |
| `QualityType getQuality() const` | `bf_get_quality` |
| `void setQuality(QualityType)` | `bf_set_quality` |
| `WidthType getWidth() const` | `bf_get_width` |
| `void setWidth(WidthType)` | `bf_set_width` |
| `[EXTRA] IntensityType getIntensity() const` | `bf_get_intensity_const` |
| `[EXTRA] const PositionType& getPosition() const` | `bf_get_position_const` |
| `[EXTRA] IntensityType& getIntensity()` | `bf_mutable_intensity` |
| `[EXTRA] PositionType& getPosition()` | `bf_mutable_position` |
| `const ChargeType& getCharge() const` | `bf_get_charge` |
| `void setCharge(const ChargeType&)` | `bf_set_charge` |
| `BaseFeature(const BaseFeature&)` | `bf_copy_constructor` |
| `BaseFeature(BaseFeature&&)` | `bf_move_constructor` |
| `BaseFeature(const Peak2D&)` | `bf_from_peak2d` |
| `BaseFeature(const RichPeak2D&)` | `bf_from_rich_peak2d` |
| `BaseFeature& operator=(const BaseFeature&)` | `bf_assignment_operator` |
| `bool operator==(const BaseFeature&) const` | `bf_equality_operator` |
| `bool operator!=(const BaseFeature&) const` | `bf_inequality_operator` |
| `[EXTRA] meta info with copy constructor` | `bf_meta_info_with_copy_constructor` |
| `[EXTRA] meta info with assignment` | `bf_meta_info_with_assignment` |
| `[QualityLess] (BaseFeature, BaseFeature)` | `bf_quality_less_feature_feature` |
| `[QualityLess] (BaseFeature, QualityType)` | `bf_quality_less_feature_value` |
| `[QualityLess] (QualityType, BaseFeature)` | `bf_quality_less_value_feature` |
| `[QualityLess] (QualityType, QualityType)` | `bf_quality_less_value_value` |
| `const PeptideIdentificationList& getPeptideIdentifications() const` | `bf_get_peptide_identifications_const` |
| `void setPeptideIdentifications(…)` | `bf_set_peptide_identifications` |
| `PeptideIdentificationList& getPeptideIdentifications()` | `bf_get_peptide_identifications_mut` |
| `AnnotationState getAnnotationState() const` | `bf_get_annotation_state` |
| `sortPeptideIdentifications()` | `bf_sort_peptide_identifications` |
| `Feature()` | `f_default_constructor` |
| `~Feature()` | `f_destructor` |
| `QualityType getOverallQuality() const` | `f_get_overall_quality` |
| `void setOverallQuality(QualityType)` | `f_set_overall_quality` |
| `QualityType getQuality(Size) const` | `f_get_quality_by_index` |
| `void setQuality(Size, QualityType)` | `f_set_quality_by_index` |
| `const vector<ConvexHull2D>& getConvexHulls() const` | `f_get_convex_hulls_const` |
| `vector<ConvexHull2D>& getConvexHulls()` | `f_get_convex_hulls_mut` |
| `void setConvexHulls(…)` | `f_set_convex_hulls` |
| `ConvexHull2D& getConvexHull() const` | `f_get_convex_hull` |
| `bool encloses(double, double) const` | `f_encloses` |
| `Feature(const Feature&)` | `f_copy_constructor` |
| `Feature(const Feature&&)` | `f_move_constructor` |
| `Feature& operator=(const Feature&)` | `f_assignment_operator` |
| `bool operator==(const Feature&) const` | `f_equality_operator` |
| `[EXTRA] operator!=` | `f_inequality_operator` |
| `[EXTRA] meta info with copy constructor` | `f_meta_info_with_copy_constructor` |
| `[EXTRA] meta info with assignment` | `f_meta_info_with_assignment` |
| `std::vector<Feature>& getSubordinates()` | `f_get_subordinates_mut` |
| `void setSubordinates(…)` | `f_set_subordinates` |
| `const std::vector<Feature>& getSubordinates() const` | `f_get_subordinates_const` |
| `applyMemberFunction(Size (Type::*)())` | `f_apply_member_function_mutable` |
| `applyMemberFunction(Size (Type::*)() const) const` | `f_apply_member_function_const` |
| `ConsensusFeature()` | `cf_default_constructor` |
| `virtual ~ConsensusFeature()` | `cf_destructor` |
| `[SizeLess] (ConsensusFeature, ConsensusFeature)` | `cf_size_less_feature_feature` |
| `[SizeLess] (ConsensusFeature, UInt64)` | `cf_size_less_feature_value` |
| `[SizeLess] (UInt64, ConsensusFeature)` | `cf_size_less_value_feature` |
| `[SizeLess] (UInt64, UInt64)` | `cf_size_less_value_value` |
| `[MapsLess] (ConsensusFeature, ConsensusFeature)` | `cf_maps_less` |
| `ConsensusFeature& operator=(const ConsensusFeature&)` | `cf_assignment_operator` |
| `ConsensusFeature(const ConsensusFeature&)` | `cf_copy_constructor` |
| `ConsensusFeature(ConsensusFeature&&)` | `cf_move_constructor` |
| `void insert(const HandleSetType&)` | `cf_insert_handle_set` |
| `void insert(UInt64, const Peak2D&, UInt64)` | `cf_insert_peak_with_element_index` |
| `void insert(UInt64, const BaseFeature&)` (line 234) | `cf_insert_map_index_and_base_feature` |
| `ConsensusFeature(const BaseFeature&)` | `cf_from_base_feature` |
| `ConsensusFeature(UInt64, const BaseFeature&)` | `cf_from_map_index_and_base_feature` |
| `[EXTRA] ConsensusFeature(UInt64, const Feature&)` | `cf_from_map_index_and_feature` |
| `ConsensusFeature(UInt64, const Peak2D&, UInt64)` | `cf_from_map_index_peak_and_element_index` |
| `[EXTRA] ConsensusFeature(UInt64, const ConsensusFeature&)` | `cf_from_map_index_and_consensus_feature` |
| `DRange<1> getIntensityRange() const` | `cf_get_intensity_range` |
| `DRange<2> getPositionRange() const` | `cf_get_position_range` |
| `const HandleSetType& getFeatures() const` | `cf_get_features` |
| `std::vector<FeatureHandle> getFeatureList() const` | `cf_get_feature_list` |
| `void insert(const ConsensusFeature&)` | `cf_insert_consensus_feature` |
| `void insert(const FeatureHandle&)` | `cf_insert_feature_handle` |
| `void insert(UInt64, const BaseFeature&)` (line 468) | `cf_insert_map_index_and_base_feature_again` |
| `void computeConsensus()` | `cf_compute_consensus` |
| `void computeMonoisotopicConsensus()` | `cf_compute_monoisotopic_consensus` |
| `void computeDechargeConsensus(const FeatureMap&, bool)` | `cf_compute_decharge_consensus` |
| `Size size() const` | `cf_size` |
| `const_iterator begin() const` | `cf_begin_const` |
| `iterator begin()` | `cf_begin_mut` |
| `const_iterator end() const` | `cf_end_const` |
| `iterator end()` | `cf_end_mut` |
| `const_reverse_iterator rbegin() const` | `cf_rbegin_const` |
| `reverse_iterator rbegin()` | `cf_rbegin_mut` |
| `const_reverse_iterator rend() const` | `cf_rend_const` |
| `reverse_iterator rend()` | `cf_rend_mut` |
| `void clear()` | `cf_clear` |
| `bool empty() const` | `cf_empty` |

Ten further tests have no source section, because the class tests never touch
the primary ID, the observation matches, the reference translation or the
ratios: `annotation_state_from_matches_ignores_legacy_identifications`,
`annotation_state_names_follow_the_source_array`,
`annotation_state_counts_only_identifications_that_have_hits`,
`sort_peptide_identifications_is_checked_and_atomic`,
`primary_id_assignment_clearing_and_checking`,
`update_id_references_translates_atomically`,
`update_all_id_references_covers_subordinates_or_changes_nothing`,
`unique_id_traversal_is_pre_order_and_bounded`,
`ratios_are_checked_and_replaced_atomically`,
`clone_with_map_index_stamps_every_attached_identification`.

Self-audit (`BaseFeature.h`, `Feature.h`, `ConsensusFeature.h`): 92 ported,
0 mapped-with-evidence, 0 mapped-without-evidence, 0 unaccounted.

### Sections whose C++ construct has no Rust form

Three sections assert something about C++ itself rather than about behaviour.
They are ported as the closest observable property, and the difference is stated
in the test:

* `BaseFeature(BaseFeature&&)` and `Feature(const Feature&&)` assert the move
  constructor is `noexcept` so `std::vector` moves instead of copying. A Rust
  move is a memcpy that cannot fail, and the moved-from value cannot be read at
  all, which is stronger than the source's "the hulls are gone" check.
* `Feature::getQuality(Size)` asserts `TEST_PRECONDITION_VIOLATED(p.getQuality(10))`.
  With named `quality_rt` / `quality_mz` fields there is no index to be out of
  range.
* The metadata sections address values by registry index (`setMetaValue(2, …)`).
  `MetaInfo` is keyed by name, so the port uses the name `"2"` and asserts the
  same copy-independence.

## Validation

```sh
cargo fmt --all -- --check
cargo clippy --locked --all-features --lib --test feature_identification -- -D warnings
cargo nextest run --locked --all-features --test feature_identification
cargo nextest run --locked --no-default-features --test feature_identification
cargo test --locked --all-features --doc
RUSTDOCFLAGS="-D warnings" cargo doc --locked --all-features --no-deps
cargo +1.85.0 check --locked --all-features --lib --test feature_identification
python3 tools/check_doc_coverage.py --report | grep feature_identification
```

The module needs no cargo features, so it is wired into the `minimum-rust` job
on the kernel no-default-features line (`.github/workflows/rust.yml:79`, the one
already carrying `--test ranges --test range_utils --test spectrum_helper`).

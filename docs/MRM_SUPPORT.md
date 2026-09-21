# MRM peak groups and transition groups

The native `kernel::mrm` module implements `KERNEL/MRMFeature.h` (with
`KERNEL/MRMFeature.cpp`) and the header-only template
`KERNEL/MRMTransitionGroup.h` in Core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. See the
[implementation](../src/kernel/mrm.rs), [tests](../tests/mrm.rs) and
[hashed provenance](../tests/data/mrm_provenance.json). One Rust module covers
both headers, because `MRMTransitionGroup` stores `MRMFeature` by value and the
two invariants worth checking span them.

The increment adds no dependency and no feature gate, and does not build or
execute C++. Neither source file carries `#pragma omp`; source and port are both
serial, so there is no parallelism gap to record here.

## The generic parameters, and why they are generic

`MRMTransitionGroup` is a C++ template over `ChromatogramType` and
`TransitionType`. Both stay generic in Rust, for different reasons.

- **`T: Transition`.** The source instantiates the template with
  `ReactionMonitoringTransition` (`ANALYSIS/MRM`), which is **not ported and is
  out of scope for this work package**. Rather than wait for it, the module
  declares a `Transition` trait carrying the seven values the group and its
  immediate consumers read — `native_id`, `precursor_mz`, `product_mz`,
  `library_intensity`, `is_detecting`, `is_identifying`, `is_quantifying` —
  verified against `ReactionMonitoringTransition.h` lines 79–238. The group
  itself only ever calls `getNativeID` (in `subset` and `subsetDependent`) and
  `getLibraryIntensity` (in `getLibraryIntensity`); the other five are carried
  because every OpenSWATH consumer reads them straight off the transitions a
  group holds, and a trait without them would push those callers back onto a
  concrete type. `SimpleTransition` is a small concrete implementation for tests
  and for callers with no assay-library type; it is **not** a port of
  `ReactionMonitoringTransition` and deliberately omits its peptide/compound
  references, decoy type, CV terms and interpretation list.
- **`C: NativeIdentified`.** The source header states that "since not all the
  functions in OpenMS will work with `MSChromatogram` data structures, this
  needs to accept also `MSSpectrum` as a type for raw data storage", so the
  chromatogram parameter must stay open. The only thing the group asks of it is
  `getNativeID`, which is the `NativeIdentified` trait; it is implemented here
  for `MSChromatogram` and `MSSpectrum`.

## Public mapping

Every public member of both headers is listed.

### `KERNEL/MRMFeature.h`

| Source member | Native counterpart | Status |
| --- | --- | --- |
| `typedef std::vector<Feature> FeatureListType` | `Vec<Feature>` used directly; no alias | not ported: a one-word alias carries nothing in Rust |
| `MRMFeature()` | `MRMFeature::new()` / `MRMFeature::default()` | ported |
| `MRMFeature(const MRMFeature& rhs)` | `Clone` | ported |
| `MRMFeature(MRMFeature&& rhs) = default` | Rust move; unconditional and non-unwinding | ported |
| `MRMFeature& operator=(const MRMFeature& rhs)` | `Clone` plus assignment | ported |
| `MRMFeature& operator=(MRMFeature&&) & = default` | Rust move assignment | ported |
| `~MRMFeature() override` | `Drop` glue; no explicit destructor needed | ported |
| *(base class)* `public Feature` | `pub feature: Feature` plus `Deref`/`DerefMut` and `From<Feature>` | ported |
| `const OpenSwath_Scores& getScores() const` | `scores(&self) -> &OpenSwathScores` | ported |
| `OpenSwath_Scores& getScores()` | `scores_mut(&mut self) -> &mut OpenSwathScores` | ported |
| `void setScores(const OpenSwath_Scores&)` | `set_scores(&mut self, OpenSwathScores)` | ported |
| `void addScore(const std::string&, double)` | `add_score(&mut self, impl Into<String>, f64) -> Result<()>` | ported |
| `Feature& getFeature(const std::string& key)` | `feature_mut(&mut self, &str) -> Result<&mut Feature>` | ported |
| `const Feature& getFeature(const std::string& key) const` | `feature(&self, &str) -> Result<&Feature>` | ported |
| `void addFeature(const Feature&, const std::string&)` | `add_feature(&mut self, Feature, impl AsRef<str>) -> Result<()>` | ported |
| `void addFeature(Feature&&, const std::string&)` | same entry point; Rust takes `Feature` by value, so there is one overload | ported |
| `const std::vector<Feature>& getFeatures() const` | `features(&self) -> &[Feature]` | ported |
| `void getFeatureIDs(std::vector<std::string>&) const` | `feature_ids(&self) -> impl ExactSizeIterator<Item = &str>` | ported |
| `void addPrecursorFeature(const Feature&, const std::string&)` | `add_precursor_feature(&mut self, Feature, impl AsRef<str>) -> Result<()>` | ported |
| `void addPrecursorFeature(Feature&&, const std::string&)` | same entry point | ported |
| `void getPrecursorFeatureIDs(std::vector<std::string>&) const` | `precursor_feature_ids(&self) -> impl ExactSizeIterator<Item = &str>` | ported |
| `Feature& getPrecursorFeature(const std::string&)` | `precursor_feature_mut(&mut self, &str) -> Result<&mut Feature>` | ported |
| `const Feature& getPrecursorFeature(const std::string&) const` | `precursor_feature(&self, &str) -> Result<&Feature>` | ported |
| `void IDScoresAsMetaValue(bool, const OpenSwath_Ind_Scores&)` | `id_scores_as_meta_value(&mut self, bool, &OpenSwathIndScores) -> Result<()>` | ported |
| *(protected)* `features_`, `precursor_features_`, `pg_scores_`, `feature_map_`, `precursor_feature_map_` | private fields of the same shape; `BTreeMap<String, usize>` for the two maps | ported as private state |
| *(removed member)* `double getScore(const std::string&)`, still present as a `NOT_TESTABLE` class-test section | `score(&self, &str) -> Option<f64>` | native |
| — | `has_feature`, `has_precursor_feature`, `precursor_features`, `add_feature_with`, `add_precursor_feature_with`, `DuplicateKeyPolicy` | native |

### `KERNEL/MRMTransitionGroup.h`

| Source member | Native counterpart | Status |
| --- | --- | --- |
| `typedef std::vector<MRMFeature> MRMFeatureListType` | `Vec<MRMFeature>` used directly | not ported: alias only |
| `typedef std::vector<TransitionType> TransitionsType` | `Vec<T>` used directly | not ported: alias only |
| `typedef typename ChromatogramType::PeakType PeakType` | — | not ported: the group never uses it, and the `NativeIdentified` bound deliberately does not reach into the raw-data type's peaks |
| `MRMTransitionGroup()` | `MRMTransitionGroup::new()` / `default()` | ported |
| `MRMTransitionGroup(const MRMTransitionGroup&)` | `Clone` | ported |
| `virtual ~MRMTransitionGroup()` | `Drop` glue | ported |
| `MRMTransitionGroup& operator=(const MRMTransitionGroup&)` | `Clone` plus assignment | ported |
| `Size size() const` | `size(&self) -> usize` | ported |
| `const std::string& getTransitionGroupID() const` | `transition_group_id(&self) -> &str` | ported |
| `void setTransitionGroupID(const std::string&)` | `set_transition_group_id(&mut self, impl Into<String>)` | ported |
| `const std::vector<TransitionType>& getTransitions() const` | `transitions(&self) -> &[T]` | ported |
| `std::vector<TransitionType>& getTransitionsMuteable()` | `transitions_mut(&mut self) -> &mut [T]` | ported with a narrower return; see native differences |
| `void addTransition(const TransitionType&, const std::string&)` | `add_transition(&mut self, T, impl AsRef<str>) -> Result<()>` | ported |
| `bool hasTransition(const std::string&) const` | `has_transition(&self, &str) -> bool` | ported |
| `const TransitionType& getTransition(const std::string&)` | `transition(&self, &str) -> Result<&T>` | ported |
| `std::vector<ChromatogramType>& getChromatograms()` | `chromatograms_mut(&mut self) -> &mut [C]` | ported with a narrower return |
| `const std::vector<ChromatogramType>& getChromatograms() const` | `chromatograms(&self) -> &[C]` | ported |
| `void addChromatogram(const ChromatogramType&, const std::string&)` | `add_chromatogram(&mut self, C, impl AsRef<str>) -> Result<()>` | ported |
| `bool hasChromatogram(const std::string&) const` | `has_chromatogram(&self, &str) -> bool` | ported |
| `ChromatogramType& getChromatogram(const std::string&)` | `chromatogram_mut(&mut self, &str) -> Result<&mut C>` | ported |
| `const ChromatogramType& getChromatogram(const std::string&) const` | `chromatogram(&self, &str) -> Result<&C>` | ported |
| `std::vector<ChromatogramType>& getPrecursorChromatograms()` | `precursor_chromatograms_mut(&mut self) -> &mut [C]` | ported with a narrower return |
| `const std::vector<ChromatogramType>& getPrecursorChromatograms() const` | `precursor_chromatograms(&self) -> &[C]` | ported |
| `void addPrecursorChromatogram(const ChromatogramType&, const std::string&)` | `add_precursor_chromatogram(&mut self, C, impl AsRef<str>) -> Result<()>` | ported |
| `bool hasPrecursorChromatogram(const std::string&) const` | `has_precursor_chromatogram(&self, &str) -> bool` | ported |
| `ChromatogramType& getPrecursorChromatogram(const std::string&)` | `precursor_chromatogram_mut(&mut self, &str) -> Result<&mut C>` | ported |
| `const ChromatogramType& getPrecursorChromatogram(const std::string&) const` | `precursor_chromatogram(&self, &str) -> Result<&C>` | ported |
| `const std::vector<MRMFeature>& getFeatures() const` | `features(&self) -> &[MRMFeature]` | ported |
| `std::vector<MRMFeature>& getFeaturesMuteable()` | `features_mut(&mut self) -> &mut Vec<MRMFeature>` | ported, full mutability kept |
| `void addFeature(const MRMFeature&)` | `add_feature(&mut self, MRMFeature) -> Result<()>` | ported |
| `void addFeature(MRMFeature&&)` | same entry point | ported |
| `bool isInternallyConsistent() const` | `is_internally_consistent(&self) -> bool` | ported, and now actually evaluated; see native differences |
| `bool chromatogramIdsMatch() const` | `chromatogram_ids_match(&self) -> bool` | ported |
| `void getLibraryIntensity(std::vector<double>&) const` | `library_intensity(&self) -> Result<Vec<f64>>` | ported |
| `MRMTransitionGroup subset(std::vector<std::string>) const` | `subset(&self, &[String]) -> Result<Self>` | ported |
| `MRMTransitionGroup subsetDependent(std::vector<std::string>) const` | `subset_dependent(&self, &[String]) -> Result<Self>` | ported |
| `const MRMFeature& getBestFeature() const` | `best_feature(&self) -> Result<&MRMFeature>` | ported |
| *(protected)* `isMappingConsistent_()` | private `is_mapping_consistent` | ported |
| *(protected)* `tr_gr_id_`, `transitions_`, `chromatograms_`, `precursor_chromatograms_`, `mrm_features_`, `chromatogram_map_`, `precursor_chromatogram_map_`, `transition_map_` | private fields of the same shape | ported as private state |
| — | `is_empty`, `MAX_ITEMS`, `MAX_BYTES`, `NativeIdentified`, `Transition`, `SimpleTransition` | native |

### Supporting types read from `ANALYSIS/OPENSWATH/OpenSwathScores.h`

`MRMFeature` stores an `OpenSwath_Scores` by value and `IDScoresAsMetaValue`
takes an `OpenSwath_Ind_Scores`, so both data structs are ported here. The
header itself is **not** claimed as ported; it belongs to the OpenSWATH analysis
port.

| Source member | Native counterpart | Status |
| --- | --- | --- |
| `struct OpenSwath_Scores` — 62 `double` data members | `OpenSwathScores`, same field names, same per-field defaults | ported |
| `OpenSwath_Scores::get_quick_lda_score`, `calculate_lda_prescore`, `calculate_lda_single_transition`, `calculate_swath_lda_prescore` | — | not ported: the bodies live in `OpenSwathScores.cpp` and belong to the OpenSWATH analysis work package, not to the kernel |
| `struct OpenSwath_Ind_Scores` — `ind_num_transitions`, `ind_transition_names` and 40 `std::vector<double>` members | `OpenSwathIndScores`, same field names | ported |
| `struct OpenSwath_Scores_Usage` — 20 `bool` switches | — | not ported: an analysis-side switchboard that neither kernel header touches |
| — | `OpenSwathIndScores::KEY_SUFFIXES` | native: the 42 metadata suffixes `IDScoresAsMetaValue` writes |

### Doxygen accounting

Both headers are lightly documented: `@brief` on the two classes and on
`getBestFeature`, `///` one-liners on most members, two `@param[in]` tags on
`addPrecursorChromatogram`, and prose paragraphs on `addTransition`,
`addChromatogram` and `addPrecursorChromatogram` about the key convention. All
of it is carried onto the corresponding Rust items, including the class-level
sentence that a consistent structure needs the same identifiers for the
chromatograms as for the transitions and the sentence about accepting
`MSSpectrum`. There are **no** `@note`, `@exception`, `@throw`, `@warning`,
`@pre`, `@see`, `@deprecated`, `@todo` or `@code` tags in either header; the
`@name` groups became rustdoc paragraphs and the source's own inline comment on
the library-intensity clamp is carried as prose.

## Preserved source conventions

- `size()` is the **chromatogram** count, not the transition or feature count.
- The three key maps are independent. The same key may name a transition, a
  fragment-ion chromatogram and a precursor chromatogram at once, and the
  feature and precursor-feature maps inside `MRMFeature` are likewise separate.
- `addTransition`, `addChromatogram` and `addPrecursorChromatogram` refuse a
  repeated key, exactly as the source's `Exception::InvalidValue`.
- `getFeatureIDs` and `getPrecursorFeatureIDs` walk a `std::map`, so both id
  lists come out in lexicographic key order; the ported iterators do the same
  from a `BTreeMap`.
- `addScore` writes a *metadata* entry, not a field of the score record. A
  repeated score name overwrites. `score()` reads the same place back.
- `getLibraryIntensity` clamps negative library intensities to zero and leaves
  everything else, including non-finite values, alone — `NaN < 0.0` is false in
  both languages, so a `NaN` is not clamped.
- `subset` selects transitions by their **own native ID**, not by the key they
  were registered under, and re-registers both transition and chromatogram under
  that native ID. The transition and chromatogram transfers are guarded
  separately, so a transition registered under some other key is dropped while a
  chromatogram stored under the native ID is still carried.
- `subset` carries **every** precursor chromatogram regardless of `tr_ids`, and
  re-keys each one by its own native ID rather than by its stored key.
- `subset` rebuilds each peak group from intensity, retention time and metadata
  only: quality, m/z, charge, width, unique ID, convex hulls, subordinates,
  peptide identifications and the `OpenSwath_Scores` record are left at their
  defaults. Its per-transition features are copied for the selected transitions
  only; all of its precursor features are copied.
- `subsetDependent` copies each peak group **whole** — that is what "dependent"
  refers to — keeps no precursor chromatogram, and requires a chromatogram under
  each selected transition's native ID.
- `getBestFeature` compares `Feature::getOverallQuality()` with a strict `>`, so
  the first of several equal maxima wins.
- `chromatogramIdsMatch` checks the fragment-ion map first and the precursor map
  second, and short-circuits on the first mismatch.
- Every `OpenSwath_Scores` default is preserved, including the five members that
  start at `-1` rather than `0`.
- `IDScoresAsMetaValue` writes all 42 keys under the `id_target_` / `id_decoy_`
  prefix even when the corresponding list is empty, and uses the source's own
  key names, which do not always match the field they come from
  (`ind_fwhm` is written as `width_at_50`, `ind_apex_position` as
  `peak_apex_position`, `ind_intensity_ratio` as `intensity_ratio_score`).

## Native differences

- **Unknown keys are errors, not silent first elements.** The source guards
  `getTransition`, `getChromatogram` and `getPrecursorChromatogram` with
  `OPENMS_PRECONDITION`, which expands to nothing unless `OPENMS_ASSERTIONS` is
  defined (`CONCEPT/Macros.h:91`). In an ordinary release build those lookups then index
  `map[key]`, which default-inserts `0`, and return the first element or read out
  of bounds on an empty list. The port returns `Error::MissingInformation`.
  `MRMFeature::getFeature` and `getPrecursorFeature` have the same shape in their
  non-const overloads, with the added twist that the lookup **mutates** the map.
- **`is_internally_consistent` can return `false`.** The source states its three
  conditions as `OPENMS_PRECONDITION`s and then unconditionally `return true`, so
  in a release build the function cannot report an inconsistent group, and in an
  assertions build a violation throws instead of returning `false`. The port
  evaluates the three conditions — equal list lengths, equal map sizes, and every
  chromatogram key naming a transition — and returns the answer. The cost is one
  pass over the chromatogram key map.
- **Duplicate feature keys are refused by default.** `MRMFeature::addFeature`
  appends and then assigns `feature_map_[key]`, so a repeated key leaves the
  previously keyed feature in the list with no key pointing at it: the data is
  still there but unreachable except by iterating. The native default is
  `Error::InvalidValue`; `add_feature_with(..., DuplicateKeyPolicy::SourceOverwrite)`
  reproduces the source exactly, and `add_precursor_feature_with` does the same
  for the precursor list. The transition-group lists always reject, as their
  source does.
- **Non-finite scores are refused.** `MetaValue` guarantees finite floating
  members, so `add_score` and `id_scores_as_meta_value` return
  `Error::InvalidValue` where the source stores a non-finite `DataValue`
  silently.
- **Mutable list accessors return slices, not vectors.**
  `getTransitionsMuteable`, the non-const `getChromatograms` and the non-const
  `getPrecursorChromatograms` hand out the `std::vector` itself, so a caller can
  push, erase or sort it and invalidate every stored index without anything
  noticing. The port returns `&mut [T]` / `&mut [C]`: fields may be edited, the
  length cannot change. Reordering a slice or editing a native ID still breaks
  the mapping, which is what `is_internally_consistent` and
  `chromatogram_ids_match` are for. `getFeaturesMuteable` keeps full `&mut Vec`
  mutability, because the feature list has no key map to desynchronise.
- **Out-parameters become return values.** `getFeatureIDs`,
  `getPrecursorFeatureIDs` and `getLibraryIntensity` fill a caller-supplied
  vector without clearing it; the port returns an iterator or a fresh `Vec`. The
  visible consequence is in `getLibraryIntensity`, whose clamp loop runs over the
  *whole* result vector and so also clamps negative entries the caller had put
  there before the call; the port clamps only its own entries.
- **`Deref` replaces public inheritance.** `MRMFeature` holds a `pub feature:
  Feature` and dereferences to it, the shape `Feature` already uses for
  `BaseFeature`. There is no virtual dispatch and no slicing.
- **`std::string` keys become `BTreeMap<String, usize>`.** The source stores
  `int` indices into a `std::vector`; `usize` removes the signed/unsigned
  conversion the source performs with `boost::numeric_cast` and casts at the call
  sites, and `BTreeMap` keeps the source's `std::map` ordering.
- **Errors replace uncaught STL exceptions.** `subset` throws
  `std::out_of_range` from `std::map::at` when a peak group has no feature for a
  selected transition, and `subsetDependent` throws the same when a selected
  transition has no chromatogram. Both are `Error::MissingInformation` here.
- **The transition type is a trait, not a class.** See the section above; this
  is the one structural divergence forced by scope rather than chosen.

## Checked boundaries and evidence

Ceilings, all checked in a preflight before anything is allocated or mutated:

| Constant | Value | Applies to |
| --- | --- | --- |
| `MRMFeature::MAX_FEATURES` | 1 000 000 | each of the two feature lists, checked per insertion |
| `MRMFeature::MAX_SCORE_VALUES` | 4 000 000 | total values across all `OpenSwathIndScores` lists in one `id_scores_as_meta_value` call |
| `MRMFeature::MAX_BYTES` | 64 MiB | cumulative payload of one `id_scores_as_meta_value` call |
| `MRMTransitionGroup::MAX_ITEMS` | 1 000 000 | each of the four lists, and the identifier list handed to `subset` / `subset_dependent` |
| `MRMTransitionGroup::MAX_BYTES` | 256 MiB | cumulative payload a subset operation may build |

Atomicity: `subset` and `subset_dependent` build a separate group and never
touch the receiver, so a rejected call leaves it byte-identical
(`subset_leaves_the_receiver_unchanged_on_error`). `id_scores_as_meta_value`
builds the entire 42-entry block before committing it in one `extend`, so a
non-finite value in the thirty-ninth list leaves the metadata untouched
(`non_finite_scores_are_refused`). Every `add_*` checks the ceiling and the key
before pushing, so a rejected insertion stores nothing.

Evidence is **tier 3 (source review)** throughout: every literal in
`tests/mrm.rs` is transcribed from the two pinned class tests, and no C++ was
built or executed. No retained C++ output exists for these headers, and there is
no TOPP tool whose reference output would exercise them — `MRMTransitionGroupPicker`
is the single direct consumer recorded in `docs/core-sdk-coverage.json` and is
not ported. The checks the port makes that the source does not — duplicate-key
rejection, the error paths the source leaves to a compiled-out precondition, the
ceilings, and `is_internally_consistent` returning `false` — are tier 4
(Rust-only invariants), derived from the source lines anchored in the manifest.

Class-test section accounting: `MRMFeature_test.cpp` 16 sections, 16 ported, 0
mapped, 0 unaccounted; `MRMTransitionGroup_test.cpp` 28 sections, 28 ported, 0
mapped, 0 unaccounted. The nine `NOT_TESTABLE` sections across the two files are
ported as real tests rather than mapped. The upstream `subsetDependent` section
calls `subset`, not `subsetDependent`; it is ported against `subset` exactly as
written, and `subset_dependent` is covered separately by native tests.

## Known gaps

- `ReactionMonitoringTransition` (`ANALYSIS/MRM`) is not ported, so no caller can
  yet build a group from a real assay library. `SimpleTransition` is a stand-in,
  not a substitute.
- The LDA scoring methods on `OpenSwath_Scores` and the whole
  `OpenSwath_Scores_Usage` switchboard are left to the OpenSWATH analysis port.
- `MRMTransitionGroupPicker`, the only direct TOPP consumer of these headers, is
  not ported, so there is no end-to-end differential test for this module yet.
- The per-list ceilings are asserted only on the two paths a test can reach
  cheaply (`subset_refuses_too_many_identifiers`,
  `id_scores_refuse_oversized_input`); the `MAX_FEATURES` and per-list
  `MAX_ITEMS` paths share the same helper but filling a list to a million
  entries is not worth a test.

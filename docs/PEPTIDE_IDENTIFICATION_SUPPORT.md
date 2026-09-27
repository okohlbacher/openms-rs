# PeptideIdentification

`identification::PeptideIdentification` ports
`src/openms/include/OpenMS/METADATA/PeptideIdentification.h` and
`src/openms/source/METADATA/PeptideIdentification.cpp` at revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The record and its accessors are in
`src/identification.rs`; the members that need run mappings, consensus maps,
ProForma or USIs, the source `empty()` and hashing are in
`src/identification/peptide_identification.rs`. The returned USI type is
[USI_SUPPORT.md](USI_SUPPORT.md). The record family (hits, evidences, protein
runs) is described in [IDENTIFICATION_SUPPORT.md](IDENTIFICATION_SUPPORT.md).

Direct TOPP consumers: DatabaseFilter, IDConflictResolver, IDExtractor,
IDFilter, IDRTCalibration, IDScoreSwitcher, IDSplitter, MultiplexResolver,
NovorAdapter, PeptideIndexer, QualityControl.

Tests: `tests/peptide_identification.rs` (17),
`tests/peptide_identification_consensus.rs` (2, `consensusxml` feature),
`tests/identification.rs` and the file-format tests that round-trip the record.

## API mapping

| C++ member | Rust | Notes |
|---|---|---|
| `using SpectrumIdentification = PeptideIdentification` | `pub type SpectrumIdentification` | The header deprecates the old name in favour of this one; the Rust type keeps its name and carries no `#[deprecated]`, which would flag every use in the crate. |
| `typedef PeptideHit HitType` | not needed | Stable Rust has no inherent associated types; the hit type is `PeptideHit`. |
| default constructor | `new()` / `Default` | No hits, empty run, higher-is-better, threshold 0, RT and m/z unknown. |
| copy/move constructors, assignments, destructor | `Clone`, moves, `Drop` | |
| `operator==`, `operator!=` | derived `PartialEq` | Unset coordinates are `None == None`, the source's NaN rule. |
| `getRT`/`setRT`/`hasRT` | field `rt: Option<f64>` | `hasRT` is `rt.is_some()`; the source quiet NaN is `None`. |
| `getMZ`/`setMZ`/`hasMZ` | field `mz: Option<f64>` | Same. |
| `getHits` (both), `insertHit` (both), `setHits` (both) | field `hits` | |
| `getSignificanceThreshold`/`setSignificanceThreshold` | field `significance_threshold` | The source stores the meta value `significance_threshold` and removes it at 0. |
| `getScoreType`/`setScoreType` | field `score_type` | |
| `isHigherScoreBetter`/`setHigherScoreBetter` | field `higher_score_better` | |
| `getIdentifier`/`setIdentifier` | field `identifier` | |
| `getExperimentLabel`/`setExperimentLabel` | `experiment_label()` / `set_experiment_label()` | Empty label ignored, as in the source (below). |
| `getSpectrumReference`/`setSpectrumReference` | `spectrum_reference()` / `set_spectrum_reference()` | Lenient source text form for non-string values. |
| `static getScoreComparator(bool)` | `PeptideIdentification::score_comparator` | Returns `fn(&PeptideHit, &PeptideHit) -> bool`: `PeptideHit::ScoreMore` (`>`) or `ScoreLess` (`<`). |
| `sort()` | `sort() -> Result<()>` | Stable; refuses invalid records (NaN scores) before moving anything. |
| `empty()` | `source_empty()` | Run, hits, threshold, score type, orientation. `is_empty()` is the Rust `len() == 0` idiom. |
| `static getReferencingHits(hits, accessions)` | `PeptideIdentification::hits_referencing(&[PeptideHit], ...)`, and `referencing_hits(&self, ...)` over own hits | Borrows instead of copies. |
| `static buildUIDsFromAllPepIDs(cmap)` | `PeptideIdentification::build_uids_from_all(&ConsensusMap) -> Result<ConsensusUidMap>` | `BTreeMap<String, Vec<(Option<usize>, usize)>>`: `None` replaces the source's `Size(-1)` for unassigned identifications; per-key vectors keep multimap insertion order. |
| `static buildUIDFromPepID(pep_id, mapping)` | `build_uid(&self, &IdentifierMSRunMapper) -> Result<String>` | `Error::MissingInformation` with the source messages. |
| `buildUSI(ms_run_name, dataset_id = "local", include_interpretation = false)` | `build_usi(&self, ms_run_name, dataset_id, include_interpretation) -> Result<Usi>` | `Usi::LOCAL_COLLECTION` is the default dataset. |
| `buildUSI(mapping, dataset_id = "local", include_interpretation = false)` | `build_usi_with_mapping(&self, &IdentifierMSRunMapper, ...)` | |
| `std::hash<PeptideIdentification>` | `impl Hash` (+ `impl Eq`) | Consistent with `==`; see below. |
| protected `id_`, `hits_`, `score_type_`, `higher_score_better_`, `mz_`, `rt_` | public fields | |
| `MetaInfoInterface` base | field `metadata` | |

Native additions: `len`, `is_empty`, `validate`, `best_hit`,
`MAX_UID_IDENTIFICATIONS`, `ConsensusIdLocation`, `ConsensusUidMap`, and `Hash`
for `PeptideHit`, `PeakAnnotation` and `AnalysisResult` (needed to hash the
hits; `PeptideHit`'s own `std::hash` belongs to `PeptideHit.h`).

Every public member of the header is mapped; the status is `complete`.

## Preserved source conventions

1. **`setExperimentLabel("")` is a no-op.** The source stores only a non-empty
   label, so an empty one leaves an earlier label in place. The port used to
   remove the entry; it now matches the executed Release build (case `labels`),
   and the assertion in `tests/identification.rs` that pinned the old
   behaviour was corrected with that evidence.
2. **Lenient reference text.** `getSpectrumReference` and `getExperimentLabel`
   stringify any value type with `StringUtils::toStr(DataValue)`: an integer
   `17`, a double `1.5` or `1.23456789125e08`, a list `[a, b]`. The port uses
   the same text form (`ListFormat` for `MetaValue`) instead of Rust's float
   display, which would give `123456789.125`.
3. **UIDs.** A run with exactly one MS run path gives `<path>|<reference>`;
   otherwise a `map_index` meta value gives `<map_index text>|<reference>` (also
   for identifiers missing from the mapping, which have no paths); otherwise
   `MissingInformation`. An empty or integer reference is written as its text.
4. **Consensus UID order.** Features in map order, then the unassigned
   identifications; keys ascend bytewise, equal keys keep that visiting order.
5. **USI construction.** The scan number comes from the spectrum reference via
   `USI::extractScanNumberFromNativeID`; without one a non-empty reference
   becomes a `nativeId` USI and an empty one the invalid default. The
   interpretation is the **first** hit (not necessarily the best: call `sort`
   first), converted with `ProForma::fromAASequence` and written losslessly
   with `/charge` unless the charge is 0 (`EM[UNIMOD:35]K/2`,
   `[UNIMOD:1]-PEPTIDEK`, `PEPTIDEK-[UNIMOD:2]/-2`, `PEPTIDE[+42.0106]K/2`).
6. **Mapped USIs** take `getPrimaryMSRunPath` (the path at `id_merge_index`,
   default 0, else the legacy `base_name`) and its basename; an empty path or an
   index beyond the run's paths without `base_name` gives the invalid USI, and a
   negative index is an error (the source throws `ConversionError`).
7. **Stable sort** with a strict comparator: equal scores keep input order.

## Native differences

* Hashing uses Rust's `Hash` with the caller's `Hasher`, not the source's
  FNV-1a/golden-ratio mixing (as for every native metadata type, see
  [METADATA_HASH_SUPPORT.md](METADATA_HASH_SUPPORT.md)); absent coordinates hash
  as their own `Option` variant where the source uses the sentinels
  `0xDEADBEEF`/`0xCAFEBABE`. The two float zeros hash alike, as they compare
  equal. `AASequence` has no `Hash`, so a hit contributes its unmodified residue
  string: hits differing only in modifications collide, equal hits never
  differ.
* `Eq` is implemented so identifications can key hash containers, as the class
  test's `unordered_set` does; it is reflexive for every value `validate`
  accepts.
* `build_uids_from_all` is preflighted against `MAX_UID_IDENTIFICATIONS` (10
  million) before any key is built, and returns nothing on error.
* `sort` validates first; the source would sort NaN scores with a comparator
  that is not a strict weak ordering.

## Evidence

Tier 1 (executed differential): `../oracle/B8b-peptide-id/driver.cpp`,
compiled with libOpenMS's Release flags against the Linux x86_64 Release
install `openms4-release-bc9cc12-c19e494-174b576` and run twice on ibminode06,
byte-identical. Cases and the tests replaying them:

| Case | What | Test |
|---|---|---|
| `uids` | all 109 multimap entries of `buildUIDsFromAllPepIDs` over the class-test `MQEvidence_2.consensusXML`, three `buildUIDFromPepID` values, 34 mapped USIs with interpretations | `tests/peptide_identification_consensus.rs` |
| `uid_errors` | both `MissingInformation` paths and their messages, int/string/double `map_index`, unknown runs, integer and empty references | `build_uid_uses_the_single_run_path_or_map_index_and_reports_missing_information` |
| `empty` | `empty()` over each contributing field and three that do not | `source_empty_consults_run_hits_threshold_score_type_and_orientation_only` |
| `labels` | label set/empty/int, references of five value types, `setSpectrumReference("")`, threshold meta value | `labels_and_references_are_stringified_as_the_source_does` |
| `sort` | stable order with ties both ways, comparator truth table | `sort_is_stable_and_the_comparator_is_strict` |
| `usi_build` | 27 references x with/without interpretation, 9 interpretation sequences | `build_usi_reads_scan_numbers_and_the_first_hit_as_proforma` |
| `usi_mapping` | single/merged/out-of-range/negative/unknown runs, `base_name` fallback | `build_usi_with_mapping_resolves_merged_runs_and_legacy_base_names` |

Tier 3: the `PeptideIdentification_test.cpp` literals for `buildUSI`,
`getReferencingHits` and `std::hash` are transcribed (tests named
`*_source_literals`, `referencing_hits_*`, `hash_*`). Of the file's 36
`START_SECTION`s, the constructor, assignment and plain accessor sections map
to derived traits and public fields; `empty`, `sort`, `getReferencingHits`,
both UID builders, `getSpectrumReference`, `std::hash`, `buildUSI` and its four
`[EXTRA]` sections are covered by the tier-1 cases and literals above.
`testSignificanceThresholdMetaValue` checks the source's meta-value storage of
the threshold, which the port replaces with a field (see the API table); its
observable value semantics are covered by the `labels` case.

The consensus map is loaded with a raised `max_work`: at this revision
`consensusxml::load` charges each modified peptide sequence about 3.8 million
work units (`AASequence::parse_with_budget`) and refuses this 186 KB class-test
file under its default 50-million ceiling, which the C++ reader does not have.
That is a reader-side limit, recorded here and in the lane report.

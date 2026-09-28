# OMS file adapter

Source: `src/openms/include/OpenMS/FORMAT/OMSFile.h` (86 lines),
`src/openms/source/FORMAT/OMSFile.cpp` (71 lines) and the two internal helpers
it delegates to, `FORMAT/OMSFileLoad.h`/`.cpp` (199 + 1,179 lines) and
`FORMAT/OMSFileStore.h`/`.cpp` (300 + 1,597 lines), at
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

Rust: `src/format/oms_file.rs` (the public class, options, JSON export and the
atomic database writer), `src/format/oms_file/store.rs` (`OMSFileStore`) and
`src/format/oms_file/load.rs` (`OMSFileLoad`). Tests: `tests/oms_file.rs`.
Provenance: `tests/data/oms_file_provenance.json`. Feature: `sqlite` (the
bundled SQLite of `rusqlite`, already pinned for sqMass).

Status: **complete**. All seven public operations of `OMSFile` and every
private helper of `OMSFileLoad` and `OMSFileStore` that they reach are ported,
for every schema version the source reads (1-5). The databases the port writes
have the same schema text and the same content as the databases the Release
build writes from the same input (tier 1, below).

---

## API mapping

| Source (`OMSFile.h`) | Port |
|---|---|
| `explicit OMSFile(LogType = NONE)` | `OmsFile::new(ProgressLogType)`, `OmsFile::default()` |
| protected `log_type_`; inherited `get/setLogType` | `OmsFile::log_type`, `OmsFile::set_log_type` |
| `store(filename, const IdentificationData&)` | `OmsFile::store_identifications` / `_with(OmsStoreOptions)` |
| `store(filename, const FeatureMap&)` | `OmsFile::store_feature_map(filename, &FeatureMap, Option<&IdentificationData>)` / `_with` |
| `store(filename, const ConsensusMap&)` | `OmsFile::store_consensus_map(filename, &ConsensusMap, Option<&IdentificationData>)` / `_with` |
| `load(filename, IdentificationData&)` | `OmsFile::load_identifications -> IdentificationData` / `_with(OmsLoadOptions)` |
| `load(filename, FeatureMap&)` | `OmsFile::load_feature_map -> OmsMapLoad<FeatureMap>` / `_with` |
| `load(filename, ConsensusMap&)` | `OmsFile::load_consensus_map -> OmsMapLoad<ConsensusMap>` / `_with` |
| `exportToJSON(filename_in, filename_out)` | `OmsFile::export_to_json`, and `OmsFile::json_text` for the text |
| (none) | `OmsFile::schema_version`: the `version` table's `OMSFile` value |

The remaining inherited `ProgressLogger` members (`startProgress`, ...) are the
crate's `ProgressLogger`; each load and store reports through one built from
the adapter's log type, with the source's labels and step counts.

The source's `FeatureMap`/`ConsensusMap` own an `IdentificationData`
(`getIdentificationData()`); the crate's maps reference a graph held beside
them (`docs/FEATURE_IDENTIFICATION_SUPPORT.md`), so the map stores take the
graph and the map loads return it in `OmsMapLoad::identifications`.

### Internal helpers

`OMSFileStore` is `store::Store`: `createTable_`, `createTableDataValue_DataType_`,
`createTableCVTerm_`, `storeCVTerm_`, `createTableMetaInfo_`, `storeMetaInfo_`,
`storeMetaInfos_`, `storeVersionAndDate_`, the twelve `store*_` functions of
identification data, `createTableMoleculeType_`, `createTableAppliedProcessingStep_`,
`storeAppliedProcessingStep_`, `storeScoredProcessingResults_`,
`createTableIdentifiedMolecule_`, `getDatabaseKey_`, `createTableParentMatches_`,
`storeParentMatches_`, `createTableBaseFeature_`, `storeBaseFeature_`,
`storeFeatures_`, `storeFeatureAndSubordinates_`, `anyFeaturePredicate_`,
`storeMapMetaData_`, `storeDataProcessing_`, `storeConsensusFeatures_`,
`storeConsensusColumnHeaders_`; the free functions `raiseDBError_`,
`execAndReset` and `execWithExceptionAndReset` are `sql_error` and the
exactly-one-row `insert`. `OMSFileLoad` is `load::Load`: `isEmpty_` (an
`Option<Statement>`), `makeDataValue_`, the twelve `load*_` functions of
identification data, `prepareQueryMetaInfo_`/`handleQueryMetaInfo_`,
`prepareQueryAppliedProcessingStep_`/`handleQueryAppliedProcessingStep_`,
`handleQueryParentMatch_`, `handleQueryPeakAnnotation_`,
`loadMapMetaDataTemplate_`/`loadMapMetaData_`, `loadDataProcessing_`,
`makeBaseFeature_`, `prepareQueriesBaseFeature_`, `loadFeatureAndSubordinates_`,
`loadFeatures_`, `loadConsensusColumnHeaders_`, `loadConsensusFeatures_`,
`exportToJSON` with `exportTableToJSON_` and `export_order_by_`.

---

## Preserved source conventions

- Tables are created only for record types that have records; database keys
  count from 1 in the graph's record order (compounds, then peptides, then
  oligonucleotides share `ID_IdentifiedMolecule`'s key range); features count
  from 0 depth-first, a consensus feature's handles taking the keys after it.
- Every `CREATE TABLE` text, column declaration and bound storage class is the
  source's, so `sqlite_master` and every row match (compared, below).
- Compound formulas are written in `EmpiricalFormula::toString()` order:
  symbols sorted as text, so isotopes such as `(2)H` precede `C`
  (`crate::format::mztab_m::source_formula_text`; the crate's `Display` orders
  by element).
- Meta values are stored as `DataValue::toString()` text with the type index
  plus one (`NULL` for an empty value); units are not stored.
- A processing step's actions, score types' metadata and a `ParentMatch`'s
  metadata are not stored, in the source as here.
- The version row records the schema version 5, `datetime('now')` and the core
  SDK version; the source's library build time is `NULL` (the crate records none).
- Loading reads schema versions 1 (no `enzyme_term_specificity`), before 3
  (`grouping_order`, `position` orders), before 4 (shared `DataValue` table) and
  before 5 (`FEAT_Feature` holds everything, `overall_quality`).
- Loading a feature sets its width through `set_width`, which records the
  `FWHM` meta value as the source's `setWidth` does; stored metadata follows.
- An unknown data processing action name is skipped (the source logs it); the
  port returns the names in `OmsMapLoad::skipped_actions`.
- `exportToJSON` orders every table without a listed order by `id`; the
  `FEAT_*` tables that have no `id` column make it fail ("no such column: id"),
  so no feature or consensus map can be exported. Ported as it is (tested).

## Native differences

Where the source loses or corrupts data, the default does not, and
`OmsStoreOptions::source()` / `OmsLoadOptions::source()` select the source:

| Source behaviour | Default | Option |
|---|---|---|
| `storeIdentifiedCompounds_` binds the compound's name to `smile` | writes the SMILES | `compound_smile_from_name` |
| ratios are stored under the last handle's key | under the consensus feature's key | `consensus_ratios_under_last_handle` |
| `BaseFeature(const FeatureHandle&)` drops each handle's unique ID | keeps it | `drop_handle_unique_ids` |
| a compound formula is written by `EmpiricalFormula::toString()`, without its charge | appends a nonzero charge (`+1`) | `drop_formula_charge` |
| string-list items keep the space after each comma; every round trip adds one | removes that one space | `keep_string_list_spaces` |
| the loaded file type is left unknown | restored from the stored name | `unknown_file_type` |
| ratios are read into a copy of the feature, so none are loaded | loaded (from the feature's key, or else the last handle's key) | `drop_consensus_ratios` |

Where the source fails, the port does not:

- A map without identified molecules cannot be stored by the source:
  `FEAT_BaseFeature` declares a foreign key into `ID_IdentifiedMolecule`,
  which then does not exist, and SQLite refuses every insert (verified on the
  Release build, oracle `F2-oms-noid`). The port declares that key only when the
  table exists, as the source does for `ID_ObservationMatch.adduct_id`.
- `storeScoredProcessingResults_` creates `ID_IdentifiedMolecule_AppliedProcessingStep`
  and `..._MetaInfo` without `IF NOT EXISTS` once per molecule kind, so the
  source fails when two of compounds, peptides and oligonucleotides carry steps
  or metadata; the port creates each table once.

Other differences:

- A store writes a temporary file beside the target and renames it over the
  target when complete; the source deletes the target first. A failed store
  leaves an existing file unchanged (tested).
- A load returns a new graph and map instead of adding to the caller's, so a
  failure leaves nothing half-loaded.
- A dangling key (a row referencing a key its table lacks) is an error; the
  source default-constructs an invalid reference. A reference to a record not
  in the graph passed to a store is an error; the source stores key 0.
- Counts and positions the source binds as `uint32_t` (truncating) are refused
  beyond 32 bits.
- The source's `setLogType` does not reach its helpers (they take the
  constructor's `log_type_`); the port's loads and stores use the current type.

## Checked boundaries

`MAX_TABLE_ROWS` (10,000,000) rows per load query, `MAX_SUBORDINATE_DEPTH`
(256) levels of subordinate features in loads and stores, and the
identification graph's own `GraphLimits`; `MAX_JSON_BYTES` and `MAX_TABLES`
for the JSON export. The source has none.

## Evidence

Tier 1, executed on the Linux x86_64 Release build
(`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`, host
ibminode06, gcc 14.4.0, each run twice with identical results; oracle
directories `oracle/F2-oms`, `oracle/F2-oms-noid` and `oracle/F2-oms-mztabm` beside the repository):

- The retained JSONExporter outputs are reproduced by `export_to_json`
  (TOPP_JSONExporter_1/2, FuzzyDiff contract).
- The Release build loaded and re-stored the retained `.oms` files
  `JSONExporter_protein.oms`, `JSONExporter_RNA.oms` (schema 3) and
  `IDMerger_6_input1/2.oms` (schema 1), and loaded and re-stored the files it
  wrote from the `OMSFile_test.cpp` fixtures (sections 1-6). The port loads the
  same inputs and stores them with the source options; the canonical dumps of
  both databases (schema SQL and every row) are identical.
- The canonical dump replaces foreign keys by the content of the referenced
  row. The source orders observations, observation matches, processing steps
  and parent groups by memory address, so its own keys change between round
  trips (observed: the protein file's second round trip renumbers them); under
  the dump its first and second round trips agree, and so does the port.
- The loaded feature and consensus maps match the Release build's loaded
  values field by field (`feature_loaded.tsv`, `consensus_loaded.tsv`), and the
  class test's assertions (2 features, 2 subordinates; 6 consensus features
  with 1 and 2 handles; every identification count and the Cl- adduct after a
  round trip) hold.

Tier 3: the native options (a consensus map with ratios, handle IDs, a file
type and a string list round-trips; the source options lose each), atomicity,
refusal of non-databases.

- The MzTabM/MzTabMFile class-test input `MzTabMFile_input_1.oms` (schema 1;
  83 features, 303 compounds with formulas, SMILES and InChI, 2 adducts, 312
  feature-level observation matches) round-trips as the Release build stores it
  (oracle `F2-oms-mztabm`; dump gzip-compressed, test needs the
  `file-compression` feature), and the MzTab-M export of the loaded map
  reproduces the retained `MzTabMFile_output_1.mztab` (`tests/mztab_m.rs`).
  The Release build's own second round trip of that file differs only in the
  string-list meta values, each item gaining one more leading space (issue 7).

Not covered by an executed fixture: consensus ratios, which none of the
upstream fixtures contain.

## OpenMS C++ issues found

1. `OMSFileStore::storeIdentifiedCompounds_` binds `compound.name` to `:smile`;
   the SMILES string is never stored (source reading, OMSFileStore.cpp).
2. `OMSFileStore::storeConsensusFeatures_` stores ratios under the last
   handle's key, and `OMSFileLoad::loadConsensusFeatures_` reads them into a
   copy after `push_back`, so consensus ratios never round-trip (source reading;
   no upstream fixture has ratios).
3. `OMSFileStore::storeConsensusFeatures_` converts handles with
   `BaseFeature(const FeatureHandle&)`, which does not copy the unique ID
   (observed: every handle row of the Release build's consensus.oms has ID 0).
4. A feature or consensus map without identification data cannot be stored
   (foreign key into the missing `ID_IdentifiedMolecule`; Release build:
   "no such table: main.ID_IdentifiedMolecule").
5. `OMSFileLoad::exportToJSON` cannot export any feature or consensus map
   (orders `FEAT_Feature`, `FEAT_FeatureHandle`, `FEAT_ConsensusRatio` by the
   missing `id`; Release build: "no such column: id"). The JSONExporter tool
   therefore fails on every feature OMS file.
6. `storeScoredProcessingResults_` creates the shared `ID_IdentifiedMolecule_*`
   tables once per molecule kind without `IF NOT EXISTS`; mixed data fails
   (source reading, not executed).
7. `makeDataValue_` splits string lists at `,` without removing the space
   `DataValue::toString` writes after it; every round trip grows each item
   (Release build: `[H:+:0.6,  Na:+:0.3]` becomes `[H:+:0.6,   Na:+:0.3]` on
   the second round trip of `MzTabMFile_input_1.oms`).
8. `loadMapMetaDataTemplate_` does not restore the loaded file type.
9. Keys of observations, matches, steps and parent groups depend on memory
   addresses (IdentificationData compares iterators), so the same data is
   written with different keys from one run or round trip to the next.

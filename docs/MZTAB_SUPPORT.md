# MzTab data model support

Native Rust port of the MzTab **data model**: the shared cell vocabulary of
`src/openms/include/OpenMS/FORMAT/MzTabBase.h` (385 header lines, 884
implementation lines) and the record structs plus document type of
`src/openms/include/OpenMS/FORMAT/MzTab.h` (889 header lines, 3,414
implementation lines), at source revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

| Artifact | Path |
|---|---|
| Rust module | `src/format/mztab.rs` |
| Integration test | `tests/mztab.rs` (60 tests) |
| Provenance manifest | `tests/data/mztab_provenance.json` |

This is stage 1 of the MzTab family. It is the data model only: it does not read
or write `.mzTab` files (`MzTabFile.h`), it does not carry the metabolomics
variant (`MzTabM.h`), and it does not export a `FeatureMap`, a `ConsensusMap` or
an identification run to MzTab. Those are named in [Deferred](#deferred), with
the exact source members each of them owns.

## What MzTab is

A tab-separated PSI reporting format. A file is a metadata section of `MTD`
key/value lines, followed by any of these table sections, each a header row plus
data rows in a fixed column order with arbitrary trailing `opt_…` columns:

| Prefix | Section | Source row struct |
|---|---|---|
| `PRT` | protein | `MzTabProteinSectionRow` |
| `PEP` | peptide | `MzTabPeptideSectionRow` |
| `PSM` | peptide-spectrum match | `MzTabPSMSectionRow` |
| `SML` | small molecule | `MzTabSmallMoleculeSectionRow` |
| `NUC` | nucleic acid (OpenMS extension) | `MzTabNucleicAcidSectionRow` |
| `OLI` | oligonucleotide (OpenMS extension) | `MzTabOligonucleotideSectionRow` |
| `OSM` | oligonucleotide-spectrum match (OpenMS extension) | `MzTabOSMSectionRow` |
| `COM` | comment line | recorded by line index in `MzTab::comment_rows` |

The cell vocabulary carries the format's real rules. A numeric column has four
textual states — a value, `null`, `NaN` and `Inf` — and all four are distinct
from an optional column that is absent altogether. Every cell type renders and
parses itself; `crate::format::mztab::MzTabCell` states that contract once so a
reader or writer can drive a column without knowing its concrete type.

## API mapping — `MzTabBase.h`

### `enum MzTabCellStateType` → `MzTabCellState`

| Source | Rust | Notes |
|---|---|---|
| `MZTAB_CELLSTATE_DEFAULT` | `MzTabCellState::Default` | |
| `MZTAB_CELLSTATE_NULL` | `MzTabCellState::Null` | `#[default]`, as the two numeric cells' default constructors |
| `MZTAB_CELLSTATE_NAN` | `MzTabCellState::NaN` | |
| `MZTAB_CELLSTATE_INF` | `MzTabCellState::Inf` | |
| `SIZE_OF_MZTAB_CELLTYPE` | **not ported**: a variant-count sentinel for sizing C arrays; Rust needs none | |
| — | `MzTabCellState::as_str` | native: the textual form of a non-`Default` state |

### `class MzTabDouble` → `MzTabDouble`

| Source | Rust | Notes |
|---|---|---|
| `MzTabDouble()` | `MzTabDouble::default` / `MzTabDouble::null` | `null` state, stored value `0.0` |
| `explicit MzTabDouble(double)` | `MzTabDouble::new`, `From<f64>` | |
| `set(const double&)` | `MzTabDouble::set` | enters `Default` |
| `get()` | `MzTabDouble::get -> Result<f64>` | `Error::MissingInformation` for `Exception::ElementNotFound` |
| `toCellString()` | `MzTabDouble::to_cell_string` | |
| `fromCellString(const std::string&)` | `MzTabDouble::from_cell_string -> Result<()>` | |
| `isNull()` / `setNull(bool)` | `MzTabDouble::is_null` / `set_null` | |
| `isNaN()` / `setNaN()` | `MzTabDouble::is_nan` / `set_nan` | |
| `isInf()` / `setInf()` | `MzTabDouble::is_inf` / `set_inf` | |
| `~MzTabDouble()` | compiler-generated drop | |
| `operator<` | `MzTabDouble::source_less` | value-only, state ignored; deliberately **not** `PartialOrd` |
| `operator==` | `MzTabDouble::source_equal` | value-only; Rust `==` also compares the state |
| — | `MzTabDouble::nan`, `inf`, `state`, `raw_value` | native constructors and accessors |

### `class MzTabDoubleList` → `MzTabDoubleList`

| Source | Rust | Notes |
|---|---|---|
| `MzTabDoubleList() = default` | `MzTabDoubleList::default` / `null` | |
| `isNull()` / `setNull(bool)` | `is_null` / `set_null` | null means empty; `set_null(false)` is ignored |
| `toCellString()` / `fromCellString()` | `to_cell_string` / `from_cell_string -> Result<()>` | `\|`-separated |
| `get()` | `MzTabDoubleList::get -> &[MzTabDouble]` | borrowed instead of a copied vector |
| `set(const std::vector<MzTabDouble>&)` | `MzTabDoubleList::set` | |
| `~MzTabDoubleList()` | compiler-generated drop | |
| — | `MzTabDoubleList::get_mut` | native |

### `class MzTabInteger` → `MzTabInteger`

| Source | Rust | Notes |
|---|---|---|
| `MzTabInteger()` | `MzTabInteger::default` / `null` | `null` state, stored value `0` |
| `explicit MzTabInteger(int)` | `MzTabInteger::new`, `From<i32>` | |
| `set(const Int&)` | `MzTabInteger::set` | `Int` is `i32` |
| `get()` | `MzTabInteger::get -> Result<i32>` | `Error::MissingInformation` |
| `toCellString()` / `fromCellString()` | `to_cell_string` / `from_cell_string -> Result<()>` | |
| `isNull` / `setNull` / `isNaN` / `setNaN` / `isInf` / `setInf` | same names in snake case | |
| `~MzTabInteger()` | compiler-generated drop | |
| — | `MzTabInteger::nan`, `inf`, `state`, `raw_value` | native |

### `class MzTabIntegerList` → `MzTabIntegerList`

Same shape as `MzTabDoubleList`, with `,` as the separator; `get`, `set`,
`is_null`, `set_null`, `to_cell_string`, `from_cell_string`, plus the native
`get_mut`.

### `class MzTabBoolean` → `MzTabBoolean`

| Source | Rust | Notes |
|---|---|---|
| `MzTabBoolean()` | `MzTabBoolean::default` / `null` | stored `-1` |
| `explicit MzTabBoolean(bool)` | `MzTabBoolean::new`, `From<bool>` | |
| `set(const bool&)` | `MzTabBoolean::set` | |
| `get()` returning `Int` | `MzTabBoolean::value -> i32` | keeps the source's `-1`-when-null reading |
| `isNull()` / `setNull(bool)` | `is_null` / `set_null` | `set_null(false)` stores `0` |
| `toCellString()` / `fromCellString()` | `to_cell_string` / `from_cell_string -> Result<()>` | |
| `~MzTabBoolean()` | compiler-generated drop | |
| — | `MzTabBoolean::as_bool -> Option<bool>` | native |

### `class MzTabString` → `MzTabString`

| Source | Rust | Notes |
|---|---|---|
| `MzTabString()` | `MzTabString::default` / `null` | |
| `explicit MzTabString(const std::string&)` | `MzTabString::from_text`, `From<&str>` | forwards to `set`, so it trims and maps `null` |
| `set(const std::string&)` | `MzTabString::set` | |
| `get()` | `MzTabString::get -> &str` | borrowed |
| `isNull()` / `setNull(bool)` | `is_null` / `set_null` | null means empty |
| `toCellString()` | `to_cell_string` | |
| `fromCellString()` | `MzTabString::from_cell_string` | infallible, as the source |
| `~MzTabString()` | compiler-generated drop | |

### `typedef MzTabOptionalColumnEntry` → `MzTabOptionalColumnEntry`

The source is `std::pair<std::string, MzTabString>`, commented "column name (not
null able), value (null able)". The Rust struct has `pub name: String` for
`first`, `pub value: MzTabString` for `second`, plus
`MzTabOptionalColumnEntry::new`.

### `class MzTabParameter` → `MzTabParameter`

| Source | Rust | Notes |
|---|---|---|
| `MzTabParameter()` | `MzTabParameter::default` / `null` | all four parts empty |
| `isNull()` / `setNull(bool)` | `is_null` / `set_null` | |
| `setCVLabel` / `getCVLabel` | `set_cv_label` / `cv_label` | |
| `setAccession` / `getAccession` | `set_accession` / `accession` | |
| `setName` / `getName` | `set_name` / `name` | |
| `setValue` / `getValue` | `set_value` / `value` | |
| `toCellString()` / `fromCellString()` | `to_cell_string` / `from_cell_string -> Result<()>` | |
| `~MzTabParameter()` | compiler-generated drop | |
| — | `MzTabParameter::from_parts`, `MzTabParameter::parse` | native constructor and parsing convenience |

The four getters carry `assert(!isNull())`, which is compiled out of a release
build and then returns the empty strings. The Rust accessors return the stored
`&str` unconditionally, which is the release behaviour.

### `class MzTabParameterList` → `MzTabParameterList`

`default`/`null`, `is_null`, `set_null`, `to_cell_string`,
`from_cell_string -> Result<()>`, `get -> &[MzTabParameter]`, `set`, native
`get_mut`. Separator `|`. A member spelled `null` is a conversion error.

### `class MzTabStringList` → `MzTabStringList`

| Source | Rust | Notes |
|---|---|---|
| `MzTabStringList()` | `default` / `null` | separator `\|` |
| `isNull()` / `setNull(bool)` | `is_null` / `set_null` | |
| `setSeparator(char)` | `MzTabStringList::set_separator` | needed for `ambiguity_members` and GO accessions, which use `,` |
| `toCellString()` / `fromCellString()` | `to_cell_string` / `from_cell_string -> Result<()>` | |
| `get()` / `set()` | `get -> &[MzTabString]` / `set` | |
| `~MzTabStringList()` | compiler-generated drop | |
| — | `MzTabStringList::separator`, `get_mut` | native |

### `class MzTabSpectraRef` → `MzTabSpectraRef`

| Source | Rust | Notes |
|---|---|---|
| `MzTabSpectraRef()` | `default` / `null` | run index `0`, empty reference |
| `isNull()` / `setNull(bool)` | `is_null` / `set_null` | null when index `< 1` **or** reference empty |
| `setMSFile(Size)` | `set_ms_file -> Result<()>` | `Error::InvalidValue` for `0` |
| `setSpecRef(const std::string&)` | `set_spec_ref -> Result<()>` | `Error::InvalidValue` for empty |
| `setSpecRefFile(const std::string&)` | `set_spec_ref_file -> Result<()>` | identical to `set_spec_ref`; the source pair differs only by a log line |
| `getSpecRef()` | `spec_ref -> &str` | |
| `getMSFile()` | `ms_file -> usize` | |
| `toCellString()` / `fromCellString()` | `to_cell_string` / `from_cell_string -> Result<()>` | |
| `~MzTabSpectraRef()` | compiler-generated drop | |
| — | `MzTabSpectraRef::new`, `resolved -> Option<(usize, &str)>` | native |

### Metadata structs

| Source | Rust | Members |
|---|---|---|
| `struct MzTabSoftwareMetaData` | `MzTabSoftwareMetaData` | `software`, `setting` |
| `struct MzTabSampleMetaData` | `MzTabSampleMetaData` | `description`, `species`, `tissue`, `cell_type`, `disease`, `custom` |
| `struct MzTabCVMetaData` | `MzTabCVMetaData` | `label`, `full_name`, `version`, `url` |
| `struct MzTabInstrumentMetaData` | `MzTabInstrumentMetaData` | `name`, `source`, `analyzer`, `detector` |
| `struct MzTabContactMetaData` | `MzTabContactMetaData` | `name`, `affiliation`, `email` |

Every member keeps its source name and type, with `std::map<Size, T>` becoming
`BTreeMap<usize, T>` so a writer emits indexed keys in index order.

### `class MzTabBase`

| Source | Rust | Notes |
|---|---|---|
| `MzTabBase() = default` | **not ported**: the class has no state | |
| `virtual ~MzTabBase()` | **not ported**: Rust has no inheritance, so no virtual destructor is needed | |
| `getOptionalColumnNames_` (protected template) | `crate::format::mztab::optional_column_names` plus the `MzTabOptionalColumns` trait | the template's only requirement on its argument is a public `opt_` member; the trait states that |

## API mapping — `MzTab.h`

### `class MzTabModification` → `MzTabModification`

| Source | Rust | Notes |
|---|---|---|
| `MzTabModification()` | `default` / `null` | |
| `isNull()` / `setNull(bool)` | `is_null` / `set_null` | null when no positions and a null identifier |
| `setPositionsAndParameters` | `set_positions_and_parameters` | `Vec<(usize, MzTabParameter)>` for `std::vector<std::pair<Size, MzTabParameter>>` |
| `getPositionsAndParameters` | `positions_and_parameters -> &[(usize, MzTabParameter)]` | |
| `setModificationIdentifier` | `set_modification_identifier` | |
| `getModOrSubstIdentifier` | `mod_or_subst_identifier -> &MzTabString` | source asserts `!isNull()`; this returns the stored cell |
| `toCellString()` | `to_cell_string -> Result<String>` | fails when positions are present and the identifier is null |
| `fromCellString()` | `from_cell_string -> Result<()>` | |
| `~MzTabModification()` | compiler-generated drop | |

### `class MzTabModificationList` → `MzTabModificationList`

`default`/`null`, `is_null`, `set_null`, `to_cell_string -> Result<String>`,
`from_cell_string -> Result<()>`, `get -> &[MzTabModification]`, `set`, native
`get_mut`. Separator `,`, with commas inside an unquoted `[…]` protected.

### Metadata structs

| Source | Rust | Members |
|---|---|---|
| `struct MzTabModificationMetaData` | `MzTabModificationMetaData` | `modification`, `site`, `position` |
| `struct MzTabAssayMetaData` | `MzTabAssayMetaData` | `quantification_reagent`, `quantification_mod`, `sample_ref`, `ms_run_ref` (`Vec<i32>`) |
| `struct MzTabMSRunMetaData` | `MzTabMSRunMetaData` | `format`, `location`, `id_format`, `fragmentation_method` |
| `struct MzTabStudyVariableMetaData` | `MzTabStudyVariableMetaData` | `assay_refs`, `sample_refs`, `description` |

### `class MzTabMetaData` → `MzTabMetaData`

`MzTabMetaData()` becomes `Default`, which sets `mz_tab_version` to `1.0.0`
(MzTab.cpp:301). All 37 public data members are ported as public fields with
their source names: `mz_tab_version`, `mz_tab_mode`, `mz_tab_type`, `mz_tab_id`,
`title`, `description`, `protein_search_engine_score`,
`peptide_search_engine_score`, `psm_search_engine_score`,
`smallmolecule_search_engine_score`, `nucleic_acid_search_engine_score`,
`oligonucleotide_search_engine_score`, `osm_search_engine_score`,
`sample_processing`, `instrument`, `software`, `false_discovery_rate`,
`publication`, `contact`, `uri`, `fixed_mod`, `variable_mod`,
`quantification_method`, `protein_quantification_unit`,
`peptide_quantification_unit`, `small_molecule_quantification_unit`, `ms_run`,
`custom`, `sample`, `assay`, `study_variable`, `cv`, `colunit_protein`,
`colunit_peptide`, `colunit_psm`, `colunit_small_molecule`.

### Section rows

Every public data member of all seven row structs is ported as a public field
with its source name in snake case. Two renames are worth stating explicitly:

| Source field | Rust field |
|---|---|
| `opt_` (all seven rows) | `opt` |
| `PSM_ID` (`MzTabPSMSectionRow`) | `psm_id` |

| Source | Rust | Notes |
|---|---|---|
| `MzTabProteinSectionRow::MzTabProteinSectionRow()` | `Default` | sets `,` separators on `go_terms` and `ambiguity_members` (MzTab.cpp:294) |
| `MzTabProteinSectionRow::RowCompare` | `MzTabProteinSectionRow::row_order -> Ordering` | accession text only |
| `MzTabPeptideSectionRow::RowCompare` | `MzTabPeptideSectionRow::row_order -> Ordering` | sequence, then accession |
| `MzTabPSMSectionRow::addPepEvidenceToRows` | `MzTabPSMSectionRow::add_pep_evidence_to_rows -> Result<()>` | fills `pre`, `post`, `start`, `end`, `accession` |
| `MzTabPSMSectionRow::RowCompare` | `MzTabPSMSectionRow::row_order -> Ordering` | sequence, run index, native identifier, accession |
| `MzTabSmallMoleculeSectionRow` | `MzTabSmallMoleculeSectionRow` | no comparator in the source, none here |
| `MzTabNucleicAcidSectionRow::RowCompare` | `MzTabNucleicAcidSectionRow::row_order -> Ordering` | accession text only |
| `MzTabOligonucleotideSectionRow::RowCompare` | `MzTabOligonucleotideSectionRow::row_order -> Result<Ordering>` | sequence, accession, start, end; `Result` because the source comparator calls `MzTabInteger::get()` |
| `MzTabOSMSectionRow::RowCompare` | `MzTabOSMSectionRow::row_order -> Ordering` | sequence, run index, native identifier |

The seven `typedef std::vector<…Row> …Rows` become the type aliases
`MzTabProteinSectionRows`, `MzTabPeptideSectionRows`, `MzTabPSMSectionRows`,
`MzTabSmallMoleculeSectionRows`, `MzTabNucleicAcidSectionRows`,
`MzTabOligonucleotideSectionRows` and `MzTabOSMSectionRows`.

### `class MzTab` → `MzTab`

| Source | Rust | Notes |
|---|---|---|
| `MzTab() = default` | `MzTab::default` | metadata already declares `mzTab-version 1.0.0` |
| `~MzTab() = default` | compiler-generated drop | |
| `getMetaData()` / `setMetaData()` | `MzTab::meta_data` / `set_meta_data`, or the `meta_data` field | |
| `getProteinSectionRows()` (both overloads) / `setProteinSectionRows()` | `protein_section_rows` / `set_protein_section_rows`, or the `protein_data` field | |
| `getPeptideSectionRows()` (both) / `setPeptideSectionRows()` | `peptide_section_rows` / `set_peptide_section_rows`, or `peptide_data` | |
| `getPSMSectionRows()` (both) / `setPSMSectionRows()` | `psm_section_rows` / `set_psm_section_rows`, or `psm_data` | |
| `getNumberOfPSMs()` | `MzTab::number_of_psms -> Result<usize>` | distinct `PSM_ID`s, not rows |
| `getSmallMoleculeSectionRows()` / `setSmallMoleculeSectionRows()` | `small_molecule_section_rows` / `set_small_molecule_section_rows`, or `small_molecule_data` | |
| `getNucleicAcidSectionRows()` / `setNucleicAcidSectionRows()` | `nucleic_acid_section_rows` / `set_nucleic_acid_section_rows`, or `nucleic_acid_data` | |
| `getOligonucleotideSectionRows()` / `setOligonucleotideSectionRows()` | `oligonucleotide_section_rows` / `set_oligonucleotide_section_rows`, or `oligonucleotide_data` | |
| `getOSMSectionRows()` / `setOSMSectionRows()` | `osm_section_rows` / `set_osm_section_rows`, or `osm_data` | |
| `setCommentRows()` / `getCommentRows()` | `set_comment_rows` / `comment_rows`, or the `comment_rows` field | |
| `setEmptyRows()` / `getEmptyRows()` | `set_empty_rows` / `empty_rows`, or the `empty_rows` field | |
| `getProteinOptionalColumnNames()` | `protein_optional_column_names -> Result<Vec<String>>` | |
| `getPeptideOptionalColumnNames()` | `peptide_optional_column_names -> Result<Vec<String>>` | |
| `getPSMOptionalColumnNames()` | `psm_optional_column_names -> Result<Vec<String>>` | |
| `getSmallMoleculeOptionalColumnNames()` | `small_molecule_optional_column_names -> Result<Vec<String>>` | |
| `getNucleicAcidOptionalColumnNames()` | `nucleic_acid_optional_column_names -> Result<Vec<String>>` | |
| `getOligonucleotideOptionalColumnNames()` | `oligonucleotide_optional_column_names -> Result<Vec<String>>` | |
| `getOSMOptionalColumnNames()` | `osm_optional_column_names -> Result<Vec<String>>` | |
| `static addMetaInfoToOptionalColumns()` | `crate::format::mztab::add_meta_info_to_optional_columns -> Result<()>` | free function; the source's `MetaInfoInterface` is `MetaInfo` |
| `static generateMzTabStringFromModifications()` | `modification_metadata` and `modification_metadata_with` | returns `ModificationMetaDataReport` so skipped names are visible |
| `static generateMzTabStringFromVariableModifications()` | `variable_modification_metadata` | |
| `static generateMzTabStringFromFixedModifications()` | `fixed_modification_metadata` | |
| `static exportFeatureMapToMzTab()` | **not ported**: needs `FeatureMap` export, [deferred](#deferred) to the exporter stage | |
| `static exportIdentificationsToMzTab()` | **not ported**: [deferred](#deferred) | |
| `static extractModificationList()` | **not ported**: needs `PeptideHit`, [deferred](#deferred) | |
| `static exportConsensusMapToMzTab()` | **not ported**: needs `ConsensusMap`, [deferred](#deferred) | |
| `class MzTab::IDMzTabStream` | **not ported**: [deferred](#deferred). Its whole public surface — the constructor, `getMetaData`, `getProteinOptionalColumnNames`, `getPeptideOptionalColumnNames`, `getPSMOptionalColumnNames`, `nextPRTRow`, `nextPEPRow`, `nextPSMRow` — belongs to the exporter stage | |
| `class MzTab::CMMzTabStream` | **not ported**: [deferred](#deferred), same public surface | |
| — | `MzTab::MAX_ROWS`, `MzTab::MAX_OPTIONAL_COLUMNS` | native resource ceilings |

The 22 protected static helpers of `MzTab` — `mapIDRunIdentifier2IDRunIndex_`,
`PSMSectionRowFromPeptideID_`, `peptideSectionRowFromConsensusFeature_`,
`peptideSectionRowFromFeature_`, `proteinSectionRowFromProteinHit_`,
`nextProteinSectionRowFromProteinGroup_`,
`nextProteinSectionRowFromIndistinguishableGroup_`, `addMSRunMetaData_`,
`mapBetweenMSFileNameAndMSRunIndex_`, `getQuantStudyVariables_`,
`getProteinScoreType_`, `getConsensusMapMetaValues_`, `getFeatureMapMetaValues_`,
`getIdentificationMetaValues_`, `getMSRunSpectrumIdentifierType_`,
`mapBetweenRunAndSearchEngines_`, `mapGroupsToProteins_`, `addSearchMetaData_`,
`mapIDRunFileIndex2MSFileIndex_`, `getSearchModifications_`,
`getModificationIdentifier_`, `checkSequenceUniqueness_` — and the two
file-static helpers `remapTargetDecoyPSMAndPeptideSection_` and
`remapTargetDecoyProteinSection_` are all **not ported**: they exist only to
build the export streams and are [deferred](#deferred) with them.

### Native additions

| Rust | Purpose |
|---|---|
| `MzTabCell` trait (`is_null`, `set_null`, `write_cell`, `read_cell`) | one contract for all thirteen cell types, so a reader or writer can drive a column generically |
| `parse_cell::<T>` | parse straight into a fresh cell; the source only has in-place `fromCellString` |
| `MzTabOptionalColumns` trait | replaces the protected `getOptionalColumnNames_` template's implicit requirement |
| `ModificationMetaDataReport` | carries the modification names the registry could not resolve, which the source only writes to a log |
| `MAX_CELL_BYTES`, `MAX_CELL_ITEMS`, `MAX_MODIFICATION_NAMES` | resource ceilings the source does not have |

## Preserved source conventions

- **A default numeric cell is `null`, not zero.** `MzTabDouble()` and
  `MzTabInteger()` start in `MZTAB_CELLSTATE_NULL` with a stored `0`/`0.0`, and
  `get()` refuses to publish it.
- **`null`, `NaN` and `Inf` are parsed case-insensitively after trimming**, and
  rendered exactly `null`, `NaN` and `Inf` (MzTabBase.cpp:645, 769).
- **Trimming is the source's four ASCII bytes only** — space, tab, newline,
  carriage return (StringUtils.h:395). U+00A0 and other Unicode whitespace are
  content and stay.
- **`MzTabString` normalises on store.** `set` trims, and the literal `null` in
  any letter case stores nothing, so `MzTabString("null")` is the null cell and
  `MzTabString("nullable")` is not (MzTabBase.cpp:426).
- **An empty payload is the null state** for `MzTabString`, `MzTabParameter`,
  all four list types and `MzTabModification`. `setNull(false)` is a no-op for
  those, because there is no value to restore; only `MzTabDouble`,
  `MzTabInteger` and `MzTabBoolean` act on it.
- **`MzTabBoolean` stores an `int` and `get()` returns `-1` when null**
  (MzTabBase.cpp:493), and `setNull(false)` stores `0`, not the previous value
  (MzTabBase.cpp:505).
- **`MzTabBoolean::fromCellString` compares the untrimmed text** to `0` and `1`
  even though it trims for the `null` test, so `" 1"` is a conversion error
  (MzTabBase.cpp:537).
- **`MzTabInteger::fromCellString` parses a double first** and then requires it
  to be integral, because mzTab files from external sources write `4.0` in
  integer columns (MzTabBase.cpp:673).
- **Numbers render with `StringUtils::toStr(double)`**: fifteen *fractional*
  digits with trailing zeros trimmed to at least one inside `[1e-2, 1e4)`, and
  shortest-round-trip scientific notation with a `+`-free, two-or-more-digit
  exponent outside it (NumericFormatting.h:42, 82, 126). Fifteen fractional
  digits is not fifteen significant digits, so `51.9678841193106` renders as
  `51.967884119310597`.
- **List separators are per type, not uniform**: `|` for `MzTabDoubleList`,
  `MzTabParameterList` and `MzTabStringList`'s default; `,` for
  `MzTabIntegerList` and `MzTabModificationList`; and `,` for the `go_terms` and
  `ambiguity_members` string lists of the protein row, because `|` occurs inside
  GO terms and accessions (MzTab.cpp:294).
- **`MzTabNucleicAcidSectionRow` keeps `|`** for the same two columns, because
  the source declares no constructor for it (MzTab.h:302). The inconsistency
  with the protein row is preserved.
- **An empty cell yields no list entries.** `StringUtils::split` clears its
  output and returns for an empty subject (StringUtils.h:601, 670), so parsing
  `""` into a list leaves it empty, which is the null state — not a
  one-element list holding an empty value.
- **`MzTabParameter` renders exactly four comma-space-separated fields in
  brackets**, quoting `name` and `value` only when they contain a comma
  *followed by a space* (MzTabBase.cpp:337).
- **The parameter scanner drops every `[` and `]`**, including inside quotes
  (MzTabBase.cpp:393), so `"a, [b]"` parses as `a, b`.
- **`MzTabParameterList` refuses a member spelled `null`** (MzTabBase.cpp:68).
- **`MzTabSpectraRef` is null when the run index is below one or the reference
  is empty** (MzTabBase.cpp:169), and parsing requires exactly two
  colon-separated fields (MzTabBase.cpp:248), so a native identifier containing
  a colon cannot be read back.
- **`MzTabModification` splits on `-` and requires exactly two fields**
  (MzTab.cpp:142); text with no `-` is a bare identifier (MzTab.cpp:131).
- **`MzTabModificationList` protects a comma only when the scanner is outside
  quotes and inside a bracket** (MzTab.cpp:263), so a comma inside a quoted
  parameter part still splits the list.
- **`MzTabModification::toCellString` refuses a null identifier** when the
  modification is not null (MzTab.cpp:101).
- **`MzTabMetaData` declares `mzTab-version 1.0.0`** on construction
  (MzTab.cpp:301).
- **`getNumberOfPSMs` counts distinct `PSM_ID`s**, because a PSM row is
  duplicated per parent protein (MzTab.cpp:402), and its `@note` records that it
  relies on `PSM_ID` being set.
- **Optional column names keep first-occurrence order**, not sorted order: the
  source collects into a vector and deduplicates with a linear find
  (MzTabBase.h:376).
- **`addMetaInfoToOptionalColumns` replaces spaces in a key with underscores
  and takes values verbatim** (MzTab.cpp:600), and still emits a column for a
  key the metadata does not carry, so every row of a section declares the same
  columns.
- **`addPepEvidenceToRows` writes the literal text `null` per evidence** for an
  unknown flanking residue or position, `-` for a terminal one, and a position
  **plus one** because MzTab counts from one (MzTab.cpp:566). Its empty-list
  early return clears `pre`, `post`, `start` and `end` but **not** `accession`
  (MzTab.cpp:512).
- **`generateMzTabStringFromModifications` skips a name the registry does not
  know and still advances its index** (MzTab.cpp:657, 661), leaving a gap in the
  returned map; the empty-list variants emit
  `[MS, MS:1002453, No fixed modifications searched, ]` and
  `[MS, MS:1002454, No variable modifications searched, ]` at index 1
  (MzTab.cpp:672, 688).
- **A modification metadata `site` of `X` means "any residue"**, because
  `ResidueModification` initialises `origin_` to `'X'`.

## Native differences

Each of these differs from the source on purpose, and each is documented at the
Rust item as well.

1. **Modification positions render as decimal digits.** MzTab.cpp:85 appends the
   `Size` position with `std::string::operator+=`, whose only viable overload
   for an integer is `operator+=(char)`, so position `3` is written as the byte
   `0x03` and position `65` as `A`. This port writes `3`, which is what the
   specification requires and what the source's own `fromCellString` reads back.
   **Any position-annotated modification cell therefore differs byte for byte
   from what the C++ produces.** See issue 1 below.
2. **A second parse replaces rather than appends.** Every source list parser
   pushes onto the existing vector without clearing it (for example
   MzTabBase.cpp:71, 870), so `fromCellString` twice into the same object
   concatenates. This port builds into a temporary and commits it, so a cell is
   a value and a refused parse leaves the previous content intact.
3. **Rust `==` on `MzTabDouble` compares the state too.** The source's
   `operator==` and `operator<` read only `value_` (MzTabBase.cpp:810, 815), so
   a `null` cell equals `MzTabDouble(0.0)` and sorts as zero. Those exact
   semantics stay available as `source_equal` and `source_less`; `PartialOrd` is
   deliberately not implemented, so nothing silently inherits them.
4. **Assertions became errors.** `MzTabSpectraRef::setMSFile(0)` and
   `setSpecRef("")` assert in a debug build and silently do nothing in a release
   build (MzTabBase.cpp:183, 192); they return `Error::InvalidValue` here. The
   `assert(!isNull())` on the parameter and spectra-ref getters has no error
   path — those accessors return the stored value, which is release behaviour.
5. **Out-of-range numbers are refused rather than cast.**
   `MzTabInteger::fromCellString` casts the parsed double to `int`
   (MzTabBase.cpp:679), which is undefined behaviour outside the `int` range;
   this checks the range. `MzTabSpectraRef::fromCellString` casts a negative run
   index to `Size` (MzTabBase.cpp:254), producing an index near `2^64`; this
   rejects it. A double literal that overflows is a conversion error, matching
   the source's `std::from_chars` reporting `result_out_of_range` — Rust's own
   `parse` would have returned infinity.
6. **`getModOrSubstIdentifier` borrows.** The source returns a copied
   `MzTabString`; this returns `&MzTabString`.
7. **Getters borrow instead of copying.** `MzTabDoubleList::get` and its
   siblings return `&[T]` where the source returns a copied `std::vector<T>`, so
   reading a long list is free. A `get_mut` is added for in-place edits.
8. **Skipped modification names are returned, not logged.** The source writes
   "Skipping unknown residue modification" to `OPENMS_LOG_WARN` and drops the
   name (MzTab.cpp:659); `ModificationMetaDataReport::skipped` carries it.
9. **The modification registry is a parameter.** `modification_metadata_with`
   takes `&ModificationsDB` instead of reaching for the source's
   `ModificationsDB::getInstance()` singleton; `modification_metadata` uses the
   crate's shared global.
10. **Optional column names return a `Result`.** The membership test uses an
    auxiliary ordered set instead of the source's quadratic linear find, and the
    row and column ceilings are checked. The produced sequence is identical.
11. **`MzTabOligonucleotideSectionRow::row_order` returns a `Result`.** The
    source's comparator calls `MzTabInteger::get()`, which throws for a null
    `start` or `end`; returning the error keeps that visible instead of
    inventing an order for an unset position.
12. **No lowercased copy of a cell is built for the state tests.** The source
    lowercases the whole cell with `std::tolower` per `unsigned char`, which is
    locale-dependent above 0x7F and can corrupt UTF-8. This compares
    ASCII-case-insensitively, which agrees on every input for which the test can
    succeed, and never rewrites non-ASCII bytes.
13. **The metadata section and every indexed key use `BTreeMap`**, so iteration
    is in index order. `std::map` is ordered too, so this matches; it is called
    out because the writer stage depends on it.
14. **`opt_` became `opt` and `PSM_ID` became `psm_id`**, and the row vectors
    and metadata are public fields rather than getter/setter pairs, because the
    source accessors only return references.
15. **No threads.** The source does not parallelise this header, so there is no
    performance gap to record here; the exporters it feeds do not use OpenMP
    either.

## Checked boundaries and evidence

### Resource ceilings

| Constant | Value | Guards |
|---|---|---|
| `MAX_CELL_BYTES` | 4 MiB | every `from_cell_string`/`read_cell`, before any split or allocation |
| `MAX_CELL_ITEMS` | 100,000 | the separated entry count of every list cell, and the evidence list of `add_pep_evidence_to_rows` |
| `MAX_MODIFICATION_NAMES` | 100,000 | the input list of the three modification-metadata generators |
| `MzTab::MAX_ROWS` | 10,000,000 | `number_of_psms` and `optional_column_names` |
| `MzTab::MAX_OPTIONAL_COLUMNS` | 100,000 | distinct optional columns per section, and `add_meta_info_to_optional_columns` |

Every ceiling is checked in a preflight before anything is allocated or mutated,
and every parser builds into a temporary and commits it, so a refusal leaves the
receiver exactly as it was. `tests/mztab.rs` asserts that for the list ceiling,
the optional-column ceiling and a refused list member.

### No panics on untrusted input

Cell text arrives from a file. No public entry point indexes, slices or does
unchecked arithmetic on it:

- Strings are trimmed with `trim_matches` over characters and split with
  `str::split`, never byte-sliced. The one place an index is taken —
  `MzTabModification::from_cell_string` locating `[` — uses `str::find`, whose
  result is a character boundary of that string, and `[` is single-byte, so both
  halves are valid slices. The comma scanners of `MzTabParameter` and
  `MzTabModificationList` accumulate characters into an owned `String` rather
  than slicing, which is also how the port avoids the source's in-place
  ASCII-bell substitution (MzTab.cpp:265).
- Every number is parsed with `str::parse` and range-checked; `+ 1` on a
  peptide position uses `checked_add`.
- `tests/mztab.rs` feeds non-ASCII text to `MzTabString`, `MzTabParameter`,
  `MzTabDouble`, `MzTabModification` and `MzTabStringList`, including a
  multi-byte character adjacent to a list separator, a fullwidth `ＮＵＬＬ` that
  must not read as the null keyword, and a U+00A0-padded value that must not be
  trimmed.

### Evidence tier

**Tier 3 (source review), with the class test transcribed.** `MzTab_test.cpp`
is 307 lines and 7 `START_SECTION`s for 5,572 lines of C++, so source review
carries most of the weight and the API table above is the primary deliverable.
No C++ was built or executed and no C++ output was retained, so this is **not** a
tier 1 differential. Transcribed C++ literals are tier-3 evidence.

Class-test accounting — 7 sections, 5 ported, 0 mapped, 2 unaccounted:

| Section | Assertion macros | Status |
|---|---|---|
| `MzTab()` | 1 | **ported** — `document_default_construction_declares_version_1_0_0`. The upstream assertion is `ptr != null_ptr`; the Rust equivalent is that a default document exists and carries `mzTab-version 1.0.0`. |
| `~MzTab()` | 0 | **ported** — `document_destruction_releases_a_populated_document`. The section only calls `delete`; drop is compiler-generated, so the test drops a populated document while an independent clone stays valid. |
| `std::vector<std::string> getPSMOptionalColumnNames() const` | 2 | **ported** — `psm_optional_column_names_from_the_upstream_two_row_fixture`, with both literals (`rows.size() == 2`, `optional_columns.size() == 5`) and every cell of the fixture. |
| `static void addMetaInfoToOptionalColumns(...)` | 7 | **ported** — `add_meta_info_to_optional_columns_matches_the_upstream_seven_assertions`, all seven literals including `[0.5, 1.4, -2.0, 0.1]`. |
| `[EXTRA] exportIdentificationsToMzTab terminates with export_all_psms on empty-hit PeptideIdentification` | 1 | **unaccounted** — it calls `MzTab::exportIdentificationsToMzTab`, which this package does not port. The regression it pins lives in `IDMzTabStream::nextPSMRow`, whose advance condition underflowed for an empty hit list; nothing in `src/format/mztab.rs` contains that loop, so there is no Rust behaviour to assert. It must be ported with the exporter stage. |
| `[EXTRA] MzTabBoolean setNull / isNull polarity` | 5 | **ported** — `boolean_set_null_polarity`, all five literals. |
| `[EXTRA] consensus-map assays use fraction-group/label grain and samples are study variables` | 18 | **unaccounted** — it calls `MzTab::exportConsensusMapToMzTab` and asserts assay/study-variable grain and abundance placement produced by `CMMzTabStream`. Those are 3,400 lines of `MzTab.cpp` plus `ConsensusMap` plumbing that this package does not port. It is above the five-macro threshold and must be ported, not mapped, when the exporter stage lands; it is listed here as unaccounted rather than claimed. |

Beyond the class test, `tests/mztab.rs` has 60 tests covering every cell type's
render/parse round trip, the null/`NaN`/`Inf` state machine, the source number
convention, every list separator, the two scanners, `MzTabSpectraRef`'s null
rule and colon restriction, the modification position and `CHEMMOD` behaviour,
all seven row structs' full column sets, every metadata struct, the document
accessors, `number_of_psms`, the seven optional-column readers, the three
modification-metadata generators, the resource ceilings and the non-ASCII cases.

## Public API for the dependent stages

The reader/writer (`MzTabFile.h`) and the metabolomics variant (`MzTabM.h`)
should build on:

- Cell vocabulary: `crate::format::mztab::{MzTabCellState, MzTabCell,
  parse_cell, MzTabDouble, MzTabDoubleList, MzTabInteger, MzTabIntegerList,
  MzTabBoolean, MzTabString, MzTabOptionalColumnEntry, MzTabParameter,
  MzTabParameterList, MzTabStringList, MzTabSpectraRef, MzTabModification,
  MzTabModificationList}`
- Metadata records: `crate::format::mztab::{MzTabMetaData,
  MzTabSoftwareMetaData, MzTabSampleMetaData, MzTabCVMetaData,
  MzTabInstrumentMetaData, MzTabContactMetaData, MzTabModificationMetaData,
  MzTabAssayMetaData, MzTabMSRunMetaData, MzTabStudyVariableMetaData}`
- Section rows and their aliases: `crate::format::mztab::{MzTabProteinSectionRow,
  MzTabPeptideSectionRow, MzTabPSMSectionRow, MzTabSmallMoleculeSectionRow,
  MzTabNucleicAcidSectionRow, MzTabOligonucleotideSectionRow,
  MzTabOSMSectionRow}` plus the seven `…Rows` aliases
- Document and helpers: `crate::format::mztab::{MzTab, MzTabOptionalColumns,
  optional_column_names, add_meta_info_to_optional_columns,
  ModificationMetaDataReport, modification_metadata, modification_metadata_with,
  variable_modification_metadata, fixed_modification_metadata}`
- Ceilings: `crate::format::mztab::{MAX_CELL_BYTES, MAX_CELL_ITEMS,
  MAX_MODIFICATION_NAMES}` and `MzTab::{MAX_ROWS, MAX_OPTIONAL_COLUMNS}`

The module is not feature-gated: it needs no optional dependency and builds
under `--no-default-features`.

## C++ defects found while porting

Reported to the integrating agent for `OpenMS_CPP_ISSUES.md`; this package does
not own that file.

1. **`MzTabModification::toCellString` writes a control byte for every
   modification position.** MzTab.cpp:85,
   `pos_param_string += pos_param_pairs_[i].first;` where `first` is a `Size`.
   The only viable `std::string::operator+=` overload for an integer is the
   `char` one, so position `3` becomes the byte `0x03`, position `10` a line
   feed and position `65` the letter `A`. Every `modifications` cell that
   carries position information is malformed, the line feed truncates the row,
   and `fromCellString` cannot read any of it back. Fix: append
   `StringUtils::toStr(pos_param_pairs_[i].first)`.
2. **A `CHEMMOD` identifier with a negative mass delta cannot be read back.**
   `getModificationIdentifier_` (MzTab.cpp:1528) writes
   `"CHEMMOD:" + toStr(r.getDiffMonoMass())`, which is negative for any loss,
   while `MzTabModification::fromCellString` (MzTab.cpp:140) splits on `-` and
   requires exactly two fields. `CHEMMOD:-18.010565` is misread as position
   `CHEMMOD:` with identifier `18.010565` and then throws on the position, and
   `8-CHEMMOD:-18.010565` yields three fields and throws. Fix: split on the
   last `-` that is not part of the identifier, or delimit the identifier
   explicitly.
3. **`MzTabModificationList::fromCellString` splits inside a quoted parameter.**
   The protection test at MzTab.cpp:263 is
   `ss[pos] == ',' && !in_quotes && in_param_bracket`, so a comma inside quotes
   is *not* protected. The loop's own comment says the opposite, and the
   worked example in the comment above it,
   `3|4[a,b,,v]|8[,,"blabla, [bla]",v],1|2|3[a,b,,v]-mod:123`, splits into three
   entries instead of two. Fix: protect a comma while `in_quotes` as well.
4. **Every list cell parser appends instead of replacing.**
   `MzTabParameterList` (MzTabBase.cpp:71), `MzTabStringList` (147),
   `MzTabIntegerList` (602), `MzTabDoubleList` (870),
   `MzTabModificationList` (229, 278) and `MzTabModification`'s
   `pos_param_pairs_` (157, 167) all push onto the existing container without
   clearing it. Parsing a second cell into a reused object silently concatenates
   — and `MzTab_test.cpp` itself reuses one row object across two rows, which is
   why the second row carries nine optional entries rather than five. Fix: clear
   the container, or build into a temporary and swap.
5. **`MzTabParameter::toCellString` does not quote a name or value containing a
   bare comma.** The test at MzTabBase.cpp:337 is `hasSubstring(name_, ", ")`,
   comma *and* space. A name like `a,b` is written unquoted and produces a cell
   with five fields that `fromCellString` rejects. Fix: test for `,`.
6. **`MzTabDouble::operator==` and `operator<` ignore the cell state.**
   MzTabBase.cpp:810, 815 compare `value_` only, so the `null` cell equals
   `MzTabDouble(0.0)` and a `NaN`-state cell sorts as zero. Any container
   keyed or sorted on `MzTabDouble` conflates absent with zero. Fix: compare
   `state_` first.
7. **`MzTabSpectraRef::setSpecRefFile` is a silent duplicate of
   `setSpecRef`.** MzTabBase.cpp:215 differs from :190 only in omitting the
   warning, despite a name that suggests it sets a file rather than a spectrum.
   Fix: remove it or give it the documented behaviour.

## Deferred

- **`MzTabFile.h`** (10,609 header bytes, 127,062 implementation bytes) is the
  reader and writer, and is the next stage. It owns `load`, `store`, the header
  row layout, the `MTD` key grammar and the `opt_` column ordering.
- **`MzTabM.h`** and **`MzTabMFile.h`** are the metabolomics variant and its
  file adapter, a separate stage with its own class test (`MzTabM_test.cpp`,
  14,453 bytes).
- **The export surface of `MzTab.h`**: `exportFeatureMapToMzTab`,
  `exportIdentificationsToMzTab`, `exportConsensusMapToMzTab`,
  `extractModificationList`, the nested `IDMzTabStream` and `CMMzTabStream`
  classes and all 22 protected static helpers. These are roughly 2,900 of
  `MzTab.cpp`'s 3,414 lines and depend on `FeatureMap`, `ConsensusMap`,
  `ProteinIdentification`, `PeptideIdentification`, `PeptideHit`,
  `ExperimentalDesign` and `IDFilter`. Two `MzTab_test.cpp` sections, including
  the 18-macro consensus-map section, belong to that stage.
- **`docs/core-sdk-coverage.json` and `docs/CORE_SDK_COMPLETION.md` are stale**
  after this package: `src/format/mztab.rs` makes `MzTab.h` a candidate, so the
  generated coverage moves one header from `unmapped` to
  `evidence_requires_review`. Regenerate with
  `python3 tools/core_sdk_coverage.py --write` after adding the ledger entries.
  Both files, and the `--write` flag, are outside this package's scope.
- **`docs/doc-coverage.json` was not rewritten.** `src/format/mztab.rs`
  measures 100.0% (220/220) and `tools/check_doc_coverage.py` reports an
  improvement, but recording the new floor touches a file outside this
  package's scope.
- **CI wiring was not added.** `tests/mztab.rs` needs no features: append
  `--test mztab` to the `minimum-rust` job's no-feature line in
  `.github/workflows/rust.yml`. The test passes under `--all-features`, under
  `--no-default-features`, and under `cargo +1.85.0`.

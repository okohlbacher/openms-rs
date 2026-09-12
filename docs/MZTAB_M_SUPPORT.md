# MzTab-M support

Native Rust port of the metabolomics profile of MzTab: the data model of
`src/openms/include/OpenMS/FORMAT/MzTabM.h` (291 header lines, 745
implementation lines) and the file adapter of
`src/openms/include/OpenMS/FORMAT/MzTabMFile.h` (99 header lines, 634
implementation lines), at source revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

| Artifact | Path |
|---|---|
| Rust module | `src/format/mztab_m.rs` (2,712 lines) |
| Integration test | `tests/mztab_m.rs` (48 tests) |
| Provenance manifest | `tests/data/mztab_m_provenance.json` |
| Fixtures | `tests/data/MzTabMFile_output_1.mztab`, `tests/data/AccurateMassSearchEngine_output1_mztabm_featureXML.mzTab` |

This is stage 2 of the MzTab family. Stage 1, `docs/MZTAB_SUPPORT.md`, ported
the shared cell vocabulary of `MzTabBase.h` and the proteomics data model of
`MzTab.h`. This stage adds the MzTab-M 2.0.0-M data model, the `FeatureMap`
exporter and the writer. Neither stage ports a **reader**; the source has none
for either flavour. `MzTabFile.h`, the proteomics writer, is still unported —
this module borrows one static helper from it, named below.

The module is registered ungated in `src/format/mod.rs`: it needs
`crate::format::mztab`, `crate::format::controlled_vocabulary`,
`crate::kernel` and `crate::identification::graph`, none of which is behind a
Cargo feature, so MzTab-M is available in the `--no-default-features` build.

## What MzTab-M is, and what it shares with MzTab

MzTab-M keeps MzTab's line grammar exactly: an `MTD` metadata section of
tab-separated key/value lines, then one header row plus data rows per table
section, in a fixed column order with arbitrary trailing `opt_…` columns. It
replaces the proteomics section set with three metabolomics sections:

| Prefix | Header | Section | Row struct |
|---|---|---|---|
| `SML` | `SMH` | small-molecule summary | `MzTabMSmallMoleculeSectionRow` |
| `SMF` | `SFH` | small-molecule feature | `MzTabMSmallMoleculeFeatureSectionRow` |
| `SME` | `SEH` | small-molecule evidence | `MzTabMSmallMoleculeEvidenceSectionRow` |

The three are a chain: an `SML` row references `SMF` rows through
`SMF_ID_REFS`, and an `SMF` row references `SME` rows through `SME_ID_REFS`.

### Shared unchanged with `MzTab.h` / `MzTabBase.h`

An auditor should be able to confirm that this module defines **no cell type
and no general metadata record of its own**. Everything in this list is used
from `crate::format::mztab` verbatim:

| Shared item | Used for |
|---|---|
| `MzTabCellState`, `MzTabCell`, `parse_cell` | the three-state `null`/`NaN`/`Inf` machinery and the generic column contract |
| `MzTabString`, `MzTabInteger`, `MzTabDouble`, `MzTabBoolean` | every scalar cell |
| `MzTabStringList`, `MzTabIntegerList`, `MzTabDoubleList`, `MzTabParameterList` | every list cell |
| `MzTabParameter` | every `[CV, accession, name, value]` cell |
| `MzTabSpectraRef` | `SME spectra_ref` |
| `MzTabOptionalColumnEntry` | every `opt_…` column |
| `MzTabOptionalColumns`, `optional_column_names` | the three `get…OptionalColumnNames` methods (replacing `MzTabBase::getOptionalColumnNames_`) |
| `add_meta_info_to_optional_columns` | `MzTabM::addMetaInfoToOptionalColumns`, whose C++ body is byte-identical to `MzTab`'s |
| `MzTabSoftwareMetaData` | `MTD software[n]`, `software[n]-setting[m]` |
| `MzTabSampleMetaData` | `MTD sample[n]-…`, all six keys |
| `MzTabInstrumentMetaData` | `MTD instrument[n]-…`, all four keys |
| `MzTabContactMetaData` | `MTD contact[n]-…`, all three keys |
| `MzTabCVMetaData` | `MTD cv[n]-…`, all four keys |
| `MzTab::MAX_ROWS`, `MzTab::MAX_OPTIONAL_COLUMNS` | the row and column ceilings, re-exported as `MzTabM::MAX_ROWS` / `MzTabM::MAX_OPTIONAL_COLUMNS` |

`MzTabModification`, `MzTabModificationList`, `MzTabModificationMetaData` and
the modification-metadata generators are **not** used: MzTab-M has no
modification column and no `fixed_mod`/`variable_mod` metadata.

### Replaced by the profile

| MzTab (`MzTab.h`) | MzTab-M (`MzTabM.h`) | Why |
|---|---|---|
| `MzTabAssayMetaData` (`quantification_reagent`, `quantification_mod[m]`, `sample_ref`, `ms_run_ref` list) | `MzTabMAssayMetaData` (`name`, `custom[m]`, `external_uri`, `sample_ref`, single `ms_run_ref`) | a metabolomics assay is one sample in one run, not a label channel |
| `MzTabMSRunMetaData` (`format`, `location`, `id_format`, `fragmentation_method` as one `Param[]`) | `MzTabMMSRunMetaData` (adds `instrument_ref`, mandatory `scan_polarity[m]`, `hash`, `hash_method`; `fragmentation_method[m]` indexed) | polarity is mandatory in the profile, and the file hash is part of its provenance model |
| `MzTabStudyVariableMetaData` (`assay_refs`, `sample_refs`, `description`) | `MzTabMStudyVariableMetaData` (drops `sample_refs`; adds `name`, `average_function`, `variation_function`, `factors`) | the profile names how replicates are summarised |
| — | `MzTabMDatabaseMetaData` | a metabolomics file declares its identification databases once and the row `identifier` columns carry the declared `prefix` |
| `MzTabMetaData` | `MzTabMMetaData` | see below |
| `MzTabSmallMoleculeSectionRow` (the proteomics file's small-molecule table) | `MzTabMSmallMoleculeSectionRow` | different column set entirely; the profile row is an `SML_ID`-keyed summary over feature rows |
| `MzTabPSMSectionRow` | `MzTabMSmallMoleculeEvidenceSectionRow` | the profile's per-identification row |
| — | `MzTabMSmallMoleculeFeatureSectionRow` | MzTab has no feature section |
| `MzTab` | `MzTabM` | three sections instead of seven |

`MzTabMMetaData` versus `MzTabMetaData`, key by key:

* **dropped**: `mzTab-mode`, `mzTab-type`, `protein_search_engine_score[n]`,
  `peptide_search_engine_score[n]`, `psm_search_engine_score[n]`,
  `smallmolecule_search_engine_score[n]`,
  `nucleic_acid_search_engine_score[n]`,
  `oligonucleotide_search_engine_score[n]`, `osm_search_engine_score[n]`,
  `false_discovery_rate`, `fixed_mod[n]`, `variable_mod[n]`,
  `protein-quantification_unit`, `peptide-quantification_unit`,
  `colunit-protein`, `colunit-peptide`, `colunit-psm`.
* **added**: `external_study_uri[n]`, `database[n]`,
  `derivatization_agent[n]`, `small_molecule_feature-quantification_unit`,
  `small_molecule-identification_reliability`, `id_confidence_measure[n]`,
  `colunit-small_molecule_feature`, `colunit-small_molecule_evidence`.
* **kept**: `mzTab-version` (default `2.0.0-M`, not `1.0.0`), `mzTab-ID`,
  `title`, `description`, `sample_processing[n]`, `instrument[n]`,
  `software[n]`, `publication[n]`, `contact[n]`, `uri[n]`,
  `quantification_method`, `sample[n]`, `ms_run[n]`, `assay[n]`,
  `study_variable[n]`, `custom[n]`, `cv[n]`,
  `small_molecule-quantification_unit`, `colunit-small_molecule`.

## API mapping — `MzTabM.h`

### `struct CompareMzTabMMatchRef`

| Source member | Rust |
|---|---|
| `bool operator()(const ObservationMatchRef& lhs, const ObservationMatchRef& rhs) const` | `compare_match_by_compound(&IdentificationData, ObservationMatchId, ObservationMatchId) -> Result<Ordering>` |

A free function rather than a functor: the Rust graph resolves a record ID
through its owning graph, so the comparison needs the graph as an argument and
can fail. The source calls
`identified_molecule_var.getIdentifiedCompoundRef()`, which throws for a
peptide or oligonucleotide match; this returns
`Error::MissingInformation` there.

### `class MzTabMAssayMetaData` → `MzTabMAssayMetaData`

| Source member | Rust |
|---|---|
| `MzTabString name` | `name` |
| `std::map<Size, MzTabParameter> custom` | `custom: BTreeMap<usize, MzTabParameter>` |
| `MzTabString external_uri` | `external_uri` |
| `MzTabInteger sample_ref` | `sample_ref` |
| `MzTabInteger ms_run_ref` | `ms_run_ref` |

### `class MzTabMMSRunMetaData` → `MzTabMMSRunMetaData`

| Source member | Rust |
|---|---|
| `MzTabString location` | `location` |
| `MzTabInteger instrument_ref` | `instrument_ref` |
| `MzTabParameter format` | `format` |
| `MzTabParameter id_format` | `id_format` — stored, and written unless `MzTabMWriteOptions::omit_ms_run_id_format` |
| `std::map<Size, MzTabParameter> fragmentation_method` | `fragmentation_method: BTreeMap<usize, MzTabParameter>` |
| `std::map<Size, MzTabParameter> scan_polarity` | `scan_polarity: BTreeMap<usize, MzTabParameter>` |
| `MzTabString hash` | `hash` |
| `MzTabParameter hash_method` | `hash_method` |

### `class MzTabMStudyVariableMetaData` → `MzTabMStudyVariableMetaData`

| Source member | Rust |
|---|---|
| `MzTabString name` | `name` |
| `std::vector<int> assay_refs` | `assay_refs: Vec<i32>` |
| `MzTabParameter average_function` | `average_function` |
| `MzTabParameter variation_function` | `variation_function` |
| `MzTabString description` | `description` |
| `MzTabParameterList factors` | `factors` |

### `class MzTabMDatabaseMetaData` → `MzTabMDatabaseMetaData`

| Source member | Rust |
|---|---|
| `MzTabParameter database` | `database` |
| `MzTabString prefix` | `prefix` |
| `MzTabString version` | `version` |
| `MzTabString uri` | `uri` |

### `class MzTabMMetaData` → `MzTabMMetaData`

| Source member | Rust |
|---|---|
| `MzTabMMetaData()` | `MzTabMMetaData::new()` and `Default`, both setting `mz_tab_version` to `2.0.0-M` |
| `MzTabString mz_tab_version` | `mz_tab_version` |
| `MzTabString mz_tab_id` | `mz_tab_id` |
| `MzTabString title` | `title` |
| `MzTabString description` | `description` |
| `std::map<Size, MzTabParameterList> sample_processing` | `sample_processing` |
| `std::map<Size, MzTabInstrumentMetaData> instrument` | `instrument` |
| `std::map<Size, MzTabSoftwareMetaData> software` | `software` |
| `std::map<Size, MzTabString> publication` | `publication` |
| `std::map<Size, MzTabContactMetaData> contact` | `contact` |
| `std::map<Size, MzTabString> uri` | `uri` |
| `std::map<Size, MzTabString> external_study_uri` | `external_study_uri` |
| `MzTabParameter quantification_method` | `quantification_method` |
| `std::map<Size, MzTabSampleMetaData> sample` | `sample` |
| `std::map<Size, MzTabMMSRunMetaData> ms_run` | `ms_run` |
| `std::map<Size, MzTabMAssayMetaData> assay` | `assay` |
| `std::map<Size, MzTabMStudyVariableMetaData> study_variable` | `study_variable` |
| `std::map<Size, MzTabParameter> custom` | `custom` |
| `std::map<Size, MzTabCVMetaData> cv` | `cv` |
| `std::map<Size, MzTabMDatabaseMetaData> database` | `database` |
| `std::map<Size, MzTabParameter> derivatization_agent` | `derivatization_agent` |
| `MzTabParameter small_molecule_quantification_unit` | `small_molecule_quantification_unit` |
| `MzTabParameter small_molecule_feature_quantification_unit` | `small_molecule_feature_quantification_unit` |
| `MzTabParameter small_molecule_identification_reliability` | `small_molecule_identification_reliability` |
| `std::map<Size, MzTabParameter> id_confidence_measure` | `id_confidence_measure` |
| `std::vector<MzTabString> colunit_small_molecule` | `colunit_small_molecule` |
| `std::vector<MzTabString> colunit_small_molecule_feature` | `colunit_small_molecule_feature` |
| `std::vector<MzTabString> colunit_small_molecule_evidence` | `colunit_small_molecule_evidence` |

Every `std::map<Size, T>` is a `BTreeMap<usize, T>`, so a writer emits indices
in numeric order rather than hash order; `std::map` already did.

### `class MzTabMSmallMoleculeSectionRow` → `MzTabMSmallMoleculeSectionRow`

| Source member | Rust | Column |
|---|---|---|
| `MzTabString sml_identifier` | `sml_identifier` | `SML_ID` |
| `MzTabStringList smf_id_refs` | `smf_id_refs` | `SMF_ID_REFS` |
| `MzTabStringList database_identifier` | `database_identifier` | `database_identifier` |
| `MzTabStringList chemical_formula` | `chemical_formula` | `chemical_formula` |
| `MzTabStringList smiles` | `smiles` | `smiles` |
| `MzTabStringList inchi` | `inchi` | `inchi` |
| `MzTabStringList chemical_name` | `chemical_name` | `chemical_name` |
| `MzTabStringList uri` | `uri` | `uri` |
| `MzTabDoubleList theoretical_neutral_mass` | `theoretical_neutral_mass` | `theoretical_neutral_mass` |
| `MzTabStringList adducts` | `adducts` | `adduct_ions` |
| `MzTabString reliability` | `reliability` | `reliability` |
| `MzTabParameter best_id_confidence_measure` | `best_id_confidence_measure` | `best_id_confidence_measure` |
| `MzTabDouble best_id_confidence_value` | `best_id_confidence_value` | `best_id_confidence_value` |
| `std::map<Size, MzTabDouble> small_molecule_abundance_assay` | `small_molecule_abundance_assay` | `abundance_assay[n]` |
| `std::map<Size, MzTabDouble> small_molecule_abundance_study_variable` | `small_molecule_abundance_study_variable` | `abundance_study_variable[n]` |
| `std::map<Size, MzTabDouble> small_molecule_abundance_variation_study_variable` | `small_molecule_abundance_variation_study_variable` | `abundance_variation_study_variable[n]` |
| `std::vector<MzTabOptionalColumnEntry> opt_` | `opt` | `opt_…` |

The Rust field is `opt`, not `opt_`: the trailing underscore is a C++ naming
convention for a member, and the field is public in both. The
`MzTabOptionalColumns` trait reaches it generically.

### `class MzTabMSmallMoleculeFeatureSectionRow` → `MzTabMSmallMoleculeFeatureSectionRow`

| Source member | Rust | Column |
|---|---|---|
| `MzTabString smf_identifier` | `smf_identifier` | `SMF_ID` |
| `MzTabStringList sme_id_refs` | `sme_id_refs` | `SME_ID_REFS` |
| `MzTabInteger sme_id_ref_ambiguity_code` | `sme_id_ref_ambiguity_code` | `SME_ID_REF_ambiguity_code` |
| `MzTabString adduct` | `adduct` | `adduct_ion` |
| `MzTabParameter isotopomer` | `isotopomer` | `isotopomer` |
| `MzTabDouble exp_mass_to_charge` | `exp_mass_to_charge` | `exp_mass_to_charge` |
| `MzTabInteger charge` | `charge` | `charge` |
| `MzTabDouble retention_time` | `retention_time` | `retention_time_in_seconds` |
| `MzTabDouble rt_start` | `rt_start` | `retention_time_in_seconds_start` |
| `MzTabDouble rt_end` | `rt_end` | `retention_time_in_seconds_end` |
| `std::map<Size, MzTabDouble> small_molecule_feature_abundance_assay` | `small_molecule_feature_abundance_assay` | `abundance_assay[n]` |
| `std::vector<MzTabOptionalColumnEntry> opt_` | `opt` | `opt_…` |

### `class MzTabMSmallMoleculeEvidenceSectionRow` → `MzTabMSmallMoleculeEvidenceSectionRow`

| Source member | Rust | Column |
|---|---|---|
| `MzTabString sme_identifier` | `sme_identifier` | `SME_ID` |
| `MzTabString evidence_input_id` | `evidence_input_id` | `evidence_input_id` |
| `MzTabString database_identifier` | `database_identifier` | `database_identifier` |
| `MzTabString chemical_formula` | `chemical_formula` | `chemical_formula` |
| `MzTabString smiles` | `smiles` | `smiles` |
| `MzTabString inchi` | `inchi` | `inchi` |
| `MzTabString chemical_name` | `chemical_name` | `chemical_name` |
| `MzTabString uri` | `uri` | `uri` |
| `MzTabParameter derivatized_form` | `derivatized_form` | `derivatized_form` |
| `MzTabString adduct` | `adduct` | `adduct_ion` |
| `MzTabDouble exp_mass_to_charge` | `exp_mass_to_charge` | `exp_mass_to_charge` |
| `MzTabInteger charge` | `charge` | `charge` |
| `MzTabDouble calc_mass_to_charge` | `calc_mass_to_charge` | `theoretical_mass_to_charge` |
| `MzTabSpectraRef spectra_ref` | `spectra_ref` | `spectra_ref` |
| `MzTabParameter identification_method` | `identification_method` | `identification_method` |
| `MzTabParameter ms_level` | `ms_level` | `ms_level` |
| `std::map<Size, MzTabDouble> id_confidence_measure` | `id_confidence_measure` | `id_confidence_measure[n]` |
| `MzTabInteger rank` | `rank` | `rank` |
| `std::vector<MzTabOptionalColumnEntry> opt_` | `opt` | `opt_…` |

### Typedefs

| Source | Rust |
|---|---|
| `typedef std::vector<MzTabMSmallMoleculeSectionRow> MzTabMSmallMoleculeSectionRows` | `type MzTabMSmallMoleculeSectionRows = Vec<MzTabMSmallMoleculeSectionRow>` |
| `typedef std::vector<MzTabMSmallMoleculeFeatureSectionRow> MzTabMSmallMoleculeFeatureSectionRows` | `type MzTabMSmallMoleculeFeatureSectionRows = Vec<…>` |
| `typedef std::vector<MzTabMSmallMoleculeEvidenceSectionRow> MzTabMSmallMoleculeEvidenceSectionRows` | `type MzTabMSmallMoleculeEvidenceSectionRows = Vec<…>` |

### `class MzTabM : public MzTabBase` → `MzTabM`

| Source member | Rust |
|---|---|
| `MzTabM() = default` | `MzTabM::new()` and `Default` |
| `~MzTabM() = default` | Rust `Drop` glue; nothing to port |
| `const MzTabMMetaData& getMetaData() const` | `meta_data()` and the public `meta_data` field |
| `void setMetaData(const MzTabMMetaData&)` | `set_meta_data` |
| `const MzTabMSmallMoleculeSectionRows& getMSmallMoleculeSectionRows() const` | `small_molecule_section_rows()` and the public `small_molecule_data` field |
| `void setMSmallMoleculeSectionRows(…)` | `set_small_molecule_section_rows` |
| `const MzTabMSmallMoleculeFeatureSectionRows& getMSmallMoleculeFeatureSectionRows() const` | `small_molecule_feature_section_rows()` / `small_molecule_feature_data` |
| `void setMSmallMoleculeFeatureSectionRows(…)` | `set_small_molecule_feature_section_rows` |
| `const MzTabMSmallMoleculeEvidenceSectionRows& getMSmallMoleculeEvidenceSectionRows() const` | `small_molecule_evidence_section_rows()` / `small_molecule_evidence_data` |
| `void setMSmallMoleculeEvidenceSectionRows(…)` | `set_small_molecule_evidence_section_rows` |
| `void setCommentRows(const std::map<Size, std::string>&)` | `set_comment_rows` |
| `const std::map<Size, std::string>& getCommentRows() const` | `comment_rows()` / `comment_rows` |
| `void setEmptyRows(const std::vector<Size>&)` | `set_empty_rows` |
| `const std::vector<Size>& getEmptyRows() const` | `empty_rows()` / `empty_rows` |
| `std::vector<std::string> getMSmallMoleculeOptionalColumnNames() const` | `small_molecule_optional_column_names() -> Result<Vec<String>>` |
| `std::vector<std::string> getMSmallMoleculeFeatureOptionalColumnNames() const` | `small_molecule_feature_optional_column_names() -> Result<Vec<String>>` |
| `std::vector<std::string> getMSmallMoleculeEvidenceOptionalColumnNames() const` | `small_molecule_evidence_optional_column_names() -> Result<Vec<String>>` |
| `static void addMetaInfoToOptionalColumns(keys, opt, id, meta)` | `MzTabM::add_meta_info_to_optional_columns`, forwarding to the shared `crate::format::mztab::add_meta_info_to_optional_columns` |
| `static MzTabM exportFeatureMapToMzTabM(const FeatureMap&)` | `MzTabM::export_feature_map(&FeatureMap, &IdentificationData)` and `MzTabM::export_feature_map_with(&FeatureMap, &IdentificationData, &ControlledVocabulary, &MzTabMExportOptions)` |
| `protected MzTabMMetaData m_meta_data_` | public field `meta_data` |
| `protected MzTabMSmallMoleculeSectionRows m_small_molecule_data_` | public field `small_molecule_data` |
| `protected MzTabMSmallMoleculeFeatureSectionRows m_small_molecule_feature_data_` | public field `small_molecule_feature_data` |
| `protected MzTabMSmallMoleculeEvidenceSectionRows m_small_molecule_evidence_data_` | public field `small_molecule_evidence_data` |
| `protected std::vector<Size> empty_rows_` | public field `empty_rows` |
| `protected std::map<Size, std::string> comment_rows_` | public field `comment_rows` |
| `protected std::vector<std::string> sml_optional_column_names_` | **not ported: dead member.** The implementation never reads or writes it; the three `get…OptionalColumnNames` methods recompute from the rows each call |
| `protected std::vector<std::string> smf_optional_column_names_` | **not ported: dead member** |
| `protected std::vector<std::string> sme_optional_column_names_` | **not ported: dead member** |
| `protected static std::string getAdductString_(const ObservationMatchRef&)` | private `export::adduct_string`; its behaviour is asserted through the exporter |
| `protected static void getFeatureMapMetaValues_(feature_map, &feature_keys, &match_keys, &compound_keys)` | private `export::meta_value_keys`, returning the three sets |
| inherited `template MzTabBase::getOptionalColumnNames_(rows)` | the shared `crate::format::mztab::optional_column_names` free function plus the `MzTabOptionalColumns` trait, implemented for all three profile rows |

`MzTabBase` itself is not reproduced as a type: its only member is that
protected template, which C++ needs a base class for and Rust does not.

## API mapping — `MzTabMFile.h`

| Source member | Rust |
|---|---|
| `class SVOutStream;` (forward declaration) | **not ported: unused.** The declaration is a leftover; neither the header nor the implementation mentions `SVOutStream` again |
| `MzTabMFile()` | `MzTabMFile::new()` and `Default` |
| `~MzTabMFile()` | Rust `Drop` glue; nothing to port |
| `void store(const std::string& filename, const MzTabM& mztab_m) const` | `MzTabMFile::store(impl AsRef<Path>, &MzTabM) -> Result<()>` |
| `protected void generateMzTabMMetaDataSection_(const MzTabMMetaData& map, StringList& sl) const` | `generate_meta_data_section(&MzTabMMetaData) -> Result<Vec<String>>` — public; the `sl` out-parameter is the return value |
| `protected std::string generateMzTabMSmallMoleculeHeader_(meta, optional_columns, size_t& n_columns) const` | `generate_small_molecule_header(&MzTabMMetaData, &[String]) -> Result<MzTabMSectionLine>` |
| `protected std::string generateMzTabMSmallMoleculeSectionRow_(row, optional_columns, size_t& n_columns) const` | `generate_small_molecule_section_row(&row, &MzTabMMetaData, &[String]) -> Result<MzTabMSectionLine>` |
| `protected std::string generateMzTabMSmallMoleculeFeatureHeader_(meta, optional_columns, size_t& n_columns) const` | `generate_small_molecule_feature_header(&MzTabMMetaData, &[String]) -> Result<MzTabMSectionLine>` |
| `protected std::string generateMzTabMSmallMoleculeFeatureSectionRow_(row, optional_columns, size_t& n_columns) const` | `generate_small_molecule_feature_section_row(&row, &MzTabMMetaData, &[String]) -> Result<MzTabMSectionLine>` |
| `protected std::string generateMzTabMSmallMoleculeEvidenceHeader_(meta, optional_columns, size_t& n_columns) const` | `generate_small_molecule_evidence_header(&MzTabMMetaData, &[String]) -> Result<MzTabMSectionLine>` |
| `protected std::string generateMzTabMSmallMoleculeEvidenceSectionRow_(row, optional_columns, size_t& n_columns) const` | `generate_small_molecule_evidence_section_row(&row, &MzTabMMetaData, &[String]) -> Result<MzTabMSectionLine>` |

Rust has no protected visibility, so the six generators are `pub`. Each
returns `MzTabMSectionLine { text, columns }`, replacing the `size_t&
n_columns` out-parameter. The three row generators take one extra argument the
source does not, `&MzTabMMetaData`, because a row's abundance and confidence
cells have to be aligned with the columns the header derived from the metadata;
see `MzTabMWriteOptions::source_row_abundance_cells`.

Two native additions: `MzTabMFile::with_options`, and
`MzTabMFile::generate_lines`, which returns the whole document's lines — the
`StringList` the source builds inside `store` and does not expose.

### Borrowed from the unported `MzTabFile.h`

| Source member | Rust |
|---|---|
| `static void MzTabFile::addOptionalColumnsToSectionRow_(column_names, column_entries, StringList& output)` | private `write::optional_cells`. `MzTabMFile` calls this static from all three row generators, so it had to be reproduced here; `MzTabFile.h` itself is not ported by this package |

Its contract is preserved exactly: one output cell per requested column name,
taken from the **first** row entry with that name, and the literal `null` when
the row has no such entry. An entry whose name is not in the requested list is
not written at all.

## Preserved source conventions

1. **`mzTab-version` defaults to `2.0.0-M`.** `MzTabMMetaData`'s constructor
   sets it; `MzTabMMetaData::default()` does the same. Asserted in
   `fill_data_structure` against the class test's own literal.
2. **Six metadata keys are written even when null.** `mzTab-version`,
   `mzTab-ID`, `quantification_method`,
   `small_molecule-quantification_unit`,
   `small_molecule_feature-quantification_unit` and
   `small_molecule-identification_reliability` are emitted unconditionally, so
   an unset one appears as `null`. `title`, `description` and the optional
   per-record keys are emitted only when their cell is not null. Those six are
   exactly the keys the specification marks mandatory. `assay[n]`,
   `assay[n]-ms_run_ref`, `study_variable[n]`,
   `study_variable[n]-assay_refs`, `study_variable[n]-description`,
   `ms_run[n]-location`, all four `cv[n]-…` keys, all four `database[n]-…`
   keys and every `id_confidence_measure[n]` are likewise unconditional per
   record.
3. **`assay[n]-ms_run_ref` is not guarded, `assay[n]-sample_ref` is.** A null
   `ms_run_ref` therefore writes the unusable spelling `ms_run[null]`, while a
   null `sample_ref` writes no line at all. Asserted in
   `mandatory_metadata_keys_are_written_even_when_null`.
4. **`study_variable[n]-assay_refs` is a `String[]` cell of `assay[k]`
   tokens**, `|`-separated, built on the fly from the `Vec<i32>`.
5. **One empty line precedes each of `SMH`, `SFH` and `SEH`**, and there is no
   trailing empty line. Both retained outputs carry exactly three empty lines.
6. **Every line gets one terminator**, including the last, because the source
   writes through `TextFile::store`. This port writes through the ported
   `crate::format::TextFile`, which has the same convention (and the same
   CRLF-on-Windows behaviour).
7. **The `SML` header spells the adduct column `adduct_ions`; `SFH` and `SEH`
   spell it `adduct_ion`.** The `SEH` column for `calc_mass_to_charge` is
   `theoretical_mass_to_charge`.
8. **`id_confidence_measure[n]` columns sit between `ms_level` and `rank`** in
   the `SEH` header and in every `SME` row. Asserted against the
   AccurateMassSearch fixture, which declares two.
9. **The source's number spelling.** Every `MzTabDouble` cell and every plain
   number the exporter concatenates renders through the crate's
   `format_float(value, true)`, which is `StringUtils::toStr(double)`: fifteen
   *fractional* digits with trailing zeros trimmed to at least one for
   `|value|` in `[1e-2, 1e4)`, and shortest-round-trip scientific notation with
   a `+`-free two-or-more-digit exponent outside it. That is why the retained
   files read `117.078979350899999`, `3.631019921875e04` and
   `-6.848277357391908e-04`, and all three are asserted byte for byte.
10. **`MzTabDouble(NaN)` renders `NaN`, not the `NaN` *state*.** The exporter
    stores the NaN that `ScoredProcessingResult::getScore` returns for an
    absent score, and the cell renders `NaN`. Asserted in
    `export_score_cells_are_nan_when_a_match_carries_no_score`.
11. **`getAdductString_` reformats `M+H;1+` to `[M+H]1+`** by splitting on the
    *first* semicolon only, and yields the literal `null` — which
    `MzTabString` stores as the null cell — when the match has no adduct.
12. **`evidence_input_id` is `mass=<mz>,rt=<rt>`**, built from the feature's
    own coordinates.
13. **`calc_mass_to_charge` is the compound's neutral monoisotopic mass**, not
    an m/z: the source does not correct for the adduct or the charge. Confirmed
    by the retained `SME` row, whose `theoretical_mass_to_charge` for
    `C5H11N1O2` is `117.078979350899999`.
14. **The exporter writes one `SML` row per `SMF` row.** OpenMS does not
    aggregate two features whose adducts imply one neutral mass, so an
    `[M+H]1+` and an `[M+Na]1+` trace of the same compound stay separate; the
    source comment says so explicitly.
15. **A feature without any identification still gets an `SMF` and an `SML`
    row**, with a null `SME_ID_REFS`, a null ambiguity code, and its adduct
    taken from the feature's own `adducts` meta value when it has one.
16. **`SME_ID_REF_ambiguity_code` is `1` when a feature row references more
    than one evidence, and null otherwise.**
17. **`software[n]` and `database[n]` indices start at 1** and are assigned as
    `map.size() + 1`, in the iteration order of the source's ordered sets.
18. **Processing softwares are visited in `(name, version)` order**, because
    the source holds them in a `std::set<ProcessingSoftware>` ordered by
    `Software::operator<`. The Rust graph iterates in registration order, so
    the port sorts explicitly; the order is observable because it fixes the
    `software[n]` indices. Asserted by
    `export_reproduces_retained_metadata_section`, which reproduces four
    `software[n]` lines in the retained order.
19. **`TOPP <name>` first, `analysis software` as the fallback.** The exporter
    looks up `"TOPP " + software.getName()` in the vocabulary and falls back to
    `analysis software` (`MS:1001456`) when that name is not a registered term.
    The retained output carries `analysis software` three times, for
    `AccurateMassSearch` and two versions of a second tool, and
    `TOPP FeatureFinderMetabo` (`MS:1002169`) once.
20. **`quantification_method` is always
    `LC-MS label-free quantitation analysis`.** The recognised branch fires for
    `FeatureFinderMetabo` and the fallback uses the very same term, so the
    output cannot differ.
21. **The quantification unit is guessed from
    `parameter: algorithm:mtd:quant_method` on `FeatureFinderMetabo`**:
    `area` → `MS1 feature area`, `median` → `median`, anything else →
    `MS1 feature maximum intensity`; absent → `MS1 feature area`. Both
    `small_molecule-quantification_unit` and
    `small_molecule_feature-quantification_unit` get the same value.
22. **`identification_method` and `ms_level` are guessed per identification
    tool**: `AccurateMassSearch` → `accurate mass` at MS level 1,
    `SiriusAdapter` → `de novo search` at level 2,
    `MetaboliteSpectralMatcher` → `TOPP SpecLibSearcher` at level 2, anything
    else → level 1 with the placeholder `[, , OpenMS TOPP, ]`, which is then
    replaced by `data processing action` because its CV label and accession are
    empty. `ms_level` defaults to level 1 when no tool performed
    identification. All four branches asserted in
    `export_identification_method_and_ms_level_follow_the_tool`; the
    `data processing action` branch is also what the retained `store` output
    shows.
23. **Scan polarity is positive only when the first adduct's name ends in
    `+`**, exactly as `MzTabM.cpp:287` (`at(size() - 1) == '+'`), so a name
    with no charge suffix — `M+H`, `M+Na`, a bare `H` — is negative in both.
    Polarity defaults to positive when the graph declares no adduct, where the
    source writes no `scan_polarity` at all even though the profile makes it
    mandatory; the source also warns, and this port does not log. An empty
    adduct name is the one divergence: `at(size() - 1)` underflows and throws
    `std::out_of_range` there, while there is no sign to read here, so the
    mandatory field is written as positive. Tests:
    `export_scan_polarity_follows_the_first_adduct_and_defaults_to_positive`,
    `export_scan_polarity_is_positive_only_for_a_name_ending_in_plus`.
24. **`ms_run[1]-location` is `file://`-prefixed and backslash-normalised.**
    Every `\` becomes `/`, and `file://` is prepended unless already present.
    The same normalisation applies to each `|`-separated entry of a search
    parameter's `database_location` meta value.
25. **The assay and study-variable names are
    `assay_<stem>` / `study_variable_<stem>`**, where `<stem>` is the input
    file's basename up to the first `.`, trimmed.
26. **`database[n]` accumulates.** The record is declared outside the loop, so
    a search parameter without a `database_location` keeps the *previous*
    parameter's URI rather than the `https://hmdb.ca/` default. Asserted in
    `export_database_uri_carries_over_from_an_earlier_search_parameter`.
27. **A `custom` database nulls the prefix and spaces the parameter name.**
    `db.database.contains("custom")` selects `[,, <name>, ]` with a null
    prefix; otherwise `[,,<name>, ]` with the database name as the prefix.
    Both render identically once parsed, because `MzTabParameter` trims each
    field.
28. **`reliability` defaults to `2`** — "putatively annotated compound" — and
    is overridden by a `reliability` meta value on any processing software,
    last one in `(name, version)` order winning.
29. **`IDConverter_trace*` keys are excluded** from the observation-match
    optional columns, because `IdentificationDataConverter` adds one per traced
    record.
30. **The port is serial.** `MzTabM.cpp` and `MzTabMFile.cpp` carry no
    `#pragma omp`, so there is no parallelism gap to record for this package:
    the source is serial here too.

## Native differences

Each item says what the source does and why this differs.

1. **`export_feature_map` takes the identification graph explicitly.** The
   source reads it from `FeatureMap::getIdentificationData()`; the Rust
   `FeatureMap` does not own an `IdentificationData`, its features carry
   owner-tagged `ObservationMatchId`s into a graph held alongside the map (see
   `docs/FEATURE_IDENTIFICATION_SUPPORT.md`). The signature is therefore
   `export_feature_map(&FeatureMap, &IdentificationData)`.
2. **The vocabulary is a parameter, not a file read.** The source calls
   `ControlledVocabulary::loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"))`,
   resolving a runtime share directory. This crate embeds its vocabularies, so
   `export_feature_map_with` takes a `&ControlledVocabulary` and
   `export_feature_map` passes `ControlledVocabulary::psi_ms()`. **That global
   holds all five pinned vocabularies in one object**, so its `name`, `label`,
   `version` and `url` identify the one loaded last and the `cv[1]` block it
   produces is generic. A caller that needs the source's `cv[1]` must load
   `psi-ms.obo` alone under the name `PSI-MS` and use
   `export_feature_map_with`; the tests do exactly that, from
   `resources/cv/psi-ms.obo`. Term lookups are unaffected.
3. **An empty identification graph is an error, not undefined behaviour.** The
   source states `OPENMS_PRECONDITION(!id_data.empty(), …)`, which a release
   build compiles out and then reads an empty graph. This returns
   `Error::MissingInformation`.
4. **A non-compound match is an error.** `getIdentifiedCompoundRef()` throws
   `Exception::IllegalArgument` for a peptide or oligonucleotide match; this
   returns `Error::MissingInformation` with the same meaning.
5. **A dangling `SME_ID_REFS` entry is an error.** The source finds the
   referenced evidence row with `std::find_if` and dereferences the iterator
   **without comparing it to `end()`**, which is undefined behaviour when the
   reference does not resolve. This returns `Error::MissingInformation`.
6. **An input file whose basename holds no `.` is not an error.**
   `String::prefix(char)` throws `Exception::ElementNotFound` when the
   delimiter is absent, so a feature map read from an extensionless file makes
   the source throw. This uses the whole trimmed basename. Asserted, together
   with a non-ASCII basename, in
   `export_assay_name_survives_a_basename_without_a_dot_and_non_ascii`.
7. **The five source log lines are not written.** The source logs an info line
   in `store` and warnings for an unrecognised tool name, a missing
   quantification method, a missing quantification unit, an unassessable
   identification method and a missing adduct. Each of those five paths still
   takes exactly the source's default; the defaults are listed under
   *Preserved source conventions* items 19-23 so nothing is lost, only the log.
8. **Four writer defects are opt-in.** `MzTabMWriteOptions::default()` writes
   the keys the specification names;
   `MzTabMWriteOptions::source()` reproduces the source bytes. They are
   defects, not preferences, so the library default is the correct output —
   following the `dta::WriteOptions::source()` pattern:
   * `source_assay_custom_key` — the source writes an assay's `custom[m]`
     under the key `ms_run[<assay index>]-custom[m]`.
   * `source_colunit_keys` — the source writes all three `colunit` families
     under `colunit_small_molecule`.
   * `source_derivatization_agent_key` — the source appends `-uri` to
     `derivatization_agent[n]`.
   * `omit_ms_run_id_format` — the source never emits `ms_run[n]-id_format`,
     silently discarding a set value.
9. **Row abundance and confidence cells are aligned with the header by
   default.** The header derives `abundance_assay[n]` from
   `MzTabMMetaData::assay`, `abundance_study_variable[n]` and
   `abundance_variation_study_variable[n]` from
   `MzTabMMetaData::study_variable`, and `id_confidence_measure[n]` from
   `MzTabMMetaData::id_confidence_measure`, while the source's row generators
   iterate the row's own maps. A row that carries fewer entries than the
   metadata declares therefore produces fewer cells than there are columns, and
   one that carries more produces cells under no column at all; the source's
   only check is an `OPENMS_POSTCONDITION` that a release build drops. The
   native default emits one cell per declared column, `null` where the row is
   silent, and **refuses** a row whose map names an index the metadata does not
   declare — non-lossy, because the alternative would silently drop the value.
   `source_row_abundance_cells` reproduces the source's cells exactly.
10. **Three export quirks are opt-in.** `MzTabMExportOptions::default()` is the
    native reading of each; `MzTabMExportOptions::source()` reproduces the C++:
    * `deduplicate_matches_by_compound` — the source collects a feature's
      matches into a `std::set` whose comparator orders by the identified
      compound's `identifier` alone, so matches agreeing on that identifier are
      *equivalent* and all but one are discarded, including matches that differ
      in adduct or score. Which one survives depends on the iteration order of
      the feature's own reference set, which is ordered by container-iterator
      address and so is not determined by the data. The native default keeps
      every match, ordered by identifier and then by graph ID; the source
      option keeps the lowest graph ID of each group, which is deterministic
      where the source is not.
    * `substitute_keys_before_lookup` — `getFeatureMapMetaValues_`
      substitutes `_` for spaces in every collected meta-value key, and the
      substituted key is then used for the *lookup* as well, so a key
      containing a space can never match and its column is always `null`. The
      native default collects the raw keys and lets only the column *name*
      carry the substitution.
    * `swap_cv_label_and_full_name` — the source writes
      `meta_cv.label = cv.name()` and `meta_cv.full_name = cv.label()`,
      crossing the two: a PSI-MS load produces `cv[1]-label PSI-MS` and
      `cv[1]-full_name MS`, putting the short namespace label in the full-name
      key. Both retained fixtures show it.
11. **Action-to-software names are sorted and deduplicated per action.** The
    source fills a `std::map<ProcessingAction, std::vector<std::string>>` in the
    iteration order of its processing-step set and reads it back twice: a
    membership test for `QUANTITATION`, and a loop over `IDENTIFICATION` whose
    last recognised tool wins. This port sorts and deduplicates the names, so
    the winner is the lexicographically last recognised identification tool
    rather than the last one the source's set happened to visit. The two agree
    whenever one tool performed the identification, which is the case the
    source comment assumes.
12. **Adducts are ordered by `(charge, rendered formula, name)`.** The source's
    `AdductCompare` orders by `(getCharge(), getEmpiricalFormula())`, where the
    formula comparison is `std::map`'s over the element table. The ordering is
    observable only through "the first adduct's name ends in `-` or not",
    i.e. one polarity bit, and all adducts of one run share it in practice.
13. **`MzTabSpectraRef` is left null for an empty observation identifier.** The
    source's `setSpecRef("")` logs a warning and keeps the previous, empty,
    value, which renders `null`; this leaves the cell at its null default.
    `MzTabSpectraRef::set_spec_ref` itself returns an error for empty text, as
    `docs/MZTAB_SUPPORT.md` records.
14. **The `size_t& n_columns` out-parameters become `MzTabMSectionLine`.**
15. **`generate_lines` exposes the whole rendered document**, which the source
    keeps inside `store`.
16. **The three dead `*_optional_column_names_` members are not ported.**
17. **`class SVOutStream;` is not ported.** Nothing uses it.
18. **`store` renders every line before creating the file**, so a formatting or
    ceiling error leaves an existing file untouched.
19. **The extension check maps to `Error::InvalidValue`.** The source throws
    `Exception::UnableToCreateFile`; the crate has no such variant and the
    established convention for this check is `InvalidValue` (see
    `src/format/idxml.rs`, `src/format/fasta.rs`). An *unknown* extension is
    accepted, because `FileHandler::hasValidExtension` accepts one too.
20. **Comment and empty rows are carried but not written.** The source stores
    them on `MzTabM` and its writer emits neither; they exist for a reader that
    does not yet exist in either language.

## Checked boundaries and evidence

### Resource ceilings

| Ceiling | Value | Guards |
|---|---|---|
| `MzTabM::MAX_ROWS` | 10,000,000 | each of the three sections, checked in `export_feature_map_with` before the loops and again as rows accumulate, and in `generate_lines` before any line is built |
| `MzTabM::MAX_OPTIONAL_COLUMNS` | 100,000 | `optional_column_names` (shared), `add_meta_info_to_optional_columns` (shared), and every header and row generator |
| `MzTabM::MAX_INDEXED_ENTRIES` | 100,000 | the *number of entries* in every indexed metadata map, checked in the metadata-section preflight. It does not bound the index values themselves: a key inserted at `usize::MAX` is written as `ms_run[18446744073709551615]-…`, which is what the source's unchecked `Size` does too, and nothing refuses it (`an_extreme_index_is_written_rather_than_refused`) |
| `MzTabMFile::MAX_LINES` | 1,000,000 | the estimated metadata line count and the total document line count, both before allocating, and the *physical* line count of the rendered document — a verbatim cell carrying a line break (source options only) turns one row into several, and those are charged too |

The metadata-section preflight computes an upper bound — eight lines is the
largest any one indexed record produces — checks every indexed map against
`MAX_INDEXED_ENTRIES`, checks the bound against `MAX_LINES`, and then
`try_reserve`s, so a metadata section that cannot fit fails before anything is
allocated. `run_mapping.rs` is the pattern.

Atomicity: `export_feature_map_with` builds the whole document in locals and
returns it only on success; `generate_lines` builds a `Vec<String>` and `store`
creates the file only after every line exists. `add_meta_info_to_optional_columns`
stages its entries and appends them only on success, so `opt` is unchanged on
error.

### No panics on untrusted input

* No indexing by a file-derived value. Section rows are iterated, never
  indexed; every metadata map is a `BTreeMap` reached through `get`/`keys`.
* **No string is byte-sliced unless this module constructed it.** The one place
  a borrowed name is split is `export::adduct_string`, which splits at a
  `str::find(';')` result — a character boundary by construction, and `;` is
  single-byte. `export::assay_stem` uses `split('.')` and
  `trim_matches`, both character-based.
* All arithmetic on counts uses `checked_add`/`saturating_add`; the two
  identifier counters are `i64` with `checked_add`.
* Non-ASCII coverage: `non_ascii_text_survives_every_cell_and_the_writer`
  drives Japanese and accented text through the metadata section, a list cell,
  an `opt_` column *name*, and a file whose own name is not ASCII;
  `export_assay_name_survives_a_basename_without_a_dot_and_non_ascii` exports
  with the input path `/data/日本語.mzML`.
* The module contains no `unsafe`, which the crate forbids anyway, and no
  `unwrap`/`expect` on anything derived from input.

### Section accounting

`MzTabM_test.cpp` — 342 lines, 4 sections. `MzTabMFile_test.cpp` — 59 lines,
3 sections. **7 sections, 7 ported, 0 mapped, 0 unaccounted.**

| Source section | Macros | Rust test | Status |
|---|---|---|---|
| `MzTabM_test.cpp` `START_SECTION(MzTabM())` | 1 | `default_constructor_declares_the_profile_version` | ported |
| `MzTabM_test.cpp` `START_SECTION(~MzTabM())` | 0 | `destructor_releases_a_populated_document` | ported |
| `MzTabM_test.cpp` `START_SECTION(Fill data structure)` | 21 | `fill_data_structure` | ported, all 21 macros and every literal the section sets |
| `MzTabM_test.cpp` `START_SECTION(MzTabM::exportFeatureMapToMzTabM(const FeatureMap&))` | 6 | `export_feature_map_to_mztab_m` | ported, see the caveat below |
| `MzTabMFile_test.cpp` `START_SECTION(MzTabMFile())` | 1 | `file_default_constructor` | ported |
| `MzTabMFile_test.cpp` `START_SECTION(~MzTabFile())` | 0 | `file_destructor` | ported |
| `MzTabMFile_test.cpp` `START_SECTION(void store(const std::string&, MzTabM&))` | 1 | `store_writes_the_retained_bytes` | ported, see the caveat below |

**The `.oms` caveat, stated plainly.** Both of the two sections that carry real
data load `MzTabMFile_input_1.oms`, an SQLite `.oms` file this crate has no
reader for (`OMSFile.h` is unported), so neither section's *input* can be fed
to the port.

* For `exportFeatureMapToMzTabM` the six asserted values — 83, 83, 312, 0, 18,
  6 — are verified against `tests/data/MzTabMFile_output_1.mztab`, which is the
  file `MzTabMFile_test.cpp` produced by storing the document *that very call*
  returned: the `SML`/`SMF`/`SME` line counts are the three section sizes and
  the `opt_` columns of the three header lines are the three optional-column
  lists. Every one of the six is therefore a value the C++ produced, not a
  transcribed literal. What the test does **not** do is re-run the exporter on
  the same input — and because its helper `document_from_retained` seeds the
  document with exactly those counts, the second half of the assertion is
  arithmetic over the retained file rather than a call into
  `MzTabM::export_feature_map`. A regression in the exporter is caught by
  `export_builds_one_summary_row_per_feature_row` and by the metadata
  differentials below, not by this section. The exporter is exercised separately, on synthetic graphs, in
  `export_reproduces_retained_metadata_section` (which reproduces all 25 `MTD`
  lines of that same retained file byte for byte),
  `export_reproduces_retained_ams_metadata_section` (all 24 `MTD` lines of the
  AccurateMassSearch output, including its two `id_confidence_measure[n]`
  entries) and nine further behavioural tests.
* For `store` the document is rebuilt from the retained file's own metadata and
  its three shortest data rows, stored through `MzTabMFile::store`, and the
  written bytes are compared line by line with the corresponding retained
  lines. `TEST_FILE_SIMILAR` in the source is a fuzzy comparison; this is
  exact.

### Evidence tier

**Tier 1 (executed differential against retained C++ output)** for the writer's
line grammar and for the exporter's metadata section. Two unmodified retained
outputs are committed:

* `tests/data/MzTabMFile_output_1.mztab` — 509 lines, the exact bytes
  `MzTabMFile_test.cpp`'s `store` section compares against. Used for: all 25
  `MTD` lines (twice — once from a hand-built metadata section, once from the
  exporter), all three header lines, `SML` row 3 and `SML` row 5 — the latter
  carrying two entries in each of its eight list cells, which pins the `|`
  separators and the `null|null` spelling of a two-entry list of null strings
  against the single-cell `null` of an empty list — `SMF` row 1 with all 18 of
  its optional columns, `SME` row 1 with all 6 of its optional columns, the
  three section sizes, the three optional-column counts, and the
  three-empty-line structure.
* `tests/data/AccurateMassSearchEngine_output1_mztabm_featureXML.mzTab` — the
  retained MzTab-M output of the AccurateMassSearch TOPP test, which is the one
  direct TOPP consumer `docs/core-sdk-coverage.json` records for
  `MzTabMFile.h`. Used for: all 24 `MTD` lines (twice), all three header lines
  including the two `id_confidence_measure[n]` columns, `SML` row 3, `SMF` row
  1, `SME` row 1, the `accurate mass` identification method, and the
  scientific-notation spellings `3.631019921875e04` and
  `-6.848277357391908e-04`.

Ten CV accessions the exporter selects are confirmed against the pinned
`resources/cv/psi-ms.obo` *and* against the retained outputs: `MS:1001456`,
`MS:1002169`, `MS:1001834`, `MS:1000130`, `MS:1001844`, `MS:1002896`,
`MS:1000207`, `MS:1000511`, `MS:1000543` and `MS:1000129`.

**Tier 3 (source review; transcribed C++ literals)** for the data model. Every
literal of `Fill data structure` is transcribed from `MzTabM_test.cpp:41-318`,
including `"313.168900000000008"`, `"156.0"` and
`"[MS, MS:1000752, TOPP Software, ]"`. Transcribed C++ literals are tier-3
evidence.

**Tier 4 (independently derived / Rust-only)** for the resource ceilings, their
atomicity, the non-ASCII inputs, the five native error paths, the option
matrices and the rectangularity invariant.

**The rectangularity guarantee, stated exactly.** With the default options the
output is a rectangle: one cell per declared column in every row, and no cell
may carry a tab, a carriage return or a line feed — such a cell is refused with
`Error::InvalidValue` before anything is written, because the source would pass
it through and split the row (defect 11 below). `MzTabMWriteOptions::source()`
turns both halves off and reproduces the source's output, corruption included.
Tests: `a_cell_carrying_a_tab_or_a_line_break_is_refused_by_default`,
`a_metadata_key_or_value_carrying_a_separator_is_refused_by_default`.

No C++ was built or executed for this package. The `.oms` input cannot be read,
so the exporter has no end-to-end tier-1 differential; that is the single
largest gap and the reason this package's ledger status for `MzTabM.h` is
`partial` rather than `complete`.

## Public API for the dependent stages

`crate::format::mztab_m::` —
`MzTabMAssayMetaData`, `MzTabMMSRunMetaData`, `MzTabMStudyVariableMetaData`,
`MzTabMDatabaseMetaData`, `MzTabMMetaData`,
`MzTabMSmallMoleculeSectionRow`, `MzTabMSmallMoleculeFeatureSectionRow`,
`MzTabMSmallMoleculeEvidenceSectionRow`,
`MzTabMSmallMoleculeSectionRows`, `MzTabMSmallMoleculeFeatureSectionRows`,
`MzTabMSmallMoleculeEvidenceSectionRows`,
`MzTabM` (with `MAX_ROWS`, `MAX_OPTIONAL_COLUMNS`, `MAX_INDEXED_ENTRIES`),
`MzTabMExportOptions`, `compare_match_by_compound`,
`MzTabMWriteOptions`, `MzTabMSectionLine`, `MzTabMFile` (with `MAX_LINES`).

The `AccurateMassSearch` TOPP tool is the one direct consumer of
`MzTabMFile.h` in the pinned source; it needs `MzTabM::export_feature_map` and
`MzTabMFile::store`, both of which this package provides.

## C++ defects found while porting

Eleven: the eight below, then three milder ones. All eleven are reproducible
on demand through the two option structs. They are proposed for
`OpenMS_CPP_ISSUES.md`, which this package may not edit.

1. **`MzTabMFile.cpp:230` writes an assay's custom parameter under an `ms_run`
   key.** `"MTD\tms_run[" + toStr(assay.first) + "]-custom[…"` inside the
   `md.assay` loop. The value is reported as an MS run's, under an index that
   need not name an existing run. Fix: `assay[`. Rust: opt-in through
   `MzTabMWriteOptions::source_assay_custom_key`.
2. **`MzTabMFile.cpp:355` and `:361` write the wrong `colunit` key.** All three
   `colunit` loops use the literal `colunit_small_molecule`, so
   `colunit_small_molecule_feature` and `colunit_small_molecule_evidence`
   become indistinguishable in the output. Fix: use each family's own key.
   Rust: `source_colunit_keys`.
3. **`MzTabMFile.cpp:331` appends `-uri` to `derivatization_agent[n]`.** The
   specification defines `derivatization_agent[n]`, and the value written is
   the agent parameter, not a URI. Rust:
   `source_derivatization_agent_key`.
4. **`MzTabMFile.cpp:190-200` never writes `ms_run[n]-id_format`.** It is the
   only member of `MzTabMMSRunMetaData` the metadata generator does not read,
   so a caller that sets it loses it. Rust: `omit_ms_run_id_format`.
5. **`MzTabMFile.cpp:422-435`, `:482` and `:546` emit row cells from the row's
   own maps while the header derives its columns from the metadata section.** A
   row with fewer abundances than the metadata declares produces a short row
   and an unparseable file; the guard is an `OPENMS_POSTCONDITION` that a
   release build drops. Rust: `source_row_abundance_cells`, with the native
   default aligning and refusing an undeclared index.
6. **`MzTabM.cpp:680-698` dereferences a `std::find_if` iterator without
   comparing it to `end()`.** An `SMF_ID_REFS` entry that no `SME` row carries
   is undefined behaviour. Rust: `Error::MissingInformation`.
7. **`MzTabM.cpp:572` discards observation matches.** The
   `std::set<ObservationMatchRef, CompareMzTabMMatchRef>` comparator orders by
   the identified compound's `identifier` alone, so matches agreeing on it are
   equivalent and all but one are dropped — including matches with different
   adducts or scores. Which one survives is decided by the address ordering of
   the feature's own reference set, so the output is not a function of the data.
   Rust: `deduplicate_matches_by_compound`, default off.
8. **`MzTabM.cpp:108/118/137` substitutes meta-value keys before they are used
   as lookup keys.** `getFeatureMapMetaValues_` applies
   `substitute(' ', '_')` to every collected key and
   `addMetaInfoToOptionalColumns` then calls `metaValueExists(key)` with the
   substituted spelling, so a meta value whose key contains a space always
   produces a `null` column. `MzTab.h`'s own callers pass raw keys and do not
   have this problem. Rust: `substitute_keys_before_lookup`, default off.

Three further points, defects of a milder kind:

9. **`MzTabM.cpp:287` indexes a `std::string_view` at `size() - 1`.** An adduct
   registered with an empty name underflows and `at()` throws
   `std::out_of_range` from inside an export. Rust: a match on the last
   character, which reproduces the `== '+'` test for every non-empty name and
   selects positive for the empty one rather than throwing.
10. **`MzTabM.cpp:349-351` crosses `cv[n]-label` and `cv[n]-full_name`.** See
    `swap_cv_label_and_full_name`. Both retained fixtures carry the swap, so
    fixing it upstream changes reference files.
11. **`MzTabMFile.cpp:626-631` writes every cell verbatim.** Each rendered row
    goes to `TextFile` unchanged, so a cell carrying a tab gains a column and
    one carrying a line break splits the row across physical lines — and text
    beginning with `MTD`, `SMH` or `SML` forges a line of that kind in the
    middle of a section. The text is file-derived: `chemical_name`, `uri`,
    `smiles` and every `opt_` value come from a featureXML or `.oms` text node,
    which may legally contain both characters. The source's own
    `OPENMS_POSTCONDITION` column check cannot see it, because it counts the
    cells it pushed rather than the tabs in the line it wrote. Rust:
    `MzTabMWriteOptions::source_verbatim_cells`, default off.

## Deferred

* **`OMSFile.h`.** Reading `MzTabMFile_input_1.oms` would close the exporter's
  tier-1 gap end to end. Until then the exporter's evidence is the metadata
  section plus behavioural tests.
* **`MzTabFile.h`.** The proteomics writer, 176 header lines and 2,930
  implementation lines. This package reproduces one of its statics,
  `addOptionalColumnsToSectionRow_`; when `MzTabFile.h` lands, that helper
  should move there and this module should call it.
* **An MzTab or MzTab-M reader.** Neither exists in the source. Every cell type
  already parses itself, so a reader is a section dispatcher over the column
  orders tabulated above.
* **The `AccurateMassSearch` TOPP tool**, the one direct consumer of
  `MzTabMFile.h`. With this package plus a `.featureXML` reader it becomes a
  tier-1 TOPP differential against
  `AccurateMassSearchEngine_output1_mztabm_featureXML.mzTab`, already committed
  here as a fixture.

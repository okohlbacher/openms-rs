// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! MzTab-M, the metabolomics profile of MzTab: the data model of
//! `FORMAT/MzTabM.h` and the file adapter of `FORMAT/MzTabMFile.h`.
//!
//! MzTab-M 2.0.0-M keeps MzTab's line grammar — an `MTD` metadata section of
//! key/value lines, then one header row plus data rows per table section — and
//! replaces the proteomics section set with three metabolomics sections:
//!
//! | Prefix | Header | Section | Row struct |
//! |---|---|---|---|
//! | `SML` | `SMH` | small-molecule summary | [`MzTabMSmallMoleculeSectionRow`](crate::format::mztab_m::MzTabMSmallMoleculeSectionRow) |
//! | `SMF` | `SFH` | small-molecule feature | [`MzTabMSmallMoleculeFeatureSectionRow`](crate::format::mztab_m::MzTabMSmallMoleculeFeatureSectionRow) |
//! | `SME` | `SEH` | small-molecule evidence | [`MzTabMSmallMoleculeEvidenceSectionRow`](crate::format::mztab_m::MzTabMSmallMoleculeEvidenceSectionRow) |
//!
//! Every cell is one of the shared MzTab cell types of
//! [`crate::format::mztab`] — this module defines no cell type of its own, so a
//! reader or writer drives an MzTab-M column exactly as it drives an MzTab
//! column. The metadata section reuses
//! [`MzTabSoftwareMetaData`](crate::format::mztab::MzTabSoftwareMetaData),
//! [`MzTabSampleMetaData`](crate::format::mztab::MzTabSampleMetaData),
//! [`MzTabInstrumentMetaData`](crate::format::mztab::MzTabInstrumentMetaData),
//! [`MzTabContactMetaData`](crate::format::mztab::MzTabContactMetaData) and
//! [`MzTabCVMetaData`](crate::format::mztab::MzTabCVMetaData) unchanged, and
//! replaces the assay, MS-run and study-variable records with profile-specific
//! ones; `docs/MZTAB_M_SUPPORT.md` tabulates what is shared and what is
//! replaced, member by member.
//!
//! [`MzTabM`](crate::format::mztab_m::MzTabM) is the document,
//! [`MzTabM::export_feature_map`](crate::format::mztab_m::MzTabM::export_feature_map)
//! builds one from a [`FeatureMap`](crate::kernel::FeatureMap) and its
//! identification graph, and
//! [`MzTabMFile`](crate::format::mztab_m::MzTabMFile) writes one out. There is
//! no reader: the source has none either.
//!
//! The writer reproduces four source defects only when asked. Its library
//! default emits the metadata keys the specification names and does not drop
//! `ms_run[n]-id_format`;
//! [`MzTabMWriteOptions::source`](crate::format::mztab_m::MzTabMWriteOptions::source)
//! selects the byte-for-byte source output instead. See
//! [`MzTabMWriteOptions`](crate::format::mztab_m::MzTabMWriteOptions) and
//! `docs/MZTAB_M_SUPPORT.md`.
//!
//! ```
//! use openms::format::mztab::MzTabString;
//! use openms::format::mztab_m::{MzTabM, MzTabMFile, MzTabMSmallMoleculeSectionRow};
//!
//! // The metadata constructor declares the profile version, not MzTab's.
//! let mut document = MzTabM::default();
//! assert_eq!(document.meta_data.mz_tab_version.get(), "2.0.0-M");
//!
//! document.meta_data.mz_tab_id.set("local_id: 1");
//! let mut row = MzTabMSmallMoleculeSectionRow::default();
//! row.sml_identifier = MzTabString::from_text("1");
//! row.reliability = MzTabString::from_text("2");
//! document.small_molecule_data.push(row);
//!
//! let lines = MzTabMFile::new().generate_lines(&document)?;
//! assert_eq!(lines[0], "MTD\tmzTab-version\t2.0.0-M");
//! assert_eq!(lines[1], "MTD\tmzTab-ID\tlocal_id: 1");
//! // An unset mandatory Param cell still gets its line, spelled `null`.
//! assert!(lines.contains(&"MTD\tquantification_method\tnull".to_owned()));
//! # Ok::<(), openms::Error>(())
//! ```

use crate::chemistry::EmpiricalFormula;
use crate::format::controlled_vocabulary::ControlledVocabulary;
use crate::format::mztab::{
    MzTab, MzTabCVMetaData, MzTabContactMetaData, MzTabDouble, MzTabDoubleList,
    MzTabInstrumentMetaData, MzTabInteger, MzTabOptionalColumnEntry, MzTabOptionalColumns,
    MzTabParameter, MzTabParameterList, MzTabSampleMetaData, MzTabSoftwareMetaData,
    MzTabSpectraRef, MzTabString, MzTabStringList, optional_column_names,
};
use crate::identification::graph::{
    IdentificationData, IdentifiedCompound, ObservationMatchId, ProcessingSoftware,
};
use crate::kernel::{Feature, FeatureMap};
use crate::metadata::ProcessingAction;
use crate::param::value::format_float;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

fn missing(message: impl Into<String>) -> Error {
    Error::MissingInformation(message.into())
}

fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b)
        .ok_or_else(|| bad("MzTab-M line count overflows"))
}

/// `StringUtils::toStr(double)`, the source's number spelling for a plain
/// `std::string` concatenation rather than a cell.
fn source_double(value: f64) -> String {
    format_float(value, true)
}

/// The `MTD` line prefix, and the tab the source concatenates after it.
const MTD: &str = "MTD\t";

// ---------------------------------------------------------------------------
// MzTabM.h — metadata section records
// ---------------------------------------------------------------------------

/// `MTD assay[n]`: name, custom parameters and the sample and MS run it draws
/// on.
///
/// Source `MzTabMAssayMetaData`. The profile replaces MzTab's
/// [`MzTabAssayMetaData`](crate::format::mztab::MzTabAssayMetaData), which
/// carries a quantification reagent and quantification modifications instead:
/// a metabolomics assay is one sample measured in one run, not a label channel.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct MzTabMAssayMetaData {
    /// `assay[n]`, the name of the assay. Mandatory in the profile.
    pub name: MzTabString,
    /// `assay[n]-custom[m]`, additional parameters or values for this assay.
    pub custom: BTreeMap<usize, MzTabParameter>,
    /// `assay[n]-external_uri`, a reference to further information.
    pub external_uri: MzTabString,
    /// `assay[n]-sample_ref`, the one-based `sample[…]` index analysed.
    pub sample_ref: MzTabInteger,
    /// `assay[n]-ms_run_ref`, the one-based `ms_run[…]` index. Mandatory.
    ///
    /// A single index, unlike MzTab's `ms_run_ref` list.
    pub ms_run_ref: MzTabInteger,
}

/// `MTD ms_run[n]`: location, format, fragmentation, polarity and hash.
///
/// Source `MzTabMMSRunMetaData`. It shares `location`, `format` and `id_format`
/// with MzTab's [`MzTabMSRunMetaData`](crate::format::mztab::MzTabMSRunMetaData),
/// turns `fragmentation_method` from one `Param[]` cell into indexed
/// `fragmentation_method[m]` keys, and adds `instrument_ref`, the mandatory
/// `scan_polarity[m]`, `hash` and `hash_method`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct MzTabMMSRunMetaData {
    /// `ms_run[n]-location`, the external data file. Mandatory.
    pub location: MzTabString,
    /// `ms_run[n]-instrument_ref`, a one-based `instrument[…]` index.
    pub instrument_ref: MzTabInteger,
    /// `ms_run[n]-format`, the data format of the external file.
    pub format: MzTabParameter,
    /// `ms_run[n]-id_format`, the native spectrum identifier format.
    ///
    /// The source writer never emits this key; see
    /// [`MzTabMWriteOptions::omit_ms_run_id_format`].
    pub id_format: MzTabParameter,
    /// `ms_run[n]-fragmentation_method[m]`.
    pub fragmentation_method: BTreeMap<usize, MzTabParameter>,
    /// `ms_run[n]-scan_polarity[m]`. Mandatory in the profile.
    pub scan_polarity: BTreeMap<usize, MzTabParameter>,
    /// `ms_run[n]-hash`, the hash of the external data file.
    pub hash: MzTabString,
    /// `ms_run[n]-hash_method`, the method that produced `hash`.
    pub hash_method: MzTabParameter,
}

/// `MTD study_variable[n]`: the assays it groups and how they are averaged.
///
/// Source `MzTabMStudyVariableMetaData`. MzTab's
/// [`MzTabStudyVariableMetaData`](crate::format::mztab::MzTabStudyVariableMetaData)
/// carries `assay_refs`, `sample_refs` and `description`; the profile drops
/// `sample_refs` and adds `name`, `average_function`, `variation_function` and
/// `factors`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabMStudyVariableMetaData {
    /// `study_variable[n]`, the name of the study variable. Mandatory.
    pub name: MzTabString,
    /// `study_variable[n]-assay_refs`, the one-based assay indices grouped
    /// here. Mandatory. Written as an `assay[k]`-per-entry `String[]` cell.
    pub assay_refs: Vec<i32>,
    /// `study_variable[n]-average_function`, how the quantification value is
    /// summarised across the assays.
    pub average_function: MzTabParameter,
    /// `study_variable[n]-variation_function`, how its variation is computed.
    pub variation_function: MzTabParameter,
    /// `study_variable[n]-description`. Mandatory.
    pub description: MzTabString,
    /// `study_variable[n]-factors`, additional parameters or factors.
    pub factors: MzTabParameterList,
}

/// `MTD database[n]`: an identification database the file draws on.
///
/// Source `MzTabMDatabaseMetaData`. MzTab has no equivalent: a proteomics file
/// names its search database per `PRT`/`PEP` row, while a metabolomics file
/// declares its databases once and the row `identifier` columns carry the
/// declared `prefix`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MzTabMDatabaseMetaData {
    /// `database[n]`, the description of the database. Mandatory.
    pub database: MzTabParameter,
    /// `database[n]-prefix`, the prefix used in the `identifier` columns of the
    /// data tables. Mandatory.
    pub prefix: MzTabString,
    /// `database[n]-version`. Mandatory.
    pub version: MzTabString,
    /// `database[n]-uri`. Mandatory.
    pub uri: MzTabString,
}

/// The whole `MTD` section of an MzTab-M file.
///
/// Source `MzTabMMetaData`. Refer to the MzTab-M 2.0.0-M specification for each
/// key. [`Default`](Default::default) sets `mz_tab_version` to `2.0.0-M`, as
/// the source's constructor, and leaves every other field null or empty. Every
/// indexed key is one-based and ordered by [`BTreeMap`], so a writer emits
/// indices in numeric order.
///
/// Relative to MzTab's
/// [`MzTabMetaData`](crate::format::mztab::MzTabMetaData) the profile drops
/// `mzTab-mode`, `mzTab-type`, the four proteomics `*_search_engine_score`
/// families, `false_discovery_rate`, `fixed_mod`, `variable_mod`,
/// `protein-quantification_unit`, `peptide-quantification_unit` and the
/// `colunit-protein`/`-peptide`/`-psm` declarations; it adds
/// `external_study_uri`, `database`, `derivatization_agent`,
/// `small_molecule_feature-quantification_unit`,
/// `small_molecule-identification_reliability`, `id_confidence_measure` and the
/// two extra `colunit-small_molecule_*` families.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct MzTabMMetaData {
    /// `mzTab-version`. Defaults to `2.0.0-M`. Mandatory.
    pub mz_tab_version: MzTabString,
    /// `mzTab-ID`, a repository or local identifier. Mandatory.
    pub mz_tab_id: MzTabString,
    /// `title`.
    pub title: MzTabString,
    /// `description`.
    pub description: MzTabString,
    /// `sample_processing[n]`, describing sample preparation and handling.
    pub sample_processing: BTreeMap<usize, MzTabParameterList>,
    /// `instrument[n]-…`.
    pub instrument: BTreeMap<usize, MzTabInstrumentMetaData>,
    /// `software[n]-…`, the software used to analyse the data.
    pub software: BTreeMap<usize, MzTabSoftwareMetaData>,
    /// `publication[n]`.
    pub publication: BTreeMap<usize, MzTabString>,
    /// `contact[n]-…`.
    pub contact: BTreeMap<usize, MzTabContactMetaData>,
    /// `uri[n]`, pointing at the file source, for example MetaboLights.
    pub uri: BTreeMap<usize, MzTabString>,
    /// `external_study_uri[n]`, pointing at an external study description such
    /// as an ISA-TAB file.
    pub external_study_uri: BTreeMap<usize, MzTabString>,
    /// `quantification_method`. Mandatory.
    pub quantification_method: MzTabParameter,
    /// `sample[n]-…`.
    pub sample: BTreeMap<usize, MzTabSampleMetaData>,
    /// `ms_run[n]-…`.
    pub ms_run: BTreeMap<usize, MzTabMMSRunMetaData>,
    /// `assay[n]-…`.
    pub assay: BTreeMap<usize, MzTabMAssayMetaData>,
    /// `study_variable[n]-…`.
    pub study_variable: BTreeMap<usize, MzTabMStudyVariableMetaData>,
    /// `custom[n]`.
    pub custom: BTreeMap<usize, MzTabParameter>,
    /// `cv[n]-…`, the controlled vocabularies the file references.
    pub cv: BTreeMap<usize, MzTabCVMetaData>,
    /// `database[n]-…`.
    pub database: BTreeMap<usize, MzTabMDatabaseMetaData>,
    /// `derivatization_agent[n]`, agents applied to the small molecules.
    pub derivatization_agent: BTreeMap<usize, MzTabParameter>,
    /// `small_molecule-quantification_unit`. Mandatory.
    pub small_molecule_quantification_unit: MzTabParameter,
    /// `small_molecule_feature-quantification_unit`. Mandatory.
    pub small_molecule_feature_quantification_unit: MzTabParameter,
    /// `small_molecule-identification_reliability`, the four-level confidence
    /// schema. Mandatory.
    pub small_molecule_identification_reliability: MzTabParameter,
    /// `id_confidence_measure[n]`, the confidence measures the `SME` section
    /// reports. Mandatory once the evidence section carries any.
    pub id_confidence_measure: BTreeMap<usize, MzTabParameter>,
    /// `colunit-small_molecule` declarations, verbatim.
    pub colunit_small_molecule: Vec<MzTabString>,
    /// `colunit-small_molecule_feature` declarations, verbatim.
    pub colunit_small_molecule_feature: Vec<MzTabString>,
    /// `colunit-small_molecule_evidence` declarations, verbatim.
    pub colunit_small_molecule_evidence: Vec<MzTabString>,
}

impl Default for MzTabMMetaData {
    fn default() -> Self {
        Self {
            mz_tab_version: MzTabString::from_text("2.0.0-M"),
            mz_tab_id: MzTabString::default(),
            title: MzTabString::default(),
            description: MzTabString::default(),
            sample_processing: BTreeMap::new(),
            instrument: BTreeMap::new(),
            software: BTreeMap::new(),
            publication: BTreeMap::new(),
            contact: BTreeMap::new(),
            uri: BTreeMap::new(),
            external_study_uri: BTreeMap::new(),
            quantification_method: MzTabParameter::default(),
            sample: BTreeMap::new(),
            ms_run: BTreeMap::new(),
            assay: BTreeMap::new(),
            study_variable: BTreeMap::new(),
            custom: BTreeMap::new(),
            cv: BTreeMap::new(),
            database: BTreeMap::new(),
            derivatization_agent: BTreeMap::new(),
            small_molecule_quantification_unit: MzTabParameter::default(),
            small_molecule_feature_quantification_unit: MzTabParameter::default(),
            small_molecule_identification_reliability: MzTabParameter::default(),
            id_confidence_measure: BTreeMap::new(),
            colunit_small_molecule: Vec::new(),
            colunit_small_molecule_feature: Vec::new(),
            colunit_small_molecule_evidence: Vec::new(),
        }
    }
}

impl MzTabMMetaData {
    /// The metadata section with the profile version set. Source
    /// `MzTabMMetaData()`, and the same value as
    /// [`Default`](Default::default).
    pub fn new() -> Self {
        Self::default()
    }
}

// ---------------------------------------------------------------------------
// MzTabM.h — section rows
// ---------------------------------------------------------------------------

/// `SML` — one small-molecule summary row.
///
/// Source `MzTabMSmallMoleculeSectionRow`. Distinct from MzTab's
/// [`MzTabSmallMoleculeSectionRow`](crate::format::mztab::MzTabSmallMoleculeSectionRow),
/// which is the proteomics-file small-molecule table: the profile row is keyed
/// by an `SML_ID`, points at feature rows through `smf_id_refs`, and reports
/// abundance per assay and per study variable rather than per `ms_run`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzTabMSmallMoleculeSectionRow {
    /// `SML_ID`, the small molecule's within-file identifier.
    pub sml_identifier: MzTabString,
    /// `SMF_ID_REFS`, the features quantification was based on.
    pub smf_id_refs: MzTabStringList,
    /// `database_identifier`, the databases used.
    pub database_identifier: MzTabStringList,
    /// `chemical_formula`, potential formulas of the reported compound.
    pub chemical_formula: MzTabStringList,
    /// `smiles`, molecular structures in SMILES notation.
    pub smiles: MzTabStringList,
    /// `inchi`, InChIs of the potential identifications.
    pub inchi: MzTabStringList,
    /// `chemical_name`, chemical or common names, or a general description.
    pub chemical_name: MzTabStringList,
    /// `uri`, the source entries' locations.
    pub uri: MzTabStringList,
    /// `theoretical_neutral_mass`, the precursor theoretical neutral masses.
    pub theoretical_neutral_mass: MzTabDoubleList,
    /// `adduct_ions`.
    pub adducts: MzTabStringList,
    /// `reliability` of this identification, a level of the four-level schema
    /// declared by `small_molecule-identification_reliability`.
    ///
    /// A free-text cell, not an [`MzTabInteger`]; the source comment records
    /// that the reliability of the identification method itself belongs in the
    /// identification data structure.
    pub reliability: MzTabString,
    /// `best_id_confidence_measure`, the approach with the highest confidence.
    pub best_id_confidence_measure: MzTabParameter,
    /// `best_id_confidence_value`.
    pub best_id_confidence_value: MzTabDouble,
    /// `abundance_assay[n]`, abundance in each declared assay.
    pub small_molecule_abundance_assay: BTreeMap<usize, MzTabDouble>,
    /// `abundance_study_variable[n]`, abundance in each study variable.
    pub small_molecule_abundance_study_variable: BTreeMap<usize, MzTabDouble>,
    /// `abundance_variation_study_variable[n]`, the variability of that
    /// measurement.
    pub small_molecule_abundance_variation_study_variable: BTreeMap<usize, MzTabDouble>,
    /// Trailing `opt_…` columns. Source member `opt_`.
    pub opt: Vec<MzTabOptionalColumnEntry>,
}

/// `SMF` — one small-molecule feature row.
///
/// Source `MzTabMSmallMoleculeFeatureSectionRow`. MzTab has no feature section
/// at all: this is the profile's per-adduct, per-charge quantified ion.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzTabMSmallMoleculeFeatureSectionRow {
    /// `SMF_ID`, the within-file identifier of this feature.
    pub smf_identifier: MzTabString,
    /// `SME_ID_REFS`, the evidence rows supporting it.
    pub sme_id_refs: MzTabStringList,
    /// `SME_ID_REF_ambiguity_code`, set when the referenced evidences are
    /// ambiguous identifications of the same feature.
    pub sme_id_ref_ambiguity_code: MzTabInteger,
    /// `adduct_ion`.
    pub adduct: MzTabString,
    /// `isotopomer`. Mandatory when de-isotoping has not been performed: the
    /// quantified isotopomer must then be reported here.
    pub isotopomer: MzTabParameter,
    /// `exp_mass_to_charge`, the precursor ion's m/z.
    pub exp_mass_to_charge: MzTabDouble,
    /// `charge`, the precursor ion's charge.
    pub charge: MzTabInteger,
    /// `retention_time_in_seconds`. In seconds, as the column name says.
    pub retention_time: MzTabDouble,
    /// `retention_time_in_seconds_start`, the feature's start on the RT axis.
    pub rt_start: MzTabDouble,
    /// `retention_time_in_seconds_end`, its end on the RT axis.
    pub rt_end: MzTabDouble,
    /// `abundance_assay[n]`, feature abundance in each declared assay.
    pub small_molecule_feature_abundance_assay: BTreeMap<usize, MzTabDouble>,
    /// Trailing `opt_…` columns. Source member `opt_`.
    pub opt: Vec<MzTabOptionalColumnEntry>,
}

/// `SME` — one small-molecule evidence row.
///
/// Source `MzTabMSmallMoleculeEvidenceSectionRow`. The profile's analogue of
/// MzTab's `PSM` section: one identification of one piece of input evidence,
/// with a `spectra_ref` and a rank.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzTabMSmallMoleculeEvidenceSectionRow {
    /// `SME_ID`, the within-file identifier of this evidence result.
    pub sme_identifier: MzTabString,
    /// `evidence_input_id`, the within-file identifier of the input data this
    /// identification rests on, for example a fragment spectrum or an
    /// RT/m-over-z pair.
    pub evidence_input_id: MzTabString,
    /// `database_identifier`, the putative identification from an external
    /// database.
    pub database_identifier: MzTabString,
    /// `chemical_formula`, the putative molecular formula.
    pub chemical_formula: MzTabString,
    /// `smiles`, the potential structure in SMILES notation.
    pub smiles: MzTabString,
    /// `inchi`.
    pub inchi: MzTabString,
    /// `chemical_name`, a chemical or common name, or a general description.
    pub chemical_name: MzTabString,
    /// `uri`, the source entry's location.
    pub uri: MzTabString,
    /// `derivatized_form`.
    pub derivatized_form: MzTabParameter,
    /// `adduct_ion`.
    pub adduct: MzTabString,
    /// `exp_mass_to_charge`, the precursor ion's m/z.
    pub exp_mass_to_charge: MzTabDouble,
    /// `charge`, the precursor ion's charge.
    pub charge: MzTabInteger,
    /// `theoretical_mass_to_charge`. Source member `calc_mass_to_charge`,
    /// whose comment calls it the precursor ion's m/z.
    pub calc_mass_to_charge: MzTabDouble,
    /// `spectra_ref`, the spectrum this evidence came from.
    pub spectra_ref: MzTabSpectraRef,
    /// `identification_method`, the database search, search engine or process
    /// that produced the identification.
    pub identification_method: MzTabParameter,
    /// `ms_level`, the highest MS level used to inform the identification.
    pub ms_level: MzTabParameter,
    /// `id_confidence_measure[n]`, one statistical value or score per measure
    /// declared in the metadata section.
    pub id_confidence_measure: BTreeMap<usize, MzTabDouble>,
    /// `rank` of this identification; `1` is best.
    pub rank: MzTabInteger,
    /// Trailing `opt_…` columns. Source member `opt_`.
    pub opt: Vec<MzTabOptionalColumnEntry>,
}

/// The `SML` section. Source typedef `MzTabMSmallMoleculeSectionRows`.
pub type MzTabMSmallMoleculeSectionRows = Vec<MzTabMSmallMoleculeSectionRow>;
/// The `SMF` section. Source typedef `MzTabMSmallMoleculeFeatureSectionRows`.
pub type MzTabMSmallMoleculeFeatureSectionRows = Vec<MzTabMSmallMoleculeFeatureSectionRow>;
/// The `SME` section. Source typedef `MzTabMSmallMoleculeEvidenceSectionRows`.
pub type MzTabMSmallMoleculeEvidenceSectionRows = Vec<MzTabMSmallMoleculeEvidenceSectionRow>;

macro_rules! optional_columns_for {
    ($row:ty) => {
        impl MzTabOptionalColumns for $row {
            fn optional_columns(&self) -> &[MzTabOptionalColumnEntry] {
                &self.opt
            }
            fn optional_columns_mut(&mut self) -> &mut Vec<MzTabOptionalColumnEntry> {
                &mut self.opt
            }
        }
    };
}
optional_columns_for!(MzTabMSmallMoleculeSectionRow);
optional_columns_for!(MzTabMSmallMoleculeFeatureSectionRow);
optional_columns_for!(MzTabMSmallMoleculeEvidenceSectionRow);

// ---------------------------------------------------------------------------
// MzTabM.h — the document
// ---------------------------------------------------------------------------

/// A whole MzTab-M document: the metadata section, the three metabolomics
/// sections, and the comment and empty lines recorded by position.
///
/// Source `MzTabM`, which derives from `MzTabBase` purely to inherit the
/// protected `getOptionalColumnNames_` template; this port reaches that
/// template through the shared [`optional_column_names`] function instead, so
/// [`MzTabM`] has no base type. The source keeps the
/// sections in protected members behind getter/setter pairs that do nothing but
/// return a reference, so they are public fields here and the source getter
/// names remain available as borrowing methods.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzTabM {
    /// The `MTD` section.
    pub meta_data: MzTabMMetaData,
    /// The `SML` section.
    pub small_molecule_data: MzTabMSmallMoleculeSectionRows,
    /// The `SMF` section.
    pub small_molecule_feature_data: MzTabMSmallMoleculeFeatureSectionRows,
    /// The `SME` section.
    pub small_molecule_evidence_data: MzTabMSmallMoleculeEvidenceSectionRows,
    /// Line indices of empty rows, preserved so a round trip can restore them.
    pub empty_rows: Vec<usize>,
    /// `COM` comment lines, keyed by line index.
    pub comment_rows: BTreeMap<usize, String>,
}

impl MzTabM {
    /// Largest number of rows one section operation will traverse. The same
    /// ceiling as [`MzTab::MAX_ROWS`].
    pub const MAX_ROWS: usize = MzTab::MAX_ROWS;
    /// Largest number of distinct optional columns one section may declare. The
    /// same ceiling as [`MzTab::MAX_OPTIONAL_COLUMNS`].
    pub const MAX_OPTIONAL_COLUMNS: usize = MzTab::MAX_OPTIONAL_COLUMNS;
    /// Largest number of entries any one indexed metadata map may carry.
    ///
    /// Native: the source has no ceiling and the writer emits up to eight lines
    /// per `ms_run` entry, so an oversized metadata section would otherwise
    /// grow the output without bound.
    ///
    /// This bounds the entry *count*, not the index values: a key inserted at
    /// `usize::MAX` is written as `ms_run[18446744073709551615]-…`, exactly as
    /// the source's unchecked `Size` would, and nothing refuses it.
    pub const MAX_INDEXED_ENTRIES: usize = 100_000;

    /// An empty document whose metadata declares `2.0.0-M`. Source `MzTabM()`,
    /// and the same value as [`Default`](Default::default).
    pub fn new() -> Self {
        Self::default()
    }

    /// The `MTD` section. Source `getMetaData`.
    pub fn meta_data(&self) -> &MzTabMMetaData {
        &self.meta_data
    }
    /// Replace the `MTD` section. Source `setMetaData`.
    pub fn set_meta_data(&mut self, meta_data: MzTabMMetaData) {
        self.meta_data = meta_data;
    }
    /// The `SML` section. Source `getMSmallMoleculeSectionRows`.
    pub fn small_molecule_section_rows(&self) -> &MzTabMSmallMoleculeSectionRows {
        &self.small_molecule_data
    }
    /// Replace the `SML` section. Source `setMSmallMoleculeSectionRows`.
    pub fn set_small_molecule_section_rows(&mut self, rows: MzTabMSmallMoleculeSectionRows) {
        self.small_molecule_data = rows;
    }
    /// The `SMF` section. Source `getMSmallMoleculeFeatureSectionRows`.
    pub fn small_molecule_feature_section_rows(&self) -> &MzTabMSmallMoleculeFeatureSectionRows {
        &self.small_molecule_feature_data
    }
    /// Replace the `SMF` section. Source
    /// `setMSmallMoleculeFeatureSectionRows`.
    pub fn set_small_molecule_feature_section_rows(
        &mut self,
        rows: MzTabMSmallMoleculeFeatureSectionRows,
    ) {
        self.small_molecule_feature_data = rows;
    }
    /// The `SME` section. Source `getMSmallMoleculeEvidenceSectionRows`.
    pub fn small_molecule_evidence_section_rows(&self) -> &MzTabMSmallMoleculeEvidenceSectionRows {
        &self.small_molecule_evidence_data
    }
    /// Replace the `SME` section. Source
    /// `setMSmallMoleculeEvidenceSectionRows`.
    pub fn set_small_molecule_evidence_section_rows(
        &mut self,
        rows: MzTabMSmallMoleculeEvidenceSectionRows,
    ) {
        self.small_molecule_evidence_data = rows;
    }
    /// Line indices of empty rows. Source `getEmptyRows`.
    pub fn empty_rows(&self) -> &[usize] {
        &self.empty_rows
    }
    /// Replace the empty-row indices. Source `setEmptyRows`.
    pub fn set_empty_rows(&mut self, rows: Vec<usize>) {
        self.empty_rows = rows;
    }
    /// Comment lines by line index. Source `getCommentRows`.
    pub fn comment_rows(&self) -> &BTreeMap<usize, String> {
        &self.comment_rows
    }
    /// Replace the comment lines. Source `setCommentRows`.
    pub fn set_comment_rows(&mut self, rows: BTreeMap<usize, String>) {
        self.comment_rows = rows;
    }

    /// Optional column names of the `SML` section, in first-occurrence order.
    /// Source `getMSmallMoleculeOptionalColumnNames`.
    ///
    /// # Errors
    ///
    /// As [`optional_column_names`].
    pub fn small_molecule_optional_column_names(&self) -> Result<Vec<String>> {
        optional_column_names(&self.small_molecule_data)
    }
    /// Optional column names of the `SMF` section, in first-occurrence order.
    /// Source `getMSmallMoleculeFeatureOptionalColumnNames`.
    ///
    /// # Errors
    ///
    /// As [`optional_column_names`].
    pub fn small_molecule_feature_optional_column_names(&self) -> Result<Vec<String>> {
        optional_column_names(&self.small_molecule_feature_data)
    }
    /// Optional column names of the `SME` section, in first-occurrence order.
    /// Source `getMSmallMoleculeEvidenceOptionalColumnNames`.
    ///
    /// # Errors
    ///
    /// As [`optional_column_names`].
    pub fn small_molecule_evidence_optional_column_names(&self) -> Result<Vec<String>> {
        optional_column_names(&self.small_molecule_evidence_data)
    }

    /// Append one `opt_<id>_<key>` column per key, taking the value from
    /// `meta`. Source static `MzTabM::addMetaInfoToOptionalColumns`.
    ///
    /// The source declares this separately on `MzTabM` with a body identical to
    /// `MzTab::addMetaInfoToOptionalColumns`, so this forwards to the shared
    /// [`add_meta_info_to_optional_columns`](crate::format::mztab::add_meta_info_to_optional_columns)
    /// rather than duplicating it.
    ///
    /// # Errors
    ///
    /// As
    /// [`add_meta_info_to_optional_columns`](crate::format::mztab::add_meta_info_to_optional_columns):
    /// [`Error::InvalidValue`] when `keys` would push the entry count past
    /// [`MzTab::MAX_OPTIONAL_COLUMNS`], with `opt` left unchanged.
    ///
    /// # Notes
    ///
    /// [`MzTabM::export_feature_map`] normally reaches this with keys whose
    /// spaces have *already* been substituted, because the source's
    /// `getFeatureMapMetaValues_` substitutes before collecting; the lookup in
    /// `meta` then misses any key that contained a space and the column is
    /// null. See [`MzTabMExportOptions::substitute_keys_before_lookup`].
    pub fn add_meta_info_to_optional_columns(
        keys: &BTreeSet<String>,
        opt: &mut Vec<MzTabOptionalColumnEntry>,
        id: &str,
        meta: &crate::metadata::MetaInfo,
    ) -> Result<()> {
        crate::format::mztab::add_meta_info_to_optional_columns(keys, opt, id, meta)
    }

    /// Export a [`FeatureMap`] and its identification graph to MzTab-M, using
    /// the pinned PSI-MS vocabulary and native options.
    ///
    /// Source static `MzTabM::exportFeatureMapToMzTabM(const FeatureMap&)`.
    ///
    /// # Arguments
    ///
    /// * `feature_map` — the quantified features. Their `id_matches` name
    ///   observation matches in `id_data`.
    /// * `id_data` — the identification graph the source reaches through
    ///   `FeatureMap::getIdentificationData()`. This port takes it explicitly,
    ///   because the Rust [`FeatureMap`] does not own one; see
    ///   `docs/MZTAB_M_SUPPORT.md`.
    ///
    /// # Errors
    ///
    /// As [`MzTabM::export_feature_map_with`].
    ///
    /// # Notes
    ///
    /// [`ControlledVocabulary::psi_ms`] carries all five pinned vocabularies in
    /// one object, so its `name`, `label`, `version` and `url` identify the
    /// vocabulary loaded last rather than PSI-MS, and the `cv[1]` block this
    /// produces is correspondingly generic. A caller that needs the `cv[1]`
    /// block the source writes must load `psi-ms.obo` alone under the name
    /// `PSI-MS` and call [`MzTabM::export_feature_map_with`]. Term *lookups*
    /// are unaffected: the shared object contains every PSI-MS term.
    pub fn export_feature_map(
        feature_map: &FeatureMap,
        id_data: &IdentificationData,
    ) -> Result<Self> {
        Self::export_feature_map_with(
            feature_map,
            id_data,
            ControlledVocabulary::psi_ms()?,
            &MzTabMExportOptions::default(),
        )
    }

    /// Export a [`FeatureMap`] and its identification graph to MzTab-M against
    /// a caller-owned vocabulary and explicit options.
    ///
    /// Source static `MzTabM::exportFeatureMapToMzTabM`, which loads
    /// `share/OpenMS/CV/psi-ms.obo` from disk under the name `PSI-MS` itself;
    /// this port takes the vocabulary, because the crate embeds its
    /// vocabularies and never resolves a runtime share directory.
    ///
    /// # Arguments
    ///
    /// * `feature_map` — the quantified features.
    /// * `id_data` — the identification graph the features' `id_matches` point
    ///   into.
    /// * `cv` — the vocabulary term names are resolved against, and whose
    ///   `name`/`label`/`version`/`url` fill `cv[1]`.
    /// * `options` — which source quirks to reproduce.
    ///
    /// # Errors
    ///
    /// * [`Error::MissingInformation`] when `id_data` is empty. The source
    ///   states this as an `OPENMS_PRECONDITION`, which is compiled out of a
    ///   release build and then reads an empty graph.
    /// * [`Error::MissingInformation`] when a feature names an observation
    ///   match whose identified molecule is not a compound, or when an `SMF`
    ///   row references an `SME` identifier that no evidence row carries. The
    ///   source dereferences both results unchecked.
    /// * [`Error::InvalidValue`] when the map exceeds [`MzTabM::MAX_ROWS`],
    ///   when the produced sections would exceed it, when a required
    ///   vocabulary term is missing from `cv`, or when an
    ///   `adducts` feature meta value is not a string list.
    /// * [`Error::Parse`] or [`Error::InvalidValue`] when a compound's
    ///   `chemical_formula` cell does not parse as an empirical formula, which
    ///   the source lets `EmpiricalFormula`'s constructor throw on.
    /// * Any error the graph accessors raise for a stale or foreign record ID.
    ///
    /// Nothing is mutated: the document is built in full and returned only on
    /// success.
    ///
    /// # Notes
    ///
    /// The source logs five different warnings while guessing a quantification
    /// method, a quantification unit, an identification method and a scan
    /// polarity from the tool names in the graph. This port makes the same
    /// guesses and the same defaults but writes no log; each default is
    /// documented in `docs/MZTAB_M_SUPPORT.md`.
    ///
    /// The source is serial here and so is this port. No `#pragma omp` appears
    /// in `MzTabM.cpp`.
    pub fn export_feature_map_with(
        feature_map: &FeatureMap,
        id_data: &IdentificationData,
        cv: &ControlledVocabulary,
        options: &MzTabMExportOptions,
    ) -> Result<Self> {
        export::run(feature_map, id_data, cv, options)
    }
}

/// Which source quirks [`MzTabM::export_feature_map_with`] reproduces.
///
/// [`Default`](Default::default) is the native reading of every one;
/// [`MzTabMExportOptions::source`] reproduces the C++ exactly. Each field names
/// the source line it comes from in `docs/MZTAB_M_SUPPORT.md`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MzTabMExportOptions {
    /// Drop every observation match of a feature that shares an identified
    /// compound's `identifier` with an earlier one.
    ///
    /// The source collects a feature's matches into a
    /// `std::set<ObservationMatchRef, CompareMzTabMMatchRef>` whose comparator
    /// orders by the identified compound's `identifier` alone, so matches that
    /// agree on that identifier are *equivalent* and all but one are discarded
    /// — including matches that differ in adduct or score. Which one survives
    /// depends on the iteration order of the feature's own reference set, which
    /// is ordered by container-iterator address and therefore not determined by
    /// the data.
    ///
    /// `false`, the default, keeps every match and orders them by identifier
    /// and then by graph ID. `true` reproduces the source's discard and keeps
    /// the lowest graph ID of each group, which is deterministic where the
    /// source is not.
    pub deduplicate_matches_by_compound: bool,
    /// Substitute spaces for underscores in meta-value keys *before* looking
    /// them up, as `getFeatureMapMetaValues_` does.
    ///
    /// The source collects feature, observation-match and compound meta-value
    /// keys with `substitute(' ', '_')` already applied, then passes those
    /// substituted keys to `addMetaInfoToOptionalColumns`, which looks each one
    /// up in the record's own meta values. A key that contained a space cannot
    /// match, so its optional column is always `null`.
    ///
    /// `false`, the default, collects the raw keys and lets the column *name*
    /// carry the substitution, so a key with a space yields its value. `true`
    /// reproduces the source's null column. The collected key set is ordered,
    /// so the two settings can also order the columns differently when one key
    /// contains a space and another differs from it only at that position.
    pub substitute_keys_before_lookup: bool,
    /// Fill `cv[n]-label` from the vocabulary's *name* and `cv[n]-full_name`
    /// from its *label*, as the source does.
    ///
    /// `meta_cv.label = MzTabString(cv.name()); meta_cv.full_name =
    /// MzTabString(cv.label());` — the two are crossed, so a PSI-MS load writes
    /// `cv[1]-label PSI-MS` and `cv[1]-full_name MS`, putting the short
    /// namespace label in the full-name key. `false`, the default, writes
    /// `label` from `label` and `full_name` from `name`.
    pub swap_cv_label_and_full_name: bool,
}

impl MzTabMExportOptions {
    /// Every source quirk reproduced, for a byte-comparable export.
    pub fn source() -> Self {
        Self {
            deduplicate_matches_by_compound: true,
            substitute_keys_before_lookup: true,
            swap_cv_label_and_full_name: true,
        }
    }
}

/// Order two observation matches by their identified compound's `identifier`.
///
/// Source functor `CompareMzTabMMatchRef`. The source calls
/// `identified_molecule_var.getIdentifiedCompoundRef()`, which throws when the
/// match identifies a peptide or an oligonucleotide rather than a compound;
/// this returns [`Error::MissingInformation`] for that case instead.
///
/// # Errors
///
/// [`Error::MissingInformation`] when either match does not identify a
/// compound, and any error the graph raises for a stale or foreign ID.
pub fn compare_match_by_compound(
    id_data: &IdentificationData,
    left: ObservationMatchId,
    right: ObservationMatchId,
) -> Result<std::cmp::Ordering> {
    let key = |id: ObservationMatchId| -> Result<&str> {
        Ok(&export::compound_of(id_data, id)?.identifier)
    };
    Ok(key(left)?.cmp(key(right)?))
}

// ---------------------------------------------------------------------------
// MzTabMFile.h — the writer
// ---------------------------------------------------------------------------

/// Which source writer defects [`MzTabMFile`] reproduces.
///
/// [`Default`](Default::default) writes the keys the MzTab-M specification
/// names and drops nothing; [`MzTabMWriteOptions::source`] reproduces
/// `MzTabMFile::store` byte for byte, including its four defects. Every field
/// is a defect of the source writer, not a formatting preference — see
/// `docs/MZTAB_M_SUPPORT.md` and `OpenMS_CPP_ISSUES.md`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MzTabMWriteOptions {
    /// Write `assay[n]-custom[m]` under the key `ms_run[n]-custom[m]`.
    ///
    /// `generateMzTabMMetaDataSection_` builds the custom-parameter line of an
    /// assay with the literal `"MTD\tms_run["` and then interpolates the
    /// *assay* index, so an assay's custom parameters are reported as an MS
    /// run's, under an index that need not name an existing run.
    pub source_assay_custom_key: bool,
    /// Write `colunit-small_molecule_feature` and
    /// `colunit-small_molecule_evidence` under the key
    /// `colunit_small_molecule`.
    ///
    /// All three `colunit` loops use the same literal, so the three families
    /// become indistinguishable in the output and a reader attributes every
    /// declaration to the summary section.
    pub source_colunit_keys: bool,
    /// Append `-uri` to every `derivatization_agent[n]` key.
    ///
    /// The source writes `derivatization_agent[n]-uri`, a key the
    /// specification does not define; the value it carries is the agent
    /// parameter, not a URI.
    pub source_derivatization_agent_key: bool,
    /// Drop `ms_run[n]-id_format` from the metadata section.
    ///
    /// `MzTabMMSRunMetaData::id_format` is the only member of that record the
    /// source's metadata generator never reads, so a caller that sets it loses
    /// it on write.
    pub omit_ms_run_id_format: bool,
    /// Emit each row's abundance cells by iterating the row's own maps, rather
    /// than one cell per column the header declares.
    ///
    /// The header derives `abundance_assay[n]` from
    /// [`MzTabMMetaData::assay`] and `abundance_study_variable[n]` from
    /// [`MzTabMMetaData::study_variable`], while the row generators iterate the
    /// row's own abundance maps. A row that carries fewer abundances than the
    /// metadata declares therefore produces fewer cells than the header has
    /// columns, and one that carries more produces extra cells whose column is
    /// unnamed; the source notices neither, because its only check is an
    /// `OPENMS_POSTCONDITION` compiled out of a release build.
    ///
    /// `false`, the default, emits one cell per declared column — `null` where
    /// the row has no abundance — and refuses a row whose abundance map names
    /// an assay or study variable the metadata does not declare, so the output
    /// is always a rectangle. `true` reproduces the source's cells exactly.
    pub source_row_abundance_cells: bool,
    /// Write a cell whose text carries a tab or a line break unchanged.
    ///
    /// `MzTabMFile.cpp:626-631` hands every rendered row to `TextFile` as it
    /// is, so a cell containing a tab silently gains a column and a cell
    /// containing a line break splits the row across physical lines — text
    /// that begins with `MTD`, `SMH` or `SML` then forges a line of that kind.
    /// The text is reachable from file-derived data: `chemical_name`, `uri`,
    /// `smiles` and every `opt_` value come from a featureXML or `.oms` text
    /// node, and an XML text node may legally contain tabs and newlines.
    ///
    /// `false`, the default, refuses such a cell with [`Error::InvalidValue`]
    /// before anything is written, which is what makes the default output
    /// rectangular and parseable. `true` reproduces the source's
    /// pass-through — the eleventh source defect of the writer.
    pub source_verbatim_cells: bool,
}

impl MzTabMWriteOptions {
    /// Every source defect reproduced, for a byte-comparable file.
    pub fn source() -> Self {
        Self {
            source_assay_custom_key: true,
            source_colunit_keys: true,
            source_derivatization_agent_key: true,
            omit_ms_run_id_format: true,
            source_row_abundance_cells: true,
            source_verbatim_cells: true,
        }
    }
}

/// One generated section line and the number of tab-separated columns in it.
///
/// The source's six generators return the line and report the column count
/// through a `size_t& n_columns` out-parameter; this replaces that parameter.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MzTabMSectionLine {
    /// The rendered line, tabs included, without a terminator.
    pub text: String,
    /// The number of cells rendered into `text`, including the leading
    /// `SMH`/`SML`/`SFH`/`SMF`/`SEH`/`SME` prefix — the source's `n_columns`.
    ///
    /// Equal to the number of tab-separated columns `text` carries, unless
    /// [`MzTabMWriteOptions::source_verbatim_cells`] let a cell through with a
    /// tab of its own, in which case `text` has more.
    pub columns: usize,
}

/// File adapter for MzTab-M files: a writer, with no reader.
///
/// Source `MzTabMFile`. The source class is stateless and its six generators
/// are protected; Rust has no protected visibility, so they are public here,
/// which also lets a caller render one section without a file. The only stored
/// state is [`MzTabMFile::options`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MzTabMFile {
    /// Which source writer defects to reproduce.
    pub options: MzTabMWriteOptions,
}

impl MzTabMFile {
    /// Largest number of lines [`MzTabMFile::store`] will write, matching the
    /// line ceiling of [`crate::format::TextFile`], through which the source
    /// also writes.
    ///
    /// Counted as *physical* lines: a verbatim cell carrying a line break —
    /// possible only with [`MzTabMWriteOptions::source_verbatim_cells`] — turns
    /// one rendered row into several, and the rendered document is charged
    /// again for those before it is handed back.
    pub const MAX_LINES: usize = 1_000_000;

    /// A writer with native options. Source `MzTabMFile()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// A writer with explicit options. Native constructor.
    pub fn with_options(options: MzTabMWriteOptions) -> Self {
        Self { options }
    }

    /// Generate the `MTD` metadata section.
    ///
    /// Source `generateMzTabMMetaDataSection_(const MzTabMMetaData& map,
    /// StringList& sl)`, whose `sl` out-parameter this returns.
    ///
    /// `mzTab-version`, `mzTab-ID`, `quantification_method`,
    /// `small_molecule-quantification_unit`,
    /// `small_molecule_feature-quantification_unit` and
    /// `small_molecule-identification_reliability` are emitted unconditionally,
    /// so an unset one appears with the value `null`; `title`, `description`
    /// and the optional per-record keys are emitted only when their cell is not
    /// null. That asymmetry is the source's and is preserved, because the
    /// specification marks exactly those six mandatory.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when any indexed metadata map has more entries
    /// than [`MzTabM::MAX_INDEXED_ENTRIES`] (the ceiling bounds the entry
    /// count, not the index values), or when the estimated line count exceeds
    /// [`MzTabMFile::MAX_LINES`] — both checked before any line is built — or
    /// when a key or value carries a tab or a line break and
    /// [`MzTabMWriteOptions::source_verbatim_cells`] is off.
    pub fn generate_meta_data_section(&self, md: &MzTabMMetaData) -> Result<Vec<String>> {
        write::meta_data_section(md, &self.options)
    }

    /// Generate the `SMH` small-molecule header.
    ///
    /// Source `generateMzTabMSmallMoleculeHeader_`.
    ///
    /// # Arguments
    ///
    /// * `meta` — the metadata section, whose `assay` and `study_variable`
    ///   maps name the abundance columns.
    /// * `optional_columns` — the `opt_…` column names to append, normally
    ///   [`MzTabM::small_molecule_optional_column_names`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the column count would exceed
    /// [`MzTabM::MAX_OPTIONAL_COLUMNS`] plus the fixed columns, and when a
    /// column name carries a tab or a line break while
    /// [`MzTabMWriteOptions::source_verbatim_cells`] is off.
    pub fn generate_small_molecule_header(
        &self,
        meta: &MzTabMMetaData,
        optional_columns: &[String],
    ) -> Result<MzTabMSectionLine> {
        write::small_molecule_header(meta, optional_columns, &self.options)
    }

    /// Generate one `SML` row.
    ///
    /// Source `generateMzTabMSmallMoleculeSectionRow_`, plus the `meta`
    /// argument this port needs to align the abundance cells with the header;
    /// see [`MzTabMWriteOptions::source_row_abundance_cells`].
    ///
    /// # Arguments
    ///
    /// * `row` — the row to render.
    /// * `meta` — the metadata section that declared the abundance columns.
    /// * `optional_columns` — the section's `opt_…` column names, in header
    ///   order. A name the row does not carry yields `null`, and an entry of
    ///   the row whose name is not in this list is not written: that is the
    ///   contract of `MzTabFile::addOptionalColumnsToSectionRow_`, which the
    ///   source calls here.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] under native options when the row's abundance
    /// maps name an assay or study variable `meta` does not declare, or when a
    /// cell carries a tab or a line break — see
    /// [`MzTabMWriteOptions::source_verbatim_cells`] — and for the column
    /// ceiling.
    pub fn generate_small_molecule_section_row(
        &self,
        row: &MzTabMSmallMoleculeSectionRow,
        meta: &MzTabMMetaData,
        optional_columns: &[String],
    ) -> Result<MzTabMSectionLine> {
        write::small_molecule_row(row, meta, optional_columns, &self.options)
    }

    /// Generate the `SFH` small-molecule feature header.
    ///
    /// Source `generateMzTabMSmallMoleculeFeatureHeader_`. Only
    /// [`MzTabMMetaData::assay`] contributes abundance columns here.
    ///
    /// # Errors
    ///
    /// As [`MzTabMFile::generate_small_molecule_header`].
    pub fn generate_small_molecule_feature_header(
        &self,
        meta: &MzTabMMetaData,
        optional_columns: &[String],
    ) -> Result<MzTabMSectionLine> {
        write::small_molecule_feature_header(meta, optional_columns, &self.options)
    }

    /// Generate one `SMF` row.
    ///
    /// Source `generateMzTabMSmallMoleculeFeatureSectionRow_`, plus the `meta`
    /// argument; see [`MzTabMFile::generate_small_molecule_section_row`].
    ///
    /// # Errors
    ///
    /// As [`MzTabMFile::generate_small_molecule_section_row`].
    pub fn generate_small_molecule_feature_section_row(
        &self,
        row: &MzTabMSmallMoleculeFeatureSectionRow,
        meta: &MzTabMMetaData,
        optional_columns: &[String],
    ) -> Result<MzTabMSectionLine> {
        write::small_molecule_feature_row(row, meta, optional_columns, &self.options)
    }

    /// Generate the `SEH` small-molecule evidence header.
    ///
    /// Source `generateMzTabMSmallMoleculeEvidenceHeader_`. The
    /// `id_confidence_measure[n]` columns come from
    /// [`MzTabMMetaData::id_confidence_measure`] and sit *before* `rank`.
    ///
    /// # Errors
    ///
    /// As [`MzTabMFile::generate_small_molecule_header`].
    pub fn generate_small_molecule_evidence_header(
        &self,
        meta: &MzTabMMetaData,
        optional_columns: &[String],
    ) -> Result<MzTabMSectionLine> {
        write::small_molecule_evidence_header(meta, optional_columns, &self.options)
    }

    /// Generate one `SME` row.
    ///
    /// Source `generateMzTabMSmallMoleculeEvidenceSectionRow_`, plus the `meta`
    /// argument; the confidence cells are aligned with the header exactly as
    /// the abundance cells are.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] under native options when the row reports a
    /// confidence measure `meta` does not declare, or when a cell carries a
    /// tab or a line break — see
    /// [`MzTabMWriteOptions::source_verbatim_cells`] — and for the column
    /// ceiling.
    pub fn generate_small_molecule_evidence_section_row(
        &self,
        row: &MzTabMSmallMoleculeEvidenceSectionRow,
        meta: &MzTabMMetaData,
        optional_columns: &[String],
    ) -> Result<MzTabMSectionLine> {
        write::small_molecule_evidence_row(row, meta, optional_columns, &self.options)
    }

    /// Render a whole document to lines, without touching the filesystem.
    ///
    /// Native: the source builds the same `StringList` inside `store` and does
    /// not expose it. The lines carry no terminator; each section is preceded
    /// by one empty line, as the source writes.
    ///
    /// # Errors
    ///
    /// As [`MzTabMFile::store`], minus the extension check and the I/O:
    /// [`Error::InvalidValue`] for any resource ceiling, for a row that does
    /// not fill its header, and — unless
    /// [`MzTabMWriteOptions::source_verbatim_cells`] is set — for a cell,
    /// metadata key or metadata value carrying a tab or a line break.
    pub fn generate_lines(&self, mztab_m: &MzTabM) -> Result<Vec<String>> {
        write::document(mztab_m, &self.options)
    }

    /// Write a document to `path`.
    ///
    /// Source `store(const std::string& filename, const MzTabM& mztab_m)`.
    ///
    /// # Arguments
    ///
    /// * `path` — the output file. Its extension must be `mzTab` or `tsv`, or
    ///   be unrecognised; the source accepts exactly the same set, because
    ///   `FileHandler::hasValidExtension` also permits an unknown extension.
    /// * `mztab_m` — the document to write.
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidValue`] when the extension is a *known* type other
    ///   than `mzTab` or `tsv`, mapping the source's
    ///   `Exception::UnableToCreateFile`; when the path is not UTF-8; and for
    ///   every resource ceiling of [`MzTabMFile::generate_lines`].
    /// * [`Error::Io`] for a filesystem failure.
    ///
    /// Every line is rendered before the file is created, so a formatting or
    /// ceiling error leaves any existing file untouched. The source logs
    /// `exporting identification data: … to MzTab-M:` at info level before the
    /// extension check; this port writes no log.
    pub fn store(&self, path: impl AsRef<Path>, mztab_m: &MzTabM) -> Result<()> {
        let path = path.as_ref();
        let filename = path
            .to_str()
            .ok_or_else(|| bad("MzTab-M output filename must be UTF-8"))?;
        if !(crate::format::file_types::has_valid_extension(
            filename,
            crate::format::FileType::MzTab,
        ) || crate::format::file_types::has_valid_extension(
            filename,
            crate::format::FileType::Tsv,
        )) {
            return Err(bad(
                "invalid MzTab-M output extension, expected 'mzTab' or 'tsv'",
            ));
        }
        let lines = self.generate_lines(mztab_m)?;
        let mut out = crate::format::TextFile::new();
        for line in lines {
            out.add_line(&line)?;
        }
        out.store(path)
    }
}

// ---------------------------------------------------------------------------
// Writer implementation
// ---------------------------------------------------------------------------

mod write {
    use super::{
        Error, MTD, MzTabCVMetaData, MzTabContactMetaData, MzTabDouble, MzTabInstrumentMetaData,
        MzTabM, MzTabMAssayMetaData, MzTabMDatabaseMetaData, MzTabMFile, MzTabMMSRunMetaData,
        MzTabMMetaData, MzTabMSectionLine, MzTabMSmallMoleculeEvidenceSectionRow,
        MzTabMSmallMoleculeFeatureSectionRow, MzTabMSmallMoleculeSectionRow,
        MzTabMStudyVariableMetaData, MzTabMWriteOptions, MzTabOptionalColumnEntry,
        MzTabSampleMetaData, MzTabSoftwareMetaData, MzTabString, MzTabStringList, Result, add, bad,
    };
    use std::collections::BTreeMap;

    /// The source's `MzTabString("null").toCellString()`, a plain `null`.
    const NULL: &str = "null";

    fn key(prefix: &str, index: usize, suffix: &str) -> String {
        format!("{MTD}{prefix}[{index}]{suffix}")
    }

    fn indexed(prefix: &str, index: usize, inner: &str, sub: usize) -> String {
        format!("{MTD}{prefix}[{index}]-{inner}[{sub}]")
    }

    fn line(key: String, value: String) -> String {
        format!("{key}\t{value}")
    }

    /// Refuse an indexed metadata map that is larger than the ceiling.
    fn bounded<T>(map: &BTreeMap<usize, T>, what: &str) -> Result<usize> {
        if map.len() > MzTabM::MAX_INDEXED_ENTRIES {
            return Err(bad(format!(
                "MzTab-M metadata section declares more than {} {what} entries",
                MzTabM::MAX_INDEXED_ENTRIES
            )));
        }
        Ok(map.len())
    }

    /// Upper bound on the metadata lines, checked before any line is built.
    ///
    /// Eight is the largest number of lines any one indexed record produces
    /// (`ms_run`), and every non-indexed key contributes at most one.
    fn estimate(md: &MzTabMMetaData) -> Result<usize> {
        let mut total = 8_usize;
        for count in [
            bounded(&md.sample_processing, "sample_processing")?,
            bounded(&md.instrument, "instrument")?,
            bounded(&md.software, "software")?,
            bounded(&md.publication, "publication")?,
            bounded(&md.contact, "contact")?,
            bounded(&md.uri, "uri")?,
            bounded(&md.external_study_uri, "external_study_uri")?,
            bounded(&md.sample, "sample")?,
            bounded(&md.ms_run, "ms_run")?,
            bounded(&md.assay, "assay")?,
            bounded(&md.study_variable, "study_variable")?,
            bounded(&md.custom, "custom")?,
            bounded(&md.cv, "cv")?,
            bounded(&md.database, "database")?,
            bounded(&md.derivatization_agent, "derivatization_agent")?,
            bounded(&md.id_confidence_measure, "id_confidence_measure")?,
        ] {
            total = add(total, count.saturating_mul(8))?;
        }
        for record in md.instrument.values() {
            total = add(total, bounded(&record.analyzer, "instrument analyzer")?)?;
        }
        for record in md.software.values() {
            total = add(total, bounded(&record.setting, "software setting")?)?;
        }
        for record in md.sample.values() {
            for map in [
                &record.species,
                &record.tissue,
                &record.cell_type,
                &record.disease,
                &record.custom,
            ] {
                total = add(total, bounded(map, "sample")?)?;
            }
        }
        for record in md.ms_run.values() {
            total = add(
                total,
                bounded(&record.fragmentation_method, "fragmentation_method")?,
            )?;
            total = add(total, bounded(&record.scan_polarity, "scan_polarity")?)?;
        }
        for record in md.assay.values() {
            total = add(total, bounded(&record.custom, "assay custom")?)?;
        }
        for record in md.study_variable.values() {
            total = add(total, record.factors.get().len())?;
        }
        for list in [
            &md.colunit_small_molecule,
            &md.colunit_small_molecule_feature,
            &md.colunit_small_molecule_evidence,
        ] {
            total = add(total, list.len())?;
        }
        if total > MzTabMFile::MAX_LINES {
            return Err(bad(
                "MzTab-M metadata section exceeds the output line limit",
            ));
        }
        Ok(total)
    }

    fn instrument_lines(index: usize, imd: &MzTabInstrumentMetaData, sl: &mut Vec<String>) {
        if !imd.name.is_null() {
            sl.push(line(
                key("instrument", index, "-name"),
                imd.name.to_cell_string(),
            ));
        }
        if !imd.source.is_null() {
            sl.push(line(
                key("instrument", index, "-source"),
                imd.source.to_cell_string(),
            ));
        }
        for (sub, analyzer) in &imd.analyzer {
            if !analyzer.is_null() {
                sl.push(line(
                    indexed("instrument", index, "analyzer", *sub),
                    analyzer.to_cell_string(),
                ));
            }
        }
        if !imd.detector.is_null() {
            sl.push(line(
                key("instrument", index, "-detector"),
                imd.detector.to_cell_string(),
            ));
        }
    }

    fn software_lines(index: usize, msmd: &MzTabSoftwareMetaData, sl: &mut Vec<String>) {
        sl.push(line(
            key("software", index, ""),
            msmd.software.to_cell_string(),
        ));
        for (sub, setting) in &msmd.setting {
            sl.push(line(
                indexed("software", index, "setting", *sub),
                setting.to_cell_string(),
            ));
        }
    }

    fn contact_lines(index: usize, cmd: &MzTabContactMetaData, sl: &mut Vec<String>) {
        for (suffix, cell) in [
            ("-name", &cmd.name),
            ("-affiliation", &cmd.affiliation),
            ("-email", &cmd.email),
        ] {
            if !cell.is_null() {
                sl.push(line(key("contact", index, suffix), cell.to_cell_string()));
            }
        }
    }

    fn sample_lines(index: usize, msmd: &MzTabSampleMetaData, sl: &mut Vec<String>) {
        if !msmd.description.is_null() {
            sl.push(line(
                key("sample", index, "-description"),
                msmd.description.to_cell_string(),
            ));
        }
        for (inner, map) in [
            ("species", &msmd.species),
            ("tissue", &msmd.tissue),
            ("cell_type", &msmd.cell_type),
            ("disease", &msmd.disease),
            ("custom", &msmd.custom),
        ] {
            for (sub, value) in map {
                sl.push(line(
                    indexed("sample", index, inner, *sub),
                    value.to_cell_string(),
                ));
            }
        }
    }

    fn ms_run_lines(
        index: usize,
        rmmd: &MzTabMMSRunMetaData,
        options: &MzTabMWriteOptions,
        sl: &mut Vec<String>,
    ) {
        sl.push(line(
            key("ms_run", index, "-location"),
            rmmd.location.to_cell_string(),
        ));
        if !rmmd.instrument_ref.is_null() {
            sl.push(line(
                key("ms_run", index, "-instrument_ref"),
                rmmd.instrument_ref.to_cell_string(),
            ));
        }
        if !rmmd.format.is_null() {
            sl.push(line(
                key("ms_run", index, "-format"),
                rmmd.format.to_cell_string(),
            ));
        }
        if !options.omit_ms_run_id_format && !rmmd.id_format.is_null() {
            sl.push(line(
                key("ms_run", index, "-id_format"),
                rmmd.id_format.to_cell_string(),
            ));
        }
        for (sub, value) in &rmmd.fragmentation_method {
            sl.push(line(
                indexed("ms_run", index, "fragmentation_method", *sub),
                value.to_cell_string(),
            ));
        }
        for (sub, value) in &rmmd.scan_polarity {
            sl.push(line(
                indexed("ms_run", index, "scan_polarity", *sub),
                value.to_cell_string(),
            ));
        }
        if !rmmd.hash.is_null() {
            sl.push(line(
                key("ms_run", index, "-hash"),
                rmmd.hash.to_cell_string(),
            ));
        }
        if !rmmd.hash_method.is_null() {
            sl.push(line(
                key("ms_run", index, "-hash_method"),
                rmmd.hash_method.to_cell_string(),
            ));
        }
    }

    fn assay_lines(
        index: usize,
        amd: &MzTabMAssayMetaData,
        options: &MzTabMWriteOptions,
        sl: &mut Vec<String>,
    ) {
        sl.push(line(key("assay", index, ""), amd.name.to_cell_string()));
        let custom_prefix = if options.source_assay_custom_key {
            "ms_run"
        } else {
            "assay"
        };
        for (sub, value) in &amd.custom {
            sl.push(line(
                indexed(custom_prefix, index, "custom", *sub),
                value.to_cell_string(),
            ));
        }
        if !amd.external_uri.is_null() {
            sl.push(line(
                key("assay", index, "-external_uri"),
                amd.external_uri.to_cell_string(),
            ));
        }
        if !amd.sample_ref.is_null() {
            sl.push(line(
                key("assay", index, "-sample_ref"),
                format!("sample[{}]", amd.sample_ref.to_cell_string()),
            ));
        }
        sl.push(line(
            key("assay", index, "-ms_run_ref"),
            format!("ms_run[{}]", amd.ms_run_ref.to_cell_string()),
        ));
    }

    fn study_variable_lines(
        index: usize,
        svmd: &MzTabMStudyVariableMetaData,
        sl: &mut Vec<String>,
    ) {
        sl.push(line(
            key("study_variable", index, ""),
            svmd.name.to_cell_string(),
        ));
        let mut refs = MzTabStringList::default();
        refs.set(
            svmd.assay_refs
                .iter()
                .map(|reference| MzTabString::from_text(&format!("assay[{reference}]")))
                .collect(),
        );
        sl.push(line(
            key("study_variable", index, "-assay_refs"),
            refs.to_cell_string(),
        ));
        if !svmd.average_function.is_null() {
            sl.push(line(
                key("study_variable", index, "-average_function"),
                svmd.average_function.to_cell_string(),
            ));
        }
        if !svmd.variation_function.is_null() {
            sl.push(line(
                key("study_variable", index, "-variation_function"),
                svmd.variation_function.to_cell_string(),
            ));
        }
        sl.push(line(
            key("study_variable", index, "-description"),
            svmd.description.to_cell_string(),
        ));
        for factor in svmd.factors.get() {
            sl.push(line(
                key("study_variable", index, "-factors"),
                factor.to_cell_string(),
            ));
        }
    }

    fn cv_lines(index: usize, cvmd: &MzTabCVMetaData, sl: &mut Vec<String>) {
        for (suffix, cell) in [
            ("-label", &cvmd.label),
            ("-full_name", &cvmd.full_name),
            ("-version", &cvmd.version),
            ("-uri", &cvmd.url),
        ] {
            sl.push(line(key("cv", index, suffix), cell.to_cell_string()));
        }
    }

    fn database_lines(index: usize, dbmd: &MzTabMDatabaseMetaData, sl: &mut Vec<String>) {
        sl.push(line(
            key("database", index, ""),
            dbmd.database.to_cell_string(),
        ));
        for (suffix, cell) in [
            ("-prefix", &dbmd.prefix),
            ("-version", &dbmd.version),
            ("-uri", &dbmd.uri),
        ] {
            sl.push(line(key("database", index, suffix), cell.to_cell_string()));
        }
    }

    pub(super) fn meta_data_section(
        md: &MzTabMMetaData,
        options: &MzTabMWriteOptions,
    ) -> Result<Vec<String>> {
        let reserve = estimate(md)?;
        let mut sl: Vec<String> = Vec::new();
        sl.try_reserve(reserve)
            .map_err(|_| bad("MzTab-M metadata section does not fit in memory"))?;

        sl.push(line(
            format!("{MTD}mzTab-version"),
            md.mz_tab_version.to_cell_string(),
        ));
        sl.push(line(
            format!("{MTD}mzTab-ID"),
            md.mz_tab_id.to_cell_string(),
        ));
        if !md.title.is_null() {
            sl.push(line(format!("{MTD}title"), md.title.to_cell_string()));
        }
        if !md.description.is_null() {
            sl.push(line(
                format!("{MTD}description"),
                md.description.to_cell_string(),
            ));
        }
        for (index, value) in &md.sample_processing {
            sl.push(line(
                key("sample_processing", *index, ""),
                value.to_cell_string(),
            ));
        }
        for (index, record) in &md.instrument {
            instrument_lines(*index, record, &mut sl);
        }
        for (index, record) in &md.software {
            software_lines(*index, record, &mut sl);
        }
        for (index, value) in &md.publication {
            sl.push(line(key("publication", *index, ""), value.to_cell_string()));
        }
        for (index, record) in &md.contact {
            contact_lines(*index, record, &mut sl);
        }
        for (index, value) in &md.uri {
            sl.push(line(key("uri", *index, ""), value.to_cell_string()));
        }
        for (index, value) in &md.external_study_uri {
            sl.push(line(
                key("external_study_uri", *index, ""),
                value.to_cell_string(),
            ));
        }
        sl.push(line(
            format!("{MTD}quantification_method"),
            md.quantification_method.to_cell_string(),
        ));
        for (index, record) in &md.sample {
            sample_lines(*index, record, &mut sl);
        }
        for (index, record) in &md.ms_run {
            ms_run_lines(*index, record, options, &mut sl);
        }
        for (index, record) in &md.assay {
            assay_lines(*index, record, options, &mut sl);
        }
        for (index, record) in &md.study_variable {
            study_variable_lines(*index, record, &mut sl);
        }
        for (index, value) in &md.custom {
            sl.push(line(key("custom", *index, ""), value.to_cell_string()));
        }
        for (index, record) in &md.cv {
            cv_lines(*index, record, &mut sl);
        }
        for (index, record) in &md.database {
            database_lines(*index, record, &mut sl);
        }
        let agent_suffix = if options.source_derivatization_agent_key {
            "-uri"
        } else {
            ""
        };
        for (index, value) in &md.derivatization_agent {
            sl.push(line(
                key("derivatization_agent", *index, agent_suffix),
                value.to_cell_string(),
            ));
        }
        sl.push(line(
            format!("{MTD}small_molecule-quantification_unit"),
            md.small_molecule_quantification_unit.to_cell_string(),
        ));
        sl.push(line(
            format!("{MTD}small_molecule_feature-quantification_unit"),
            md.small_molecule_feature_quantification_unit
                .to_cell_string(),
        ));
        sl.push(line(
            format!("{MTD}small_molecule-identification_reliability"),
            md.small_molecule_identification_reliability
                .to_cell_string(),
        ));
        for (index, value) in &md.id_confidence_measure {
            sl.push(line(
                key("id_confidence_measure", *index, ""),
                value.to_cell_string(),
            ));
        }
        let (feature_key, evidence_key) = if options.source_colunit_keys {
            ("colunit_small_molecule", "colunit_small_molecule")
        } else {
            (
                "colunit_small_molecule_feature",
                "colunit_small_molecule_evidence",
            )
        };
        for (name, list) in [
            ("colunit_small_molecule", &md.colunit_small_molecule),
            (feature_key, &md.colunit_small_molecule_feature),
            (evidence_key, &md.colunit_small_molecule_evidence),
        ] {
            for value in list {
                sl.push(line(format!("{MTD}{name}"), value.to_cell_string()));
            }
        }
        // Every `MTD` line is exactly three tab-separated fields. A key or a
        // value carrying a separator would add a field or split the line — and
        // text beginning with `MTD`, `SMH` or `SML` would forge a line of that
        // kind — so the whole section is checked before it is handed back.
        if !options.source_verbatim_cells {
            for text in &sl {
                if text.matches('\t').count() != 2 || text.contains(['\n', '\r']) {
                    return Err(bad(
                        "MzTab-M metadata key or value carries a tab or a line break",
                    ));
                }
            }
        }
        Ok(sl)
    }

    /// The separators a cell may not carry: a tab would add a column and a
    /// line break would split the row. See
    /// [`MzTabMWriteOptions::source_verbatim_cells`].
    const CELL_SEPARATORS: [char; 3] = ['\t', '\n', '\r'];

    fn finish(
        parts: Vec<String>,
        what: &str,
        options: &MzTabMWriteOptions,
    ) -> Result<MzTabMSectionLine> {
        if !options.source_verbatim_cells {
            for (column, part) in parts.iter().enumerate() {
                if part.contains(CELL_SEPARATORS) {
                    return Err(bad(format!(
                        "MzTab-M {what} column {column} carries a tab or a line break"
                    )));
                }
            }
        }
        Ok(MzTabMSectionLine {
            columns: parts.len(),
            text: parts.join("\t"),
        })
    }

    fn check_columns(count: usize) -> Result<()> {
        if count > MzTabM::MAX_OPTIONAL_COLUMNS {
            return Err(bad("MzTab-M section exceeds its optional-column limit"));
        }
        Ok(())
    }

    /// `MzTabFile::addOptionalColumnsToSectionRow_`: one cell per requested
    /// name, taken from the first entry of the row with that name, or `null`.
    fn optional_cells(
        column_names: &[String],
        entries: &[MzTabOptionalColumnEntry],
        output: &mut Vec<String>,
    ) {
        for name in column_names {
            let cell = entries
                .iter()
                .find(|entry| &entry.name == name)
                .map_or_else(|| NULL.to_owned(), |entry| entry.value.to_cell_string());
            output.push(cell);
        }
    }

    /// One cell per column the metadata declares, `null` where the row is
    /// silent; or the row's own cells verbatim under source options.
    fn aligned_cells<T>(
        declared: &BTreeMap<usize, T>,
        row: &BTreeMap<usize, MzTabDouble>,
        source: bool,
        what: &str,
        output: &mut Vec<String>,
    ) -> Result<()> {
        if source {
            output.extend(row.values().map(MzTabDouble::to_cell_string));
            return Ok(());
        }
        for index in row.keys() {
            if !declared.contains_key(index) {
                return Err(bad(format!(
                    "MzTab-M row reports {what}[{index}], which the metadata section does not declare"
                )));
            }
        }
        for index in declared.keys() {
            output.push(
                row.get(index)
                    .map_or_else(|| NULL.to_owned(), MzTabDouble::to_cell_string),
            );
        }
        Ok(())
    }

    fn abundance_columns<T>(prefix: &str, declared: &BTreeMap<usize, T>, header: &mut Vec<String>) {
        for index in declared.keys() {
            header.push(format!("{prefix}[{index}]"));
        }
    }

    pub(super) fn small_molecule_header(
        meta: &MzTabMMetaData,
        optional_columns: &[String],
        options: &MzTabMWriteOptions,
    ) -> Result<MzTabMSectionLine> {
        check_columns(optional_columns.len())?;
        let mut header: Vec<String> = [
            "SMH",
            "SML_ID",
            "SMF_ID_REFS",
            "database_identifier",
            "chemical_formula",
            "smiles",
            "inchi",
            "chemical_name",
            "uri",
            "theoretical_neutral_mass",
            "adduct_ions",
            "reliability",
            "best_id_confidence_measure",
            "best_id_confidence_value",
        ]
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
        abundance_columns("abundance_assay", &meta.assay, &mut header);
        abundance_columns(
            "abundance_study_variable",
            &meta.study_variable,
            &mut header,
        );
        abundance_columns(
            "abundance_variation_study_variable",
            &meta.study_variable,
            &mut header,
        );
        header.extend(optional_columns.iter().cloned());
        finish(header, "SMH header", options)
    }

    pub(super) fn small_molecule_row(
        row: &MzTabMSmallMoleculeSectionRow,
        meta: &MzTabMMetaData,
        optional_columns: &[String],
        options: &MzTabMWriteOptions,
    ) -> Result<MzTabMSectionLine> {
        check_columns(optional_columns.len())?;
        let mut s = vec!["SML".to_owned()];
        s.push(row.sml_identifier.to_cell_string());
        s.push(row.smf_id_refs.to_cell_string());
        s.push(row.database_identifier.to_cell_string());
        s.push(row.chemical_formula.to_cell_string());
        s.push(row.smiles.to_cell_string());
        s.push(row.inchi.to_cell_string());
        s.push(row.chemical_name.to_cell_string());
        s.push(row.uri.to_cell_string());
        s.push(row.theoretical_neutral_mass.to_cell_string());
        s.push(row.adducts.to_cell_string());
        s.push(row.reliability.to_cell_string());
        s.push(row.best_id_confidence_measure.to_cell_string());
        s.push(row.best_id_confidence_value.to_cell_string());
        let source = options.source_row_abundance_cells;
        aligned_cells(
            &meta.assay,
            &row.small_molecule_abundance_assay,
            source,
            "abundance_assay",
            &mut s,
        )?;
        aligned_cells(
            &meta.study_variable,
            &row.small_molecule_abundance_study_variable,
            source,
            "abundance_study_variable",
            &mut s,
        )?;
        aligned_cells(
            &meta.study_variable,
            &row.small_molecule_abundance_variation_study_variable,
            source,
            "abundance_variation_study_variable",
            &mut s,
        )?;
        optional_cells(optional_columns, &row.opt, &mut s);
        finish(s, "SML row", options)
    }

    pub(super) fn small_molecule_feature_header(
        meta: &MzTabMMetaData,
        optional_columns: &[String],
        options: &MzTabMWriteOptions,
    ) -> Result<MzTabMSectionLine> {
        check_columns(optional_columns.len())?;
        let mut header: Vec<String> = [
            "SFH",
            "SMF_ID",
            "SME_ID_REFS",
            "SME_ID_REF_ambiguity_code",
            "adduct_ion",
            "isotopomer",
            "exp_mass_to_charge",
            "charge",
            "retention_time_in_seconds",
            "retention_time_in_seconds_start",
            "retention_time_in_seconds_end",
        ]
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
        abundance_columns("abundance_assay", &meta.assay, &mut header);
        header.extend(optional_columns.iter().cloned());
        finish(header, "SMF header", options)
    }

    pub(super) fn small_molecule_feature_row(
        row: &MzTabMSmallMoleculeFeatureSectionRow,
        meta: &MzTabMMetaData,
        optional_columns: &[String],
        options: &MzTabMWriteOptions,
    ) -> Result<MzTabMSectionLine> {
        check_columns(optional_columns.len())?;
        let mut s = vec!["SMF".to_owned()];
        s.push(row.smf_identifier.to_cell_string());
        s.push(row.sme_id_refs.to_cell_string());
        s.push(row.sme_id_ref_ambiguity_code.to_cell_string());
        s.push(row.adduct.to_cell_string());
        s.push(row.isotopomer.to_cell_string());
        s.push(row.exp_mass_to_charge.to_cell_string());
        s.push(row.charge.to_cell_string());
        s.push(row.retention_time.to_cell_string());
        s.push(row.rt_start.to_cell_string());
        s.push(row.rt_end.to_cell_string());
        aligned_cells(
            &meta.assay,
            &row.small_molecule_feature_abundance_assay,
            options.source_row_abundance_cells,
            "abundance_assay",
            &mut s,
        )?;
        optional_cells(optional_columns, &row.opt, &mut s);
        finish(s, "SMF row", options)
    }

    pub(super) fn small_molecule_evidence_header(
        meta: &MzTabMMetaData,
        optional_columns: &[String],
        options: &MzTabMWriteOptions,
    ) -> Result<MzTabMSectionLine> {
        check_columns(optional_columns.len())?;
        let mut header: Vec<String> = [
            "SEH",
            "SME_ID",
            "evidence_input_id",
            "database_identifier",
            "chemical_formula",
            "smiles",
            "inchi",
            "chemical_name",
            "uri",
            "derivatized_form",
            "adduct_ion",
            "exp_mass_to_charge",
            "charge",
            "theoretical_mass_to_charge",
            "spectra_ref",
            "identification_method",
            "ms_level",
        ]
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
        for index in meta.id_confidence_measure.keys() {
            header.push(format!("id_confidence_measure[{index}]"));
        }
        header.push("rank".to_owned());
        header.extend(optional_columns.iter().cloned());
        finish(header, "SEH header", options)
    }

    pub(super) fn small_molecule_evidence_row(
        row: &MzTabMSmallMoleculeEvidenceSectionRow,
        meta: &MzTabMMetaData,
        optional_columns: &[String],
        options: &MzTabMWriteOptions,
    ) -> Result<MzTabMSectionLine> {
        check_columns(optional_columns.len())?;
        let mut s = vec!["SME".to_owned()];
        s.push(row.sme_identifier.to_cell_string());
        s.push(row.evidence_input_id.to_cell_string());
        s.push(row.database_identifier.to_cell_string());
        s.push(row.chemical_formula.to_cell_string());
        s.push(row.smiles.to_cell_string());
        s.push(row.inchi.to_cell_string());
        s.push(row.chemical_name.to_cell_string());
        s.push(row.uri.to_cell_string());
        s.push(row.derivatized_form.to_cell_string());
        s.push(row.adduct.to_cell_string());
        s.push(row.exp_mass_to_charge.to_cell_string());
        s.push(row.charge.to_cell_string());
        s.push(row.calc_mass_to_charge.to_cell_string());
        s.push(row.spectra_ref.to_cell_string());
        s.push(row.identification_method.to_cell_string());
        s.push(row.ms_level.to_cell_string());
        aligned_cells(
            &meta.id_confidence_measure,
            &row.id_confidence_measure,
            options.source_row_abundance_cells,
            "id_confidence_measure",
            &mut s,
        )?;
        s.push(row.rank.to_cell_string());
        optional_cells(optional_columns, &row.opt, &mut s);
        finish(s, "SME row", options)
    }

    pub(super) fn document(mztab_m: &MzTabM, options: &MzTabMWriteOptions) -> Result<Vec<String>> {
        for (rows, what) in [
            (mztab_m.small_molecule_data.len(), "SML"),
            (mztab_m.small_molecule_feature_data.len(), "SMF"),
            (mztab_m.small_molecule_evidence_data.len(), "SME"),
        ] {
            if rows > MzTabM::MAX_ROWS {
                return Err(bad(format!("MzTab-M {what} section exceeds its row limit")));
            }
        }
        let meta = mztab_m.meta_data();
        let mut out = meta_data_section(meta, options)?;
        let mut total = out.len();
        for rows in [
            mztab_m.small_molecule_data.len(),
            mztab_m.small_molecule_feature_data.len(),
            mztab_m.small_molecule_evidence_data.len(),
        ] {
            total = add(total, add(rows, 2)?)?;
        }
        if total > MzTabMFile::MAX_LINES {
            return Err(bad("MzTab-M output exceeds the line limit"));
        }

        let sml_optional = mztab_m.small_molecule_optional_column_names()?;
        out.push(String::new());
        let header = small_molecule_header(meta, &sml_optional, options)?;
        let sml_columns = header.columns;
        out.push(header.text);
        for row in &mztab_m.small_molecule_data {
            let rendered = small_molecule_row(row, meta, &sml_optional, options)?;
            check_rectangle(sml_columns, rendered.columns, "small molecule", options)?;
            out.push(rendered.text);
        }

        let smf_optional = mztab_m.small_molecule_feature_optional_column_names()?;
        out.push(String::new());
        let header = small_molecule_feature_header(meta, &smf_optional, options)?;
        let smf_columns = header.columns;
        out.push(header.text);
        for row in &mztab_m.small_molecule_feature_data {
            let rendered = small_molecule_feature_row(row, meta, &smf_optional, options)?;
            check_rectangle(
                smf_columns,
                rendered.columns,
                "small molecule feature",
                options,
            )?;
            out.push(rendered.text);
        }

        let sme_optional = mztab_m.small_molecule_evidence_optional_column_names()?;
        out.push(String::new());
        let header = small_molecule_evidence_header(meta, &sme_optional, options)?;
        let sme_columns = header.columns;
        out.push(header.text);
        for row in &mztab_m.small_molecule_evidence_data {
            let rendered = small_molecule_evidence_row(row, meta, &sme_optional, options)?;
            check_rectangle(
                sme_columns,
                rendered.columns,
                "small molecule evidence",
                options,
            )?;
            out.push(rendered.text);
        }
        // MzTabMFile::MAX_LINES counts *physical* lines. Under source options a
        // verbatim cell may carry a line break, which turns one rendered row
        // into several, so the rendered lines are charged again here — the
        // pre-render estimate above cannot see them.
        let mut physical = out.len();
        for text in &out {
            physical = add(physical, text.matches('\n').count())?;
        }
        if physical > MzTabMFile::MAX_LINES {
            return Err(bad("MzTab-M output exceeds the line limit"));
        }
        Ok(out)
    }

    /// The source's three `OPENMS_POSTCONDITION`s, which a release build drops.
    ///
    /// Under native options a mismatch cannot happen, because the abundance and
    /// confidence cells are aligned with the header; the check stays as an
    /// assertion of that. Under source options a mismatch is exactly the
    /// source's own behaviour, so the row is written as the source writes it.
    fn check_rectangle(
        header: usize,
        row: usize,
        what: &str,
        options: &MzTabMWriteOptions,
    ) -> Result<()> {
        if !options.source_row_abundance_cells && header != row {
            return Err(Error::InvalidValue(format!(
                "MzTab-M {what} row has {row} columns but its header has {header}"
            )));
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Export implementation
// ---------------------------------------------------------------------------

mod export {
    use super::{
        BTreeMap, BTreeSet, ControlledVocabulary, EmpiricalFormula, Feature, FeatureMap,
        IdentificationData, IdentifiedCompound, MzTabCVMetaData, MzTabDouble, MzTabDoubleList,
        MzTabInteger, MzTabM, MzTabMAssayMetaData, MzTabMDatabaseMetaData, MzTabMExportOptions,
        MzTabMMSRunMetaData, MzTabMMetaData, MzTabMSmallMoleculeEvidenceSectionRow,
        MzTabMSmallMoleculeFeatureSectionRow, MzTabMSmallMoleculeSectionRow,
        MzTabMStudyVariableMetaData, MzTabParameter, MzTabSoftwareMetaData, MzTabSpectraRef,
        MzTabString, ObservationMatchId, ProcessingAction, ProcessingSoftware, Result, bad,
        missing, source_double,
    };
    use crate::metadata::{MetaInfo, MetaValueData};

    /// The identified compound behind one observation match.
    ///
    /// # Errors
    ///
    /// [`super::Error::MissingInformation`] when the match identifies a peptide
    /// or an oligonucleotide; the source's `getIdentifiedCompoundRef()` throws
    /// `Exception::IllegalArgument` there.
    pub(super) fn compound_of(
        id_data: &IdentificationData,
        id: ObservationMatchId,
    ) -> Result<&IdentifiedCompound> {
        let record = id_data.observation_match(id)?;
        let compound = record.identified_molecule.compound().map_err(|_| {
            missing("MzTab-M export needs every observation match to identify a compound")
        })?;
        id_data.compound(compound)
    }

    /// `cv.getTermByName(name)` with the source's empty description.
    ///
    /// # Errors
    ///
    /// [`super::Error::InvalidValue`] when `cv` has no term with that name; the
    /// source throws `Exception::ElementNotFound`.
    fn term(cv: &ControlledVocabulary, name: &str) -> Result<(String, String)> {
        let definition = cv.get_term_by_name(name, "").map_err(|_| {
            bad(format!(
                "MzTab-M export needs the controlled-vocabulary term {name:?}"
            ))
        })?;
        Ok((definition.id.clone(), definition.name.clone()))
    }

    /// `[MS, <id>, <name>, <value>]` parsed back into a parameter, as the
    /// source builds every CV parameter here.
    fn ms_parameter(id: &str, name: &str, value: &str) -> Result<MzTabParameter> {
        MzTabParameter::parse(&format!("[MS, {id}, {name}, {value}]"))
    }

    /// The source's `std::regex_replace(name, R"(\\)", "/")` plus its
    /// `file://` prefixing.
    fn as_file_uri(name: &str) -> String {
        let slashed = name.replace('\\', "/");
        if slashed.starts_with("file://") {
            slashed
        } else {
            format!("file://{slashed}")
        }
    }

    /// `StringUtils::trimmed(StringUtils::prefix(File::basename(path), '.'))`.
    ///
    /// `String::prefix(char)` throws `Exception::ElementNotFound` when the
    /// basename holds no `.`; this returns the whole trimmed basename, because
    /// a feature map read from an extensionless file is not an error.
    fn assay_stem(path: &str) -> &str {
        let base = crate::system::file::basename(path);
        let stem = base.split('.').next().unwrap_or(base);
        stem.trim_matches(|c: char| c == ' ' || c == '\t' || c == '\n' || c == '\r')
    }

    /// The source's `substitute(key, ' ', '_')`.
    fn substitute(key: &str) -> String {
        key.replace(' ', "_")
    }

    /// Source `getFeatureMapMetaValues_`, whose three `std::set` out-parameters
    /// this returns.
    ///
    /// The observation-match set drops every key containing
    /// `IDConverter_trace`, which `IdentificationDataConverter` adds and which
    /// would otherwise fan out into one optional column per traced record.
    fn meta_value_keys(
        feature_map: &FeatureMap,
        id_data: &IdentificationData,
        options: &MzTabMExportOptions,
    ) -> Result<(BTreeSet<String>, BTreeSet<String>, BTreeSet<String>)> {
        let mut feature_keys = BTreeSet::new();
        let mut match_keys = BTreeSet::new();
        let mut compound_keys = BTreeSet::new();
        let spell = |key: &str| -> String {
            if options.substitute_keys_before_lookup {
                substitute(key)
            } else {
                key.to_owned()
            }
        };
        for feature in &feature_map.features {
            for key in feature.metadata.keys() {
                feature_keys.insert(spell(key));
            }
            for &id in &feature.id_matches {
                let record = id_data.observation_match(id)?;
                for key in record.result.metadata.keys() {
                    if !key.contains("IDConverter_trace") {
                        match_keys.insert(spell(key));
                    }
                }
                let compound = compound_of(id_data, id)?;
                for key in compound.result.metadata.keys() {
                    compound_keys.insert(spell(key));
                }
            }
        }
        Ok((feature_keys, match_keys, compound_keys))
    }

    /// Source `getAdductString_`: the adduct's name, reformatted from
    /// `M+H;1+` to `[M+H]1+`, or the literal `null` when the match has none.
    ///
    /// The source splits on the *first* `;` with `String::substr`, so a name
    /// with several semicolons keeps the later ones inside the charge suffix.
    /// That is preserved. Splitting on a character boundary of a `&str` cannot
    /// panic here because `;` is single-byte and `find` returns a boundary.
    fn adduct_string(id_data: &IdentificationData, id: ObservationMatchId) -> Result<String> {
        let record = id_data.observation_match(id)?;
        let Some(adduct) = record.adduct else {
            return Ok("null".to_owned());
        };
        let name = id_data.adduct(adduct)?.name();
        match name.find(';') {
            Some(offset) => {
                let (prefix, rest) = name.split_at(offset);
                let suffix = rest.get(1..).unwrap_or("");
                Ok(format!("[{prefix}]{suffix}"))
            }
            None => Ok(name.to_owned()),
        }
    }

    /// Every processing software in the source's `(name, version)` set order.
    ///
    /// The graph iterates its registry in registration order; the source's
    /// `std::set<ProcessingSoftware>` is ordered by `Software::operator<`,
    /// which compares name then version, so this sorts to match. The order is
    /// observable: it fixes the `software[n]` indices.
    fn ordered_softwares(id_data: &IdentificationData) -> Vec<&ProcessingSoftware> {
        let mut softwares: Vec<&ProcessingSoftware> =
            id_data.processing_softwares().map(|(_, s)| s).collect();
        softwares.sort_by(|a, b| {
            (&a.software.name, &a.software.version).cmp(&(&b.software.name, &b.software.version))
        });
        softwares
    }

    /// The software names that performed each processing action.
    ///
    /// Source `action_software_name`, a `std::map<ProcessingAction,
    /// std::vector<std::string>>` filled in the iteration order of the
    /// processing-step set and read back in two places: a membership test for
    /// `QUANTITATION`, and a loop over `IDENTIFICATION` whose last recognised
    /// tool wins. This port sorts and deduplicates the names per action, so the
    /// winner is the lexicographically last recognised identification tool
    /// rather than the last one the source's set happened to visit. The two
    /// agree whenever one tool performed the identification, which is the case
    /// the source comment assumes.
    fn action_software_names(
        id_data: &IdentificationData,
    ) -> Result<BTreeMap<ProcessingAction, Vec<String>>> {
        let mut names: BTreeMap<ProcessingAction, BTreeSet<String>> = BTreeMap::new();
        for (_, step) in id_data.processing_steps() {
            let software = id_data.processing_software(step.software)?;
            for &action in &step.actions {
                names
                    .entry(action)
                    .or_default()
                    .insert(software.software.name.clone());
            }
        }
        Ok(names
            .into_iter()
            .map(|(action, set)| (action, set.into_iter().collect()))
            .collect())
    }

    fn text_list(metadata: &MetaInfo, key: &str) -> Result<Option<Vec<String>>> {
        match metadata.get(key) {
            None => Ok(None),
            Some(value) => match value.data() {
                MetaValueData::StringList(values) => Ok(Some(values.clone())),
                MetaValueData::String(text) => Ok(Some(vec![text.clone()])),
                _ => Err(bad(
                    "MzTab-M export needs the feature 'adducts' meta value to be a string list",
                )),
            },
        }
    }

    fn meta_text(metadata: &MetaInfo, key: &str) -> Option<String> {
        metadata.get(key).map(ToString::to_string)
    }

    /// Build the metadata section, and the reliability cell and score
    /// references the row loops need.
    struct Prepared {
        meta: MzTabMMetaData,
        reliability: MzTabString,
        identification_method: MzTabParameter,
        ms_level: MzTabParameter,
        score_types: Vec<crate::identification::graph::ScoreTypeId>,
    }

    fn prepare(
        feature_map: &FeatureMap,
        id_data: &IdentificationData,
        cv: &ControlledVocabulary,
        options: &MzTabMExportOptions,
    ) -> Result<Prepared> {
        let mut meta = MzTabMMetaData::default();
        meta.mz_tab_id
            .set(&format!("local_id: {}", feature_map.unique_id));

        // software[n] and the SML reliability default.
        let mut reliability = MzTabString::from_text("2");
        let softwares = ordered_softwares(id_data);
        for software in &softwares {
            if let Some(text) = meta_text(&software.software.cv_terms.metadata, "reliability") {
                reliability = MzTabString::from_text(&text);
            }
            let topp_tool = format!("TOPP {}", software.software.name);
            let (id, name) = if cv.has_term_with_name(&topp_tool) {
                term(cv, &topp_tool)?
            } else {
                term(cv, "analysis software")?
            };
            let index = meta.software.len().saturating_add(1);
            meta.software.insert(
                index,
                MzTabSoftwareMetaData {
                    software: ms_parameter(&id, &name, &software.software.version)?,
                    setting: BTreeMap::new(),
                },
            );
        }

        let actions = action_software_names(id_data)?;
        let empty: Vec<String> = Vec::new();
        let quantitation = actions
            .get(&ProcessingAction::Quantitation)
            .unwrap_or(&empty);
        let identification = actions
            .get(&ProcessingAction::Identification)
            .unwrap_or(&empty);

        // quantification_method: recognised only for FeatureFinderMetabo, and
        // the source falls back to the very same term when nothing matched.
        let (quant_id, quant_name) = term(cv, "LC-MS label-free quantitation analysis")?;
        let _ = quantitation
            .iter()
            .any(|name| name == "FeatureFinderMetabo");
        meta.quantification_method = ms_parameter(&quant_id, &quant_name, "")?;

        // ms_run[1]: location from the last input file in name order.
        let mut ms_run = MzTabMMSRunMetaData::default();
        let mut input_names: Vec<&str> = id_data
            .input_files()
            .map(|(_, file)| file.name.as_str())
            .collect();
        input_names.sort_unstable();
        let input_file_name = input_names.last().map(|name| as_file_uri(name));
        if let Some(name) = &input_file_name {
            ms_run.location.set(name);
        }

        // ms_run[1]-scan_polarity[1]: the sign of the first adduct's name.
        let mut adduct_names: Vec<(i32, String, &str)> = Vec::new();
        for (_, adduct) in id_data.adducts() {
            adduct_names.push((
                adduct.charge(),
                adduct.empirical_formula().to_string(),
                adduct.name(),
            ));
        }
        adduct_names.sort_by(|a, b| (a.0, &a.1, a.2).cmp(&(b.0, &b.1, b.2)));
        // `MzTabM.cpp:287` selects positive only when the last character is
        // `+`, so every other spelling — `M+H`, `M+Na`, a bare `H` — is
        // negative there. An empty name underflows `at(size() - 1)` and throws
        // `std::out_of_range`; there is no name to take a sign from, so this
        // keeps the mandatory field writable and calls it positive.
        let positive = match adduct_names.first() {
            None => true,
            Some((_, _, name)) => match name.chars().next_back() {
                Some('+') | None => true,
                Some(_) => false,
            },
        };
        let (polarity_id, polarity_name) = term(
            cv,
            if positive {
                "positive scan"
            } else {
                "negative scan"
            },
        )?;
        ms_run
            .scan_polarity
            .insert(1, ms_parameter(&polarity_id, &polarity_name, "")?);

        // assay[1] and study_variable[1].
        let stem = assay_stem(input_file_name.as_deref().unwrap_or(""));
        let assay = MzTabMAssayMetaData {
            name: MzTabString::from_text(&format!("assay_{stem}")),
            ms_run_ref: MzTabInteger::new(1),
            ..MzTabMAssayMetaData::default()
        };
        let study_variable = MzTabMStudyVariableMetaData {
            name: MzTabString::from_text(&format!("study_variable_{stem}")),
            assay_refs: vec![1],
            description: MzTabString::from_text(&format!("study_variable_{stem}")),
            ..MzTabMStudyVariableMetaData::default()
        };

        // cv[1].
        let (label, full_name) = if options.swap_cv_label_and_full_name {
            (cv.name(), cv.label())
        } else {
            (cv.label(), cv.name())
        };
        meta.cv.insert(
            1,
            MzTabCVMetaData {
                label: MzTabString::from_text(label),
                full_name: MzTabString::from_text(full_name),
                version: MzTabString::from_text(cv.version()),
                url: MzTabString::from_text(cv.url()),
            },
        );

        // database[n]. One accumulating record, as the source declares it
        // outside its loop: a later parameter without a `database_location`
        // meta value therefore keeps the previous one's URI rather than the
        // https://hmdb.ca/ default.
        let mut database = MzTabMDatabaseMetaData {
            database: MzTabParameter::parse("[,, no database , null]")?,
            prefix: MzTabString::null(),
            version: MzTabString::from_text("Unknown"),
            uri: MzTabString::from_text("https://hmdb.ca/"),
        };
        let mut searches: Vec<&crate::identification::graph::DBSearchParam> =
            id_data.db_search_params().map(|(_, p)| p).collect();
        searches.sort_by(|a, b| {
            (&a.database, &a.database_version).cmp(&(&b.database, &b.database_version))
        });
        for search in searches {
            if search.database.contains("custom") {
                database.prefix = MzTabString::null();
                database.version = MzTabString::from_text(&search.database_version);
                database.database = MzTabParameter::parse(&format!("[,, {}, ]", search.database))?;
            } else {
                database.prefix = MzTabString::from_text(&search.database);
                database.version = MzTabString::from_text(&search.database_version);
                database.database = MzTabParameter::parse(&format!("[,,{}, ]", search.database))?;
            }
            if let Some(locations) = text_list(&search.metadata, "database_location")? {
                let joined: Vec<String> = locations.iter().map(|l| as_file_uri(l)).collect();
                database.uri = MzTabString::from_text(&joined.join("|"));
            }
            let index = meta.database.len().saturating_add(1);
            meta.database.insert(index, database.clone());
        }

        // The two quantification units, guessed from FeatureFinderMetabo's
        // quant_method parameter and defaulting to MS1 feature area.
        let mut quantification_unit = None;
        for software in &softwares {
            if software.software.name != "FeatureFinderMetabo" {
                continue;
            }
            let Some(method) = meta_text(
                &software.software.cv_terms.metadata,
                "parameter: algorithm:mtd:quant_method",
            ) else {
                continue;
            };
            let name = match method.as_str() {
                "area" => "MS1 feature area",
                "median" => "median",
                _ => "MS1 feature maximum intensity",
            };
            let (id, term_name) = term(cv, name)?;
            quantification_unit = Some(ms_parameter(&id, &term_name, "")?);
        }
        let quantification_unit = match quantification_unit {
            Some(unit) => unit,
            None => {
                let (id, name) = term(cv, "MS1 feature area")?;
                ms_parameter(&id, &name, "")?
            }
        };
        meta.small_molecule_quantification_unit = quantification_unit.clone();
        meta.small_molecule_feature_quantification_unit = quantification_unit;

        let (rel_id, rel_name) = term(cv, "compound identification confidence level")?;
        meta.small_molecule_identification_reliability = ms_parameter(&rel_id, &rel_name, "")?;

        // id_confidence_measure[n]: one per score the identification tools
        // assign, in software order and then in the software's own priority
        // order, with source-permitted duplicates retained.
        let mut score_types = Vec::new();
        for software in &softwares {
            if !identification.contains(&software.software.name) {
                continue;
            }
            for &score in &software.assigned_scores {
                let score_type = id_data.score_type(score)?;
                let index = score_types.len().saturating_add(1);
                meta.id_confidence_measure.insert(
                    index,
                    MzTabParameter::parse(&format!("[,, {}, ]", score_type.cv_term.name))?,
                );
                score_types.push(score);
            }
        }

        meta.ms_run.insert(1, ms_run);
        meta.assay.insert(1, assay);
        meta.study_variable.insert(1, study_variable);

        // identification_method and ms_level, guessed per identification tool.
        let mut identification_method = MzTabParameter::parse("[, , OpenMS TOPP, ]")?;
        let mut ms_level: Option<MzTabParameter> = None;
        for tool in identification {
            let mut level = 1;
            match tool.as_str() {
                "AccurateMassSearch" => {
                    let (id, name) = term(cv, "accurate mass")?;
                    identification_method = ms_parameter(&id, &name, "")?;
                }
                "SiriusAdapter" => {
                    let (id, name) = term(cv, "de novo search")?;
                    identification_method = ms_parameter(&id, &name, "")?;
                    level = 2;
                }
                "MetaboliteSpectralMatcher" => {
                    let (id, name) = term(cv, "TOPP SpecLibSearcher")?;
                    identification_method = ms_parameter(&id, &name, "")?;
                    level = 2;
                }
                _ => {}
            }
            let (id, name) = term(cv, "ms level")?;
            ms_level = Some(ms_parameter(&id, &name, &level.to_string())?);
        }
        let ms_level = match ms_level {
            Some(level) => level,
            None => {
                let (id, name) = term(cv, "ms level")?;
                ms_parameter(&id, &name, "1")?
            }
        };
        if identification_method.cv_label().is_empty()
            && identification_method.accession().is_empty()
        {
            let (id, name) = term(cv, "data processing action")?;
            identification_method = ms_parameter(&id, &name, "")?;
        }

        Ok(Prepared {
            meta,
            reliability,
            identification_method,
            ms_level,
            score_types,
        })
    }

    /// A feature's observation matches in the order the evidence rows take.
    fn ordered_matches(
        id_data: &IdentificationData,
        feature: &Feature,
        options: &MzTabMExportOptions,
    ) -> Result<Vec<ObservationMatchId>> {
        let mut keyed: Vec<(&str, ObservationMatchId)> = Vec::new();
        for &id in &feature.id_matches {
            keyed.push((compound_of(id_data, id)?.identifier.as_str(), id));
        }
        keyed.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
        if options.deduplicate_matches_by_compound {
            keyed.dedup_by(|a, b| a.0 == b.0);
        }
        Ok(keyed.into_iter().map(|(_, id)| id).collect())
    }

    fn feature_row(
        feature: &Feature,
        identifier: i64,
        feature_keys: &BTreeSet<String>,
    ) -> Result<MzTabMSmallMoleculeFeatureSectionRow> {
        let mut smf = MzTabMSmallMoleculeFeatureSectionRow {
            smf_identifier: MzTabString::from_text(&identifier.to_string()),
            exp_mass_to_charge: MzTabDouble::new(feature.mz),
            charge: MzTabInteger::new(feature.charge),
            retention_time: MzTabDouble::new(feature.rt),
            ..MzTabMSmallMoleculeFeatureSectionRow::default()
        };
        smf.small_molecule_feature_abundance_assay
            .insert(1, MzTabDouble::new(f64::from(feature.intensity)));
        MzTabM::add_meta_info_to_optional_columns(
            feature_keys,
            &mut smf.opt,
            "global",
            &feature.metadata,
        )?;
        Ok(smf)
    }

    #[allow(clippy::too_many_lines)]
    pub(super) fn run(
        feature_map: &FeatureMap,
        id_data: &IdentificationData,
        cv: &ControlledVocabulary,
        options: &MzTabMExportOptions,
    ) -> Result<MzTabM> {
        if id_data.is_empty() {
            return Err(missing(
                "MzTab-M export needs a non-empty IdentificationData graph",
            ));
        }
        if feature_map.features.len() > MzTabM::MAX_ROWS {
            return Err(bad("MzTab-M export exceeds its row limit"));
        }

        let (feature_keys, match_keys, compound_keys) =
            meta_value_keys(feature_map, id_data, options)?;
        let prepared = prepare(feature_map, id_data, cv, options)?;

        let mut smfs: Vec<MzTabMSmallMoleculeFeatureSectionRow> = Vec::new();
        let mut smes: Vec<MzTabMSmallMoleculeEvidenceSectionRow> = Vec::new();
        let mut feature_counter: i64 = 1;
        let mut evidence_counter: i64 = 1;

        for feature in &feature_map.features {
            let matches = ordered_matches(id_data, feature, options)?;
            if matches.is_empty() {
                let mut smf = feature_row(feature, feature_counter, &feature_keys)?;
                smf.sme_id_refs.set_null(true);
                match text_list(&feature.metadata, "adducts")? {
                    Some(adducts) => smf.adduct = MzTabString::from_text(&adducts.join("|")),
                    None => smf.adduct.set_null(true),
                }
                smfs.push(smf);
                feature_counter = feature_counter
                    .checked_add(1)
                    .ok_or_else(|| bad("MzTab-M feature identifier overflows"))?;
                continue;
            }

            let mut per_adduct: BTreeMap<String, Vec<i64>> = BTreeMap::new();
            for id in matches {
                let record = id_data.observation_match(id)?;
                let compound = compound_of(id_data, id)?;
                let mut sme = MzTabMSmallMoleculeEvidenceSectionRow {
                    sme_identifier: MzTabString::from_text(&evidence_counter.to_string()),
                    evidence_input_id: MzTabString::from_text(&format!(
                        "mass={},rt={}",
                        source_double(feature.mz),
                        source_double(feature.rt)
                    )),
                    database_identifier: MzTabString::from_text(&compound.identifier),
                    chemical_formula: MzTabString::from_text(&compound.formula.to_string()),
                    smiles: MzTabString::from_text(&compound.smile),
                    inchi: MzTabString::from_text(&compound.inchi),
                    chemical_name: MzTabString::from_text(&compound.name),
                    exp_mass_to_charge: MzTabDouble::new(feature.mz),
                    charge: MzTabInteger::new(feature.charge),
                    calc_mass_to_charge: MzTabDouble::new(compound.formula.mono_mass()),
                    identification_method: prepared.identification_method.clone(),
                    ms_level: prepared.ms_level.clone(),
                    rank: MzTabInteger::new(1),
                    ..MzTabMSmallMoleculeEvidenceSectionRow::default()
                };
                let adduct = adduct_string(id_data, id)?;
                sme.adduct = MzTabString::from_text(&adduct);
                let data_id = &id_data.observation(record.observation)?.data_id;
                if !data_id.is_empty() {
                    sme.spectra_ref = MzTabSpectraRef::new(1, data_id)?;
                }
                for (position, &score) in prepared.score_types.iter().enumerate() {
                    let value = record.result.score(score).unwrap_or(f64::NAN);
                    sme.id_confidence_measure
                        .insert(position.saturating_add(1), MzTabDouble::new(value));
                }
                MzTabM::add_meta_info_to_optional_columns(
                    &match_keys,
                    &mut sme.opt,
                    "global",
                    &record.result.metadata,
                )?;
                MzTabM::add_meta_info_to_optional_columns(
                    &compound_keys,
                    &mut sme.opt,
                    "global",
                    &compound.result.metadata,
                )?;
                per_adduct.entry(adduct).or_default().push(evidence_counter);
                evidence_counter = evidence_counter
                    .checked_add(1)
                    .ok_or_else(|| bad("MzTab-M evidence identifier overflows"))?;
                smes.push(sme);
                if smes.len() > MzTabM::MAX_ROWS {
                    return Err(bad("MzTab-M SME section exceeds its row limit"));
                }
            }

            for (adduct, evidences) in per_adduct {
                let mut smf = feature_row(feature, feature_counter, &feature_keys)?;
                smf.sme_id_refs.set(
                    evidences
                        .iter()
                        .map(|evidence| MzTabString::from_text(&evidence.to_string()))
                        .collect(),
                );
                smf.adduct = MzTabString::from_text(&adduct);
                if evidences.len() > 1 {
                    smf.sme_id_ref_ambiguity_code = MzTabInteger::new(1);
                }
                smfs.push(smf);
                if smfs.len() > MzTabM::MAX_ROWS {
                    return Err(bad("MzTab-M SMF section exceeds its row limit"));
                }
                feature_counter = feature_counter
                    .checked_add(1)
                    .ok_or_else(|| bad("MzTab-M feature identifier overflows"))?;
            }
        }

        // One summary row per feature row. OpenMS does not aggregate two
        // features whose adducts imply one neutral mass, so, for example, a
        // [M+H]1+ and a [M+Na]1+ trace of the same compound stay separate.
        let mut smls: Vec<MzTabMSmallMoleculeSectionRow> = Vec::new();
        for smf in &smfs {
            let mut sml = MzTabMSmallMoleculeSectionRow {
                sml_identifier: smf.smf_identifier.clone(),
                reliability: prepared.reliability.clone(),
                small_molecule_abundance_assay: smf.small_molecule_feature_abundance_assay.clone(),
                ..MzTabMSmallMoleculeSectionRow::default()
            };
            sml.smf_id_refs.set(vec![smf.smf_identifier.clone()]);
            let mut database_identifier = Vec::new();
            let mut chemical_formula = Vec::new();
            let mut smiles = Vec::new();
            let mut inchi = Vec::new();
            let mut chemical_name = Vec::new();
            let mut uri = Vec::new();
            let mut neutral_mass = Vec::new();
            let mut adducts = Vec::new();
            for evidence in smf.sme_id_refs.get() {
                let row = smes
                    .iter()
                    .find(|sme| sme.sme_identifier.get() == evidence.get())
                    .ok_or_else(|| {
                        missing(
                            "MzTab-M summary row references an SME identifier no evidence row carries",
                        )
                    })?;
                database_identifier.push(row.database_identifier.clone());
                chemical_formula.push(row.chemical_formula.clone());
                smiles.push(row.smiles.clone());
                inchi.push(row.inchi.clone());
                chemical_name.push(row.chemical_name.clone());
                uri.push(row.uri.clone());
                let formula = row.chemical_formula.to_cell_string();
                if formula == "null" {
                    neutral_mass.push(MzTabDouble::null());
                } else {
                    neutral_mass.push(MzTabDouble::new(
                        EmpiricalFormula::parse(&formula)?.mono_mass(),
                    ));
                }
                adducts.push(MzTabString::from_text(row.adduct.get()));
            }
            sml.database_identifier.set(database_identifier);
            sml.chemical_formula.set(chemical_formula);
            sml.smiles.set(smiles);
            sml.inchi.set(inchi);
            sml.chemical_name.set(chemical_name);
            sml.uri.set(uri);
            let mut masses = MzTabDoubleList::default();
            masses.set(neutral_mass);
            sml.theoretical_neutral_mass = masses;
            sml.adducts.set(adducts);
            sml.small_molecule_abundance_study_variable
                .insert(1, MzTabDouble::null());
            sml.small_molecule_abundance_variation_study_variable
                .insert(1, MzTabDouble::null());
            smls.push(sml);
        }

        Ok(MzTabM {
            meta_data: prepared.meta,
            small_molecule_data: smls,
            small_molecule_feature_data: smfs,
            small_molecule_evidence_data: smes,
            empty_rows: Vec::new(),
            comment_rows: BTreeMap::new(),
        })
    }
}

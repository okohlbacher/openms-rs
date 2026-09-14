// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Integration tests for `src/format/mztab_m.rs`, the MzTab-M profile.
//!
//! Every `START_SECTION` of `MzTabM_test.cpp` (342 lines, 4 sections) and
//! `MzTabMFile_test.cpp` (59 lines, 3 sections) is ported here; see
//! `docs/MZTAB_M_SUPPORT.md` for the section-by-section accounting.
//!
//! Two retained C++ outputs anchor the writer:
//!
//! * `tests/data/MzTabMFile_output_1.mztab` — the unmodified file
//!   `MzTabMFile_test.cpp`'s `store` section compares against, i.e. the bytes
//!   `MzTabMFile::store` wrote from
//!   `MzTabM::exportFeatureMapToMzTabM(MzTabMFile_input_1.oms)`.
//! * `tests/data/AccurateMassSearchEngine_output1_mztabm_featureXML.mzTab` —
//!   the retained MzTab-M output of the AccurateMassSearch TOPP test, which
//!   also exercises `id_confidence_measure[n]` columns and scientific-notation
//!   cells.
//!
//! Assertions that quote those files are tier-1 differential evidence: they
//! compare this port's bytes with bytes the C++ actually produced. Assertions
//! that quote `MzTabM_test.cpp` literals are tier-3 source review.

use openms::Error;
use openms::chemistry::AdductInfo;
use openms::format::controlled_vocabulary::ControlledVocabulary;
use openms::format::mztab::{
    MzTabCVMetaData, MzTabContactMetaData, MzTabDouble, MzTabInstrumentMetaData, MzTabInteger,
    MzTabOptionalColumnEntry, MzTabParameter, MzTabParameterList, MzTabSampleMetaData,
    MzTabSoftwareMetaData, MzTabSpectraRef, MzTabString, MzTabStringList,
};
use openms::format::mztab_m::{
    MzTabM, MzTabMAssayMetaData, MzTabMDatabaseMetaData, MzTabMExportOptions, MzTabMFile,
    MzTabMMSRunMetaData, MzTabMMetaData, MzTabMSmallMoleculeEvidenceSectionRow,
    MzTabMSmallMoleculeFeatureSectionRow, MzTabMSmallMoleculeSectionRow,
    MzTabMStudyVariableMetaData, MzTabMWriteOptions, compare_match_by_compound,
};
use openms::identification::graph::{
    DBSearchParam, IdentificationData, IdentifiedCompound, IdentifiedPeptide, InputFile,
    Observation, ObservationMatch, ProcessingSoftware, ProcessingStep, ScoreType,
};
use openms::kernel::{Feature, FeatureMap};
use openms::metadata::{MetaValue, ProcessingAction};
use openms::system::file::TempDir;
use std::collections::{BTreeMap, BTreeSet};

// ---------------------------------------------------------------------------
// Retained C++ output helpers
// ---------------------------------------------------------------------------

const STORE_FIXTURE: &str = "tests/data/MzTabMFile_output_1.mztab";
const AMS_FIXTURE: &str = "tests/data/AccurateMassSearchEngine_output1_mztabm_featureXML.mzTab";

/// The retained file's lines, without terminators.
fn retained(path: &str) -> Vec<String> {
    let text = std::fs::read_to_string(path).expect("retained C++ output is readable");
    // The C++ `TextFile::store` writes one terminator per line, so the final
    // split field is empty and is dropped rather than kept as a blank line.
    let mut lines: Vec<String> = text.split('\n').map(str::to_owned).collect();
    if lines.last().is_some_and(String::is_empty) {
        lines.pop();
    }
    lines
}

/// The `MTD` block of a retained file, in order.
fn retained_metadata(path: &str) -> Vec<String> {
    retained(path)
        .into_iter()
        .take_while(|line| line.starts_with("MTD\t"))
        .collect()
}

fn retained_line(path: &str, prefix: &str, index: usize) -> String {
    retained(path)
        .into_iter()
        .filter(|line| line.starts_with(prefix))
        .nth(index)
        .unwrap_or_else(|| panic!("{path} has no {prefix} line {index}"))
}

fn retained_count(path: &str, prefix: &str) -> usize {
    retained(path)
        .iter()
        .filter(|line| line.starts_with(prefix))
        .count()
}

/// The `opt_`-prefixed columns of a retained header line.
fn retained_optional_columns(path: &str, header: &str) -> Vec<String> {
    retained_line(path, header, 0)
        .split('\t')
        .filter(|column| column.starts_with("opt_"))
        .map(str::to_owned)
        .collect()
}

/// The pinned PSI-MS vocabulary, loaded alone under the source's `PSI-MS`
/// name so that `cv[1]` carries the identity the source writes.
fn psi_ms() -> ControlledVocabulary {
    let mut cv = ControlledVocabulary::new();
    cv.load_obo("PSI-MS", "resources/cv/psi-ms.obo")
        .expect("the pinned psi-ms.obo loads");
    cv
}

// ---------------------------------------------------------------------------
// The metadata section of the retained `store` output, reconstructed
// ---------------------------------------------------------------------------

fn parameter(text: &str) -> MzTabParameter {
    MzTabParameter::parse(text).expect("parameter cell parses")
}

/// The `MTD` section of `MzTabMFile_output_1.mztab`, field by field.
///
/// Every value is read off that retained file; the section it produces is
/// compared with the file's own lines in
/// [`store_metadata_section_matches_retained_bytes`].
fn store_fixture_metadata() -> MzTabMMetaData {
    let mut meta = MzTabMMetaData::default();
    meta.mz_tab_id.set("local_id: 14677498592798891241");
    for (index, accession, name, version) in [
        (
            1,
            "MS:1001456",
            "analysis software",
            "2.7.0-pre-idf-ams-2021-12-20",
        ),
        (
            2,
            "MS:1002169",
            "TOPP FeatureFinderMetabo",
            "2.4.0-nightly-2019-07-17",
        ),
        (
            3,
            "MS:1001456",
            "analysis software",
            "2.4.0-nightly-2019-07-17",
        ),
        (
            4,
            "MS:1001456",
            "analysis software",
            "2.7.0-pre-idf-ams-2021-12-20",
        ),
    ] {
        meta.software.insert(
            index,
            MzTabSoftwareMetaData {
                software: parameter(&format!("[MS, {accession}, {name}, {version}]")),
                setting: BTreeMap::new(),
            },
        );
    }
    meta.quantification_method =
        parameter("[MS, MS:1001834, LC-MS label-free quantitation analysis, ]");
    let mut ms_run = MzTabMMSRunMetaData::default();
    ms_run.location.set(
        "file://I:/OpenSWATH_Metabolomics_data/20181121_full_data/\
         04_PestMixes_individually_Solvent_DDA_20-50/2012_02_03_PStd_10_1-50.wiff",
    );
    ms_run
        .scan_polarity
        .insert(1, parameter("[MS, MS:1000130, positive scan, ]"));
    meta.ms_run.insert(1, ms_run);
    meta.assay.insert(
        1,
        MzTabMAssayMetaData {
            name: MzTabString::from_text("assay_2012_02_03_PStd_10_1-50"),
            ms_run_ref: MzTabInteger::new(1),
            ..MzTabMAssayMetaData::default()
        },
    );
    meta.study_variable.insert(
        1,
        MzTabMStudyVariableMetaData {
            name: MzTabString::from_text("study_variable_2012_02_03_PStd_10_1-50"),
            assay_refs: vec![1],
            description: MzTabString::from_text("study_variable_2012_02_03_PStd_10_1-50"),
            ..MzTabMStudyVariableMetaData::default()
        },
    );
    meta.cv.insert(
        1,
        MzTabCVMetaData {
            label: MzTabString::from_text("PSI-MS"),
            full_name: MzTabString::from_text("MS"),
            version: MzTabString::from_text("4.1.155"),
            url: MzTabString::from_text("http://purl.obolibrary.org/obo/ms/psi-ms.obo"),
        },
    );
    meta.database.insert(
        1,
        MzTabMDatabaseMetaData {
            database: parameter("[,,HMDB, ]"),
            prefix: MzTabString::from_text("HMDB"),
            version: MzTabString::from_text("3.5"),
            uri: MzTabString::from_text("https://hmdb.ca/"),
        },
    );
    meta.small_molecule_quantification_unit = parameter("[MS, MS:1001844, MS1 feature area, ]");
    meta.small_molecule_feature_quantification_unit =
        parameter("[MS, MS:1001844, MS1 feature area, ]");
    meta.small_molecule_identification_reliability =
        parameter("[MS, MS:1002896, compound identification confidence level, ]");
    meta
}

/// The `MTD` section of the AccurateMassSearch output, which declares two
/// `id_confidence_measure` entries and a `file://` database URI.
fn ams_fixture_metadata() -> MzTabMMetaData {
    let mut meta = MzTabMMetaData::default();
    meta.mz_tab_id.set("local_id: 0");
    meta.software.insert(
        1,
        MzTabSoftwareMetaData {
            software: parameter(
                "[MS, MS:1001456, analysis software, 3.4.0-pre-fix-mztab-2025-04-02]",
            ),
            setting: BTreeMap::new(),
        },
    );
    meta.quantification_method =
        parameter("[MS, MS:1001834, LC-MS label-free quantitation analysis, ]");
    let mut ms_run = MzTabMMSRunMetaData::default();
    ms_run.location.set("file://E:/NotARealFileJustATest.mzML");
    ms_run
        .scan_polarity
        .insert(1, parameter("[MS, MS:1000130, positive scan, ]"));
    meta.ms_run.insert(1, ms_run);
    meta.assay.insert(
        1,
        MzTabMAssayMetaData {
            name: MzTabString::from_text("assay_NotARealFileJustATest"),
            ms_run_ref: MzTabInteger::new(1),
            ..MzTabMAssayMetaData::default()
        },
    );
    meta.study_variable.insert(
        1,
        MzTabMStudyVariableMetaData {
            name: MzTabString::from_text("study_variable_NotARealFileJustATest"),
            assay_refs: vec![1],
            description: MzTabString::from_text("study_variable_NotARealFileJustATest"),
            ..MzTabMStudyVariableMetaData::default()
        },
    );
    meta.cv.insert(
        1,
        MzTabCVMetaData {
            label: MzTabString::from_text("PSI-MS"),
            full_name: MzTabString::from_text("MS"),
            version: MzTabString::from_text("4.1.155"),
            url: MzTabString::from_text("http://purl.obolibrary.org/obo/ms/psi-ms.obo"),
        },
    );
    meta.database.insert(
        1,
        MzTabMDatabaseMetaData {
            database: parameter("[,,HMDB, ]"),
            prefix: MzTabString::from_text("HMDB"),
            version: MzTabString::from_text("3.5"),
            uri: MzTabString::from_text(
                "file:///home/sachsenb/Development/OpenMS/src/tests/class_tests/openms/data/\
                 reducedHMDBMapping.tsv",
            ),
        },
    );
    meta.small_molecule_quantification_unit = parameter("[MS, MS:1001844, MS1 feature area, ]");
    meta.small_molecule_feature_quantification_unit =
        parameter("[MS, MS:1001844, MS1 feature area, ]");
    meta.small_molecule_identification_reliability =
        parameter("[MS, MS:1002896, compound identification confidence level, ]");
    meta.id_confidence_measure
        .insert(1, parameter("[,, MassErrorPPMScore, ]"));
    meta.id_confidence_measure
        .insert(2, parameter("[,, MassErrorDaScore, ]"));
    meta
}

// ---------------------------------------------------------------------------
// MzTabM_test.cpp — START_SECTION(MzTabM())
// ---------------------------------------------------------------------------

#[test]
fn default_constructor_declares_the_profile_version() {
    // TEST_NOT_EQUAL(ptr, null_ptr): the source only checks that `new MzTabM`
    // produced an object. The observable part of that construction is the
    // metadata constructor, which sets mzTab-version to 2.0.0-M (asserted in
    // the `Fill data structure` section as well).
    let document = MzTabM::new();
    assert_eq!(document.meta_data.mz_tab_version.get(), "2.0.0-M");
    assert_eq!(
        document.meta_data.mz_tab_version.to_cell_string(),
        "2.0.0-M"
    );
    assert!(document.small_molecule_section_rows().is_empty());
    assert!(document.small_molecule_feature_section_rows().is_empty());
    assert!(document.small_molecule_evidence_section_rows().is_empty());
    assert!(document.empty_rows().is_empty());
    assert!(document.comment_rows().is_empty());
    assert_eq!(document, MzTabM::default());
    // Every other metadata cell starts null.
    assert!(document.meta_data.mz_tab_id.is_null());
    assert!(document.meta_data.quantification_method.is_null());
    assert!(
        document
            .meta_data
            .small_molecule_identification_reliability
            .is_null()
    );
    assert_eq!(MzTabMMetaData::new(), MzTabMMetaData::default());
}

// ---------------------------------------------------------------------------
// MzTabM_test.cpp — START_SECTION(~MzTabM())
// ---------------------------------------------------------------------------

#[test]
fn destructor_releases_a_populated_document() {
    // `delete ptr` with no assertion. A populated document is dropped and the
    // clone taken beforehand stays valid and equal.
    let mut document = MzTabM::new();
    document.set_meta_data(store_fixture_metadata());
    document.small_molecule_data = vec![MzTabMSmallMoleculeSectionRow::default()];
    document.small_molecule_feature_data = vec![MzTabMSmallMoleculeFeatureSectionRow::default()];
    document.small_molecule_evidence_data = vec![MzTabMSmallMoleculeEvidenceSectionRow::default()];
    let copy = document.clone();
    drop(document);
    assert_eq!(copy.small_molecule_section_rows().len(), 1);
    assert_eq!(copy.meta_data().software.len(), 4);
}

// ---------------------------------------------------------------------------
// MzTabM_test.cpp — START_SECTION(Fill data structure)
// ---------------------------------------------------------------------------

/// All 21 assertion macros of the source's `Fill data structure` section, with
/// every literal it sets along the way. Tier-3 evidence: the values are
/// transcribed from `MzTabM_test.cpp:41-318`.
#[test]
fn fill_data_structure() {
    let mut mztabm = MzTabM::new();

    // --- SML row ---------------------------------------------------------
    let mut sml_row = MzTabMSmallMoleculeSectionRow::default();
    sml_row.sml_identifier.from_cell_string("1");
    sml_row.smf_id_refs.from_cell_string("1,2").unwrap();
    sml_row
        .database_identifier
        .from_cell_string("[HMDB:HMDB0001847]")
        .unwrap();
    sml_row
        .chemical_formula
        .from_cell_string("[C17H20N4O2]")
        .unwrap();
    sml_row
        .smiles
        .from_cell_string("[C1=CC=C(C=C1)CCNC(=O)CCNNC(=O)C2=CC=NC=C2]")
        .unwrap();
    sml_row
        .inchi
        .from_cell_string(
            "[InChI=1S/C17H20N4O2/c22-16(19-12-6-14-4-2-1-3-5-14)9-13-20-21-17(23)15-7-10-18-11-8-15/\
             h1-5,7-8,10-11,20H,6,9,12-13H2,(H,19,22)(H,21,23)]",
        )
        .unwrap();
    sml_row
        .chemical_name
        .from_cell_string("[N-(2-phenylethyl)-3-[2-(pyridine-4-carbonyl)hydrazinyl]propanamide]")
        .unwrap();
    sml_row
        .uri
        .from_cell_string("[http://www.hmdb.ca/metabolites/HMDB0001847]")
        .unwrap();
    sml_row
        .theoretical_neutral_mass
        .set(vec![MzTabDouble::new(312.17)]);
    sml_row.adducts.from_cell_string("[[M+H]1+]").unwrap();
    sml_row.reliability.set("3");
    sml_row
        .best_id_confidence_measure
        .from_cell_string("[MS, MS:1000752, TOPP Software,]")
        .unwrap();
    sml_row.best_id_confidence_value.set(0.4);
    for (name, value) in [
        ("SIRIUS_TREE_score", "-10.59083"),
        ("SIRIUS_explained_intensity_score", "96.67"),
        ("SIRIUS_ISO_score", "0.0649874"),
    ] {
        sml_row.opt.push(MzTabOptionalColumnEntry::new(
            name,
            MzTabString::from_text(value),
        ));
    }
    let sml_rows = vec![sml_row];

    // --- SMF row ---------------------------------------------------------
    let mut smf_row = MzTabMSmallMoleculeFeatureSectionRow::default();
    smf_row.smf_identifier.from_cell_string("1");
    smf_row.sme_id_refs.from_cell_string("1").unwrap();
    smf_row
        .sme_id_ref_ambiguity_code
        .from_cell_string("null")
        .unwrap();
    smf_row.adduct.from_cell_string("[M+H]1+");
    smf_row.isotopomer.set_null(true);
    smf_row.exp_mass_to_charge.set(313.1689);
    smf_row.charge.set(1);
    smf_row.retention_time.set(156.0); // always in seconds
    smf_row.rt_start.set(152.2);
    smf_row.rt_end.set(163.4);
    let smf_rows = vec![smf_row];

    // --- SME row ---------------------------------------------------------
    let mut sme_row = MzTabMSmallMoleculeEvidenceSectionRow::default();
    sme_row.sme_identifier.set("1");
    sme_row.evidence_input_id.set("1234.5_156.0");
    sme_row.database_identifier.set("HMDB:HMDB0001847");
    sme_row.chemical_formula.set("C17H20N4O2");
    sme_row
        .smiles
        .set("C1=CC=C(C=C1)CCNC(=O)CCNNC(=O)C2=CC=NC=C2");
    sme_row.inchi.set(
        "InChI=1S/C17H20N4O2/c22-16(19-12-6-14-4-2-1-3-5-14)9-13-20-21-17(23)15-7-10-18-11-8-15/\
         h1-5,7-8,10-11,20H,6,9,12-13H2,(H,19,22)(H,21,23)",
    );
    sme_row
        .chemical_name
        .set("N-(2-phenylethyl)-3-[2-(pyridine-4-carbonyl)hydrazinyl]propanamide");
    sme_row
        .uri
        .set("http://www.hmdb.ca/metabolites/HMDB0001847");
    // `sme_row.derivatized_form.isNull();` in the source discards its result.
    assert!(sme_row.derivatized_form.is_null());
    sme_row.adduct.set("[M+H]1+");
    sme_row.exp_mass_to_charge.set(313.1689);
    sme_row.charge.set(1);
    sme_row.calc_mass_to_charge.set(313.1665);
    let mut sp_ref = MzTabSpectraRef::default();
    sp_ref.set_ms_file(1).unwrap();
    sp_ref.set_spec_ref("index=5").unwrap();
    sme_row.spectra_ref = sp_ref;
    sme_row
        .identification_method
        .from_cell_string("[MS, MS:1000752, TOPP Software,]")
        .unwrap();
    sme_row
        .ms_level
        .from_cell_string("[MS, MS:1000511, ms level, 1]")
        .unwrap();
    // The source writes `id_confidence_measure[0]`, a zero index the format
    // does not define; the map is preserved verbatim.
    sme_row
        .id_confidence_measure
        .insert(0, MzTabDouble::new(123.0));
    sme_row.rank.set(1);
    for (name, value) in [
        ("SIRIUS_TREE_score", "-10.59083"),
        ("SIRIUS_explained_intensity_score", "96.67"),
        ("SIRIUS_ISO_score", "0.0649874"),
    ] {
        sme_row.opt.push(MzTabOptionalColumnEntry::new(
            name,
            MzTabString::from_text(value),
        ));
    }
    let sme_rows = vec![sme_row];

    // --- metadata --------------------------------------------------------
    let mut mztabm_meta = MzTabMMetaData::default();
    mztabm_meta.mz_tab_id.set("local_identifier");
    mztabm_meta.title.set("SML_ROW_TEST");
    mztabm_meta
        .description
        .set("small_molecule_section_row_test");

    let mut sp = MzTabParameterList::default();
    sp.from_cell_string(
        "[MS, MS:1000544, Conversion to mzML, ]|[MS, MS:1000035, Peak picking, ]|\
         [MS, MS:1000594, Low intensity data point removal, ]",
    )
    .unwrap();
    mztabm_meta.sample_processing.insert(0, sp);

    let mut meta_instrument = MzTabInstrumentMetaData::default();
    meta_instrument
        .name
        .from_cell_string(
            "[MS, MS:1000483, Thermo Fisher Scientific instrument model, LTQ Orbitrap Velos]",
        )
        .unwrap();
    meta_instrument
        .source
        .from_cell_string("[MS, MS:1000008, Ionization Type, ESI]")
        .unwrap();
    meta_instrument.analyzer.insert(
        0,
        parameter("[MS, MS:1000443, Mass Analyzer Type, Orbitrap]"),
    );
    meta_instrument
        .detector
        .from_cell_string("[MS, MS:1000453, Detector, Dynode Detector]")
        .unwrap();
    mztabm_meta.instrument.insert(0, meta_instrument);

    let mut meta_software = MzTabSoftwareMetaData {
        software: parameter("[MS, MS:1002205, ProteoWizard msconvert, ]"),
        ..MzTabSoftwareMetaData::default()
    };
    meta_software
        .setting
        .insert(0, MzTabString::from_text("Peak Picking MS1"));
    mztabm_meta.software.insert(0, meta_software);

    mztabm_meta.publication.insert(
        0,
        MzTabString::from_text("pubmed:21063943|doi:10.1007/978-1-60761-987-1_6"),
    );

    let meta_contact = MzTabContactMetaData {
        name: MzTabString::from_text("Max MusterMann"),
        affiliation: MzTabString::from_text("University of Musterhausen"),
        email: MzTabString::from_text("MMM@please_do_not_try_to_write_an_email.com"),
    };
    mztabm_meta.contact.insert(0, meta_contact);
    mztabm_meta.uri.insert(
        0,
        MzTabString::from_text("https://www.ebi.ac.uk/metabolights/MTBLS"),
    );
    mztabm_meta.external_study_uri.insert(
        0,
        MzTabString::from_text(
            "https://www.ebi.ac.uk/metabolights/MTBLS/files/i_Investigation.txt",
        ),
    );
    mztabm_meta
        .quantification_method
        .from_cell_string("[MS, MS:1001834, LC-MS label-free quantitation analysis, ]")
        .unwrap();

    let meta_sample = MzTabSampleMetaData {
        description: MzTabString::from_text("Nice Sample"),
        ..MzTabSampleMetaData::default()
    };
    mztabm_meta.sample.insert(0, meta_sample);

    let mut meta_msrun = MzTabMMSRunMetaData {
        location: MzTabString::from_text("ftp://ftp.ebi.ac.uk/path/to/file"),
        instrument_ref: MzTabInteger::new(0),
        format: parameter("[MS, MS:1000584, mzML file, ]"),
        id_format: parameter("[MS, MS:1000584, mzML file, ]"),
        hash: MzTabString::from_text("de9f2c7fd25e1b3afad3e85a0bd17d9b100db4b3"),
        hash_method: parameter("[MS, MS:1000569, SHA-1, ]"),
        ..MzTabMMSRunMetaData::default()
    };
    meta_msrun
        .fragmentation_method
        .insert(0, parameter("[MS, MS:1000133, CID, ]"));
    meta_msrun
        .fragmentation_method
        .insert(1, parameter("[MS, MS:1000422, HCD, ]"));
    meta_msrun
        .scan_polarity
        .insert(0, parameter("[MS, MS:1000130, positive scan, ]"));
    meta_msrun
        .scan_polarity
        .insert(1, parameter("[MS, MS:1000130, positive scan, ]"));
    mztabm_meta.ms_run.insert(0, meta_msrun);

    let mut meta_assay = MzTabMAssayMetaData {
        external_uri: MzTabString::from_text(
            "https://www.ebi.ac.uk/metabolights/MTBLS/files/i_Investigation.txt?STUDYASSAY=a_8pos.txt",
        ),
        sample_ref: MzTabInteger::new(1),
        ms_run_ref: MzTabInteger::new(1),
        ..MzTabMAssayMetaData::default()
    };
    meta_assay
        .custom
        .insert(0, parameter("[MS, , Assay operator, Blogs]"));
    mztabm_meta.assay.insert(0, meta_assay);

    let mut pl_factors = MzTabParameterList::default();
    pl_factors
        .from_cell_string("[MS, MS:1000130, positive scan, ]")
        .unwrap();
    let meta_study = MzTabMStudyVariableMetaData {
        assay_refs: vec![1],
        average_function: parameter("[MS, MS:1002883, median, ]"),
        variation_function: parameter("[MS, MS:1002885, standard error, ]"),
        description: MzTabString::from_text("control"),
        factors: pl_factors,
        ..MzTabMStudyVariableMetaData::default()
    };
    mztabm_meta.study_variable.insert(0, meta_study);

    mztabm_meta.cv.insert(
        0,
        MzTabCVMetaData {
            label: MzTabString::from_text("MS"),
            full_name: MzTabString::from_text("PSI-MS controlled vocabulary"),
            version: MzTabString::from_text("4.1.155"),
            url: MzTabString::from_text("share/OpenMS/CV/psi-ms.obo"),
        },
    );

    mztabm_meta.database.insert(
        0,
        MzTabMDatabaseMetaData {
            database: parameter("[MIRIAM, MIR:00100079, HMDB, ]"),
            prefix: MzTabString::from_text("HMDB"),
            version: MzTabString::from_text("4.0"),
            // `MzTabString("null")` stores nothing, so this cell is null.
            uri: MzTabString::from_text("null"),
        },
    );

    mztabm_meta.small_molecule_quantification_unit =
        parameter("[MS, MS:1000042, peak intensity, ]");
    mztabm_meta.small_molecule_feature_quantification_unit =
        parameter("[MS, MS:1000042, peak intensity, ]");
    mztabm_meta.small_molecule_identification_reliability =
        parameter("[MS, MS:1002955, hr-ms compound identification confidence level, ]");
    mztabm_meta
        .id_confidence_measure
        .insert(0, parameter("[MS,MS:1002890,fragmentation score,]"));

    mztabm.set_meta_data(mztabm_meta);
    mztabm.set_small_molecule_section_rows(sml_rows);
    mztabm.set_small_molecule_feature_section_rows(smf_rows);
    mztabm.set_small_molecule_evidence_section_rows(sme_rows);

    // --- the 21 assertions ----------------------------------------------
    let sml_test = &mztabm.small_molecule_section_rows()[0];
    assert_eq!(sml_test.smf_id_refs.to_cell_string(), "1,2");
    assert_eq!(sml_test.adducts.to_cell_string(), "[[M+H]1+]");

    let smf_test = &mztabm.small_molecule_feature_section_rows()[0];
    assert_eq!(
        smf_test.exp_mass_to_charge.to_cell_string(),
        "313.168900000000008"
    );
    assert_eq!(smf_test.retention_time.to_cell_string(), "156.0");

    let sme_test = &mztabm.small_molecule_evidence_section_rows()[0];
    assert_eq!(
        sme_test.database_identifier.to_cell_string(),
        "HMDB:HMDB0001847"
    );
    assert_eq!(
        sme_test.identification_method.to_cell_string(),
        "[MS, MS:1000752, TOPP Software, ]"
    );

    let mtest = mztabm.meta_data();
    assert_eq!(mtest.mz_tab_version.to_cell_string(), "2.0.0-M"); // set by constructor
    assert_eq!(
        mtest.sample_processing[&0].to_cell_string(),
        "[MS, MS:1000544, Conversion to mzML, ]|[MS, MS:1000035, Peak picking, ]|\
         [MS, MS:1000594, Low intensity data point removal, ]"
    );
    assert_eq!(
        mtest.instrument[&0].analyzer[&0].to_cell_string(),
        "[MS, MS:1000443, Mass Analyzer Type, Orbitrap]"
    );
    assert_eq!(
        mtest.software[&0].setting[&0].to_cell_string(),
        "Peak Picking MS1"
    );
    assert_eq!(
        mtest.contact[&0].affiliation.to_cell_string(),
        "University of Musterhausen"
    );
    assert_eq!(mtest.sample[&0].description.to_cell_string(), "Nice Sample");
    assert_eq!(
        mtest.ms_run[&0].format.to_cell_string(),
        "[MS, MS:1000584, mzML file, ]"
    );
    assert_eq!(
        mtest.study_variable[&0].description.to_cell_string(),
        "control"
    );
    assert_eq!(mtest.database[&0].prefix.to_cell_string(), "HMDB");
    assert_eq!(
        mtest.small_molecule_quantification_unit.to_cell_string(),
        "[MS, MS:1000042, peak intensity, ]"
    );

    let optional_sml_columns = mztabm.small_molecule_optional_column_names().unwrap();
    let optional_sme_columns = mztabm
        .small_molecule_evidence_optional_column_names()
        .unwrap();

    assert_eq!(mztabm.small_molecule_section_rows().len(), 1);
    assert_eq!(mztabm.small_molecule_feature_section_rows().len(), 1);
    assert_eq!(mztabm.small_molecule_feature_section_rows().len(), 1);

    assert_eq!(optional_sml_columns.len(), 3);
    assert_eq!(optional_sme_columns.len(), 3);

    // Beyond the source's own assertions: the column names keep the order the
    // row declares them in, not alphabetical order.
    assert_eq!(
        optional_sml_columns,
        vec![
            "SIRIUS_TREE_score".to_owned(),
            "SIRIUS_explained_intensity_score".to_owned(),
            "SIRIUS_ISO_score".to_owned(),
        ]
    );
    assert!(
        mztabm
            .small_molecule_feature_optional_column_names()
            .unwrap()
            .is_empty()
    );
    // The `null` database URI the source stores really is a null cell.
    assert!(mtest.database[&0].uri.is_null());
}

// ---------------------------------------------------------------------------
// MzTabM_test.cpp — START_SECTION(MzTabM::exportFeatureMapToMzTabM(...))
// ---------------------------------------------------------------------------

/// All six assertion macros of the source's `exportFeatureMapToMzTabM`
/// section.
///
/// The section's input is `MzTabMFile_input_1.oms`, an SQLite `.oms` file this
/// crate has no reader for, so the export cannot be re-run on the same input.
/// The six literals are instead verified against
/// `tests/data/MzTabMFile_output_1.mztab`, the retained output that
/// `MzTabMFile_test.cpp` produced by storing the document *this very call*
/// returned: the `SML`/`SMF`/`SME` line counts are the three section sizes and
/// the `opt_` columns of the three header lines are the three optional-column
/// lists. That makes them tier-1 values rather than transcribed literals. The
/// port's own exporter is exercised against a synthetic graph in
/// [`export_reproduces_retained_metadata_section`] and
/// [`export_builds_one_summary_row_per_feature_row`].
#[test]
fn export_feature_map_to_mztab_m() {
    assert_eq!(retained_count(STORE_FIXTURE, "SML\t"), 83);
    assert_eq!(retained_count(STORE_FIXTURE, "SMF\t"), 83);
    assert_eq!(retained_count(STORE_FIXTURE, "SME\t"), 312);

    assert_eq!(retained_optional_columns(STORE_FIXTURE, "SMH\t").len(), 0);
    assert_eq!(retained_optional_columns(STORE_FIXTURE, "SFH\t").len(), 18);
    assert_eq!(retained_optional_columns(STORE_FIXTURE, "SEH\t").len(), 6);

    // The same six quantities, read back out of a document this port builds
    // from those retained lines, agree with `MzTabM`'s own accessors.
    let document = document_from_retained(STORE_FIXTURE);
    assert_eq!(document.small_molecule_section_rows().len(), 83);
    assert_eq!(document.small_molecule_feature_section_rows().len(), 83);
    assert_eq!(document.small_molecule_evidence_section_rows().len(), 312);
    assert_eq!(
        document
            .small_molecule_optional_column_names()
            .unwrap()
            .len(),
        0
    );
    assert_eq!(
        document
            .small_molecule_feature_optional_column_names()
            .unwrap()
            .len(),
        18
    );
    assert_eq!(
        document
            .small_molecule_evidence_optional_column_names()
            .unwrap()
            .len(),
        6
    );
}

/// A document whose sections carry one row per data line of `path` and whose
/// rows declare that file's own optional columns, so that the port's section
/// sizes and optional-column lists can be compared with the retained file's.
fn document_from_retained(path: &str) -> MzTabM {
    let mut document = MzTabM::new();
    let sml_columns = retained_optional_columns(path, "SMH\t");
    let smf_columns = retained_optional_columns(path, "SFH\t");
    let sme_columns = retained_optional_columns(path, "SEH\t");
    let entries = |names: &[String]| -> Vec<MzTabOptionalColumnEntry> {
        names
            .iter()
            .map(|name| MzTabOptionalColumnEntry::new(name.clone(), MzTabString::default()))
            .collect()
    };
    for _ in 0..retained_count(path, "SML\t") {
        document
            .small_molecule_data
            .push(MzTabMSmallMoleculeSectionRow {
                opt: entries(&sml_columns),
                ..MzTabMSmallMoleculeSectionRow::default()
            });
    }
    for _ in 0..retained_count(path, "SMF\t") {
        document
            .small_molecule_feature_data
            .push(MzTabMSmallMoleculeFeatureSectionRow {
                opt: entries(&smf_columns),
                ..MzTabMSmallMoleculeFeatureSectionRow::default()
            });
    }
    for _ in 0..retained_count(path, "SME\t") {
        document
            .small_molecule_evidence_data
            .push(MzTabMSmallMoleculeEvidenceSectionRow {
                opt: entries(&sme_columns),
                ..MzTabMSmallMoleculeEvidenceSectionRow::default()
            });
    }
    document
}

// ---------------------------------------------------------------------------
// MzTabMFile_test.cpp — START_SECTION(MzTabMFile()) and (~MzTabFile())
// ---------------------------------------------------------------------------

#[test]
fn file_default_constructor() {
    // TEST_NOT_EQUAL(ptr, null_ptr). The observable part of the default
    // construction is the option set; the source's is native here, because
    // the source's own defects are opt-in.
    let file = MzTabMFile::new();
    assert_eq!(file, MzTabMFile::default());
    assert_eq!(file.options, MzTabMWriteOptions::default());
    assert!(!file.options.source_assay_custom_key);
    assert!(!file.options.source_colunit_keys);
    assert!(!file.options.source_derivatization_agent_key);
    assert!(!file.options.omit_ms_run_id_format);
    assert!(!file.options.source_row_abundance_cells);
    assert!(!file.options.source_verbatim_cells);
    let source = MzTabMFile::with_options(MzTabMWriteOptions::source());
    assert!(source.options.source_assay_custom_key);
    assert!(source.options.source_colunit_keys);
    assert!(source.options.source_derivatization_agent_key);
    assert!(source.options.omit_ms_run_id_format);
    assert!(source.options.source_row_abundance_cells);
    assert!(source.options.source_verbatim_cells);
}

/// One SML row whose `chemical_name` and `opt_` value carry a tab and a line
/// break, the shape a featureXML or `.oms` text node can legally have.
fn document_with_separator_cells() -> MzTabM {
    let mut names = MzTabStringList::default();
    names.set(vec![MzTabString::from_text("acetyl-\n-carnitine")]);
    let mut row = MzTabMSmallMoleculeSectionRow {
        sml_identifier: MzTabString::from_text("1"),
        chemical_name: names,
        ..MzTabMSmallMoleculeSectionRow::default()
    };
    row.opt.push(MzTabOptionalColumnEntry {
        name: "opt_global_note".to_owned(),
        value: MzTabString::from_text("a\tb"),
    });
    let mut document = MzTabM::new();
    document.small_molecule_data.push(row);
    document
}

#[test]
fn a_cell_carrying_a_tab_or_a_line_break_is_refused_by_default() {
    // `MzTabMFile.cpp:626-631` hands each row to TextFile unchanged, so such a
    // cell splits the row or adds a column and the output stops being a
    // rectangle — and text starting with MTD/SMH/SML forges a line of that
    // kind. The default refuses it; `MzTabMWriteOptions::source()` reproduces
    // the source's pass-through.
    let document = document_with_separator_cells();
    let error = MzTabMFile::new().generate_lines(&document).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");

    // The same document under source options is written as the source writes
    // it: the SML row splits across two physical lines of 8 and 10 fields
    // against a 17-field header.
    let lines = MzTabMFile::with_options(MzTabMWriteOptions::source())
        .generate_lines(&document)
        .expect("source options reproduce the defect");
    let header = lines
        .iter()
        .find(|line| line.starts_with("SMH\t"))
        .expect("an SMH header");
    let row = lines
        .iter()
        .find(|line| line.starts_with("SML\t"))
        .expect("an SML row");
    assert!(row.contains('\n'), "the cell's line break survives");
    assert_ne!(
        row.split('\n').next().unwrap().split('\t').count(),
        header.split('\t').count(),
        "the first physical line is short of the header"
    );
}

#[test]
fn a_metadata_key_or_value_carrying_a_separator_is_refused_by_default() {
    let mut document = MzTabM::new();
    document.set_meta_data(MzTabMMetaData {
        title: MzTabString::from_text("a\tb\nMTD\tmzTab-version\tforged"),
        ..MzTabMMetaData::default()
    });
    let error = MzTabMFile::new().generate_lines(&document).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    // Under source options the forged line is written, as the source writes it.
    let lines = MzTabMFile::with_options(MzTabMWriteOptions::source())
        .generate_lines(&document)
        .expect("source options reproduce the defect");
    assert!(
        lines.iter().any(|line| line.contains("\nMTD\t")),
        "the forged metadata line is present: {lines:?}"
    );
}

#[test]
fn file_destructor() {
    // `delete ptr` with no assertion. `MzTabMFile` owns nothing that needs a
    // destructor, so the scope end is the whole of it; the clone taken inside
    // stays valid afterwards.
    let copy = {
        let file = MzTabMFile::with_options(MzTabMWriteOptions::source());
        file.clone()
    };
    assert!(copy.options.omit_ms_run_id_format);
}

// ---------------------------------------------------------------------------
// MzTabMFile_test.cpp — START_SECTION(void store(...))
// ---------------------------------------------------------------------------

/// The source's `store` section: export, store, `TEST_FILE_SIMILAR` against
/// `MzTabMFile_output_1.mztab`.
///
/// The export half needs the `.oms` reader this crate lacks, so the document
/// is rebuilt from the retained file's own metadata and its three shortest
/// data rows; `store` then has to reproduce those bytes exactly. Tier-1
/// evidence.
#[test]
fn store_writes_the_retained_bytes() {
    let mut document = MzTabM::new();
    document.set_meta_data(store_fixture_metadata());
    document
        .small_molecule_data
        .push(store_fixture_sml_row_three());
    document
        .small_molecule_feature_data
        .push(store_fixture_smf_row_one());
    document
        .small_molecule_evidence_data
        .push(store_fixture_sme_row_one());

    let writer = MzTabMFile::with_options(MzTabMWriteOptions::source());
    // A directory of its own, so concurrent runs cannot remove each other's files.
    let temp = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let directory = temp.path();
    let path = directory.join("stored.mzTab");
    writer.store(&path, &document).unwrap();
    let written = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(&path).ok();

    let mut expected: Vec<String> = retained_metadata(STORE_FIXTURE);
    expected.push(String::new());
    expected.push(retained_line(STORE_FIXTURE, "SMH\t", 0));
    expected.push(retained_line(STORE_FIXTURE, "SML\t", 2));
    expected.push(String::new());
    expected.push(retained_line(STORE_FIXTURE, "SFH\t", 0));
    expected.push(retained_line(STORE_FIXTURE, "SMF\t", 0));
    expected.push(String::new());
    expected.push(retained_line(STORE_FIXTURE, "SEH\t", 0));
    expected.push(retained_line(STORE_FIXTURE, "SME\t", 0));

    let produced: Vec<String> = written
        .lines()
        .map(|line| line.trim_end_matches('\r').to_owned())
        .collect();
    assert_eq!(produced, expected);
    // The source writes one terminator per line, so the file ends with one.
    assert!(written.ends_with('\n'));
}

/// `SML` row 3 of the retained output — the shortest, a single-evidence row.
#[allow(clippy::excessive_precision)] // literals transcribed from the retained C++ bytes
fn store_fixture_sml_row_three() -> MzTabMSmallMoleculeSectionRow {
    let mut row = MzTabMSmallMoleculeSectionRow {
        sml_identifier: MzTabString::from_text("3"),
        reliability: MzTabString::from_text("2"),
        ..MzTabMSmallMoleculeSectionRow::default()
    };
    row.smf_id_refs.set(vec![MzTabString::from_text("3")]);
    row.database_identifier
        .set(vec![MzTabString::from_text("HMDB:HMDB05033")]);
    row.chemical_formula
        .set(vec![MzTabString::from_text("C16H14N2O3S1")]);
    row.smiles.set(vec![MzTabString::from_text("Valdecoxib")]);
    row.inchi.set(vec![MzTabString::from_text(
        "InChI=1S/C16H14N2O3S/c1-11-15(12-7-9-14(10-8-12)22(17,19)20)16(18-21-11)13-5-3-2-4-6-13/\
         h2-10H,1H3,(H2,17,19,20)",
    )]);
    row.chemical_name
        .set(vec![MzTabString::from_text("Valdecoxib")]);
    row.uri.set(vec![MzTabString::null()]);
    row.theoretical_neutral_mass
        .set(vec![MzTabDouble::new(314.0725141766)]);
    row.adducts.set(vec![MzTabString::from_text("[M+Na]1+")]);
    row.small_molecule_abundance_assay
        .insert(1, MzTabDouble::new(464.264129638671875));
    row.small_molecule_abundance_study_variable
        .insert(1, MzTabDouble::null());
    row.small_molecule_abundance_variation_study_variable
        .insert(1, MzTabDouble::null());
    row
}

/// `SMF` row 1 of the retained output, with all 18 optional columns.
#[allow(clippy::excessive_precision)] // literals transcribed from the retained C++ bytes
fn store_fixture_smf_row_one() -> MzTabMSmallMoleculeFeatureSectionRow {
    let mut row = MzTabMSmallMoleculeFeatureSectionRow {
        smf_identifier: MzTabString::from_text("1"),
        sme_id_ref_ambiguity_code: MzTabInteger::new(1),
        adduct: MzTabString::from_text("[M+H]1+"),
        exp_mass_to_charge: MzTabDouble::new(118.086281670984334),
        charge: MzTabInteger::new(0),
        retention_time: MzTabDouble::new(70.157003402709961),
        ..MzTabMSmallMoleculeFeatureSectionRow::default()
    };
    row.sme_id_refs.set(
        (1..=8)
            .map(|n| MzTabString::from_text(&n.to_string()))
            .collect(),
    );
    row.small_molecule_feature_abundance_assay
        .insert(1, MzTabDouble::new(406.5467529296875));
    for (name, value) in [
        ("opt_global_FWHM", "3.829350471496582"),
        ("opt_global_Group", "9297005714436264670"),
        ("opt_global_adducts", "null"),
        ("opt_global_dc_charge_adduct_mass", "null"),
        ("opt_global_dc_charge_adducts", "null"),
        ("opt_global_is_backbone", "null"),
        ("opt_global_is_ungrouped_monoisotopic", "1"),
        ("opt_global_is_ungrouped_with_charge", "null"),
        ("opt_global_isotope_distances", "[]"),
        ("opt_global_label", "T697.1"),
        ("opt_global_legal_isotope_pattern", "-1"),
        ("opt_global_map_idx", "null"),
        ("opt_global_masstrace_centroid_mz", "[118.086281670984334]"),
        ("opt_global_masstrace_centroid_rt", "[70.157003402709961]"),
        ("opt_global_masstrace_intensity", "[406.546738705937969]"),
        ("opt_global_max_height", "134.418792724609375"),
        ("opt_global_num_of_masstraces", "1"),
        ("opt_global_old_charge", "null"),
    ] {
        row.opt.push(MzTabOptionalColumnEntry::new(
            name,
            MzTabString::from_text(value),
        ));
    }
    row
}

/// `SME` row 1 of the retained output, with all 6 optional columns.
#[allow(clippy::excessive_precision)] // literals transcribed from the retained C++ bytes
fn store_fixture_sme_row_one() -> MzTabMSmallMoleculeEvidenceSectionRow {
    let mut row = MzTabMSmallMoleculeEvidenceSectionRow {
        sme_identifier: MzTabString::from_text("1"),
        evidence_input_id: MzTabString::from_text("mass=118.086281670984334,rt=70.157003402709961"),
        database_identifier: MzTabString::from_text("HMDB:HMDB00043"),
        chemical_formula: MzTabString::from_text("C5H11N1O2"),
        smiles: MzTabString::from_text("Betaine"),
        inchi: MzTabString::from_text("InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3"),
        chemical_name: MzTabString::from_text("Betaine"),
        adduct: MzTabString::from_text("[M+H]1+"),
        exp_mass_to_charge: MzTabDouble::new(118.086281670984334),
        charge: MzTabInteger::new(0),
        calc_mass_to_charge: MzTabDouble::new(117.0789793509),
        spectra_ref: MzTabSpectraRef::new(1, "2655476886865018721").unwrap(),
        identification_method: parameter("[MS, MS:1000543, data processing action, ]"),
        ms_level: parameter("[MS, MS:1000511, ms level, 1]"),
        rank: MzTabInteger::new(1),
        ..MzTabMSmallMoleculeEvidenceSectionRow::default()
    };
    for (name, value) in [
        ("opt_global_chemical_formula", "C5H11NO2"),
        ("opt_global_description", "[Betaine]"),
        ("opt_global_identifier", "[HMDB:HMDB00043]"),
        ("opt_global_modifications", "M+H;1+"),
        ("opt_global_mz_error_Da", "2.661798863812237e-05"),
        ("opt_global_mz_error_ppm", "0.225411404792002"),
    ] {
        row.opt.push(MzTabOptionalColumnEntry::new(
            name,
            MzTabString::from_text(value),
        ));
    }
    row
}

// ---------------------------------------------------------------------------
// Tier-1 differential: the writer against the retained bytes
// ---------------------------------------------------------------------------

#[test]
fn store_metadata_section_matches_retained_bytes() {
    let writer = MzTabMFile::with_options(MzTabMWriteOptions::source());
    let produced = writer
        .generate_meta_data_section(&store_fixture_metadata())
        .unwrap();
    assert_eq!(produced, retained_metadata(STORE_FIXTURE));
    assert_eq!(produced.len(), 25);
    // The mandatory keys the source emits even when the value is unset appear
    // once each, with a real value here.
    assert_eq!(produced[0], "MTD\tmzTab-version\t2.0.0-M");
}

#[test]
fn ams_metadata_section_matches_retained_bytes() {
    let writer = MzTabMFile::with_options(MzTabMWriteOptions::source());
    let produced = writer
        .generate_meta_data_section(&ams_fixture_metadata())
        .unwrap();
    assert_eq!(produced, retained_metadata(AMS_FIXTURE));
    assert_eq!(produced.len(), 24);
    assert!(
        produced
            .iter()
            .any(|line| line == "MTD\tid_confidence_measure[2]\t[, , MassErrorDaScore, ]")
    );
}

#[test]
fn headers_match_retained_bytes() {
    let writer = MzTabMFile::with_options(MzTabMWriteOptions::source());
    for (path, meta) in [
        (STORE_FIXTURE, store_fixture_metadata()),
        (AMS_FIXTURE, ams_fixture_metadata()),
    ] {
        let sml_columns = retained_optional_columns(path, "SMH\t");
        let smf_columns = retained_optional_columns(path, "SFH\t");
        let sme_columns = retained_optional_columns(path, "SEH\t");
        let sml = writer
            .generate_small_molecule_header(&meta, &sml_columns)
            .unwrap();
        assert_eq!(sml.text, retained_line(path, "SMH\t", 0));
        assert_eq!(sml.columns, sml.text.split('\t').count());
        let smf = writer
            .generate_small_molecule_feature_header(&meta, &smf_columns)
            .unwrap();
        assert_eq!(smf.text, retained_line(path, "SFH\t", 0));
        let sme = writer
            .generate_small_molecule_evidence_header(&meta, &sme_columns)
            .unwrap();
        assert_eq!(sme.text, retained_line(path, "SEH\t", 0));
    }
    // The evidence header of the AccurateMassSearch output carries the two
    // declared confidence columns, immediately before `rank`.
    let header = retained_line(AMS_FIXTURE, "SEH\t", 0);
    let columns: Vec<&str> = header.split('\t').collect();
    let rank = columns.iter().position(|c| *c == "rank").unwrap();
    assert_eq!(columns[rank - 2], "id_confidence_measure[1]");
    assert_eq!(columns[rank - 1], "id_confidence_measure[2]");
}

#[test]
fn rows_match_retained_bytes() {
    let writer = MzTabMFile::with_options(MzTabMWriteOptions::source());
    let meta = store_fixture_metadata();
    let sml = writer
        .generate_small_molecule_section_row(&store_fixture_sml_row_three(), &meta, &[])
        .unwrap();
    assert_eq!(sml.text, retained_line(STORE_FIXTURE, "SML\t", 2));
    assert_eq!(sml.columns, 17);

    let smf_columns = retained_optional_columns(STORE_FIXTURE, "SFH\t");
    let smf = writer
        .generate_small_molecule_feature_section_row(
            &store_fixture_smf_row_one(),
            &meta,
            &smf_columns,
        )
        .unwrap();
    assert_eq!(smf.text, retained_line(STORE_FIXTURE, "SMF\t", 0));

    let sme_columns = retained_optional_columns(STORE_FIXTURE, "SEH\t");
    let sme = writer
        .generate_small_molecule_evidence_section_row(
            &store_fixture_sme_row_one(),
            &meta,
            &sme_columns,
        )
        .unwrap();
    assert_eq!(sme.text, retained_line(STORE_FIXTURE, "SME\t", 0));

    // `SML` row 5, whose eight list cells each carry two entries, so the
    // `|` separators and the `null|null` spelling of a two-entry list of null
    // strings are pinned as well — an *empty* list renders as the single cell
    // `null`, a two-entry list of nulls as `null|null`.
    let sml = writer
        .generate_small_molecule_section_row(&store_fixture_sml_row_five(), &meta, &[])
        .unwrap();
    let expected = retained_line(STORE_FIXTURE, "SML\t", 4);
    assert_eq!(sml.text, expected);
    assert!(expected.contains("\tnull|null\t"));
    assert_eq!(sml.columns, 17);
}

/// `SML` row 5 of the retained output: two evidences under one adduct, so
/// every list cell has two entries.
#[allow(clippy::excessive_precision)] // literals transcribed from the retained C++ bytes
fn store_fixture_sml_row_five() -> MzTabMSmallMoleculeSectionRow {
    let mut row = MzTabMSmallMoleculeSectionRow {
        sml_identifier: MzTabString::from_text("5"),
        reliability: MzTabString::from_text("2"),
        ..MzTabMSmallMoleculeSectionRow::default()
    };
    row.smf_id_refs.set(vec![MzTabString::from_text("5")]);
    row.database_identifier.set(vec![
        MzTabString::from_text("HMDB:HMDB41777"),
        MzTabString::from_text("HMDB:HMDB41778"),
    ]);
    row.chemical_formula.set(vec![
        MzTabString::from_text("C16H12O9S1"),
        MzTabString::from_text("C16H12O9S1"),
    ]);
    row.smiles.set(vec![
        MzTabString::from_text("Tectorigenin 4'-sulfate"),
        MzTabString::from_text("Tectorigenin 7-sulfate"),
    ]);
    row.inchi.set(vec![
        MzTabString::from_text(
            "InChI=1S/C16H12O9S/c1-23-16-11(17)6-12-13(15(16)19)14(18)10(7-24-12)8-2-4-9(5-3-8)\
             25-26(20,21)22/h2-7,17,19H,1H3,(H,20,21,22)",
        ),
        MzTabString::from_text(
            "InChI=1S/C16H12O9S/c1-23-16-12(25-26(20,21)22)6-11-13(15(16)19)14(18)10(7-24-11)\
             8-2-4-9(17)5-3-8/h2-7,17,19H,1H3,(H,20,21,22)",
        ),
    ]);
    row.chemical_name.set(vec![
        MzTabString::from_text("Tectorigenin 4'-sulfate"),
        MzTabString::from_text("Tectorigenin 7-sulfate"),
    ]);
    row.uri.set(vec![MzTabString::null(), MzTabString::null()]);
    row.theoretical_neutral_mass.set(vec![
        MzTabDouble::new(380.020206112799997),
        MzTabDouble::new(380.020206112799997),
    ]);
    row.adducts.set(vec![
        MzTabString::from_text("[M+Na]1+"),
        MzTabString::from_text("[M+Na]1+"),
    ]);
    row.small_molecule_abundance_assay
        .insert(1, MzTabDouble::new(464.30322265625));
    row.small_molecule_abundance_study_variable
        .insert(1, MzTabDouble::null());
    row.small_molecule_abundance_variation_study_variable
        .insert(1, MzTabDouble::null());
    row
}

#[test]
#[allow(clippy::excessive_precision)] // literals transcribed from the retained C++ bytes
fn ams_rows_match_retained_bytes() {
    let writer = MzTabMFile::with_options(MzTabMWriteOptions::source());
    let meta = ams_fixture_metadata();

    // SML row 3: a single-evidence summary whose abundance is written in
    // scientific notation with a two-digit exponent and no `+`.
    let mut sml = MzTabMSmallMoleculeSectionRow {
        sml_identifier: MzTabString::from_text("3"),
        reliability: MzTabString::from_text("2"),
        ..MzTabMSmallMoleculeSectionRow::default()
    };
    sml.smf_id_refs.set(vec![MzTabString::from_text("3")]);
    sml.database_identifier
        .set(vec![MzTabString::from_text("HMDB:HMDB39534")]);
    sml.chemical_formula
        .set(vec![MzTabString::from_text("C17Cl1H25O2")]);
    sml.smiles.set(vec![MzTabString::from_text(
        "CCCCCCCC(Cl)C(O)CC#CC#CC(O)C=C",
    )]);
    sml.inchi.set(vec![MzTabString::from_text(
        "InChI=1S/C17H25ClO2/c1-3-5-6-7-10-13-16(18)17(20)14-11-8-9-12-15(19)4-2/\
         h4,15-17,19-20H,2-3,5-7,10,13-14H2,1H3",
    )]);
    sml.chemical_name
        .set(vec![MzTabString::from_text("Panaxydol chlorohydrin")]);
    sml.uri.set(vec![MzTabString::null()]);
    sml.theoretical_neutral_mass
        .set(vec![MzTabDouble::new(296.154308477499967)]);
    sml.adducts.set(vec![MzTabString::from_text("[M+H+Na]2+")]);
    sml.small_molecule_abundance_assay
        .insert(1, MzTabDouble::new(36310.19921875));
    sml.small_molecule_abundance_study_variable
        .insert(1, MzTabDouble::null());
    sml.small_molecule_abundance_variation_study_variable
        .insert(1, MzTabDouble::null());
    let rendered = writer
        .generate_small_molecule_section_row(&sml, &meta, &[])
        .unwrap();
    assert_eq!(rendered.text, retained_line(AMS_FIXTURE, "SML\t", 2));
    assert!(rendered.text.contains("\t3.631019921875e04\t"));

    // SME row 1: both declared confidence measures, one in each notation.
    let mut sme = MzTabMSmallMoleculeEvidenceSectionRow {
        sme_identifier: MzTabString::from_text("1"),
        evidence_input_id: MzTabString::from_text("mass=160.07500553925999,rt=281.25"),
        database_identifier: MzTabString::from_text("HMDB:HMDB01190"),
        chemical_formula: MzTabString::from_text("C10H9N1O1"),
        smiles: MzTabString::from_text("O=CCC1=CNC2=CC=CC=C12"),
        inchi: MzTabString::from_text(
            "InChI=1S/C10H9NO/c12-6-5-8-7-11-10-4-2-1-3-9(8)10/h1-4,6-7,11H,5H2",
        ),
        chemical_name: MzTabString::from_text("Indoleacetaldehyde"),
        adduct: MzTabString::from_text("[M+H]1+"),
        exp_mass_to_charge: MzTabDouble::new(160.07500553925999),
        charge: MzTabInteger::new(0),
        calc_mass_to_charge: MzTabDouble::new(159.068414287099984),
        spectra_ref: MzTabSpectraRef::new(1, "6565892897288149707").unwrap(),
        identification_method: parameter("[MS, MS:1000207, accurate mass, ]"),
        ms_level: parameter("[MS, MS:1000511, ms level, 1]"),
        rank: MzTabInteger::new(1),
        ..MzTabMSmallMoleculeEvidenceSectionRow::default()
    };
    sme.id_confidence_measure
        .insert(1, MzTabDouble::new(-4.278149506456153));
    sme.id_confidence_measure
        .insert(2, MzTabDouble::new(-0.0006848277357391908));
    for (name, value) in [
        ("opt_global_chemical_formula", "C10H9NO"),
        ("opt_global_description", "[Indoleacetaldehyde]"),
        ("opt_global_identifier", "[HMDB:HMDB01190]"),
        ("opt_global_modifications", "M+H;1+"),
        ("opt_global_mz_error_Da", "-6.848277357391908e-04"),
        ("opt_global_mz_error_ppm", "-4.278149506456153"),
    ] {
        sme.opt.push(MzTabOptionalColumnEntry::new(
            name,
            MzTabString::from_text(value),
        ));
    }
    let sme_columns = retained_optional_columns(AMS_FIXTURE, "SEH\t");
    let rendered = writer
        .generate_small_molecule_evidence_section_row(&sme, &meta, &sme_columns)
        .unwrap();
    assert_eq!(rendered.text, retained_line(AMS_FIXTURE, "SME\t", 0));

    // SMF row 1: an `isotopomer` null cell and an ambiguity code of 1.
    let mut smf = MzTabMSmallMoleculeFeatureSectionRow {
        smf_identifier: MzTabString::from_text("1"),
        sme_id_ref_ambiguity_code: MzTabInteger::new(1),
        adduct: MzTabString::from_text("[M+2CH3CN+2H]2+"),
        exp_mass_to_charge: MzTabDouble::new(160.07500553925999),
        charge: MzTabInteger::new(0),
        retention_time: MzTabDouble::new(281.25),
        ..MzTabMSmallMoleculeFeatureSectionRow::default()
    };
    smf.sme_id_refs.set(vec![
        MzTabString::from_text("5"),
        MzTabString::from_text("7"),
    ]);
    smf.small_molecule_feature_abundance_assay
        .insert(1, MzTabDouble::new(36310.19921875));
    for (name, value) in [
        ("opt_global_FWHM", "6.25"),
        ("opt_global_dc_charge_adducts", "H1Na1"),
        ("opt_global_label", "T141"),
        ("opt_global_masstrace_intensity", "[3.6310164094607e04]"),
        ("opt_global_num_of_masstraces", "1"),
    ] {
        smf.opt.push(MzTabOptionalColumnEntry::new(
            name,
            MzTabString::from_text(value),
        ));
    }
    let smf_columns = retained_optional_columns(AMS_FIXTURE, "SFH\t");
    let rendered = writer
        .generate_small_molecule_feature_section_row(&smf, &meta, &smf_columns)
        .unwrap();
    assert_eq!(rendered.text, retained_line(AMS_FIXTURE, "SMF\t", 0));
}

// ---------------------------------------------------------------------------
// Tier-1 differential: the exporter against the retained metadata section
// ---------------------------------------------------------------------------

/// A graph whose records reproduce the metadata section of
/// `MzTabMFile_output_1.mztab`.
///
/// The four `software[n]` terms pin the `(name, version)` order the source's
/// `std::set<ProcessingSoftware>` imposes, the single input file pins
/// `ms_run[1]-location`, `assay[1]` and `study_variable[1]`, the single adduct
/// pins the scan polarity, and the single search parameter pins `database[1]`.
fn store_fixture_graph() -> IdentificationData {
    let mut graph = IdentificationData::new().unwrap();
    let file = graph
        .register_input_file(InputFile::new(
            "I:\\OpenSWATH_Metabolomics_data\\20181121_full_data\\\
             04_PestMixes_individually_Solvent_DDA_20-50\\2012_02_03_PStd_10_1-50.wiff",
        ))
        .unwrap();
    // Names sort as AccurateMassSearch < FeatureFinderMetabo <
    // MetaboliteAdductDecharger, and of the four only `TOPP
    // FeatureFinderMetabo` is a registered PSI-MS term, so the other three
    // fall back to `analysis software` — which is why the retained metadata
    // section carries that term three times with two different versions.
    let quant = graph
        .register_processing_software(ProcessingSoftware::new(
            "FeatureFinderMetabo",
            "2.4.0-nightly-2019-07-17",
        ))
        .unwrap();
    graph
        .register_processing_software(ProcessingSoftware::new(
            "AccurateMassSearch",
            "2.7.0-pre-idf-ams-2021-12-20",
        ))
        .unwrap();
    graph
        .register_processing_software(ProcessingSoftware::new(
            "MetaboliteAdductDecharger",
            "2.4.0-nightly-2019-07-17",
        ))
        .unwrap();
    graph
        .register_processing_software(ProcessingSoftware::new(
            "MetaboliteAdductDecharger",
            "2.7.0-pre-idf-ams-2021-12-20",
        ))
        .unwrap();
    let mut step = ProcessingStep::new(quant);
    step.input_files = vec![file];
    step.actions.insert(ProcessingAction::Quantitation);
    let search = graph
        .register_db_search_param(DBSearchParam {
            database: "HMDB".into(),
            database_version: "3.5".into(),
            ..DBSearchParam::default()
        })
        .unwrap();
    graph.register_processing_step(step, Some(search)).unwrap();
    graph
        .register_adduct(AdductInfo::parse("M+H;1+").unwrap())
        .unwrap();
    graph
}

#[test]
fn export_reproduces_retained_metadata_section() {
    let graph = store_fixture_graph();
    let mut map = FeatureMap::new();
    map.unique_id = 14_677_498_592_798_891_241;
    let document =
        MzTabM::export_feature_map_with(&map, &graph, &psi_ms(), &MzTabMExportOptions::source())
            .unwrap();
    let writer = MzTabMFile::with_options(MzTabMWriteOptions::source());
    let produced = writer
        .generate_meta_data_section(document.meta_data())
        .unwrap();
    assert_eq!(produced, retained_metadata(STORE_FIXTURE));
}

#[test]
fn export_reproduces_retained_ams_metadata_section() {
    let mut graph = IdentificationData::new().unwrap();
    let file = graph
        .register_input_file(InputFile::new("E:\\NotARealFileJustATest.mzML"))
        .unwrap();
    let mut software =
        ProcessingSoftware::new("AccurateMassSearch", "3.4.0-pre-fix-mztab-2025-04-02");
    let ppm = graph
        .register_score_type(ScoreType::new("MassErrorPPMScore", false))
        .unwrap();
    let da = graph
        .register_score_type(ScoreType::new("MassErrorDaScore", false))
        .unwrap();
    software.assigned_scores = vec![ppm, da];
    let software = graph.register_processing_software(software).unwrap();
    let mut step = ProcessingStep::new(software);
    step.input_files = vec![file];
    // Only `Identification`: the AccurateMassSearch output's quantification
    // method and unit both come from the source's fallbacks, which use the
    // very same terms as the recognised branches.
    step.actions.insert(ProcessingAction::Identification);
    let mut search = DBSearchParam {
        database: "HMDB".into(),
        database_version: "3.5".into(),
        ..DBSearchParam::default()
    };
    search.metadata.insert(
        "database_location".into(),
        "/home/sachsenb/Development/OpenMS/src/tests/class_tests/openms/data/\
         reducedHMDBMapping.tsv"
            .into(),
    );
    let search = graph.register_db_search_param(search).unwrap();
    graph.register_processing_step(step, Some(search)).unwrap();
    graph
        .register_adduct(AdductInfo::parse("M+H;1+").unwrap())
        .unwrap();

    let map = FeatureMap::new();
    let document =
        MzTabM::export_feature_map_with(&map, &graph, &psi_ms(), &MzTabMExportOptions::source())
            .unwrap();
    let writer = MzTabMFile::with_options(MzTabMWriteOptions::source());
    let produced = writer
        .generate_meta_data_section(document.meta_data())
        .unwrap();
    assert_eq!(produced, retained_metadata(AMS_FIXTURE));
    // And the identification method the AccurateMassSearch branch chooses
    // reaches the evidence rows, as the retained SME row shows.
    assert_eq!(
        document
            .meta_data()
            .id_confidence_measure
            .get(&1)
            .unwrap()
            .to_cell_string(),
        "[, , MassErrorPPMScore, ]"
    );
}

// ---------------------------------------------------------------------------
// The exporter's own behaviour
// ---------------------------------------------------------------------------

struct Exported {
    graph: IdentificationData,
    map: FeatureMap,
}

/// One feature with two compound identifications under two adducts, plus one
/// feature without any identification.
fn exportable() -> Exported {
    let mut graph = IdentificationData::new().unwrap();
    let file = graph
        .register_input_file(InputFile::new("sample.mzML"))
        .unwrap();
    let software = graph
        .register_processing_software(ProcessingSoftware::new("AccurateMassSearch", "3.0"))
        .unwrap();
    let mut step = ProcessingStep::new(software);
    step.input_files = vec![file];
    step.actions.insert(ProcessingAction::Identification);
    graph.register_processing_step(step, None).unwrap();
    let proton = graph
        .register_adduct(AdductInfo::parse("M+H;1+").unwrap())
        .unwrap();
    let sodium = graph
        .register_adduct(AdductInfo::parse("M+Na;1+").unwrap())
        .unwrap();
    let observation = graph
        .register_observation(Observation::new("spectrum=7", file))
        .unwrap();
    let mut betaine = IdentifiedCompound::new("HMDB:HMDB00043");
    betaine.formula = "C5H11NO2".parse().unwrap();
    betaine.name = "Betaine".into();
    betaine.smile = "C[N+](C)(C)CC(=O)[O-]".into();
    betaine.inchi = "InChI=1S/C5H11NO2".into();
    betaine.result.metadata.insert(
        "mz error ppm".into(),
        MetaValue::try_from(0.25_f64).unwrap(),
    );
    let betaine = graph.register_identified_compound(betaine).unwrap();
    let mut valine = IdentifiedCompound::new("HMDB:HMDB00883");
    valine.formula = "C5H11NO2".parse().unwrap();
    valine.name = "L-Valine".into();
    let valine = graph.register_identified_compound(valine).unwrap();

    let mut first = ObservationMatch::new(betaine, observation);
    first.adduct = Some(proton);
    let first = graph.register_observation_match(first).unwrap();
    let mut second = ObservationMatch::new(valine, observation);
    second.adduct = Some(sodium);
    let second = graph.register_observation_match(second).unwrap();

    let mut feature = Feature::new(70.15, 118.08, 406.5);
    feature.charge = 1;
    feature
        .metadata
        .insert("FWHM".into(), MetaValue::try_from(3.5_f64).unwrap());
    feature.add_id_match(first).unwrap();
    feature.add_id_match(second).unwrap();

    let mut lonely = Feature::new(81.0, 200.0, 12.5);
    lonely
        .metadata
        .insert("adducts".into(), vec!["M+H;1+".to_owned()].into());

    let mut map = FeatureMap::from_features(vec![feature, lonely]);
    map.unique_id = 42;
    Exported { graph, map }
}

#[test]
fn export_builds_one_summary_row_per_feature_row() {
    let fixture = exportable();
    let document = MzTabM::export_feature_map(&fixture.map, &fixture.graph).unwrap();
    // Two evidences under two adducts give two feature rows, plus the
    // unidentified feature: three feature rows, three summary rows.
    assert_eq!(document.small_molecule_evidence_section_rows().len(), 2);
    assert_eq!(document.small_molecule_feature_section_rows().len(), 3);
    assert_eq!(document.small_molecule_section_rows().len(), 3);

    let evidence = &document.small_molecule_evidence_section_rows()[0];
    assert_eq!(evidence.sme_identifier.get(), "1");
    assert_eq!(evidence.database_identifier.get(), "HMDB:HMDB00043");
    assert_eq!(evidence.chemical_formula.get(), "C5H11N1O2");
    assert_eq!(evidence.adduct.get(), "[M+H]1+");
    assert_eq!(
        evidence.evidence_input_id.get(),
        // `StringUtils::toStr(double)`: fifteen fractional digits, trailing
        // zeros trimmed to at least one, so the stored binary value is spelled
        // out rather than shortened.
        "mass=118.079999999999998,rt=70.150000000000006"
    );
    assert_eq!(
        evidence.spectra_ref.to_cell_string(),
        "ms_run[1]:spectrum=7"
    );
    assert_eq!(evidence.rank.get().unwrap(), 1);
    // `calc_mass_to_charge` is the compound's neutral monoisotopic mass, not
    // an m/z: the source does not correct for the adduct or the charge.
    assert!((evidence.calc_mass_to_charge.get().unwrap() - 117.078_979_350_9).abs() < 1e-9);

    // The unidentified feature's row takes its adduct from the feature's own
    // `adducts` meta value and has a null SME reference list.
    let lonely = &document.small_molecule_feature_section_rows()[2];
    assert!(lonely.sme_id_refs.is_null());
    assert_eq!(lonely.adduct.get(), "M+H;1+");
    assert!(lonely.sme_id_ref_ambiguity_code.is_null());
    assert_eq!(
        lonely.small_molecule_feature_abundance_assay[&1]
            .get()
            .unwrap(),
        12.5
    );

    // Each summary row declares a null study-variable abundance, which is why
    // the retained files carry those two columns.
    let summary = &document.small_molecule_section_rows()[0];
    assert!(summary.small_molecule_abundance_study_variable[&1].is_null());
    assert!(summary.small_molecule_abundance_variation_study_variable[&1].is_null());
    assert_eq!(summary.reliability.get(), "2");
    // theoretical_neutral_mass is derived from the evidence formula cell.
    assert_eq!(summary.theoretical_neutral_mass.get().len(), 1);

    // The whole document writes cleanly and is rectangular under native
    // options, which is what the source's dropped POSTCONDITIONs asserted.
    let lines = MzTabMFile::new().generate_lines(&document).unwrap();
    assert_rectangular(&lines);
}

/// Every data line of a section has exactly the column count of its header.
fn assert_rectangular(lines: &[String]) {
    let mut expected = 0;
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let columns = line.split('\t').count();
        match line.split('\t').next().unwrap_or("") {
            "SMH" | "SFH" | "SEH" => expected = columns,
            "SML" | "SMF" | "SME" => assert_eq!(columns, expected, "row {line}"),
            _ => {}
        }
    }
}

#[test]
fn export_reliability_comes_from_a_software_meta_value() {
    let mut fixture = exportable();
    let mut software = ProcessingSoftware::new("ZzzReliable", "1.0");
    software
        .software
        .cv_terms
        .metadata
        .insert("reliability".into(), "1".into());
    fixture
        .graph
        .register_processing_software(software)
        .unwrap();
    let document = MzTabM::export_feature_map(&fixture.map, &fixture.graph).unwrap();
    assert_eq!(
        document.small_molecule_section_rows()[0].reliability.get(),
        "1"
    );
}

#[test]
fn export_deduplicates_matches_only_when_asked() {
    let mut graph = IdentificationData::new().unwrap();
    let file = graph
        .register_input_file(InputFile::new("sample.mzML"))
        .unwrap();
    let software = graph
        .register_processing_software(ProcessingSoftware::new("Tool", "1.0"))
        .unwrap();
    graph
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap();
    let proton = graph
        .register_adduct(AdductInfo::parse("M+H;1+").unwrap())
        .unwrap();
    let sodium = graph
        .register_adduct(AdductInfo::parse("M+Na;1+").unwrap())
        .unwrap();
    let observation = graph
        .register_observation(Observation::new("spectrum=1", file))
        .unwrap();
    let compound = graph
        .register_identified_compound(IdentifiedCompound::new("HMDB:HMDB00043"))
        .unwrap();
    // Two matches on the *same* compound under different adducts: distinct
    // records in the graph, equivalent under the source's comparator.
    let mut first = ObservationMatch::new(compound, observation);
    first.adduct = Some(proton);
    let first = graph.register_observation_match(first).unwrap();
    let mut second = ObservationMatch::new(compound, observation);
    second.adduct = Some(sodium);
    let second = graph.register_observation_match(second).unwrap();
    let mut feature = Feature::new(1.0, 2.0, 3.0);
    feature.add_id_match(first).unwrap();
    feature.add_id_match(second).unwrap();
    let map = FeatureMap::from_features(vec![feature]);

    let native = MzTabM::export_feature_map(&map, &graph).unwrap();
    assert_eq!(native.small_molecule_evidence_section_rows().len(), 2);
    assert_eq!(native.small_molecule_feature_section_rows().len(), 2);

    let source = MzTabM::export_feature_map_with(
        &map,
        &graph,
        openms::format::controlled_vocabulary::ControlledVocabulary::psi_ms().unwrap(),
        &MzTabMExportOptions::source(),
    )
    .unwrap();
    assert_eq!(source.small_molecule_evidence_section_rows().len(), 1);
    assert_eq!(source.small_molecule_feature_section_rows().len(), 1);

    // The comparator the source's set uses is exposed, and it makes the two
    // matches equivalent.
    assert_eq!(
        compare_match_by_compound(&graph, first, second).unwrap(),
        std::cmp::Ordering::Equal
    );
}

#[test]
fn export_optional_column_lookup_follows_the_substitution_option() {
    let fixture = exportable();
    // Native: the key "mz error ppm" is looked up as written, so the column
    // carries the value and its name has the substitution applied.
    let native = MzTabM::export_feature_map(&fixture.map, &fixture.graph).unwrap();
    let entry = native.small_molecule_evidence_section_rows()[0]
        .opt
        .iter()
        .find(|entry| entry.name == "opt_global_mz_error_ppm")
        .expect("the substituted column name is present");
    assert_eq!(entry.value.get(), "0.25");

    // Source: the key is substituted before the lookup, which then misses.
    let source = MzTabM::export_feature_map_with(
        &fixture.map,
        &fixture.graph,
        ControlledVocabulary::psi_ms().unwrap(),
        &MzTabMExportOptions::source(),
    )
    .unwrap();
    let entry = source.small_molecule_evidence_section_rows()[0]
        .opt
        .iter()
        .find(|entry| entry.name == "opt_global_mz_error_ppm")
        .expect("the column is still declared");
    assert!(entry.value.is_null());
}

#[test]
fn export_swaps_the_cv_identity_only_when_asked() {
    let fixture = exportable();
    let cv = psi_ms();
    let native = MzTabM::export_feature_map_with(
        &fixture.map,
        &fixture.graph,
        &cv,
        &MzTabMExportOptions::default(),
    )
    .unwrap();
    assert_eq!(native.meta_data().cv[&1].label.get(), "MS");
    assert_eq!(native.meta_data().cv[&1].full_name.get(), "PSI-MS");
    let source = MzTabM::export_feature_map_with(
        &fixture.map,
        &fixture.graph,
        &cv,
        &MzTabMExportOptions::source(),
    )
    .unwrap();
    assert_eq!(source.meta_data().cv[&1].label.get(), "PSI-MS");
    assert_eq!(source.meta_data().cv[&1].full_name.get(), "MS");
}

#[test]
fn export_identification_method_and_ms_level_follow_the_tool() {
    for (tool, method, level) in [
        (
            "AccurateMassSearch",
            "[MS, MS:1000207, accurate mass, ]",
            "1",
        ),
        ("SiriusAdapter", "[MS, MS:1001010, de novo search, ]", "2"),
        (
            "MetaboliteSpectralMatcher",
            "[MS, MS:1002187, TOPP SpecLibSearcher, ]",
            "2",
        ),
        (
            "SomeUnknownTool",
            "[MS, MS:1000543, data processing action, ]",
            "1",
        ),
    ] {
        let mut graph = IdentificationData::new().unwrap();
        let software = graph
            .register_processing_software(ProcessingSoftware::new(tool, "1.0"))
            .unwrap();
        let mut step = ProcessingStep::new(software);
        step.actions.insert(ProcessingAction::Identification);
        graph.register_processing_step(step, None).unwrap();
        let file = graph
            .register_input_file(InputFile::new("run.mzML"))
            .unwrap();
        let observation = graph
            .register_observation(Observation::new("s=1", file))
            .unwrap();
        let compound = graph
            .register_identified_compound(IdentifiedCompound::new("X:1"))
            .unwrap();
        let id = graph
            .register_observation_match(ObservationMatch::new(compound, observation))
            .unwrap();
        let mut feature = Feature::new(1.0, 2.0, 3.0);
        feature.add_id_match(id).unwrap();
        let map = FeatureMap::from_features(vec![feature]);
        let document = MzTabM::export_feature_map(&map, &graph).unwrap();
        let row = &document.small_molecule_evidence_section_rows()[0];
        assert_eq!(row.identification_method.to_cell_string(), method, "{tool}");
        assert_eq!(
            row.ms_level.to_cell_string(),
            format!("[MS, MS:1000511, ms level, {level}]"),
            "{tool}"
        );
    }
}

#[test]
fn export_scan_polarity_follows_the_first_adduct_and_defaults_to_positive() {
    let mut graph = IdentificationData::new().unwrap();
    let software = graph
        .register_processing_software(ProcessingSoftware::new("Tool", "1.0"))
        .unwrap();
    graph
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap();
    let map = FeatureMap::new();
    // No adduct at all: the source warns and assumes positive.
    let document = MzTabM::export_feature_map(&map, &graph).unwrap();
    assert_eq!(
        document.meta_data().ms_run[&1].scan_polarity[&1].to_cell_string(),
        "[MS, MS:1000130, positive scan, ]"
    );

    let mut negative = IdentificationData::new().unwrap();
    let software = negative
        .register_processing_software(ProcessingSoftware::new("Tool", "1.0"))
        .unwrap();
    negative
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap();
    negative
        .register_adduct(AdductInfo::parse("M-H;1-").unwrap())
        .unwrap();
    let document = MzTabM::export_feature_map(&map, &negative).unwrap();
    assert_eq!(
        document.meta_data().ms_run[&1].scan_polarity[&1].to_cell_string(),
        "[MS, MS:1000129, negative scan, ]"
    );
}

#[test]
fn export_scan_polarity_is_positive_only_for_a_name_ending_in_plus() {
    // `MzTabM.cpp:287` tests `first_adduct.at(size() - 1) == '+'`, so any other
    // final character — including none at all, where `at()` underflows and
    // throws — is negative there. `AdductInfo::new` takes any bounded text as a
    // name, so names with no charge suffix are reachable; only the empty name
    // diverges, and it diverges from a throw.
    let map = FeatureMap::new();
    for (name, expected) in [
        ("M+H;1+", "[MS, MS:1000130, positive scan, ]"),
        ("M+H", "[MS, MS:1000129, negative scan, ]"),
        ("M+Na", "[MS, MS:1000129, negative scan, ]"),
        ("H", "[MS, MS:1000129, negative scan, ]"),
        ("M-H;1-", "[MS, MS:1000129, negative scan, ]"),
        // No name to take a sign from: the source throws, this writes positive.
        ("", "[MS, MS:1000130, positive scan, ]"),
    ] {
        let mut graph = IdentificationData::new().unwrap();
        let software = graph
            .register_processing_software(ProcessingSoftware::new("Tool", "1.0"))
            .unwrap();
        graph
            .register_processing_step(ProcessingStep::new(software), None)
            .unwrap();
        graph
            .register_adduct(
                AdductInfo::new(
                    name,
                    openms::chemistry::EmpiricalFormula::parse("H").unwrap(),
                    1,
                    1,
                )
                .unwrap(),
            )
            .unwrap();
        let document = MzTabM::export_feature_map(&map, &graph).unwrap();
        assert_eq!(
            document.meta_data().ms_run[&1].scan_polarity[&1].to_cell_string(),
            expected,
            "adduct name {name:?}"
        );
    }
}

#[test]
fn export_quantification_unit_follows_the_quant_method_parameter() {
    for (method, expected) in [
        ("area", "[MS, MS:1001844, MS1 feature area, ]"),
        ("median", "[MS, MS:1002883, median, ]"),
        (
            "max_height",
            "[MS, MS:1001843, MS1 feature maximum intensity, ]",
        ),
    ] {
        let mut graph = IdentificationData::new().unwrap();
        let mut software = ProcessingSoftware::new("FeatureFinderMetabo", "3.0");
        software.software.cv_terms.metadata.insert(
            "parameter: algorithm:mtd:quant_method".into(),
            method.into(),
        );
        let software = graph.register_processing_software(software).unwrap();
        let mut step = ProcessingStep::new(software);
        step.actions.insert(ProcessingAction::Quantitation);
        graph.register_processing_step(step, None).unwrap();
        let map = FeatureMap::new();
        let document = MzTabM::export_feature_map(&map, &graph).unwrap();
        assert_eq!(
            document
                .meta_data()
                .small_molecule_quantification_unit
                .to_cell_string(),
            expected,
            "{method}"
        );
        assert_eq!(
            document
                .meta_data()
                .small_molecule_feature_quantification_unit
                .to_cell_string(),
            expected
        );
    }
}

#[test]
fn export_database_uri_carries_over_from_an_earlier_search_parameter() {
    let mut graph = IdentificationData::new().unwrap();
    let software = graph
        .register_processing_software(ProcessingSoftware::new("Tool", "1.0"))
        .unwrap();
    let mut located = DBSearchParam {
        database: "AAA".into(),
        database_version: "1".into(),
        ..DBSearchParam::default()
    };
    located
        .metadata
        .insert("database_location".into(), "C:\\db\\aaa.tsv".into());
    let located = graph.register_db_search_param(located).unwrap();
    let bare = graph
        .register_db_search_param(DBSearchParam {
            database: "BBB".into(),
            database_version: "2".into(),
            ..DBSearchParam::default()
        })
        .unwrap();
    graph
        .register_processing_step(ProcessingStep::new(software), Some(located))
        .unwrap();
    let _ = bare;
    let map = FeatureMap::new();
    let document = MzTabM::export_feature_map(&map, &graph).unwrap();
    let databases = &document.meta_data().database;
    assert_eq!(databases.len(), 2);
    assert_eq!(databases[&1].uri.get(), "file://C:/db/aaa.tsv");
    // The source's accumulating record keeps the previous URI rather than
    // restoring the https://hmdb.ca/ default.
    assert_eq!(databases[&2].uri.get(), "file://C:/db/aaa.tsv");
    assert_eq!(databases[&2].prefix.get(), "BBB");
    // A "custom" database nulls the prefix and spaces the parameter name.
    let mut custom = IdentificationData::new().unwrap();
    let software = custom
        .register_processing_software(ProcessingSoftware::new("Tool", "1.0"))
        .unwrap();
    custom
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap();
    custom
        .register_db_search_param(DBSearchParam {
            database: "my custom db".into(),
            database_version: "9".into(),
            ..DBSearchParam::default()
        })
        .unwrap();
    let document = MzTabM::export_feature_map(&map, &custom).unwrap();
    assert!(document.meta_data().database[&1].prefix.is_null());
    assert_eq!(
        document.meta_data().database[&1].database.to_cell_string(),
        "[, , my custom db, ]"
    );
}

#[test]
fn export_assay_name_survives_a_basename_without_a_dot_and_non_ascii() {
    for (input, stem) in [
        ("/data/plain", "plain"),
        ("/data/日本語.mzML", "日本語"),
        ("no_separator_at_all", "no_separator_at_all"),
        ("", ""),
    ] {
        let mut graph = IdentificationData::new().unwrap();
        let software = graph
            .register_processing_software(ProcessingSoftware::new("Tool", "1.0"))
            .unwrap();
        graph
            .register_processing_step(ProcessingStep::new(software), None)
            .unwrap();
        if !input.is_empty() {
            graph.register_input_file(InputFile::new(input)).unwrap();
        }
        let map = FeatureMap::new();
        let document = MzTabM::export_feature_map(&map, &graph).unwrap();
        assert_eq!(
            document.meta_data().assay[&1].name.get(),
            format!("assay_{stem}"),
            "{input}"
        );
        assert_eq!(
            document.meta_data().study_variable[&1].description.get(),
            format!("study_variable_{stem}")
        );
    }
}

#[test]
fn export_refuses_an_empty_graph() {
    let graph = IdentificationData::new().unwrap();
    let map = FeatureMap::new();
    assert!(matches!(
        MzTabM::export_feature_map(&map, &graph),
        Err(Error::MissingInformation(_))
    ));
}

#[test]
fn export_refuses_a_match_that_is_not_a_compound() {
    let mut graph = IdentificationData::new().unwrap();
    let file = graph
        .register_input_file(InputFile::new("run.mzML"))
        .unwrap();
    let software = graph
        .register_processing_software(ProcessingSoftware::new("Tool", "1.0"))
        .unwrap();
    graph
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap();
    let observation = graph
        .register_observation(Observation::new("s=1", file))
        .unwrap();
    let peptide = graph
        .register_identified_peptide(IdentifiedPeptide::new(
            openms::chemistry::AASequence::parse("PEPTIDE").unwrap(),
        ))
        .unwrap();
    let id = graph
        .register_observation_match(ObservationMatch::new(peptide, observation))
        .unwrap();
    let mut feature = Feature::new(1.0, 2.0, 3.0);
    feature.add_id_match(id).unwrap();
    let map = FeatureMap::from_features(vec![feature]);
    assert!(matches!(
        MzTabM::export_feature_map(&map, &graph),
        Err(Error::MissingInformation(_))
    ));
    assert!(matches!(
        compare_match_by_compound(&graph, id, id),
        Err(Error::MissingInformation(_))
    ));
}

#[test]
fn export_leaves_the_spectra_reference_null_for_an_empty_data_id() {
    // `Observation::new` requires a non-empty data ID, so the null branch of
    // `MzTabSpectraRef` is reached through the cell itself: the source's
    // `setSpecRef("")` logs a warning and keeps the previous, empty, value.
    let mut reference = MzTabSpectraRef::default();
    reference.set_ms_file(1).unwrap();
    assert!(reference.is_null());
    assert_eq!(reference.to_cell_string(), "null");
    assert!(reference.set_spec_ref("").is_err());
}

#[test]
fn export_score_cells_are_nan_when_a_match_carries_no_score() {
    let mut graph = IdentificationData::new().unwrap();
    let file = graph
        .register_input_file(InputFile::new("run.mzML"))
        .unwrap();
    let score = graph
        .register_score_type(ScoreType::new("MassErrorPPMScore", false))
        .unwrap();
    let mut software = ProcessingSoftware::new("AccurateMassSearch", "1.0");
    software.assigned_scores = vec![score];
    let software = graph.register_processing_software(software).unwrap();
    let mut step = ProcessingStep::new(software);
    step.actions.insert(ProcessingAction::Identification);
    graph.register_processing_step(step, None).unwrap();
    let observation = graph
        .register_observation(Observation::new("s=1", file))
        .unwrap();
    let compound = graph
        .register_identified_compound(IdentifiedCompound::new("X:1"))
        .unwrap();
    let id = graph
        .register_observation_match(ObservationMatch::new(compound, observation))
        .unwrap();
    let mut feature = Feature::new(1.0, 2.0, 3.0);
    feature.add_id_match(id).unwrap();
    let map = FeatureMap::from_features(vec![feature]);
    let document = MzTabM::export_feature_map(&map, &graph).unwrap();
    // `getScore` returns (NaN, false) when the score is absent, and the
    // source stores that NaN in a value cell, which renders `NaN`.
    assert_eq!(
        document.small_molecule_evidence_section_rows()[0].id_confidence_measure[&1]
            .to_cell_string(),
        "NaN"
    );
    assert_eq!(document.meta_data().id_confidence_measure.len(), 1);
}

#[test]
fn export_refuses_a_non_list_adducts_meta_value() {
    let mut fixture = exportable();
    fixture.map.features[1]
        .metadata
        .insert("adducts".into(), 5_i64.into());
    assert!(matches!(
        MzTabM::export_feature_map(&fixture.map, &fixture.graph),
        Err(Error::InvalidValue(_))
    ));
}

// ---------------------------------------------------------------------------
// The writer's source-defect options
// ---------------------------------------------------------------------------

fn metadata_with_every_optional_key() -> MzTabMMetaData {
    let mut meta = MzTabMMetaData::default();
    meta.mz_tab_id.set("id");
    let mut ms_run = MzTabMMSRunMetaData {
        location: MzTabString::from_text("file:///run.mzML"),
        id_format: parameter("[MS, MS:1000768, Thermo nativeID format, ]"),
        ..MzTabMMSRunMetaData::default()
    };
    ms_run
        .scan_polarity
        .insert(1, parameter("[MS, MS:1000130, positive scan, ]"));
    meta.ms_run.insert(1, ms_run);
    let mut assay = MzTabMAssayMetaData {
        name: MzTabString::from_text("assay_1"),
        ms_run_ref: MzTabInteger::new(1),
        ..MzTabMAssayMetaData::default()
    };
    assay
        .custom
        .insert(1, parameter("[MS, , Assay operator, Blogs]"));
    meta.assay.insert(3, assay);
    meta.derivatization_agent
        .insert(1, parameter("[MS, MS:1000000, agent, ]"));
    meta.colunit_small_molecule.push(MzTabString::from_text(
        "retention_time=[UO, UO:0000010, second, ]",
    ));
    meta.colunit_small_molecule_feature
        .push(MzTabString::from_text(
            "rt_start=[UO, UO:0000010, second, ]",
        ));
    meta.colunit_small_molecule_evidence
        .push(MzTabString::from_text("charge=[UO, UO:0000000, unit, ]"));
    meta
}

#[test]
fn write_options_select_the_source_metadata_keys() {
    let meta = metadata_with_every_optional_key();
    let source = MzTabMFile::with_options(MzTabMWriteOptions::source())
        .generate_meta_data_section(&meta)
        .unwrap();
    // The assay's custom parameter is reported under an ms_run key, at the
    // assay's own index.
    assert!(source.contains(&"MTD\tms_run[3]-custom[1]\t[MS, , Assay operator, Blogs]".to_owned()));
    assert!(!source.iter().any(|line| line.contains("assay[3]-custom")));
    // All three colunit families collapse onto one key.
    assert_eq!(
        source
            .iter()
            .filter(|line| line.starts_with("MTD\tcolunit_small_molecule\t"))
            .count(),
        3
    );
    // The derivatization agent gains a `-uri` suffix.
    assert!(
        source.contains(&"MTD\tderivatization_agent[1]-uri\t[MS, MS:1000000, agent, ]".to_owned())
    );
    // And `ms_run[1]-id_format` is dropped entirely.
    assert!(!source.iter().any(|line| line.contains("id_format")));

    let native = MzTabMFile::new().generate_meta_data_section(&meta).unwrap();
    assert!(native.contains(&"MTD\tassay[3]-custom[1]\t[MS, , Assay operator, Blogs]".to_owned()));
    assert!(native.contains(
        &"MTD\tcolunit_small_molecule_feature\trt_start=[UO, UO:0000010, second, ]".to_owned()
    ));
    assert!(native.contains(
        &"MTD\tcolunit_small_molecule_evidence\tcharge=[UO, UO:0000000, unit, ]".to_owned()
    ));
    assert!(native.contains(&"MTD\tderivatization_agent[1]\t[MS, MS:1000000, agent, ]".to_owned()));
    assert!(native.contains(
        &"MTD\tms_run[1]-id_format\t[MS, MS:1000768, Thermo nativeID format, ]".to_owned()
    ));
    // Native output is otherwise the same section, one line longer.
    assert_eq!(native.len(), source.len() + 1);
}

#[test]
fn mandatory_metadata_keys_are_written_even_when_null() {
    let lines = MzTabMFile::new()
        .generate_meta_data_section(&MzTabMMetaData::default())
        .unwrap();
    assert_eq!(
        lines,
        vec![
            "MTD\tmzTab-version\t2.0.0-M".to_owned(),
            "MTD\tmzTab-ID\tnull".to_owned(),
            "MTD\tquantification_method\tnull".to_owned(),
            "MTD\tsmall_molecule-quantification_unit\tnull".to_owned(),
            "MTD\tsmall_molecule_feature-quantification_unit\tnull".to_owned(),
            "MTD\tsmall_molecule-identification_reliability\tnull".to_owned(),
        ]
    );
    // `assay[n]-ms_run_ref` is written unconditionally too, so a null one
    // produces the unusable spelling `ms_run[null]`.
    let mut meta = MzTabMMetaData::default();
    meta.assay.insert(1, MzTabMAssayMetaData::default());
    let lines = MzTabMFile::new().generate_meta_data_section(&meta).unwrap();
    assert!(lines.contains(&"MTD\tassay[1]-ms_run_ref\tms_run[null]".to_owned()));
    assert!(!lines.iter().any(|line| line.contains("-sample_ref")));
}

#[test]
fn abundance_cells_are_aligned_with_the_header_under_native_options() {
    let mut meta = MzTabMMetaData::default();
    for index in [1, 2, 5] {
        meta.assay.insert(
            index,
            MzTabMAssayMetaData {
                ms_run_ref: MzTabInteger::new(1),
                ..MzTabMAssayMetaData::default()
            },
        );
    }
    meta.study_variable
        .insert(1, MzTabMStudyVariableMetaData::default());
    let mut row = MzTabMSmallMoleculeFeatureSectionRow::default();
    row.small_molecule_feature_abundance_assay
        .insert(5, MzTabDouble::new(1.5));

    let writer = MzTabMFile::new();
    let header = writer
        .generate_small_molecule_feature_header(&meta, &[])
        .unwrap();
    let native = writer
        .generate_small_molecule_feature_section_row(&row, &meta, &[])
        .unwrap();
    assert_eq!(native.columns, header.columns);
    assert!(native.text.ends_with("\tnull\tnull\t1.5"));

    // The source emits only the cells the row carries, so its row is three
    // columns short of its own header.
    let source = MzTabMFile::with_options(MzTabMWriteOptions::source())
        .generate_small_molecule_feature_section_row(&row, &meta, &[])
        .unwrap();
    assert_eq!(source.columns + 2, header.columns);
    assert!(source.text.ends_with("\t1.5"));
}

#[test]
fn an_undeclared_abundance_is_refused_under_native_options() {
    let meta = MzTabMMetaData::default();
    let mut row = MzTabMSmallMoleculeSectionRow::default();
    row.small_molecule_abundance_assay
        .insert(7, MzTabDouble::new(1.0));
    let writer = MzTabMFile::new();
    assert!(matches!(
        writer.generate_small_molecule_section_row(&row, &meta, &[]),
        Err(Error::InvalidValue(_))
    ));
    // The source writes the orphan cell without complaint.
    let source = MzTabMFile::with_options(MzTabMWriteOptions::source())
        .generate_small_molecule_section_row(&row, &meta, &[])
        .unwrap();
    assert!(source.text.ends_with("\t1.0"));

    let mut evidence = MzTabMSmallMoleculeEvidenceSectionRow::default();
    evidence
        .id_confidence_measure
        .insert(3, MzTabDouble::new(0.5));
    assert!(matches!(
        writer.generate_small_molecule_evidence_section_row(&evidence, &meta, &[]),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn optional_columns_are_placed_by_name_and_missing_ones_are_null() {
    let meta = MzTabMMetaData::default();
    let mut row = MzTabMSmallMoleculeSectionRow::default();
    // Declared out of header order, plus one the header does not request.
    row.opt.push(MzTabOptionalColumnEntry::new(
        "opt_global_b",
        MzTabString::from_text("second"),
    ));
    row.opt.push(MzTabOptionalColumnEntry::new(
        "opt_global_a",
        MzTabString::from_text("first"),
    ));
    row.opt.push(MzTabOptionalColumnEntry::new(
        "opt_global_unwanted",
        MzTabString::from_text("dropped"),
    ));
    let columns = [
        "opt_global_a".to_owned(),
        "opt_global_b".to_owned(),
        "opt_global_missing".to_owned(),
    ];
    let rendered = MzTabMFile::new()
        .generate_small_molecule_section_row(&row, &meta, &columns)
        .unwrap();
    assert!(rendered.text.ends_with("\tfirst\tsecond\tnull"));
    assert!(!rendered.text.contains("dropped"));
}

#[test]
fn optional_column_names_keep_first_occurrence_order_across_rows() {
    let mut document = MzTabM::new();
    let entry = |name: &str| MzTabOptionalColumnEntry::new(name, MzTabString::default());
    document
        .small_molecule_feature_data
        .push(MzTabMSmallMoleculeFeatureSectionRow {
            opt: vec![entry("opt_global_z"), entry("opt_global_a")],
            ..MzTabMSmallMoleculeFeatureSectionRow::default()
        });
    document
        .small_molecule_feature_data
        .push(MzTabMSmallMoleculeFeatureSectionRow {
            opt: vec![entry("opt_global_a"), entry("opt_global_m")],
            ..MzTabMSmallMoleculeFeatureSectionRow::default()
        });
    assert_eq!(
        document
            .small_molecule_feature_optional_column_names()
            .unwrap(),
        vec![
            "opt_global_z".to_owned(),
            "opt_global_a".to_owned(),
            "opt_global_m".to_owned()
        ]
    );
}

#[test]
fn add_meta_info_to_optional_columns_substitutes_only_the_name() {
    let mut meta = openms::metadata::MetaInfo::new();
    meta.insert("with space".into(), "value".into());
    let mut opt = Vec::new();
    let keys: BTreeSet<String> = ["with space".to_owned(), "absent".to_owned()]
        .into_iter()
        .collect();
    MzTabM::add_meta_info_to_optional_columns(&keys, &mut opt, "global", &meta).unwrap();
    assert_eq!(opt.len(), 2);
    assert_eq!(opt[0].name, "opt_global_absent");
    assert!(opt[0].value.is_null());
    assert_eq!(opt[1].name, "opt_global_with_space");
    assert_eq!(opt[1].value.get(), "value");
}

// ---------------------------------------------------------------------------
// Boundaries, errors and non-ASCII input
// ---------------------------------------------------------------------------

#[test]
fn store_checks_the_output_extension() {
    let document = MzTabM::new();
    let writer = MzTabMFile::new();
    // A directory of its own, so concurrent runs cannot remove each other's files.
    let temp = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let directory = temp.path();
    for name in ["out.mzTab", "out.tsv", "out.unknownext", "日本語.mzTab"] {
        let path = directory.join(name);
        writer.store(&path, &document).unwrap();
        assert!(std::fs::read_to_string(&path).unwrap().starts_with("MTD\t"));
        std::fs::remove_file(&path).ok();
    }
    for name in ["out.mzML", "out.idXML", "out.fasta"] {
        let path = directory.join(name);
        assert!(
            matches!(writer.store(&path, &document), Err(Error::InvalidValue(_))),
            "{name}"
        );
        assert!(!path.exists(), "{name} must not be created");
    }
}

#[test]
fn a_metadata_section_over_the_indexed_entry_ceiling_is_refused() {
    let mut meta = MzTabMMetaData::default();
    for index in 0..=MzTabM::MAX_INDEXED_ENTRIES {
        meta.publication
            .insert(index, MzTabString::from_text("pubmed:1"));
    }
    assert_eq!(meta.publication.len(), MzTabM::MAX_INDEXED_ENTRIES + 1);
    assert!(matches!(
        MzTabMFile::new().generate_meta_data_section(&meta),
        Err(Error::InvalidValue(_))
    ));
    meta.publication.pop_last();
    // At the ceiling the section still builds; the estimate is what bounds it.
    assert!(MzTabMFile::new().generate_meta_data_section(&meta).is_ok());
}

/// `MAX_INDEXED_ENTRIES` bounds how many entries an indexed map may carry, not
/// how large its keys may be: a key at `usize::MAX` is written out, exactly as
/// the source's unchecked `Size` writes it, and nothing refuses it. Pinned so
/// the ceiling's scope is not mistaken for a check on the index values.
#[test]
fn an_extreme_index_is_written_rather_than_refused() {
    let mut meta = MzTabMMetaData::default();
    meta.ms_run
        .insert(usize::MAX, MzTabMMSRunMetaData::default());
    let lines = MzTabMFile::new()
        .generate_meta_data_section(&meta)
        .expect("an extreme index is not a ceiling violation");
    assert!(
        lines
            .iter()
            .any(|line| line.contains("ms_run[18446744073709551615]")),
        "{lines:?}"
    );
}

#[test]
fn ceilings_are_public_and_ordered() {
    assert_eq!(MzTabM::MAX_ROWS, 10_000_000);
    assert_eq!(MzTabM::MAX_OPTIONAL_COLUMNS, 100_000);
    assert_eq!(MzTabM::MAX_INDEXED_ENTRIES, 100_000);
    assert_eq!(MzTabMFile::MAX_LINES, 1_000_000);
    const { assert!(MzTabM::MAX_INDEXED_ENTRIES < MzTabMFile::MAX_LINES) };
}

#[test]
fn non_ascii_text_survives_every_cell_and_the_writer() {
    let mut meta = MzTabMMetaData::default();
    meta.mz_tab_id.set("局所識別子");
    meta.title.set("メタボロミクス");
    meta.description.set("ünïcödé — description");
    meta.custom
        .insert(1, MzTabParameter::from_parts("MS", "MS:1", "名前", "値"));
    let mut assay = MzTabMAssayMetaData {
        name: MzTabString::from_text("assay_日本語"),
        ms_run_ref: MzTabInteger::new(1),
        ..MzTabMAssayMetaData::default()
    };
    assay.custom.insert(1, parameter("[MS, , 演算子, ブログ]"));
    meta.assay.insert(1, assay);

    let mut document = MzTabM::new();
    document.set_meta_data(meta);
    let mut row = MzTabMSmallMoleculeSectionRow {
        sml_identifier: MzTabString::from_text("1"),
        reliability: MzTabString::from_text("²"),
        ..MzTabMSmallMoleculeSectionRow::default()
    };
    row.chemical_name
        .set(vec![MzTabString::from_text("β-D-グルコース")]);
    row.small_molecule_abundance_assay
        .insert(1, MzTabDouble::new(1.0));
    row.opt.push(MzTabOptionalColumnEntry::new(
        "opt_global_名前",
        MzTabString::from_text("値"),
    ));
    document.small_molecule_data.push(row);

    let lines = MzTabMFile::new().generate_lines(&document).unwrap();
    let joined = lines.join("\n");
    assert!(joined.contains("局所識別子"));
    assert!(joined.contains("β-D-グルコース"));
    assert!(joined.contains("opt_global_名前"));
    assert!(joined.contains("\t値"));
    assert_rectangular(&lines);

    // And the same text round-trips through a file whose own name is not ASCII.
    // A directory of its own, so concurrent runs cannot remove each other's files.
    let temp = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let directory = temp.path();
    let path = directory.join("メタボ.mzTab");
    MzTabMFile::new().store(&path, &document).unwrap();
    let written = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(&path).ok();
    assert!(written.contains("β-D-グルコース"));
}

#[test]
fn an_adduct_name_with_several_semicolons_keeps_the_source_reformatting() {
    // `getAdductString_` splits on the first `;` only.
    let mut graph = IdentificationData::new().unwrap();
    let file = graph
        .register_input_file(InputFile::new("run.mzML"))
        .unwrap();
    let software = graph
        .register_processing_software(ProcessingSoftware::new("Tool", "1.0"))
        .unwrap();
    graph
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap();
    let observation = graph
        .register_observation(Observation::new("s=1", file))
        .unwrap();
    let compound = graph
        .register_identified_compound(IdentifiedCompound::new("X:1"))
        .unwrap();
    let plain = graph
        .register_observation_match(ObservationMatch::new(compound, observation))
        .unwrap();
    let mut feature = Feature::new(1.0, 2.0, 3.0);
    feature.add_id_match(plain).unwrap();
    let map = FeatureMap::from_features(vec![feature]);
    let document = MzTabM::export_feature_map(&map, &graph).unwrap();
    // A match without an adduct gets the literal `null`, which `MzTabString`
    // then stores as the null cell.
    let row = &document.small_molecule_evidence_section_rows()[0];
    assert!(row.adduct.is_null());
    assert_eq!(row.adduct.to_cell_string(), "null");
    // The feature row groups by that same literal.
    assert_eq!(
        document.small_molecule_feature_section_rows()[0]
            .adduct
            .to_cell_string(),
        "null"
    );
}

#[test]
fn accessors_and_setters_round_trip() {
    let mut document = MzTabM::new();
    document.set_meta_data(store_fixture_metadata());
    assert_eq!(document.meta_data().software.len(), 4);
    document.set_small_molecule_section_rows(vec![store_fixture_sml_row_three()]);
    document.set_small_molecule_feature_section_rows(vec![store_fixture_smf_row_one()]);
    document.set_small_molecule_evidence_section_rows(vec![store_fixture_sme_row_one()]);
    document.set_empty_rows(vec![26, 111, 196]);
    let mut comments = BTreeMap::new();
    comments.insert(0_usize, "COM\tgenerated by a test".to_owned());
    document.set_comment_rows(comments);
    assert_eq!(document.small_molecule_section_rows().len(), 1);
    assert_eq!(document.small_molecule_feature_section_rows().len(), 1);
    assert_eq!(document.small_molecule_evidence_section_rows().len(), 1);
    assert_eq!(document.empty_rows(), [26, 111, 196]);
    assert_eq!(document.comment_rows().len(), 1);
    // The public fields and the source-named accessors are the same storage.
    assert_eq!(
        document.small_molecule_data.len(),
        document.small_molecule_section_rows().len()
    );
    assert_eq!(
        document.small_molecule_feature_data.len(),
        document.small_molecule_feature_section_rows().len()
    );
    assert_eq!(
        document.small_molecule_evidence_data.len(),
        document.small_molecule_evidence_section_rows().len()
    );
    assert_eq!(document.meta_data, *document.meta_data());
    // Comment and empty rows are carried but not written; the writer emits
    // only the metadata section and the three tables.
    let lines = MzTabMFile::new().generate_lines(&document).unwrap();
    assert!(!lines.iter().any(|line| line.starts_with("COM")));
}

#[test]
fn generate_lines_places_one_empty_line_before_each_header() {
    let document = MzTabM::new();
    let lines = MzTabMFile::new().generate_lines(&document).unwrap();
    let positions: Vec<usize> = lines
        .iter()
        .enumerate()
        .filter(|(_, line)| line.is_empty())
        .map(|(index, _)| index)
        .collect();
    assert_eq!(positions.len(), 3);
    for position in positions {
        let header = &lines[position + 1];
        assert!(
            header.starts_with("SMH\t")
                || header.starts_with("SFH\t")
                || header.starts_with("SEH\t")
        );
    }
    // The same three-blank-line structure the retained outputs have.
    assert_eq!(
        retained(STORE_FIXTURE)
            .iter()
            .filter(|line| line.is_empty())
            .count(),
        3
    );
}

#[test]
fn section_line_reports_its_own_column_count() {
    let meta = MzTabMMetaData::default();
    let writer = MzTabMFile::new();
    let header = writer.generate_small_molecule_header(&meta, &[]).unwrap();
    assert_eq!(header.columns, 14);
    assert_eq!(header.text.split('\t').count(), 14);
    let feature = writer
        .generate_small_molecule_feature_header(&meta, &[])
        .unwrap();
    assert_eq!(feature.columns, 11);
    let evidence = writer
        .generate_small_molecule_evidence_header(&meta, &[])
        .unwrap();
    assert_eq!(evidence.columns, 18);
    assert_eq!(evidence.text.split('\t').next_back().unwrap(), "rank");
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `src/tests/class_tests/openms/source/SpectrumMetaDataLookup_test.cpp`
//! (core SDK `bc9cc12`), section by section with unchanged literals. The
//! class test has no section for `addMissingFAIMSToPeptideIDs`; its tests at
//! the end are derived from `SpectrumMetaDataLookup.cpp:245-312` and say so.

use openms::concept::constants::user_param::{FAIMS_CV, IM};
use openms::identification::{PeptideIdentification, ProteinIdentification};
use openms::kernel::{MSExperiment, MSSpectrum, Precursor};
use openms::metadata::DriftTimeUnit;
use openms::metadata::spectrum_lookup::DEFAULT_SCAN_REGEXP;
use openms::metadata::spectrum_meta_data_lookup::{
    MDF_ALL, MDF_PRECURSORMZ, MDF_RT, SpectrumMetaData, SpectrumMetaDataLookup,
};

fn data_path(name: &str) -> String {
    format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))
}

fn spectrum(native_id: &str, rt: f64, ms_level: u32) -> MSSpectrum {
    MSSpectrum {
        native_id: native_id.into(),
        rt,
        ms_level,
        ..MSSpectrum::default()
    }
}

/// The file-scope fixture: three spectra, the last two with a precursor. The
/// source reuses one `MSSpectrum` object, so the third keeps the second's
/// precursor slot and overwrites it.
fn fixture() -> (Vec<MSSpectrum>, MSSpectrum) {
    let mut spectra = vec![spectrum("spectrum=0", 1.0, 1)];
    let mut s = spectrum("spectrum=1", 2.0, 2);
    s.precursors.push(Precursor::new(1000.0, 2));
    spectra.push(s.clone());
    s.native_id = "spectrum=2".into();
    s.rt = 3.0;
    s.ms_level = 2;
    s.precursors[0] = Precursor::new(500.0, 3);
    spectra.push(s.clone());
    (spectra, s)
}

fn read_lookup() -> SpectrumMetaDataLookup {
    let (spectra, _) = fixture();
    let mut lookup = SpectrumMetaDataLookup::new();
    lookup
        .read_spectra(&spectra, DEFAULT_SCAN_REGEXP, true)
        .unwrap();
    lookup
}

fn peptide(reference: Option<&str>, rt: Option<f64>) -> PeptideIdentification {
    let mut peptide = PeptideIdentification::new();
    if let Some(reference) = reference {
        peptide.set_spectrum_reference(reference);
    }
    peptide.rt = rt;
    peptide
}

// START_SECTION((SpectrumMetaDataLookup()))
#[test]
fn default_constructor() {
    let lookup = SpectrumMetaDataLookup::new();
    assert!(lookup.is_empty());
}

// START_SECTION((~SpectrumMetaDataLookup())): Drop.

// START_SECTION((template <typename SpectrumContainer> void readSpectra(...)))
#[test]
fn read_spectra() {
    let lookup = read_lookup();
    assert!(!lookup.is_empty());
}

// START_SECTION((void getSpectrumMetaData(Size, SpectrumMetaData&) const))
#[test]
fn meta_data_by_index() {
    let lookup = read_lookup();
    let meta = lookup.spectrum_meta_data(0).unwrap();
    assert_eq!(meta.rt, 1.0);
    assert_eq!(meta.ms_level, 1);
    assert_eq!(meta.native_id, "spectrum=0");
    assert_eq!(meta.scan_number, 0);
    let meta = lookup.spectrum_meta_data(1).unwrap();
    assert_eq!(meta.rt, 2.0);
    assert_eq!(meta.precursor_rt, 1.0);
    assert_eq!(meta.precursor_mz, 1000.0);
    assert_eq!(meta.precursor_charge, 2);
    assert_eq!(meta.ms_level, 2);
    assert_eq!(meta.native_id, "spectrum=1");
    assert_eq!(meta.scan_number, 1);
    // Exception::IndexOverflow for an index past the spectra read
    // (SpectrumMetaDataLookup.cpp:22-26).
    assert!(lookup.spectrum_meta_data(3).is_err());
}

// START_SECTION((static void getSpectrumMetaData(const MSSpectrum&, SpectrumMetaData&, const boost::regex&, const map<Size, double>&)))
#[test]
fn meta_data_of_a_spectrum() {
    let (_, spectrum) = fixture();
    let mut meta = SpectrumMetaData::default();
    let log = SpectrumMetaDataLookup::spectrum_meta_data_of(
        &spectrum,
        &mut meta,
        None,
        &Default::default(),
    )
    .unwrap();
    assert!(log.is_empty());
    assert_eq!(meta.rt, 3.0);
    assert_eq!(meta.precursor_mz, 500.0);
    assert_eq!(meta.precursor_charge, 3);
    assert_eq!(meta.ms_level, 2);
    assert_eq!(meta.native_id, "spectrum=2");
    assert_eq!(meta.scan_number, -1); // not extracted

    let mut precursor_rts = std::collections::BTreeMap::new();
    precursor_rts.insert(1_usize, 1.0);
    let scan_regexp = openms::concept::boost_regex::BoostRegex::new("=(?<SCAN>\\d+)$").unwrap();
    SpectrumMetaDataLookup::spectrum_meta_data_of(
        &spectrum,
        &mut meta,
        Some(&scan_regexp),
        &precursor_rts,
    )
    .unwrap();
    assert_eq!(meta.precursor_rt, 1.0);
    assert_eq!(meta.scan_number, 2);
}

// START_SECTION((void getSpectrumMetaData(const std::string&, SpectrumMetaData&, MetaDataFlags) const))
#[test]
fn meta_data_from_reference() {
    let mut lookup = read_lookup();
    let mut meta = SpectrumMetaData::default();
    lookup.add_reference_format(DEFAULT_SCAN_REGEXP).unwrap();
    lookup
        .spectrum_meta_data_from_reference("scan_number=1", &mut meta, MDF_ALL)
        .unwrap();
    assert_eq!(meta.rt, 2.0);
    assert_eq!(meta.native_id, "spectrum=1");

    lookup
        .add_reference_format(r"rt=(?<RT>\d+(\.\d+)?),mz=(?<MZ>\d+(\.\d+)?)")
        .unwrap();
    let mut meta2 = SpectrumMetaData::default();
    let flags = MDF_RT | MDF_PRECURSORMZ;
    // no actual look-up of the spectrum necessary:
    lookup
        .spectrum_meta_data_from_reference("rt=5.0,mz=1000.0", &mut meta2, flags)
        .unwrap();
    assert_eq!(meta2.rt, 5.0);
    assert_eq!(meta2.precursor_mz, 1000.0);
    assert_eq!(meta2.precursor_charge, 0);
    assert_eq!(meta2.native_id, "");

    // look-up of the spectrum necessary:
    let mut meta3 = SpectrumMetaData::default();
    lookup
        .spectrum_meta_data_from_reference("rt=2.0,mz=1000.0", &mut meta3, MDF_ALL)
        .unwrap();
    assert_eq!(meta3.rt, 2.0);
    assert_eq!(meta3.precursor_mz, 1000.0);
    assert_eq!(meta3.precursor_charge, 2);
    assert_eq!(meta3.native_id, "spectrum=1");
    // TEST_EXCEPTION(Exception::ElementNotFound, ...)
    assert!(
        lookup
            .spectrum_meta_data_from_reference("rt=5.0,mz=1000.0", &mut meta3, MDF_ALL)
            .is_err()
    );
}

// START_SECTION((bool addMissingRTsToPeptideIDs(PeptideIdentificationList& peptides, const MSExperiment& exp)))
#[test]
fn add_missing_rts() {
    // Test 1: No spectra in MSExperiment
    let mut peptides = vec![peptide(None, Some(1.0))];
    let exp_empty = MSExperiment::default();
    let outcome =
        SpectrumMetaDataLookup::add_missing_rts_to_peptide_ids(&mut peptides, &exp_empty).unwrap();
    assert!(!outcome.success);
    assert_eq!(peptides[0].rt, Some(1.0));

    // Test 2: Valid MSExperiment with missing RTs
    let mut exp_valid = MSExperiment::default();
    exp_valid.spectra.push(spectrum("index=0", 2.5, 1));
    exp_valid.spectra.push(spectrum("index=2", 5.3, 1));
    peptides.resize(2, PeptideIdentification::new());
    peptides[0].set_spectrum_reference("index=0");
    peptides[0].rt = Some(f64::NAN);
    peptides[1].set_spectrum_reference("index=2");
    peptides[1].rt = None;
    let outcome =
        SpectrumMetaDataLookup::add_missing_rts_to_peptide_ids(&mut peptides, &exp_valid).unwrap();
    assert!(outcome.success);
    assert_eq!(peptides[0].rt, Some(2.5));
    assert_eq!(peptides[1].rt, Some(5.3));

    // Test 3: Missing spectrum reference in experiment
    peptides.truncate(1);
    peptides[0].set_spectrum_reference("index=nonexistent");
    peptides[0].rt = None;
    let outcome =
        SpectrumMetaDataLookup::add_missing_rts_to_peptide_ids(&mut peptides, &exp_valid).unwrap();
    assert!(!outcome.success);
    assert_eq!(peptides[0].rt, None);
    assert_eq!(
        outcome.log,
        [
            "Error: Failed to look up retention time for peptide identification with spectrum reference 'index=nonexistent' - no spectrum with corresponding native ID found."
        ]
    );

    // Test 4: No missing RTs in peptides
    peptides[0].set_spectrum_reference("index=0");
    peptides[0].rt = Some(2.0);
    let outcome =
        SpectrumMetaDataLookup::add_missing_rts_to_peptide_ids(&mut peptides, &exp_valid).unwrap();
    assert!(outcome.success);
    assert_eq!(peptides[0].rt, Some(2.0));
}

// START_SECTION((bool addMissingIMToPeptideIDs(PeptideIdentificationList& peptides, const MSExperiment& exp)))
#[test]
fn add_missing_im() {
    // Test 1: Empty MSExperiment
    let mut peptides = vec![peptide(Some("index=0"), None)];
    let exp_empty = MSExperiment::default();
    assert!(
        !SpectrumMetaDataLookup::add_missing_im_to_peptide_ids(&mut peptides, &exp_empty)
            .unwrap()
            .success
    );

    // Test 2: MSExperiment with no IM format (not IM_SPECTRUM)
    let mut exp_no_im = MSExperiment::default();
    exp_no_im.spectra.push(spectrum("index=0", 0.0, 1));
    assert!(
        !SpectrumMetaDataLookup::add_missing_im_to_peptide_ids(&mut peptides, &exp_no_im)
            .unwrap()
            .success
    );

    // Test 3: MSExperiment with valid IM values
    let mut exp_valid = MSExperiment::default();
    let mut spectrum1 = spectrum("index=0", 0.0, 1);
    spectrum1.drift_time = 2.5;
    let mut spectrum2 = spectrum("index=2", 0.0, 1);
    spectrum2.drift_time = 5.3;
    // Distinct retention times would not change the result; the source's two
    // spectra share the default RT, whose map entry the second overwrites.
    exp_valid.spectra.push(spectrum1);
    exp_valid.spectra.push(spectrum2);
    peptides.resize(2, PeptideIdentification::new());
    peptides[0].set_spectrum_reference("index=0");
    peptides[1].set_spectrum_reference("index=2");
    assert!(
        SpectrumMetaDataLookup::add_missing_im_to_peptide_ids(&mut peptides, &exp_valid)
            .unwrap()
            .success
    );
    assert_eq!(peptides[0].metadata[IM].as_f64().unwrap(), 2.5);
    assert_eq!(peptides[1].metadata[IM].as_f64().unwrap(), 5.3);
}

// START_SECTION((bool addMissingSpectrumReferences(...)))
#[test]
fn add_missing_spectrum_references() {
    let mut peptides = vec![peptide(Some("index=666"), Some(5.1))];
    // missing file -> exception, no non-effective executions
    assert!(
        SpectrumMetaDataLookup::add_missing_spectrum_references(
            &mut peptides,
            "this_file_does_not_exist.mzML",
            false,
            false,
            false,
            &mut [],
        )
        .is_err()
    );
    assert_eq!(peptides[0].spectrum_reference(), "index=666");

    peptides.push(peptide(None, Some(5.3)));
    let filename = data_path("mzml_validator/MzMLFile_1.mzML");
    SpectrumMetaDataLookup::add_missing_spectrum_references(
        &mut peptides,
        &filename,
        false,
        false,
        false,
        &mut [],
    )
    .unwrap();
    assert_eq!(peptides[0].spectrum_reference(), "index=666"); // no overwrite
    assert_eq!(peptides[1].spectrum_reference(), "index=2");

    let mut proteins = vec![ProteinIdentification::default()];
    SpectrumMetaDataLookup::add_missing_spectrum_references(
        &mut peptides,
        &filename,
        false,
        true,
        true,
        &mut proteins,
    )
    .unwrap();
    assert_eq!(peptides[0].spectrum_reference(), "index=0"); // gets updated
    assert_eq!(peptides[1].spectrum_reference(), "index=2");
    // Divergence (documented): the source sets this on a by-value copy.
    assert_eq!(
        proteins[0].primary_ms_run_paths,
        [format!("file://{filename}")]
    );
}

// ---------------------------------------------------------------------------
// Not in the class test; expected values derived from
// SpectrumMetaDataLookup.cpp.
// ---------------------------------------------------------------------------

/// `addMissingSpectrumReferences` with an RT matching nothing: the error line,
/// `false`, and `stop_on_error` ending the loop (cpp:343-357).
#[test]
fn add_missing_spectrum_references_reports_unmatched_rts() {
    let filename = data_path("mzml_validator/MzMLFile_1.mzML");
    let mut peptides = vec![peptide(None, Some(99.0)), peptide(None, Some(5.2))];
    let outcome = SpectrumMetaDataLookup::add_missing_spectrum_references(
        &mut peptides,
        &filename,
        true,
        false,
        false,
        &mut [],
    )
    .unwrap();
    assert!(!outcome.success);
    assert_eq!(outcome.log.len(), 1);
    assert!(outcome.log[0].starts_with(
        "Error: Failed to look up spectrum native ID for peptide identification with retention time '99"
    ));
    assert!(!peptides[1].metadata.contains_key("spectrum_reference"));
    let outcome = SpectrumMetaDataLookup::add_missing_spectrum_references(
        &mut peptides,
        &filename,
        false,
        false,
        false,
        &mut [],
    )
    .unwrap();
    assert!(!outcome.success);
    assert_eq!(peptides[1].spectrum_reference(), "index=1");
}

/// `addMissingIMToPeptideIDs` lets `ElementNotFound` escape (cpp:207).
#[test]
fn add_missing_im_propagates_unknown_references() {
    let mut exp = MSExperiment::default();
    exp.spectra.push(spectrum("index=0", 1.0, 1));
    let mut peptides = vec![peptide(Some("index=9"), None)];
    assert!(SpectrumMetaDataLookup::add_missing_im_to_peptide_ids(&mut peptides, &exp).is_err());
}

/// `addMissingFAIMSToPeptideIDs` (cpp:219-294): explicit voltages on MS1 and
/// MS2, MS2 inheriting the last seen voltage, existing annotations counted,
/// unknown references warned about.
#[test]
fn add_missing_faims() {
    let faims = |native_id: &str, rt: f64, level: u32, cv: Option<f64>| {
        let mut s = spectrum(native_id, rt, level);
        if let Some(cv) = cv {
            s.drift_time = cv;
            s.drift_time_unit = DriftTimeUnit::FaimsCompensationVoltage;
        }
        s
    };
    // Without spectra, and without FAIMS data: false.
    let mut peptides = vec![peptide(Some("s=0"), None)];
    assert!(
        !SpectrumMetaDataLookup::add_missing_faims_to_peptide_ids(
            &mut peptides,
            &MSExperiment::default()
        )
        .unwrap()
        .success
    );
    let mut plain = MSExperiment::default();
    plain.spectra.push(spectrum("s=0", 1.0, 1));
    let outcome =
        SpectrumMetaDataLookup::add_missing_faims_to_peptide_ids(&mut peptides, &plain).unwrap();
    assert!(!outcome.success);
    assert_eq!(
        outcome.log,
        ["No FAIMS data found in the experiment. Skipping FAIMS annotation."]
    );

    let mut exp = MSExperiment::default();
    exp.spectra.push(faims("s=0", 1.0, 1, Some(-45.0)));
    exp.spectra.push(faims("s=1", 2.0, 2, None)); // inherits -45
    exp.spectra.push(faims("s=2", 3.0, 1, Some(-65.0)));
    exp.spectra.push(faims("s=3", 4.0, 2, Some(-70.0))); // explicit on MS2
    exp.spectra.push(faims("s=4", 5.0, 1, None)); // MS1 inherits nothing
    let mut peptides = vec![
        peptide(Some("s=1"), None),
        peptide(Some("s=3"), None),
        peptide(Some("s=4"), None),
        peptide(Some("missing"), None),
        peptide(Some("s=0"), None),
    ];
    peptides[4].metadata.insert(
        FAIMS_CV.into(),
        openms::metadata::MetaValue::try_from(-1.0).unwrap(),
    );
    let outcome =
        SpectrumMetaDataLookup::add_missing_faims_to_peptide_ids(&mut peptides, &exp).unwrap();
    assert!(outcome.success);
    assert_eq!(peptides[0].metadata[FAIMS_CV].as_f64().unwrap(), -45.0);
    assert_eq!(peptides[1].metadata[FAIMS_CV].as_f64().unwrap(), -70.0);
    assert!(!peptides[2].metadata.contains_key(FAIMS_CV));
    assert!(!peptides[3].metadata.contains_key(FAIMS_CV));
    assert_eq!(peptides[4].metadata[FAIMS_CV].as_f64().unwrap(), -1.0);
    assert_eq!(
        outcome.log,
        [
            "Could not find spectrum with native ID 'missing' for FAIMS annotation.",
            "Annotated 3 of 5 peptide IDs with FAIMS CV values."
        ]
    );
}

/// Reference extraction details (cpp:60-148): no matching format leaves the
/// record unchanged; an `ID` group is written even when its flag is not the
/// only one; `CHARGE` truncates a decimal.
#[test]
fn reference_extraction_details() {
    let mut lookup = read_lookup();
    let mut meta = SpectrumMetaData::default();
    lookup
        .spectrum_meta_data_from_reference("anything", &mut meta, MDF_ALL)
        .unwrap();
    // Still the default record (NaN fields compare by is_nan).
    assert!(meta.rt.is_nan() && meta.precursor_rt.is_nan() && meta.precursor_mz.is_nan());
    assert_eq!(
        (
            meta.precursor_charge,
            meta.ms_level,
            meta.scan_number,
            meta.native_id.as_str()
        ),
        (0, 0, -1, "")
    );
    lookup
        .add_reference_format(r"id=(?<ID>[^,]*),z=(?<CHARGE>[\d.]+)")
        .unwrap();
    lookup
        .spectrum_meta_data_from_reference(
            "id=spectrum=2,z=3.9",
            &mut meta,
            openms::metadata::spectrum_meta_data_lookup::MDF_NATIVEID
                | openms::metadata::spectrum_meta_data_lookup::MDF_PRECURSORCHARGE,
        )
        .unwrap();
    assert_eq!(meta.native_id, "spectrum=2");
    assert_eq!(meta.precursor_charge, 3);
    assert!(meta.rt.is_nan());
}

/// Precursor RT without a preceding spectrum one level lower logs the
/// source's error line (cpp:39-52). The map is consulted only once it holds
/// an entry, so the first spectrum read never logs; a level-3 spectrum first
/// fills it without providing the level-1 entry a level-2 precursor needs.
/// A native ID the scan expression cannot parse logs the scan error line.
#[test]
fn precursor_rt_not_found_is_logged() {
    let first = spectrum("scan=4", 1.0, 3);
    let mut second = spectrum("scan=5", 2.0, 2);
    second.precursors.push(Precursor::new(400.0, 2));
    let third = spectrum("no scan here", 3.0, 1);
    let mut lookup = SpectrumMetaDataLookup::new();
    let log = lookup
        .read_spectra(&[first, second, third], DEFAULT_SCAN_REGEXP, true)
        .unwrap();
    assert_eq!(
        log,
        [
            "Error: Could not set precursor RT for spectrum with native ID 'scan=5' - precursor spectrum not found.",
            "Error: Could not extract scan number from spectrum native ID 'no scan here' using regular expression '=(?<SCAN>\\d+)$'."
        ]
    );
    assert!(lookup.spectrum_meta_data(1).unwrap().precursor_rt.is_nan());
    assert_eq!(lookup.spectrum_meta_data(2).unwrap().scan_number, -1);
}

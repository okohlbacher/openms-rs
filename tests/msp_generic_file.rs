// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/msp_generic.rs`, the port of
//! `FORMAT/MSPGenericFile.h`.
//!
//! Every `START_SECTION` of the upstream `MSPGenericFile_test.cpp` (core
//! `bc9cc12`) is reproduced with its literals, including the friend-class
//! section for the private `addSpectrumToLibrary`. The fixture is the upstream
//! class-test file `MSPGenericFile_input.msp` (CRLF line ends), recorded in
//! `tests/data/msp_generic_file_provenance.json`.

use openms::Error;
use openms::format::msp_generic::{IS_VALID, MSPGenericFile, NumberFormat, StoreOptions};
use openms::kernel::{MSExperiment, MSSpectrum, Peak1D};
use openms::metadata::{MetaValue, MetaValueData};
use openms::param::ParamValue;
use openms::system::file::TempDir;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/msp_generic_file")
        .join(name)
}

fn meta(spectrum: &MSSpectrum, key: &str) -> String {
    spectrum.metadata[key].to_string()
}

/// The source's `TEST_REAL_SIMILAR` with its default tolerances.
fn similar(a: f64, b: f64) -> bool {
    if (a - b).abs() <= 1e-5 {
        return true;
    }
    if a == 0.0 || b == 0.0 {
        return false;
    }
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    (1.0 / (1.0 + 1e-5)..=1.0 + 1e-5).contains(&ratio)
}

// START_SECTION(MSPGenericFile()) / START_SECTION(~MSPGenericFile())
#[test]
fn construction_and_parameters() {
    let file = MSPGenericFile::new();
    assert_eq!(file.synonyms_separator(), "|");
    let defaults = MSPGenericFile::defaults().unwrap();
    assert_eq!(
        defaults.value("synonyms_separator").unwrap(),
        &ParamValue::String("|".into())
    );
    let mut param = defaults.clone();
    param
        .set_value(
            "synonyms_separator",
            ParamValue::String("/".into()),
            "",
            &[],
        )
        .unwrap();
    let file = MSPGenericFile::from_param(&param).unwrap();
    assert_eq!(file.synonyms_separator(), "/");
    assert_eq!(
        file.to_param()
            .unwrap()
            .value("synonyms_separator")
            .unwrap(),
        &ParamValue::String("/".into())
    );
    drop(file);
}

// START_SECTION(void load(const std::string& filename, MSExperiment& experiment) const)
#[test]
fn load_matches_the_source_section() {
    let mut msp = MSPGenericFile::new();
    let experiment = msp.load(data("MSPGenericFile_input.msp")).unwrap();
    let spectra = &experiment.spectra;
    assert_eq!(spectra.len(), 3);

    let s1 = &spectra[0];
    assert_eq!(s1.peaks.len(), 14);
    assert_eq!(s1.name, "name1 of first");
    assert_eq!(meta(s1, "Synon"), "name2 of 1st|name3 of firsttt");
    assert_eq!(meta(s1, "Formula"), "A11B22C333");
    assert_eq!(meta(s1, "MW"), "156");
    assert_eq!(meta(s1, "CAS#"), "0123-45-6");
    assert_eq!(meta(s1, "NIST#"), "654321");
    assert_eq!(meta(s1, "DB#"), "1");
    assert_eq!(meta(s1, "Comments"), "Some comment");
    assert_eq!(meta(s1, "Num Peaks"), "14");
    assert_eq!(s1.peaks[0], Peak1D::new(27.0, 29.0));
    assert_eq!(s1.peaks[5], Peak1D::new(60.0, 41.0));
    assert_eq!(s1.peaks[10], Peak1D::new(90.0, 168.0));
    assert_eq!(s1.peaks[13], Peak1D::new(105.0, 36.0));

    let s2 = &spectra[1];
    assert_eq!(s2.peaks.len(), 15);
    assert_eq!(s2.name, "name1 of second");
    assert_eq!(meta(s2, "Synon"), "name2 of 2nd|name3 of seconddd");
    assert_eq!(meta(s2, "Formula"), "A44B55C666");
    assert_eq!(meta(s2, "MW"), "589");
    assert_eq!(meta(s2, "CAS#"), "3210-45-6");
    assert_eq!(meta(s2, "NIST#"), "789564");
    assert_eq!(meta(s2, "DB#"), "2");
    assert_eq!(meta(s2, "Comments"), "Some other comment");
    assert_eq!(meta(s2, "Num Peaks"), "15");
    assert_eq!(s2.peaks[0], Peak1D::new(27.0, 29.0));
    assert_eq!(s2.peaks[5], Peak1D::new(260.0, 41.0));
    assert_eq!(s2.peaks[10], Peak1D::new(290.0, 168.0));
    assert_eq!(s2.peaks[14], Peak1D::new(310.0, 20.0));

    let s3 = &spectra[2];
    assert_eq!(s3.peaks.len(), 16);
    assert_eq!(s3.name, "name1 of third");
    assert_eq!(meta(s3, "Synon"), "name2 of 3rd|name3 of thirddd");
    assert_eq!(meta(s3, "Formula"), "A12B12C123");
    assert_eq!(meta(s3, "MW"), "562");
    assert_eq!(meta(s3, "CAS#"), "4210-47-4");
    assert_eq!(meta(s3, "NIST#"), "749514");
    assert_eq!(meta(s3, "DB#"), "3");
    assert!(!s3.metadata.contains_key("Comments")); // this spectrum has no comment
    assert_eq!(meta(s3, "Num Peaks"), "16");
    assert_eq!(s3.peaks[0], Peak1D::new(27.0, 29.0));
    assert_eq!(s3.peaks[5], Peak1D::new(260.0, 41.0));
    assert_eq!(s3.peaks[10], Peak1D::new(290.0, 168.0));
    assert_eq!(s3.peaks[14], Peak1D::new(310.0, 20.0));
    assert_eq!(s3.peaks[15], Peak1D::new(111.0, 44.0));

    // CCS metadata: a plain number is stored as a float under MSM_CCS.
    let ccs = msp
        .read("Name: ccs_test\nCCS: 123.45\nNum Peaks: 1\n100:100\n".as_bytes())
        .unwrap();
    assert_eq!(ccs.spectra.len(), 1);
    match ccs.spectra[0].metadata["CCS"].data() {
        MetaValueData::Float(value) => assert!(similar(*value, 123.45)),
        other => panic!("{other:?}"),
    }
    // A non-numeric value is kept verbatim as a string.
    let ccs = msp
        .read("Name: ccs_text\nCCS: not_a_number\nNum Peaks: 1\n100:100\n".as_bytes())
        .unwrap();
    assert_eq!(ccs.spectra.len(), 1);
    assert_eq!(meta(&ccs.spectra[0], "CCS"), "not_a_number");
}

fn spectrum(name: &str, rt: f64, peaks: &[(f64, f32)]) -> MSSpectrum {
    let mut spec = MSSpectrum {
        name: name.into(),
        rt,
        ms_level: 1,
        ..MSSpectrum::default()
    };
    for &(mz, intensity) in peaks {
        spec.peaks.push(Peak1D::new(mz, intensity));
    }
    spec
}

fn set(spec: &mut MSSpectrum, key: &str, value: &str) {
    spec.metadata.insert(key.into(), MetaValue::from(value));
}

// START_SECTION(void store(const std::string& filename, const MSExperiment& library) const)
#[test]
fn store_matches_the_source_section() {
    for options in [StoreOptions::default(), StoreOptions::source()] {
        let msp = MSPGenericFile::new();
        let mut exp = MSExperiment::new();

        let mut spec = spectrum(
            "first spectrum",
            11.1,
            &[(1.0, 1.50), (2.0, 2.5), (3.0, 3.5)],
        );
        set(&mut spec, "Synon", "first1|first2|first3");
        set(&mut spec, "CAS#", "0123-45-6");
        set(&mut spec, "NIST#", "654321");
        exp.spectra.push(spec);

        let mut spec = spectrum(
            "second spectrum",
            22.2,
            &[
                (11.0, 11.50),
                (12.0, 12.5),
                (13.0, 13.5),
                (14.0, 14.5),
                (15.0, 15.5),
                (16.0, 16.5),
            ],
        );
        set(&mut spec, "Synon", "second1");
        set(&mut spec, "CAS#", "0123-45-2");
        set(&mut spec, "NIST#", "654322");
        set(&mut spec, "other_metadata1", "value1");
        set(&mut spec, "other_metadata2", "value2");
        exp.spectra.push(spec);

        let mut spec = spectrum(
            "third spectrum",
            33.3,
            &[(101.0, 101.50), (102.0, 102.5), (103.0, 103.5)],
        );
        set(&mut spec, "CAS#", "0123-45-3");
        set(&mut spec, "NIST#", "654323");
        exp.spectra.push(spec);

        let dir = TempDir::new(false).unwrap();
        let output_filepath = dir.path().join("out.msp");
        msp.store_with_options(&output_filepath, &exp, options)
            .unwrap();

        // read back created file
        let exp_test = MSPGenericFile::new().load(&output_filepath).unwrap();
        let spectra = &exp_test.spectra;
        assert_eq!(spectra.len(), 3);

        let s1 = &spectra[0];
        assert_eq!(s1.name, "first spectrum");
        assert_eq!(meta(s1, "Synon"), "first1|first2|first3");
        assert_eq!(meta(s1, "CAS#"), "0123-45-6");
        assert_eq!(meta(s1, "NIST#"), "654321");
        assert_eq!(s1.peaks.len(), 3);
        for (peak, (mz, intensity)) in s1.peaks.iter().zip([(1.0, 1.5), (2.0, 2.5), (3.0, 3.5)]) {
            assert!(similar(peak.mz, mz));
            assert!(similar(f64::from(peak.intensity), intensity));
        }

        let s2 = &spectra[1];
        assert_eq!(s2.name, "second spectrum");
        assert_eq!(meta(s2, "Synon"), "second1");
        assert_eq!(meta(s2, "CAS#"), "0123-45-2");
        assert_eq!(meta(s2, "NIST#"), "654322");
        assert_eq!(s2.peaks.len(), 6);
        for (peak, (mz, intensity)) in s2.peaks.iter().zip([
            (11.0, 11.5),
            (12.0, 12.5),
            (13.0, 13.5),
            (14.0, 14.5),
            (15.0, 15.5),
            (16.0, 16.5),
        ]) {
            assert!(similar(peak.mz, mz));
            assert!(similar(f64::from(peak.intensity), intensity));
        }

        let s3 = &spectra[2];
        assert_eq!(s3.name, "third spectrum");
        assert!(!s3.metadata.contains_key("Synon"));
        assert_eq!(meta(s3, "CAS#"), "0123-45-3");
        assert_eq!(meta(s3, "NIST#"), "654323");
        assert_eq!(s3.peaks.len(), 3);
        for (peak, (mz, intensity)) in
            s3.peaks
                .iter()
                .zip([(101.0, 101.5), (102.0, 102.5), (103.0, 103.5)])
        {
            assert!(similar(peak.mz, mz));
            assert!(similar(f64::from(peak.intensity), intensity));
        }

        // test invalid spectrum (no name).
        let mut invalid_exp = MSExperiment::new();
        let mut invalid_spec = spectrum("", 11.1, &[(1.0, 1.50), (2.0, 2.5), (3.0, 3.5)]);
        set(&mut invalid_spec, "Synon", "first1|first2|first3");
        set(&mut invalid_spec, "CAS#", "0123-45-6");
        set(&mut invalid_spec, "NIST#", "654321");
        invalid_exp.spectra.push(invalid_spec);
        let invalid_path = dir.path().join("invalid.msp");
        assert!(matches!(
            msp.store_with_options(&invalid_path, &invalid_exp, options),
            Err(Error::MissingInformation(_))
        ));
    }
}

// START_SECTION(void addSpectrumToLibrary(MSSpectrum& spectrum, MSExperiment& library))
#[test]
fn add_spectrum_to_library_matches_the_source_section() {
    let mut msp = MSPGenericFile::new();
    let mut lib = MSExperiment::new();

    let mut spec = MSSpectrum {
        name: String::new(), // empty name
        ..MSSpectrum::default()
    };
    spec.metadata.insert(IS_VALID.into(), MetaValue::from(1i64));
    assert!(matches!(
        msp.add_spectrum_to_library(&mut spec, &mut lib),
        Err(Error::MissingInformation(_))
    ));
    assert_eq!(lib.len(), 0);

    spec.name = "foo".into(); // Num Peaks still absent!
    assert!(matches!(
        msp.add_spectrum_to_library(&mut spec, &mut lib),
        Err(Error::MissingInformation(_))
    ));
    assert_eq!(lib.len(), 0);

    set(&mut spec, "Num Peaks", "2");
    // Num Peaks is set but raw data points have not been added
    assert!(matches!(
        msp.add_spectrum_to_library(&mut spec, &mut lib),
        Err(Error::Parse { .. })
    ));
    assert_eq!(lib.len(), 0);

    spec.peaks.push(Peak1D::new(1.0, 2.0));
    spec.peaks.push(Peak1D::new(3.0, 4.0)); // now the spectrum is valid
    msp.add_spectrum_to_library(&mut spec, &mut lib).unwrap();
    assert_eq!(lib.len(), 1);

    spec.name = "bar".into();
    spec.metadata.insert(IS_VALID.into(), MetaValue::from(1i64));
    msp.add_spectrum_to_library(&mut spec, &mut lib).unwrap();
    assert_eq!(lib.len(), 2);

    spec.metadata.insert(IS_VALID.into(), MetaValue::from(1i64));
    msp.add_spectrum_to_library(&mut spec, &mut lib).unwrap(); // duplicate, won't be added
    assert_eq!(lib.len(), 2);

    spec.metadata.insert(IS_VALID.into(), MetaValue::from(0i64));
    spec.name = "not a duplicate".into();
    msp.add_spectrum_to_library(&mut spec, &mut lib).unwrap();
    assert_eq!(lib.len(), 2);

    // Native observation of the source's in-place update: the first add set
    // the caller's retention time from -1 to 0 (the library size then), so the
    // second copy kept 0 rather than taking index 1.
    assert_eq!(lib.spectra[0].rt, 0.0);
    assert_eq!(lib.spectra[1].rt, 0.0);
    assert!(!lib.spectra[0].metadata.contains_key(IS_VALID));
}

// ----------------------------------------------------------------------------
// Native regressions derived from MSPGenericFile.cpp
// ----------------------------------------------------------------------------

fn read(text: &str) -> openms::Result<MSExperiment> {
    MSPGenericFile::new().read(text.as_bytes())
}

#[test]
fn line_order_and_special_fields() {
    let lib = read(
        "Comment before any name: dropped\n\
         12 5\n\
         NAME: A\n\
         SYNONYMS: first\n\
         synon: second\n\
         RETENTIONTIME: 7.5 min\n\
         PRECURSORMZ: 301.2\n\
         INCHIKEY: KEY\n\
         SMILES: CCO\r\n\
         FORMULA: C2H6O\n\
         PRECURSORTYPE: [M+H]+\n\
         Num peaks: 2\n\
         a: b: c\n\
         1.5:2 3\t4;\n",
    )
    .unwrap();
    assert_eq!(lib.spectra.len(), 1);
    let s = &lib.spectra[0];
    assert_eq!(s.name, "A");
    assert_eq!(meta(s, "Metabolite_Name"), "A");
    assert_eq!(meta(s, "Synon"), "first|second");
    assert_eq!(s.rt, 7.5);
    assert_eq!(s.precursors.len(), 1);
    assert_eq!(s.precursors[0].mz, 301.2);
    assert_eq!(meta(s, "Inchi_String"), "KEY");
    // Without no_mod_s, `.` also takes the CR before the LF.
    assert_eq!(meta(s, "SMILES_String"), "CCO\r");
    assert_eq!(meta(s, "Sum_Formula"), "C2H6O");
    assert_eq!(meta(s, "Precursor_Ion"), "[M+H]+");
    assert_eq!(meta(s, "a: b"), "c");
    assert_eq!(s.peaks, vec![Peak1D::new(1.5, 2.0), Peak1D::new(3.0, 4.0)]);
    assert!(!s.metadata.contains_key("Comment before any name"));
}

#[test]
fn duplicates_missing_retention_times_and_errors() {
    let lib =
        read("Name: A\nNum Peaks: 1\n1 1\nName: A\nNum Peaks: 1\n2 2\nName: B\nNum Peaks: 0\n")
            .unwrap();
    assert_eq!(lib.spectra.len(), 2);
    assert_eq!(lib.spectra[0].peaks[0].mz, 1.0);
    // Retention times default to the library index.
    assert_eq!(lib.spectra[0].rt, 0.0);
    assert_eq!(lib.spectra[1].rt, 1.0);
    // A point count mismatch aborts the load.
    assert!(matches!(
        read("Name: A\nNum Peaks: 2\n1 1\n"),
        Err(Error::Parse { .. })
    ));
    // A peak line without a point is the source's std::stod failure.
    assert!(matches!(
        read("Name: A\nNum Peaks: 0\n5\n"),
        Err(Error::Parse { line: 3, .. })
    ));
    assert!(matches!(
        read("Name: A\nRetention Time: soon\nNum Peaks: 0\n"),
        Err(Error::Parse { line: 2, .. })
    ));
}

#[test]
fn round_trip_default_is_exact_and_source_format_is_six_digits() {
    let mut spec = spectrum("x", 1234.5678901, &[(123.456789012, 1_500_000.0)]);
    set(&mut spec, "Num Peaks", "ignored on store");
    let mut exp = MSExperiment::new();
    exp.spectra.push(spec);
    let msp = MSPGenericFile::new();

    let mut native = Vec::new();
    msp.write_with_options(&mut native, &exp, StoreOptions::default())
        .unwrap();
    let native = String::from_utf8(native).unwrap();
    assert_eq!(
        native,
        "Name: x\nRetention Time: 1234.5678901\nNum Peaks: 1\n123.456789012:1500000 \n\n"
    );
    let back = read(&native).unwrap();
    assert_eq!(back.spectra[0].rt, 1234.5678901);
    assert_eq!(
        back.spectra[0].peaks,
        vec![Peak1D::new(123.456789012, 1_500_000.0)]
    );

    let mut source = Vec::new();
    msp.write_with_options(
        &mut source,
        &exp,
        StoreOptions {
            numbers: NumberFormat::Source,
        },
    )
    .unwrap();
    assert_eq!(
        String::from_utf8(source).unwrap(),
        "Name: x\nRetention Time: 1234.57\nNum Peaks: 1\n123.457:1.5e+06 \n\n"
    );
}

#[test]
fn separators_split_synonyms_on_store_and_join_on_load() {
    let mut param = MSPGenericFile::defaults().unwrap();
    param
        .set_value(
            "synonyms_separator",
            ParamValue::String("||".into()),
            "",
            &[],
        )
        .unwrap();
    let mut msp = MSPGenericFile::from_param(&param).unwrap();
    let lib = msp
        .read("Name: A\nSynon: s1\nSynon: s2\nNum Peaks: 1\n1 1\n".as_bytes())
        .unwrap();
    assert_eq!(meta(&lib.spectra[0], "Synon"), "s1||s2");
    let mut out = Vec::new();
    msp.write_with_options(&mut out, &lib, StoreOptions::default())
        .unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(text.contains("Synon: s1\nSynon: s2\n"), "{text}");
}

#[test]
fn line_breaks_in_stored_text_are_refused() {
    let mut spec = spectrum("x", 1.0, &[(1.0, 1.0)]);
    set(&mut spec, "note", "two\nlines");
    let mut exp = MSExperiment::new();
    exp.spectra.push(spec);
    assert!(matches!(
        MSPGenericFile::new().write_with_options(Vec::new(), &exp, StoreOptions::default()),
        Err(Error::InvalidValue(_))
    ));
}

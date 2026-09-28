// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tests for `src/format/flash_deconv_spectrum_file.rs` and
//! `src/format/flash_deconv_feature_file.rs`, the header writers of
//! `FORMAT/FLASHDeconvSpectrumFile.h` and `FORMAT/FLASHDeconvFeatureFile.h`.
//!
//! The header sections of the upstream class tests (core `bc9cc12`) are
//! reproduced with their literals. Each header is also compared in full with
//! the text transcribed from the source's stream statements
//! (`FLASHDeconvSpectrumFile.cpp:236-324`, `FLASHDeconvFeatureFile.cpp:18-44`),
//! tier 3.

use openms::format::flash_deconv_feature_file::FLASHDeconvFeatureFile;
use openms::format::flash_deconv_spectrum_file::FLASHDeconvSpectrumFile;
use openms::param::{Param, ParamValue};

fn spectrum_header(ms_level: u32, detail: bool, decoy: bool) -> String {
    let mut out = Vec::new();
    FLASHDeconvSpectrumFile::write_deconvolved_masses_header(&mut out, ms_level, detail, decoy)
        .unwrap();
    String::from_utf8(out).unwrap()
}

// FLASHDeconvSpectrumFile_test.cpp: START_SECTION(static void
// writeDeconvolvedMassesHeader(std::ostream& os, uint ms_level, bool detail,
// bool report_decoy))
#[test]
fn deconvolved_masses_header_matches_the_class_test() {
    let header = spectrum_header(1, false, false);
    for column in [
        "Index",
        "FileName",
        "ScanNum",
        "RetentionTime",
        "MonoisotopicMass",
        "AverageMass",
        "SumIntensity",
        "IsotopeCosine",
        "Qscore",
    ] {
        assert!(header.contains(column), "{column}");
    }
    assert!(!header.contains("PrecursorScanNum"));
    assert!(!header.contains("PeakMZs"));
    assert!(!header.contains("TargetDecoyType"));

    let header = spectrum_header(1, true, false);
    for column in [
        "PeakMZs",
        "PeakIntensities",
        "PeakCharges",
        "PerChargeIntensity",
        "PerIsotopeIntensity",
    ] {
        assert!(header.contains(column), "{column}");
    }

    let header = spectrum_header(1, false, true);
    assert!(header.contains("TargetDecoyType"));
    assert!(header.contains("Qvalue"));

    let header = spectrum_header(2, false, false);
    for column in [
        "PrecursorScanNum",
        "PrecursorMz",
        "PrecursorCharge",
        "PrecursorMonoisotopicMass",
    ] {
        assert!(header.contains(column), "{column}");
    }

    let header = spectrum_header(2, true, true);
    for column in [
        "PrecursorScanNum",
        "PeakMZs",
        "TargetDecoyType",
        "PrecursorQvalue",
    ] {
        assert!(header.contains(column), "{column}");
    }
}

// The eight headers in full, transcribed from FLASHDeconvSpectrumFile.cpp:236-324.
#[test]
fn deconvolved_masses_headers_match_the_source_text() {
    let lead = "Index\tFileName\tScanNum\tFeatureIndex\t";
    let mass = "RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\tSumIntensity\tMinCharge\tMaxCharge\t";
    let detail = "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\tNoisePeakMZs\tNoisePeakIntensities\tNoisePeakCharges\tNoisePeakMasses\tNoisePeakIsotopeIndices\tNoisePeakPPMErrors\t";
    let precursor = "PrecursorScanNum\tPrecursorFeatureIndex\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorSNR\tPrecursorMonoisotopicMass\tPrecursorQscore\tPrecursorQscore2D\t";
    let scores = "IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQscore\tQscore2D\t";
    let tail = "PerChargeIntensity\tPerIsotopeIntensity\n";

    // detail, MS1 (lines 240-257)
    assert_eq!(
        spectrum_header(1, true, false),
        format!("{lead}{mass}{detail}{scores}{tail}")
    );
    assert_eq!(
        spectrum_header(1, true, true),
        format!("{lead}TargetDecoyType\t{mass}{detail}{scores}Qvalue\t{tail}")
    );
    // detail, MS2 (lines 258-280)
    assert_eq!(
        spectrum_header(2, true, false),
        format!("{lead}{mass}{detail}{precursor}{scores}{tail}")
    );
    assert_eq!(
        spectrum_header(3, true, true),
        format!(
            "{lead}TargetDecoyType\t{mass}{detail}{precursor}PrecursorQvalue\t{scores}Qvalue\t{tail}"
        )
    );
    // summary, MS1 (lines 284-300): no tab after Qvalue.
    assert_eq!(
        spectrum_header(1, false, false),
        format!("{lead}{mass}PeakCount\t{scores}\n")
    );
    assert_eq!(
        spectrum_header(1, false, true),
        format!("{lead}TargetDecoyType\t{mass}PeakCount\t{scores}Qvalue\n")
    );
    // summary, MS2 (lines 301-322)
    assert_eq!(
        spectrum_header(2, false, false),
        format!("{lead}{mass}PeakCount\t{precursor}{scores}\n")
    );
    assert_eq!(
        spectrum_header(0, false, true),
        format!(
            "{lead}TargetDecoyType\t{mass}PeakCount\t{precursor}PrecursorQvalue\t{scores}Qvalue\t\n"
        )
    );
}

// START_SECTION(static void writeTopFDHeader(std::ostream& os, const Param& param))
#[test]
fn top_fd_header_matches_the_class_test() {
    let mut params = Param::new();
    params
        .set_value("test_param", ParamValue::Integer(42), "", &[])
        .unwrap();
    params
        .set_value("another_param", ParamValue::String("value".into()), "", &[])
        .unwrap();
    let mut out = Vec::new();
    FLASHDeconvSpectrumFile::write_top_fd_header(&mut out, &params).unwrap();
    let header = String::from_utf8(out).unwrap();
    for text in [
        "#FLASHDeconv generated msalign file",
        "Parameters",
        "test_param",
        "42",
        "another_param",
        "value",
    ] {
        assert!(header.contains(text), "{text}");
    }
    let rule = "####################### Parameters ######################\n";
    assert_eq!(
        header,
        format!(
            "#FLASHDeconv generated msalign file\n{rule}#test_param: 42\n#another_param: value\n{rule}"
        )
    );

    // Entry names are written without their section path.
    let mut nested = Param::new();
    nested
        .set_value("algorithm:tol", ParamValue::Float(10.0), "", &[])
        .unwrap();
    let mut out = Vec::new();
    FLASHDeconvSpectrumFile::write_top_fd_header(&mut out, &nested).unwrap();
    assert!(String::from_utf8(out).unwrap().contains("\n#tol: 10\n"));
}

fn feature_header(decoy: bool) -> String {
    let mut out = Vec::new();
    FLASHDeconvFeatureFile::write_header(&mut out, decoy).unwrap();
    String::from_utf8(out).unwrap()
}

fn top_fd_feature_header(ms_level: u32) -> String {
    let mut out = Vec::new();
    FLASHDeconvFeatureFile::write_top_fd_feature_header(&mut out, ms_level).unwrap();
    String::from_utf8(out).unwrap()
}

// FLASHDeconvFeatureFile_test.cpp: START_SECTION(static void
// writeHeader(std::ostream& os, bool report_decoy = false))
#[test]
fn feature_header_matches_the_class_test() {
    let header = feature_header(false);
    for column in [
        "FeatureIndex",
        "FileName",
        "MSLevel",
        "MonoisotopicMass",
        "AverageMass",
        "SumIntensity",
        "PerChargeIntensity",
        "PerIsotopeIntensity",
        "\n",
    ] {
        assert!(header.contains(column), "{column}");
    }
    assert!(!header.contains("IsDecoy"));
    let header = feature_header(true);
    assert!(header.contains("IsDecoy"));
    assert!(header.contains("FeatureIndex"));

    // FLASHDeconvFeatureFile.cpp:18-27 in full.
    let rest = "\tMonoisotopicMass\tAverageMass\tMassCount\tStartRetentionTime\tEndRetentionTime\tRetentionTimeDuration\tApexRetentionTime\tSumIntensity\tMaxIntensity\tFeatureQuantity\tMinCharge\tMaxCharge\tChargeCount\tIsotopeCosineScore\tQscore2D\tPerChargeIntensity\tPerIsotopeIntensity\n";
    assert_eq!(
        feature_header(false),
        format!("FeatureIndex\tFileName\tMSLevel{rest}")
    );
    assert_eq!(
        feature_header(true),
        format!("FeatureIndex\tFileName\tMSLevel\tIsDecoy{rest}")
    );
}

// START_SECTION(static void writeTopFDFeatureHeader(std::ostream& os, uint ms_level))
#[test]
fn top_fd_feature_header_matches_the_class_test() {
    let header = top_fd_feature_header(1);
    for column in [
        "File_name",
        "Feature_ID",
        "Mass",
        "Intensity",
        "Min_time",
        "Max_time",
        "Apex_time",
        "EC_score",
    ] {
        assert!(header.contains(column), "{column}");
    }
    assert!(!header.contains("Spectrum_ID"));
    let header = top_fd_feature_header(2);
    for column in [
        "File_name",
        "Spectrum_ID",
        "Precursor_charge",
        "Precursor_intensity",
        "Fraction_feature_ID",
    ] {
        assert!(header.contains(column), "{column}");
    }

    // FLASHDeconvFeatureFile.cpp:35-42 in full (the MS1 and other-level branches).
    assert_eq!(
        top_fd_feature_header(1),
        "File_name\tFraction_ID\tFeature_ID\tMass\tIntensity\tMin_time\tMax_time\tMin_scan\tMax_scan\tMin_charge\tMax_charge\tApex_time\tApex_scan\tApex_intensity\tRep_charge\tRep_average_mz\tEnvelope_num\tEC_score\n"
    );
    assert_eq!(
        top_fd_feature_header(2),
        "File_name\tFraction_ID\tSpectrum_ID\tScans\tMS_one_ID\tMS_one_scans\tFraction_feature_ID\tFraction_feature_intensity\tFraction_feature_score\tFraction_feature_min_time\tFraction_feature_max_time\tFraction_feature_apex_time\tPrecursor_monoisotopic_mz\tPrecursor_average_mz\tPrecursor_charge\tPrecursor_intensity\n"
    );
    assert_eq!(top_fd_feature_header(0), top_fd_feature_header(2));
}

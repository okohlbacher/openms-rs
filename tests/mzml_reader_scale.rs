#![cfg(feature = "mzml")]
// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Size-derived mzML reader ceilings, and source-compatible reading of the
//! timestamp sentinels ProteoWizard and Boost write.
//!
//! Two blockers of the OpenMS4 benchmark are covered here.
//!
//! 1. `run/@startTimeStamp="-infinity"`, which ProteoWizard writes for a vendor
//!    file without an acquisition date, made every TOPP tool exit 6 on the
//!    PXD001819 `50amol_R1.mzML` benchmark input.
//!    `data/mzml_reader_scale/cpp_timestamp_oracle.tsv` records what the C++
//!    Release build (core `bc9cc12`, install prefix
//!    `openms4-release-bc9cc12-c19e494-174b576`) did with that attribute and
//!    with `infinity`, `not-a-date-time`, an empty value and a valid control on
//!    `data/mzml_reader_scale/pxd001819_50amol_r1_first3.mzML`: `FileInfo` and
//!    `FileConverter` both exit 0, log `DateTime conversion error of "<text>"`
//!    as a non-fatal error, and write the file back without `startTimeStamp`.
//!    The source is lenient unconditionally, so the reader is too by default;
//!    a caller can opt out with `source_invalid_timestamps: false`.
//!    The driver is `../oracle/mzml-reader-scale/datetime_sentinel_cpp.sh`;
//!    hashes are in `data/mzml_reader_scale_provenance.json`.
//!
//! 2. The reader's fixed cumulative ceilings rejected every benchmark input
//!    from 0.5 to 2.3 GB. They are now size-derived: `floor + rate *
//!    consumed`, where counts that need their own start tag (records, binary
//!    arrays, parameter groups) grow once per group of consumed bytes. The
//!    synthetic documents below are large enough that the former fixed floors
//!    reject them, which `InputScaling::fixed` reproduces, while amplifying
//!    documents of the same size stay rejected, and the per-array ceiling stays
//!    absolute. The `#[ignore]`d tests at the end read the real benchmark
//!    inputs on the HPC nodes.

use base64::{Engine, engine::general_purpose::STANDARD};
use flate2::{Compression, write::ZlibEncoder};
use openms::{
    Error, MSExperiment,
    concept::log_stream::{LogLevel, LogSink, with_thread_local_log},
    data_structures::DateTime,
    format::{
        FileHandler, FileType, PeakFileOptions,
        mzml::{self, Allowance, InputScaling, LoadOptions, ReadOptions},
    },
    kernel::NumericRange,
};
use std::{
    io::{self, BufReader, Cursor, Write},
    path::PathBuf,
    sync::{Arc, Mutex},
};

/// The three-spectrum slice of the PXD001819 benchmark input the C++ oracle ran
/// on; `startTimeStamp="-infinity"` is ProteoWizard's own value.
const VELOS_FIRST3: &str = include_str!("data/mzml_reader_scale/pxd001819_50amol_r1_first3.mzML");
/// The executed C++ Release results for each timestamp case.
const CPP_ORACLE: &str = include_str!("data/mzml_reader_scale/cpp_timestamp_oracle.tsv");

/// The `startTimeStamp` spelling of one oracle case.
fn case(name: &str) -> (String, String) {
    let row = CPP_ORACLE
        .lines()
        .skip(1)
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .find(|row| row[0] == name)
        .unwrap_or_else(|| panic!("no oracle case {name}"));
    let stamp = if row[1] == "(empty)" { "" } else { row[1] };
    let mut xml = VELOS_FIRST3.replace(
        "startTimeStamp=\"-infinity\"",
        &format!("startTimeStamp=\"{stamp}\""),
    );
    if row[2] != "-" {
        let inserted = format!(
            "<cvParam cvRef=\"MS\" accession=\"MS:1000747\" name=\"completion time\" value=\"{}\"/>",
            row[2]
        );
        let anchor = "<processingMethod order=\"1\" softwareRef=\"pwiz\">";
        assert!(xml.contains(anchor));
        xml = xml.replacen(anchor, &format!("{anchor}{inserted}"), 1);
    }
    (xml, row[3].to_owned())
}

fn read(xml: &str, options: &ReadOptions) -> openms::Result<MSExperiment> {
    mzml::read_with_options(Cursor::new(xml.as_bytes()), options)
}

/// The `startTimeStamp` attribute of a written run element, or `None`.
fn written_start_time(experiment: &MSExperiment) -> Option<String> {
    let mut out = Vec::new();
    mzml::write(&mut out, experiment).unwrap();
    let text = String::from_utf8(out).unwrap();
    let run = text
        .lines()
        .find(|line| line.trim_start().starts_with("<run "))
        .expect("written run element");
    run.split("startTimeStamp=\"")
        .nth(1)
        .map(|rest| rest.split('"').next().unwrap().to_owned())
}

/// The options a caller passes to keep an unparseable timestamp an error, which
/// the executed C++ never does.
fn strict_timestamps() -> ReadOptions {
    ReadOptions {
        source_invalid_timestamps: false,
        ..ReadOptions::default()
    }
}

#[test]
fn timestamp_sentinels_follow_the_executed_cpp_result_by_default() {
    discard_warnings();
    let source = ReadOptions::source();
    assert!(source.source_invalid_timestamps && source.source_dangling_references);
    // The executed C++ is unconditionally lenient here, so the default is too:
    // `FileHandler::load_experiment` and every TOPP tool behind it read the
    // PXD001819 input whose `startTimeStamp` is `-infinity`.
    assert!(ReadOptions::default().source_invalid_timestamps);
    assert!(!ReadOptions::default().source_dangling_references);
    for name in [
        "minus_infinity",
        "infinity",
        "not_a_date_time",
        "empty",
        "valid_control",
        "completion_not_a_date_time",
    ] {
        let (xml, written) = case(name);
        let experiment = read(&xml, &source).unwrap_or_else(|e| panic!("{name}: {e}"));
        // The default reads every case exactly as the source option does.
        assert_eq!(
            read(&xml, &ReadOptions::default()).unwrap_or_else(|e| panic!("{name}: {e}")),
            experiment,
            "{name}"
        );
        // C++ writes the run element back with, or without, startTimeStamp.
        let expected = (written != "absent").then(|| written.clone());
        assert_eq!(written_start_time(&experiment), expected, "{name}");
        assert_eq!(
            experiment.settings.date_time.is_valid(),
            expected.is_some(),
            "{name}"
        );
        // Source `asDateTime_` keeps no raw text for a rejected timestamp, and
        // the C++ writer omits `mzml_start_time_stamp` from the run userParams.
        assert!(
            !experiment
                .settings
                .metadata
                .contains_key("mzml_start_time_stamp"),
            "{name}"
        );
        // Every spectrum carries the file's three pwiz processing methods; the
        // completion time the C++ run dropped is unset in the second one.
        let processing = &experiment.spectra[0].data_processing;
        assert_eq!(processing.len(), 3, "{name}");
        assert!(
            processing.iter().all(|step| step.completion_time.is_none()),
            "{name}"
        );
    }
    // The dropped completion time is the only difference between those two
    // documents, as it is for the C++ outputs (one MS:1000747 each, the
    // converter's own).
    assert_eq!(
        read(&case("completion_not_a_date_time").0, &source).unwrap(),
        read(&case("minus_infinity").0, &source).unwrap()
    );
    assert_eq!(
        read(&case("valid_control").0, &source)
            .unwrap()
            .settings
            .date_time,
        DateTime::parse("2014-01-01T00:00:00").unwrap()
    );
}

#[test]
fn a_caller_that_opts_out_still_refuses_an_unparseable_timestamp() {
    for name in [
        "minus_infinity",
        "infinity",
        "not_a_date_time",
        "empty",
        "completion_not_a_date_time",
    ] {
        let (xml, _) = case(name);
        let error = read(&xml, &strict_timestamps()).unwrap_err();
        assert!(
            matches!(&error, Error::InvalidValue(text)
                if text.contains("invalid DateTime input or calendar fields")),
            "{name}: {error}"
        );
    }
    // The valid control reads under every policy, with the same result.
    let (xml, _) = case("valid_control");
    assert_eq!(
        read(&xml, &strict_timestamps()).unwrap(),
        read(&xml, &ReadOptions::source()).unwrap()
    );
    assert_eq!(
        read(&xml, &ReadOptions::default()).unwrap(),
        read(&xml, &ReadOptions::source()).unwrap()
    );
}

#[test]
fn each_dropped_timestamp_warns_once_and_valid_or_empty_ones_stay_silent() {
    let capture = Capture::default();
    let sink = LogSink::new(capture.clone());
    with_thread_local_log(LogLevel::Warn, |log| {
        log.remove_all_streams()?;
        log.set_color(None);
        log.insert(&sink)
    })
    .unwrap();
    let take = || {
        with_thread_local_log(LogLevel::Warn, |log| log.clear_cache()).unwrap();
        String::from_utf8(std::mem::take(&mut *capture.0.lock().unwrap())).unwrap()
    };
    read(&case("minus_infinity").0, &ReadOptions::source()).unwrap();
    assert_eq!(
        take(),
        "Warning: mzML run startTimeStamp '-infinity' is not a date-time; \
         source-compatible reading leaves it unset.\n"
    );
    read(
        &case("completion_not_a_date_time").0,
        &ReadOptions::source(),
    )
    .unwrap();
    let warnings = take();
    let lines: Vec<_> = warnings.lines().collect();
    assert_eq!(lines.len(), 2, "{lines:?}");
    assert!(
        lines
            .iter()
            .any(|line| line.contains("processingMethod completion time 'not-a-date-time'")),
        "{lines:?}"
    );
    assert!(
        lines
            .iter()
            .any(|line| line.contains("run startTimeStamp '-infinity'")),
        "{lines:?}"
    );
    // Neither an empty nor a valid timestamp is a dropped value.
    for name in ["empty", "valid_control"] {
        read(&case(name).0, &ReadOptions::source()).unwrap();
        assert_eq!(take(), "", "{name}");
    }
    // A read that opted out fails instead of warning.
    assert!(read(&case("minus_infinity").0, &strict_timestamps()).is_err());
    assert_eq!(take(), "");
}

/// Send this test thread's warnings to a discarding sink.
fn discard_warnings() {
    with_thread_local_log(LogLevel::Warn, |log| {
        log.remove_all_streams()?;
        log.insert(&LogSink::new(io::sink()))
    })
    .unwrap();
}

#[derive(Clone, Default)]
struct Capture(Arc<Mutex<Vec<u8>>>);
impl Write for Capture {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// How a synthetic document encodes its binary arrays.
#[derive(Clone, Copy, PartialEq)]
enum Encoding {
    /// Uncompressed 64-bit m/z and 32-bit intensity arrays.
    Plain,
    /// The same arrays, zlib compressed.
    Zlib,
    /// Numpress linear m/z and short-logged-float intensity arrays.
    Numpress,
}

/// Base64 of `bytes`, zlib compressed when `compress`.
fn encode(bytes: &[u8], compress: bool) -> String {
    if !compress {
        return STANDARD.encode(bytes);
    }
    let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(bytes).unwrap();
    STANDARD.encode(encoder.finish().unwrap())
}

/// One binary array element of a synthetic spectrum.
fn array(kind: &str, terms: &str, text: &str) -> String {
    format!(
        "<binaryDataArray encodedLength=\"{}\">{terms}\
         <cvParam cvRef=\"MS\" accession=\"{kind}\" name=\"array\"/><binary>{text}</binary>\
         </binaryDataArray>",
        text.len()
    )
}

/// A synthetic mzML document: `spectra` spectra of `peaks` points each, with
/// the parameter payload a converter writes.
///
/// Returns the document and the number of binary bytes behind its base64
/// payload, which the Numpress work allowance is measured against.
fn synthetic(spectra: usize, peaks: usize, encoding: Encoding) -> (String, usize) {
    let mut xml = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\
         <mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">\
         <cvList count=\"1\"><cv id=\"MS\" fullName=\"PSI-MS\" URI=\"https://purl.obolibrary.org/obo/ms.obo\"/></cvList>\
         <fileDescription><fileContent/></fileDescription>\
         <referenceableParamGroupList count=\"1\"><referenceableParamGroup id=\"common\">\
         <cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\"/>\
         <cvParam cvRef=\"MS\" accession=\"MS:1000130\" name=\"positive scan\"/>\
         </referenceableParamGroup></referenceableParamGroupList>\
         <softwareList count=\"1\"><software id=\"sw\" version=\"1.0\"/></softwareList>\
         <instrumentConfigurationList count=\"1\"><instrumentConfiguration id=\"ic\"/></instrumentConfigurationList>\
         <dataProcessingList count=\"1\"><dataProcessing id=\"dp\">\
         <processingMethod order=\"0\" softwareRef=\"sw\">\
         <cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\"/>\
         </processingMethod></dataProcessing></dataProcessingList>\
         <run id=\"run\" defaultInstrumentConfigurationRef=\"ic\" startTimeStamp=\"2016-11-18T23:31:16\">",
    );
    xml += &format!("<spectrumList count=\"{spectra}\" defaultDataProcessingRef=\"dp\">");
    let mz: Vec<f64> = (0..peaks).map(|i| 300.0 + i as f64 * 0.5).collect();
    let intensity: Vec<f64> = (0..peaks).map(|i| 100.0 + (i % 17) as f64).collect();
    let (mz_terms, mz_text, intensity_terms, intensity_text) = match encoding {
        Encoding::Plain | Encoding::Zlib => {
            let zlib = encoding == Encoding::Zlib;
            let compression = if zlib {
                "<cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\"/>"
            } else {
                "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>"
            };
            let mz_bytes: Vec<u8> = mz.iter().flat_map(|v| v.to_le_bytes()).collect();
            let intensity_bytes: Vec<u8> = intensity
                .iter()
                .flat_map(|v| (*v as f32).to_le_bytes())
                .collect();
            (
                format!(
                    "<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>{compression}"
                ),
                encode(&mz_bytes, zlib),
                format!(
                    "<cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\"/>{compression}"
                ),
                encode(&intensity_bytes, zlib),
            )
        }
        Encoding::Numpress => {
            let fixed = openms::format::numpress::optimal_linear_fixed_point(&mz).unwrap();
            let mz_bytes = openms::format::numpress::encode_linear(&mz, fixed).unwrap();
            let slof = openms::format::numpress::optimal_slof_fixed_point(&intensity).unwrap();
            let intensity_bytes = openms::format::numpress::encode_slof(&intensity, slof).unwrap();
            (
                "<cvParam cvRef=\"MS\" accession=\"MS:1002312\" name=\"MS-Numpress linear\"/>"
                    .to_owned(),
                encode(&mz_bytes, false),
                "<cvParam cvRef=\"MS\" accession=\"MS:1002314\" name=\"MS-Numpress slof\"/>"
                    .to_owned(),
                encode(&intensity_bytes, false),
            )
        }
    };
    let payload = (mz_text.len() + intensity_text.len()) / 4 * 3;
    for index in 0..spectra {
        xml += &format!(
            "<spectrum id=\"controllerType=0 controllerNumber=1 scan={}\" index=\"{index}\" defaultArrayLength=\"{peaks}\">\
             <referenceableParamGroupRef ref=\"common\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"1\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000504\" name=\"base peak m/z\" value=\"400.5\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000505\" name=\"base peak intensity\" value=\"1000.0\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000285\" name=\"total ion current\" value=\"12345.0\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000528\" name=\"lowest observed m/z\" value=\"300.0\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000527\" name=\"highest observed m/z\" value=\"1500.0\"/>\
             <userParam name=\"filter string\" value=\"FTMS + p NSI Full ms [300.00-1500.00]\"/>\
             <userParam name=\"preset scan configuration\" value=\"1\"/>\
             <scanList count=\"1\"><cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\"/>\
             <scan><cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"{}\" unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"/>\
             <scanWindowList count=\"1\"><scanWindow>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000501\" name=\"scan window lower limit\" value=\"300.0\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000500\" name=\"scan window upper limit\" value=\"1500.0\"/>\
             </scanWindow></scanWindowList></scan></scanList>\
             <binaryDataArrayList count=\"2\">{}{}</binaryDataArrayList></spectrum>",
            index + 1,
            index as f64 * 0.5,
            array("MS:1000514", &mz_terms, &mz_text),
            array("MS:1000515", &intensity_terms, &intensity_text),
        );
    }
    xml += "</spectrumList></run></mzML>";
    (xml, payload * spectra)
}

/// The scientific options `FileHandler::load_experiment_with_options` passes for
/// FeatureFinderCentroided: MS1 only, positive intensities, sorted by m/z.
fn feature_finder_load() -> LoadOptions {
    let mut scientific = PeakFileOptions::default();
    scientific.set_ms_levels(&[1]).unwrap();
    scientific.set_intensity_range(NumericRange {
        min: f64::MIN_POSITIVE,
        max: f64::MAX,
    });
    LoadOptions {
        scientific,
        ..Default::default()
    }
}

#[test]
fn a_realistic_document_beyond_every_former_fixed_ceiling_reads_with_the_defaults() {
    let (xml, _) = synthetic(20_000, 3, Encoding::Plain);
    // A few MB of XML: the size a converter writes for 20,000 sparse spectra.
    assert!((10..40).contains(&(xml.len() / (1 << 20))), "{}", xml.len());
    let experiment = read(&xml, &ReadOptions::default()).unwrap();
    assert_eq!(experiment.spectra.len(), 20_000);
    assert_eq!(experiment.spectra[19_999].peaks.len(), 3);
    assert_eq!(experiment.spectra[19_999].peaks[0].mz, 300.0);
    // The former fixed ceilings reject the same document; only the parameter
    // storage estimate is reached, at about 7 KB of charge per spectrum.
    let fixed = ReadOptions {
        scaling: InputScaling::default().fixed(),
        ..Default::default()
    };
    let error = read(&xml, &fixed).unwrap_err();
    assert!(
        error.to_string().contains("parameter bytes exceed"),
        "{error}"
    );
    // Both the counter and the FeatureFinderCentroided load path agree.
    assert_eq!(
        mzml::read_size_with_options(
            Cursor::new(xml.as_bytes()),
            &PeakFileOptions::default(),
            &ReadOptions::default()
        )
        .unwrap()
        .spectra,
        20_000
    );
    assert!(
        mzml::read_size_with_options(
            Cursor::new(xml.as_bytes()),
            &PeakFileOptions::default(),
            &fixed
        )
        .is_err()
    );
    let selected = mzml::read_with_load_options(
        Cursor::new(xml.as_bytes()),
        &feature_finder_load(),
        &ReadOptions::default(),
    )
    .unwrap();
    assert_eq!(selected.spectra.len(), 20_000);
}

#[test]
fn compressed_and_numpress_documents_read_at_the_same_scale() {
    for (encoding, spectra, peaks) in [
        (Encoding::Zlib, 4_000, 200),
        (Encoding::Numpress, 6_500, 500),
    ] {
        let (xml, payload) = synthetic(spectra, peaks, encoding);
        let experiment = read(&xml, &ReadOptions::default()).unwrap();
        assert_eq!(experiment.spectra.len(), spectra);
        assert_eq!(experiment.spectra[spectra - 1].peaks.len(), peaks);
        assert!((experiment.spectra[0].peaks[7].mz - 303.5).abs() < 1e-4);
        if encoding == Encoding::Numpress {
            // `numpress_coder::decode_text` charges 64 work units per encoded
            // byte, so this document costs more than the coder's whole former
            // per-document allowance of 500,000,000 work units; the allowance
            // is now per array.
            assert!(payload * 64 > 500_000_000, "{payload}");
        }
    }
}

#[test]
fn a_small_document_that_declares_or_amplifies_huge_quantities_is_still_refused() {
    let (prefix, _) = synthetic(1, 3, Encoding::Plain);
    // A record declaring a billion points, with a three-point payload.
    let huge = prefix.replace(
        "defaultArrayLength=\"3\"",
        "defaultArrayLength=\"1000000000\"",
    );
    let error = read(&huge, &ReadOptions::default()).unwrap_err();
    assert!(error.to_string().contains("peak count exceeds"), "{error}");
    // A zlib payload that inflates 4,000-fold, declared as 50 million points.
    let zeros = vec![0u8; 400 * (1 << 20)];
    let mut encoder = ZlibEncoder::new(Vec::new(), Compression::best());
    encoder.write_all(&zeros).unwrap();
    let bomb = STANDARD.encode(encoder.finish().unwrap());
    assert!(bomb.len() < (1 << 20), "{}", bomb.len());
    let (mut xml, _) = synthetic(1, 3, Encoding::Zlib);
    let start = xml.find("<binaryDataArray ").unwrap();
    let end = xml.find("</binaryDataArrayList>").unwrap();
    xml.replace_range(
        start..end,
        &format!(
            "{}{}",
            array(
                "MS:1000514",
                "<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>\
                 <cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\"/>",
                &bomb
            ),
            array(
                "MS:1000515",
                "<cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\"/>\
                 <cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\"/>",
                &bomb
            ),
        ),
    );
    let xml = xml.replace(
        "defaultArrayLength=\"3\"",
        "defaultArrayLength=\"50000000\"",
    );
    let error = read(&xml, &ReadOptions::default()).unwrap_err();
    assert!(error.to_string().contains("peak count exceeds"), "{error}");
    // The same bomb below the peak allowance is stopped by the decoded-byte
    // allowance derived from the document's own size.
    let xml = xml.replace(
        "defaultArrayLength=\"50000000\"",
        "defaultArrayLength=\"9000000\"",
    );
    let error = read(&xml, &ReadOptions::default()).unwrap_err();
    assert!(
        error.to_string().contains("byte limit")
            || error.to_string().contains("element limit exceeded"),
        "{error}"
    );
}

#[test]
fn parameter_group_amplification_is_bounded_by_the_consumed_bytes() {
    // A parameter-group reference is 40 bytes of XML and is charged for every
    // parameter it expands to. A group of 20,000 parameters referenced by every
    // spectrum therefore charges about 17,000 storage bytes per input byte,
    // 66 times the 256 the default allowance grants, so the read fails after
    // about a megabyte however large the document is.
    let mut group = String::new();
    for index in 0..20_000 {
        group += &format!("<userParam name=\"p{index}\" value=\"v\"/>");
    }
    let (prefix, _) = synthetic(1_000, 3, Encoding::Plain);
    let xml = prefix.replace(
        "<cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\"/>\
         <cvParam cvRef=\"MS\" accession=\"MS:1000130\" name=\"positive scan\"/>",
        &group,
    );
    assert!(xml.len() > 2 * (1 << 20), "{}", xml.len());
    let error = read(&xml, &ReadOptions::default()).unwrap_err();
    assert!(
        error.to_string().contains("parameter bytes exceed")
            || error.to_string().contains("parameter count exceeds"),
        "{error}"
    );
    // The same 1,000 spectra without the amplifying group read normally.
    assert_eq!(
        read(&prefix, &ReadOptions::default())
            .unwrap()
            .spectra
            .len(),
        1_000
    );
    // Deep nesting is refused whatever the allowances say.
    let deep = prefix.replace(
        "<scanWindowList count=\"1\">",
        &"<scanWindowList count=\"1\">".repeat(200),
    );
    assert!(read(&deep, &ReadOptions::default()).is_err());
}

#[test]
fn explicit_absolute_ceilings_still_win_over_the_size_derived_allowances() {
    let (xml, _) = synthetic(200, 3, Encoding::Plain);
    assert!(read(&xml, &ReadOptions::default()).is_ok());
    for (options, expected) in [
        (
            ReadOptions {
                max_param_bytes: 10_000,
                ..Default::default()
            },
            "parameter bytes exceed",
        ),
        (
            ReadOptions {
                max_total_peaks: 10,
                ..Default::default()
            },
            "peak count exceeds",
        ),
        (
            ReadOptions {
                max_records: 10,
                ..Default::default()
            },
            "record count exceeds",
        ),
        (
            ReadOptions {
                max_total_array_elements: 10,
                ..Default::default()
            },
            "element limit exceeded",
        ),
    ] {
        let error = read(&xml, &options)
            .err()
            .unwrap_or_else(|| panic!("expected {expected}"));
        assert!(error.to_string().contains(expected), "{error}");
    }
    // A short XML ceiling truncates the stream, which fails as malformed XML
    // or, for a document that ends inside the allowance, with the ceiling's own
    // message; both are refusals of the same explicit ceiling.
    let error = read(
        &xml,
        &ReadOptions {
            max_xml_bytes: xml.len() as u64 - 1,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("XML exceeds configured byte limit"),
        "{error}"
    );
    // A caller can also tighten the size-derived side alone.
    let strict = ReadOptions {
        scaling: InputScaling {
            array_bytes: Allowance::new(0, 0),
            ..InputScaling::default()
        },
        ..Default::default()
    };
    let error = read(&xml, &strict).unwrap_err();
    assert!(error.to_string().contains("byte limit"), "{error}");
    assert_eq!(Allowance::new(7, 3).after(5), 22);
}

#[test]
fn the_per_array_ceiling_does_not_grow_with_the_document() {
    // 64 MiB per array, 8 million f64 values. The largest single binary array
    // in any benchmark input is the `UK222.mzML` TIC chromatogram, 53,824
    // elements and 431 KB decoded (measured on ibminode06, see
    // `docs/MZML_READER_SCALE_SUPPORT.md`), so the ceiling is 155x the largest
    // real array and reading one costs no more than the array itself.
    assert_eq!(ReadOptions::default().max_array_bytes, 64 * (1 << 20));
    let (real, _) = synthetic(1, 53_824, Encoding::Plain);
    assert_eq!(
        read(&real, &ReadOptions::default()).unwrap().spectra[0]
            .peaks
            .len(),
        53_824
    );
    // A small document whose two zlib arrays each inflate to 80 MB stays below
    // the cumulative peak, element and byte allowances a document of any size
    // is granted, and is refused by the per-array ceiling alone.
    let zeros = vec![0u8; 80 * (1 << 20)];
    let mut encoder = ZlibEncoder::new(Vec::new(), Compression::best());
    encoder.write_all(&zeros).unwrap();
    let bomb = STANDARD.encode(encoder.finish().unwrap());
    let (mut xml, _) = synthetic(1, 3, Encoding::Zlib);
    let start = xml.find("<binaryDataArray ").unwrap();
    let end = xml.find("</binaryDataArrayList>").unwrap();
    let terms = |kind: &str| {
        format!(
            "<cvParam cvRef=\"MS\" accession=\"{kind}\" name=\"float\"/>\
             <cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\"/>"
        )
    };
    xml.replace_range(
        start..end,
        &format!(
            "{}{}",
            array("MS:1000514", &terms("MS:1000523"), &bomb),
            array("MS:1000515", &terms("MS:1000523"), &bomb),
        ),
    );
    let xml = xml.replace("defaultArrayLength=\"3\"", "defaultArrayLength=\"9999999\"");
    assert!(xml.len() < 1 << 20, "{}", xml.len());
    let scaling = InputScaling::default();
    assert!(scaling.peaks.after(0) >= 9_999_999);
    assert!(scaling.array_elements.after(0) >= 2 * 9_999_999);
    assert!(scaling.array_bytes.after(0) >= 2 * 8 * 9_999_999);
    let error = read(&xml, &ReadOptions::default()).unwrap_err();
    assert!(error.to_string().contains("byte limit"), "{error}");
}

#[test]
fn record_array_and_group_counts_are_bounded_by_the_consumed_bytes() {
    let scaling = InputScaling::default();
    // The floors are the former fixed ceilings, which already hold every real
    // file: the largest benchmark input declares 40,857 records, 81,714 arrays
    // and one parameter group.
    assert_eq!(scaling.records.after(0), 1_000_000);
    assert_eq!(scaling.arrays.after(0), 1_000_000);
    assert_eq!(scaling.param_groups.after(0), 100_000);
    // One record per 512 consumed bytes, one array per 256, one group per 4 KiB.
    assert_eq!(scaling.records.after(3 * 512), 1_000_003);
    assert_eq!(scaling.arrays.after(3 * 256), 1_000_003);
    assert_eq!(scaling.param_groups.after(3 * 4_096), 100_003);
    // The densest benchmark input writes one record per 3,951 bytes and one
    // array per 3,559, so the rates are 7.7x and 13.9x the measured density,
    // while 2,000,000 empty records in 473 MiB of XML (236 bytes each) is
    // refused.
    assert!(scaling.records.after(473 * (1 << 20)) < 2_000_000);
    assert!(scaling.records.after(2_317_975_830) > 80 * 40_857);
    assert!(scaling.arrays.after(2_317_975_830) > 80 * 81_714);
    // A tightened allowance refuses a document the defaults read, at the count
    // it tightens.
    let (xml, _) = synthetic(2_000, 3, Encoding::Plain);
    assert_eq!(
        read(&xml, &ReadOptions::default()).unwrap().spectra.len(),
        2_000
    );
    for (tightened, expected) in [
        (
            InputScaling {
                records: Allowance::every(10, 4_096),
                ..InputScaling::default()
            },
            "record count exceeds",
        ),
        (
            InputScaling {
                arrays: Allowance::every(10, 4_096),
                ..InputScaling::default()
            },
            "binary array count limit exceeded",
        ),
        (
            InputScaling {
                param_groups: Allowance::fixed(0),
                ..InputScaling::default()
            },
            "parameter group count",
        ),
    ] {
        let error = read(
            &xml,
            &ReadOptions {
                scaling: tightened,
                ..Default::default()
            },
        )
        .unwrap_err();
        assert!(error.to_string().contains(expected), "{expected}: {error}");
    }
    // The counter agrees with the reader on the same document.
    assert_eq!(
        mzml::read_size_with_options(
            Cursor::new(xml.as_bytes()),
            &PeakFileOptions::default(),
            &ReadOptions::default()
        )
        .unwrap()
        .spectra,
        2_000
    );
    assert!(
        mzml::read_size_with_options(
            Cursor::new(xml.as_bytes()),
            &PeakFileOptions::default(),
            &ReadOptions {
                scaling: InputScaling {
                    records: Allowance::every(10, 4_096),
                    ..InputScaling::default()
                },
                ..Default::default()
            }
        )
        .is_err()
    );
}

// ---------------------------------------------------------------------------
// HPC scale tests. They read the staged benchmark inputs under
// /ceph/ibmi/abi/oliver/bench/openms4/inputs (MANIFEST.json), which exist only
// on the IBMI nodes, and are therefore `#[ignore]`d. Run one per process on
// ibminode06 to measure peak RSS:
//
//   /usr/bin/time -v <test binary> --ignored --exact hpc_profile_hr_qe_silac_uk222
//
// `OPENMS_BENCH_INPUTS` overrides the directory (the node-local staging copy
// holds the same files with the same names).
// ---------------------------------------------------------------------------

/// The staged path of a benchmark input.
fn bench_input(dataset: &str, name: &str) -> PathBuf {
    match std::env::var("OPENMS_BENCH_INPUTS") {
        Ok(directory) => PathBuf::from(directory).join(name),
        Err(_) => PathBuf::from("/ceph/ibmi/abi/oliver/bench/openms4/inputs")
            .join(dataset)
            .join(name),
    }
}

/// Read one benchmark input the way a TOPP tool does and check its spectrum
/// count against `MANIFEST.json`.
fn hpc_read(dataset: &str, name: &str, spectra: usize) {
    discard_warnings();
    let path = bench_input(dataset, name);
    let file = std::fs::File::open(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let start = std::time::Instant::now();
    let experiment = mzml::read_with_options(BufReader::new(file), &ReadOptions::source()).unwrap();
    let peaks: usize = experiment.spectra.iter().map(|s| s.peaks.len()).sum();
    eprintln!(
        "{name}: {} spectra, {} chromatograms, {peaks} peaks, {:.1} s",
        experiment.spectra.len(),
        experiment.chromatograms.len(),
        start.elapsed().as_secs_f64()
    );
    assert_eq!(experiment.spectra.len(), spectra);
    // The counter and the FeatureFinderCentroided load path read the same file.
    assert_eq!(
        mzml::load_size_with_options(&path, &PeakFileOptions::default(), &ReadOptions::source())
            .unwrap()
            .spectra,
        spectra
    );
}

#[test]
#[ignore = "reads /ceph/ibmi/abi/oliver/bench/openms4/inputs on the IBMI nodes"]
fn hpc_profile_hr_qe_silac_uk222() {
    hpc_read("profile_hr_qe_silac_uk222", "UK222.mzML", 40_856);
}

#[test]
#[ignore = "reads /ceph/ibmi/abi/oliver/bench/openms4/inputs on the IBMI nodes"]
fn hpc_profile_hr_ltqorbitrapxl_ecoli() {
    hpc_read(
        "profile_hr_ltqorbitrapxl_ecoli",
        "20100219_SvNa_SA_Ecoli_preccorrected.mzML",
        34_894,
    );
}

#[test]
#[ignore = "reads /ceph/ibmi/abi/oliver/bench/openms4/inputs on the IBMI nodes"]
fn hpc_centroid_lcms_velos_pxd001819() {
    // The `startTimeStamp="-infinity"` input of the smoke benchmark; zlib.
    hpc_read(
        "centroid_lcms_velos_pxd001819_50amol_r1",
        "50amol_R1.mzML",
        43_745,
    );
}

#[test]
#[ignore = "reads /ceph/ibmi/abi/oliver/bench/openms4/inputs on the IBMI nodes"]
fn hpc_centroid_lcms_qe_silac_uk222_picked() {
    hpc_read(
        "centroid_lcms_qe_silac_uk222_picked",
        "UK222_picked.mzML",
        40_856,
    );
}

#[test]
#[ignore = "reads /ceph/ibmi/abi/oliver/bench/openms4/inputs on the IBMI nodes"]
fn hpc_sanity_inputs_and_the_tool_load_path() {
    discard_warnings();
    for (dataset, name, spectra) in [
        (
            "sanity_profile_velos_phospho_int20000_filtered",
            "20120210_v_QRLi_S1_Phospho_240min_Exp_1_int20000_filtered.mzML",
            730,
        ),
        (
            "sanity_centroid_ltqorbitrapxl_zeitz_small",
            "Zeitz_SIP_13-II_020_small_knubbel.mzML",
            1_929,
        ),
    ] {
        hpc_read(dataset, name, spectra);
        // These two carry valid timestamps, so the strict tool path reads them.
        let path = bench_input(dataset, name);
        assert_eq!(
            FileHandler::load_experiment(&path, &[FileType::MzMl])
                .unwrap()
                .spectra
                .len(),
            spectra
        );
    }
}

#[test]
#[ignore = "reads /ceph/ibmi/abi/oliver/bench/openms4/inputs on the IBMI nodes"]
fn hpc_feature_finder_load_path_on_the_largest_input() {
    discard_warnings();
    let path = bench_input("profile_hr_qe_silac_uk222", "UK222.mzML");
    let start = std::time::Instant::now();
    let experiment = mzml::read_with_load_options(
        BufReader::new(std::fs::File::open(&path).unwrap()),
        &feature_finder_load(),
        &ReadOptions::source(),
    )
    .unwrap();
    // MS1 only, as MANIFEST.json records for this run.
    eprintln!(
        "UK222 MS1 selection: {} spectra, {:.1} s",
        experiment.spectra.len(),
        start.elapsed().as_secs_f64()
    );
    assert_eq!(experiment.spectra.len(), 6_911);
    assert!(experiment.spectra.iter().all(|s| s.ms_level == 1));
}

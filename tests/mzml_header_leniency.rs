#![cfg(feature = "mzml")]
// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source-compatible reading of dangling mzML header references (decision D10).
//!
//! Source `MzMLHandler` resolves `softwareRef` and every data-processing
//! reference through `std::map::operator[]`, so a reference that names no
//! definition becomes an empty `Software` or an empty processing history.
//! `mzml::ReadOptions::source_dangling_references` selects that behaviour; the
//! default stays strict.
//!
//! Evidence: `data/mzml_header_leniency/p2_mzml_leniency_oracle.tsv` is the
//! stdout of the product-sdk oracle driver `../oracle/p2-mzml-leniency`
//! (Debug, core 4fdec46, tier 1 executed differential) on the upstream
//! `PeakPickerHiRes_5_input.mzML` (test-data 0cb15f2) and four synthetic cases.
//! Hashes are in `data/mzml_header_leniency_provenance.json`.

use openms::{
    Error, MSChromatogram, MSExperiment, MSSpectrum, Result,
    concept::log_stream::{LogColor, LogLevel, LogSink, with_thread_local_log},
    format::{
        PeakFileOptions,
        mzml::{self, LoadOptions, ReadOptions, TransformOptions},
    },
    interfaces::MSDataConsumer,
    kernel::DataArray,
    metadata::{DataProcessing, ExperimentalSettings},
};
use std::{
    io::{self, Cursor, Write},
    ops::ControlFlow,
    path::PathBuf,
    sync::{Arc, Mutex},
};

const ORACLE: &str = include_str!("data/mzml_header_leniency/p2_mzml_leniency_oracle.tsv");
const PEAK_PICKER_HI_RES_5: &str =
    include_str!("data/mzml_header_leniency/PeakPickerHiRes_5_input.mzML");
const RECORD_PROCESSING: &str =
    include_str!("data/mzml_header_leniency/record_processing_dangling.mzML");

/// The oracle's cases, in the order `run.sh` passes them to every mode.
const CASES: [(&str, &str); 5] = [
    ("PeakPickerHiRes_5_input", PEAK_PICKER_HI_RES_5),
    (
        "instrument_software_dangling",
        include_str!("data/mzml_header_leniency/instrument_software_dangling.mzML"),
    ),
    (
        "processing_method_software_dangling",
        include_str!("data/mzml_header_leniency/processing_method_software_dangling.mzML"),
    ),
    ("record_processing_dangling", RECORD_PROCESSING),
    (
        "forward_software_reference",
        include_str!("data/mzml_header_leniency/forward_software_reference.mzML"),
    ),
];

fn source() -> ReadOptions {
    ReadOptions {
        source_dangling_references: true,
        ..Default::default()
    }
}

/// `FileHandler::load_experiment_with_options` hands the mzML reader exactly
/// these load options for default `PeakFileOptions`.
fn file_handler_load() -> LoadOptions {
    LoadOptions {
        scientific: PeakFileOptions::default(),
        ..Default::default()
    }
}

/// Send this test thread's warning route to a discarding sink, so the expected
/// dangling-reference warnings do not flood the test output. Only
/// `each_distinct_dangling_reference_warns_once_per_read_and_strict_reads_stay_silent`
/// inspects them.
fn discard_warnings() {
    with_thread_local_log(LogLevel::Warn, |log| {
        log.remove_all_streams()?;
        log.insert(&LogSink::new(io::sink()))
    })
    .unwrap();
}

fn parse_message(result: Result<impl Sized>) -> String {
    match result {
        Err(Error::Parse { message, .. }) => message,
        Err(other) => panic!("expected a parse error, got {other}"),
        Ok(_) => panic!("expected a parse error, got a result"),
    }
}

// Rendering in the oracle driver's line format.

fn text(value: &str) -> &str {
    if value.is_empty() { "<empty>" } else { value }
}
fn settings(out: &mut Vec<String>, mode: &str, label: &str, s: &ExperimentalSettings) {
    let i = &s.instrument;
    out.push(format!(
        "{mode}_instrument\t{label}\tname\t{}\tsoftware\t{}\t{}",
        text(&i.name),
        text(&i.software.name),
        text(&i.software.version)
    ));
}
fn processing(
    out: &mut Vec<String>,
    mode: &str,
    label: &str,
    record: &str,
    id: &str,
    entries: &[Arc<DataProcessing>],
) {
    out.push(format!(
        "{mode}_{record}\t{label}\t{id}\tprocessing\t{}",
        entries.len()
    ));
    for (k, entry) in entries.iter().enumerate() {
        let actions = entry
            .actions
            .iter()
            .map(|action| action.name())
            .collect::<Vec<_>>()
            .join(",");
        out.push(format!(
            "{mode}_{record}_processing\t{label}\t{id}\t{k}\t{}\t{}\t{}",
            text(&entry.software.name),
            text(&entry.software.version),
            if actions.is_empty() { "-" } else { &actions }
        ));
    }
}
fn arrays(out: &mut Vec<String>, mode: &str, label: &str, id: &str, arrays: &[DataArray<f32>]) {
    for array in arrays {
        out.push(format!(
            "{mode}_float_array\t{label}\t{id}\t{}\tprocessing\t{}",
            array.name,
            array.data_processing.len()
        ));
    }
}
fn spectrum(out: &mut Vec<String>, mode: &str, label: &str, s: &MSSpectrum) {
    processing(
        out,
        mode,
        label,
        "spectrum",
        &s.native_id,
        &s.data_processing,
    );
    arrays(out, mode, label, &s.native_id, &s.float_data_arrays);
}
fn chromatogram(out: &mut Vec<String>, mode: &str, label: &str, c: &MSChromatogram) {
    processing(
        out,
        mode,
        label,
        "chromatogram",
        &c.native_id,
        &c.data_processing,
    );
    arrays(out, mode, label, &c.native_id, &c.float_data_arrays);
}
fn experiment(out: &mut Vec<String>, mode: &str, label: &str, e: &MSExperiment) {
    out.push(format!(
        "{mode}\t{label}\tspectra\t{}\tchromatograms\t{}",
        e.spectra.len(),
        e.chromatograms.len()
    ));
    settings(out, mode, label, &e.settings);
    for s in &e.spectra {
        spectrum(out, mode, label, s);
    }
    for c in &e.chromatograms {
        chromatogram(out, mode, label, c);
    }
}

struct Recorder<'a> {
    label: &'a str,
    lines: Vec<String>,
}
impl MSDataConsumer for Recorder<'_> {
    fn set_expected_size(&mut self, spectra: usize, chromatograms: usize) -> Result<()> {
        self.lines.push(format!(
            "transform\t{}\tspectra\t{spectra}\tchromatograms\t{chromatograms}",
            self.label
        ));
        Ok(())
    }
    fn set_experimental_settings(&mut self, s: &ExperimentalSettings) -> Result<()> {
        settings(&mut self.lines, "transform", self.label, s);
        Ok(())
    }
    fn consume_spectrum(&mut self, s: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        spectrum(&mut self.lines, "transform", self.label, s);
        Ok(ControlFlow::Continue(()))
    }
    fn consume_chromatogram(&mut self, c: &mut MSChromatogram) -> Result<ControlFlow<()>> {
        chromatogram(&mut self.lines, "transform", self.label, c);
        Ok(ControlFlow::Continue(()))
    }
}

fn load_lines(label: &str, xml: &str, options: &ReadOptions) -> Result<Vec<String>> {
    let e = mzml::read_with_load_options(Cursor::new(xml), &file_handler_load(), options)?;
    let mut out = Vec::new();
    experiment(&mut out, "load", label, &e);
    Ok(out)
}
fn metadata_lines(label: &str, xml: &str, options: &ReadOptions) -> Result<Vec<String>> {
    let mut load = file_handler_load();
    load.scientific.metadata_only = true;
    let e = mzml::read_with_load_options(Cursor::new(xml), &load, options)?;
    assert_eq!(
        mzml::read_metadata_with_options(Cursor::new(xml), options)?,
        e.settings
    );
    let mut out = Vec::new();
    experiment(&mut out, "metadata", label, &e);
    Ok(out)
}
fn transform_lines(label: &str, xml: &str, options: &ReadOptions) -> Result<Vec<String>> {
    let mut recorder = Recorder {
        label,
        lines: Vec::new(),
    };
    let transform = TransformOptions {
        read: *options,
        ..Default::default()
    };
    mzml::transform_from(
        || Ok(Cursor::new(xml.as_bytes())),
        &mut recorder,
        &transform,
    )?;
    Ok(recorder.lines)
}

#[test]
fn source_option_reproduces_the_executed_cpp_oracle() {
    discard_warnings();
    let mut rust = Vec::new();
    for (label, xml) in CASES {
        rust.extend(load_lines(label, xml, &source()).unwrap());
    }
    for (label, xml) in CASES {
        rust.extend(metadata_lines(label, xml, &source()).unwrap());
    }
    for (label, xml) in CASES {
        rust.extend(transform_lines(label, xml, &source()).unwrap());
    }
    let oracle = ORACLE.lines().collect::<Vec<_>>();
    assert_eq!(rust.len(), oracle.len());
    for (line, (rust, oracle)) in rust.iter().zip(&oracle).enumerate() {
        assert_eq!(rust, oracle, "oracle line {}", line + 1);
    }
}

#[test]
fn strict_default_keeps_every_existing_error() {
    assert!(!ReadOptions::default().source_dangling_references);
    let strict = ReadOptions::default();
    for (label, xml, message) in [
        (
            "PeakPickerHiRes_5_input",
            PEAK_PICKER_HI_RES_5,
            "unresolved softwareRef",
        ),
        (
            "instrument_software_dangling",
            CASES[1].1,
            "unresolved softwareRef",
        ),
        (
            "processing_method_software_dangling",
            CASES[2].1,
            "unresolved softwareRef",
        ),
        (
            "record_processing_dangling",
            RECORD_PROCESSING,
            "unresolved dataProcessingRef",
        ),
        (
            "forward_software_reference",
            CASES[4].1,
            "unresolved softwareRef",
        ),
    ] {
        assert_eq!(
            parse_message(mzml::read(Cursor::new(xml))),
            message,
            "{label}"
        );
        assert_eq!(
            parse_message(load_lines(label, xml, &strict)),
            message,
            "{label}"
        );
        assert_eq!(
            parse_message(transform_lines(label, xml, &strict)),
            message,
            "{label}"
        );
    }
    // With the software reference removed, the dangling default processing
    // reference of the chromatogram list is the next strict error, also on the
    // metadata-only path, which resolves it before stopping.
    let without_software = PEAK_PICKER_HI_RES_5.replace("<softwareRef ref=\"so_in_0\" />", "");
    assert_ne!(without_software, PEAK_PICKER_HI_RES_5);
    assert_eq!(
        parse_message(mzml::read(Cursor::new(&without_software))),
        "unresolved dataProcessingRef"
    );
    assert_eq!(
        parse_message(mzml::read_metadata(Cursor::new(&without_software))),
        "unresolved dataProcessingRef"
    );
}

#[test]
fn source_option_drops_only_the_dangling_references_of_peak_picker_hi_res_5() {
    discard_warnings();
    // Removing the two dangling references by hand gives a document the strict
    // reader accepts. The source option must produce exactly that experiment,
    // on the stream, load-option and compressed-path readers alike.
    const SOFTWARE_REF: &str = "<softwareRef ref=\"so_in_0\" />";
    const DEFAULT_PROCESSING: &str = " defaultDataProcessingRef=\"dp_sp_0\"";
    assert_eq!(PEAK_PICKER_HI_RES_5.matches(SOFTWARE_REF).count(), 1);
    assert_eq!(PEAK_PICKER_HI_RES_5.matches(DEFAULT_PROCESSING).count(), 1);
    // Neither ID is defined anywhere in the upstream fixture.
    assert!(!PEAK_PICKER_HI_RES_5.contains("id=\"so_in_0\""));
    assert!(!PEAK_PICKER_HI_RES_5.contains("id=\"dp_sp_0\""));
    let repaired = PEAK_PICKER_HI_RES_5
        .replace(SOFTWARE_REF, "")
        .replace(DEFAULT_PROCESSING, "");
    let expected = mzml::read(Cursor::new(&repaired)).unwrap();
    assert_eq!(expected.chromatograms.len(), 3);
    let lenient = mzml::read_with_options(Cursor::new(PEAK_PICKER_HI_RES_5), &source()).unwrap();
    assert_eq!(lenient, expected);
    let loaded = mzml::read_with_load_options(
        Cursor::new(PEAK_PICKER_HI_RES_5),
        &file_handler_load(),
        &source(),
    )
    .unwrap();
    assert_eq!(
        loaded,
        mzml::read_with_load_options(Cursor::new(&repaired), &file_handler_load(), &source())
            .unwrap()
    );
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/mzml_header_leniency/PeakPickerHiRes_5_input.mzML");
    let from_path = mzml::load_with_options(&path, &file_handler_load(), &source()).unwrap();
    assert_eq!(from_path.spectra, loaded.spectra);
    assert_eq!(from_path.chromatograms, loaded.chromatograms);
    assert_eq!(from_path.settings.instrument, loaded.settings.instrument);
    assert!(mzml::load_with_options(&path, &file_handler_load(), &ReadOptions::default()).is_err());
}

#[test]
fn source_option_leaves_other_references_and_malformed_ids_strict() {
    discard_warnings();
    const HEADER: &str = r#"<fileDescription><fileContent/><sourceFileList count="1"><sourceFile id="sf" name="x.raw" location="file:///data"/></sourceFileList></fileDescription>
<sampleList count="1"><sample id="sa" name="s"/></sampleList>
<softwareList count="1"><software id="sw" version="1"><cvParam accession="MS:1000799" name="custom unreleased software tool" value="tool"/></software></softwareList>
<instrumentConfigurationList count="1"><instrumentConfiguration id="ic"/></instrumentConfigurationList>
<dataProcessingList count="1"><dataProcessing id="dp"><processingMethod order="0" softwareRef="sw"><cvParam accession="MS:1000035" name="peak picking"/></processingMethod></dataProcessing></dataProcessingList>"#;
    let document = |run: &str, spectrum: &str| {
        format!(
            r#"<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">{HEADER}<run id="r" {run}><spectrumList count="1" defaultDataProcessingRef="dp">{spectrum}</spectrumList></run></mzML>"#
        )
    };
    let valid = document(
        r#"sampleRef="sa" defaultInstrumentConfigurationRef="ic""#,
        r#"<spectrum id="a" defaultArrayLength="0" sourceFileRef="sf"/>"#,
    );
    assert!(mzml::read(Cursor::new(&valid)).is_ok());
    let scan = |id: &str| {
        valid.replacen(
            r#"defaultArrayLength="0" sourceFileRef="sf"/>"#,
            &format!(
                r#"defaultArrayLength="0"><scanList count="1"><scan instrumentConfigurationRef="{id}"/></scanList></spectrum>"#
            ),
            1,
        )
    };
    assert!(mzml::read_with_options(Cursor::new(scan("ic")), &source()).is_ok());
    for (changed, message) in [
        (
            valid.replacen(r#"sourceFileRef="sf""#, r#"sourceFileRef="absent""#, 1),
            "unresolved sourceFileRef",
        ),
        (
            valid.replacen(r#"sampleRef="sa""#, r#"sampleRef="absent""#, 1),
            "unresolved sampleRef",
        ),
        (
            valid.replacen(
                r#"defaultInstrumentConfigurationRef="ic""#,
                r#"defaultInstrumentConfigurationRef="absent""#,
                1,
            ),
            "unresolved defaultInstrumentConfigurationRef",
        ),
        (scan("absent"), "unresolved scan instrumentConfigurationRef"),
        (
            valid.replacen(
                r#"defaultArrayLength="0" sourceFileRef="sf"/>"#,
                r#"defaultArrayLength="0"><referenceableParamGroupRef ref="absent"/></spectrum>"#,
                1,
            ),
            "unknown parameter group absent",
        ),
        (
            valid.replacen(
                r#"defaultDataProcessingRef="dp""#,
                r#"defaultDataProcessingRef="1malformed""#,
                1,
            ),
            "invalid parameter group ID/IDREF",
        ),
        (
            valid.replacen(r#"softwareRef="sw""#, r#"softwareRef="""#, 1),
            "invalid parameter group ID/IDREF",
        ),
    ] {
        assert_ne!(changed, valid, "{message}");
        assert_eq!(
            parse_message(mzml::read_with_options(Cursor::new(&changed), &source())),
            message
        );
    }
}

#[test]
fn dangling_processing_references_in_every_position_become_empty_histories() {
    discard_warnings();
    let e = mzml::read_with_options(Cursor::new(RECORD_PROCESSING), &source()).unwrap();
    let spectra = e
        .spectra
        .iter()
        .map(|s| s.data_processing.len())
        .collect::<Vec<_>>();
    assert_eq!(spectra, [0, 0, 1, 0]);
    let aux = &e.spectra[3].float_data_arrays;
    assert_eq!(aux.len(), 1);
    assert_eq!(aux[0].name, "aux");
    assert!(aux[0].data_processing.is_empty());
    assert_eq!(e.spectra[3].peaks.len(), 2);
    let chromatograms = e
        .chromatograms
        .iter()
        .map(|c| c.data_processing.len())
        .collect::<Vec<_>>();
    assert_eq!(chromatograms, [1, 0]);
    // Resolved references keep sharing the definition's handles.
    assert!(Arc::ptr_eq(
        &e.spectra[2].data_processing[0],
        &e.chromatograms[0].data_processing[0]
    ));
}

#[test]
fn placeholder_marker_lookup_neither_panics_nor_accepts_a_dangling_software() {
    discard_warnings();
    let empty = MSExperiment {
        spectra: vec![MSSpectrum {
            native_id: "scan=1".into(),
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut xml = Vec::new();
    mzml::write(&mut xml, &empty).unwrap();
    let xml = String::from_utf8(xml).unwrap();
    assert!(xml.contains("openms-rust:empty-processing-history"));
    let start = xml.find("softwareRef=\"").unwrap() + "softwareRef=\"".len();
    // A whitespace-padded reference names the placeholder software after ID
    // normalization. The marker lookup indexed the raw attribute and panicked;
    // it now uses the normalized ID, under both policies.
    let padded = format!("{} {}", &xml[..start], &xml[start..]);
    for options in [ReadOptions::default(), source()] {
        assert_eq!(
            mzml::read_with_options(Cursor::new(&padded), &options).unwrap(),
            empty
        );
    }
    // A dangling reference cannot name the exact placeholder payload, so the
    // marker is refused rather than silently normalized.
    let end = start + xml[start..].find('"').unwrap();
    let dangling = format!("{}absent{}", &xml[..start], &xml[end..]);
    assert_eq!(
        parse_message(mzml::read_with_options(Cursor::new(&dangling), &source())),
        "processing marker contains additional or discarded XML payload"
    );
    assert_eq!(
        parse_message(mzml::read(Cursor::new(&dangling))),
        "unresolved softwareRef"
    );
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

#[test]
fn each_distinct_dangling_reference_warns_once_per_read_and_strict_reads_stay_silent() {
    // The warning route of this test thread only; other tests keep stderr.
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

    mzml::read_with_options(Cursor::new(PEAK_PICKER_HI_RES_5), &source()).unwrap();
    assert_eq!(
        take(),
        "Warning: mzML softwareRef 'so_in_0' names no definition; source-compatible reading uses empty software.\n\
         Warning: mzML dataProcessingRef 'dp_sp_0' names no definition; source-compatible reading uses an empty processing history.\n"
    );

    // Three spectra name the same dangling ID, and the list default is a
    // second one: two lines, no repetition summary from the log cache.
    let repeated = RECORD_PROCESSING.replace(
        "<spectrum id=\"scan=3\" index=\"2\" defaultArrayLength=\"0\" dataProcessingRef=\"dp\"/>",
        "<spectrum id=\"scan=3\" index=\"2\" defaultArrayLength=\"0\" dataProcessingRef=\"absent_explicit\"/>\
         <spectrum id=\"scan=5\" index=\"4\" defaultArrayLength=\"0\" dataProcessingRef=\"absent_explicit\"/>",
    );
    assert_ne!(repeated, RECORD_PROCESSING);
    let e = mzml::read_with_options(Cursor::new(&repeated), &source()).unwrap();
    assert!(e.spectra.iter().all(|s| s.data_processing.is_empty()));
    let warnings = take();
    let lines = warnings.lines().collect::<Vec<_>>();
    assert_eq!(
        lines
            .iter()
            .map(|line| line.split('\'').nth(1).unwrap())
            .collect::<Vec<_>>(),
        [
            "absent_default",
            "absent_explicit",
            "absent_array",
            "absent_aux",
            "absent_chrom"
        ],
        "{warnings}"
    );

    // Each read starts with a fresh record, and a strict read warns nothing.
    mzml::read_with_options(Cursor::new(PEAK_PICKER_HI_RES_5), &source()).unwrap();
    assert_eq!(take().lines().count(), 2);
    assert!(mzml::read(Cursor::new(PEAK_PICKER_HI_RES_5)).is_err());
    assert_eq!(take(), "");

    with_thread_local_log(LogLevel::Warn, |log| {
        log.remove_all_streams()?;
        log.set_color(Some(LogColor::Yellow));
        log.insert(&LogSink::stderr())
    })
    .unwrap();
}

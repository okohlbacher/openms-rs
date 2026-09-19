// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports every `START_SECTION` of `MSDataWritingConsumer_test.cpp` (pinned
//! revision bc9cc12, 11 sections) plus native boundary tests.
//!
//! The upstream test supplies no expected values at all: ten of its eleven
//! sections contain only `// TODO`, the eleventh constructs an abstract class
//! with a constructor it does not have, and `executables.cmake:239` has the
//! whole test commented out so it is never built. Every expectation here is
//! therefore derived from the header and `.cpp` (tier 3) or from mzML itself
//! (tier 4); see `tests/data/ms_data_writing_consumer_provenance.json`.

#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::ms_data_writing_consumer::{
    CountPolicy, MSDataWritingConsumer, MSDataWritingProcessor, NoopMSDataWritingConsumer,
    PlainMSDataWritingConsumer, PlainProcessor, ReferencePolicy, WritingLimits,
};
use openms::format::mzml::{self, WriteOptions};
use openms::interfaces::MSDataConsumer;
use openms::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D,
};
use openms::metadata::{
    DataProcessing, ExperimentalSettings, ProcessingAction, Software, SourceFile,
};
use openms::system::file::TempDir;
use std::sync::Arc;

fn spectrum(native_id: &str, rt: f64, mz: f64) -> MSSpectrum {
    MSSpectrum {
        native_id: native_id.to_owned(),
        rt,
        ms_level: 1,
        peaks: vec![
            Peak1D {
                mz,
                intensity: 10.0,
            },
            Peak1D {
                mz: mz + 1.0,
                intensity: 20.0,
            },
        ],
        ..Default::default()
    }
}

fn chromatogram(native_id: &str, rt: f64) -> MSChromatogram {
    MSChromatogram {
        native_id: native_id.to_owned(),
        peaks: vec![
            ChromatogramPeak { rt, intensity: 3.0 },
            ChromatogramPeak {
                rt: rt + 1.0,
                intensity: 4.0,
            },
        ],
        ..Default::default()
    }
}

fn written(consumer: PlainMSDataWritingConsumer<Vec<u8>>) -> String {
    String::from_utf8(consumer.finish().expect("finish")).expect("UTF-8 mzML")
}

// ---------------------------------------------------------------------------
// START_SECTION(MSDataWritingConsumer())
// The one upstream section with an assertion macro (TEST_NOT_EQUAL on the
// pointer). It cannot compile: the class is abstract and has no default
// constructor. The constructed state the source's initialiser list sets is
// asserted here instead.
// ---------------------------------------------------------------------------
#[test]
fn a_new_consumer_starts_in_the_source_initial_state() {
    let consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    assert!(!consumer.started_writing());
    assert!(!consumer.writing_spectra());
    assert!(!consumer.writing_chromatograms());
    assert_eq!(consumer.spectra_written(), 0);
    assert_eq!(consumer.chromatograms_written(), 0);
    assert_eq!(consumer.expected_size(), (0, 0));
    assert!(consumer.additional_data_processing().is_none());
    assert_eq!(consumer.limits(), WritingLimits::default());
    assert_eq!(consumer.count_policy(), CountPolicy::Checked);
    assert_eq!(consumer.processor(), &PlainProcessor);
    assert_eq!(consumer.settings(), &ExperimentalSettings::default());
    assert!(!consumer.write_options().zlib_compression);
}

// ---------------------------------------------------------------------------
// START_SECTION(~MSDataWritingConsumer())
// START_SECTION((virtual ~MSDataWritingConsumer()))
// Two upstream sections, both `// TODO`. The destructor calls doCleanup_,
// which writes a footer only when writing actually started.
// ---------------------------------------------------------------------------
#[test]
fn a_consumer_that_received_nothing_writes_nothing() {
    let consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    assert_eq!(consumer.finish().unwrap(), Vec::<u8>::new());
}

#[test]
fn cleanup_closes_the_open_list_and_the_document() {
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new())
        .with_count_policy(CountPolicy::SourceInconsistent);
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    let text = written(consumer);
    // `doCleanup_` closes the open list, then the document, then hands over to
    // `MzMLHandlerHelper::writeFooter_`, whose inherited `PeakFileOptions` has
    // `write_index_` true: the index, its offset and the checksum follow, and
    // the `indexedmzML` wrapper the header opened is closed last.
    assert!(
        text.contains("</spectrumList>\n</run></mzML>\n<indexList count=\"1\">\n"),
        "{text}"
    );
    assert!(text.contains("<index name=\"spectrum\">\n"), "{text}");
    assert!(text.contains("<offset idRef=\"scan=1\">"), "{text}");
    assert!(
        text.ends_with("</fileChecksum>\n</indexedmzML>\n"),
        "{text}"
    );

    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new())
        .with_count_policy(CountPolicy::SourceInconsistent);
    consumer
        .consume_chromatogram(&mut chromatogram("chrom=1", 1.0))
        .unwrap();
    let text = written(consumer);
    assert!(
        text.contains("</chromatogramList>\n</run></mzML>\n<indexList count=\"1\">\n"),
        "{text}"
    );
    assert!(text.contains("<index name=\"chromatogram\">\n"), "{text}");
    assert!(text.contains("<offset idRef=\"chrom=1\">"), "{text}");
    assert!(
        text.ends_with("</fileChecksum>\n</indexedmzML>\n"),
        "{text}"
    );

    // Dropping without finishing leaves the document unclosed, deliberately.
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    drop(consumer);
}

// ---------------------------------------------------------------------------
// START_SECTION((virtual void setExperimentalSettings(const ExperimentalSettings &exp)))
// Upstream `// TODO`.
// ---------------------------------------------------------------------------
#[test]
fn experimental_settings_reach_the_header_and_freeze_with_it() {
    let settings = ExperimentalSettings {
        source_files: vec![SourceFile {
            name: "raw.RAW".into(),
            path: "file://.".into(),
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new())
        .with_count_policy(CountPolicy::SourceInconsistent);
    consumer.set_experimental_settings(&settings).unwrap();
    assert_eq!(consumer.settings().source_files.len(), 1);
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    // Once the header exists the settings cannot change, because the header
    // was already rendered from them.
    let error = consumer
        .set_experimental_settings(&ExperimentalSettings::default())
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    let text = written(consumer);
    assert!(text.contains("raw.RAW"), "{text}");
    assert!(text.contains("<sourceFileList"), "{text}");
}

// ---------------------------------------------------------------------------
// START_SECTION((virtual void setExpectedSize(Size expectedSpectra, Size expectedChromatograms)))
// Upstream `// TODO`. The counts go into the two list `count` attributes and
// the source enforces nothing.
// ---------------------------------------------------------------------------
#[test]
fn expected_sizes_become_the_list_counts() {
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(2, 1).unwrap();
    assert_eq!(consumer.expected_size(), (2, 1));
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    consumer
        .consume_spectrum(&mut spectrum("scan=2", 2.0, 200.0))
        .unwrap();
    // A list tag is open, so the announcement can no longer change.
    let error = consumer.set_expected_size(5, 5).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    consumer
        .consume_chromatogram(&mut chromatogram("chrom=1", 1.0))
        .unwrap();
    let text = written(consumer);
    assert!(text.contains("<spectrumList count=\"2\""), "{text}");
    assert!(text.contains("<chromatogramList count=\"1\""), "{text}");
}

#[test]
fn a_wrong_announcement_is_reported_rather_than_shipped_silently() {
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(7, 0).unwrap();
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    let error = consumer.finish().unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");

    // The source behaviour is available explicitly: the file is written with
    // the wrong count and no error is raised.
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new())
        .with_count_policy(CountPolicy::SourceInconsistent);
    consumer.set_expected_size(7, 0).unwrap();
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    let text = written(consumer);
    assert!(text.contains("<spectrumList count=\"7\""), "{text}");
    // The declared count disagrees with the one spectrum present, exactly the
    // inconsistency the source's class note warns about.
    assert_eq!(text.matches("<spectrum id=").count(), 1);
}

// ---------------------------------------------------------------------------
// START_SECTION((virtual void consumeSpectrum(SpectrumType &s)))
// Upstream `// TODO`.
// ---------------------------------------------------------------------------
#[test]
fn spectra_are_streamed_and_read_back_unchanged() {
    let mut source = MSExperiment::new();
    source.spectra.push(spectrum("scan=1", 1.0, 100.0));
    source.spectra.push(spectrum("scan=2", 2.5, 200.0));
    source.spectra.push(spectrum("scan=3", 3.5, 300.0));

    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(source.spectra.len(), 0).unwrap();
    for record in &source.spectra {
        let mut copy = record.clone();
        consumer.consume_spectrum(&mut copy).unwrap();
        // The source copies the record, so the caller's value is untouched.
        assert_eq!(&copy, record);
    }
    assert_eq!(consumer.spectra_written(), 3);
    let text = written(consumer);

    // Every record carries its position, not the zero of its own render.
    assert!(text.contains("index=\"0\""), "{text}");
    assert!(text.contains("index=\"1\""), "{text}");
    assert!(text.contains("index=\"2\""), "{text}");
    let loaded = mzml::read(text.as_bytes()).unwrap();
    assert_eq!(loaded.spectra.len(), 3);
    assert_eq!(loaded.spectra, source.spectra);

    // Byte-identical to the **indexed** whole-document writer, which is what
    // `MzMLFile::store` and this consumer both are: the record blocks come from
    // the same encoder, the header differs only in the list `count` this test
    // announces correctly, and the index entries, `indexListOffset` and SHA-1
    // `fileChecksum` fall out of the same byte positions. This is the strongest
    // statement the streaming path can make - the file a caller gets from
    // streaming is the file it would have got from holding the experiment.
    let mut whole = Vec::new();
    mzml::write(&mut whole, &source).unwrap();
    assert_eq!(text, String::from_utf8(whole).unwrap());
    assert!(text.contains("<indexedmzML "), "{text}");
}

#[test]
fn spectra_after_chromatograms_are_refused() {
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(0, 1).unwrap();
    consumer
        .consume_chromatogram(&mut chromatogram("chrom=1", 1.0))
        .unwrap();
    let error = consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    assert_eq!(consumer.spectra_written(), 0);
    let text = written(consumer);
    assert!(!text.contains("<spectrumList"), "{text}");
}

#[test]
fn an_empty_native_id_is_filled_in_with_the_record_position() {
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(2, 1).unwrap();
    consumer
        .consume_spectrum(&mut MSSpectrum::default())
        .unwrap();
    consumer
        .consume_spectrum(&mut MSSpectrum::default())
        .unwrap();
    consumer
        .consume_chromatogram(&mut MSChromatogram::default())
        .unwrap();
    let text = written(consumer);
    assert!(
        text.contains("<spectrum id=\"index=0\" index=\"0\""),
        "{text}"
    );
    assert!(
        text.contains("<spectrum id=\"index=1\" index=\"1\""),
        "{text}"
    );
    assert!(
        text.contains("<chromatogram id=\"chromatogram=0\" index=\"0\""),
        "{text}"
    );
    mzml::read(text.as_bytes()).unwrap();
}

#[test]
fn a_repeated_native_identifier_is_refused() {
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(2, 0).unwrap();
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    let error = consumer
        .consume_spectrum(&mut spectrum("scan=1", 2.0, 200.0))
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    assert_eq!(consumer.spectra_written(), 1);
}

#[test]
fn a_non_ascii_native_identifier_is_written_and_read_back() {
    let id = "scan=\u{65e5}\u{672c}\u{8a9e}";
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(1, 0).unwrap();
    consumer
        .consume_spectrum(&mut spectrum(id, 1.0, 100.0))
        .unwrap();
    let text = written(consumer);
    let loaded = mzml::read(text.as_bytes()).unwrap();
    assert_eq!(loaded.spectra.len(), 1);
    assert_eq!(loaded.spectra[0].native_id, id);
}

// ---------------------------------------------------------------------------
// START_SECTION((virtual void consumeChromatogram(ChromatogramType &c)))
// Upstream `// TODO`.
// ---------------------------------------------------------------------------
#[test]
fn a_chromatogram_closes_an_open_spectrum_list() {
    let mut source = MSExperiment::new();
    source.spectra.push(spectrum("scan=1", 1.0, 100.0));
    source.chromatograms.push(chromatogram("chrom=1", 1.0));
    source.chromatograms.push(chromatogram("chrom=2", 5.0));

    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(1, 2).unwrap();
    let mut record = source.spectra[0].clone();
    consumer.consume_spectrum(&mut record).unwrap();
    assert!(consumer.writing_spectra());
    for record in &source.chromatograms {
        let mut copy = record.clone();
        consumer.consume_chromatogram(&mut copy).unwrap();
        assert_eq!(&copy, record);
    }
    assert!(!consumer.writing_spectra());
    assert!(consumer.writing_chromatograms());
    assert_eq!(consumer.chromatograms_written(), 2);
    let text = written(consumer);

    let spectrum_list_close = text.find("</spectrumList>").expect("list closed");
    let chromatogram_list_open = text.find("<chromatogramList").expect("list opened");
    assert!(spectrum_list_close < chromatogram_list_open, "{text}");
    let loaded = mzml::read(text.as_bytes()).unwrap();
    assert_eq!(loaded.spectra, source.spectra);
    assert_eq!(loaded.chromatograms, source.chromatograms);
}

// ---------------------------------------------------------------------------
// START_SECTION((MSDataWritingConsumer(std::string filename)))
// Upstream `// TODO`. The source's only constructor; it truncates the file.
// ---------------------------------------------------------------------------
#[test]
fn the_filename_constructor_creates_and_truncates() {
    // A directory of its own. Under a fixed name in the shared temporary
    // directory, a concurrent run of this binary removed the file between
    // `finish` and `load`.
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("out.mzML");
    std::fs::write(&path, b"stale contents that must not survive").unwrap();

    let mut consumer = PlainMSDataWritingConsumer::create_plain(&path).unwrap();
    consumer.set_expected_size(1, 0).unwrap();
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    consumer.finish().unwrap();
    let loaded = mzml::load(&path).unwrap();
    assert_eq!(loaded.spectra.len(), 1);
    assert_eq!(loaded.spectra[0].native_id, "scan=1");

    let error =
        MSDataWritingConsumer::create(dir.path().join("missing").join("out.mzML"), PlainProcessor)
            .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error}");
}

// ---------------------------------------------------------------------------
// START_SECTION((virtual void addDataProcessing(DataProcessing d)))
// Upstream `// TODO`. The entry is appended to every record written after it.
// ---------------------------------------------------------------------------
#[test]
fn an_added_data_processing_entry_reaches_every_record() {
    let processing = DataProcessing {
        actions: [ProcessingAction::Smoothing].into_iter().collect(),
        ..Default::default()
    };
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(2, 1).unwrap();
    consumer.add_data_processing(processing.clone()).unwrap();
    assert_eq!(
        consumer.additional_data_processing().map(|p| p.as_ref()),
        Some(&processing)
    );
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    // A second call replaces the entry rather than adding another, but the
    // header has been written by now, so it is refused.
    let error = consumer
        .add_data_processing(processing.clone())
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    consumer
        .consume_spectrum(&mut spectrum("scan=2", 2.0, 200.0))
        .unwrap();
    consumer
        .consume_chromatogram(&mut chromatogram("chrom=1", 1.0))
        .unwrap();
    let text = written(consumer);

    let loaded = mzml::read(text.as_bytes()).unwrap();
    for record in &loaded.spectra {
        assert_eq!(record.data_processing.len(), 1);
        assert_eq!(record.data_processing[0].actions, processing.actions);
    }
    assert_eq!(loaded.chromatograms.len(), 1);
    assert_eq!(loaded.chromatograms[0].data_processing.len(), 1);

    // Replacing the entry before the header is written does take effect.
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(1, 0).unwrap();
    consumer
        .add_data_processing(DataProcessing::default())
        .unwrap();
    consumer.add_data_processing(processing.clone()).unwrap();
    assert_eq!(
        consumer.additional_data_processing().map(|p| p.as_ref()),
        Some(&processing)
    );
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    let loaded = mzml::read(written(consumer).as_bytes()).unwrap();
    assert_eq!(loaded.spectra[0].data_processing.len(), 1);
}

// ---------------------------------------------------------------------------
// START_SECTION((virtual Size getNrSpectraWritten()))
// START_SECTION((virtual Size getNrChromatogramsWritten()))
// Both upstream `// TODO`.
// ---------------------------------------------------------------------------
#[test]
fn the_counters_follow_the_records_actually_written() {
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(2, 2).unwrap();
    assert_eq!(consumer.spectra_written(), 0);
    assert_eq!(consumer.chromatograms_written(), 0);
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    assert_eq!(consumer.spectra_written(), 1);
    consumer
        .consume_spectrum(&mut spectrum("scan=2", 2.0, 200.0))
        .unwrap();
    assert_eq!(consumer.spectra_written(), 2);
    assert_eq!(consumer.chromatograms_written(), 0);
    consumer
        .consume_chromatogram(&mut chromatogram("chrom=1", 1.0))
        .unwrap();
    consumer
        .consume_chromatogram(&mut chromatogram("chrom=2", 2.0))
        .unwrap();
    assert_eq!(consumer.spectra_written(), 2);
    assert_eq!(consumer.chromatograms_written(), 2);
    consumer.finish().unwrap();
}

// ---------------------------------------------------------------------------
// Native tests beyond the upstream sections.
// ---------------------------------------------------------------------------

/// The template-method hook, as a derived class would define it.
struct Halving {
    spectra: usize,
    chromatograms: usize,
}

impl MSDataWritingProcessor for Halving {
    fn process_spectrum(&mut self, spectrum: &mut MSSpectrum) -> openms::Result<()> {
        self.spectra += 1;
        for peak in &mut spectrum.peaks {
            peak.intensity /= 2.0;
        }
        Ok(())
    }
    fn process_chromatogram(&mut self, chromatogram: &mut MSChromatogram) -> openms::Result<()> {
        self.chromatograms += 1;
        chromatogram.native_id.push_str("_processed");
        Ok(())
    }
}

#[test]
fn the_processor_hook_transforms_the_copy_not_the_caller_record() {
    let mut consumer = MSDataWritingConsumer::new(
        Vec::new(),
        Halving {
            spectra: 0,
            chromatograms: 0,
        },
    );
    consumer.set_expected_size(1, 1).unwrap();
    let mut record = spectrum("scan=1", 1.0, 100.0);
    consumer.consume_spectrum(&mut record).unwrap();
    assert_eq!(record.peaks[0].intensity, 10.0);
    let mut chrom = chromatogram("chrom=1", 1.0);
    consumer.consume_chromatogram(&mut chrom).unwrap();
    assert_eq!(chrom.native_id, "chrom=1");
    assert_eq!(consumer.processor().spectra, 1);
    assert_eq!(consumer.processor().chromatograms, 1);
    let text = String::from_utf8(consumer.finish().unwrap()).unwrap();
    let loaded = mzml::read(text.as_bytes()).unwrap();
    assert_eq!(loaded.spectra[0].peaks[0].intensity, 5.0);
    assert_eq!(loaded.spectra[0].peaks[1].intensity, 10.0);
    assert_eq!(loaded.chromatograms[0].native_id, "chrom=1_processed");
}

/// A hook that refuses, which the source's `void` hooks cannot express.
struct Refusing;

impl MSDataWritingProcessor for Refusing {
    fn process_spectrum(&mut self, _spectrum: &mut MSSpectrum) -> openms::Result<()> {
        Err(Error::Unsupported("no spectra here".into()))
    }
    fn process_chromatogram(&mut self, _chromatogram: &mut MSChromatogram) -> openms::Result<()> {
        Err(Error::Unsupported("no chromatograms here".into()))
    }
}

#[test]
fn a_refusing_processor_aborts_the_record_before_anything_is_written() {
    let mut consumer = MSDataWritingConsumer::new(Vec::new(), Refusing);
    let error = consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error}");
    let error = consumer
        .consume_chromatogram(&mut chromatogram("chrom=1", 1.0))
        .unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error}");
    assert!(!consumer.started_writing());
    assert_eq!(consumer.finish().unwrap(), Vec::<u8>::new());
}

#[test]
fn a_record_referencing_an_undeclared_header_element_is_refused() {
    // Only the first record contributes to the header, so a second record
    // with its own source file has nothing to point at. The source emits a
    // dangling sourceFileRef instead; see MzMLHandler.cpp:5254.
    let mut first = spectrum("scan=1", 1.0, 100.0);
    let mut second = spectrum("scan=2", 2.0, 200.0);
    second.source_file = SourceFile {
        name: "other.RAW".into(),
        path: "file://.".into(),
        ..Default::default()
    };
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(2, 0).unwrap();
    consumer.consume_spectrum(&mut first).unwrap();
    let error = consumer.consume_spectrum(&mut second).unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error}");
    assert_eq!(consumer.spectra_written(), 1);

    // A differing data-processing history is refused for the same reason.
    let mut third = spectrum("scan=3", 3.0, 300.0);
    third.data_processing = vec![Arc::new(DataProcessing {
        actions: [ProcessingAction::Smoothing].into_iter().collect(),
        ..Default::default()
    })];
    let error = consumer.consume_spectrum(&mut third).unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error}");

    // The same source file on every record is fine: the header declares it.
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(2, 0).unwrap();
    let mut a = spectrum("scan=1", 1.0, 100.0);
    a.source_file = second.source_file.clone();
    let mut b = spectrum("scan=2", 2.0, 200.0);
    b.source_file = second.source_file.clone();
    consumer.consume_spectrum(&mut a).unwrap();
    consumer.consume_spectrum(&mut b).unwrap();
    let text = written(consumer);
    let loaded = mzml::read(text.as_bytes()).unwrap();
    assert_eq!(loaded.spectra.len(), 2);
    assert_eq!(loaded.spectra[0].source_file.name, "other.RAW");
}

/// A second software with a name and a version, for the histories below.
fn history(name: &str, version: &str) -> Vec<Arc<DataProcessing>> {
    vec![Arc::new(DataProcessing {
        software: Software {
            name: name.to_owned(),
            version: version.to_owned(),
            ..Default::default()
        },
        actions: [ProcessingAction::Smoothing].into_iter().collect(),
        ..Default::default()
    })]
}

/// Two histories that differ only in their software render the same
/// `dataProcessingList`, and must still not share one header.
///
/// The rendered `processingMethod` names its software by the history's
/// *position* (`so_dp_<history>_<method>`), which is zero for the single
/// record of every per-record render, so the `dataProcessingList` text of two
/// such records is identical and only the `softwareList` differs. Without that
/// list in the comparison the second record would be written silently under
/// the first record's software — the reader would then report a
/// `PeakPickerHiRes 1.0` step for data a `PeakPickerHiRes 2.0` produced.
#[test]
fn two_histories_differing_only_in_their_software_are_not_one_header() {
    let mut first = spectrum("scan=1", 1.0, 100.0);
    first.data_processing = history("PeakPickerHiRes", "1.0");
    let mut second = spectrum("scan=2", 2.0, 200.0);
    second.data_processing = history("PeakPickerHiRes", "2.0");

    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(2, 0).unwrap();
    consumer.consume_spectrum(&mut first).unwrap();
    let error = consumer.consume_spectrum(&mut second).unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error}");
    assert_eq!(consumer.spectra_written(), 1);
}

// ---------------------------------------------------------------------------
// ReferencePolicy: the source's dangling header references.
// ---------------------------------------------------------------------------

#[test]
fn the_reference_policy_defaults_to_refusing() {
    let consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    assert_eq!(consumer.reference_policy(), ReferencePolicy::Checked);
    let consumer = consumer.with_reference_policy(ReferencePolicy::SourceDangling);
    assert_eq!(consumer.reference_policy(), ReferencePolicy::SourceDangling);
}

/// Under the source policy a record whose history the header does not declare
/// is written with the reference the source writes: `dp_sp_<index>`, numbered
/// by the record's position in the stream and naming nothing.
///
/// `MzMLHandler.cpp:5258-5272` with `dps_` holding the one entry
/// `writeHeader_` filled it with: the search for a matching entry fails, so
/// `dp_ref_num` keeps the record's own index. Executed against the C++ Release
/// build on the `refs` fixture (`../oracle/p4-lowmemory`,
/// `logs/closediff2_06.log`, case `refs_low`), where records 1 and 2 of five
/// come out as `dataProcessingRef="dp_sp_1"` and `"dp_sp_2"`.
#[test]
fn the_source_policy_writes_the_dangling_processing_reference() {
    let mut first = spectrum("scan=1", 1.0, 100.0);
    first.data_processing = history("PeakPickerHiRes", "1.0");
    let mut second = spectrum("scan=2", 2.0, 200.0);
    second.data_processing = history("PeakPickerHiRes", "2.0");
    let mut third = spectrum("scan=3", 3.0, 300.0);
    third.data_processing = first.data_processing.clone();

    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new())
        .with_reference_policy(ReferencePolicy::SourceDangling);
    consumer.set_expected_size(3, 0).unwrap();
    consumer.consume_spectrum(&mut first).unwrap();
    consumer.consume_spectrum(&mut second).unwrap();
    consumer.consume_spectrum(&mut third).unwrap();
    assert_eq!(consumer.spectra_written(), 3);
    let text = written(consumer);
    let tags: Vec<&str> = text
        .match_indices("<spectrum id=")
        .map(|(at, _)| {
            let rest = &text[at..];
            &rest[..rest.find('>').expect("a start tag")]
        })
        .collect();
    assert_eq!(tags.len(), 3, "{text:.600}");
    // The first record IS the header, so it needs no reference of its own.
    assert!(!tags[0].contains("dataProcessingRef="), "{}", tags[0]);
    assert!(
        tags[1].ends_with(" dataProcessingRef=\"dp_sp_1\""),
        "{}",
        tags[1]
    );
    // The third record's history is the header's again, so the source writes
    // no reference and it inherits the list's default.
    assert!(!tags[2].contains("dataProcessingRef="), "{}", tags[2]);
    // The identifier dangles, which is the whole point, and cannot collide
    // with anything this writer declares.
    assert!(
        !text.contains("<dataProcessing id=\"dp_sp_1\""),
        "{text:.600}"
    );
    assert!(text.contains("<dataProcessing id=\"dp_00000000000000000000\""));
}

/// The source renumbers a `sourceFileRef` by the record's position for **every**
/// record after the first that carries one, whether or not the source file is
/// the first record's (`MzMLHandler.cpp:5252-5255`, which does not look at
/// `dps` or at the header at all).
///
/// Executed on the `refs` fixture: the C++ low-memory output carries
/// `sourceFileRef="sf_sp_0"` … `"sf_sp_4"` over five records against a header
/// declaring one record source file. Record 3's source file is not the first
/// record's and record 4's is, and both get `sf_sp_3` and `sf_sp_4` just the
/// same.
#[test]
fn the_source_policy_renumbers_every_later_source_file_reference() {
    let one = SourceFile {
        name: "part_one.RAW".into(),
        path: "file://.".into(),
        ..Default::default()
    };
    let two = SourceFile {
        name: "part_two.RAW".into(),
        path: "file://.".into(),
        ..Default::default()
    };
    let mut first = spectrum("scan=1", 1.0, 100.0);
    first.source_file = one.clone();
    let mut second = spectrum("scan=2", 2.0, 200.0);
    second.source_file = one;
    let mut third = spectrum("scan=3", 3.0, 300.0);
    third.source_file = two;

    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new())
        .with_reference_policy(ReferencePolicy::SourceDangling);
    consumer.set_expected_size(3, 0).unwrap();
    consumer.consume_spectrum(&mut first).unwrap();
    consumer.consume_spectrum(&mut second).unwrap();
    consumer.consume_spectrum(&mut third).unwrap();
    let text = written(consumer);
    assert!(text.contains(" sourceFileRef=\"sf_00000000000000000000\""));
    assert!(text.contains(" sourceFileRef=\"sf_sp_1\""), "{text:.600}");
    assert!(text.contains(" sourceFileRef=\"sf_sp_2\""), "{text:.600}");
    assert!(!text.contains("<sourceFile id=\"sf_sp_"), "{text:.600}");
}

/// `SourceDangling` decides "differs from the first record's" by content,
/// where the source decides it by pointer — this pins where they part.
///
/// The source compares `spec.getDataProcessing() != dps[0]`
/// (`MzMLHandler.cpp:5258`) over
/// `std::vector<std::shared_ptr<const DataProcessing>>`
/// (`SpectrumSettings.h:165`), which is pointer identity, so two histories
/// that are equal in every field but held in different objects count as
/// different and the second record gets a dangling `dp_sp_<s>`. This port has
/// no pointer identity in its model: it compares the text each history renders
/// into the declaration blocks, so the second record gets no
/// `dataProcessingRef` and inherits the list's `defaultDataProcessingRef`.
///
/// Measured against the executed C++ on `ibminode06`
/// (`../oracle/integ-w7/dupdp_06.sh`): on the committed `refs` fixture with
/// `dp_sp_1`'s `softwareRef` repointed at `so_dp_0`, so that `dp_sp_0` and
/// `dp_sp_1` render identically, the C++ low-memory output is byte-identical
/// to its output on the unmodified fixture and still dangles `dp_sp_1` and
/// `dp_sp_2`, while this port writes neither. Recorded as a divergence in
/// native difference 12 of `docs/TOPP_PEAK_PICKER_HI_RES_SUPPORT.md`; this
/// test fails loudly if the decision ever becomes pointer-like.
#[test]
fn a_history_equal_to_the_headers_by_content_is_not_renumbered() {
    let mut first = spectrum("scan=1", 1.0, 100.0);
    first.data_processing = history("PeakPickerHiRes", "1.0");
    // Separately constructed, equal in every field: a distinct object holding
    // the same content, which is exactly what the source's pointer comparison
    // calls different.
    let mut second = spectrum("scan=2", 2.0, 200.0);
    second.data_processing = history("PeakPickerHiRes", "1.0");
    assert_eq!(first.data_processing, second.data_processing);
    assert!(!Arc::ptr_eq(
        &first.data_processing[0],
        &second.data_processing[0]
    ));
    let mut third = spectrum("scan=3", 3.0, 300.0);
    third.data_processing = history("PeakPickerHiRes", "2.0");

    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new())
        .with_reference_policy(ReferencePolicy::SourceDangling);
    consumer.set_expected_size(3, 0).unwrap();
    consumer.consume_spectrum(&mut first).unwrap();
    consumer.consume_spectrum(&mut second).unwrap();
    consumer.consume_spectrum(&mut third).unwrap();
    let text = written(consumer);

    // The second record renders the header's own declaration blocks, so this
    // port writes it no reference at all where the source would dangle
    // `dp_sp_1`.
    assert!(
        !text.contains("dataProcessingRef=\"dp_sp_1\""),
        "{text:.900}"
    );
    // The third record genuinely differs and is renumbered, as the source does.
    assert!(
        text.contains("dataProcessingRef=\"dp_sp_2\""),
        "{text:.900}"
    );
    assert!(!text.contains("<dataProcessing id=\"dp_sp_"), "{text:.900}");
}

/// A chromatogram carries neither reference in the source
/// (`MzMLHandler.cpp:5879` writes `id`, `index` and `defaultArrayLength` and
/// nothing else), so a chromatogram whose history the header does not declare
/// is written unchanged and inherits the list's default — the information is
/// lost, silently, on both sides.
#[test]
fn the_source_policy_writes_a_chromatogram_without_any_reference() {
    let mut first = chromatogram("chrom=1", 1.0);
    first.data_processing = history("PeakPickerHiRes", "1.0");
    let mut second = chromatogram("chrom=2", 2.0);
    second.data_processing = history("PeakPickerHiRes", "2.0");

    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new())
        .with_reference_policy(ReferencePolicy::SourceDangling);
    consumer.set_expected_size(0, 2).unwrap();
    consumer.consume_chromatogram(&mut first).unwrap();
    consumer.consume_chromatogram(&mut second).unwrap();
    assert_eq!(consumer.chromatograms_written(), 2);
    let text = written(consumer);
    assert!(!text.contains("dp_sp_"), "{text:.600}");
    assert!(!text.contains("sf_sp_"), "{text:.600}");
    let loaded = mzml::read(text.as_bytes()).unwrap();
    assert_eq!(loaded.chromatograms.len(), 2);
    // Both now report the header's history, which is the first one's.
    assert_eq!(loaded.chromatograms[1].data_processing.len(), 1);
    assert_eq!(
        loaded.chromatograms[1].data_processing[0].software.version,
        "1.0"
    );
}

/// A binary array's own `dataProcessingRef` is renumbered into the source's
/// array namespace, so it dangles instead of resolving to the wrong entry.
///
/// The per-record render numbers each record's array histories from one again,
/// so leaving this writer's `dp_<1>` in place would name whatever the *first*
/// record's render declared at that position — a reference that resolves, to
/// the wrong processing. The source writes `dp_sp_<s>_bi_<m>` there
/// (`MzMLHandler.cpp:5567`, `:5597`, `:5806`), which names nothing from the
/// second record on.
#[test]
fn the_source_policy_renumbers_an_array_history_reference() {
    let smoothing = Arc::new(DataProcessing {
        actions: [ProcessingAction::Smoothing].into_iter().collect(),
        ..Default::default()
    });
    let calibration = Arc::new(DataProcessing {
        actions: [ProcessingAction::MzCalibration].into_iter().collect(),
        ..Default::default()
    });
    let with_array = |id: &str, history: Vec<Arc<DataProcessing>>| {
        let mut record = spectrum(id, 1.0, 100.0);
        record.float_data_arrays.push(DataArray {
            name: "signal to noise".into(),
            data: vec![1.0, 2.0],
            metadata: Default::default(),
            data_processing: history,
        });
        record
    };

    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new())
        .with_reference_policy(ReferencePolicy::SourceDangling);
    consumer.set_expected_size(2, 0).unwrap();
    consumer
        .consume_spectrum(&mut with_array("scan=1", vec![Arc::clone(&smoothing)]))
        .unwrap();
    consumer
        .consume_spectrum(&mut with_array("scan=2", vec![Arc::clone(&calibration)]))
        .unwrap();
    let text = written(consumer);

    // The first record's array reference is the one the header declares.
    let declared = "dp_00000000000000000001";
    assert!(
        text.contains(&format!("<dataProcessing id=\"{declared}\"")),
        "{text:.900}"
    );
    assert_eq!(
        text.matches(&format!("dataProcessingRef=\"{declared}\""))
            .count(),
        1,
        "only the first record may point at the header's array history\n{text:.900}"
    );
    // The second record's names nothing.
    assert!(
        text.contains("dataProcessingRef=\"dp_sp_1_bi_0\""),
        "{text:.900}"
    );
    assert!(!text.contains("<dataProcessing id=\"dp_sp_1_bi_0\""));
}

/// The policy changes nothing when every record fits the header, which is the
/// ordinary case: the bytes are the same under either policy.
#[test]
fn the_source_policy_is_inert_when_every_record_fits_the_header() {
    let write = |policy| {
        let mut consumer =
            PlainMSDataWritingConsumer::plain(Vec::new()).with_reference_policy(policy);
        consumer.set_expected_size(2, 1).unwrap();
        consumer
            .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
            .unwrap();
        consumer
            .consume_spectrum(&mut spectrum("scan=2", 2.0, 200.0))
            .unwrap();
        consumer
            .consume_chromatogram(&mut chromatogram("chrom=1", 1.0))
            .unwrap();
        written(consumer)
    };
    assert_eq!(
        write(ReferencePolicy::Checked),
        write(ReferencePolicy::SourceDangling)
    );
}

#[test]
fn an_auxiliary_array_history_is_part_of_the_declarations() {
    let processing = Arc::new(DataProcessing {
        actions: [ProcessingAction::Smoothing].into_iter().collect(),
        ..Default::default()
    });
    let with_array = |id: &str, history: Vec<Arc<DataProcessing>>| {
        let mut record = spectrum(id, 1.0, 100.0);
        record.float_data_arrays.push(DataArray {
            name: "signal to noise".into(),
            data: vec![1.0, 2.0],
            metadata: Default::default(),
            data_processing: history,
        });
        record
    };

    // A later record whose auxiliary array carries its own data processing
    // needs a dataProcessingList the header does not have.
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(2, 0).unwrap();
    consumer
        .consume_spectrum(&mut with_array("scan=1", Vec::new()))
        .unwrap();
    let error = consumer
        .consume_spectrum(&mut with_array("scan=2", vec![Arc::clone(&processing)]))
        .unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error}");

    // The same array history on every record is fine.
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(2, 0).unwrap();
    consumer
        .consume_spectrum(&mut with_array("scan=1", vec![Arc::clone(&processing)]))
        .unwrap();
    consumer
        .consume_spectrum(&mut with_array("scan=2", vec![Arc::clone(&processing)]))
        .unwrap();
    let loaded = mzml::read(written(consumer).as_bytes()).unwrap();
    assert_eq!(loaded.spectra.len(), 2);
    for record in &loaded.spectra {
        assert_eq!(record.float_data_arrays.len(), 1);
        assert_eq!(record.float_data_arrays[0].data_processing.len(), 1);
    }
}

#[test]
fn records_may_differ_in_everything_the_header_does_not_index() {
    // fileContent legitimately differs between an MS1 and an MS2 record and is
    // deliberately not part of the declaration comparison; only the lists a
    // record's references are numbered against have to match.
    let mut first = spectrum("scan=1", 1.0, 100.0);
    first.spectrum_type = openms::kernel::SpectrumType::Profile;
    let mut second = spectrum("scan=2", 2.0, 200.0);
    second.ms_level = 2;
    second.spectrum_type = openms::kernel::SpectrumType::Centroid;
    second.precursors = vec![openms::kernel::Precursor {
        mz: 500.0,
        charge: 2,
        ..Default::default()
    }];
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_expected_size(2, 0).unwrap();
    consumer.consume_spectrum(&mut first).unwrap();
    consumer.consume_spectrum(&mut second).unwrap();
    let loaded = mzml::read(written(consumer).as_bytes()).unwrap();
    assert_eq!(loaded.spectra.len(), 2);
    assert_eq!(loaded.spectra[0].ms_level, 1);
    assert_eq!(loaded.spectra[1].ms_level, 2);
    assert_eq!(loaded.spectra[1].precursors.len(), 1);
}

#[test]
fn a_record_source_file_is_numbered_after_the_settings_source_files() {
    // The header's sourceFileList begins with the experimental settings' own
    // source files, so a record's sourceFileRef is offset by their count.
    // Rendering a later record without the settings would number it from zero
    // and silently point it at the run's source file instead of its own.
    let settings = ExperimentalSettings {
        source_files: vec![
            SourceFile {
                name: "run-a.RAW".into(),
                path: "file://.".into(),
                ..Default::default()
            },
            SourceFile {
                name: "run-b.RAW".into(),
                path: "file://.".into(),
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    let record_source = SourceFile {
        name: "record.RAW".into(),
        path: "file://.".into(),
        ..Default::default()
    };
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    consumer.set_experimental_settings(&settings).unwrap();
    consumer.set_expected_size(2, 0).unwrap();
    for (index, id) in ["scan=1", "scan=2"].iter().enumerate() {
        let mut record = spectrum(id, 1.0 + index as f64, 100.0);
        record.source_file = record_source.clone();
        consumer.consume_spectrum(&mut record).unwrap();
    }
    let text = written(consumer);
    assert_eq!(
        text.matches("sourceFileRef=\"sf_00000000000000000002\"")
            .count(),
        2,
        "{text}"
    );
    let loaded = mzml::read(text.as_bytes()).unwrap();
    // The reader returns the whole sourceFileList as run-level source files,
    // record entries included, so three come back from two settings entries
    // plus the shared record entry. What matters here is that every record
    // resolves to its own source file and not to a run-level one.
    let names: Vec<&str> = loaded
        .settings
        .source_files
        .iter()
        .map(|file| file.name.as_str())
        .collect();
    assert_eq!(names, ["run-a.RAW", "run-b.RAW", "record.RAW"]);
    for record in &loaded.spectra {
        assert_eq!(record.source_file.name, "record.RAW");
    }
}

#[test]
fn the_record_count_ceilings_are_checked_before_rendering() {
    let limits = WritingLimits {
        max_spectra: 1,
        max_chromatograms: 1,
        ..WritingLimits::default()
    };
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new()).with_limits(limits);
    consumer.set_expected_size(1, 1).unwrap();
    let error = consumer.set_expected_size(2, 0).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    let error = consumer
        .consume_spectrum(&mut spectrum("scan=2", 2.0, 200.0))
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    consumer
        .consume_chromatogram(&mut chromatogram("chrom=1", 1.0))
        .unwrap();
    let error = consumer
        .consume_chromatogram(&mut chromatogram("chrom=2", 2.0))
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    consumer.finish().unwrap();
}

#[test]
fn the_record_byte_ceiling_stops_an_oversized_render() {
    let limits = WritingLimits {
        max_record_bytes: 64,
        ..WritingLimits::default()
    };
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new()).with_limits(limits);
    let error = consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error}");
    assert!(!consumer.started_writing());
    assert_eq!(consumer.finish().unwrap(), Vec::<u8>::new());
}

#[test]
fn the_native_identifier_ceiling_is_charged() {
    let limits = WritingLimits {
        max_native_id_bytes: 8,
        ..WritingLimits::default()
    };
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new())
        .with_limits(limits)
        .with_count_policy(CountPolicy::SourceInconsistent);
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    let error = consumer
        .consume_spectrum(&mut spectrum("scan=22", 2.0, 200.0))
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error}");
    consumer.finish().unwrap();
}

#[test]
fn compression_options_reach_the_streamed_records() {
    let mut consumer =
        PlainMSDataWritingConsumer::plain(Vec::new()).with_write_options(WriteOptions {
            zlib_compression: true,
        });
    assert!(consumer.write_options().zlib_compression);
    consumer.set_expected_size(1, 0).unwrap();
    consumer
        .consume_spectrum(&mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    let text = written(consumer);
    assert!(text.contains("MS:1000574"), "{text}");
    let loaded = mzml::read(text.as_bytes()).unwrap();
    assert_eq!(loaded.spectra[0].peaks.len(), 2);
}

#[test]
fn the_consumer_can_be_driven_through_the_streaming_interface() {
    // The reader's transform operation is the consumer's intended caller.
    let mut source = MSExperiment::new();
    source.spectra.push(spectrum("scan=1", 1.0, 100.0));
    source.spectra.push(spectrum("scan=2", 2.0, 200.0));
    source.chromatograms.push(chromatogram("chrom=1", 1.0));
    let mut input = Vec::new();
    mzml::write(&mut input, &source).unwrap();

    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    let report = mzml::transform_from(
        || Ok(std::io::BufReader::new(std::io::Cursor::new(input.clone()))),
        &mut consumer,
        &mzml::TransformOptions::default(),
    )
    .unwrap();
    assert!(!report.stopped);
    assert_eq!(consumer.spectra_written(), 2);
    assert_eq!(consumer.chromatograms_written(), 1);
    let text = written(consumer);
    let loaded = mzml::read(text.as_bytes()).unwrap();
    assert_eq!(loaded.spectra, source.spectra);
    assert_eq!(loaded.chromatograms, source.chromatograms);
}

#[test]
fn the_interface_methods_forward_to_the_inherent_ones() {
    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new());
    MSDataConsumer::set_expected_size(&mut consumer, 1, 0).unwrap();
    MSDataConsumer::set_experimental_settings(&mut consumer, &ExperimentalSettings::default())
        .unwrap();
    assert_eq!(consumer.expected_size(), (1, 0));
    let flow = MSDataConsumer::consume_spectrum(&mut consumer, &mut spectrum("scan=1", 1.0, 100.0))
        .unwrap();
    assert_eq!(flow, std::ops::ControlFlow::Continue(()));
    assert_eq!(consumer.spectra_written(), 1);
    consumer.finish().unwrap();
}

// ---------------------------------------------------------------------------
// NoopMSDataWritingConsumer
// ---------------------------------------------------------------------------
#[test]
fn the_noop_consumer_accepts_everything_and_touches_no_file() {
    let mut consumer = NoopMSDataWritingConsumer::new();
    assert_eq!(consumer, NoopMSDataWritingConsumer::default());
    consumer.set_expected_size(9, 9).unwrap();
    consumer
        .set_experimental_settings(&ExperimentalSettings::default())
        .unwrap();
    let mut record = spectrum("scan=1", 1.0, 100.0);
    assert_eq!(
        consumer.consume_spectrum(&mut record).unwrap(),
        std::ops::ControlFlow::Continue(())
    );
    let mut chrom = chromatogram("chrom=1", 1.0);
    assert_eq!(
        consumer.consume_chromatogram(&mut chrom).unwrap(),
        std::ops::ControlFlow::Continue(())
    );
    assert_eq!(consumer.spectra_written(), 1);
    assert_eq!(consumer.chromatograms_written(), 1);
    // Records are untouched, and a duplicate identifier is accepted because
    // nothing is written.
    assert_eq!(record, spectrum("scan=1", 1.0, 100.0));
    assert_eq!(
        consumer.consume_spectrum(&mut record).unwrap(),
        std::ops::ControlFlow::Continue(())
    );
    assert_eq!(consumer.spectra_written(), 2);
}

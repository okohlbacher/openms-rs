// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ported from `IndexedMzMLFile_test.cpp`, the class test of
//! `Internal::IndexedMzMLHandler`, plus native bounded-work coverage.
//!
//! Every literal taken from the upstream test or from the unmodified fixture is
//! transcribed source review (tier 3); the synthetic documents and the
//! resource-limit expectations are independently derived (tier 4).

#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::indexed_mzml_handler::{IndexedMzMLHandler, RecordKind, RecordReadLimits};
use openms::format::mzml;
use openms::format::peak_options::PeakFileOptions;
use openms::kernel::NumericRange;
use openms::system::file::TempDir;
use std::path::PathBuf;

const SCAN_1: &str = "controllerType=0 controllerNumber=1 scan=1";
const SCAN_2: &str = "controllerType=0 controllerNumber=1 scan=2";

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// The unmodified upstream `IndexedmzMLFile_1.mzML`, as the class test uses it.
fn upstream() -> PathBuf {
    data("indexed_mzml/IndexedmzMLFile_1.mzML")
}

fn open_upstream() -> IndexedMzMLHandler {
    IndexedMzMLHandler::open(upstream()).expect("upstream indexed fixture opens")
}

fn temp(bytes: &[u8]) -> (TempDir, PathBuf) {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("input.mzML");
    std::fs::write(&path, bytes).unwrap();
    (dir, path)
}

// ---------------------------------------------------------------------------
// Synthetic indexed documents, used for the offset-ordering branch that the
// upstream fixture cannot reach and for hostile index entries.
// ---------------------------------------------------------------------------

const SPECTRUM: &str = concat!(
    "<spectrum id=\"scan=1\" index=\"0\" defaultArrayLength=\"2\">",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"1\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000128\" name=\"profile spectrum\"/>",
    "<scanList count=\"1\">",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\"/>",
    "<scan><cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"12.5\"",
    " unitAccession=\"UO:0000010\" unitName=\"second\" unitCvRef=\"UO\"/></scan>",
    "</scanList>",
    "<binaryDataArrayList count=\"2\">",
    "<binaryDataArray encodedLength=\"24\">",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\"",
    " unitAccession=\"MS:1000040\" unitName=\"m/z\" unitCvRef=\"MS\"/>",
    "<binary>AAAAAAAAWUAAAAAAAABpQA==</binary>",
    "</binaryDataArray>",
    "<binaryDataArray encodedLength=\"24\">",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\"",
    " unitAccession=\"MS:1000131\" unitName=\"number of detector counts\" unitCvRef=\"MS\"/>",
    "<binary>AAAAAAAAJEAAAAAAAAA0QA==</binary>",
    "</binaryDataArray>",
    "</binaryDataArrayList>",
    "</spectrum>",
);

const CHROMATOGRAM: &str = concat!(
    "<chromatogram id=\"TIC\" index=\"0\" defaultArrayLength=\"2\">",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000235\" name=\"total ion current chromatogram\"/>",
    "<binaryDataArrayList count=\"2\">",
    "<binaryDataArray encodedLength=\"24\">",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000595\" name=\"time array\"",
    " unitAccession=\"UO:0000010\" unitName=\"second\" unitCvRef=\"UO\"/>",
    "<binary>AAAAAAAA8D8AAAAAAAAAQA==</binary>",
    "</binaryDataArray>",
    "<binaryDataArray encodedLength=\"24\">",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\"",
    " unitAccession=\"MS:1000131\" unitName=\"number of detector counts\" unitCvRef=\"MS\"/>",
    "<binary>AAAAAAAAJEAAAAAAAAA0QA==</binary>",
    "</binaryDataArray>",
    "</binaryDataArrayList>",
    "</chromatogram>",
);

const PROLOGUE: &str = concat!(
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",
    "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\">\n",
    "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">\n",
    "<cvList count=\"1\"><cv id=\"MS\" fullName=\"PSI-MS\" version=\"4.1\"",
    " URI=\"https://example.invalid/psi-ms.obo\"/></cvList>\n",
    "<fileDescription><fileContent/></fileDescription>\n",
    "<softwareList count=\"1\"><software id=\"sw\" version=\"1\"/></softwareList>\n",
    "<instrumentConfigurationList count=\"1\"><instrumentConfiguration id=\"ic\">",
    "</instrumentConfiguration></instrumentConfigurationList>\n",
    "<dataProcessingList count=\"1\"><dataProcessing id=\"dp\">",
    "<processingMethod order=\"1\" softwareRef=\"sw\"/></dataProcessing></dataProcessingList>\n",
    "<run id=\"run0\" defaultInstrumentConfigurationRef=\"ic\">\n",
);

fn spectrum_xml(id: &str) -> String {
    SPECTRUM.replace("id=\"scan=1\"", &format!("id=\"{id}\""))
}

/// The mzML half of a synthetic indexed document: `spectra` spectra named
/// `scan=1..`, one `TIC` chromatogram, in the requested list order.
fn body_with(spectra: usize, chromatograms_first: bool) -> String {
    let records: String = (1..=spectra)
        .map(|n| format!("{}\n", spectrum_xml(&format!("scan={n}"))))
        .collect();
    let spectrum_list = format!(
        "<spectrumList count=\"{spectra}\" defaultDataProcessingRef=\"dp\">\n{records}</spectrumList>\n"
    );
    let chromatogram_list = format!(
        "<chromatogramList count=\"1\" defaultDataProcessingRef=\"dp\">\n{CHROMATOGRAM}\n</chromatogramList>\n"
    );
    let lists = if chromatograms_first {
        format!("{chromatogram_list}{spectrum_list}")
    } else {
        format!("{spectrum_list}{chromatogram_list}")
    };
    format!("{PROLOGUE}{lists}</run>\n</mzML>\n")
}

fn body(chromatograms_first: bool) -> String {
    body_with(1, chromatograms_first)
}

/// Byte offsets of every `<spectrum` and `<chromatogram` start tag in `body`.
fn true_offsets(body: &str) -> (Vec<u64>, Vec<u64>) {
    let find = |needle: &str| {
        body.match_indices(needle)
            .map(|(at, _)| at as u64)
            .collect::<Vec<_>>()
    };
    (find("<spectrum id="), find("<chromatogram id="))
}

/// Append a footer index naming the given entries. `indexListOffset` always
/// points at the `<indexList>` element that follows the mzML body.
fn assemble(body: &str, spectra: &[(&str, u64)], chromatograms: &[(&str, u64)]) -> Vec<u8> {
    let section = |name: &str, entries: &[(&str, u64)]| {
        let offsets: String = entries
            .iter()
            .map(|(id, at)| format!("<offset idRef=\"{id}\">{at}</offset>"))
            .collect();
        format!("<index name=\"{name}\">{offsets}</index>")
    };
    let footer = format!(
        "<indexList count=\"2\">{}{}</indexList><indexListOffset>{}</indexListOffset>\
         <fileChecksum>0</fileChecksum>\n</indexedmzML>\n",
        section("spectrum", spectra),
        section("chromatogram", chromatograms),
        body.len(),
    );
    format!("{body}{footer}").into_bytes()
}

/// The one-spectrum document with a truthful index.
fn synthetic(chromatograms_first: bool) -> Vec<u8> {
    let body = body(chromatograms_first);
    let (spectra, chromatograms) = true_offsets(&body);
    assemble(
        &body,
        &[("scan=1", spectra[0])],
        &[("TIC", chromatograms[0])],
    )
}

// ---------------------------------------------------------------------------
// START_SECTION((IndexedMzMLHandler(std::string filename)))
// ---------------------------------------------------------------------------

#[test]
fn constructor_from_filename_parses_the_index() {
    // The upstream section only checks that `new` returned non-null; the Rust
    // equivalent of "the object exists" is that `open` produced `Ok`, so the
    // index contents are asserted here instead.
    let handler = open_upstream();
    assert_eq!(handler.path(), upstream().as_path());
    assert_eq!(handler.index_list_offset(), 667742);
    assert_eq!(handler.offset(RecordKind::Spectrum, 0), Some(24146));
    assert_eq!(handler.offset(RecordKind::Spectrum, 1), Some(345745));
    assert_eq!(handler.offset(RecordKind::Chromatogram, 0), Some(665563));
    assert_eq!(handler.native_id(RecordKind::Spectrum, 0), Some(SCAN_1));
    assert_eq!(handler.native_id(RecordKind::Chromatogram, 0), Some("TIC"));
    assert!(handler.spectra_before_chromatograms());
}

// ---------------------------------------------------------------------------
// START_SECTION((~IndexedMzMLHandler()))
// ---------------------------------------------------------------------------

#[test]
fn dropping_a_handler_releases_the_file() {
    // The source destructor is `= default`; the owned `std::ifstream` closes.
    // The Rust equivalent is `File`'s drop, observable because the temporary
    // directory can be removed and reopened afterwards.
    let (dir, path) = temp(&synthetic(false));
    let handler = IndexedMzMLHandler::open(&path).unwrap();
    assert_eq!(handler.spectrum_count(), 1);
    drop(handler);
    let reopened = IndexedMzMLHandler::open(&path).unwrap();
    assert_eq!(reopened.spectrum_count(), 1);
    drop(reopened);
    drop(dir);
    assert!(!path.exists());
}

// ---------------------------------------------------------------------------
// START_SECTION((IndexedMzMLHandler()))
// ---------------------------------------------------------------------------

#[test]
fn there_is_no_unopened_handler() {
    // The source default constructor leaves `parsing_success_` false and its own
    // documentation calls any retrieval on that object invalid. This port has no
    // such state: the only constructor parses the index and reports failure.
    let error =
        IndexedMzMLHandler::open(data("indexed_mzml_handler/fileDoesNotExist")).unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error:?}");
}

// ---------------------------------------------------------------------------
// START_SECTION((IndexedMzMLHandler(const IndexedMzMLHandler &source)))
// ---------------------------------------------------------------------------

#[test]
fn a_second_handler_on_one_file_reads_the_same_data() {
    // The source copy constructor deliberately reopens the file rather than
    // copying the stream. It is not ported: a handler owns a seek position, so
    // an independent reader is an independent `open`. The upstream assertions
    // are reproduced across two such handlers.
    let mut first = open_upstream();
    let mut second = open_upstream();

    assert_eq!(first.spectrum_count(), second.spectrum_count());
    assert_eq!(first.chromatogram_count(), second.chromatogram_count());
    assert_eq!(first.spectrum_count(), 2);
    assert_eq!(first.chromatogram_count(), 1);

    for index in 0..2 {
        let a = first.spectrum(index).unwrap().unwrap();
        let b = second.spectrum(index).unwrap().unwrap();
        let mz = |s: &openms::kernel::MSSpectrum| s.peaks.iter().map(|p| p.mz).collect::<Vec<_>>();
        let intensity = |s: &openms::kernel::MSSpectrum| {
            s.peaks.iter().map(|p| p.intensity).collect::<Vec<_>>()
        };
        assert_eq!(mz(&a), mz(&b));
        assert_eq!(intensity(&a), intensity(&b));
    }

    let a = first.chromatogram(0).unwrap().unwrap();
    let b = second.chromatogram(0).unwrap().unwrap();
    assert_eq!(
        a.peaks.iter().map(|p| p.rt).collect::<Vec<_>>(),
        b.peaks.iter().map(|p| p.rt).collect::<Vec<_>>()
    );
    assert_eq!(
        a.peaks.iter().map(|p| p.intensity).collect::<Vec<_>>(),
        b.peaks.iter().map(|p| p.intensity).collect::<Vec<_>>()
    );

    // The source copy constructor omits both native-id maps from its initialiser
    // list, so every `getMSSpectrumByNativeId` on a copy throws. Here the second
    // handler resolves identifiers exactly like the first.
    assert_eq!(second.index_of(RecordKind::Spectrum, SCAN_2), Some(1));
    assert_eq!(second.index_of(RecordKind::Chromatogram, "TIC"), Some(0));
    assert!(second.spectrum_by_native_id(SCAN_1).unwrap().is_some());
}

// ---------------------------------------------------------------------------
// START_SECTION((bool getParsingSuccess() const))
// ---------------------------------------------------------------------------

#[test]
fn parsing_success_is_the_result_of_open() {
    // Missing file: the source throws FileNotFound out of openFile and leaves
    // parsing_success_ false.
    let missing = IndexedMzMLHandler::open(data("indexed_mzml_handler/fileDoesNotExist"));
    assert!(matches!(missing, Err(Error::Io(_))));

    // A plain, non-indexed mzML: the source records parsing_success_ == false.
    let plain = IndexedMzMLHandler::open(data("mzml_upstream_minimal.mzML"));
    assert!(
        matches!(&plain, Err(Error::Parse { message, .. }) if message.contains("indexListOffset")),
        "{plain:?}"
    );

    // The indexed fixture: parsing_success_ == true.
    let handler = open_upstream();
    assert_eq!(handler.spectrum_count(), 2);
    assert_eq!(handler.chromatogram_count(), 1);
    assert!(!handler.is_empty());
}

// ---------------------------------------------------------------------------
// START_SECTION((void openFile(std::string filename)))
// ---------------------------------------------------------------------------

#[test]
fn opening_replaces_rather_than_accumulates() {
    // `parseFooter_` never clears the offset vectors, so a second `openFile` on
    // one source object appends: `getNrSpectra` reports four after opening the
    // two-spectrum fixture twice. `open` here always yields a fresh handler.
    let first = open_upstream();
    assert_eq!(first.spectrum_count(), 2);
    let second = open_upstream();
    assert_eq!(second.spectrum_count(), 2);
    assert_eq!(second.chromatogram_count(), 1);
    assert_eq!(second.offset(RecordKind::Spectrum, 2), None);
    assert!(IndexedMzMLHandler::open(data("mzml_upstream_minimal.mzML")).is_err());
}

// ---------------------------------------------------------------------------
// START_SECTION((size_t getNrSpectra() const))
// ---------------------------------------------------------------------------

#[test]
fn spectrum_count_matches_the_index() {
    let handler = open_upstream();
    assert_eq!(handler.spectrum_count(), 2);
    assert_eq!(handler.count(RecordKind::Spectrum), 2);
}

// ---------------------------------------------------------------------------
// START_SECTION((size_t getNrChromatograms() const))
// ---------------------------------------------------------------------------

#[test]
fn chromatogram_count_matches_the_index() {
    let handler = open_upstream();
    assert_eq!(handler.chromatogram_count(), 1);
    assert_eq!(handler.count(RecordKind::Chromatogram), 1);
}

// ---------------------------------------------------------------------------
// START_SECTION((OpenMS::Interfaces::SpectrumPtr getSpectrumById(int id)))
// ---------------------------------------------------------------------------

#[test]
fn spectrum_arrays_match_a_whole_file_load() {
    // The upstream section compares the handler against `MzMLFile().load` of the
    // same file. The two `Interfaces::Spectrum` arrays become the single peak
    // vector here, so both upstream size assertions land on `peaks.len()`.
    let experiment = mzml::load(upstream()).unwrap();
    let mut handler = open_upstream();
    assert_eq!(handler.spectrum_count(), experiment.spectra.len());

    let spectrum = handler.spectrum(0).unwrap().unwrap();
    assert_eq!(spectrum.peaks.len(), experiment.spectra[0].peaks.len());
    assert_eq!(spectrum.peaks.len(), 19914);
    assert_eq!(
        spectrum.peaks.iter().map(|p| p.mz).collect::<Vec<_>>(),
        experiment.spectra[0]
            .peaks
            .iter()
            .map(|p| p.mz)
            .collect::<Vec<_>>()
    );

    // Source `getSpectrumById(-1)` throws IllegalArgument. An index is `usize`
    // here, so the negative case cannot be expressed; the upper bound remains.
    let error = handler.spectrum(handler.spectrum_count() + 1).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert!(handler.spectrum(2).is_err());

    // Source: retrieval on an object whose parsing failed throws ParseError.
    // Unreachable here, because such an object cannot be constructed.
    assert!(IndexedMzMLHandler::open(data("mzml_upstream_minimal.mzML")).is_err());
}

// ---------------------------------------------------------------------------
// START_SECTION((OpenMS::MSSpectrum getMSSpectrumById(int id)))
// ---------------------------------------------------------------------------

#[test]
fn spectrum_by_index_carries_full_record_metadata() {
    let experiment = mzml::load(upstream()).unwrap();
    let mut handler = open_upstream();
    assert_eq!(handler.spectrum_count(), experiment.spectra.len());

    let spectrum = handler.spectrum(0).unwrap().unwrap();
    assert_eq!(spectrum.peaks.len(), experiment.spectra[0].peaks.len());
    // The upstream test comments out the native-id assertion; it holds here.
    assert_eq!(spectrum.native_id, experiment.spectra[0].native_id);
    assert_eq!(spectrum.native_id, SCAN_1);
    // Source `domParseString_` reads only `id`, `defaultArrayLength` and the
    // binaryDataArray elements, so MS level and RT come back default-constructed.
    // Routing the record through the whole-file reader keeps them.
    assert_eq!(spectrum.ms_level, 1);
    assert!((spectrum.rt - 0.2961).abs() < 1e-9, "{}", spectrum.rt);
    assert_eq!(spectrum.rt, experiment.spectra[0].rt);

    let second = handler.spectrum(1).unwrap().unwrap();
    assert_eq!(second.native_id, SCAN_2);
    assert_eq!(second.peaks.len(), 19800);
    assert!((second.rt - 0.4738).abs() < 1e-9, "{}", second.rt);

    let error = handler.spectrum(handler.spectrum_count() + 1).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
}

// ---------------------------------------------------------------------------
// START_SECTION((void getMSSpectrumByNativeId(std::string id, MSSpectrum& s)))
// ---------------------------------------------------------------------------

#[test]
fn spectrum_by_native_id_resolves_through_the_index() {
    let experiment = mzml::load(upstream()).unwrap();
    let mut handler = open_upstream();
    assert_eq!(handler.spectrum_count(), experiment.spectra.len());

    let spectrum = handler.spectrum_by_native_id(SCAN_1).unwrap().unwrap();
    assert_eq!(spectrum.peaks.len(), experiment.spectra[0].peaks.len());
    assert_eq!(spectrum.native_id, experiment.spectra[0].native_id);

    // Source: an unknown native id throws IllegalArgument.
    let error = handler.spectrum_by_native_id("TEST").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert_eq!(handler.index_of(RecordKind::Spectrum, "TEST"), None);
    assert_eq!(handler.index_of(RecordKind::Spectrum, SCAN_2), Some(1));
    // The same exception on a handler whose parsing failed: not constructible.
    assert!(IndexedMzMLHandler::open(data("mzml_upstream_minimal.mzML")).is_err());
}

// ---------------------------------------------------------------------------
// START_SECTION((OpenMS::Interfaces::ChromatogramPtr getChromatogramById(int)))
// ---------------------------------------------------------------------------

#[test]
fn chromatogram_arrays_match_a_whole_file_load() {
    let experiment = mzml::load(upstream()).unwrap();
    let mut handler = open_upstream();
    assert_eq!(handler.chromatogram_count(), experiment.chromatograms.len());

    let chromatogram = handler.chromatogram(0).unwrap().unwrap();
    assert_eq!(
        chromatogram.peaks.len(),
        experiment.chromatograms[0].peaks.len()
    );
    assert_eq!(chromatogram.peaks.len(), 48);
}

// ---------------------------------------------------------------------------
// START_SECTION((OpenMS::MSChromatogram getMSChromatogramById(int id)))
// ---------------------------------------------------------------------------

#[test]
fn chromatogram_by_index_carries_its_native_id() {
    let experiment = mzml::load(upstream()).unwrap();
    let mut handler = open_upstream();
    assert_eq!(handler.chromatogram_count(), experiment.chromatograms.len());

    let chromatogram = handler.chromatogram(0).unwrap().unwrap();
    assert_eq!(
        chromatogram.peaks.len(),
        experiment.chromatograms[0].peaks.len()
    );
    assert_eq!(
        chromatogram.native_id,
        experiment.chromatograms[0].native_id
    );
    assert_eq!(chromatogram.native_id, "TIC");
    assert_eq!(
        chromatogram.peaks.iter().map(|p| p.rt).collect::<Vec<_>>(),
        experiment.chromatograms[0]
            .peaks
            .iter()
            .map(|p| p.rt)
            .collect::<Vec<_>>()
    );
}

// ---------------------------------------------------------------------------
// START_SECTION((void getMSChromatogramByNativeId(std::string, MSChromatogram&)))
// ---------------------------------------------------------------------------

#[test]
fn chromatogram_by_native_id_resolves_through_the_index() {
    let experiment = mzml::load(upstream()).unwrap();
    let mut handler = open_upstream();
    assert_eq!(handler.chromatogram_count(), experiment.chromatograms.len());

    let chromatogram = handler.chromatogram_by_native_id("TIC").unwrap().unwrap();
    assert_eq!(
        chromatogram.peaks.len(),
        experiment.chromatograms[0].peaks.len()
    );
    assert_eq!(
        chromatogram.native_id,
        experiment.chromatograms[0].native_id
    );

    let error = handler.chromatogram_by_native_id("TEST").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert_eq!(handler.index_of(RecordKind::Chromatogram, "TEST"), None);
    assert!(IndexedMzMLHandler::open(data("mzml_upstream_minimal.mzML")).is_err());
}

// ---------------------------------------------------------------------------
// START_SECTION(([EXTRA] load broken file)) -- indexListOffset 2^64
// ---------------------------------------------------------------------------

#[test]
fn footer_offset_beyond_sixty_three_bits_is_rejected() {
    // `IndexedmzMLFile_2_broken.mzML` carries 18446744073709551616 in
    // indexListOffset. The source throws ConversionError because it does not fit
    // a long long; the decoder here rejects anything above 2^63-1 identically.
    let error =
        IndexedMzMLHandler::open(data("indexed_mzml_handler/IndexedmzMLFile_2_broken.mzML"))
            .unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
}

// ---------------------------------------------------------------------------
// START_SECTION(([EXTRA] load broken file)) -- indexListOffset 2^63-1
// ---------------------------------------------------------------------------

#[test]
fn footer_offset_past_end_of_file_is_rejected() {
    // `IndexedmzMLFile_3_broken.mzML` carries 9223372036854775807, which fits a
    // 64-bit streampos but is far past the 3202-byte file, so the source's
    // parsing_success_ stays false.
    let error =
        IndexedMzMLHandler::open(data("indexed_mzml_handler/IndexedmzMLFile_3_broken.mzML"))
            .unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
}

// ---------------------------------------------------------------------------
// Native coverage beyond the class test.
// ---------------------------------------------------------------------------

#[test]
fn the_last_record_is_trimmed_at_its_own_closing_tag() {
    // Spectrum 1 is the last spectrum, so the source reads [345745, 665563) and
    // hands `</spectrum>\n\t\t</spectrumList>\n\t\t<chromatogramList ...>` to a
    // DOM parser that was constructed without an error handler and therefore
    // discards the resulting well-formedness errors. The chromatogram is the last
    // record of the file, so its range runs to <indexList> and carries
    // `</chromatogramList></run></mzML>`.
    let mut handler = open_upstream();

    let spectrum = handler.record_xml(RecordKind::Spectrum, 1).unwrap();
    assert!(spectrum.starts_with(b"<spectrum id=\"controllerType=0"));
    assert!(spectrum.ends_with(b"</spectrum>"));
    assert!(spectrum.len() < (665563 - 345745));
    assert!(!spectrum.windows(15).any(|w| w == b"</spectrumList>"));

    let chromatogram = handler.record_xml(RecordKind::Chromatogram, 0).unwrap();
    assert!(chromatogram.starts_with(b"<chromatogram id=\"TIC\""));
    assert!(chromatogram.ends_with(b"</chromatogram>"));
    assert!(chromatogram.len() < (667742 - 665563));
    assert!(!chromatogram.windows(6).any(|w| w == b"</run>"));

    // A record that is not the last one ends at the next offset, and trimming
    // removes only the inter-record whitespace.
    let first = handler.record_xml(RecordKind::Spectrum, 0).unwrap();
    assert!(first.ends_with(b"</spectrum>"));
    assert!(first.len() < (345745 - 24146));
}

#[test]
fn chromatograms_before_spectra_bound_the_last_record_correctly() {
    // `spectra_before_chroms_` is false only when both sections are populated and
    // the first chromatogram precedes the first spectrum. The upstream fixture
    // cannot reach that branch, so this uses a synthetic document.
    let (_dir, path) = temp(&synthetic(true));
    let mut handler = IndexedMzMLHandler::open(&path).unwrap();
    assert!(!handler.spectra_before_chromatograms());
    assert_eq!(handler.spectrum_count(), 1);
    assert_eq!(handler.chromatogram_count(), 1);

    // The last chromatogram now ends at the first spectrum, and the last spectrum
    // ends at <indexList>.
    let chromatogram = handler.record_xml(RecordKind::Chromatogram, 0).unwrap();
    assert!(chromatogram.ends_with(b"</chromatogram>"));
    let spectrum = handler.record_xml(RecordKind::Spectrum, 0).unwrap();
    assert!(spectrum.ends_with(b"</spectrum>"));

    let decoded = handler.spectrum(0).unwrap().unwrap();
    assert_eq!(decoded.native_id, "scan=1");
    assert_eq!(decoded.peaks.len(), 2);
    assert_eq!(decoded.peaks[0].mz, 100.0);
    assert_eq!(decoded.peaks[1].mz, 200.0);
    assert_eq!(decoded.rt, 12.5);
    let decoded = handler.chromatogram(0).unwrap().unwrap();
    assert_eq!(decoded.native_id, "TIC");
    assert_eq!(decoded.peaks.len(), 2);
    assert_eq!(decoded.peaks[0].rt, 1.0);
    assert_eq!(decoded.peaks[1].intensity, 20.0);

    // The same document with spectra first decodes identically.
    let (_dir, path) = temp(&synthetic(false));
    let mut handler = IndexedMzMLHandler::open(&path).unwrap();
    assert!(handler.spectra_before_chromatograms());
    assert_eq!(handler.spectrum(0).unwrap().unwrap().peaks.len(), 2);
    assert_eq!(handler.chromatogram(0).unwrap().unwrap().peaks.len(), 2);
}

#[test]
fn the_synthetic_document_is_a_valid_whole_file_mzml() {
    // Guards the synthetic fixtures above: the handler must not be the only
    // reader that accepts them.
    let (_dir, path) = temp(&synthetic(false));
    let experiment = mzml::load(&path).unwrap();
    assert_eq!(experiment.spectra.len(), 1);
    assert_eq!(experiment.chromatograms.len(), 1);
    assert_eq!(experiment.spectra[0].peaks.len(), 2);
    assert_eq!(experiment.chromatograms[0].peaks.len(), 2);
}

#[test]
fn decreasing_index_offsets_are_rejected_before_allocating() {
    // The source computes `endidx - startidx` as a signed streampos and passes it
    // to `new char[readl + 1]`; a decreasing pair makes that length negative.
    let body = body_with(2, false);
    let (spectra, chromatograms) = true_offsets(&body);
    assert!(spectra[0] < spectra[1]);

    let truthful = assemble(
        &body,
        &[("scan=1", spectra[0]), ("scan=2", spectra[1])],
        &[("TIC", chromatograms[0])],
    );
    let (_dir, path) = temp(&truthful);
    let mut handler = IndexedMzMLHandler::open(&path).unwrap();
    assert_eq!(handler.spectrum_count(), 2);
    assert!(handler.spectrum(0).unwrap().is_some());

    // The same file with the two spectrum offsets exchanged: the range of the
    // first entry now runs backwards.
    let reversed = assemble(
        &body,
        &[("scan=2", spectra[1]), ("scan=1", spectra[0])],
        &[("TIC", chromatograms[0])],
    );
    let (_dir, path) = temp(&reversed);
    let mut handler = IndexedMzMLHandler::open(&path).unwrap();
    let error = handler.spectrum(0).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("do not increase")),
        "{error:?}"
    );
    assert!(handler.record_xml(RecordKind::Spectrum, 0).is_err());
}

#[test]
fn offsets_past_the_end_of_the_file_are_rejected() {
    let body = body(false);
    let (_, chromatograms) = true_offsets(&body);
    let bytes = assemble(
        &body,
        &[("scan=1", u64::from(u32::MAX))],
        &[("TIC", chromatograms[0])],
    );
    let (_dir, path) = temp(&bytes);
    let error = IndexedMzMLHandler::open(&path).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("past the end")),
        "{error:?}"
    );
}

#[test]
fn an_offset_that_is_not_a_record_start_is_rejected() {
    // The source seeks to whatever the index claims and feeds the bytes to a
    // parser whose errors are discarded, yielding an empty record.
    let body = body(false);
    let (spectra, chromatograms) = true_offsets(&body);
    let bytes = assemble(
        &body,
        &[("scan=1", spectra[0] + 3)],
        &[("TIC", chromatograms[0])],
    );
    let (_dir, path) = temp(&bytes);
    let mut handler = IndexedMzMLHandler::open(&path).unwrap();
    let error = handler.spectrum(0).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("expected element")),
        "{error:?}"
    );
}

#[test]
fn an_oversized_record_range_is_refused_before_reading() {
    let (_dir, path) = temp(&synthetic(false));
    let limits = RecordReadLimits {
        max_record_bytes: 16,
        ..RecordReadLimits::default()
    };
    let mut handler =
        IndexedMzMLHandler::open_with_limits(&path, limits, mzml::ReadOptions::default()).unwrap();
    assert_eq!(handler.limits().max_record_bytes, 16);
    let error = handler.spectrum(0).unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(message) if message.contains("exceeds the configured limit")),
        "{error:?}"
    );
    // The failure leaves the handler usable: the index is untouched.
    assert_eq!(handler.spectrum_count(), 1);
    assert_eq!(handler.native_id(RecordKind::Spectrum, 0), Some("scan=1"));
}

#[test]
fn a_tiny_header_ceiling_is_refused_at_open() {
    let (_dir, path) = temp(&synthetic(false));
    let limits = RecordReadLimits {
        max_header_bytes: 8,
        ..RecordReadLimits::default()
    };
    let error = IndexedMzMLHandler::open_with_limits(&path, limits, mzml::ReadOptions::default())
        .unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
}

#[test]
fn a_record_count_ceiling_is_refused_at_open() {
    let (_dir, path) = temp(&synthetic(false));
    let limits = RecordReadLimits {
        max_records: 0,
        ..RecordReadLimits::default()
    };
    let error = IndexedMzMLHandler::open_with_limits(&path, limits, mzml::ReadOptions::default())
        .unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(message) if message.contains("record count")),
        "{error:?}"
    );
    let limits = RecordReadLimits {
        max_native_id_bytes: 1,
        ..RecordReadLimits::default()
    };
    let error = IndexedMzMLHandler::open_with_limits(&path, limits, mzml::ReadOptions::default())
        .unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(message) if message.contains("native identifiers")),
        "{error:?}"
    );
}

#[test]
fn a_narrow_list_tag_window_is_reported_rather_than_guessed() {
    let (_dir, path) = temp(&synthetic(false));
    let limits = RecordReadLimits {
        max_list_tag_bytes: 4,
        ..RecordReadLimits::default()
    };
    let error = IndexedMzMLHandler::open_with_limits(&path, limits, mzml::ReadOptions::default())
        .unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("search window")),
        "{error:?}"
    );
}

#[test]
fn a_mismatched_index_identifier_is_detected_after_decoding() {
    // The index claims the spectrum offset carries the chromatogram's id. The
    // source never compares the two.
    let body = body(false);
    let (spectra, chromatograms) = true_offsets(&body);
    let bytes = assemble(
        &body,
        &[("wrong", spectra[0])],
        &[("TIC", chromatograms[0])],
    );
    let (_dir, path) = temp(&bytes);
    let mut handler = IndexedMzMLHandler::open(&path).unwrap();
    assert_eq!(handler.native_id(RecordKind::Spectrum, 0), Some("wrong"));
    let error = handler.spectrum(0).unwrap_err();
    assert!(
        matches!(&error, Error::Parse { message, .. } if message.contains("does not match")),
        "{error:?}"
    );
}

#[test]
fn peak_file_options_filter_whole_records_and_individual_peaks() {
    let (_dir, path) = temp(&synthetic(false));
    let mut handler = IndexedMzMLHandler::open(&path).unwrap();

    // MS level filter excludes the record entirely: no peak decoding happens.
    let mut options = PeakFileOptions::new();
    options.set_ms_levels(&[2]).unwrap();
    handler.set_options(options);
    assert!(handler.spectrum(0).unwrap().is_none());
    assert!(handler.options().has_ms_levels());

    // RT range filter, checked against the record's own scan start time.
    handler.set_options(PeakFileOptions::new());
    handler
        .options_mut()
        .set_rt_range(NumericRange { min: 0.0, max: 1.0 });
    assert!(handler.spectrum(0).unwrap().is_none());
    handler.options_mut().set_rt_range(NumericRange {
        min: 0.0,
        max: 100.0,
    });
    assert_eq!(handler.spectrum(0).unwrap().unwrap().peaks.len(), 2);

    // m/z range selects peaks and keeps the record.
    handler.set_options(PeakFileOptions::new());
    handler.options_mut().set_mz_range(NumericRange {
        min: 150.0,
        max: 500.0,
    });
    let spectrum = handler.spectrum(0).unwrap().unwrap();
    assert_eq!(spectrum.peaks.len(), 1);
    assert_eq!(spectrum.peaks[0].mz, 200.0);

    // Intensity range likewise.
    handler.set_options(PeakFileOptions::new());
    handler.options_mut().set_intensity_range(NumericRange {
        min: 15.0,
        max: 100.0,
    });
    let spectrum = handler.spectrum(0).unwrap().unwrap();
    assert_eq!(spectrum.peaks.len(), 1);
    assert_eq!(spectrum.peaks[0].intensity, 20.0);

    // Chromatograms can be skipped wholesale.
    handler.set_options(PeakFileOptions::new());
    handler.options_mut().skip_chromatograms = true;
    assert!(handler.chromatogram(0).unwrap().is_none());

    // metadata_only decodes no records at all.
    handler.set_options(PeakFileOptions::new());
    handler.options_mut().metadata_only = true;
    assert!(handler.spectrum(0).unwrap().is_none());
    assert!(handler.chromatogram(0).unwrap().is_none());

    // fill_data keeps the record and drops its peaks.
    handler.set_options(PeakFileOptions::new());
    handler.options_mut().fill_data = false;
    let spectrum = handler.spectrum(0).unwrap().unwrap();
    assert_eq!(spectrum.native_id, "scan=1");
    assert!(spectrum.peaks.is_empty());
}

#[test]
fn a_repeated_native_id_resolves_to_its_first_index() {
    // Source `parseFooter_` uses `unordered_map::emplace`, which does not
    // overwrite; the ordered map here keeps the first entry for the same reason.
    let body = body(false);
    let (spectra, chromatograms) = true_offsets(&body);
    let bytes = assemble(
        &body,
        &[("dup", spectra[0]), ("dup", spectra[0])],
        &[("TIC", chromatograms[0])],
    );
    let (_dir, path) = temp(&bytes);
    let handler = IndexedMzMLHandler::open(&path).unwrap();
    assert_eq!(handler.spectrum_count(), 2);
    assert_eq!(handler.index_of(RecordKind::Spectrum, "dup"), Some(0));
    assert_eq!(handler.native_id(RecordKind::Spectrum, 1), Some("dup"));
}

#[test]
fn record_kind_names_the_mzml_elements() {
    assert_eq!(RecordKind::Spectrum.element(), "spectrum");
    assert_eq!(RecordKind::Spectrum.list_element(), "spectrumList");
    assert_eq!(RecordKind::Chromatogram.element(), "chromatogram");
    assert_eq!(RecordKind::Chromatogram.list_element(), "chromatogramList");
    assert_ne!(RecordKind::Spectrum, RecordKind::Chromatogram);
}

#[test]
fn default_limits_are_explicit_and_reported() {
    let limits = RecordReadLimits::default();
    assert_eq!(limits.max_record_bytes, 256 << 20);
    assert_eq!(limits.max_header_bytes, 16 << 20);
    assert_eq!(limits.max_list_tag_bytes, 64 << 10);
    assert_eq!(limits.max_records, 1_000_000);
    assert_eq!(limits.max_native_id_bytes, 64 << 20);
    assert_eq!(limits.max_depth, 64);
    let handler = open_upstream();
    assert_eq!(handler.limits(), limits);
    assert_eq!(
        handler.read_options().max_records,
        mzml::ReadOptions::default().max_records
    );
    assert_eq!(handler.options(), &PeakFileOptions::default());
}

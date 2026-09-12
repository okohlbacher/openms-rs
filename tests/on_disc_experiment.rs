// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ported from `OnDiscMSExperiment_test.cpp`, the class test of
//! `OnDiscMSExperiment` / `OnDiscPeakMap`, plus native bounded-work coverage.
//!
//! Every literal taken from the upstream test or from its unmodified fixtures is
//! transcribed source review (tier 3); the resource-limit and atomicity
//! expectations are independently derived (tier 4). No C++ was executed.
//!
//! Both fixtures the upstream sections open are committed unmodified and used
//! as they are; they are recorded in
//! `tests/data/on_disc_experiment_provenance.json`:
//!
//! - the non-failure sections open `IndexedmzMLFile_1.mzML`, at
//!   `tests/data/indexed_mzml/IndexedmzMLFile_1.mzML`;
//! - the failure sections open the non-indexed `MzMLFile_1.mzML`, at
//!   `tests/data/mzml_validator/MzMLFile_1.mzML`.
//!
//! One derived fixture:
//! `a_duplicate_native_identifier_resolves_to_the_first_record` writes a copy of
//! `IndexedmzMLFile_1.mzML` into a temporary directory with the second
//! spectrum's native identifier rewritten to the first's. The rewrite is
//! byte-length preserving, so every `indexList` offset in the fixture stays
//! correct, and the helper asserts that; no upstream section covers a duplicate
//! identifier, so the input has to be made.

#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::indexed_mzml_handler::RecordReadLimits;
use openms::format::peak_options::PeakFileOptions;
use openms::kernel::NumericRange;
use openms::kernel::on_disc_experiment::{OnDiscLimits, OnDiscMSExperiment, OnDiscPeakMap};
use openms::system::file::TempDir;
use std::path::{Path, PathBuf};

const SCAN_1: &str = "controllerType=0 controllerNumber=1 scan=1";
const SCAN_2: &str = "controllerType=0 controllerNumber=1 scan=2";

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// The unmodified upstream `IndexedmzMLFile_1.mzML`: 2 spectra, 1 chromatogram.
fn indexed() -> PathBuf {
    data("indexed_mzml/IndexedmzMLFile_1.mzML")
}

/// Upstream `MzMLFile_4_indexed.mzML`: 4 spectra, index 1 is MS2 with a
/// precursor at m/z 5.55.
fn indexed_with_precursor() -> PathBuf {
    data("mzml_validator/MzMLFile_4_indexed.mzML")
}

/// The upstream `MzMLFile_1.mzML`: plain mzML with no `indexListOffset` footer,
/// which is what the failure sections need. Committed unmodified.
fn not_indexed() -> PathBuf {
    data("mzml_validator/MzMLFile_1.mzML")
}

fn opened(path: impl AsRef<Path>, skip_metadata: bool) -> OnDiscPeakMap {
    let mut experiment = OnDiscPeakMap::new();
    let parsed = experiment
        .open_file(path, skip_metadata)
        .expect("the fixture is readable mzML");
    assert!(parsed, "the fixture carries a usable index");
    experiment
}

/// `tmp` of the upstream sections: metadata loaded.
fn with_metadata() -> OnDiscPeakMap {
    opened(indexed(), false)
}

/// `tmp2` of the upstream sections: `openFile(..., true)`.
fn without_metadata() -> OnDiscPeakMap {
    opened(indexed(), true)
}

/// `failed` of the upstream sections: a readable file with no index.
fn failed_open() -> OnDiscPeakMap {
    let mut experiment = OnDiscPeakMap::new();
    let parsed = experiment.open_file(not_indexed(), false).unwrap();
    assert!(!parsed);
    experiment
}

fn range(min: f64, max: f64) -> NumericRange {
    NumericRange { min, max }
}

// ---------------------------------------------------------------------------
// START_SECTION((OnDiscMSExperiment()))
// ---------------------------------------------------------------------------

#[test]
fn default_construction_holds_no_file() {
    // The upstream section only checks that `new OnDiscPeakMap()` is not null.
    // The Rust equivalent is that the default value exists and is inert: the
    // source's own constructor documentation says "use openFile to open a file".
    let experiment = OnDiscMSExperiment::new();
    assert_eq!(experiment, OnDiscMSExperiment::default());
    assert_eq!(experiment.path(), Path::new(""));
    assert!(!experiment.is_indexed());
    assert!(experiment.is_empty());
    assert_eq!(experiment.len(), 0);
    assert_eq!(experiment.spectrum_count(), 0);
    assert_eq!(experiment.chromatogram_count(), 0);
    assert!(experiment.metadata().is_none());
    assert!(experiment.experimental_settings().is_none());
    assert!(!experiment.is_sorted_by_rt());
    assert!(!experiment.skip_xml_checks());
    assert!(!experiment.options().has_filters());
}

// ---------------------------------------------------------------------------
// START_SECTION((~OnDiscMSExperiment()))
// ---------------------------------------------------------------------------

#[test]
fn dropping_an_experiment_releases_the_file() {
    // The source destructor is implicit; the owned `std::ifstream` inside the
    // handler closes. Observable here because the temporary directory holding
    // the copy can be removed afterwards.
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("copy.mzML");
    std::fs::copy(indexed(), &path).unwrap();
    let experiment = OnDiscMSExperiment::open(&path).unwrap();
    assert_eq!(experiment.len(), 2);
    drop(experiment);
    let reopened = OnDiscMSExperiment::open(&path).unwrap();
    assert_eq!(reopened.len(), 2);
    drop(reopened);
    drop(dir);
    assert!(!path.exists());
}

// ---------------------------------------------------------------------------
// START_SECTION((OnDiscMSExperiment(const OnDiscMSExperiment& filename)))
// ---------------------------------------------------------------------------

#[test]
fn a_reopened_clone_sees_the_same_file() {
    // The source copy constructor copies the filename, the handler (which
    // reopens the file), the metadata pointer and the options. `try_clone` does
    // the same; it can fail, so it is not `Clone`.
    let source = with_metadata();
    let copy = source.try_clone().unwrap();
    let (mine, theirs) = (
        &copy.experimental_settings().unwrap().instrument,
        &source.experimental_settings().unwrap().instrument,
    );
    assert_eq!(mine.name, theirs.name);
    assert_eq!(mine.vendor, theirs.vendor);
    assert_eq!(mine.model, theirs.model);
    assert_eq!(mine.mass_analyzers.len(), theirs.mass_analyzers.len());
    assert_eq!(copy.len(), source.len());
    assert_eq!(copy, source);

    // The source copy constructor omits both native-identifier caches from its
    // initialiser list, and so does the handler's, so every by-native-id lookup
    // on a source copy raises Exception::IllegalArgument. Here the caches are
    // ordinary owned state and the copy resolves identifiers.
    let mut copy = copy;
    assert_eq!(
        copy.spectrum_by_native_id(SCAN_1).unwrap().peaks.len(),
        19914
    );
}

// ---------------------------------------------------------------------------
// START_SECTION((OnDiscMSExperiment(const std::string& filename)))
//
// Commented out in the upstream test (`OnDiscMSExperiment_test.cpp:53-58`)
// because no such constructor exists. `OnDiscMSExperiment::open` is this port's
// native equivalent, and the commented assertion `tmp.size() == 2` is what it
// checks.
// ---------------------------------------------------------------------------

#[test]
fn open_is_the_constructor_from_a_filename() {
    let experiment = OnDiscMSExperiment::open(indexed()).unwrap();
    assert_eq!(experiment.len(), 2);
    assert_eq!(experiment.path(), indexed());

    // A file with no index is the source's `false` return; this constructor
    // refuses it instead of yielding an empty experiment.
    let error = OnDiscMSExperiment::open(not_indexed()).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
    let error = OnDiscMSExperiment::open(data("indexed_mzml/does_not_exist.mzML")).unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error:?}");
}

// ---------------------------------------------------------------------------
// START_SECTION((bool operator== (const OnDiscMSExperiment& rhs) const))
// ---------------------------------------------------------------------------

#[test]
fn equality_compares_the_file_and_the_metadata() {
    let with = with_metadata();
    let without = without_metadata();
    let same = with_metadata();
    let failed = failed_open();

    assert_eq!(with, same);
    // The source falls back to comparing the two metadata pointers when either
    // is null, so a skip-metadata experiment never equals one that loaded it.
    assert_ne!(without, same);
    assert_eq!(without, without_metadata());
    assert_eq!(
        with.experimental_settings().unwrap(),
        same.experimental_settings().unwrap()
    );
    assert_ne!(with, failed);
}

// ---------------------------------------------------------------------------
// START_SECTION((bool operator!= (const OnDiscMSExperiment& rhs) const))
// ---------------------------------------------------------------------------

#[test]
fn inequality_is_the_negation_of_equality() {
    let with = with_metadata();
    let without = without_metadata();
    let same = with_metadata();
    let failed = failed_open();

    assert!(!(with != same));
    assert!(without != same);
    assert!(with != failed);
}

// ---------------------------------------------------------------------------
// START_SECTION((bool openFile(const std::string& filename, bool skipMetaData = false)))
// ---------------------------------------------------------------------------

#[test]
fn open_file_reports_whether_the_index_parsed() {
    let mut experiment = OnDiscPeakMap::new();
    assert!(experiment.open_file(indexed(), false).unwrap());
    assert!(experiment.open_file(indexed(), true).unwrap());

    let mut other = OnDiscPeakMap::new();
    assert!(other.open_file(indexed(), false).unwrap());
    assert!(other.open_file(indexed(), true).unwrap());

    let mut failed = OnDiscPeakMap::new();
    assert!(!failed.open_file(not_indexed(), false).unwrap());
    assert!(!failed.open_file(not_indexed(), true).unwrap());

    // Reopening replaces. The source's handler appends to its offset vectors
    // and never clears its native-id maps, so two successful `openFile` calls
    // on one object report the sum of both files' counts; the upstream section
    // does call `openFile` twice but its own repetitions are of the same file
    // after a failure, so the defect is invisible there.
    assert_eq!(experiment.len(), 2);
    assert_eq!(experiment.chromatogram_count(), 1);
    assert!(experiment.metadata().is_none());

    // An empty filename is the source's skip: no index parse, no metadata load.
    let mut empty = OnDiscPeakMap::new();
    assert!(!empty.open_file("", false).unwrap());
    assert!(empty.metadata().is_none());
    assert!(empty.is_empty());
}

#[test]
fn a_failed_open_leaves_the_previous_file_in_place() {
    // Native (tier 4). The source assigns `filename_` before loading the
    // metadata, so a throwing metadata load leaves the object naming a file
    // whose metadata it does not have. Everything here is committed at the end.
    let mut experiment = with_metadata();
    let error = experiment
        .open_file(data("indexed_mzml/does_not_exist.mzML"), false)
        .unwrap_err();
    assert!(matches!(error, Error::Io(_)), "{error:?}");
    assert_eq!(experiment.path(), indexed());
    assert_eq!(experiment.len(), 2);
    assert_eq!(experiment.spectrum(0).unwrap().peaks.len(), 19914);
}

// ---------------------------------------------------------------------------
// START_SECTION((bool isSortedByRT() const))
// ---------------------------------------------------------------------------

#[test]
fn sorted_by_rt_reads_the_metadata_only() {
    assert!(with_metadata().is_sorted_by_rt());
    // The metadata RTs are 0.2961 s and 0.4738 s, in that order.
    let experiment = with_metadata();
    let rts: Vec<f64> = experiment
        .metadata()
        .unwrap()
        .spectra
        .iter()
        .map(|spectrum| spectrum.rt)
        .collect();
    assert_eq!(rts, vec![0.2961, 0.4738]);
    // Without metadata the source cannot tell, and answers false.
    assert!(!without_metadata().is_sorted_by_rt());
}

// ---------------------------------------------------------------------------
// START_SECTION((Size size() const))
// ---------------------------------------------------------------------------

#[test]
fn size_is_the_indexed_spectrum_count() {
    assert_eq!(with_metadata().len(), 2);
    assert_eq!(without_metadata().len(), 2);
    assert_eq!(failed_open().len(), 0);
}

// ---------------------------------------------------------------------------
// START_SECTION((bool empty() const))
// ---------------------------------------------------------------------------

#[test]
fn empty_asks_only_about_spectra() {
    assert!(!with_metadata().is_empty());
    assert!(!without_metadata().is_empty());
    assert!(failed_open().is_empty());
}

// ---------------------------------------------------------------------------
// START_SECTION((Size getNrSpectra() const))
// ---------------------------------------------------------------------------

#[test]
fn spectrum_count_matches_the_index() {
    assert_eq!(with_metadata().spectrum_count(), 2);
    assert_eq!(without_metadata().spectrum_count(), 2);
    assert_eq!(failed_open().spectrum_count(), 0);
    // The failed open still loaded the metadata, as the source does, but the
    // counts come from the index and are therefore zero. `MzMLFile_1.mzML`
    // declares 4 spectra and 2 chromatograms, so the two answers are
    // distinguishable: the metadata is present and the index is not.
    let failed = failed_open();
    let meta = failed
        .metadata()
        .expect("a failed index parse still loads the metadata");
    assert_eq!(meta.spectra.len(), 4);
    assert_eq!(meta.chromatograms.len(), 2);
}

// ---------------------------------------------------------------------------
// START_SECTION((Size getNrChromatograms() const))
// ---------------------------------------------------------------------------

#[test]
fn chromatogram_count_matches_the_index() {
    assert_eq!(with_metadata().chromatogram_count(), 1);
    assert_eq!(without_metadata().chromatogram_count(), 1);
    assert_eq!(failed_open().chromatogram_count(), 0);
}

// ---------------------------------------------------------------------------
// START_SECTION((std::shared_ptr<const ExperimentalSettings> getExperimentalSettings() const))
// ---------------------------------------------------------------------------

#[test]
fn experimental_settings_come_from_the_metadata() {
    let experiment = with_metadata();
    let settings = experiment.experimental_settings().unwrap();
    assert_eq!(settings.instrument.name, "LTQ FT");
    assert_eq!(settings.instrument.mass_analyzers.len(), 1);

    // The source returns a null shared_ptr after `openFile(..., true)`.
    assert!(without_metadata().experimental_settings().is_none());
}

// ---------------------------------------------------------------------------
// START_SECTION((MSSpectrum operator[] (Size n) const))
// ---------------------------------------------------------------------------

#[test]
fn indexing_is_an_alias_for_spectrum() {
    // `operator[]` forwards to `getSpectrum`; `Index` cannot express it, because
    // the value is produced on demand and the fetch needs `&mut self`.
    let mut experiment = with_metadata();
    assert!(!experiment.is_empty());
    let spectrum = experiment.spectrum(0).unwrap();
    assert!(!spectrum.peaks.is_empty());
    assert_eq!(spectrum.peaks.len(), 19914);

    let mut experiment = without_metadata();
    assert!(!experiment.is_empty());
    let spectrum = experiment.spectrum(0).unwrap();
    assert!(!spectrum.peaks.is_empty());
    assert_eq!(spectrum.peaks.len(), 19914);
}

// ---------------------------------------------------------------------------
// START_SECTION((MSSpectrum getSpectrum(Size id)))
// ---------------------------------------------------------------------------

#[test]
fn spectrum_by_index_merges_peaks_into_the_metadata() {
    let mut experiment = with_metadata();
    let spectrum = experiment.spectrum(0).unwrap();
    assert_eq!(spectrum.peaks.len(), 19914);
    // The source merges the decoded arrays into the metadata record, so its
    // fields survive.
    assert_eq!(spectrum.native_id, SCAN_1);
    assert_eq!(spectrum.ms_level, 1);
    assert_eq!(spectrum.rt, 0.2961);

    let mut experiment = without_metadata();
    let first = experiment.spectrum(0).unwrap();
    assert_eq!(first.peaks.len(), 19914);
    let second = experiment.spectrum(1).unwrap();
    assert_eq!(second.peaks.len(), 19800);
    assert_eq!(second.native_id, SCAN_2);

    // Out of range. The source indexes the metadata vector with
    // `operator[]`, which is unchecked, and the offset vector through the
    // handler, which throws.
    let error = experiment.spectrum(2).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    let mut experiment = with_metadata();
    let error = experiment.spectrum(2).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");

    // No index means no peak data, even when the metadata describes spectra.
    let mut failed = failed_open();
    assert!(failed.spectrum(0).is_err());
}

// ---------------------------------------------------------------------------
// START_SECTION(OpenMS::Interfaces::SpectrumPtr getSpectrumById(Size id))
//
// The crate has no OpenSWATH `Interfaces` layer, so the pointer overload is not
// ported; its two arrays are the one peak vector here, and both upstream size
// assertions therefore land on `peaks.len()`.
// ---------------------------------------------------------------------------

#[test]
fn spectrum_arrays_are_the_peak_vector() {
    let mut experiment = with_metadata();
    let spectrum = experiment.spectrum(0).unwrap();
    assert!(!spectrum.peaks.is_empty());
    assert_eq!(spectrum.peaks.iter().filter(|p| p.mz > 0.0).count(), 19914);
    assert_eq!(spectrum.peaks.len(), 19914);

    let mut experiment = without_metadata();
    let spectrum = experiment.spectrum(0).unwrap();
    assert!(!spectrum.peaks.is_empty());
    assert_eq!(spectrum.peaks.len(), 19914);
}

// ---------------------------------------------------------------------------
// START_SECTION((MSChromatogram getChromatogram(Size id)))
// ---------------------------------------------------------------------------

#[test]
fn chromatogram_by_index_merges_points_into_the_metadata() {
    let mut experiment = with_metadata();
    assert_eq!(experiment.chromatogram_count(), 1);
    assert!(!experiment.is_empty());
    let chromatogram = experiment.chromatogram(0).unwrap();
    assert!(!chromatogram.peaks.is_empty());
    assert_eq!(chromatogram.peaks.len(), 48);
    assert_eq!(chromatogram.native_id, "TIC");

    let mut experiment = without_metadata();
    assert_eq!(experiment.chromatogram_count(), 1);
    assert!(!experiment.is_empty());
    let chromatogram = experiment.chromatogram(0).unwrap();
    assert!(!chromatogram.peaks.is_empty());
    assert_eq!(chromatogram.peaks.len(), 48);

    let error = experiment.chromatogram(1).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
}

// ---------------------------------------------------------------------------
// START_SECTION(OpenMS::Interfaces::ChromatogramPtr getChromatogramById(Size id))
//
// As the spectrum pointer overload: not ported, its time and intensity arrays
// are the one point vector.
// ---------------------------------------------------------------------------

#[test]
fn chromatogram_arrays_are_the_point_vector() {
    let mut experiment = with_metadata();
    let chromatogram = experiment.chromatogram(0).unwrap();
    assert!(!chromatogram.peaks.is_empty());
    assert_eq!(chromatogram.peaks.len(), 48);
    assert_eq!(
        chromatogram.peaks.iter().filter(|p| p.rt >= 0.0).count(),
        48
    );

    let mut experiment = without_metadata();
    let chromatogram = experiment.chromatogram(0).unwrap();
    assert!(!chromatogram.peaks.is_empty());
    assert_eq!(chromatogram.peaks.len(), 48);
}

// ---------------------------------------------------------------------------
// START_SECTION(MSChromatogram getChromatogramByNativeId(const std::string& id))
// ---------------------------------------------------------------------------

#[test]
fn chromatogram_by_native_id_resolves_through_the_metadata() {
    let mut experiment = with_metadata();
    assert!(!experiment.is_empty());
    let chromatogram = experiment.chromatogram_by_native_id("TIC").unwrap();
    assert!(!chromatogram.peaks.is_empty());
    assert_eq!(chromatogram.peaks.len(), 48);
    let error = experiment.chromatogram_by_native_id("TIK").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");

    // Without metadata the source goes straight to the handler, whose own
    // identifier map raises the same Exception::IllegalArgument.
    let mut experiment = without_metadata();
    assert!(!experiment.is_empty());
    let chromatogram = experiment.chromatogram_by_native_id("TIC").unwrap();
    assert!(!chromatogram.peaks.is_empty());
    assert_eq!(chromatogram.peaks.len(), 48);
    let error = experiment.chromatogram_by_native_id("TIK").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
}

// ---------------------------------------------------------------------------
// START_SECTION(MSMSSpectrum getSpectrumByNativeId(const std::string& id))
// ---------------------------------------------------------------------------

#[test]
fn spectrum_by_native_id_resolves_through_the_metadata() {
    let mut experiment = with_metadata();
    assert!(!experiment.is_empty());
    let first = experiment.spectrum_by_native_id(SCAN_1).unwrap();
    assert!(!first.peaks.is_empty());
    assert_eq!(first.peaks.len(), 19914);
    let second = experiment.spectrum_by_native_id(SCAN_2).unwrap();
    assert!(!second.peaks.is_empty());
    assert_eq!(second.peaks.len(), 19800);
    let error = experiment.spectrum_by_native_id("TIK").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");

    let mut experiment = without_metadata();
    assert!(!experiment.is_empty());
    let first = experiment.spectrum_by_native_id(SCAN_1).unwrap();
    assert!(!first.peaks.is_empty());
    assert_eq!(first.peaks.len(), 19914);
    let error = experiment.spectrum_by_native_id("TIK").unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
}

#[test]
fn by_native_id_ignores_the_options() {
    // Native (tier 4), documenting a source inconsistency rather than repairing
    // it: `getSpectrumByNativeId` and `getChromatogramByNativeId` never consult
    // `options_`, while `getSpectrum` and `getChromatogram` do.
    let mut experiment = with_metadata();
    experiment.options_mut().set_mz_range(range(400.0, 600.0));
    let filtered = experiment.spectrum(0).unwrap();
    let unfiltered = experiment.spectrum_by_native_id(SCAN_1).unwrap();
    assert!(filtered.peaks.len() < unfiltered.peaks.len());
    assert_eq!(unfiltered.peaks.len(), 19914);

    experiment.options_mut().set_rt_range(range(0.0, 0.1));
    let filtered = experiment.chromatogram(0).unwrap();
    let unfiltered = experiment.chromatogram_by_native_id("TIC").unwrap();
    assert!(filtered.peaks.len() < unfiltered.peaks.len());
    assert_eq!(unfiltered.peaks.len(), 48);
}

// ---------------------------------------------------------------------------
// START_SECTION((PeakFileOptions& getOptions()))
// ---------------------------------------------------------------------------

#[test]
fn mutable_options_persist() {
    let mut experiment = OnDiscPeakMap::new();
    {
        let options = experiment.options();
        assert!(!options.has_mz_range());
        assert!(!options.has_rt_range());
        assert!(!options.has_intensity_range());
    }
    experiment.options_mut().set_mz_range(range(400.0, 600.0));
    assert!(experiment.options().has_mz_range());
    assert_eq!(experiment.options().mz_range().min, 400.0);
    assert_eq!(experiment.options().mz_range().max, 600.0);
}

// ---------------------------------------------------------------------------
// START_SECTION((const PeakFileOptions& getOptions() const))
// ---------------------------------------------------------------------------

#[test]
fn shared_options_read_back_the_same_values() {
    let mut experiment = OnDiscPeakMap::new();
    experiment.options_mut().set_mz_range(range(400.0, 600.0));
    let shared: &OnDiscPeakMap = &experiment;
    let options = shared.options();
    assert!(options.has_mz_range());
    assert_eq!(options.mz_range().min, 400.0);
}

// ---------------------------------------------------------------------------
// START_SECTION((void setOptions(const PeakFileOptions& options)))
// ---------------------------------------------------------------------------

#[test]
fn set_options_replaces_the_whole_value() {
    let mut experiment = OnDiscPeakMap::new();
    let mut options = PeakFileOptions::new();
    options.set_mz_range(range(400.0, 600.0));
    options.set_intensity_range(range(100.0, 1000.0));

    experiment.set_options(options);
    assert!(experiment.options().has_mz_range());
    assert!(experiment.options().has_intensity_range());
    assert_eq!(experiment.options().mz_range().min, 400.0);
    assert_eq!(experiment.options().intensity_range().min, 100.0);

    // Native: `skip_xml_checks` on these options is inert, as in the source,
    // which forwards the flag only through `setSkipXMLChecks`.
    let mut options = PeakFileOptions::new();
    options.skip_xml_checks = true;
    experiment.set_options(options);
    assert!(!experiment.skip_xml_checks());
}

// ---------------------------------------------------------------------------
// START_SECTION(([EXTRA] Test m/z range filtering on getSpectrum))
// ---------------------------------------------------------------------------

#[test]
fn mz_range_selects_peaks_after_loading() {
    let mut experiment = with_metadata();
    let unfiltered = experiment.spectrum(0).unwrap();
    assert_eq!(unfiltered.peaks.len(), 19914);

    experiment.options_mut().set_mz_range(range(400.0, 600.0));
    let filtered = experiment.spectrum(0).unwrap();
    assert!(filtered.peaks.len() < unfiltered.peaks.len());
    assert!(!filtered.peaks.is_empty());
    for peak in &filtered.peaks {
        assert!(peak.mz >= 400.0);
        // `DRange::encloses` is half-open, so 600.0 itself is excluded.
        assert!(peak.mz < 600.0);
    }

    // Divergence: the source builds the filtered spectrum by assigning only the
    // `SpectrumSettings` base, which does not carry retention time, MS level,
    // drift time, name or the auxiliary arrays, so all of those are lost on this
    // path (`OnDiscMSExperiment.cpp:178`). This port selects peaks in place, so
    // the record keeps its metadata and its aligned annotation arrays.
    assert_eq!(filtered.rt, 0.2961);
    assert_eq!(filtered.ms_level, 1);
    assert_eq!(filtered.native_id, SCAN_1);
}

// ---------------------------------------------------------------------------
// START_SECTION(([EXTRA] Test intensity range filtering on getSpectrum))
// ---------------------------------------------------------------------------

#[test]
fn intensity_range_selects_peaks_after_loading() {
    let mut experiment = with_metadata();
    let unfiltered_size = experiment.spectrum(0).unwrap().peaks.len();

    experiment
        .options_mut()
        .set_intensity_range(range(1000.0, 1_000_000.0));
    let filtered = experiment.spectrum(0).unwrap();
    assert!(filtered.peaks.len() < unfiltered_size);
    for peak in &filtered.peaks {
        assert!(f64::from(peak.intensity) >= 1000.0);
        assert!(f64::from(peak.intensity) < 1_000_000.0);
    }
}

// ---------------------------------------------------------------------------
// START_SECTION(([EXTRA] Test copy constructor copies options))
// ---------------------------------------------------------------------------

#[test]
fn a_reopened_clone_copies_the_options() {
    let mut experiment = OnDiscPeakMap::new();
    experiment.options_mut().set_mz_range(range(400.0, 600.0));
    assert!(experiment.open_file(indexed(), false).unwrap());

    let copy = experiment.try_clone().unwrap();
    assert!(copy.options().has_mz_range());
    assert_eq!(copy.options().mz_range().min, 400.0);
    assert_eq!(copy.options().mz_range().max, 600.0);
}

// ---------------------------------------------------------------------------
// START_SECTION(([EXTRA] Test RT range filter skips loading peak data))
// ---------------------------------------------------------------------------

#[test]
fn rt_range_filter_skips_loading_peak_data() {
    let mut experiment = with_metadata();
    let first = experiment.spectrum(0).unwrap();
    let rt = first.rt;
    assert_eq!(first.peaks.len(), 19914);

    experiment
        .options_mut()
        .set_rt_range(range(rt + 1000.0, rt + 2000.0));
    let filtered = experiment.spectrum(0).unwrap();
    assert!(filtered.peaks.is_empty());
    assert_eq!(filtered.rt, rt);
    assert_eq!(filtered.ms_level, 1);
    assert_eq!(filtered.native_id, SCAN_1);

    // Without metadata the source cannot test RT before the read and applies no
    // RT filter at all to a spectrum, so the peaks arrive unfiltered.
    let mut experiment = without_metadata();
    experiment
        .options_mut()
        .set_rt_range(range(rt + 1000.0, rt + 2000.0));
    assert_eq!(experiment.spectrum(0).unwrap().peaks.len(), 19914);
}

// ---------------------------------------------------------------------------
// START_SECTION(([EXTRA] Test MS level filter skips loading peak data))
// ---------------------------------------------------------------------------

#[test]
fn ms_level_filter_skips_loading_peak_data() {
    let mut experiment = with_metadata();
    let first = experiment.spectrum(0).unwrap();
    let ms_level = first.ms_level;
    assert_eq!(ms_level, 1);
    assert_eq!(first.peaks.len(), 19914);

    experiment
        .options_mut()
        .set_ms_levels(&[i32::try_from(ms_level).unwrap() + 1])
        .unwrap();
    let filtered = experiment.spectrum(0).unwrap();
    assert!(filtered.peaks.is_empty());
    assert_eq!(filtered.ms_level, ms_level);

    // The matching level loads the peaks again.
    experiment
        .options_mut()
        .set_ms_levels(&[i32::try_from(ms_level).unwrap()])
        .unwrap();
    assert_eq!(experiment.spectrum(0).unwrap().peaks.len(), 19914);
}

// ---------------------------------------------------------------------------
// START_SECTION(([EXTRA] Test precursor m/z range filter skips loading peak data for MS2))
// ---------------------------------------------------------------------------

#[test]
fn precursor_mz_range_filter_skips_loading_peak_data() {
    // `MzMLFile_4_indexed.mzML`: MS1 at index 0, MS2 at index 1 with precursor
    // m/z 5.55, MS1 at indices 2 and 3.
    let mut experiment = opened(indexed_with_precursor(), false);
    let ms2 = experiment.spectrum(1).unwrap();
    assert_eq!(ms2.ms_level, 2);
    let unfiltered_size = ms2.peaks.len();
    assert!(unfiltered_size > 0);
    assert_eq!(unfiltered_size, 10);

    assert!(!ms2.precursors.is_empty());
    assert_eq!(ms2.precursors[0].mz, 5.55);

    experiment
        .options_mut()
        .set_precursor_mz_range(range(10.0, 20.0));
    let filtered = experiment.spectrum(1).unwrap();
    assert!(filtered.peaks.is_empty());
    assert_eq!(filtered.ms_level, 2);
    assert!(!filtered.precursors.is_empty());

    experiment
        .options_mut()
        .set_precursor_mz_range(range(5.0, 6.0));
    let passing = experiment.spectrum(1).unwrap();
    assert_eq!(passing.peaks.len(), unfiltered_size);
    assert_eq!(passing.ms_level, 2);

    // A spectrum with no precursor is not tested against the range at all.
    let ms1 = experiment.spectrum(0).unwrap();
    assert_eq!(ms1.ms_level, 1);
    assert_eq!(ms1.peaks.len(), 15);
    experiment
        .options_mut()
        .set_precursor_mz_range(range(10.0, 20.0));
    assert_eq!(experiment.spectrum(0).unwrap().peaks.len(), 15);
}

// ---------------------------------------------------------------------------
// Native coverage: materialisation, bounded work, skip_xml_checks, caches.
// ---------------------------------------------------------------------------

#[test]
fn load_experiment_materialises_every_record() {
    // The loop of `IndexedMzMLFileLoader::store`, collected into an experiment.
    let mut experiment = with_metadata();
    let loaded = experiment.load_experiment().unwrap();
    assert_eq!(loaded.spectra.len(), 2);
    assert_eq!(loaded.chromatograms.len(), 1);
    assert_eq!(loaded.spectra[0].peaks.len(), 19914);
    assert_eq!(loaded.spectra[1].peaks.len(), 19800);
    assert_eq!(loaded.chromatograms[0].peaks.len(), 48);
    assert_eq!(loaded.settings.instrument.name, "LTQ FT");
    assert!(loaded.is_sorted(false));

    // Filters apply as they do to a single fetch, and an excluded record keeps
    // its slot so the indices still match the file.
    experiment.options_mut().set_ms_levels(&[2]).unwrap();
    let filtered = experiment.load_experiment().unwrap();
    assert_eq!(filtered.spectra.len(), 2);
    assert!(filtered.spectra.iter().all(|s| s.peaks.is_empty()));
    assert_eq!(filtered.spectra[1].rt, 0.4738);

    // Without metadata the settings stay default, where the source's own
    // `store` loop dereferences the null metadata pointer instead.
    let mut experiment = without_metadata();
    let loaded = experiment.load_experiment().unwrap();
    assert_eq!(loaded.spectra.len(), 2);
    assert_eq!(loaded.settings.instrument.name, "");
}

#[test]
fn materialisation_ceilings_are_checked() {
    let limits = OnDiscLimits {
        max_materialized_points: 19_913,
        ..OnDiscLimits::default()
    };
    let mut experiment = OnDiscMSExperiment::with_limits(limits);
    assert!(experiment.open_file(indexed(), false).unwrap());
    let error = experiment.load_experiment().unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    // The rejection leaves the experiment usable.
    assert_eq!(experiment.spectrum(0).unwrap().peaks.len(), 19914);

    let limits = OnDiscLimits {
        max_records: 2,
        ..OnDiscLimits::default()
    };
    let mut experiment = OnDiscMSExperiment::with_limits(limits);
    // Skipping the metadata gets past the identifier-cache ceiling, which the
    // same two-record limit would already refuse, so the materialisation
    // ceiling is the one under test: two spectra plus one chromatogram.
    assert!(experiment.open_file(indexed(), true).unwrap());
    let error = experiment.load_experiment().unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert_eq!(experiment.spectrum(0).unwrap().peaks.len(), 19914);
}

#[test]
fn metadata_ceilings_are_checked_before_the_caches_are_built() {
    let limits = OnDiscLimits {
        max_records: 1,
        ..OnDiscLimits::default()
    };
    let mut experiment = OnDiscMSExperiment::with_limits(limits);
    let error = experiment.open_file(indexed(), false).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    // Nothing was committed.
    assert_eq!(experiment.path(), Path::new(""));
    assert!(experiment.metadata().is_none());

    let limits = OnDiscLimits {
        max_native_id_bytes: 8,
        ..OnDiscLimits::default()
    };
    let mut experiment = OnDiscMSExperiment::with_limits(limits);
    let error = experiment.open_file(indexed(), false).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert!(!experiment.is_indexed());

    // Skipping the metadata skips both ceilings, because neither cache is built.
    let mut experiment = OnDiscMSExperiment::with_limits(limits);
    assert!(experiment.open_file(indexed(), true).unwrap());
    assert_eq!(experiment.len(), 2);
}

#[test]
fn record_ceilings_are_carried_to_the_handler() {
    let limits = OnDiscLimits {
        record: RecordReadLimits {
            max_record_bytes: 1024,
            ..RecordReadLimits::default()
        },
        ..OnDiscLimits::default()
    };
    let mut experiment = OnDiscMSExperiment::with_limits(limits);
    assert!(experiment.open_file(indexed(), false).unwrap());
    let error = experiment.spectrum(0).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    // The index and the metadata are untouched by the rejected fetch.
    assert_eq!(experiment.len(), 2);
    assert_eq!(experiment.metadata().unwrap().spectra.len(), 2);
}

#[test]
fn skip_xml_checks_reaches_the_decoder_and_survives_reopening() {
    let mut experiment = with_metadata();
    assert!(!experiment.skip_xml_checks());
    experiment.set_skip_xml_checks(true);
    assert!(experiment.skip_xml_checks());
    // The flag suppresses the Base64 whitespace strip; this fixture's arrays
    // carry no whitespace, so the decoded content is identical either way.
    assert_eq!(experiment.spectrum(0).unwrap().peaks.len(), 19914);

    // The source keeps the flag on the handler, which `openFile` never resets.
    assert!(experiment.open_file(indexed(), false).unwrap());
    assert!(experiment.skip_xml_checks());
    assert_eq!(experiment.spectrum(0).unwrap().peaks.len(), 19914);

    // A clone carries it too, as the source's handler copy constructor does.
    let copy = experiment.try_clone().unwrap();
    assert!(copy.skip_xml_checks());
}

#[test]
fn the_metadata_is_loaded_without_peak_data() {
    // Source `loadMetaData_` sets `fillData=false` and no other filter, so the
    // metadata indices correspond one-to-one with the index sections.
    let experiment = with_metadata();
    let meta = experiment.metadata().unwrap();
    assert_eq!(meta.spectra.len(), 2);
    assert_eq!(meta.chromatograms.len(), 1);
    assert!(meta.spectra.iter().all(|s| s.peaks.is_empty()));
    assert!(meta.chromatograms.iter().all(|c| c.peaks.is_empty()));
    assert_eq!(meta.spectra[0].native_id, SCAN_1);
    assert_eq!(meta.spectra[1].native_id, SCAN_2);
    assert_eq!(meta.chromatograms[0].native_id, "TIC");

    // The facade's own filters never reach the metadata load: the source
    // deliberately loads everything so that the indices stay valid.
    let mut experiment = OnDiscPeakMap::new();
    experiment.options_mut().set_ms_levels(&[2]).unwrap();
    experiment.options_mut().set_rt_range(range(100.0, 200.0));
    assert!(experiment.open_file(indexed(), false).unwrap());
    assert_eq!(experiment.metadata().unwrap().spectra.len(), 2);
    assert_eq!(experiment.len(), 2);
}

#[test]
fn half_open_range_endpoints_follow_drange_encloses() {
    // `DRange<1>::encloses` rejects a value equal to the maximum
    // (`DATASTRUCTURES/DRange.h:158`), which is why an RT range that starts at
    // the spectrum's own RT keeps it and one that ends there does not.
    let mut experiment = with_metadata();
    experiment.options_mut().set_rt_range(range(0.2961, 1.0));
    assert_eq!(experiment.spectrum(0).unwrap().peaks.len(), 19914);

    experiment.options_mut().set_rt_range(range(0.0, 0.2961));
    assert!(experiment.spectrum(0).unwrap().peaks.is_empty());
}

/// A copy of `IndexedmzMLFile_1.mzML` in `dir` whose second spectrum carries
/// the *first* spectrum's native identifier, in both the record and the index.
///
/// The identifier appears exactly twice in the fixture, once on the
/// `<spectrum>` element and once as the `idRef` of its `<indexList>` entry, and
/// `scan=1` and `scan=2` are the same length, so rewriting both leaves every
/// byte offset in the index correct. Nothing else about the file changes:
/// spectrum 0 still has 19914 peaks and spectrum 1 still has 19800.
fn duplicate_native_id_fixture(dir: &TempDir) -> PathBuf {
    let text = std::fs::read_to_string(indexed()).expect("the fixture is UTF-8");
    assert_eq!(
        text.matches(SCAN_2).count(),
        2,
        "the rewrite expects exactly the <spectrum> element and its <offset>"
    );
    assert_eq!(
        SCAN_1.len(),
        SCAN_2.len(),
        "the rewrite must preserve length"
    );
    let rewritten = text.replace(SCAN_2, SCAN_1);
    assert_eq!(rewritten.len(), text.len(), "byte offsets must not move");
    let path = dir.path().join("duplicate_native_id.mzML");
    std::fs::write(&path, rewritten).unwrap();
    path
}

#[test]
fn a_duplicate_native_identifier_resolves_to_the_first_record() {
    // Source `getMetaSpectrumById_` fills its map with
    // `unordered_map::emplace`, which does not overwrite, so for two records
    // sharing an identifier the first wins. No upstream section covers that:
    // `IndexedmzMLFile_1.mzML`'s identifiers are unique, so the input is derived
    // here, and the derivation shows where the rule is and is not observable.
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let duplicated = duplicate_native_id_fixture(&dir);

    // Loading the metadata is refused outright: mzML's schema makes record
    // identifiers unique per kind (`KEY_SPECTRUM_ID`, `KEY_CHROMATOGRAM_ID`)
    // and this crate's reader enforces that, keyed by `(tag, id)`
    // (`src/format/mzml.rs`, "empty or duplicate record id"). So the facade's
    // *own* metadata cache never sees a duplicate through `open_file`; its
    // first-wins `entry().or_insert_with()` is defensive, not a path a file can
    // reach.
    // The refusal is atomic, as `a_failed_open_leaves_the_previous_file_in_place`
    // asserts for the general case.
    let mut refused = OnDiscPeakMap::new();
    let error = refused.open_file(&duplicated, false).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
    assert!(refused.metadata().is_none());
    assert_eq!(refused.path(), Path::new(""));

    // Where first-wins *is* observable is the index: `skip_metadata` loads no
    // metadata, so the identifier goes straight to the handler's own
    // `<indexList>` map, which is built with `or_insert` for the same reason.
    // Both offsets are still in the index, and the duplicated identifier
    // resolves to the first of them - spectrum 0, with 19914 peaks, not
    // spectrum 1 with 19800.
    let mut experiment = OnDiscPeakMap::new();
    assert!(experiment.open_file(&duplicated, true).unwrap());
    assert_eq!(experiment.len(), 2, "both records are still indexed");
    let resolved = experiment.spectrum_by_native_id(SCAN_1).unwrap();
    assert_eq!(
        resolved.peaks.len(),
        19914,
        "the first index entry won, not the second"
    );

    // The identifier the rewrite overwrote is gone from the index entirely.
    let error = experiment.spectrum_by_native_id(SCAN_2).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");

    // First-wins loses no data: it only decides what one identifier resolves
    // to, and both records stay reachable by index.
    assert_eq!(experiment.spectrum(0).unwrap().peaks.len(), 19914);
    assert_eq!(experiment.spectrum(1).unwrap().peaks.len(), 19800);
}

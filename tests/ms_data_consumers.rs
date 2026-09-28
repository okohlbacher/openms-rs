// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The in-memory MS data consumers of `FORMAT/DATAACCESS`:
//! `NoopMSDataConsumer`, `MSDataStoringConsumer`, `MSDataChainingConsumer`,
//! `MSDataTransformingConsumer` and `MSDataAggregatingConsumer`.
//!
//! Every section of the five upstream class tests
//! (`src/tests/class_tests/openms/source/<Class>_test.cpp` at core `bc9cc12`)
//! is transcribed below with its line numbers: source review, tier 3. The
//! sections that open `MzMLFile_1.mzML` use the unmodified upstream fixture at
//! `tests/data/mzml_validator/MzMLFile_1.mzML`. The `native_*` tests are
//! independently derived (tier 4): the source has no soft stop and no
//! fallible callbacks, so their behaviour is this port's own contract.
//!
//! Two adaptations, both documented at their use:
//!
//! - `MSSpectrum` has no comment field in this port (docs/MZXML_SUPPORT.md,
//!   item 16), so the class tests' `setComment` is carried by a metadata entry
//!   `comment`. What the sections check — that a stored or merged spectrum
//!   keeps the metadata of the one it came from — is the same for any
//!   metadata.
//! - `MSDataAggregatingConsumer` merges with `SpectrumAddition::addUpSpectra`
//!   upstream, which is not ported. Its sections run with [`add_up_on_grid`], a
//!   test double that adds intensities on a shared m/z grid. The literals
//!   17/120/230 are the class test's; any merge that adds intensities on the
//!   fixture's shared grid reaches them, so they check the consumer's grouping
//!   and metadata plumbing, not `SpectrumAddition`.

#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::ms_data_aggregating_consumer::{MSDataAggregatingConsumer, RT_TOLERANCE};
use openms::format::ms_data_chaining_consumer::MSDataChainingConsumer;
use openms::format::ms_data_storing_consumer::MSDataStoringConsumer;
use openms::format::ms_data_transforming_consumer::MSDataTransformingConsumer;
use openms::format::mzml;
use openms::format::noop_ms_data_consumer::NoopMSDataConsumer;
use openms::interfaces::MSDataConsumer;
use openms::kernel::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};
use openms::metadata::{ExperimentalSettings, MetaValue};
use std::fs::File;
use std::io::BufReader;
use std::ops::ControlFlow;

/// The upstream `MzMLFile_1.mzML`, loaded as `MzMLFile().load(...)` does.
fn mzml_file_1() -> MSExperiment {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/mzml_validator/MzMLFile_1.mzML"
    );
    mzml::read(BufReader::new(File::open(path).unwrap())).unwrap()
}

/// The class tests' `setComment`, as a metadata entry.
fn set_comment(spectrum: &mut MSSpectrum, comment: &str) {
    spectrum
        .metadata
        .insert("comment".into(), MetaValue::from(comment));
}

fn comment(spectrum: &MSSpectrum) -> MetaValue {
    spectrum.metadata["comment"].clone()
}

fn named(name: &str, rt: f64) -> MSSpectrum {
    MSSpectrum {
        name: name.into(),
        rt,
        ..MSSpectrum::default()
    }
}

/// Test double for `SpectrumAddition::addUpSpectra(list, -1, true)`: adds the
/// intensities of peaks at identical m/z, in first-seen m/z order.
fn add_up_on_grid(spectra: &[MSSpectrum]) -> openms::Result<MSSpectrum> {
    let mut merged = MSSpectrum::default();
    for spectrum in spectra {
        for peak in &spectrum.peaks {
            match merged.peaks.iter_mut().find(|p| p.mz == peak.mz) {
                Some(existing) => existing.intensity += peak.intensity,
                None => merged.peaks.push(*peak),
            }
        }
    }
    Ok(merged)
}

fn continued(flow: openms::Result<ControlFlow<()>>) {
    assert_eq!(flow.unwrap(), ControlFlow::Continue(()));
}

// ---------------------------------------------------------------------------
// NoopMSDataConsumer_test.cpp
// ---------------------------------------------------------------------------

#[test]
fn noop_constructs() {
    // NoopMSDataConsumer_test.cpp:25-36: construction and destruction. The
    // four member sections (:38-60) are empty `// TODO` bodies upstream.
    let _consumer = NoopMSDataConsumer::new();
    let _default = NoopMSDataConsumer;
}

#[test]
fn native_noop_accepts_everything_and_changes_nothing() {
    let experiment = mzml_file_1();
    let mut consumer = NoopMSDataConsumer::new();
    consumer.set_expected_size(usize::MAX, usize::MAX).unwrap();
    consumer
        .set_experimental_settings(&experiment.settings)
        .unwrap();
    let mut spectrum = experiment.spectra[0].clone();
    continued(consumer.consume_spectrum(&mut spectrum));
    assert_eq!(spectrum, experiment.spectra[0]);
    let mut chromatogram = experiment.chromatograms[0].clone();
    continued(consumer.consume_chromatogram(&mut chromatogram));
    assert_eq!(chromatogram, experiment.chromatograms[0]);
}

// ---------------------------------------------------------------------------
// MSDataStoringConsumer_test.cpp
// ---------------------------------------------------------------------------

#[test]
fn storing_constructs() {
    // MSDataStoringConsumer_test.cpp:31-38.
    let consumer = MSDataStoringConsumer::new();
    assert!(consumer.data().spectra.is_empty());
}

#[test]
fn storing_consume_spectrum() {
    // MSDataStoringConsumer_test.cpp:40-72.
    let mut consumer = MSDataStoringConsumer::new();
    let mut s = MSSpectrum::default();
    for (name, comm, rt) in [
        ("spec1", "comm1", 5.0),
        ("spec2", "comm2", 15.0),
        ("spec3", "comm3", 25.0),
    ] {
        s.name = name.into();
        set_comment(&mut s, comm);
        s.rt = rt;
        continued(consumer.consume_spectrum(&mut s));
    }
    assert_eq!(consumer.data().spectra.len(), 3);
    assert_eq!(consumer.data().chromatograms.len(), 0);
    assert_eq!(consumer.data().spectra[0].name, "spec1");
    assert_eq!(consumer.data().spectra[1].name, "spec2");
    assert_eq!(consumer.data().spectra[2].name, "spec3");
    assert_eq!(comment(&consumer.data().spectra[0]), MetaValue::from("comm1"));
    assert_eq!(comment(&consumer.data().spectra[1]), MetaValue::from("comm2"));
    assert_eq!(comment(&consumer.data().spectra[2]), MetaValue::from("comm3"));
}

#[test]
fn storing_consume_chromatogram() {
    // MSDataStoringConsumer_test.cpp:74-88.
    let mut consumer = MSDataStoringConsumer::new();
    let mut c = MSChromatogram {
        native_id: "testid".into(),
        ..MSChromatogram::default()
    };
    continued(consumer.consume_chromatogram(&mut c));
    assert_eq!(consumer.data().spectra.len(), 0);
    assert_eq!(consumer.data().chromatograms.len(), 1);
    assert_eq!(consumer.data().chromatograms[0].native_id, "testid");
}

#[test]
fn storing_set_experimental_settings_keeps_the_data() {
    // MSDataStoringConsumer_test.cpp:94-120: settings set after two records
    // replace the settings and keep both records.
    let mut consumer = MSDataStoringConsumer::new();
    consumer.set_expected_size(1, 1).unwrap();
    let mut c = MSChromatogram {
        native_id: "testid".into(),
        ..MSChromatogram::default()
    };
    continued(consumer.consume_chromatogram(&mut c));
    let mut spec = named("spec1", 5.0);
    continued(consumer.consume_spectrum(&mut spec));
    let settings = ExperimentalSettings {
        comment: "mySettings".into(),
        ..ExperimentalSettings::default()
    };
    consumer.set_experimental_settings(&settings).unwrap();
    assert_eq!(consumer.data().spectra.len(), 1);
    assert_eq!(consumer.data().chromatograms.len(), 1);
    assert_eq!(consumer.data().settings.comment, "mySettings");
}

#[test]
fn native_storing_expected_size_is_a_bounded_hint() {
    // An expected size is a reservation only; a count taken from an untrusted
    // file must not become an allocation of that size.
    let mut consumer = MSDataStoringConsumer::new();
    consumer.set_expected_size(usize::MAX, usize::MAX).unwrap();
    let mut s = named("spec1", 1.0);
    continued(consumer.consume_spectrum(&mut s));
    // The caller's record is unchanged: the consumer stores a copy.
    assert_eq!(s, named("spec1", 1.0));
    let experiment = consumer.into_data();
    assert_eq!(experiment.spectra, vec![named("spec1", 1.0)]);
    assert_eq!(experiment.sql_run_id, 0);
}

// ---------------------------------------------------------------------------
// MSDataTransformingConsumer_test.cpp
// ---------------------------------------------------------------------------

#[test]
fn transforming_constructs() {
    // MSDataTransformingConsumer_test.cpp:35-42.
    let consumer = MSDataTransformingConsumer::new();
    assert!(!consumer.has_spectra_processing_func());
    assert!(!consumer.has_chromatogram_processing_func());
    assert!(!consumer.has_experimental_settings_func());
}

#[test]
fn transforming_without_functions_changes_nothing() {
    // MSDataTransformingConsumer_test.cpp:44-72.
    let expc = mzml_file_1();
    let mut consumer = MSDataTransformingConsumer::new();

    let mut exp = expc.clone();
    assert!(!exp.spectra.is_empty());
    let first_spectrum = exp.spectra[0].clone();
    consumer.set_expected_size(2, 0).unwrap();
    continued(consumer.consume_spectrum(&mut exp.spectra[0]));
    assert_eq!(first_spectrum, exp.spectra[0]);

    let mut exp = expc.clone();
    assert!(!exp.chromatograms.is_empty());
    let first_chromatogram = exp.chromatograms[0].clone();
    consumer.set_expected_size(0, 1).unwrap();
    continued(consumer.consume_chromatogram(&mut exp.chromatograms[0]));
    assert_eq!(first_chromatogram, exp.chromatograms[0]);

    // :78-88: setExperimentalSettings without a function, NOT_TESTABLE there.
    consumer.set_expected_size(2, 0).unwrap();
    consumer
        .set_experimental_settings(&ExperimentalSettings::default())
        .unwrap();
}

#[test]
fn transforming_set_spectra_processing_func() {
    // MSDataTransformingConsumer_test.cpp:90-113.
    let mut exp = mzml_file_1();
    exp.spectra[0].sort_by_position().unwrap();
    let first_spectrum = exp.spectra[0].clone();
    let mut consumer = MSDataTransformingConsumer::new();
    consumer.set_expected_size(2, 0).unwrap();
    consumer.set_spectra_processing_func(|s: &mut MSSpectrum| s.sort_by_intensity(false));
    continued(consumer.consume_spectrum(&mut exp.spectra[0]));
    assert_ne!(first_spectrum, exp.spectra[0]);
    assert!(first_spectrum.is_sorted());
    assert!(!exp.spectra[0].is_sorted());
}

#[test]
fn transforming_set_chromatogram_processing_func() {
    // MSDataTransformingConsumer_test.cpp:115-138.
    let mut exp = mzml_file_1();
    exp.chromatograms[0].sort_by_position().unwrap();
    let first_chromatogram = exp.chromatograms[0].clone();
    let mut consumer = MSDataTransformingConsumer::new();
    consumer.set_expected_size(0, 1).unwrap();
    consumer
        .set_chromatogram_processing_func(|c: &mut MSChromatogram| c.sort_by_intensity(false));
    continued(consumer.consume_chromatogram(&mut exp.chromatograms[0]));
    assert_ne!(first_chromatogram, exp.chromatograms[0]);
    assert!(first_chromatogram.is_sorted());
    assert!(!exp.chromatograms[0].is_sorted());
}

#[test]
fn native_transforming_functions_capture_state_fail_and_clear() {
    let exp = mzml_file_1();
    let mut seen_spectra = 0usize;
    let mut seen_chromatograms = 0usize;
    let mut seen_settings = Vec::new();
    {
        let mut consumer = MSDataTransformingConsumer::new();
        consumer.set_spectra_processing_func(|_: &mut MSSpectrum| {
            seen_spectra += 1;
            Ok(())
        });
        consumer.set_chromatogram_processing_func(|_: &mut MSChromatogram| {
            seen_chromatograms += 1;
            Ok(())
        });
        // The header calls setExperimentalSettings "ignored"; the
        // implementation calls the settings function (MSDataTransformingConsumer.cpp).
        consumer.set_experimental_settings_func(|settings: &ExperimentalSettings| {
            seen_settings.push(settings.comment.clone());
            Ok(())
        });
        let settings = ExperimentalSettings {
            comment: "run".into(),
            ..ExperimentalSettings::default()
        };
        consumer.set_experimental_settings(&settings).unwrap();
        let mut s = exp.spectra[0].clone();
        continued(consumer.consume_spectrum(&mut s));
        continued(consumer.consume_spectrum(&mut s));
        let mut c = exp.chromatograms[0].clone();
        continued(consumer.consume_chromatogram(&mut c));

        // A failing function fails the consumption.
        consumer.set_spectra_processing_func(|_: &mut MSSpectrum| {
            Err(Error::InvalidValue("refused".into()))
        });
        assert!(matches!(
            consumer.consume_spectrum(&mut s),
            Err(Error::InvalidValue(_))
        ));

        // Clearing is the source's nullptr: records pass unchanged again.
        consumer.clear_spectra_processing_func();
        consumer.clear_chromatogram_processing_func();
        consumer.clear_experimental_settings_func();
        assert!(!consumer.has_spectra_processing_func());
        continued(consumer.consume_spectrum(&mut s));
        continued(consumer.consume_chromatogram(&mut c));
        consumer.set_experimental_settings(&settings).unwrap();
    }
    assert_eq!(seen_spectra, 2);
    assert_eq!(seen_chromatograms, 1);
    assert_eq!(seen_settings, vec!["run".to_string()]);
}

// ---------------------------------------------------------------------------
// MSDataChainingConsumer_test.cpp
// ---------------------------------------------------------------------------

#[test]
fn chaining_constructs() {
    // MSDataChainingConsumer_test.cpp:37-51.
    let chain = MSDataChainingConsumer::new();
    assert!(chain.is_empty());
    let chain = MSDataChainingConsumer::from_consumers(Vec::new());
    assert_eq!(chain.len(), 0);
}

#[test]
fn chaining_consume_spectrum_through_noops() {
    // MSDataChainingConsumer_test.cpp:53-77.
    let (mut a, mut b, mut c) = (
        NoopMSDataConsumer::new(),
        NoopMSDataConsumer::new(),
        NoopMSDataConsumer::new(),
    );
    let mut chain = MSDataChainingConsumer::from_consumers(vec![&mut a, &mut b, &mut c]);
    let mut exp = mzml_file_1();
    assert!(!exp.spectra.is_empty());
    let first_spectrum = exp.spectra[0].clone();
    chain.set_expected_size(2, 0).unwrap();
    continued(chain.consume_spectrum(&mut exp.spectra[0]));
    assert_eq!(first_spectrum, exp.spectra[0]);
}

#[test]
fn chaining_consume_spectrum_through_a_transformation() {
    // MSDataChainingConsumer_test.cpp:79-122.
    let mut transforming = MSDataTransformingConsumer::new();
    transforming.set_expected_size(2, 0).unwrap();
    transforming.set_spectra_processing_func(|s: &mut MSSpectrum| s.sort_by_intensity(false));
    let (mut a, mut b) = (NoopMSDataConsumer::new(), NoopMSDataConsumer::new());
    let mut exp = mzml_file_1();
    assert!(!exp.spectra.is_empty());
    let first_spectrum = exp.spectra[0].clone();
    {
        let steps: Vec<&mut dyn MSDataConsumer> = vec![&mut a, &mut transforming, &mut b];
        let mut chain = MSDataChainingConsumer::from_consumers(steps);
        chain.set_expected_size(2, 0).unwrap();
        continued(chain.consume_spectrum(&mut exp.spectra[0]));
    }
    assert_ne!(first_spectrum, exp.spectra[0]);
    assert!(first_spectrum.is_sorted());
    assert!(!exp.spectra[0].is_sorted());

    // :109-114: the chain did not own the transforming consumer, which still
    // works once the chain is gone.
    continued(transforming.consume_spectrum(&mut exp.spectra[0]));
    assert!(first_spectrum.is_sorted());
    assert!(!exp.spectra[0].is_sorted());
}

#[test]
fn chaining_consume_chromatogram_through_noops() {
    // MSDataChainingConsumer_test.cpp:124-148.
    let (mut a, mut b, mut c) = (
        NoopMSDataConsumer::new(),
        NoopMSDataConsumer::new(),
        NoopMSDataConsumer::new(),
    );
    let mut chain = MSDataChainingConsumer::from_consumers(vec![&mut a, &mut b, &mut c]);
    let mut exp = mzml_file_1();
    assert!(!exp.chromatograms.is_empty());
    let first_chromatogram = exp.chromatograms[0].clone();
    chain.set_expected_size(0, 1).unwrap();
    continued(chain.consume_chromatogram(&mut exp.chromatograms[0]));
    assert_eq!(first_chromatogram, exp.chromatograms[0]);
}

#[test]
fn chaining_consume_chromatogram_through_a_transformation() {
    // MSDataChainingConsumer_test.cpp:150-185.
    let mut transforming = MSDataTransformingConsumer::new();
    transforming.set_expected_size(2, 0).unwrap();
    transforming
        .set_chromatogram_processing_func(|c: &mut MSChromatogram| c.sort_by_intensity(false));
    let (mut a, mut b) = (NoopMSDataConsumer::new(), NoopMSDataConsumer::new());
    let steps: Vec<&mut dyn MSDataConsumer> = vec![&mut a, &mut transforming, &mut b];
    let mut chain = MSDataChainingConsumer::from_consumers(steps);
    let mut exp = mzml_file_1();
    assert!(!exp.chromatograms.is_empty());
    let first_chromatogram = exp.chromatograms[0].clone();
    chain.set_expected_size(0, 1).unwrap();
    continued(chain.consume_chromatogram(&mut exp.chromatograms[0]));
    assert_ne!(first_chromatogram, exp.chromatograms[0]);
    assert!(first_chromatogram.is_sorted());
    assert!(!exp.chromatograms[0].is_sorted());
}

#[test]
fn chaining_set_experimental_settings() {
    // MSDataChainingConsumer_test.cpp:191-202.
    let mut chain = MSDataChainingConsumer::new();
    chain.set_expected_size(2, 0).unwrap();
    chain
        .set_experimental_settings(&ExperimentalSettings::default())
        .unwrap();
}

#[test]
fn chaining_append_consumer() {
    // MSDataChainingConsumer_test.cpp:204-239.
    let mut transforming = MSDataTransformingConsumer::new();
    transforming.set_expected_size(2, 0).unwrap();
    transforming.set_spectra_processing_func(|s: &mut MSSpectrum| s.sort_by_intensity(false));
    let (mut a, mut b) = (NoopMSDataConsumer::new(), NoopMSDataConsumer::new());
    let mut chain = MSDataChainingConsumer::from_consumers(vec![&mut a, &mut b]);
    chain.append_consumer(&mut transforming);
    assert_eq!(chain.len(), 3);
    let mut exp = mzml_file_1();
    assert!(!exp.spectra.is_empty());
    let first_spectrum = exp.spectra[0].clone();
    chain.set_expected_size(2, 0).unwrap();
    continued(chain.consume_spectrum(&mut exp.spectra[0]));
    assert_ne!(first_spectrum, exp.spectra[0]);
    assert!(first_spectrum.is_sorted());
    assert!(!exp.spectra[0].is_sorted());
}

#[test]
fn native_chaining_forwards_settings_and_sizes_in_order() {
    let mut first = MSDataStoringConsumer::new();
    let mut second = MSDataStoringConsumer::new();
    let settings = ExperimentalSettings {
        comment: "chained".into(),
        ..ExperimentalSettings::default()
    };
    {
        let mut chain = MSDataChainingConsumer::from_consumers(vec![&mut first, &mut second]);
        chain.set_expected_size(1, 1).unwrap();
        chain.set_experimental_settings(&settings).unwrap();
        let mut s = named("spec1", 5.0);
        continued(chain.consume_spectrum(&mut s));
    }
    assert_eq!(first.data().settings.comment, "chained");
    assert_eq!(second.data().settings.comment, "chained");
    assert_eq!(first.data().spectra, second.data().spectra);
}

#[test]
fn native_chaining_stops_at_the_first_break_or_error() {
    // A consumer that refuses spectra named "stop" with Break and spectra
    // named "fail" with an error.
    let mut gate = MSDataTransformingConsumer::new();
    let mut after = MSDataStoringConsumer::new();
    let mut before = MSDataStoringConsumer::new();
    struct Gate;
    impl MSDataConsumer for Gate {
        fn set_expected_size(&mut self, _: usize, _: usize) -> openms::Result<()> {
            Err(Error::InvalidValue("size".into()))
        }
        fn set_experimental_settings(&mut self, _: &ExperimentalSettings) -> openms::Result<()> {
            Ok(())
        }
        fn consume_spectrum(&mut self, s: &mut MSSpectrum) -> openms::Result<ControlFlow<()>> {
            match s.name.as_str() {
                "stop" => Ok(ControlFlow::Break(())),
                "fail" => Err(Error::InvalidValue("fail".into())),
                _ => Ok(ControlFlow::Continue(())),
            }
        }
        fn consume_chromatogram(
            &mut self,
            _: &mut MSChromatogram,
        ) -> openms::Result<ControlFlow<()>> {
            Ok(ControlFlow::Break(()))
        }
    }
    gate.set_spectra_processing_func(|s: &mut MSSpectrum| {
        s.name.push('!');
        Ok(())
    });
    let mut stop = Gate;
    {
        let mut chain = MSDataChainingConsumer::from_consumers(vec![
            &mut before,
            &mut stop,
            &mut gate,
            &mut after,
        ]);
        let mut s = named("go", 1.0);
        continued(chain.consume_spectrum(&mut s));
        assert_eq!(s.name, "go!");
        let mut s = named("stop", 2.0);
        assert_eq!(chain.consume_spectrum(&mut s).unwrap(), ControlFlow::Break(()));
        assert_eq!(s.name, "stop");
        let mut s = named("fail", 3.0);
        assert!(chain.consume_spectrum(&mut s).is_err());
        let mut c = MSChromatogram::default();
        assert_eq!(chain.consume_chromatogram(&mut c).unwrap(), ControlFlow::Break(()));
        assert!(chain.set_expected_size(1, 1).is_err());
    }
    // Consumers before the stopping one saw every record; those after it only
    // the record nobody stopped.
    let names: Vec<_> = before.data().spectra.iter().map(|s| s.name.clone()).collect();
    assert_eq!(names, ["go", "stop", "fail"]);
    let names: Vec<_> = after.data().spectra.iter().map(|s| s.name.clone()).collect();
    assert_eq!(names, ["go!"]);
    assert_eq!(before.data().chromatograms.len(), 1);
    assert!(after.data().chromatograms.is_empty());
}

// ---------------------------------------------------------------------------
// MSDataAggregatingConsumer_test.cpp
// ---------------------------------------------------------------------------

#[test]
fn aggregating_no_adding_up() {
    // MSDataAggregatingConsumer_test.cpp:42-70.
    let mut storage = MSDataStoringConsumer::new();
    {
        let mut agg = MSDataAggregatingConsumer::new(&mut storage, add_up_on_grid);
        let mut s = MSSpectrum::default();
        for (name, rt) in [("spec1", 5.0), ("spec2", 15.0), ("spec3", 25.0)] {
            s.name = name.into();
            s.rt = rt;
            continued(agg.consume_spectrum(&mut s));
        }
        // "note how we can / have to destroy the aggregate consumer to ensure
        // it flushes the data": here, `finish`.
        continued(agg.finish());
    }
    assert_eq!(storage.data().spectra.len(), 3);
    assert_eq!(storage.data().chromatograms.len(), 0);
    assert_eq!(storage.data().spectra[0].name, "spec1");
    assert_eq!(storage.data().spectra[1].name, "spec2");
    assert_eq!(storage.data().spectra[2].name, "spec3");
}

#[test]
fn aggregating_empty_spectra() {
    // MSDataAggregatingConsumer_test.cpp:72-115.
    let mut storage = MSDataStoringConsumer::new();
    {
        let mut agg = MSDataAggregatingConsumer::new(&mut storage, add_up_on_grid);
        let mut s = MSSpectrum::default();
        for (name, comm, rt) in [
            ("spec1", "comm1", 5.0),
            ("spec2", "comm2", 5.0),
            ("spec3", "comm3", 25.0),
            ("spec4", "comm4", 25.0),
            ("spec5", "comm5", 35.0),
        ] {
            s.name = name.into();
            set_comment(&mut s, comm);
            s.rt = rt;
            continued(agg.consume_spectrum(&mut s));
        }
        // Dropping flushes, as the source's destructor does.
    }
    let data = storage.data();
    assert_eq!(data.spectra.len(), 3);
    assert_eq!(data.chromatograms.len(), 0);
    assert_eq!(data.spectra[0].name, "spec1");
    assert_eq!(data.spectra[1].name, "spec3");
    assert_eq!(data.spectra[2].name, "spec5");
    assert_eq!(comment(&data.spectra[0]), MetaValue::from("comm1"));
    assert_eq!(comment(&data.spectra[1]), MetaValue::from("comm3"));
    assert_eq!(comment(&data.spectra[2]), MetaValue::from("comm5"));
}

#[test]
fn aggregating_full_spectra() {
    // MSDataAggregatingConsumer_test.cpp:117-185, with the merge double.
    let mut storage = MSDataStoringConsumer::new();
    {
        let mut agg = MSDataAggregatingConsumer::new(&mut storage, add_up_on_grid);
        let mut s = MSSpectrum::default();
        s.name = "spec1".into();
        set_comment(&mut s, "comm1");
        s.rt = 5.0;
        s.peaks = vec![
            Peak1D::new(5.0, 7.0),
            Peak1D::new(10.0, 20.0),
            Peak1D::new(15.0, 30.0),
        ];
        continued(agg.consume_spectrum(&mut s));
        s.clear(true);
        s.name = "spec2".into();
        set_comment(&mut s, "comm2");
        s.rt = 5.0;
        s.peaks = vec![
            Peak1D::new(5.0, 10.0),
            Peak1D::new(10.0, 100.0),
            Peak1D::new(15.0, 200.0),
        ];
        continued(agg.consume_spectrum(&mut s));
        for (name, comm, rt) in [
            ("spec3", "comm3", 25.0),
            ("spec4", "comm4", 25.0),
            ("spec5", "comm5", 35.0),
        ] {
            s.clear(true);
            s.name = name.into();
            set_comment(&mut s, comm);
            s.rt = rt;
            continued(agg.consume_spectrum(&mut s));
        }
        s.clear(true);
        continued(agg.finish());
    }
    let data = storage.data();
    assert_eq!(data.spectra.len(), 3);
    assert_eq!(data.chromatograms.len(), 0);
    assert_eq!(data.spectra[0].name, "spec1");
    assert_eq!(data.spectra[1].name, "spec3");
    assert_eq!(data.spectra[2].name, "spec5");
    assert_eq!(comment(&data.spectra[0]), MetaValue::from("comm1"));
    assert_eq!(comment(&data.spectra[1]), MetaValue::from("comm3"));
    assert_eq!(comment(&data.spectra[2]), MetaValue::from("comm5"));
    let snew = &data.spectra[0];
    assert_eq!(snew.peaks.len(), 3);
    // TEST_REAL_SIMILAR, the class test's default relative tolerance.
    let similar = |a: f64, b: f64| (a - b).abs() <= 1e-5 * b.abs().max(1.0);
    assert!(similar(snew.peaks[0].mz, 5.0));
    assert!(similar(f64::from(snew.peaks[0].intensity), 17.0));
    assert!(similar(snew.peaks[1].mz, 10.0));
    assert!(similar(f64::from(snew.peaks[1].intensity), 120.0));
    assert!(similar(snew.peaks[2].mz, 15.0));
    assert!(similar(f64::from(snew.peaks[2].intensity), 230.0));
}

#[test]
fn aggregating_consume_chromatogram() {
    // MSDataAggregatingConsumer_test.cpp:189-206.
    let mut storage = MSDataStoringConsumer::new();
    {
        let mut agg = MSDataAggregatingConsumer::new(&mut storage, add_up_on_grid);
        let mut c = MSChromatogram {
            native_id: "testid".into(),
            ..MSChromatogram::default()
        };
        continued(agg.consume_chromatogram(&mut c));
    }
    assert_eq!(storage.data().spectra.len(), 0);
    assert_eq!(storage.data().chromatograms.len(), 1);
    assert_eq!(storage.data().chromatograms[0].native_id, "testid");
}

#[test]
fn native_aggregating_does_not_forward_sizes_or_settings() {
    // Both overrides are empty in the source header.
    let mut storage = MSDataStoringConsumer::new();
    {
        let mut agg = MSDataAggregatingConsumer::new(&mut storage, add_up_on_grid);
        agg.set_expected_size(10, 10).unwrap();
        let settings = ExperimentalSettings {
            comment: "not forwarded".into(),
            ..ExperimentalSettings::default()
        };
        agg.set_experimental_settings(&settings).unwrap();
        continued(agg.finish());
    }
    assert_eq!(storage.data().settings, ExperimentalSettings::default());
    assert!(storage.data().spectra.is_empty());
}

#[test]
fn native_aggregating_groups_by_the_previous_retention_time() {
    // The tolerance is to the previous spectrum, so a slow drift chains one
    // group; a NaN never joins one.
    let mut storage = MSDataStoringConsumer::new();
    let mut group_sizes = Vec::new();
    {
        let mut agg = MSDataAggregatingConsumer::new(&mut storage, |spectra: &[MSSpectrum]| {
            group_sizes.push(spectra.len());
            add_up_on_grid(spectra)
        });
        let step = RT_TOLERANCE * 0.6;
        for (name, rt) in [
            ("a", 5.0),
            ("b", 5.0 + step),
            ("c", 5.0 + 2.0 * step),
            // Two tolerances past "c": one tolerance exactly would sit on the
            // strict boundary, where the rounding of `5.0 + x` decides.
            ("d", 5.0 + 2.0 * step + 2.0 * RT_TOLERANCE),
            ("e", f64::NAN),
            ("f", f64::NAN),
        ] {
            let mut s = named(name, rt);
            continued(agg.consume_spectrum(&mut s));
            // The caller's spectrum is not changed.
            assert_eq!(s.name, name);
        }
        assert_eq!(agg.pending(), 1);
        continued(agg.finish());
    }
    assert_eq!(group_sizes, [3, 1, 1, 1]);
    let names: Vec<_> = storage.data().spectra.iter().map(|s| s.name.clone()).collect();
    assert_eq!(names, ["a", "d", "e", "f"]);
}

#[test]
fn native_aggregating_errors_leave_the_group_and_break_stops() {
    let mut storage = MSDataStoringConsumer::new();
    let mut fail = true;
    {
        let mut agg = MSDataAggregatingConsumer::new(&mut storage, |spectra: &[MSSpectrum]| {
            if fail {
                fail = false;
                return Err(Error::InvalidValue("merge".into()));
            }
            add_up_on_grid(spectra)
        });
        continued(agg.consume_spectrum(&mut named("a", 1.0)));
        // The merge of group "a" fails when "b" arrives: "b" is not taken.
        assert!(agg.consume_spectrum(&mut named("b", 2.0)).is_err());
        assert_eq!(agg.pending(), 1);
        continued(agg.consume_spectrum(&mut named("b", 2.0)));
        continued(agg.finish());
    }
    let names: Vec<_> = storage.data().spectra.iter().map(|s| s.name.clone()).collect();
    assert_eq!(names, ["a", "b"]);

    // A next consumer that stops: the group it received was delivered, the
    // spectrum that triggered the flush is not kept, and later records are
    // refused with Break.
    let mut received = 0usize;
    let mut stopper = MSDataTransformingConsumer::new();
    stopper.set_spectra_processing_func(|_: &mut MSSpectrum| {
        received += 1;
        Ok(())
    });
    struct StopAfterFirst<'a, 'b>(&'a mut MSDataTransformingConsumer<'b>);
    impl MSDataConsumer for StopAfterFirst<'_, '_> {
        fn set_expected_size(&mut self, s: usize, c: usize) -> openms::Result<()> {
            self.0.set_expected_size(s, c)
        }
        fn set_experimental_settings(&mut self, e: &ExperimentalSettings) -> openms::Result<()> {
            self.0.set_experimental_settings(e)
        }
        fn consume_spectrum(&mut self, s: &mut MSSpectrum) -> openms::Result<ControlFlow<()>> {
            let _ = self.0.consume_spectrum(s)?;
            Ok(ControlFlow::Break(()))
        }
        fn consume_chromatogram(
            &mut self,
            c: &mut MSChromatogram,
        ) -> openms::Result<ControlFlow<()>> {
            self.0.consume_chromatogram(c)
        }
    }
    {
        let mut next = StopAfterFirst(&mut stopper);
        let mut agg = MSDataAggregatingConsumer::new(&mut next, add_up_on_grid);
        continued(agg.consume_spectrum(&mut named("a", 1.0)));
        assert_eq!(
            agg.consume_spectrum(&mut named("b", 2.0)).unwrap(),
            ControlFlow::Break(())
        );
        assert_eq!(agg.pending(), 0);
        assert_eq!(
            agg.consume_spectrum(&mut named("c", 3.0)).unwrap(),
            ControlFlow::Break(())
        );
        assert_eq!(
            agg.consume_chromatogram(&mut MSChromatogram::default())
                .unwrap(),
            ControlFlow::Break(())
        );
        continued(agg.finish());
    }
    drop(stopper);
    assert_eq!(received, 1);
}

#[test]
fn native_aggregating_chromatograms_pass_through_immediately() {
    let mut storage = MSDataStoringConsumer::new();
    {
        let mut agg = MSDataAggregatingConsumer::new(&mut storage, add_up_on_grid);
        continued(agg.consume_spectrum(&mut named("a", 1.0)));
        let mut c = MSChromatogram {
            native_id: "chrom".into(),
            peaks: vec![ChromatogramPeak { rt: 1.0, intensity: 2.0 }],
            ..MSChromatogram::default()
        };
        continued(agg.consume_chromatogram(&mut c));
        assert_eq!(agg.pending(), 1);
        continued(agg.finish());
    }
    // The chromatogram arrived before the spectrum group was flushed.
    assert_eq!(storage.data().chromatograms.len(), 1);
    assert_eq!(storage.data().spectra.len(), 1);
}

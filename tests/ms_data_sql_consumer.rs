// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FORMAT/DATAACCESS/MSDataSqlConsumer`.
//!
//! The source has no class test for this consumer, so every expectation here is
//! independently derived (tier 4) from its implementation at core `bc9cc12`
//! (`MSDataSqlConsumer.cpp`), cited at each use. The central invariant: the
//! consumer is the streaming form of a whole-experiment sqMass write, so what
//! it stores for an experiment reads back record for record as what
//! `MzMLSqliteHandler::write_experiment` stores for the same experiment, for
//! every buffer size. The input is the unmodified upstream
//! `MzMLSqliteHandler_1.mzML`, already recorded in
//! `tests/data/mzml_sqlite_handler_provenance.json`.

#![cfg(feature = "sqmass")]

use openms::Error;
use openms::format::ms_data_sql_consumer::{
    DEFAULT_BUFFER_SIZE, DEFAULT_LINEAR_MASS_ACC, MSDataSqlConsumer,
};
use openms::format::mzml_sqlite_handler::{LossPolicy, MzMLSqliteHandler};
use openms::identification::PeptideIdentification;
use openms::interfaces::MSDataConsumer;
use openms::kernel::{MSExperiment, MSSpectrum, Peak1D};
use openms::metadata::Polarity;
use openms::system::file::TempDir;
use std::ops::ControlFlow;
use std::path::{Path, PathBuf};

fn original() -> MSExperiment {
    openms::format::mzml::load(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/mzml_sqlite_handler/MzMLSqliteHandler_1.mzML"),
    )
    .unwrap()
}

fn temporary(name: &str) -> (TempDir, PathBuf) {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join(name);
    (dir, path)
}

fn read_back(path: &Path) -> MSExperiment {
    MzMLSqliteHandler::new(path, 0)
        .read_experiment(false)
        .unwrap()
}

/// What a whole-experiment write stores for `input`.
fn reference(input: &MSExperiment, full_meta: bool) -> MSExperiment {
    let (_dir, path) = temporary("reference.sqMass");
    let mut handler = MzMLSqliteHandler::new(&path, 7);
    handler
        .set_config(full_meta, false, DEFAULT_LINEAR_MASS_ACC, 500)
        .unwrap();
    handler.set_loss_policy(LossPolicy::Source);
    handler.create_tables().unwrap();
    handler.write_experiment(input).unwrap();
    read_back(&path)
}

fn consume_all(consumer: &mut MSDataSqlConsumer, input: &MSExperiment) {
    let mut records = input.clone();
    for spectrum in &mut records.spectra {
        assert_eq!(
            consumer.consume_spectrum(spectrum).unwrap(),
            ControlFlow::Continue(())
        );
    }
    for chromatogram in &mut records.chromatograms {
        assert_eq!(
            consumer.consume_chromatogram(chromatogram).unwrap(),
            ControlFlow::Continue(())
        );
    }
}

fn continued(flow: openms::Result<ControlFlow<()>>) {
    assert_eq!(flow.unwrap(), ControlFlow::Continue(()));
}

#[test]
fn defaults_are_the_source_constructor_defaults() {
    // MSDataSqlConsumer.h: flush_after = 500, full_meta = true,
    // lossy_compression = false, linear_mass_acc = 1e-4.
    assert_eq!(DEFAULT_BUFFER_SIZE, 500);
    assert_eq!(DEFAULT_LINEAR_MASS_ACC, 1e-4);
    let (_dir, path) = temporary("defaults.sqMass");
    let consumer = MSDataSqlConsumer::new(&path, 3).unwrap();
    assert_eq!(consumer.flush_after(), 500);
    assert!(consumer.full_meta());
    assert!(!consumer.wrote_any_run());
    assert_eq!(consumer.filename(), path.as_path());
    let config = consumer.handler().config();
    assert!(config.write_full_meta);
    assert!(!config.use_lossy_compression);
    assert_eq!(config.linear_abs_mass_acc, 1e-4);
    assert_eq!(config.sql_batch_size, 500);
    // The constructor creates the tables (MSDataSqlConsumer.cpp:26).
    assert!(path.exists());
    consumer.finish().unwrap();
}

#[test]
fn streamed_records_read_back_as_a_whole_experiment_write() {
    let input = original();
    for full_meta in [true, false] {
        let expected = reference(&input, full_meta);
        for buffer in [1, 2, DEFAULT_BUFFER_SIZE] {
            let (_dir, path) = temporary("streamed.sqMass");
            let mut consumer = MSDataSqlConsumer::with_config(
                &path,
                7,
                buffer,
                full_meta,
                false,
                DEFAULT_LINEAR_MASS_ACC,
            )
            .unwrap();
            consumer.set_loss_policy(LossPolicy::Source);
            consume_all(&mut consumer, &input);
            consumer.finish().unwrap();
            let actual = read_back(&path);
            assert_eq!(
                actual.spectra, expected.spectra,
                "full_meta {full_meta}, buffer {buffer}"
            );
            assert_eq!(
                actual.chromatograms, expected.chromatograms,
                "full_meta {full_meta}, buffer {buffer}"
            );
            assert_eq!(MzMLSqliteHandler::new(&path, 0).run_id().unwrap(), 7);
        }
    }
}

#[test]
fn consuming_clears_the_callers_peaks_and_keeps_its_metadata() {
    // MSDataSqlConsumer.cpp:74-77: push a copy, then `s.clear(false)`, and the
    // cleared record is what `peak_meta_` keeps.
    let input = original();
    let (_dir, path) = temporary("cleared.sqMass");
    let mut consumer = MSDataSqlConsumer::new(&path, 1).unwrap();
    consumer.set_loss_policy(LossPolicy::Source);
    let mut spectrum = input.spectra[0].clone();
    assert!(!spectrum.peaks.is_empty());
    continued(consumer.consume_spectrum(&mut spectrum));
    assert!(spectrum.peaks.is_empty());
    assert_eq!(spectrum.native_id, input.spectra[0].native_id);
    assert_eq!(spectrum.rt, input.spectra[0].rt);
    assert_eq!(consumer.buffered_spectra().len(), 1);
    assert_eq!(consumer.buffered_spectra()[0].peaks, input.spectra[0].peaks);
    assert_eq!(consumer.peak_meta().spectra, [spectrum.clone()]);

    let mut chromatogram = input.chromatograms[0].clone();
    assert!(!chromatogram.peaks.is_empty());
    continued(consumer.consume_chromatogram(&mut chromatogram));
    assert!(chromatogram.peaks.is_empty());
    assert_eq!(chromatogram.native_id, input.chromatograms[0].native_id);
    assert_eq!(consumer.peak_meta().chromatograms, [chromatogram.clone()]);
    consumer.finish().unwrap();

    // Without full metadata nothing is kept for the snapshot
    // (MSDataSqlConsumer.cpp:75, `if (full_meta_)`).
    let (_dir, path) = temporary("lean.sqMass");
    let mut consumer =
        MSDataSqlConsumer::with_config(&path, 1, 10, false, false, DEFAULT_LINEAR_MASS_ACC)
            .unwrap();
    consumer.set_loss_policy(LossPolicy::Source);
    let mut spectrum = input.spectra[0].clone();
    continued(consumer.consume_spectrum(&mut spectrum));
    assert!(consumer.peak_meta().spectra.is_empty());
    consumer.finish().unwrap();
}

#[test]
fn a_full_buffer_is_written_at_once() {
    // MSDataSqlConsumer.cpp:78-81: flush when `size() >= flush_after_`.
    let input = original();
    assert_eq!(input.spectra.len(), 2);
    let (_dir, path) = temporary("buffer.sqMass");
    let mut consumer =
        MSDataSqlConsumer::with_config(&path, 1, 2, true, false, DEFAULT_LINEAR_MASS_ACC).unwrap();
    consumer.set_loss_policy(LossPolicy::Source);
    let reader = MzMLSqliteHandler::new(&path, 0);
    let mut first = input.spectra[0].clone();
    continued(consumer.consume_spectrum(&mut first));
    assert_eq!(reader.nr_spectra().unwrap(), 0);
    assert_eq!(consumer.buffered_spectra().len(), 1);
    let mut second = input.spectra[1].clone();
    continued(consumer.consume_spectrum(&mut second));
    assert_eq!(reader.nr_spectra().unwrap(), 2);
    assert!(consumer.buffered_spectra().is_empty());

    // An explicit flush writes a partial buffer; spectra and chromatograms
    // have separate buffers (MSDataSqlConsumer.cpp:50-64).
    let mut chromatogram = input.chromatograms[0].clone();
    continued(consumer.consume_chromatogram(&mut chromatogram));
    assert_eq!(reader.nr_chromatograms().unwrap(), 0);
    consumer.flush().unwrap();
    assert_eq!(reader.nr_chromatograms().unwrap(), 1);
    assert!(consumer.buffered_chromatograms().is_empty());
    consumer.finish().unwrap();
    assert_eq!(reader.nr_spectra().unwrap(), 2);
    assert_eq!(reader.nr_chromatograms().unwrap(), 1);
}

#[test]
fn sizes_and_settings_are_ignored() {
    // MSDataSqlConsumer.cpp:98-100: both overrides are empty.
    let (_dir, path) = temporary("ignored.sqMass");
    let mut consumer = MSDataSqlConsumer::new(&path, 1).unwrap();
    consumer.set_expected_size(10, 10).unwrap();
    let mut settings = original().settings;
    settings.comment = "not stored".into();
    consumer.set_experimental_settings(&settings).unwrap();
    assert_eq!(consumer.peak_meta().settings.comment, "");
    consumer.finish().unwrap();
}

#[test]
fn finish_writes_the_run_once_and_names_the_sqmass_file() {
    // MSDataSqlConsumer.cpp:29-39: the destructor sets the loaded file path to
    // the sqMass file name and writes the run unless one was written.
    let (_dir, path) = temporary("named.sqMass");
    let consumer =
        MSDataSqlConsumer::with_config(&path, 9, 5, false, false, DEFAULT_LINEAR_MASS_ACC).unwrap();
    consumer.finish().unwrap();
    let read = read_back(&path);
    assert_eq!(
        read.settings.document.loaded_file_path,
        path.to_string_lossy()
    );
    assert_eq!(MzMLSqliteHandler::new(&path, 0).run_id().unwrap(), 9);
}

#[test]
fn add_run_replaces_the_run_written_at_the_end() {
    // MSDataSqlConsumer.cpp:42-49: `addRun` sets the handler's run id, writes a
    // run naming `filename`, and the destructor then writes none.
    let (_dir, path) = temporary("added.sqMass");
    let mut consumer =
        MSDataSqlConsumer::with_config(&path, 9, 5, false, false, DEFAULT_LINEAR_MASS_ACC).unwrap();
    consumer.add_run("other.mzML", 11).unwrap();
    assert!(consumer.wrote_any_run());
    let mut spectrum = MSSpectrum {
        native_id: "scan=1".into(),
        rt: 1.0,
        peaks: vec![Peak1D::new(100.0, 1.0)],
        ..MSSpectrum::default()
    };
    // SQL-only storage keeps positive or negative polarity only, so the
    // default loss policy refuses the default (unknown) one.
    spectrum.instrument_settings.polarity = Polarity::Positive;
    continued(consumer.consume_spectrum(&mut spectrum));
    consumer.finish().unwrap();
    let reader = MzMLSqliteHandler::new(&path, 0);
    assert_eq!(reader.run_id().unwrap(), 11);
    assert_eq!(reader.nr_spectra().unwrap(), 1);
    assert_eq!(
        read_back(&path).settings.document.loaded_file_path,
        "other.mzML"
    );
}

#[test]
fn set_run_id_masks_bit_63() {
    // MSDataSqlConsumer.cpp:51-54 forwards to the handler, which clears the
    // sign bit of the SQL identifier.
    let (_dir, path) = temporary("masked.sqMass");
    let mut consumer =
        MSDataSqlConsumer::with_config(&path, 0, 5, false, false, DEFAULT_LINEAR_MASS_ACC).unwrap();
    consumer.set_run_id(u64::MAX);
    consumer.finish().unwrap();
    assert_eq!(
        MzMLSqliteHandler::new(&path, 0).run_id().unwrap(),
        i64::MAX as u64
    );
}

#[test]
fn dropping_without_finish_still_ends_the_file() {
    // The source's destructor flushes and writes the run; a drop does the same.
    let input = original();
    let (_dir, path) = temporary("dropped.sqMass");
    {
        let mut consumer = MSDataSqlConsumer::new(&path, 5).unwrap();
        consumer.set_loss_policy(LossPolicy::Source);
        consume_all(&mut consumer, &input);
    }
    let read = read_back(&path);
    assert_eq!(read.spectra.len(), input.spectra.len());
    assert_eq!(read.chromatograms.len(), input.chromatograms.len());
    assert_eq!(MzMLSqliteHandler::new(&path, 0).run_id().unwrap(), 5);
}

#[test]
fn native_invalid_configuration_is_refused() {
    // The source passes an `int` batch size and the accuracy to its handler
    // unchecked.
    let (_dir, path) = temporary("invalid.sqMass");
    for (buffer, accuracy) in [(0, 1e-4), (5, f64::NAN), (5, f64::INFINITY)] {
        let error =
            MSDataSqlConsumer::with_config(&path, 1, buffer, true, false, accuracy).unwrap_err();
        assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    }
}

#[test]
fn native_full_metadata_refusals_happen_before_the_record_is_taken() {
    let (_dir, path) = temporary("refused.sqMass");
    let mut consumer = MSDataSqlConsumer::new(&path, 1).unwrap();
    assert_eq!(consumer.loss_policy(), LossPolicy::Reject);

    // An empty native identifier cannot be matched to its snapshot record.
    let mut spectrum = MSSpectrum {
        rt: 1.0,
        peaks: vec![Peak1D::new(100.0, 1.0)],
        ..MSSpectrum::default()
    };
    let error = consumer.consume_spectrum(&mut spectrum).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    assert_eq!(spectrum.peaks.len(), 1);
    assert!(consumer.buffered_spectra().is_empty());
    assert!(consumer.peak_meta().spectra.is_empty());

    // Peptide identifications have no place in the mzML snapshot, under
    // either policy.
    spectrum.native_id = "scan=1".into();
    spectrum
        .peptide_identifications
        .push(PeptideIdentification::default());
    for policy in [LossPolicy::Reject, LossPolicy::Source] {
        consumer.set_loss_policy(policy);
        let error = consumer.consume_spectrum(&mut spectrum).unwrap_err();
        assert!(matches!(error, Error::Unsupported(_)), "{error:?}");
        assert_eq!(spectrum.peaks.len(), 1);
    }
    assert!(consumer.buffered_spectra().is_empty());

    // Auxiliary arrays are refused under the default policy only.
    spectrum.peptide_identifications.clear();
    spectrum.float_data_arrays.push(Default::default());
    consumer.set_loss_policy(LossPolicy::Reject);
    let error = consumer.consume_spectrum(&mut spectrum).unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error:?}");
    consumer.set_loss_policy(LossPolicy::Source);
    continued(consumer.consume_spectrum(&mut spectrum));
    assert_eq!(consumer.buffered_spectra().len(), 1);
    consumer.finish().unwrap();
}

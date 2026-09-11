// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml")]
use openms::{
    MSChromatogram, MSExperiment, MSSpectrum, Result,
    format::mzml::{self, LoadOptions, ReadOptions, TransformOptions},
    interfaces::MSDataConsumer,
    metadata::ExperimentalSettings,
};
use std::{
    io::{Cursor, Write},
    ops::ControlFlow,
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
fn cv(id: &str, value: &str) -> String {
    format!("<cvParam cvRef=\"MS\" accession=\"{id}\" name=\"fixture\" value=\"{value}\"/>")
}
fn array(role: &str, compression: &str, payload: &str, declared: usize) -> String {
    format!(
        "<binaryDataArray encodedLength=\"{declared}\">{}{}{}<binary>{payload}</binary></binaryDataArray>",
        cv("MS:1000523", ""),
        cv(compression, ""),
        cv(role, "")
    )
}
fn document(compression: &str, payload: &str, declared: usize) -> String {
    let intensity = array("MS:1000515", compression, payload, declared);
    let mz = array("MS:1000514", compression, payload, declared);
    let time = array("MS:1000595", compression, payload, declared).replace("accession=\"MS:1000595\"", "accession=\"MS:1000595\" unitAccession=\"UO:0000010\" unitCvRef=\"UO\" unitName=\"second\"");
    format!(
        "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\"><fileDescription><fileContent/></fileDescription><softwareList count=\"1\"><software id=\"sw\" version=\"1\">{}</software></softwareList><instrumentConfigurationList count=\"1\"><instrumentConfiguration id=\"ic\"/></instrumentConfigurationList><dataProcessingList count=\"1\"><dataProcessing id=\"dp\"><processingMethod order=\"0\" softwareRef=\"sw\">{}</processingMethod></dataProcessing></dataProcessingList><run id=\"r\" defaultInstrumentConfigurationRef=\"ic\"><spectrumList count=\"1\" defaultDataProcessingRef=\"dp\"><spectrum id=\"scan=1\" index=\"0\" defaultArrayLength=\"4\">{}<scanList count=\"1\"><scan><cvParam accession=\"MS:1000016\" name=\"scan start time\" value=\"1\" unitAccession=\"UO:0000010\"/></scan></scanList><binaryDataArrayList count=\"2\">{mz}{intensity}</binaryDataArrayList></spectrum></spectrumList><chromatogramList count=\"1\" defaultDataProcessingRef=\"dp\"><chromatogram id=\"trace\" index=\"0\" defaultArrayLength=\"4\"><binaryDataArrayList count=\"2\">{time}{intensity}</binaryDataArrayList></chromatogram></chromatogramList></run></mzML>",
        cv("MS:1000799", "test"),
        cv("MS:1000544", ""),
        cv("MS:1000511", "1")
    )
}
fn cases() -> Vec<(&'static str, &'static str)> {
    include_str!("data/mzml_normalization/payloads.tsv")
        .lines()
        .skip(1)
        .map(|row| {
            let fields: Vec<_> = row.split('\t').collect();
            (fields[2], fields[3])
        })
        .collect()
}
fn options(skip: bool) -> LoadOptions {
    let mut options = LoadOptions::default();
    options.scientific.skip_xml_checks = skip;
    options
}
fn read(xml: &str, options: &LoadOptions) -> Result<MSExperiment> {
    mzml::read_with_load_options(Cursor::new(xml), options, &ReadOptions::default())
}
fn whitespace(text: &str, separator: &str) -> String {
    let mut result = separator.to_string();
    for chunk in text.as_bytes().chunks(4) {
        result.push_str(std::str::from_utf8(chunk).unwrap());
        result.push_str(separator);
    }
    result
}
#[derive(Default)]
struct Log {
    setup: Vec<&'static str>,
    spectra: Vec<MSSpectrum>,
    chromatograms: Vec<MSChromatogram>,
}
impl MSDataConsumer for Log {
    fn set_expected_size(&mut self, spectra: usize, chromatograms: usize) -> Result<()> {
        assert_eq!((spectra, chromatograms), (1, 1));
        self.setup.push("size");
        Ok(())
    }
    fn set_experimental_settings(&mut self, _: &ExperimentalSettings) -> Result<()> {
        self.setup.push("settings");
        Ok(())
    }
    fn consume_spectrum(&mut self, s: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        self.spectra.push(s.clone());
        Ok(ControlFlow::Continue(()))
    }
    fn consume_chromatogram(&mut self, c: &mut MSChromatogram) -> Result<ControlFlow<()>> {
        self.chromatograms.push(c.clone());
        Ok(ControlFlow::Continue(()))
    }
}
#[test]
fn contiguous_ordinary_zlib_and_all_numpress_modes_are_identical_under_both_settings() {
    assert_eq!(cases().len(), 8);
    for (compression, payload) in cases() {
        let xml = document(compression, payload, payload.len());
        let expected = read(&xml, &options(false)).unwrap();
        assert_eq!(read(&xml, &options(true)).unwrap(), expected);
        assert_eq!(expected.spectra[0].len(), 4);
        assert_eq!(expected.chromatograms[0].len(), 4);
        // Independent source input values, allowing the source lossy codecs.
        for (p, value) in expected.spectra[0]
            .peaks
            .iter()
            .zip([100.0, 200.0, 300.00005, 400.00010])
        {
            assert!((p.mz - value).abs() <= 0.0001 * value);
        }
    }
}
#[test]
fn only_four_xml_whitespaces_are_normalized_and_explicit_skip_rejects_them() {
    for (compression, payload) in cases() {
        let expected = read(
            &document(compression, payload, payload.len()),
            &options(false),
        )
        .unwrap();
        for separator in [" ", "\t", "\n", "\r", " \t\r\n"] {
            let padded = whitespace(payload, separator);
            let xml = document(compression, &padded, payload.len());
            assert_eq!(read(&xml, &options(false)).unwrap(), expected);
            assert_eq!(mzml::read(Cursor::new(&xml)).unwrap(), expected);
            assert!(read(&xml, &options(true)).is_err());
            // Even when declared length includes whitespace, strict Base64 may
            // not silently accept the source SIMD decoder's malformed alphabet.
            let advertised = document(compression, &padded, padded.len());
            assert!(read(&advertised, &options(true)).is_err());
        }
    }
}
#[test]
fn malformed_alphabet_padding_and_xml_stay_errors_cpp055() {
    for (compression, payload) in cases() {
        for prefix in ["!", "_", "-", "~", "=", "é"] {
            let bad = format!("{prefix}{}", &payload[1..]);
            let xml = document(compression, &bad, bad.len());
            for skip in [false, true] {
                assert!(read(&xml, &options(skip)).is_err());
            }
        }
        let xml = document(compression, payload, payload.len());
        for bad in [
            xml.replacen("</binary>", "</wrong>", 1),
            xml.replace("id=\"scan=1\"", "id=\"scan<1\""),
            xml.replace("<binary>", "<binary>\u{b}"),
            xml.replace("<binary>", "<binary>\u{c}"),
        ] {
            for skip in [false, true] {
                assert!(read(&bad, &options(skip)).is_err());
            }
        }
    }
}
#[test]
fn skipped_population_metadata_count_and_setup_never_decode_binary() {
    for skip in [false, true] {
        let xml = document("MS:1000576", "!!!!", 4);
        let mut load = options(skip);
        load.scientific.fill_data = false;
        let limits = ReadOptions {
            max_array_bytes: 0,
            max_total_array_bytes: 0,
            max_total_array_elements: 0,
            ..Default::default()
        };
        let result = mzml::read_with_load_options(Cursor::new(&xml), &load, &limits).unwrap();
        assert!(result.spectra[0].is_empty() && result.chromatograms[0].is_empty());
        let count =
            mzml::read_size_with_options(Cursor::new(&xml), &load.scientific, &limits).unwrap();
        assert_eq!((count.spectra, count.chromatograms), (1, 1));
        let mut consumer = Log::default();
        let mut transform = TransformOptions {
            load: load.clone(),
            read: limits,
            ..Default::default()
        };
        let report =
            mzml::transform_from(|| Ok(Cursor::new(&xml)), &mut consumer, &transform).unwrap();
        assert_eq!(consumer.setup, ["size", "settings"]);
        assert_eq!(
            (report.delivered.spectra, report.delivered.chromatograms),
            (1, 1)
        );
        assert!(consumer.spectra[0].is_empty() && consumer.chromatograms[0].is_empty());
        transform.load.scientific.fill_data = true;
        transform.read = ReadOptions::default();
        let mut failed = Log::default();
        assert!(mzml::transform_from(|| Ok(Cursor::new(&xml)), &mut failed, &transform).is_err());
        assert_eq!(failed.setup, ["size", "settings"]); // setup completed without decoding
        assert!(failed.spectra.is_empty() && failed.chromatograms.is_empty());
        load.scientific.fill_data = true;
        load.scientific.metadata_only = true;
        let result = mzml::read_with_load_options(Cursor::new(&xml), &load, &limits).unwrap();
        assert!(result.spectra.is_empty() && result.chromatograms.is_empty());
        transform.load = load;
        transform.read = limits;
        let mut metadata = Log::default();
        mzml::transform_from(|| Ok(Cursor::new(&xml)), &mut metadata, &transform).unwrap();
        assert_eq!(metadata.setup, ["size", "settings"]);
        assert!(metadata.spectra.is_empty() && metadata.chromatograms.is_empty());
    }
}
struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        static N: AtomicU64 = AtomicU64::new(0);
        let p = std::env::temp_dir().join(format!(
            "openms-normalize-{}-{}",
            std::process::id(),
            N.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&p).unwrap();
        Self(p)
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
#[test]
fn path_and_consumer_modes_share_normalization_and_atomic_retention() {
    let dir = Directory::new();
    for (compression, payload) in cases().into_iter().step_by(3) {
        let padded = whitespace(payload, "\n");
        for spaced in [false, true] {
            let text = if spaced { &padded } else { payload };
            let xml = document(compression, text, payload.len());
            for suffix in ["mzML", "mzML.gz", "mzML.bz2"] {
                let path = dir.0.join(format!("input.{suffix}"));
                if suffix.ends_with("gz") {
                    let mut writer =
                        flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
                    writer.write_all(xml.as_bytes()).unwrap();
                    std::fs::write(&path, writer.finish().unwrap()).unwrap();
                } else if suffix.ends_with("bz2") {
                    let mut writer =
                        bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
                    writer.write_all(xml.as_bytes()).unwrap();
                    std::fs::write(&path, writer.finish().unwrap()).unwrap();
                } else {
                    std::fs::write(&path, &xml).unwrap();
                }
                for skip in [false, true] {
                    let load = options(skip);
                    let parsed = mzml::load_with_options(&path, &load, &Default::default());
                    let mut log = Log::default();
                    let mut retained = MSExperiment {
                        spectra: vec![MSSpectrum {
                            name: "old".into(),
                            ..Default::default()
                        }],
                        ..Default::default()
                    };
                    let old = retained.clone();
                    let pointer = retained.spectra.as_ptr();
                    let transform = TransformOptions {
                        load: load.clone(),
                        ..Default::default()
                    };
                    let report = mzml::transform_into_with_options(
                        &path,
                        &mut log,
                        &mut retained,
                        &transform,
                    );
                    assert_eq!(log.setup, ["size", "settings"]);
                    if spaced && skip {
                        assert!(parsed.is_err() && report.is_err());
                        assert_eq!(retained, old);
                        assert_eq!(retained.spectra.as_ptr(), pointer);
                    } else {
                        let expected = parsed.unwrap();
                        assert_eq!(report.unwrap().delivered.spectra, 1);
                        assert_eq!(retained.spectra[1], expected.spectra[0]);
                        assert_eq!(log.chromatograms, expected.chromatograms);
                    }
                    assert_eq!(
                        mzml::load_size_with_options(&path, &load.scientific, &Default::default())
                            .unwrap()
                            .spectra,
                        1
                    );
                }
            }
        }
    }
}
#[test]
fn declared_lengths_counts_and_cumulative_resource_guards_remain_active() {
    let (compression, payload) = cases()[0];
    let xml = document(compression, payload, payload.len());
    for skip in [false, true] {
        let load = options(skip);
        for bad in [
            document(compression, payload, payload.len() - 1),
            xml.replace("defaultArrayLength=\"4\"", "defaultArrayLength=\"3\""),
        ] {
            assert!(read(&bad, &load).is_err());
        }
        for limits in [
            ReadOptions {
                max_xml_bytes: 10,
                ..Default::default()
            },
            ReadOptions {
                max_array_bytes: 1,
                ..Default::default()
            },
            ReadOptions {
                max_total_array_bytes: 64,
                ..Default::default()
            },
            ReadOptions {
                max_total_array_elements: 8,
                ..Default::default()
            },
            ReadOptions {
                max_total_arrays: 3,
                ..Default::default()
            },
        ] {
            assert!(mzml::read_with_load_options(Cursor::new(&xml), &load, &limits).is_err());
        }
        let mut no_work = load;
        no_work.max_selection_work = 0;
        assert!(read(&xml, &no_work).is_err());
        let mut destination = MSExperiment {
            spectra: vec![MSSpectrum {
                name: "old".into(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let before = destination.clone();
        let mut log = Log::default();
        let opt = TransformOptions {
            load: options(skip),
            read: ReadOptions {
                max_total_array_elements: 8,
                ..Default::default()
            },
            ..Default::default()
        };
        assert!(
            mzml::transform_from_into(|| Ok(Cursor::new(&xml)), &mut log, &mut destination, &opt)
                .is_err()
        );
        assert_eq!(destination, before);
    }
}

#[test]
fn original_transform_skip_flag_literals_retain_four_spectra_forty_peaks_tic350() {
    let source = include_str!("data/mzml_consumer/source.mzML");
    for skip in [false, true] {
        for retaining in [false, true] {
            let mut load = options(skip);
            load.scientific.fill_data = true;
            load.scientific.max_data_pool_size = 100;
            load.scientific.always_append_data = false;
            let transform = TransformOptions {
                skip_full_count: true,
                skip_first_pass: true,
                load,
                ..Default::default()
            };
            let mut log = Log::default();
            let mut destination = MSExperiment::default();
            let report = if retaining {
                mzml::transform_from_into(
                    || Ok(Cursor::new(source)),
                    &mut log,
                    &mut destination,
                    &transform,
                )
            } else {
                mzml::transform_from(|| Ok(Cursor::new(source)), &mut log, &transform)
            }
            .unwrap();
            assert!(log.setup.is_empty());
            assert_eq!(report.delivered.spectra, 4);
            assert_eq!(log.spectra.iter().map(|s| s.len()).sum::<usize>(), 40);
            assert_eq!(
                log.spectra
                    .iter()
                    .flat_map(|s| &s.peaks)
                    .map(|p| f64::from(p.intensity))
                    .sum::<f64>(),
                350.
            );
            assert_eq!(destination.spectra.len(), if retaining { 4 } else { 0 });
        }
    }
}

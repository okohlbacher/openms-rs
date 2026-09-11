// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml")]

use openms::format::mzml::{self, LoadOptions, ReadOptions, WriteOptions};
use openms::kernel::{DataArray, NumericRange};
use openms::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};
use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

const SOURCE: &[u8] = include_bytes!("data/mzml_load_source_projection.mzML");
struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "openms-mzml-paths-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
    fn write(&self, name: &str, bytes: &[u8]) -> PathBuf {
        let path = self.path(name);
        std::fs::write(&path, bytes).unwrap();
        path
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn signal() -> MSExperiment {
    let mut spectrum = MSSpectrum::from_peaks(vec![
        Peak1D::new(300., 3.),
        Peak1D::new(100., 1.),
        Peak1D::new(200., 2.),
    ]);
    spectrum.rt = 10.;
    spectrum.native_id = "scan=1".into();
    spectrum
        .integer_data_arrays
        .push(DataArray::new("label_ids", vec![30, 10, 20]));
    let mut chromatogram = MSChromatogram::from_peaks(vec![
        ChromatogramPeak::new(3., 30.),
        ChromatogramPeak::new(1., 10.),
        ChromatogramPeak::new(2., 20.),
    ]);
    chromatogram.native_id = "tic".into();
    MSExperiment {
        spectra: vec![spectrum],
        chromatograms: vec![chromatogram],
        ..Default::default()
    }
}
fn encoded() -> Vec<u8> {
    let mut xml = Vec::new();
    mzml::write(&mut xml, &signal()).unwrap();
    xml
}

#[test]
fn file_loading_applies_scientific_defaults_and_keeps_stream_order_compatibility() {
    let directory = Directory::new();
    let xml = encoded();
    let stream = mzml::read(xml.as_slice()).unwrap();
    assert_eq!(stream.spectra[0].peaks[0].mz, 300.);
    for suffix in [
        "sample.mzML",
        "sample.mzML.gz",
        "sample.unknown",
        "sample.idXML",
    ] {
        let path = directory.write(suffix, &xml);
        let loaded = mzml::load(&path).unwrap();
        assert_eq!(
            loaded.spectra[0]
                .peaks
                .iter()
                .map(|p| p.mz)
                .collect::<Vec<_>>(),
            [100., 200., 300.]
        );
        assert_eq!(loaded.spectra[0].integer_data_arrays[0].data, [10, 20, 30]);
        assert_eq!(
            loaded.chromatograms[0]
                .peaks
                .iter()
                .map(|p| p.rt)
                .collect::<Vec<_>>(),
            [1., 2., 3.]
        );
    }
}

#[test]
fn source_projection_filters_through_paths_and_atomic_replacement() {
    let directory = Directory::new();
    let path = directory.write("source.mzML", SOURCE);
    let mut options = LoadOptions::default();
    options.scientific.set_rt_range(NumericRange {
        min: 5.15,
        max: 5.35,
    });
    let filtered = mzml::load_with_options(&path, &options, &ReadOptions::default()).unwrap();
    assert_eq!(
        filtered.spectra.iter().map(|s| s.rt).collect::<Vec<_>>(),
        [5.2, 5.3]
    );
    let mut destination = signal();
    mzml::load_into_with_options(&path, &mut destination, &options, &ReadOptions::default())
        .unwrap();
    assert_eq!(destination, filtered);
    let before = destination.clone();
    std::fs::write(&path, b"<mzML><broken></mzML>").unwrap();
    assert!(mzml::load_into(&path, &mut destination).is_err());
    assert_eq!(destination, before);
    assert!(mzml::load_into(directory.path("missing.mzML"), &mut destination).is_err());
    assert_eq!(destination, before);
    std::fs::write(&path, SOURCE).unwrap();
    assert!(
        mzml::load_into_with_options(
            &path,
            &mut destination,
            &options,
            &ReadOptions {
                max_xml_bytes: 64,
                ..Default::default()
            }
        )
        .is_err()
    );
    assert_eq!(destination, before);
    options.scientific.metadata_only = true;
    // The historical scientific projection omits the required default processing
    // reference. Use a complete document for the newly supported header-only path.
    let header_path = directory.path("headers.mzML");
    mzml::store(&header_path, &signal()).unwrap();
    let full = mzml::load(&header_path).unwrap();
    let metadata =
        mzml::load_with_options(&header_path, &options, &ReadOptions::default()).unwrap();
    assert!(metadata.spectra.is_empty() && metadata.chromatograms.is_empty());
    assert_eq!(metadata.settings, full.settings);
    assert!(matches!(
        mzml::load_with_options(
            directory.path("missing.mzML"),
            &options,
            &ReadOptions::default()
        ),
        Err(openms::Error::Io(_))
    ));
    options.scientific.metadata_only = false;
    options.scientific.fill_data = false;
    assert!(matches!(
        mzml::load_with_options(
            directory.path("missing.mzML"),
            &options,
            &ReadOptions::default()
        ),
        Err(openms::Error::Unsupported(_))
    ));
}

#[test]
fn compressed_input_uses_magic_and_limits_decoded_xml() {
    let directory = Directory::new();
    let expected =
        mzml::read_with_load_options(SOURCE, &LoadOptions::default(), &ReadOptions::default())
            .unwrap();
    let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gzip.write_all(SOURCE).unwrap();
    let mut bzip = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
    bzip.write_all(SOURCE).unwrap();
    for (name, bytes) in [
        ("gzip.data", gzip.finish().unwrap()),
        ("bzip.mzML.gz", bzip.finish().unwrap()),
    ] {
        let path = directory.write(name, &bytes);
        assert_eq!(mzml::load(&path).unwrap(), expected);
        assert!(
            mzml::load_with_options(
                &path,
                &LoadOptions::default(),
                &ReadOptions {
                    max_xml_bytes: 64,
                    ..Default::default()
                }
            )
            .is_err()
        );
        std::fs::write(&path, &bytes[..bytes.len() / 2]).unwrap();
        let mut destination = signal();
        let before = destination.clone();
        assert!(mzml::load_into(&path, &mut destination).is_err());
        assert_eq!(destination, before);
    }
}

#[test]
fn outer_compression_and_binary_zlib_are_independent_and_failures_preserve_files() {
    let directory = Directory::new();
    for binary in [false, true] {
        for (name, magic) in [
            ("output.mzML", b"<?xml".as_slice()),
            ("output.gz", b"\x1f\x8b".as_slice()),
            ("output.bz2", b"BZh".as_slice()),
            ("output.extension", b"<?xml".as_slice()),
        ] {
            let path = directory.path(name);
            mzml::store_with_options(
                &path,
                &signal(),
                &WriteOptions {
                    zlib_compression: binary,
                },
            )
            .unwrap();
            let bytes = std::fs::read(&path).unwrap();
            assert!(bytes.starts_with(magic));
            let loaded = mzml::load(&path).unwrap();
            assert_eq!(loaded.spectra[0].integer_data_arrays[0].data, [10, 20, 30]);
            let mut invalid = signal();
            invalid.spectra[0].peaks[1].mz = f64::NAN;
            assert!(mzml::store(&path, &invalid).is_err());
            assert_eq!(std::fs::read(&path).unwrap(), bytes);
        }
    }
    let zip = directory.write("output.zip", b"existing");
    assert!(matches!(
        mzml::store(&zip, &signal()),
        Err(openms::Error::Unsupported(_))
    ));
    assert_eq!(std::fs::read(&zip).unwrap(), b"existing");
    assert_eq!(std::fs::read_dir(&directory.0).unwrap().count(), 5);
}

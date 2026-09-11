// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::{file_handler::FileHandler, file_types::FileType, mzml};
use openms::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D,
};
use openms::metadata::{
    Acquisition, ChromatogramType, DataProcessing, Product, ScanMode, ScanWindow,
};
use std::{io::Cursor, sync::Arc};

fn base() -> MSExperiment {
    let mut spectrum = MSSpectrum::from_peaks(vec![Peak1D::new(123., 45.)]);
    spectrum.native_id = "scan=1".into();
    spectrum.precursors.push(openms::Precursor::new(500., 2));
    let mut chromatogram = MSChromatogram::from_peaks(vec![ChromatogramPeak::new(10., 11.)]);
    chromatogram.native_id = "transition".into();
    chromatogram.product.mz = 250.;
    MSExperiment {
        spectra: vec![spectrum],
        chromatograms: vec![chromatogram],
        ..Default::default()
    }
}
fn rejected(experiment: &MSExperiment, owner: &str) {
    for mode in 0..4 {
        let mut bytes = b"preserved".to_vec();
        let result = match mode {
            0 => mzml::write(&mut bytes, experiment),
            1 => mzml::write_with_options(
                &mut bytes,
                experiment,
                &mzml::WriteOptions {
                    zlib_compression: true,
                },
            ),
            2 => mzml::write_with_numpress(
                &mut bytes,
                experiment,
                &mzml::NumpressWriteOptions::default(),
            )
            .map(|_| ()),
            _ => FileHandler::write_experiment(&mut bytes, experiment, FileType::MzMl),
        };
        assert!(
            matches!(result, Err(Error::Unsupported(_))),
            "{owner}, mode{mode}: {result:?}"
        );
        assert!(result.unwrap_err().to_string().contains(owner));
        assert_eq!(bytes, b"preserved");
    }
}
#[test]
fn default_source_settings_and_existing_precursor_product_transport_remain_supported() {
    // Source InstrumentSettings defaults: Unknown scan and polarity, no zoom or
    // windows; AcquisitionInfo defaults to empty method/vector/metadata.
    let experiment = base();
    assert!(!experiment.spectra[0].has_acquisition_settings());
    assert!(!experiment.chromatograms[0].has_acquisition_settings());
    assert_eq!(
        experiment.chromatograms[0].chromatogram_type,
        ChromatogramType::Mass
    );
    for mode in 0..2 {
        let mut bytes = Vec::new();
        if mode == 0 {
            mzml::write(&mut bytes, &experiment).unwrap();
        } else {
            mzml::write_with_numpress(
                &mut bytes,
                &experiment,
                &mzml::NumpressWriteOptions::default(),
            )
            .unwrap();
        }
        let read = mzml::read(Cursor::new(bytes)).unwrap();
        assert!(!read.spectra[0].has_acquisition_settings());
        assert!(!read.chromatograms[0].has_acquisition_settings());
        assert_eq!(read.spectra[0].precursors, experiment.spectra[0].precursors);
        assert_eq!(
            read.chromatograms[0].product,
            experiment.chromatograms[0].product
        );
        assert_eq!(read, experiment);
    }
    // Source equality considers -0 equal to zero; the native writer is stricter
    // because silently defaulting this attached field would erase its IEEE sign.
    for owner in ["spectrum", "chromatogram"] {
        let mut e = base();
        if owner == "spectrum" {
            e.spectra[0].source_file.size_mb = -0.0;
        } else {
            e.chromatograms[0].source_file.size_mb = -0.0;
        }
        rejected(&e, "source-file");
    }
}
#[test]
fn unrepresented_spectrum_acquisition_categories_reject_before_output() {
    let mut e = base();
    e.spectra[0]
        .instrument_settings
        .metadata
        .insert("empty".into(), "".into());
    rejected(&e, "spectrum");
    for arbitrary_cv in [false, true] {
        let mut e = base();
        if arbitrary_cv {
            e.spectra[0]
                .source_file
                .cv_terms
                .add(openms::metadata::CVTerm::new(
                    "MS:1000563",
                    "Thermo RAW format",
                    "MS",
                ))
                .unwrap();
        } else {
            e.spectra[0].source_file.size_mb = 1.;
        }
        rejected(&e, "source-file");
    }
}

#[test]
fn chromatogram_settings_and_unknown_type_remain_rejected() {
    for case in [0, 1, 4, 5] {
        let mut e = base();
        let c = &mut e.chromatograms[0];
        match case {
            0 => c.instrument_settings.scan_mode = ScanMode::MassSpectrum,
            1 => c.acquisition_info.acquisitions.push(Acquisition::default()),
            4 => {
                c.instrument_settings
                    .metadata
                    .insert("empty".into(), "".into());
            }
            _ => {
                c.acquisition_info
                    .metadata
                    .insert("empty".into(), "".into());
            }
        }
        rejected(&e, "chromatogram");
    }
    let mut e = base();
    e.chromatograms[0].chromatogram_type = ChromatogramType::Unknown;
    rejected(&e, "chromatogram");
}
#[test]
fn invalid_spectrum_products_error_before_output() {
    let mut e = base();
    e.spectra[0].products.push(Product {
        mz: f64::NAN,
        ..Default::default()
    });
    for numpress in [false, true] {
        let mut bytes = b"preserved".to_vec();
        let result = if numpress {
            mzml::write_with_numpress(&mut bytes, &e, &Default::default()).map(|_| ())
        } else {
            mzml::write(&mut bytes, &e)
        };
        assert!(result.is_err());
        assert_eq!(bytes, b"preserved");
    }
}
#[test]
fn guards_precede_deep_acquisition_validation_and_preserve_owned_shared_handles() {
    let mut e = base();
    let shared = Arc::new(DataProcessing::default());
    e.spectra[0].data_processing.push(Arc::clone(&shared));
    e.spectra[0]
        .instrument_settings
        .metadata
        .insert("unsupported".into(), "x".into());
    e.spectra[0]
        .instrument_settings
        .scan_windows
        .push(ScanWindow {
            begin: f64::NAN,
            end: 0.,
            ..Default::default()
        });
    assert!(e.validate().is_err());
    rejected(&e, "spectrum");
    assert_eq!(Arc::strong_count(&shared), 2);
    let mut e = base();
    e.chromatograms[0].source_file.size_mb = f32::NAN;
    assert!(e.validate().is_err());
    rejected(&e, "source-file");
}
#[test]
fn unrepresentable_array_description_types_fail_before_output() {
    for owner in 0..2 {
        let mut e = base();
        let mut array = DataArray::new("described", vec![1.]);
        array
            .metadata
            .insert("list".into(), vec!["x".to_string()].into());
        // A nonempty processing vector remains meaningful even when its record
        // is otherwise default; descriptions cannot disappear through XML.
        array
            .data_processing
            .push(Arc::new(DataProcessing::default()));
        if owner == 0 {
            e.spectra[0].float_data_arrays.push(array);
        } else {
            e.chromatograms[0].float_data_arrays.push(array);
        }
        rejected(&e, "Empty/list metadata");
    }
}
#[test]
fn path_and_dispatch_failures_preserve_existing_files() {
    let mut e = base();
    e.spectra[0].source_file.size_mb = 1.;
    for suffix in ["mzML", "mzML.gz", "mzML.bz2"] {
        let path = std::env::temp_dir().join(format!(
            "openms-acquisition-guard-{}.{}",
            std::process::id(),
            suffix
        ));
        std::fs::write(&path, b"existing").unwrap();
        assert!(matches!(mzml::store(&path, &e), Err(Error::Unsupported(_))));
        assert_eq!(std::fs::read(&path).unwrap(), b"existing");
        assert!(matches!(
            FileHandler::store_experiment(&path, &e, None),
            Err(Error::Unsupported(_))
        ));
        assert_eq!(std::fs::read(&path).unwrap(), b"existing");
        std::fs::remove_file(path).unwrap();
    }
}

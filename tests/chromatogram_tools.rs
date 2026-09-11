// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::kernel::{
    ChromatogramConversionLimits, ChromatogramPeak, ChromatogramTools, DataArray, MSChromatogram,
    MSExperiment, MSSpectrum, Peak1D, Precursor,
};
use openms::metadata::{Acquisition, ChromatogramType, DataProcessing, Product, ScanMode};
use std::sync::Arc;
fn spectrum(mode: ScanMode, rt: f64, precursors: &[f64], peaks: &[(f64, f32)]) -> MSSpectrum {
    let mut s = MSSpectrum {
        rt,
        precursors: precursors.iter().map(|&mz| Precursor::new(mz, 0)).collect(),
        peaks: peaks.iter().map(|&(m, i)| Peak1D::new(m, i)).collect(),
        ..Default::default()
    };
    s.instrument_settings.scan_mode = mode;
    s
}
fn chrom(precursor: f64, product: f64, kind: ChromatogramType, rt: &[f64]) -> MSChromatogram {
    MSChromatogram {
        precursor: Precursor::new(precursor, 0),
        product: Product {
            mz: product,
            ..Default::default()
        },
        chromatogram_type: kind,
        peaks: rt.iter().map(|&r| ChromatogramPeak::new(r, 0.)).collect(),
        ..Default::default()
    }
}
fn annotate(s: &mut MSSpectrum) {
    s.instrument_settings
        .metadata
        .insert("device".into(), "kept".into());
    s.instrument_settings.zoom_scan = true;
    s.acquisition_info.method_of_combination = "sum".into();
    s.acquisition_info.acquisitions.push(Acquisition {
        identifier: "scan x".into(),
        ..Default::default()
    });
    s.source_file.name = "input.raw".into();
    s.source_file
        .cv_terms
        .metadata
        .insert("origin".into(), 42_i64.into());
    s.products.push(Product {
        mz: 999.,
        ..Default::default()
    });
    s.data_processing.push(Arc::new(DataProcessing::default()));
}
#[test]
fn source_chromatograms_to_four_ms2_scans_literal() {
    let mut e = MSExperiment {
        chromatograms: vec![
            chrom(
                100.1,
                200.1,
                ChromatogramType::SelectedReactionMonitoring,
                &[0.1, 0.2],
            ),
            chrom(
                100.2,
                200.2,
                ChromatogramType::SelectedReactionMonitoring,
                &[0.2, 0.2],
            ),
        ],
        ..Default::default()
    };
    let r = ChromatogramTools::default()
        .convert_chromatograms_to_spectra(&mut e)
        .unwrap();
    assert_eq!(r.added_spectra, 4);
    assert_eq!(r.removed_chromatograms.len(), 2);
    assert!(e.chromatograms.is_empty());
    assert_eq!(
        e.spectra.iter().map(|s| s.rt).collect::<Vec<_>>(),
        [0.1, 0.2, 0.2, 0.2]
    );
    assert_eq!(e.spectra[0].peaks[0].mz, 200.1);
    assert_eq!(e.spectra[0].precursors[0].mz, 100.1);
    assert_eq!(e.spectra[2].products[0].mz, 200.2);
    for s in &e.spectra {
        assert_eq!(s.ms_level, 2);
        assert_eq!(
            s.instrument_settings.scan_mode,
            ScanMode::SelectedReactionMonitoring
        );
    }
}
#[test]
fn source_four_srm_scans_group_to_two_chromatograms_and_remove_only_srm() {
    let scans = vec![
        spectrum(
            ScanMode::SelectedReactionMonitoring,
            0.1,
            &[500.1],
            &[(100.1, 20_000_000.)],
        ),
        spectrum(
            ScanMode::SelectedReactionMonitoring,
            0.3,
            &[500.2],
            &[(100.2, 30_000_000.)],
        ),
        spectrum(
            ScanMode::SelectedReactionMonitoring,
            0.4,
            &[500.1],
            &[(100.1, 40_000_000.)],
        ),
        spectrum(
            ScanMode::SelectedReactionMonitoring,
            0.5,
            &[500.2],
            &[(100.2, 50_000_000.)],
        ),
        spectrum(ScanMode::MassSpectrum, -1., &[], &[]),
    ];
    for remove in [false, true] {
        let mut e = MSExperiment {
            spectra: scans.clone(),
            ..Default::default()
        };
        let report = ChromatogramTools::default()
            .convert_spectra_to_chromatograms(&mut e, remove, false)
            .unwrap();
        assert_eq!(report.added_chromatograms, 2);
        assert_eq!(e.spectra.len(), if remove { 1 } else { 5 });
        assert_eq!(report.removed_spectra.len(), if remove { 4 } else { 0 });
        assert_eq!(
            e.chromatograms[0].peaks,
            vec![
                ChromatogramPeak::new(0.1, 20_000_000.),
                ChromatogramPeak::new(0.4, 40_000_000.)
            ]
        );
        assert_eq!(e.chromatograms[1].precursor.mz, 500.2);
        assert_eq!(e.chromatograms[1].product.mz, 100.2);
    }
}
#[test]
fn conversion_copies_only_source_fields_and_returns_original_owned_payload() {
    let mut s = spectrum(
        ScanMode::SelectedReactionMonitoring,
        2.,
        &[500.],
        &[(200., 4.), (100., 3.)],
    );
    annotate(&mut s);
    s.native_id = "scan=source".into();
    s.metadata
        .insert("not propagated".into(), "still owned".into());
    s.precursors[0]
        .cv_terms
        .metadata
        .insert("precursor label".into(), "kept".into());
    s.float_data_arrays
        .push(DataArray::new("discard from output", vec![5., 6.]));
    let source = s.clone();
    let buffer = s.peaks.as_ptr();
    let mut e = MSExperiment {
        spectra: vec![s],
        ..Default::default()
    };
    let report = ChromatogramTools::default()
        .convert_spectra_to_chromatograms(&mut e, true, false)
        .unwrap();
    assert_eq!(report.removed_spectra[0], source);
    assert_eq!(report.removed_spectra[0].peaks.as_ptr(), buffer);
    for c in &e.chromatograms {
        assert_eq!(c.instrument_settings, source.instrument_settings);
        assert_eq!(c.acquisition_info, source.acquisition_info);
        assert_eq!(c.source_file, source.source_file);
        assert_eq!(c.precursor, source.precursors[0]);
        assert!(c.data_processing.is_empty());
        assert!(c.metadata.is_empty());
        assert!(c.float_data_arrays.is_empty());
        assert_eq!(c.native_id, "chromatogram=scan=source");
        assert_eq!(
            c.product,
            Product {
                mz: c.product.mz,
                ..Default::default()
            }
        );
    }
    assert_eq!(e.chromatograms[0].product.mz, 100.);
    assert_eq!(e.chromatograms[1].product.mz, 200.);
    e.chromatograms[0]
        .product
        .cv_terms
        .metadata
        .insert("product".into(), "kept".into());
    e.chromatograms[0].chromatogram_type = ChromatogramType::SelectedIonMonitoring;
    e.chromatograms[0]
        .data_processing
        .push(Arc::new(DataProcessing::default()));
    let c = e.chromatograms[0].clone();
    let prior = e.chromatograms[0].peaks.as_ptr();
    let report = ChromatogramTools::default()
        .convert_chromatograms_to_spectra(&mut e)
        .unwrap();
    assert_eq!(report.removed_chromatograms[0], c);
    assert_eq!(report.removed_chromatograms[0].peaks.as_ptr(), prior);
    assert_eq!(e.spectra[0].products[0], c.product);
    assert_eq!(e.spectra[0].source_file, c.source_file);
    assert_eq!(
        e.spectra[0].instrument_settings.scan_mode,
        ScanMode::SelectedIonMonitoring
    );
    assert!(e.spectra[0].native_id.is_empty());
    assert!(e.spectra[0].data_processing.is_empty());
}
#[test]
fn forced_xic_order_first_payload_repeated_keys_and_removal_quirks() {
    let mut first = spectrum(ScanMode::MassSpectrum, 3., &[], &[(200., 2.), (100., 1.)]);
    first.source_file.name = "first".into();
    let mut next = spectrum(
        ScanMode::MassSpectrum,
        1.,
        &[8., 9.],
        &[(100., 4.), (100., 5.)],
    );
    next.source_file.name = "later".into();
    let invalid_srm = spectrum(ScanMode::SelectedReactionMonitoring, 7., &[], &[]);
    let normal = spectrum(ScanMode::SelectedReactionMonitoring, 2., &[5.], &[(9., 6.)]);
    let mut e = MSExperiment {
        spectra: vec![first, next, invalid_srm, normal],
        chromatograms: vec![MSChromatogram {
            native_id: "prefix".into(),
            ..Default::default()
        }],
        ..Default::default()
    };
    let r = ChromatogramTools::default()
        .convert_spectra_to_chromatograms(&mut e, true, true)
        .unwrap();
    assert_eq!(r.added_chromatograms, 3);
    assert_eq!(r.removed_spectra.len(), 2);
    assert_eq!(e.spectra.len(), 2);
    assert_eq!(e.chromatograms[0].native_id, "prefix");
    let x = &e.chromatograms[1];
    assert_eq!(x.precursor.mz, 100.);
    assert_eq!(x.product.mz, 0.);
    assert_eq!(x.chromatogram_type, ChromatogramType::Mass);
    assert!(x.native_id.is_empty());
    assert_eq!(x.source_file.name, "first");
    assert_eq!(
        x.peaks,
        vec![
            ChromatogramPeak::new(3., 1.),
            ChromatogramPeak::new(1., 4.),
            ChromatogramPeak::new(1., 5.)
        ]
    );
    assert_eq!(
        x.precursor.cv_terms.metadata["description"]
            .as_str()
            .unwrap(),
        "XIC @ 100.0"
    );
    assert_eq!(e.chromatograms[2].precursor.mz, 200.);
    assert_eq!(
        e.chromatograms[3].chromatogram_type,
        ChromatogramType::SelectedReactionMonitoring
    );
}
#[test]
fn skipped_srm_is_removed_even_without_conversion_and_forced_one_precursor_is_srm() {
    let a = spectrum(
        ScanMode::SelectedReactionMonitoring,
        0.,
        &[1., 2.],
        &[(3., 4.)],
    );
    let b = spectrum(ScanMode::MassSpectrum, 5., &[9.], &[(10., 11.)]);
    let mut e = MSExperiment {
        spectra: vec![a, b],
        ..Default::default()
    };
    let report = ChromatogramTools::default()
        .convert_spectra_to_chromatograms(&mut e, true, false)
        .unwrap();
    assert_eq!(report.skipped_spectra, 1);
    assert_eq!(report.removed_spectra.len(), 1);
    assert!(e.chromatograms.is_empty());
    ChromatogramTools::default()
        .convert_spectra_to_chromatograms(&mut e, true, true)
        .unwrap();
    assert_eq!(e.spectra.len(), 1);
    assert_eq!(
        e.chromatograms[0].chromatogram_type,
        ChromatogramType::SelectedReactionMonitoring
    );
    assert_eq!(
        e.chromatograms[0].instrument_settings.scan_mode,
        ScanMode::MassSpectrum
    );
}
#[test]
fn exact_signed_zero_groups_retain_first_key_and_encounter_order() {
    let mut e = MSExperiment {
        spectra: vec![
            spectrum(
                ScanMode::SelectedReactionMonitoring,
                2.,
                &[-0.],
                &[(-0., 1.)],
            ),
            spectrum(
                ScanMode::SelectedReactionMonitoring,
                -1.,
                &[0.],
                &[(0., 2.)],
            ),
        ],
        ..Default::default()
    };
    ChromatogramTools::default()
        .convert_spectra_to_chromatograms(&mut e, false, false)
        .unwrap();
    assert_eq!(e.chromatograms.len(), 1);
    let c = &e.chromatograms[0];
    assert!(c.precursor.mz.is_sign_negative());
    assert!(c.product.mz.is_sign_negative());
    assert_eq!(c.peaks[1].rt, -1.);
}
#[test]
fn late_errors_and_resource_limits_are_atomic_and_skip_unused_payload() {
    let make = || MSExperiment {
        chromatograms: vec![chrom(1., 2., ChromatogramType::Mass, &[1., 2.])],
        ..Default::default()
    };
    for limit in [
        ChromatogramConversionLimits {
            max_output_records: 1,
            ..Default::default()
        },
        ChromatogramConversionLimits {
            max_points: 1,
            ..Default::default()
        },
        ChromatogramConversionLimits {
            max_work: 0,
            ..Default::default()
        },
        ChromatogramConversionLimits {
            max_bytes: 1,
            ..Default::default()
        },
    ] {
        let mut e = make();
        let old = e.clone();
        assert!(
            ChromatogramTools { limits: limit }
                .convert_chromatograms_to_spectra(&mut e)
                .is_err()
        );
        assert_eq!(e, old);
    }
    let mut e = make();
    e.chromatograms[0].peaks[1].rt = f64::NAN;
    let buffer = e.chromatograms.as_ptr();
    assert!(
        ChromatogramTools::default()
            .convert_chromatograms_to_spectra(&mut e)
            .is_err()
    );
    assert_eq!(e.chromatograms.as_ptr(), buffer);
    assert!(e.spectra.is_empty());
    let mut e = MSExperiment {
        spectra: vec![
            spectrum(ScanMode::SelectedReactionMonitoring, 1., &[5.], &[(2., 3.)]),
            spectrum(
                ScanMode::SelectedReactionMonitoring,
                2.,
                &[6.],
                &[(f64::NAN, 4.)],
            ),
        ],
        ..Default::default()
    };
    assert!(
        ChromatogramTools::default()
            .convert_spectra_to_chromatograms(&mut e, true, false)
            .is_err()
    );
    assert_eq!(e.spectra.len(), 2);
    assert!(e.chromatograms.is_empty());
    e.spectra[1].instrument_settings.scan_mode = ScanMode::MassSpectrum;
    ChromatogramTools::default()
        .convert_spectra_to_chromatograms(&mut e, false, false)
        .unwrap();
    assert_eq!(e.chromatograms.len(), 1);
}
#[test]
fn repeated_acquisition_copies_share_work_and_byte_limits() {
    let mut c = chrom(1., 2., ChromatogramType::Mass, &[1.]);
    c.source_file.name = "x".repeat(4096);
    let make = |c: MSChromatogram| MSExperiment {
        chromatograms: vec![c],
        ..Default::default()
    };
    let tools = ChromatogramTools {
        limits: ChromatogramConversionLimits {
            max_work: 7000,
            ..Default::default()
        },
    };
    assert!(
        tools
            .convert_chromatograms_to_spectra(&mut make(c.clone()))
            .is_ok()
    );
    c.peaks.push(ChromatogramPeak::new(2., 0.));
    let mut e = make(c);
    let original = e.clone();
    assert!(tools.convert_chromatograms_to_spectra(&mut e).is_err());
    assert_eq!(e, original);
}
#[test]
fn attached_settings_value_equality_shared_processing_clear_select_and_swap() {
    let mut s = spectrum(ScanMode::MassSpectrum, 3., &[], &[(1., 2.), (3., 4.)]);
    annotate(&mut s);
    s.float_data_arrays
        .push(DataArray::new("aligned", vec![5., 6.]));
    let clone = s.clone();
    assert!(Arc::ptr_eq(
        &clone.data_processing[0],
        &s.data_processing[0]
    ));
    let mut same = clone.clone();
    same.data_processing[0] = Arc::new(DataProcessing::default());
    assert_eq!(same, clone);
    s.select(&[1]).unwrap();
    assert_eq!(s.source_file, clone.source_file);
    assert_eq!(s.products, clone.products);
    let mut other = MSSpectrum::default();
    std::mem::swap(&mut s, &mut other);
    assert!(s.data_processing.is_empty());
    assert_eq!(other.data_processing, clone.data_processing);
    other.clear(false);
    assert!(other.peaks.is_empty());
    assert!(other.float_data_arrays.is_empty());
    assert_eq!(other.source_file, clone.source_file);
    other.clear(true);
    assert_eq!(other, MSSpectrum::default());
    let mut c = chrom(1., 2., ChromatogramType::SelectedIonMonitoring, &[3.]);
    c.source_file.name = "source".into();
    c.data_processing = clone.data_processing;
    c.clear(false);
    assert_eq!(c.chromatogram_type, ChromatogramType::SelectedIonMonitoring);
    assert!(!c.data_processing.is_empty());
    c.clear(true);
    assert_eq!(c, MSChromatogram::default());
}
#[test]
fn empty_conversions_preserve_unrelated_buffer_and_metadata() {
    let mut e = MSExperiment {
        spectra: vec![spectrum(ScanMode::MassSpectrum, 0., &[], &[(1., 2.)])],
        chromatograms: vec![MSChromatogram::default()],
        ..Default::default()
    };
    e.metadata.insert("experiment".into(), "unchanged".into());
    let ptr = e.spectra.as_ptr();
    let report = ChromatogramTools::default()
        .convert_chromatograms_to_spectra(&mut e)
        .unwrap();
    assert_eq!(e.spectra.as_ptr(), ptr);
    assert_eq!(report.removed_chromatograms.len(), 1);
    let report = ChromatogramTools::default()
        .convert_spectra_to_chromatograms(&mut e, false, false)
        .unwrap();
    assert_eq!(report.added_chromatograms, 0);
    assert_eq!(e.spectra.as_ptr(), ptr);
    assert_eq!(e.metadata["experiment"], "unchanged");
}

#[test]
fn two_dimensional_replacement_and_tic_leave_acquisition_ownership_explicit() {
    use openms::kernel::Peak2D;
    let mut s = spectrum(ScanMode::MassSpectrum, 2., &[], &[(10., 4.)]);
    annotate(&mut s);
    let old_source = s.source_file.clone();
    let handle = s.data_processing[0].clone();
    let mut e = MSExperiment {
        spectra: vec![s],
        chromatograms: vec![chrom(
            1.,
            2.,
            ChromatogramType::SelectedReactionMonitoring,
            &[3.],
        )],
        ..Default::default()
    };
    let tic = e.calculate_tic(1);
    assert!(!tic.has_acquisition_settings());
    assert_eq!(tic.chromatogram_type, ChromatogramType::Mass);
    let previous = e.set_2d_data(&[Peak2D::new(5., 20., 6.)]).unwrap();
    assert_eq!(previous.spectra[0].source_file, old_source);
    assert!(Arc::ptr_eq(
        &previous.spectra[0].data_processing[0],
        &handle
    ));
    assert_eq!(
        previous.chromatograms[0].chromatogram_type,
        ChromatogramType::SelectedReactionMonitoring
    );
    assert!(e.chromatograms.is_empty());
    assert!(!e.spectra[0].has_acquisition_settings());
}
#[test]
#[allow(clippy::excessive_precision)] // Preserve the source EMG fixture literals.
fn theoretical_append_and_emg_copy_preserve_acquisition_fields() {
    use openms::analysis::emg::EmgGradientDescent;
    use openms::chemistry::theoretical::TheoreticalSpectrumGenerator;
    let mut s = MSSpectrum::default();
    annotate(&mut s);
    let source = s.source_file.clone();
    let handle = s.data_processing[0].clone();
    TheoreticalSpectrumGenerator::default()
        .append_to(&mut s, &"PEPTIDE".parse().unwrap(), 1, 1, None)
        .unwrap();
    assert_eq!(s.source_file, source);
    assert!(Arc::ptr_eq(&s.data_processing[0], &handle));
    assert_eq!(s.products[0].mz, 999.);
    let mut input = spectrum(
        ScanMode::MassSpectrum,
        2.,
        &[],
        // Reuse the finite source cutoff example from tests/emg.rs; a broad
        // symmetric toy curve is outside this fitter's stable initial shape.
        &[
            (15.34253311, 3.48297429),
            (15.35624981, 15.54384613),
            (15.36995029, 50.31319046),
            (15.38366699, 151.8971405),
            (15.39736652, 411.25631714),
            (15.41156673, 946.44311523),
            (15.42574978, 1642.56152344),
            (15.44018364, 2118.89526367),
            (15.45436668, 2055.13647461),
            (15.46856689, 1665.13232422),
            (15.48274994, 1275.53015137),
            (15.49695015, 1009.70056152),
        ],
    );
    annotate(&mut input);
    let fitter = EmgGradientDescent {
        max_iterations: 1,
        compute_additional_points: false,
        ..Default::default()
    };
    let fit = fitter.fit_spectrum(&input, None, None).unwrap();
    assert_eq!(fit.spectrum.instrument_settings, input.instrument_settings);
    assert_eq!(fit.spectrum.source_file, input.source_file);
    assert_eq!(fit.spectrum.products, input.products);
    assert!(Arc::ptr_eq(
        &fit.spectrum.data_processing[0],
        &input.data_processing[0]
    ));
    let mut c = MSChromatogram {
        peaks: input
            .peaks
            .iter()
            .map(|p| ChromatogramPeak::new(p.mz, p.intensity))
            .collect(),
        instrument_settings: input.instrument_settings,
        source_file: input.source_file,
        data_processing: input.data_processing,
        chromatogram_type: ChromatogramType::SelectedIonMonitoring,
        ..Default::default()
    };
    c.acquisition_info.method_of_combination = "sum".into();
    let fit = fitter.fit_chromatogram(&c, None, None).unwrap();
    assert_eq!(fit.chromatogram.chromatogram_type, c.chromatogram_type);
    assert_eq!(fit.chromatogram.acquisition_info, c.acquisition_info);
    assert!(Arc::ptr_eq(
        &fit.chromatogram.data_processing[0],
        &c.data_processing[0]
    ));
}

#[test]
fn empty_descriptors_in_second_conversion_pass_consume_shared_work() {
    let tools = ChromatogramTools {
        limits: ChromatogramConversionLimits {
            max_work: 50,
            ..Default::default()
        },
    };
    let populated = chrom(1., 2., ChromatogramType::Mass, &[3.]);
    let mut small = MSExperiment {
        chromatograms: vec![populated.clone()],
        ..Default::default()
    };
    assert!(tools.convert_chromatograms_to_spectra(&mut small).is_ok());
    let mut many = MSExperiment {
        chromatograms: vec![MSChromatogram::default(); 31],
        ..Default::default()
    };
    many.chromatograms.push(populated);
    let old = many.clone();
    assert!(tools.convert_chromatograms_to_spectra(&mut many).is_err());
    assert_eq!(many, old);
}

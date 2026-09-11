use openms::kernel::{
    ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, SpectrumType,
};
use openms::metadata::{
    Acquisition, ChromatogramType, DataProcessing, Polarity, Product, ScanMode, ScanWindow,
};
use openms::processing::chromatogram::{ChromatogramSmoothing, PeakPickerChromatogram};
use openms::processing::deisotoping::{AveragineDeisotoper, Deisotoper};
use openms::processing::iterative::PeakPickerIterative;
use openms::processing::peak_picking::PeakPickerHiRes;
use openms::processing::smoothing::{GaussFilter, GaussianWidth, SavitzkyGolayFilter};
use openms::processing::window_mower::WindowMower;
use openms::processing::{Normalizer, SpectrumFilter, ThresholdMower};
use std::sync::Arc;

fn spectrum() -> MSSpectrum {
    let mut value = MSSpectrum {
        rt: 42.0,
        spectrum_type: SpectrumType::Profile,
        peaks: [0., 0., 1., 4., 10., 4., 1., 0., 0.]
            .iter()
            .enumerate()
            .map(|(i, &y)| Peak1D::new(100.0 + i as f64 * 0.01, y))
            .collect(),
        ..Default::default()
    };
    value.instrument_settings.scan_mode = ScanMode::MassSpectrum;
    value.instrument_settings.zoom_scan = true;
    value.instrument_settings.polarity = Polarity::Positive;
    value
        .instrument_settings
        .metadata
        .insert("device".into(), "kept".into());
    value.instrument_settings.scan_windows.push(ScanWindow {
        begin: 50.0,
        end: 200.0,
        ..Default::default()
    });
    value.acquisition_info.method_of_combination = "sum".into();
    value
        .acquisition_info
        .metadata
        .insert("combined".into(), "true".into());
    let mut scan = Acquisition {
        identifier: "scan=7".into(),
        ..Default::default()
    };
    scan.metadata.insert("injection".into(), 7_i64.into());
    value.acquisition_info.acquisitions.push(scan);
    value.source_file.name = "input.raw".into();
    value.source_file.path = "/runs".into();
    value.source_file.size_mb = 1.25;
    value.source_file.file_type = "raw".into();
    value
        .source_file
        .cv_terms
        .metadata
        .insert("origin".into(), "instrument".into());
    let mut product = Product {
        mz: 200.0,
        isolation_window_lower_offset: 0.4,
        isolation_window_upper_offset: 0.5,
        ..Default::default()
    };
    product
        .cv_terms
        .metadata
        .insert("transition".into(), 3_i64.into());
    value.products.push(product);
    let mut processing = DataProcessing::default();
    processing.software.name = "decode".into();
    processing.software.version = "1".into();
    processing
        .metadata
        .insert("reason".into(), "fixture".into());
    value.data_processing.push(Arc::new(processing));
    value
}
fn chromatogram(source: &MSSpectrum) -> MSChromatogram {
    MSChromatogram {
        peaks: source
            .peaks
            .iter()
            .enumerate()
            .map(|(i, p)| ChromatogramPeak::new(i as f64, p.intensity))
            .collect(),
        instrument_settings: source.instrument_settings.clone(),
        acquisition_info: source.acquisition_info.clone(),
        source_file: source.source_file.clone(),
        data_processing: source.data_processing.clone(),
        product: source.products[0].clone(),
        chromatogram_type: ChromatogramType::SelectedReactionMonitoring,
        ..Default::default()
    }
}
fn assert_spectrum(expected: &MSSpectrum, actual: &MSSpectrum) {
    assert_eq!(actual.instrument_settings, expected.instrument_settings);
    assert_eq!(actual.acquisition_info, expected.acquisition_info);
    assert_eq!(actual.source_file, expected.source_file);
    assert_eq!(actual.products, expected.products);
    assert_eq!(actual.data_processing.len(), expected.data_processing.len());
    assert!(
        actual
            .data_processing
            .iter()
            .zip(&expected.data_processing)
            .all(|(a, b)| Arc::ptr_eq(a, b))
    );
}
fn assert_chromatogram(expected: &MSChromatogram, actual: &MSChromatogram) {
    assert_eq!(actual.instrument_settings, expected.instrument_settings);
    assert_eq!(actual.acquisition_info, expected.acquisition_info);
    assert_eq!(actual.source_file, expected.source_file);
    assert_eq!(actual.product, expected.product);
    assert_eq!(actual.chromatogram_type, expected.chromatogram_type);
    assert_eq!(actual.data_processing.len(), expected.data_processing.len());
    assert!(
        actual
            .data_processing
            .iter()
            .zip(&expected.data_processing)
            .all(|(a, b)| Arc::ptr_eq(a, b))
    );
}
fn hires() -> PeakPickerHiRes {
    PeakPickerHiRes {
        signal_to_noise: 0.0,
        ..Default::default()
    }
}
fn iterative() -> PeakPickerIterative {
    PeakPickerIterative {
        signal_to_noise: 0.0,
        iterations: 1,
        ..Default::default()
    }
}

#[test]
fn direct_clone_based_spectrum_outputs_preserve_full_settings_and_shared_processing() {
    let input = spectrum();
    let mut outputs = vec![
        WindowMower::default().filtered_spectrum(&input).unwrap(),
        Deisotoper::default().deisotope(&input).unwrap().spectrum,
        AveragineDeisotoper::default()
            .deisotope(&input)
            .unwrap()
            .spectrum,
        hires().pick_spectrum(&input).unwrap().spectrum,
        iterative().pick_spectrum(&input).unwrap().picked.spectrum,
    ];
    for output in &outputs {
        assert_spectrum(&input, output);
    }
    outputs[0].source_file.name.push_str("-edited");
    outputs[0].acquisition_info.acquisitions[0]
        .identifier
        .push_str("-edited");
    assert_eq!(input.source_file.name, "input.raw");
    assert_eq!(outputs[1].source_file.name, "input.raw");
    assert_eq!(input.acquisition_info.acquisitions[0].identifier, "scan=7");
    let mut empty = input.clone();
    empty.peaks.clear();
    assert_spectrum(
        &empty,
        &Deisotoper::default().deisotope(&empty).unwrap().spectrum,
    );
    assert_spectrum(
        &empty,
        &AveragineDeisotoper::default()
            .deisotope(&empty)
            .unwrap()
            .spectrum,
    );
    assert_spectrum(
        &empty,
        &iterative().pick_spectrum(&empty).unwrap().picked.spectrum,
    );
}

#[test]
fn bundled_experiment_filters_preserve_each_records_acquisition_ownership() {
    let source = spectrum();
    let chrom = chromatogram(&source);
    let input = MSExperiment {
        spectra: vec![source.clone(), source.clone()],
        chromatograms: vec![chrom.clone()],
        ..Default::default()
    };
    let filters: Vec<Box<dyn SpectrumFilter>> = vec![
        Box::new(Normalizer::default()),
        Box::new(ThresholdMower::default()),
        Box::new(WindowMower::default()),
        Box::new(Deisotoper::default()),
        Box::new(AveragineDeisotoper::default()),
        Box::new(hires()),
        Box::new(iterative()),
        Box::new(GaussFilter::new(GaussianWidth::Absolute(0.05)).unwrap()),
        Box::new(SavitzkyGolayFilter::new(3, 1).unwrap()),
    ];
    for filter in filters {
        let mut output = input.clone();
        filter.filter_experiment(&mut output).unwrap();
        assert_eq!(output.spectra.len(), 2);
        assert_eq!(output.chromatograms.len(), 1);
        for value in &output.spectra {
            assert_spectrum(&source, value);
        }
        assert_chromatogram(&chrom, &output.chromatograms[0]);
    }
}

#[test]
fn nested_chromatogram_picker_preserves_both_returned_records_and_type() {
    let input = chromatogram(&spectrum());
    let picked = PeakPickerChromatogram {
        signal_to_noise: 0.0,
        seed_signal_to_noise: 0.0,
        smoothing: ChromatogramSmoothing::SavitzkyGolay {
            frame_length: 3,
            polynomial_order: 1,
        },
        ..Default::default()
    }
    .pick_chromatogram(&input)
    .unwrap();
    assert_chromatogram(&input, &picked.smoothed);
    assert_chromatogram(&input, &picked.picked.chromatogram);
    assert_chromatogram(
        &input,
        &hires().pick_chromatogram(&input).unwrap().chromatogram,
    );
}

#[test]
fn later_filter_error_preserves_original_record_ownership_and_metadata() {
    let good = spectrum();
    let mut bad = good.clone();
    bad.peaks[2].intensity = -1.0;
    let mut experiment = MSExperiment {
        spectra: vec![good, bad],
        ..Default::default()
    };
    let before = experiment.clone();
    let address = experiment.spectra.as_ptr();
    assert!(
        Deisotoper::default()
            .filter_experiment(&mut experiment)
            .is_err()
    );
    assert_eq!(experiment, before);
    assert_eq!(experiment.spectra.as_ptr(), address);
    assert_spectrum(&before.spectra[0], &experiment.spectra[0]);
    assert_spectrum(&before.spectra[1], &experiment.spectra[1]);
}

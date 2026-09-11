use openms::chemistry::{
    AASequence, NASequence, NucleicAcidSpectrumGenerator, TheoreticalSpectrumGenerator,
};
use openms::kernel::{DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, SummaryLimits};
use openms::metadata::{DataProcessing, MetaValue, MetaValueData, Unit};
use openms::processing::peak_picking::PeakPickerHiRes;
use std::sync::Arc;

fn describe<T>(array: &mut DataArray<T>) -> Arc<DataProcessing> {
    array.metadata.insert(
        "calibration".into(),
        MetaValue::new(MetaValueData::Float(2.))
            .unwrap()
            .with_unit(Unit::new("UO:0000010", "second", "UO").unwrap())
            .unwrap(),
    );
    let processing = Arc::new(DataProcessing::default());
    array.data_processing.push(processing.clone());
    processing
}

#[test]
fn array_description_defaults_equality_validation_and_shared_processing() {
    let mut a = DataArray::new("score", vec![1.]);
    assert!(!a.has_description_metadata());
    let processing = describe(&mut a);
    let b = a.clone();
    assert_eq!(a, b);
    assert!(Arc::ptr_eq(&a.data_processing[0], &processing));
    assert!(Arc::ptr_eq(&b.data_processing[0], &processing));
    a.metadata.insert("another".into(), "value".into());
    assert_ne!(a, b);
    a.validate_description().unwrap();
}

#[test]
fn spectrum_and_chromatogram_selection_preserve_full_descriptions() {
    let mut spectrum = MSSpectrum::from_peaks(vec![
        Peak1D::new(3., 3.),
        Peak1D::new(1., 1.),
        Peak1D::new(2., 2.),
    ]);
    let mut values = DataArray::new("identity", vec![3, 1, 2]);
    let processing = describe(&mut values);
    let metadata = values.metadata.clone();
    spectrum.integer_data_arrays.push(values);
    spectrum.sort_by_position().unwrap();
    spectrum.select(&[2, 0]).unwrap();
    assert_eq!(spectrum.integer_data_arrays[0].data, [3, 1]);
    assert_eq!(spectrum.integer_data_arrays[0].metadata, metadata);
    assert!(Arc::ptr_eq(
        &spectrum.integer_data_arrays[0].data_processing[0],
        &processing
    ));
    let mut chrom = MSChromatogram::from_peaks(vec![
        openms::ChromatogramPeak::new(3., 3.),
        openms::ChromatogramPeak::new(1., 1.),
    ]);
    chrom.integer_data_arrays = spectrum.integer_data_arrays;
    chrom.sort_by_position().unwrap();
    assert_eq!(chrom.integer_data_arrays[0].metadata, metadata);
}

#[test]
fn theoretical_append_preserves_existing_annotation_descriptions() {
    let mut spectrum = MSSpectrum::new();
    spectrum
        .string_data_arrays
        .push(DataArray::new("IonNames", Vec::new()));
    let processing = describe(&mut spectrum.string_data_arrays[0]);
    let before = spectrum.string_data_arrays[0].metadata.clone();
    TheoreticalSpectrumGenerator::default()
        .append_to(
            &mut spectrum,
            &AASequence::parse("PEPTIDE").unwrap(),
            1,
            1,
            None,
        )
        .unwrap();
    assert!(!spectrum.is_empty());
    assert_eq!(spectrum.string_data_arrays[0].metadata, before);
    assert!(Arc::ptr_eq(
        &spectrum.string_data_arrays[0].data_processing[0],
        &processing
    ));
    let mut rna = MSSpectrum::new();
    rna.string_data_arrays
        .push(DataArray::new("IonNames", Vec::new()));
    let processing = describe(&mut rna.string_data_arrays[0]);
    let before = rna.string_data_arrays[0].metadata.clone();
    NucleicAcidSpectrumGenerator {
        add_metainfo: true,
        ..Default::default()
    }
    .append_to(&mut rna, &NASequence::parse("ACGU").unwrap(), 1, 1)
    .unwrap();
    assert!(!rna.is_empty());
    assert_eq!(rna.string_data_arrays[0].metadata, before);
    assert!(Arc::ptr_eq(
        &rna.string_data_arrays[0].data_processing[0],
        &processing
    ));
}

#[test]
fn retained_peak_picker_mobility_array_keeps_its_description() {
    let mut input = MSSpectrum::from_peaks(
        [100., 100.03, 100.06, 100.07, 100.08]
            .into_iter()
            .zip([200., 250., 450., 250., 200.])
            .map(|(mz, intensity)| Peak1D::new(mz, intensity))
            .collect(),
    );
    let mut array = DataArray::new("Ion Mobility", vec![1., 1.1, 1.2, 1.3, 1.4]);
    let processing = describe(&mut array);
    let before = array.metadata.clone();
    input.float_data_arrays.push(array);
    let result = PeakPickerHiRes {
        allow_missing_flank: true,
        ..Default::default()
    }
    .pick_spectrum(&input)
    .unwrap();
    assert_eq!(result.spectrum.len(), 1);
    let a = &result.spectrum.float_data_arrays[0];
    assert_eq!(a.metadata, before);
    assert!(Arc::ptr_eq(&a.data_processing[0], &processing));
}

#[test]
fn bounded_array_clear_accounts_metadata_before_any_removal() {
    let mut experiment = MSExperiment::new();
    let mut s = MSSpectrum::new();
    let mut array = DataArray::new("", Vec::<i32>::new());
    describe(&mut array);
    s.integer_data_arrays.push(array);
    experiment.spectra.push(s);
    let before = experiment.clone();
    let limits = SummaryLimits {
        max_bytes: 1,
        ..Default::default()
    };
    assert!(
        experiment
            .clear_meta_data_arrays_with_limits(limits)
            .is_err()
    );
    assert_eq!(experiment, before);
    assert!(experiment.clear_meta_data_arrays().unwrap());
    assert!(experiment.spectra[0].integer_data_arrays.is_empty());
}

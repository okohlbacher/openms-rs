// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source literals and independent arithmetic; no C++ execution.
//! Provenance: data/spectrum_annotation_provenance.json.

use openms::chemistry::{
    AASequence, SpectrumAnnotator, TheoreticalIonSeries, TheoreticalSpectrumGenerator as Generator,
    ion_naming,
};
use openms::comparison::{SpectrumAlignment, Tolerance};
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::kernel::{DataArray, MSSpectrum, Peak1D, Precursor};
use openms::metadata::MetaValueData;

const INPUT: &str = include_str!("data/spectrum_annotation_source_peaks.tsv");
const SCALARS: &str = include_str!("data/spectrum_annotation_source_statistics.tsv");
const DERIVED: &str = include_str!("data/spectrum_annotation_derived_masses.tsv");

fn rows(text: &str) -> impl Iterator<Item = Vec<&str>> {
    text.lines().skip(1).map(|line| line.split('\t').collect())
}
fn f64_bits(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).unwrap())
}
fn f32_bits(text: &str) -> f32 {
    f32::from_bits(u32::from_str_radix(text, 16).unwrap())
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance}"
    );
}
fn hit(text: &str, charge: i32) -> PeptideHit {
    PeptideHit::new(1.0, 1, charge, AASequence::parse(text).unwrap()).unwrap()
}
fn identification(hit: PeptideHit) -> PeptideIdentification {
    PeptideIdentification {
        hits: vec![hit],
        ..Default::default()
    }
}
fn generator() -> Generator {
    Generator {
        add_metainfo: true,
        ..Default::default()
    }
}
fn alignment(tolerance: Tolerance) -> SpectrumAlignment {
    SpectrumAlignment {
        tolerance,
        ..Default::default()
    }
}
fn source_spectrum(unmatched: bool) -> MSSpectrum {
    MSSpectrum {
        ms_level: 2,
        peaks: rows(INPUT)
            .filter(|row| unmatched || row[6] == "matched")
            .map(|row| {
                assert_eq!(row[1].parse::<f64>().unwrap(), f64_bits(row[2]));
                assert_eq!(row[3].parse::<f32>().unwrap(), f32_bits(row[4]));
                Peak1D::new(f64_bits(row[2]), f32_bits(row[4]))
            })
            .collect(),
        ..Default::default()
    }
}
fn source_labels() -> Vec<String> {
    let mut rows: Vec<_> = rows(INPUT).filter(|r| r[6] == "matched").collect();
    rows.sort_by(|a, b| f64_bits(a[2]).total_cmp(&f64_bits(b[2])));
    rows.into_iter().map(|r| r[5].to_owned()).collect()
}
fn number(id: &PeptideIdentification, key: &str) -> f64 {
    id.hits[0].metadata[key].as_f64().unwrap()
}
fn text<'a>(id: &'a PeptideIdentification, key: &str) -> &'a str {
    id.hits[0].metadata[key].as_str().unwrap()
}

#[test]
fn source_eleven_peak_annotations_and_thirteen_peak_unmatched_variant() {
    let annotator = SpectrumAnnotator::default();
    let mut spectrum = source_spectrum(false);
    spectrum
        .float_data_arrays
        .push(DataArray::new("old float", vec![1.0; 11]));
    spectrum
        .integer_data_arrays
        .push(DataArray::new("old int", vec![7; 11]));
    spectrum
        .string_data_arrays
        .push(DataArray::new("old text", vec!["old".into(); 11]));
    spectrum
        .metadata
        .insert("retained".into(), "metadata".into());
    let align = alignment(Tolerance::Absolute(0.1));
    annotator
        .annotate_matches(&mut spectrum, &hit("IFSQVGK", 2), &generator(), &align)
        .unwrap();
    assert!(spectrum.is_sorted());
    assert_eq!(spectrum.string_data_arrays.len(), 1);
    assert_eq!(spectrum.integer_data_arrays.len(), 1);
    assert_eq!(spectrum.float_data_arrays.len(), 1);
    assert_eq!(spectrum.string_data_arrays[0].name, "IonNames");
    assert_eq!(spectrum.float_data_arrays[0].name, "IonMatchError");
    assert_eq!(spectrum.integer_data_arrays[0].name, "Charges");
    assert_eq!(spectrum.string_data_arrays[0].data, source_labels());
    assert_eq!(spectrum.integer_data_arrays[0].data, vec![1; 11]);
    assert_eq!(
        spectrum.metadata["fragment_mass_tolerance"]
            .as_f64()
            .unwrap(),
        0.1
    );
    assert_eq!(
        spectrum.metadata["fragment_mass_tolerance_ppm"]
            .as_i64()
            .unwrap(),
        0
    );
    assert_eq!(spectrum.metadata["retained"].as_str().unwrap(), "metadata");
    for (actual, row) in spectrum.float_data_arrays[0].data.iter().zip(rows(DERIVED)) {
        near(f64::from(*actual), f64::from(f32_bits(row[3])), 1e-9);
    }
    for include_unmatched in [false, true] {
        let input = source_spectrum(true);
        let before = input.clone();
        let mut target = hit("IFSQVGK", 2);
        annotator
            .add_peak_annotations(&mut target, &input, &generator(), &align, include_unmatched)
            .unwrap();
        assert_eq!(input, before);
        assert_eq!(
            target.peak_annotations.len(),
            if include_unmatched { 13 } else { 11 }
        );
        let labels: Vec<_> = target
            .peak_annotations
            .iter()
            .filter(|p| !p.annotation.is_empty())
            .map(|p| p.annotation.clone())
            .collect();
        assert_eq!(labels, source_labels());
        for peak in target.peak_annotations {
            assert_eq!(
                peak.intensity,
                if peak.annotation.is_empty() {
                    0.5
                } else {
                    f64::from(1.1_f32)
                }
            );
            assert_eq!(peak.charge, if peak.annotation.is_empty() { 0 } else { 1 });
        }
    }
}

#[test]
fn source_statistics_keep_f32_intensities_and_positive_independent_mse() {
    let mut spectrum = source_spectrum(false);
    let mut id = identification(hit("IFSQVGK", 2));
    SpectrumAnnotator::default()
        .add_ion_match_statistics(
            &mut id,
            &mut spectrum,
            &generator(),
            &alignment(Tolerance::Absolute(0.1)),
        )
        .unwrap();
    let mut scalar_count = 0;
    for row in rows(SCALARS) {
        match row[1] {
            "string" => assert_eq!(text(&id, row[0]), row[2]),
            // C++ bool promotes to DataValue's int constructor, not a string.
            "bool" => assert_eq!(
                id.hits[0].metadata[row[0]].data(),
                &MetaValueData::Integer(i64::from(row[2].parse::<bool>().unwrap()))
            ),
            _ => {
                let expected: f64 = row[2].parse().unwrap();
                let tolerance = match row[0] {
                    "sum_intensity" | "matched_intensity" => 1e-6,
                    "topN_MSEfragmenterror" => 1e-6, // Source zero is an approximate assertion.
                    "NTermIonCurrentRatio" | "CTermIonCurrentRatio" => 1e-6,
                    _ if expected.fract() == 0.0 => 0.0,
                    _ => 1e-7,
                };
                near(number(&id, row[0]), expected, tolerance);
            }
        }
        scalar_count += 1;
    }
    assert_eq!(scalar_count, 17);
    assert_eq!(number(&id, "sum_intensity"), 11.0 * f64::from(1.1_f32));
    assert_eq!(
        number(&id, "matched_intensity"),
        number(&id, "sum_intensity")
    );
    let errors: Vec<f64> = rows(DERIVED).map(|r| f64::from(f32_bits(r[3]))).collect();
    let top: Vec<_> = errors.iter().rev().take(7).copied().collect();
    let mean = top.iter().sum::<f64>() / 7.0;
    let mse = top.iter().map(|x| x * x).sum::<f64>() / 7.0;
    let sd = (top.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / 6.0).sqrt();
    let mut sorted = errors;
    sorted.sort_by(f64::total_cmp);
    for (key, expected) in [
        ("topN_meanfragmenterror", mean),
        ("topN_MSEfragmenterror", mse),
        ("topN_stddevfragmenterror", sd),
        ("median_fragment_error", sorted[5]),
        ("IQR_fragment_error", sorted[7] - sorted[2]),
    ] {
        near(number(&id, key), expected, 1e-9);
    }
    assert!(number(&id, "topN_MSEfragmenterror") > 0.0);
}

fn controlled_trace(count: usize) -> MSSpectrum {
    let errors = [0.125, 0.25, 0.5, 1.0];
    MSSpectrum {
        ms_level: 2,
        peaks: rows(DERIVED)
            .take(count)
            .enumerate()
            .map(|(i, row)| Peak1D::new(f64_bits(row[2]) + errors[i], (1_u32 << i) as f32))
            .collect(),
        ..Default::default()
    }
}
fn controlled_stats(
    count: usize,
    options: &SpectrumAnnotator,
) -> (PeptideIdentification, MSSpectrum) {
    let mut id = identification(hit("IFSQVGK", 1));
    let mut spectrum = controlled_trace(count);
    options
        .add_ion_match_statistics(
            &mut id,
            &mut spectrum,
            &generator(),
            &alignment(Tolerance::Absolute(1.1)),
        )
        .unwrap();
    (id, spectrum)
}

#[test]
fn exact_binary_errors_pin_top_n_padding_sample_sd_and_both_signal_ratios() {
    let options = SpectrumAnnotator {
        list_of_ions_matched: false,
        fragment_error_statistics: false,
        ..Default::default()
    };
    let (id, spectrum) = controlled_stats(4, &options);
    assert_eq!(text(&id, "matched_ions"), "y1+,y2+,b2+,y3+");
    assert_eq!(
        spectrum.float_data_arrays[0].data,
        vec![0.125, 0.25, 0.5, 1.0]
    );
    near(number(&id, "topN_meanfragmenterror"), 1.875 / 7.0, 1e-15);
    near(number(&id, "topN_MSEfragmenterror"), 1.328125 / 7.0, 1e-15);
    near(
        number(&id, "topN_stddevfragmenterror"),
        ((1.328125 - 1.875_f64.powi(2) / 7.0) / 6.0).sqrt(),
        1e-15,
    );
    assert_eq!(number(&id, "median_fragment_error"), 0.375);
    assert_eq!(number(&id, "IQR_fragment_error"), 0.75);
    assert_eq!(number(&id, "sn_by_matched_intensity"), 0.0);
    assert_eq!(number(&id, "sn_by_median_intensity"), 4.0);
    assert_eq!(number(&id, "NTermIonCurrentRatio"), 4.0 / 15.0);
    assert_eq!(number(&id, "CTermIonCurrentRatio"), 11.0 / 15.0);
}

#[test]
fn small_list_quartiles_are_a_checked_native_extension_and_n_one_is_atomic_error() {
    // C++ nth_element subranges are invalid for n=1,2,3. These expected values
    // represent the documented intended sorted indices, not executed C++ goldens.
    for (n, median, iqr) in [(1, 0.125, 0.0), (2, 0.1875, 0.125), (3, 0.25, 0.125)] {
        let (id, _) = controlled_stats(n, &SpectrumAnnotator::default());
        assert_eq!(number(&id, "median_fragment_error"), median);
        assert_eq!(number(&id, "IQR_fragment_error"), iqr);
    }
    let mut id = identification(hit("IFSQVGK", 1));
    let mut spectrum = controlled_trace(4);
    let before = (id.clone(), spectrum.clone());
    let annotator = SpectrumAnnotator {
        top_n_fragment_errors: 1,
        ..Default::default()
    };
    assert!(
        annotator
            .add_ion_match_statistics(
                &mut id,
                &mut spectrum,
                &generator(),
                &alignment(Tolerance::Absolute(1.1))
            )
            .is_err()
    );
    assert_eq!((id, spectrum), before);
}

#[test]
fn repeated_ppm_targets_overwrite_arrays_but_retain_matched_only_annotations() {
    let generator = Generator {
        ion_series: vec![TheoreticalIonSeries::B],
        add_first_prefix_ion: true,
        add_metainfo: true,
        ..Default::default()
    };
    let annotator = SpectrumAnnotator::default();
    let align = alignment(Tolerance::Ppm(1_000_000.0));
    let spectrum = MSSpectrum {
        peaks: vec![Peak1D::new(50.0, 2.0)],
        ..Default::default()
    };
    let mut annotated = spectrum.clone();
    annotator
        .annotate_matches(&mut annotated, &hit("AG", 2), &generator, &align)
        .unwrap();
    assert_eq!(annotated.string_data_arrays[0].data, vec!["b1+"]);
    assert_eq!(annotated.integer_data_arrays[0].data, vec![1]);
    for include in [false, true] {
        let mut target = hit("AG", 2);
        annotator
            .add_peak_annotations(&mut target, &spectrum, &generator, &align, include)
            .unwrap();
        let expected = if include {
            vec!["b1+"]
        } else {
            vec!["b1++", "b1+"]
        };
        assert_eq!(
            target
                .peak_annotations
                .iter()
                .map(|p| p.annotation.as_str())
                .collect::<Vec<_>>(),
            expected
        );
        assert!(
            target
                .peak_annotations
                .iter()
                .all(|p| p.mz == 50.0 && p.intensity == 2.0)
        );
    }
}

#[test]
fn source_precursor_check_uses_raw_ppm_value_as_daltons_and_restores_mz_order() {
    let mut spectrum = source_spectrum(false);
    for (i, peak) in spectrum.peaks.iter_mut().enumerate() {
        peak.intensity = (11 - i) as f32;
    }
    spectrum.precursors.push(Precursor {
        mz: 150.0,
        ..Default::default()
    });
    let mut id = identification(hit("IFSQVGK", 2));
    SpectrumAnnotator::default()
        .add_ion_match_statistics(
            &mut id,
            &mut spectrum,
            &generator(),
            &alignment(Tolerance::Ppm(5.0)),
        )
        .unwrap();
    assert!(spectrum.is_sorted());
    let value = id.hits[0].metadata["precursor_in_ms2"].data();
    assert_eq!(value, &MetaValueData::Integer(1));
    assert_eq!(
        spectrum.metadata["fragment_mass_tolerance_ppm"]
            .as_i64()
            .unwrap(),
        1
    );
    assert!((spectrum.peaks[0].mz - 150.0).abs() > 150.0 * 5e-6);
    assert_eq!(
        id.metadata["fragment_match_tolerance"].as_f64().unwrap(),
        5.0
    );
}

#[test]
fn ion_name_source_literals_and_conservative_fallbacks_are_distinct_from_statistics_grammar() {
    for (name, expected) in [
        ("y4-H2O1^2/-1", 2),
        ("y3/1.2e-05", 0),
        ("y5++\nsome user comment", 2),
        ("note ^2text", 0),
        ("y3^x/--", -2),
        ("y3^x/-2", 0),
        ("y3^2^0/+9", 0),
        ("y3+-", -1),
        ("y3^2147483648", 0),
        ("y3^-2147483648", i32::MIN),
    ] {
        assert_eq!(ion_naming::charge_from_name(name), expected, "{name}");
    }
    for charge in [i32::MIN, -20, -9, -1, 1, 8, 9, 20, i32::MAX] {
        assert_eq!(
            ion_naming::charge_from_name(&format!("y5{}", ion_naming::charge_suffix(charge))),
            charge
        );
    }
    assert_eq!(
        ion_naming::with_charge("y3\r\nλ comment", 2).unwrap(),
        "y3++\r\nλ comment"
    );
    assert_eq!(ion_naming::with_charge("y3^2/-1", 3).unwrap(), "y3^2/-1");
    for (name, ordinal) in [
        ("y3-H2O1+", 3),
        ("y3+12", 3),
        ("[alpha|ci$y3]++", 0),
        ("iY+U-H3PO4+", 0),
        ("y000000003^2", 3),
        ("y0000000003^2", 0),
        ("A999999999", 999_999_999),
    ] {
        assert_eq!(ion_naming::ordinal_from_name(name), ordinal, "{name}");
    }
}

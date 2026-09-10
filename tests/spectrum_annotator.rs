// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source: pinned SpectrumAnnotator.cpp and focused native boundary regressions.

use openms::chemistry::spectrum_annotator::{
    MAX_ANNOTATION_HITS, MAX_ANNOTATION_PEAKS, MAX_ANNOTATION_TOP_N,
};
use openms::chemistry::{AASequence, SpectrumAnnotator, TheoreticalSpectrumGenerator};
use openms::comparison::{SpectrumAlignment, Tolerance};
use openms::identification::{PeakAnnotation, PeptideHit, PeptideIdentification};
use openms::kernel::{DataArray, MSSpectrum, Peak1D, Precursor};
use openms::metadata::MetaValue;
use std::collections::BTreeMap;

fn hit() -> PeptideHit {
    PeptideHit {
        sequence: AASequence::parse("PEPTIDE").unwrap(),
        charge: 1,
        ..Default::default()
    }
}
fn generator() -> TheoreticalSpectrumGenerator {
    TheoreticalSpectrumGenerator {
        add_metainfo: true,
        ..Default::default()
    }
}
fn alignment() -> SpectrumAlignment {
    SpectrumAlignment {
        tolerance: Tolerance::Absolute(0.1),
        ..Default::default()
    }
}
fn measured() -> MSSpectrum {
    let mut spectrum = generator().generate(&hit().sequence, 1, 1, None).unwrap();
    spectrum.float_data_arrays.clear();
    spectrum.integer_data_arrays.clear();
    spectrum.string_data_arrays.clear();
    spectrum.precursors.clear();
    spectrum
}
fn identification() -> PeptideIdentification {
    PeptideIdentification {
        hits: vec![hit()],
        ..Default::default()
    }
}
fn without_ratios() -> SpectrumAnnotator {
    SpectrumAnnotator {
        sn_statistics: false,
        terminal_series_match_ratio: false,
        ..Default::default()
    }
}
fn number(id: &PeptideIdentification, name: &str) -> f64 {
    id.hits[0].metadata[name].as_f64().unwrap()
}

#[test]
fn array_replacement_is_aligned_and_does_not_validate_unrelated_state() {
    let original = generator().generate(&hit().sequence, 1, 1, None).unwrap();
    let expected: BTreeMap<_, _> = original
        .peaks
        .iter()
        .zip(&original.string_data_arrays[0].data)
        .map(|(peak, name)| (peak.mz.to_bits(), name.clone()))
        .collect();
    let mut spectrum = measured();
    spectrum.peaks.reverse();
    // These old arrays will all be discarded, not reordered or deep-cloned.
    spectrum.float_data_arrays = vec![DataArray::new("profile", vec![f32::NAN])];
    spectrum.integer_data_arrays = vec![DataArray::new("old", vec![])];
    spectrum.string_data_arrays = vec![DataArray::new("old", vec!["retired".into()])];
    spectrum.rt = f64::NAN;
    spectrum.ms_level = 0;
    spectrum.native_id = "scan=25".into();
    spectrum.metadata.insert("keep".into(), "unrelated".into());
    spectrum.precursors.push(Precursor {
        activation_energy: f64::NAN,
        ..Default::default()
    });
    spectrum
        .peptide_identifications
        .push(PeptideIdentification {
            hits: vec![PeptideHit {
                score: f64::NAN,
                ..hit()
            }],
            ..Default::default()
        });
    SpectrumAnnotator::default()
        .annotate_matches(&mut spectrum, &hit(), &generator(), &alignment())
        .unwrap();
    assert!(spectrum.is_sorted());
    assert_eq!(spectrum.float_data_arrays[0].name, "IonMatchError");
    assert_eq!(spectrum.integer_data_arrays[0].name, "Charges");
    assert_eq!(spectrum.string_data_arrays[0].name, "IonNames");
    assert_eq!(spectrum.float_data_arrays.len(), 1);
    assert_eq!(spectrum.integer_data_arrays.len(), 1);
    assert_eq!(spectrum.string_data_arrays.len(), 1);
    for (i, peak) in spectrum.peaks.iter().enumerate() {
        assert_eq!(
            spectrum.string_data_arrays[0].data[i],
            expected[&peak.mz.to_bits()]
        );
        assert_eq!(spectrum.integer_data_arrays[0].data[i], 1);
        assert_eq!(spectrum.float_data_arrays[0].data[i], 0.0);
    }
    assert!(spectrum.rt.is_nan());
    assert_eq!(spectrum.ms_level, 0);
    assert!(spectrum.precursors[0].activation_energy.is_nan());
    assert!(spectrum.peptide_identifications[0].hits[0].score.is_nan());
    assert_eq!(spectrum.native_id, "scan=25");
    assert_eq!(spectrum.metadata["keep"], "unrelated");
    assert_eq!(spectrum.metadata["fragment_mass_tolerance"], "0.1");
    assert_eq!(spectrum.metadata["fragment_mass_tolerance_ppm"], "0");
}

#[test]
fn absent_theoretical_metadata_follows_each_source_entrypoint() {
    let generator = TheoreticalSpectrumGenerator::default();
    let annotator = SpectrumAnnotator::default();
    let mut spectrum = measured();
    let saved = spectrum.clone();
    assert!(
        annotator
            .annotate_matches(&mut spectrum, &hit(), &generator, &alignment())
            .is_err()
    );
    assert_eq!(spectrum, saved);
    let mut candidate = hit();
    candidate.score = f64::NAN; // Scoring and old annotations are not consumed.
    candidate.peak_annotations.push(PeakAnnotation {
        mz: f64::NAN,
        ..Default::default()
    });
    annotator
        .add_peak_annotations(&mut candidate, &spectrum, &generator, &alignment(), false)
        .unwrap();
    assert_eq!(candidate.peak_annotations.len(), spectrum.len());
    assert!(
        candidate
            .peak_annotations
            .iter()
            .all(|a| a.annotation.is_empty() && a.charge == 0)
    );
    assert!(candidate.score.is_nan());
    assert_eq!(spectrum, saved);
    // No match means source never indexes the absent arrays.
    spectrum.peaks = vec![Peak1D::new(0.0, 1.0)];
    annotator
        .annotate_matches(&mut spectrum, &hit(), &generator, &alignment())
        .unwrap();
    assert_eq!(spectrum.string_data_arrays[0].data, [""]);
    spectrum.peaks.clear();
    annotator
        .annotate_matches(&mut spectrum, &hit(), &generator, &alignment())
        .unwrap();
    assert!(spectrum.string_data_arrays[0].data.is_empty());
}

#[test]
fn source_noop_flags_and_disabled_statistics_preserve_existing_keys() {
    let mut spectrum = measured();
    let mut id = identification();
    let annotator = SpectrumAnnotator {
        list_of_ions_matched: false,
        fragment_error_statistics: false,
        ..Default::default()
    };
    annotator
        .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
        .unwrap();
    assert!(
        !id.hits[0].metadata["matched_ions"]
            .as_str()
            .unwrap()
            .is_empty()
    );
    assert_eq!(number(&id, "topN_meanfragmenterror"), 0.0);
    let mut id = identification();
    id.hits[0]
        .metadata
        .insert("matched_ions".into(), "previous".into());
    id.hits[0]
        .metadata
        .insert("topN_meanfragmenterror".into(), 42_i32.into());
    let disabled = SpectrumAnnotator {
        basic_statistics: false,
        top_n_fragment_errors: 0,
        max_series: false,
        ..Default::default()
    };
    disabled
        .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
        .unwrap();
    assert_eq!(
        id.hits[0].metadata["matched_ions"].as_str().unwrap(),
        "previous"
    );
    assert_eq!(number(&id, "topN_meanfragmenterror"), 42.0);
    assert!(!id.hits[0].metadata.contains_key("peak_number"));
    assert!(!id.hits[0].metadata.contains_key("max_series_type"));
}

#[test]
fn statistics_last_hit_controls_annotations_and_precursor_presence_controls_sort() {
    let mut spectrum = measured();
    let n = spectrum.len();
    for (i, peak) in spectrum.peaks.iter_mut().enumerate() {
        peak.intensity = (n - i) as f32;
    }
    let mut id = identification();
    SpectrumAnnotator::default()
        .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
        .unwrap();
    assert!(spectrum.peaks.windows(2).all(|p| p[0].mz > p[1].mz));
    assert_eq!(
        id.metadata["fragment_match_tolerance"].as_f64().unwrap(),
        0.1
    );
    let mut second = hit();
    second.sequence = AASequence::parse("AAAAAAA").unwrap();
    id.hits.push(second);
    spectrum.precursors.push(Precursor::new(1e9, 2));
    without_ratios()
        .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
        .unwrap();
    assert!(spectrum.is_sorted());
    assert_eq!(
        id.hits[1].metadata["matched_ion_number"].as_i64().unwrap(),
        0
    );
    assert!(
        spectrum.string_data_arrays[0]
            .data
            .iter()
            .all(String::is_empty)
    );
    assert!(id.hits[0].metadata["matched_ion_number"].as_i64().unwrap() > 0);
    assert_eq!(
        id.hits[1].metadata["precursor_in_ms2"],
        MetaValue::from(0_i32)
    );
}

#[test]
fn a_late_failure_rolls_back_both_inputs_and_all_three_entrypoints_are_atomic() {
    let mut spectrum = measured();
    spectrum.metadata.insert("keep".into(), "original".into());
    let mut id = identification();
    id.hits[0].metadata.insert("keep".into(), "original".into());
    id.hits.push(PeptideHit { charge: 0, ..hit() });
    let saved_spectrum = spectrum.clone();
    let saved_id = id.clone();
    assert!(
        SpectrumAnnotator::default()
            .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
            .is_err()
    );
    assert_eq!(id, saved_id);
    assert_eq!(spectrum, saved_spectrum);
    for charge in [0, -1, i32::MIN] {
        let mut candidate = PeptideHit { charge, ..hit() };
        let saved_hit = candidate.clone();
        assert!(
            SpectrumAnnotator::default()
                .annotate_matches(&mut spectrum, &candidate, &generator(), &alignment())
                .is_err()
        );
        assert!(
            SpectrumAnnotator::default()
                .add_peak_annotations(&mut candidate, &spectrum, &generator(), &alignment(), true)
                .is_err()
        );
        assert_eq!(candidate, saved_hit);
        assert_eq!(spectrum, saved_spectrum);
    }
    // A positive source charge is clamped before narrowing into the generator API.
    SpectrumAnnotator::default()
        .annotate_matches(
            &mut spectrum,
            &PeptideHit {
                charge: i32::MAX,
                ..hit()
            },
            &generator(),
            &alignment(),
        )
        .unwrap();
}

#[test]
fn finite_signed_ratios_and_f32_zero_group_fallbacks_are_preserved() {
    let mut spectrum = measured();
    spectrum.peaks.truncate(1);
    spectrum.peaks[0].intensity = 1.0;
    spectrum.peaks.push(Peak1D::new(0.0, -2.0));
    let mut id = identification();
    SpectrumAnnotator::default()
        .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
        .unwrap();
    assert_eq!(number(&id, "sn_by_matched_intensity"), -0.5);
    assert_eq!(number(&id, "sn_by_median_intensity"), -0.5);
    // The source's f32 median overflows to infinity, leaving no signal group.
    // Its explicit group-empty fallback still produces a finite zero statistic.
    let mut spectrum = measured();
    spectrum.peaks.truncate(2);
    for peak in &mut spectrum.peaks {
        peak.intensity = f32::MAX;
    }
    let mut id = identification();
    SpectrumAnnotator::default()
        .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
        .unwrap();
    assert_eq!(number(&id, "sn_by_matched_intensity"), 0.0);
    assert_eq!(number(&id, "sn_by_median_intensity"), 0.0);
    assert_eq!(number(&id, "sum_intensity"), 2.0 * f64::from(f32::MAX));
}

#[test]
fn undefined_enabled_ratios_and_float_error_overflow_are_checked() {
    let mut spectrum = measured();
    spectrum.peaks.truncate(1);
    spectrum.peaks.push(Peak1D::new(0.0, 0.0));
    let mut id = identification();
    let saved = (id.clone(), spectrum.clone());
    assert!(
        SpectrumAnnotator::default()
            .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
            .unwrap_err()
            .to_string()
            .contains("sn_by_matched_intensity")
    );
    assert_eq!((id, spectrum), saved);
    let mut spectrum = MSSpectrum {
        peaks: vec![Peak1D::new(1e100, 1.0)],
        ..Default::default()
    };
    let saved = spectrum.clone();
    let broad = SpectrumAlignment {
        tolerance: Tolerance::Absolute(1e100),
        ..Default::default()
    };
    assert!(
        SpectrumAnnotator::default()
            .annotate_matches(&mut spectrum, &hit(), &generator(), &broad)
            .unwrap_err()
            .to_string()
            .contains("f32")
    );
    assert_eq!(spectrum, saved);
}

#[test]
fn top_n_one_is_valid_only_for_empty_matched_error_lists() {
    let config = SpectrumAnnotator {
        top_n_fragment_errors: 1,
        ..without_ratios()
    };
    let mut spectrum = measured();
    spectrum.peaks.truncate(1);
    let mut id = identification();
    let saved = (id.clone(), spectrum.clone());
    assert!(
        config
            .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
            .is_err()
    );
    assert_eq!((id.clone(), spectrum), saved);
    let mut spectrum = MSSpectrum {
        peaks: vec![Peak1D::new(0.0, 1.0)],
        ..Default::default()
    };
    config
        .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
        .unwrap();
    for key in [
        "median_fragment_error",
        "IQR_fragment_error",
        "topN_meanfragmenterror",
        "topN_MSEfragmenterror",
        "topN_stddevfragmenterror",
    ] {
        assert_eq!(number(&id, key), 0.0);
    }
}

#[test]
fn statistics_empty_input_shortcuts_and_scope_limits_precede_unused_work() {
    let invalid = SpectrumAnnotator {
        top_n_fragment_errors: usize::MAX,
        ..Default::default()
    };
    let mut empty = MSSpectrum::default();
    let mut id = PeptideIdentification {
        hits: vec![PeptideHit::default(); MAX_ANNOTATION_HITS + 1],
        ..Default::default()
    };
    invalid
        .add_ion_match_statistics(&mut id, &mut empty, &generator(), &alignment())
        .unwrap();
    let mut spectrum = MSSpectrum {
        peaks: vec![Peak1D::new(f64::NAN, f32::NAN)],
        ..Default::default()
    };
    let mut no_hits = PeptideIdentification::default();
    invalid
        .add_ion_match_statistics(&mut no_hits, &mut spectrum, &generator(), &alignment())
        .unwrap();
    assert!(spectrum.peaks[0].mz.is_nan());
    let mut spectrum = measured();
    assert!(
        SpectrumAnnotator::default()
            .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
            .unwrap_err()
            .to_string()
            .contains("hit limit")
    );
    let mut id = identification();
    let excessive = SpectrumAnnotator {
        top_n_fragment_errors: MAX_ANNOTATION_TOP_N + 1,
        ..Default::default()
    };
    assert!(
        excessive
            .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
            .unwrap_err()
            .to_string()
            .contains("top-N limit")
    );
    let mut too_large = MSSpectrum {
        peaks: vec![Peak1D::new(0.0, 0.0); MAX_ANNOTATION_PEAKS + 1],
        ..Default::default()
    };
    assert!(
        SpectrumAnnotator::default()
            .annotate_matches(&mut too_large, &hit(), &generator(), &alignment())
            .unwrap_err()
            .to_string()
            .contains("peak limit")
    );
    assert_eq!(too_large.len(), MAX_ANNOTATION_PEAKS + 1);
}

#[test]
fn consumed_precursor_mz_is_validated_without_reading_other_acquisition_fields() {
    let mut spectrum = measured();
    spectrum.precursors.push(Precursor {
        mz: spectrum.peaks[0].mz,
        activation_energy: f64::NAN,
        ..Default::default()
    });
    let mut id = identification();
    SpectrumAnnotator::default()
        .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
        .unwrap();
    assert_eq!(
        id.hits[0].metadata["precursor_in_ms2"],
        MetaValue::from(1_i32)
    );
    assert!(spectrum.precursors[0].activation_energy.is_nan());
    spectrum.precursors[0].mz = f64::NAN;
    id.hits[0]
        .metadata
        .insert("keep".into(), MetaValue::from("original"));
    let old_peaks = spectrum.peaks.clone();
    let old_arrays = spectrum.string_data_arrays.clone();
    let old_id = id.clone();
    assert!(
        SpectrumAnnotator::default()
            .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
            .unwrap_err()
            .to_string()
            .contains("precursor m/z")
    );
    assert_eq!(id, old_id);
    assert_eq!(spectrum.peaks, old_peaks);
    assert_eq!(spectrum.string_data_arrays, old_arrays);
}

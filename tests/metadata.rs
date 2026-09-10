// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::kernel::{NumericRange, Precursor, SpectrumType};
use openms::metadata::*;
use std::collections::BTreeMap;

#[test]
fn data_value_types_empty_strings_lists_and_boolean_conventions() {
    let empty = MetaValue::default();
    assert!(empty.is_empty());
    assert_eq!(empty.to_string(), "");
    assert!(empty.as_str().is_err());
    assert!(empty.as_i64().is_err());
    assert!(empty.as_f64().is_err());
    let string = MetaValue::from("");
    assert!(!string.is_empty());
    assert_eq!(string.as_str().unwrap(), "");
    let value = MetaValue::from(5i64);
    assert_eq!(value.as_i64().unwrap(), 5);
    assert_eq!(value.as_f64().unwrap(), 5.0);
    assert_eq!(value.to_string(), "5");
    assert!(value.as_str().is_err());
    assert!(value.to_bool().is_err());
    let value = MetaValue::try_from(47.11).unwrap();
    assert_eq!(value.as_f64().unwrap(), 47.11);
    assert!(value.as_i64().is_err());
    let values = MetaValue::from(vec![
        "test string".to_owned(),
        "string2".to_owned(),
        "last string".to_owned(),
    ]);
    assert_eq!(values.to_string(), "[test string, string2, last string]");
    assert_eq!(values.as_string_list().unwrap().len(), 3);
    let values = MetaValue::from(vec![1i64, 2, 3, 4, 5]);
    assert_eq!(values.to_string(), "[1, 2, 3, 4, 5]");
    assert_eq!(values.as_integer_list().unwrap(), [1, 2, 3, 4, 5]);
    let values = MetaValue::try_from(vec![1.2, 47.11]).unwrap();
    assert_eq!(values.as_float_list().unwrap(), [1.2, 47.11]);
    assert!(MetaValue::from("true").to_bool().unwrap());
    assert!(!MetaValue::from("false").to_bool().unwrap());
    for value in ["TRUE", "False", "1", "yes", ""] {
        assert!(MetaValue::from(value).to_bool().is_err());
    }
}

#[test]
fn numeric_validation_precision_and_equality_are_explicit() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(MetaValue::try_from(value).is_err());
        assert!(MetaValue::new(MetaValueData::Float(value)).is_err());
        assert!(MetaValue::try_from(vec![1.0, value]).is_err());
    }
    assert!(MetaValue::try_from(u64::MAX).is_err());
    assert_eq!(
        MetaValue::try_from(i64::MAX as u64)
            .unwrap()
            .as_i64()
            .unwrap(),
        i64::MAX
    );
    assert_ne!(MetaValue::from(1), MetaValue::try_from(1.0).unwrap());
    assert_ne!(
        MetaValue::try_from(1.0).unwrap(),
        MetaValue::try_from(1.0000001).unwrap()
    );
    assert_eq!(
        MetaValue::try_from(-0.0).unwrap(),
        MetaValue::try_from(0.0).unwrap()
    );
    assert!(
        MetaValue::from(vec![] as Vec<String>)
            .as_float_list()
            .is_err()
    );
    let value = MetaValue::from(9_007_199_254_740_993i64);
    assert_eq!(value.as_f64().unwrap(), 9_007_199_254_740_992.0); // Documented f64 conversion.
}

#[test]
fn units_are_owned_and_participate_in_equality() {
    let seconds = Unit::new("UO:0000010", "second", "UO").unwrap();
    assert_eq!(seconds.accession(), "UO:0000010");
    assert_eq!(seconds.name(), "second");
    assert_eq!(seconds.cv_ref(), "UO");
    assert_eq!(
        Unit::from_ontology(UnitOntology::Unit, 10).accession(),
        "UO:0000010"
    );
    assert_eq!(
        Unit::from_ontology(UnitOntology::MassSpectrometry, 1000040).accession(),
        "MS:1000040"
    );
    let plain = MetaValue::try_from(47.11).unwrap();
    let value = plain.clone().with_unit(seconds.clone()).unwrap();
    assert_eq!(value.unit(), Some(&seconds));
    assert_ne!(value, plain);
    assert_eq!(value.to_string(), "47.11");
    assert_eq!(value.without_unit(), plain);
    assert!(Unit::new("", "", "UO").is_err());
    assert!(Unit::new("UO: 1", "", "UO").is_err());
}

#[test]
fn legacy_string_maps_round_trip_without_guessing_types() {
    let old = BTreeMap::from([
        ("score".into(), "47.11".into()),
        ("scan".into(), "123".into()),
    ]);
    let typed = meta_from_strings(&old);
    assert_eq!(typed["scan"].as_str().unwrap(), "123");
    assert!(typed["scan"].as_i64().is_err());
    assert_eq!(meta_to_strings(&typed).unwrap(), old);
    let mut typed = typed;
    typed.insert("scan".into(), 123.into());
    assert!(meta_to_strings(&typed).is_err());
    assert_eq!(meta_to_strings_lossy(&typed), old);
    typed.insert(
        "scan".into(),
        MetaValue::from("123")
            .with_unit(Unit::from_ontology(UnitOntology::Unit, 10))
            .unwrap(),
    );
    assert!(meta_to_strings(&typed).is_err());
}

#[test]
fn metadata_merges_overwrite_like_source_and_reject_conflicts_atomically() {
    let mut left = MetaInfo::from([("one".into(), 1.into()), ("two".into(), 2.into())]);
    let right = MetaInfo::from([("one".into(), 3.into()), ("three".into(), 4.into())]);
    let before = left.clone();
    assert!(merge_meta(&mut left, &right, MetaMergePolicy::RejectConflicts).is_err());
    assert_eq!(left, before);
    merge_meta(&mut left, &right, MetaMergePolicy::KeepExisting).unwrap();
    assert_eq!(left["one"].as_i64().unwrap(), 1);
    assert_eq!(left["three"].as_i64().unwrap(), 4);
    merge_meta(&mut left, &right, MetaMergePolicy::Overwrite).unwrap();
    assert_eq!(left["one"].as_i64().unwrap(), 3);
    merge_meta(&mut left, &right, MetaMergePolicy::RejectConflicts).unwrap();
    validate_meta(&left).unwrap();
}

fn term(accession: &str, value: &str) -> CVTerm {
    CVTerm {
        value: value.into(),
        ..CVTerm::new(accession, "term name", "MS")
    }
}

#[test]
fn cv_terms_distinguish_absent_and_empty_values_with_independent_units() {
    let mut t = CVTerm::default();
    assert!(!t.has_value());
    assert!(!t.has_unit());
    t.value = MetaValue::from("");
    assert!(t.has_value());
    t.value = MetaValue::default()
        .with_unit(Unit::from_ontology(UnitOntology::Unit, 10))
        .unwrap();
    assert!(!t.has_value());
    assert!(t.has_unit());
    assert!(t.validate().is_err());
    t.accession = "MS:1000016".into();
    t.validate().unwrap();
}

#[test]
fn source_cv_append_replace_consume_and_empty_key_semantics() {
    let mut list = CVTermList::new();
    list.metadata.insert("meta".into(), "keep".into());
    assert!(list.is_empty());
    list.add_terms(&[term("my_accession", "3.0"), term("my_accession2", "2.0")])
        .unwrap();
    list.add_terms(&[term("my_accession", "1.0")]).unwrap();
    assert_eq!(list.get("my_accession").unwrap().len(), 2);
    list.replace(term("my_accession", "2.0")).unwrap();
    assert_eq!(
        list.get("my_accession").unwrap()[0].value.as_str().unwrap(),
        "2.0"
    );
    list.replace_accession(
        "my_accession",
        vec![term("my_accession", "3.0"), term("my_accession", "2.0")],
    )
    .unwrap();
    let mut other = CVTermList::new();
    other.add(term("my_accession", "4.0")).unwrap();
    other.metadata.insert("meta".into(), "ignored".into());
    list.consume(&other).unwrap();
    assert_eq!(list.get("my_accession").unwrap().len(), 3);
    assert_eq!(list.metadata["meta"].as_str().unwrap(), "keep");
    list.replace_all(BTreeMap::from([("empty".into(), vec![])]))
        .unwrap();
    assert!(list.contains("empty"));
    assert!(!list.is_empty());
    assert_eq!(list.get("empty"), Some([].as_slice()));
    list.remove("empty");
    assert!(list.is_empty());
}

#[test]
fn cv_list_invalid_changes_are_atomic_and_cannot_break_accession_index() {
    let mut list = CVTermList::new();
    list.add(term("a", "one")).unwrap();
    let before = list.clone();
    assert!(
        list.add_terms(&[term("b", "two"), CVTerm::default()])
            .is_err()
    );
    assert_eq!(list, before);
    assert!(list.replace_accession("a", vec![term("b", "bad")]).is_err());
    assert_eq!(list, before);
    assert!(
        list.replace_all(BTreeMap::from([("a".into(), vec![term("b", "bad")])]))
            .is_err()
    );
    assert_eq!(list, before);
    assert!(list.replace_accession("", vec![]).is_err());
    assert_eq!(list, before);
}

#[test]
fn source_enum_names_order_and_defaults() {
    assert_eq!(ActivationMethod::ALL.len(), 19);
    assert_eq!(
        ActivationMethod::ALL
            .iter()
            .map(|m| m.short_name())
            .collect::<Vec<_>>(),
        vec![
            "CID", "PSD", "PD", "SID", "BIRD", "ECD", "IMD", "SORI", "HCID", "LCID", "PHD", "ETD",
            "ETciD", "EThcD", "PQD", "TRAP", "HCD", "INSOURCE", "LIFT"
        ]
    );
    for &method in ActivationMethod::ALL {
        assert_eq!(method.name().parse::<ActivationMethod>().unwrap(), method);
        assert_eq!(
            method.short_name().parse::<ActivationMethod>().unwrap(),
            method
        );
    }
    assert_eq!(
        "HCD".parse::<ActivationMethod>().unwrap().name(),
        "beam-type collision-induced dissociation"
    );
    assert!("hcd".parse::<ActivationMethod>().is_err());
    assert!(
        "SIZE_OF_ACTIVATIONMETHOD"
            .parse::<ActivationMethod>()
            .is_err()
    );
    assert_eq!(ScanMode::ALL.len(), 15);
    assert_eq!(ScanMode::Ms1Spectrum.name(), "MS1Spectrum");
    assert_eq!(ScanMode::MsnSpectrum.name(), "MSnSpectrum");
    for &mode in ScanMode::ALL {
        assert_eq!(mode.to_string().parse::<ScanMode>().unwrap(), mode);
    }
    assert_eq!(ProcessingAction::ALL.len(), 22);
    assert_eq!(
        ProcessingAction::ConversionMzML.name(),
        "Conversion to mzML format"
    );
    for &action in ProcessingAction::ALL {
        assert_eq!(
            action.to_string().parse::<ProcessingAction>().unwrap(),
            action
        );
    }
    assert_eq!(ChromatogramType::default(), ChromatogramType::Mass);
    assert_eq!(
        ChromatogramType::ALL.last(),
        Some(&ChromatogramType::Unknown)
    );
    assert_eq!(DriftTimeUnit::default().name(), "<NONE>");
    assert_eq!(Polarity::default().name(), "unknown");
    assert_eq!(ChecksumType::Sha1.name(), "SHA-1");
}

#[test]
fn source_acquisition_defaults_are_complete_and_valid() {
    let precursor = PrecursorInfo::default();
    precursor.validate().unwrap();
    assert_eq!(precursor.peak, Precursor::default());
    assert_eq!(precursor.drift_time, None);
    assert_eq!(precursor.drift_time_unit, DriftTimeUnit::None);
    assert!(precursor.activation_methods.is_empty());
    assert_eq!(precursor.activation_energy, 0.0);
    assert_eq!(
        precursor.isolation_window().unwrap(),
        NumericRange { min: 0.0, max: 0.0 }
    );
    let instrument = InstrumentSettings::default();
    instrument.validate().unwrap();
    assert_eq!(instrument.scan_mode, ScanMode::Unknown);
    assert!(!instrument.zoom_scan);
    assert!(instrument.scan_windows.is_empty());
    SourceFile::default().validate().unwrap();
    Product::default().validate().unwrap();
    Acquisition::default().validate().unwrap();
    AcquisitionInfo::default().validate().unwrap();
    DataProcessing::default().validate().unwrap();
    let settings = SpectrumSettings::default();
    settings.validate().unwrap();
    assert_eq!(settings.spectrum_type, SpectrumType::Unknown);
    assert_eq!(settings.ion_mobility_format, IonMobilityFormat::Unknown);
    assert_eq!(
        settings.ion_mobility_peak_type,
        IonMobilityPeakType::Unknown
    );
    let settings = ChromatogramSettings::default();
    settings.validate().unwrap();
    assert_eq!(settings.chromatogram_type, ChromatogramType::Mass);
}

#[test]
fn precursor_windows_activation_and_source_mass_convention() {
    let mut info = PrecursorInfo::from(Precursor::new(500.5, 2));
    info.activation_energy = 47.11;
    info.isolation_window_lower_offset = 22.8;
    info.isolation_window_upper_offset = 22.7;
    info.activation_methods.extend([
        ActivationMethod::Hcd,
        ActivationMethod::Cid,
        ActivationMethod::Cid,
    ]);
    assert_eq!(info.activation_methods.len(), 2);
    info.possible_charge_states = vec![2, 3, 2, -1, 0];
    info.drift_time = Some(-45.0);
    info.drift_time_unit = DriftTimeUnit::FaimsCompensationVoltage;
    info.validate().unwrap();
    assert_eq!(
        info.isolation_window().unwrap(),
        NumericRange {
            min: 477.7,
            max: 523.2
        }
    );
    let expected = 2.0 * 500.5 - 2.0 * openms::chemistry::PROTON_MASS_U;
    assert_eq!(info.uncharged_mass().unwrap(), expected);
    info.peak.charge = 0;
    assert_eq!(info.uncharged_mass().unwrap(), expected);
    info.peak.charge = -2;
    assert_eq!(info.uncharged_mass().unwrap(), -expected);
    assert_eq!(info.possible_charge_states, [2, 3, 2, -1, 0]);
    info.isolation_window_lower_offset = -1.0;
    assert!(info.validate().is_err());
    info.isolation_window_lower_offset = 0.0;
    info.drift_time = Some(f64::NAN);
    assert!(info.validate().is_err());
    info.drift_time = None;
    info.peak.mz = f64::MAX;
    info.isolation_window_upper_offset = f64::MAX;
    assert!(info.isolation_window().is_err());
    assert!(info.uncharged_mass().is_err());
}

#[test]
fn product_and_scan_windows_are_inclusive_and_finite() {
    let product = Product {
        mz: 400.0,
        isolation_window_lower_offset: 1.5,
        isolation_window_upper_offset: 2.0,
        ..Product::default()
    };
    assert_eq!(
        product.isolation_window().unwrap(),
        NumericRange {
            min: 398.5,
            max: 402.0
        }
    );
    let window = ScanWindow::new(50.0, 1500.0).unwrap();
    for mz in [50.0, 1000.0, 1500.0] {
        assert!(window.contains(mz).unwrap());
    }
    assert!(!window.contains(1500.01).unwrap());
    assert!(window.contains(f64::NAN).is_err());
    assert!(ScanWindow::new(2.0, 1.0).is_err());
    assert!(ScanWindow::new(0.0, f64::INFINITY).is_err());
}

#[test]
fn source_file_checksums_sizes_and_acquisition_values() {
    let mut file = SourceFile {
        name: "test.mzML".into(),
        path: "/data".into(),
        size_mb: 47.11,
        checksum_type: ChecksumType::Sha1,
        checksum: "0123456789abcdef0123456789abcdef01234567".into(),
        ..SourceFile::default()
    };
    file.validate().unwrap();
    file.checksum.pop();
    assert!(file.validate().is_err());
    file.checksum_type = ChecksumType::Unknown;
    file.validate().unwrap();
    file.size_mb = -1.0;
    assert!(file.validate().is_err());
    let info = AcquisitionInfo {
        acquisitions: vec![Acquisition {
            identifier: "acq_1".into(),
            metadata: MetaInfo::from([("scan".into(), 42.into())]),
        }],
        method_of_combination: "sum".into(),
        ..AcquisitionInfo::default()
    };
    info.validate().unwrap();
    assert_eq!(info.acquisitions[0].metadata["scan"].as_i64().unwrap(), 42);
}

#[test]
fn completion_time_calendar_validation_does_not_guess_timezones() {
    for text in [
        "2024-02-29 23:59:59",
        "2000-02-29T00:00:00",
        "0001-01-01 00:00:00",
    ] {
        let stamp = text.parse::<CompletionTime>().unwrap();
        assert_eq!(stamp.to_string().parse::<CompletionTime>().unwrap(), stamp);
    }
    for text in [
        "1900-02-29 00:00:00",
        "2025-02-29 00:00:00",
        "2024-02-30 00:00:00",
        "2024-00-01 00:00:00",
        "2024-01-01 24:00:00",
        "2024-01-01 00:00:60",
        "0000-01-01 00:00:00",
        "2024-01-01 00:00:00Z",
        "２０２４-01-01 00:00:00",
        "2024-01-01 0+:00:00",
    ] {
        assert!(text.parse::<CompletionTime>().is_err(), "{text}");
    }
}

#[test]
fn spectrum_settings_unify_preserves_source_overwrite_append_rules() {
    let mut left = SpectrumSettings {
        spectrum_type: SpectrumType::Profile,
        native_id: "left".into(),
        comment: "Original Comment".into(),
        ..SpectrumSettings::default()
    };
    left.metadata.insert("1".into(), "will be gone".into());
    left.metadata
        .insert("2".into(), "will be still present".into());
    left.precursors.push(Precursor::new(1.0, 2).into());
    left.products.push(Product {
        mz: 1.0,
        ..Product::default()
    });
    left.data_processing.push(DataProcessing {
        software: Software {
            name: "org_software".into(),
            ..Software::default()
        },
        ..DataProcessing::default()
    });
    left.instrument_settings.polarity = Polarity::Positive;
    left.acquisition_info.method_of_combination = "keep".into();
    left.source_file.name = "left.mzML".into();
    let mut right = SpectrumSettings {
        spectrum_type: SpectrumType::Profile,
        native_id: "right".into(),
        comment: "Appended to org Commment".into(),
        ..SpectrumSettings::default()
    };
    right
        .metadata
        .insert("1".into(), "will overwrite org comment".into());
    right.precursors.push(Precursor::new(2.0, 3).into());
    right.products.push(Product {
        mz: 2.0,
        ..Product::default()
    });
    right.data_processing.push(DataProcessing {
        software: Software {
            name: "appended_software".into(),
            ..Software::default()
        },
        ..DataProcessing::default()
    });
    left.unify(&right).unwrap();
    assert_eq!(
        left.metadata["1"].as_str().unwrap(),
        "will overwrite org comment"
    );
    assert_eq!(
        left.metadata["2"].as_str().unwrap(),
        "will be still present"
    );
    assert_eq!(left.comment, "Original CommentAppended to org Commment");
    assert_eq!(
        left.precursors
            .iter()
            .map(|p| p.peak.mz)
            .collect::<Vec<_>>(),
        vec![1.0, 2.0]
    );
    assert_eq!(
        left.products.iter().map(|p| p.mz).collect::<Vec<_>>(),
        vec![1.0, 2.0]
    );
    assert_eq!(left.data_processing[1].software.name, "appended_software");
    assert_eq!(left.spectrum_type, SpectrumType::Profile);
    assert_eq!(left.native_id, "left");
    assert_eq!(left.instrument_settings.polarity, Polarity::Positive);
    assert_eq!(left.acquisition_info.method_of_combination, "keep");
    assert_eq!(left.source_file.name, "left.mzML");
    left.unify(&SpectrumSettings {
        spectrum_type: SpectrumType::Centroid,
        ..SpectrumSettings::default()
    })
    .unwrap();
    assert_eq!(left.spectrum_type, SpectrumType::Unknown);
}

#[test]
fn settings_unify_invalid_nested_values_is_atomic_and_equality_includes_mobility() {
    let mut left = SpectrumSettings::default();
    let before = left.clone();
    let mut right = SpectrumSettings {
        comment: "must not append".into(),
        ..SpectrumSettings::default()
    };
    right.products.push(Product {
        isolation_window_upper_offset: -1.0,
        ..Product::default()
    });
    assert!(left.unify(&right).is_err());
    assert_eq!(left, before);
    let mut changed = before.clone();
    changed.ion_mobility_format = IonMobilityFormat::PerPeak;
    assert_ne!(changed, before);
    let mut changed = before.clone();
    changed.ion_mobility_peak_type = IonMobilityPeakType::Centroid;
    assert_ne!(changed, before);
    let mut chromatogram = ChromatogramSettings::default();
    chromatogram.source_file.size_mb = f32::INFINITY;
    assert!(chromatogram.validate().is_err());
}

#[test]
fn native_settings_identification_attachments_append_and_validate() {
    let mut left = SpectrumSettings::default();
    let mut right = SpectrumSettings::default();
    right
        .peptide_identifications
        .push(openms::identification::PeptideIdentification::default());
    left.unify(&right).unwrap();
    assert_eq!(left.peptide_identifications.len(), 1);
    right.peptide_identifications[0].rt = Some(f64::NAN);
    let before = left.clone();
    assert!(left.unify(&right).is_err());
    assert_eq!(left, before);
}

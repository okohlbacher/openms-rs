use openms::data_structures::DateTime;
use openms::{
    kernel::{MSExperiment, MSSpectrum, Peak1D, Peak2D},
    metadata::*,
};

fn complete() -> ExperimentalSettings {
    let mut s = ExperimentalSettings::new();
    s.document.identifier = "lsid".into();
    s.document.loaded_file_path = "/source/file.mzML".into();
    s.document.loaded_file_type = openms::format::FileType::MzMl;
    s.sample.name = "bla2".into();
    s.sample.organism = "source organism".into();
    s.sample.subsamples.push(Sample {
        name: "child".into(),
        mass: -2.5,
        ..Default::default()
    });
    s.source_files.push(SourceFile {
        name: "raw".into(),
        size_mb: -1.,
        ..Default::default()
    });
    s.contacts = vec![
        ContactPerson {
            first_name: "bla17".into(),
            last_name: "blubb17".into(),
            ..Default::default()
        },
        ContactPerson {
            first_name: "bla18".into(),
            last_name: "blubb18".into(),
            ..Default::default()
        },
    ];
    s.instrument.name = "bla".into();
    s.instrument.ion_sources.push(IonSource {
        order: -3,
        ..Default::default()
    });
    s.instrument.mass_analyzers.push(MassAnalyzer {
        accuracy: -2.,
        ..Default::default()
    });
    s.instrument.ion_detectors.push(IonDetector {
        order: 4,
        ..Default::default()
    });
    s.instrument_configurations.insert(
        "astral".into(),
        Instrument {
            name: "Orbitrap Astral".into(),
            ..Default::default()
        },
    );
    s.hplc.flux = 5;
    s.hplc.gradient.add_eluent("A").unwrap();
    s.hplc.gradient.add_timepoint(0).unwrap();
    s.date_time = "02/07/2006 01:02:03".parse().unwrap();
    s.comment = "bla".into();
    s.fraction_identifier = "bla2".into();
    s.metadata.insert("label".into(), "label".into());
    s
}
fn experiment(settings: ExperimentalSettings) -> MSExperiment {
    MSExperiment {
        spectra: vec![MSSpectrum::from_peaks(vec![
            Peak1D::new(100., 1.),
            Peak1D::new(101., 3.),
            Peak1D::new(102., 1.),
        ])],
        settings,
        ..Default::default()
    }
}
#[test]
fn complete_source_defaults_literals_copy_and_diagnostic_stream() {
    let empty = ExperimentalSettings::new();
    assert_eq!(empty.date_time.get(), "0000-00-00 00:00:00");
    assert_eq!(empty.hplc.temperature, 21);
    assert!(empty.source_files.is_empty() && empty.contacts.is_empty());
    assert!(!empty.has_transport_metadata());
    assert_eq!(
        empty.to_string(),
        "-- EXPERIMENTALSETTINGS BEGIN --\n-- EXPERIMENTALSETTINGS END --\n"
    );
    let original = complete();
    let mut copy = original.checked_clone().unwrap();
    assert_eq!(copy, original);
    assert_eq!(copy.date_time.get(), "2006-02-07 01:02:03");
    assert_eq!(copy.hplc.flux, 5);
    assert_eq!(copy.contacts[0].first_name, "bla17");
    assert_eq!(copy.contacts[1].last_name, "blubb18");
    assert_eq!(copy.metadata["label"].as_str().unwrap(), "label");
    copy.instrument_configurations
        .get_mut("astral")
        .unwrap()
        .name = "Q Exactive".into();
    assert_ne!(copy, original);
    assert_eq!(
        original.instrument_configurations["astral"].name,
        "Orbitrap Astral"
    );
    copy = ExperimentalSettings::default();
    assert_eq!(copy, empty);
}
#[test]
fn source_equality_includes_every_owned_field_except_document_provenance() {
    let empty = ExperimentalSettings::new();
    let changes: Vec<fn(&mut ExperimentalSettings)> = vec![
        |s| s.sample.name = "a".into(),
        |s| s.source_files.push(SourceFile::default()),
        |s| s.contacts.push(ContactPerson::default()),
        |s| s.instrument.name = "a".into(),
        |s| {
            s.instrument_configurations
                .insert("a".into(), Instrument::default());
        },
        |s| s.hplc.flux = 5,
        |s| s.date_time = "2006-02-07 01:02:03".parse().unwrap(),
        |s| s.comment = "a".into(),
        |s| s.fraction_identifier = "a".into(),
        |s| s.document.identifier = "a".into(),
        |s| {
            s.metadata.insert("a".into(), 1.into());
        },
    ];
    for change in changes {
        let mut changed = empty.clone();
        change(&mut changed);
        assert_ne!(changed, empty);
        assert!(changed.has_transport_metadata());
    }
    let mut provenance = empty.clone();
    provenance.document.loaded_file_path = "different".into();
    provenance.document.loaded_file_type = openms::format::FileType::Mgf;
    assert_eq!(provenance, empty);
    assert!(!provenance.has_transport_metadata());
    let cloned = provenance.checked_clone().unwrap();
    assert_eq!(cloned.document.loaded_file_path, "different");
    assert_eq!(
        cloned.document.loaded_file_type,
        openms::format::FileType::Mgf
    );
}
#[test]
fn checked_copy_preserves_float_bits_partial_datetime_and_stale_gradient() {
    let mut s = complete();
    s.sample.mass = f64::from_bits(0x7ff8_0000_0000_1234);
    s.instrument.mass_analyzers[0].accuracy = -0.;
    s.date_time = DateTime::default();
    s.date_time.set_time("01:02:03").unwrap();
    s.hplc.gradient.clear_eluents();
    s.validate().unwrap();
    let cloned = s.checked_clone().unwrap();
    assert_eq!(cloned.sample.mass.to_bits(), s.sample.mass.to_bits());
    assert_eq!(
        cloned.instrument.mass_analyzers[0].accuracy.to_bits(),
        (-0.0f64).to_bits()
    );
    assert_eq!(cloned.date_time, s.date_time);
    assert_eq!(
        cloned.hplc.gradient.percentages(),
        s.hplc.gradient.percentages()
    );
    assert_eq!(cloned.source_files[0].size_mb, -1.);
}
#[test]
fn recursive_samples_and_complete_payload_obey_checked_limits_before_copy() {
    let mut s = complete();
    for limits in [
        ExperimentalSettingsLimits {
            max_records: 1,
            ..Default::default()
        },
        ExperimentalSettingsLimits {
            max_work: 0,
            ..Default::default()
        },
        ExperimentalSettingsLimits {
            max_bytes: 0,
            ..Default::default()
        },
        ExperimentalSettingsLimits {
            max_sample_depth: 1,
            ..Default::default()
        },
    ] {
        assert!(s.checked_clone_with_limits(limits).is_err());
    }
    let mut sample = Sample::default();
    for _ in 0..64 {
        sample = Sample {
            subsamples: vec![sample],
            ..Default::default()
        };
    }
    s.sample = sample;
    assert!(s.validate().is_err());
    let extended = ExperimentalSettingsLimits {
        max_sample_depth: usize::MAX,
        ..Default::default()
    };
    assert!(s.checked_clone_with_limits(extended).is_err());
    let mut wide = ExperimentalSettings::default();
    wide.instrument
        .software
        .cv_terms
        .metadata
        .insert("large".into(), "x".repeat(4096).into());
    assert!(
        wide.checked_clone_with_limits(ExperimentalSettingsLimits {
            max_bytes: 1024,
            ..Default::default()
        })
        .is_err()
    );
}
#[test]
fn one_metadata_owner_clear_swap_and_2d_replacement_preserve_source_semantics() {
    let mut e = experiment(complete());
    let before = e.settings.clone();
    e.clear(false);
    assert_eq!(e.settings, before);
    e.clear(true);
    assert_eq!(e.settings, ExperimentalSettings::default());
    assert_eq!(e.settings.hplc.temperature, 21);
    let mut e = experiment(complete());
    let label_ptr = e.settings.metadata["label"].as_str().unwrap().as_ptr();
    let previous = e.set_2d_data(&[Peak2D::new(1., 2., 3.)]).unwrap();
    assert_eq!(
        previous.settings.metadata["label"]
            .as_str()
            .unwrap()
            .as_ptr(),
        label_ptr
    );
    assert_eq!(e.settings, ExperimentalSettings::default());
    let mut other = previous;
    std::mem::swap(&mut e, &mut other);
    assert_eq!(e.settings, before);
    assert_eq!(other.spectra[0].peaks[0].mz, 2.);
}
#[test]
fn content_filters_and_whole_experiment_clones_preserve_complete_settings() {
    use openms::processing::{
        Normalizer, SpectrumFilter,
        iterative::PeakPickerIterative,
        peak_picking::PeakPickerHiRes,
        smoothing::{GaussFilter, SavitzkyGolayFilter},
    };
    let mut e = experiment(complete());
    let original = e.settings.clone();
    Normalizer::default().filter_experiment(&mut e).unwrap();
    assert_eq!(e.settings, original);
    GaussFilter::default().filter_experiment(&mut e).unwrap();
    assert_eq!(e.settings, original);
    SavitzkyGolayFilter::default()
        .filter_experiment(&mut e)
        .unwrap();
    assert_eq!(e.settings, original);
    assert_eq!(
        PeakPickerHiRes::default()
            .pick_experiment(&e)
            .unwrap()
            .experiment
            .settings,
        original
    );
    assert_eq!(
        PeakPickerIterative::default()
            .pick_experiment(&e)
            .unwrap()
            .experiment
            .settings,
        original
    );
    let mut deep = Sample::default();
    for _ in 0..64 {
        deep = Sample {
            subsamples: vec![deep],
            ..Default::default()
        };
    }
    e.settings.sample = deep;
    let peak_ptr = e.spectra[0].peaks.as_ptr();
    assert!(GaussFilter::default().filter_experiment(&mut e).is_err());
    assert_eq!(e.spectra[0].peaks.as_ptr(), peak_ptr);
}

#[cfg(feature = "mzml")]
#[test]
fn typed_run_metadata_and_units_roundtrip_with_both_mzml_writers() {
    use openms::format::mzml;
    let mut e = experiment(ExperimentalSettings::new());
    let seconds = Unit::new("UO:0000010", "second", "UO").unwrap();
    e.settings
        .metadata
        .insert("string".into(), "<&>µ\n\t\r\"".into());
    e.settings
        .metadata
        .insert("integer".into(), i64::MAX.into());
    e.settings.metadata.insert(
        "float".into(),
        MetaValue::try_from(-0.0f64)
            .unwrap()
            .with_unit(seconds)
            .unwrap(),
    );
    e.settings.document.loaded_file_path = "local-only".into();
    e.settings.document.loaded_file_type = openms::format::FileType::Mgf;
    for compressed in [false, true] {
        let mut xml = Vec::new();
        if compressed {
            mzml::write_with_numpress(&mut xml, &e, &Default::default()).unwrap();
        } else {
            mzml::write(&mut xml, &e).unwrap();
        }
        let parsed = mzml::read(xml.as_slice()).unwrap();
        assert_eq!(parsed.settings, e.settings);
        assert_eq!(
            parsed.settings.metadata["float"]
                .as_f64()
                .unwrap()
                .to_bits(),
            (-0.0f64).to_bits()
        );
        assert!(parsed.settings.document.loaded_file_path.is_empty());
    }
}
#[cfg(feature = "mzml")]
#[test]
fn unrepresented_headers_empty_and_list_metadata_fail_before_output() {
    use openms::format::mzml;
    let mut cases = vec![complete()];
    let mut negative_zero = ExperimentalSettings::new();
    negative_zero.sample.mass = -0.;
    cases.push(negative_zero);
    let mut stale = ExperimentalSettings::new();
    stale.hplc.gradient.add_eluent("A").unwrap();
    stale.hplc.gradient.clear_eluents();
    cases.push(stale);
    for value in [
        MetaValue::default(),
        MetaValue::from(vec!["a".to_string()]),
        MetaValue::from(vec![1i64]),
        MetaValue::try_from(vec![1f64]).unwrap(),
    ] {
        let mut s = ExperimentalSettings::new();
        s.metadata.insert("typed".into(), value);
        cases.push(s);
    }
    for settings in cases {
        let e = experiment(settings);
        for compressed in [false, true] {
            let mut output = b"unchanged".to_vec();
            let failed = if compressed {
                mzml::write_with_numpress(&mut output, &e, &Default::default()).map(|_| ())
            } else {
                mzml::write(&mut output, &e)
            };
            assert!(matches!(failed, Err(openms::Error::Unsupported(_))));
            assert_eq!(output, b"unchanged");
        }
    }
    let mut e = experiment(ExperimentalSettings::new());
    e.settings.sample.mass = -0.;
    e.spectra[0].peaks[0].intensity = f32::NAN;
    for compressed in [false, true] {
        let mut output = Vec::new();
        let failed = if compressed {
            mzml::write_with_numpress(&mut output, &e, &Default::default()).map(|_| ())
        } else {
            mzml::write(&mut output, &e)
        };
        assert!(matches!(failed, Err(openms::Error::Unsupported(_))));
        assert!(output.is_empty());
    }
}
#[test]
fn peak_only_transports_reject_settings_and_explicit_tic_projection_ignores_them() {
    use openms::format::{FileHandler, FileType, dta2d};
    let e = experiment(complete());
    for kind in [FileType::Mgf, FileType::Dta, FileType::Dta2d, FileType::Ms2] {
        let mut bytes = b"old".to_vec();
        assert!(FileHandler::write_experiment(&mut bytes, &e, kind).is_err());
        assert_eq!(bytes, b"old");
    }
    let mut bytes = Vec::new();
    dta2d::write_tic(&mut bytes, &e).unwrap();
    assert!(std::str::from_utf8(&bytes).unwrap().contains("\t0\t5"));
}
#[cfg(feature = "mzml")]
#[test]
fn path_loads_populate_provenance_without_blocking_rewrite_or_changing_source_equality() {
    use openms::format::{FileHandler, FileType, mzml};
    let path =
        std::env::temp_dir().join(format!("openms-settings-{}-load.mzML", std::process::id()));
    let input = experiment(ExperimentalSettings::default());
    mzml::store(&path, &input).unwrap();
    for loaded in [
        mzml::load(&path).unwrap(),
        FileHandler::load_experiment(&path, &[FileType::MzMl]).unwrap(),
    ] {
        assert_eq!(loaded.settings, input.settings);
        assert_eq!(
            loaded.settings.document.loaded_file_path,
            path.to_str().unwrap()
        );
        assert_eq!(loaded.settings.document.loaded_file_type, FileType::MzMl);
        mzml::write(Vec::new(), &loaded).unwrap();
    }
    std::fs::remove_file(path).unwrap();
}

#[test]
fn every_nested_owned_payload_is_precharged_before_checked_clone() {
    let budget = ExperimentalSettingsLimits {
        max_bytes: 32_000,
        ..Default::default()
    };
    ExperimentalSettings::default()
        .checked_clone_with_limits(budget)
        .unwrap();
    let changes: Vec<fn(&mut ExperimentalSettings)> = vec![
        |s| s.document.loaded_file_path = "x".repeat(64_000),
        |s| s.sample.organism = "x".repeat(64_000),
        |s| {
            s.sample.subsamples.push(Sample {
                comment: "x".repeat(64_000),
                ..Default::default()
            })
        },
        |s| {
            s.contacts.push(ContactPerson {
                address: "x".repeat(64_000),
                ..Default::default()
            })
        },
        |s| {
            s.source_files.push(SourceFile {
                native_id_type_accession: "x".repeat(64_000),
                ..Default::default()
            })
        },
        |s| s.instrument.customizations = "x".repeat(64_000),
        |s| {
            s.instrument_configurations
                .insert("x".repeat(64_000), Instrument::default());
        },
        |s| {
            s.instrument
                .software
                .cv_terms
                .metadata
                .insert("k".into(), MetaValue::from(vec!["x".repeat(64_000)]));
        },
        |s| {
            s.instrument.ion_sources.push(IonSource {
                metadata: [("k".into(), "x".repeat(64_000).into())].into(),
                ..Default::default()
            })
        },
        |s| {
            s.instrument.mass_analyzers.push(MassAnalyzer {
                metadata: [("k".into(), "x".repeat(64_000).into())].into(),
                ..Default::default()
            })
        },
        |s| {
            s.instrument.ion_detectors.push(IonDetector {
                metadata: [("k".into(), "x".repeat(64_000).into())].into(),
                ..Default::default()
            })
        },
        |s| s.hplc.gradient.add_eluent(&"x".repeat(64_000)).unwrap(),
        |s| {
            s.metadata.insert(
                "k".into(),
                MetaValue::from(1)
                    .with_unit(Unit::new("UO:0000010", "x".repeat(64_000), "UO").unwrap())
                    .unwrap(),
            );
        },
    ];
    for change in changes {
        let mut settings = ExperimentalSettings::default();
        change(&mut settings);
        assert!(settings.checked_clone_with_limits(budget).is_err());
    }
}

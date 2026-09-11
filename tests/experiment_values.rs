// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::metadata::*;
use std::collections::{BTreeMap, hash_map::DefaultHasher};
use std::hash::{Hash, Hasher};

fn hash(value: &impl Hash) -> u64 {
    let mut state = DefaultHasher::new();
    value.hash(&mut state);
    state.finish()
}
fn meta() -> MetaInfo {
    BTreeMap::from([("label".into(), "label".into())])
}

#[test]
fn all_154_source_literal_names_in_order_and_exact_parsing() {
    let mut counts = BTreeMap::<&str, usize>::new();
    let mut lengths = BTreeMap::new();
    for row in include_str!("data/experiment_value_names.tsv")
        .lines()
        .skip(1)
    {
        let fields: Vec<_> = row.split('\t').collect();
        assert_eq!(fields.len(), 6);
        let (family, index, label) = (fields[0], fields[1].parse::<usize>().unwrap(), fields[3]);
        let count = counts.entry(family).or_default();
        assert_eq!(index, *count);
        *count += 1;
        macro_rules! check {
            ($ty:ty) => {{
                let value = <$ty>::ALL[index];
                assert_eq!(value.name(), label, "{}:{}", fields[4], fields[5]);
                assert_eq!(label.parse::<$ty>().unwrap(), value);
                assert_eq!(value.to_string(), label);
                assert!(format!(" {label}").parse::<$ty>().is_err());
                assert!(format!("{label}\n").parse::<$ty>().is_err());
                assert!(fields[2].parse::<$ty>().is_err() || fields[2] == label);
                <$ty>::ALL.len()
            }};
        }
        let length = match family {
            "SampleState" => check!(SampleState),
            "InletType" => check!(InletType),
            "IonizationMethod" => check!(IonizationMethod),
            "Polarity" => check!(Polarity),
            "AnalyzerType" => check!(AnalyzerType),
            "ResolutionMethod" => check!(ResolutionMethod),
            "ResolutionType" => check!(ResolutionType),
            "ScanDirection" => check!(ScanDirection),
            "ScanLaw" => check!(ScanLaw),
            "ReflectronState" => check!(ReflectronState),
            "DetectorType" => check!(DetectorType),
            "DetectorAcquisitionMode" => check!(DetectorAcquisitionMode),
            "IonOpticsType" => check!(IonOpticsType),
            _ => panic!("unknown fixture family"),
        };
        lengths.insert(family, length);
    }
    assert_eq!(counts, lengths);
    assert_eq!(counts.values().sum::<usize>(), 154);
    assert_eq!(
        counts.values().copied().collect::<Vec<_>>(),
        vec![15, 5, 22, 20, 12, 52, 3, 4, 4, 3, 7, 3, 4]
    );
    assert!("Linear".parse::<ScanLaw>().is_err());
    assert!("Membrane separator".parse::<InletType>().is_err());
    assert!("unknown".parse::<SampleState>().is_err());
    assert!("Unknown".parse::<Polarity>().is_err());
    assert!("SIZE_OF_TYPE".parse::<DetectorType>().is_err());
}

#[test]
fn complete_source_defaults() {
    let s = Sample::default();
    assert!(
        s.name.is_empty() && s.organism.is_empty() && s.number.is_empty() && s.comment.is_empty()
    );
    assert_eq!(s.state, SampleState::Unknown);
    assert_eq!(
        [
            s.mass.to_bits(),
            s.volume.to_bits(),
            s.concentration.to_bits()
        ],
        [0; 3]
    );
    assert!(s.subsamples.is_empty() && s.metadata.is_empty());
    let source = IonSource::default();
    assert_eq!(source.inlet_type, InletType::Unknown);
    assert_eq!(source.ionization_method, IonizationMethod::Unknown);
    assert_eq!(source.polarity, Polarity::Unknown);
    assert_eq!(source.order, 0);
    assert!(source.metadata.is_empty());
    let analyzer = MassAnalyzer::default();
    assert_eq!(analyzer.analyzer_type, AnalyzerType::Unknown);
    assert_eq!(analyzer.resolution_method, ResolutionMethod::Unknown);
    assert_eq!(analyzer.resolution_type, ResolutionType::Unknown);
    assert_eq!(analyzer.scan_direction, ScanDirection::Unknown);
    assert_eq!(analyzer.scan_law, ScanLaw::Unknown);
    assert_eq!(analyzer.reflectron_state, ReflectronState::Unknown);
    assert_eq!(
        [
            analyzer.resolution,
            analyzer.accuracy,
            analyzer.scan_rate,
            analyzer.scan_time,
            analyzer.tof_total_path_length,
            analyzer.isolation_width,
            analyzer.magnetic_field_strength
        ]
        .map(f64::to_bits),
        [0; 7]
    );
    assert_eq!((analyzer.final_ms_exponent, analyzer.order), (0, 0));
    assert!(analyzer.metadata.is_empty());
    let detector = IonDetector::default();
    assert_eq!(detector.detector_type, DetectorType::Unknown);
    assert_eq!(detector.acquisition_mode, DetectorAcquisitionMode::Unknown);
    assert_eq!(
        [
            detector.resolution.to_bits(),
            detector.adc_sampling_frequency.to_bits()
        ],
        [0; 2]
    );
    assert_eq!(detector.order, 0);
    assert!(detector.metadata.is_empty());
    let instrument = Instrument::default();
    assert!(
        instrument.name.is_empty()
            && instrument.vendor.is_empty()
            && instrument.model.is_empty()
            && instrument.customizations.is_empty()
    );
    assert!(
        instrument.ion_sources.is_empty()
            && instrument.mass_analyzers.is_empty()
            && instrument.ion_detectors.is_empty()
    );
    assert_eq!(instrument.software, Software::default());
    assert_eq!(instrument.ion_optics, IonOpticsType::Unknown);
    assert!(instrument.metadata.is_empty());
}

#[test]
fn source_sample_literals_order_and_owned_copy() {
    // Sample_test.cpp:191–234, with its explicit ordered subsamples case165–179.
    let original = Sample {
        name: "TTEST".into(),
        organism: "TTEST2".into(),
        number: "Sample4711".into(),
        comment: "Sample Description".into(),
        state: SampleState::Liquid,
        mass: 4711.2,
        volume: 4711.3,
        concentration: 4711.4,
        metadata: BTreeMap::from([("label".into(), "horse".into())]),
        subsamples: vec![
            Sample {
                name: "2".into(),
                ..Default::default()
            },
            Sample {
                name: "3".into(),
                ..Default::default()
            },
        ],
    };
    let mut copied = original.clone();
    assert_eq!(copied, original);
    assert_eq!(
        (copied.mass, copied.volume, copied.concentration),
        (4711.2, 4711.3, 4711.4)
    );
    assert_eq!(
        (
            &*copied.name,
            &*copied.organism,
            &*copied.number,
            &*copied.comment
        ),
        ("TTEST", "TTEST2", "Sample4711", "Sample Description")
    );
    assert_eq!(
        copied
            .subsamples
            .iter()
            .map(|s| s.name.as_str())
            .collect::<Vec<_>>(),
        ["2", "3"]
    );
    copied.subsamples[0].subsamples.push(Sample {
        name: "nested Δ".into(),
        metadata: meta(),
        ..Default::default()
    });
    copied.metadata.clear();
    assert_ne!(copied, original);
    assert!(original.subsamples[0].subsamples.is_empty());
    assert_eq!(original.metadata["label"].as_str().unwrap(), "horse");
    let mut moved = copied;
    let mut empty = Sample::default();
    std::mem::swap(&mut moved, &mut empty);
    assert_eq!(moved, Sample::default());
    assert_ne!(empty, moved);
}

fn source_analyzer() -> MassAnalyzer {
    // MassAnalyzer_test.cpp:205–245: all fields are distinct source literals.
    MassAnalyzer {
        analyzer_type: AnalyzerType::Quadrupole,
        accuracy: 47.11,
        final_ms_exponent: 47,
        isolation_width: 47.12,
        magnetic_field_strength: 47.13,
        reflectron_state: ReflectronState::On,
        resolution: 47.14,
        resolution_method: ResolutionMethod::Fwhm,
        resolution_type: ResolutionType::Constant,
        scan_direction: ScanDirection::Up,
        scan_law: ScanLaw::Linear,
        scan_rate: 47.15,
        scan_time: 47.16,
        tof_total_path_length: 47.17,
        metadata: meta(),
        order: 45,
    }
}

#[test]
fn source_component_literals_copy_and_hash_contract() {
    // IonSource_test.cpp:82–95 and IonDetector_test.cpp:92–107.
    let source = IonSource {
        inlet_type: InletType::Direct,
        ionization_method: IonizationMethod::Esi,
        polarity: Polarity::Positive,
        order: 45,
        metadata: meta(),
    };
    assert_eq!(source.clone(), source);
    assert_eq!(hash(&source), hash(&source.clone()));
    let mut changed = source.clone();
    changed.metadata.clear();
    assert_ne!(changed, source);
    assert_ne!(hash(&changed), hash(&source));
    let analyzer = source_analyzer();
    assert_eq!(analyzer.clone(), analyzer);
    assert_eq!(
        [
            analyzer.accuracy,
            analyzer.isolation_width,
            analyzer.magnetic_field_strength,
            analyzer.resolution,
            analyzer.scan_rate,
            analyzer.scan_time,
            analyzer.tof_total_path_length
        ],
        [47.11, 47.12, 47.13, 47.14, 47.15, 47.16, 47.17]
    );
    assert_eq!((analyzer.final_ms_exponent, analyzer.order), (47, 45));
    let mut changed = analyzer.clone();
    changed.metadata.clear();
    assert_ne!(changed, analyzer);
    assert_eq!(hash(&changed), hash(&analyzer)); // Source specialization omits metadata.
    let detector = IonDetector {
        resolution: 47.11,
        adc_sampling_frequency: 47.21,
        acquisition_mode: DetectorAcquisitionMode::PulseCounting,
        detector_type: DetectorType::ElectronMultiplier,
        metadata: meta(),
        order: 45,
    };
    assert_eq!(detector.clone(), detector);
    assert_eq!(
        (
            detector.resolution,
            detector.adc_sampling_frequency,
            detector.order
        ),
        (47.11, 47.21, 45)
    );
    let mut changed = detector.clone();
    changed.metadata.clear();
    assert_ne!(changed, detector);
    assert_eq!(hash(&changed), hash(&detector));
}

#[test]
fn source_instrument_literals_family_order_software_and_metadata_copy() {
    // Instrument_test.cpp:127–149,213–247; independent unsorted order edge.
    let instrument = Instrument {
        name: "Name".into(),
        vendor: "Vendor".into(),
        model: "Model".into(),
        customizations: "Customizations".into(),
        ion_sources: vec![
            IonSource {
                order: 20,
                ..Default::default()
            },
            IonSource {
                order: -1,
                ..Default::default()
            },
        ],
        mass_analyzers: vec![
            MassAnalyzer {
                scan_time: 47.11,
                ..Default::default()
            },
            MassAnalyzer {
                scan_time: 47.12,
                ..Default::default()
            },
        ],
        ion_detectors: vec![IonDetector::default()],
        software: Software {
            name: "sn".into(),
            ..Default::default()
        },
        ion_optics: IonOpticsType::Reflectron,
        metadata: meta(),
    };
    let mut copy = instrument.clone();
    assert_eq!(copy, instrument);
    assert_eq!(
        copy.mass_analyzers
            .iter()
            .map(|x| x.scan_time)
            .collect::<Vec<_>>(),
        [47.11, 47.12]
    );
    assert_eq!(
        copy.ion_sources.iter().map(|x| x.order).collect::<Vec<_>>(),
        [20, -1]
    );
    assert_eq!(copy.software.name, "sn");
    copy.software
        .cv_terms
        .metadata
        .insert("changed".into(), 1i32.into());
    copy.ion_sources.swap(0, 1);
    copy.mass_analyzers[0]
        .metadata
        .insert("changed".into(), 2i32.into());
    assert_ne!(copy, instrument);
    assert!(instrument.software.cv_terms.metadata.is_empty());
    assert!(instrument.mass_analyzers[0].metadata.is_empty());
    assert_eq!(instrument.ion_sources[0].order, 20);
}

#[test]
fn every_owned_field_participates_in_equality() {
    macro_rules! changes {
        ($ty:ty; $($field:ident = $value:expr),+ $(,)?) => {{
            let original = <$ty>::default();
            $(let mut changed = original.clone(); changed.$field = $value;
              assert_ne!(original, changed, stringify!($field));)+
        }};
    }
    changes!(Sample; name="x".into(), organism="x".into(), number="x".into(), comment="x".into(), state=SampleState::Gas,
        mass=1.0, volume=1.0, concentration=1.0, subsamples=vec![Sample::default()], metadata=meta());
    changes!(IonSource; inlet_type=InletType::Direct, ionization_method=IonizationMethod::Esi, polarity=Polarity::Negative, order=i32::MIN, metadata=meta());
    changes!(MassAnalyzer; analyzer_type=AnalyzerType::Lit, resolution_method=ResolutionMethod::Baseline, resolution_type=ResolutionType::Proportional,
        scan_direction=ScanDirection::Down, scan_law=ScanLaw::Quadratic, reflectron_state=ReflectronState::None,
        resolution=1.0, accuracy=1.0, scan_rate=1.0, scan_time=1.0, tof_total_path_length=1.0,
        isolation_width=1.0, magnetic_field_strength=1.0, final_ms_exponent=i32::MIN, order=i32::MAX, metadata=meta());
    changes!(IonDetector; detector_type=DetectorType::Channeltron, acquisition_mode=DetectorAcquisitionMode::Tdc,
        resolution=1.0, adc_sampling_frequency=1.0, order=i32::MIN, metadata=meta());
    changes!(Instrument; name="x".into(), vendor="x".into(), model="x".into(), customizations="x".into(),
        ion_sources=vec![IonSource::default()], mass_analyzers=vec![MassAnalyzer::default()], ion_detectors=vec![IonDetector::default()],
        software=Software { name: "x".into(), ..Default::default() }, ion_optics=IonOpticsType::Reflectron, metadata=meta());
}

#[test]
fn ieee_scalar_storage_equality_and_zero_hashing_are_preserved() {
    // Independent domain cases: source setters assign directly, equality uses ==.
    let cases = [
        0.0,
        -0.0,
        -2.5,
        f64::MIN,
        f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::from_bits(0x7ff8_0000_0000_0042),
        f64::from_bits(0xfff0_0000_0000_0001),
    ];
    for x in cases {
        let sample = Sample {
            mass: x,
            volume: x,
            concentration: x,
            ..Default::default()
        };
        let copied = sample.clone();
        assert_eq!(
            [
                copied.mass.to_bits(),
                copied.volume.to_bits(),
                copied.concentration.to_bits()
            ],
            [x.to_bits(); 3]
        );
        assert_eq!(sample == copied, !x.is_nan());
        let analyzer = MassAnalyzer {
            resolution: x,
            accuracy: x,
            scan_rate: x,
            scan_time: x,
            tof_total_path_length: x,
            isolation_width: x,
            magnetic_field_strength: x,
            ..Default::default()
        };
        let copied = analyzer.clone();
        assert_eq!(
            [
                copied.resolution,
                copied.accuracy,
                copied.scan_rate,
                copied.scan_time,
                copied.tof_total_path_length,
                copied.isolation_width,
                copied.magnetic_field_strength
            ]
            .map(f64::to_bits),
            [x.to_bits(); 7]
        );
        assert_eq!(analyzer == copied, !x.is_nan());
        let detector = IonDetector {
            resolution: x,
            adc_sampling_frequency: x,
            ..Default::default()
        };
        assert_eq!(
            [
                detector.clone().resolution.to_bits(),
                detector.clone().adc_sampling_frequency.to_bits()
            ],
            [x.to_bits(); 2]
        );
        assert_eq!(detector == detector.clone(), !x.is_nan());
        let instrument = Instrument {
            mass_analyzers: vec![analyzer],
            ion_detectors: vec![detector],
            ..Default::default()
        };
        assert_eq!(instrument == instrument.clone(), !x.is_nan());
    }
    assert_eq!(
        hash(&MassAnalyzer::default()),
        hash(&MassAnalyzer {
            resolution: -0.0,
            ..Default::default()
        })
    );
    assert_eq!(
        hash(&IonDetector::default()),
        hash(&IonDetector {
            resolution: -0.0,
            ..Default::default()
        })
    );
    let nested = Sample {
        subsamples: vec![Sample {
            mass: f64::NAN,
            ..Default::default()
        }],
        ..Default::default()
    };
    assert_ne!(nested, nested.clone());
}

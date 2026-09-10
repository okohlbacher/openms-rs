// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Source: TheoreticalSpectrumGenerator.cpp/test.cpp, pinned OpenMS4-core 7c029e8.

use openms::Error;
use openms::chemistry::{
    AASequence, TheoreticalIsotopeModel, TheoreticalSpectrumGenerator as Generator,
};
use openms::kernel::{DataArray, MSSpectrum, Peak1D, SpectrumType};
use openms::metadata::ActivationMethod;

fn peptide(text: &str) -> AASequence {
    text.parse().unwrap()
}
fn names(spectrum: &MSSpectrum) -> &[String] {
    &spectrum.string_data_arrays[0].data
}
fn immonium_settings() -> Generator {
    Generator {
        ion_series: Vec::new(),
        add_abundant_immonium_ions: true,
        add_metainfo: true,
        sort_by_position: false,
        ..Generator::default()
    }
}

#[test]
fn source_cid_and_hcid_factory_mass_goldens() {
    for (method, sequence, charge, expected) in [
        (
            ActivationMethod::Cid,
            "HFYLWCP",
            1,
            vec![
                116.0706, 219.0797, 285.1346, 405.1591, 448.1979, 518.2431, 561.2819, 681.3064,
                747.3613, 828.3749, 850.3704,
            ],
        ),
        (
            ActivationMethod::Hcid,
            "PEP",
            3,
            vec![
                58.5389, 100.0574, 114.0549, 116.0706, 123.0602, 199.1077, 227.1026, 245.1131,
            ],
        ),
    ] {
        let spectrum =
            Generator::generate_for_activation(method, &peptide(sequence), charge).unwrap();
        assert_eq!(spectrum.len(), expected.len());
        for (peak, mz) in spectrum.peaks.iter().zip(expected) {
            assert!((peak.mz - mz).abs() < 0.0001, "{} != {mz}", peak.mz);
            assert_eq!(peak.intensity, 1.0);
        }
        assert_eq!(spectrum.ms_level, 2);
        assert_eq!(spectrum.spectrum_type, SpectrumType::Centroid);
        assert!(spectrum.string_data_arrays.is_empty());
        assert!(spectrum.integer_data_arrays.is_empty());
    }
}

#[test]
fn all_activation_variants_and_source_precursor_inference_are_explicit() {
    let p = peptide("HFYLWCP");
    for &method in ActivationMethod::ALL {
        let count = match method {
            ActivationMethod::Cid => 11,
            ActivationMethod::Hcid | ActivationMethod::Hcd => 16,
            ActivationMethod::Ecd | ActivationMethod::Etd => 17,
            ActivationMethod::Etcid | ActivationMethod::Ethcd => 45,
            _ => {
                assert!(matches!(
                    Generator::generate_for_activation(method, &p, 2),
                    Err(Error::Unsupported(_))
                ));
                continue;
            }
        };
        let low = Generator::generate_for_activation(method, &p, 2).unwrap();
        assert_eq!(low.len(), count);
        assert_eq!(low.precursors[0].charge, 2);
        assert_eq!(low.precursors[0].mz, p.mz(2).unwrap());
        for charge in [0, 1] {
            assert_eq!(
                Generator::generate_for_activation(method, &p, charge).unwrap(),
                low
            );
        }
        let high = Generator::generate_for_activation(method, &p, 3).unwrap();
        assert_eq!(high.len(), 2 * count);
        assert_eq!(high.precursors[0].charge, 3);
        assert_eq!(high.precursors[0].mz, p.mz(3).unwrap());
        assert_eq!(
            Generator::generate_for_activation(method, &p, u16::MAX).unwrap(),
            high
        );
    }
}

#[test]
fn factory_preserves_empty_and_checked_unresolved_peptide_behavior() {
    assert_eq!(
        Generator::generate_for_activation(ActivationMethod::Cid, &AASequence::default(), 0)
            .unwrap(),
        MSSpectrum::default()
    );
    assert!(matches!(
        Generator::generate_for_activation(ActivationMethod::Cid, &peptide("AXA"), 2),
        Err(Error::Unsupported(_))
    ));
    assert!(Generator::generate_for_activation(ActivationMethod::Etd, &peptide("A"), 2).is_err());
    assert!(matches!(
        Generator::generate_for_activation(ActivationMethod::Sori, &AASequence::default(), 0),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn fixed_immonium_peaks_have_source_order_unit_intensity_and_charge_one() {
    let mut generator = immonium_settings();
    generator.intensities.b = 29.0;
    generator.relative_loss_intensity = 0.0;
    generator.isotope_model = TheoreticalIsotopeModel::Coarse { max_peaks: 3 };
    let spectrum = generator
        .generate(&peptide("WHYLPCFWHYLPCF"), 2, 3, None)
        .unwrap();
    let masses = [
        70.0656, 76.0221, 86.09698, 110.0718, 120.0813, 136.0762, 159.0922,
    ];
    assert_eq!(spectrum.len(), 7);
    for (peak, expected) in spectrum.peaks.iter().zip(masses) {
        assert_eq!(peak.mz, expected);
        assert_eq!(peak.intensity, 1.0);
    }
    assert_eq!(
        names(&spectrum),
        ["iP+", "iC+", "iL/I+", "iH+", "iF+", "iY+", "iW+"]
    );
    assert_eq!(spectrum.integer_data_arrays[0].data, [1; 7]);
    assert_eq!(spectrum.precursors[0].charge, 4);
    spectrum.validate().unwrap();
    generator.add_metainfo = false;
    let without_annotations = generator.generate(&peptide("WHYLPCF"), 2, 3, None).unwrap();
    assert_eq!(without_annotations.peaks, spectrum.peaks);
    assert!(without_annotations.string_data_arrays.is_empty());
    assert!(without_annotations.integer_data_arrays.is_empty());
}

#[test]
fn immonium_eligibility_uses_unmodified_residues_but_ignores_terminal_slots() {
    let generator = immonium_settings();
    for p in ["I", "A", "C(Carbamidomethyl)L[+12.3456789]"] {
        assert!(
            generator
                .generate(&peptide(p), 1, 1, None)
                .unwrap()
                .is_empty()
        );
    }
    let free_and_modified = generator
        .generate(&peptide("C(Carbamidomethyl)CL[+12.3456789]L"), 1, 1, None)
        .unwrap();
    assert_eq!(names(&free_and_modified), ["iC+", "iL/I+"]);
    let terminal = generator
        .generate(&peptide(".(Acetyl)HFYLWCP.(Amidated)"), 1, 1, None)
        .unwrap();
    assert_eq!(terminal.len(), 7);
    let h = generator.generate(&peptide("H"), 1, 1, None).unwrap();
    assert_eq!(names(&h), ["iH+"]);
    let mut disabled = generator.clone();
    disabled.add_abundant_immonium_ions = false;
    assert!(
        disabled
            .generate(&peptide("HFYLWCP"), 1, 1, None)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn immonium_append_preserves_aligned_annotations_and_is_atomic_on_other_arrays() {
    let generator = immonium_settings();
    let mut spectrum = MSSpectrum {
        peaks: vec![Peak1D::new(900.0, 2.0)],
        string_data_arrays: vec![DataArray::new("IonNames", vec!["existing".into()])],
        integer_data_arrays: vec![DataArray::new("Charges", vec![3])],
        name: "retained".into(),
        ..MSSpectrum::default()
    };
    generator
        .append_to(&mut spectrum, &peptide("H"), 2, 2, None)
        .unwrap();
    assert_eq!(spectrum.name, "retained");
    assert_eq!(names(&spectrum), ["existing", "iH+"]);
    assert_eq!(spectrum.integer_data_arrays[0].data, [3, 1]);
    assert_eq!(
        spectrum.peaks,
        [Peak1D::new(900.0, 2.0), Peak1D::new(110.0718, 1.0)]
    );
    spectrum.validate().unwrap();
    spectrum
        .float_data_arrays
        .push(DataArray::new("profile", vec![1.0, 2.0]));
    let before = spectrum.clone();
    assert!(
        generator
            .append_to(&mut spectrum, &peptide("H"), 1, 1, None)
            .is_err()
    );
    assert_eq!(spectrum, before);
}

use openms::chemistry::{
    IMSElement, IMSIsotopeDistribution as Distribution, IMSIsotopeOptions as Options,
    IMSIsotopePeak as Peak,
};

fn hydrogen() -> Distribution {
    Distribution::from_peaks(
        vec![
            Peak {
                mass: 0.0078250319,
                abundance: 0.999885,
            },
            Peak {
                mass: 0.01410178,
                abundance: 0.000115,
            },
            Peak {
                mass: 0.01604927,
                abundance: 0.0,
            },
        ],
        1,
    )
    .unwrap()
}

#[test]
fn source_constructors_and_hydrogen_oxygen_mass_literals() {
    let default = IMSElement::default();
    assert_eq!(default.name(), "");
    assert_eq!(default.sequence(), "");
    assert!(default.isotope_distribution().is_empty());
    assert!(default.mass(0).is_err());
    let h = IMSElement::from_distribution("H", hydrogen()).unwrap();
    assert_eq!(h.name(), "H");
    assert_eq!(h.sequence(), "H");
    assert_eq!(h.nominal_mass(), 1);
    assert_eq!(h.mass(0).unwrap(), 1.0078250319);
    assert_eq!(h.mass(1).unwrap(), 2.01410178);
    assert_eq!(h.mass(2).unwrap(), 3.01604927);
    let oxygen = IMSElement::from_mass("O", 15.9994).unwrap();
    assert_eq!(oxygen.mass(0).unwrap(), 15.9994);
    assert_eq!(oxygen.nominal_mass(), 0);
    let oxygen_nominal = IMSElement::new("O", 16).unwrap();
    assert_eq!(oxygen_nominal.nominal_mass(), 16);
    assert!(oxygen_nominal.mass(0).is_err());
}

#[test]
fn source_independent_labels_replacement_and_complete_value_equality() {
    let mut h = IMSElement::from_distribution("H", hydrogen()).unwrap();
    let original = h.clone();
    h.set_name("D").unwrap();
    assert_eq!(h.sequence(), "H");
    assert_ne!(h, original);
    h.set_name("H").unwrap();
    h.set_sequence("H2").unwrap();
    assert_eq!(h.name(), "H");
    h.set_sequence("H").unwrap();
    assert_eq!(h, original);
    let mut peaks = hydrogen().peaks().to_vec();
    peaks.push(Peak {
        mass: 0.03604927,
        abundance: 0.0,
    });
    h.set_isotope_distribution(Distribution::from_peaks(peaks, 1).unwrap());
    assert_ne!(h, original); // hidden zero-abundance tail is part of equality
    assert_eq!(h.average_mass().unwrap(), original.average_mass().unwrap());
}

#[test]
fn source_electron_constant_and_full_signed_charge_range() {
    let h = IMSElement::from_distribution("H", hydrogen()).unwrap();
    assert_eq!(IMSElement::ELECTRON_MASS_IN_U, 0.00054858);
    for electrons in [i32::MIN, -2, 0, 1, 2, i32::MAX] {
        assert_eq!(
            h.ion_mass(electrons).unwrap(),
            1.0078250319 - f64::from(electrons) * 0.00054858
        );
    }
    assert!(IMSElement::default().ion_mass(1).is_err());
    assert_eq!(
        h.average_mass().unwrap(),
        1.0078250319 * 0.999885 + 2.01410178 * 0.000115
    );
}

#[test]
fn exact_labels_and_classic_stream_layout() {
    let h = IMSElement::from_distribution("H", hydrogen()).unwrap();
    assert_eq!(
        h.to_text(Options::default()).unwrap(),
        "name:\tH\nsequence:\tH\nisotope distribution:\n\n"
    );
    assert_eq!(
        h.to_text(Options {
            size: 3,
            abundances_sum_error: f64::NAN
        })
        .unwrap(),
        "name:\tH\nsequence:\tH\nisotope distribution:\n1.00783 0.999885\n2.0141 0.000115\n3.01605 0\n\n"
    );
    let strange = IMSElement::new("λ\n\0", 0).unwrap();
    assert_eq!(strange.name(), "λ\n\0");
}

#[test]
fn invalid_values_and_bounded_label_changes_are_atomic() {
    assert!(IMSElement::from_mass("H", f64::NAN).is_err());
    let excessive = "x".repeat(1024 * 1024 + 1);
    assert!(IMSElement::new(&excessive, 0).is_err());
    let mut element = IMSElement::new("H", 1).unwrap();
    let before = element.clone();
    assert!(element.set_name(&excessive).is_err());
    assert_eq!(element, before);
    assert!(element.set_sequence(&excessive).is_err());
    assert_eq!(element, before);
}

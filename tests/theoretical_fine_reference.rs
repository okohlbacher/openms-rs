// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Literal fine-spectrum goldens from TheoreticalSpectrumGenerator_test.cpp
//! at OpenMS4-core 7c029e8, plus independent small-formula multinomial oracles.
//! These checks do not invoke C++ or use the native fine generator as an oracle.

use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationRecord, ModificationsDB, ResidueModification,
    TheoreticalIonSeries as Ion, TheoreticalIsotopeModel as Isotopes,
    TheoreticalSpectrumGenerator as Generator,
};
use openms::kernel::{DataArray, MSSpectrum, Peak1D, Precursor};

fn peptide(text: &str) -> AASequence {
    text.parse().unwrap()
}
fn generator(unexplained_probability: f64) -> Generator {
    Generator {
        isotope_model: Isotopes::Fine {
            unexplained_probability,
        },
        add_metainfo: true,
        ..Default::default()
    }
}
fn names(spectrum: &MSSpectrum) -> &[String] {
    &spectrum
        .string_data_arrays
        .iter()
        .find(|array| array.name == "IonNames")
        .unwrap()
        .data
}
fn charges(spectrum: &MSSpectrum) -> &[i32] {
    &spectrum
        .integer_data_arrays
        .iter()
        .find(|array| array.name == "Charges")
        .unwrap()
        .data
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.13} != {expected:.13} (tolerance {tolerance})"
    );
}

#[test]
fn literal_source_fine_y_ladder_goldens_and_unexplained_probability() {
    let mut g = generator(0.05);
    g.ion_series = vec![Ion::Y];
    // Pinned class-test lines 681–705. Its four approximate monoisotopic
    // literals retain the source's 0.001 tolerance; the others have more digits.
    let expected = [
        78.54206,
        79.04424117545,
        107.05279,
        107.5549732233,
        185.10335,
        185.6023689147,
        185.6055289147,
        263.15390,
        263.6529246061,
        263.6560846061,
    ];
    let s = g.generate(&peptide("ARRGH"), 2, 2, None).unwrap();
    assert_eq!(s.len(), expected.len());
    for (peak, mass) in s.peaks.iter().zip(expected) {
        near(peak.mz, mass, 0.001);
    }
    assert_eq!(
        names(&s),
        ["y1", "y1", "y2", "y2", "y3", "y3", "y3", "y4", "y4", "y4"]
    );
    assert_eq!(charges(&s), [2; 10]);

    g.isotope_model = Isotopes::Fine {
        unexplained_probability: 0.20,
    };
    let s = g.generate(&peptide("ARRGH"), 2, 2, None).unwrap();
    assert_eq!(s.len(), 5);
    for (peak, mass) in
        s.peaks
            .iter()
            .zip([78.54206, 107.05279, 185.10335, 263.15390, 263.6560846061])
    {
        near(peak.mz, mass, 0.001);
    }
}

#[test]
fn literal_source_fine_losses_and_precursors() {
    let mut g = generator(0.05);
    g.ion_series = vec![Ion::Y];
    g.add_losses = true;
    let s = g.generate(&peptide("ARRGH"), 1, 2, None).unwrap();
    assert_eq!(s.len(), 50);
    for ((peak, mass), intensity) in s
        .peaks
        .iter()
        .zip([78.5426, 79.0442, 107.0532, 107.5549])
        .zip([0.921514, 0.0598011, 0.896088, 0.0775347])
    {
        near(peak.mz, mass, 0.001);
        near(f64::from(peak.intensity), intensity, 0.000001);
    }
    for (peak, mass) in s.peaks[s.len() - 5..]
        .iter()
        .zip([509.271, 509.277, 525.301, 526.298, 526.304])
    {
        near(peak.mz, mass, 0.001);
    }
    assert!(names(&s).iter().any(|name| name == "y3-C1H2N1O1++"));

    g.ion_series.clear();
    g.add_precursor_peaks = true;
    let s = g.generate(&peptide("ARRGH"), 2, 2, None).unwrap();
    assert_eq!(s.len(), 12);
    const PROTON: f64 = 1.007_276_466_771;
    near(s.peaks[0].mz, (578.32698 + PROTON) / 2.0, 0.001);
    near(s.peaks[1].mz, (579.31100 + PROTON) / 2.0, 0.001);
    near(s.peaks[11].mz, (598.344_813_339_43 + PROTON) / 2.0, 0.001);
    assert!(
        names(&s)
            .iter()
            .all(|name| ["[M+2H]++", "[M+2H-H2O]++", "[M+2H-NH3]++"].contains(&name.as_str()))
    );
    assert_eq!(charges(&s), [2; 12]);
}

fn choose(n: u32, k: u32) -> f64 {
    (0..k).fold(1.0, |value, j| value * f64::from(n - j) / f64::from(j + 1))
}
fn binary_probability(n: u32, heavy: u32, light_p: f32, heavy_p: f32) -> f64 {
    choose(n, heavy)
        * f64::from(light_p).powi((n - heavy) as i32)
        * f64::from(heavy_p).powi(heavy as i32)
}

// Enumerate the entire small C/H/N/O Cartesian product using closed-form
// binomial/trinomial probabilities, unlike production's lazy search. All masses
// and probabilities are literal ElementDB.cpp inputs. Source ElementDB stores
// abundances as f32, and IsoSpecWrapper stores each resulting probability as f32
// before its coverage trim. No normalization is applied at either step.
fn chno_envelope(counts: [u32; 4], unexplained: f64) -> Vec<(f64, f32)> {
    let [c, h, n, o] = counts;
    assert!(c <= 10 && h <= 20 && n <= 4 && o <= 5);
    let mut all = Vec::new();
    for c13 in 0..=c {
        for h2 in 0..=h {
            for n15 in 0..=n {
                for o17 in 0..=o {
                    for o18 in 0..=o - o17 {
                        let mass = f64::from(c - c13) * 12.0
                            + f64::from(c13) * 13.003_355
                            + f64::from(h - h2) * 1.007_825_031_9
                            + f64::from(h2) * 2.014_101_78
                            + f64::from(n - n15) * 14.003_074
                            + f64::from(n15) * 15.000_109
                            + f64::from(o - o17 - o18) * 15.994_915
                            + f64::from(o17) * 16.999_132
                            + f64::from(o18) * 17.999_169;
                        let probability = binary_probability(c, c13, 0.9893, 0.0107)
                            * binary_probability(h, h2, 0.999885, 0.000115)
                            * binary_probability(n, n15, 0.99632, 0.00368)
                            * choose(o, o17)
                            * choose(o - o17, o18)
                            * f64::from(0.99757_f32).powi((o - o17 - o18) as i32)
                            * f64::from(0.00038_f32).powi(o17 as i32)
                            * f64::from(0.00205_f32).powi(o18 as i32);
                        all.push((mass, probability as f32));
                    }
                }
            }
        }
    }
    all.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.total_cmp(&b.0)));
    let mut sum = 0.0;
    let mut selected = Vec::new();
    for peak in all {
        if sum >= 1.0 - unexplained {
            break;
        }
        sum += f64::from(peak.1);
        selected.push(peak);
    }
    assert!(sum >= 1.0 - unexplained);
    selected.sort_by(|a, b| a.0.total_cmp(&b.0));
    selected
}

fn assert_group(
    spectrum: &MSSpectrum,
    name: &str,
    charge: i32,
    counts: [u32; 4],
    unexplained: f64,
    intensity: f64,
) {
    let actual: Vec<_> = spectrum
        .peaks
        .iter()
        .zip(names(spectrum))
        .zip(charges(spectrum))
        .filter(|((_, label), z)| label.as_str() == name && **z == charge)
        .map(|((peak, _), _)| peak)
        .collect();
    let expected = chno_envelope(counts, unexplained);
    assert_eq!(actual.len(), expected.len(), "{name}/{charge}");
    for (peak, (mass, probability)) in actual.into_iter().zip(expected) {
        near(peak.mz, mass / f64::from(charge), 1e-9);
        assert_eq!(
            peak.intensity.to_bits(),
            ((intensity * f64::from(probability)) as f32).to_bits(),
            "{name}/{charge} at {}",
            peak.mz
        );
    }
}

#[test]
fn explicit_natural_hydrogen_and_charge_division_match_an_independent_full_oracle() {
    let mut g = generator(0.0001);
    g.ion_series = vec![Ion::B];
    g.add_first_prefix_ion = true;
    g.intensities.b = 0.7;
    let s = g.generate(&peptide("AA"), 1, 2, Some(4)).unwrap();
    // Alanine b1 is C3H5NO. Each charge adds natural neutral H, whose heavy
    // isotope participates, and only then is the full mass divided by charge.
    assert_group(&s, "b1", 1, [3, 6, 1, 1], 0.0001, f64::from(0.7_f32));
    assert_group(&s, "b1", 2, [3, 7, 1, 1], 0.0001, f64::from(0.7_f32));
    let two: Vec<_> = s
        .peaks
        .iter()
        .zip(charges(&s))
        .filter(|(_, charge)| **charge == 2)
        .map(|(peak, _)| peak)
        .collect();
    // Nitrogen-15 and carbon-13 remain separate lines within the nominal M+1.
    let nitrogen_shift = (15.000_109 - 14.003_074) / 2.0;
    let carbon_shift = (13.003_355 - 12.0) / 2.0;
    assert!(
        two.iter()
            .any(|peak| (peak.mz - two[0].mz - nitrogen_shift).abs() < 1e-9)
    );
    assert!(
        two.iter()
            .any(|peak| (peak.mz - two[0].mz - carbon_shift).abs() < 1e-9)
    );
    assert_eq!(s.precursors[0].charge, 4);
    near(s.precursors[0].mz, peptide("AA").mz(4).unwrap(), 1e-12);
    assert!(names(&s).iter().all(|name| name == "b1"));
}

#[test]
fn terminal_formulas_loss_names_and_f32_probabilities_retain_source_arithmetic() {
    let p = peptide(".(Acetyl)SG.(Amidated)");
    let mut g = generator(0.005);
    g.add_first_prefix_ion = true;
    g.add_losses = true;
    g.add_precursor_peaks = true;
    g.sort_by_position = false;
    g.intensities.b = 0.7;
    g.relative_loss_intensity = 0.9;
    let s = g.generate(&p, 1, 1, None).unwrap();
    // Explicit compositions derive from Ser=C3H7NO3, Gly=C2H5NO2,
    // peptide condensation, Acetyl=C2H2O, Amidated=HNO-1, and neutral H.
    for (name, formula, intensity) in [
        ("b1", [5, 8, 1, 3], f64::from(0.7_f32)),
        ("y1", [2, 7, 2, 1], 1.0),
        (
            "b1-H2O1+",
            [5, 6, 1, 2],
            f64::from(0.7_f32) * f64::from(0.9_f32),
        ),
        ("[M+H]+", [7, 14, 3, 4], 1.0),
        ("[M+H-H2O]+", [7, 12, 3, 3], 1.0),
        ("[M+H-NH3]+", [7, 11, 2, 4], 1.0),
    ] {
        assert_group(&s, name, 1, formula, 0.005, intensity);
    }
    let p0 = chno_envelope([5, 6, 1, 2], 0.005)[0].1;
    assert_ne!(
        (f64::from(0.7_f32) * f64::from(0.9_f32) * f64::from(p0)) as f32,
        (f64::from(0.7_f32 * 0.9_f32) * f64::from(p0)) as f32,
        "chosen source loss distinguishes an extra intermediate f32 rounding"
    );
    assert!(!names(&s).iter().any(|name| name.starts_with("y1-")));
    assert!(names(&s).iter().all(|name| !name.contains("Acetyl")));
}

#[test]
fn internal_and_immonium_peaks_do_not_acquire_fine_envelopes() {
    let mut g = generator(0.0001);
    g.ion_series.clear();
    g.add_internal_fragments = true;
    g.add_abundant_immonium_ions = true;
    let p = peptide(".(Acetyl)APGGA.(Amidated)");
    let fine = g.generate(&p, 1, 2, None).unwrap();
    g.isotope_model = Isotopes::None;
    let plain = g.generate(&p, 1, 2, None).unwrap();
    assert_eq!(fine, plain);
    assert_eq!(fine.len(), 9); // four internal fragments per charge, one P ion.
    assert!(names(&fine).iter().any(|name| name == "iP+"));
    // A formula-free, mass-tagged internal residue stays usable when the only
    // enabled peak families are the source's independent internal/immonium ones.
    g.isotope_model = Isotopes::Fine {
        unexplained_probability: 0.0001,
    };
    assert!(
        g.generate(&peptide("AGX[111.23456789]GA"), 1, 1, None)
            .is_ok()
    );
}

#[test]
fn late_precursor_envelope_limit_keeps_complete_append_target_unchanged() {
    let delta: EmpiricalFormula = "F1000000".parse().unwrap();
    let record = ResidueModification::from_record(ModificationRecord {
        name: "LargeFormula".into(),
        full_name: "Fine isotope atom-limit boundary".into(),
        origin: Some('A'),
        diff_mono_mass: delta.mono_mass(),
        diff_average_mass: delta.average_mass(),
        diff_formula: delta,
        ..Default::default()
    })
    .unwrap();
    let registry = ModificationsDB::from_records(vec![record]).unwrap();
    let p = AASequence::parse_with_registry("GGA(LargeFormula)", &registry).unwrap();
    let mut g = generator(0.05);
    g.ion_series = vec![Ion::B];
    g.add_first_prefix_ion = true;
    assert!(!g.generate(&p, 1, 1, None).unwrap().is_empty());
    // The b1/b2 envelopes fit first; only the later complete precursor exceeds
    // the fine atom cap. Existing annotations, metadata and precursor survive.
    g.add_precursor_peaks = true;
    let mut target = MSSpectrum {
        peaks: vec![Peak1D::new(42.0, 3.0)],
        name: "original".into(),
        precursors: vec![Precursor::new(123.0, 2)],
        string_data_arrays: vec![DataArray::new("IonNames", vec!["observed".into()])],
        integer_data_arrays: vec![DataArray::new("Charges", vec![0])],
        ..Default::default()
    };
    target.metadata.insert("sample".into(), "unchanged".into());
    let before = target.clone();
    let error = g.append_to(&mut target, &p, 1, 1, None).unwrap_err();
    assert!(error.to_string().contains("fine isotope atom limit"));
    assert_eq!(target, before);
}

#[test]
fn fine_work_is_cumulative_across_envelopes_before_any_append_mutation() {
    let delta: EmpiricalFormula = "C990000".parse().unwrap();
    let record = ResidueModification::from_record(ModificationRecord {
        name: "LargeCarbon".into(),
        full_name: "Cumulative fine isotope setup boundary".into(),
        origin: Some('G'),
        diff_mono_mass: delta.mono_mass(),
        diff_average_mass: delta.average_mass(),
        diff_formula: delta,
        ..Default::default()
    })
    .unwrap();
    let registry = ModificationsDB::from_records(vec![record]).unwrap();
    let p =
        AASequence::parse_with_registry(&format!("G(LargeCarbon){}", "G".repeat(105)), &registry)
            .unwrap();
    let before_peptide = p.clone();
    let mut g = generator(0.999999);
    g.ion_series = vec![Ion::B];
    g.add_first_prefix_ion = true;
    // Each envelope needs only its most probable configuration and stays below
    // one million atoms. The factorial setup alone costs >990,000 work units
    // per prefix; 105 prefixes exceed the shared 100-million allowance while
    // requesting merely 105 output peaks. Resetting the budget per envelope
    // incorrectly succeeds. A single prefix demonstrates individually valid work.
    let one = g.generate(&p.prefix(2).unwrap(), 1, 1, None).unwrap();
    assert_eq!(one.len(), 1);
    let mut target = MSSpectrum {
        peaks: vec![Peak1D::new(42.0, 3.0)],
        name: "original".into(),
        precursors: vec![Precursor::new(123.0, 2)],
        string_data_arrays: vec![DataArray::new("IonNames", vec!["observed".into()])],
        integer_data_arrays: vec![DataArray::new("Charges", vec![0])],
        ..Default::default()
    };
    target.metadata.insert("sample".into(), "unchanged".into());
    let before = target.clone();
    let error = g.append_to(&mut target, &p, 1, 1, None).unwrap_err();
    assert!(error.to_string().contains("fine isotope work limit"));
    assert_eq!(target, before);
    assert_eq!(p, before_peptide);
}

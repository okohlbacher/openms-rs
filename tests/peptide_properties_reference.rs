// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Independent source tables, literal assertions and mathematical references.
//! See data/peptide_properties_provenance.json. No C++ execution.

use openms::chemistry::{
    AAIndex, AAIndexScale, AASequence, HydrophobicityProfile as Hydro, HydrophobicityScale,
    IsoelectricPoint, ProteomicsPkaScale,
};
use std::collections::BTreeMap;

const AA: &str = include_str!("data/peptide_properties_aaindex_gb.tsv");
const HYDRO: &str = include_str!("data/peptide_properties_hydrophobicity.tsv");
const PKA: &str = include_str!("data/peptide_properties_pka.tsv");
const SCALARS: &str = include_str!("data/peptide_properties_source_scalar_literals.tsv");
const GB_DECIMAL: &str = include_str!("data/peptide_properties_gb_decimal.tsv");
const CANONICAL: &str = "ACDEFGHIKLMNPQRSTVWY";

fn rows(text: &str) -> impl Iterator<Item = Vec<&str>> {
    text.lines().skip(1).map(|line| line.split('\t').collect())
}
fn sequence(text: &str) -> AASequence {
    text.parse().unwrap()
}
fn bits(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).unwrap())
}
fn literal(decimal: &str, encoded: &str) -> f64 {
    let value = bits(encoded);
    assert_eq!(decimal.parse::<f64>().unwrap().to_bits(), value.to_bits());
    value
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; absolute tolerance {tolerance}"
    );
}
fn ulps(actual: f64, expected: f64, maximum: u64) {
    assert!(actual.is_finite() && expected.is_finite() && actual > 0.0 && expected > 0.0);
    assert!(
        actual.to_bits().abs_diff(expected.to_bits()) <= maximum,
        "{actual:.17e} != {expected:.17e}; maximum {maximum} ULPs"
    );
}
fn pka_scale(name: &str) -> ProteomicsPkaScale {
    match name {
        "lehninger" => ProteomicsPkaScale::Lehninger,
        "emboss" => ProteomicsPkaScale::Emboss,
        "sillero" => ProteomicsPkaScale::Sillero,
        "bjellqvist" => ProteomicsPkaScale::Bjellqvist,
        _ => panic!("unknown fixture scale {name}"),
    }
}

#[test]
fn all_public_aaindex_literals_and_unusual_indicator_memberships() {
    let mut count = 0;
    for row in rows(AA).filter(|r| r[0].starts_with("get")) {
        let scale = AAIndexScale::ALL
            .into_iter()
            .find(|scale| scale.accession() == &row[0][3..])
            .unwrap();
        let expected = literal(row[2], row[3]);
        assert_eq!(
            scale
                .value(row[1].chars().next().unwrap())
                .unwrap()
                .to_bits(),
            expected.to_bits()
        );
        count += 1;
    }
    assert_eq!(count, 200);
    assert_eq!(AAIndexScale::ALL.len(), 10);
    for residue in ('A'..='Z').chain(['a', '?', 'é']) {
        for (actual, members) in [
            (AAIndex::aliphatic(residue), "AGFIMLPV"),
            (AAIndex::acidic(residue), "DE"),
            (AAIndex::basic(residue), "KRHW"),
            (AAIndex::polar(residue), "STYHCNQW"),
        ] {
            assert_eq!(actual, if members.contains(residue) { 1.0 } else { 0.0 });
        }
        for scale in AAIndexScale::ALL {
            assert_eq!(scale.value(residue).is_ok(), CANONICAL.contains(residue));
        }
    }
}

#[test]
fn all_hydrophobicity_cells_include_explicit_sentinel_rejections() {
    let mut counts = [0, 0];
    assert_eq!(HydrophobicityScale::ALL.len(), 7);
    for row in rows(HYDRO) {
        let scale = HydrophobicityScale::ALL[row[0].parse::<usize>().unwrap()];
        let residue = row[2].chars().next().unwrap();
        let expected = literal(row[3], row[4]);
        if row[5] == "true" {
            assert_eq!(scale.value(residue).unwrap().to_bits(), expected.to_bits());
            counts[0] += 1;
        } else {
            assert_eq!(expected, 999.0);
            assert!(scale.value(residue).is_err());
            counts[1] += 1;
        }
    }
    assert_eq!(counts, [140, 42]);
    for scale in HydrophobicityScale::ALL {
        for residue in ['a', '?', 'é'] {
            assert!(scale.value(residue).is_err());
        }
    }
}

fn gb_tables() -> BTreeMap<(&'static str, char), f64> {
    let table: BTreeMap<_, _> = rows(AA)
        .filter(|row| !row[0].starts_with("get"))
        .map(|row| {
            (
                (row[0], row[1].chars().next().unwrap()),
                literal(row[2], row[3]),
            )
        })
        .collect();
    assert_eq!(table.len(), 62);
    table
}
fn direct_gb(text: &str, temperature: f64, table: &BTreeMap<(&str, char), f64>) -> f64 {
    // Original source constants and split pairing, independently read from fixtures.
    let rt = ((6.022_136_7e23 * 1.380_657e-23) / 1000.0) * temperature;
    let residues: Vec<_> = text.chars().collect();
    let mut total = 0.0;
    for split in 0..=residues.len() {
        let left = if split == 0 { '>' } else { residues[split - 1] };
        let right = residues.get(split).copied().unwrap_or('<');
        let mut contribution =
            ((table[&("GBleft_", left)] + table[&("GBdeltaright_", right)]) / rt).exp();
        if split > 0 && split < residues.len() {
            contribution += (table[&("GBsidechain_", right)] / rt).exp();
        }
        total += contribution;
    }
    rt * total.ln() / 2.0_f64.ln()
}

#[test]
fn all_private_gb_entries_are_exercised_by_independent_site_expressions() {
    let table = gb_tables();
    // At 100000 K, small-energy terms remain observable beside the strongest site.
    for first in CANONICAL.chars() {
        let text = first.to_string();
        ulps(
            AAIndex::calculate_gb(&sequence(&text), 100_000.0).unwrap(),
            direct_gb(&text, 100_000.0, &table),
            4,
        );
        for second in CANONICAL.chars() {
            let text = format!("{first}{second}");
            ulps(
                AAIndex::calculate_gb(&sequence(&text), 100_000.0).unwrap(),
                direct_gb(&text, 100_000.0, &table),
                4,
            );
        }
    }
    for text in ["AK", "KA", "ARR", "ALEGDEK", "GTVVTGR", "EHVLLAR"] {
        ulps(
            AAIndex::calculate_gb(&sequence(text), 500.0).unwrap(),
            direct_gb(text, 500.0, &table),
            4,
        );
    }
}

#[test]
fn gas_basicity_matches_independent_ninety_digit_arithmetic() {
    let mut count = 0;
    for row in rows(GB_DECIMAL) {
        let temperature = literal(row[1], row[2]);
        let expected = literal(row[3], row[4]);
        ulps(
            AAIndex::calculate_gb(&sequence(row[0]), temperature).unwrap(),
            expected,
            4,
        );
        count += 1;
    }
    assert_eq!(count, 100);
}

fn concentration_charge(pka: f64, ph: f64, acidic: bool) -> f64 {
    // Concentration-ratio form, algebraically independent of the implementation's
    // 1/(1+10^(difference)) expression. Test pH values keep both powers finite.
    let proton = 10.0_f64.powf(-ph);
    let constant = 10.0_f64.powf(-pka);
    if acidic {
        -constant / (constant + proton)
    } else {
        proton / (constant + proton)
    }
}
fn base_sidechain(scale: &str, residue: &str, ph: f64) -> f64 {
    rows(PKA)
        .find(|row| row[0] == scale && row[1] == "base" && row[2] == residue)
        .map_or(0.0, |row| {
            concentration_charge(literal(row[4], row[5]), ph, "DECY".contains(residue))
        })
}

#[test]
fn every_pka_constant_has_an_independently_isolated_charge_check() {
    let mut count = 0;
    for row in rows(PKA) {
        let pka = literal(row[4], row[5]);
        let scales: Vec<_> = if row[0] == "all" {
            vec!["lehninger", "emboss", "sillero", "bjellqvist"]
        } else {
            vec![row[0]]
        };
        for name in scales {
            let calc = IsoelectricPoint {
                scale: pka_scale(name),
                ..Default::default()
            };
            for ph in [pka - 0.75, pka, pka + 0.75] {
                let terminal = row[2] == "nterm" || row[2] == "cterm";
                let code = if terminal {
                    if row[3].is_empty() { "X" } else { row[3] }
                } else {
                    row[2]
                };
                let mut peptide = sequence(code);
                if row[2] != "nterm" {
                    peptide.set_n_terminal_mass_tag("+0.123456789").unwrap();
                }
                if row[2] != "cterm" {
                    peptide.set_c_terminal_mass_tag("+0.123456789").unwrap();
                }
                let acidic = row[2] == "cterm" || "DECYU".contains(row[2]);
                let expected = concentration_charge(pka, ph, acidic)
                    + if terminal {
                        base_sidechain(name, code, ph)
                    } else {
                        0.0
                    };
                near(calc.compute_charge(&peptide, ph).unwrap(), expected, 2e-14);
            }
        }
        count += 1;
    }
    assert_eq!(count, 59);
}

fn source_pi(text: &str, scale: ProteomicsPkaScale, tolerance: f64) -> f64 {
    IsoelectricPoint { scale, tolerance }
        .compute_pi(&sequence(text))
        .unwrap()
}

#[test]
fn all_thirty_upstream_scalar_assertions_remain_separate_from_derived_goldens() {
    let kd = HydrophobicityScale::KyteDoolittle;
    let le = ProteomicsPkaScale::Lehninger;
    let bj = ProteomicsPkaScale::Bjellqvist;
    let em = ProteomicsPkaScale::Emboss;
    let mut count = 0;
    for row in rows(SCALARS) {
        let expected = literal(row[2], row[3]);
        let line: u32 = row[5].parse().unwrap();
        let actual = match (row[0], line) {
            ("AAIndex", 37) => AAIndex::calculate_gb(&sequence("ALEGDEK"), 500.0).unwrap(),
            ("AAIndex", 38) => AAIndex::calculate_gb(&sequence("GTVVTGR"), 500.0).unwrap(),
            ("AAIndex", 39) => AAIndex::calculate_gb(&sequence("EHVLLAR"), 500.0).unwrap(),
            ("IsoelectricPoint", 77) => source_pi("A", le, 1e-4),
            ("IsoelectricPoint", 83) => source_pi("K", le, 1e-4),
            ("IsoelectricPoint", 89) => source_pi("D", le, 1e-4),
            ("IsoelectricPoint", 97 | 178) => {
                let (text, scale) = if line == 97 {
                    ("PEPTIDE", le)
                } else {
                    ("DEKRPEPTIDE", ProteomicsPkaScale::Sillero)
                };
                let calc = IsoelectricPoint {
                    scale,
                    ..Default::default()
                };
                let peptide = sequence(text);
                calc.compute_charge(&peptide, calc.compute_pi(&peptide).unwrap())
                    .unwrap()
            }
            ("IsoelectricPoint", 114) => source_pi("U", le, 1e-4),
            ("IsoelectricPoint", 118) => source_pi(&"R".repeat(100), le, 1e-4),
            ("IsoelectricPoint", 125) => source_pi("A", bj, 1e-4),
            ("IsoelectricPoint", 132) => source_pi("P", bj, 1e-4),
            ("IsoelectricPoint", 149) => source_pi("ACDEFGHIK", em, 1e-6),
            ("IsoelectricPoint", 150) => source_pi("PEPTIDER", em, 1e-6),
            ("IsoelectricPoint", 151) => source_pi("SAMPLER", em, 1e-6),
            ("IsoelectricPoint", 152) => source_pi("YGGFMR", em, 1e-6),
            ("IsoelectricPoint", 153) => source_pi("DEEEAANK", em, 1e-6),
            ("HydrophobicityProfile", 48) => Hydro::compute_gravy(&sequence("ACDE")).unwrap(),
            ("HydrophobicityProfile", 64..=67) => {
                Hydro::compute_profile(&sequence("ACDE"), kd).unwrap()[(line - 64) as usize]
            }
            ("HydrophobicityProfile", 68) => {
                Hydro::compute_profile(&sequence("ACDE"), HydrophobicityScale::Eisenberg).unwrap()
                    [0]
            }
            ("HydrophobicityProfile", 79..=81) => {
                Hydro::compute_windowed_profile(&sequence("ACDEF"), 3, kd).unwrap()
                    [(line - 79) as usize]
            }
            ("HydrophobicityProfile", 82) => {
                Hydro::compute_windowed_profile(&sequence("ACDEF"), 6, kd).unwrap()[0]
            }
            ("HydrophobicityProfile", 94..=96) => {
                Hydro::compute_hydrophobic_moment(&sequence("ACDEF"), 3, 100.0).unwrap()
                    [(line - 94) as usize]
            }
            other => panic!("unmapped source assertion {other:?}"),
        };
        // Missing ClassTest relative defaults are not inferred. These tolerances
        // explicitly cover rounded source literals and default pI interval width.
        let tolerance = match row[0] {
            "AAIndex" => 0.01,
            "IsoelectricPoint" if expected == 0.0 => 1e-3,
            "IsoelectricPoint" => 5e-4,
            _ => 1e-8,
        };
        near(actual, expected, tolerance);
        count += 1;
    }
    assert_eq!(count, 30);
}

#[test]
fn independent_normalized_eisenberg_moments_and_mass_unavailable_annotations() {
    let eisenberg: BTreeMap<_, _> = rows(HYDRO)
        .filter(|row| row[0] == "1" && row[5] == "true")
        .map(|row| (row[2].chars().next().unwrap(), literal(row[3], row[4])))
        .collect();
    let text = "FLIGKWVRP";
    let values: Vec<_> = text.chars().map(|code| eisenberg[&code]).collect();
    for (window, degrees) in [(3, 100.0_f64), (5, 160.0), (4, 0.0), (3, -100.0)] {
        let expected: Vec<_> = values
            .windows(window)
            .map(|slice| {
                let phase = degrees * std::f64::consts::PI / 180.0;
                let (sine, cosine) = slice.iter().enumerate().fold((0.0, 0.0), |(s, c), (i, h)| {
                    (
                        s + h * (phase * i as f64).sin(),
                        c + h * (phase * i as f64).cos(),
                    )
                });
                (sine * sine + cosine * cosine).sqrt() / window as f64
            })
            .collect();
        let actual = Hydro::compute_hydrophobic_moment(&sequence(text), window, degrees).unwrap();
        assert_eq!(actual.len(), expected.len());
        for (actual, expected) in actual.into_iter().zip(expected) {
            near(actual, expected, 2e-14);
        }
    }
    let mut annotated = sequence(text);
    annotated.set_mass_tag(4, "+0.123456789").unwrap();
    annotated.set_n_terminal_mass_tag("+0.123456789").unwrap();
    annotated.set_c_terminal_mass_tag("+0.123456789").unwrap();
    assert!(annotated.formula().is_err());
    for scale in HydrophobicityScale::ALL {
        assert_eq!(
            Hydro::compute_profile(&annotated, scale).unwrap(),
            Hydro::compute_profile(&sequence(text), scale).unwrap()
        );
    }
    assert_eq!(
        AAIndex::calculate_gb(&annotated, 500.0).unwrap(),
        AAIndex::calculate_gb(&sequence(text), 500.0).unwrap()
    );
    let mut unknown = sequence("BZXJ");
    unknown.set_n_terminal_mass_tag("+0.123456789").unwrap();
    unknown.set_c_terminal_mass_tag("+0.123456789").unwrap();
    assert!(unknown.mono_mass().is_err());
    assert_eq!(
        IsoelectricPoint::default().compute_pi(&unknown).unwrap(),
        0.0
    );
    assert_eq!(
        IsoelectricPoint::default()
            .compute_charge(&unknown, 7.0)
            .unwrap(),
        0.0
    );
}

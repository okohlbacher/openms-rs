// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{
    AAIndex, AASequence, HydrophobicityProfile as Hydro, HydrophobicityScale, IsoelectricPoint,
    ModifiedPeptideGenerator, ProteaseDigestion,
};

fn seq(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-12,
        "{actual:.16} != {expected:.16}"
    );
}

#[test]
fn digestion_retains_only_original_terminal_blocks_and_properties_follow_parent_codes() {
    let protein = seq("(Acetyl)AC(Carbamidomethyl)M(Oxidation)KAGHIK(Amidated)");
    let original = protein.clone();
    let products = ProteaseDigestion::default().digest(&protein).unwrap();
    assert_eq!(products.len(), 2);
    assert_eq!((products[0].start, products[0].end), (0, 4));
    assert_eq!((products[1].start, products[1].end), (4, 9));
    let first = &products[0].sequence;
    let second = &products[1].sequence;
    assert!(first.n_terminal_modification().is_some());
    assert!(first.c_terminal_modification().is_none());
    assert!(second.n_terminal_modification().is_none());
    assert!(second.c_terminal_modification().is_some());

    // Independent Henderson-Hasselbalch terms for the surviving groups at pH7:
    // first: C terminus, Cys, Lys; second: new N terminus, His, Lys.
    let acidic = |pka: f64| -1.0 / (1.0 + 10.0_f64.powf(pka - 7.0));
    let basic = |pka: f64| 1.0 / (1.0 + 10.0_f64.powf(7.0 - pka));
    let calculator = IsoelectricPoint::default();
    near(
        calculator.compute_charge(first, 7.0).unwrap(),
        acidic(2.34) + acidic(8.18) + basic(10.53),
    );
    near(
        calculator.compute_charge(second, 7.0).unwrap(),
        basic(9.69) + basic(6.0) + basic(10.53),
    );
    for (peptide, expected) in [
        (first, vec![1.8, 2.5, 1.9, -3.9]),
        (second, vec![1.8, -0.4, -3.2, 4.5, -3.9]),
    ] {
        assert_eq!(
            Hydro::compute_profile(peptide, HydrophobicityScale::KyteDoolittle).unwrap(),
            expected
        );
        near(
            Hydro::compute_gravy(peptide).unwrap(),
            expected.iter().sum::<f64>() / expected.len() as f64,
        );
        let parent = seq(peptide.as_str());
        assert_eq!(
            AAIndex::calculate_gb(peptide, 500.0).unwrap(),
            AAIndex::calculate_gb(&parent, 500.0).unwrap()
        );
    }
    assert_eq!(protein, original);
}

#[test]
fn formula_free_annotations_do_not_require_invented_compositions_for_properties() {
    let peptide = seq("AC[+0.123456789]MK");
    let parent = seq("ACMK");
    let before = peptide.clone();
    assert!(peptide.formula().is_err());
    let calculator = IsoelectricPoint::default();
    assert_eq!(
        calculator.compute_pi(&peptide).unwrap(),
        calculator.compute_pi(&parent).unwrap()
    );
    for scale in HydrophobicityScale::ALL {
        assert_eq!(
            Hydro::compute_profile(&peptide, scale).unwrap(),
            Hydro::compute_profile(&parent, scale).unwrap()
        );
        assert_eq!(
            Hydro::compute_windowed_profile(&peptide, 3, scale).unwrap(),
            Hydro::compute_windowed_profile(&parent, 3, scale).unwrap()
        );
    }
    assert_eq!(
        Hydro::compute_hydrophobic_moment(&peptide, 3, 100.0).unwrap(),
        Hydro::compute_hydrophobic_moment(&parent, 3, 100.0).unwrap()
    );
    for temperature in [100.0, 500.0] {
        assert_eq!(
            AAIndex::calculate_gb(&peptide, temperature).unwrap(),
            AAIndex::calculate_gb(&parent, temperature).unwrap()
        );
    }
    assert_eq!(peptide, before);
    let unknown = seq("X");
    assert!(unknown.mono_mass().is_err());
    assert_eq!(
        calculator.compute_pi(&unknown).unwrap(),
        calculator.compute_pi(&seq("A")).unwrap()
    );
    assert!(Hydro::compute_gravy(&unknown).is_err());
    assert!(AAIndex::calculate_gb(&unknown, 500.0).is_err());
}

#[test]
fn source_modification_placement_distinguishes_residue_annotations_from_terminal_blocks() {
    let modifications =
        ModifiedPeptideGenerator::get_modifications(&["Gln->pyro-Glu (N-term Q)"]).unwrap();
    let generator = ModifiedPeptideGenerator::default();
    let parent = seq("Q");
    let residue_variants = generator
        .variable_modifications(&modifications, &parent, 1, false)
        .unwrap();
    let terminal_variants = generator
        .variable_modifications(&modifications, &parent, 2, false)
        .unwrap();
    let residue = residue_variants
        .iter()
        .find(|p| {
            p.n_terminal_modification().is_none() && p.residue_modification(0).unwrap().is_some()
        })
        .unwrap();
    let terminal = terminal_variants
        .iter()
        .find(|p| p.n_terminal_modification().is_some())
        .unwrap();
    let calculator = IsoelectricPoint::default();
    assert_eq!(
        calculator.compute_pi(residue).unwrap(),
        calculator.compute_pi(&parent).unwrap()
    );
    assert_eq!(calculator.compute_pi(terminal).unwrap(), 0.0);
    near(
        calculator.compute_charge(terminal, 7.0).unwrap(),
        -1.0 / (1.0 + 10.0_f64.powf(2.34 - 7.0)),
    );
    assert_eq!(
        Hydro::compute_gravy(residue).unwrap(),
        Hydro::compute_gravy(terminal).unwrap()
    );
    assert_eq!(
        AAIndex::calculate_gb(residue, 500.0).unwrap(),
        AAIndex::calculate_gb(terminal, 500.0).unwrap()
    );
}

#[cfg(feature = "idxml")]
#[test]
fn derived_properties_and_annotated_sequences_roundtrip_as_typed_identification_metadata() {
    use openms::format::idxml::{self, IdXmlDocument};
    use openms::identification::{PeptideHit, PeptideIdentification, ProteinIdentification};
    use openms::metadata::MetaValue;
    let peptide = seq("(Acetyl)AC(Carbamidomethyl)M(Oxidation)K");
    let calculator = IsoelectricPoint::default();
    let mut hit = PeptideHit::new(1.0, 1, 2, peptide.clone()).unwrap();
    for (name, value) in [
        (
            "model:Lehninger:pI",
            calculator.compute_pi(&peptide).unwrap(),
        ),
        (
            "model:Lehninger:charge_pH7",
            calculator.compute_charge(&peptide, 7.0).unwrap(),
        ),
        (
            "model:KyteDoolittle:GRAVY",
            Hydro::compute_gravy(&peptide).unwrap(),
        ),
        (
            "model:GB:100K",
            AAIndex::calculate_gb(&peptide, 100.0).unwrap(),
        ),
    ] {
        hit.metadata
            .insert(name.into(), MetaValue::try_from(value).unwrap());
    }
    let document = IdXmlDocument {
        protein_identifications: vec![ProteinIdentification {
            identifier: "properties".into(),
            date_time: Some("2026-09-10T12:00:00".into()),
            search_engine: "property demonstration".into(),
            ..Default::default()
        }],
        peptide_identifications: vec![PeptideIdentification {
            identifier: "properties".into(),
            hits: vec![hit],
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut output = Vec::new();
    idxml::write(&mut output, &document).unwrap();
    assert_eq!(idxml::read(output.as_slice()).unwrap(), document);
}

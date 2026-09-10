// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::MSSpectrum;
use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationDefinitionsSet, ModificationMatchOptions,
    ModifiedPeptideGenerator, ProteaseDigestion, TheoreticalSpectrumGenerator,
};
use std::collections::BTreeMap;

fn variants() -> Vec<AASequence> {
    let products = ProteaseDigestion::default()
        .digest(&AASequence::parse("ACMMKAGHIK").unwrap())
        .unwrap();
    assert_eq!((products[0].start, products[0].end), (0, 5));
    assert_eq!(products[1].sequence.as_str(), "AGHIK");
    let generator = ModifiedPeptideGenerator::default();
    let fixed = ModifiedPeptideGenerator::get_modifications(&["Carbamidomethyl (C)"]).unwrap();
    let variable = ModifiedPeptideGenerator::get_modifications(&["Oxidation (M)"]).unwrap();
    let mut peptide = products[0].sequence.clone();
    generator
        .apply_fixed_modifications(&fixed, &mut peptide)
        .unwrap();
    generator
        .variable_modifications(&variable, &peptide, 2, true)
        .unwrap()
}
fn named(spectrum: &MSSpectrum) -> BTreeMap<String, f64> {
    let names = &spectrum
        .string_data_arrays
        .iter()
        .find(|a| a.name == "IonNames")
        .unwrap()
        .data;
    names
        .iter()
        .cloned()
        .zip(spectrum.peaks.iter().map(|p| p.mz))
        .collect()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-9,
        "{actual:.12} != {expected:.12}"
    );
}

#[test]
fn digested_variants_retain_fixed_chemistry_and_independent_fragment_shifts() {
    let peptides = variants();
    let expected = [
        "AC(Carbamidomethyl)MMK",
        "AC(Carbamidomethyl)MM(Oxidation)K",
        "AC(Carbamidomethyl)M(Oxidation)MK",
        "AC(Carbamidomethyl)M(Oxidation)M(Oxidation)K",
    ];
    assert_eq!(peptides.len(), expected.len());
    let generator = TheoreticalSpectrumGenerator {
        add_metainfo: true,
        add_first_prefix_ion: true,
        ..Default::default()
    };
    let original = AASequence::parse("ACMMK").unwrap();
    let fixed_formula = original
        .formula()
        .unwrap()
        .checked_add(&EmpiricalFormula::parse("C2H3N1O1").unwrap())
        .unwrap();
    let oxygen = EmpiricalFormula::parse("O1").unwrap();
    let base_peaks = named(&generator.generate(&peptides[0], 1, 2, Some(2)).unwrap());
    for (i, peptide) in peptides.iter().enumerate() {
        assert_eq!(peptide, &AASequence::parse(expected[i]).unwrap());
        let sites: &[usize] = match i {
            0 => &[],
            1 => &[3],
            2 => &[2],
            _ => &[2, 3],
        };
        let formula = fixed_formula
            .checked_add(&oxygen.checked_scale(sites.len() as i32).unwrap())
            .unwrap();
        assert_eq!(peptide.formula().unwrap(), formula);
        near(peptide.mono_mass().unwrap(), formula.mono_mass());
        near(
            peptide.mz(2).unwrap(),
            (formula.mono_mass() + 2.0 * openms::chemistry::PROTON_MASS_U) / 2.0,
        );
        let spectrum = generator.generate(peptide, 1, 2, Some(2)).unwrap();
        spectrum.validate().unwrap();
        let peaks = named(&spectrum);
        assert_eq!(peaks.len(), base_peaks.len());
        for charge in 1..=2 {
            for ordinal in 1..peptide.len() {
                for series in ['b', 'y'] {
                    let included = sites
                        .iter()
                        .filter(|&&site| {
                            if series == 'b' {
                                site < ordinal
                            } else {
                                site >= peptide.len() - ordinal
                            }
                        })
                        .count();
                    let name = format!("{series}{ordinal}{}", "+".repeat(charge));
                    near(
                        peaks[&name],
                        base_peaks[&name] + included as f64 * oxygen.mono_mass() / charge as f64,
                    );
                }
            }
        }
    }
}

fn identifications(
    peptides: Vec<AASequence>,
) -> Vec<openms::identification::PeptideIdentification> {
    use openms::identification::{PeptideHit, PeptideIdentification};
    peptides
        .into_iter()
        .enumerate()
        .map(|(i, sequence)| PeptideIdentification {
            identifier: "modification-workflow".into(),
            rt: Some(i as f64),
            mz: Some(sequence.mz(2).unwrap()),
            hits: vec![PeptideHit::new(0.0, 0, 2, sequence).unwrap()],
            ..Default::default()
        })
        .collect()
}

#[test]
fn generated_evidence_recovers_fixed_and_variable_search_definitions() {
    let ids = identifications(variants());
    let definitions =
        ModificationDefinitionsSet::from_names(&["Carbamidomethyl (C)"], &["Oxidation (M)"])
            .unwrap();
    let mut inferred = ModificationDefinitionsSet::default();
    inferred.infer_from_peptides(&ids).unwrap();
    assert_eq!(inferred.fixed_names(), definitions.fixed_names());
    assert_eq!(inferred.variable_names(), definitions.variable_names());
    for id in &ids {
        assert!(inferred.is_compatible(&id.hits[0].sequence).unwrap());
    }
    let matches = inferred
        .find_matches(
            15.994915,
            &ModificationMatchOptions {
                residue: "M".into(),
                consider_fixed: false,
                tolerance: 0.000001,
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].definition.modification_name(), "Oxidation (M)");
    assert!(!matches[0].definition.fixed);
    assert!(matches[0].mass_error <= 0.000001);
}

#[cfg(feature = "idxml")]
fn document(peptides: Vec<AASequence>) -> openms::format::idxml::IdXmlDocument {
    use openms::identification::ProteinIdentification;
    openms::format::idxml::IdXmlDocument {
        protein_identifications: vec![ProteinIdentification {
            identifier: "modification-workflow".into(),
            date_time: Some("2026-09-10T12:00:00".into()),
            search_engine: "variant demonstration".into(),
            ..Default::default()
        }],
        peptide_identifications: identifications(peptides),
        ..Default::default()
    }
}

#[cfg(feature = "idxml")]
#[test]
fn generated_variants_roundtrip_through_identification_xml() {
    use openms::format::idxml;
    let mut doc = document(variants());
    let search = &mut doc.protein_identifications[0].search_parameters;
    let definitions =
        ModificationDefinitionsSet::from_names(&["Carbamidomethyl (C)"], &["Oxidation (M)"])
            .unwrap();
    search.fixed_modifications = definitions.fixed_names().into_iter().collect();
    search.variable_modifications = definitions.variable_names().into_iter().collect();
    let mut encoded = Vec::new();
    idxml::write(&mut encoded, &doc).unwrap();
    assert_eq!(idxml::read(std::io::Cursor::new(encoded)).unwrap(), doc);
}

#[cfg(feature = "idxml")]
#[test]
fn source_terminal_placements_fail_before_idxml_can_relocate_or_lose_them() {
    use openms::format::idxml;
    use std::io::Write;
    #[derive(Default)]
    struct ObservedWriter {
        writes: usize,
        flushes: usize,
    }
    impl Write for ObservedWriter {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.writes += 1;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            self.flushes += 1;
            Ok(())
        }
    }
    let generator = ModifiedPeptideGenerator::default();
    let pyro = ModifiedPeptideGenerator::get_modifications(&["Gln->pyro-Glu (N-term Q)"]).unwrap();
    // The max-one source path stores this terminal record on the residue.
    // The general path can place it on a terminus of a nonmatching residue.
    for (peptide, maximum) in [("Q", 1), ("A", 2)] {
        let generated = generator
            .variable_modifications(&pyro, &AASequence::parse(peptide).unwrap(), maximum, false)
            .unwrap();
        assert!(!generated.is_empty());
        let doc = document(generated);
        let before = doc.clone();
        let mut writer = ObservedWriter::default();
        let error = idxml::write(&mut writer, &doc).unwrap_err();
        assert!(
            error.to_string().contains("modification placement"),
            "{error}"
        );
        assert_eq!((writer.writes, writer.flushes), (0, 0));
        assert_eq!(doc, before);
    }
}

#[test]
fn custom_mass_only_variants_propagate_to_fragments_without_inventing_composition() {
    use openms::Peak1D;
    use openms::chemistry::{ModificationsDB, TheoreticalIsotopeModel};
    let database = ModificationsDB::from_tsv(
        "1000001\tLabMass\tLaboratory mass-only modification\tM\tanywhere\t12.5\t12.6\t\t0\tOther\t\n"
    ).unwrap();
    let original = AASequence::parse("AMK").unwrap();
    let peptides = ModifiedPeptideGenerator::default()
        .variable_modifications(&[database.entries()[0].clone()], &original, 1, false)
        .unwrap();
    drop(database);
    assert_eq!(peptides.len(), 1);
    let modified = &peptides[0];
    near(
        modified.mono_mass().unwrap(),
        original.mono_mass().unwrap() + 12.5,
    );
    assert!(modified.formula().is_err());
    assert!(modified.average_mass().is_err());
    assert!(
        modified
            .residue_modification(1)
            .unwrap()
            .unwrap()
            .diff_formula()
            .is_err()
    );
    let generator = TheoreticalSpectrumGenerator {
        add_first_prefix_ion: true,
        add_metainfo: true,
        ..Default::default()
    };
    let base = named(&generator.generate(&original, 1, 2, Some(2)).unwrap());
    let actual = named(&generator.generate(modified, 1, 2, Some(2)).unwrap());
    for charge in 1..=2 {
        for series in ['b', 'y'] {
            let suffix = "+".repeat(charge);
            for ordinal in 1..=2 {
                let key = format!("{series}{ordinal}{suffix}");
                near(
                    actual[&key],
                    base[&key]
                        + if ordinal == 2 {
                            12.5 / charge as f64
                        } else {
                            0.0
                        },
                );
            }
        }
    }
    let mut output = MSSpectrum::from_peaks(vec![Peak1D::new(80.0, 1.0)]);
    let before = output.clone();
    let isotopes = TheoreticalSpectrumGenerator {
        isotope_model: TheoreticalIsotopeModel::Coarse { max_peaks: 3 },
        ..Default::default()
    };
    assert!(
        isotopes
            .append_to(&mut output, modified, 1, 2, Some(2))
            .is_err()
    );
    assert_eq!(output, before);
    #[cfg(feature = "idxml")]
    {
        // The private custom registry is not implicitly registered by idXML.
        let doc = document(peptides);
        let mut bytes = Vec::new();
        assert!(openms::format::idxml::write(&mut bytes, &doc).is_err());
        assert!(bytes.is_empty());
    }
}

// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Independent integration oracles traced from OpenMS4-core 7c029e8cdba6abab503708ecdd56f6ab55e38ce4.
// Source paths below are relative to that repository; SHA-256:
// src/openms/source/CHEMISTRY/TheoreticalSpectrumGenerator.cpp
// 763d7e37041d3bed8f0b9fb122e43b5d6890fc699753e0f877ee3c50f2698325
// src/tests/class_tests/openms/source/AASequence_test.cpp
// d40acc86f1b46211db7b66c613c0b0156ee1c0917021f6c0941c50c9494cb95e
// src/openms/source/CHEMISTRY/AASequence.cpp
// 8307ab416dedff8473f36c26e4259382eed7aaa7f82711c581c9c0832641b693
// src/openms/source/CHEMISTRY/Residue.cpp
// 88524929645ad4ff2e1a9ce0ece281942da37f0e70986e72b4b8e1a73cacad83
// share/OpenMS/CHEMISTRY/Enzymes.xml
// 2f160f3ec32db6257cb4eee43fd594b48cb36398a17b7af21bbba76bc8b16995
// src/openms/source/ANALYSIS/ID/PeptideIndexing.cpp
// ccf2f8ec617c2272ca779d619bb3a29d4de936284d0c8e4820ad37430034d4a5

use openms::analysis::peptide_indexing::{DecoyRule, MissingDecoyAction, PeptideIndexing};
use openms::chemistry::{
    AASequence, ProteaseDigestion, TheoreticalIsotopeModel, TheoreticalSpectrumGenerator,
};
use openms::format::fasta::FASTAEntry;
use openms::identification::{
    FlankingResidue, PeptideHit, PeptideIdentification, ProteinIdentification, TargetDecoyType,
};
use openms::kernel::{MSSpectrum, Peak1D};
use std::collections::BTreeMap;

fn sequence(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-9,
        "{actual:.12} != {expected:.12}"
    );
}
fn annotated(spectrum: &MSSpectrum) -> BTreeMap<String, f64> {
    let names = &spectrum
        .string_data_arrays
        .iter()
        .find(|array| array.name == "IonNames")
        .unwrap()
        .data;
    let values: BTreeMap<_, _> = names
        .iter()
        .cloned()
        .zip(spectrum.peaks.iter().map(|p| p.mz))
        .collect();
    assert_eq!(values.len(), spectrum.len(), "annotations must be unique");
    values
}

#[test]
fn ambiguous_sequences_digest_and_regain_chemistry_after_slicing() {
    // Source parsing does not calculate mass. The actual pinned Trypsin regex
    // is (?<=[KRX])(?!P), so X itself is also a cleavage residue.
    let protein = sequence("BZXKACDRXAR");
    assert!(protein.mono_mass().is_err());
    assert!(protein.average_mass().is_err());
    assert!(protein.formula().is_err());
    let products = ProteaseDigestion::default().digest(&protein).unwrap();
    assert_eq!(
        products
            .iter()
            .map(|p| (p.sequence.as_str(), p.start, p.end))
            .collect::<Vec<_>>(),
        [
            ("BZX", 0, 3),
            ("K", 3, 4),
            ("ACDR", 4, 8),
            ("X", 8, 9),
            ("AR", 9, 11)
        ]
    );
    assert!(products[0].sequence.mono_mass().is_err());
    assert_eq!(products[2].sequence, sequence("ACDR"));
    assert!(products[2].sequence.mono_mass().unwrap() > 0.0);
    assert!(products[3].sequence.mono_mass().is_err());
    assert_eq!(protein.subsequence(4..8).unwrap(), sequence("ACDR"));
}

#[test]
fn ambiguous_identifications_index_without_a_peptide_mass() {
    // Literal B/Z/X equality spends no ambiguity allowance. PeptideIndexing's
    // N-terminal methionine clipping retains transformed coordinates 1..=4.
    let mut runs = vec![ProteinIdentification {
        identifier: "ambiguous-run".into(),
        ..Default::default()
    }];
    let mut peptides = vec![PeptideIdentification {
        identifier: "ambiguous-run".into(),
        hits: vec![PeptideHit {
            sequence: sequence("BZXK"),
            score: 7.5,
            charge: 2,
            ..Default::default()
        }],
        ..Default::default()
    }];
    let fasta = vec![FASTAEntry {
        identifier: "target".into(),
        sequence: "MBZXK".into(),
        description: "literal ambiguous protein".into(),
    }];
    let report = PeptideIndexing {
        decoy_rule: DecoyRule::Prefix("DECOY_".into()),
        missing_decoy_action: MissingDecoyAction::Silent,
        max_ambiguities: 0,
        max_mismatches: 0,
        write_protein_sequence: true,
        ..Default::default()
    }
    .run(&fasta, &mut runs, &mut peptides)
    .unwrap();
    assert_eq!((report.target_hits, report.evidence_count), (1, 1));
    let hit = &peptides[0].hits[0];
    assert_eq!(hit.sequence.as_str(), "BZXK");
    assert!(hit.sequence.mono_mass().is_err());
    assert_eq!(hit.target_decoy_type().unwrap(), TargetDecoyType::Target);
    assert_eq!(hit.evidences[0].start, Some(1));
    assert_eq!(hit.evidences[0].end, Some(4));
    assert_eq!(hit.evidences[0].aa_before, FlankingResidue::Residue('M'));
    assert_eq!(hit.evidences[0].aa_after, FlankingResidue::CTerminus);
    assert_eq!(runs[0].hits[0].sequence, "MBZXK");
    peptides[0].validate().unwrap();
    runs[0].validate().unwrap();
}

#[test]
fn source_absolute_mass_tag_shifts_only_retaining_fragments_and_neighbor_losses() {
    // AASequence_test.cpp: PEPTX[999]IDE = PEPTIDE + 999 Da.
    // TSG.cpp's non-isotope paths add residue masses and subtract declared
    // losses. Insertion before I changes prefix/suffix ordinal independently.
    let generator = TheoreticalSpectrumGenerator {
        add_metainfo: true,
        add_first_prefix_ion: true,
        add_losses: true,
        ..Default::default()
    };
    let base = sequence("PEPTIDE");
    let tagged = sequence("PEPTX[999]IDE");
    near(
        tagged.mono_mass().unwrap(),
        base.mono_mass().unwrap() + 999.0,
    );
    assert!(tagged.formula().is_err());
    assert!(tagged.average_mass().is_err());
    let base_peaks = annotated(&generator.generate(&base, 1, 2, None).unwrap());
    let spectrum = generator.generate(&tagged, 1, 2, None).unwrap();
    let actual = annotated(&spectrum);
    let mut expected = BTreeMap::new();
    for charge in 1..=2 {
        let suffix = "+".repeat(charge);
        for (series, unchanged) in [("b", 4), ("y", 3)] {
            for ordinal in 1..=7 {
                let shifted = ordinal > unchanged;
                let original = ordinal - usize::from(shifted);
                for loss in ["", "-H2O1"] {
                    if let Some(&mass) =
                        base_peaks.get(&format!("{series}{original}{loss}{suffix}"))
                    {
                        expected.insert(
                            format!("{series}{ordinal}{loss}{suffix}"),
                            mass + if shifted { 999.0 / charge as f64 } else { 0.0 },
                        );
                    }
                }
            }
        }
    }
    assert_eq!(actual.len(), expected.len());
    for (name, mass) in expected {
        near(actual[&name], mass);
    }
    assert_eq!(spectrum.precursors[0].charge, 3);
    near(spectrum.precursors[0].mz, tagged.mz(3).unwrap());
    assert_eq!(tagged.suffix(3).unwrap(), sequence("IDE"));
    spectrum.validate().unwrap();
}

#[test]
fn anonymous_residue_modification_replaces_original_neutral_losses() {
    // Residue::setModification clears base losses. A mass-only tag has no
    // declared loss formula, so modified S must not keep unmodified S's water.
    let generator = TheoreticalSpectrumGenerator {
        add_metainfo: true,
        add_first_prefix_ion: true,
        add_losses: true,
        ..Default::default()
    };
    let ordinary = annotated(&generator.generate(&sequence("ASA"), 1, 1, None).unwrap());
    assert!(ordinary.contains_key("b2-H2O1+"));
    assert!(ordinary.contains_key("y2-H2O1+"));
    let tagged = sequence("AS[+0.001]A");
    let peaks = annotated(&generator.generate(&tagged, 1, 1, None).unwrap());
    assert_eq!(peaks.len(), 4);
    assert!(peaks.keys().all(|name| !name.contains('-')));
    near(peaks["b1+"], ordinary["b1+"]);
    near(peaks["y1+"], ordinary["y1+"]);
    near(peaks["b2+"], ordinary["b2+"] + 0.001);
    near(peaks["y2+"], ordinary["y2+"] + 0.001);
}

#[test]
fn small_complementary_fragment_survives_a_large_finite_tag() {
    // A suffix that does not contain the tag has its own known mass. Computing
    // it as total minus a huge prefix loses all of that mass to cancellation.
    // TSG.cpp's suffix loop accumulates from the C terminus independently.
    let peptide = sequence(&format!("X[1{}]A", "0".repeat(200)));
    let ions = peptide.fragment_ions(2).unwrap();
    for charge in 1..=2 {
        let ion = ions
            .iter()
            .find(|ion| ion.series == openms::chemistry::IonSeries::Y && ion.charge == charge)
            .unwrap();
        near(ion.mz, sequence("A").mz(i32::from(charge)).unwrap());
    }
}

#[test]
fn cancelling_terminal_mass_tags_preserve_the_peptide_mass() {
    // AASequence.cpp applies the terminal deltas before summing residue masses.
    // Equal opposite tags must not erase the known residue sum by cancellation.
    let mass = format!("1{}", "0".repeat(200));
    let peptide = sequence(&format!(".[+{mass}]ACDK.[-{mass}]"));
    near(
        peptide.mono_mass().unwrap(),
        sequence("ACDK").mono_mass().unwrap(),
    );
}

#[test]
fn unknown_chemistry_and_formula_dependent_generation_fail_atomically() {
    let mut existing = MSSpectrum {
        peaks: vec![Peak1D::new(90.0, 12.0)],
        ..Default::default()
    };
    existing.metadata.insert("keep".into(), "unchanged".into());
    let before = existing.clone();
    for bare in ["ABK", "AZK", "AXK", "(Acetyl)AXK"] {
        let generator = TheoreticalSpectrumGenerator::default();
        assert!(
            generator
                .append_to(&mut existing, &sequence(bare), 1, 2, None)
                .is_err()
        );
        assert_eq!(existing, before, "failed append with {bare}");
    }
    for generator in [
        TheoreticalSpectrumGenerator {
            isotope_model: TheoreticalIsotopeModel::Coarse { max_peaks: 3 },
            ..Default::default()
        },
        TheoreticalSpectrumGenerator {
            add_precursor_peaks: true,
            ..Default::default()
        },
    ] {
        let tagged = sequence("PEPTX[999]IDE");
        assert!(generator.generate(&tagged, 1, 2, None).is_err());
        assert!(
            generator
                .append_to(&mut existing, &tagged, 1, 2, None)
                .is_err()
        );
        assert_eq!(existing, before);
    }
}

#[cfg(feature = "idxml")]
#[test]
fn idxml_preserves_ambiguous_and_mass_only_sequence_identity() {
    use openms::format::idxml::{self, IdXmlDocument};
    use std::io::Cursor;

    // Independent schema-shaped input. Reading this string exercises the
    // parser directly; it does not depend on the Rust writer's chosen syntax.
    let xml = r#"<IdXML version="1.5">
<SearchParameters id="SP" db="synthetic" db_version="1" taxonomy="" mass_type="monoisotopic" charges="2" enzyme="Trypsin" missed_cleavages="0" precursor_peak_tolerance="0.01" peak_mass_tolerance="0.01"/>
<IdentificationRun date="2026-09-10T12:00:00" search_engine="synthetic" search_engine_version="1" search_parameters_ref="SP">
<ProteinIdentification score_type="score" higher_score_better="true" significance_threshold="0"/>
<PeptideIdentification score_type="score" higher_score_better="true" significance_threshold="0">
<PeptideHit score="3" sequence="BZXK" charge="2"/>
<PeptideHit score="2" sequence="PEPTX[999]IDE" charge="2"/>
<PeptideHit score="1" sequence=".[+1234.56789]ACDK.[+987.654321]" charge="2"/>
</PeptideIdentification>
</IdentificationRun></IdXML>"#;
    let mut document = idxml::read(Cursor::new(xml)).unwrap();
    let hits = &document.peptide_identifications[0].hits;
    assert_eq!(hits[0].sequence.as_str(), "BZXK");
    assert!(hits[0].sequence.mono_mass().is_err());
    near(
        hits[1].sequence.mono_mass().unwrap(),
        sequence("PEPTIDE").mono_mass().unwrap() + 999.0,
    );
    assert!(hits[1].sequence.formula().is_err());
    assert!(hits[2].sequence.average_mass().is_err());
    // Native numeric attachment is position-stable: removing A must not
    // reinterpret the formerly internal Q tag as N-terminal pyro-Glu.
    let sliced = sequence("AQ[111]AR").subsequence(1..4).unwrap();
    assert!(sliced.n_terminal_modification().is_none());
    assert!(sliced.residue_modification(0).unwrap().is_some());
    assert!(sliced.formula().is_err());
    document.peptide_identifications[0].hits.push(PeptideHit {
        sequence: sliced,
        score: 0.5,
        charge: 2,
        ..Default::default()
    });
    let mut encoded = Vec::new();
    idxml::write(&mut encoded, &document).unwrap();
    let roundtrip: IdXmlDocument = idxml::read(Cursor::new(encoded)).unwrap();
    assert_eq!(document, roundtrip);
}

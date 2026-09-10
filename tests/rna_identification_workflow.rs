// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{
    NAFragmentType, NASequence, NucleicAcidSpectrumGenerator, RNaseDigestion, Ribonucleotide,
    RibonucleotideDB, RibonucleotideRecord,
};
use openms::identification::graph::*;

fn parent(accession: &str, sequence: &str) -> ParentSequence {
    let mut parent = ParentSequence::new(accession);
    parent.molecule_type = MoleculeType::RNA;
    parent.sequence = sequence.into();
    parent
}

#[test]
fn shared_products_keep_both_parents_processing_history_and_fragment_chemistry() {
    let mut graph = IdentificationData::new().unwrap();
    let first = graph
        .register_parent_sequence(parent("first", "AGUC"))
        .unwrap();
    let second = graph
        .register_parent_sequence(parent("second", "UC"))
        .unwrap();
    let software = graph
        .register_processing_software(ProcessingSoftware::new("RNase_T1", "1"))
        .unwrap();
    let step = graph
        .register_processing_step(ProcessingStep::new(software), None)
        .unwrap();
    graph.set_current_processing_step(step).unwrap();
    RNaseDigestion::default()
        .digest_identification_data(&mut graph)
        .unwrap();
    assert_eq!(graph.oligo_count(), 2);
    let (id, oligo) = graph
        .oligos()
        .find(|(_, o)| o.sequence.to_string() == "UC")
        .unwrap();
    assert_eq!(
        oligo.parent_matches[&first]
            .iter()
            .map(|m| (m.start_pos, m.end_pos))
            .collect::<Vec<_>>(),
        [(Some(2), Some(3))]
    );
    assert_eq!(
        oligo.parent_matches[&second]
            .iter()
            .map(|m| (m.start_pos, m.end_pos))
            .collect::<Vec<_>>(),
        [(Some(0), Some(1))]
    );
    assert_eq!(oligo.result.steps_and_scores[0].processing_step, Some(step));
    let generator = NucleicAcidSpectrumGenerator {
        add_metainfo: true,
        add_first_prefix_ion: true,
        ..Default::default()
    };
    let spectrum = generator.generate(&oligo.sequence, -1, -1).unwrap();
    assert_eq!(spectrum.string_data_arrays[0].data.len(), 2);
    let before = oligo.sequence.mono_mass(NAFragmentType::Full, -1).unwrap();
    let (mut copied, translator) = graph.try_clone_with_translation().unwrap();
    let copied_id = translator.oligo(id).unwrap();
    assert!(copied.oligo(id).is_err());
    drop(graph);
    let retained = &copied.oligo(copied_id).unwrap().sequence;
    assert_eq!(
        retained
            .mono_mass(NAFragmentType::Full, -1)
            .unwrap()
            .to_bits(),
        before.to_bits()
    );
    assert_eq!(
        generator.generate(retained, -1, -1).unwrap().peaks,
        spectrum.peaks
    );
    copied.calculate_coverages(true).unwrap();
    assert_eq!(
        copied
            .parent(translator.parent(first).unwrap())
            .unwrap()
            .coverage,
        1.0
    );
    assert_eq!(
        copied
            .parent(translator.parent(second).unwrap())
            .unwrap()
            .coverage,
        1.0
    );
}

#[test]
fn all_parents_share_product_and_residue_limits_even_when_products_deduplicate() {
    for product_limit in [true, false] {
        let mut graph = IdentificationData::new().unwrap();
        let first = graph.register_parent_sequence(parent("a", "G")).unwrap();
        graph.register_parent_sequence(parent("b", "G")).unwrap();
        let mut digestion = RNaseDigestion::default();
        if product_limit {
            digestion.max_products = 1;
        } else {
            digestion.max_residues = 1;
        }
        assert_eq!(
            digestion
                .digest(&NASequence::parse("G").unwrap())
                .unwrap()
                .len(),
            1
        );
        assert!(digestion.digest_identification_data(&mut graph).is_err());
        assert_eq!(graph.oligo_count(), 0);
        assert_eq!(graph.parent(first).unwrap().sequence, "G");
        assert_eq!(graph.parent_count(), 2);
    }
}

fn custom_registry(code: &str) -> RibonucleotideDB {
    let mut records = RibonucleotideDB::global()
        .entries()
        .iter()
        .map(|r| r.as_ref().clone())
        .collect::<Vec<_>>();
    records.push(
        Ribonucleotide::from_record(RibonucleotideRecord {
            code: code.into(),
            origin: 'G',
            formula: "C10H13N5O4S".parse().unwrap(),
            mono_mass: 299.068825,
            average_mass: 299.308,
            ..Default::default()
        })
        .unwrap(),
    );
    RibonucleotideDB::from_records(records).unwrap()
}

#[test]
fn custom_chemistry_survives_registry_drop_and_flanks_use_raw_code_initials() {
    let registry = custom_registry("m1G");
    let mut graph = IdentificationData::new().unwrap();
    let parent = graph
        .register_parent_sequence(parent("custom", "A[m1G]UC"))
        .unwrap();
    RNaseDigestion::default()
        .digest_identification_data_with_registry(&mut graph, &registry)
        .unwrap();
    drop(registry);
    let (_, suffix) = graph
        .oligos()
        .find(|(_, o)| o.sequence.to_string() == "UC")
        .unwrap();
    let matched = suffix.parent_matches[&parent].first().unwrap();
    assert_eq!((matched.start_pos, matched.end_pos), (Some(2), Some(3)));
    assert_eq!(
        (&*matched.left_neighbor, &*matched.right_neighbor),
        ("m", "]")
    );
    let (_, prefix) = graph
        .oligos()
        .find(|(_, o)| o.sequence.to_string() == "A[m1G]p")
        .unwrap();
    assert_eq!(
        prefix.sequence.residues()[1].formula(),
        &"C10H13N5O4S".parse().unwrap()
    );
    assert_eq!(prefix.sequence.residues()[1].mono_mass(), 299.068825);
}

#[test]
fn unrepresentable_byte_flank_rolls_back_earlier_registrations() {
    let registry = custom_registry("éG");
    let mut graph = IdentificationData::new().unwrap();
    let first = graph.register_parent_sequence(parent("a", "G")).unwrap();
    graph
        .register_parent_sequence(parent("b", "A[éG]UC"))
        .unwrap();
    let error = RNaseDigestion::default()
        .digest_identification_data_with_registry(&mut graph, &registry)
        .unwrap_err();
    assert!(matches!(error, openms::Error::Unsupported(_)));
    assert_eq!(graph.oligo_count(), 0);
    assert_eq!(graph.parent(first).unwrap().sequence, "G");
}

#[test]
fn parent_list_cannot_leave_parsing_an_unbounded_allocation_allowance() {
    let mut graph = IdentificationData::new().unwrap();
    let id = graph
        .register_parent_sequence(parent("large", &"A".repeat(100_000)))
        .unwrap();
    let mut digestion = RNaseDigestion::new("RNase_T1").unwrap();
    digestion.max_output_bytes = std::mem::size_of::<ParentId>();
    assert!(digestion.digest_identification_data(&mut graph).is_err());
    assert_eq!(graph.oligo_count(), 0);
    assert_eq!(graph.parent(id).unwrap().sequence.len(), 100_000);
}

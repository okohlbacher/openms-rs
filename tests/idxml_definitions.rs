// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "idxml")]
use openms::chemistry::{
    AASequence, ModificationRecord, ModificationsDB, ResidueModification, TermSpecificity,
};
use openms::format::{
    idxml::{self, IdXmlDocument},
    modification_definitions as defs,
};
use openms::identification::{PeptideHit, PeptideIdentification, ProteinIdentification};
fn custom(
    name: &str,
    origin: Option<char>,
    term: TermSpecificity,
    delta: f64,
) -> ResidueModification {
    ResidueModification::from_record(ModificationRecord {
        name: name.into(),
        full_name: format!("Full {name}"),
        origin,
        term_specificity: term,
        diff_mono_mass: delta,
        diff_average_mass: delta + 0.01,
        ..Default::default()
    })
    .unwrap()
}
fn document() -> IdXmlDocument {
    let db = ModificationsDB::from_records(vec![
        custom("lab|end;\\one", None, TermSpecificity::NTerm, 8.5),
        custom("lab-mid", Some('M'), TermSpecificity::Anywhere, 12.5),
        custom("search-only", Some('S'), TermSpecificity::Anywhere, 7.5),
    ])
    .unwrap();
    let sequence =
        AASequence::parse_with_registry("(lab|end;\\one)AM(lab-mid)K[999]", &db).unwrap();
    let mut run = ProteinIdentification {
        identifier: "run".into(),
        date_time: Some("2026-09-10T12:00:00".into()),
        ..Default::default()
    };
    run.search_parameters
        .variable_modifications
        .push("search-only (S)".into());
    // Supply only the search-space definition explicitly; hit definitions must
    // be recovered from their owned Arc handles after this registry is dropped.
    defs::attach(
        &mut run.search_parameters,
        &[db.get_modification_handle("search-only (S)", None, None)
            .unwrap()],
    )
    .unwrap();
    IdXmlDocument {
        protein_identifications: vec![run],
        peptide_identifications: vec![PeptideIdentification {
            identifier: "run".into(),
            hits: vec![PeptideHit {
                sequence,
                ..Default::default()
            }],
            ..Default::default()
        }],
        ..Default::default()
    }
}
#[test]
fn custom_definitions_make_default_idxml_roundtrip_self_contained() {
    let input = document();
    let before = input.clone();
    let mut bytes = Vec::new();
    idxml::write(&mut bytes, &input).unwrap();
    assert_eq!(input, before);
    let restored = idxml::read(bytes.as_slice()).unwrap();
    assert_eq!(
        restored.peptide_identifications,
        input.peptide_identifications
    );
    let text = restored.protein_identifications[0]
        .search_parameters
        .metadata[defs::METADATA_KEY]
        .as_str()
        .unwrap();
    assert_eq!(
        ResidueModification::split_definition_records(text)
            .unwrap()
            .len(),
        3
    );
    let mut empty = ModificationsDB::default();
    assert_eq!(defs::register_from(text, &mut empty).unwrap().registered, 3);
    let second = idxml::read_with_registry(bytes.as_slice(), &Default::default(), &empty).unwrap();
    assert_eq!(restored, second);
    let mut again = Vec::new();
    idxml::write(&mut again, &restored).unwrap();
    assert_eq!(bytes, again);
}
#[test]
fn definitions_are_resolved_before_hits_and_do_not_mutate_the_global_registry() {
    let input = document();
    let mut bytes = Vec::new();
    idxml::write(&mut bytes, &input).unwrap();
    assert!(
        ModificationsDB::global()
            .find("lab-mid", None, None)
            .is_empty()
    );
    let restored = idxml::read(bytes.as_slice()).unwrap();
    drop(bytes);
    assert!(
        ModificationsDB::global()
            .find("lab-mid", None, None)
            .is_empty()
    );
    let original = input.peptide_identifications[0].hits[0]
        .sequence
        .mono_mass()
        .unwrap();
    assert_eq!(
        restored.peptide_identifications[0].hits[0]
            .sequence
            .mono_mass()
            .unwrap(),
        original
    );
    assert!(
        restored.peptide_identifications[0].hits[0]
            .sequence
            .formula()
            .is_err()
    );
}
#[test]
fn strict_custom_conflict_and_nonportable_definition_leave_writer_untouched() {
    let mut input = document();
    let mut bytes = vec![99];
    let conflicting = ModificationsDB::from_records(vec![custom(
        "lab-mid",
        Some('M'),
        TermSpecificity::Anywhere,
        99.0,
    )])
    .unwrap();
    assert!(
        idxml::write_with_registry(&mut bytes, &input, &Default::default(), &conflicting).is_err()
    );
    assert_eq!(bytes, [99]);
    input.protein_identifications[0]
        .search_parameters
        .metadata
        .insert(defs::METADATA_KEY.into(), "broken record".into());
    assert!(idxml::write(&mut bytes, &input).is_err());
    assert_eq!(bytes, [99]);
}
#[test]
fn shared_xml_decodes_latin1_utf16_and_rejects_missing_attribute_separators() {
    let input = IdXmlDocument {
        document_id: "é".into(),
        protein_identifications: vec![ProteinIdentification {
            identifier: "run".into(),
            date_time: Some("2026-09-10T12:00:00".into()),
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut output = Vec::new();
    idxml::write(&mut output, &input).unwrap();
    let xml = String::from_utf8(output).unwrap();
    let latin = xml
        .replace("UTF-8", "ISO-8859-1")
        .chars()
        .map(|c| c as u8)
        .collect::<Vec<_>>();
    assert_eq!(idxml::read(latin.as_slice()).unwrap(), input);
    let utf = xml.replace("UTF-8", "UTF-16");
    for little in [false, true] {
        let mut bytes = if little {
            vec![0xff, 0xfe]
        } else {
            vec![0xfe, 0xff]
        };
        for unit in utf.encode_utf16() {
            bytes.extend(if little {
                unit.to_le_bytes()
            } else {
                unit.to_be_bytes()
            });
        }
        assert_eq!(idxml::read(bytes.as_slice()).unwrap(), input);
    }
    assert!(idxml::read(xml.replace(" xmlns:xsi=", "xmlns:xsi=").as_bytes()).is_err());
}

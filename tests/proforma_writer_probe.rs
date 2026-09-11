// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::chemistry::proforma::*;

fn scored(tag: ModificationTag, value: f64) -> Modification {
    Modification {
        alternatives: vec![(
            tag,
            Some(Label {
                label_type: LabelType::Ambiguous,
                identifier: "g1".into(),
                score: Some(value),
            }),
        )],
        ..Default::default()
    }
}
fn chain(modifications: Vec<Modification>) -> Peptidoform {
    Peptidoform {
        sequence: vec![SequenceSection::Element(SequenceElement {
            amino_acid: 'X',
            modifications,
        })],
        ..Default::default()
    }
}
fn info(text: &str) -> ModificationTag {
    ModificationTag::InfoTag(InfoTag { text: text.into() })
}

#[test]
fn compiled_source_writer_matches_all_160_numeric_and_stream_cases() {
    let mut cases = 0;
    for line in include_str!("data/proforma_writer_reference.tsv").lines() {
        let fields: Vec<_> = line.splitn(4, '\t').collect();
        let value = f64::from_bits(u64::from_str_radix(fields[0], 16).unwrap());
        let mode = if fields[1] == "0" {
            WriteMode::Lossless
        } else {
            WriteMode::Canonical
        };
        let kind: u8 = fields[2].parse().unwrap();
        let mut mods = Vec::new();
        if kind == 1 {
            mods.push(scored(info("before"), value));
        }
        mods.push(scored(
            ModificationTag::MassDelta(MassDelta {
                source: MassDeltaSource::Obs,
                mass: value,
                original_text: if kind == 2 {
                    "+001.2300".into()
                } else {
                    String::new()
                },
            }),
            value,
        ));
        mods.push(scored(info("after"), value));
        let mut p = chain(mods);
        let text = if kind == 3 {
            p.charge = Some(ChargeState::Simple(2));
            let mut second = chain(vec![scored(info("fresh"), value)]);
            second.charge = Some(ChargeState::Simple(-1));
            PeptidoformIon {
                name: Some("omitted".into()),
                chains: vec![p, second],
                charge: Some(ChargeState::Simple(3)),
                is_chimeric: true,
            }
            .to_text(mode)
            .unwrap()
        } else {
            p.to_text(mode).unwrap()
        };
        assert_eq!(
            text, fields[3],
            "bits={} mode={mode:?} kind={kind}",
            fields[0]
        );
        cases += 1;
    }
    assert_eq!(cases, 160);
}

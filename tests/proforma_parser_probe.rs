// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::chemistry::proforma::*;

fn unhex(text: &str) -> String {
    assert_eq!(text.len() % 2, 0);
    let bytes = (0..text.len())
        .step_by(2)
        .map(|index| u8::from_str_radix(&text[index..index + 2], 16).unwrap())
        .collect();
    String::from_utf8(bytes).unwrap()
}

#[test]
fn compiled_source_parser_matches_476_outputs_and_error_payloads() {
    let mut count = 0;
    let (mut successes, mut failures) = (0, 0);
    for line in include_str!("data/proforma_parser_reference.tsv").lines() {
        let fields: Vec<_> = line.split('\t').collect();
        let input = unhex(fields[0]);
        let actual = if fields[1] == "0" {
            Peptidoform::parse(&input).map(|p| {
                (
                    p.to_text(WriteMode::Lossless).unwrap(),
                    p.to_text(WriteMode::Canonical).unwrap(),
                )
            })
        } else {
            PeptidoformIon::parse(&input).map(|p| {
                (
                    p.to_text(WriteMode::Lossless).unwrap(),
                    p.to_text(WriteMode::Canonical).unwrap(),
                )
            })
        };
        if fields[2] == "ok" {
            let actual = actual.unwrap_or_else(|e| panic!("C++ accepted {input:?}: {e}"));
            assert_eq!(actual, (unhex(fields[3]), unhex(fields[4])), "{input:?}");
            successes += 1;
        } else {
            let error = match actual {
                Err(ParseFailure::Syntax(error)) => error,
                other => panic!("C++ rejected {input:?}: got {other:?}"),
            };
            assert_eq!(
                error.code() as usize,
                fields[3].parse::<usize>().unwrap(),
                "{input:?}"
            );
            assert_eq!(
                error.position(),
                fields[4].parse::<usize>().unwrap(),
                "{input:?}"
            );
            assert_eq!(error.message(), unhex(fields[5]), "{input:?}");
            failures += 1;
        }
        count += 1;
    }
    assert_eq!(count, 476);
    assert_eq!(successes + failures, count);
    assert!(successes > 300 && failures > 50);
}

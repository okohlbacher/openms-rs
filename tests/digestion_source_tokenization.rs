// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The source's own digestion arithmetic against the Linux x86_64 Release
//! build: `ProteaseDigestion::digest` (pair overload), `peptideCount` and
//! `EnzymaticDigestion::countInternalCleavageSites`, through
//! [`ProteaseDigestion::digest_ranges_source`],
//! [`ProteaseDigestion::peptide_count_source`] and
//! [`ProteaseDigestion::count_internal_cleavage_sites_source`].
//!
//! `data/digestion_source_oracle.tsv` is `../oracle/integ-B-digest-source`'s
//! `results/records.tsv` verbatim (ibminode06, the Release install
//! `openms4-release-bc9cc12-c19e494-174b576`, conda g++ 14.4.0, `-O3 -DNDEBUG`,
//! every enzyme in its own process, two byte-identical runs): eight enzymes
//! (Asp-N, Asp-N_ambic and Lys-N, whose expressions are lookahead-only;
//! Trypsin, Trypsin/P and Chymotrypsin; `unspecific cleavage`; `no cleavage`),
//! sixteen sequences including the empty one, specificities full, semi and
//! none, zero, one, two and five missed cleavages and four length filters.
//! Nothing in this file is derived from Rust output.

use openms::chemistry::{DigestionSpecificity, ProductValidation, Protease, ProteaseDigestion};

const ORACLE: &str = include_str!("data/digestion_source_oracle.tsv");

fn records(kind: &str) -> Vec<Vec<&'static str>> {
    ORACLE
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|fields| fields[0] == kind)
        .collect()
}

fn digestion(enzyme: &str, spec: &str, mc: &str, min: &str, max: &str) -> ProteaseDigestion {
    let max: usize = max.parse().unwrap();
    ProteaseDigestion {
        enzyme: Protease::from_name(enzyme).unwrap(),
        specificity: spec.parse::<DigestionSpecificity>().unwrap(),
        missed_cleavages: mc.parse().unwrap(),
        min_length: min.parse().unwrap(),
        // The source's 0 is "unbounded".
        max_length: (max != 0).then_some(max),
        ..ProteaseDigestion::default()
    }
}

fn pairs(text: &str) -> Vec<(usize, usize)> {
    if text.is_empty() {
        return Vec::new();
    }
    text.split(',')
        .map(|pair| {
            let (a, b) = pair.split_once('-').unwrap();
            (a.parse().unwrap(), b.parse().unwrap())
        })
        .collect()
}

#[test]
fn count_internal_cleavage_sites_is_the_token_count_minus_one() {
    let rows = records("count");
    assert_eq!(rows.len(), 8 * 16);
    let mut leading = 0;
    let mut wrapped = 0;
    for row in &rows {
        let (enzyme, sequence, expected) = (row[1], row[2], row[3]);
        let digestion = digestion(enzyme, "full", "0", "1", "0");
        let got = digestion
            .count_internal_cleavage_sites_source(sequence)
            .unwrap();
        assert_eq!(got.to_string(), expected, "{enzyme} {sequence:?}");
        if got == usize::MAX {
            wrapped += 1;
        }
        // The native count differs exactly where the source's tokenisation does.
        let native = digestion.count_internal_cleavage_sites(sequence).unwrap();
        if got != native && got != usize::MAX {
            assert_eq!(got, native + 1, "{enzyme} {sequence:?}");
            leading += 1;
        }
    }
    // Asp-N, Asp-N_ambic and Lys-N on sequences that start with a matching
    // residue; the empty sequence wraps for every enzyme but `no cleavage`.
    assert!(leading > 0, "{leading}");
    assert_eq!(wrapped, 7, "{wrapped}");
}

#[test]
fn digest_matches_the_release_build_in_products_order_and_discarded_count() {
    let rows = records("digest");
    assert_eq!(rows.len(), 8 * 3 * 4 * 16 * 4);
    let (mut compared, mut refused) = (0, 0);
    for row in &rows {
        let (enzyme, spec, mc, min, max, sequence) =
            (row[1], row[2], row[3], row[4], row[5], row[6]);
        let label = format!("{enzyme} {spec} mc={mc} len={min}..{max} {sequence:?}");
        let result = digestion(enzyme, spec, mc, min, max).digest_ranges_source(sequence);
        if row[7].starts_with("ERR:") {
            let error = result.expect_err(&label);
            if spec == "none" {
                // ProteaseDigestion.cpp:93-98: Exception::InvalidValue.
                assert!(
                    error.to_string().contains(
                        "Specificity value set on current ProteaseDigestion object is not \
                         supported by ProteaseDigestion::digest()."
                    ),
                    "{label}: {error}"
                );
            } else {
                // EnzymaticDigestion.cpp:121-126: fewer than two positions.
                assert!(sequence.is_empty(), "{label}");
                assert!(
                    error
                        .to_string()
                        .contains("Too few cleavage positions - at least sequence start and end positions are required."),
                    "{label}: {error}"
                );
            }
            refused += 1;
            continue;
        }
        let got = result.unwrap_or_else(|error| panic!("{label}: {error}"));
        assert_eq!(got.discarded.to_string(), row[7], "{label}: discarded");
        assert_eq!(got.ranges, pairs(row[8]), "{label}: products");
        compared += 1;
    }
    // 6,144 records: 2,048 of specificity none, which `digest` refuses, and 112
    // semi-specific digests of the empty sequence, which it refuses too.
    assert_eq!((compared, refused), (3984, 2160));
}

#[test]
fn peptide_count_matches_the_release_build() {
    let rows = records("pepcnt");
    assert_eq!(rows.len(), 8 * 3 * 4 * 16);
    let mut differs_from_native = 0;
    for row in &rows {
        let (enzyme, spec, mc, sequence, expected) = (row[1], row[2], row[3], row[4], row[5]);
        let label = format!("{enzyme} {spec} mc={mc} {sequence:?}");
        let digestion = digestion(enzyme, spec, mc, "1", "0");
        let got = digestion.peptide_count_source(sequence).unwrap();
        assert_eq!(got.to_string(), expected, "{label}");
        // The native count follows the generation settings; the source's does
        // not look at the specificity or the length filter.
        if digestion.peptide_count_unmodified(sequence).ok() != Some(got) {
            differs_from_native += 1;
        }
    }
    assert!(differs_from_native > 0);
}

#[test]
fn digest_unmodified_pair_overload_matches_the_release_build() {
    // EnzymaticDigestion::digestUnmodified returns 0 and no product before it
    // tokenises when the sequence is shorter than min_length, enumerates every
    // substring for `unspecific cleavage` and specificity none, and otherwise
    // is digest's loops; the source writes (start, length).
    let (mut digested, mut enumerated, mut short) = (0, 0, 0);
    let rows = records("unmod");
    assert_eq!(rows.len(), 8 * 3 * 4 * 16 * 4);
    for row in rows {
        let (enzyme, spec, mc, min, max, sequence) =
            (row[1], row[2], row[3], row[4], row[5], row[6]);
        let label = format!("{enzyme} {spec} mc={mc} len={min}..{max} {sequence:?}");
        let got = digestion(enzyme, spec, mc, min, max)
            .digest_unmodified_source(sequence)
            .unwrap_or_else(|error| panic!("{label}: {error}"));
        assert_eq!(got.discarded.to_string(), row[7], "{label}: discarded");
        let converted: Vec<(usize, usize)> = got
            .ranges
            .iter()
            .map(|&(start, end)| (start, end - start))
            .collect();
        assert_eq!(converted, pairs(row[8]), "{label}: products");
        let minimum: usize = min.parse().unwrap();
        if sequence.len() < minimum {
            short += 1;
        } else if spec == "none" || enzyme == "unspecific cleavage" {
            enumerated += 1;
        } else {
            digested += 1;
        }
    }
    assert!(digested > 2000 && enumerated > 1000 && short > 500);
}

#[test]
fn unspecific_cleavage_orders_by_length_where_the_native_digest_orders_by_start() {
    // ProteaseDigestion.cpp:106: the missed-cleavage ceiling is the largest
    // Size, so the loops emit every product of one length before the next.
    let unspecific = ProteaseDigestion {
        enzyme: Protease::UnspecificCleavage,
        ..ProteaseDigestion::default()
    };
    let source = unspecific.digest_ranges_source("PEPT").unwrap();
    assert_eq!(
        source.ranges,
        [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (0, 2),
            (1, 3),
            (2, 4),
            (0, 3),
            (1, 4),
            (0, 4)
        ]
    );
    let mut native: Vec<(usize, usize)> = unspecific
        .digest_ranges("PEPT")
        .unwrap()
        .iter()
        .map(|p| (p.start, p.end))
        .collect();
    assert_ne!(native, source.ranges);
    native.sort_by_key(|&(start, end)| (end - start, start));
    assert_eq!(native, source.ranges);
}

#[test]
fn is_valid_product_matches_the_release_build_outside_the_sequence_too() {
    // isValidProduct_ over pos = -1..=n+1 and length = -1..=n+1 for two
    // sequences, with and without the missed-cleavage check, the N-terminal
    // methionine allowance and the random D|P allowance. A cell is the answer
    // and the warning the call wrote to OPENMS_LOG_WARN; XX is a cell the
    // driver did not run, where the source is undefined.
    let rows = records("valid");
    assert!(rows.len() > 500, "{}", rows.len());
    let (mut cells, mut skipped) = (0, 0);
    let mut kinds = [0_usize; 4];
    let mut accepted = 0_usize;
    for row in rows {
        let (enzyme, spec, mc, sequence) = (row[1], row[2], row[3], row[4]);
        let (ignore, nterm, asp) = (row[5] == "1", row[6] == "1", row[7] == "1");
        let digestion = digestion(enzyme, spec, mc, "1", "0");
        let options = ProductValidation {
            ignore_missed_cleavages: ignore,
            allow_nterm_protein_cleavage: nterm,
            allow_random_asp_pro_cleavage: asp,
        };
        let n = sequence.len() as i32;
        let width = (n + 3) as usize;
        let text = row[8].as_bytes();
        assert_eq!(text.len(), 2 * width * width, "{row:?}");
        for (index, cell) in text.chunks(2).enumerate() {
            let pos = (index / width) as i32 - 1;
            let length = (index % width) as i32 - 1;
            let label = format!(
                "{enzyme} {spec} mc={mc} {sequence} ignore={ignore} nterm={nterm} asp={asp} pos={pos} len={length}"
            );
            if cell == b"XX" {
                skipped += 1;
                continue;
            }
            let got = digestion
                .is_valid_product_source(sequence, pos, length, options)
                .unwrap_or_else(|error| panic!("{label}: {error}"));
            assert_eq!(got.valid, cell[0] == b'1', "{label}");
            // The warning's numbers are those after the allowance.
            let (adjusted_pos, adjusted_length) = if nterm && pos <= 2 && sequence.starts_with('M')
            {
                (0, length + pos)
            } else {
                (pos, length)
            };
            let expected: Vec<String> = match cell[1] {
                b'-' => Vec::new(),
                b's' => vec![format!(
                    "Warning: start of fragment ({adjusted_pos}) is beyond end of sequence '{sequence}'!"
                )],
                b'e' => vec![format!(
                    "Warning: end of fragment ({}) is beyond end of sequence '{sequence}'!",
                    adjusted_pos + adjusted_length
                )],
                b'z' => vec!["Warning: fragment and sequence must not be empty!".to_owned()],
                other => panic!("{label}: unexpected warning kind {}", other as char),
            };
            assert_eq!(got.warnings, expected, "{label}");
            kinds[match cell[1] {
                b'-' => 0,
                b's' => 1,
                b'e' => 2,
                _ => 3,
            }] += 1;
            if got.valid {
                accepted += 1;
            }
            cells += 1;
        }
    }
    assert!(kinds.iter().all(|&count| count > 100), "{kinds:?}");
    assert!(accepted > 1000, "{accepted}");
    assert!(cells > 20_000 && skipped > 1_000, "{cells} {skipped}");
}

#[test]
fn the_source_applies_the_methionine_allowance_before_it_looks_at_the_length() {
    // ProteaseDigestion.h: "Regard peptide as n-terminal of protein if it
    // starts only at pos=1 or 2 and protein starts with 'M'". isValidProduct_
    // resets (pos 1, length 0) to (0, 1) before it tests the length, so
    // that fragment is the valid first token of MDKPR under Asp-N (a cut
    // before the D); the native range API cannot express a length of zero and
    // answers false for the empty range 1..1.
    let asp_n = digestion("Asp-N", "full", "0", "1", "0");
    let options = ProductValidation {
        allow_nterm_protein_cleavage: true,
        ..ProductValidation::default()
    };
    let result = asp_n
        .is_valid_product_source("MDKPR", 1, 0, options)
        .unwrap();
    assert!(result.valid);
    assert!(result.warnings.is_empty());
    assert!(
        !asp_n
            .is_valid_product_unmodified("MDKPR", 1..1, options)
            .unwrap()
    );
    // Without the allowance a length of zero is the source's empty fragment.
    let plain = asp_n
        .is_valid_product_source("MDKPR", 1, 0, ProductValidation::default())
        .unwrap();
    assert!(!plain.valid);
    assert_eq!(
        plain.warnings,
        ["Warning: fragment and sequence must not be empty!"]
    );
}

#[test]
fn specificity_values_names_and_lookup_are_the_source_enumeration() {
    // EnzymaticDigestion.h:41-49 (values) and EnzymaticDigestion.cpp:20 (names)
    // at bc9cc12; tier 3.
    let table = [
        (DigestionSpecificity::None, 0, "none"),
        (DigestionSpecificity::Semi, 1, "semi"),
        (DigestionSpecificity::Full, 2, "full"),
        (DigestionSpecificity::Unknown, 3, "unknown"),
        (DigestionSpecificity::NoCTerm, 8, "no-cterm"),
        (DigestionSpecificity::NoNTerm, 9, "no-nterm"),
    ];
    for (value, number, name) in table {
        assert_eq!(value.value(), number);
        assert_eq!(value.name(), name);
        assert_eq!(DigestionSpecificity::from_value(number), Some(value));
        assert_eq!(DigestionSpecificity::by_name(name), value);
    }
    assert_eq!(
        DigestionSpecificity::NAMES_OF_SPECIFICITY,
        [
            "none", "semi", "full", "unknown", "unknown", "unknown", "unknown", "unknown",
            "no-cterm", "no-nterm"
        ]
    );
    for number in [4, 5, 6, 7, 10, 255] {
        assert_eq!(DigestionSpecificity::from_value(number), None);
    }
    // getSpecificityByName answers SPEC_UNKNOWN for a name it does not know.
    assert_eq!(
        DigestionSpecificity::by_name("Full"),
        DigestionSpecificity::Unknown
    );
    assert_eq!(
        DigestionSpecificity::by_name(""),
        DigestionSpecificity::Unknown
    );
    // FromStr stays the three configuration names.
    assert!("unknown".parse::<DigestionSpecificity>().is_err());
    assert_eq!(openms::chemistry::digestion::NO_CLEAVAGE, "no cleavage");
    assert_eq!(
        openms::chemistry::digestion::UNSPECIFIC_CLEAVAGE,
        "unspecific cleavage"
    );
}

#[test]
fn the_values_digest_does_not_support_are_refused_and_validation_treats_them_as_full() {
    // ProteaseDigestion.cpp:93-98 refuses what is neither full nor semi;
    // EnzymaticDigestion::isValidProduct_ (:211-303) tests `== SPEC_NONE`,
    // `== SPEC_FULL` for no cleavage and `== SPEC_SEMI`, so every other value
    // is the full test.
    let protein: openms::chemistry::AASequence = "DFPIANGERDK".parse().unwrap();
    let text = protein.as_str().to_owned();
    for spec in [
        DigestionSpecificity::Unknown,
        DigestionSpecificity::NoCTerm,
        DigestionSpecificity::NoNTerm,
    ] {
        let digestion = ProteaseDigestion {
            enzyme: Protease::AspN,
            specificity: spec,
            ..ProteaseDigestion::default()
        };
        assert!(digestion.digest(&protein).is_err(), "{spec}");
        assert!(digestion.digest_ranges_source(&text).is_err(), "{spec}");
        let full = ProteaseDigestion {
            specificity: DigestionSpecificity::Full,
            ..digestion.clone()
        };
        for start in 0..text.len() {
            for end in start + 1..=text.len() {
                let options = openms::chemistry::ProductValidation::default();
                assert_eq!(
                    digestion
                        .is_valid_product_unmodified(&text, start..end, options)
                        .unwrap(),
                    full.is_valid_product_unmodified(&text, start..end, options)
                        .unwrap(),
                    "{spec} {start}..{end}"
                );
            }
        }
    }
}

#[test]
fn filter_by_missed_cleavages_hands_the_filter_the_source_count_as_an_int() {
    // The oracle's count rows: Asp-N on DFPIANGERDK is 2 (the leading site is
    // counted), Trypsin on the empty sequence is Size(0) - 1.
    let asp_n = digestion("Asp-N", "full", "0", "1", "0");
    let mut seen = None;
    let filtered = asp_n
        .filter_by_missed_cleavages("DFPIANGERDK", |count| {
            seen = Some(count);
            count > 1
        })
        .unwrap();
    assert_eq!(seen, Some(2));
    assert!(filtered);
    let trypsin = digestion("Trypsin", "full", "0", "1", "0");
    let mut seen = None;
    let filtered = trypsin
        .filter_by_missed_cleavages("", |count| {
            seen = Some(count);
            count > 1
        })
        .unwrap();
    // 0 - 1 in Size arithmetic, converted to Int.
    assert_eq!(seen, Some(-1));
    assert!(!filtered);
}

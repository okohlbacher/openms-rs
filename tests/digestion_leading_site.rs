// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `ProteaseDigestion` and `EnzymaticDigestion` against the Linux x86_64
//! Release build, with and without a cleavage site at the protein's first
//! residue.
//!
//! Asp-N `(?=[DBX])`, Asp-N_ambic `(?=[DBEZX])` and Lys-N `(?=[KX])` are
//! lookahead-only (`share/OpenMS/CHEMISTRY/Enzymes.xml:43`, `:67`, `:135` at
//! core `bc9cc12`), so they match with zero width at position 0 of a protein
//! that starts with one of those residues; Trypsin `(?<=[KRX])(?!P)` (`:6`)
//! cannot. `EnzymaticDigestion::tokenize_` (`EnzymaticDigestion.cpp:92-116`)
//! then pushes position 0 twice. The source's own `isValidProduct` is
//! unaffected (it only counts sites strictly inside a product), but
//! `countInternalCleavageSites` (`:195-198`) counts one more site,
//! `ProteaseDigestion::digest` (`ProteaseDigestion.cpp:87-158`) discards an
//! empty product, repeats the first product as a "missed-cleavage" one and,
//! semi-specifically, loses products its `isValidProduct` accepts, and
//! `peptideCount` (`:48-69`) counts the empty product. The port keeps each
//! boundary once (`docs/DIGESTION_SUPPORT.md`, "Leading cleavage site"), so
//! these tests pin what it shares with the Release build and exactly where it
//! differs.
//!
//! `data/digestion_leading_site_oracle.tsv` is `../oracle/integ-B-digest-leading`'s
//! `results/records.tsv` (ibminode06, the Release install
//! `openms4-release-bc9cc12-c19e494-174b576`, conda g++ 14.4.0, `-O3
//! -DNDEBUG`, two byte-identical runs), with the `isValidProduct` answers
//! packed one line per case by that oracle's `make_fixture.py`.

use openms::chemistry::{ProductValidation, Protease, ProteaseDigestion};
use std::collections::BTreeSet;

const ORACLE: &str = include_str!("data/digestion_leading_site_oracle.tsv");

type Range = (usize, usize);

fn records(kind: &str) -> Vec<Vec<&'static str>> {
    ORACLE
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|fields| fields[0] == kind)
        .collect()
}

/// Whether the source's expression matches with zero width at position 0.
fn leading_site(enzyme: &str, sequence: &str) -> bool {
    let class: &[u8] = match enzyme {
        "Asp-N" => b"DBX",
        "Asp-N_ambic" => b"DBEZX",
        "Lys-N" => b"KX",
        _ => b"",
    };
    sequence
        .as_bytes()
        .first()
        .is_some_and(|first| class.contains(first))
}

fn digestion(enzyme: &str, specificity: &str, missed_cleavages: &str) -> ProteaseDigestion {
    ProteaseDigestion {
        enzyme: Protease::from_name(enzyme).unwrap(),
        specificity: specificity.parse().unwrap(),
        missed_cleavages: missed_cleavages.parse().unwrap(),
        ..ProteaseDigestion::default()
    }
}

/// `start-end` pairs (`digest`) or `start-length` pairs (`digestUnmodified`).
fn pairs(text: &str) -> Vec<Range> {
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

fn without_repeats(ranges: &[Range]) -> Vec<Range> {
    let mut seen = BTreeSet::new();
    ranges.iter().copied().filter(|r| seen.insert(*r)).collect()
}

fn port_ranges(digestion: &ProteaseDigestion, sequence: &str) -> Vec<Range> {
    digestion
        .digest_ranges(sequence)
        .unwrap()
        .into_iter()
        .map(|p| (p.start, p.end))
        .collect()
}

/// The ranges the Release build's `isValidProduct(.., ignore_missed_cleavages
/// = false)` accepts, from the packed `valid` records.
fn release_valid_set(enzyme: &str, specificity: &str, mc: &str, sequence: &str) -> BTreeSet<Range> {
    let row = records("valid")
        .into_iter()
        .find(|f| {
            f[1] == enzyme && f[2] == specificity && f[3] == mc && f[4] == sequence && f[5] == "0"
        })
        .unwrap();
    let n = sequence.len();
    let mut answers = row[6].bytes();
    let mut set = BTreeSet::new();
    for pos in 0..n {
        for len in 1..=n - pos {
            if answers.next().unwrap() == b'1' {
                set.insert((pos, pos + len));
            }
        }
    }
    assert!(answers.next().is_none());
    set
}

#[test]
fn fixture_holds_every_driver_case() {
    assert_eq!(records("count").len(), 40);
    assert_eq!(records("digest").len(), 240);
    assert_eq!(records("unmod").len(), 240);
    assert_eq!(records("pepcnt").len(), 240);
    assert_eq!(records("valid").len(), 480);
    let leading = records("count")
        .iter()
        .filter(|f| leading_site(f[1], f[2]))
        .count();
    assert_eq!(leading, 11);
}

#[test]
fn product_validity_matches_the_release_build() {
    let mut checked = 0_usize;
    for f in records("valid") {
        let d = digestion(f[1], f[2], f[3]);
        let sequence = f[4];
        let options = ProductValidation {
            ignore_missed_cleavages: f[5] == "1",
            ..ProductValidation::default()
        };
        let mut answers = f[6].bytes();
        for pos in 0..sequence.len() {
            for len in 1..=sequence.len() - pos {
                let release = answers.next().unwrap() == b'1';
                assert_eq!(
                    d.is_valid_product_unmodified(sequence, pos..pos + len, options)
                        .unwrap(),
                    release,
                    "{f:?} pos {pos} len {len}"
                );
                checked += 1;
            }
        }
        assert!(answers.next().is_none());
    }
    assert_eq!(checked, 17_376);
}

#[test]
fn internal_cleavage_count_leaves_out_only_the_leading_zero_width_site() {
    for f in records("count") {
        let d = digestion(f[1], "full", "0");
        let release: usize = f[3].parse().unwrap();
        let port = d.count_internal_cleavage_sites(f[2]).unwrap();
        assert_eq!(
            port + usize::from(leading_site(f[1], f[2])),
            release,
            "{f:?}"
        );
    }
}

#[test]
fn full_digestion_is_the_release_order_with_repeats_removed() {
    for f in records("digest").into_iter().filter(|f| f[2] == "full") {
        let (enzyme, sequence) = (f[1], f[4]);
        let release = pairs(f[6]);
        let discarded: usize = f[5].parse().unwrap();
        let leading = leading_site(enzyme, sequence);
        // The source discards exactly one product for length: the empty one
        // between the two copies of position 0.
        assert_eq!(discarded, usize::from(leading), "{f:?}");
        let expected = without_repeats(&release);
        if !leading {
            assert_eq!(expected, release, "{f:?}");
        }
        assert_eq!(
            port_ranges(&digestion(enzyme, "full", f[3]), sequence),
            expected,
            "{f:?}"
        );
    }
}

#[test]
fn digestion_is_the_release_validity_set() {
    for f in records("digest") {
        let (enzyme, specificity, mc, sequence) = (f[1], f[2], f[3], f[4]);
        let release = pairs(f[6]);
        let port = port_ranges(&digestion(enzyme, specificity, mc), sequence);
        let port_set: BTreeSet<Range> = port.iter().copied().collect();
        assert_eq!(port_set.len(), port.len(), "{f:?}");
        assert_eq!(
            port_set,
            release_valid_set(enzyme, specificity, mc, sequence),
            "{f:?}"
        );
        if leading_site(enzyme, sequence) {
            // The source's product list is a subset: semi-specifically it
            // loses products its own isValidProduct accepts.
            assert!(release.iter().all(|r| port_set.contains(r)), "{f:?}");
        } else {
            assert_eq!(port, release, "{f:?}");
        }
    }
}

#[test]
fn unmodified_digestion_matches_the_release_build() {
    let digests = records("digest");
    for (f, d) in records("unmod").into_iter().zip(digests) {
        assert_eq!(f[1..5], d[1..5]);
        let release: Vec<Range> = pairs(f[6])
            .into_iter()
            .map(|(start, length)| (start, start + length))
            .collect();
        // Both source paths (digestUnmodified and ProteaseDigestion::digest)
        // give the same products and discarded count.
        assert_eq!(release, pairs(d[6]), "{f:?}");
        assert_eq!(f[5], d[5], "{f:?}");
        let sequence = f[4];
        let port = digestion(f[1], f[2], f[3])
            .digest_unmodified(sequence)
            .unwrap();
        if f[2] == "full" || !leading_site(f[1], sequence) {
            let expected: Vec<&str> = without_repeats(&release)
                .into_iter()
                .map(|(start, end)| &sequence[start..end])
                .collect();
            assert_eq!(port, expected, "{f:?}");
        }
    }
}

#[test]
fn peptide_count_counts_the_ports_products() {
    for f in records("pepcnt") {
        let (enzyme, specificity, mc, sequence) = (f[1], f[2], f[3], f[4]);
        let release: usize = f[5].parse().unwrap();
        let d = digestion(enzyme, specificity, mc);
        let protein = sequence.parse().unwrap();
        let port = d.peptide_count(&protein).unwrap();
        assert_eq!(port, port_ranges(&d, sequence).len(), "{f:?}");
        let digest = records("digest")
            .into_iter()
            .find(|g| g[1..5] == f[1..5])
            .unwrap();
        let discarded: usize = digest[5].parse().unwrap();
        if specificity == "full" {
            // The source's count is its product list plus the discarded empty
            // product, repeats included.
            assert_eq!(release, pairs(digest[6]).len() + discarded, "{f:?}");
            if !leading_site(enzyme, sequence) {
                assert_eq!(port, release, "{f:?}");
            }
        } else {
            // The source ignores semi-specificity (docs/DIGESTION_SUPPORT.md).
            let full = records("pepcnt")
                .into_iter()
                .find(|g| g[1] == enzyme && g[2] == "full" && g[3] == mc && g[4] == sequence)
                .unwrap();
            assert_eq!(f[5], full[5], "{f:?}");
        }
    }
}

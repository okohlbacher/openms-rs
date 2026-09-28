// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Tier 1: the chemistry value members lane A5a closed, against the executed
//! C++ Release build.
//!
//! Every expected value is a record of `../oracle/A5a-chem/results/records.tsv`
//! (sha256 `8de8db89310d5c776c0b908936d04e8394a338fb817f49f2e4bea50736f59263`),
//! written by `driver.cpp` (sha256 `939f86c7...0664`) built with
//! `x86_64-conda-linux-gnu-g++ 14.4.0` and libOpenMS's core flags against the
//! Linux x86_64 Release install `openms4-release-bc9cc12-c19e494-174b576`
//! (libOpenMS.so sha256 `abd4fc99...6d6c`) and run twice on ibminode06; both
//! runs agreed byte for byte. Hashes are the source's `std::hash`
//! specialisations; doubles are compared by their bits.

use openms::chemistry::element_db::{atomic_numbers, get_element, names, symbols};
use openms::chemistry::isotopes::{
    CoarseIsotopePatternGenerator as Generator, CoarseMassMode, IsotopeDistribution, IsotopePeak,
    ProbabilityPrecision,
};
use openms::chemistry::{
    EmpiricalFormula, Ribonucleotide, RibonucleotideRecord, RibonucleotideTermSpecificity,
};

fn ef(text: &str) -> EmpiricalFormula {
    EmpiricalFormula::parse(text).unwrap()
}

fn bits(value: f64) -> u64 {
    value.to_bits()
}

#[test]
fn element_db_sizes_and_lookups() {
    // ElementDB.names.size / symbols.size / atomic_numbers.size
    assert_eq!(names().len(), 367);
    assert_eq!(symbols().len(), 369);
    assert_eq!(atomic_numbers().len(), 84);
    for (query, expected) in [
        ("C", "Carbon C 6"),
        ("Carbon", "Carbon C 6"),
        ("(13)C", "(13)Carbon (13)C 6"),
        ("(13)Carbon", "(13)Carbon (13)C 6"),
        ("D", "(2)Hydrogen (2)H 1"),
        ("T", "(3)Hydrogen (3)H 1"),
        ("(2)H", "(2)Hydrogen (2)H 1"),
        ("(3)Hydrogen", "(3)Hydrogen (3)H 1"),
    ] {
        let found = get_element(query).unwrap();
        assert_eq!(
            format!(
                "{} {} {}",
                found.name(),
                found.symbol(),
                found.atomic_number()
            ),
            expected,
            "{query}"
        );
    }
    assert!(get_element("Xx").is_none());
    // The SDK builds iridium from rhenium's tables (CPP-249), so it has
    // (185)Ir where the port, with the declared iridium table, has (191)Ir.
    assert!(get_element("(191)Ir").is_some());
    assert!(get_element("(185)Ir").is_none());
}

#[test]
fn element_stream_text_and_hash() {
    for (symbol, text, hash) in [
        (
            "H",
            "Hydrogen H 1 1.00794 1.00783 1.0078250319=99.9885% 2.01410178=0.0115%",
            0xca09fec9c39f3305_u64,
        ),
        (
            "C",
            "Carbon C 6 12.0107 12 12.0=98.93% 13.003355000000001=1.07%",
            0xd576291dcc8c6af0,
        ),
        (
            "O",
            "Oxygen O 8 15.9994 15.9949 15.994915000000001=99.757% 16.999131999999999=0.038% 17.999168999999998=0.205%",
            0x4fd02469f85c982c,
        ),
        (
            "S",
            "Sulfur S 16 32.0661 31.9721 31.972070729999999=94.93% 32.971457999999998=0.76% 33.967866999999998=4.29% 35.967081=0.02%",
            0x02605cd955bcea84,
        ),
        (
            "(13)C",
            "(13)Carbon (13)C 6 13.0034 13.0034 13.003355000000001=100%",
            0x9f29bb17d72a64e8,
        ),
        (
            "Pt",
            "Platinum Pt 78 195.051 194.965 191.961038000000002=0.782% 193.962680299999988=32.967% 194.964791100000014=33.832% 195.964951500000012=25.242% 197.967893000000004=7.163%",
            0xd834a9845ebf251f,
        ),
        (
            "Re",
            "Rhenium Re 75 186.207 186.956 184.952955000000003=37.4% 186.95575310000001=62.6%",
            0x6a03104a950ce90a,
        ),
        (
            "U",
            "Uranium U 92 238.029 238.051 234.040950000000009=0.0054% 235.043927999999994=0.7204% 238.050790000000006=99.2742%",
            0x5a3df0d8da4e5218,
        ),
    ] {
        let element = get_element(symbol).unwrap();
        assert_eq!(element.source_stream_text(), text, "{symbol}");
        assert_eq!(element.source_hash64(), hash, "{symbol}");
    }
    for (symbol, mono, average) in [
        ("H", 0x3ff0200d24019c4a_u64, 0x3ff020867bd8f13c_u64),
        ("C", 0x4028000000000000, 0x4028057f2cf9e6bf),
        ("O", 0x402ffd657fb69985, 0x402fffb20df4d74e),
        ("S", 0x403ff8d9a09abfb4, 0x4040087576f06a61),
        ("(13)C", 0x402a01b7bf1e8e61, 0x402a01b7bf1e8e61),
        ("Pt", 0x40685edf9195bf1a, 0x406861a38a78fe4d),
        ("U", 0x406dc1a012599ed8, 0x406dc0ecd8ee7c8f),
    ] {
        let element = get_element(symbol).unwrap();
        assert_eq!(bits(element.mono_mass()), mono, "{symbol} mono");
        assert_eq!(bits(element.average_mass()), average, "{symbol} average");
    }
    let (h, c, c13, s) = (
        get_element("H").unwrap(),
        get_element("C").unwrap(),
        get_element("(13)C").unwrap(),
        get_element("S").unwrap(),
    );
    assert!(h < c);
    // `*s < *c` is false in the SDK, as is `*(13)C < *C`.
    assert!(s >= c);
    assert!(c < c13);
    assert!(c13 >= c);
}

#[test]
fn empirical_formula_text_hash_and_isotope_number() {
    // (formula, toString, operator<<, std::hash, calculateTheoreticalIsotopesNumber)
    for (text, to_string, streamed, hash, number) in [
        ("C2H5", "C2H5", "C2H5", 0xbde37b11a9de0a2c_u64, Some(10.0)),
        ("CH4+", "C1H4", "CH4+", 0xe54214f2a647147c, Some(12.0)),
        ("C2+3", "C2", "C2+3", 0xd208804ec185128d, Some(1.0)),
        ("C2-", "C2", "C2-", 0xbe7bdc230b6e8b78, Some(1.0)),
        ("H-2O", "H-2O1", "HO", 0xf397c9137907bbd0, None),
        (
            "C2(13)C",
            "(13)C1C2",
            "(13)CC2",
            0x22e1600ac3cbde1f,
            Some(1.0),
        ),
        ("C5D5", "(2)H5C5", "(2)H5C5", 0xb2ac83827c18be77, Some(50.0)),
        (
            "C6H12O6",
            "C6H12O6",
            "C6H12O6",
            0x84867025d9133fbe,
            Some(19_800.0),
        ),
        ("", "", "", 0xeb5cf0391d18900a, Some(1.0)),
        ("O2", "O2", "O2", 0x72cde47971d9d8cb, None),
        ("C-2", "C-2", "C", 0x1c13f5af6be7c3e8, None),
        ("(13)C4", "(13)C4", "(13)C4", 0xccb8c1185ca66834, Some(4.0)),
        ("C", "C1", "C", 0xf3dcb7cc1c569c12, Some(2.0)),
    ] {
        let formula = if text.is_empty() {
            EmpiricalFormula::default()
        } else {
            ef(text)
        };
        assert_eq!(formula.to_source_string(), to_string, "{text}");
        assert_eq!(formula.source_stream_text(), streamed, "{text}");
        assert_eq!(formula.source_hash64(), hash, "{text}");
        match number {
            Some(value) => assert_eq!(
                formula.theoretical_isotopes_number().unwrap(),
                value,
                "{text}"
            ),
            // O2: Boost's domain error for k > n. H-2O and C-2: the SDK
            // converts the negative count to unsigned and returns 2.767e19
            // and 9.223e18; the port refuses a negative count.
            None => assert!(formula.theoretical_isotopes_number().is_err(), "{text}"),
        }
    }
    let charged = ef("C2").with_charge(-2);
    assert_eq!(charged.source_stream_text(), "C2--2");
    assert_eq!(charged.source_hash64(), 0x517fdc367ebf3ccb);

    // getLightestIsotopeWeight for single-element formulas, where the SDK's
    // element order cannot matter.
    for (text, expected) in [
        ("C2+3", 0x403b05969c906f85_u64),
        ("C2-", 0x4036fe2321253029),
        ("O2", 0x403ffd657fb69985),
        ("C-2", 0xc038000000000000),
        ("(13)C4", 0x404a01b7bf1e8e61),
        ("C", 0x4028000000000000),
    ] {
        assert_eq!(bits(ef(text).lightest_isotope_mass()), expected, "{text}");
    }
    assert_eq!(bits(EmpiricalFormula::default().lightest_isotope_mass()), 0);
    // Multi-element sums follow the SDK's address order; within a few ulps.
    for (text, expected) in [
        ("C2H5", 0x403d0a041b4080d7_u64),
        ("CH4+", 0x403109e027db36ea),
        ("C6H12O6", 0x406682074b44a039),
        ("C636H1007N167O178S5PtH-2", 0x40cbb3259a79a7d0),
    ] {
        let expected = f64::from_bits(expected);
        let actual = ef(text).lightest_isotope_mass();
        assert!(
            (actual - expected).abs() <= expected * 4.0 * f64::EPSILON,
            "{text}: {actual}"
        );
    }

    let mut adduct = ef("C2+");
    adduct.add_charge_adduct(1, "Na+").unwrap();
    assert_eq!(adduct.to_source_string(), "C2Na1");
    assert_eq!(adduct.charge(), 0);
    assert_eq!(
        EmpiricalFormula::water(3).unwrap().to_source_string(),
        "H6O3"
    );
    // The SDK's water(0) keeps H0 and O0 entries ("H0O0", not empty); the port
    // stores no zero counts.
    assert!(EmpiricalFormula::water(0).unwrap().is_empty());
    let c13 = get_element("(13)C").unwrap();
    let labelled = EmpiricalFormula::from_element(2, c13, 0);
    assert_eq!(labelled.to_source_string(), "(13)C2");
    assert_eq!(labelled, ef("(13)C2"));
    assert_eq!(labelled.number_of(c13), 2);
    assert_eq!(labelled.number_of(get_element("C").unwrap()), 0);
}

fn ribo(name: &str, code: &str, formula: &str, mono: f64, average: f64) -> RibonucleotideRecord {
    RibonucleotideRecord {
        name: name.into(),
        code: code.into(),
        new_code: code.into(),
        html_code: code.into(),
        formula: ef(formula),
        origin: 'A',
        mono_mass: mono,
        average_mass: average,
        term_specificity: RibonucleotideTermSpecificity::Anywhere,
        baseloss_formula: ef("C5H10O5"),
    }
}

#[test]
fn ribonucleotide_stream_text_and_hash() {
    let a = Ribonucleotide::from_record(ribo("adenosine", "A", "C10H13N5O4", 267.0968, 267.24))
        .unwrap();
    let pa = Ribonucleotide::from_record(ribo(
        "adenosine monophosphate",
        "pA",
        "C10H14N5O7P",
        347.063,
        347.22,
    ))
    .unwrap();
    let mut five_record = ribo("adenosine", "A", "C10H13N5O4", 267.0968, 267.24);
    five_record.term_specificity = RibonucleotideTermSpecificity::FivePrime;
    let five = Ribonucleotide::from_record(five_record).unwrap();
    let default = Ribonucleotide::default();
    assert_eq!(
        a.source_stream_text(),
        "Ribonucleotide 'A' (adenosine, C10H13N5O4)"
    );
    assert_eq!(
        pa.source_stream_text(),
        "Ribonucleotide 'pA' (adenosine monophosphate, C10H14N5O7P)"
    );
    assert_eq!(
        default.source_stream_text(),
        "Ribonucleotide '.' (unknown ribonucleotide, )"
    );
    assert_eq!(a.source_hash64(), 0x10c7dfd8ca4cd2de);
    assert_eq!(pa.source_hash64(), 0x03c43990fc7e3721);
    assert_eq!(five.source_hash64(), 0xb0ca63d1603b1e2d);
    assert_eq!(default.source_hash64(), 0x8d546d5fa808057d);
}

fn distribution(peaks: &[(f64, f64)]) -> IsotopeDistribution {
    IsotopeDistribution::from_peaks(
        peaks
            .iter()
            .map(|&(mass, probability)| IsotopePeak { mass, probability })
            .collect(),
    )
    .unwrap()
}

fn pairs(distribution: &IsotopeDistribution) -> Vec<(f64, f64)> {
    distribution
        .peaks()
        .iter()
        .map(|p| (p.mass, p.probability))
        .collect()
}

#[test]
fn isotope_distribution_merge_hash_and_order() {
    type Case = (
        &'static [(f64, f64)],
        f64,
        f64,
        Option<&'static [(f64, f64)]>,
    );
    let cases: [Case; 7] = [
        (
            &[(100.0, 0.5), (100.4, 0.25), (101.0, 0.125), (101.5, 0.125)],
            1.0,
            0.0,
            Some(&[(100.0, 0.75), (100.75, 0.125)]),
        ),
        (
            &[(101.0, 0.125), (100.0, 0.5), (101.5, 0.125), (100.4, 0.25)],
            1.0,
            0.0,
            Some(&[(100.0, 0.75), (100.75, 0.125)]),
        ),
        (
            &[(100.0, 0.5), (101.0, 0.25), (101.25, 0.125), (102.0, 0.125)],
            1.0,
            0.0,
            Some(&[(100.0, 0.5), (101.0, 0.375)]),
        ),
        (
            &[
                (99.0, 0.0625),
                (100.0, 0.5),
                (100.5, 0.25),
                (103.0, 0.25),
                (104.0, 0.0625),
            ],
            1.5,
            0.1,
            Some(&[(100.0, 0.75)]),
        ),
        (&[(100.0, 1.0)], 1.0, 0.0, Some(&[])),
        // "New Isotope Distribution has more points than the old one."
        (&[(100.0, 0.5), (100.5, 0.5)], 0.1, 0.0, None),
        (
            &[(100.0, 0.5), (100.2, 0.25), (100.4, 0.125), (103.0, 0.125)],
            1.0,
            0.0,
            Some(&[(100.0, 0.875), (0.0, 0.0), (0.0, 0.0)]),
        ),
    ];
    for (input, resolution, min_prob, expected) in cases {
        let mut merged = distribution(input);
        match expected {
            Some(peaks) => {
                merged.merge(resolution, min_prob).unwrap();
                assert_eq!(pairs(&merged), peaks, "{input:?}");
            }
            None => assert!(merged.merge(resolution, min_prob).is_err()),
        }
    }
    assert_eq!(
        IsotopeDistribution::default().source_hash64(),
        0xa8d8756188585c3d
    );
    assert_eq!(
        IsotopeDistribution::empty().source_hash64(),
        0xa8c7f832281a39c5
    );
    assert_eq!(
        distribution(&[(1.0, 0.5), (2.0, 0.25)]).source_hash64(),
        0x2c59cc9903467746
    );
    // The float intensity: 0.1 narrows to 0.1f.
    assert_eq!(
        distribution(&[(1.0, 0.1)]).source_hash64(),
        0x0ebdc28399a152c4
    );
    assert!(distribution(&[(5.0, 1.0)]) < distribution(&[(1.0, 0.5), (2.0, 0.25)]));
    assert!(distribution(&[(1.0, 0.8), (3.0, 0.1)]) < distribution(&[(1.0, 0.9), (2.0, 0.1)]));
}

#[test]
fn sulfur_fragment_estimate_matches_the_sdk_bits() {
    let generator = Generator::new(None, CoarseMassMode::Approximate)
        .unwrap()
        .with_precision(ProbabilityPrecision::SourceSingle);
    // (sulfur, raw [(mass bits, weight bits)], renormalized weight bits)
    type Case = (u32, [(u64, u64); 3], [u64; 3]);
    let expected: [Case; 4] = [
        (
            0,
            [
                (0x40590ae391a1bedd, 0x3f66976ec0000000),
                (0x40594b1a88d76777, 0x3f626012c0000000),
                (0x40598b51800d1010, 0x3f66976ec0000000),
            ],
            [0x3fd6bf9ea0000000, 0x3fd280c2a0000000, 0x3fd6bf9ea0000000],
        ),
        (
            1,
            [
                (0x4058ff15de87edf2, 0x3f65c472c0000000),
                (0x40593f4cd5bd968c, 0x3f5fc804c0000000),
                (0x40597f83ccf33f26, 0x3fa55fa720000000),
            ],
            [0x3fad5b5540000000, 0x3fa56e5ba0000000, 0x3fecd36500000000],
        ),
        (
            2,
            [
                (0x4059026f470dad28, 0x3f650d0600000000),
                (0x405942a63e4355c2, 0x3f5acc2860000000),
                (0x405982dd3578fe5c, 0x3fb39266a0000000),
            ],
            [0x3fa04fd500000000, 0x3f94c3c7e0000000, 0x3fee54e480000000],
        ),
        (
            3,
            [
                (0x4058fca40ab4298c, 0x3f64806f60000000),
                (0x40593cdb01e9d226, 0x3f505ed4a0000000),
                (0x40597d11f91f7ac0, 0x3fbc73c100000000),
            ],
            [0x3f965a94e0000000, 0x3f81d96ea0000000, 0x3fef05c5a0000000],
        ),
    ];
    for (sulfur, raw, renormalized) in expected {
        let mut pattern = generator
            .estimate_fragment_from_peptide_weight_and_sulfur(200.0, sulfur, 100.0, sulfur, &[2])
            .unwrap();
        let actual: Vec<(u64, u64)> = pattern
            .peaks()
            .iter()
            .map(|p| (bits(p.mass), bits(p.probability)))
            .collect();
        assert_eq!(actual, raw, "sulfur {sulfur} raw");
        pattern
            .renormalize_with(ProbabilityPrecision::SourceSingle)
            .unwrap();
        let weights: Vec<u64> = pattern
            .peaks()
            .iter()
            .map(|p| bits(p.probability))
            .collect();
        assert_eq!(weights, renormalized, "sulfur {sulfur} renormalized");
    }
}

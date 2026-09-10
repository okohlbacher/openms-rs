// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::ion_naming::{
    MAX_ION_NAME_BYTES, MAX_REPEATED_SIGNS, charge_from_name, charge_suffix, ordinal_from_name,
    with_charge,
};

// Literal assertions from OpenMS4-core rev7c029e8 IonNaming_test.cpp. They remain
// explicit here so these tests work without access to the C++ source checkout.
#[test]
fn source_charge_suffix_literals_include_extreme_and_compact_charges() {
    assert_eq!(MAX_REPEATED_SIGNS, 8);
    for (charge, expected) in [
        (0, ""),
        (1, "+"),
        (2, "++"),
        (5, "+++++"),
        (8, "++++++++"),
        (-1, "-"),
        (-3, "---"),
        (9, "+9"),
        (-12, "-12"),
        (2_000_000_000, "+2000000000"),
        (i32::MIN, "-2147483648"),
    ] {
        let suffix = charge_suffix(charge);
        assert_eq!(suffix, expected);
        assert!(suffix.len() <= 11);
    }
}

#[test]
fn source_charge_parsing_literals_cover_all_three_forms() {
    for (name, expected) in [
        ("", 0),
        ("y5", 0),
        ("[alpha|ci$y3]", 0),
        ("a2-B", 0),
        ("my note", 0),
        ("y1-H2O1", 0),
        ("b3-H3N1", 0),
        ("y5+", 1),
        ("b3++", 2),
        ("y1-H2O1+", 1),
        ("[alpha|ci$y3]++", 2),
        ("c1-", -1),
        ("a3-B-", -1),
        ("w12--", -2),
        ("y5+3", 3),
        ("y5-3", -3),
        ("b2^2", 2),
        ("y3-H2O^2", 2),
        ("y4-H2O1^2/3.2ppm", 2),
        ("y3^2*0.75", 2),
        ("c1^-1", -1),
        ("y4-H2O1^2/-1", 2),
        ("y4^2/-1", 2),
        ("note ^ see docs", 0),
        ("y5++\nsome user comment", 2),
        ("y5\nsome user comment", 0),
        ("y5+99999999999999999999", 0),
        ("y5^99999999999999999999", 0),
        ("y5+2147483648", 0),
        ("y5-2147483649", 0),
        ("y3/1.2e-05", 0),
        ("y3/-1", 0),
        ("y3*0.75", 0),
        ("y3^1/1.2e-05", 1),
        ("y3^2/-1", 2),
        ("note ^2text", 0),
        ("y4^2abc", 0),
    ] {
        assert_eq!(charge_from_name(name), expected, "{name:?}");
    }
}

#[test]
fn source_with_charge_literals_preserve_existing_priority_and_free_text() {
    for (name, charge, expected) in [
        ("y3+", 1, "y3+"),
        ("b2++", 2, "b2++"),
        ("c1-", -1, "c1-"),
        ("y3^2", 2, "y3^2"),
        ("y1-H2O1+", 1, "y1-H2O1+"),
        ("[alpha|ci$y3]", 2, "[alpha|ci$y3]++"),
        ("y3", 1, "y3+"),
        ("w1", -1, "w1-"),
        ("y3+", 2, "y3+"),
        ("y3", 0, "y3"),
        ("", 0, ""),
        ("y3\nmy comment", 1, "y3+\nmy comment"),
        ("y3\r\nmy comment", 2, "y3++\r\nmy comment"),
        ("y3+\nmy comment", 1, "y3+\nmy comment"),
    ] {
        assert_eq!(with_charge(name, charge).unwrap(), expected);
    }
}

#[test]
fn source_ordinal_literals_do_not_merge_charge_or_loss_digits() {
    for (name, expected) in [
        ("y3", 3),
        ("y12++", 12),
        ("b7+", 7),
        ("y3-H2O1+", 3),
        ("y3+12", 3),
        ("y3^2", 3),
        ("", 0),
        ("y", 0),
        ("[alpha|ci$y3]++", 0),
        ("iY+U-H3PO4+", 0),
        ("+3", 0),
        ("-2", 0),
        ("12", 0),
        ("y99999999999999", 0),
    ] {
        assert_eq!(ordinal_from_name(name), expected, "{name:?}");
    }
}

#[test]
fn caret_precedence_zero_tokens_and_malformed_fallback_match_source_branches() {
    for (name, expected) in [
        ("y3^+2", 2),
        ("y3^2^3", 3),
        ("y3^2^oops", 0),
        ("y3^2^oops--", -2),
        ("y3^2^oops+4", 4),
        ("y3^x/--", -2),
        ("y3^x/-2", 0),
        ("y3^x*++", 2),
        ("y3^x*+2", 0),
        ("y3^", 0),
        ("y3^+", 1),
        ("y3^-", -1),
        ("y3^--", -2),
        ("y3^2abc+3", 3),
        ("y3^2abc/--", -2),
        ("y3^2abc/-2", 0),
        ("y3^2147483648/-1", 0),
        ("y3^2147483648/--", -2),
        ("y3^0/-2", 0),
        ("y3^0/--", 0),
        ("y3^-0*++", 0),
        ("y3^2/--", 2),
        ("y3^2*++", 2),
        ("y3++---", -3),
        ("y3--++", 2),
    ] {
        assert_eq!(charge_from_name(name), expected, "{name:?}");
    }
    // Zero caret charge is unknown, yet remains the highest-priority token:
    // source appends a suffix without rewriting it or claiming parse idempotence.
    let named = with_charge("y3^0/--", 2).unwrap();
    assert_eq!(named, "y3^0/--++");
    assert_eq!(charge_from_name(&named), 0);
}

#[test]
fn numeric_lengths_and_i32_boundaries_are_checked_without_normalizing_tokens() {
    for (token, expected) in [
        ("+2147483647", i32::MAX),
        ("-2147483648", i32::MIN),
        ("+2147483648", 0),
        ("-2147483649", 0),
        ("+9999999999", 0),
        ("-9999999999", 0),
        ("+0000000001", 1),
        ("-0000000001", -1),
        ("+00000000001", 0),
        ("-00000000001", 0),
        ("+0000000000", 0),
    ] {
        for prefix in ["y3", "y3^"] {
            assert_eq!(charge_from_name(&format!("{prefix}{token}")), expected);
        }
    }
    assert_eq!(ordinal_from_name("Y999999999+"), 999_999_999);
    assert_eq!(ordinal_from_name("z000000001^2"), 1);
    assert_eq!(ordinal_from_name("a0000000001+"), 0);
    assert_eq!(ordinal_from_name("A4294967295"), 0); // The source limit is nine digits.
}

#[test]
fn ascii_grammar_preserves_unicode_nuls_and_exact_first_line_boundaries() {
    for (name, expected) in [
        ("α3++", 2),
        ("🧬y3^-2/注記", -2),
        ("y3^２", 0),
        ("y3+٣", 0),
        ("y3−2", 0),
        ("y3＋2", 0),
        ("y3\0+2", 2),
        ("y3+\rcomment--", 1),
        ("y3\rcomment^2", 0),
        ("y3-\r\ncomment++", -1),
        ("y3\n\rcomment++", 0),
        ("\n++", 0),
        ("\r^2", 0),
    ] {
        assert_eq!(charge_from_name(name), expected, "{name:?}");
    }
    for (name, expected) in [
        ("α12++", 0),
        ("é12", 0),
        ("y２", 0),
        ("y٣", 0),
        ("y12🧬", 12),
        ("a12\0", 12),
        ("Y12\r\nannotation^8", 12),
        ("y\n12", 0),
        ("\ny12", 0),
    ] {
        assert_eq!(ordinal_from_name(name), expected, "{name:?}");
    }
    for line_ending in ["\r", "\n", "\r\n", "\n\r"] {
        let name = format!("α3🧬{line_ending}free\0text^9");
        assert_eq!(
            with_charge(&name, -2).unwrap(),
            format!("α3🧬--{line_ending}free\0text^9")
        );
    }
    assert_eq!(with_charge("\r\ncomment", 1).unwrap(), "+\r\ncomment");
}

#[test]
fn emitted_suffixes_round_trip_for_small_extreme_and_distributed_charges() {
    for charge in (-4096..=4096).chain([i32::MIN, i32::MAX, 2_000_000_000, -2_000_000_000]) {
        let name = format!("y5{}", charge_suffix(charge));
        assert_eq!(charge_from_name(&name), charge);
        assert_eq!(ordinal_from_name(&name), 5);
        assert_eq!(with_charge(&name, charge).unwrap(), name);
    }
    // Deterministic full-width values include both signs and large magnitudes;
    // widening arithmetic is independently exercised beyond the source loop.
    let mut bits = 0_u32;
    for _ in 0..10_000 {
        bits = bits.wrapping_add(0x9e37_79b9);
        let charge = i32::from_ne_bytes(bits.to_ne_bytes());
        let name = with_charge("b19\r\nfree text", charge).unwrap();
        assert_eq!(charge_from_name(&name), charge);
        assert_eq!(with_charge(&name, charge.wrapping_add(1)).unwrap(), name);
    }
}

#[test]
fn parsers_handle_long_tokens_while_copy_limits_include_unchanged_and_multiline_names() {
    let digits = "0".repeat(MAX_ION_NAME_BYTES + 1);
    assert_eq!(charge_from_name(&format!("y3+{digits}")), 0);
    assert_eq!(charge_from_name(&format!("y3^{digits}/-2")), 0);
    assert_eq!(charge_from_name(&format!("y3^{digits}/---")), -3);
    assert_eq!(ordinal_from_name(&format!("y{digits}")), 0);
    let signs = "+".repeat(MAX_ION_NAME_BYTES + 1);
    assert_eq!(charge_from_name(&signs), (MAX_ION_NAME_BYTES + 1) as i32);

    let at_limit = "A".repeat(MAX_ION_NAME_BYTES);
    assert_eq!(with_charge(&at_limit, 0).unwrap().len(), MAX_ION_NAME_BYTES);
    assert!(with_charge(&at_limit, 1).is_err());
    assert!(with_charge(&digits, 0).is_err());
    assert!(with_charge(&signs, 1).is_err()); // Existing charge does not bypass cap.
    let fits = "α".repeat((MAX_ION_NAME_BYTES - 2) / 2);
    assert_eq!(with_charge(&fits, 2).unwrap().len(), MAX_ION_NAME_BYTES);
    assert!(with_charge(&fits, 3).is_err());
    let multiline = format!("y3\r\n{}", "A".repeat(MAX_ION_NAME_BYTES - 5));
    assert_eq!(
        with_charge(&multiline, 1).unwrap().len(),
        MAX_ION_NAME_BYTES
    );
    assert!(with_charge(&multiline, 2).is_err());
}

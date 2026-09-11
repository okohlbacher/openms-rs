// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::MassDecomposition as MD;
use openms::chemistry::mass_decomposition::{
    MAX_MASS_DECOMPOSITION_INPUT_BYTES, MAX_MASS_DECOMPOSITION_OUTPUT_BYTES,
};
use std::cmp::Ordering;

#[test]
fn source_construction_copy_assignment_and_compact_literals() -> openms::Result<()> {
    let empty = MD::new();
    assert_eq!(empty.number_of_max_aa(), 0);
    assert_eq!(empty.to_text()?, "");
    assert_eq!(empty.to_expanded_string()?, "");
    for (text, max) in [("C3", 3), ("C3 M4", 4), ("C3 M4 S200", 200)] {
        let md: MD = text.parse()?;
        assert_eq!(md.number_of_max_aa(), max);
        assert_eq!(md.to_text()?, text);
        let cloned = md.clone();
        assert_eq!(cloned, md);
        let mut assigned = MD::new();
        assigned.clone_from(&md);
        assert_eq!(assigned, md);
        assert_eq!(assigned.number_of_max_aa(), max);
    }
    Ok(())
}

#[test]
fn source_addition_maximum_and_expanded_literals() -> openms::Result<()> {
    let (c, m, s) = (MD::parse("C3")?, MD::parse("M4")?, MD::parse("S200")?);
    let mut result = MD::new();
    for (rhs, text, max) in [(&c, "C3", 3), (&m, "C3 M4", 4), (&s, "C3 M4 S200", 200)] {
        result.checked_add_assign(rhs)?;
        assert_eq!(result.to_text()?, text);
        assert_eq!(result.number_of_max_aa(), max);
    }
    assert_eq!(MD::new().checked_add(&c)?, c);
    assert_eq!(c.checked_add(&m)?.to_text()?, "C3 M4");
    assert_eq!(c.checked_add(&m)?.checked_add(&s)?, result);
    assert_eq!(c.to_expanded_string()?, "CCC");
    assert_eq!(c.checked_add(&m)?.to_expanded_string()?, "CCCMMMM");
    assert_eq!(c.checked_add(&c)?.to_text()?, "C6");
    assert_eq!(c.to_text()?, "C3");
    Ok(())
}

#[test]
fn source_order_and_string_equality_literals() -> openms::Result<()> {
    let (mut md, c, m) = (MD::new(), MD::parse("C3")?, MD::parse("M4")?);
    assert!(!md.equals_text(&c.to_text()?)?);
    md.checked_add_assign(&c)?;
    assert!(!m.source_less(&c));
    assert!(md.source_less(&m));
    md.checked_add_assign(&m)?;
    assert!(md.source_less(&m));
    md = m.clone();
    assert!(md.equals_text(&m.to_text()?)?);
    let s = MD::parse("S200")?;
    md = m.checked_add(&s)?;
    assert!(!md.equals_text(&s.to_text()?)?);
    assert!(MD::parse("A2")?.source_less(&MD::parse("A10")?));
    Ok(())
}

#[test]
fn source_tag_and_compatibility_literals() -> openms::Result<()> {
    let (empty, c, cm, cms) = (
        MD::new(),
        MD::parse("C3")?,
        MD::parse("C3 M4")?,
        MD::parse("C3 M4 S200")?,
    );
    assert!(!empty.contains_tag("C")?);
    assert!(!empty.contains_tag("CCC")?);
    assert!(c.contains_tag("CCC")?);
    assert!(!c.contains_tag("CCCC")?);
    assert!(cm.contains_tag("CMC")?);
    assert!(cms.contains_tag("CCCSSMSSSSSSSSSSSSSSM")?);
    assert!(empty.compatible(&MD::parse("")?));
    assert!(!empty.compatible(&MD::parse("C1")?));
    assert!(c.compatible(&MD::parse("C1")?));
    assert!(cm.compatible(&MD::parse("C2 M4")?));
    assert!(!cm.compatible(&MD::parse("C2 M5")?));
    for rhs in ["C3 S200", "C3 M4", "S2", "M4 S200"] {
        assert!(cms.compatible(&MD::parse(rhs)?));
    }
    Ok(())
}

#[test]
fn duplicate_history_affects_equality_but_not_map_order_and_right_cache_is_ignored()
-> openms::Result<()> {
    let historical = MD::parse("A9 A1")?;
    let plain = MD::parse("A1")?;
    assert_eq!(historical.to_text()?, "A1");
    assert_eq!(historical.to_expanded_string()?, "A");
    assert_eq!(historical.number_of_max_aa(), 9);
    assert_ne!(historical, plain);
    assert_eq!(historical.source_cmp(&plain), Ordering::Equal);
    assert!(!historical.source_less(&plain));
    assert!(historical.equals_text("A9 A1")?);
    assert!(!historical.equals_text("A1")?);
    assert_eq!(MD::new().checked_add(&historical)?, plain);
    assert_eq!(historical.checked_add(&MD::new())?, historical);
    assert_eq!(historical.checked_add(&plain)?.number_of_max_aa(), 9);
    assert_eq!(plain.checked_add(&historical)?.number_of_max_aa(), 2);
    assert_eq!(
        historical
            .checked_add(&MD::parse("A20")?)?
            .number_of_max_aa(),
        21
    );
    Ok(())
}

#[test]
fn source_zero_keys_suffix_trim_and_non_amino_acid_byte_symbols() -> openms::Result<()> {
    let zero = MD::parse("Z0")?;
    assert_eq!(zero.to_text()?, "Z0");
    assert_eq!(zero.to_expanded_string()?, "");
    assert!(!MD::new().compatible(&zero));
    assert!(zero.compatible(&MD::new()));
    assert!(zero.contains_tag("")?);
    assert!(!zero.contains_tag("Z")?);
    assert!(MD::parse("A2 B1")?.contains_tag("ABA")?);
    assert!(MD::parse("A2 B1")?.contains_tag("BAA")?);
    assert!(!MD::parse("A2")?.contains_tag("Å")?);
    assert_eq!(MD::parse(" B2 A1 \t(annotation)")?.to_text()?, "A1 B2");
    assert_eq!(MD::parse("(anything, even malformed counts)")?, MD::new());
    assert_eq!(MD::parse("A+2 B-0")?.to_text()?, "A2 B0");
    assert_eq!(MD::parse("12 +3 a1")?.to_text()?, "+3 12 a1");
    assert_eq!(MD::parse("12 +3 a1")?.to_expanded_string()?, "+++11a");
    assert_eq!(MD::parse("\0\t2")?.to_expanded_string()?, "\0\0");
    // Source trims the complete compact representation, even a leading symbol.
    assert_eq!(MD::parse("\t2 A1")?.to_text()?, "2 A1");
    assert_eq!(MD::parse("\t2 A1")?.to_expanded_string()?, "\t\tA");
    Ok(())
}

#[test]
fn malformed_counts_and_source_i32_parse_boundary_are_checked() -> openms::Result<()> {
    for text in [
        " ",
        " A1",
        "A1 ",
        "A1  B1",
        "A",
        "A-1",
        "A+-1",
        "A++1",
        "A1.0",
        "A2147483648",
        "A-2147483649",
        "A1\tB2",
        "Å2",
    ] {
        assert!(MD::parse(text).is_err(), "{text:?}");
        assert!(MD::new().equals_text(text).is_err(), "{text:?}");
    }
    let md = MD::parse("A2147483647")?.checked_add(&MD::parse("A1")?)?;
    assert_eq!(md.to_text()?, "A2147483648");
    assert_eq!(md.number_of_max_aa(), 2_147_483_648);
    assert!(md.equals_text(&md.to_text()?).is_err()); // Source constructor still narrows to i32.
    Ok(())
}

#[test]
fn unsigned_addition_overflow_is_atomic_even_after_an_earlier_valid_key() -> openms::Result<()> {
    let one = MD::parse("Z1")?;
    let mut limit = one.clone();
    for _ in 1..usize::BITS {
        limit = limit.checked_add(&limit)?.checked_add(&one)?;
    }
    assert_eq!(limit.number_of_max_aa(), usize::MAX);
    assert_eq!(limit.to_text()?, format!("Z{}", usize::MAX));
    let saved = limit.clone();
    assert!(limit.checked_add_assign(&MD::parse("A1 Z1")?).is_err());
    assert_eq!(limit, saved);
    assert!(limit.checked_add(&one).is_err());
    assert_eq!(limit, saved);
    assert!(
        limit
            .checked_add(&MD::parse("A1")?)?
            .to_expanded_string()
            .is_err()
    ); // total length overflow
    Ok(())
}

#[test]
fn input_tag_work_and_expanded_output_limits_precede_large_work() -> openms::Result<()> {
    let allowed = format!("A1 ({}", "x".repeat(MAX_MASS_DECOMPOSITION_INPUT_BYTES - 4));
    assert_eq!(MD::parse(&allowed)?, MD::parse("A1")?);
    assert!(MD::parse(&(allowed + "x")).is_err());
    let oversize_tag = "A".repeat(MAX_MASS_DECOMPOSITION_INPUT_BYTES + 1);
    assert!(
        MD::parse("A2147483647")?
            .contains_tag(&oversize_tag)
            .is_err()
    );
    let large = MD::parse(&format!("A{}", MAX_MASS_DECOMPOSITION_OUTPUT_BYTES + 1))?;
    assert!(large.to_expanded_string().is_err());
    assert_eq!(
        large.to_text()?,
        format!("A{}", MAX_MASS_DECOMPOSITION_OUTPUT_BYTES + 1)
    );
    let mut tokens = vec!["A2147483647".into()];
    for byte in b'!'..=b'~' {
        if !matches!(byte, b'(' | b'A') {
            tokens.push(format!("{}0", char::from(byte)));
        }
    }
    let many_keys = MD::parse(&tokens.join(" "))?;
    let large_tag = "A".repeat(MAX_MASS_DECOMPOSITION_INPUT_BYTES);
    assert!(many_keys.contains_tag(&large_tag).is_err());
    assert!(many_keys.contains_tag("AAA")?);
    Ok(())
}

#[test]
fn source_plus_and_plus_assign_differ_for_later_new_keys() -> openms::Result<()> {
    // cpp operator+ compares each new count with original-left maximum; +=
    // compares with the running maximum. This is a finite observable difference.
    let left = MD::parse("A1")?;
    let rhs = MD::parse("B10 C5")?;
    let added = left.checked_add(&rhs)?;
    let mut assigned = left.clone();
    assigned.checked_add_assign(&rhs)?;
    assert_eq!(added.to_text()?, "A1 B10 C5");
    assert_eq!(added.number_of_max_aa(), 5);
    assert!(!added.equals_text(&added.to_text()?)?);
    assert_eq!(assigned.number_of_max_aa(), 10);
    assert_ne!(added, assigned);
    assert_eq!(added.source_cmp(&assigned), Ordering::Equal);
    // A new key can also lower a maximum produced by an earlier existing key.
    let combined = MD::parse("A10")?.checked_add(&MD::parse("A100 B11")?)?;
    assert_eq!(combined.to_text()?, "A110 B11");
    assert_eq!(combined.number_of_max_aa(), 11);
    // The original left maximum remains the comparison threshold throughout.
    assert_eq!(
        MD::parse("A10")?
            .checked_add(&MD::parse("B20 C9")?)?
            .number_of_max_aa(),
        20
    );
    Ok(())
}

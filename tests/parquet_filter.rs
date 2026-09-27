// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `ParquetFilter` and `ParquetFilterBuilder`.
//!
//! Every section of the pinned class test
//! `src/tests/class_tests/openms/source/ParquetFilter_test.cpp` (core
//! `bc9cc12`) is transcribed here with its literals (evidence tier 3); each
//! test names the section it ports. The module needs no Arrow type, so this
//! binary runs in the default build.

use openms::format::parquet_filter::{
    AND, ColumnType, Condition, FilterExpression, OR, ParquetFilter, ParquetFilterBuilder,
};

/// Section `ParquetFilter()`.
#[test]
fn a_fresh_filter_has_no_conditions() {
    let f = ParquetFilter::new();
    assert!(f.is_empty());
    assert!(f.expression().is_empty());
    assert_eq!(f.expression().conditions.len(), 0);
    assert_eq!(f.expression().connectors.len(), 0);
    assert_eq!(ParquetFilter::default(), f);
}

/// Section `eq(const std::string&, Int64)`.
#[test]
fn eq_int_records_an_integer_equality() {
    let mut f = ParquetFilter::new();
    f.eq_int("PRECURSOR_ID", 1318);
    assert!(!f.is_empty());
    let c = &f.expression().conditions[0];
    assert_eq!(c.column, "PRECURSOR_ID");
    assert_eq!(c.op, "=");
    assert_eq!(c.column_type, ColumnType::Int);
    assert_eq!(c.values.len(), 1);
    assert_eq!(c.values[0], "1318");
}

/// Section `ne(const std::string&, Int64)`.
#[test]
fn ne_int_keeps_the_sign() {
    let mut f = ParquetFilter::new();
    f.ne_int("ID", -5);
    assert_eq!(f.expression().conditions[0].op, "!=");
    assert_eq!(f.expression().conditions[0].values[0], "-5");
    assert_eq!(f.expression().conditions[0].column_type, ColumnType::Int);
}

/// Sections `lt`, `le`, `gt`, `ge` for `Int64`.
#[test]
fn integer_comparison_operators() {
    for (op, expected) in [("lt", "<"), ("le", "<="), ("gt", ">"), ("ge", ">=")] {
        let mut f = ParquetFilter::new();
        match op {
            "lt" => f.lt_int("ID", 10),
            "le" => f.le_int("ID", 10),
            "gt" => f.gt_int("ID", 10),
            _ => f.ge_int("ID", 10),
        };
        assert_eq!(f.expression().conditions[0].op, expected, "{op}");
        assert_eq!(f.expression().conditions[0].values, ["10"], "{op}");
        assert_eq!(f.expression().conditions[0].column_type, ColumnType::Int, "{op}");
    }
}

/// Section `in(const std::string&, const std::vector<Int64>&)`.
#[test]
fn in_int_records_every_value() {
    let mut f = ParquetFilter::new();
    f.in_int("ID", &[1, 2, 3]);
    let c = &f.expression().conditions[0];
    assert_eq!(c.op, "IN");
    assert_eq!(c.column_type, ColumnType::Int);
    assert_eq!(c.values.len(), 3);
    assert_eq!(c.values[0], "1");
    assert_eq!(c.values[2], "3");
}

/// Section `eq(const std::string&, const std::string&)`.
#[test]
fn eq_str_records_a_string_equality() {
    let mut f = ParquetFilter::new();
    f.eq_str("ANNOTATION", "y3^1");
    let c = &f.expression().conditions[0];
    assert_eq!(c.op, "=");
    assert_eq!(c.column_type, ColumnType::String);
    assert_eq!(c.values[0], "y3^1");
}

/// Sections `ne`, `lt`, `le`, `gt`, `ge` for `std::string`.
#[test]
fn string_comparison_operators() {
    let mut f = ParquetFilter::new();
    f.ne_str("ANNOTATION", "b2");
    assert_eq!(f.expression().conditions[0].op, "!=");
    assert_eq!(f.expression().conditions[0].column_type, ColumnType::String);
    for (op, expected) in [("lt", "<"), ("le", "<="), ("gt", ">"), ("ge", ">=")] {
        let mut f = ParquetFilter::new();
        match op {
            "lt" => f.lt_str("NAME", "m"),
            "le" => f.le_str("NAME", "m"),
            "gt" => f.gt_str("NAME", "m"),
            _ => f.ge_str("NAME", "m"),
        };
        assert_eq!(f.expression().conditions[0].op, expected, "{op}");
        assert_eq!(f.expression().conditions[0].column_type, ColumnType::String, "{op}");
        assert_eq!(f.expression().conditions[0].values, ["m"], "{op}");
    }
}

/// Section `in(const std::string&, const std::vector<std::string>&)`.
#[test]
fn in_str_records_every_value() {
    let mut f = ParquetFilter::new();
    f.in_str("NAME", &["a", "b"]);
    let c = &f.expression().conditions[0];
    assert_eq!(c.op, "IN");
    assert_eq!(c.column_type, ColumnType::String);
    assert_eq!(c.values.len(), 2);
    assert_eq!(c.values[1], "b");
}

/// Section `andNext()`: AND is the default connector and `and_next` is
/// equivalent to it.
#[test]
fn and_is_the_default_connector() {
    let mut f = ParquetFilter::new();
    f.eq_int("A", 1).eq_int("B", 2);
    assert_eq!(f.expression().conditions.len(), 2);
    assert_eq!(f.expression().connectors.len(), 1);
    assert_eq!(f.expression().connectors[0], "AND");
    let mut g = ParquetFilter::new();
    g.eq_int("A", 1).and_next().eq_int("B", 2);
    assert_eq!(g.expression().connectors[0], "AND");
    assert_eq!(AND, "AND");
}

/// Section `orNext()`: OR applies to the next condition only and then
/// resets to AND.
#[test]
fn or_next_applies_to_one_connector() {
    let mut f = ParquetFilter::new();
    f.eq_int("A", 1).or_next().eq_int("B", 2).eq_int("C", 3);
    assert_eq!(f.expression().conditions.len(), 3);
    assert_eq!(f.expression().connectors.len(), 2);
    assert_eq!(f.expression().connectors[0], "OR");
    assert_eq!(f.expression().connectors[1], "AND");
    assert_eq!(OR, "OR");
}

/// Derived from `addCondition_` (tier 4): a connector chosen before the first
/// condition is consumed by nothing, and the first condition resets it.
#[test]
fn a_connector_before_the_first_condition_is_dropped() {
    let mut f = ParquetFilter::new();
    f.or_next().eq_int("A", 1).eq_int("B", 2);
    assert_eq!(f.expression().connectors, ["AND"]);
}

/// Section `expression()`.
#[test]
fn expression_exposes_conditions_and_connectors() {
    let mut f = ParquetFilter::new();
    f.eq_int("A", 1).or_next().eq_str("B", "x");
    let e = f.expression();
    assert_eq!(e.conditions.len(), 2);
    assert_eq!(e.conditions[0].column, "A");
    assert_eq!(e.conditions[1].column, "B");
    assert_eq!(e.connectors[0], "OR");
}

/// Section `empty()`.
#[test]
fn empty_follows_the_conditions() {
    let mut f = ParquetFilter::new();
    assert!(f.is_empty());
    f.eq_int("A", 1);
    assert!(!f.is_empty());
}

/// Section `merge(lhs, rhs, connector = "AND")`.
#[test]
fn merge_joins_two_expressions() {
    let mut a = ParquetFilter::new();
    a.eq_int("A", 1).eq_int("B", 2);
    let mut b = ParquetFilter::new();
    b.eq_int("C", 3);
    let m = ParquetFilter::merge(a.expression(), b.expression(), "OR");
    assert_eq!(m.conditions.len(), 3);
    assert_eq!(m.connectors.len(), 2);
    assert_eq!(m.connectors[0], "AND");
    assert_eq!(m.connectors[1], "OR");
    assert_eq!(m.conditions[2].column, "C");
    let empty = FilterExpression::default();
    let m2 = ParquetFilter::merge(&empty, b.expression(), AND);
    assert_eq!(m2.conditions.len(), 1);
    assert_eq!(m2.connectors.len(), 0);
    let m3 = ParquetFilter::merge(a.expression(), &empty, AND);
    assert_eq!(m3.conditions.len(), 2);
    assert_eq!(m3.connectors.len(), 1);
}

/// Derived from `merge` (tier 4): an empty side returns the other side
/// unchanged, connectors included, even when the empty side carries
/// connectors of its own.
#[test]
fn merge_with_an_empty_side_ignores_its_connectors() {
    let mut b = ParquetFilter::new();
    b.eq_int("C", 3);
    let odd = FilterExpression {
        conditions: Vec::new(),
        connectors: vec!["OR".to_owned()],
    };
    assert_eq!(ParquetFilter::merge(&odd, b.expression(), "AND"), *b.expression());
    assert_eq!(ParquetFilter::merge(b.expression(), &odd, "AND"), *b.expression());
    assert!(odd.is_empty());
}

/// Builder section `eq(const std::string&, Int64)`.
#[test]
fn builder_forwards_to_the_filter() {
    let mut fb = ParquetFilterBuilder::new();
    fb.eq_int("PRECURSOR_ID", 1318).and_next().eq_str("ANNOTATION", "y3^1");
    assert!(!fb.is_empty());
    let pf = fb.filter();
    assert_eq!(pf.expression().conditions.len(), 2);
    assert_eq!(pf.expression().conditions[0].column, "PRECURSOR_ID");
    assert_eq!(pf.expression().conditions[0].column_type, ColumnType::Int);
    assert_eq!(pf.expression().conditions[1].column, "ANNOTATION");
    assert_eq!(pf.expression().conditions[1].column_type, ColumnType::String);
    assert_eq!(pf.expression().connectors[0], "AND");
}

/// Builder sections `ne`, `lt`, `le`, `gt`, `ge`, `in` for both value types.
#[test]
fn builder_operators() {
    let op = |build: &dyn Fn(&mut ParquetFilterBuilder)| {
        let mut fb = ParquetFilterBuilder::new();
        build(&mut fb);
        fb.filter().expression().conditions[0].clone()
    };
    assert_eq!(op(&|fb| { fb.ne_int("ID", 5); }).op, "!=");
    assert_eq!(op(&|fb| { fb.lt_int("ID", 5); }).op, "<");
    assert_eq!(op(&|fb| { fb.le_int("ID", 5); }).op, "<=");
    assert_eq!(op(&|fb| { fb.gt_int("ID", 5); }).op, ">");
    assert_eq!(op(&|fb| { fb.ge_int("ID", 5); }).op, ">=");
    let c = op(&|fb| {
        fb.in_int("ID", &[7, 8]);
    });
    assert_eq!(c.op, "IN");
    assert_eq!(c.values.len(), 2);
    let c = op(&|fb| {
        fb.eq_str("NAME", "x");
    });
    assert_eq!(c.column_type, ColumnType::String);
    assert_eq!(c.op, "=");
    assert_eq!(op(&|fb| { fb.ne_str("NAME", "x"); }).op, "!=");
    assert_eq!(op(&|fb| { fb.lt_str("NAME", "x"); }).op, "<");
    assert_eq!(op(&|fb| { fb.le_str("NAME", "x"); }).op, "<=");
    assert_eq!(op(&|fb| { fb.gt_str("NAME", "x"); }).op, ">");
    assert_eq!(op(&|fb| { fb.ge_str("NAME", "x"); }).op, ">=");
    let c = op(&|fb| {
        fb.in_str("NAME", &["a", "b", "c"]);
    });
    assert_eq!(c.op, "IN");
    assert_eq!(c.column_type, ColumnType::String);
    assert_eq!(c.values.len(), 3);
}

/// Builder sections `andNext()`, `orNext()`, `filter()` and `empty()`.
#[test]
fn builder_connectors_filter_and_empty() {
    let mut fb = ParquetFilterBuilder::new();
    fb.eq_int("A", 1).and_next().eq_int("B", 2);
    assert_eq!(fb.filter().expression().connectors[0], "AND");
    let mut fb = ParquetFilterBuilder::new();
    fb.eq_int("A", 1).or_next().eq_int("B", 2);
    assert_eq!(fb.filter().expression().connectors[0], "OR");
    let mut fb = ParquetFilterBuilder::new();
    assert!(fb.is_empty());
    fb.eq_int("A", 1);
    assert_eq!(fb.filter().expression().conditions.len(), 1);
    assert!(!fb.is_empty());
}

/// `Condition`'s default type is INT, as the source member initialiser.
#[test]
fn a_default_condition_is_an_integer_condition() {
    assert_eq!(Condition::default().column_type, ColumnType::Int);
}

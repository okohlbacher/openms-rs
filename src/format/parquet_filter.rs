// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Typed filter expressions for Parquet-backed readers.
//!
//! Ports `FORMAT/ParquetFilter.h` and `ParquetFilter.cpp`: the column type
//! tag, one condition, the flat condition/connector expression, the
//! [`ParquetFilter`] builder and its [`ParquetFilterBuilder`] wrapper. See
//! `docs/ARROW_PARQUET_SUPPORT.md` for the full API mapping.
//!
//! The header describes an expression that a reader pushes down to an Arrow
//! Dataset scanner or evaluates in memory. Neither happens here: the source
//! class only *records* conditions, and so does this port. Evaluation belongs
//! to the readers that consume an expression (`XICParquetFile` and its
//! siblings), which are separate headers.
//!
//! The module needs no Arrow type, so unlike the rest of the Arrow/Parquet
//! family it is not behind the `arrow` feature.
//!
//! Values are stored as text, as in the source: an integer condition keeps
//! the decimal spelling `StringUtils::toStr(Int64)` produces, which for an
//! `i64` is exactly Rust's `Display`.
//!
//! ```
//! use openms::format::parquet_filter::{ColumnType, ParquetFilter};
//!
//! let mut filter = ParquetFilter::new();
//! filter.eq_int("PRECURSOR_ID", 1318).and_next().eq_str("ANNOTATION", "y3^1");
//! let expression = filter.expression();
//! assert_eq!(expression.conditions.len(), 2);
//! assert_eq!(expression.connectors, ["AND"]);
//! assert_eq!(expression.conditions[1].column_type, ColumnType::String);
//! ```

/// The connector [`ParquetFilter::and_next`] selects, and the default
/// between two consecutive conditions.
pub const AND: &str = "AND";
/// The connector [`ParquetFilter::or_next`] selects.
pub const OR: &str = "OR";

/// Column type for typed Parquet filters.
///
/// Ports `OpenMS::ColumnType`. The source enumerators are `INT` and `STRING`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ColumnType {
    /// An integer column; the source `ColumnType::INT` and a condition's
    /// default.
    #[default]
    Int,
    /// A string column; the source `ColumnType::STRING`.
    String,
}

/// A single filter condition: column, operator and values.
///
/// Ports `OpenMS::Condition`. The operator is one of `=`, `!=`, `<`, `<=`,
/// `>`, `>=` or `IN` when a [`ParquetFilter`] built it; the fields are public,
/// as in the source, so a caller can store anything.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct Condition {
    /// The column name, as the source `column`.
    pub column: String,
    /// The comparison operator, as the source `op`.
    pub op: String,
    /// The compared values as text, as the source `values`: one for a
    /// comparison, any number (including none) for `IN`.
    pub values: Vec<String>,
    /// How the values are to be read, as the source `type` (default
    /// [`ColumnType::Int`]). Renamed because `type` is a Rust keyword.
    pub column_type: ColumnType,
}

/// A left-to-right chain of conditions joined by connectors.
///
/// Ports `OpenMS::FilterExpression`. A well-formed expression has one
/// connector fewer than it has conditions; `connectors[i]` joins
/// `conditions[i]` and `conditions[i + 1]`. The fields are public, as in the
/// source, and nothing checks that shape.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct FilterExpression {
    /// The conditions, in the order they were added.
    pub conditions: Vec<Condition>,
    /// The connectors, `"AND"` or `"OR"` when a [`ParquetFilter`] built them.
    pub connectors: Vec<String>,
}

impl FilterExpression {
    /// Whether the expression has no conditions (the source `empty()`).
    ///
    /// Connectors are not consulted, exactly as in the source: an expression
    /// with connectors but no conditions is empty.
    pub fn is_empty(&self) -> bool {
        self.conditions.is_empty()
    }
}

/// Typed filter builder for Parquet-backed datasets.
///
/// Ports `OpenMS::ParquetFilter`. Conditions are combined left to right.
/// [`and_next`](Self::and_next) and [`or_next`](Self::or_next) choose how the
/// *next* condition is connected to the previous one; after a condition is
/// added the choice resets to AND, so `or_next` affects exactly one
/// connector. A connector chosen before the first condition is consumed by
/// nothing and reset by it.
///
/// The source overloads each operator for `Int64` and `std::string` values;
/// the overload set becomes the `_int` and `_str` method pairs.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ParquetFilter {
    expr: FilterExpression,
    next_connector: String,
}

impl Default for ParquetFilter {
    fn default() -> Self {
        Self {
            expr: FilterExpression::default(),
            next_connector: AND.to_owned(),
        }
    }
}

impl ParquetFilter {
    /// An empty filter whose next connector is AND.
    pub fn new() -> Self {
        Self::default()
    }

    /// Combine the next condition with logical AND.
    pub fn and_next(&mut self) -> &mut Self {
        AND.clone_into(&mut self.next_connector);
        self
    }

    /// Combine the next condition with logical OR.
    pub fn or_next(&mut self) -> &mut Self {
        OR.clone_into(&mut self.next_connector);
        self
    }

    /// Add an equality condition (`=`) for an integer column.
    pub fn eq_int(&mut self, column: &str, value: i64) -> &mut Self {
        self.add_condition(column, "=", ColumnType::Int, vec![value.to_string()])
    }

    /// Add an inequality condition (`!=`) for an integer column.
    pub fn ne_int(&mut self, column: &str, value: i64) -> &mut Self {
        self.add_condition(column, "!=", ColumnType::Int, vec![value.to_string()])
    }

    /// Add a less-than condition (`<`) for an integer column.
    pub fn lt_int(&mut self, column: &str, value: i64) -> &mut Self {
        self.add_condition(column, "<", ColumnType::Int, vec![value.to_string()])
    }

    /// Add a less-than-or-equal condition (`<=`) for an integer column.
    pub fn le_int(&mut self, column: &str, value: i64) -> &mut Self {
        self.add_condition(column, "<=", ColumnType::Int, vec![value.to_string()])
    }

    /// Add a greater-than condition (`>`) for an integer column.
    pub fn gt_int(&mut self, column: &str, value: i64) -> &mut Self {
        self.add_condition(column, ">", ColumnType::Int, vec![value.to_string()])
    }

    /// Add a greater-than-or-equal condition (`>=`) for an integer column.
    pub fn ge_int(&mut self, column: &str, value: i64) -> &mut Self {
        self.add_condition(column, ">=", ColumnType::Int, vec![value.to_string()])
    }

    /// Add an `IN` condition for an integer column.
    ///
    /// An empty `values` is recorded as an `IN` with no values, as in the
    /// source; what it selects is the evaluating reader's decision.
    pub fn in_int(&mut self, column: &str, values: &[i64]) -> &mut Self {
        let values = values.iter().map(i64::to_string).collect();
        self.add_condition(column, "IN", ColumnType::Int, values)
    }

    /// Add an equality condition (`=`) for a string column.
    pub fn eq_str(&mut self, column: &str, value: &str) -> &mut Self {
        self.add_condition(column, "=", ColumnType::String, vec![value.to_owned()])
    }

    /// Add an inequality condition (`!=`) for a string column.
    pub fn ne_str(&mut self, column: &str, value: &str) -> &mut Self {
        self.add_condition(column, "!=", ColumnType::String, vec![value.to_owned()])
    }

    /// Add a less-than condition (`<`) for a string column.
    pub fn lt_str(&mut self, column: &str, value: &str) -> &mut Self {
        self.add_condition(column, "<", ColumnType::String, vec![value.to_owned()])
    }

    /// Add a less-than-or-equal condition (`<=`) for a string column.
    pub fn le_str(&mut self, column: &str, value: &str) -> &mut Self {
        self.add_condition(column, "<=", ColumnType::String, vec![value.to_owned()])
    }

    /// Add a greater-than condition (`>`) for a string column.
    pub fn gt_str(&mut self, column: &str, value: &str) -> &mut Self {
        self.add_condition(column, ">", ColumnType::String, vec![value.to_owned()])
    }

    /// Add a greater-than-or-equal condition (`>=`) for a string column.
    pub fn ge_str(&mut self, column: &str, value: &str) -> &mut Self {
        self.add_condition(column, ">=", ColumnType::String, vec![value.to_owned()])
    }

    /// Add an `IN` condition for a string column.
    pub fn in_str<S: AsRef<str>>(&mut self, column: &str, values: &[S]) -> &mut Self {
        let values = values.iter().map(|v| v.as_ref().to_owned()).collect();
        self.add_condition(column, "IN", ColumnType::String, values)
    }

    /// The constructed filter expression.
    pub fn expression(&self) -> &FilterExpression {
        &self.expr
    }

    /// Whether the filter has no conditions (the source `empty()`).
    pub fn is_empty(&self) -> bool {
        self.expr.is_empty()
    }

    /// Merge two expressions with a connector.
    ///
    /// The result is `lhs`'s conditions and connectors, then `connector`,
    /// then `rhs`'s conditions and connectors. When either side has no
    /// conditions the other side is returned unchanged and `connector` is not
    /// used. The source's default connector is `"AND"` ([`AND`]); it is not
    /// validated, so any text is stored, exactly as in the source.
    pub fn merge(
        lhs: &FilterExpression,
        rhs: &FilterExpression,
        connector: &str,
    ) -> FilterExpression {
        if lhs.conditions.is_empty() {
            return rhs.clone();
        }
        if rhs.conditions.is_empty() {
            return lhs.clone();
        }
        let mut combined = FilterExpression {
            conditions: Vec::with_capacity(lhs.conditions.len() + rhs.conditions.len()),
            connectors: Vec::with_capacity(lhs.connectors.len() + rhs.connectors.len() + 1),
        };
        combined.conditions.extend(lhs.conditions.iter().cloned());
        combined.connectors.extend(lhs.connectors.iter().cloned());
        combined.connectors.push(connector.to_owned());
        combined.conditions.extend(rhs.conditions.iter().cloned());
        combined.connectors.extend(rhs.connectors.iter().cloned());
        combined
    }

    /// The source's private `addCondition_`: push the pending connector when
    /// a condition already exists, append the condition, reset to AND.
    fn add_condition(
        &mut self,
        column: &str,
        op: &str,
        column_type: ColumnType,
        values: Vec<String>,
    ) -> &mut Self {
        if !self.expr.conditions.is_empty() {
            self.expr.connectors.push(self.next_connector.clone());
        }
        self.expr.conditions.push(Condition {
            column: column.to_owned(),
            op: op.to_owned(),
            values,
            column_type,
        });
        AND.clone_into(&mut self.next_connector);
        self
    }
}

/// Fluent builder for [`ParquetFilter`] objects.
///
/// Ports `OpenMS::ParquetFilterBuilder`: a wrapper that forwards every
/// condition-building call to the [`ParquetFilter`] it owns.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct ParquetFilterBuilder {
    filter: ParquetFilter,
}

impl ParquetFilterBuilder {
    /// An empty builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Combine the next condition with logical AND.
    pub fn and_next(&mut self) -> &mut Self {
        self.filter.and_next();
        self
    }

    /// Combine the next condition with logical OR.
    pub fn or_next(&mut self) -> &mut Self {
        self.filter.or_next();
        self
    }

    /// Add an equality condition (`=`) for an integer column.
    pub fn eq_int(&mut self, column: &str, value: i64) -> &mut Self {
        self.filter.eq_int(column, value);
        self
    }

    /// Add an inequality condition (`!=`) for an integer column.
    pub fn ne_int(&mut self, column: &str, value: i64) -> &mut Self {
        self.filter.ne_int(column, value);
        self
    }

    /// Add a less-than condition (`<`) for an integer column.
    pub fn lt_int(&mut self, column: &str, value: i64) -> &mut Self {
        self.filter.lt_int(column, value);
        self
    }

    /// Add a less-than-or-equal condition (`<=`) for an integer column.
    pub fn le_int(&mut self, column: &str, value: i64) -> &mut Self {
        self.filter.le_int(column, value);
        self
    }

    /// Add a greater-than condition (`>`) for an integer column.
    pub fn gt_int(&mut self, column: &str, value: i64) -> &mut Self {
        self.filter.gt_int(column, value);
        self
    }

    /// Add a greater-than-or-equal condition (`>=`) for an integer column.
    pub fn ge_int(&mut self, column: &str, value: i64) -> &mut Self {
        self.filter.ge_int(column, value);
        self
    }

    /// Add an `IN` condition for an integer column.
    pub fn in_int(&mut self, column: &str, values: &[i64]) -> &mut Self {
        self.filter.in_int(column, values);
        self
    }

    /// Add an equality condition (`=`) for a string column.
    pub fn eq_str(&mut self, column: &str, value: &str) -> &mut Self {
        self.filter.eq_str(column, value);
        self
    }

    /// Add an inequality condition (`!=`) for a string column.
    pub fn ne_str(&mut self, column: &str, value: &str) -> &mut Self {
        self.filter.ne_str(column, value);
        self
    }

    /// Add a less-than condition (`<`) for a string column.
    pub fn lt_str(&mut self, column: &str, value: &str) -> &mut Self {
        self.filter.lt_str(column, value);
        self
    }

    /// Add a less-than-or-equal condition (`<=`) for a string column.
    pub fn le_str(&mut self, column: &str, value: &str) -> &mut Self {
        self.filter.le_str(column, value);
        self
    }

    /// Add a greater-than condition (`>`) for a string column.
    pub fn gt_str(&mut self, column: &str, value: &str) -> &mut Self {
        self.filter.gt_str(column, value);
        self
    }

    /// Add a greater-than-or-equal condition (`>=`) for a string column.
    pub fn ge_str(&mut self, column: &str, value: &str) -> &mut Self {
        self.filter.ge_str(column, value);
        self
    }

    /// Add an `IN` condition for a string column.
    pub fn in_str<S: AsRef<str>>(&mut self, column: &str, values: &[S]) -> &mut Self {
        self.filter.in_str(column, values);
        self
    }

    /// The built filter, to pass to a Parquet-backed reader.
    pub fn filter(&self) -> &ParquetFilter {
        &self.filter
    }

    /// Whether the filter has no conditions.
    pub fn is_empty(&self) -> bool {
        self.filter.is_empty()
    }
}

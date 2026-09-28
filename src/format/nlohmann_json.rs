// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! JSON text as the `nlohmann::json` library the source links writes it.
//!
//! Two writers need the library's exact text: the OMS JSON export
//! (`OMSFileLoad::exportToJSON`, `dump(4)`) and the mzQC writer
//! (`MzQCFile::store`, an `ordered_json` with `dump(2)`). This module holds
//! the number and string formatting both share, and a small ordered value
//! type with the library's pretty printer. `serde_json` is not used for this:
//! its number text differs (`1e-5` against `1e-05`, `100` against `100.0` for
//! an integral double) and its maps are sorted unless a crate-wide feature is
//! switched on.

/// A number as `nlohmann::json` prints a double: the shortest digits that
/// round-trip, in fixed notation for decimal exponents in `(-4, 15]` (with
/// `.0` appended to an integral value) and in scientific notation otherwise,
/// the exponent signed and at least two digits long. Non-finite values print
/// `null`, as the library does.
pub(crate) fn number(value: f64) -> String {
    if !value.is_finite() {
        return "null".into();
    }
    if value == 0.0 {
        return if value.is_sign_negative() {
            "-0.0"
        } else {
            "0.0"
        }
        .into();
    }
    // Rust's `{:e}` is the shortest round-trip representation.
    let scientific = format!("{:e}", value.abs());
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    let k = digits.len() as i32;
    // `n` is the position of the decimal point relative to the digit string.
    let n = exponent + 1;
    let sign = if value.is_sign_negative() { "-" } else { "" };
    let body = if k <= n && n <= 15 {
        format!("{digits}{}.0", "0".repeat((n - k) as usize))
    } else if 0 < n && n <= 15 {
        format!("{}.{}", &digits[..n as usize], &digits[n as usize..])
    } else if -4 < n && n <= 0 {
        format!("0.{}{digits}", "0".repeat((-n) as usize))
    } else {
        let exponent = n - 1;
        let exponent_text = format!(
            "{}{:02}",
            if exponent < 0 { '-' } else { '+' },
            exponent.abs()
        );
        if k == 1 {
            format!("{digits}e{exponent_text}")
        } else {
            format!("{}.{}e{exponent_text}", &digits[..1], &digits[1..])
        }
    };
    format!("{sign}{body}")
}

/// A JSON string literal as `nlohmann::json::dump` escapes it: `"` and `\`,
/// the short escapes for backspace, form feed, newline, carriage return and
/// tab, `\u00xx` (lower-case hex) for the other control characters, and every
/// other character verbatim (the library's default `ensure_ascii = false`).
pub(crate) fn string_literal(text: &str) -> String {
    let mut escaped = String::with_capacity(text.len() + 2);
    escaped.push('"');
    for c in text.chars() {
        match c {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\u{08}' => escaped.push_str("\\b"),
            '\u{0c}' => escaped.push_str("\\f"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            c if (c as u32) < 0x20 => escaped.push_str(&format!("\\u{:04x}", c as u32)),
            c => escaped.push(c),
        }
    }
    escaped.push('"');
    escaped
}

/// An ordered JSON value, `nlohmann::ordered_json`: object members keep
/// their insertion order.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Json {
    /// `null`, also what a default-constructed library value dumps as.
    Null,
    /// An unsigned integer, the library's `number_unsigned`.
    Unsigned(u64),
    /// A double, the library's `number_float`; a `float` is widened first.
    Float(f64),
    /// A string.
    String(String),
    /// An array.
    Array(Vec<Json>),
    /// An object, members in insertion order.
    Object(Vec<(String, Json)>),
}

impl Json {
    /// A string value.
    pub(crate) fn string(text: impl Into<String>) -> Self {
        Self::String(text.into())
    }
    /// An object from `(name, value)` members in order.
    pub(crate) fn object<const N: usize>(members: [(&str, Json); N]) -> Self {
        Self::Object(
            members
                .into_iter()
                .map(|(name, value)| (name.to_owned(), value))
                .collect(),
        )
    }
    /// The text `dump(indent)` writes: pretty-printed with `indent` spaces
    /// per level, `": "` between a member's name and value, one element or
    /// member per line, `[]`/`{}` for an empty array/object, and no final
    /// newline.
    pub(crate) fn dump(&self, indent: usize) -> String {
        let mut out = String::new();
        self.write(&mut out, indent, 0);
        out
    }
    fn write(&self, out: &mut String, indent: usize, level: usize) {
        match self {
            Self::Null => out.push_str("null"),
            Self::Unsigned(value) => out.push_str(&value.to_string()),
            Self::Float(value) => out.push_str(&number(*value)),
            Self::String(text) => out.push_str(&string_literal(text)),
            Self::Array(values) => {
                if values.is_empty() {
                    out.push_str("[]");
                    return;
                }
                out.push_str("[\n");
                for (position, value) in values.iter().enumerate() {
                    if position > 0 {
                        out.push_str(",\n");
                    }
                    push_indent(out, indent * (level + 1));
                    value.write(out, indent, level + 1);
                }
                out.push('\n');
                push_indent(out, indent * level);
                out.push(']');
            }
            Self::Object(members) => {
                if members.is_empty() {
                    out.push_str("{}");
                    return;
                }
                out.push_str("{\n");
                for (position, (name, value)) in members.iter().enumerate() {
                    if position > 0 {
                        out.push_str(",\n");
                    }
                    push_indent(out, indent * (level + 1));
                    out.push_str(&string_literal(name));
                    out.push_str(": ");
                    value.write(out, indent, level + 1);
                }
                out.push('\n');
                push_indent(out, indent * level);
                out.push('}');
            }
        }
    }
}

fn push_indent(out: &mut String, width: usize) {
    out.extend(std::iter::repeat_n(' ', width));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_follow_the_nlohmann_format() {
        assert_eq!(number(0.3), "0.3");
        assert_eq!(number(1.0), "1.0");
        assert_eq!(number(-2.5), "-2.5");
        assert_eq!(number(1234.5678), "1234.5678");
        assert_eq!(number(0.0001), "0.0001");
        assert_eq!(number(0.00001), "1e-05");
        assert_eq!(number(1.5e-7), "1.5e-07");
        assert_eq!(number(1e14), "100000000000000.0");
        assert_eq!(number(1e15), "1e+15");
        assert_eq!(number(1.25e16), "1.25e+16");
        assert_eq!(number(123_456_789_012_345_680.0), "1.2345678901234568e+17");
        assert_eq!(number(f64::NAN), "null");
        assert_eq!(number(0.0), "0.0");
        // A float widened to double, as the mzQC writer stores them.
        assert_eq!(number(f64::from(91.666_664_f32)), "91.66666412353516");
    }

    #[test]
    fn strings_escape_like_nlohmann() {
        assert_eq!(
            string_literal("a\"b\\c\nd\u{01}é/"),
            "\"a\\\"b\\\\c\\nd\\u0001é/\""
        );
    }

    #[test]
    fn dump_pretty_prints_like_nlohmann() {
        let value = Json::object([
            ("a", Json::Array(vec![Json::Unsigned(1), Json::Float(-1.0)])),
            ("b", Json::Array(Vec::new())),
            ("c", Json::Object(Vec::new())),
            ("d", Json::object([("e", Json::Null)])),
        ]);
        assert_eq!(
            value.dump(2),
            "{\n  \"a\": [\n    1,\n    -1.0\n  ],\n  \"b\": [],\n  \"c\": {},\n  \"d\": {\n    \"e\": null\n  }\n}"
        );
    }
}

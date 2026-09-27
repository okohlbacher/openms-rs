// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The string helpers of `PEFFFile.cpp`'s anonymous namespace
//! (`PEFFFile.cpp:29-503`), and the `StringUtils` conversions PEFF relies on.

use super::PEFFCustomKeyDef;
use crate::{Error, Result};

/// Source `StringUtils::trim` and the unquoted `StringUtils::split`, shared
/// with the other text readers.
pub(super) use crate::data_structures::string_utils::{split, trim};

/// Source `StringUtils::toInt32` (`StringUtils.cpp:136-166`) through the shared
/// [`crate::data_structures::string_utils::to_int32`]; the source's
/// `Exception::ConversionError` is reported as [`Error::Parse`] with its message.
pub(super) fn to_int32(text: &str) -> Result<i32> {
    crate::data_structures::string_utils::to_int32(text)
        .map_err(|message| Error::Parse { line: 0, message })
}

/// A `Size` field assigned from `StringUtils::toInt32`: a negative value is
/// converted as C++ converts it, modulo 2^64.
pub(super) fn to_size(text: &str) -> Result<usize> {
    Ok(to_int32(text)? as i64 as u64 as usize)
}

/// Source `peffUnescape`: `\\`, `\|`, `\(` and `\)` become the character
/// they escape; any other backslash stays.
pub(super) fn peff_unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '\\' {
            if let Some(&next) = chars.peek() {
                if matches!(next, '\\' | '|' | '(' | ')') {
                    out.push(next);
                    chars.next();
                    continue;
                }
            }
        }
        out.push(c);
    }
    out
}

/// Source `escape_peff`: prefix `\`, `|`, `(` and `)` with a backslash.
pub(super) fn escape_peff(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        if matches!(c, '\\' | '|' | '(' | ')') {
            out.push('\\');
        }
        out.push(c);
    }
    out
}

/// Source `splitByPipeEscapeAware`: split on unescaped `|`, un-escaping
/// `\\`, `\|`, `\(` and `\)` in every field. Always at least one field.
pub(super) fn split_by_pipe_escape_aware(text: &str) -> Vec<String> {
    let mut result = Vec::new();
    let mut current = String::new();
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '\\' {
            if let Some(&next) = chars.peek() {
                if matches!(next, '|' | '\\' | '(' | ')') {
                    current.push(next);
                    chars.next();
                    continue;
                }
            }
        }
        if c == '|' {
            result.push(std::mem::take(&mut current));
        } else {
            current.push(c);
        }
    }
    result.push(current);
    result
}

/// Source `extractAnnotationId`: an all-digit prefix before the first `:` is
/// the annotation identifier, and the rest the value; otherwise there is no
/// identifier ([`super::ANNOTATION_ID_NOT_SET`]) and the value is the field.
///
/// # Errors
///
/// [`Error::Parse`] when the digits exceed `int32`, the source's
/// `StringUtils::toInt32` `ConversionError`.
pub(super) fn extract_annotation_id(field: &str) -> Result<(u32, &str)> {
    if let Some((prefix, rest)) = field.split_once(':') {
        if !prefix.is_empty() && prefix.bytes().all(|b| b.is_ascii_digit()) {
            // Digits only, so the int32 is non-negative and fits u32.
            return Ok((to_int32(prefix)? as u32, rest));
        }
    }
    Ok((super::ANNOTATION_ID_NOT_SET, field))
}

/// Source `isXsdInteger_`.
pub(super) fn is_xsd_integer(text: &str) -> bool {
    let digits = text.strip_prefix(['+', '-']).unwrap_or(text);
    !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit())
}

/// Source `isUnsignedInteger_`.
pub(super) fn is_unsigned_integer(text: &str) -> bool {
    !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit())
}

/// Source `isXsdDecimal_`: an optional sign, digits with at most one `.`,
/// and at least one digit.
fn is_xsd_decimal(text: &str) -> bool {
    let body = text.strip_prefix(['+', '-']).unwrap_or(text);
    if body.is_empty() {
        return false;
    }
    let mut saw_digit = false;
    let mut saw_dot = false;
    for b in body.bytes() {
        if b.is_ascii_digit() {
            saw_digit = true;
        } else if b == b'.' && !saw_dot {
            saw_dot = true;
        } else {
            return false;
        }
    }
    saw_digit
}

/// Source `isXsdBoolean_`.
fn is_xsd_boolean(text: &str) -> bool {
    matches!(
        text.to_ascii_lowercase().as_str(),
        "true" | "false" | "1" | "0"
    )
}

fn digit(bytes: &[u8], index: usize) -> bool {
    bytes.get(index).is_some_and(u8::is_ascii_digit)
}

/// The optional `Z` or `+hh:mm` / `-hh:mm` zone at `at`, reaching the end.
fn timezone(bytes: &[u8], at: usize) -> bool {
    if at == bytes.len() {
        return true;
    }
    if bytes[at] == b'Z' {
        return at + 1 == bytes.len();
    }
    if bytes[at] != b'+' && bytes[at] != b'-' {
        return false;
    }
    at + 6 == bytes.len()
        && digit(bytes, at + 1)
        && digit(bytes, at + 2)
        && bytes[at + 3] == b':'
        && digit(bytes, at + 4)
        && digit(bytes, at + 5)
}

/// Source `isXsdDate_`: `[-]YYYY-MM-DD` with an optional zone; ranges are
/// deliberately not checked.
fn is_xsd_date(text: &str) -> bool {
    let bytes = text.as_bytes();
    if bytes.len() < 10 {
        return false;
    }
    let i = usize::from(bytes[0] == b'-');
    if i + 10 > bytes.len() {
        return false;
    }
    (0..4).all(|k| digit(bytes, i + k))
        && bytes[i + 4] == b'-'
        && digit(bytes, i + 5)
        && digit(bytes, i + 6)
        && bytes[i + 7] == b'-'
        && digit(bytes, i + 8)
        && digit(bytes, i + 9)
        && timezone(bytes, i + 10)
}

/// Source `isXsdTime_`: `hh:mm:ss`, optional fraction and zone; ranges are
/// deliberately not checked.
fn is_xsd_time(text: &str) -> bool {
    let bytes = text.as_bytes();
    if bytes.len() < 8 {
        return false;
    }
    if !(digit(bytes, 0)
        && digit(bytes, 1)
        && bytes[2] == b':'
        && digit(bytes, 3)
        && digit(bytes, 4)
        && bytes[5] == b':'
        && digit(bytes, 6)
        && digit(bytes, 7))
    {
        return false;
    }
    let mut position = 8;
    if bytes.get(position) == Some(&b'.') {
        position += 1;
        let start = position;
        while digit(bytes, position) {
            position += 1;
        }
        if position == start {
            return false;
        }
    }
    timezone(bytes, position)
}

/// Source `validateCustomKeyFieldType_`: whether `value` fits the declared
/// field type. Unknown types and the empty type always fit.
pub(super) fn validate_custom_key_field_type(field_type: &str, value: &str) -> bool {
    let declared = trim(field_type);
    let lower = declared.to_ascii_lowercase();
    match lower.as_str() {
        "" | "string" => true,
        "integer" => is_xsd_integer(value),
        "decimal" => is_xsd_decimal(value),
        "boolean" => is_xsd_boolean(value),
        "date" => is_xsd_date(value),
        "time" => is_xsd_time(value),
        _ if lower.starts_with("enumeration(") && declared.ends_with(')') => {
            let (Some(left), Some(right)) = (declared.find('('), declared.rfind(')')) else {
                return false;
            };
            if right <= left + 1 {
                return false;
            }
            split(&declared[left + 1..right], '|')
                .into_iter()
                .any(|option| trim(option) == value)
        }
        _ => true,
    }
}

/// Source `splitCustomKeyDefSegments_`: split on `|` outside double quotes
/// and outside parentheses, keeping escapes verbatim.
fn split_custom_key_def_segments(text: &str) -> Vec<String> {
    let mut result = Vec::new();
    let mut current = String::new();
    let mut in_quotes = false;
    let mut depth = 0usize;
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '\\' {
            if let Some(next) = chars.next() {
                current.push(c);
                current.push(next);
                continue;
            }
        }
        if c == '"' {
            in_quotes = !in_quotes;
            current.push(c);
            continue;
        }
        if !in_quotes {
            if c == '(' {
                depth += 1;
            } else if c == ')' && depth > 0 {
                depth -= 1;
            } else if c == '|' && depth == 0 {
                result.push(std::mem::take(&mut current));
                continue;
            }
        }
        current.push(c);
    }
    result.push(current);
    result
}

/// Source `unquoteCustomKeyDefValue_`: trim; a value wrapped in double quotes
/// loses them and has `\"` un-escaped.
fn unquote_custom_key_def_value(value: &str) -> String {
    let trimmed = trim(value);
    if trimmed.len() >= 2 && trimmed.starts_with('"') && trimmed.ends_with('"') {
        let inner = &trimmed[1..trimmed.len() - 1];
        return inner.replace("\\\"", "\"");
    }
    trimmed.to_owned()
}

/// Source `splitCommaSeparatedList_`: trimmed, non-empty items.
fn split_comma_separated_list(value: &str) -> Vec<String> {
    split(value, ',')
        .into_iter()
        .map(trim)
        .filter(|item| !item.is_empty())
        .map(str::to_owned)
        .collect()
}

/// Source `parseCustomKeyDef_`, with its `OPENMS_LOG_WARN` lines appended to
/// `warnings`.
pub(super) fn parse_custom_key_def(value: &str, warnings: &mut Vec<String>) -> PEFFCustomKeyDef {
    let mut definition = PEFFCustomKeyDef::default();
    let trimmed = trim(value);
    if trimmed.len() >= 2 && trimmed.starts_with('(') && trimmed.ends_with(')') {
        for segment in split_custom_key_def_segments(&trimmed[1..trimmed.len() - 1]) {
            let segment = trim(&segment);
            if segment.is_empty() {
                continue;
            }
            let Some((key, value)) = segment.split_once('=') else {
                warnings.push(format!(
                    "PEFF: CustomKeyDef segment is missing '=': '{segment}'."
                ));
                continue;
            };
            let (key, value) = (trim(key), trim(value));
            match key {
                "KeyName" => definition.key_name = value.to_owned(),
                "Description" => definition.description = unquote_custom_key_def_value(value),
                "ConceptCURIE" => definition.concept_curie = value.to_owned(),
                "RegExp" => definition.regexp = unquote_custom_key_def_value(value),
                "FieldNames" => definition.field_names = split_comma_separated_list(value),
                "FieldTypes" => definition.field_types = split_comma_separated_list(value),
                _ => warnings.push(format!("PEFF: Unknown CustomKeyDef attribute '{key}'.")),
            }
        }
    } else if let Some((name, description)) = trimmed.split_once(':') {
        // The historical KeyName:description form.
        definition.key_name = name.to_owned();
        definition.description = description.to_owned();
    } else {
        definition.key_name = trimmed.to_owned();
    }
    if definition.key_name.is_empty() {
        warnings.push("PEFF: CustomKeyDef is missing required KeyName.".to_owned());
    }
    definition
}

/// Source `escapeCustomKeyDefQuotedValue_`: `"` becomes `\"`.
pub(super) fn escape_custom_key_def_quoted_value(value: &str) -> String {
    value.replace('"', "\\\"")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn to_int32_follows_string_utils() {
        assert_eq!(to_int32(" 42 ").unwrap(), 42);
        assert_eq!(to_int32("+7").unwrap(), 7);
        assert_eq!(to_int32("-7").unwrap(), -7);
        assert_eq!(to_int32("+-7").unwrap(), -7);
        for bad in ["", " ", "x", "7x", "2147483648", "++7", "7 8"] {
            assert!(to_int32(bad).is_err(), "{bad:?}");
        }
        assert_eq!(to_size("-1").unwrap(), usize::MAX);
    }

    #[test]
    fn escapes_invert_each_other_on_reserved_characters() {
        let text = r"a\b|c(d)e";
        assert_eq!(escape_peff(text), r"a\\b\|c\(d\)e");
        assert_eq!(peff_unescape(&escape_peff(text)), text);
        assert_eq!(peff_unescape(r"\x\"), r"\x\");
        assert_eq!(
            split_by_pipe_escape_aware(r"1|MOD:1|a\|b|t"),
            ["1", "MOD:1", "a|b", "t"]
        );
    }

    #[test]
    fn annotation_ids_need_a_digit_prefix() {
        assert_eq!(extract_annotation_id("3:10").unwrap(), (3, "10"));
        assert_eq!(
            extract_annotation_id("x:10").unwrap(),
            (super::super::ANNOTATION_ID_NOT_SET, "x:10")
        );
        assert_eq!(
            extract_annotation_id(":10").unwrap(),
            (super::super::ANNOTATION_ID_NOT_SET, ":10")
        );
        assert!(extract_annotation_id("99999999999:1").is_err());
    }

    #[test]
    fn field_types_follow_the_source_validators() {
        assert!(validate_custom_key_field_type("integer", "-12"));
        assert!(!validate_custom_key_field_type("integer", "1.5"));
        assert!(validate_custom_key_field_type(" Decimal ", ".5"));
        assert!(!validate_custom_key_field_type("decimal", "."));
        assert!(validate_custom_key_field_type("boolean", "TRUE"));
        assert!(validate_custom_key_field_type("date", "2024-01-31+01:00"));
        assert!(!validate_custom_key_field_type("date", "2024-1-31"));
        assert!(validate_custom_key_field_type("time", "12:30:00.5Z"));
        assert!(!validate_custom_key_field_type("time", "12:30:00."));
        assert!(validate_custom_key_field_type("enumeration(a|b)", "b"));
        assert!(!validate_custom_key_field_type("enumeration(a|b)", "c"));
        assert!(validate_custom_key_field_type("unknown", "anything"));
    }
}

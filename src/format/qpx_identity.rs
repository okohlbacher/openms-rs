// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Opaque QPX row identities (`feature_id`, `psm_id`, `pg_id`).
//!
//! Ports the `QPXIdentity` namespace (`FORMAT/QPXIdentity.h` and
//! `QPXIdentity.cpp`). See `docs/QPX_SUPPORT.md` for the API mapping.
//!
//! Every QPX view carries a mandatory `int64` identity column that is its
//! primary key. The value is opaque -- it is not meant to be parsed or
//! reversed, only compared -- and is derived deterministically from a
//! footer-declared `identity_composite` of ordinary columns.
//!
//! **These functions reproduce `qpx.core.data.identity` byte for byte.** That
//! is a hard requirement, not a nicety: `qpxc` re-derives `feature_id` for
//! identified rows from the declared composite but never rewrites the
//! `psm.feature_id` that points at it, so an OpenMS collection whose ids were
//! invented rather than derived would turn into one with dangling
//! cross-references the moment it was converted (qpx's dataset validation
//! reports them as `dangling_feature_id`).
//!
//! The derivation is BLAKE2b with an 8-byte digest over the canonical JSON
//! encoding of the composite, read as a big-endian signed 64-bit integer.
//! Negative values are normal and carry no meaning. BLAKE2b-64 is **not**
//! BLAKE2b-512 truncated: the digest length is part of the parameter block and
//! therefore of the initial state. The source vendors an RFC 7693
//! implementation; this port uses the RustCrypto `blake2` crate with a
//! compile-time 8-byte output, which sets the same parameter block (see
//! `docs/THIRD_PARTY_CRATE_DECISIONS.md`).
//!
//! Identity is meaningful **within a file only**. Two QPX files must not be
//! joined on `feature_id` alone, and a test reference must not pin id values
//! across files: the feature composite contains `rt` and `observed_mz`, so one
//! ULP of platform drift changes the whole id rather than one digit of it.
//!
//! **Evidence.** Every expected value in `tests/qpx_identity.rs` is transcribed
//! from the upstream class test `QPXIdentity_test.cpp`, whose values were
//! produced by the reference implementation `qpx.core.data.identity`, not by
//! running the C++.

use blake2::digest::consts::U8;
use blake2::{Blake2b, Digest};
use std::collections::HashMap;
use std::fmt::Write as _;

/// The feature-to-PSM edge of one QPX collection, as `psm_id -> feature_id`.
///
/// Ports `QPXIdentity::FeatureLinks` (`std::unordered_map<Int64, Int64>`).
/// The consensus-map exporter collects it and the psm exporter fills
/// `psm.feature_id` from it. Keying on the *id* rather than on a pointer or a
/// row index is what keeps the two directions consistent: both views derive
/// the key from the same four persisted columns, so they agree by construction
/// rather than by both happening to walk the identifications in the same
/// order.
pub type FeatureLinks = HashMap<i64, i64>;

/// One component of an identity composite, in its persisted Arrow type.
///
/// Ports `QPXIdentity::Value`, a `std::variant` of the source's `Null`,
/// `std::string`, `Int64`, `Float32`, `std::vector<Int32>` and
/// `std::vector<std::string>`.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    /// A composite component that is JSON `null` (a null column value); the
    /// source `QPXIdentity::Null`.
    Null,
    /// A JSON string.
    String(String),
    /// A JSON integer.
    Int64(i64),
    /// A `float32` component, a JSON number; the source `QPXIdentity::Float32`.
    ///
    /// Distinct from a double so a caller cannot pass an unrounded value by
    /// accident: QPX hashes the value *as persisted*, so a column stored as
    /// `float32` must be narrowed to `f32` before it is encoded, whatever
    /// precision it was computed at.
    Float32(f32),
    /// A JSON array of integers (`scan`).
    Int32List(Vec<i32>),
    /// A JSON array of strings (`grouped_runs`, `pg_accessions`).
    StringList(Vec<String>),
}

/// Footer `identity_composite` for the feature view.
///
/// Ports `QPXIdentity::FEATURE_COMPOSITE`. These strings go into the Parquet
/// footer and `qpxc` reads them back to re-derive the ids, so they name the
/// columns in the order the values are hashed in.
pub const FEATURE_COMPOSITE: &str = "run_file_name,peptidoform,charge,rt,scan,observed_mz";
/// Footer `identity_composite` for the psm view; the source
/// `QPXIdentity::PSM_COMPOSITE`.
pub const PSM_COMPOSITE: &str = "run_file_name,scan,peptidoform,charge";
/// Footer `identity_composite` for the pg view; the source
/// `QPXIdentity::PG_COMPOSITE`.
pub const PG_COMPOSITE: &str = "pg_accessions,grouped_runs,label";

/// Render one `float32` the way Python's `repr` writes it inside JSON.
///
/// Ports `QPXIdentity::formatFloat`. QPX encodes composites with
/// `json.dumps`, which formats floats with `float.__repr__`: the shortest
/// decimal that round-trips, switching to exponential notation when the
/// decimal point would fall at position `<= -4` or `> 16`, and appending
/// `".0"` to a value that would otherwise look like an integer. The exponent
/// is always signed and at least two digits.
///
/// The value is widened to `f64` first, exactly as Arrow's `cast` +
/// `to_pylist` does on the reading side, so `0.1f32` renders as
/// `"0.10000000149011612"`. Other examples: `"1000000000000000.0"`,
/// `"1e+16"`, `"9.999e-05"`, `"-0.0"`. Non-finite values are written bare, as
/// `json.dumps` writes them: `"NaN"`, `"Infinity"`, `"-Infinity"`.
///
/// The source reads the shortest digits off `std::to_chars` in scientific
/// form; this reads them off Rust's `{:e}` formatting of the `f64`, which is
/// also the shortest round-tripping representation.
pub fn format_float(value: f32) -> String {
    let widened = f64::from(value);
    if widened.is_nan() {
        return "NaN".to_owned();
    }
    if widened.is_infinite() {
        return if widened < 0.0 {
            "-Infinity"
        } else {
            "Infinity"
        }
        .to_owned();
    }

    let mut out = String::new();
    if widened.is_sign_negative() {
        out.push('-'); // -0.0 keeps its sign, as repr() does
    }
    let magnitude = widened.abs();
    if magnitude == 0.0 {
        out.push_str("0.0");
        return out;
    }

    // Shortest round-tripping digits in scientific form, e.g. "1.0000000149011612e-1".
    let scientific = format!("{magnitude:e}");
    let (mantissa, exponent) = scientific
        .split_once('e')
        .unwrap_or((scientific.as_str(), "0"));
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    let exponent: i32 = exponent.parse().unwrap_or(0);

    // The value is 0.<digits> * 10^decimal_point; Python's format_float_short
    // switches on this quantity.
    let decimal_point = exponent + 1;
    let digit_count = i32::try_from(digits.len()).unwrap_or(i32::MAX);

    if decimal_point <= -4 || decimal_point > 16 {
        out.push_str(&digits[..1]);
        if digits.len() > 1 {
            out.push('.');
            out.push_str(&digits[1..]);
        }
        let printed = decimal_point - 1;
        let _ = write!(
            out,
            "e{}{:02}",
            if printed < 0 { '-' } else { '+' },
            printed.unsigned_abs()
        );
    } else if decimal_point <= 0 {
        out.push_str("0.");
        for _ in 0..(-decimal_point) {
            out.push('0');
        }
        out.push_str(&digits);
    } else if decimal_point >= digit_count {
        out.push_str(&digits);
        for _ in 0..(decimal_point - digit_count) {
            out.push('0');
        }
        out.push_str(".0"); // Py_DTSF_ADD_DOT_0: a float never prints as a bare integer
    } else {
        let split = usize::try_from(decimal_point).unwrap_or(0);
        out.push_str(&digits[..split]);
        out.push('.');
        out.push_str(&digits[split..]);
    }
    out
}

/// Append one `\uXXXX` escape, lower-case hex as Python writes it.
fn push_escaped_code_unit(out: &mut String, unit: u32) {
    let _ = write!(out, "\\u{unit:04x}");
}

/// Encode one string as `json.dumps(..., ensure_ascii=True)` does.
///
/// Printable ASCII (space through `~`) passes through except `"` and `\`;
/// the JSON short escapes are used for `\n`, `\r`, `\t`, backspace and form
/// feed; every other character, DEL included, becomes `\uXXXX`, outside the
/// BMP as a surrogate pair.
fn encode_string(text: &str, out: &mut String) {
    out.push('"');
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            ' '..='~' => out.push(c),
            _ => {
                let code_point = u32::from(c);
                if code_point <= 0xFFFF {
                    push_escaped_code_unit(out, code_point);
                } else {
                    let rest = code_point - 0x1_0000;
                    push_escaped_code_unit(out, 0xD800 + (rest >> 10));
                    push_escaped_code_unit(out, 0xDC00 + (rest & 0x3FF));
                }
            }
        }
    }
    out.push('"');
}

fn encode_int32_list(values: &[i32], out: &mut String) {
    out.push('[');
    for (index, value) in values.iter().enumerate() {
        if index != 0 {
            out.push(',');
        }
        let _ = write!(out, "{value}");
    }
    out.push(']');
}

fn encode_string_list<S: AsRef<str>>(values: &[S], out: &mut String) {
    out.push('[');
    for (index, value) in values.iter().enumerate() {
        if index != 0 {
            out.push(',');
        }
        encode_string(value.as_ref(), out);
    }
    out.push(']');
}

/// Encode a set-valued string list: deduplicated and ordered by each
/// element's JSON encoding.
///
/// qpx builds a dict keyed by each element's encoding and walks the sorted
/// keys (`qpx/core/data/identity.py`), so equal elements collapse to one and
/// the survivors are ordered by their encoded form -- a set, not just a
/// sorted list.
fn encode_string_set<S: AsRef<str>>(values: &[S], out: &mut String) {
    let mut encoded: Vec<String> = values
        .iter()
        .map(|value| {
            let mut element = String::new();
            encode_string(value.as_ref(), &mut element);
            element
        })
        .collect();
    encoded.sort_unstable();
    encoded.dedup();
    out.push('[');
    out.push_str(&encoded.join(","));
    out.push(']');
}

fn encode_value(value: &Value, out: &mut String) {
    match value {
        Value::Null => out.push_str("null"),
        Value::String(text) => encode_string(text, out),
        Value::Int64(number) => {
            let _ = write!(out, "{number}");
        }
        Value::Float32(number) => out.push_str(&format_float(*number)),
        Value::Int32List(values) => encode_int32_list(values, out),
        Value::StringList(values) => encode_string_list(values, out),
    }
}

/// Canonical JSON encoding of an identity composite.
///
/// Ports `QPXIdentity::canonical`: `json.dumps` with compact separators and
/// ASCII-escaped strings -- the exact bytes qpx hashes, e.g.
/// `["BSA1_F1",[2311],"PEPTIDER",2]`.
///
/// `unordered_list_indices` names the positions whose list value is a *set*;
/// a string list there is deduplicated and sorted by its elements' JSON
/// encoding before it is written. QPX uses this for `grouped_runs` and
/// `pg_accessions`. Ordered lists such as `scan` keep their order, and an
/// integer list is never reordered even at an unordered position (the source
/// applies the set treatment to string lists only).
pub fn canonical(values: &[Value], unordered_list_indices: &[usize]) -> String {
    let mut out = String::from("[");
    for (index, value) in values.iter().enumerate() {
        if index != 0 {
            out.push(',');
        }
        match value {
            Value::StringList(list) if unordered_list_indices.contains(&index) => {
                encode_string_set(list, &mut out);
            }
            _ => encode_value(value, &mut out),
        }
    }
    out.push(']');
    out
}

/// BLAKE2b with an 8-byte digest, read as a signed big-endian integer.
fn blake2b64(message: &str) -> i64 {
    let digest = Blake2b::<U8>::digest(message.as_bytes());
    let mut bytes = [0u8; 8];
    bytes.copy_from_slice(&digest);
    i64::from_be_bytes(bytes)
}

/// Derive an opaque signed 64-bit identity from a composite.
///
/// Ports `QPXIdentity::deriveId`: BLAKE2b-64 of [`canonical`]`(values,
/// unordered_list_indices)`. Uniqueness is not guaranteed by construction but
/// by the primary-key uniqueness check of
/// [`QpxValueValidation`](crate::format::qpx_value_validation::QpxValueValidation),
/// which reports a collision rather than losing a row to it.
pub fn derive_id(values: &[Value], unordered_list_indices: &[usize]) -> i64 {
    blake2b64(&canonical(values, unordered_list_indices))
}

/// `feature_id` from a feature row's persisted column values.
///
/// Ports `QPXIdentity::featureId`. Composite ([`FEATURE_COMPOSITE`]):
/// `(run_file_name, peptidoform, charge, rt, scan, observed_mz)`. The floats
/// are what make an *unidentified* feature row identifiable at all -- such a
/// row has an empty peptidoform and an empty scan, and its run and charge are
/// shared with dozens of others.
///
/// # Arguments
///
/// * `run_file_name` -- bare run name (no directory, no extension)
/// * `peptidoform` -- ProForma notation; empty for an unidentified feature
/// * `charge` -- charge state as persisted (`int16`)
/// * `rt` -- retention time in seconds as persisted (`float32`), or `None`
///   when the column is null; a null is a key value of its own
/// * `scan` -- scan components; empty for an unidentified feature
/// * `observed_mz` -- experimental m/z as persisted (`float32`)
pub fn feature_id(
    run_file_name: &str,
    peptidoform: &str,
    charge: i64,
    rt: Option<f32>,
    scan: &[i32],
    observed_mz: f32,
) -> i64 {
    let mut out = String::from("[");
    encode_string(run_file_name, &mut out);
    out.push(',');
    encode_string(peptidoform, &mut out);
    let _ = write!(out, ",{charge},");
    match rt {
        Some(rt) => out.push_str(&format_float(rt)),
        None => out.push_str("null"),
    }
    out.push(',');
    encode_int32_list(scan, &mut out);
    out.push(',');
    out.push_str(&format_float(observed_mz));
    out.push(']');
    blake2b64(&out)
}

/// `psm_id` from a psm row's persisted column values.
///
/// Ports `QPXIdentity::psmId`. Composite ([`PSM_COMPOSITE`]):
/// `(run_file_name, scan, peptidoform, charge)` -- no floats, so a psm
/// identity is stable across platforms. `scan` keeps its order (a psm
/// identity is not a set); `charge` is the precursor charge as persisted
/// (`int16`).
pub fn psm_id(run_file_name: &str, scan: &[i32], peptidoform: &str, charge: i64) -> i64 {
    let mut out = String::from("[");
    encode_string(run_file_name, &mut out);
    out.push(',');
    encode_int32_list(scan, &mut out);
    out.push(',');
    encode_string(peptidoform, &mut out);
    let _ = write!(out, ",{charge}]");
    blake2b64(&out)
}

/// `pg_id` from a protein-group row's persisted column values.
///
/// Ports `QPXIdentity::pgId`. Composite ([`PG_COMPOSITE`]):
/// `(pg_accessions, grouped_runs, label)`. Both lists are treated as sets --
/// deduplicated and sorted by their JSON encoding -- because neither the order
/// in which a group lists its members nor the order of the runs aggregated
/// into one quantity is part of the group's identity. qpx picks the unordered
/// positions by field name, which for this composite are positions 0 and 1.
///
/// The FULL membership keys the id, not the leading protein alone: two
/// distinct groups that share a leader (`P1;P2` and `P1;P3`) would otherwise
/// derive the same `pg_id`, and since the id is the view's primary key that
/// collision would be a refused export.
///
/// `label` is the channel label, or `None` for an identification-only group;
/// a null label is distinct from an empty one.
pub fn pg_id<A: AsRef<str>, R: AsRef<str>>(
    pg_accessions: &[A],
    grouped_runs: &[R],
    label: Option<&str>,
) -> i64 {
    let mut out = String::from("[");
    encode_string_set(pg_accessions, &mut out);
    out.push(',');
    encode_string_set(grouped_runs, &mut out);
    out.push(',');
    match label {
        Some(label) => encode_string(label, &mut out),
        None => out.push_str("null"),
    }
    out.push(']');
    blake2b64(&out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The per-view helpers encode inline for speed; they must produce the
    /// bytes the generic encoder produces for the same composite.
    #[test]
    fn helpers_agree_with_the_generic_encoder() {
        let feature = derive_id(
            &[
                Value::String("run".into()),
                Value::String("PEP".into()),
                Value::Int64(-3),
                Value::Null,
                Value::Int32List(vec![5, 1]),
                Value::Float32(2.5),
            ],
            &[],
        );
        assert_eq!(feature_id("run", "PEP", -3, None, &[5, 1], 2.5), feature);

        let psm = derive_id(
            &[
                Value::String("r\u{e4}".into()),
                Value::Int32List(vec![]),
                Value::String("A\"B".into()),
                Value::Int64(0),
            ],
            &[],
        );
        assert_eq!(psm_id("r\u{e4}", &[], "A\"B", 0), psm);

        let pg = derive_id(
            &[
                Value::StringList(vec!["b".into(), "a".into(), "b".into()]),
                Value::StringList(vec!["r2".into(), "r1".into()]),
                Value::String("LFQ".into()),
            ],
            &[0, 1],
        );
        assert_eq!(pg_id(&["b", "a", "b"], &["r2", "r1"], Some("LFQ")), pg);
    }

    #[test]
    fn integer_lists_are_never_reordered() {
        assert_eq!(
            canonical(&[Value::Int32List(vec![3, 1, 2])], &[0]),
            "[[3,1,2]]"
        );
    }
}

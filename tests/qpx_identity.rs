// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `QPXIdentity` (`FORMAT/QPXIdentity.h`).
//!
//! Every expected value is transcribed from the upstream class test
//! `src/tests/class_tests/openms/source/QPXIdentity_test.cpp` at core
//! `bc9cc12`. Its header states that those values were produced by the
//! reference implementation `qpx.core.data.identity` (bigbio/qpx), not by
//! running the C++, so they are independent of both implementations.
#![cfg(feature = "arrow")]
// The float literals are the class test's, verbatim: several are written with
// more digits than a float holds, which is the point of the narrowing checks.
#![allow(clippy::excessive_precision, clippy::inconsistent_digit_grouping)]

use openms::format::qpx_identity::{
    FEATURE_COMPOSITE, PG_COMPOSITE, PSM_COMPOSITE, Value, canonical, derive_id, feature_id,
    format_float, pg_id, psm_id,
};

fn s(text: &str) -> Value {
    Value::String(text.to_owned())
}

fn strings(values: &[&str]) -> Value {
    Value::StringList(values.iter().map(|v| (*v).to_owned()).collect())
}

/// `START_SECTION((std::string formatFloat(float value)))`
#[test]
fn format_float_follows_python_repr() {
    assert_eq!(format_float(0.0), "0.0");
    assert_eq!(format_float(-0.0), "-0.0");
    assert_eq!(format_float(1.0), "1.0");
    assert_eq!(format_float(-1.0), "-1.0");
    assert_eq!(format_float(123456.0), "123456.0");

    // Widened to double first: 0.1f is 0.10000000149011612.
    assert_eq!(format_float(0.1), "0.10000000149011612");
    assert_eq!(format_float(549.8573), "549.8572998046875");
    assert_eq!(format_float(1534.9807), "1534.980712890625");

    // The fixed/exponential switch on both sides of each boundary.
    assert_eq!(format_float(1125899906842624.0), "1125899906842624.0"); // 2^50, decpt 16
    assert_eq!(format_float(18014398509481984.0), "1.8014398509481984e+16"); // 2^54, decpt 17
    assert_eq!(format_float(0.0001220703125), "0.0001220703125"); // 2^-13, decpt -3
    assert_eq!(format_float(0.00006103515625), "6.103515625e-05"); // 2^-14, decpt -4
    assert_eq!(format_float(1e-4), "9.999999747378752e-05");

    // Subnormal and maximum float32.
    assert_eq!(format_float(1e-45), "1.401298464324817e-45");
    assert_eq!(format_float(3.4028235e38), "3.4028234663852886e+38");

    // Written bare, as json.dumps does.
    assert_eq!(format_float(f32::INFINITY), "Infinity");
    assert_eq!(format_float(f32::NEG_INFINITY), "-Infinity");
    assert_eq!(format_float(f32::NAN), "NaN");
}

/// `START_SECTION((std::string canonical(...)))`
#[test]
fn canonical_is_compact_ascii_json() {
    assert_eq!(
        canonical(
            &[
                s("BSA1_F1"),
                Value::Int32List(vec![2311]),
                s("PEPTIDER"),
                Value::Int64(2)
            ],
            &[]
        ),
        "[\"BSA1_F1\",[2311],\"PEPTIDER\",2]"
    );
    assert_eq!(canonical(&[Value::Int64(1), Value::Int64(2)], &[]), "[1,2]");

    assert_eq!(canonical(&[Value::Null], &[]), "[null]");
    assert_ne!(canonical(&[Value::Null], &[]), canonical(&[s("")], &[]));

    assert_eq!(canonical(&[s("B\u{e4}cker")], &[]), "[\"B\\u00e4cker\"]");
    assert_eq!(canonical(&[s("\u{1F600}")], &[]), "[\"\\ud83d\\ude00\"]");
    assert_eq!(canonical(&[s("a\"b\\c\nd")], &[]), "[\"a\\\"b\\\\c\\nd\"]");
    assert_eq!(canonical(&[s("a/b\u{7f}")], &[]), "[\"a/b\\u007f\"]");

    assert_ne!(
        canonical(&[s("a,b"), s("c")], &[]),
        canonical(&[s("a"), s("b,c")], &[])
    );

    let forward = strings(&["BSA1_F1", "BSA1_F2"]);
    let reversed = strings(&["BSA1_F2", "BSA1_F1"]);
    assert_eq!(
        canonical(std::slice::from_ref(&reversed), &[0]),
        "[[\"BSA1_F1\",\"BSA1_F2\"]]"
    );
    assert_eq!(
        canonical(std::slice::from_ref(&forward), &[0]),
        canonical(std::slice::from_ref(&reversed), &[0])
    );
    assert_ne!(
        canonical(std::slice::from_ref(&forward), &[]),
        canonical(std::slice::from_ref(&reversed), &[])
    );
    assert_eq!(
        canonical(&[Value::Int32List(vec![3, 1, 2])], &[]),
        "[[3,1,2]]"
    );
}

/// `START_SECTION((Int64 deriveId(...)))`
#[test]
fn derive_id_matches_qpx() {
    assert_eq!(
        derive_id(
            &[
                s("BSA1_F1"),
                Value::Int32List(vec![2311]),
                s("PEPTIDER"),
                Value::Int64(2)
            ],
            &[]
        ),
        249030645178071364
    );
    assert_eq!(
        derive_id(
            &[
                s("B\u{e4}cker/prote\u{ef}n"),
                Value::Int32List(vec![1, 2, 3]),
                s("PEP[Phospho]TIDE"),
                Value::Int64(3)
            ],
            &[]
        ),
        -3403684637348881565
    );
    assert_eq!(
        derive_id(
            &[
                s("run,with]brackets"),
                Value::Int32List(vec![]),
                s(""),
                Value::Int64(0)
            ],
            &[]
        ),
        260970448792327776
    );
}

/// `START_SECTION((Int64 featureId(...)))`
#[test]
fn feature_id_matches_qpx() {
    assert_eq!(
        feature_id("BSA1_F1", "PEPTIDER", 2, Some(1534.9807), &[2311], 549.8573),
        -9185367007480417157
    );
    assert_eq!(
        feature_id("BSA1_F1", "", 2, Some(1534.9807), &[], 549.8573),
        6679072350175538872
    );

    let base = feature_id("BSA1_F1", "PEPTIDER", 2, Some(1534.9807), &[2311], 549.8573);
    assert_ne!(
        feature_id("BSA1_F2", "PEPTIDER", 2, Some(1534.9807), &[2311], 549.8573),
        base
    );
    assert_ne!(
        feature_id("BSA1_F1", "PEPTIDEK", 2, Some(1534.9807), &[2311], 549.8573),
        base
    );
    assert_ne!(
        feature_id("BSA1_F1", "PEPTIDER", 3, Some(1534.9807), &[2311], 549.8573),
        base
    );
    assert_ne!(
        feature_id("BSA1_F1", "PEPTIDER", 2, Some(1534.9808), &[2311], 549.8573),
        base
    );
    assert_ne!(
        feature_id("BSA1_F1", "PEPTIDER", 2, Some(1534.9807), &[2312], 549.8573),
        base
    );
    assert_ne!(
        feature_id("BSA1_F1", "PEPTIDER", 2, Some(1534.9807), &[2311], 549.8574),
        base
    );
    assert_ne!(
        feature_id("BSA1_F1", "PEPTIDER", 2, None, &[2311], 549.8573),
        base
    );
    // A double that narrows to the same float32 gives the same id.
    let narrowed = 1534.98071289062500001_f64 as f32;
    assert_eq!(
        feature_id("BSA1_F1", "PEPTIDER", 2, Some(narrowed), &[2311], 549.8573),
        base
    );
}

/// `START_SECTION((Int64 psmId(...)))`
#[test]
fn psm_id_matches_qpx() {
    assert_eq!(
        psm_id("BSA1_F1", &[2311], "PEPTIDER", 2),
        249030645178071364
    );
    assert_ne!(
        psm_id("BSA1_F1", &[1, 2], "PEPTIDER", 2),
        psm_id("BSA1_F1", &[2, 1], "PEPTIDER", 2)
    );
    assert_ne!(
        psm_id("BSA1_F1", &[2311], "PEPTIDER", 3),
        psm_id("BSA1_F1", &[2311], "PEPTIDER", 2)
    );
}

/// `START_SECTION((Int64 pgId(...)))`
#[test]
fn pg_id_matches_qpx() {
    let none: &[&str] = &[];
    assert_eq!(
        pg_id(&["P02769"], &["BSA1_F1", "BSA1_F2"], Some("LFQ")),
        4537562256413832418
    );
    assert_eq!(
        pg_id(&["P02769"], &["BSA1_F2", "BSA1_F1"], Some("LFQ")),
        4537562256413832418
    );
    assert_eq!(
        pg_id(&["P02769", "P00761"], &["BSA1_F1"], Some("LFQ")),
        7169929772723875882
    );
    assert_eq!(
        pg_id(&["P00761", "P02769"], &["BSA1_F1"], Some("LFQ")),
        7169929772723875882
    );
    assert_eq!(
        pg_id(&["P02769", "P02769"], &["BSA1_F1"], Some("LFQ")),
        -3038322759200455733
    );
    assert_eq!(
        pg_id(&["P02769"], &["BSA1_F1"], Some("LFQ")),
        -3038322759200455733
    );
    assert_ne!(
        pg_id(&["P1", "P2"], &["r"], Some("LFQ")),
        pg_id(&["P1", "P3"], &["r"], Some("LFQ"))
    );
    assert_eq!(pg_id(&["P02769"], &["BSA1_F1"], None), -6310942379590911755);
    assert_ne!(
        pg_id(&["P02769"], &["BSA1_F1"], Some("")),
        pg_id(&["P02769"], &["BSA1_F1"], None)
    );
    assert_ne!(
        pg_id(&["P02769"], &["BSA1_F1"], Some("LFQ")),
        pg_id(&["P02769"], &["BSA1_F1", "BSA1_F2"], Some("LFQ"))
    );
    assert_eq!(pg_id(none, &["BSA1_F1"], Some("LFQ")), -1719918802529113715);
}

/// `START_SECTION([EXTRA] identity_composite footer declarations)`
#[test]
fn composite_declarations() {
    assert_eq!(
        FEATURE_COMPOSITE,
        "run_file_name,peptidoform,charge,rt,scan,observed_mz"
    );
    assert_eq!(PSM_COMPOSITE, "run_file_name,scan,peptidoform,charge");
    assert_eq!(PG_COMPOSITE, "pg_accessions,grouped_runs,label");
}

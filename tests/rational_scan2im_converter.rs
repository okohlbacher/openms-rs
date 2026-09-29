//! `RationalScan2ImConverter` against the five `RationalScan2ImConverter`
//! sections of the upstream class test `BrukerTimsFile_test.cpp` (core pin
//! bc9cc12, lines 129-285), transcribed with the same coefficients, inputs and
//! assertions (evidence tier 3: class-test literals). `TEST_REAL_SIMILAR` uses
//! the ClassTest defaults: relative 1e-5, absolute 1e-5.
//!
//! The remaining tests are native: the empty-map refusal, the fallback segment,
//! opentims' linear converter, and (`sqlite`) `tryCreateRationalConverter` on
//! synthetic `analysis.tdf` databases. Their expectations are derived from the
//! source's control flow and formulas (tier 4), not from C++ output.

use openms::format::rational_scan2im_converter::{
    Coefficients, LinearScan2ImConverter, RationalScan2ImConverter, Scan2InvIonMobilityConverter,
};
use std::collections::BTreeMap;

fn similar(actual: f64, expected: f64) -> bool {
    if actual == expected || (actual - expected).abs() <= 1e-5 {
        return true;
    }
    if actual == 0.0 || expected == 0.0 {
        return false;
    }
    let ratio = actual / expected;
    let ratio = if ratio < 1.0 { 1.0 / ratio } else { ratio };
    ratio > 0.0 && ratio <= 1.0 + 1e-5
}

/// The opentims ModelType-2 coefficients every section uses.
fn opentims_coefficients() -> Coefficients {
    Coefficients {
        c0: 1.0,
        c1: 917.0,
        c2: 213.5998,
        c3: 75.81729,
        c4: 33.0,
        c5: 1.0,
        c6: -0.009065829,
        c7: 135.4364,
        c8: 13.32608,
        c9: 1663.341,
    }
}

fn single(c: Coefficients) -> RationalScan2ImConverter {
    RationalScan2ImConverter::new(BTreeMap::from([(1, c)]), vec![0, 1]).unwrap()
}

#[test]
fn forward_conversion() {
    let converter = single(opentims_coefficients());
    let result = converter.convert(1, &[500.0]);
    assert!(similar(result[0], 1.070429), "{}", result[0]);
    let result2 = converter.convert_indices(1, &[500]);
    assert!(similar(result2[0], 1.070429), "{}", result2[0]);

    let results = converter.convert(1, &[0.0, 250.0, 500.0, 916.0]);
    assert!(similar(results[0], 1.63883), "{}", results[0]);
    assert!(results.iter().all(|&r| r > 0.0));
    for i in 1..4 {
        assert!(results[i] < results[i - 1]);
    }
}

#[test]
fn round_trip_via_inverse_convert() {
    let converter = single(opentims_coefficients());
    let inv_k0 = converter.convert(1, &[500.0]);
    assert_eq!(converter.inverse_convert(1, &inv_k0), vec![500]);
    for s in (10u32..900).step_by(50) {
        let ik0 = converter.convert(1, &[f64::from(s)]);
        let back = converter.inverse_convert(1, &ik0)[0];
        assert!(back >= s - 1 && back <= s + 1, "{s} -> {back}");
    }
}

#[test]
fn per_frame_calibration() {
    let c2 = Coefficients {
        c0: 1.0,
        c1: 917.0,
        c2: 220.0,
        c3: 80.0,
        c4: 33.0,
        c5: 1.0,
        c6: -0.01,
        c7: 140.0,
        c8: 13.0,
        c9: 1660.0,
    };
    let converter = RationalScan2ImConverter::new(
        BTreeMap::from([(1, opentims_coefficients()), (2, c2)]),
        vec![0, 1, 2],
    )
    .unwrap();
    let f1 = converter.convert(1, &[500.0])[0];
    let f2 = converter.convert(2, &[500.0])[0];
    assert_ne!(f1, f2);
    assert!(f1 > 0.0);
    assert!(f2 > 0.0);
}

#[test]
fn description() {
    let desc = single(opentims_coefficients()).description();
    assert!(desc.contains("RationalScan2ImConverter"));
    assert!(desc.contains('1'));
}

#[test]
fn singularity_edge_cases() {
    let mut zero_slope = opentims_coefficients();
    zero_slope.c3 = zero_slope.c2;
    assert!(single(zero_slope).convert(1, &[500.0])[0].is_finite());

    let mut v_zero = opentims_coefficients();
    v_zero.c2 = 0.0;
    v_zero.c3 = 0.0;
    assert!(single(v_zero).convert(1, &[500.0])[0].is_finite());
}

/// Native: an empty calibration map is refused (the source dereferences
/// `begin()` of the empty map on the first conversion).
#[test]
fn empty_calibrations_are_refused() {
    assert!(RationalScan2ImConverter::new(BTreeMap::new(), vec![0]).is_err());
}

/// Out-of-range frames use the first segment, as the source's
/// `getCalibration` does; here "first" is the smallest Id.
#[test]
fn out_of_range_frames_use_the_first_segment() {
    let mut other = opentims_coefficients();
    other.c7 = 140.0;
    let converter = RationalScan2ImConverter::new(
        BTreeMap::from([(7, other), (3, opentims_coefficients())]),
        vec![0, 7],
    )
    .unwrap();
    let first = converter.convert(0, &[500.0]);
    assert_eq!(first, converter.convert(99, &[500.0]));
    assert!(similar(first[0], 1.070429));
    assert_ne!(converter.convert(1, &[500.0]), first);
}

/// opentims `OpenSourceScan2ImConverter`: `1/K0 = im_max + (im_min - im_max)
/// / scan_max * scan` (independently derived; values chosen exact in binary).
#[test]
fn linear_converter() {
    let linear = LinearScan2ImConverter::new(0.5, 1.5, 1024);
    assert_eq!(linear.intercept(), 1.5);
    assert_eq!(linear.slope(), -1.0 / 1024.0);
    assert_eq!(linear.convert(1, &[0.0, 1024.0]), vec![1.5, 0.5]);
    assert_eq!(linear.inverse_convert(1, &[1.5, 0.5]), vec![0, 1024]);
}

#[cfg(feature = "sqlite")]
mod tdf {
    use super::*;
    use openms::format::rational_scan2im_converter::try_create_rational_converter;
    use rusqlite::Connection;
    use std::path::PathBuf;

    fn tdf(name: &str, sql: &str) -> PathBuf {
        let dir =
            std::env::temp_dir().join(format!("openms-rational-{}-{name}.d", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let db = Connection::open(dir.join("analysis.tdf")).unwrap();
        db.execute_batch(sql).unwrap();
        dir
    }

    const CAL: &str = "CREATE TABLE TimsCalibration (Id INTEGER PRIMARY KEY, ModelType INTEGER, \
        C0 REAL, C1 REAL, C2 REAL, C3 REAL, C4 REAL, C5 REAL, C6 REAL, C7 REAL, C8 REAL, C9 REAL);";

    #[test]
    fn reads_calibrations_and_frames() {
        let dir = tdf(
            "ok",
            &format!(
                "{CAL}
                INSERT INTO TimsCalibration VALUES (1, 2, 1.0, 917.0, 213.5998, 75.81729, 33.0, 1.0, -0.009065829, 135.4364, 13.32608, 1663.341);
                CREATE TABLE Frames (Id INTEGER PRIMARY KEY, TimsCalibration INTEGER);
                INSERT INTO Frames VALUES (1, 1), (3, 1);"
            ),
        );
        let converter = try_create_rational_converter(&dir).unwrap();
        assert_eq!(converter.frame_to_cal(), &[0, 1, 0, 1]);
        assert_eq!(converter.calibrations()[&1], opentims_coefficients());
        assert!(similar(converter.convert(3, &[500.0])[0], 1.070429));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn falls_back_on_every_source_condition() {
        let one = "INSERT INTO TimsCalibration VALUES (1, 2, 0,1,0,0,0,0,0,0,0,0);";
        let frames = "CREATE TABLE Frames (Id INTEGER, TimsCalibration INTEGER);";
        let cases = [
            ("missing", "CREATE TABLE Other (x);".to_string()),
            (
                "model",
                format!(
                    "{CAL} INSERT INTO TimsCalibration VALUES (1, 1, 0,0,0,0,0,0,0,0,0,0);
                    {frames} INSERT INTO Frames VALUES (1, 1);"
                ),
            ),
            ("empty", format!("{CAL} {frames}")),
            (
                "oldframes",
                format!(
                    "{CAL} {one} CREATE TABLE Frames (Id INTEGER); INSERT INTO Frames VALUES (1);"
                ),
            ),
            (
                "null",
                format!("{CAL} {one} {frames} INSERT INTO Frames VALUES (1, NULL);"),
            ),
            (
                "unknown",
                format!("{CAL} {one} {frames} INSERT INTO Frames VALUES (1, 5);"),
            ),
            ("noframes", format!("{CAL} {one} {frames}")),
            (
                "ceiling",
                format!("{CAL} {one} {frames} INSERT INTO Frames VALUES (4294967295, 1);"),
            ),
        ];
        for (name, sql) in cases {
            let dir = tdf(name, &sql);
            assert!(try_create_rational_converter(&dir).is_none(), "{name}");
            let _ = std::fs::remove_dir_all(dir);
        }
        assert!(try_create_rational_converter("/nonexistent/path.d").is_none());
    }
}

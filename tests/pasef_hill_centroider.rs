//! `PASEFHillCentroider::centroidFrame` against every section of the upstream
//! class test `PASEFHillCentroider_test.cpp` (core pin bc9cc12), transcribed
//! section by section with the same inputs, parameters and assertions
//! (evidence tier 3: class-test literals). `TEST_REAL_SIMILAR` uses the
//! ClassTest defaults: relative 1e-5, absolute 1e-5.

use openms::format::pasef_hill_centroider::{Centroid, Params, centroid_frame};

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

macro_rules! assert_similar {
    ($actual:expr, $expected:expr) => {{
        let (a, e): (f64, f64) = ($actual, $expected);
        assert!(similar(a, e), "{} = {a}, expected {e}", stringify!($actual));
    }};
}

/// `buildFlat`: one IM value per scan, `(m/z, intensity)` peaks.
fn flat(scans: &[(f64, &[(f64, f64)])]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (mut mz, mut intensity, mut im) = (Vec::new(), Vec::new(), Vec::new());
    for &(scan_im, peaks) in scans {
        for &(m, i) in peaks {
            mz.push(m);
            intensity.push(i);
            im.push(scan_im);
        }
    }
    (mz, intensity, im)
}

/// `defaultParams()` of the class test.
fn default_params() -> Params {
    Params {
        mz_tol_ppm: 20.0,
        valley_factor: 1.3,
        min_hill_length: 2,
        min_intensity: 0.0,
        im_group_tolerance: 0.0,
        ..Params::default()
    }
}

fn run(scans: &[(f64, &[(f64, f64)])], params: &Params) -> Vec<Centroid> {
    let (mz, intensity, im) = flat(scans);
    centroid_frame(&mz, &intensity, &im, params, None).unwrap()
}

#[test]
fn source_defaults() {
    let p = Params::default();
    assert_eq!(p.mz_tol_ppm, 8.0);
    assert_eq!(p.valley_factor, 1.3);
    assert_eq!(p.min_hill_length, 1);
    assert_eq!(p.min_intensity, 0.0);
    assert_eq!(p.im_group_tolerance, 0.0);
    assert_eq!(p.max_scan_gap, 0);
}

#[test]
fn empty_input_returns_empty_output() {
    assert!(
        centroid_frame(&[], &[], &[], &default_params(), None)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn isolated_peak_across_n_im_scans_yields_one_centroid() {
    let result = run(
        &[
            (1.00, &[(500.0, 100.0)]),
            (1.05, &[(500.0, 400.0)]),
            (1.10, &[(500.0, 800.0)]),
            (1.15, &[(500.0, 400.0)]),
            (1.20, &[(500.0, 100.0)]),
        ],
        &default_params(),
    );
    assert_eq!(result.len(), 1);
    assert_similar!(result[0].mz, 500.0);
    assert_similar!(result[0].intensity, 100.0 + 400.0 + 800.0 + 400.0 + 100.0);
    assert_similar!(result[0].im_apex, 1.10);
    assert_similar!(result[0].im_lower, 1.00);
    assert_similar!(result[0].im_upper, 1.20);
    assert_eq!(result[0].length, 5);
}

#[test]
fn two_co_mz_ions_separated_by_im_valley_are_split() {
    let result = run(
        &[
            (1.00, &[(500.0, 200.0)]),
            (1.05, &[(500.0, 700.0)]),
            (1.10, &[(500.0, 200.0)]),
            (1.15, &[(500.0, 50.0)]),
            (1.20, &[(500.0, 200.0)]),
            (1.25, &[(500.0, 700.0)]),
            (1.30, &[(500.0, 200.0)]),
        ],
        &default_params(),
    );
    assert_eq!(result.len(), 2);
    assert_similar!(result[0].mz, 500.0);
    assert_similar!(result[1].mz, 500.0);
    assert!(result[0].im_apex < 1.15);
    assert!(result[1].im_apex > 1.15);
}

#[test]
fn single_im_scan_blip_is_rejected_by_min_hill_length_2() {
    let result = run(
        &[
            (1.00, &[(500.0, 100.0)]),
            (1.05, &[(500.0, 400.0)]),
            (1.10, &[(500.0, 800.0)]),
            (1.15, &[(500.0, 400.0)]),
            (1.20, &[(500.0, 100.0), (600.0, 999.0)]),
        ],
        &default_params(),
    );
    assert_eq!(result.len(), 1);
    assert_similar!(result[0].mz, 500.0);
}

#[test]
fn min_hill_length_1_keeps_single_scan_peaks() {
    let params = Params {
        min_hill_length: 1,
        ..default_params()
    };
    let result = run(
        &[
            (1.00, &[(500.0, 100.0)]),
            (1.05, &[(500.0, 400.0)]),
            (1.10, &[(500.0, 800.0), (600.0, 999.0)]),
        ],
        &params,
    );
    assert_eq!(result.len(), 2);
    assert_similar!(result[0].mz, 500.0);
    assert_similar!(result[1].mz, 600.0);
}

fn gap_case(
    intensity: &[f64],
    im: &[f64],
    scans: &[u32],
    max_scan_gap: usize,
    min_hill_length: usize,
) -> Vec<Centroid> {
    let mz = vec![500.0; intensity.len()];
    let params = Params {
        max_scan_gap,
        min_hill_length,
        ..default_params()
    };
    centroid_frame(&mz, intensity, im, &params, Some(scans)).unwrap()
}

#[test]
fn max_scan_gap_0_splits_a_hill_broken_by_one_empty_scan() {
    let result = gap_case(
        &[100.0, 400.0, 400.0, 100.0],
        &[1.00, 1.05, 1.15, 1.20],
        &[0, 1, 3, 4],
        0,
        2,
    );
    assert_eq!(result.len(), 2);
}

#[test]
fn max_scan_gap_1_bridges_a_single_empty_im_scan() {
    let result = gap_case(
        &[100.0, 400.0, 400.0, 100.0],
        &[1.00, 1.05, 1.15, 1.20],
        &[0, 1, 3, 4],
        1,
        2,
    );
    assert_eq!(result.len(), 1);
    assert_similar!(result[0].mz, 500.0);
    assert_eq!(result[0].length, 4);
}

#[test]
fn max_scan_gap_1_does_not_let_a_length_1_hill_bridge_a_gap() {
    let result = gap_case(&[100.0, 200.0], &[1.00, 1.10], &[0, 2], 1, 2);
    assert_eq!(result.len(), 0);
}

#[test]
fn max_scan_gap_1_lets_a_paired_hill_bridge_a_subsequent_gap() {
    let result = gap_case(
        &[100.0, 400.0, 200.0],
        &[1.00, 1.05, 1.15],
        &[0, 1, 3],
        1,
        2,
    );
    assert_eq!(result.len(), 1);
    assert_eq!(result[0].length, 3);
}

#[test]
fn max_scan_gap_1_does_not_bridge_two_consecutive_empty_scans() {
    let result = gap_case(
        &[100.0, 400.0, 400.0, 100.0],
        &[1.00, 1.05, 1.20, 1.25],
        &[0, 1, 4, 5],
        1,
        2,
    );
    assert_eq!(result.len(), 2);
}

#[test]
fn two_close_mz_ions_outside_tolerance_are_kept_separate() {
    let result = run(
        &[
            (1.00, &[(500.0, 200.0), (500.1, 150.0)]),
            (1.05, &[(500.0, 700.0), (500.1, 500.0)]),
            (1.10, &[(500.0, 200.0), (500.1, 150.0)]),
        ],
        &default_params(),
    );
    assert_eq!(result.len(), 2);
    assert!(result[0].mz < result[1].mz);
    assert_similar!(result[0].mz, 500.0);
    assert_similar!(result[1].mz, 500.1);
}

#[test]
fn determinism_same_input_produces_identical_output() {
    let scans: &[(f64, &[(f64, f64)])] = &[
        (1.00, &[(500.0, 100.0), (600.0, 50.0)]),
        (1.05, &[(500.0, 400.0), (600.0, 200.0)]),
        (1.10, &[(500.0, 800.0), (600.0, 400.0)]),
        (1.15, &[(500.0, 400.0), (600.0, 200.0)]),
        (1.20, &[(500.0, 100.0), (600.0, 50.0)]),
    ];
    let r1 = run(scans, &default_params());
    let r2 = run(scans, &default_params());
    assert_eq!(r1.len(), r2.len());
    for (a, b) in r1.iter().zip(&r2) {
        assert_similar!(a.mz, b.mz);
        assert_similar!(a.intensity, b.intensity);
        assert_similar!(a.im_apex, b.im_apex);
        assert_eq!(a.length, b.length);
    }
}

#[test]
fn min_intensity_filter_drops_low_intensity_input_peaks() {
    let params = Params {
        min_intensity: 1.0,
        ..default_params()
    };
    let result = run(
        &[
            (1.00, &[(500.0, 100.0), (600.0, 0.5)]),
            (1.05, &[(500.0, 400.0), (600.0, 0.5)]),
            (1.10, &[(500.0, 800.0), (600.0, 0.5)]),
            (1.15, &[(500.0, 400.0)]),
            (1.20, &[(500.0, 100.0)]),
        ],
        &params,
    );
    assert_eq!(result.len(), 1);
    assert_similar!(result[0].mz, 500.0);
}

#[test]
fn no_valley_split_when_intensities_form_a_smooth_single_peak() {
    let result = run(
        &[
            (1.00, &[(500.0, 100.0)]),
            (1.05, &[(500.0, 300.0)]),
            (1.10, &[(500.0, 600.0)]),
            (1.15, &[(500.0, 1000.0)]),
            (1.20, &[(500.0, 600.0)]),
            (1.25, &[(500.0, 300.0)]),
            (1.30, &[(500.0, 100.0)]),
        ],
        &default_params(),
    );
    assert_eq!(result.len(), 1);
    assert_eq!(result[0].length, 7);
}

#[test]
fn output_is_sorted_by_mz_ascending() {
    let result = run(
        &[
            (1.00, &[(700.0, 50.0), (500.0, 100.0), (600.0, 80.0)]),
            (1.05, &[(700.0, 50.0), (500.0, 400.0), (600.0, 80.0)]),
            (1.10, &[(700.0, 50.0), (500.0, 800.0), (600.0, 80.0)]),
        ],
        &default_params(),
    );
    assert_eq!(result.len(), 3);
    for pair in result.windows(2) {
        assert!(pair[1].mz >= pair[0].mz);
    }
}

/// Parse a C99 `%a` / Python `float.hex()` hexadecimal float: `[-]0x1.hhhp±d`,
/// `0x0p+0`, or a subnormal `0x0.hhhp-1022`.
fn hex_f64(text: &str) -> f64 {
    let (negative, rest) = match text.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, text),
    };
    let rest = rest.strip_prefix("0x").expect("hex float");
    let (mantissa, exponent) = rest.split_once('p').expect("hex exponent");
    let exponent: i32 = exponent.parse().expect("exponent");
    let (lead, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let lead = u64::from_str_radix(lead, 16).expect("lead digit");
    assert!(lead <= 1 && fraction.len() <= 13, "{text}");
    let bits_fraction = if fraction.is_empty() {
        0
    } else {
        u64::from_str_radix(fraction, 16).expect("fraction") << (4 * (13 - fraction.len()))
    };
    let magnitude = if lead == 0 && bits_fraction == 0 {
        0.0
    } else if lead == 0 {
        assert_eq!(exponent, -1022, "{text}");
        f64::from_bits(bits_fraction)
    } else {
        let biased = u64::try_from(exponent + 1023).expect("normal exponent");
        f64::from_bits((biased << 52) | bits_fraction)
    };
    if negative { -magnitude } else { magnitude }
}

/// Run one case of `tests/data/pasef_hill_centroider_oracle/` and compare
/// every field of every centroid bit for bit with the Release build's output.
fn oracle_case(name: &str) -> Result<(), String> {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/pasef_hill_centroider_oracle");
    let input = std::fs::read_to_string(dir.join(format!("pasef_{name}.txt"))).unwrap();
    let mut lines = input.lines();
    let header: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
    let params = Params {
        mz_tol_ppm: hex_f64(header[0]),
        valley_factor: hex_f64(header[1]),
        min_hill_length: header[2].parse().unwrap(),
        min_intensity: hex_f64(header[3]),
        im_group_tolerance: hex_f64(header[4]),
        max_scan_gap: header[5].parse().unwrap(),
    };
    let use_scan = header[6] == "1";
    let (mut mz, mut intensity, mut im, mut scan) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for line in lines.filter(|line| !line.is_empty()) {
        let fields: Vec<&str> = line.split_whitespace().collect();
        mz.push(hex_f64(fields[0]));
        intensity.push(hex_f64(fields[1]));
        im.push(hex_f64(fields[2]));
        scan.push(fields[3].parse::<u32>().unwrap());
    }
    let actual = centroid_frame(&mz, &intensity, &im, &params, use_scan.then_some(&scan[..]))
        .map_err(|e| e.to_string())?;

    let expected = std::fs::read_to_string(dir.join(format!("pasef_{name}.tsv"))).unwrap();
    let mut expected = expected.lines();
    let count: usize = expected
        .next()
        .unwrap()
        .strip_prefix("n\t")
        .unwrap()
        .parse()
        .unwrap();
    if actual.len() != count {
        return Err(format!(
            "{name}: {} centroids, oracle {count}",
            actual.len()
        ));
    }
    for (index, (centroid, line)) in actual.iter().zip(expected).enumerate() {
        let f: Vec<&str> = line.split('\t').collect();
        let want: Vec<f64> = f[..7].iter().map(|field| hex_f64(field)).collect();
        let got = [
            centroid.mz,
            centroid.intensity,
            centroid.im_apex,
            centroid.im_lower,
            centroid.im_upper,
            centroid.mz_lower,
            centroid.mz_upper,
        ];
        for (field, (g, w)) in got.iter().zip(want).enumerate() {
            if g.to_bits() != w.to_bits() {
                return Err(format!(
                    "{name}: centroid {index} field {field}: {g} != oracle {w}"
                ));
            }
        }
        let length: usize = f[7].parse().unwrap();
        if centroid.length != length {
            return Err(format!(
                "{name}: centroid {index} length {} != oracle {length}",
                centroid.length
            ));
        }
    }
    Ok(())
}

/// Tier 1: synthetic TIMS-PASEF frames (`../oracle/F5-xmass/gen_cases.py`,
/// seed 20260928, continuous intensities so no two peaks of a scan tie) run
/// through `Internal::PASEFHillCentroider::centroidFrame` of the Linux x86_64
/// Release install `openms4-release-bc9cc12-c19e494-174b576` on ibminode06,
/// twice, with identical output. Cases cover the MS1 and MS2 defaults of
/// `BrukerTimsFile::Config`, `max_scan_gap` 1 and 2 with dropped scans,
/// `min_intensity`, and the IM-value grouping path without scan indices, with
/// exact and with tolerant grouping.
#[test]
fn matches_the_release_build_bit_for_bit() {
    for name in [
        "ms1_default",
        "ms2_hill2",
        "gap1",
        "gap2_minint",
        "no_scan_idx",
        "im_tol",
    ] {
        oracle_case(name).unwrap();
    }
}

/// The same oracle on integer (detector-count) intensities, where peaks of
/// one scan can tie. The source sorts with `std::sort`, whose order among
/// equal keys is unspecified; the port sorts stably. This case measures
/// whether that difference is observable on tie-rich input.
#[test]
fn integer_intensities_match_the_release_build() {
    oracle_case("int_ties").unwrap();
}

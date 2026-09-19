// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Differential cover for the two peak pickers that build their own
//! `SignalToNoiseEstimatorMedian`: `PeakPickerChromatogram` (source
//! `ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp:171`, `snt_.init(chromatogram)`)
//! and `PeakPickerIterative` (source `PROCESSING/CENTROIDING/PeakPickerIterative.h:321`,
//! `snt.init(input)`).
//!
//! Every expected value here was measured on the Linux x86_64 Release build
//! `openms4-release-bc9cc12-c19e494-174b576`; the drivers, the case generator
//! and the two identical repeat runs are under `oracle/picker-consumers/`, and
//! `tests/data/picker_consumers_provenance.json` records the source hashes.
//! Nothing in this file is derived from Rust output.

use openms::Error;
use openms::kernel::{ChromatogramPeak, MSChromatogram, MSSpectrum, Peak1D};
use openms::processing::chromatogram::{
    ChromatogramPickingMethod, ChromatogramSmoothing, PeakPickerChromatogram,
};
use openms::processing::iterative::PeakPickerIterative;
use openms::processing::peak_picking::{PickingCompatibility, SignalToNoiseEstimatorMedian};
use std::collections::BTreeMap;

const DATA: &str = "tests/data/picker_consumers";

// ---------------------------------------------------------------- fixtures

fn parse_scalar(token: &str) -> f64 {
    if let Some(hex) = token.strip_prefix("0x") {
        f64::from_bits(u64::from_str_radix(hex, 16).expect("f64 bits"))
    } else if let Some(hex) = token.strip_prefix("f32:") {
        f64::from(f32::from_bits(
            u32::from_str_radix(hex, 16).expect("f32 bits"),
        ))
    } else {
        match token {
            "nan" => f64::NAN,
            "inf" => f64::INFINITY,
            "-inf" => f64::NEG_INFINITY,
            other => other.parse().expect("decimal"),
        }
    }
}

/// One case's points, exactly as the C++ probe read them.
fn points(name: &str) -> Vec<(f64, f32)> {
    let path = format!("{DATA}/cases/{name}.tsv");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    text.lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let mut it = line.split('\t');
            let x = parse_scalar(it.next().expect("position"));
            let y = parse_scalar(it.next().expect("intensity")) as f32;
            (x, y)
        })
        .collect()
}

fn chromatogram(name: &str) -> MSChromatogram {
    MSChromatogram {
        peaks: points(name)
            .into_iter()
            .map(|(x, y)| ChromatogramPeak::new(x, y))
            .collect(),
        ..Default::default()
    }
}

fn spectrum(name: &str) -> MSSpectrum {
    MSSpectrum {
        peaks: points(name)
            .into_iter()
            .map(|(x, y)| Peak1D::new(x, y))
            .collect(),
        ..Default::default()
    }
}

/// One `case` row of `snt_oracle.tsv` with its measured ratios.
struct SntCase {
    kind: String,
    data: String,
    win_len: f64,
    bin_count: usize,
    write_log: bool,
    status: String,
    ratios: Vec<u64>,
}

fn snt_oracle() -> BTreeMap<String, SntCase> {
    let path = format!("{DATA}/snt_oracle.tsv");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut cases: BTreeMap<String, SntCase> = BTreeMap::new();
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "case" => {
                cases.insert(
                    f[1].to_string(),
                    SntCase {
                        kind: f[2].to_string(),
                        data: f[3].to_string(),
                        win_len: parse_scalar(f[4]),
                        bin_count: f[5].parse().expect("bin_count"),
                        write_log: f[6] == "true",
                        status: f[7].to_string(),
                        ratios: Vec::new(),
                    },
                );
            }
            "sn" => cases
                .get_mut(f[1])
                .expect("sn before case")
                .ratios
                .push(u64::from_str_radix(&f[3][2..], 16).expect("f64 bits")),
            other => panic!("unknown snt_oracle row {other}"),
        }
    }
    assert!(!cases.is_empty());
    cases
}

/// One `case` row of `pick_oracle.tsv` with its measured output.
struct PickCase {
    which: String,
    data: String,
    params: BTreeMap<String, String>,
    status: String,
    arrays: BTreeMap<String, Vec<u32>>,
    out: Vec<(u64, u32)>,
}
impl PickCase {
    fn number(&self, key: &str, fallback: f64) -> f64 {
        self.params.get(key).map_or(fallback, |v| parse_scalar(v))
    }
    fn flag(&self, key: &str, fallback: bool) -> bool {
        self.params.get(key).map_or(fallback, |v| v == "true")
    }
    fn array(&self, name: &str) -> &[u32] {
        self.arrays
            .get(name)
            .unwrap_or_else(|| panic!("no {name} array"))
    }
}

fn pick_oracle() -> BTreeMap<String, PickCase> {
    let path = format!("{DATA}/pick_oracle.tsv");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut cases: BTreeMap<String, PickCase> = BTreeMap::new();
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "case" => {
                let mut params = BTreeMap::new();
                if f[4] != "-" {
                    for entry in f[4].split(' ').filter(|e| !e.is_empty()) {
                        let (k, v) = entry.split_once('=').expect("key=value");
                        params.insert(k.to_string(), v.to_string());
                    }
                }
                cases.insert(
                    f[1].to_string(),
                    PickCase {
                        which: f[2].to_string(),
                        data: f[3].to_string(),
                        params,
                        status: f[5].to_string(),
                        arrays: BTreeMap::new(),
                        out: Vec::new(),
                    },
                );
            }
            "smooth" => {}
            "out" => cases.get_mut(f[1]).expect("out before case").out.push((
                u64::from_str_radix(&f[3][2..], 16).expect("f64 bits"),
                u32::from_str_radix(&f[4][2..], 16).expect("f32 bits"),
            )),
            "arr" => cases
                .get_mut(f[1])
                .expect("arr before case")
                .arrays
                .entry(f[2].to_string())
                .or_default()
                .push(u32::from_str_radix(&f[4][2..], 16).expect("f32 bits")),
            other => panic!("unknown pick_oracle row {other}"),
        }
    }
    assert!(!cases.is_empty());
    cases
}

/// The estimator each consumer builds, with that consumer's own parameters.
fn consumer_estimator(
    win_len: f64,
    bin_count: usize,
    write_log: bool,
) -> SignalToNoiseEstimatorMedian {
    SignalToNoiseEstimatorMedian {
        window_length: win_len,
        bin_count,
        write_log_messages: write_log,
        ..Default::default()
    }
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}
fn f32_bits(values: &[f32]) -> Vec<u32> {
    values.iter().map(|v| v.to_bits()).collect()
}
fn array<'a>(arrays: &'a [openms::kernel::DataArray<f32>], name: &str) -> &'a [f32] {
    &arrays
        .iter()
        .find(|a| a.name == name)
        .unwrap_or_else(|| panic!("no {name} array"))
        .data
}

// ------------------------------------------------- the estimator both pickers build

/// Every configuration either consumer reaches `init()` with, over inputs the
/// native profile refuses: duplicate, unsorted, negative and non-finite
/// samples, and the `win_len`/`bin_count` values the source's own parameter
/// restrictions accept or refuse.
#[test]
fn consumer_noise_estimates_match_the_release_build() {
    let cases = snt_oracle();
    let mut checked = 0;
    for (name, case) in &cases {
        let estimator = consumer_estimator(case.win_len, case.bin_count, case.write_log);
        // Both consumers hand the estimator a plain source object.
        let profile = PickingCompatibility::source();
        let result = match case.kind.as_str() {
            "chrom" => estimator.estimate_chromatogram(&chromatogram(&case.data), &profile),
            "spec" => estimator.estimate_spectrum(&spectrum(&case.data), &profile),
            other => panic!("unknown kind {other}"),
        };
        match case.status.as_str() {
            "ok" => {
                let estimates = result
                    .unwrap_or_else(|e| panic!("{name}: the Release build computed this: {e:?}"));
                assert_eq!(
                    bits(&estimates.signal_to_noise),
                    case.ratios,
                    "{name}: signal-to-noise bits differ from the Release build"
                );
            }
            "InvalidParameter" => {
                // `Param::checkDefaults` rejects the value before estimation:
                // `win_len` below its `setMinFloat(1.0)` (which `-inf` also
                // fails while NaN and `+inf` pass) and `bin_count` below its
                // `setMinInt(3)`.
                assert!(
                    matches!(result, Err(Error::InvalidValue(_))),
                    "{name}: the Release build threw Exception::InvalidParameter, got {result:?}"
                );
            }
            other => panic!("unexpected status {other} for {name}"),
        }
        checked += 1;
    }
    assert_eq!(checked, cases.len());
    assert!(checked >= 36, "the fixture lost cases: {checked}");
}

/// The estimator reads each `f32` intensity through the widening the Release
/// build's `cvtss2sd` performs, so a NaN keeps its sign and payload and becomes
/// quiet.
///
/// Reading the peaks is what makes that definite. Rust does not specify which
/// NaN a `f32`-to-`f64` widening produces, and the compiler may fold one at
/// build time, so `f64::from` on a copied slice happens to agree on the two
/// hosts this was checked on (x86_64 and aarch64) without being required to.
/// `estimate_peaks` goes through `x86::widen`, which is defined for every
/// payload, so these bit patterns are a contract rather than a coincidence.
#[test]
fn non_finite_intensities_widen_as_the_release_build_widens_them() {
    let case = &snt_oracle()["ppi_nonfinite"];
    let estimator = consumer_estimator(case.win_len, case.bin_count, case.write_log);
    let estimates = estimator
        .estimate_spectrum(&spectrum(&case.data), &PickingCompatibility::source())
        .expect("the source value domain accepts non-finite samples");
    assert_eq!(bits(&estimates.signal_to_noise), case.ratios);
    // The three inputs whose payload the widening has to carry: a signaling
    // NaN with payload 1, a quiet NaN with payload 1, and a negative default
    // quiet NaN. `(bits & 0x003f_ffff) << 29`, quieted, keeps the first two
    // apart from the default NaN the third produces.
    assert_eq!(
        estimates.signal_to_noise[7].to_bits(),
        0x7ff8_0000_2000_0000
    );
    assert_eq!(
        estimates.signal_to_noise[11].to_bits(),
        0x7ff8_0000_2000_0000
    );
    assert_eq!(
        estimates.signal_to_noise[13].to_bits(),
        0xfff8_0000_0000_0000
    );
    assert_eq!(
        estimates.signal_to_noise[17].to_bits(),
        0x7ff0_0000_0000_0000
    );
    assert_eq!(
        estimates.signal_to_noise[19].to_bits(),
        0xfff0_0000_0000_0000
    );
}

/// The native profile refuses exactly the samples and parameters the source
/// computes with, so the two profiles have to disagree on these cases.
#[test]
fn the_native_profile_refuses_what_the_source_computes() {
    for name in [
        "ppc_neg",
        "ppi_neg",
        "ppc_dup",
        "ppi_unsorted",
        "ppi_nonfinite",
    ] {
        let case = &snt_oracle()[name];
        assert_eq!(case.status, "ok", "{name} must be an accepted source case");
        let estimator = consumer_estimator(case.win_len, case.bin_count, case.write_log);
        let native = PickingCompatibility::default();
        let refused = match case.kind.as_str() {
            "chrom" => estimator.estimate_chromatogram(&chromatogram(&case.data), &native),
            _ => estimator.estimate_spectrum(&spectrum(&case.data), &native),
        };
        assert!(refused.is_err(), "{name}: the native profile must refuse");
    }
    // A NaN or infinite window passes `setMinFloat("win_len", 1.0)` because
    // `NaN < 1` is false and the restriction sets no upper bound.
    for name in ["ppc_winnan", "ppc_wininf", "ppi_winnan", "ppi_wininf"] {
        let case = &snt_oracle()[name];
        assert_eq!(case.status, "ok");
        let estimator = consumer_estimator(case.win_len, case.bin_count, case.write_log);
        assert!(
            estimator
                .estimate_spectrum(&spectrum(&case.data), &PickingCompatibility::default())
                .is_err(),
            "{name}: the native profile must refuse a non-finite window"
        );
    }
}

// ------------------------------------------------------- PeakPickerChromatogram

fn chromatogram_picker(case: &PickCase) -> PeakPickerChromatogram {
    let peak_width = case.number("peak_width", -1.0);
    PeakPickerChromatogram {
        method: if case.params.get("method").map(String::as_str) == Some("legacy") {
            ChromatogramPickingMethod::Legacy
        } else {
            ChromatogramPickingMethod::Corrected
        },
        smoothing: if case.flag("use_gauss", true) {
            ChromatogramSmoothing::Gaussian {
                width: case.number("gauss_width", 50.0),
            }
        } else {
            ChromatogramSmoothing::SavitzkyGolay {
                frame_length: case.number("sgolay_frame_length", 15.0) as usize,
                polynomial_order: case.number("sgolay_polynomial_order", 3.0) as usize,
            }
        },
        peak_width: (peak_width > 0.0).then_some(peak_width),
        signal_to_noise: case.number("signal_to_noise", 1.0),
        noise_estimator: consumer_estimator(
            case.number("sn_win_len", 1000.0),
            case.number("sn_bin_count", 30.0) as usize,
            false,
        ),
        report_sn: case.flag("report_sn", false),
        remove_overlapping_peaks: case.flag("remove_overlapping_peaks", false),
        ..Default::default()
    }
}

/// The boundary signal of the legacy method is the caller's chromatogram, so a
/// baseline-subtracted one puts negative intensities straight into
/// `snt_.init`. The source picks it; the native profile refuses it and names
/// the flag that lifts the refusal.
#[test]
fn a_negative_chromatogram_is_refused_natively_and_reproduced_in_source_mode() {
    let cases = pick_oracle();
    for name in ["ppc_legacy_neg", "ppc_legacy_negbase"] {
        let case = &cases[name];
        assert_eq!((case.which.as_str(), case.status.as_str()), ("ppc", "ok"));
        let mut picker = chromatogram_picker(case);
        let input = chromatogram(&case.data);

        let refused = picker
            .pick_chromatogram(&input)
            .expect_err("native refusal");
        match &refused {
            Error::InvalidValue(message) => assert!(
                message.contains("allow_negative_intensities"),
                "{name}: the refusal must name the flag, got {message}"
            ),
            other => panic!("{name}: unexpected {other:?}"),
        }

        picker.compatibility = PickingCompatibility::source();
        let picked = picker
            .pick_chromatogram(&input)
            .unwrap_or_else(|e| panic!("{name}: the Release build picked this: {e:?}"));
        let arrays = &picked.picked.chromatogram.float_data_arrays;
        assert_eq!(
            f32_bits(array(arrays, "SN")),
            case.array("SN"),
            "{name}: apex signal-to-noise differs from the Release build"
        );
        assert_eq!(
            f32_bits(array(arrays, "leftWidth")),
            case.array("leftWidth")
        );
        assert_eq!(
            f32_bits(array(arrays, "rightWidth")),
            case.array("rightWidth")
        );
        assert_eq!(
            f32_bits(array(arrays, "IntegratedIntensity")),
            case.array("IntegratedIntensity")
        );
    }
}

/// `MSChromatogram::isSorted` accepts equal retention times, so the source
/// reaches both the seed picker and `snt_.init` with duplicates.
#[test]
fn duplicate_retention_times_are_refused_natively_and_reproduced_in_source_mode() {
    let cases = pick_oracle();
    let case = &cases["ppc_gauss_dup"];
    assert_eq!(case.status, "ok");
    let mut picker = chromatogram_picker(case);
    let input = chromatogram(&case.data);

    match picker
        .pick_chromatogram(&input)
        .expect_err("native refusal")
    {
        Error::InvalidValue(message) => assert!(
            message.contains("allow_duplicate_positions"),
            "the refusal must name the flag, got {message}"
        ),
        other => panic!("unexpected {other:?}"),
    }

    picker.compatibility = PickingCompatibility::source();
    let picked = picker.pick_chromatogram(&input).expect("source mode picks");
    let arrays = &picked.picked.chromatogram.float_data_arrays;
    assert_eq!(f32_bits(array(arrays, "SN")), case.array("SN"));
    assert_eq!(
        f32_bits(array(arrays, "leftWidth")),
        case.array("leftWidth")
    );
    assert_eq!(
        f32_bits(array(arrays, "rightWidth")),
        case.array("rightWidth")
    );
}

/// `pickChromatogram` throws `Exception::IllegalArgument` for a chromatogram
/// that is not sorted by position, so `allow_unsorted_positions` must not lift
/// this port's refusal: the source refuses here too.
#[test]
fn an_unsorted_chromatogram_stays_refused_in_both_profiles() {
    let case = &pick_oracle()["ppc_sg_unsorted"];
    assert_eq!(case.status, "IllegalArgument");
    let input = chromatogram(&case.data);
    for profile in [
        PickingCompatibility::default(),
        PickingCompatibility::source(),
    ] {
        let picker = PeakPickerChromatogram {
            compatibility: profile,
            ..chromatogram_picker(case)
        };
        assert!(matches!(
            picker.pick_chromatogram(&input),
            Err(Error::UnsortedData)
        ));
    }
}

/// The chromatogram picker's internal estimator is a plain source object, not
/// one that follows the picker's profile: the caller's chromatogram has already
/// been validated above and, under `corrected`, the signal the estimator reads
/// is the smoothed one this picker produced. A `win_len` of NaN or `+inf`
/// passes the source's parameter restriction and the Release build picks with
/// it, while the native median-noise profile refuses both.
#[test]
fn the_chromatogram_estimator_always_uses_the_source_profile() {
    let cases = pick_oracle();
    for name in ["ppc_winnan_pick", "ppc_wininf_pick"] {
        let case = &cases[name];
        assert_eq!(case.status, "ok");
        let picker = chromatogram_picker(case);
        assert_eq!(picker.compatibility, PickingCompatibility::default());
        assert!(
            !picker.noise_estimator.window_length.is_finite(),
            "{name}: the case must carry a non-finite window"
        );
        let picked = picker
            .pick_chromatogram(&chromatogram(&case.data))
            .unwrap_or_else(|e| panic!("{name}: the Release build picked this: {e:?}"));
        let arrays = &picked.picked.chromatogram.float_data_arrays;
        assert_eq!(
            f32_bits(array(arrays, "SN")),
            case.array("SN"),
            "{name}: apex signal-to-noise differs from the Release build"
        );
        assert_eq!(
            f32_bits(array(arrays, "leftWidth")),
            case.array("leftWidth")
        );
        assert_eq!(
            f32_bits(array(arrays, "rightWidth")),
            case.array("rightWidth")
        );
    }
}

/// The values the source's own restrictions refuse stay refused in both
/// profiles: `Param::checkDefaults` throws before `init` ever runs.
#[test]
fn chromatogram_window_and_bin_count_below_the_source_minimum_stay_refused() {
    let input = chromatogram("clean");
    for (window, bins) in [(0.5, 30usize), (f64::NEG_INFINITY, 30), (1000.0, 1)] {
        for profile in [
            PickingCompatibility::default(),
            PickingCompatibility::source(),
        ] {
            let picker = PeakPickerChromatogram {
                noise_estimator: consumer_estimator(window, bins, false),
                report_sn: true,
                compatibility: profile,
                ..Default::default()
            };
            assert!(
                picker.pick_chromatogram(&input).is_err(),
                "win_len {window} bin_count {bins} must be refused"
            );
        }
    }
}

// --------------------------------------------------------- PeakPickerIterative

fn iterative_picker(case: &PickCase) -> PeakPickerIterative {
    PeakPickerIterative {
        signal_to_noise: case.number("signal_to_noise_", 1.0),
        peak_width: case.number("peak_width", 0.0),
        spacing_difference: case.number("spacing_difference", 1.5),
        noise_estimator: consumer_estimator(
            case.number("sn_win_len_", 20.0),
            case.number("sn_bin_count_", 30.0) as usize,
            true,
        ),
        iterations: case.number("nr_iterations_", 5.0) as usize,
        check_width_internally: case.flag("check_width_internally", false),
        ..Default::default()
    }
}

/// `snt.init(input)` reads the caller's raw spectrum, so the iterative picker's
/// estimate follows this picker's own profile, as `PeakPickerHiRes` does.
#[test]
fn iterative_negative_intensities_are_refused_natively_and_picked_in_source_mode() {
    let case = &pick_oracle()["ppi_neg"];
    assert_eq!((case.which.as_str(), case.status.as_str()), ("ppi", "ok"));
    let mut picker = iterative_picker(case);
    let input = spectrum(&case.data);

    match picker.pick_spectrum(&input).expect_err("native refusal") {
        Error::InvalidValue(message) => assert!(
            message.contains("allow_negative_intensities"),
            "the refusal must name the flag, got {message}"
        ),
        other => panic!("unexpected {other:?}"),
    }

    picker.compatibility = PickingCompatibility::source();
    let picked = picker
        .pick_spectrum(&input)
        .unwrap_or_else(|e| panic!("the Release build picked this: {e:?}"));
    let peaks: Vec<(u64, u32)> = picked
        .picked
        .spectrum
        .peaks
        .iter()
        .map(|p| (p.mz.to_bits(), p.intensity.to_bits()))
        .collect();
    assert_eq!(
        peaks, case.out,
        "picked peaks differ from the Release build"
    );
    let arrays = &picked.picked.spectrum.float_data_arrays;
    assert_eq!(
        f32_bits(array(arrays, "IntegratedIntensity")),
        case.array("IntegratedIntensity")
    );
    assert_eq!(
        f32_bits(array(arrays, "leftWidth")),
        case.array("leftWidth")
    );
    assert_eq!(
        f32_bits(array(arrays, "rightWidth")),
        case.array("rightWidth")
    );
}

/// A NaN or infinite `win_len` reaches the estimator through this picker's own
/// profile, so the native profile refuses it and the source profile reproduces
/// the Release build.
#[test]
fn the_iterative_estimator_follows_the_picker_profile() {
    let cases = pick_oracle();
    for name in ["ppi_winnan", "ppi_wininf"] {
        let case = &cases[name];
        assert_eq!(case.status, "ok");
        let mut picker = iterative_picker(case);
        assert!(!picker.noise_estimator.window_length.is_finite());
        assert!(
            picker.pick_spectrum(&spectrum(&case.data)).is_err(),
            "{name}: the native profile must refuse a non-finite window"
        );

        picker.compatibility = PickingCompatibility::source();
        let picked = picker
            .pick_spectrum(&spectrum(&case.data))
            .unwrap_or_else(|e| panic!("{name}: the Release build picked this: {e:?}"));
        let peaks: Vec<(u64, u32)> = picked
            .picked
            .spectrum
            .peaks
            .iter()
            .map(|p| (p.mz.to_bits(), p.intensity.to_bits()))
            .collect();
        assert_eq!(peaks, case.out, "{name}: picked peaks differ");
        let arrays = &picked.picked.spectrum.float_data_arrays;
        assert_eq!(
            f32_bits(array(arrays, "IntegratedIntensity")),
            case.array("IntegratedIntensity")
        );
        assert_eq!(
            f32_bits(array(arrays, "leftWidth")),
            case.array("leftWidth")
        );
        assert_eq!(
            f32_bits(array(arrays, "rightWidth")),
            case.array("rightWidth")
        );
    }
}

/// `snt.setParameters` throws for these before `init` runs, so both profiles
/// refuse them.
#[test]
fn iterative_window_and_bin_count_below_the_source_minimum_stay_refused() {
    let input = spectrum("clean");
    for (window, bins) in [
        (0.5, 30usize),
        (f64::NEG_INFINITY, 30),
        (20.0, 1),
        (20.0, 2),
    ] {
        for profile in [
            PickingCompatibility::default(),
            PickingCompatibility::source(),
        ] {
            let picker = PeakPickerIterative {
                noise_estimator: consumer_estimator(window, bins, true),
                compatibility: profile,
                ..Default::default()
            };
            assert!(
                picker.pick_spectrum(&input).is_err(),
                "win_len {window} bin_count {bins} must be refused"
            );
        }
    }
}

/// The source reaches `snt.init` with duplicate and decreasing m/z, but
/// `pickRecenterPeaks_` then collects each peak's support in a
/// `std::map<double, double>` keyed by m/z and bounds it with
/// `begin()`/`rbegin()` and `std::fabs` spacings, which this port has not
/// ported. Both refusals therefore stand in both profiles rather than
/// returning different peaks under a flag that claims source behaviour.
#[test]
fn iterative_duplicate_and_unsorted_positions_stay_refused_in_both_profiles() {
    for (data, expected_unsorted) in [("dup", false), ("unsorted", true)] {
        let input = spectrum(data);
        for profile in [
            PickingCompatibility::default(),
            PickingCompatibility::source(),
        ] {
            let picker = PeakPickerIterative {
                compatibility: profile,
                ..Default::default()
            };
            let error = picker.pick_spectrum(&input).expect_err("refused");
            if expected_unsorted {
                assert!(matches!(error, Error::UnsortedData), "{data}: {error:?}");
            } else {
                assert!(
                    matches!(&error, Error::InvalidValue(m) if m.contains("distinct profile coordinates")),
                    "{data}: {error:?}"
                );
            }
        }
    }
    // The Release build accepts both, which is why these stay tracked.
    let cases = pick_oracle();
    assert_eq!(cases["ppi_dup"].status, "ok");
    assert_eq!(cases["ppi_unsorted"].status, "ok");
}

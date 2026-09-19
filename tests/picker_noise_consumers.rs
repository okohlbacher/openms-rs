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
use openms::processing::peak_picking::{
    NoiseHistogramRange, PickingCompatibility, SignalToNoiseEstimatorMedian,
};
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
    /// `-1` leaves the automatic standard-deviation range, which is the only
    /// range either picker can reach; a nonnegative value is the source's
    /// manual upper end with `auto_mode = -1`.
    max_intensity: i32,
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
                        max_intensity: f[7].parse().expect("max_intensity"),
                        status: f[8].to_string(),
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
    /// The smoothed trace the chromatogram picker produced, empty for `ppi`.
    smooth: Vec<(u64, u32)>,
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
                        smooth: Vec::new(),
                    },
                );
            }
            "smooth" => cases
                .get_mut(f[1])
                .expect("smooth before case")
                .smooth
                .push((
                    u64::from_str_radix(&f[3][2..], 16).expect("f64 bits"),
                    u32::from_str_radix(&f[4][2..], 16).expect("f32 bits"),
                )),
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

/// The estimator of one `snt_oracle.tsv` case, including the manual histogram
/// range that only the two `*_bigmax` cases carry.
fn snt_estimator(case: &SntCase) -> SignalToNoiseEstimatorMedian {
    let mut estimator = consumer_estimator(case.win_len, case.bin_count, case.write_log);
    if case.max_intensity >= 0 {
        estimator.histogram_range = NoiseHistogramRange::Manual {
            max_intensity: f64::from(case.max_intensity),
        };
        estimator.range_parameters.max_intensity = case.max_intensity;
    }
    estimator
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
        let estimator = snt_estimator(case);
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
    assert!(checked >= 42, "the fixture lost cases: {checked}");
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

/// Neither the source nor this picker checks the sign of a retention time, so
/// a chromatogram that starts before zero picks in both profiles. Recorded
/// because the iterative picker below does refuse the matching spectrum.
#[test]
fn negative_retention_times_pick_in_both_profiles() {
    let case = &pick_oracle()["ppc_gauss_negmz"];
    assert_eq!(case.status, "ok");
    for profile in [
        PickingCompatibility::default(),
        PickingCompatibility::source(),
    ] {
        let picker = PeakPickerChromatogram {
            compatibility: profile,
            ..chromatogram_picker(case)
        };
        let picked = picker
            .pick_chromatogram(&chromatogram(&case.data))
            .expect("negative retention times are accepted");
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
}

/// The Release build picks a spectrum with negative m/z, and `PeakPickerHiRes`
/// does not check the sign of a position either, but this picker refuses it in
/// both profiles and no `PickingCompatibility` flag covers positions. The
/// refusal is therefore pinned here together with what the source does, so the
/// gap stays visible rather than being mistaken for source behaviour.
#[test]
fn negative_mz_stays_refused_by_the_iterative_picker_although_the_source_picks_it() {
    let case = &pick_oracle()["ppi_negmz"];
    assert_eq!(case.status, "ok");
    assert_eq!(case.out.len(), 1, "the Release build returns one peak");
    // m/z 7.9908..., integrated 6723.68 over leftWidth 4.0 and rightWidth 12.0.
    assert_eq!(case.out[0], (0x401f_f751_a000_0000, 0x45d2_1d70));
    assert_eq!(case.array("leftWidth"), [0x4080_0000]);
    assert_eq!(case.array("rightWidth"), [0x4140_0000]);

    let input = spectrum(&case.data);
    for profile in [
        PickingCompatibility::default(),
        PickingCompatibility::source(),
    ] {
        let picker = PeakPickerIterative {
            compatibility: profile,
            ..iterative_picker(case)
        };
        assert!(
            matches!(picker.pick_spectrum(&input),
                Err(Error::InvalidValue(m)) if m.contains("nonnegative m/z")),
            "negative m/z must stay refused"
        );
    }
}

// ------------------------------------ every measured case, every measured column

/// Which profile a recorded refusal stands in. No case is refused in the source
/// profile and accepted natively — `PickingCompatibility::source` sets every
/// flag, so it accepts everything the native default accepts — and the replay
/// below fails if that ever stops holding.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum RefusedIn {
    /// Accepted once `compatibility` is `PickingCompatibility::source()`.
    NativeOnly,
    /// Refused in both profiles: a gap between this port and the source that no
    /// flag covers yet. Each of these is listed in the branch's `not_done`.
    BothProfiles,
}
use RefusedIn::{BothProfiles, NativeOnly};

/// Every `ok` case of `pick_oracle.tsv` this port refuses, the profile the
/// refusal stands in, and the text the message must carry.
///
/// The replay fails if a case outside this table is refused, if a case in it is
/// accepted where the table says it is refused, or if an entry here is never
/// reached — so the table cannot drift away from what the port does.
const REFUSED: &[(&str, RefusedIn, &str)] = &[
    // Negative intensities: no source check exists, so these are native-only.
    ("ppc_gauss_neg", NativeOnly, "allow_negative_intensities"),
    ("ppc_gauss_neg_ov", NativeOnly, "allow_negative_intensities"),
    ("ppc_gauss_neg_pw", NativeOnly, "allow_negative_intensities"),
    ("ppc_gauss2_neg", NativeOnly, "allow_negative_intensities"),
    ("ppc_gauss8_neg", NativeOnly, "allow_negative_intensities"),
    ("ppc_gauss16_neg", NativeOnly, "allow_negative_intensities"),
    ("ppc_sg_neg", NativeOnly, "allow_negative_intensities"),
    ("ppc_legacy_neg", NativeOnly, "allow_negative_intensities"),
    (
        "ppc_winnan_legacy",
        NativeOnly,
        "allow_negative_intensities",
    ),
    (
        "ppc_gauss_negbase",
        NativeOnly,
        "allow_negative_intensities",
    ),
    (
        "ppc_gauss_negbase_sn0",
        NativeOnly,
        "allow_negative_intensities",
    ),
    (
        "ppc_legacy_negbase",
        NativeOnly,
        "allow_negative_intensities",
    ),
    ("ppc_sg_negbase", NativeOnly, "allow_negative_intensities"),
    ("ppc_gauss_allneg", NativeOnly, "allow_negative_intensities"),
    (
        "ppc_legacy_allneg",
        NativeOnly,
        "allow_negative_intensities",
    ),
    ("ppi_neg", NativeOnly, "allow_negative_intensities"),
    ("ppi_negbase", NativeOnly, "allow_negative_intensities"),
    ("ppi_allneg", NativeOnly, "allow_negative_intensities"),
    ("ppi_neg_sn0", NativeOnly, "allow_negative_intensities"),
    ("ppi_negbase_sn0", NativeOnly, "allow_negative_intensities"),
    ("ppi_allneg_sn0", NativeOnly, "allow_negative_intensities"),
    // Duplicate retention times, which `MSChromatogram::isSorted` lets through.
    ("ppc_sg_dup", NativeOnly, "allow_duplicate_positions"),
    ("ppc_gauss_dup", NativeOnly, "allow_duplicate_positions"),
    // Both at once; the intensity check runs first.
    ("ppc_gauss_dupneg", NativeOnly, "allow_negative_intensities"),
    // A `win_len` the source's `setMinFloat(1.0)` lets through, refused by the
    // native median-noise profile. Only the iterative picker's estimate follows
    // the picker's profile, which is why no `ppc_win*` case appears here.
    ("ppi_winnan", NativeOnly, "invalid median-noise options"),
    ("ppi_wininf", NativeOnly, "invalid median-noise options"),
    // --- refused in both profiles: recorded gaps, all in `not_done` ---
    // `pickRecenterPeaks_` keys each peak's support on `std::map<double,double>`;
    // see `PeakPickerIterative::compatibility`.
    ("ppi_dup", BothProfiles, "distinct profile coordinates"),
    ("ppi_unsorted", BothProfiles, "UnsortedData"),
    // No `PickingCompatibility` flag covers the sign of a position.
    ("ppi_negmz", BothProfiles, "nonnegative m/z"),
    // The kernel validators refuse a non-finite sample before either picker
    // runs, and `PeakPickerHiRes::validate_points` refuses one again behind
    // them; no flag lifts either.
    ("ppi_nfint", BothProfiles, "peak intensity must be finite"),
    (
        "ppi_nonfinite",
        BothProfiles,
        "peak intensity must be finite",
    ),
    (
        "ppc_gauss_nfint",
        BothProfiles,
        "chromatogram intensity must be finite",
    ),
];

fn refusal(name: &str, native: bool) -> Option<&'static str> {
    REFUSED
        .iter()
        .find(|(n, r, _)| *n == name && (native || *r == BothProfiles))
        .map(|(_, _, message)| *message)
}

/// Replay every case the Release build picked, in both profiles, asserting
/// every column the fixture carries: the picked peaks, the smoothed trace, and
/// each of `IntegratedIntensity`, `leftWidth`, `rightWidth` and `SN`.
///
/// The per-case tests above pin the individual behaviours and their refusal
/// messages; this one exists so that no measured case sits in the fixture
/// unasserted, and so that a case this port refuses has to be a listed gap.
#[test]
fn every_release_picked_case_matches_in_every_profile_that_accepts_it() {
    let cases = pick_oracle();
    let mut reached: Vec<&str> = Vec::new();
    let mut accepted = 0usize;
    let mut refused = 0usize;
    for (name, case) in &cases {
        if case.status != "ok" {
            continue;
        }
        for native in [true, false] {
            let profile = if native {
                PickingCompatibility::default()
            } else {
                PickingCompatibility::source()
            };
            let tag = if native { "native" } else { "source" };
            let expected = refusal(name, native);
            let (peaks, arrays, smoothed) = match case.which.as_str() {
                "ppc" => {
                    let picker = PeakPickerChromatogram {
                        compatibility: profile,
                        ..chromatogram_picker(case)
                    };
                    match picker.pick_chromatogram(&chromatogram(&case.data)) {
                        Ok(picked) => {
                            let peaks = picked
                                .picked
                                .chromatogram
                                .peaks
                                .iter()
                                .map(|p| (p.rt.to_bits(), p.intensity.to_bits()))
                                .collect::<Vec<_>>();
                            let smoothed = picked
                                .smoothed
                                .peaks
                                .iter()
                                .map(|p| (p.rt.to_bits(), p.intensity.to_bits()))
                                .collect::<Vec<_>>();
                            (
                                peaks,
                                picked.picked.chromatogram.float_data_arrays.clone(),
                                Some(smoothed),
                            )
                        }
                        Err(error) => {
                            let message = expected.unwrap_or_else(|| {
                                panic!("{name}/{tag}: the Release build picked this, got {error:?}")
                            });
                            assert!(
                                format!("{error:?}").contains(message),
                                "{name}/{tag}: refusal must mention {message}, got {error:?}"
                            );
                            reached.push(name);
                            refused += 1;
                            continue;
                        }
                    }
                }
                "ppi" => {
                    let picker = PeakPickerIterative {
                        compatibility: profile,
                        ..iterative_picker(case)
                    };
                    match picker.pick_spectrum(&spectrum(&case.data)) {
                        Ok(picked) => {
                            let peaks = picked
                                .picked
                                .spectrum
                                .peaks
                                .iter()
                                .map(|p| (p.mz.to_bits(), p.intensity.to_bits()))
                                .collect::<Vec<_>>();
                            (
                                peaks,
                                picked.picked.spectrum.float_data_arrays.clone(),
                                None,
                            )
                        }
                        Err(error) => {
                            let message = expected.unwrap_or_else(|| {
                                panic!("{name}/{tag}: the Release build picked this, got {error:?}")
                            });
                            assert!(
                                format!("{error:?}").contains(message),
                                "{name}/{tag}: refusal must mention {message}, got {error:?}"
                            );
                            reached.push(name);
                            refused += 1;
                            continue;
                        }
                    }
                }
                other => panic!("unknown consumer {other}"),
            };
            assert!(
                expected.is_none(),
                "{name}/{tag}: listed as refused ({}), but the port accepted it",
                expected.unwrap_or_default()
            );
            assert_eq!(
                peaks, case.out,
                "{name}/{tag}: picked peaks differ from the Release build"
            );
            if let Some(smoothed) = smoothed {
                assert_eq!(
                    smoothed, case.smooth,
                    "{name}/{tag}: the smoothed trace differs from the Release build"
                );
            }
            // Six measured cases pick nothing at all — the Release build runs
            // them without throwing and returns an empty record. The port then
            // has to return the three (or four) named arrays empty, not absent.
            if case.out.is_empty() {
                assert!(
                    case.arrays.is_empty(),
                    "{name}: no peaks but arrays in the fixture"
                );
                for data in &arrays {
                    assert!(
                        data.data.is_empty(),
                        "{name}/{tag}: {} must be empty when no peak is picked",
                        data.name
                    );
                }
            } else {
                assert!(
                    !case.arrays.is_empty(),
                    "{name}: the fixture carries no output array"
                );
            }
            for (array_name, expected_bits) in &case.arrays {
                assert_eq!(
                    &f32_bits(array(&arrays, array_name)),
                    expected_bits,
                    "{name}/{tag}: {array_name} differs from the Release build"
                );
            }
            accepted += 1;
        }
    }
    for (name, _, _) in REFUSED {
        assert!(
            reached.contains(name),
            "{name} is listed as refused but the replay never refused it"
        );
    }
    let picked_by_source = cases.values().filter(|c| c.status == "ok").count();
    assert_eq!(
        accepted + refused,
        2 * picked_by_source,
        "every case the Release build accepted must be replayed in both profiles"
    );
    assert!(
        accepted >= 56 && refused >= 38,
        "the fixture lost cases: {accepted} accepted runs, {refused} refusals"
    );
}

/// Every case the Release build itself threw on, and the refusal this port
/// answers with. Two of them are refused a step earlier than the source
/// throws, which is the point of recording the text rather than the kind.
const EXCEPTIONS: &[(&str, &str)] = &[
    // `pickChromatogram` (`PeakPickerChromatogram.cpp:68-72`) throws
    // `IllegalArgument`, "Chromatogram must be sorted by position", because
    // `MSChromatogram::isSorted` compares `prev.getRT() > next.getRT()` and the
    // `+inf` retention time at index 29 is greater than the finite one after
    // it. Here `MSChromatogram::validate` refuses the same chromatogram one
    // step earlier, for its non-finite intensities.
    ("ppc_sg_nonfinite", "chromatogram intensity must be finite"),
    ("ppc_sg_unsorted", "UnsortedData"),
    // `Param::checkDefaults` throws before `init` runs: `setMinInt("bin_count",
    // 3)` and `setMinFloat("win_len", 1.0)`.
    ("ppi_bins1", "invalid median-noise options"),
    ("ppi_winsmall", "invalid median-noise options"),
    // `signal_to_noise_` carries no restriction of its own
    // (`PeakPickerIterative.h:92`) and is copied into `PeakPickerHiRes`'s
    // `signal_to_noise`, which carries `setMinFloat(0.0)`
    // (`PeakPickerHiRes.cpp:31-32`), so the Release build aborts naming a class
    // the caller never mentioned. This port refuses the option itself. CPP-348.
    ("ppi_snneg", "invalid iterative picker options"),
];

/// The cases the Release build itself refused stay refused here, in both
/// profiles, with the refusal `EXCEPTIONS` records.
#[test]
fn every_release_exception_case_stays_refused_in_both_profiles() {
    let cases = pick_oracle();
    let mut checked = 0usize;
    for (name, case) in &cases {
        if case.status == "ok" {
            continue;
        }
        for profile in [
            PickingCompatibility::default(),
            PickingCompatibility::source(),
        ] {
            let error = match case.which.as_str() {
                "ppc" => PeakPickerChromatogram {
                    compatibility: profile,
                    ..chromatogram_picker(case)
                }
                .pick_chromatogram(&chromatogram(&case.data))
                .err(),
                _ => PeakPickerIterative {
                    compatibility: profile,
                    ..iterative_picker(case)
                }
                .pick_spectrum(&spectrum(&case.data))
                .err(),
            };
            let error =
                error.unwrap_or_else(|| panic!("{name}: the Release build threw, port accepted"));
            let (_, expected) = EXCEPTIONS
                .iter()
                .find(|(n, _)| *n == name)
                .unwrap_or_else(|| {
                    panic!("{name}: no recorded refusal for a case the source threw on")
                });
            assert!(
                format!("{error:?}").contains(expected),
                "{name}: refusal must mention {expected}, got {error:?}"
            );
            checked += 1;
        }
    }
    assert_eq!(
        checked,
        2 * EXCEPTIONS.len(),
        "every recorded exception case must be replayed in both profiles"
    );
}

/// `signal_to_noise_ = 0.0` is what `PeakPickerIterative.h:92` documents as the
/// way to turn the S/N gate off, and `:185`, `:205` and `:315` honour it by
/// skipping `snt.init` and every S/N break. A candidate's support then extends
/// across the negative samples, `weighted_mz /= integrated_intensity` (`:228`)
/// divides by a negative sum, and `:231-234` stores both the quotient and that
/// negative sum. The Release build does all of this; `allow_negative_intensities`
/// is what makes the branch reachable, so the refusal it used to carry belongs
/// to the native profile alone.
#[test]
fn a_negative_integrated_intensity_divides_as_the_source_divides_by_it() {
    let cases = pick_oracle();
    for name in ["ppi_neg_sn0", "ppi_negbase_sn0", "ppi_allneg_sn0"] {
        let case = &cases[name];
        assert_eq!((case.which.as_str(), case.status.as_str()), ("ppi", "ok"));
        assert_eq!(
            case.number("signal_to_noise_", 1.0),
            0.0,
            "{name}: the case must turn the S/N gate off"
        );
        // Without at least one negative sum this case would not reach the
        // branch, and the test below would pass for the wrong reason.
        let negative = case
            .array("IntegratedIntensity")
            .iter()
            .filter(|bits| *bits & 0x8000_0000 != 0)
            .count();
        assert!(
            negative > 0,
            "{name}: the Release build stored no negative integrated intensity"
        );

        let mut picker = iterative_picker(case);
        match picker.pick_spectrum(&spectrum(&case.data)) {
            Err(Error::InvalidValue(message)) => assert!(
                message.contains("allow_negative_intensities"),
                "{name}: the native refusal must name the flag, got {message}"
            ),
            other => panic!("{name}: the native profile must refuse, got {other:?}"),
        }

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
            case.array("IntegratedIntensity"),
            "{name}: integrated intensities differ from the Release build"
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

/// Non-finite samples are refused by the kernel validators, before either
/// picker's own checks and before `PeakPickerHiRes::validate_points` refuses
/// them again behind those. The Release build accepts all three — it picks
/// peaks from the two spectra and returns an empty chromatogram for the third —
/// so the gap is pinned here together with what it produced, as the negative-m/z
/// and duplicate gaps are.
#[test]
fn non_finite_samples_stay_refused_by_both_pickers_although_the_release_build_picks_them() {
    let cases = pick_oracle();

    // `nfint` and `nonfinite` carry the same five non-finite intensities; the
    // latter also carries a non-finite position, which the chromatogram side
    // rejects as unsorted first, so only the spectrum side reaches it here.
    for name in ["ppi_nfint", "ppi_nonfinite"] {
        let case = &cases[name];
        assert_eq!(case.status, "ok", "{name}: the Release build picked this");
        assert!(
            !case.out.is_empty(),
            "{name}: the Release build returned peaks"
        );
        for profile in [
            PickingCompatibility::default(),
            PickingCompatibility::source(),
        ] {
            let picker = PeakPickerIterative {
                compatibility: profile,
                ..iterative_picker(case)
            };
            assert!(
                matches!(picker.pick_spectrum(&spectrum(&case.data)),
                    Err(Error::InvalidValue(m)) if m == "peak intensity must be finite"),
                "{name}: a non-finite intensity must stay refused"
            );
        }
    }

    let case = &cases["ppc_gauss_nfint"];
    // The Release build runs this one without throwing and returns an empty
    // chromatogram: the non-finite intensities propagate through the Gaussian
    // smoother, so no seed survives. Recorded so the gap is measured rather
    // than asserted to exist.
    assert_eq!(case.status, "ok");
    assert!(
        case.out.is_empty() && case.arrays.is_empty(),
        "the Release build returned peaks for ppc_gauss_nfint"
    );
    for profile in [
        PickingCompatibility::default(),
        PickingCompatibility::source(),
    ] {
        let picker = PeakPickerChromatogram {
            compatibility: profile,
            ..chromatogram_picker(case)
        };
        assert!(
            matches!(picker.pick_chromatogram(&chromatogram(&case.data)),
                Err(Error::InvalidValue(m)) if m == "chromatogram intensity must be finite"),
            "a non-finite intensity must stay refused"
        );
    }
}

/// The chromatogram picker's estimate runs under `PickingCompatibility::source`
/// whatever the picker's own profile says, which also selects the Linux x86-64
/// Release build's bin-index conversion: a histogram quotient outside `int`
/// range becomes `INT_MIN` and lands in bin `0`, where clamping before
/// truncation would put it in the last bin.
///
/// Reaching that needs a hand-set histogram range, which the source's picker
/// never sets (`PeakPickerChromatogram.cpp:408-412` configures `win_len`,
/// `bin_count` and `write_log_messages` only). The Rust field exposes the whole
/// estimator, so the configuration exists here; the numbers it produces are the
/// Release build's own, measured at the estimator as `ppc_bigmax`.
#[test]
fn the_chromatogram_estimator_bins_out_of_range_quotients_as_the_release_build_does() {
    let snt = snt_oracle();
    let case = &snt["ppc_bigmax"];
    assert_eq!((case.kind.as_str(), case.status.as_str()), ("chrom", "ok"));
    assert_eq!(case.max_intensity, 10);
    let input = chromatogram(&case.data);
    assert!(
        input.peaks.iter().any(|p| f64::from(p.intensity)
            > f64::from(u32::MAX) * f64::from(case.max_intensity) / case.bin_count as f64),
        "the case must carry an intensity whose quotient leaves int range"
    );

    // What the Release build's estimator answers over exactly these samples.
    let estimator = snt_estimator(case);
    let source = estimator
        .estimate_chromatogram(&input, &PickingCompatibility::source())
        .expect("the source profile estimates this");
    assert_eq!(
        bits(&source.signal_to_noise),
        case.ratios,
        "the source profile must reproduce the Release build"
    );
    // The native profile accepts every sample of this chromatogram — finite,
    // nonnegative, strictly increasing retention times — and answers
    // differently only because it clamps the quotient before truncating it.
    let native = estimator
        .estimate_chromatogram(&input, &PickingCompatibility::default())
        .expect("the native profile accepts this chromatogram in full");
    assert_ne!(
        bits(&native.signal_to_noise),
        case.ratios,
        "the two bin-index conversions must disagree here, or this pins nothing"
    );

    // The picker at its default profile reports the source values. `legacy`
    // makes the boundary signal the caller's chromatogram, so the estimate runs
    // over exactly the samples the fixture was measured over.
    let picker = PeakPickerChromatogram {
        method: ChromatogramPickingMethod::Legacy,
        noise_estimator: estimator,
        report_sn: true,
        ..Default::default()
    };
    assert_eq!(picker.compatibility, PickingCompatibility::default());
    let picked = picker
        .pick_chromatogram(&input)
        .expect("the chromatogram is accepted at the default profile");
    let reported = f32_bits(array(&picked.picked.chromatogram.float_data_arrays, "SN"));
    assert!(!reported.is_empty(), "the picker reported no apex S/N");
    let from_source: Vec<u32> = source
        .signal_to_noise
        .iter()
        .map(|v| (*v as f32).to_bits())
        .collect();
    let from_native: Vec<u32> = native
        .signal_to_noise
        .iter()
        .map(|v| (*v as f32).to_bits())
        .collect();
    for value in &reported {
        assert!(
            from_source.contains(value),
            "reported apex S/N {value:#010x} is not a Release-measured ratio"
        );
    }
    assert!(
        reported.iter().any(|v| !from_native.contains(v)),
        "the reported ratios are also native ones, so this does not discriminate"
    );
}

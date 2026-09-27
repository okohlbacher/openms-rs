// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! DEBT-04: the inputs the port used to refuse because "the order is
//! unspecified", against what the Linux x86_64 Release build of OpenMS
//! `bc9cc12` does with them (tier 1, executed differential).
//!
//! Lead decision D16 puts reproducing libstdc++'s own `std::sort`,
//! `std::stable_sort`, `std::is_sorted`, `std::lower_bound` and red-black-tree
//! `lower_bound` results in scope. Every expected value below is read from
//! `tests/data/a3a_source_order/oracle.tsv`, produced by
//! `../oracle/A3a-source-order/a3a_driver.cpp` on ibminode06 against the
//! Release install `openms4-release-bc9cc12-c19e494-174b576` (two invocations
//! of two fresh-process runs each, byte-identical); the inputs are
//! `tests/data/a3a_source_order/cases.tsv`, shared byte for byte with the
//! driver. `tests/data/a3a_source_order_provenance.json` records the hashes.
//!
//! Doubles are compared by their bits, floats likewise: a NaN's sign and
//! payload included, since every NaN here is either an input or made by
//! `crate::math::x86_64`'s x86_64 rules.

use openms::analysis::feature_finder_picked::debug::AbortReasons;
use openms::analysis::feature_finder_picked::helper_structs::Seed;
use openms::kernel::faims_helper::FaimsHelper;
use openms::kernel::range_utils::{IsInIsolationWindow, SpectrumPredicate};
use openms::kernel::{
    ChromatogramPeak, ChromatogramTools, DataArray, MSChromatogram, MSExperiment, MSSpectrum,
    Peak1D, Precursor,
};
use openms::kernel::{Feature, FeatureMap};
use openms::math::statistic_functions::{compute_rank, rank_correlation_coefficient};
use openms::metadata::{ChromatogramType, DriftTimeUnit, Product, ScanMode};
use openms::processing::feature_overlap_filter::{
    CentroidTolerances, FeatureOverlapFilter, FeatureOverlapMode,
};
use openms::processing::iterative::PeakPickerIterative;
use openms::processing::peak_picking::PickingCompatibility;
use std::collections::BTreeMap;

const CASES: &str = include_str!("data/a3a_source_order/cases.tsv");
const ORACLE: &str = include_str!("data/a3a_source_order/oracle.tsv");

/// A case: its mode and its fields.
struct Case {
    mode: String,
    fields: Vec<String>,
}

fn cases() -> BTreeMap<String, Case> {
    let mut out = BTreeMap::new();
    for line in CASES.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t').map(str::to_owned);
        let name = fields.next().expect("name");
        let mode = fields.next().expect("mode");
        out.insert(
            name,
            Case {
                mode,
                fields: fields.collect(),
            },
        );
    }
    out
}

/// The records of one case: key and values, in the driver's order.
fn oracle(name: &str) -> Vec<(String, Vec<String>)> {
    ORACLE
        .lines()
        .filter_map(|line| {
            let mut fields = line.split('\t');
            (fields.next() == Some(name)).then(|| {
                let key = fields.next().expect("key").to_owned();
                (key, fields.map(str::to_owned).collect())
            })
        })
        .collect()
}

/// The values of the only record `key` of a case.
fn record(name: &str, key: &str) -> Vec<String> {
    let found: Vec<Vec<String>> = oracle(name)
        .into_iter()
        .filter(|(k, _)| k == key)
        .map(|(_, v)| v)
        .collect();
    assert_eq!(found.len(), 1, "{name}: one '{key}' record");
    found.into_iter().next().expect("one")
}

/// Every case of `mode` must have run to completion in the Release build.
fn names(mode: &str) -> Vec<String> {
    let names: Vec<String> = cases()
        .into_iter()
        .filter(|(_, case)| case.mode == mode)
        .map(|(name, _)| name)
        .collect();
    assert!(!names.is_empty(), "no {mode} cases");
    for name in &names {
        assert_eq!(record(name, "status"), ["ok"], "{name}");
        assert_eq!(record(name, "exit"), ["0"], "{name}");
    }
    names
}

fn list(field: &str) -> Vec<&str> {
    if field == "-" || field.is_empty() {
        Vec::new()
    } else {
        field.split(',').collect()
    }
}

fn f64_token(token: &str) -> f64 {
    let hex = token.strip_prefix("0x").expect("0x token");
    f64::from_bits(u64::from_str_radix(hex, 16).expect("hex"))
}

fn f32_token(token: &str) -> f32 {
    let hex = token.strip_prefix("f32:").expect("f32: token");
    f32::from_bits(u32::from_str_radix(hex, 16).expect("hex"))
}

fn doubles(field: &str) -> Vec<f64> {
    list(field).into_iter().map(f64_token).collect()
}

fn floats(field: &str) -> Vec<f32> {
    list(field).into_iter().map(f32_token).collect()
}

fn ints(field: &str) -> Vec<i64> {
    list(field)
        .into_iter()
        .map(|token| token.parse().expect("integer"))
        .collect()
}

fn h64(value: f64) -> String {
    format!("0x{:016x}", value.to_bits())
}

fn h32(value: f32) -> String {
    format!("f32:{:08x}", value.to_bits())
}

fn h64s(values: impl IntoIterator<Item = f64>) -> Vec<String> {
    values.into_iter().map(h64).collect()
}

fn h32s(values: impl IntoIterator<Item = f32>) -> Vec<String> {
    values.into_iter().map(h32).collect()
}

// --- MATH/StatisticFunctions.h: computeRank and rankCorrelationCoefficient ---

/// `computeRank` (`StatisticFunctions.h:815-863`) sorts `(index, value)`
/// pairs with a lambda on `.second` (`:829-830`). With a NaN among the values
/// the introsort's permutation decides which original index gets which rank,
/// and this is that row: the ranks of NaN-bearing vectors of 3 to 64 values,
/// ties, signed zeros and infinities included.
#[test]
fn compute_rank_reproduces_the_release_build_on_nan_keys() {
    let all = cases();
    let mut with_nan = 0;
    for name in names("rank") {
        let mut values = doubles(&all[&name].fields[0]);
        with_nan += usize::from(values.iter().any(|v| v.is_nan()));
        compute_rank(&mut values).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!(h64s(values), record(&name, "rank"), "{name}");
    }
    assert!(with_nan >= 10, "the NaN rows are the point of this test");
}

#[test]
fn rank_correlation_reproduces_the_release_build_on_nan_keys() {
    let all = cases();
    for name in names("rankcorr") {
        let a = doubles(&all[&name].fields[0]);
        let b = doubles(&all[&name].fields[1]);
        let rho = rank_correlation_coefficient(&a, &b).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!([h64(rho)], record(&name, "rho").as_slice(), "{name}");
    }
}

// --- KERNEL/RangeUtils.h: IsInIsolationWindow ---

#[test]
fn isolation_window_sorts_and_searches_as_the_release_build() {
    let all = cases();
    for name in names("isolation") {
        let case = &all[&name];
        let mz = doubles(&case.fields[0]);
        let reverse = ints(&case.fields[1])[0] != 0;
        let probes = doubles(&case.fields[2]);
        let predicate =
            IsInIsolationWindow::new(mz, reverse).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        // The driver replays the constructor's `std::sort` of a
        // `std::vector<double>` (`RangeUtils.h:634`) on a copy.
        assert_eq!(
            h64s(predicate.mz_values().iter().copied()),
            record(&name, "sorted"),
            "{name}"
        );
        let keep: Vec<String> = probes
            .chunks_exact(3)
            .map(|probe| {
                let mut precursor = Precursor::new(probe[0], 0);
                precursor.isolation_window_lower_offset = probe[1];
                precursor.isolation_window_upper_offset = probe[2];
                let spectrum = MSSpectrum {
                    ms_level: 2,
                    precursors: vec![precursor],
                    ..MSSpectrum::default()
                };
                if predicate.keep(&spectrum) { "1" } else { "0" }.to_owned()
            })
            .collect();
        assert_eq!(keep, record(&name, "keep"), "{name}");
    }
}

// --- KERNEL/MSSpectrum: isSortedByIM and sortByIonMobility ---

#[test]
fn ion_mobility_sort_reproduces_the_release_build() {
    let all = cases();
    for name in names("mobility") {
        let case = &all[&name];
        let mz = doubles(&case.fields[0]);
        let intensity = floats(&case.fields[1]);
        let im = floats(&case.fields[2]);
        let extra = floats(&case.fields[3]);
        let mut spectrum = MSSpectrum {
            peaks: mz
                .iter()
                .zip(&intensity)
                .map(|(&m, &i)| Peak1D::new(m, i))
                .collect(),
            float_data_arrays: vec![DataArray::new("mean ion mobility array", im)],
            ..MSSpectrum::default()
        };
        if !extra.is_empty() {
            spectrum
                .float_data_arrays
                .push(DataArray::new("extra", extra));
        }
        let sorted = spectrum
            .is_sorted_by_im()
            .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!(
            [if sorted { "1" } else { "0" }],
            record(&name, "is_sorted_by_im").as_slice(),
            "{name}"
        );
        spectrum
            .sort_by_ion_mobility()
            .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!(
            h64s(spectrum.peaks.iter().map(|p| p.mz)),
            record(&name, "sorted_mz"),
            "{name}"
        );
        assert_eq!(
            h32s(spectrum.peaks.iter().map(|p| p.intensity)),
            record(&name, "sorted_int"),
            "{name}"
        );
        let arrays: Vec<Vec<String>> = oracle(&name)
            .into_iter()
            .filter(|(key, _)| key == "sorted_farr")
            .map(|(_, values)| values)
            .collect();
        assert_eq!(arrays.len(), spectrum.float_data_arrays.len(), "{name}");
        for (array, expected) in spectrum.float_data_arrays.iter().zip(&arrays) {
            let mut got = vec![array.name.clone()];
            got.extend(h32s(array.data.iter().copied()));
            assert_eq!(&got, expected, "{name}");
        }
        let after = spectrum.is_sorted_by_im().unwrap();
        assert_eq!(
            [if after { "1" } else { "0" }],
            record(&name, "is_sorted_by_im_after").as_slice(),
            "{name}"
        );
    }
}

// --- IONMOBILITY/FAIMSHelper: getCompensationVoltages ---

#[test]
fn faims_voltages_reproduce_the_release_build_on_nan() {
    let all = cases();
    for name in names("faims") {
        let case = &all[&name];
        let drift = doubles(&case.fields[0]);
        let units = ints(&case.fields[1]);
        let experiment = MSExperiment {
            spectra: drift
                .iter()
                .zip(&units)
                .map(|(&drift_time, &unit)| MSSpectrum {
                    drift_time,
                    drift_time_unit: if unit != 0 {
                        DriftTimeUnit::FaimsCompensationVoltage
                    } else {
                        DriftTimeUnit::None
                    },
                    ..MSSpectrum::default()
                })
                .collect(),
            ..MSExperiment::default()
        };
        let cvs = FaimsHelper::get_compensation_voltages(&experiment)
            .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!(h64s(cvs.values()), record(&name, "cvs"), "{name}");
        let warned_in_cpp = oracle(&name).iter().any(|(key, values)| {
            key == "stderr"
                && values
                    .iter()
                    .any(|v| v == FaimsHelper::MISSING_VOLTAGE_WARNING)
        });
        assert_eq!(!cvs.warnings.is_empty(), warned_in_cpp, "{name}");
    }
}

// --- FEATUREFINDER: abort_reasons_[seed] = reason ---

#[test]
fn abort_reason_map_reproduces_std_map_on_nan_intensities() {
    let all = cases();
    for name in names("abortmap") {
        let case = &all[&name];
        let intensity = floats(&case.fields[0]);
        let spectrum = ints(&case.fields[1]);
        let peak = ints(&case.fields[2]);
        let reason = ints(&case.fields[3]);
        let mut map = AbortReasons::new();
        for i in 0..intensity.len() {
            let seed = Seed::new(spectrum[i] as usize, peak[i] as usize, intensity[i]);
            map.insert(seed, &format!("r{}", reason[i]))
                .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        }
        let got: Vec<Vec<String>> = map
            .iter()
            .map(|(seed, reason)| {
                vec![
                    seed.spectrum.to_string(),
                    seed.peak.to_string(),
                    h32(seed.intensity),
                    reason.to_owned(),
                ]
            })
            .collect();
        let expected: Vec<Vec<String>> = oracle(&name)
            .into_iter()
            .filter(|(key, _)| key == "entry")
            .map(|(_, values)| values)
            .collect();
        assert_eq!(got, expected, "{name}");
        assert_eq!(map.len(), expected.len(), "{name}");
    }
}

// --- KERNEL/ChromatogramTools ---

fn chromatogram_type_code(kind: ChromatogramType) -> String {
    match kind {
        ChromatogramType::Mass => "0",
        ChromatogramType::SelectedIonMonitoring => "4",
        ChromatogramType::SelectedReactionMonitoring => "5",
        other => panic!("unexpected {other:?}"),
    }
    .to_owned()
}

fn scan_mode_code(mode: ScanMode) -> String {
    match mode {
        ScanMode::Unknown => "0",
        ScanMode::SelectedIonMonitoring => "4",
        ScanMode::SelectedReactionMonitoring => "5",
        other => panic!("unexpected {other:?}"),
    }
    .to_owned()
}

#[test]
fn spectra_to_chromatograms_group_nan_keys_as_the_release_build() {
    let all = cases();
    for name in names("s2c") {
        let f = &all[&name].fields;
        let srm = ints(&f[0]);
        let rt = doubles(&f[1]);
        let precursor_counts = ints(&f[2]);
        let precursor_mz = doubles(&f[3]);
        let peak_counts = ints(&f[4]);
        let peak_mz = doubles(&f[5]);
        let peak_intensity = floats(&f[6]);
        let flags = ints(&f[7]);
        let mut experiment = MSExperiment::default();
        let (mut next_precursor, mut next_peak) = (0, 0);
        for i in 0..srm.len() {
            let mut spectrum = MSSpectrum {
                rt: rt[i],
                native_id: format!("s{i}"),
                ..MSSpectrum::default()
            };
            if srm[i] != 0 {
                spectrum.instrument_settings.scan_mode = ScanMode::SelectedReactionMonitoring;
            }
            for _ in 0..precursor_counts[i] {
                spectrum
                    .precursors
                    .push(Precursor::new(precursor_mz[next_precursor], 0));
                next_precursor += 1;
            }
            for _ in 0..peak_counts[i] {
                spectrum
                    .peaks
                    .push(Peak1D::new(peak_mz[next_peak], peak_intensity[next_peak]));
                next_peak += 1;
            }
            experiment.spectra.push(spectrum);
        }
        ChromatogramTools::default()
            .convert_spectra_to_chromatograms(&mut experiment, flags[0] != 0, flags[1] != 0)
            .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let mut got = Vec::new();
        for chromatogram in &experiment.chromatograms {
            let description = chromatogram
                .precursor
                .cv_terms
                .metadata
                .get("description")
                .map_or_else(
                    || "-".to_owned(),
                    |value| value.as_str().expect("text").to_owned(),
                );
            let native_id = if chromatogram.native_id.is_empty() {
                "-".to_owned()
            } else {
                chromatogram.native_id.clone()
            };
            got.push((
                "chrom".to_owned(),
                vec![
                    chromatogram_type_code(chromatogram.chromatogram_type),
                    h64(chromatogram.precursor.mz),
                    h64(chromatogram.product.mz),
                    native_id,
                    description,
                ],
            ));
            got.push((
                "chrom_rt".to_owned(),
                h64s(chromatogram.peaks.iter().map(|p| p.rt)),
            ));
            got.push((
                "chrom_int".to_owned(),
                h32s(chromatogram.peaks.iter().map(|p| p.intensity)),
            ));
        }
        got.push((
            "spectra".to_owned(),
            experiment
                .spectra
                .iter()
                .map(|s| s.native_id.clone())
                .collect(),
        ));
        let expected: Vec<(String, Vec<String>)> = oracle(&name)
            .into_iter()
            .filter(|(key, _)| key.starts_with("chrom") || key == "spectra")
            .collect();
        assert_eq!(got, expected, "{name}");
    }
}

#[test]
fn chromatograms_to_spectra_copy_nonfinite_values_as_the_release_build() {
    let all = cases();
    for name in names("c2s") {
        let f = &all[&name].fields;
        let kinds = list(&f[0]);
        let product = doubles(&f[1]);
        let precursor = doubles(&f[2]);
        let counts = ints(&f[3]);
        let rt = doubles(&f[4]);
        let intensity = floats(&f[5]);
        let mut experiment = MSExperiment::default();
        let mut next = 0;
        for i in 0..kinds.len() {
            let mut chromatogram = MSChromatogram {
                chromatogram_type: match kinds[i] {
                    "srm" => ChromatogramType::SelectedReactionMonitoring,
                    "sim" => ChromatogramType::SelectedIonMonitoring,
                    _ => ChromatogramType::Mass,
                },
                product: Product {
                    mz: product[i],
                    ..Default::default()
                },
                precursor: Precursor::new(precursor[i], 0),
                ..MSChromatogram::default()
            };
            for _ in 0..counts[i] {
                chromatogram
                    .peaks
                    .push(ChromatogramPeak::new(rt[next], intensity[next]));
                next += 1;
            }
            experiment.chromatograms.push(chromatogram);
        }
        ChromatogramTools::default()
            .convert_chromatograms_to_spectra(&mut experiment)
            .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let mut got: Vec<(String, Vec<String>)> = experiment
            .spectra
            .iter()
            .map(|s| {
                (
                    "spectrum".to_owned(),
                    vec![
                        h64(s.rt),
                        s.ms_level.to_string(),
                        scan_mode_code(s.instrument_settings.scan_mode),
                        h64(s.precursors[0].mz),
                        h64(s.products[0].mz),
                        h64(s.peaks[0].mz),
                        h32(s.peaks[0].intensity),
                    ],
                )
            })
            .collect();
        got.push((
            "chromatograms_left".to_owned(),
            vec![experiment.chromatograms.len().to_string()],
        ));
        let expected: Vec<(String, Vec<String>)> = oracle(&name)
            .into_iter()
            .filter(|(key, _)| key == "spectrum" || key == "chromatograms_left")
            .collect();
        assert_eq!(got, expected, "{name}");
    }
}

// --- PROCESSING/FEATURE/FeatureOverlapFilter: std::stable_sort ---

/// The driver filters with the source's default comparator
/// `left.getOverallQuality() > right.getOverallQuality()` on features far
/// apart, so its output order is the `std::stable_sort` order of
/// `FeatureOverlapFilter.cpp:141`. A `Feature` here cannot hold a NaN quality
/// (`Feature::validate`), so the port's features carry a finite quality and
/// the comparator reads the case's qualities by unique ID: the same
/// comparisons, with the same outcomes, in the order the sort makes them.
#[test]
fn feature_overlap_filter_sorts_as_the_release_builds_stable_sort() {
    let all = cases();
    for name in names("fof") {
        let f = &all[&name].fields;
        let quality = doubles(&f[0]);
        let rt = doubles(&f[1]);
        let mz = doubles(&f[2]);
        let features: Vec<Feature> = (0..quality.len())
            .map(|i| {
                let mut feature = Feature::new(rt[i], mz[i], 1.0);
                feature.base.charge = 2;
                feature.base.unique_id = i as u64 + 1;
                feature
            })
            .collect();
        let mut map = FeatureMap::from_features(features);
        // The source's QualityType is `float`: its setter narrows each value.
        let key = |feature: &Feature| quality[feature.base.unique_id as usize - 1] as f32;
        FeatureOverlapFilter::filter_with_mode(
            &mut map,
            |left, right| key(left) > key(right),
            |_, _| true,
            FeatureOverlapMode::CentroidBased,
            &CentroidTolerances::default(),
        )
        .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let order: Vec<String> = map
            .features
            .iter()
            .map(|feature| feature.base.unique_id.to_string())
            .collect();
        assert_eq!(order, record(&name, "order"), "{name}");
    }
}

// --- PROCESSING/CENTROIDING/PeakPickerIterative ---

/// `PeakPickerIterative::pick` in the source profile: the introsort of the
/// seeds by descending apex intensity (`PeakPickerIterative.h:311`), zero
/// integrated intensities that divide to an infinite or NaN centroid, and the
/// final `sortByPosition` (`:366`) over such centroids. The `ppi_tie_*` rows
/// hold more than 16 equally intense seeds whose centroids all come out equal
/// or NaN, so the output order, read through `leftWidth`, is the order the
/// introsort left them in.
#[test]
fn iterative_picker_reproduces_the_release_build_order_and_zero_sums() {
    let all = cases();
    let mut nonfinite = 0;
    for name in names("ppi") {
        let f = &all[&name].fields;
        let mz = doubles(&f[0]);
        let intensity = floats(&f[1]);
        let mut picker = PeakPickerIterative {
            compatibility: PickingCompatibility::source(),
            ..PeakPickerIterative::default()
        };
        for pair in list(&f[2].replace(';', ",")) {
            let (key, value) = pair.split_once('=').expect("key=value");
            let value: f64 = value.parse().expect("number");
            match key {
                "signal_to_noise_" => picker.signal_to_noise = value,
                "peak_width" => picker.peak_width = value,
                other => panic!("{name}: unexpected parameter {other}"),
            }
        }
        let input = MSSpectrum {
            peaks: mz
                .iter()
                .zip(&intensity)
                .map(|(&m, &i)| Peak1D::new(m, i))
                .collect(),
            ..MSSpectrum::default()
        };
        let picked = picker
            .pick_spectrum(&input)
            .unwrap_or_else(|e| panic!("{name}: the Release build picked this: {e:?}"))
            .picked
            .spectrum;
        nonfinite += picked.peaks.iter().filter(|p| !p.mz.is_finite()).count();
        assert_eq!(
            h64s(picked.peaks.iter().map(|p| p.mz)),
            record(&name, "out_mz"),
            "{name}"
        );
        assert_eq!(
            h32s(picked.peaks.iter().map(|p| p.intensity)),
            record(&name, "out_int"),
            "{name}"
        );
        let arrays: Vec<Vec<String>> = oracle(&name)
            .into_iter()
            .filter(|(key, _)| key == "out_farr")
            .map(|(_, values)| values)
            .collect();
        assert_eq!(arrays.len(), picked.float_data_arrays.len(), "{name}");
        for (array, expected) in picked.float_data_arrays.iter().zip(&arrays) {
            let mut got = vec![array.name.clone()];
            got.extend(h32s(array.data.iter().copied()));
            assert_eq!(&got, expected, "{name}");
        }
    }
    assert!(
        nonfinite >= 60,
        "the zero-sum rows are the point of this test"
    );
}

/// The `ppi_tie_*` rows discriminate: more than 16 seeds share one apex
/// intensity, and the Release build's output is not in their stable order.
#[test]
fn the_tie_rows_carry_more_than_sixteen_equal_seeds_out_of_stable_order() {
    for name in ["ppi_tie_ninf_20", "ppi_tie_pinf_20"] {
        let seeds = record(name, "seeds_int");
        let tied = seeds.iter().filter(|s| **s == seeds[1]).count();
        assert!(tied > 16, "{name}: {tied} tied seeds");
        let widths: Vec<Vec<String>> = oracle(name)
            .into_iter()
            .filter(|(key, values)| key == "out_farr" && values[0] == "leftWidth")
            .map(|(_, values)| values[1..].to_vec())
            .collect();
        let left: Vec<f32> = widths[0].iter().map(|t| f32_token(t)).collect();
        assert!(
            left.windows(2).any(|pair| pair[0] > pair[1]),
            "{name}: the equal centroids keep the introsort's order, not the input order"
        );
    }
}

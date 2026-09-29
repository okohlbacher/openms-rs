// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `TheoreticalSpectrumGenerator` configured through its `DefaultParamHandler`
//! contract, against the Linux x86-64 Release build.
//!
//! Tier 1: `tests/data/param_defaults/tsg_intensity/<case>.tsv` is the stdout
//! of `../oracle/B6-tsg-intensity` (`b6_tsg <case>`) on ibminode06 against the
//! Release install `openms4-release-bc9cc12-c19e494-174b576`, run twice in
//! fresh processes with identical output: a `peaks<TAB>n` line, then one line
//! per peak with the m/z (`%.17g`), the `float` intensity's bits, the ion name
//! and the charge. The source's intensity parameters are `double` members
//! (`TheoreticalSpectrumGenerator.h:141-153`) multiplied in double precision
//! and narrowed once into `Peak1D`; `b_intensity` 0.7 and
//! `relative_loss_intensity` 0.9 are not representable in `f32`, so these
//! cases fail for a port that holds them as `f32`.

use openms::chemistry::{AASequence, TheoreticalSpectrumGenerator};
use openms::kernel::MSSpectrum;
use openms::param::{DefaultParamHandling, Param, ParamValue};
use std::path::PathBuf;

fn fixture(name: &str) -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/param_defaults/tsg_intensity")
        .join(name);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The C++ peaks: m/z, intensity bits, ion name, charge.
fn oracle(case: &str) -> Vec<(f64, u32, String, i32)> {
    let text = fixture(&format!("{case}.tsv"));
    let mut lines = text.lines();
    let count: usize = lines
        .next()
        .and_then(|line| line.strip_prefix("peaks\t"))
        .unwrap()
        .parse()
        .unwrap();
    let peaks: Vec<_> = lines
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), 4, "{line}");
            (
                fields[0].parse().unwrap(),
                u32::from_str_radix(fields[1], 16).unwrap(),
                fields[2].to_owned(),
                fields[3].parse().unwrap(),
            )
        })
        .collect();
    assert_eq!(peaks.len(), count, "{case}");
    peaks
}

/// A generator after `setParameters` with `entries` over the defaults, as the
/// driver builds it.
fn generator(entries: &[(&str, ParamValue)]) -> TheoreticalSpectrumGenerator {
    let mut param = TheoreticalSpectrumGenerator::defaults().unwrap();
    for (key, value) in entries {
        param.replace_value(key, value.clone()).unwrap();
    }
    let (generator, warnings) =
        TheoreticalSpectrumGenerator::from_param_with_warnings(&param).unwrap();
    assert!(warnings.is_empty(), "{warnings:?}");
    generator
}

fn yes() -> ParamValue {
    ParamValue::String("true".into())
}

fn assert_spectrum(case: &str, spectrum: &MSSpectrum) {
    let expected = oracle(case);
    let names = &spectrum
        .string_data_arrays
        .iter()
        .find(|array| array.name == "IonNames")
        .unwrap()
        .data;
    let charges = &spectrum
        .integer_data_arrays
        .iter()
        .find(|array| array.name == "Charges")
        .unwrap()
        .data;
    assert_eq!(spectrum.peaks.len(), expected.len(), "{case}");
    let mut worst_ulps = 0;
    for (i, (mz, bits, name, charge)) in expected.iter().enumerate() {
        let peak = &spectrum.peaks[i];
        assert_eq!(
            (names[i].as_str(), charges[i]),
            (name.as_str(), *charge),
            "{case} peak {i}"
        );
        // The subject of these cases: bit-identical intensities.
        assert_eq!(
            peak.intensity.to_bits(),
            *bits,
            "{case} peak {i} {name}: {} != {}",
            peak.intensity,
            f32::from_bits(*bits)
        );
        // The port's m/z arithmetic sums masses in its own order and differs
        // from the Release build in the last bits (reported in the lane
        // report); this is the 1e-9 the TheoreticalSpectrumGenerator
        // reference tests already hold m/z to.
        assert!(
            (peak.mz - mz).abs() <= 1e-9,
            "{case} peak {i} {name}: {} != {mz}",
            peak.mz
        );
        worst_ulps = worst_ulps.max(peak.mz.to_bits().abs_diff(mz.to_bits()));
    }
    eprintln!("{case}: largest m/z difference {worst_ulps} ulp");
}

#[test]
fn fine_envelopes_scale_by_the_double_ion_intensity() {
    let g = generator(&[
        ("isotope_model", ParamValue::String("fine".into())),
        ("max_isotope_probability", ParamValue::Float(0.0001)),
        ("add_first_prefix_ion", yes()),
        ("add_metainfo", yes()),
        ("add_y_ions", ParamValue::String("false".into())),
        ("b_intensity", ParamValue::Float(0.7)),
    ]);
    let peptide: AASequence = "AA".parse().unwrap();
    assert_spectrum("fine_AA", &g.generate(&peptide, 1, 2, Some(4)).unwrap());
}

#[test]
fn fine_losses_multiply_both_doubles_before_one_narrowing() {
    let g = generator(&[
        ("isotope_model", ParamValue::String("fine".into())),
        ("max_isotope_probability", ParamValue::Float(0.005)),
        ("add_first_prefix_ion", yes()),
        ("add_losses", yes()),
        ("add_precursor_peaks", yes()),
        ("sort_by_position", ParamValue::String("false".into())),
        ("add_metainfo", yes()),
        ("b_intensity", ParamValue::Float(0.7)),
        ("relative_loss_intensity", ParamValue::Float(0.9)),
    ]);
    let peptide: AASequence = ".(Acetyl)SG.(Amidated)".parse().unwrap();
    assert_spectrum("fine_SG", &g.generate(&peptide, 1, 1, None).unwrap());
}

#[test]
fn plain_losses_narrow_the_double_product() {
    let g = generator(&[
        ("add_losses", yes()),
        ("add_metainfo", yes()),
        ("add_a_ions", yes()),
        ("b_intensity", ParamValue::Float(0.7)),
        ("y_intensity", ParamValue::Float(0.3)),
        ("a_intensity", ParamValue::Float(0.1)),
        ("relative_loss_intensity", ParamValue::Float(0.9)),
    ]);
    let peptide: AASequence = "PEPTIDERK".parse().unwrap();
    assert_spectrum("plain_losses", &g.generate(&peptide, 1, 2, None).unwrap());
}

/// The parameter tree the generator reports is the one it was built from,
/// with the values as given (no `f32` rounding on the way through).
#[test]
fn configured_intensities_round_trip_exactly() {
    let g = generator(&[
        ("b_intensity", ParamValue::Float(0.7)),
        ("relative_loss_intensity", ParamValue::Float(0.9)),
        ("precursor_H2O_intensity", ParamValue::Float(0.3)),
    ]);
    let param: Param = g.to_param().unwrap();
    for (key, value) in [
        ("b_intensity", 0.7),
        ("relative_loss_intensity", 0.9),
        ("precursor_H2O_intensity", 0.3),
    ] {
        assert_eq!(
            param.value(key).unwrap(),
            &ParamValue::Float(value),
            "{key}"
        );
    }
}

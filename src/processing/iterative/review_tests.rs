// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Independently derived source oracles; no upstream numerical tests exist.

use super::*;

fn fixture(name: &str, text: &str) -> MSSpectrum {
    MSSpectrum::from_peaks(
        text.lines()
            .skip(1)
            .filter_map(|line| {
                let fields: Vec<_> = line.split('\t').collect();
                (fields[0] == name).then(|| {
                    Peak1D::new(
                        fields[2].parse().unwrap(),
                        f32::from_bits(u32::from_str_radix(fields[3], 16).unwrap()),
                    )
                })
            })
            .collect(),
    )
}

#[test]
fn independent_seed_recenter_extension_and_suppression_oracles() {
    for row in include_str!("../../../tests/data/iterative_options.tsv")
        .lines()
        .skip(1)
    {
        let config: Vec<_> = row.split('\t').collect();
        let name = config[0];
        let input = fixture(
            name,
            include_str!("../../../tests/data/iterative_traces.tsv"),
        );
        let seeds = fixture(
            name,
            include_str!("../../../tests/data/iterative_seeds.tsv"),
        );
        let picker = PeakPickerIterative {
            iterations: config[1].parse().unwrap(),
            peak_width: config[2].parse().unwrap(),
            spacing_difference: config[3].parse().unwrap(),
            check_width_internally: config[4].parse().unwrap(),
            signal_to_noise: config[5].parse().unwrap(),
            ..Default::default()
        };
        let expected: Vec<Vec<_>> = include_str!("../../../tests/data/iterative_expected.tsv")
            .lines()
            .skip(1)
            .map(|l| l.split('\t').collect::<Vec<_>>())
            .filter(|v| v[0] == name)
            .collect();
        let result = picker
            .refine_with_seeds(&input, &seeds)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let output = &result.picked.spectrum;
        assert_eq!(output.len(), expected.len(), "{name}");
        assert_eq!(result.regions.len(), expected.len());
        assert_eq!(result.picked.boundaries.len(), expected.len());
        for (i, expected) in expected.iter().enumerate() {
            assert_eq!(expected[1].parse::<usize>().unwrap(), i);
            let region = &result.regions[i];
            assert_eq!(
                region.seed_index,
                expected[2].parse::<usize>().unwrap(),
                "{name}"
            );
            assert_eq!(
                region.initial_center_index,
                expected[3].parse::<usize>().unwrap(),
                "{name}"
            );
            assert_eq!(
                region.center_index,
                expected[4].parse::<usize>().unwrap(),
                "{name}"
            );
            let left = expected[5].parse::<usize>().unwrap();
            let right = expected[6].parse::<usize>().unwrap();
            assert_eq!(
                (region.left_index, region.right_index),
                (left, right),
                "{name}"
            );
            let mz = f32::from_bits(u32::from_str_radix(expected[8], 16).unwrap());
            let intensity = f32::from_bits(u32::from_str_radix(expected[10], 16).unwrap());
            assert_eq!(
                output.peaks[i],
                Peak1D::new(f64::from(mz), intensity),
                "{name}"
            );
            assert_eq!(
                result.picked.boundaries[i],
                PeakBoundary {
                    min: input.peaks[left].mz,
                    max: input.peaks[right].mz
                }
            );
            assert_eq!(output.float_data_arrays[0].data[i], intensity);
            assert_eq!(
                output.float_data_arrays[1].data[i],
                input.peaks[left].mz as f32
            );
            assert_eq!(
                output.float_data_arrays[2].data[i],
                input.peaks[right].mz as f32
            );
        }
        assert_eq!(
            output
                .float_data_arrays
                .iter()
                .map(|a| a.name.as_str())
                .collect::<Vec<_>>(),
            ["IntegratedIntensity", "leftWidth", "rightWidth"]
        );
        output.validate().unwrap();
    }
}

#[test]
fn equal_seed_priorities_use_the_documented_stable_native_order() {
    let input = fixture(
        "one_seed_per_raw_index",
        include_str!("../../../tests/data/iterative_traces.tsv"),
    );
    let mut seeds = fixture(
        "one_seed_per_raw_index",
        include_str!("../../../tests/data/iterative_seeds.tsv"),
    );
    for peak in &mut seeds.peaks {
        peak.intensity = 10.;
    }
    let picker = PeakPickerIterative {
        signal_to_noise: 0.,
        spacing_difference: 1.,
        iterations: 1,
        ..Default::default()
    };
    let result = picker.refine_with_seeds(&input, &seeds).unwrap();
    // Source std::sort does not specify equal-key order. Stable seed order is a
    // deterministic native policy, not an asserted upstream numerical result.
    assert_eq!(result.regions.len(), 1);
    assert_eq!(result.regions[0].seed_index, 0);
    assert_eq!(result.regions[0].initial_center_index, 3);
    assert_eq!(result.picked.spectrum.peaks, [Peak1D::new(103., 18.)]);
}

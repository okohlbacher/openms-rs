use super::*;

fn from_hex(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).unwrap())
}

#[test]
fn decisions_and_labels_match_executed_libsvm() {
    let mut count = 0;
    let mut labels = [[0; 2]; 2];
    for line in include_str!("../../../tests/data/metabo_predictor_libsvm.tsv").lines() {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 7);
        let (model, data, index) = match fields[0] {
            "2" => (Model::Noise2, &NOISE2, 0),
            "5" => (Model::Noise5, &NOISE5, 1),
            _ => panic!("invalid oracle model"),
        };
        let raw = std::array::from_fn(|i| from_hex(fields[i + 1]));
        let expected = from_hex(fields[5]);
        let mut work = 100_000;
        let actual = decision(data, raw, &mut work).unwrap();
        // libm implementations can differ in their final bits. The fixture
        // retains the exact observed binary64 values, with a portable bound.
        assert!(
            (actual - expected).abs() <= 1e-10 + expected.abs() * 1e-11,
            "model {} case {count}: {actual} != {expected}",
            fields[0]
        );
        let accepted = fields[6] == "2";
        assert_eq!(
            predict(model, raw, &mut work, &mut 0).unwrap(),
            accepted,
            "model {} case {count}",
            fields[0]
        );
        labels[index][usize::from(accepted)] += 1;
        count += 1;
    }
    assert_eq!(count, 488);
    for model in labels {
        assert!(model[0] > 20 && model[1] > 20);
    }
}

#[test]
fn all_projected_coefficients_features_and_scales_retain_source_bits() {
    let models = [
        (
            &NOISE2,
            include_str!(
                "../../../resources/metabolite_isotope_models/MetaboliteIsoModelNoised2.svm"
            ),
            include_str!(
                "../../../resources/metabolite_isotope_models/MetaboliteIsoModelNoised2.scale"
            ),
            548,
        ),
        (
            &NOISE5,
            include_str!(
                "../../../resources/metabolite_isotope_models/MetaboliteIsoModelNoised5.svm"
            ),
            include_str!(
                "../../../resources/metabolite_isotope_models/MetaboliteIsoModelNoised5.scale"
            ),
            999,
        ),
    ];
    for (model, source, scale, count) in models {
        assert_eq!(model.support.len(), count);
        let (header, rows) = source.split_once("\nSV\n").unwrap();
        for (key, expected) in [("gamma ", model.gamma), ("rho ", model.rho)] {
            let literal: f64 = header
                .lines()
                .find_map(|row| row.strip_prefix(key))
                .unwrap()
                .parse()
                .unwrap();
            assert_eq!(literal.to_bits(), expected.to_bits());
        }
        for (record, line) in model.support.iter().zip(rows.lines()) {
            let fields: Vec<_> = line.split_whitespace().collect();
            assert_eq!(fields.len(), 5);
            assert_eq!(
                fields[0].parse::<f64>().unwrap().to_bits(),
                record.coefficient.to_bits()
            );
            for (i, field) in fields[1..].iter().enumerate() {
                let (index, value) = field.split_once(':').unwrap();
                assert_eq!(index.parse::<usize>().unwrap(), i + 1);
                assert_eq!(
                    value.parse::<f64>().unwrap().to_bits(),
                    record.features[i].to_bits()
                );
            }
        }
        assert_eq!(rows.lines().count(), count);
        assert_eq!(scale.lines().count(), 4);
        for (i, row) in scale.lines().enumerate() {
            let values: Vec<f64> = row.split_whitespace().map(|v| v.parse().unwrap()).collect();
            assert_eq!(values[0].to_bits(), model.centers[i].to_bits());
            assert_eq!(values[1].to_bits(), model.scales[i].to_bits());
        }
    }
}

#[test]
fn shared_work_is_precharged_and_prediction_needs_no_heap_budget() {
    let costs = [(Model::Noise2, 17_552), (Model::Noise5, 31_984)];
    for (model, cost) in costs {
        let mut work = cost * 2;
        let mut bytes = 0;
        for expected in [cost, 0] {
            predict(model, [500., 0.3, 0.1, 0.], &mut work, &mut bytes).unwrap();
            assert_eq!(work, expected);
            assert_eq!(bytes, 0);
        }
        assert!(predict(model, [0.; 4], &mut work, &mut bytes).is_err());
        assert_eq!(work, 0);
        work = cost - 1;
        assert!(predict(model, [0.; 4], &mut work, &mut bytes).is_err());
        assert_eq!(work, cost - 1);
    }
}

#[test]
fn nonfinite_raw_values_fail_but_finite_rbf_overflow_keeps_source_limit() {
    for (kind, model) in [(Model::Noise2, &NOISE2), (Model::Noise5, &NOISE5)] {
        for position in 0..4 {
            for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let mut raw = [0.; 4];
                raw[position] = value;
                assert!(predict(kind, raw, &mut 100_000, &mut 0).is_err());
            }
        }
        for raw in [[f64::MAX; 4], [-f64::MAX; 4], [500., 1e300, 0., 0.]] {
            assert_eq!(decision(model, raw, &mut 100_000).unwrap(), -model.rho);
            assert!(!predict(kind, raw, &mut 100_000, &mut 0).unwrap());
        }
    }
}

#[test]
fn literal_row_order_and_strict_zero_boundary_are_retained() {
    static ROWS: [SupportVector; 3] = [
        SupportVector {
            coefficient: 1e16,
            features: [0.; 4],
        },
        SupportVector {
            coefficient: -1e16,
            features: [0.; 4],
        },
        SupportVector {
            coefficient: 1.,
            features: [0.; 4],
        },
    ];
    static REORDERED: [SupportVector; 3] = [
        SupportVector {
            coefficient: 1e16,
            features: [0.; 4],
        },
        SupportVector {
            coefficient: 1.,
            features: [0.; 4],
        },
        SupportVector {
            coefficient: -1e16,
            features: [0.; 4],
        },
    ];
    let mut model = ModelData {
        gamma: 2.,
        rho: 1.,
        centers: [0.; 4],
        scales: [1.; 4],
        support: &ROWS,
    };
    let value = decision(&model, [0.; 4], &mut 1000).unwrap();
    assert_eq!(value, 0.);
    assert!(value <= 0.); // Source zero votes for label 1, not acceptance label 2.
    model.rho = 0.5;
    assert_eq!(decision(&model, [0.; 4], &mut 1000).unwrap(), 0.5);
    model.support = &REORDERED;
    assert_eq!(decision(&model, [0.; 4], &mut 1000).unwrap(), -0.5);
}

#[test]
fn scaling_is_center_then_divide_including_absent_isotope_zeros() {
    static ROWS: [SupportVector; 1] = [SupportVector {
        coefficient: 4.,
        features: [0.5, -0.5, -0.5, -0.5],
    }];
    let model = ModelData {
        gamma: 1.,
        rho: 1.,
        centers: [2.; 4],
        scales: [4.; 4],
        support: &ROWS,
    };
    assert_eq!(decision(&model, [4., 0., 0., 0.], &mut 1000).unwrap(), 3.);
    let shifted = decision(&model, [8., 0., 0., 0.], &mut 1000).unwrap();
    assert_eq!(shifted, 4. * (-1f64).exp() - 1.);
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Every `START_SECTION` of `FAIMSHelper_test.cpp` at Core SDK `bc9cc12`, the
//! executed product-sdk oracle for both `FAIMSHelper` functions, and the native
//! boundaries of `src/kernel/faims_helper.rs`.
//!
//! Evidence labels:
//!
//! - Expectations citing a `FAIMSHelper_test.cpp` line are transcribed
//!   class-test literals (tier 3, source review).
//! - `tests/data/faims_helper/oracle_cases.tsv` is oracle-generated (tier 1
//!   executed differential): the unmodified product-sdk `libOpenMS` (Debug,
//!   core 4fdec46, identical to `bc9cc12` for every file involved) run by
//!   `../oracle/pte-faims-helper/driver.cpp`. The upstream class tests do not
//!   exercise `filterPeptidesByFAIMSCV` at all, so its only executed evidence
//!   is this oracle.
//! - The NaN refusal, the parameter refusals and the ceilings are native
//!   (tier 4). Where the source gives a result for input this port refuses,
//!   the test asserts the refusal *and* the recorded source result, so the
//!   divergence stays visible.

use openms::concept::constants::user_param::FAIMS_CV;
use openms::identification::PeptideIdentification;
use openms::kernel::faims_helper::{CompensationVoltage, FaimsHelper};
use openms::metadata::{DriftTimeUnit, ImTypes, MetaValue, MetaValueData, to_drift_time_unit};
use openms::{Error, MSExperiment, MSSpectrum};
use std::hash::{DefaultHasher, Hash, Hasher};

const ORACLE: &str = include_str!("data/faims_helper/oracle_cases.tsv");

/// The oracle records of one kind, split into tab-separated fields.
fn records(kind: &str) -> Vec<Vec<&'static str>> {
    ORACLE
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|fields| fields[0] == kind)
        .collect()
}

/// Parses the C `printf("%a")` spelling the oracle writes (`-0x1.9p+5`,
/// `-0x0p+0`, `inf`, `nan`). Every mantissa has at most 53 significant bits,
/// so the conversion is exact.
fn hex_f64(text: &str) -> f64 {
    let (negative, body) = match text.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, text),
    };
    let magnitude = match body {
        "inf" => f64::INFINITY,
        "nan" => f64::NAN,
        _ => {
            let body = body.strip_prefix("0x").expect("hexadecimal float");
            let (mantissa, exponent) = body.split_once('p').expect("binary exponent");
            let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
            let mut bits: u64 = 0;
            for digit in whole.chars().chain(fraction.chars()) {
                bits = bits * 16 + u64::from(digit.to_digit(16).expect("hexadecimal digit"));
            }
            let shift = exponent.parse::<i32>().expect("exponent") - 4 * fraction.len() as i32;
            let mut value = bits as f64;
            for _ in 0..shift.unsigned_abs() {
                value = if shift > 0 { value * 2.0 } else { value / 2.0 };
            }
            value
        }
    };
    if negative { -magnitude } else { magnitude }
}

/// The bit patterns of a `-`-or-comma-separated oracle voltage list, so that
/// the sign of zero is compared too.
fn oracle_bits(list: &str) -> Vec<u64> {
    if list == "-" {
        return Vec::new();
    }
    list.split(',')
        .map(|text| hex_f64(text).to_bits())
        .collect()
}

fn spectrum(unit: DriftTimeUnit, drift_time: f64) -> MSSpectrum {
    MSSpectrum {
        drift_time,
        drift_time_unit: unit,
        ..MSSpectrum::default()
    }
}

fn experiment(spectra: Vec<MSSpectrum>) -> MSExperiment {
    MSExperiment {
        spectra,
        ..MSExperiment::default()
    }
}

fn voltage_bits(experiment: &MSExperiment) -> Vec<u64> {
    FaimsHelper::get_compensation_voltages(experiment)
        .unwrap()
        .values()
        .map(f64::to_bits)
        .collect()
}

#[test]
fn constructor_and_destructor_sections() {
    // FAIMSHelper_test.cpp:32-39 allocates and deletes an instance of the
    // stateless class; the Rust type is a zero-sized unit struct.
    let constructed = FaimsHelper;
    assert_eq!(constructed, FaimsHelper);
    assert_eq!(std::mem::size_of::<FaimsHelper>(), 0);
    assert!(!std::mem::needs_drop::<FaimsHelper>());
}

#[cfg(feature = "mzml")]
fn load_im_faims_test() -> MSExperiment {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/faims_helper/IM_FAIMS_test.mzML");
    let file = std::fs::File::open(path).unwrap();
    openms::format::mzml::read(std::io::BufReader::new(file)).unwrap()
}

#[cfg(feature = "mzml")]
fn assert_class_literals_and_oracle_voltages(experiment: &MSExperiment) {
    let cvs = FaimsHelper::get_compensation_voltages(experiment).unwrap();
    // FAIMSHelper_test.cpp:54-57.
    assert_eq!(cvs.voltages.len(), 3);
    for volts in [-65.0, -55.0, -45.0] {
        assert!(
            cvs.voltages
                .contains(&CompensationVoltage::new(volts).unwrap()),
            "{volts}"
        );
    }
    assert!(cvs.warnings.is_empty());
    let recorded = records("faims_file_cvs");
    assert_eq!(recorded.len(), 1);
    assert_eq!(recorded[0][1], "3");
    assert_eq!(
        cvs.values().map(f64::to_bits).collect::<Vec<_>>(),
        oracle_bits(recorded[0][2])
    );
}

#[cfg(feature = "mzml")]
#[test]
fn get_compensation_voltages_section_through_the_mzml_reader() {
    // FAIMSHelper_test.cpp:43-59 exactly: load IM_FAIMS_test.mzML, whose CVs
    // are spectrum-level cvParams, and query the voltages.
    let experiment = load_im_faims_test();
    assert_eq!(experiment.spectra.len(), 19);
    let rows = records("faims_file_spectrum");
    assert_eq!(rows.len(), 19);
    for (spectrum, row) in experiment.spectra.iter().zip(&rows) {
        assert_eq!(spectrum.native_id, row[2]);
        assert_eq!(spectrum.drift_time.to_bits(), hex_f64(row[4]).to_bits());
        assert_eq!(
            spectrum.drift_time_unit,
            to_drift_time_unit(row[5]).unwrap()
        );
    }
    assert_class_literals_and_oracle_voltages(&experiment);
}

#[cfg(feature = "mzml")]
#[test]
fn get_compensation_voltages_section_on_the_values_the_cpp_reader_produced() {
    // FAIMSHelper_test.cpp:43-59 with the reader step split off: the Rust
    // reader loads the 19 spectra, and each spectrum receives the drift time
    // and unit the C++ MzMLFile produced for it (the oracle's per-spectrum
    // records). This checks the helper on the real file's voltage sequence
    // independently of A3-FORMAT-IO's reader change; the test above, re-enabled
    // since A3-FORMAT-IO merged, checks the whole path through the Rust reader.
    let mut experiment = load_im_faims_test();
    assert_eq!(experiment.spectra.len(), 19); // FAIMSHelper_test.cpp:50
    assert_eq!(records("faims_file_spectra")[0][1], "19");
    let rows = records("faims_file_spectrum");
    assert_eq!(rows.len(), 19);
    for (index, (spectrum, row)) in experiment.spectra.iter_mut().zip(&rows).enumerate() {
        assert_eq!(row[1], index.to_string());
        assert_eq!(spectrum.native_id, row[2]);
        assert_eq!(spectrum.ms_level.to_string(), row[3]);
        spectrum.drift_time = hex_f64(row[4]);
        spectrum.drift_time_unit = to_drift_time_unit(row[5]).unwrap();
    }
    assert!(
        rows.iter().all(|row| row[5] == "FAIMS_CV"),
        "every spectrum of the file is a FAIMS spectrum"
    );
    assert_class_literals_and_oracle_voltages(&experiment);
}

#[test]
fn get_compensation_voltages_detects_faims_beyond_the_first_spectrum_and_ignores_the_sentinel_section()
 {
    // FAIMSHelper_test.cpp:61-87.
    let experiment = experiment(vec![
        spectrum(DriftTimeUnit::Millisecond, 12.3),
        spectrum(DriftTimeUnit::FaimsCompensationVoltage, -50.0),
        spectrum(
            DriftTimeUnit::FaimsCompensationVoltage,
            ImTypes::DRIFTTIME_NOT_SET,
        ),
    ]);
    let cvs = FaimsHelper::get_compensation_voltages(&experiment).unwrap();
    assert_eq!(cvs.voltages.len(), 1);
    assert!(
        cvs.voltages
            .contains(&CompensationVoltage::new(-50.0).unwrap())
    );
    // The source logs this warning (FAIMSHelper.cpp:50-53); the oracle's stderr
    // shows it right after this case's marker.
    assert_eq!(cvs.warnings, [FaimsHelper::MISSING_VOLTAGE_WARNING]);
}

#[test]
fn get_compensation_voltages_returns_empty_for_non_faims_section() {
    // FAIMSHelper_test.cpp:89-104.
    let experiment = experiment(vec![
        spectrum(DriftTimeUnit::Millisecond, 1.0),
        spectrum(DriftTimeUnit::Millisecond, 2.0),
    ]);
    let cvs = FaimsHelper::get_compensation_voltages(&experiment).unwrap();
    assert!(cvs.voltages.is_empty());
    assert!(cvs.warnings.is_empty());
}

#[test]
fn every_executed_voltage_case_matches_the_oracle_or_is_refused_for_nan() {
    let rows = records("faims_cvs");
    assert_eq!(rows.len(), 16);
    let mut warned = Vec::new();
    let mut refused = Vec::new();
    for row in &rows {
        let name = row[1];
        let inputs: Vec<(DriftTimeUnit, f64)> = if row[4] == "-" {
            Vec::new()
        } else {
            row[4]
                .split(',')
                .map(|item| {
                    let (unit, value) = item.rsplit_once(':').expect("unit:value");
                    (to_drift_time_unit(unit).unwrap(), hex_f64(value))
                })
                .collect()
        };
        let experiment = experiment(
            inputs
                .iter()
                .map(|&(unit, value)| spectrum(unit, value))
                .collect(),
        );
        let first_faims_nan = inputs.iter().position(|&(unit, value)| {
            unit == DriftTimeUnit::FaimsCompensationVoltage && value.is_nan()
        });
        let result = FaimsHelper::get_compensation_voltages(&experiment);
        if let Some(index) = first_faims_nan {
            match result {
                Err(Error::InvalidValue(message)) => {
                    assert!(message.contains(&format!("spectrum {index} ")), "{message}");
                }
                other => panic!("{name}: expected a NaN refusal, got {other:?}"),
            }
            refused.push((name, row[2], row[3]));
            continue;
        }
        let cvs = result.unwrap();
        assert_eq!(cvs.voltages.len().to_string(), row[2], "{name}");
        assert_eq!(
            cvs.values().map(f64::to_bits).collect::<Vec<_>>(),
            oracle_bits(row[3]),
            "{name}"
        );
        if !cvs.warnings.is_empty() {
            assert_eq!(cvs.warnings, [FaimsHelper::MISSING_VOLTAGE_WARNING]);
            warned.push(name);
        }
    }
    assert_eq!(
        warned,
        ["class_section_beyond_first_and_sentinel", "only_sentinel"]
    );
    // What the source does with a NaN voltage depends on where it is inserted:
    // as the first element it becomes the tree root, every later voltage
    // compares equal to it and is dropped, and erase(DRIFTTIME_NOT_SET) then
    // removes the NaN itself (and logs the missing-voltage warning, the third
    // occurrence the oracle stderr counts). Anywhere else the NaN is dropped.
    assert_eq!(
        refused,
        [
            ("nan_first", "0", "-"),
            ("nan_middle", "2", "-0x1.9p+5,-0x1.4p+5"),
            ("nan_last", "2", "-0x1.9p+5,-0x1.4p+5"),
        ]
    );
}

#[test]
fn voltages_order_numerically_and_not_by_bit_pattern() {
    let ordered: Vec<CompensationVoltage> = [
        f64::NEG_INFINITY,
        -65.0,
        -55.0,
        -45.0,
        -0.0,
        0.5,
        10.0,
        f64::INFINITY,
    ]
    .into_iter()
    .map(|volts| CompensationVoltage::new(volts).unwrap())
    .collect();
    assert!(ordered.windows(2).all(|pair| pair[0] < pair[1]));
    // Why the key is not f64::to_bits: sign-magnitude patterns order negative
    // voltages by magnitude, after all positive ones.
    assert!((-45.0f64).to_bits() < (-65.0f64).to_bits());
    assert!(10.0f64.to_bits() < (-45.0f64).to_bits());

    let shuffled = experiment(
        [-45.0, 10.0, -65.0, 0.5, -55.0]
            .into_iter()
            .map(|volts| spectrum(DriftTimeUnit::FaimsCompensationVoltage, volts))
            .collect(),
    );
    assert_eq!(
        FaimsHelper::get_compensation_voltages(&shuffled)
            .unwrap()
            .values()
            .collect::<Vec<_>>(),
        [-65.0, -55.0, -45.0, 0.5, 10.0]
    );
}

#[test]
fn signed_zeros_are_one_voltage_that_keeps_the_first_sign() {
    let negative = CompensationVoltage::new(-0.0).unwrap();
    let positive = CompensationVoltage::new(0.0).unwrap();
    assert_eq!(negative, positive);
    let hash = |voltage: CompensationVoltage| {
        let mut hasher = DefaultHasher::new();
        voltage.hash(&mut hasher);
        hasher.finish()
    };
    assert_eq!(hash(negative), hash(positive));
    assert!(negative.volts().is_sign_negative());
    assert!(positive.volts().is_sign_positive());

    let unit = DriftTimeUnit::FaimsCompensationVoltage;
    assert_eq!(
        voltage_bits(&experiment(vec![spectrum(unit, -0.0), spectrum(unit, 0.0)])),
        [(-0.0f64).to_bits()]
    );
    assert_eq!(
        voltage_bits(&experiment(vec![spectrum(unit, 0.0), spectrum(unit, -0.0)])),
        [0.0f64.to_bits()]
    );
    // No tolerance: adjacent doubles are distinct voltages.
    let nearly = f64::from_bits((-50.0f64).to_bits() - 1);
    assert_eq!(
        voltage_bits(&experiment(vec![
            spectrum(unit, -50.0),
            spectrum(unit, nearly)
        ])),
        [(-50.0f64).to_bits(), nearly.to_bits()]
    );
}

#[test]
fn compensation_voltage_refuses_nan_only() {
    assert!(matches!(
        CompensationVoltage::new(f64::NAN),
        Err(Error::InvalidValue(_))
    ));
    assert!(CompensationVoltage::try_from(f64::NAN).is_err());
    for volts in [f64::NEG_INFINITY, -65.0, 0.0, f64::INFINITY] {
        let voltage = CompensationVoltage::try_from(volts).unwrap();
        assert_eq!(f64::from(voltage).to_bits(), volts.to_bits());
    }
}

#[test]
fn nan_drift_times_matter_only_on_faims_spectra() {
    // nan_on_non_faims in the oracle: a NaN with another unit is ignored.
    let experiment = experiment(vec![
        spectrum(DriftTimeUnit::Millisecond, f64::NAN),
        spectrum(DriftTimeUnit::None, f64::NAN),
        spectrum(DriftTimeUnit::FaimsCompensationVoltage, -40.0),
    ]);
    assert_eq!(voltage_bits(&experiment), [(-40.0f64).to_bits()]);
    let sentinel_elsewhere = MSExperiment {
        spectra: vec![
            spectrum(DriftTimeUnit::Millisecond, ImTypes::DRIFTTIME_NOT_SET),
            spectrum(DriftTimeUnit::FaimsCompensationVoltage, -50.0),
        ],
        ..MSExperiment::default()
    };
    let cvs = FaimsHelper::get_compensation_voltages(&sentinel_elsewhere).unwrap();
    assert_eq!(cvs.values().collect::<Vec<_>>(), [-50.0]);
    assert!(cvs.warnings.is_empty());
}

#[test]
fn constants_carry_the_source_values_and_the_native_ceilings() {
    assert_eq!(FaimsHelper::DEFAULT_CV_TOLERANCE, 0.01);
    assert_eq!(
        FaimsHelper::MISSING_VOLTAGE_WARNING,
        "Warning: FAIMS compensation voltage is missing for at least one spectrum!"
    );
    assert_eq!(FaimsHelper::MAX_SPECTRA, ImTypes::MAX_SPECTRA);
    assert_eq!(FaimsHelper::MAX_SPECTRA, 100_000_000);
    assert_eq!(FaimsHelper::MAX_PEPTIDE_IDENTIFICATIONS, 100_000_000);
    assert_eq!(FAIMS_CV, "FAIMS_CV");
    let empty = FaimsHelper::get_compensation_voltages(&MSExperiment::new()).unwrap();
    assert!(empty.voltages.is_empty() && empty.warnings.is_empty());
}

fn identification(id: &str, key: &str, value: Option<MetaValue>) -> PeptideIdentification {
    let mut peptide = PeptideIdentification::new();
    peptide.identifier = id.to_owned();
    if let Some(value) = value {
        peptide.metadata.insert(key.to_owned(), value);
    }
    peptide
}

fn float(value: f64) -> Option<MetaValue> {
    Some(MetaValue::try_from(value).unwrap())
}

/// The nine identifications the oracle driver built, in its order.
fn oracle_identifications() -> Vec<PeptideIdentification> {
    vec![
        identification("p0_double_-45", FAIMS_CV, float(-45.0)),
        identification("p1_unannotated", FAIMS_CV, None),
        identification("p2_double_-44.995", FAIMS_CV, float(-44.995)),
        identification("p3_double_-44.99", FAIMS_CV, float(-44.99)),
        identification("p4_int_-45", FAIMS_CV, Some(MetaValue::from(-45_i64))),
        identification("p5_double_-55", FAIMS_CV, float(-55.0)),
        identification("p6_double_-44.5", FAIMS_CV, float(-44.5)),
        identification("p7_double_+45", FAIMS_CV, float(45.0)),
        identification("p8_other_key_only", "FAIMS", float(-45.0)),
    ]
}

fn identifiers(peptides: &[PeptideIdentification]) -> String {
    if peptides.is_empty() {
        return "-".to_owned();
    }
    peptides
        .iter()
        .map(|peptide| peptide.identifier.as_str())
        .collect::<Vec<_>>()
        .join(",")
}

#[test]
fn filter_matches_every_executed_oracle_case_or_refuses_the_silent_ones() {
    let rows = records("faims_filter");
    assert_eq!(rows.len(), 15);
    assert_eq!(records("faims_filter_empty_value_exists")[0][1], "1");
    let mut refused = Vec::new();
    let mut infinite_targets = Vec::new();
    for row in &rows {
        let name = row[1];
        let target = hex_f64(row[2]);
        let tolerance = hex_f64(row[3]);
        let input = match name {
            "empty_input" => Vec::new(),
            "empty_datavalue_annotation" => vec![
                identification("e0_unannotated", FAIMS_CV, None),
                identification(
                    "e1_empty_datavalue",
                    FAIMS_CV,
                    Some(MetaValue::new(MetaValueData::Empty).unwrap()),
                ),
            ],
            _ => oracle_identifications(),
        };
        let result = FaimsHelper::filter_peptides_by_faims_cv(&input, target, tolerance);
        if name == "empty_datavalue_annotation" {
            // The source throws Exception::ConversionError converting the empty
            // DataValue to double.
            assert_eq!(row[4], "exception:ConversionError");
            match result {
                Err(Error::InvalidValue(message)) => {
                    assert!(message.contains("peptide identification 1 "), "{message}");
                }
                other => panic!("expected a refusal, got {other:?}"),
            }
            continue;
        }
        assert_eq!(row[4], "ok", "{name}");
        // An infinite target is not silent in this sense: it is a voltage
        // get_compensation_voltages can return, and it is filtered exactly as
        // the source filters it.
        let silent = target.is_nan() || tolerance.is_nan() || tolerance <= 0.0;
        if silent {
            assert!(
                matches!(result, Err(Error::InvalidValue(_))),
                "{name}: {result:?}"
            );
            // The source result the refusal replaces: its strict comparison can
            // never succeed, so only the unannotated identifications survive.
            assert_eq!(row[5], "p1_unannotated,p8_other_key_only", "{name}");
            refused.push(name);
        } else {
            if target.is_infinite() {
                infinite_targets.push(name);
            }
            assert_eq!(identifiers(&result.unwrap()), row[5], "{name}");
        }
    }
    assert_eq!(
        refused,
        [
            "tolerance_zero",
            "tolerance_negative",
            "tolerance_nan",
            "target_nan",
            "target_nan_tolerance_infinite"
        ]
    );
    assert_eq!(
        infinite_targets,
        [
            "target_positive_infinity",
            "target_negative_infinity",
            "target_positive_infinity_tolerance_infinite",
            "target_negative_infinity_tolerance_infinite"
        ]
    );
}

#[test]
fn filter_keeps_the_strict_boundary_order_and_input_unchanged() {
    let input = oracle_identifications();
    let before = input.clone();
    // |-44.5 - -45| is exactly 0.5 and is dropped; -44.99 differs from -45 by
    // 0.00999999999999801 in binary arithmetic and is kept at the default.
    let half = FaimsHelper::filter_peptides_by_faims_cv(&input, -45.0, 0.5).unwrap();
    assert!(!half.iter().any(|p| p.identifier == "p6_double_-44.5"));
    let default =
        FaimsHelper::filter_peptides_by_faims_cv(&input, -45.0, FaimsHelper::DEFAULT_CV_TOLERANCE)
            .unwrap();
    assert_eq!(
        identifiers(&default),
        "p0_double_-45,p1_unannotated,p2_double_-44.995,p3_double_-44.99,p4_int_-45,p8_other_key_only"
    );
    assert_eq!(default[0], input[0]);
    assert_eq!(input, before);
    let everything =
        FaimsHelper::filter_peptides_by_faims_cv(&input, -45.0, f64::INFINITY).unwrap();
    assert_eq!(everything, input);
}

#[test]
fn filter_refuses_non_numeric_annotations_and_accepts_infinite_targets() {
    let unannotated = identification("u", FAIMS_CV, None);
    for value in [
        MetaValue::from("-45"),
        MetaValue::new(MetaValueData::FloatList(vec![-45.0])).unwrap(),
        MetaValue::new(MetaValueData::IntegerList(vec![-45])).unwrap(),
    ] {
        let input = vec![
            unannotated.clone(),
            identification("bad", FAIMS_CV, Some(value)),
        ];
        match FaimsHelper::filter_peptides_by_faims_cv(&input, -45.0, 0.01) {
            Err(Error::InvalidValue(message)) => {
                assert!(message.contains("peptide identification 1 "), "{message}");
            }
            other => panic!("expected a refusal, got {other:?}"),
        }
    }
    // FAIMSHelper.cpp has no target check. The executed source keeps only the
    // unannotated identifications for an infinite target (oracle
    // target_positive_infinity, target_negative_infinity), because
    // |cv - target| is infinite and never less than the tolerance.
    for (target, case) in [
        (f64::INFINITY, "target_positive_infinity"),
        (f64::NEG_INFINITY, "target_negative_infinity"),
    ] {
        let row = records("faims_filter")
            .into_iter()
            .find(|row| row[1] == case)
            .expect("oracle case");
        assert_eq!(hex_f64(row[2]).to_bits(), target.to_bits(), "{case}");
        assert_eq!(hex_f64(row[3]).to_bits(), 0.01f64.to_bits(), "{case}");
        assert_eq!(row[4], "ok", "{case}");
        assert_eq!(row[5], "p1_unannotated,p8_other_key_only", "{case}");
        assert_eq!(
            FaimsHelper::filter_peptides_by_faims_cv(
                std::slice::from_ref(&unannotated),
                target,
                0.01
            )
            .unwrap(),
            std::slice::from_ref(&unannotated)
        );
        let kept =
            FaimsHelper::filter_peptides_by_faims_cv(&oracle_identifications(), target, 0.01)
                .unwrap();
        assert_eq!(identifiers(&kept), row[5], "{case}");
    }
    let kept =
        FaimsHelper::filter_peptides_by_faims_cv(std::slice::from_ref(&unannotated), -45.0, 1e-300)
            .unwrap();
    assert_eq!(kept, [unannotated]);
}

#[test]
fn an_infinite_target_keeps_no_annotated_identification_not_even_an_infinite_one() {
    // The oracle filtered {unannotated, +inf, -inf, -45} by +inf and by -inf at
    // an infinite tolerance, and C++ kept only the unannotated identification:
    // the same infinity gives |inf - inf| = NaN, the opposite infinity and -45
    // give inf, and neither is less than inf. The two infinite annotations
    // cannot be built here, because MetaValue refuses non-finite floats, so
    // the port runs the two representable identifications and must keep what
    // C++ kept: nothing C++ would keep is lost.
    assert!(MetaValue::try_from(f64::INFINITY).is_err());
    assert!(MetaValue::try_from(f64::NEG_INFINITY).is_err());
    let representable = [
        identification("i0_unannotated", FAIMS_CV, None),
        identification("i3_double_-45", FAIMS_CV, float(-45.0)),
    ];
    let rows = records("faims_filter_infinite_annotation");
    assert_eq!(rows.len(), 2);
    for (row, (case, target)) in rows.iter().zip([
        ("target_positive_infinity", f64::INFINITY),
        ("target_negative_infinity", f64::NEG_INFINITY),
    ]) {
        assert_eq!(row[1], case);
        assert_eq!(hex_f64(row[2]).to_bits(), target.to_bits(), "{case}");
        let tolerance = hex_f64(row[3]);
        assert_eq!(tolerance, f64::INFINITY, "{case}");
        assert_eq!(row[4], "ok", "{case}");
        assert_eq!(row[5], "i0_unannotated", "{case}");
        let kept =
            FaimsHelper::filter_peptides_by_faims_cv(&representable, target, tolerance).unwrap();
        assert_eq!(identifiers(&kept), row[5], "{case}");
    }
}

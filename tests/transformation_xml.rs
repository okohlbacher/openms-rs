// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FORMAT/TransformationXMLFile.h`: TrafoXML reading and writing.
//!
//! Expectations come from `TransformationXMLFile_test.cpp` and its four
//! retained `.trafoXML` fixtures. See
//! `tests/data/transformation_xml_provenance.json`.

#![cfg(any(feature = "featurexml", feature = "consensusxml"))]

use openms::Error;
use openms::analysis::transformations::{
    DataPoint, Extrapolation, Interpolation, InterpolationOptions, LinearOptions, ModelConfig,
    TransformationDescription, WeightFunction,
};
use openms::format::transformation_xml::{
    self as trafo, ReadOptions, SCHEMA_LOCATION, TransformationRecord, VERSION, WriteOptions,
};
use openms::param::ParamValue;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// The class test writes pi and e out to nineteen digits as its slope and
/// intercept; the f64 constants are those literals.
use std::f64::consts::{E, PI};

/// `TEST_REAL_SIMILAR`'s default relative tolerance.
const REAL_SIMILAR: f64 = 1e-5;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn similar(got: f64, want: f64) {
    assert!(
        (got - want).abs() <= REAL_SIMILAR * want.abs().max(1.0),
        "{got} is not similar to {want}"
    );
}

// START_SECTION((TransformationXMLFile())) — the source constructor fixes the
// handler version and the schema it validates against. This port has no handler
// object; the two constants and the option defaults carry that state.
#[test]
fn the_constructor_state_is_version_1_1_and_its_schema() {
    assert_eq!(VERSION, "1.1");
    assert!(
        SCHEMA_LOCATION.ends_with("/SCHEMAS/TrafoXML_1_1.xsd"),
        "{SCHEMA_LOCATION}"
    );
    // The source `load` defaults to fitting the named model.
    assert!(ReadOptions::default().fit_model);
    assert_eq!(WriteOptions::default().max_records, 1_000_000);
}

// START_SECTION([EXTRA] static bool isValid(const std::string& filename)) —
// files 1, 2 and 4 validate against the schema and file 3 does not. This port
// has no XSD validator; its structural reader stands in, and reaches the same
// verdict on all four files. File 3 carries a second `</Transformation>` close
// tag, so it is not even well-formed XML.
#[test]
fn the_three_schema_valid_fixtures_read_and_the_invalid_one_does_not() {
    for name in [
        "transformation_xml_1.trafoXML",
        "transformation_xml_2.trafoXML",
        "transformation_xml_4.trafoXML",
    ] {
        trafo::load(data(name)).unwrap_or_else(|e| panic!("{name} should read: {e}"));
    }
    let error = trafo::load(data("transformation_xml_3.trafoXML")).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
}

// START_SECTION(void load(const std::string & filename,
//   TransformationDescription & transformation, bool fit_model=true))
#[test]
fn loading_fits_the_named_model_from_the_file() {
    // File 1 names the "none" model and carries no parameter.
    let record = trafo::load_record_with_options(
        data("transformation_xml_1.trafoXML"),
        &ReadOptions::default(),
    )
    .unwrap();
    assert_eq!(record.model_type, "none");
    assert!(record.parameters.is_empty());
    let description = trafo::load(data("transformation_xml_1.trafoXML")).unwrap();
    assert_eq!(description.model().name(), "none");

    // File 2 names "linear" with slope and intercept. The fixture stores them
    // at six significant digits, which is why the class test compares them
    // with TEST_REAL_SIMILAR rather than exactly.
    let record = trafo::load_record_with_options(
        data("transformation_xml_2.trafoXML"),
        &ReadOptions::default(),
    )
    .unwrap();
    assert_eq!(record.model_type, "linear");
    assert_eq!(record.parameters.len(), 2);
    similar(record.parameters["slope"].to_f64().unwrap(), PI);
    similar(record.parameters["intercept"].to_f64().unwrap(), E);
    let description = trafo::load(data("transformation_xml_2.trafoXML")).unwrap();
    assert_eq!(description.model().name(), "linear");
    // Derived from the file's own coefficients: y = slope * x + intercept.
    // The fixture rounds pi and e to six significant digits, which is inside
    // the tolerance the class test itself uses.
    similar(description.apply(1.0).unwrap(), PI + E);

    // File 4 names "interpolated" with a linear interpolation type and three
    // pairs.
    let record = trafo::load_record_with_options(
        data("transformation_xml_4.trafoXML"),
        &ReadOptions::default(),
    )
    .unwrap();
    assert_eq!(record.model_type, "interpolated");
    assert_eq!(
        record.parameters["interpolation_type"],
        ParamValue::String("linear".into())
    );
    let description = trafo::load(data("transformation_xml_4.trafoXML")).unwrap();
    assert_eq!(description.model().name(), "interpolated");
    let points = description.data_points();
    assert_eq!(points.len(), 3);
    similar(points[0].x, 1.2);
    similar(points[1].x, 2.2);
    similar(points[2].x, 3.2);
    similar(points[0].y, 5.2);
    similar(points[1].y, 6.25);
    similar(points[2].y, 7.3);
    // The fitted model reproduces its own anchors.
    similar(description.apply(2.2).unwrap(), 6.25);

    // Not performing the model fit leaves the pairs and no model, because the
    // source's setDataPoints resets the model even when it was "identity".
    let unfitted = trafo::load_with_options(
        data("transformation_xml_4.trafoXML"),
        &ReadOptions {
            fit_model: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(unfitted.model().name(), "none");
    assert_eq!(unfitted.data_points().len(), 3);
    let unfitted = trafo::load_with_options(
        data("transformation_xml_2.trafoXML"),
        &ReadOptions {
            fit_model: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(unfitted.model().name(), "none");
    assert!(matches!(unfitted.model_config(), ModelConfig::None));
}

// START_SECTION(void store(std::string filename,
//   const TransformationDescription& transformation))
#[test]
fn storing_and_reloading_preserves_the_model_and_its_parameters() {
    let directory = std::env::temp_dir().join(format!("openms-trafo-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();

    // "none": no parameter is written, and the reload reports no parameter.
    let mut description = TransformationDescription::default();
    description.fit_model(ModelConfig::None).unwrap();
    let none_path = directory.join("none.trafoXML");
    trafo::store(&none_path, &description).unwrap();
    let record = trafo::load_record_with_options(&none_path, &ReadOptions::default()).unwrap();
    assert_eq!(record.model_type, "none");
    assert!(record.parameters.is_empty());
    assert_eq!(trafo::load(&none_path).unwrap().model().name(), "none");
    // The source exercises apply() on a stored "none" model; it is the identity.
    similar(description.apply(234255132.43212).unwrap(), 234255132.43212);

    // "linear" from explicit coefficients: slope and intercept round-trip at
    // full precision, unlike the pre-3.0 fixture's six digits.
    let mut description = TransformationDescription::default();
    description
        .fit_model(ModelConfig::Linear(LinearOptions {
            coefficients: Some(openms::analysis::transformations::LinearCoefficients {
                slope: PI,
                intercept: E,
            }),
            ..Default::default()
        }))
        .unwrap();
    let linear_path = directory.join("linear.trafoXML");
    trafo::store(&linear_path, &description).unwrap();
    let record = trafo::load_record_with_options(&linear_path, &ReadOptions::default()).unwrap();
    assert_eq!(record.model_type, "linear");
    assert_eq!(record.parameters.len(), 2);
    similar(record.parameters["slope"].to_f64().unwrap(), PI);
    similar(record.parameters["intercept"].to_f64().unwrap(), E);
    let reloaded = trafo::load(&linear_path).unwrap();
    assert_eq!(reloaded.model().name(), "linear");
    similar(
        reloaded.apply(234255132.43212).unwrap(),
        description.apply(234255132.43212).unwrap(),
    );

    // "interpolated" with three pairs: the reload carries both the
    // interpolation type that was set and the extrapolation default the source
    // fills in, so two parameters and three pairs.
    let pairs = vec![
        DataPoint::new(1.2, 5.2),
        DataPoint::new(2.2, 6.25),
        DataPoint::new(3.2, 7.3),
    ];
    let mut description = TransformationDescription::new(pairs).unwrap();
    description
        .fit_model(ModelConfig::Interpolated(InterpolationOptions {
            interpolation: Interpolation::Linear,
            ..Default::default()
        }))
        .unwrap();
    let pairs_path = directory.join("pairs.trafoXML");
    trafo::store(&pairs_path, &description).unwrap();
    let record = trafo::load_record_with_options(&pairs_path, &ReadOptions::default()).unwrap();
    assert_eq!(record.model_type, "interpolated");
    assert_eq!(record.parameters.len(), 2);
    assert_eq!(
        record.parameters["interpolation_type"],
        ParamValue::String("linear".into())
    );
    assert_eq!(
        record.parameters["extrapolation_type"],
        ParamValue::String("two-point-linear".into())
    );
    let reloaded = trafo::load(&pairs_path).unwrap();
    assert_eq!(reloaded.model().name(), "interpolated");
    let points = reloaded.data_points();
    assert_eq!(points.len(), 3);
    similar(points[0].x, 1.2);
    similar(points[1].x, 2.2);
    similar(points[2].x, 3.2);
    similar(points[0].y, 5.2);
    similar(points[1].y, 6.25);
    similar(points[2].y, 7.3);
    similar(
        reloaded.apply(234255132.43212).unwrap(),
        description.apply(234255132.43212).unwrap(),
    );

    // TEST_EXCEPTION(Exception::IllegalArgument,
    //   trafo.fitModel("mumble_pfrwoarpfz")): an unknown model name is refused
    // when the record is resolved rather than when it is read.
    let record = TransformationRecord {
        model_type: "mumble_pfrwoarpfz".into(),
        parameters: BTreeMap::new(),
        data: Vec::new(),
    };
    let error = record.model_config().unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
    // b_spline is a model the source implements and this port does not.
    let record = TransformationRecord {
        model_type: "b_spline".into(),
        ..Default::default()
    };
    assert!(matches!(
        record.model_config().unwrap_err(),
        Error::Unsupported(_)
    ));

    std::fs::remove_dir_all(&directory).unwrap();
}

/// The source refuses to write a transformation with an empty model name.
#[test]
fn an_empty_model_name_is_refused() {
    let mut bytes = Vec::new();
    let error = trafo::write_record_with_options(
        &mut bytes,
        &TransformationRecord::default(),
        &WriteOptions::default(),
    )
    .unwrap_err();
    assert!(
        format!("{error}").contains("empty name"),
        "unexpected error: {error}"
    );
    assert!(bytes.is_empty());
}

/// Notes on pairs survive a round trip, including non-ASCII and XML-special
/// text, which the source escapes with `writeXMLEscape` on writing only.
#[test]
fn pair_notes_round_trip_through_escaping() {
    let mut description = TransformationDescription::new(vec![
        DataPoint::with_note(1.2, 5.2, "日本語 & <anchor>"),
        DataPoint::new(2.2, 6.25),
        DataPoint::with_note(3.2, 7.3, "\"quoted\""),
    ])
    .unwrap();
    description.fit_model(ModelConfig::None).unwrap();
    let mut bytes = Vec::new();
    trafo::write(&mut bytes, &description).unwrap();
    let record =
        trafo::read_record_with_options(bytes.as_slice(), &ReadOptions::default()).unwrap();
    assert_eq!(record.data[0].note, "日本語 & <anchor>");
    assert_eq!(record.data[1].note, "");
    assert_eq!(record.data[2].note, "\"quoted\"");
}

/// Coordinate weights are written and read back, so a weighted linear model
/// survives a round trip instead of silently losing its weighting.
#[test]
fn coordinate_weights_round_trip() {
    let mut description = TransformationDescription::new(vec![
        DataPoint::new(1.0, 2.0),
        DataPoint::new(2.0, 4.0),
        DataPoint::new(3.0, 6.0),
    ])
    .unwrap();
    description
        .fit_model(ModelConfig::Linear(LinearOptions {
            x_weight: openms::analysis::transformations::CoordinateWeight {
                function: WeightFunction::Log,
                min: 1e-15,
                max: 1e15,
            },
            ..Default::default()
        }))
        .unwrap();
    let mut bytes = Vec::new();
    trafo::write(&mut bytes, &description).unwrap();
    let text = String::from_utf8(bytes.clone()).unwrap();
    assert!(text.contains("name=\"x_weight\""), "{text}");
    assert!(text.contains("value=\"ln(x)\""), "{text}");
    let record =
        trafo::read_record_with_options(bytes.as_slice(), &ReadOptions::default()).unwrap();
    let ModelConfig::Linear(options) = record.model_config().unwrap() else {
        panic!("expected a linear model");
    };
    assert_eq!(options.x_weight.function, WeightFunction::Log);
    assert_eq!(options.y_weight.function, WeightFunction::Identity);
    // An invalid weight name is refused, as the source's model constructor is.
    let mut record = record;
    record
        .parameters
        .insert("x_weight".into(), ParamValue::String("1/z".into()));
    assert!(matches!(
        record.model_config().unwrap_err(),
        Error::InvalidValue(_)
    ));
}

/// A document version above the parser's is refused rather than warned about.
#[test]
fn a_newer_document_version_is_refused() {
    let text = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <TrafoXML version=\"9.9\">\n<Transformation name=\"none\"/>\n</TrafoXML>\n";
    let error = trafo::read(text.as_bytes()).unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error:?}");
}

/// An unknown element inside `<Transformation>` is refused; the source logs it
/// at debug level and drops it.
#[test]
fn an_unknown_element_is_refused() {
    let text = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <TrafoXML version=\"1.1\">\n<Transformation name=\"none\">\n\
         <Mystery value=\"1\"/>\n</Transformation>\n</TrafoXML>\n";
    let error = trafo::read(text.as_bytes()).unwrap_err();
    assert!(matches!(error, Error::Unsupported(_)), "{error:?}");
}

/// Extrapolation names round-trip through the file.
#[test]
fn extrapolation_names_round_trip() {
    for (extrapolation, name) in [
        (Extrapolation::TwoPointLinear, "two-point-linear"),
        (Extrapolation::FourPointLinear, "four-point-linear"),
        (Extrapolation::GlobalLinear, "global-linear"),
    ] {
        let mut description = TransformationDescription::new(vec![
            DataPoint::new(1.0, 2.0),
            DataPoint::new(2.0, 4.0),
            DataPoint::new(3.0, 6.0),
        ])
        .unwrap();
        description
            .fit_model(ModelConfig::Interpolated(InterpolationOptions {
                interpolation: Interpolation::Linear,
                extrapolation,
                ..Default::default()
            }))
            .unwrap();
        let mut bytes = Vec::new();
        trafo::write(&mut bytes, &description).unwrap();
        let record =
            trafo::read_record_with_options(bytes.as_slice(), &ReadOptions::default()).unwrap();
        assert_eq!(
            record.parameters["extrapolation_type"],
            ParamValue::String(name.into())
        );
        let ModelConfig::Interpolated(options) = record.model_config().unwrap() else {
            panic!("expected an interpolated model");
        };
        assert_eq!(options.extrapolation, extrapolation);
    }
}

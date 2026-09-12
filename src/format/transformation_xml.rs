// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded native TrafoXML 1.1 interchange for retention-time transformations.
//!
//! Ports `FORMAT/TransformationXMLFile.h`. A TrafoXML document holds one
//! `<Transformation>`: the fitted model's name, the parameters it was fitted
//! with, and the `(from, to)` pairs the fit came from. Reading rebuilds a
//! [`TransformationDescription`](crate::analysis::transformations::TransformationDescription),
//! optionally refitting the named model; writing serialises one.
//!
//! See `docs/TRANSFORMATION_XML_SUPPORT.md` for the API mapping, the preserved
//! source conventions and the native differences.

use super::identification_xml::{self as xml, Node};
use crate::analysis::transformations::{
    CoordinateWeight, DataPoint, Extrapolation, Interpolation, InterpolationOptions,
    LinearCoefficients, LinearOptions, LowessOptions, ModelConfig, TransformationDescription,
    TransformationModel, WeightFunction,
};
use crate::param::ParamValue;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::{BufRead, Write};
use std::path::Path;

/// Document version this port writes, and the highest it accepts.
pub const VERSION: &str = "1.1";
/// Schema location written into the root element, as in the source.
pub const SCHEMA_LOCATION: &str =
    "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/TrafoXML_1_1.xsd";

fn bad(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}
fn unsupported(message: impl Into<String>) -> Error {
    Error::Unsupported(message.into())
}
fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
fn finite(value: f64, label: &str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad(format!("TrafoXML {label} must be finite")))
    }
}
/// `StringUtils::toStr(double)`: the source writes attribute values at full
/// precision, so the port reuses the crate's reproduction of that formatter.
fn text(value: f64) -> String {
    crate::param::value::format_float(value, true)
}

/// Reading limits and the `fit_model` switch of `TransformationXMLFile::load`.
#[derive(Clone, Copy, Debug)]
pub struct ReadOptions {
    /// Maximum undecoded input bytes.
    pub max_xml_bytes: u64,
    /// Maximum XML elements in the document, `<Pair>` entries included.
    pub max_records: usize,
    /// Maximum entries in one decoded list.
    pub max_list_items: usize,
    /// Maximum cumulative decoded-tree allocation.
    pub max_payload_bytes: usize,
    /// Maximum cumulative parser work units.
    pub max_work: usize,
    /// Maximum `<Pair>` elements kept, checked before any pair is stored.
    pub max_data_points: usize,
    /// Refit the named model after reading, as the source `fit_model` argument.
    /// With `false` the returned description keeps the pairs and no model, which
    /// is what the source leaves behind because `setDataPoints` resets the model.
    pub fit_model: bool,
}
impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            max_xml_bytes: 64 * 1024 * 1024,
            max_records: 1_000_000,
            max_list_items: 1_000_000,
            max_payload_bytes: 256 * 1024 * 1024,
            max_work: 50_000_000,
            max_data_points: 1_000_000,
            fit_model: true,
        }
    }
}
impl ReadOptions {
    fn xml(&self) -> xml::ReadOptions {
        xml::ReadOptions {
            max_xml_bytes: self.max_xml_bytes,
            max_records: self.max_records,
            max_list_items: self.max_list_items,
        }
    }
    fn validate(&self) -> Result<()> {
        if self.max_xml_bytes == 0
            || self.max_records == 0
            || self.max_list_items == 0
            || self.max_payload_bytes == 0
            || self.max_work == 0
            || self.max_data_points == 0
        {
            return Err(invalid("TrafoXML limits must be positive"));
        }
        Ok(())
    }
}

/// Writing limits.
#[derive(Clone, Copy, Debug)]
pub struct WriteOptions {
    /// Maximum serialised output bytes.
    pub max_xml_bytes: usize,
    /// Maximum XML elements written.
    pub max_records: usize,
}
impl Default for WriteOptions {
    fn default() -> Self {
        Self {
            max_xml_bytes: 64 * 1024 * 1024,
            max_records: 1_000_000,
        }
    }
}

/// One TrafoXML document as written: the model name, the parameters recorded
/// for it, and the fitted pairs.
///
/// The source keeps these three in `TransformationXMLFile`'s protected
/// `model_type_`, `params_` and `data_` members and hands them to
/// `TransformationDescription::fitModel`. They are public here because the
/// untyped `(name, Param)` pair is exactly what the file contains, and because
/// the source's `fit_model=false` path returns a description that has the pairs
/// but has thrown the model name away.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TransformationRecord {
    /// Value of the `<Transformation name="...">` attribute. `none` and
    /// `identity` name the two models that carry no parameters.
    pub model_type: String,
    /// `<Param>` entries in document order, keyed by name. Only the source's
    /// three written types occur when reading: `int`, `float` and `string`.
    pub parameters: BTreeMap<String, ParamValue>,
    /// `<Pair>` entries, with the optional `note` attribute preserved.
    pub data: Vec<DataPoint>,
}

impl TransformationRecord {
    /// The typed model configuration this record names.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] for `b_spline`, which the source
    /// implements and this port does not, and [`Error::InvalidValue`] for a
    /// name that is not a source model type at all. The source
    /// `TransformationDescription::fitModel` throws `Exception::IllegalArgument`
    /// for both, so an unknown name is an error here as it is there.
    ///
    /// Returns [`Error::Unsupported`] when a parameter has a type this model
    /// cannot accept, and [`Error::InvalidValue`] when a weight name is not one
    /// of the source's `x`/`1/x`/`1/x2`/`ln(x)` set. The source
    /// `TransformationModel` constructor throws `Exception::InvalidParameter`
    /// for an invalid weight.
    pub fn model_config(&self) -> Result<ModelConfig> {
        match self.model_type.as_str() {
            "none" => Ok(ModelConfig::None),
            "identity" => Ok(ModelConfig::Identity),
            "linear" => Ok(ModelConfig::Linear(self.linear_options()?)),
            "interpolated" => Ok(ModelConfig::Interpolated(self.interpolation_options()?)),
            "lowess" => Ok(ModelConfig::Lowess(self.lowess_options()?)),
            "b_spline" => Err(unsupported(
                "b_spline transformation model is not implemented in this port",
            )),
            other => Err(invalid(format!(
                "unknown transformation model type {other:?}"
            ))),
        }
    }

    /// Rebuild a description from this record.
    ///
    /// With `fit_model` the named model is fitted to the pairs, as
    /// `TransformationXMLFile::load` does by default. Without it the
    /// description carries the pairs and no model, which is the state the
    /// source leaves after `setDataPoints`; `getModelType` then reports `none`.
    ///
    /// # Errors
    ///
    /// Propagates [`Self::model_config`] and any error from fitting, for
    /// instance [`Error::InvalidValue`] when `linear` has neither pairs nor an
    /// explicit slope and intercept.
    pub fn description(&self, fit_model: bool) -> Result<TransformationDescription> {
        let mut description = TransformationDescription::new(self.data.clone())?;
        if fit_model {
            description.fit_model(self.model_config()?)?;
        }
        Ok(description)
    }

    /// The record a description serialises to.
    ///
    /// Named after the fitted model, with the parameters that model actually
    /// holds. See the module's support document for why this is a smaller set
    /// than the source writes on its data-fitted `linear` path.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] if a stored coordinate is not finite.
    pub fn from_description(description: &TransformationDescription) -> Result<Self> {
        let mut parameters = BTreeMap::new();
        let model_type = match description.model() {
            TransformationModel::None => "none",
            TransformationModel::Identity => "identity",
            TransformationModel::Linear(model) => {
                let coefficients = model.coefficients();
                parameters.insert(
                    "slope".to_owned(),
                    ParamValue::Float(finite(coefficients.slope, "slope")?),
                );
                parameters.insert(
                    "intercept".to_owned(),
                    ParamValue::Float(finite(coefficients.intercept, "intercept")?),
                );
                let options = model.options();
                insert_weights(&mut parameters, options.x_weight, options.y_weight)?;
                "linear"
            }
            TransformationModel::Interpolated(model) => {
                let options = model.options();
                parameters.insert(
                    "interpolation_type".to_owned(),
                    ParamValue::String(interpolation_name(options.interpolation).to_owned()),
                );
                parameters.insert(
                    "extrapolation_type".to_owned(),
                    ParamValue::String(extrapolation_name(options.extrapolation).to_owned()),
                );
                "interpolated"
            }
            TransformationModel::Lowess(model) => {
                let options = model.options();
                parameters.insert(
                    "span".to_owned(),
                    ParamValue::Float(finite(options.span, "span")?),
                );
                parameters.insert(
                    "num_iterations".to_owned(),
                    ParamValue::Integer(
                        i64::try_from(options.iterations)
                            .map_err(|_| invalid("lowess iteration count exceeds i64"))?,
                    ),
                );
                parameters.insert(
                    "delta".to_owned(),
                    ParamValue::Float(finite(options.delta.unwrap_or(-1.0), "delta")?),
                );
                parameters.insert(
                    "interpolation_type".to_owned(),
                    ParamValue::String(interpolation_name(options.interpolation).to_owned()),
                );
                parameters.insert(
                    "extrapolation_type".to_owned(),
                    ParamValue::String(extrapolation_name(options.extrapolation).to_owned()),
                );
                "lowess"
            }
        };
        for point in description.data_points() {
            finite(point.x, "pair source")?;
            finite(point.y, "pair target")?;
        }
        Ok(Self {
            model_type: model_type.to_owned(),
            parameters,
            data: description.data_points().to_vec(),
        })
    }

    fn float(&self, key: &str) -> Result<Option<f64>> {
        match self.parameters.get(key) {
            None | Some(ParamValue::Empty) => Ok(None),
            Some(value) => Ok(Some(finite(
                value
                    .to_f64()
                    .map_err(|_| unsupported(format!("TrafoXML parameter {key} is not numeric")))?,
                key,
            )?)),
        }
    }
    fn integer(&self, key: &str) -> Result<Option<i64>> {
        match self.parameters.get(key) {
            None | Some(ParamValue::Empty) => Ok(None),
            Some(value) => Ok(Some(value.to_i64().map_err(|_| {
                unsupported(format!("TrafoXML parameter {key} is not an integer"))
            })?)),
        }
    }
    fn string(&self, key: &str) -> Result<Option<&str>> {
        match self.parameters.get(key) {
            None | Some(ParamValue::Empty) => Ok(None),
            Some(ParamValue::String(value)) => Ok(Some(value)),
            Some(_) => Err(unsupported(format!(
                "TrafoXML parameter {key} is not a string"
            ))),
        }
    }
    fn weight(&self, key: &str, identity: &str) -> Result<CoordinateWeight> {
        let name = self.string(key)?.unwrap_or(identity);
        // A TrafoXML written before OpenMS 3.0 stores an empty weight name to
        // mean "unweighted"; the source base model maps "" to identity too.
        let function = match name {
            "" => WeightFunction::Identity,
            _ if name == identity => WeightFunction::Identity,
            "ln(x)" | "ln(y)" => WeightFunction::Log,
            "1/x" | "1/y" => WeightFunction::Reciprocal,
            "1/x2" | "1/y2" => WeightFunction::ReciprocalSquared,
            other => {
                return Err(invalid(format!(
                    "{other:?} is not a valid {key} for a linear transformation"
                )));
            }
        };
        let prefix = if identity == "x" { "x" } else { "y" };
        Ok(CoordinateWeight {
            function,
            min: self.float(&format!("{prefix}_datum_min"))?.unwrap_or(1e-15),
            max: self.float(&format!("{prefix}_datum_max"))?.unwrap_or(1e15),
        })
    }
    fn linear_options(&self) -> Result<LinearOptions> {
        // The source uses explicit coefficients only when there are no pairs.
        let coefficients = match (self.float("slope")?, self.float("intercept")?) {
            (Some(slope), Some(intercept)) => Some(LinearCoefficients { slope, intercept }),
            _ => None,
        };
        Ok(LinearOptions {
            x_weight: self.weight("x_weight", "x")?,
            y_weight: self.weight("y_weight", "y")?,
            coefficients,
            ..Default::default()
        })
    }
    fn interpolation(&self) -> Result<Interpolation> {
        match self.string("interpolation_type")? {
            None => Ok(Interpolation::default()),
            Some("linear") => Ok(Interpolation::Linear),
            Some("cspline") => Ok(Interpolation::CubicSpline),
            Some("polynomial") | Some("akima") => Err(unsupported(
                "polynomial and akima interpolation are not implemented in this port",
            )),
            Some(other) => Err(invalid(format!("unknown interpolation_type {other:?}"))),
        }
    }
    fn extrapolation(&self) -> Result<Extrapolation> {
        match self.string("extrapolation_type")? {
            None => Ok(Extrapolation::default()),
            Some("two-point-linear") => Ok(Extrapolation::TwoPointLinear),
            Some("four-point-linear") => Ok(Extrapolation::FourPointLinear),
            Some("global-linear") => Ok(Extrapolation::GlobalLinear),
            Some(other) => Err(invalid(format!("unknown extrapolation_type {other:?}"))),
        }
    }
    fn interpolation_options(&self) -> Result<InterpolationOptions> {
        Ok(InterpolationOptions {
            interpolation: self.interpolation()?,
            extrapolation: self.extrapolation()?,
            ..Default::default()
        })
    }
    fn lowess_options(&self) -> Result<LowessOptions> {
        let defaults = LowessOptions::default();
        let delta = self.float("delta")?.unwrap_or(-1.0);
        Ok(LowessOptions {
            span: self.float("span")?.unwrap_or(defaults.span),
            iterations: match self.integer("num_iterations")? {
                Some(value) => usize::try_from(value)
                    .map_err(|_| invalid("lowess num_iterations must be non-negative"))?,
                None => defaults.iterations,
            },
            // The source treats a negative delta as "choose automatically".
            delta: if delta < 0.0 { None } else { Some(delta) },
            interpolation: self.interpolation()?,
            extrapolation: self.extrapolation()?,
            ..defaults
        })
    }
}

fn insert_weights(
    parameters: &mut BTreeMap<String, ParamValue>,
    x_weight: CoordinateWeight,
    y_weight: CoordinateWeight,
) -> Result<()> {
    // The source writes the weight parameters only when its Param happens to
    // carry them. They are written here whenever they are not the identity, so
    // a weighted model round-trips instead of silently losing its weighting.
    for (weight, names) in [
        (x_weight, ["x_weight", "x_datum_min", "x_datum_max"]),
        (y_weight, ["y_weight", "y_datum_min", "y_datum_max"]),
    ] {
        if weight.function == WeightFunction::Identity {
            continue;
        }
        let axis = if names[0] == "x_weight" { 'x' } else { 'y' };
        let name = match weight.function {
            WeightFunction::Identity => unreachable!("identity weights are skipped above"),
            WeightFunction::Log => format!("ln({axis})"),
            WeightFunction::Reciprocal => format!("1/{axis}"),
            WeightFunction::ReciprocalSquared => format!("1/{axis}2"),
        };
        parameters.insert(names[0].to_owned(), ParamValue::String(name));
        parameters.insert(
            names[1].to_owned(),
            ParamValue::Float(finite(weight.min, names[1])?),
        );
        parameters.insert(
            names[2].to_owned(),
            ParamValue::Float(finite(weight.max, names[2])?),
        );
    }
    Ok(())
}
fn interpolation_name(value: Interpolation) -> &'static str {
    match value {
        Interpolation::Linear => "linear",
        Interpolation::CubicSpline => "cspline",
    }
}
fn extrapolation_name(value: Extrapolation) -> &'static str {
    match value {
        Extrapolation::TwoPointLinear => "two-point-linear",
        Extrapolation::FourPointLinear => "four-point-linear",
        Extrapolation::GlobalLinear => "global-linear",
    }
}

/// Read a TrafoXML document and refit its model, as `load` does by default.
///
/// # Errors
///
/// As [`read_with_options`].
pub fn read(reader: impl BufRead) -> Result<TransformationDescription> {
    read_with_options(reader, &ReadOptions::default())
}

/// Read a TrafoXML document under explicit limits.
///
/// # Errors
///
/// Returns [`Error::Parse`] for XML that is not well formed, for a root that is
/// not `TrafoXML`, for a missing required attribute and for a non-finite
/// coordinate; [`Error::Unsupported`] for an unknown element or attribute, for a
/// `<Param type>` this port does not represent, and for a document version
/// above [`VERSION`]; [`Error::InvalidValue`] for zero limits and for an unknown
/// model or parameter value. The source instead logs and continues for an
/// unknown element, an unsupported `<Param type>` and a too-new version, which
/// drops data silently; see the support document.
pub fn read_with_options(
    reader: impl BufRead,
    options: &ReadOptions,
) -> Result<TransformationDescription> {
    read_record_with_options(reader, options)?.description(options.fit_model)
}

/// Read a TrafoXML document without interpreting the model name.
///
/// This is the surface the source keeps private. It is public here because it
/// is the only way to see the model name of a document whose model this port
/// cannot fit, and because writing back what was read requires it.
///
/// # Errors
///
/// As [`read_with_options`], except that no model is resolved or fitted, so the
/// unknown-model and fitting errors cannot occur.
pub fn read_record_with_options(
    reader: impl BufRead,
    options: &ReadOptions,
) -> Result<TransformationRecord> {
    options.validate()?;
    let mut work = options.max_work;
    let mut bytes = options.max_payload_bytes;
    let root = xml::parse_xml_with_budget(reader, &options.xml(), 8, None, &mut work, &mut bytes)?;
    if root.name != "TrafoXML" {
        return Err(bad("expected TrafoXML root"));
    }
    root.check(
        &["version", "xmlns:xsi", "xsi:noNamespaceSchemaLocation"],
        &["Transformation"],
    )?;
    if root.optional("xsi:noNamespaceSchemaLocation").is_some()
        && root.optional("xmlns:xsi").is_none()
    {
        return Err(bad("unbound xsi attribute prefix"));
    }
    if root
        .optional("xmlns:xsi")
        .is_some_and(|value| value != "http://www.w3.org/2001/XMLSchema-instance")
    {
        return Err(bad("invalid xsi namespace"));
    }
    // The source compares the file version numerically and only warns when the
    // file is newer than the parser. This rejects, because "undefined program
    // behavior" — the source's own words — is not an outcome a parser may pick.
    let version = root.optional("version").unwrap_or("1.0");
    let numeric = version
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
        .ok_or_else(|| bad(format!("invalid TrafoXML version {version:?}")))?;
    if !(1.0..=1.1).contains(&numeric) {
        return Err(unsupported(format!(
            "TrafoXML version {version} is outside 1.0 through {VERSION}"
        )));
    }

    let mut transformations = root
        .children
        .iter()
        .filter(|node| node.name == "Transformation");
    let transformation = transformations
        .next()
        .ok_or_else(|| bad("TrafoXML requires a Transformation"))?;
    if transformations.next().is_some() {
        return Err(bad("TrafoXML holds more than one Transformation"));
    }
    transformation.check(&["name"], &["Param", "Pairs"])?;
    let model_type = transformation.get("name")?.to_owned();

    let mut parameters: BTreeMap<String, ParamValue> = BTreeMap::new();
    let mut data = Vec::new();
    let mut pairs_seen = false;
    for node in &transformation.children {
        match node.name.as_str() {
            "Param" => {
                node.check(&["type", "name", "value"], &[])?;
                let name = node.get("name")?;
                let raw = node.get("value")?;
                let value = match node.get("type")? {
                    "int" => ParamValue::Integer(
                        raw.trim()
                            .parse::<i64>()
                            .map_err(|_| bad(format!("invalid int parameter {name}={raw:?}")))?,
                    ),
                    "float" => ParamValue::Float(finite(
                        raw.trim()
                            .parse::<f64>()
                            .map_err(|_| bad(format!("invalid float parameter {name}={raw:?}")))?,
                        name,
                    )?),
                    "string" => ParamValue::String(raw.to_owned()),
                    other => {
                        return Err(unsupported(format!(
                            "TrafoXML parameter type {other:?} for {name}"
                        )));
                    }
                };
                if parameters.insert(name.to_owned(), value).is_some() {
                    return Err(bad(format!("duplicate TrafoXML parameter {name}")));
                }
            }
            "Pairs" => {
                if pairs_seen {
                    return Err(bad("TrafoXML holds more than one Pairs list"));
                }
                pairs_seen = true;
                node.check(&["count"], &["Pair"])?;
                // The source only `reserve`s this count and never compares it
                // with the number of children, so a disagreement is accepted
                // here too; it is checked against the ceiling before reserving.
                let declared = node
                    .get("count")?
                    .trim()
                    .parse::<usize>()
                    .map_err(|_| bad("invalid Pairs count"))?;
                if declared > options.max_data_points
                    || node.children.len() > options.max_data_points
                {
                    return Err(bad("TrafoXML pair count exceeds its limit"));
                }
                data.reserve_exact(node.children.len());
                for pair in &node.children {
                    pair.check(&["from", "to", "note"], &[])?;
                    let from = finite(
                        pair.get("from")?
                            .trim()
                            .parse::<f64>()
                            .map_err(|_| bad("invalid Pair from"))?,
                        "pair source",
                    )?;
                    let to = finite(
                        pair.get("to")?
                            .trim()
                            .parse::<f64>()
                            .map_err(|_| bad("invalid Pair to"))?,
                        "pair target",
                    )?;
                    data.push(DataPoint::with_note(
                        from,
                        to,
                        pair.optional("note").unwrap_or(""),
                    ));
                }
            }
            other => return Err(unsupported(format!("{other} inside Transformation"))),
        }
    }
    Ok(TransformationRecord {
        model_type,
        parameters,
        data,
    })
}

/// Read a TrafoXML file from `path`, decompressing gzip or bzip2 content.
///
/// # Errors
///
/// As [`read_with_options`], plus [`Error::Io`] when the file cannot be opened.
pub fn load(path: impl AsRef<Path>) -> Result<TransformationDescription> {
    load_with_options(path, &ReadOptions::default())
}

/// Read a TrafoXML file under explicit limits.
///
/// # Errors
///
/// As [`load`].
pub fn load_with_options(
    path: impl AsRef<Path>,
    options: &ReadOptions,
) -> Result<TransformationDescription> {
    read_with_options(super::path_io::open(path.as_ref())?, options)
}

/// Read a TrafoXML file without interpreting the model name.
///
/// # Errors
///
/// As [`read_record_with_options`], plus [`Error::Io`] when the file cannot be
/// opened.
pub fn load_record_with_options(
    path: impl AsRef<Path>,
    options: &ReadOptions,
) -> Result<TransformationRecord> {
    read_record_with_options(super::path_io::open(path.as_ref())?, options)
}

/// Serialise a description as TrafoXML 1.1.
///
/// # Errors
///
/// As [`write_with_options`].
pub fn write(writer: impl Write, description: &TransformationDescription) -> Result<()> {
    write_with_options(writer, description, &WriteOptions::default())
}

/// Serialise a description under explicit limits.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] for a non-finite coordinate or a lowess
/// iteration count that does not fit `i64`, [`Error::Parse`] when the byte or
/// element ceiling is reached, and [`Error::Io`] from the writer.
pub fn write_with_options(
    writer: impl Write,
    description: &TransformationDescription,
    options: &WriteOptions,
) -> Result<()> {
    write_record_with_options(
        writer,
        &TransformationRecord::from_description(description)?,
        options,
    )
}

/// Serialise a record as TrafoXML 1.1.
///
/// Parameters are written in the source's three types: `int` for an integer,
/// `float` for a float and `string` for a string. A list value is written as
/// `string` carrying its bracketed text, which is what the source does. An
/// empty value is skipped, as in the source.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `model_type` is empty — the source
/// throws `Exception::IllegalArgument` with "will not write a transformation
/// with empty name" — or when a parameter name is not an XML name;
/// [`Error::Parse`] when a ceiling is reached; [`Error::Io`] from the writer.
pub fn write_record_with_options(
    mut writer: impl Write,
    record: &TransformationRecord,
    options: &WriteOptions,
) -> Result<()> {
    if options.max_xml_bytes == 0 || options.max_records == 0 {
        return Err(invalid("TrafoXML limits must be positive"));
    }
    if record.model_type.is_empty() {
        return Err(invalid("will not write a transformation with empty name"));
    }
    let mut root = Node::new("TrafoXML");
    root.attr("version", VERSION);
    root.attr("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    root.attr("xsi:noNamespaceSchemaLocation", SCHEMA_LOCATION);
    let mut transformation = Node::new("Transformation");
    transformation.attr("name", &record.model_type);
    for (name, value) in &record.parameters {
        let (kind, rendered) = match value {
            ParamValue::Empty => continue,
            ParamValue::Integer(number) => ("int", number.to_string()),
            ParamValue::Float(number) => ("float", text(finite(*number, name)?)),
            ParamValue::String(_)
            | ParamValue::StringList(_)
            | ParamValue::IntegerList(_)
            | ParamValue::FloatList(_) => ("string", value.to_text(true)?),
        };
        let mut node = Node::new("Param");
        node.attr("type", kind);
        node.attr("name", name);
        node.attr("value", rendered);
        transformation.children.push(node);
    }
    if !record.data.is_empty() {
        let mut pairs = Node::new("Pairs");
        pairs.attr("count", record.data.len());
        for point in &record.data {
            let mut node = Node::new("Pair");
            node.attr("from", text(finite(point.x, "pair source")?));
            node.attr("to", text(finite(point.y, "pair target")?));
            if !point.note.is_empty() {
                node.attr("note", &point.note);
            }
            pairs.children.push(node);
        }
        transformation.children.push(pairs);
    }
    root.children.push(transformation);
    let bytes = xml::render(
        &root,
        &xml::WriteOptions {
            max_xml_bytes: options.max_xml_bytes,
            max_records: options.max_records,
        },
    )?;
    writer.write_all(&bytes)?;
    writer.flush()?;
    Ok(())
}

/// Write a description to `path`, publishing the file only once complete.
///
/// # Errors
///
/// As [`write_with_options`], plus [`Error::Io`] when the output cannot be
/// created; the source throws `Exception::UnableToCreateFile` there.
pub fn store(path: impl AsRef<Path>, description: &TransformationDescription) -> Result<()> {
    store_with_options(path, description, &WriteOptions::default())
}

/// Write a description to `path` under explicit limits.
///
/// # Errors
///
/// As [`store`].
pub fn store_with_options(
    path: impl AsRef<Path>,
    description: &TransformationDescription,
    options: &WriteOptions,
) -> Result<()> {
    let record = TransformationRecord::from_description(description)?;
    store_record_with_options(path, &record, options)
}

/// Write a record to `path` under explicit limits.
///
/// # Errors
///
/// As [`write_record_with_options`], plus [`Error::Io`] when the output cannot
/// be created.
pub fn store_record_with_options(
    path: impl AsRef<Path>,
    record: &TransformationRecord,
    options: &WriteOptions,
) -> Result<()> {
    let mut bytes = Vec::new();
    write_record_with_options(&mut bytes, record, options)?;
    super::path_io::store(path.as_ref(), &bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = concat!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",
        "<TrafoXML version=\"1.0\" xsi:noNamespaceSchemaLocation=\"x\"",
        " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n"
    );

    fn document(body: &str) -> String {
        format!("{HEADER}{body}</TrafoXML>\n")
    }

    #[test]
    fn unknown_parameter_type_is_an_error_where_the_source_only_logs() {
        let text = document(
            "<Transformation name=\"linear\">\
             <Param type=\"bogus\" name=\"slope\" value=\"1\"/>\
             </Transformation>\n",
        );
        let error = read(text.as_bytes()).unwrap_err();
        assert!(matches!(error, Error::Unsupported(_)), "{error:?}");
    }

    #[test]
    fn notes_survive_a_round_trip_including_non_ascii_text() {
        let record = TransformationRecord {
            model_type: "none".into(),
            parameters: BTreeMap::new(),
            data: vec![DataPoint::with_note(1.0, 2.0, "日本語 <&\"'>")],
        };
        let mut bytes = Vec::new();
        write_record_with_options(&mut bytes, &record, &WriteOptions::default()).unwrap();
        let back = read_record_with_options(bytes.as_slice(), &ReadOptions::default()).unwrap();
        assert_eq!(back, record);
    }

    #[test]
    fn pair_ceiling_is_checked_before_any_pair_is_stored() {
        let text = document(
            "<Transformation name=\"none\"><Pairs count=\"2\">\
             <Pair from=\"1\" to=\"2\"/><Pair from=\"3\" to=\"4\"/>\
             </Pairs></Transformation>\n",
        );
        let options = ReadOptions {
            max_data_points: 1,
            ..Default::default()
        };
        assert!(read_record_with_options(text.as_bytes(), &options).is_err());
    }
}

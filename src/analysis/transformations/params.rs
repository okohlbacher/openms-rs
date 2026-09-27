// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The `Param` surface of the transformation models: each model's
//! `getDefaultParameters`, the constructors' reading of their `Param`, the
//! model parameters `getParameters` reports, and
//! `TransformationModelDefaults::getDefaults`.
//!
//! The typed option structs stay the primary API. `defaults`, `from_param` and
//! `to_param` translate between them and the source's parameter layout, as
//! `processing::peak_picking::PeakPickerHiRes` does for an algorithm. The
//! models are not `DefaultParamHandler`s in the source: a model constructor
//! fills missing values from its defaults (`Param::setDefaults`) and neither
//! checks restrictions nor warns about unknown names, and neither does
//! `from_param`. Only what a value must be converted to is checked, where the
//! source's `ParamValue` conversions throw `Exception::ConversionError`.

use super::lowess::{CrossValidationMetric, LowessAutoSpan};
use super::{
    BSplineExtrapolation, BSplineOptions, CoordinateWeight, Extrapolation, Interpolation,
    InterpolationOptions, LinearCoefficients, LinearOptions, LowessOptions, ModelConfig,
    ModelWeights, WeightFunction, bad,
};
use crate::Result;
use crate::data_structures::list;
use crate::param::{Param, ParamEntry, ParamNode, ParamValue};
use crate::processing::spline::BoundaryCondition;

/// The model names `TransformationDescription::getModelTypes` lists, in its
/// order. `none` and `identity` do not count as models there.
pub const MODEL_TYPES: [&str; 4] = ["linear", "b_spline", "interpolated", "lowess"];

/// Valid names of the x weighting (source `getValidXWeights`); `x` disables it.
pub const VALID_X_WEIGHTS: [&str; 4] = ["1/x", "1/x2", "ln(x)", "x"];

/// Valid names of the y weighting (source `getValidYWeights`); `y` disables it.
pub const VALID_Y_WEIGHTS: [&str; 4] = ["1/y", "1/y2", "ln(y)", "y"];

fn entry(name: &str, value: ParamValue, description: &str) -> ParamEntry {
    ParamEntry {
        name: name.into(),
        description: description.into(),
        value,
        ..ParamEntry::default()
    }
}

fn float_entry(name: &str, value: f64, description: &str) -> ParamEntry {
    entry(name, ParamValue::Float(value), description)
}

fn int_entry(name: &str, value: i64, description: &str) -> ParamEntry {
    entry(name, ParamValue::Integer(value), description)
}

fn string_entry(name: &str, value: &str, description: &str, valid: &[&str]) -> ParamEntry {
    let mut e = entry(name, ParamValue::String(value.into()), description);
    e.valid_strings = valid.iter().map(|v| (*v).to_owned()).collect();
    e
}

fn tree(entries: Vec<ParamEntry>) -> Result<Param> {
    Param::from_root(ParamNode {
        name: String::new(),
        description: String::new(),
        entries,
        nodes: Vec::new(),
    })
}

/// `params` with missing values filled from `defaults` (source
/// `Param::setDefaults`, which every model constructor calls).
fn with_defaults(params: &Param, defaults: &Param) -> Result<Param> {
    let mut filled = params.checked_clone()?;
    filled.set_defaults(defaults, "", false)?;
    Ok(filled)
}

fn get<'a>(params: &'a Param, key: &str) -> Result<Option<&'a ParamValue>> {
    if params.exists(key)? {
        Ok(Some(params.value(key)?))
    } else {
        Ok(None)
    }
}

/// `(double)param.getValue(key)`: an integer converts, a string or an empty
/// value does not.
fn float(params: &Param, key: &str) -> Result<f64> {
    params
        .value(key)?
        .to_f64()
        .map_err(|_| bad(&format!("transformation parameter '{key}' is not a number")))
}

/// `(int)param.getValue(key)`: only an integer converts.
fn int(params: &Param, key: &str) -> Result<i64> {
    params.value(key)?.to_i64().map_err(|_| {
        bad(&format!(
            "transformation parameter '{key}' is not an integer"
        ))
    })
}

/// `(std::string)param.getValue(key)`: only a string converts.
fn string<'a>(params: &'a Param, key: &str) -> Result<&'a str> {
    params
        .value(key)?
        .as_str()
        .map_err(|_| bad(&format!("transformation parameter '{key}' is not a string")))
}

impl WeightFunction {
    /// The source name of this weighting on the x axis (`x`, `1/x`, `1/x2`,
    /// `ln(x)`).
    pub fn x_name(self) -> &'static str {
        match self {
            Self::Identity => "x",
            Self::Reciprocal => "1/x",
            Self::ReciprocalSquared => "1/x2",
            Self::Log => "ln(x)",
        }
    }

    /// The source name of this weighting on the y axis (`y`, `1/y`, `1/y2`,
    /// `ln(y)`).
    pub fn y_name(self) -> &'static str {
        match self {
            Self::Identity => "y",
            Self::Reciprocal => "1/y",
            Self::ReciprocalSquared => "1/y2",
            Self::Log => "ln(y)",
        }
    }

    /// The weighting a source name selects on the x axis (`y == false`) or
    /// the y axis, as the base `TransformationModel` constructor reads it.
    ///
    /// An empty name means no weighting, which TrafoXML files written before
    /// OpenMS 3.0 store for an unweighted model.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for a name that is not in
    /// [`VALID_X_WEIGHTS`] (respectively [`VALID_Y_WEIGHTS`]), the source's
    /// `Exception::InvalidParameter`. The axis matters: `ln(y)` is not a valid
    /// x weight.
    pub fn from_name(name: &str, y: bool) -> Result<Self> {
        let valid = if y { VALID_Y_WEIGHTS } else { VALID_X_WEIGHTS };
        if name.is_empty() {
            return Ok(Self::Identity);
        }
        match valid.iter().position(|v| *v == name) {
            Some(0) => Ok(Self::Reciprocal),
            Some(1) => Ok(Self::ReciprocalSquared),
            Some(2) => Ok(Self::Log),
            Some(_) => Ok(Self::Identity),
            None => Err(bad(&format!(
                "Value '{name}' is not a valid weight parameter for {} values.",
                if y { "y" } else { "x" }
            ))),
        }
    }
}

impl ModelWeights {
    /// Read the weighting as the base `TransformationModel(DataPoints,
    /// Param)` constructor does: `x_weight` and `y_weight` (missing or empty:
    /// no weighting) and the four datum bounds (missing: `1e-15` and `1e15`).
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for a weight name that is not
    /// valid for its axis (the source's `Exception::InvalidParameter`), a
    /// non-string weight, or a non-numeric bound (the source's
    /// `ParamValue` conversion throws `Exception::ConversionError`).
    pub fn from_param(params: &Param) -> Result<Self> {
        let weight = |axis: char| -> Result<CoordinateWeight> {
            let y = axis == 'y';
            let name = match get(params, &format!("{axis}_weight"))? {
                None => "",
                Some(value) => value
                    .as_str()
                    .map_err(|_| bad(&format!("'{axis}_weight' is not a string")))?,
            };
            let bound = |key: &str, fallback: f64| -> Result<f64> {
                if params.exists(key)? {
                    float(params, key)
                } else {
                    Ok(fallback)
                }
            };
            Ok(CoordinateWeight {
                function: WeightFunction::from_name(name, y)?,
                min: bound(&format!("{axis}_datum_min"), 1e-15)?,
                max: bound(&format!("{axis}_datum_max"), 1e15)?,
            })
        };
        Ok(Self {
            x: weight('x')?,
            y: weight('y')?,
        })
    }
}

/// Whether `weight` is one of `valid` (source `checkValidWeight`, which also
/// logs "weight ... is not supported." for an invalid name).
pub fn check_valid_weight(weight: &str, valid: &[&str]) -> bool {
    valid.contains(&weight)
}

/// Clamp a datum into `[datum_min, datum_max]` (source `checkDatumRange`).
///
/// A value at or above the maximum becomes the maximum; otherwise a value at
/// or below the minimum becomes the minimum. The source logs both cases at
/// info level. With crossed bounds the maximum test wins, as in the source.
pub fn check_datum_range(datum: f64, datum_min: f64, datum_max: f64) -> f64 {
    if datum >= datum_max {
        datum_max
    } else if datum <= datum_min {
        datum_min
    } else {
        datum
    }
}

/// Weight a datum by a source weight name (source `weightDatum`).
///
/// `ln(x)`/`ln(y)` take the natural logarithm, `1/x`/`1/y` the reciprocal of
/// the absolute value, `1/x2`/`1/y2` the reciprocal of the square, and
/// `x`/`y` leave the datum unchanged. Any other name leaves it unchanged too;
/// the source logs "weight ... not supported" then. The arithmetic is not
/// checked, as in the source: `ln` of a non-positive datum is NaN or `-inf`.
pub fn weight_datum(datum: f64, weight: &str) -> f64 {
    match weight {
        "ln(x)" | "ln(y)" => datum.ln(),
        "1/x" | "1/y" => 1.0 / datum.abs(),
        "1/x2" | "1/y2" => 1.0 / (datum * datum),
        _ => datum,
    }
}

/// Undo [`weight_datum`] (source `unWeightDatum`): `exp` for the logarithm,
/// `1/|d|` for the reciprocal, `sqrt(1/|d|)` for the reciprocal square, and
/// unchanged otherwise.
pub fn unweight_datum(datum: f64, weight: &str) -> f64 {
    match weight {
        "ln(x)" | "ln(y)" => datum.exp(),
        "1/x" | "1/y" => 1.0 / datum.abs(),
        "1/x2" | "1/y2" => (1.0 / datum.abs()).sqrt(),
        _ => datum,
    }
}

/// Every model's parameters in one tree, with a `type` selection (source
/// `TransformationModelDefaults::getDefaults`), as the alignment tools and
/// algorithms offer them.
///
/// `type` holds `default_model`. Its valid strings are `linear`, `b_spline`,
/// `lowess` and `interpolated`, with `default_model` prepended when it is not
/// one of them (the tools use `none`). The sections `linear:`, `b_spline:`,
/// `lowess:` and `interpolated:` hold each model's defaults and are described
/// as `Parameters for '<model>' model`.
///
/// # Errors
///
/// Returns [`crate::Error::InvalidValue`] when `default_model` contains a
/// comma, which a valid-string list cannot hold.
pub fn transformation_model_defaults(default_model: &str) -> Result<Param> {
    let mut models = vec!["linear", "b_spline", "lowess", "interpolated"];
    if !models.contains(&default_model) {
        models.insert(0, default_model);
    }
    let mut params = Param::new();
    params.set_value(
        "type",
        ParamValue::String(default_model.into()),
        "Type of model",
        &[],
    )?;
    params.set_valid_strings(
        "type",
        &models.iter().map(|m| (*m).to_owned()).collect::<Vec<_>>(),
    )?;
    for (name, defaults) in [
        ("linear", LinearOptions::defaults()?),
        ("b_spline", BSplineOptions::defaults()?),
        ("lowess", LowessOptions::defaults()?),
        ("interpolated", InterpolationOptions::defaults()?),
    ] {
        params.insert(&format!("{name}:"), &defaults)?;
        params.set_section_description(name, &format!("Parameters for '{name}' model"))?;
    }
    Ok(params)
}

impl LinearOptions {
    /// The source defaults (`TransformationModelLinear::getDefaultParameters`):
    /// `symmetric_regression`, `x_weight`, `y_weight`, `x_datum_min`,
    /// `x_datum_max`, `y_datum_min` and `y_datum_max`, with their descriptions
    /// and valid strings, in declaration order.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which the constant tree cannot
    /// trigger.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// Read the options the way the linear model's constructor reads its
    /// `Param`.
    ///
    /// A missing weight means no weighting and a missing bound takes the
    /// default (`1e-15`, `1e15`). `slope` and `intercept` become
    /// [`LinearOptions::coefficients`] only when both exist; like the source,
    /// they are then used only if there are no data points to fit.
    /// `symmetric_regression` is not read, because the source reads it into a
    /// member it never uses.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for an invalid weight name (see
    /// [`WeightFunction::from_name`]) or a value that does not convert: a
    /// non-string weight, or a non-numeric bound, slope or intercept.
    pub fn from_param(params: &Param) -> Result<Self> {
        let weights = ModelWeights::from_param(params)?;
        let coefficients = if params.exists("slope")? && params.exists("intercept")? {
            Some(LinearCoefficients {
                slope: float(params, "slope")?,
                intercept: float(params, "intercept")?,
            })
        } else {
            None
        };
        Ok(Self {
            x_weight: weights.x,
            y_weight: weights.y,
            coefficients,
            ..Self::default()
        })
    }

    /// The options in the layout of [`LinearOptions::defaults`], carrying the
    /// current weights and bounds, followed by `slope` and `intercept` when
    /// explicit coefficients are set.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut entries = vec![
            string_entry(
                "symmetric_regression",
                "false",
                "Perform linear regression on 'y - x' vs. 'y + x', instead of on 'y' vs. 'x'.",
                &["true", "false"],
            ),
            string_entry(
                "x_weight",
                self.x_weight.function.x_name(),
                "Weight x values",
                &VALID_X_WEIGHTS,
            ),
            string_entry(
                "y_weight",
                self.y_weight.function.y_name(),
                "Weight y values",
                &VALID_Y_WEIGHTS,
            ),
            float_entry("x_datum_min", self.x_weight.min, "Minimum x value"),
            float_entry("x_datum_max", self.x_weight.max, "Maximum x value"),
            float_entry("y_datum_min", self.y_weight.min, "Minimum y value"),
            float_entry("y_datum_max", self.y_weight.max, "Maximum y value"),
        ];
        if let Some(coefficients) = self.coefficients {
            entries.push(float_entry("slope", coefficients.slope, ""));
            entries.push(float_entry("intercept", coefficients.intercept, ""));
        }
        tree(entries)
    }
}

impl Interpolation {
    /// The source parameter value naming this interpolation.
    pub fn name(self) -> &'static str {
        match self {
            Self::Linear => "linear",
            Self::CubicSpline => "cspline",
            Self::Akima => "akima",
        }
    }

    /// The interpolation a source name selects.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for any other name, the source's
    /// `Exception::IllegalArgument` ("unknown/unsupported interpolation
    /// type").
    pub fn from_name(name: &str) -> Result<Self> {
        match name {
            "linear" => Ok(Self::Linear),
            "cspline" => Ok(Self::CubicSpline),
            "akima" => Ok(Self::Akima),
            _ => Err(bad(&format!(
                "unknown/unsupported interpolation type '{name}'"
            ))),
        }
    }
}

impl Extrapolation {
    /// The source parameter value naming this extrapolation.
    pub fn name(self) -> &'static str {
        match self {
            Self::TwoPointLinear => "two-point-linear",
            Self::FourPointLinear => "four-point-linear",
            Self::GlobalLinear => "global-linear",
        }
    }

    /// The extrapolation a source name selects.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for any other name, the source's
    /// `Exception::IllegalArgument` ("unknown/unsupported extrapolation
    /// type").
    pub fn from_name(name: &str) -> Result<Self> {
        match name {
            "two-point-linear" => Ok(Self::TwoPointLinear),
            "four-point-linear" => Ok(Self::FourPointLinear),
            "global-linear" => Ok(Self::GlobalLinear),
            _ => Err(bad(&format!(
                "unknown/unsupported extrapolation type '{name}'"
            ))),
        }
    }
}

const INTERPOLATED_EXTRAPOLATION_DESCRIPTION: &str = "Type of extrapolation to apply: two-point-linear: use the first and last data point to build a single linear model, four-point-linear: build two linear models on both ends using the first two / last two points, global-linear: use all points to build a single linear model. Note that global-linear may not be continuous at the border.";

impl InterpolationOptions {
    /// The source defaults
    /// (`TransformationModelInterpolated::getDefaultParameters`):
    /// `interpolation_type` (`cspline`) and `extrapolation_type`
    /// (`two-point-linear`) with their valid strings.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// Read the options the way the interpolated model's constructor does,
    /// missing values taking their defaults.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for a non-string or unknown
    /// type name (the source's `Exception::IllegalArgument`).
    pub fn from_param(params: &Param) -> Result<Self> {
        let params = with_defaults(params, &Self::defaults()?)?;
        Ok(Self {
            interpolation: Interpolation::from_name(string(&params, "interpolation_type")?)?,
            extrapolation: Extrapolation::from_name(string(&params, "extrapolation_type")?)?,
            ..Self::default()
        })
    }

    /// The options in the layout of [`InterpolationOptions::defaults`].
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        tree(vec![
            string_entry(
                "interpolation_type",
                self.interpolation.name(),
                "Type of interpolation to apply.",
                &["linear", "cspline", "akima"],
            ),
            string_entry(
                "extrapolation_type",
                self.extrapolation.name(),
                INTERPOLATED_EXTRAPOLATION_DESCRIPTION,
                &["two-point-linear", "four-point-linear", "global-linear"],
            ),
        ])
    }
}

impl BSplineOptions {
    /// The source defaults (`TransformationModelBSpline::getDefaultParameters`):
    /// `wavelength` (0, minimum 0), `num_nodes` (5, minimum 0), `extrapolate`
    /// (`linear`) and `boundary_condition` (2, from 0 to 2), with their
    /// descriptions and restrictions.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// Read the options the way the B-spline model's constructor does,
    /// missing values taking their defaults.
    ///
    /// `extrapolate` follows [`BSplineExtrapolation::from_name`], so an
    /// unknown name means linear extrapolation, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for a value that does not
    /// convert, a negative `num_nodes` (the source's conversion to `Size`
    /// throws `Exception::ConversionError`), or a `boundary_condition` other
    /// than 0, 1 or 2, which the source casts to its enum unchecked.
    pub fn from_param(params: &Param) -> Result<Self> {
        let params = with_defaults(params, &Self::defaults()?)?;
        let num_nodes = usize::try_from(int(&params, "num_nodes")?)
            .map_err(|_| bad("B-spline 'num_nodes' must not be negative"))?;
        let boundary_condition = match int(&params, "boundary_condition")? {
            0 => BoundaryCondition::ZeroEndpoints,
            1 => BoundaryCondition::ZeroFirst,
            2 => BoundaryCondition::ZeroSecond,
            _ => return Err(bad("B-spline 'boundary_condition' must be 0, 1 or 2")),
        };
        Ok(Self {
            wavelength: float(&params, "wavelength")?,
            num_nodes,
            extrapolation: BSplineExtrapolation::from_name(string(&params, "extrapolate")?),
            boundary_condition,
            ..Self::default()
        })
    }

    /// The options in the layout of [`BSplineOptions::defaults`].
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] when `num_nodes` exceeds the
    /// source's 64-bit integer parameter.
    pub fn to_param(&self) -> Result<Param> {
        let mut wavelength = float_entry(
            "wavelength",
            self.wavelength,
            "Determines the amount of smoothing by setting the number of nodes for the B-spline. The number is chosen so that the spline approximates a low-pass filter with this cutoff wavelength. The wavelength is given in the same units as the data; a higher value means more smoothing. '0' sets the number of nodes to twice the number of input points.",
        );
        wavelength.min_float = 0.0;
        let mut nodes = int_entry(
            "num_nodes",
            i64::try_from(self.num_nodes).map_err(|_| bad("B-spline 'num_nodes' is too large"))?,
            "Number of nodes for B-spline fitting. Overrides 'wavelength' if set (to two or greater). A lower value means more smoothing.",
        );
        nodes.min_int = 0;
        let extrapolate = string_entry(
            "extrapolate",
            self.extrapolation.name(),
            "Method to use for extrapolation beyond the original data range. 'linear': Linear extrapolation using the slope of the B-spline at the corresponding endpoint. 'b_spline': Use the B-spline (as for interpolation). 'constant': Use the constant value of the B-spline at the corresponding endpoint. 'global_linear': Use a linear fit through the data (which will most probably introduce discontinuities at the ends of the data range).",
            &["linear", "b_spline", "constant", "global_linear"],
        );
        let mut boundary = int_entry(
            "boundary_condition",
            self.boundary_condition as i64,
            "Boundary condition at B-spline endpoints: 0 (value zero), 1 (first derivative zero) or 2 (second derivative zero)",
        );
        boundary.min_int = 0;
        boundary.max_int = 2;
        tree(vec![wavelength, nodes, extrapolate, boundary])
    }
}

impl LowessOptions {
    /// The source defaults (`TransformationModelLowess::getDefaultParameters`),
    /// in declaration order: `span`, the automatic span selection
    /// (`auto_span`, `auto_span_min`, `auto_span_max`, `auto_min_neighbors`,
    /// `auto_k_folds`, `auto_metric`, `auto_span_grid`), `num_iterations`,
    /// `delta` (-1: automatic), `interpolation_type` (`cspline`) and
    /// `extrapolation_type` (`four-point-linear`), with their descriptions and
    /// restrictions.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param_with(None)
    }

    /// Read the options the way the LOWESS model's constructor does, missing
    /// values taking their defaults, together with the automatic span
    /// selection when `auto_span` is `true`.
    ///
    /// A negative `delta` selects the automatic distance (`None`). The source
    /// converts every automatic-selection parameter whether or not the
    /// selection is enabled, so a wrong type or an unknown `auto_metric` is an
    /// error even with `auto_span` set to `false`; so it is here.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for a value that does not
    /// convert, a negative `num_iterations`, an `auto_span` other than `true`
    /// or `false` (the source's `ParamValue::toBool` throws), an unknown
    /// `auto_metric`, an `auto_span_grid` entry that is not a number (the
    /// source's `ListUtils::create` throws `Exception::ConversionError`), or
    /// an unknown interpolation or extrapolation name.
    pub fn from_param(params: &Param) -> Result<(Self, Option<LowessAutoSpan>)> {
        let params = with_defaults(params, &Self::defaults()?)?;
        let delta = float(&params, "delta")?;
        let options = Self {
            span: float(&params, "span")?,
            iterations: usize::try_from(int(&params, "num_iterations")?)
                .map_err(|_| bad("LOWESS 'num_iterations' must not be negative"))?,
            delta: if delta < 0.0 { None } else { Some(delta) },
            interpolation: Interpolation::from_name(string(&params, "interpolation_type")?)?,
            extrapolation: Extrapolation::from_name(string(&params, "extrapolation_type")?)?,
            ..Self::default()
        };
        let enabled = params
            .value("auto_span")?
            .to_bool()
            .map_err(|_| bad("LOWESS 'auto_span' must be 'true' or 'false'"))?;
        let narrow = |key: &str| -> Result<i32> {
            i32::try_from(int(&params, key)?)
                .map_err(|_| bad(&format!("LOWESS '{key}' does not fit a 32-bit integer")))
        };
        let grid_text = string(&params, "auto_span_grid")?;
        let grid = if grid_text.is_empty() {
            Vec::new()
        } else {
            list::create::<f64>(grid_text, b',')
                .map_err(|_| bad(&format!("Could not convert '{grid_text}' to span values")))?
        };
        let auto = LowessAutoSpan {
            span_min: float(&params, "auto_span_min")?,
            span_max: float(&params, "auto_span_max")?,
            min_neighbors: narrow("auto_min_neighbors")?,
            k_folds: narrow("auto_k_folds")?,
            metric: CrossValidationMetric::from_name(string(&params, "auto_metric")?)?,
            grid,
            ..LowessAutoSpan::default()
        };
        Ok((options, enabled.then_some(auto)))
    }

    /// The options in the layout of [`LowessOptions::defaults`]; the
    /// automatic-selection parameters take their defaults and `auto_span` is
    /// `false`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] when `iterations` exceeds the
    /// source's 64-bit integer parameter.
    pub fn to_param(&self) -> Result<Param> {
        self.to_param_with(None)
    }

    /// As [`LowessOptions::to_param`], with the automatic-selection
    /// parameters taken from `auto` and `auto_span` set to `true` when it is
    /// given.
    ///
    /// # Errors
    ///
    /// As [`LowessOptions::to_param`].
    pub fn to_param_with(&self, auto: Option<&LowessAutoSpan>) -> Result<Param> {
        let defaults = LowessAutoSpan::default();
        let selection = auto.unwrap_or(&defaults);
        let mut span = float_entry(
            "span",
            self.span,
            "Fraction of datapoints (f) to use for each local regression (determines the amount of smoothing). Choosing this parameter in the range .2 to .8 usually results in a good fit.",
        );
        span.min_float = 0.0;
        span.max_float = 1.0;
        let auto_span = string_entry(
            "auto_span",
            if auto.is_some() { "true" } else { "false" },
            "If true, or if 'span' is 0, automatically select LOWESS span by cross-validation.",
            &["true", "false"],
        );
        let mut span_min = float_entry(
            "auto_span_min",
            selection.span_min,
            "Lower bound for auto-selected span.",
        );
        span_min.min_float = 0.001;
        let mut span_max = float_entry(
            "auto_span_max",
            selection.span_max,
            "Upper bound for auto-selected span.",
        );
        span_max.max_float = 0.99;
        let mut neighbors = int_entry(
            "auto_min_neighbors",
            i64::from(selection.min_neighbors),
            "Minimum number of neighbors (span*n) enforced in auto mode.",
        );
        neighbors.min_int = 3;
        let mut folds = int_entry(
            "auto_k_folds",
            i64::from(selection.k_folds),
            "K-folds for CV when n>50 (else LOO is used).",
        );
        folds.min_int = 2;
        let metric = string_entry(
            "auto_metric",
            selection.metric.name(),
            "Metric for CV selection: one of {'p90','p95','p99','rmse','mae'}.",
            &["p90", "p95", "p99", "rmse", "mae"],
        );
        let grid = entry(
            "auto_span_grid",
            ParamValue::String(
                selection
                    .grid
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
            ),
            "Optional explicit grid of span candidates in (0,1]. Comma-separated list, e.g. '0.2,0.3,0.5'.  If empty, a default grid is used.",
        );
        let mut iterations = int_entry(
            "num_iterations",
            i64::try_from(self.iterations)
                .map_err(|_| bad("LOWESS 'num_iterations' is too large"))?,
            "Number of robustifying iterations for lowess fitting.",
        );
        iterations.min_int = 0;
        let delta = float_entry(
            "delta",
            self.delta.unwrap_or(-1.0),
            "Nonnegative parameter which may be used to save computations (recommended value is 0.01 of the range of the input, e.g. for data ranging from 1000 seconds to 2000 seconds, it could be set to 10). Setting a negative value will automatically do this.",
        );
        let interpolation = string_entry(
            "interpolation_type",
            self.interpolation.name(),
            "Method to use for interpolation between datapoints computed by lowess. 'linear': Linear interpolation. 'cspline': Use the cubic spline for interpolation. 'akima': Use an akima spline for interpolation",
            &["linear", "cspline", "akima"],
        );
        let extrapolation = string_entry(
            "extrapolation_type",
            self.extrapolation.name(),
            "Method to use for extrapolation outside the data range. 'two-point-linear': Uses a line through the first and last point to extrapolate. 'four-point-linear': Uses a line through the first and second point to extrapolate in front and and a line through the last and second-to-last point in the end. 'global-linear': Uses a linear regression to fit a line through all data points and use it for interpolation.",
            &["two-point-linear", "four-point-linear", "global-linear"],
        );
        tree(vec![
            span,
            auto_span,
            span_min,
            span_max,
            neighbors,
            folds,
            metric,
            grid,
            iterations,
            delta,
            interpolation,
            extrapolation,
        ])
    }
}

impl ModelConfig {
    /// The configuration `TransformationDescription::fitModel(model_type,
    /// params)` fits: the model's constructor reading `params`.
    ///
    /// `none` and `identity` ignore `params`, as the source's base model does.
    /// `lowess` with `auto_span` set to `true` yields
    /// [`ModelConfig::LowessAuto`].
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for an unknown model type (the
    /// source's `Exception::IllegalArgument`, "unknown model type") and for
    /// the parameter errors of each model's `from_param`.
    pub fn from_param(model_type: &str, params: &Param) -> Result<Self> {
        match model_type {
            "none" => Ok(Self::None),
            "identity" => Ok(Self::Identity),
            "linear" => Ok(Self::Linear(LinearOptions::from_param(params)?)),
            "b_spline" => Ok(Self::BSpline(BSplineOptions::from_param(params)?)),
            "interpolated" => Ok(Self::Interpolated(InterpolationOptions::from_param(
                params,
            )?)),
            "lowess" => Ok(match LowessOptions::from_param(params)? {
                (options, Some(auto)) => Self::LowessAuto(options, auto),
                (options, None) => Self::Lowess(options),
            }),
            other => Err(bad(&format!("unknown model type '{other}'"))),
        }
    }

    /// The default parameters of a model type (the static
    /// `getDefaultParameters` of its class); empty for `none` and `identity`,
    /// whose base class clears the tree.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for an unknown model type.
    pub fn default_parameters(model_type: &str) -> Result<Param> {
        match model_type {
            "none" | "identity" => Ok(Param::new()),
            "linear" => LinearOptions::defaults(),
            "b_spline" => BSplineOptions::defaults(),
            "interpolated" => InterpolationOptions::defaults(),
            "lowess" => LowessOptions::defaults(),
            other => Err(bad(&format!("unknown model type '{other}'"))),
        }
    }

    /// The model type this configuration fits, as `getModelType` names it.
    pub fn model_type(&self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Identity => "identity",
            Self::Linear(_) => "linear",
            Self::BSpline(_) => "b_spline",
            Self::Interpolated(_) => "interpolated",
            Self::Lowess(_) | Self::LowessAuto(..) => "lowess",
        }
    }

    /// The parameters a model built from this configuration reports: the
    /// model's default layout with the configured values; empty for `none`
    /// and `identity`.
    ///
    /// # Errors
    ///
    /// As the options' `to_param`.
    pub fn to_param(&self) -> Result<Param> {
        match self {
            Self::None | Self::Identity => Ok(Param::new()),
            Self::Linear(options) => options.to_param(),
            Self::BSpline(options) => options.to_param(),
            Self::Interpolated(options) => options.to_param(),
            Self::Lowess(options) => options.to_param(),
            Self::LowessAuto(options, auto) => options.to_param_with(Some(auto)),
        }
    }
}

/// The parameters the source model keeps after construction from `params`
/// (`params_`, which `getParameters` and `getModelParameters` return).
///
/// - `none` and `identity`: empty; the base model's default constructor
///   ignores the `Param`.
/// - `linear` with explicit `slope` and `intercept` and no data: `params`
///   unchanged.
/// - `linear` otherwise: `params` with the linear defaults filled in, and the
///   fitted `slope` and `intercept` set (without descriptions, as
///   `Param::setValue` stores them).
/// - `b_spline`, `interpolated`, `lowess`: `params` with the model's defaults
///   filled in; `lowess` after an automatic selection also carries the
///   selected `span` and `auto_span` set to `false`.
pub(super) fn stored_parameters(
    model_type: &str,
    params: &Param,
    data_given: bool,
    fitted: &super::TransformationModel,
) -> Result<Param> {
    match model_type {
        "none" | "identity" => Ok(Param::new()),
        "linear" => {
            if !data_given && params.exists("slope")? && params.exists("intercept")? {
                return params.checked_clone();
            }
            let mut stored = with_defaults(params, &LinearOptions::defaults()?)?;
            if let super::TransformationModel::Linear(model) = fitted {
                let coefficients = model.coefficients();
                stored.set_value("slope", ParamValue::Float(coefficients.slope), "", &[])?;
                stored.set_value(
                    "intercept",
                    ParamValue::Float(coefficients.intercept),
                    "",
                    &[],
                )?;
            }
            Ok(stored)
        }
        "lowess" => {
            let mut stored = with_defaults(params, &LowessOptions::defaults()?)?;
            if stored.value("auto_span")?.to_bool().unwrap_or(false) {
                if let super::TransformationModel::Lowess(model) = fitted {
                    stored.set_value("span", ParamValue::Float(model.options().span), "", &[])?;
                    stored.set_value("auto_span", ParamValue::String("false".into()), "", &[])?;
                }
            }
            Ok(stored)
        }
        other => with_defaults(params, &ModelConfig::default_parameters(other)?),
    }
}

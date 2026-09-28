// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Coordinate transformations fitted to annotated anchors.
//!
//! Ports `ANALYSIS/MAPMATCHING/TransformationDescription.h`,
//! `TransformationModel.h`, `TransformationModelLinear.h`,
//! `TransformationModelInterpolated.h`, `TransformationModelBSpline.h`,
//! `TransformationModelLowess.h` and `TransformationModelDefaults.h`.
//!
//! A [`TransformationDescription`] stores data points (coordinate pairs) from
//! which a [`TransformationModel`] is estimated; applying the transformation
//! to a coordinate evaluates the model there. The source's models are a
//! virtual base class with one subclass per model; here they are the closed
//! enum [`TransformationModel`] fitted from a typed [`ModelConfig`]. The
//! source's `Param` interface is kept beside the typed one:
//! [`TransformationDescription::fit_model_with_param`] fits a model named by
//! string from a `Param` exactly as `fitModel(model_type, params)` does, and
//! [`TransformationDescription::model_parameters`] returns what
//! `getModelParameters` returns.
//!
//! Inversion of a description refits swapped anchors and is generally
//! approximate. See `docs/TRANSFORMATIONS_SUPPORT.md`.

mod akima;
mod b_spline;
mod lowess;
mod params;
pub use b_spline::{BSplineExtrapolation, BSplineModel, BSplineOptions};
pub use lowess::{
    CrossValidationMetric, DEFAULT_SPAN_GRID, LowessAutoSpan, LowessModel, LowessOptions,
    SPAN_SCORE_TIE_TOLERANCE, build_span_grid, lowess, score_residuals,
};
pub use params::{
    MODEL_TYPES, VALID_X_WEIGHTS, VALID_Y_WEIGHTS, check_datum_range, check_valid_weight,
    transformation_model_defaults, unweight_datum, weight_datum,
};

use crate::concept::fuzzy_string_comparator::format_g;
use crate::kernel::NumericRange;
use crate::math::libstdcxx;
use crate::param::{Param, ParamValue};
use crate::processing::peak_picking::CubicSpline2d;
use crate::{Error, Result};
use std::io::Write;

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("transformation arithmetic or input is not finite"))
    }
}
fn validate_points(data: &[DataPoint], max_points: usize) -> Result<()> {
    if max_points == 0 || data.len() > max_points {
        return Err(bad("transformation point limit exceeded or zero"));
    }
    for p in data {
        finite(p.x)?;
        finite(p.y)?;
    }
    Ok(())
}

/// Coordinate pair with an optional annotation (source
/// `TransformationModel::DataPoint`, whose `first` and `second` are `x` and
/// `y` here).
///
/// Ordering and equality compare `x`, then `y`, then `note`, as the source's
/// `std::tie` operators do.
#[derive(Clone, Debug, Default, PartialEq, PartialOrd)]
pub struct DataPoint {
    /// Coordinate in the source (x) space.
    pub x: f64,
    /// Coordinate in the target (y) space.
    pub y: f64,
    /// Free-text annotation, such as the peptide the pair came from.
    pub note: String,
}
impl DataPoint {
    /// A pair without annotation.
    pub fn new(x: f64, y: f64) -> Self {
        Self {
            x,
            y,
            note: String::new(),
        }
    }
    /// An annotated pair.
    pub fn with_note(x: f64, y: f64, note: impl Into<String>) -> Self {
        Self {
            x,
            y,
            note: note.into(),
        }
    }
}
impl From<(f64, f64)> for DataPoint {
    /// The source's `DataPoint(const std::pair<double, double>&)`, without a note.
    fn from((x, y): (f64, f64)) -> Self {
        Self::new(x, y)
    }
}

/// A coordinate transformation applied to x or y before a linear fit.
///
/// The source calls these "weights" (`x_weight`, `y_weight`, `weightData`);
/// they are coordinate transformations, not per-observation statistical
/// regression weights. Their source names are [`WeightFunction::x_name`] and
/// [`WeightFunction::y_name`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum WeightFunction {
    /// `x` / `y`: unchanged.
    #[default]
    Identity,
    /// `ln(x)` / `ln(y)`: the natural logarithm.
    Log,
    /// `1/x` / `1/y`: the reciprocal of the absolute value.
    Reciprocal,
    /// `1/x2` / `1/y2`: the reciprocal of the square.
    ReciprocalSquared,
}
/// A [`WeightFunction`] with the datum range training values are clamped to
/// (source `x_weight`, `x_datum_min`, `x_datum_max` and the y equivalents).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoordinateWeight {
    /// The transformation.
    pub function: WeightFunction,
    /// Lower bound of the datum range (default `1e-15`).
    pub min: f64,
    /// Upper bound of the datum range (default `1e15`).
    pub max: f64,
}
impl Default for CoordinateWeight {
    fn default() -> Self {
        Self {
            function: WeightFunction::Identity,
            min: 1e-15,
            max: 1e15,
        }
    }
}
impl CoordinateWeight {
    fn validate(self) -> Result<()> {
        finite(self.min)?;
        finite(self.max)?;
        if self.min > self.max {
            return Err(bad("coordinate weight minimum exceeds maximum"));
        }
        Ok(())
    }
    /// Transform a training value (source `weightData` for one datum): clamp
    /// it into the datum range (source `checkDatumRange`), then transform it.
    /// The value is clamped only for a nonidentity weight, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for crossed or non-finite bounds, or a
    /// non-finite value or result. The source does not check the bounds and
    /// lets `ln` of a non-positive value produce `-inf` or NaN.
    pub fn transform_training(self, value: f64) -> Result<f64> {
        self.validate()?;
        finite(value)?;
        self.transform(if self.function == WeightFunction::Identity {
            value
        } else {
            value.clamp(self.min, self.max)
        })
    }
    /// Transform a value without clamping (source `weightDatum`), as the
    /// linear model does when it is evaluated; evaluation intentionally does
    /// not clamp coordinates, following the source.
    ///
    /// # Errors
    ///
    /// As [`CoordinateWeight::transform_training`].
    pub fn transform(self, value: f64) -> Result<f64> {
        self.validate()?;
        finite(value)?;
        finite(match self.function {
            WeightFunction::Identity => value,
            WeightFunction::Log => value.ln(),
            WeightFunction::Reciprocal => 1.0 / value.abs(),
            WeightFunction::ReciprocalSquared => 1.0 / finite(value * value)?,
        })
    }
    /// Undo [`CoordinateWeight::transform`] (source `unWeightDatum`).
    /// Reciprocal transformations discard the sign, so their inverse yields
    /// the non-negative branch.
    ///
    /// # Errors
    ///
    /// As [`CoordinateWeight::transform_training`].
    pub fn untransform(self, value: f64) -> Result<f64> {
        self.validate()?;
        finite(value)?;
        finite(match self.function {
            WeightFunction::Identity => value,
            WeightFunction::Log => value.exp(),
            WeightFunction::Reciprocal => 1.0 / value.abs(),
            WeightFunction::ReciprocalSquared => (1.0 / value.abs()).sqrt(),
        })
    }
}
/// The coordinate weighting state of the source's base `TransformationModel`
/// (`x_weight_`, `x_datum_min_`, `x_datum_max_` and the y equivalents), which
/// its `(DataPoints, Param)` constructor reads; see
/// [`ModelWeights::from_param`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ModelWeights {
    /// Weighting of the x values.
    pub x: CoordinateWeight,
    /// Weighting of the y values.
    pub y: CoordinateWeight,
}
impl ModelWeights {
    /// Whether any axis is weighted (source `weighting_`).
    pub fn is_weighting(&self) -> bool {
        self.x.function != WeightFunction::Identity || self.y.function != WeightFunction::Identity
    }
    /// Weight data points in place (source `weightData`): for each axis whose
    /// weight is not the identity, clamp each coordinate into its datum range
    /// (`checkDatumRange`) and transform it (`weightDatum`).
    ///
    /// # Errors
    ///
    /// As [`CoordinateWeight::transform_training`]: a non-finite result, which
    /// the source stores (`ln(0)` is `-inf`), is refused, and `data` is left
    /// unchanged then.
    pub fn weight_data(&self, data: &mut [DataPoint]) -> Result<()> {
        let weighted = data
            .iter()
            .map(|p| {
                Ok((
                    self.x.transform_training(p.x)?,
                    self.y.transform_training(p.y)?,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        for (point, (x, y)) in data.iter_mut().zip(weighted) {
            point.x = x;
            point.y = y;
        }
        Ok(())
    }
    /// Undo [`ModelWeights::weight_data`] in place (source `unWeightData`),
    /// without clamping.
    ///
    /// # Errors
    ///
    /// As [`CoordinateWeight::untransform`]; `data` is unchanged on error.
    pub fn unweight_data(&self, data: &mut [DataPoint]) -> Result<()> {
        let unweighted = data
            .iter()
            .map(|p| Ok((self.x.untransform(p.x)?, self.y.untransform(p.y)?)))
            .collect::<Result<Vec<_>>>()?;
        for (point, (x, y)) in data.iter_mut().zip(unweighted) {
            point.x = x;
            point.y = y;
        }
        Ok(())
    }
}
/// Slope and intercept of a [`LinearModel`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearCoefficients {
    /// The slope.
    pub slope: f64,
    /// The intercept.
    pub intercept: f64,
}
/// Options of the linear model (source `TransformationModelLinear`
/// parameters; see [`LinearOptions::defaults`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearOptions {
    /// Weighting of the x values.
    pub x_weight: CoordinateWeight,
    /// Weighting of the y values.
    pub y_weight: CoordinateWeight,
    /// Explicit coefficients (source `slope` and `intercept`). Used only when
    /// no anchors are supplied, as in the source.
    pub coefficients: Option<LinearCoefficients>,
    /// Maximum number of anchors (native limit).
    pub max_points: usize,
}
impl Default for LinearOptions {
    fn default() -> Self {
        Self {
            x_weight: Default::default(),
            y_weight: Default::default(),
            coefficients: None,
            max_points: 1_000_000,
        }
    }
}
/// A linear model `f(x) = slope * x + intercept`, fitted by least squares in
/// the weighted coordinates (source `TransformationModelLinear`).
#[derive(Clone, Debug)]
pub struct LinearModel {
    coefficients: LinearCoefficients,
    options: LinearOptions,
}
impl LinearModel {
    /// Fit the model (the source constructor).
    ///
    /// With anchors, the coordinates are weighted and one anchor gives slope
    /// one and the observed offset, two anchors their exact line, and three or
    /// more the least-squares line from centred covariance and variance.
    /// Without anchors the explicit coefficients are used.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when there are neither anchors nor
    /// coefficients (the source's `Exception::IllegalArgument`, "no data
    /// points for 'linear' model"), when the x values do not vary (the
    /// source's `Exception::UnableToFit`), for invalid weight bounds, or for
    /// non-finite input or arithmetic.
    pub fn fit(data: &[DataPoint], options: LinearOptions) -> Result<Self> {
        validate_points(data, options.max_points)?;
        options.x_weight.validate()?;
        options.y_weight.validate()?;
        if let Some(c) = options.coefficients {
            finite(c.slope)?;
            finite(c.intercept)?;
        }
        let coefficients = if data.is_empty() {
            options
                .coefficients
                .ok_or_else(|| bad("linear fit requires anchors or explicit coefficients"))?
        } else {
            let weighted: Vec<_> = data
                .iter()
                .map(|p| {
                    Ok((
                        options.x_weight.transform_training(p.x)?,
                        options.y_weight.transform_training(p.y)?,
                    ))
                })
                .collect::<Result<_>>()?;
            let (slope, intercept) = if weighted.len() == 1 {
                (1.0, weighted[0].1 - weighted[0].0)
            } else if weighted.len() == 2 {
                let slope =
                    finite(weighted[1].1 - weighted[0].1)? / finite(weighted[1].0 - weighted[0].0)?;
                (slope, weighted[0].1 - slope * weighted[0].0)
            } else {
                // `gte::ApprHeightLine2::FitIndexed`: the mean is the running
                // sum times the reciprocal of the count, as GTE's
                // `Vector::operator/=` divides; `sum / count` can differ in
                // the last bit.
                let reciprocal = 1.0 / weighted.len() as f64;
                let mut sx = 0.0;
                let mut sy = 0.0;
                for &(x, y) in &weighted {
                    sx += x;
                    sy += y;
                }
                let mx = finite(sx * reciprocal)?;
                let my = finite(sy * reciprocal)?;
                let mut xx = 0.0;
                let mut xy = 0.0;
                for &(x, y) in &weighted {
                    let dx = x - mx;
                    xx += dx * dx;
                    xy += dx * (y - my);
                }
                finite(xx)?;
                finite(xy)?;
                if xx <= 0.0 {
                    return Err(bad("linear regression has no x variance"));
                }
                let slope = xy / xx;
                (slope, -slope * mx + my)
            };
            LinearCoefficients {
                slope: finite(slope)?,
                intercept: finite(intercept)?,
            }
        };
        Ok(Self {
            coefficients,
            options,
        })
    }
    /// A model with explicit coefficients and no weighting (the source
    /// constructor with `slope` and `intercept` parameters and no data).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a non-finite coefficient.
    pub fn from_coefficients(slope: f64, intercept: f64) -> Result<Self> {
        Self::fit(
            &[],
            LinearOptions {
                coefficients: Some(LinearCoefficients { slope, intercept }),
                ..Default::default()
            },
        )
    }
    /// The fitted slope and intercept (the first two outputs of the source's
    /// eight-output `getParameters`).
    pub fn coefficients(&self) -> LinearCoefficients {
        self.coefficients
    }
    /// The weights and bounds (the other six outputs of `getParameters`).
    pub fn options(&self) -> LinearOptions {
        self.options
    }
    /// Evaluate the model (source `evaluate`): weight x, apply the line, and
    /// undo the y weighting.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a non-finite input or result.
    pub fn apply(&self, x: f64) -> Result<f64> {
        let x = self.options.x_weight.transform(x)?;
        self.options.y_weight.untransform(finite(
            self.coefficients.slope * x + self.coefficients.intercept,
        )?)
    }
    /// The model's weighting.
    pub fn weights(&self) -> ModelWeights {
        ModelWeights {
            x: self.options.x_weight,
            y: self.options.y_weight,
        }
    }
    /// Weight data points in place with this model's weights (source
    /// `weightData`); see [`ModelWeights::weight_data`].
    ///
    /// # Errors
    ///
    /// As [`ModelWeights::weight_data`].
    pub fn weight_data(&self, data: &mut [DataPoint]) -> Result<()> {
        self.weights().weight_data(data)
    }
    /// Undo [`LinearModel::weight_data`] in place (source `unWeightData`).
    ///
    /// # Errors
    ///
    /// As [`ModelWeights::unweight_data`].
    pub fn unweight_data(&self, data: &mut [DataPoint]) -> Result<()> {
        self.weights().unweight_data(data)
    }
    /// Algebraic inverse of the fitted model, swapping its coordinate weights
    /// and bounds (source `TransformationModelLinear::invert`).
    /// Reciprocal weights discard signs and therefore invert only their branch.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a zero slope (the source's
    /// `Exception::DivisionByZero`) or a non-finite result.
    pub fn inverse(&self) -> Result<Self> {
        if self.coefficients.slope == 0.0 {
            return Err(bad("zero-slope linear model has no inverse"));
        }
        Self::fit(
            &[],
            LinearOptions {
                coefficients: Some(LinearCoefficients {
                    slope: finite(1.0 / self.coefficients.slope)?,
                    intercept: finite(-self.coefficients.intercept / self.coefficients.slope)?,
                }),
                x_weight: self.options.y_weight,
                y_weight: self.options.x_weight,
                ..self.options
            },
        )
    }
}
/// Interpolation between data points (source `interpolation_type`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Interpolation {
    /// `linear`: linearly between neighbouring points.
    Linear,
    /// `cspline`: a natural cubic spline. The default.
    #[default]
    CubicSpline,
    /// `akima`: an Akima spline, less affected by outliers.
    Akima,
}
/// Extrapolation beyond the data range (source `extrapolation_type`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Extrapolation {
    /// `two-point-linear`: one line through the first and last averaged
    /// anchor. The interpolated model's default.
    #[default]
    TwoPointLinear,
    /// `four-point-linear`: a line through the first two averaged anchors in
    /// front and one through the last two at the end, which may approximate
    /// non-linear data better.
    FourPointLinear,
    /// `global-linear`: a least-squares line through all anchors, including
    /// repeated x values. The source notes it may not be continuous with the
    /// interpolation at the border.
    GlobalLinear,
}
/// Options of the interpolated model (see [`InterpolationOptions::defaults`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InterpolationOptions {
    /// Interpolation between the anchors.
    pub interpolation: Interpolation,
    /// Extrapolation beyond them.
    pub extrapolation: Extrapolation,
    /// Maximum number of anchors (native limit).
    pub max_points: usize,
}
impl Default for InterpolationOptions {
    fn default() -> Self {
        Self {
            interpolation: Default::default(),
            extrapolation: Default::default(),
            max_points: 1_000_000,
        }
    }
}
/// An interpolation model with linear extrapolation (source
/// `TransformationModelInterpolated`).
///
/// The source's public nested `Interpolator` class, an abstract interface its
/// `.cpp` implements three times, is the closed [`Interpolation`] choice here.
#[derive(Clone, Debug)]
pub struct InterpolatedModel {
    x: Vec<f64>,
    y: Vec<f64>,
    spline: Option<CubicSpline2d>,
    akima: Option<akima::Akima>,
    front: LinearModel,
    back: LinearModel,
    options: InterpolationOptions,
}
impl InterpolatedModel {
    /// Fit the model (the source constructor from `DataPoints`).
    ///
    /// The anchors are sorted by x, and y is averaged over repeated x, in the
    /// anchors' order, as the source's `std::map` preprocessing does. At least
    /// three distinct x values are required, even for linear interpolation.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for fewer than three distinct x values
    /// (the source's `Exception::IllegalArgument`), a failed spline or
    /// extrapolation fit, or non-finite input or arithmetic.
    pub fn fit(data: &[DataPoint], options: InterpolationOptions) -> Result<Self> {
        validate_points(data, options.max_points)?;
        let (x, y) = averaged(data.iter().map(|p| (p.x, p.y)).collect())?;
        let global = data.to_vec();
        Self::build(x, y, options, || Ok(global))
    }

    /// Fit the model from plain pairs (the source constructor
    /// `TransformationModelInterpolated(const std::vector<std::pair<double,
    /// double>>&, const Param&, bool preprocess)`), reproducing that
    /// overload's behaviour, which differs from [`InterpolatedModel::fit`] in
    /// two ways the source does not document:
    ///
    /// - With `preprocess`, anchors are sorted and averaged as in `fit`, but
    ///   `global-linear` extrapolation fits its line to as many `(0, 0)` points
    ///   as there are distinct x values, followed by the averaged anchors: the
    ///   source resizes the vector and then appends to it.
    /// - Without `preprocess`, the interpolation runs on as many `(0, 0)`
    ///   points as there are pairs, followed by the pairs as given, for the
    ///   same reason; linear interpolation then searches that sequence with
    ///   `std::upper_bound` whether or not it is sorted, as the source does.
    ///
    /// Both are recorded as defect candidates in the support document.
    ///
    /// # Errors
    ///
    /// As [`InterpolatedModel::fit`]. Without `preprocess`, a sequence that is
    /// not strictly increasing is refused by the cubic and Akima splines,
    /// where the source's spline divides by a zero interval (cubic) or throws
    /// `std::runtime_error` (Akima), and fewer than three points are refused.
    pub fn fit_pairs(
        pairs: &[(f64, f64)],
        options: InterpolationOptions,
        preprocess: bool,
    ) -> Result<Self> {
        if options.max_points == 0 || pairs.len() > options.max_points / 2 {
            return Err(bad("transformation point limit exceeded or zero"));
        }
        for &(x, y) in pairs {
            finite(x)?;
            finite(y)?;
        }
        if preprocess {
            let (x, y) = averaged(pairs.to_vec())?;
            let padded: Vec<DataPoint> = std::iter::repeat_n(DataPoint::default(), x.len())
                .chain(x.iter().zip(&y).map(|(&x, &y)| DataPoint::new(x, y)))
                .collect();
            Self::build(x, y, options, || Ok(padded))
        } else {
            let n = pairs.len();
            let x: Vec<f64> = std::iter::repeat_n(0.0, n)
                .chain(pairs.iter().map(|p| p.0))
                .collect();
            let y: Vec<f64> = std::iter::repeat_n(0.0, n)
                .chain(pairs.iter().map(|p| p.1))
                .collect();
            if x.len() < 3 {
                return Err(bad("interpolation requires at least three points"));
            }
            let global: Vec<DataPoint> = x
                .iter()
                .zip(&y)
                .map(|(&x, &y)| DataPoint::new(x, y))
                .collect();
            Self::build(x, y, options, || Ok(global))
        }
    }

    fn build(
        x: Vec<f64>,
        y: Vec<f64>,
        options: InterpolationOptions,
        global: impl FnOnce() -> Result<Vec<DataPoint>>,
    ) -> Result<Self> {
        for pair in x.windows(2) {
            finite(pair[1] - pair[0])?;
        }
        let spline = match options.interpolation {
            Interpolation::CubicSpline => {
                Some(CubicSpline2d::with_max_points(&x, &y, options.max_points)?)
            }
            _ => None,
        };
        let akima = match options.interpolation {
            Interpolation::Akima => Some(akima::Akima::new(&x, &y)?),
            _ => None,
        };
        let line = |a: usize, b: usize| {
            LinearModel::fit(
                &[DataPoint::new(x[a], y[a]), DataPoint::new(x[b], y[b])],
                LinearOptions::default(),
            )
        };
        let (front, back) = match options.extrapolation {
            Extrapolation::TwoPointLinear => {
                let m = line(0, x.len() - 1)?;
                (m.clone(), m)
            }
            Extrapolation::FourPointLinear => (line(0, 1)?, line(x.len() - 2, x.len() - 1)?),
            Extrapolation::GlobalLinear => {
                let m = LinearModel::fit(
                    &global()?,
                    LinearOptions {
                        max_points: options.max_points,
                        ..Default::default()
                    },
                )?;
                (m.clone(), m)
            }
        };
        Ok(Self {
            x,
            y,
            spline,
            akima,
            front,
            back,
            options,
        })
    }
    /// The options the model was fitted with.
    pub fn options(&self) -> InterpolationOptions {
        self.options
    }
    /// The points interpolated between, sorted by x with repeated x averaged.
    pub fn knots(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        self.x.iter().copied().zip(self.y.iter().copied())
    }
    /// Evaluate the model (source `evaluate`): extrapolate below the first
    /// and above the last knot, interpolate in between.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a non-finite input or result.
    pub fn apply(&self, value: f64) -> Result<f64> {
        finite(value)?;
        if value < self.x[0] {
            return self.front.apply(value);
        }
        if value > self.x[self.x.len() - 1] {
            return self.back.apply(value);
        }
        if let Some(s) = &self.spline {
            return s.eval(value);
        }
        if let Some(a) = &self.akima {
            return finite(a.eval(value));
        }
        // LinearInterpolator::eval: std::upper_bound, then the last y at the
        // end, and a chord otherwise.
        let right = libstdcxx::upper_bound(&self.x, |&x| value < x);
        if right == self.x.len() {
            return Ok(self.y[self.y.len() - 1]);
        }
        // `value >= x[0]` here, so the halving search never returns the first
        // position; the source would read before the range if it did.
        let Some(left) = right.checked_sub(1) else {
            return Err(bad("linear interpolation found no interval"));
        };
        finite(
            self.y[left]
                + (self.y[right] - self.y[left]) * (value - self.x[left])
                    / (self.x[right] - self.x[left]),
        )
    }
}

/// Sort pairs by x and average y over repeated x (source
/// `preprocessDataPoints_`), requiring three distinct x values.
fn averaged(mut sorted: Vec<(f64, f64)>) -> Result<(Vec<f64>, Vec<f64>)> {
    sorted.sort_by(|a, b| a.0.total_cmp(&b.0));
    let (mut x, mut y) = (Vec::new(), Vec::new());
    let mut i = 0;
    while i < sorted.len() {
        let mut j = i;
        // std::accumulate(begin, end, 0.0)
        let mut sum = 0.0;
        while j < sorted.len() && sorted[j].0 == sorted[i].0 {
            sum += sorted[j].1;
            j += 1;
        }
        x.push(sorted[i].0);
        y.push(finite(sum / (j - i) as f64)?);
        i = j;
    }
    if x.len() < 3 {
        return Err(bad(
            "Cubic spline model needs at least 3 data points (with unique x values)",
        ));
    }
    Ok((x, y))
}

/// Which model to fit, with its options.
#[derive(Clone, Debug)]
pub enum ModelConfig {
    /// `none`: no model, the identity.
    None,
    /// `identity`: the identity, marked as optimal (for reference files); a
    /// description holding it ignores later fits.
    Identity,
    /// `linear`.
    Linear(LinearOptions),
    /// `b_spline`.
    BSpline(BSplineOptions),
    /// `interpolated`.
    Interpolated(InterpolationOptions),
    /// `lowess` with a fixed span.
    Lowess(LowessOptions),
    /// `lowess` with the span selected by cross-validation (`auto_span`); the
    /// options' span is replaced by the selected one when the model is
    /// fitted.
    LowessAuto(LowessOptions, LowessAutoSpan),
}
/// A fitted transformation model.
#[derive(Clone, Debug)]
pub enum TransformationModel {
    /// No model: the identity (the source's base `TransformationModel`).
    None,
    /// The identity, fitted on purpose.
    Identity,
    /// A linear model.
    Linear(LinearModel),
    /// A B-spline model.
    BSpline(BSplineModel),
    /// An interpolated model.
    Interpolated(InterpolatedModel),
    /// A LOWESS model.
    Lowess(LowessModel),
}
impl TransformationModel {
    /// Evaluate the model at `value` (source `evaluate`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a non-finite input or result.
    pub fn apply(&self, value: f64) -> Result<f64> {
        match self {
            Self::None | Self::Identity => finite(value),
            Self::Linear(m) => m.apply(value),
            Self::BSpline(m) => m.apply(value),
            Self::Interpolated(m) => m.apply(value),
            Self::Lowess(m) => m.apply(value),
        }
    }
    /// The model type (source `TransformationDescription::getModelType`).
    pub fn name(&self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Identity => "identity",
            Self::Linear(_) => "linear",
            Self::BSpline(_) => "b_spline",
            Self::Interpolated(_) => "interpolated",
            Self::Lowess(_) => "lowess",
        }
    }
    /// Weight data points with the model's weights (source `weightData`).
    ///
    /// Only the linear model carries weights. The other source models are
    /// built with the base class's default constructor, which leaves the
    /// weight members uninitialized, so `weightData` on them reads an
    /// indeterminate flag; they evaluate without weighting, and this leaves
    /// the data unchanged for them.
    ///
    /// # Errors
    ///
    /// As [`LinearModel::weight_data`].
    pub fn weight_data(&self, data: &mut [DataPoint]) -> Result<()> {
        match self {
            Self::Linear(m) => m.weight_data(data),
            _ => Ok(()),
        }
    }
    /// Undo [`TransformationModel::weight_data`] (source `unWeightData`).
    ///
    /// # Errors
    ///
    /// As [`LinearModel::unweight_data`].
    pub fn unweight_data(&self, data: &mut [DataPoint]) -> Result<()> {
        match self {
            Self::Linear(m) => m.unweight_data(data),
            _ => Ok(()),
        }
    }
    /// Fit a configuration, returning the model and the configuration with an
    /// automatic LOWESS span resolved to the selected one.
    fn fit(data: &[DataPoint], config: &ModelConfig) -> Result<(Self, ModelConfig)> {
        Ok(match config {
            ModelConfig::None => (Self::None, ModelConfig::None),
            ModelConfig::Identity => (Self::Identity, ModelConfig::Identity),
            ModelConfig::Linear(o) => (Self::Linear(LinearModel::fit(data, *o)?), config.clone()),
            ModelConfig::BSpline(o) => {
                (Self::BSpline(BSplineModel::fit(data, *o)?), config.clone())
            }
            ModelConfig::Interpolated(o) => (
                Self::Interpolated(InterpolatedModel::fit(data, *o)?),
                config.clone(),
            ),
            ModelConfig::Lowess(o) => (Self::Lowess(LowessModel::fit(data, *o)?), config.clone()),
            ModelConfig::LowessAuto(o, auto) => {
                let (model, _) = LowessModel::fit_auto(data, *o, auto)?;
                let resolved = ModelConfig::Lowess(model.options());
                (Self::Lowess(model), resolved)
            }
        })
    }
    /// The parameters this model reports for a typed configuration: the
    /// configuration's [`ModelConfig::to_param`], and for a fitted linear
    /// model its slope and intercept.
    fn synthesized_parameters(&self, config: &ModelConfig, data_given: bool) -> Result<Param> {
        let mut params = config.to_param()?;
        if let Self::Linear(model) = self {
            if data_given {
                let c = model.coefficients();
                params.set_value("slope", ParamValue::Float(c.slope), "", &[])?;
                params.set_value("intercept", ParamValue::Float(c.intercept), "", &[])?;
            }
        }
        Ok(params)
    }
}
/// Anchors and the transformation model fitted to them (source
/// `TransformationDescription`).
///
/// The description stores data points; the model stores parameters. The
/// same data can therefore be modelled with different models, and a model
/// keeps whatever representation of the data it needs.
///
/// A clone copies the fitted model. The source's copy constructor and
/// assignment refit the model from its reported parameters instead, which
/// gives the same model where it succeeds; it fails in the source for an
/// inverted linear model without data, whose reported weight names (`y` as
/// the x weight) its own constructor rejects.
#[derive(Clone, Debug)]
pub struct TransformationDescription {
    data: Vec<DataPoint>,
    model: TransformationModel,
    config: ModelConfig,
    parameters: Param,
    max_points: usize,
}
impl Default for TransformationDescription {
    fn default() -> Self {
        Self {
            data: Vec::new(),
            model: TransformationModel::None,
            config: ModelConfig::None,
            parameters: Param::new(),
            max_points: 1_000_000,
        }
    }
}
impl TransformationDescription {
    /// A description of these anchors with no model (model type `none`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a non-finite coordinate or more
    /// than a million anchors.
    pub fn new(data: Vec<DataPoint>) -> Result<Self> {
        Self::with_max_points(data, 1_000_000)
    }
    /// As [`TransformationDescription::new`], with a native anchor limit that
    /// also bounds batch application.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a non-finite coordinate, more than
    /// `max_points` anchors, or a zero limit.
    pub fn with_max_points(data: Vec<DataPoint>, max_points: usize) -> Result<Self> {
        validate_points(&data, max_points)?;
        Ok(Self {
            data,
            max_points,
            ..Default::default()
        })
    }
    /// The anchors (source `getDataPoints`).
    pub fn data_points(&self) -> &[DataPoint] {
        &self.data
    }
    /// The fitted model.
    pub fn model(&self) -> &TransformationModel {
        &self.model
    }
    /// The type of the fitted model (source `getModelType`): `none`,
    /// `identity`, `linear`, `b_spline`, `interpolated` or `lowess`.
    pub fn model_type(&self) -> &'static str {
        self.model.name()
    }
    /// The model types a description can fit (source static `getModelTypes`):
    /// `linear`, `b_spline`, `interpolated`, `lowess`. `none` and `identity`
    /// do not count.
    pub fn model_types() -> [&'static str; 4] {
        MODEL_TYPES
    }
    /// The configuration the model was fitted with. An automatically selected
    /// LOWESS span is resolved to the selected value.
    pub fn model_config(&self) -> &ModelConfig {
        &self.config
    }
    /// The model's parameters (source `getModelParameters`).
    ///
    /// After [`TransformationDescription::fit_model_with_param`] these are
    /// what the source model keeps: the given parameters with the model's
    /// defaults filled in (linear: plus the fitted `slope` and `intercept`,
    /// unless explicit coefficients were used without data; LOWESS with
    /// `auto_span`: plus the selected `span` and `auto_span` reset to
    /// `false`), and nothing for `none` and `identity`. After
    /// [`TransformationDescription::fit_model`] they are the typed
    /// configuration in the same layout.
    pub fn model_parameters(&self) -> &Param {
        &self.parameters
    }
    /// Replace the anchors (source `setDataPoints`). This resets the model to
    /// `none`, including an `identity` model.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a non-finite coordinate or too
    /// many anchors; the description is unchanged then.
    pub fn set_data_points(&mut self, data: Vec<DataPoint>) -> Result<()> {
        validate_points(&data, self.max_points)?;
        self.data = data;
        self.model = TransformationModel::None;
        self.config = ModelConfig::None;
        self.parameters = Param::new();
        Ok(())
    }
    /// Replace the anchors with plain pairs (the source's backwards-compatible
    /// `setDataPoints(const std::vector<std::pair<double, double>>&)`).
    ///
    /// # Errors
    ///
    /// As [`TransformationDescription::set_data_points`].
    pub fn set_data_point_pairs(&mut self, pairs: &[(f64, f64)]) -> Result<()> {
        self.set_data_points(pairs.iter().copied().map(DataPoint::from).collect())
    }
    /// Fit a model to the anchors (source `fitModel`).
    ///
    /// An `identity` model is never replaced: fit requests are ignored until
    /// the anchors are reset, as in the source.
    ///
    /// # Errors
    ///
    /// Returns the model's fitting error; the previous model is kept then.
    pub fn fit_model(&mut self, config: ModelConfig) -> Result<()> {
        if matches!(self.model, TransformationModel::Identity) {
            return Ok(());
        }
        let (model, resolved) = TransformationModel::fit(&self.data, &config)?;
        let parameters = model.synthesized_parameters(&resolved, !self.data.is_empty())?;
        self.model = model;
        self.config = resolved;
        self.parameters = parameters;
        Ok(())
    }
    /// Fit a model named by `model_type` from a `Param` (source
    /// `fitModel(model_type, params)`).
    ///
    /// The model reads `params` as its source constructor does (see
    /// [`ModelConfig::from_param`]), and
    /// [`TransformationDescription::model_parameters`] then returns what the
    /// source model keeps. An `identity` model is never replaced.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an unknown model type (the source's
    /// `Exception::IllegalArgument`), a parameter error, or the model's
    /// fitting error; the previous model is kept then.
    pub fn fit_model_with_param(&mut self, model_type: &str, params: &Param) -> Result<()> {
        if matches!(self.model, TransformationModel::Identity) {
            return Ok(());
        }
        let config = ModelConfig::from_param(model_type, params)?;
        let (model, resolved) = TransformationModel::fit(&self.data, &config)?;
        let parameters =
            params::stored_parameters(model_type, params, !self.data.is_empty(), &model)?;
        self.model = model;
        self.config = resolved;
        self.parameters = parameters;
        Ok(())
    }
    /// Apply the transformation (source `apply`); the value is returned
    /// unchanged when no model was fitted.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a non-finite input or result.
    pub fn apply(&self, value: f64) -> Result<f64> {
        self.model.apply(value)
    }
    /// Apply the transformation to every value atomically; a failure leaves
    /// every value unchanged.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for more values than the anchor limit
    /// or the first failed application.
    pub fn apply_values(&self, values: &mut [f64]) -> Result<()> {
        if values.len() > self.max_points {
            return Err(bad("transformation batch exceeds point limit"));
        }
        let transformed: Vec<_> = values
            .iter()
            .map(|&v| self.apply(v))
            .collect::<Result<_>>()?;
        values.copy_from_slice(&transformed);
        Ok(())
    }
    /// An approximate inverse (source `invert`, without mutating this
    /// description): the anchors are swapped, keeping their notes, and the
    /// same model is refitted to them. A linear model without anchors is
    /// inverted algebraically instead, as the source's special case does.
    ///
    /// The reported [`model_parameters`](Self::model_parameters) follow the
    /// source: a refitted linear model updates `slope` and `intercept`, and
    /// the algebraic inverse also swaps the weights and datum bounds.
    ///
    /// # Errors
    ///
    /// Returns the refit's error, or [`Error::InvalidValue`] for a zero slope
    /// in the algebraic case. The source swaps the anchors before refitting
    /// and so leaves them swapped when the refit throws; see
    /// [`TransformationDescription::invert`].
    pub fn inverse(&self) -> Result<Self> {
        let data = self
            .data
            .iter()
            .map(|p| DataPoint::with_note(p.y, p.x, p.note.clone()))
            .collect();
        let mut result = Self::with_max_points(data, self.max_points)?;
        if self.data.is_empty() {
            if let TransformationModel::Linear(model) = &self.model {
                let inverse = model.inverse()?;
                result.parameters = inverted_linear_parameters(&self.parameters, model, &inverse)?;
                result.config = ModelConfig::Linear(inverse.options());
                result.model = TransformationModel::Linear(inverse);
                return Ok(result);
            }
        }
        if matches!(self.model, TransformationModel::Identity) {
            // The source's fitModel returns early for an identity model, so
            // the swapped anchors keep it.
            result.model = TransformationModel::Identity;
            result.config = ModelConfig::Identity;
            result.parameters = self.parameters.checked_clone()?;
            return Ok(result);
        }
        let (model, resolved) = TransformationModel::fit(&result.data, &self.config)?;
        let mut parameters = self.parameters.checked_clone()?;
        if let TransformationModel::Linear(linear) = &model {
            let c = linear.coefficients();
            parameters.set_value("slope", ParamValue::Float(c.slope), "", &[])?;
            parameters.set_value("intercept", ParamValue::Float(c.intercept), "", &[])?;
        }
        result.model = model;
        result.config = resolved;
        result.parameters = parameters;
        Ok(result)
    }
    /// Replace this description by its [`inverse`](Self::inverse), only if
    /// that succeeds (source `invert`).
    ///
    /// # Errors
    ///
    /// As [`TransformationDescription::inverse`]; the description is
    /// unchanged on error, where the source leaves its anchors swapped.
    pub fn invert(&mut self) -> Result<()> {
        *self = self.inverse()?;
        Ok(())
    }
    /// Absolute deviations `|x - y|` between the anchor pairs, or `|f(x) - y|`
    /// with `apply`, optionally sorted (source `getDeviations`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a failed application or a
    /// non-finite deviation.
    pub fn deviations(&self, apply: bool, sort: bool) -> Result<Vec<f64>> {
        let mut result: Vec<_> = self
            .data
            .iter()
            .map(|p| finite((if apply { self.apply(p.x)? } else { p.x } - p.y).abs()))
            .collect::<Result<_>>()?;
        if sort {
            result.sort_by(f64::total_cmp);
        }
        Ok(result)
    }
    /// Data ranges and deviation percentiles before and after the
    /// transformation (source `getStatistics`).
    ///
    /// The percentiles are 100, 99, 95, 90, 75, 50 and 25, each taken at the
    /// sorted deviation of rank `p / 100 * n - 1`, the source's discrete
    /// formula. For small `n` that rank is negative, which the source converts
    /// to `Size` (undefined behaviour); it is clamped to the first deviation
    /// here. Without anchors the result is empty.
    ///
    /// # Errors
    ///
    /// As [`TransformationDescription::deviations`].
    pub fn statistics(&self) -> Result<TransformationStatistics> {
        let mut stats = TransformationStatistics::default();
        if self.data.is_empty() {
            return Ok(stats);
        }
        let range = |values: Vec<f64>| NumericRange {
            min: values.iter().copied().fold(f64::INFINITY, f64::min),
            max: values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        };
        stats.x_range = Some(range(self.data.iter().map(|p| p.x).collect()));
        stats.y_range = Some(range(self.data.iter().map(|p| p.y).collect()));
        let before = self.deviations(false, true)?;
        let after = self.deviations(true, true)?;
        for percent in [100, 99, 95, 90, 75, 50, 25] {
            let index =
                (f64::from(percent) / 100.0 * self.data.len() as f64 - 1.0).max(0.0) as usize;
            stats.percentiles.push(DeviationPercentile {
                percent,
                before: before[index],
                after: after[index],
            });
        }
        Ok(stats)
    }
    /// Write summary statistics (source `printSummary`), formatted as the
    /// source's `std::ostream` does by default (six significant digits).
    ///
    /// The deviations after the transformation are printed only when a model
    /// other than `none` or `identity` is fitted.
    ///
    /// # Errors
    ///
    /// As [`TransformationDescription::statistics`], and
    /// [`Error::Io`] for a failed write.
    pub fn print_summary(&self, out: &mut impl Write) -> Result<()> {
        let stats = self.statistics()?;
        writeln!(
            out,
            "Number of data points (x/y pairs): {}",
            self.data.len()
        )?;
        let (Some(x), Some(y)) = (stats.x_range, stats.y_range) else {
            return Ok(());
        };
        writeln!(
            out,
            "Data range (x): {} to {}\nData range (y): {} to {}",
            format_g(x.min),
            format_g(x.max),
            format_g(y.min),
            format_g(y.max)
        )?;
        let no_model = matches!(
            self.model,
            TransformationModel::None | TransformationModel::Identity
        );
        writeln!(
            out,
            "Summary of x/y deviations{}:",
            if no_model {
                ""
            } else {
                " before transformation"
            }
        )?;
        for p in &stats.percentiles {
            writeln!(
                out,
                "- {:>3}% of data points within (+/-){}",
                p.percent,
                format_g(p.before)
            )?;
        }
        if no_model {
            writeln!(out)?;
            return Ok(());
        }
        writeln!(
            out,
            "Summary of x/y deviations after applying '{}' transformation:",
            self.model_type()
        )?;
        for p in &stats.percentiles {
            writeln!(
                out,
                "- {:>3}% of data points within (+/-){}",
                p.percent,
                format_g(p.after)
            )?;
        }
        writeln!(out)?;
        Ok(())
    }
    /// A residual-based extraction window (source `estimateWindow`).
    ///
    /// The absolute residuals are `|f(x_i) - y_i|` in transformed units, or,
    /// with [`WindowOptions::inverse`], `|x_i - f^-1(y_i)|` in original units
    /// after refitting the inverse as [`TransformationDescription::inverse`]
    /// does. The window is an adaptive quantile that is robust to sparse
    /// outliers but permissive when tails are dense:
    ///
    /// 1. the Tukey upper fence `UF = Q3 + 1.5 * IQR`;
    /// 2. `h_raw` is the interpolated `quantile` of the residuals, `h_rob` that
    ///    of the residuals capped at `UF`;
    /// 3. with `tail` the fraction above `UF`, `h = (1 - w) * h_rob + w * h_raw`
    ///    where `w` rises linearly from 0 at a tail of 1 % to 1 at 10 %.
    ///
    /// With fewer than four residuals or a zero IQR, the raw quantile is used.
    /// The result is `2h` for a full window, otherwise `h`, times the padding
    /// factor; zero without anchors.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a quantile outside `(0, 1]`, a
    /// negative or non-finite padding factor, a failed inverse, or a
    /// non-finite residual (the source silently drops those).
    pub fn estimate_window(&self, options: WindowOptions) -> Result<f64> {
        if !options.quantile.is_finite()
            || options.quantile <= 0.0
            || options.quantile > 1.0
            || !options.padding_factor.is_finite()
            || options.padding_factor < 0.0
        {
            return Err(bad("invalid residual-window quantile or padding"));
        }
        let inverse;
        let source = if options.inverse {
            inverse = self.inverse()?;
            &inverse
        } else {
            self
        };
        let residuals = source.deviations(true, true)?;
        if residuals.is_empty() {
            return Ok(0.0);
        }
        let raw = quantile(&residuals, options.quantile);
        let mut robust = raw;
        let mut weight = 0.0;
        if residuals.len() >= 4 {
            let q1 = quantile(&residuals, 0.25);
            let q3 = quantile(&residuals, 0.75);
            let iqr = q3 - q1;
            if iqr > 0.0 {
                let fence = finite(q3 + 1.5 * iqr)?;
                let tail = residuals.iter().filter(|&&v| v > fence).count() as f64
                    / residuals.len() as f64;
                let winsorized: Vec<_> = residuals.iter().map(|&v| v.min(fence)).collect();
                robust = quantile(&winsorized, options.quantile);
                weight = ((tail - 0.01) / 0.09).clamp(0.0, 1.0);
            }
        }
        finite(
            ((1.0 - weight) * robust + weight * raw)
                * if options.full_window { 2.0 } else { 1.0 }
                * options.padding_factor,
        )
    }
}

/// The source's `TransformationModelLinear::invert` parameter update: the
/// new slope and intercept, and the x and y weight names and datum bounds
/// swapped, as the member values (with their defaults) were.
fn inverted_linear_parameters(
    parameters: &Param,
    forward: &LinearModel,
    inverse: &LinearModel,
) -> Result<Param> {
    let mut updated = parameters.checked_clone()?;
    let name = |key: &str, fallback: &str| -> Result<String> {
        if parameters.exists(key)? {
            let value = parameters.value(key)?.as_str()?;
            if !value.is_empty() {
                return Ok(value.to_owned());
            }
        }
        Ok(fallback.to_owned())
    };
    let forward_options = forward.options();
    let c = inverse.coefficients();
    let x_weight = name("y_weight", forward_options.y_weight.function.y_name())?;
    let y_weight = name("x_weight", forward_options.x_weight.function.x_name())?;
    updated.set_value("slope", ParamValue::Float(c.slope), "", &[])?;
    updated.set_value("intercept", ParamValue::Float(c.intercept), "", &[])?;
    updated.set_value("x_weight", ParamValue::String(x_weight), "", &[])?;
    updated.set_value("y_weight", ParamValue::String(y_weight), "", &[])?;
    for (key, value) in [
        ("x_datum_min", forward_options.y_weight.min),
        ("x_datum_max", forward_options.y_weight.max),
        ("y_datum_min", forward_options.x_weight.min),
        ("y_datum_max", forward_options.x_weight.max),
    ] {
        updated.set_value(key, ParamValue::Float(value), "", &[])?;
    }
    Ok(updated)
}

/// Summary statistics of a description (source
/// `TransformationDescription::TransformationStatistics`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TransformationStatistics {
    /// Smallest and largest x (source `xmin`, `xmax`); `None` without anchors,
    /// where the source reports zeros.
    pub x_range: Option<NumericRange>,
    /// Smallest and largest y (source `ymin`, `ymax`).
    pub y_range: Option<NumericRange>,
    /// Deviation percentiles in the source's `percents` order (100, 99, 95,
    /// 90, 75, 50, 25); the source's `percentiles_before` and
    /// `percentiles_after` maps.
    pub percentiles: Vec<DeviationPercentile>,
}
/// One deviation percentile before and after the transformation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DeviationPercentile {
    /// The percentile.
    pub percent: u32,
    /// Deviation within which this percentage of anchors lies before.
    pub before: f64,
    /// And after applying the model; equal to `before` without a model.
    pub after: f64,
}
/// Options of [`TransformationDescription::estimate_window`], with the source
/// defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindowOptions {
    /// Quantile of the absolute residuals, in `(0, 1]` (default 0.99).
    pub quantile: f64,
    /// Measure the residuals in original x units through the inverse (default
    /// `true`); in transformed y units otherwise.
    pub inverse: bool,
    /// Return the full width, twice the half-width (default `true`).
    pub full_window: bool,
    /// Factor applied to the result (default 1).
    pub padding_factor: f64,
}
impl Default for WindowOptions {
    fn default() -> Self {
        Self {
            quantile: 0.99,
            inverse: true,
            full_window: true,
            padding_factor: 1.0,
        }
    }
}
fn quantile(sorted: &[f64], q: f64) -> f64 {
    let position = q * (sorted.len() - 1) as f64;
    let i = position.floor() as usize;
    let fraction = position - i as f64;
    if fraction == 0.0 {
        sorted[i]
    } else {
        (1.0 - fraction) * sorted[i] + fraction * sorted[i + 1]
    }
}

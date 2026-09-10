// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Coordinate transformations fitted to annotated anchors.
//! Inversion of a description refits swapped anchors and is generally approximate.

mod lowess;
pub use lowess::{LowessModel, LowessOptions, lowess};

use crate::kernel::NumericRange;
use crate::processing::peak_picking::CubicSpline2d;
use crate::{Error, Result};

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

/// Coordinate correspondence, with optional source annotation.
#[derive(Clone, Debug, Default, PartialEq, PartialOrd)]
pub struct DataPoint {
    pub x: f64,
    pub y: f64,
    pub note: String,
}
impl DataPoint {
    pub fn new(x: f64, y: f64) -> Self {
        Self {
            x,
            y,
            note: String::new(),
        }
    }
    pub fn with_note(x: f64, y: f64, note: impl Into<String>) -> Self {
        Self {
            x,
            y,
            note: note.into(),
        }
    }
}
impl From<(f64, f64)> for DataPoint {
    fn from((x, y): (f64, f64)) -> Self {
        Self::new(x, y)
    }
}

/// The source calls these coordinate transformations "weights"; they are not
/// per-observation statistical regression weights.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum WeightFunction {
    #[default]
    Identity,
    Log,
    Reciprocal,
    ReciprocalSquared,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoordinateWeight {
    pub function: WeightFunction,
    pub min: f64,
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
    /// Training values are clamped only for a nonidentity coordinate weight.
    pub fn transform_training(self, value: f64) -> Result<f64> {
        self.validate()?;
        finite(value)?;
        self.transform(if self.function == WeightFunction::Identity {
            value
        } else {
            value.clamp(self.min, self.max)
        })
    }
    /// Evaluation intentionally does not clamp coordinates, following OpenMS.
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
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearCoefficients {
    pub slope: f64,
    pub intercept: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearOptions {
    pub x_weight: CoordinateWeight,
    pub y_weight: CoordinateWeight,
    /// Used only when no anchors are supplied, as in the source.
    pub coefficients: Option<LinearCoefficients>,
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
#[derive(Clone, Debug)]
pub struct LinearModel {
    coefficients: LinearCoefficients,
    options: LinearOptions,
}
impl LinearModel {
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
                let count = weighted.len() as f64;
                let mx = finite(weighted.iter().map(|p| p.0).sum::<f64>() / count)?;
                let my = finite(weighted.iter().map(|p| p.1).sum::<f64>() / count)?;
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
    pub fn from_coefficients(slope: f64, intercept: f64) -> Result<Self> {
        Self::fit(
            &[],
            LinearOptions {
                coefficients: Some(LinearCoefficients { slope, intercept }),
                ..Default::default()
            },
        )
    }
    pub fn coefficients(&self) -> LinearCoefficients {
        self.coefficients
    }
    pub fn options(&self) -> LinearOptions {
        self.options
    }
    pub fn apply(&self, x: f64) -> Result<f64> {
        let x = self.options.x_weight.transform(x)?;
        self.options.y_weight.untransform(finite(
            self.coefficients.slope * x + self.coefficients.intercept,
        )?)
    }
    /// Algebraic inverse of the fitted model, swapping its coordinate weights.
    /// Reciprocal weights discard signs and therefore invert only their branch.
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
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Interpolation {
    Linear,
    #[default]
    CubicSpline,
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Extrapolation {
    /// One chord through the first and last averaged anchor.
    #[default]
    TwoPointLinear,
    /// Separate lines through the first two and last two averaged anchors.
    FourPointLinear,
    /// Least squares over the original anchors, including repeated x values.
    GlobalLinear,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InterpolationOptions {
    pub interpolation: Interpolation,
    pub extrapolation: Extrapolation,
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
#[derive(Clone, Debug)]
pub struct InterpolatedModel {
    x: Vec<f64>,
    y: Vec<f64>,
    spline: Option<CubicSpline2d>,
    front: LinearModel,
    back: LinearModel,
    options: InterpolationOptions,
}
impl InterpolatedModel {
    /// Sort anchors by x and average y at repeated x. At least three distinct
    /// x values are required even for linear interpolation, following OpenMS.
    pub fn fit(data: &[DataPoint], options: InterpolationOptions) -> Result<Self> {
        validate_points(data, options.max_points)?;
        let mut sorted: Vec<_> = data.iter().map(|p| (p.x, p.y)).collect();
        sorted.sort_by(|a, b| a.0.total_cmp(&b.0));
        let (mut x, mut y) = (Vec::new(), Vec::new());
        let mut i = 0;
        while i < sorted.len() {
            let mut j = i + 1;
            let mut sum = sorted[i].1;
            while j < sorted.len() && sorted[j].0 == sorted[i].0 {
                sum += sorted[j].1;
                j += 1;
            }
            x.push(sorted[i].0);
            y.push(finite(sum / (j - i) as f64)?);
            i = j;
        }
        if x.len() < 3 {
            return Err(bad("interpolation requires three distinct x coordinates"));
        }
        for pair in x.windows(2) {
            finite(pair[1] - pair[0])?;
        }
        let spline = match options.interpolation {
            Interpolation::Linear => None,
            Interpolation::CubicSpline => {
                Some(CubicSpline2d::with_max_points(&x, &y, options.max_points)?)
            }
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
                    data,
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
            front,
            back,
            options,
        })
    }
    pub fn options(&self) -> InterpolationOptions {
        self.options
    }
    pub fn knots(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        self.x.iter().copied().zip(self.y.iter().copied())
    }
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
        let right = self.x.partition_point(|&x| x <= value);
        if right == self.x.len() {
            return Ok(self.y[self.y.len() - 1]);
        }
        let left = right - 1;
        finite(
            self.y[left]
                + (self.y[right] - self.y[left]) * (value - self.x[left])
                    / (self.x[right] - self.x[left]),
        )
    }
}

#[derive(Clone, Debug)]
pub enum ModelConfig {
    None,
    Identity,
    Linear(LinearOptions),
    Interpolated(InterpolationOptions),
    Lowess(LowessOptions),
}
#[derive(Clone, Debug)]
pub enum TransformationModel {
    None,
    Identity,
    Linear(LinearModel),
    Interpolated(InterpolatedModel),
    Lowess(LowessModel),
}
impl TransformationModel {
    pub fn apply(&self, value: f64) -> Result<f64> {
        match self {
            Self::None | Self::Identity => finite(value),
            Self::Linear(m) => m.apply(value),
            Self::Interpolated(m) => m.apply(value),
            Self::Lowess(m) => m.apply(value),
        }
    }
    pub fn name(&self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Identity => "identity",
            Self::Linear(_) => "linear",
            Self::Interpolated(_) => "interpolated",
            Self::Lowess(_) => "lowess",
        }
    }
    fn fit(data: &[DataPoint], config: &ModelConfig) -> Result<Self> {
        Ok(match config {
            ModelConfig::None => Self::None,
            ModelConfig::Identity => Self::Identity,
            ModelConfig::Linear(o) => Self::Linear(LinearModel::fit(data, *o)?),
            ModelConfig::Interpolated(o) => Self::Interpolated(InterpolatedModel::fit(data, *o)?),
            ModelConfig::Lowess(o) => Self::Lowess(LowessModel::fit(data, *o)?),
        })
    }
}
#[derive(Clone, Debug)]
pub struct TransformationDescription {
    data: Vec<DataPoint>,
    model: TransformationModel,
    config: ModelConfig,
    max_points: usize,
}
impl Default for TransformationDescription {
    fn default() -> Self {
        Self {
            data: Vec::new(),
            model: TransformationModel::None,
            config: ModelConfig::None,
            max_points: 1_000_000,
        }
    }
}
impl TransformationDescription {
    pub fn new(data: Vec<DataPoint>) -> Result<Self> {
        Self::with_max_points(data, 1_000_000)
    }
    pub fn with_max_points(data: Vec<DataPoint>, max_points: usize) -> Result<Self> {
        validate_points(&data, max_points)?;
        Ok(Self {
            data,
            max_points,
            ..Default::default()
        })
    }
    pub fn data_points(&self) -> &[DataPoint] {
        &self.data
    }
    pub fn model(&self) -> &TransformationModel {
        &self.model
    }
    pub fn model_config(&self) -> &ModelConfig {
        &self.config
    }
    /// Replacing anchors resets the model, including the identity refit lock.
    pub fn set_data_points(&mut self, data: Vec<DataPoint>) -> Result<()> {
        validate_points(&data, self.max_points)?;
        self.data = data;
        self.model = TransformationModel::None;
        self.config = ModelConfig::None;
        Ok(())
    }
    /// Identity models ignore refitting requests until data are reset.
    pub fn fit_model(&mut self, config: ModelConfig) -> Result<()> {
        if matches!(self.model, TransformationModel::Identity) {
            return Ok(());
        }
        let model = TransformationModel::fit(&self.data, &config)?;
        self.model = model;
        self.config = config;
        Ok(())
    }
    pub fn apply(&self, value: f64) -> Result<f64> {
        self.model.apply(value)
    }
    /// Atomic batch application; failures do not overwrite any input value.
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
    /// Swap source/target anchors, preserve notes, and refit with the same options.
    /// With no anchors, a linear model is inverted algebraically instead.
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
                result.config = ModelConfig::Linear(inverse.options());
                result.model = TransformationModel::Linear(inverse);
                return Ok(result);
            }
        }
        result.fit_model(self.config.clone())?;
        Ok(result)
    }
    pub fn invert(&mut self) -> Result<()> {
        *self = self.inverse()?;
        Ok(())
    }
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
    /// Adaptive residual quantile in transformed units, or original units after
    /// source-style refitting of the inverse. Defaults are in WindowOptions.
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
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TransformationStatistics {
    pub x_range: Option<NumericRange>,
    pub y_range: Option<NumericRange>,
    pub percentiles: Vec<DeviationPercentile>,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DeviationPercentile {
    pub percent: u32,
    pub before: f64,
    pub after: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindowOptions {
    pub quantile: f64,
    pub inverse: bool,
    pub full_window: bool,
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

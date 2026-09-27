// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The B-spline transformation model, `TransformationModelBSpline`.
//!
//! Ports `ANALYSIS/MAPMATCHING/TransformationModelBSpline.h` and its `.cpp`.
//! The smoothing itself is [`BSpline2d`], this crate's port of
//! `MATH/MISC/BSpline2d`; this module adds what the model class adds: the
//! data-range check on the cutoff wavelength and the four extrapolation
//! methods beyond the data range. See `docs/TRANSFORMATIONS_SUPPORT.md`.

use super::{DataPoint, LinearModel, LinearOptions, bad, finite, validate_points};
use crate::Result;
use crate::processing::spline::{BSpline2d, BoundaryCondition};

/// Method used beyond the data range (source parameter `extrapolate`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum BSplineExtrapolation {
    /// `linear`: continue with the slope of the B-spline at the corresponding
    /// endpoint. The default.
    #[default]
    Linear,
    /// `b_spline`: use the B-spline itself, as for interpolation.
    BSpline,
    /// `constant`: use the value of the B-spline at the corresponding endpoint.
    Constant,
    /// `global_linear`: use a least-squares line through all data points,
    /// anchored at the endpoints; the source warns this most probably
    /// introduces discontinuities at the ends of the data range.
    GlobalLinear,
}

impl BSplineExtrapolation {
    /// The source parameter value naming this method.
    pub fn name(self) -> &'static str {
        match self {
            Self::Linear => "linear",
            Self::BSpline => "b_spline",
            Self::Constant => "constant",
            Self::GlobalLinear => "global_linear",
        }
    }

    /// The method a source `extrapolate` value selects.
    ///
    /// Source `TransformationModelBSpline` compares against `b_spline`,
    /// `global_linear` and `constant` and takes linear extrapolation for
    /// **every other string**, so an unknown name is not an error there. This
    /// reproduces that: the parameter's valid-string restriction is the only
    /// place an unknown name is refused, and a caller that validates against
    /// [`BSplineOptions::defaults`] gets that refusal.
    pub fn from_name(name: &str) -> Self {
        match name {
            "b_spline" => Self::BSpline,
            "global_linear" => Self::GlobalLinear,
            "constant" => Self::Constant,
            _ => Self::Linear,
        }
    }
}

/// Options of the B-spline model (source `getDefaultParameters`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BSplineOptions {
    /// Cutoff wavelength in the units of the data (source `wavelength`,
    /// default 0, minimum 0). It sets the node count so that the spline
    /// approximates a low-pass filter with this cutoff; a higher value means
    /// more smoothing. The source documents 0 as "twice the number of input
    /// points"; see [`BSpline2d::with_options`] for what zero actually does.
    pub wavelength: f64,
    /// Number of nodes (source `num_nodes`, default 5, minimum 0). Two or
    /// more overrides `wavelength`; a lower value means more smoothing.
    pub num_nodes: usize,
    /// Method used beyond the data range (source `extrapolate`).
    pub extrapolation: BSplineExtrapolation,
    /// Condition at the spline endpoints (source `boundary_condition`, 0, 1
    /// or 2, default 2: second derivative zero).
    pub boundary_condition: BoundaryCondition,
    /// Maximum number of data points (native limit).
    pub max_points: usize,
}

impl Default for BSplineOptions {
    fn default() -> Self {
        Self {
            wavelength: 0.0,
            num_nodes: 5,
            extrapolation: BSplineExtrapolation::Linear,
            boundary_condition: BoundaryCondition::ZeroSecond,
            max_points: 1_000_000,
        }
    }
}

/// A B-spline smoothing model with extrapolation beyond the data range
/// (source `TransformationModelBSpline`).
#[derive(Clone, Debug)]
pub struct BSplineModel {
    spline: BSpline2d,
    xmin: f64,
    xmax: f64,
    offset_min: f64,
    offset_max: f64,
    slope_min: f64,
    slope_max: f64,
    options: BSplineOptions,
}

impl BSplineModel {
    /// Fit the spline to the data points (the source constructor).
    ///
    /// The points need not be sorted; their smallest and largest x set the
    /// node grid and the extrapolation origins. For `global_linear`, the
    /// line is the ordinary least-squares [`LinearModel`] over all points,
    /// evaluated at the two ends so that extrapolation starts from the line
    /// rather than from the spline.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for fewer than two points (the
    /// source's `Exception::IllegalArgument`, "'b_spline' model requires more
    /// data"), for a `wavelength` larger than the data range (the source's
    /// `Exception::InvalidValue`), when the spline cannot be fitted (the
    /// source's `Exception::UnableToFit`; see [`BSpline2d::with_options`] for
    /// the conditions), for a failed `global_linear` fit, or for non-finite
    /// input or arithmetic.
    pub fn fit(data: &[DataPoint], options: BSplineOptions) -> Result<Self> {
        validate_points(data, options.max_points)?;
        if data.len() < 2 {
            return Err(bad("'b_spline' model requires more data"));
        }
        if !options.wavelength.is_finite() || options.wavelength < 0.0 {
            return Err(bad("B-spline wavelength must be finite and non-negative"));
        }
        let x: Vec<f64> = data.iter().map(|p| p.x).collect();
        let y: Vec<f64> = data.iter().map(|p| p.y).collect();
        let (mut xmin, mut xmax) = (x[0], x[0]);
        for &value in &x {
            if value < xmin {
                xmin = value;
            } else if value > xmax {
                xmax = value;
            }
        }
        if options.wavelength > finite(xmax - xmin)? {
            return Err(bad(
                "B-spline 'wavelength' can't be larger than the data range",
            ));
        }
        let spline = BSpline2d::with_options(
            &x,
            &y,
            options.wavelength,
            options.boundary_condition,
            options.num_nodes,
        )
        .map_err(|error| bad(&format!("Unable to fit B-spline to data points: {error}")))?;
        if !spline.ok() {
            return Err(bad("Unable to fit B-spline to data points"));
        }
        let (mut offset_min, mut offset_max, mut slope_min, mut slope_max) = (0.0, 0.0, 0.0, 0.0);
        match options.extrapolation {
            BSplineExtrapolation::BSpline => {}
            BSplineExtrapolation::GlobalLinear => {
                let line = LinearModel::fit(
                    data,
                    LinearOptions {
                        max_points: options.max_points,
                        ..LinearOptions::default()
                    },
                )?;
                slope_min = line.coefficients().slope;
                slope_max = slope_min;
                // Extrapolation to either side takes xmin / xmax as its origin.
                offset_min = line.apply(xmin)?;
                offset_max = line.apply(xmax)?;
            }
            BSplineExtrapolation::Linear | BSplineExtrapolation::Constant => {
                offset_min = spline.eval(xmin)?;
                offset_max = spline.eval(xmax)?;
                if options.extrapolation == BSplineExtrapolation::Linear {
                    slope_min = spline.derivative(xmin)?;
                    slope_max = spline.derivative(xmax)?;
                }
            }
        }
        Ok(Self {
            spline,
            xmin,
            xmax,
            offset_min,
            offset_max,
            slope_min,
            slope_max,
            options,
        })
    }

    /// The options the model was fitted with.
    pub fn options(&self) -> BSplineOptions {
        self.options
    }

    /// The smallest and largest x of the fitted data, where extrapolation
    /// takes over.
    pub fn data_range(&self) -> (f64, f64) {
        (self.xmin, self.xmax)
    }

    /// The fitted spline.
    pub fn spline(&self) -> &BSpline2d {
        &self.spline
    }

    /// Evaluate the model at `value` (source `evaluate`).
    ///
    /// Below the smallest and above the largest fitted x the chosen
    /// extrapolation applies, unless it is [`BSplineExtrapolation::BSpline`].
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for a non-finite `value` or
    /// result.
    pub fn apply(&self, value: f64) -> Result<f64> {
        finite(value)?;
        let spline_ends = self.options.extrapolation == BSplineExtrapolation::BSpline;
        if value < self.xmin && !spline_ends {
            return finite(
                if self.options.extrapolation == BSplineExtrapolation::Constant {
                    self.offset_min
                } else {
                    self.offset_min - self.slope_min * (self.xmin - value)
                },
            );
        }
        if value > self.xmax && !spline_ends {
            return finite(
                if self.options.extrapolation == BSplineExtrapolation::Constant {
                    self.offset_max
                } else {
                    self.offset_max + self.slope_max * (value - self.xmax)
                },
            );
        }
        self.spline.eval(value)
    }
}

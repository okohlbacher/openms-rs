// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The distance between two features or consensus features that feature
//! grouping minimises.
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/FeatureDistance.h`. See
//! `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! The distance is
//!
//! ```text
//! (w_RT  * (|RT_1  - RT_2|  / dRT_max)^p_RT
//!  + w_MZ  * (|MZ_1  - MZ_2|  / dMZ_max)^p_MZ
//!  + w_int * (|int_1 - int_2| / int_max)^p_int) / (w_RT + w_MZ + w_int)
//! ```
//!
//! which lies between zero and one unless a constraint is violated. `dRT_max`
//! and `dMZ_max` are the `distance_RT:max_difference` and
//! `distance_MZ:max_difference` parameters and normalise the observed
//! differences; `int_max` is the constructor's `max_intensity`, not a user
//! parameter, because it depends on the data. `p_X` is `distance_X:exponent`,
//! `w_X` is `distance_X:weight`; the intensity weight defaults to zero, so
//! intensity is ignored unless asked for, though `int_max` should still be set
//! sensibly.
//!
//! When an absolute difference exceeds its maximum, the result depends on
//! `force_constraints`: without it, that dimension may exceed one and the pair
//! is flagged invalid; with it, the distance is [`FeatureDistance::INFINITY`](crate::analysis::mapmatching::feature_distance::FeatureDistance::INFINITY).
//!
//! By default two features pair only when their charges are equal or at least
//! one is unknown (`0`); otherwise the distance is infinite. `ignore_charge`
//! lifts that. Peptide identifications are not taken into account: they are
//! stored in a form unsuitable for rapid comparison.

use super::params::{add_float, add_string, bad, flag, flag_text, float, string};
use crate::chemistry::EmpiricalFormula;
use crate::concept::constants::user_param::DC_CHARGE_ADDUCTS;
use crate::error::Result;
use crate::kernel::features::BaseFeature;
use crate::param::{DefaultParamHandler, Param};

/// Parameter-handler name of the source class.
pub const FEATURE_DISTANCE_NAME: &str = "FeatureDistance";

/// One distance component: RT, m/z or intensity (the user-settable part of
/// source `DistanceParams_`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DistanceComponent {
    /// Normalising maximum difference. For intensity it is the constructor's
    /// `max_intensity`, or its `log10(x + 1)` under log transformation, as the
    /// source writes into `distance_intensity:max_difference`.
    pub max_difference: f64,
    /// Exponent the normalised difference is raised to. `1` and `2` take the
    /// source's multiplication shortcuts; anything else uses `powf`.
    pub exponent: f64,
    /// Weight of the component.
    pub weight: f64,
}

impl DistanceComponent {
    /// The source's `relevant`: a non-zero weight and a non-zero exponent. An
    /// irrelevant component has its weight treated as zero, in the distance
    /// and in the total weight.
    fn effective_weight(&self) -> f64 {
        if self.weight != 0.0 && self.exponent != 0.0 {
            self.weight
        } else {
            0.0
        }
    }

    /// Source `distance_`: the weighted, normalised, exponentiated difference.
    /// `norm_factor` is `1 / max_difference`, multiplied rather than divided,
    /// as in the source.
    fn distance(&self, diff: f64, norm_factor: f64) -> f64 {
        let weight = self.effective_weight();
        if self.exponent == 1.0 {
            diff * norm_factor * weight
        } else if self.exponent == 2.0 {
            let tmp = diff * norm_factor;
            tmp * tmp * weight
        } else {
            (diff * norm_factor).powf(self.exponent) * weight
        }
    }
}

/// Feature distance functor (source `FeatureDistance`).
#[derive(Clone, Debug, PartialEq)]
pub struct FeatureDistance {
    /// RT component (`distance_RT:*`, defaults 100 s, exponent 1, weight 1).
    pub rt: DistanceComponent,
    /// m/z component (`distance_MZ:*`, defaults 0.3, exponent 2, weight 1).
    pub mz: DistanceComponent,
    /// Whether the m/z maximum is in ppm of the *left* feature's m/z
    /// (`distance_MZ:unit` = `ppm`) rather than in Da.
    pub mz_ppm: bool,
    /// Intensity exponent (`distance_intensity:exponent`, default 1).
    pub intensity_exponent: f64,
    /// Intensity weight (`distance_intensity:weight`, default 0).
    pub intensity_weight: f64,
    /// Compare `log10(intensity + 1)` instead of intensities
    /// (`distance_intensity:log_transform`).
    pub log_transform: bool,
    /// Pair irrespective of charge (`ignore_charge`, default false).
    pub ignore_charge: bool,
    /// Pair irrespective of the `dc_charge_adducts` annotation
    /// (`ignore_adduct`, default true).
    pub ignore_adduct: bool,
    /// Maximum feature intensity, for normalisation (constructor argument).
    pub max_intensity: f64,
    /// Return [`Self::INFINITY`] when a maximum difference is exceeded
    /// (constructor argument).
    pub force_constraints: bool,
}

impl Default for FeatureDistance {
    /// Default parameters with `max_intensity` 1 and no forced constraints, as
    /// the source's default constructor.
    fn default() -> Self {
        Self::new(1.0, false)
    }
}

impl FeatureDistance {
    /// Returned when a maximum difference is exceeded under forced
    /// constraints, or when charges or adducts do not match (source
    /// `infinity`).
    pub const INFINITY: f64 = f64::INFINITY;

    /// Default parameters with the given normalising intensity and constraint
    /// mode (source constructor).
    pub fn new(max_intensity: f64, force_constraints: bool) -> Self {
        Self {
            rt: DistanceComponent {
                max_difference: 100.0,
                exponent: 1.0,
                weight: 1.0,
            },
            mz: DistanceComponent {
                max_difference: 0.3,
                exponent: 2.0,
                weight: 1.0,
            },
            mz_ppm: false,
            intensity_exponent: 1.0,
            intensity_weight: 0.0,
            log_transform: false,
            ignore_charge: false,
            ignore_adduct: true,
            max_intensity,
            force_constraints,
        }
    }

    /// The source defaults tree.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters in the layout of [`Self::defaults`], carrying the
    /// current values.
    ///
    /// The source also stores the derived `distance_intensity:max_difference`
    /// (from `max_intensity`) in its current parameters; it is not a default
    /// and not user-settable, so it is left out here and reported by
    /// [`Self::intensity_component`] instead.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut p = Param::new();
        add_float(
            &mut p,
            "distance_RT:max_difference",
            self.rt.max_difference,
            "Never pair features with a larger RT distance (in seconds).",
            false,
            Some(0.0),
            None,
        )?;
        add_float(
            &mut p,
            "distance_RT:exponent",
            self.rt.exponent,
            "Normalized RT differences ([0-1], relative to 'max_difference') are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)",
            true,
            Some(0.0),
            None,
        )?;
        add_float(
            &mut p,
            "distance_RT:weight",
            self.rt.weight,
            "Final RT distances are weighted by this factor",
            true,
            Some(0.0),
            None,
        )?;
        p.set_section_description("distance_RT", "Distance component based on RT differences")?;
        add_float(
            &mut p,
            "distance_MZ:max_difference",
            self.mz.max_difference,
            "Never pair features with larger m/z distance (unit defined by 'unit')",
            false,
            Some(0.0),
            None,
        )?;
        add_string(
            &mut p,
            "distance_MZ:unit",
            if self.mz_ppm { "ppm" } else { "Da" },
            "Unit of the 'max_difference' parameter",
            false,
            &["Da", "ppm"],
        )?;
        add_float(
            &mut p,
            "distance_MZ:exponent",
            self.mz.exponent,
            "Normalized ([0-1], relative to 'max_difference') m/z differences are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)",
            true,
            Some(0.0),
            None,
        )?;
        add_float(
            &mut p,
            "distance_MZ:weight",
            self.mz.weight,
            "Final m/z distances are weighted by this factor",
            true,
            Some(0.0),
            None,
        )?;
        p.set_section_description("distance_MZ", "Distance component based on m/z differences")?;
        add_float(
            &mut p,
            "distance_intensity:exponent",
            self.intensity_exponent,
            "Differences in relative intensity ([0-1]) are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)",
            true,
            Some(0.0),
            None,
        )?;
        add_float(
            &mut p,
            "distance_intensity:weight",
            self.intensity_weight,
            "Final intensity distances are weighted by this factor",
            true,
            Some(0.0),
            None,
        )?;
        add_string(
            &mut p,
            "distance_intensity:log_transform",
            if self.log_transform {
                "enabled"
            } else {
                "disabled"
            },
            "Log-transform intensities? If disabled, d = |int_f2 - int_f1| / int_max. If enabled, d = |log(int_f2 + 1) - log(int_f1 + 1)| / log(int_max + 1))",
            true,
            &["enabled", "disabled"],
        )?;
        p.set_section_description("distance_intensity", "Distance component based on differences in relative intensity (usually relative to highest peak in the whole data set)")?;
        add_string(
            &mut p,
            "ignore_charge",
            flag_text(self.ignore_charge),
            "false [default]: pairing requires equal charge state (or at least one unknown charge '0'); true: Pairing irrespective of charge state",
            false,
            &["true", "false"],
        )?;
        add_string(
            &mut p,
            "ignore_adduct",
            flag_text(self.ignore_adduct),
            "true [default]: pairing requires equal adducts (or at least one without adduct annotation); true: Pairing irrespective of adducts",
            false,
            &["true", "false"],
        )?;
        Ok(p)
    }

    /// Build from parameters with `max_intensity` 1 and no forced constraints,
    /// as the default constructor followed by `setParameters`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::error::Error::InvalidValue`] when a parameter has the
    /// wrong type or violates its restriction (source
    /// `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`Self::from_param`], also returning the unknown-parameter warnings
    /// the source logs.
    ///
    /// # Errors
    ///
    /// As [`Self::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        Self::with_param(1.0, false, parameters)
    }

    /// Build from parameters with the given constructor arguments.
    ///
    /// # Errors
    ///
    /// As [`Self::from_param`].
    pub fn with_param(
        max_intensity: f64,
        force_constraints: bool,
        parameters: &Param,
    ) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(FEATURE_DISTANCE_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |p| {
            Self::from_complete_param(max_intensity, force_constraints, p)
        })
    }

    /// Read a validated tree whose every default is present.
    pub(crate) fn from_complete_param(
        max_intensity: f64,
        force_constraints: bool,
        p: &Param,
    ) -> Result<Self> {
        let component = |what: &str| -> Result<DistanceComponent> {
            Ok(DistanceComponent {
                max_difference: float(p, &format!("distance_{what}:max_difference"))?,
                exponent: float(p, &format!("distance_{what}:exponent"))?,
                weight: float(p, &format!("distance_{what}:weight"))?,
            })
        };
        let unit = string(p, "distance_MZ:unit")?;
        let log_transform = match string(p, "distance_intensity:log_transform")?.as_str() {
            "enabled" => true,
            "disabled" => false,
            other => {
                return Err(bad(format!(
                    "invalid distance_intensity:log_transform '{other}'"
                )));
            }
        };
        Ok(Self {
            rt: component("RT")?,
            mz: component("MZ")?,
            mz_ppm: unit == "ppm",
            intensity_exponent: float(p, "distance_intensity:exponent")?,
            intensity_weight: float(p, "distance_intensity:weight")?,
            log_transform,
            ignore_charge: flag(p, "ignore_charge")?,
            ignore_adduct: flag(p, "ignore_adduct")?,
            max_intensity,
            force_constraints,
        })
    }

    /// The intensity component as the source derives it: its maximum
    /// difference is `max_intensity`, or `log10(max_intensity + 1)` under log
    /// transformation.
    pub fn intensity_component(&self) -> DistanceComponent {
        DistanceComponent {
            max_difference: if self.log_transform {
                linear2log(self.max_intensity)
            } else {
                self.max_intensity
            },
            exponent: self.intensity_exponent,
            weight: self.intensity_weight,
        }
    }

    /// Check the constraints and compute the distance between two features
    /// (source `operator()`).
    ///
    /// Returns whether the constraints were satisfied and the distance, which
    /// is [`Self::INFINITY`] when charges differ (both known) and
    /// `ignore_charge` is off, when both features carry a `dc_charge_adducts`
    /// annotation with different element compositions and `ignore_adduct` is
    /// off, or when a maximum difference is exceeded under
    /// `force_constraints`. An exceeded maximum without `force_constraints`
    /// only makes the first value `false`.
    ///
    /// With a ppm m/z unit the maximum is `max_difference * left_mz * 1e-6`:
    /// relative to the **left** feature, so the distance is not symmetric.
    /// The source stores that per-call normalisation in a member, which makes
    /// its `operator()` non-`const`; here it is local.
    ///
    /// The terms are summed as `(rt + mz) + intensity` and multiplied by the
    /// reciprocal of the effective total weight, as in the source.
    ///
    /// # Errors
    ///
    /// Returns the [`EmpiricalFormula::parse`] error when an adduct annotation
    /// is compared and does not parse, or is not a string; the source throws
    /// `Exception::ParseError` from the `EmpiricalFormula` constructor.
    pub fn distance(&self, left: &BaseFeature, right: &BaseFeature) -> Result<(bool, f64)> {
        if !self.ignore_charge {
            let (l, r) = (left.charge, right.charge);
            if l != r && l != 0 && r != 0 {
                return Ok((false, Self::INFINITY));
            }
        }
        if !self.ignore_adduct {
            if let (Some(l), Some(r)) = (
                left.metadata.get(DC_CHARGE_ADDUCTS),
                right.metadata.get(DC_CHARGE_ADDUCTS),
            ) {
                if EmpiricalFormula::parse(l.as_str()?)? != EmpiricalFormula::parse(r.as_str()?)? {
                    return Ok((false, Self::INFINITY));
                }
            }
        }

        let mut valid = true;
        let mut dist_mz = (left.mz - right.mz).abs();
        let mut max_diff_mz = self.mz.max_difference;
        let mut mz_norm = 1.0 / self.mz.max_difference;
        if self.mz_ppm {
            max_diff_mz *= left.mz * 1e-6;
            mz_norm = 1.0 / max_diff_mz;
        }
        if dist_mz > max_diff_mz {
            if self.force_constraints {
                return Ok((false, Self::INFINITY));
            }
            valid = false;
        }

        let mut dist_rt = (left.rt - right.rt).abs();
        if dist_rt > self.rt.max_difference {
            if self.force_constraints {
                return Ok((false, Self::INFINITY));
            }
            valid = false;
        }

        dist_rt = self.rt.distance(dist_rt, 1.0 / self.rt.max_difference);
        dist_mz = self.mz.distance(dist_mz, mz_norm);

        let intensity = self.intensity_component();
        let mut dist_intensity = 0.0;
        if intensity.effective_weight() != 0.0 {
            dist_intensity = if self.log_transform {
                (linear2log(f64::from(left.intensity)) - linear2log(f64::from(right.intensity)))
                    .abs()
            } else {
                // `fabs(float - float)`: the difference is taken in `f32`.
                f64::from((left.intensity - right.intensity).abs())
            };
            dist_intensity = intensity.distance(dist_intensity, 1.0 / intensity.max_difference);
        }

        let total_weight_reciprocal = 1.0
            / (self.rt.effective_weight()
                + self.mz.effective_weight()
                + intensity.effective_weight());
        let dist = (dist_rt + dist_mz + dist_intensity) * total_weight_reciprocal;
        Ok((valid, dist))
    }
}

/// Source `Math::linear2log`: `log10(x + 1)`.
fn linear2log(x: f64) -> f64 {
    (x + 1.0).log10()
}

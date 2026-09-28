// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Fitting (and validating) elution models for features.
//!
//! Ports `FEATUREFINDER/ElutionModelFitter.h` and its `.cpp`, the helper that
//! `FeatureFinderMetaboIdent` (and `FeatureFinderIdentification`) use to fit a
//! Gaussian or EGH retention-time model to the mass traces of every feature,
//! validate it, and replace the feature intensity by the model area. See
//! `docs/FEATUREFINDER_MODELS_SUPPORT.md`.
//!
//! The trace fitters are the picked feature finder's
//! [`GaussTraceFitter`] and [`EGHTraceFitter`]; the imputation line is
//! [`crate::analysis::transformations::LinearModel`], the port of
//! `TransformationModelLinear`.

use crate::Result;
use crate::analysis::feature_finder_picked::egh_trace_fitter::EGHTraceFitter;
use crate::analysis::feature_finder_picked::gauss_trace_fitter::GaussTraceFitter;
use crate::analysis::feature_finder_picked::helper_structs::{MassTrace, MassTraces, TracePeak};
use crate::analysis::feature_finder_picked::trace_fitter::{
    TraceFitter, TraceFitterParams, UNABLE_TO_FIT_FINAL_SET,
};
use crate::analysis::transformations::{DataPoint, LinearModel, LinearOptions};
use crate::error::Error;
use crate::kernel::features::{Feature, FeatureMap};
use crate::math::statistic_functions::median;
use crate::metadata::MetaValue;
use crate::param::{Param, ParamValue};

use super::base_model::{complete_parameters, float_value, set_default, set_valid};

/// Status text of a valid model (meta value `model_status`).
pub const STATUS_VALID: &str = "0 (valid)";

fn missing(message: &str) -> Error {
    Error::MissingInformation(message.to_owned())
}

/// A float meta value as the source stores it, NaN and infinities included.
fn float_meta(value: f64) -> MetaValue {
    MetaValue::source_float(value)
}

/// A numeric meta value as the source's `double(getMetaValue(key))` reads it.
fn meta_f64(feature: &Feature, key: &str) -> Result<f64> {
    feature
        .base
        .metadata
        .get(key)
        .ok_or_else(|| {
            Error::MissingInformation(format!(
                "feature {} has no meta value '{key}'",
                feature.base.unique_id
            ))
        })?
        .as_f64()
}

fn status_of(feature: &Feature) -> Option<&str> {
    feature
        .base
        .metadata
        .get("model_status")
        .and_then(|v| v.as_str().ok())
}

/// What a [`ElutionModelFitter::fit_elution_models`] run reports, in place of
/// the source's log lines.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ElutionModelFitReport {
    /// Features whose model is valid (source log "successes").
    pub successes: usize,
    /// Features whose model is not valid ("failures").
    pub failures: usize,
    /// The imputation line `log(area) = slope * log(raw_intensity) +
    /// intercept`, when imputation ran.
    pub imputation: Option<(f64, f64)>,
    /// The source's `OPENMS_LOG_ERROR` lines for fits that threw
    /// `UnableToFit`: `Error fitting model to feature '<id>': <name> - <what>`.
    pub errors: Vec<String>,
    /// Whether every feature was removed because no model was valid.
    pub cleared: bool,
}

/// Fits elution models to features and validates them (source
/// `ElutionModelFitter`, a `DefaultParamHandler`).
#[derive(Clone, Debug, PartialEq)]
pub struct ElutionModelFitter {
    /// Parameter `asymmetric` (default false): fit the exponential-Gaussian
    /// hybrid instead of a Gaussian.
    pub asymmetric: bool,
    /// Parameter `add_zeros` (default 0.2, at least 0): weight of the two
    /// zero-intensity points added outside the feature range; 0 disables them.
    pub add_zeros: f64,
    /// Parameter `unweighted_fit` (default false): do not weight mass traces by
    /// their theoretical intensities.
    pub unweighted_fit: bool,
    /// Parameter `no_imputation` (default false): set a failed feature's
    /// intensity to zero instead of imputing it.
    pub no_imputation: bool,
    /// Parameter `each_trace` (default false): also fit a model to each mass
    /// trace (stored on the subordinate features).
    pub each_trace: bool,
    /// Parameter `check:min_area` (default 1, at least 0).
    pub min_area: f64,
    /// Parameter `check:boundaries` (default 0.5, 0 to 1).
    pub boundaries: f64,
    /// Parameter `check:width` (default 10, at least 0; 0 disables).
    pub width: f64,
    /// Parameter `check:asymmetry` (default 10, at least 0; 0 disables; EGH
    /// only).
    pub asymmetry: f64,
}

impl Default for ElutionModelFitter {
    /// The source defaults.
    fn default() -> Self {
        Self {
            asymmetric: false,
            add_zeros: 0.2,
            unweighted_fit: false,
            no_imputation: false,
            each_trace: false,
            min_area: 1.0,
            boundaries: 0.5,
            width: 10.0,
            asymmetry: 10.0,
        }
    }
}

fn flag(value: bool) -> ParamValue {
    ParamValue::String(if value { "true" } else { "false" }.to_owned())
}

/// Either trace fitter, as the source's `TraceFitter*`.
enum Fitter {
    Gauss(GaussTraceFitter),
    Egh(EGHTraceFitter),
}

impl Fitter {
    fn as_dyn(&self) -> &dyn TraceFitter {
        match self {
            Self::Gauss(f) => f,
            Self::Egh(f) => f,
        }
    }
    fn as_dyn_mut(&mut self) -> &mut dyn TraceFitter {
        match self {
            Self::Gauss(f) => f,
            Self::Egh(f) => f,
        }
    }
}

impl ElutionModelFitter {
    /// The handler name (source `DefaultParamHandler("ElutionModelFitter")`).
    pub const NAME: &'static str = "ElutionModelFitter";

    /// A fitter with the default parameters (source constructor).
    pub fn new() -> Self {
        Self::default()
    }

    /// The default parameters (source `getDefaults`), with the source's
    /// descriptions, restrictions and `advanced` tags.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The current parameters (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut p = Param::new();
        set_default(
            &mut p,
            "asymmetric",
            flag(self.asymmetric),
            "Fit an asymmetric (exponential-Gaussian hybrid) model? By default a symmetric (Gaussian) model is used.",
            false,
        )?;
        set_valid(&mut p, "asymmetric", &["true", "false"])?;
        set_default(
            &mut p,
            "add_zeros",
            ParamValue::Float(self.add_zeros),
            "Add zero-intensity points outside the feature range to constrain the model fit. This parameter sets the weight given to these points during model fitting; '0' to disable.",
            true,
        )?;
        p.set_min_float("add_zeros", 0.0)?;
        set_default(
            &mut p,
            "unweighted_fit",
            flag(self.unweighted_fit),
            "Suppress weighting of mass traces according to theoretical intensities when fitting elution models",
            true,
        )?;
        set_valid(&mut p, "unweighted_fit", &["true", "false"])?;
        set_default(
            &mut p,
            "no_imputation",
            flag(self.no_imputation),
            "If fitting the elution model fails for a feature, set its intensity to zero instead of imputing a value from the initial intensity estimate",
            true,
        )?;
        set_valid(&mut p, "no_imputation", &["true", "false"])?;
        set_default(
            &mut p,
            "each_trace",
            flag(self.each_trace),
            "Fit elution model to each individual mass trace",
            true,
        )?;
        set_valid(&mut p, "each_trace", &["true", "false"])?;
        set_default(
            &mut p,
            "check:min_area",
            ParamValue::Float(self.min_area),
            "Lower bound for the area under the curve of a valid elution model",
            true,
        )?;
        p.set_min_float("check:min_area", 0.0)?;
        set_default(
            &mut p,
            "check:boundaries",
            ParamValue::Float(self.boundaries),
            "Time points corresponding to this fraction of the elution model height have to be within the data region used for model fitting",
            true,
        )?;
        p.set_min_float("check:boundaries", 0.0)?;
        p.set_max_float("check:boundaries", 1.0)?;
        set_default(
            &mut p,
            "check:width",
            ParamValue::Float(self.width),
            "Upper limit for acceptable widths of elution models (Gaussian or EGH), expressed in terms of modified (median-based) z-scores. '0' to disable. Not applied to individual mass traces (parameter 'each_trace').",
            true,
        )?;
        p.set_min_float("check:width", 0.0)?;
        set_default(
            &mut p,
            "check:asymmetry",
            ParamValue::Float(self.asymmetry),
            "Upper limit for acceptable asymmetry of elution models (EGH only), expressed in terms of modified (median-based) z-scores. '0' to disable. Not applied to individual mass traces (parameter 'each_trace').",
            true,
        )?;
        p.set_min_float("check:asymmetry", 0.0)?;
        p.set_section_description(
            "check",
            "Parameters for checking the validity of elution models (and rejecting them if necessary)",
        )?;
        Ok(p)
    }

    /// Replace the parameters, taking defaults for the missing ones (source
    /// `setParameters`); returns the unknown-parameter warnings.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a wrong type or a value outside its
    /// restriction; the fitter is then unchanged.
    pub fn set_parameters(&mut self, param: &Param) -> Result<Vec<String>> {
        let (c, warnings) = complete_parameters(Self::NAME, Self::defaults()?, param)?;
        let staged = Self {
            asymmetric: c.value("asymmetric")?.to_bool()?,
            add_zeros: float_value(&c, "add_zeros")?,
            unweighted_fit: c.value("unweighted_fit")?.to_bool()?,
            no_imputation: c.value("no_imputation")?.to_bool()?,
            each_trace: c.value("each_trace")?.to_bool()?,
            min_area: float_value(&c, "check:min_area")?,
            boundaries: float_value(&c, "check:boundaries")?,
            width: float_value(&c, "check:width")?,
            asymmetry: float_value(&c, "check:asymmetry")?,
        };
        *self = staged;
        Ok(warnings)
    }

    /// A default fitter after [`ElutionModelFitter::set_parameters`].
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param(param: &Param) -> Result<Self> {
        let mut fitter = Self::default();
        fitter.set_parameters(param)?;
        Ok(fitter)
    }

    /// Mean relative error of the fitted model over the traces' peaks within
    /// the model's RT bounds and the first trace's RT range: `sum(|model *
    /// theoretical_int - intensity| / model) / sum(theoretical_int)` (source
    /// `calculateFitQuality_`).
    fn calculate_fit_quality(fitter: &dyn TraceFitter, traces: &MassTraces) -> f64 {
        let Some(first) = traces.get(0) else {
            return f64::NAN;
        };
        let (Some(first_peak), Some(last_peak)) = (first.peaks.first(), first.peaks.last()) else {
            return f64::NAN;
        };
        let rt_start = fitter.lower_rt_bound().max(first_peak.rt);
        let rt_end = fitter.upper_rt_bound().min(last_peak.rt);
        let mut mre = 0.0;
        let mut total_weights = 0.0;
        for trace in traces.iter() {
            for peak in &trace.peaks {
                let rt = peak.rt;
                if rt >= rt_start && rt <= rt_end {
                    let model_value = fitter.value(rt);
                    let diff =
                        (model_value * trace.theoretical_int - f64::from(peak.intensity)).abs();
                    mre += diff / model_value;
                    total_weights += trace.theoretical_int;
                }
            }
        }
        mre / total_weights
    }

    /// Fit one model and record its parameters, fit quality and validity on
    /// `feature` (source `fitAndValidateModel_`).
    fn fit_and_validate(
        &self,
        fitter: &mut Fitter,
        traces: &MassTraces,
        feature: &mut Feature,
        region: (f64, f64),
        report: &mut ElutionModelFitReport,
    ) -> Result<()> {
        let (region_start, region_end) = region;
        let mut fit_success = true;
        if let Err(error) = fitter.as_dyn_mut().fit(traces) {
            match &error {
                Error::InvalidValue(message) if message.starts_with(UNABLE_TO_FIT_FINAL_SET) => {
                    report.errors.push(format!(
                        "Error fitting model to feature '{}': {}",
                        feature.base.unique_id,
                        message.replacen(": ", " - ", 1)
                    ));
                    fit_success = false;
                }
                _ => return Err(error),
            }
        }
        let f = fitter.as_dyn();
        let center = f.center();
        let height = f.height();
        let meta = &mut feature.base.metadata;
        meta.insert("model_height".into(), float_meta(height));
        meta.insert("model_FWHM".into(), float_meta(f.fwhm()));
        meta.insert("model_center".into(), float_meta(center));
        meta.insert("model_lower".into(), float_meta(f.lower_rt_bound()));
        meta.insert("model_upper".into(), float_meta(f.upper_rt_bound()));
        match fitter {
            Fitter::Egh(egh) => {
                let sigma = egh.sigma();
                let tau = egh.tau();
                meta.insert("model_EGH_tau".into(), float_meta(tau));
                meta.insert("model_EGH_sigma".into(), float_meta(sigma));
                // "see implementation of EGHTraceFitter::getArea"
                meta.insert(
                    "model_width".into(),
                    float_meta(sigma * 0.6266571 + tau.abs()),
                );
                meta.insert("model_asymmetry".into(), float_meta(tau.abs() / sigma));
            }
            Fitter::Gauss(gauss) => {
                let sigma = gauss.sigma();
                meta.insert("model_Gauss_sigma".into(), float_meta(sigma));
                meta.insert("model_width".into(), float_meta(sigma));
            }
        }
        let f = fitter.as_dyn();
        let mre = if fit_success {
            Self::calculate_fit_quality(f, traces)
        } else {
            -1.0
        };
        let meta = &mut feature.base.metadata;
        meta.insert("model_error".into(), float_meta(mre));
        let area = f.area();
        meta.insert("model_area".into(), float_meta(area));
        let status = if area.is_nan() || area <= self.min_area {
            "1 (invalid area)"
        } else if center <= region_start || center >= region_end {
            "2 (center out of bounds)"
        } else if f.value(region_start) > self.boundaries * height {
            "3 (left side out of bounds)"
        } else if f.value(region_end) > self.boundaries * height {
            "4 (right side out of bounds)"
        } else {
            STATUS_VALID
        };
        meta.insert("model_status".into(), MetaValue::from(status));
        Ok(())
    }

    /// Fit elution models to all features and validate them (source
    /// `fitElutionModels`).
    ///
    /// The source states its assumptions as unchecked: every feature has the
    /// meta values `leftWidth` and `rightWidth` (the RT region), subordinates
    /// (one per mass trace or transition) with a meta value
    /// `isotope_probability` and one convex hull each, and the hull points'
    /// y coordinates are intensities. Points of zero intensity are not used.
    ///
    /// Per feature, the traces (and, with `add_zeros`, one trace of two zero
    /// points 20% of the region beyond either end, weighted `add_zeros`) are
    /// fitted with a Gaussian or EGH trace fitter - weighted by theoretical
    /// intensity unless `unweighted_fit` - and the model's height, FWHM,
    /// center, bounds, sigma (and tau, width, asymmetry), mean relative
    /// error, area and status are stored as meta values. If no model is valid
    /// every feature is removed. Otherwise valid models whose width (and, for
    /// EGH, asymmetry) has a modified z-score beyond the limits are rejected;
    /// every feature gets `raw_intensity`; valid ones take the model area as
    /// intensity, failed ones zero or - by default - an imputed
    /// `exp(slope * log(raw) + intercept)` from a linear fit of `log(area)`
    /// on `log(raw)` over the valid ones.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] for no features, a feature without
    /// subordinates or its first subordinate without a hull (the source's
    /// exceptions), and for a missing `leftWidth`, `rightWidth` or
    /// `isotope_probability` (the source's `ConversionError` on an empty
    /// value); errors of the trace fitters other than `UnableToFit`, which is
    /// caught and reported as the source logs it; and a refused imputation
    /// line (a non-finite logarithm, which the source carries along as
    /// infinity). `features` is unchanged on every error.
    pub fn fit_elution_models(&self, features: &mut FeatureMap) -> Result<ElutionModelFitReport> {
        if features.features.is_empty() {
            return Err(missing("No features provided."));
        }
        let mut report = ElutionModelFitReport::default();
        let mut staged = features.clone();
        // Source: the fitter's defaults, with `weighted` "true" unless
        // `unweighted_fit`.
        let params = TraceFitterParams {
            weighted: !self.unweighted_fit,
            ..TraceFitterParams::default()
        };
        let mut fitter = if self.asymmetric {
            let mut f = EGHTraceFitter::new();
            f.set_parameters(params);
            Fitter::Egh(f)
        } else {
            let mut f = GaussTraceFitter::new();
            f.set_parameters(params);
            Fitter::Gauss(f)
        };

        for feat in &mut staged.features {
            let region = (meta_f64(feat, "leftWidth")?, meta_f64(feat, "rightWidth")?);
            let Some(first_sub) = feat.subordinates.first() else {
                return Err(missing(
                    "No subordinate features for mass traces available.",
                ));
            };
            if first_sub.convex_hulls.is_empty() {
                return Err(missing(
                    "No hull points for mass trace in subordinate feature available.",
                ));
            }
            let first_mz = first_sub.base.mz;
            let mut traces = MassTraces::default();
            traces.reserve(feat.subordinates.len() + usize::from(self.add_zeros > 0.0))?;
            for sub in &mut feat.subordinates {
                let hull = sub.convex_hulls.first().ok_or_else(|| {
                    missing("No hull points for mass trace in subordinate feature available.")
                })?;
                let mut trace = MassTrace::default();
                for (k, point) in hull.hull_points().iter().enumerate() {
                    let intensity = point.mz;
                    // "only use non-zero intensities for fitting"
                    if intensity > 0.0 {
                        trace.peaks.push(TracePeak::new(
                            0,
                            k,
                            point.rt,
                            sub.base.mz,
                            intensity as f32,
                        ));
                    }
                }
                trace.update_maximum();
                if trace.peaks.is_empty() {
                    continue;
                }
                if self.each_trace {
                    let mut temp = MassTraces::default();
                    let mut single = trace.clone();
                    single.theoretical_int = 1.0;
                    temp.push(single);
                    temp.max_trace = 0;
                    self.fit_and_validate(&mut fitter, &temp, sub, region, &mut report)?;
                }
                trace.theoretical_int = meta_f64(sub, "isotope_probability")?;
                traces.push(trace);
            }
            // "find the trace with maximal intensity"
            let mut max_trace = 0usize;
            let mut max_intensity = 0.0f64;
            for (i, trace) in traces.iter().enumerate() {
                if let Some(peak) = trace.max_peak {
                    if f64::from(peak.intensity) > max_intensity {
                        max_trace = i;
                        max_intensity = f64::from(peak.intensity);
                    }
                }
            }
            traces.max_trace = max_trace;
            traces.baseline = 0.0;
            if self.add_zeros > 0.0 {
                let (region_start, region_end) = region;
                let offset = 0.2 * (region_start - region_end);
                let trace = MassTrace {
                    theoretical_int: self.add_zeros,
                    max_rt: 0.0,
                    peaks: vec![
                        TracePeak::new(0, 0, region_start - offset, first_mz, 0.0),
                        TracePeak::new(0, 0, region_end + offset, first_mz, 0.0),
                    ],
                    ..MassTrace::default()
                };
                traces.push(trace);
            }
            self.fit_and_validate(&mut fitter, &traces, feat, region, &mut report)?;
        }

        // "check if fit worked for at least one feature"
        let has_valid = staged
            .features
            .iter()
            .any(|f| status_of(f) == Some(STATUS_VALID));
        if !has_valid {
            staged.features.clear();
            report.cleared = true;
            *features = staged;
            return Ok(report);
        }

        if self.width > 0.0 {
            Self::reject_outliers(
                &mut staged,
                "model_width",
                self.width,
                "5 (width too large)",
                "6 (width too small)",
            )?;
        }
        let asym_limit = if self.asymmetric { self.asymmetry } else { 0.0 };
        if asym_limit > 0.0 {
            Self::reject_outliers(
                &mut staged,
                "model_asymmetry",
                asym_limit,
                "7 (asymmetry too high)",
                "8 (asymmetry too low)",
            )?;
        }

        let impute = !self.no_imputation;
        let mut quant_values = Vec::new();
        let mut failed = Vec::new();
        for (index, feat) in staged.features.iter_mut().enumerate() {
            let raw_intensity = f64::from(feat.base.intensity);
            feat.base
                .metadata
                .insert("raw_intensity".into(), float_meta(raw_intensity));
            let valid = status_of(feat).is_some_and(|s| s.starts_with('0'));
            if !valid {
                if impute {
                    failed.push(index);
                } else {
                    feat.base.intensity = 0.0;
                }
                report.failures += 1;
            } else {
                let area = meta_f64(feat, "model_area")?;
                if impute {
                    quant_values.push(DataPoint::new(raw_intensity.ln(), area.ln()));
                }
                feat.base.intensity = area as f32;
                report.successes += 1;
            }
        }
        if impute {
            let model = LinearModel::fit(&quant_values, LinearOptions::default())?;
            let c = model.coefficients();
            report.imputation = Some((c.slope, c.intercept));
            for index in failed {
                let feat = &mut staged.features[index];
                let x = f64::from(feat.base.intensity).ln();
                feat.base.intensity = (c.slope * x + c.intercept).exp() as f32;
            }
        }
        *features = staged;
        Ok(report)
    }

    /// Reject valid models whose `key` has a modified z-score `(v - median) /
    /// (1.4826 * MAD)` beyond `limit` on either side.
    fn reject_outliers(
        features: &mut FeatureMap,
        key: &str,
        limit: f64,
        too_large: &str,
        too_small: &str,
    ) -> Result<()> {
        let mut values = Vec::new();
        for f in &features.features {
            if status_of(f) == Some(STATUS_VALID) {
                values.push(meta_f64(f, key)?);
            }
        }
        let center = median(&mut values.clone())?;
        let mut abs_diffs: Vec<f64> = values.iter().map(|v| (v - center).abs()).collect();
        // "median absolute deviation (constant factor to approximate std. dev.)"
        let mad = 1.4826 * median(&mut abs_diffs)?;
        for f in &mut features.features {
            if status_of(f) == Some(STATUS_VALID) {
                let z = (meta_f64(f, key)? - center) / mad;
                if z > limit {
                    f.base
                        .metadata
                        .insert("model_status".into(), MetaValue::from(too_large));
                } else if z < -limit {
                    f.base
                        .metadata
                        .insert("model_status".into(), MetaValue::from(too_small));
                }
            }
        }
        Ok(())
    }
}

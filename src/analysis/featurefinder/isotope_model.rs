// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Widened averagine isotope patterns approximated by linear interpolation.
//!
//! Ports `FEATUREFINDER/IsotopeModel.h` and
//! `FEATUREFINDER/ExtendedIsotopeModel.h` with their `.cpp` files. See
//! `docs/FEATUREFINDER_MODELS_SUPPORT.md`.
//!
//! Both models build an averagine formula from a mass, take its coarse isotope
//! distribution, trim and renormalize it, "stretch" it onto the interpolation
//! grid (isotope `i` moves to grid index `i * distance / step / charge`,
//! rounded), and fold a peak shape - a Gaussian, or for [`IsotopeModel`] also a
//! Lorentzian - onto every isotope. They model raw (profile) data; the bare
//! distribution is `CoarseIsotopePatternGenerator`'s.
//!
//! The isotope distribution is computed in the source's binary32 `Peak1D`
//! arithmetic ([`ProbabilityPrecision::SourceSingle`]), because every weight
//! the models use is a `float` intensity there.

use crate::Result;
use crate::chemistry::EmpiricalFormula;
use crate::chemistry::isotopes::{
    CoarseIsotopePatternGenerator, CoarseMassMode, IsotopeDistribution, ProbabilityPrecision,
};
use crate::error::Error;
use crate::math::basic_statistics::BasicStatistics;
use crate::param::{Param, ParamValue};
use statrs::distribution::{Cauchy, Continuous};

use super::base_model::{
    InterpolationCore, InterpolationModel, MAX_SAMPLES, complete_parameters, float_value,
    interpolation_base_model, set_default, set_valid,
};

fn bad(message: String) -> Error {
    Error::InvalidValue(message)
}

/// The averagine atoms per dalton, in the source's `Averagines` enum order
/// `C, H, N, O, S`: the `float` literals of the source defaults.
const AVERAGINE_DEFAULTS: [f32; 5] = [
    0.044_439_89,
    0.069_815_72,
    0.012_217_73,
    0.013_293_99,
    0.000_375_25,
];
/// The element symbols in `Averagines` order.
const AVERAGINE_KEYS: [&str; 5] = ["C", "H", "N", "O", "S"];

/// The averagine formula of `mass` daltons: each element count is
/// `Int(0.5 + mass * atoms_per_dalton)`, and elements with a zero count are
/// left out, in the order `C H N O S` (source `IsotopeModel::getFormula` and
/// the start of `ExtendedIsotopeModel::setSamples`).
fn averagine_formula(mass: f64, averagine: &[f64; 5]) -> Result<EmpiricalFormula> {
    let mut form = String::new();
    // The source computes C, N, O, H, S and appends C, H, N, O, S.
    for (symbol, per_dalton) in AVERAGINE_KEYS.iter().zip(averagine.iter()) {
        let count = 0.5 + mass * per_dalton;
        // Source `Int(...)`: truncation, undefined beyond the 32-bit range.
        if !(count > f64::from(i32::MIN) - 1.0 && count < f64::from(i32::MAX) + 1.0) {
            return Err(bad(format!(
                "the averagine {symbol} count of mass {mass} does not fit a 32-bit integer"
            )));
        }
        let count = count as i32;
        if count != 0 {
            form.push_str(symbol);
            form.push_str(&count.to_string());
        }
    }
    EmpiricalFormula::parse(&form)
}

/// The coarse isotope distribution of `formula` with at most `max_isotope`
/// peaks (0 for all), trimmed on the right at `cutoff` and renormalized, in
/// the source's binary32 arithmetic.
fn trimmed_distribution(
    formula: &EmpiricalFormula,
    max_isotope: i32,
    mass_mode: CoarseMassMode,
    cutoff: f64,
) -> Result<IsotopeDistribution> {
    // Source `CoarseIsotopePatternGenerator(Size max_isotope, ...)`: 0 means all.
    let max_peaks = match max_isotope {
        0 => None,
        m if m > 0 => Some(m as usize),
        m => {
            return Err(bad(format!(
                "isotope:maximum {m} is negative; the source converts it to a huge unsigned count"
            )));
        }
    };
    let generator = CoarseIsotopePatternGenerator::new(max_peaks, mass_mode)?
        .with_precision(ProbabilityPrecision::SourceSingle);
    let mut distribution = generator.run(formula)?;
    distribution.trim_right(cutoff)?;
    distribution.renormalize_with(ProbabilityPrecision::SourceSingle)?;
    Ok(distribution)
}

/// Move isotope `i` of `isotopes` to grid index
/// `size_t(i * distance / step / charge + 0.5)` of a vector of
/// `size_t((n - 1) * distance / step + 1.6)` entries, zeroing its old slot,
/// from the last isotope down to isotope 1, as the source's stretch loop.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when there is no isotope, the stretched
/// length is not finite or exceeds [`MAX_SAMPLES`], or an index falls outside
/// it; the source's `size_t` conversion and indexing are undefined there.
fn stretch(isotopes: &[f64], distance: f64, step: f64, charge: f64) -> Result<Vec<f64>> {
    let n = isotopes.len();
    if n == 0 {
        return Err(bad(
            "the trimmed isotope distribution is empty; the source's size - 1 wraps around".into(),
        ));
    }
    let length = (n - 1) as f64 * distance / step + 1.6;
    if !(length >= 0.0 && length < MAX_SAMPLES as f64) {
        return Err(bad(format!(
            "the stretched isotope pattern needs {length} grid entries (at most {MAX_SAMPLES})"
        )));
    }
    let length = length as usize;
    if length < n {
        return Err(bad(format!(
            "the interpolation step {step} exceeds the isotope distance {distance}; the source \
             reads beyond its shrunken pattern"
        )));
    }
    let mut exact = vec![0.0; length];
    exact[..n].copy_from_slice(isotopes);
    for i in (1..n).rev() {
        let index = i as f64 * distance / step / charge + 0.5;
        if !(index >= 0.0 && index < length as f64) {
            return Err(bad(format!(
                "isotope {i} maps to grid index {index}, outside the pattern of {length}"
            )));
        }
        let value = exact[i];
        exact[index as usize] = value;
        exact[i] = 0.0;
    }
    Ok(exact)
}

/// The grid `-width, -width + step, ...` up to `width`, accumulated by
/// repeated addition as the source's `coord += interpolation_step_` loop.
fn shape_grid(width: f64, step: f64) -> Result<Vec<f64>> {
    let mut coordinates = Vec::new();
    let mut coord = -width;
    // A NaN width gives no point, as the source's loop condition.
    if coord.is_nan() || width.is_nan() || coord > width {
        return Ok(coordinates);
    }
    let estimate = 2.0 * width / step;
    if !(step > 0.0 && estimate.is_finite() && estimate < MAX_SAMPLES as f64) {
        return Err(bad(format!(
            "the peak shape of half width {width} at step {step} needs more than {MAX_SAMPLES} \
             points or never ends"
        )));
    }
    while coord <= width {
        if coordinates.len() >= MAX_SAMPLES {
            return Err(bad(format!(
                "the peak shape needs more than {MAX_SAMPLES} points"
            )));
        }
        coordinates.push(coord);
        coord += step;
    }
    Ok(coordinates)
}

/// A normal peak shape of standard deviation `stdev` over `±4 stdev`, sampled
/// at `step` and normalized to sum 1 (`BasicStatistics::normalApproximation`
/// with sum 1 and mean 0). Returns the half width and the shape.
fn gaussian_shape(stdev: f64, step: f64) -> Result<(f64, Vec<f64>)> {
    // "MAGIC alert, num stdev for smooth table for normal distribution".
    let width = stdev * 4.0;
    let coordinates = shape_grid(width, step)?;
    let mut normal = BasicStatistics::new();
    normal.set_sum(1.0);
    normal.set_mean(0.0);
    normal.set_variance(stdev * stdev);
    Ok((width, normal.normal_approximation_at(&coordinates)?))
}

/// Fold `shape` onto every non-zero entry of `pattern`, into
/// `min(pattern.len() + shape.len() - 1, trunc(2 * width / step * max_isotope + 1))`
/// entries, looping backwards "because then the small products tend to come
/// first".
fn fold(
    pattern: &[f64],
    shape: &[f64],
    width: f64,
    step: f64,
    max_isotope: i32,
) -> Result<Vec<f64>> {
    let full = pattern.len() as f64 + shape.len() as f64 - 1.0;
    let limit = (2.0 * width / step * f64::from(max_isotope) + 1.0).trunc();
    let r_max = full.min(limit);
    if !(r_max >= 0.0 && r_max < MAX_SAMPLES as f64) {
        return Err(bad(format!(
            "the folded isotope pattern needs {r_max} entries (0 to {MAX_SAMPLES})"
        )));
    }
    let r_max = r_max as usize;
    let mut result = vec![0.0; r_max];
    for i in (0..pattern.len()).rev() {
        if pattern[i] == 0.0 {
            continue;
        }
        let j_end = r_max.saturating_sub(i).min(shape.len());
        for j in (0..j_end).rev() {
            result[i + j] += pattern[i] * shape[j];
        }
    }
    Ok(result)
}

/// The peak shape [`IsotopeModel`] widens each isotope with (parameter
/// `isotope:mode:mode`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IsotopePeakShape {
    /// `Gaussian` (the default): a normal distribution of standard deviation
    /// `isotope:mode:GaussianSD`, sampled over four standard deviations on
    /// either side.
    Gaussian,
    /// `Lorentzian`: a Cauchy distribution of full width at half maximum
    /// `isotope:mode:LorentzFWHM`, sampled over eight FWHM on either side.
    Lorentzian,
}

impl IsotopePeakShape {
    fn as_str(self) -> &'static str {
        match self {
            Self::Gaussian => "Gaussian",
            Self::Lorentzian => "Lorentzian",
        }
    }
}

/// Averagine isotope distribution, widened by a Gaussian or Lorentzian peak
/// shape and approximated using linear interpolation (source `IsotopeModel`).
///
/// The model can sample raw peaks, depending on the positions queried. For the
/// distribution without widening use `EmpiricalFormula`'s isotope distribution
/// for a sum formula, or `CoarseIsotopePatternGenerator`'s peptide-weight
/// estimate for averagine.
///
/// Unlike the other models, changing the parameters does not resample: call
/// [`IsotopeModel::set_samples_from_formula`], usually with
/// [`IsotopeModel::formula`].
#[derive(Clone, Debug, PartialEq)]
pub struct IsotopeModel {
    core: InterpolationCore,
    /// Parameter `isotope:mode:GaussianSD`, default 0.1.
    isotope_stdev: f64,
    /// Parameter `isotope:mode:LorentzFWHM`, default 0.3.
    isotope_lorentz_fwhm: f64,
    /// Parameter `isotope:mode:mode`, default `Gaussian`.
    mode: IsotopePeakShape,
    /// Parameter `charge`, default 1 (the source member starts at 0 until the
    /// constructor applies the defaults).
    charge: u32,
    /// Parameter `statistics:mean`, default 0: the centroid m/z.
    mean: f64,
    /// The monoisotopic m/z of the last sampling (source `monoisotopic_mz_`,
    /// initially 0).
    monoisotopic_mz: f64,
    /// Parameters `averagines:C` ... `averagines:S`.
    averagine: [f64; 5],
    /// Parameter `isotope:maximum`, default 100.
    max_isotope: i32,
    /// Parameter `isotope:trim_right_cutoff`, default 0.001.
    trim_right_cutoff: f64,
    /// Parameter `isotope:distance`, default 1.000495.
    isotope_distance: f64,
    /// The distribution of the last sampling, without widening.
    isotope_distribution: IsotopeDistribution,
}

impl IsotopeModel {
    /// The handler name (source `setName("IsotopeModel")`).
    pub const NAME: &'static str = "IsotopeModel";

    fn default_state() -> Self {
        Self {
            core: InterpolationCore::default(),
            isotope_stdev: 0.1,
            isotope_lorentz_fwhm: 0.3,
            mode: IsotopePeakShape::Gaussian,
            charge: 1,
            mean: 0.0,
            monoisotopic_mz: 0.0,
            averagine: AVERAGINE_DEFAULTS.map(f64::from),
            max_isotope: 100,
            trim_right_cutoff: 0.001,
            isotope_distance: 1.000_495,
            isotope_distribution: IsotopeDistribution::default(),
        }
    }

    /// A model with the default parameters and no samples (source default
    /// constructor; unlike the other models its `updateMembers_` does not
    /// sample).
    ///
    /// # Errors
    ///
    /// Infallible in practice; `Result` for uniformity with the other models.
    pub fn new() -> Result<Self> {
        Ok(Self::default_state())
    }

    /// The default parameters (source `getDefaults`): the three shared ones,
    /// then the advanced `averagines:C` 0.04443989, `averagines:H`
    /// 0.06981572, `averagines:N` 0.01221773, `averagines:O` 0.01329399 and
    /// `averagines:S` 0.00037525 (each the `f64` of the source's `float`
    /// literal), `isotope:trim_right_cutoff` 0.001, `isotope:maximum` 100,
    /// `isotope:distance` 1.000495, `isotope:mode:mode` `Gaussian` (or
    /// `Lorentzian`), `isotope:mode:LorentzFWHM` 0.3,
    /// `isotope:mode:GaussianSD` 0.1, `charge` 1 and `statistics:mean` 0.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default_state().to_param_impl()
    }

    /// A default model after `set_parameters` with `param`.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param(param: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(param)?.0)
    }

    /// As [`IsotopeModel::from_param`], also returning the unknown-parameter
    /// warnings.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param_with_warnings(param: &Param) -> Result<(Self, Vec<String>)> {
        let mut model = Self::new()?;
        let warnings = model.set_parameters_impl(param)?;
        Ok((model, warnings))
    }

    /// The charge (source `getCharge`).
    pub fn charge(&self) -> u32 {
        self.charge
    }

    /// The interpolation's offset (source `getOffset`).
    pub fn offset(&self) -> f64 {
        self.core.interpolation.offset()
    }

    /// The peak shape.
    pub fn peak_shape(&self) -> IsotopePeakShape {
        self.mode
    }

    /// The isotope distribution, without widening, of the last
    /// [`IsotopeModel::set_samples_from_formula`]: useful for the number of
    /// isotopes the model contains and their positions (source
    /// `getIsotopeDistribution`). Before any sampling it is the source's
    /// default single peak `(0, 1)`.
    pub fn isotope_distribution(&self) -> &IsotopeDistribution {
        &self.isotope_distribution
    }

    /// The averagine peptide formula of mass `statistics:mean * charge`
    /// (source `getFormula`); set the mean and charge through the parameters.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when an element count does not fit the
    /// source's 32-bit integer, or the formula does not parse (a negative
    /// mass gives negative counts).
    pub fn formula(&self) -> Result<EmpiricalFormula> {
        averagine_formula(self.mean * f64::from(self.charge), &self.averagine)
    }

    /// Sample the widened isotope distribution of `formula` (source
    /// `setSamples(const EmpiricalFormula&)`).
    ///
    /// The distribution is trimmed at `isotope:trim_right_cutoff` and
    /// renormalized; its mean isotope offset, `sum(p_i * i) * distance /
    /// charge` (each product in `float`, as the source's `float * Int`), is
    /// subtracted from `statistics:mean` to give the monoisotopic m/z, which
    /// the interpolation maps the pattern's first isotope onto. The result is
    /// scaled so that `step * sum` equals the scaling factor.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a zero charge, a negative
    /// `isotope:maximum`, a peak shape of zero width, a step larger than the
    /// isotope distance, a pattern beyond [`MAX_SAMPLES`] and a distribution
    /// the generator refuses; the source divides by zero, converts to `size_t`
    /// out of range or reads past its pattern in those cases (a zero-width
    /// Lorentzian throws `std::domain_error` from Boost). The model is then
    /// unchanged.
    pub fn set_samples_from_formula(&mut self, formula: &EmpiricalFormula) -> Result<()> {
        if self.charge == 0 {
            return Err(bad(
                "an isotope model needs a non-zero charge; the source divides by it".into(),
            ));
        }
        let charge = f64::from(self.charge);
        let step = self.core.interpolation_step;
        let distribution = trimmed_distribution(
            formula,
            self.max_isotope,
            CoarseMassMode::Nominal,
            self.trim_right_cutoff,
        )?;
        let mut isotopes = Vec::with_capacity(distribution.len());
        let mut isotopes_mean = 0.0f64;
        for (cnt, peak) in distribution.peaks().iter().enumerate() {
            isotopes.push(peak.probability);
            // Source `peak.getIntensity() * cnt`: float times Int, in float.
            isotopes_mean += f64::from(peak.probability as f32 * cnt as f32);
        }
        isotopes_mean *= self.isotope_distance / charge;
        let exact = stretch(&isotopes, self.isotope_distance, step, charge)?;
        let (width, shape) = match self.mode {
            IsotopePeakShape::Gaussian => gaussian_shape(self.isotope_stdev, step)?,
            IsotopePeakShape::Lorentzian => {
                // "Lorentzian has infinite support, but we need to stop
                // sampling at some point: 8*FWHM". Cauchy uses the HWHM.
                let width = self.isotope_lorentz_fwhm * 8.0;
                let cauchy = Cauchy::new(0.0, self.isotope_lorentz_fwhm / 2.0).map_err(|e| {
                    bad(format!(
                        "isotope:mode:LorentzFWHM {} gives no Cauchy distribution: {e}",
                        self.isotope_lorentz_fwhm
                    ))
                })?;
                let shape = shape_grid(width, step)?
                    .into_iter()
                    .map(|coord| cauchy.pdf(coord))
                    .collect();
                (width, shape)
            }
        };
        let mut result = fold(&exact, &shape, width, step, self.max_isotope)?;
        let monoisotopic_mz = self.mean - isotopes_mean;
        let sum: f64 = result.iter().fold(0.0, |acc, &v| acc + v);
        let factor = self.core.scaling / (step * sum);
        for value in &mut result {
            *value *= factor;
        }
        self.isotope_distribution = distribution;
        self.monoisotopic_mz = monoisotopic_mz;
        self.core.interpolation.set_data(result);
        self.core
            .interpolation
            .set_mapping(step, width / step, monoisotopic_mz);
        Ok(())
    }

    fn to_param_impl(&self) -> Result<Param> {
        let mut p = Param::new();
        self.core.write_param(&mut p)?;
        for (i, symbol) in AVERAGINE_KEYS.iter().enumerate() {
            set_default(
                &mut p,
                &format!("averagines:{symbol}"),
                ParamValue::Float(self.averagine[i]),
                &format!("Number of {symbol} atoms per Dalton of mass."),
                true,
            )?;
        }
        set_default(
            &mut p,
            "isotope:trim_right_cutoff",
            ParamValue::Float(self.trim_right_cutoff),
            "Cutoff in averagine distribution, trailing isotopes below this relative intensity are not considered.",
            true,
        )?;
        set_default(
            &mut p,
            "isotope:maximum",
            ParamValue::Integer(i64::from(self.max_isotope)),
            "Maximum isotopic rank to be considered.",
            true,
        )?;
        set_default(
            &mut p,
            "isotope:distance",
            ParamValue::Float(self.isotope_distance),
            "Distance between consecutive isotopic peaks.",
            true,
        )?;
        set_default(
            &mut p,
            "isotope:mode:mode",
            ParamValue::String(self.mode.as_str().to_owned()),
            "Peak Shape used around each isotope peak.",
            true,
        )?;
        set_valid(&mut p, "isotope:mode:mode", &["Gaussian", "Lorentzian"])?;
        set_default(
            &mut p,
            "isotope:mode:LorentzFWHM",
            ParamValue::Float(self.isotope_lorentz_fwhm),
            "Full width of the Lorentzian (Cauchy) function applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.",
            true,
        )?;
        set_default(
            &mut p,
            "isotope:mode:GaussianSD",
            ParamValue::Float(self.isotope_stdev),
            "Standard deviation of Gaussian applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.",
            true,
        )?;
        set_default(
            &mut p,
            "charge",
            ParamValue::Integer(i64::from(self.charge)),
            "Charge state of the model.",
            true,
        )?;
        set_default(
            &mut p,
            "statistics:mean",
            ParamValue::Float(self.mean),
            "Centroid m/z (as opposed to monoisotopic m/z).",
            true,
        )?;
        Ok(p)
    }

    fn set_parameters_impl(&mut self, param: &Param) -> Result<Vec<String>> {
        let (complete, warnings) = complete_parameters(Self::NAME, Self::defaults()?, param)?;
        let mut staged = self.clone();
        staged.core.read_param(&complete)?;
        staged.charge = complete.value("charge")?.to_u32()?;
        staged.isotope_stdev = float_value(&complete, "isotope:mode:GaussianSD")?;
        staged.isotope_lorentz_fwhm = float_value(&complete, "isotope:mode:LorentzFWHM")?;
        staged.mean = float_value(&complete, "statistics:mean")?;
        staged.max_isotope = complete.value("isotope:maximum")?.to_i32()?;
        staged.trim_right_cutoff = float_value(&complete, "isotope:trim_right_cutoff")?;
        staged.isotope_distance = float_value(&complete, "isotope:distance")?;
        staged.mode = match complete.value("isotope:mode:mode")?.as_str()? {
            "Lorentzian" => IsotopePeakShape::Lorentzian,
            _ => IsotopePeakShape::Gaussian,
        };
        for (i, symbol) in AVERAGINE_KEYS.iter().enumerate() {
            staged.averagine[i] = float_value(&complete, &format!("averagines:{symbol}"))?;
        }
        *self = staged;
        Ok(warnings)
    }
}

interpolation_base_model!(IsotopeModel);

impl InterpolationModel for IsotopeModel {
    fn core(&self) -> &InterpolationCore {
        &self.core
    }
    fn core_mut(&mut self) -> &mut InterpolationCore {
        &mut self.core
    }

    /// The monoisotopic m/z of the last sampling (0 before any): an m/z, not
    /// necessarily the monoisotopic mass.
    fn center(&self) -> Result<f64> {
        Ok(self.monoisotopic_mz)
    }

    /// Shift the model without recomputing it: `statistics:mean` and the
    /// monoisotopic m/z move by `offset - offset()`.
    ///
    /// The source notes that this leaves a discrepancy, minor for shifts of
    /// one or two standard deviations but significant beyond; resample after
    /// [`crate::analysis::featurefinder::base_model::BaseModel::set_parameters`]
    /// for large shifts.
    fn set_offset(&mut self, offset: f64) {
        let diff = offset - self.core.interpolation.offset();
        self.mean += diff;
        self.monoisotopic_mz += diff;
        self.core.interpolation.set_offset(offset);
    }

    /// The source inherits the base class's `setSamples()`, which throws
    /// `Exception::NotImplemented`; this returns [`Error::Unsupported`]. Use
    /// [`IsotopeModel::set_samples_from_formula`].
    fn set_samples(&mut self) -> Result<()> {
        Err(Error::Unsupported(
            "IsotopeModel samples from a formula: use set_samples_from_formula (the source's \
             setSamples() throws NotImplemented)"
                .into(),
        ))
    }
}

/// Averagine isotope distribution of a given monoisotopic m/z, widened by a
/// Gaussian and approximated using linear interpolation (source
/// `ExtendedIsotopeModel`).
///
/// Unlike [`IsotopeModel`] it is placed by its monoisotopic m/z
/// (`isotope:monoisotopic_mz`), always uses a Gaussian, and resamples whenever
/// its parameters change.
#[derive(Clone, Debug, PartialEq)]
pub struct ExtendedIsotopeModel {
    core: InterpolationCore,
    /// Parameter `isotope:stdev`, default 0.1.
    isotope_stdev: f64,
    /// Parameter `charge`, default 1.
    charge: u32,
    /// Parameter `isotope:monoisotopic_mz`, default 1.
    monoisotopic_mz: f64,
    /// Parameters `averagines:C` ... `averagines:S`.
    averagine: [f64; 5],
    /// Parameter `isotope:maximum`, default 100.
    max_isotope: i32,
    /// Parameter `isotope:trim_right_cutoff`, default 0.001.
    trim_right_cutoff: f64,
    /// Parameter `isotope:distance`, default 1.000495.
    isotope_distance: f64,
}

impl ExtendedIsotopeModel {
    /// The handler name (source `setName("ExtendedIsotopeModel")`).
    pub const NAME: &'static str = "ExtendedIsotopeModel";

    fn default_state() -> Self {
        Self {
            core: InterpolationCore::default(),
            isotope_stdev: 0.1,
            charge: 1,
            monoisotopic_mz: 1.0,
            averagine: AVERAGINE_DEFAULTS.map(f64::from),
            max_isotope: 100,
            trim_right_cutoff: 0.001,
            isotope_distance: 1.000_495,
        }
    }

    /// A model with the default parameters, sampled (source default
    /// constructor). The defaults describe a 1 Da averagine, whose formula is
    /// empty, so the pattern is a single Gaussian.
    ///
    /// # Errors
    ///
    /// As [`InterpolationModel::set_samples`]; the defaults succeed.
    pub fn new() -> Result<Self> {
        let mut model = Self::default_state();
        model.set_samples()?;
        Ok(model)
    }

    /// The default parameters (source `getDefaults`): the three shared ones,
    /// then the advanced averagines as for [`IsotopeModel::defaults`],
    /// `isotope:trim_right_cutoff` 0.001, `isotope:maximum` 100,
    /// `isotope:distance` 1.000495, `isotope:stdev` 0.1, `charge` 1 and
    /// `isotope:monoisotopic_mz` 1.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default_state().to_param_impl()
    }

    /// A default model after `set_parameters` with `param`.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param(param: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(param)?.0)
    }

    /// As [`ExtendedIsotopeModel::from_param`], also returning the
    /// unknown-parameter warnings.
    ///
    /// # Errors
    ///
    /// As `set_parameters`.
    pub fn from_param_with_warnings(param: &Param) -> Result<(Self, Vec<String>)> {
        let mut model = Self::new()?;
        let warnings = model.set_parameters_impl(param)?;
        Ok((model, warnings))
    }

    /// The charge (source `getCharge`).
    pub fn charge(&self) -> u32 {
        self.charge
    }

    /// The interpolation's offset (source `getOffset`).
    pub fn offset(&self) -> f64 {
        self.core.interpolation.offset()
    }

    fn to_param_impl(&self) -> Result<Param> {
        let mut p = Param::new();
        self.core.write_param(&mut p)?;
        for (i, symbol) in AVERAGINE_KEYS.iter().enumerate() {
            set_default(
                &mut p,
                &format!("averagines:{symbol}"),
                ParamValue::Float(self.averagine[i]),
                &format!("Number of {symbol} atoms per Dalton of mass."),
                true,
            )?;
        }
        set_default(
            &mut p,
            "isotope:trim_right_cutoff",
            ParamValue::Float(self.trim_right_cutoff),
            "Cutoff in averagine distribution, trailing isotopes below this relative intensity are not considered.",
            true,
        )?;
        set_default(
            &mut p,
            "isotope:maximum",
            ParamValue::Integer(i64::from(self.max_isotope)),
            "Maximum isotopic rank to be considered.",
            true,
        )?;
        set_default(
            &mut p,
            "isotope:distance",
            ParamValue::Float(self.isotope_distance),
            "Distance between consecutive isotopic peaks.",
            true,
        )?;
        set_default(
            &mut p,
            "isotope:stdev",
            ParamValue::Float(self.isotope_stdev),
            "Standard deviation of gaussian applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.",
            true,
        )?;
        set_default(
            &mut p,
            "charge",
            ParamValue::Integer(i64::from(self.charge)),
            "Charge state of the model.",
            true,
        )?;
        set_default(
            &mut p,
            "isotope:monoisotopic_mz",
            ParamValue::Float(self.monoisotopic_mz),
            "Monoisotopic m/z of the model.",
            true,
        )?;
        Ok(p)
    }

    fn set_parameters_impl(&mut self, param: &Param) -> Result<Vec<String>> {
        let (complete, warnings) = complete_parameters(Self::NAME, Self::defaults()?, param)?;
        let mut staged = self.clone();
        staged.core.read_param(&complete)?;
        staged.charge = complete.value("charge")?.to_u32()?;
        staged.isotope_stdev = float_value(&complete, "isotope:stdev")?;
        staged.monoisotopic_mz = float_value(&complete, "isotope:monoisotopic_mz")?;
        staged.max_isotope = complete.value("isotope:maximum")?.to_i32()?;
        staged.trim_right_cutoff = float_value(&complete, "isotope:trim_right_cutoff")?;
        staged.isotope_distance = float_value(&complete, "isotope:distance")?;
        for (i, symbol) in AVERAGINE_KEYS.iter().enumerate() {
            staged.averagine[i] = float_value(&complete, &format!("averagines:{symbol}"))?;
        }
        staged.set_samples()?;
        *self = staged;
        Ok(warnings)
    }
}

interpolation_base_model!(ExtendedIsotopeModel);

impl InterpolationModel for ExtendedIsotopeModel {
    fn core(&self) -> &InterpolationCore {
        &self.core
    }
    fn core_mut(&mut self) -> &mut InterpolationCore {
        &mut self.core
    }

    /// The monoisotopic m/z.
    fn center(&self) -> Result<f64> {
        Ok(self.monoisotopic_mz)
    }

    /// Shift the model without recomputing it: the monoisotopic m/z moves by
    /// `offset - offset()`. The source's note on [`IsotopeModel`]'s shift
    /// applies here too.
    fn set_offset(&mut self, offset: f64) {
        let diff = offset - self.core.interpolation.offset();
        self.monoisotopic_mz += diff;
        self.core.interpolation.set_offset(offset);
    }

    /// Sample the averagine pattern of mass `monoisotopic_mz * charge`,
    /// widened by a Gaussian of `isotope:stdev` over four standard deviations,
    /// with the pattern's first isotope mapped onto the monoisotopic m/z, and
    /// scaled so that the rectangle integral equals the scaling factor.
    ///
    /// # Errors
    ///
    /// As [`IsotopeModel::set_samples_from_formula`], and an averagine count
    /// beyond 32 bits. The model is then unchanged.
    fn set_samples(&mut self) -> Result<()> {
        if self.charge == 0 {
            return Err(bad(
                "an isotope model needs a non-zero charge; the source divides by it".into(),
            ));
        }
        let charge = f64::from(self.charge);
        let step = self.core.interpolation_step;
        let formula = averagine_formula(self.monoisotopic_mz * charge, &self.averagine)?;
        let distribution = trimmed_distribution(
            &formula,
            self.max_isotope,
            CoarseMassMode::Approximate,
            self.trim_right_cutoff,
        )?;
        let isotopes: Vec<f64> = distribution.peaks().iter().map(|p| p.probability).collect();
        let exact = stretch(&isotopes, self.isotope_distance, step, charge)?;
        let (width, shape) = gaussian_shape(self.isotope_stdev, step)?;
        let mut result = fold(&exact, &shape, width, step, self.max_isotope)?;
        let sum: f64 = result.iter().fold(0.0, |acc, &v| acc + v);
        let factor = self.core.scaling / step / sum;
        for value in &mut result {
            *value *= factor;
        }
        self.core.interpolation.set_data(result);
        self.core
            .interpolation
            .set_mapping(step, width / step, self.monoisotopic_mz);
        Ok(())
    }
}

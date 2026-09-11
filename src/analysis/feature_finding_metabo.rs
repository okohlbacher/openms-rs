// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Source mass-trace isotope assembly with atomic owned feature publication.

mod predictor;
use super::feature_hypothesis::{FeatureHypothesis, FeatureHypothesisLimits};
use crate::{
    Error, Result,
    chemistry::{AveragineComposition, CoarseIsotopePatternGenerator},
    concept::{
        progress_logger::{ProgressLogType, ProgressLogger},
        unique_id::UniqueIdGenerator,
    },
    kernel::{Feature, FeatureMap, MSChromatogram, MassTrace, MassTraceQuantMethod},
    metadata::MetaValue,
};
use std::{collections::BTreeSet, mem::size_of};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum IsotopeFilteringModel {
    Metabolites2,
    #[default]
    Metabolites5,
    Peptides,
    None,
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum IsotopeRtOverlapReference {
    #[default]
    Longer,
    Shorter,
}

/// All nineteen source configuration fields. Finite negative RT/IM ranges and
/// inverted charge ranges retain their source branch behavior. A negative m/z
/// range is rejected when a positive charge consumes the source unsigned offset
/// conversion. chrom_fwhm is
/// stored by the source but never used by this algorithm.
#[derive(Clone, Debug, PartialEq)]
pub struct FeatureFindingMetaboOptions {
    pub local_rt_range: f64,
    pub local_im_range: f64,
    pub local_mz_range: f64,
    pub charge_lower_bound: usize,
    pub charge_upper_bound: usize,
    pub chrom_fwhm: f64,
    pub report_summed_intensities: bool,
    pub enable_rt_filtering: bool,
    pub isotope_rt_overlap_reference: IsotopeRtOverlapReference,
    pub min_isotope_rt_overlap: f64,
    pub isotope_filtering_model: IsotopeFilteringModel,
    pub mz_scoring_13c: bool,
    pub use_smoothed_intensities: bool,
    pub report_smoothed_intensities: bool,
    pub report_convex_hulls: bool,
    pub report_chromatograms: bool,
    pub remove_single_traces: bool,
    pub mz_scoring_by_elements: bool,
    pub elements: String,
}
impl Default for FeatureFindingMetaboOptions {
    fn default() -> Self {
        Self {
            local_rt_range: 10.0,
            local_im_range: 0.02,
            local_mz_range: 6.5,
            charge_lower_bound: 1,
            charge_upper_bound: 3,
            chrom_fwhm: 5.0,
            report_summed_intensities: false,
            enable_rt_filtering: true,
            isotope_rt_overlap_reference: IsotopeRtOverlapReference::Longer,
            min_isotope_rt_overlap: 0.7,
            isotope_filtering_model: IsotopeFilteringModel::Metabolites5,
            mz_scoring_13c: false,
            use_smoothed_intensities: true,
            report_smoothed_intensities: true,
            report_convex_hulls: false,
            report_chromatograms: false,
            remove_single_traces: false,
            mz_scoring_by_elements: false,
            elements: "CHNOPS".into(),
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FeatureFindingMetaboLimits {
    pub max_traces: usize,
    pub max_peaks: usize,
    pub max_hypotheses: usize,
    pub max_features: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for FeatureFindingMetaboLimits {
    fn default() -> Self {
        Self {
            max_traces: 1_000_000,
            max_peaks: 10_000_000,
            max_hypotheses: 1_000_000,
            max_features: 1_000_000,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}
/// Explicit equivalents of source configuration warnings.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FeatureFindingMetaboDiagnostics {
    pub smoothed_reporting_disabled: bool,
    pub isotope_filtering_disabled_by_elements: bool,
}
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FeatureFindingMetaboOutput {
    pub features: FeatureMap,
    /// Acceptance order, which can differ from the final feature-map m/z order.
    pub chromatograms: Vec<Vec<MSChromatogram>>,
    pub diagnostics: FeatureFindingMetaboDiagnostics,
}
#[derive(Clone)]
pub struct FeatureFindingMetabo {
    options: FeatureFindingMetaboOptions,
    isotope_masses: Vec<Vec<f64>>,
    effective_model: IsotopeFilteringModel,
    has_im_data: bool,
    pub limits: FeatureFindingMetaboLimits,
    pub logger: ProgressLogger,
}
impl FeatureFindingMetabo {
    pub fn new() -> Result<Self> {
        Self::with_options(FeatureFindingMetaboOptions::default())
    }
    pub fn with_options(options: FeatureFindingMetaboOptions) -> Result<Self> {
        let limits = FeatureFindingMetaboLimits::default();
        let mut work = Work::new(limits);
        let isotope_masses = prepare_options(&options, &mut work)?;
        let mut logger = ProgressLogger::new();
        logger.set_log_type(ProgressLogType::Cmd);
        Ok(Self {
            effective_model: options.isotope_filtering_model,
            options,
            isotope_masses,
            has_im_data: false,
            limits,
            logger,
        })
    }
    pub fn options(&self) -> &FeatureFindingMetaboOptions {
        &self.options
    }
    /// Checks/prepares a new alphabet before replacing any configuration state.
    pub fn set_options(&mut self, options: FeatureFindingMetaboOptions) -> Result<()> {
        let isotope_masses = prepare_options(&options, &mut Work::new(self.limits))?;
        self.effective_model = options.isotope_filtering_model;
        self.options = options;
        self.isotope_masses = isotope_masses;
        Ok(())
    }
    pub fn effective_isotope_filtering_model(&self) -> IsotopeFilteringModel {
        self.effective_model
    }
    pub fn contains_im_data(&self) -> bool {
        self.has_im_data
    }
    /// Borrowed input is sorted only after complete success; RNG, flags and owned
    /// output publish together. Progress side effects cannot be rolled back.
    pub fn run(
        &mut self,
        traces: &mut [MassTrace],
        rng: &mut UniqueIdGenerator,
    ) -> Result<FeatureFindingMetaboOutput> {
        let mut work = Work::new(self.limits);
        let effective = if self.options.mz_scoring_by_elements {
            IsotopeFilteringModel::None
        } else {
            self.effective_model
        };
        let diagnostics = FeatureFindingMetaboDiagnostics {
            smoothed_reporting_disabled: self.options.report_smoothed_intensities
                && !self.options.use_smoothed_intensities,
            isotope_filtering_disabled_by_elements: self.options.mz_scoring_by_elements
                && self.effective_model != IsotopeFilteringModel::None,
        };
        if traces.is_empty() {
            self.effective_model = effective;
            return Ok(FeatureFindingMetaboOutput {
                diagnostics,
                ..Default::default()
            });
        }
        work.inputs(traces)?;
        work.allocate::<UniqueIdGenerator>(1)?;
        work.scan(size_of::<UniqueIdGenerator>() / 8 + 1)?;
        let mut staged_rng = rng.clone();
        let mut order = work.vector(traces.len())?;
        order.extend(0..traces.len());
        work.sort(order.len())?;
        order.sort_unstable_by(|&a, &b| {
            traces[a]
                .centroid_mz()
                .partial_cmp(&traces[b].centroid_mz())
                .unwrap()
                .then(a.cmp(&b))
        });
        let has_im = traces.iter().any(MassTrace::contains_im_data);
        let count = i64::try_from(traces.len()).map_err(|_| resource())?;
        self.logger
            .start_progress(0, count, "assembling mass traces to features")?;
        let hypotheses = find_hypotheses(
            traces,
            &order,
            &self.options,
            &self.isotope_masses,
            effective,
            has_im,
            &mut self.logger,
            &mut work,
        );
        let end = self.logger.end_progress(0);
        let hypotheses = hypotheses?;
        end?;
        let mut output = accept(
            traces,
            hypotheses,
            &self.options,
            effective,
            has_im,
            &mut staged_rng,
            &mut work,
        )?;
        output.diagnostics = diagnostics;
        let mut inverse = work.vector(order.len())?;
        inverse.resize(order.len(), 0usize);
        work.scan(mul(order.len(), 3)?)?;
        for (new, &old) in order.iter().enumerate() {
            inverse[old] = new;
        }
        for i in 0..inverse.len() {
            while inverse[i] != i {
                let j = inverse[i];
                traces.swap(i, j);
                inverse.swap(i, j);
            }
        }
        *rng = staged_rng;
        self.effective_model = effective;
        self.has_im_data = has_im;
        Ok(output)
    }
    /// Returns previous output ownership without examining arbitrary metadata.
    pub fn run_into(
        &mut self,
        traces: &mut [MassTrace],
        rng: &mut UniqueIdGenerator,
        features: &mut FeatureMap,
        chromatograms: &mut Vec<Vec<MSChromatogram>>,
    ) -> Result<(
        FeatureMap,
        Vec<Vec<MSChromatogram>>,
        FeatureFindingMetaboDiagnostics,
    )> {
        let result = self.run(traces, rng)?;
        Ok((
            std::mem::replace(features, result.features),
            std::mem::replace(chromatograms, result.chromatograms),
            result.diagnostics,
        ))
    }
}
fn prepare_options(
    options: &FeatureFindingMetaboOptions,
    work: &mut Work,
) -> Result<Vec<Vec<f64>>> {
    if !options.min_isotope_rt_overlap.is_finite()
        || !(0.0..=1.0).contains(&options.min_isotope_rt_overlap)
    {
        return Err(bad("min_isotope_rt_overlap must be between zero and one"));
    }
    crate::chemistry::metabo_elements::isotope_masses(
        &options.elements,
        &mut work.remaining,
        &mut work.bytes,
    )
}
struct Hypothesis {
    members: Vec<usize>,
    score: f64,
    charge: i64,
}
impl Hypothesis {
    fn copy(&self, work: &mut Work) -> Result<Self> {
        let mut members = work.vector(self.members.len())?;
        members.extend_from_slice(&self.members);
        Ok(Self {
            members,
            score: self.score,
            charge: self.charge,
        })
    }
}
#[allow(clippy::too_many_arguments)]
fn find_hypotheses(
    traces: &[MassTrace],
    order: &[usize],
    options: &FeatureFindingMetaboOptions,
    alphabet: &[Vec<f64>],
    model: IsotopeFilteringModel,
    has_im: bool,
    logger: &mut ProgressLogger,
    work: &mut Work,
) -> Result<Vec<Hypothesis>> {
    let mut total = 0.0;
    let mut intensities = work.vector(traces.len())?;
    intensities.resize(traces.len(), 0.0);
    for &i in order {
        let value = work.intensity(&traces[i], options.use_smoothed_intensities)?;
        intensities[i] = value;
        total = finite(total + value)?;
    }
    let mut hypotheses = Vec::new();
    for (progress, &i) in order.iter().enumerate() {
        logger.set_progress(progress as i64)?;
        let mono = &traces[i];
        let mut candidates = Vec::new();
        work.push(&mut candidates, i)?;
        for &j in &order[progress + 1..] {
            work.scan(1)?;
            if finite((traces[j].centroid_mz() - mono.centroid_mz()).abs())?
                > finite(options.local_mz_range)?
            {
                break;
            }
            if finite((traces[j].centroid_rt() - mono.centroid_rt()).abs())?
                <= finite(options.local_rt_range)?
                && (!has_im
                    || finite((traces[j].centroid_im() - mono.centroid_im()).abs())?
                        <= finite(options.local_im_range)?)
            {
                work.push(&mut candidates, j)?;
            }
        }
        let mut members = work.vector(1)?;
        members.push(i);
        let singleton = Hypothesis {
            members,
            score: finite(intensities[i] / total)?,
            charge: 0,
        };
        let saved = singleton.copy(work)?;
        work.hypothesis(&mut hypotheses, saved)?;
        for charge in options.charge_lower_bound..=options.charge_upper_bound {
            work.scan(1)?;
            let mut current = singleton.copy(work)?;
            let count = finite((charge as f64 * finite(options.local_mz_range)?).floor())?;
            if count < 0.0 || count > i32::MAX as f64 {
                return Err(bad("isotope offset is outside the source signed32 domain"));
            }
            let mut last = 0;
            for offset in 1..=count as usize {
                work.scan(1)?;
                let window = isotope_window(alphabet, offset, work)?;
                let mut best = 0.0;
                let mut best_index = 0;
                for (index, &j) in candidates.iter().enumerate().skip(last + 1) {
                    work.scan(1)?;
                    let rt = score_rt(mono, &traces[j], options, work)?;
                    let mz = score_mz(mono, &traces[j], offset, charge, window, options)?;
                    let intensity = if model == IsotopeFilteringModel::Peptides {
                        let mut values = work.vector(add(current.members.len(), 1)?)?;
                        for &member in &current.members {
                            values.push(work.intensity(&traces[member], false)?);
                        }
                        values.push(intensities[j]);
                        averagine_score(
                            &values,
                            finite(traces[j].centroid_mz() * charge as f64)?,
                            work,
                        )?
                    } else {
                        1.0
                    };
                    let score = if rt > 0.0 && mz > 0.0 && intensity > 0.0 {
                        finite((rt.ln() + mz.ln() + intensity.ln()).exp())?
                    } else {
                        0.0
                    };
                    if score > best {
                        best = score;
                        best_index = index;
                    }
                }
                if best <= 0.0 {
                    break;
                }
                let selected = candidates[best_index];
                work.push(&mut current.members, selected)?;
                current.score =
                    finite(current.score + finite(intensities[selected] * best / total)?)?;
                current.charge =
                    i64::try_from(charge).map_err(|_| bad("hypothesis charge exceeds signed64"))?;
                last = best_index;
                let saved = current.copy(work)?;
                work.hypothesis(&mut hypotheses, saved)?;
            }
        }
    }
    work.sort(hypotheses.len())?;
    // Stable serial seed/charge/prefix order breaks otherwise unspecified ties.
    work.allocate::<Hypothesis>(hypotheses.len())?;
    hypotheses.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
    Ok(hypotheses)
}
fn isotope_window(
    alphabet: &[Vec<f64>],
    offset: usize,
    work: &mut Work,
) -> Result<Option<(f64, f64)>> {
    let mut bounds: Option<(f64, f64)> = None;
    work.scan(alphabet.len())?;
    for masses in alphabet {
        for &mass in masses.iter().skip(1) {
            work.scan(1)?;
            let gap = (mass.round() as i32) - (masses[0].round() as i32);
            if gap <= 0 {
                return Err(bad("isotope table has nonpositive integer mass gap"));
            }
            if gap as usize > offset {
                break;
            }
            let defect = finite((mass - masses[0] - gap as f64) * (offset / gap as usize) as f64)?;
            bounds = Some(bounds.map_or((defect, defect), |(low, high)| {
                (low.min(defect), high.max(defect))
            }));
        }
    }
    Ok(bounds.map(|(low, high)| (offset as f64 + low, offset as f64 + high)))
}
fn variance(sd: f64) -> Result<f64> {
    finite((2.0 * sd.ln()).exp())
}
fn score_mz(
    a: &MassTrace,
    b: &MassTrace,
    offset: usize,
    charge: usize,
    window: Option<(f64, f64)>,
    options: &FeatureFindingMetaboOptions,
) -> Result<f64> {
    let distance = finite((b.centroid_mz() - a.centroid_mz()).abs())?;
    let variances = finite(variance(a.centroid_sd())? + variance(b.centroid_sd())?)?;
    let charge = charge as f64;
    if options.mz_scoring_by_elements {
        let Some((low, high)) = window else {
            return Ok(0.0);
        };
        let low = finite(low / charge)?;
        let high = finite(high / charge)?;
        let sigma = variances.sqrt();
        let deviation = finite(3.0 * sigma)?;
        if distance > low && distance < high {
            return Ok(1.0);
        }
        if distance > low - deviation && distance < high + deviation {
            let exponent = finite(if distance < low {
                (low - distance) / sigma
            } else {
                (distance - high) / sigma
            })?;
            return finite((-0.5 * exponent * exponent).exp());
        }
        Ok(0.0)
    } else {
        let offset = offset as f64;
        let mean = finite(if options.mz_scoring_13c {
            crate::constants::C13C12_MASSDIFF_U * offset / charge
        } else {
            (1.000857 * offset + 0.001091) / charge
        })?;
        let sd = finite((0.0016633 * offset - 0.0004751) / charge)?;
        let sigma = finite((variance(sd)? + variances).sqrt())?;
        if distance < mean + 3.0 * sigma && distance > mean - 3.0 * sigma {
            let exponent = finite((distance - mean) / sigma)?;
            finite((-0.5 * exponent * exponent).exp())
        } else {
            Ok(0.0)
        }
    }
}
fn score_rt(
    a: &MassTrace,
    b: &MassTrace,
    options: &FeatureFindingMetaboOptions,
    work: &mut Work,
) -> Result<f64> {
    if !options.enable_rt_filtering {
        return Ok(1.0);
    }
    let mut values = Vec::new();
    for trace in [a, b] {
        let (left, right) = trace.fwhm_borders();
        let peaks = trace
            .peaks()
            .get(left..=right)
            .ok_or_else(|| bad("empty trace has no source FWHM endpoint"))?;
        work.scan(peaks.len())?;
        for peak in peaks {
            work.push(
                &mut values,
                (finite(peak.rt())?, finite(f64::from(peak.intensity))?),
            )?;
        }
    }
    work.sort(values.len())?;
    work.allocate::<(f64, f64)>(values.len())?;
    values.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut first = None;
    let mut last = 0.0;
    let mut index = 0;
    work.scan(values.len())?;
    while index < values.len() {
        let mut end = index + 1;
        while end < values.len() && values[end].0 == values[index].0 {
            end += 1;
        }
        if end - index == 2 {
            work.push(&mut x, values[index].1)?;
            work.push(&mut y, values[index + 1].1)?;
            first.get_or_insert(values[index].0);
            last = values[index].0;
        }
        index = end;
    }
    let overlap = finite(first.map_or(0.0, |start| (last - start).abs()))?;
    let denominator = if options.isotope_rt_overlap_reference == IsotopeRtOverlapReference::Longer {
        a.fwhm().max(b.fwhm())
    } else {
        a.fwhm().min(b.fwhm())
    };
    let proportion = if denominator > 0.0 {
        finite(overlap / denominator)?
    } else {
        0.0
    };
    if proportion < options.min_isotope_rt_overlap {
        return Ok(0.0);
    }
    if options.isotope_rt_overlap_reference == IsotopeRtOverlapReference::Shorter {
        let (long, short) = if a.fwhm() >= b.fwhm() { (a, b) } else { (b, a) };
        work.scan(add(long.len(), 1)?)?;
        let apex = long.find_max_by_int_peak(options.use_smoothed_intensities)?;
        let (left, right) = short.fwhm_borders();
        let rt = finite(long[apex].rt())?;
        if rt < finite(short[left].rt())? || rt > finite(short[right].rt())? {
            return Ok(0.0);
        }
    }
    cosine(&x, &y, work)
}
fn cosine(x: &[f64], y: &[f64], work: &mut Work) -> Result<f64> {
    if x.len() != y.len() {
        return Ok(0.0);
    }
    work.scan(x.len())?;
    let (mut dot, mut left, mut right) = (0.0, 0.0, 0.0);
    for (&a, &b) in x.iter().zip(y) {
        dot = finite(dot + a * b)?;
        left = finite(left + a * a)?;
        right = finite(right + b * b)?;
    }
    let denominator = finite(left.sqrt() * right.sqrt())?;
    if denominator > 0.0 {
        finite(dot / denominator)
    } else {
        Ok(0.0)
    }
}
fn averagine_score(intensities: &[f64], mass: f64, work: &mut Work) -> Result<f64> {
    let n = intensities.len();
    // Five peptide elements, each at most31 exponent bits, <=2 convolutions/bit.
    // Preserve existing f64 isotope probabilities; no f32-intermediate parity claim.
    work.scan(add(1024, mul(400, mul(n, n)?)?)?)?;
    work.allocate::<u8>(add(4096, mul(n, 6400)?)?)?;
    let formula = AveragineComposition::PEPTIDE
        .estimate_average_mass(mass)?
        .formula;
    let distribution =
        CoarseIsotopePatternGenerator::new(Some(n), crate::chemistry::CoarseMassMode::Approximate)?
            .run(&formula)?;
    if distribution.peaks().len() < n {
        return Err(bad("averagine pattern has fewer peaks than hypothesis"));
    }
    let mut max = 0.0f64;
    let mut theoretical = 0.0f64;
    for (i, &intensity) in intensities.iter().enumerate() {
        max = max.max(intensity);
        theoretical = theoretical.max(distribution.peaks()[i].probability);
    }
    let mut x = work.vector(n)?;
    let mut y = work.vector(n)?;
    for (i, &intensity) in intensities.iter().enumerate() {
        x.push(finite(distribution.peaks()[i].probability / theoretical)?);
        y.push(finite(intensity / max)?);
    }
    cosine(&x, &y, work)
}
#[allow(clippy::too_many_arguments)]
fn accept(
    traces: &[MassTrace],
    hypotheses: Vec<Hypothesis>,
    options: &FeatureFindingMetaboOptions,
    model: IsotopeFilteringModel,
    has_im: bool,
    rng: &mut UniqueIdGenerator,
    work: &mut Work,
) -> Result<FeatureFindingMetaboOutput> {
    let mut result = FeatureFindingMetaboOutput::default();
    let mut excluded = BTreeSet::<&str>::new();
    let report_smooth = options.use_smoothed_intensities && options.report_smoothed_intensities;
    for hypothesis in hypotheses {
        work.scan(hypothesis.members.len())?;
        let mut collision = false;
        for &index in &hypothesis.members {
            let label = traces[index].label();
            work.compare_key(label, excluded.len())?;
            if excluded.contains(label) {
                collision = true;
                break;
            }
        }
        if collision {
            continue;
        }
        let mut legal = -1i32;
        if hypothesis.members.len() > 1 {
            let predictor_model = match model {
                IsotopeFilteringModel::Metabolites2 => Some(predictor::Model::Noise2),
                IsotopeFilteringModel::Metabolites5 => Some(predictor::Model::Noise5),
                _ => None,
            };
            if let Some(model) = predictor_model {
                let mut values = work.vector(hypothesis.members.len())?;
                for &index in &hypothesis.members {
                    values.push(work.intensity(&traces[index], options.use_smoothed_intensities)?);
                }
                let mass =
                    finite(traces[hypothesis.members[0]].centroid_mz() * hypothesis.charge as f64)?
                        .min(1000.0);
                let mut features = [mass, 0.0, 0.0, 0.0];
                for i in 1..values.len().min(4) {
                    features[i] = finite(values[i] / values[0])?;
                }
                legal = i32::from(predictor::predict(
                    model,
                    features,
                    &mut work.remaining,
                    &mut work.bytes,
                )?);
            }
        }
        if legal == 0 || (options.remove_single_traces && hypothesis.charge == 0) {
            continue;
        }
        if result.features.features.len() >= work.limits.max_features {
            return Err(resource());
        }
        let mut feature = make_feature(
            traces,
            &hypothesis,
            options,
            report_smooth,
            has_im,
            legal,
            work,
        )?;
        work.scan(400)?; // one MT19937-64 draw, including a possible state twist
        feature.unique_id = rng.get_unique_id();
        if options.report_chromatograms && feature.intensity != 0.0 {
            let group = make_chromatograms(traces, &hypothesis, feature.unique_id, work)?;
            work.push(&mut result.chromatograms, group)?;
        }
        work.push(&mut result.features.features, feature)?;
        for &index in &hypothesis.members {
            let label = traces[index].label();
            work.compare_key(label, excluded.len())?;
            if !excluded.contains(label) {
                work.allocate::<u8>(512)?;
                excluded.insert(label);
            }
        }
    }
    work.scan(400)?;
    result.features.unique_id = rng.get_unique_id();
    work.sort(result.features.features.len())?;
    work.allocate::<Feature>(result.features.features.len())?;
    result
        .features
        .features
        .sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
    Ok(result)
}
fn borrowed_hypothesis<'a>(
    traces: &'a [MassTrace],
    hypothesis: &Hypothesis,
    work: &mut Work,
) -> Result<FeatureHypothesis<'a>> {
    // Each pointer may be copied once at insertion and during geometric growth.
    work.scan(mul(hypothesis.members.len(), 4)?)?;
    work.allocate::<&MassTrace>(mul(hypothesis.members.len(), 4)?)?;
    let mut result = FeatureHypothesis::with_limits(FeatureHypothesisLimits {
        max_traces: work.limits.max_traces,
        max_peaks: work.limits.max_peaks,
        max_work: work.limits.max_work,
        max_bytes: work.limits.max_bytes,
    });
    for &index in &hypothesis.members {
        result.add_mass_trace(&traces[index])?;
    }
    result.set_score(hypothesis.score);
    result.set_charge(hypothesis.charge);
    Ok(result)
}
#[allow(clippy::too_many_arguments)]
fn make_feature(
    traces: &[MassTrace],
    hypothesis: &Hypothesis,
    options: &FeatureFindingMetaboOptions,
    smooth: bool,
    has_im: bool,
    legal: i32,
    work: &mut Work,
) -> Result<Feature> {
    let first = &traces[hypothesis.members[0]];
    let intensity = if options.report_summed_intensities {
        let mut sum = 0.0;
        for &index in &hypothesis.members {
            sum = finite(sum + work.intensity(&traces[index], smooth)?)?;
        }
        sum
    } else {
        work.intensity(first, smooth)?
    };
    let mut feature = Feature::new(first.centroid_rt(), first.centroid_mz(), narrow(intensity)?);
    feature.quality = narrow(hypothesis.score)?;
    feature.charge =
        i32::try_from(hypothesis.charge).map_err(|_| bad("feature charge exceeds signed32"))?;
    // Up to eleven entries plus FWHM, including sparse-tree minimum allocation.
    work.allocate::<u8>(add(1024, mul(24 * 12, size_of::<(String, MetaValue)>())?)?)?;
    work.scan(512)?;
    feature.set_width(narrow(first.fwhm())?)?;
    let borrowed = borrowed_hypothesis(traces, hypothesis, work)?;
    let mut label_bytes = hypothesis.members.len().saturating_sub(1);
    for &index in &hypothesis.members {
        label_bytes = add(label_bytes, traces[index].label().len())?;
    }
    work.scan(mul(hypothesis.members.len(), 2)?)?;
    work.allocate::<u8>(label_bytes)?;
    set_meta(
        &mut feature,
        "label",
        MetaValue::from(borrowed.label()?),
        work,
    )?;
    let mut apex = 0.0f64;
    let mut intensities = work.vector(hypothesis.members.len())?;
    let mut rts = work.vector(hypothesis.members.len())?;
    let mut mzs = work.vector(hypothesis.members.len())?;
    let mut ims = if has_im {
        work.vector(hypothesis.members.len())?
    } else {
        Vec::new()
    };
    let mut distances = work.vector(hypothesis.members.len().saturating_sub(1))?;
    let mut previous = None;
    for &index in &hypothesis.members {
        let trace = &traces[index];
        work.scan(if smooth {
            trace.smoothed_intensities().len()
        } else {
            trace.len()
        })?;
        apex = apex.max(trace.max_intensity(smooth)?);
        intensities.push(work.intensity(trace, smooth)?);
        rts.push(trace.centroid_rt());
        mzs.push(trace.centroid_mz());
        if has_im {
            ims.push(trace.centroid_im());
        }
        if let Some(mz) = previous {
            distances.push(finite(trace.centroid_mz() - mz)?);
        }
        previous = Some(trace.centroid_mz());
    }
    // MetaValue validates each transferred float-list payload, without cloning it.
    work.scan(mul(hypothesis.members.len(), if has_im { 5 } else { 4 })?)?;
    set_meta(&mut feature, "max_height", MetaValue::try_from(apex)?, work)?;
    set_meta(
        &mut feature,
        "num_of_masstraces",
        MetaValue::from(i64::try_from(intensities.len()).map_err(|_| resource())?),
        work,
    )?;
    set_meta(
        &mut feature,
        "masstrace_intensity",
        MetaValue::try_from(intensities)?,
        work,
    )?;
    set_meta(
        &mut feature,
        "masstrace_centroid_rt",
        MetaValue::try_from(rts)?,
        work,
    )?;
    set_meta(
        &mut feature,
        "masstrace_centroid_mz",
        MetaValue::try_from(mzs)?,
        work,
    )?;
    if has_im {
        set_meta(
            &mut feature,
            "masstrace_centroid_im",
            MetaValue::try_from(ims)?,
            work,
        )?;
    }
    set_meta(
        &mut feature,
        "isotope_distances",
        MetaValue::try_from(distances)?,
        work,
    )?;
    set_meta(
        &mut feature,
        "legal_isotope_pattern",
        MetaValue::from(legal),
        work,
    )?;
    if options.report_convex_hulls {
        work.allocate::<crate::kernel::ConvexHull2D>(hypothesis.members.len())?;
        work.scan(hypothesis.members.len())?;
        for &index in &hypothesis.members {
            let n = traces[index].len();
            work.sort(n)?;
            work.scan(mul(n, 4)?)?;
            work.allocate::<[f64; 16]>(n)?;
        }
        feature.convex_hulls = borrowed.convex_hulls()?;
    }
    Ok(feature)
}
fn make_chromatograms(
    traces: &[MassTrace],
    hypothesis: &Hypothesis,
    id: u64,
    work: &mut Work,
) -> Result<Vec<MSChromatogram>> {
    let borrowed = borrowed_hypothesis(traces, hypothesis, work)?;
    work.allocate::<MSChromatogram>(hypothesis.members.len())?;
    work.allocate::<u8>(64)?;
    work.scan(hypothesis.members.len())?;
    for &index in &hypothesis.members {
        let n = traces[index].len();
        work.scan(mul(n, 2)?)?;
        work.sort(n)?;
        work.allocate::<crate::kernel::ChromatogramPeak>(mul(n, 2)?)?;
        work.allocate::<u8>(add(4096, mul(24, size_of::<(String, MetaValue)>())?)?)?;
    }
    borrowed.chromatograms(id)
}
fn set_meta(feature: &mut Feature, key: &str, value: MetaValue, work: &mut Work) -> Result<()> {
    work.scan(mul(key.len(), 8)?)?;
    let key = work.text(key)?;
    feature.metadata.insert(key, value);
    Ok(())
}
struct Work {
    limits: FeatureFindingMetaboLimits,
    remaining: usize,
    bytes: usize,
}
impl Work {
    fn new(limits: FeatureFindingMetaboLimits) -> Self {
        Self {
            remaining: limits.max_work,
            bytes: limits.max_bytes,
            limits,
        }
    }
    fn scan(&mut self, n: usize) -> Result<()> {
        self.remaining = self.remaining.checked_sub(n).ok_or_else(resource)?;
        Ok(())
    }
    fn allocate<T>(&mut self, n: usize) -> Result<()> {
        self.scan(n)?;
        self.bytes = self
            .bytes
            .checked_sub(mul(n, size_of::<T>())?)
            .ok_or_else(resource)?;
        Ok(())
    }
    fn vector<T>(&mut self, n: usize) -> Result<Vec<T>> {
        self.allocate::<T>(n)?;
        let mut result = Vec::new();
        result.try_reserve_exact(n).map_err(|_| resource())?;
        Ok(result)
    }
    fn push<T>(&mut self, values: &mut Vec<T>, value: T) -> Result<()> {
        self.scan(1)?;
        if values.len() == values.capacity() {
            let capacity = values.capacity().saturating_mul(2).max(1);
            self.allocate::<T>(capacity)?;
            self.scan(values.len())?;
            values
                .try_reserve_exact(capacity - values.len())
                .map_err(|_| resource())?;
        }
        values.push(value);
        Ok(())
    }
    fn hypothesis(&mut self, values: &mut Vec<Hypothesis>, value: Hypothesis) -> Result<()> {
        if values.len() >= self.limits.max_hypotheses {
            return Err(resource());
        }
        self.push(values, value)
    }
    fn sort(&mut self, n: usize) -> Result<()> {
        if n > 1 {
            self.scan(mul(
                mul(n, (usize::BITS - n.leading_zeros()) as usize)?,
                32,
            )?)?;
        }
        Ok(())
    }
    fn text(&mut self, text: &str) -> Result<String> {
        self.allocate::<u8>(text.len())?;
        let mut result = String::new();
        result
            .try_reserve_exact(text.len())
            .map_err(|_| resource())?;
        result.push_str(text);
        Ok(result)
    }
    fn compare_key(&mut self, key: &str, items: usize) -> Result<()> {
        self.scan(mul(
            add(key.len(), 1)?,
            mul((usize::BITS - items.leading_zeros()) as usize + 1, 12)?,
        )?)
    }
    fn inputs(&mut self, traces: &[MassTrace]) -> Result<()> {
        if traces.len() > self.limits.max_traces {
            return Err(resource());
        }
        self.scan(mul(traces.len(), 8)?)?;
        let mut peaks = 0;
        for trace in traces {
            peaks = add(peaks, trace.len())?;
            if peaks > self.limits.max_peaks {
                return Err(resource());
            }
            finite(trace.centroid_mz())?;
            finite(trace.centroid_rt())?;
            finite(trace.centroid_im())?;
            self.scan(trace.label().len())?;
        }
        Ok(())
    }
    fn intensity(&mut self, trace: &MassTrace, smoothed: bool) -> Result<f64> {
        match trace.quant_method() {
            MassTraceQuantMethod::Area => {
                let (left, right) = trace.fwhm_borders();
                if (left, right) != (0, 0) {
                    self.scan(right - left + 1)?;
                }
            }
            MassTraceQuantMethod::Median => {
                self.scan(trace.len())?;
                if trace.len() > 1 {
                    self.sort(trace.len())?;
                    self.allocate::<f64>(trace.len())?;
                }
            }
            MassTraceQuantMethod::MaxHeight => self.scan(if smoothed {
                trace.smoothed_intensities().len()
            } else {
                trace.len()
            })?,
        }
        trace.intensity(smoothed)
    }
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(resource)
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(resource)
}
fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn resource() -> Error {
    bad("FeatureFindingMetabo resource limit or allocation failure")
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("nonfinite FeatureFindingMetabo data or arithmetic"))
    }
}
fn narrow(value: f64) -> Result<f32> {
    let result = finite(value)? as f32;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(bad("FeatureFindingMetabo output exceeds f32"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::Peak2D;
    fn trace(mz: f64, rts: &[f64], ints: &[f32]) -> MassTrace {
        let mut trace = MassTrace::from_peaks(
            rts.iter()
                .zip(ints)
                .map(|(&r, &i)| Peak2D::new(r, mz, i))
                .collect(),
        )
        .unwrap();
        trace.update_median_mz().unwrap();
        trace.set_quant_method(MassTraceQuantMethod::MaxHeight);
        trace
    }
    fn work() -> Work {
        Work::new(FeatureFindingMetaboLimits::default())
    }
    #[test]
    fn mean_and_element_windows_have_strict_edges_and_separate_gaussian_limits() {
        let a = trace(0., &[0.], &[1.]);
        let options = FeatureFindingMetaboOptions::default();
        let mean = 1.000857 + 0.001091;
        let sd: f64 = 0.0016633 - 0.0004751;
        let sigma = (2. * sd.ln()).exp().sqrt();
        for (mz, expected) in [
            (mean, 1.),
            (mean - 3. * sigma, 0.),
            (mean + 3. * sigma, 0.),
            (mean + sigma, (-0.5f64).exp()),
        ] {
            let b = trace(mz, &[0.], &[1.]);
            let actual = score_mz(&a, &b, 1, 1, None, &options).unwrap();
            assert!((actual - expected).abs() < 1e-12, "{actual} {expected}");
        }
        let options = FeatureFindingMetaboOptions {
            mz_scoring_by_elements: true,
            ..options
        };
        for (mz, expected) in [(1., 0.), (1.5, 1.), (2., 0.)] {
            assert_eq!(
                score_mz(&a, &trace(mz, &[0.], &[1.]), 1, 1, Some((1., 2.)), &options).unwrap(),
                expected
            );
        }
        assert_eq!(
            score_mz(&a, &trace(1.5, &[0.], &[1.]), 1, 1, None, &options).unwrap(),
            0.
        );
        assert_eq!(variance(0.).unwrap(), 0.);
        assert!(variance(-1.).is_err());
    }
    #[test]
    fn isotope_offsets_use_integer_division_and_labeled_or_zero_count_alphabets() {
        assert_eq!(
            isotope_window(&[vec![10., 12.25]], 3, &mut work()).unwrap(),
            Some((3.25, 3.25))
        );
        assert_eq!(
            isotope_window(&[vec![10., 12.25]], 1, &mut work()).unwrap(),
            None
        );
        assert_eq!(
            isotope_window(&[vec![10., 11.125], vec![20., 21.25]], 2, &mut work()).unwrap(),
            Some((2.25, 2.5))
        );
        let mut work = work();
        let empty = crate::chemistry::metabo_elements::isotope_masses(
            "C0+2",
            &mut work.remaining,
            &mut work.bytes,
        )
        .unwrap();
        assert!(empty.is_empty());
        let labeled = crate::chemistry::metabo_elements::isotope_masses(
            "(13)C2",
            &mut work.remaining,
            &mut work.bytes,
        )
        .unwrap();
        assert_eq!(labeled.len(), 1);
        assert_eq!(labeled[0].len(), 1);
        let negative = crate::chemistry::metabo_elements::isotope_masses(
            "C-2H1",
            &mut work.remaining,
            &mut work.bytes,
        )
        .unwrap();
        assert_eq!(negative.len(), 2);
    }
    #[test]
    fn rt_duplicate_groups_require_exactly_two_entries_even_from_one_trace() {
        let mut a = trace(100., &[-0., 0., 1.], &[1., 2., 1.]);
        a.estimate_fwhm(false).unwrap();
        let mut b = trace(101., &[3., 4., 5.], &[1., 2., 1.]);
        b.estimate_fwhm(false).unwrap();
        let options = FeatureFindingMetaboOptions {
            min_isotope_rt_overlap: 0.,
            ..Default::default()
        };
        // The source map has exactly two zero-RT entries, both from a: x=[1], y=[2].
        assert_eq!(score_rt(&a, &b, &options, &mut work()).unwrap(), 1.);
        b.peaks_mut()[0].set_rt(0.); // a third zero entry removes the only pair
        assert_eq!(score_rt(&a, &b, &options, &mut work()).unwrap(), 0.);
    }
    #[test]
    fn shorter_overlap_checks_longer_apex_and_uses_configured_smoothing() {
        let mut a = trace(100., &[0., 1., 2., 3., 4.], &[1., 5., 10., 5., 1.]);
        a.estimate_fwhm(false).unwrap();
        a.set_smoothed_intensities(&[100., 1., 1., 1., 1.]).unwrap();
        let mut b = trace(101., &[1., 2., 3.], &[1., 10., 1.]);
        b.estimate_fwhm(false).unwrap();
        let mut options = FeatureFindingMetaboOptions {
            isotope_rt_overlap_reference: IsotopeRtOverlapReference::Shorter,
            use_smoothed_intensities: false,
            ..Default::default()
        };
        assert!(score_rt(&a, &b, &options, &mut work()).unwrap() > 0.);
        options.use_smoothed_intensities = true;
        assert_eq!(score_rt(&a, &b, &options, &mut work()).unwrap(), 0.);
        options.enable_rt_filtering = false;
        assert_eq!(
            score_rt(&MassTrace::new(), &MassTrace::new(), &options, &mut work()).unwrap(),
            1.
        );
    }
    #[test]
    fn cosine_independent_dot_norm_and_averagine_zero_boundaries() {
        let actual = cosine(&[1., 2., 3.], &[4., 5., 6.], &mut work()).unwrap();
        assert_eq!(actual, 32.0 / (14.0f64.sqrt() * 77.0f64.sqrt()));
        assert_eq!(cosine(&[], &[], &mut work()).unwrap(), 0.);
        assert_eq!(cosine(&[0.], &[1.], &mut work()).unwrap(), 0.);
        assert_eq!(cosine(&[1.], &[], &mut work()).unwrap(), 0.);
        assert!(averagine_score(&[0., 0.], 500., &mut work()).is_err());
        let formula = AveragineComposition::PEPTIDE
            .estimate_average_mass(1000.)
            .unwrap()
            .formula;
        let distribution = CoarseIsotopePatternGenerator::new(
            Some(3),
            crate::chemistry::CoarseMassMode::Approximate,
        )
        .unwrap()
        .run(&formula)
        .unwrap();
        assert!((distribution.peaks()[0].probability - 0.586906).abs() < 2e-6); // upstream existing coarse literal
        let values: Vec<_> = distribution.peaks().iter().map(|p| p.probability).collect();
        assert!((averagine_score(&values, 1000., &mut work()).unwrap() - 1.).abs() < 1e-14);
    }
    #[test]
    fn peptide_score_matches_independent_source_f32_oracle_with_declared_precision_boundary() {
        // Python source-expression oracle, all24 active element traversal orders
        // bounded in tests/data/metabo_peptide_reference.json; not Rust output.
        for (values, expected) in [
            (vec![10., 10.], 0.9572966255102681),
            (vec![10., 5.], 0.9995951520171115),
            (vec![10., 10., 10.], 0.8576672937856934),
            (vec![10., 5., 2.], 0.9991489830493052),
        ] {
            let actual = averagine_score(&values, 1000., &mut work()).unwrap();
            assert!((actual - expected).abs() <= 1e-7, "{actual} {expected}");
        }
    }
}

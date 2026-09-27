// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Map alignment based on peptide identifications.
//!
//! Ports `ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h` and its
//! `.cpp`. See `docs/MAP_ALIGNMENT_SUPPORT.md`.

use super::superimposer::{entry, tree};
use crate::analysis::transformations::{DataPoint, ModelConfig, TransformationDescription};
use crate::concept::progress_logger::ProgressReporter;
use crate::identification::PeptideIdentification;
use crate::identification::graph::{IdentificationData, ScoreTypeId};
use crate::kernel::{BaseFeature, ConsensusMap, FeatureMap};
use crate::math::statistic_functions::median;
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{Error, Result};
use std::collections::BTreeMap;

/// Parameter handler name of the source class.
pub const MAP_ALIGNMENT_ALGORITHM_IDENTIFICATION_NAME: &str = "MapAlignmentAlgorithmIdentification";

/// Retention times per peptide sequence (source `SeqToList`).
pub type SeqToList = BTreeMap<String, Vec<f64>>;
/// One representative retention time per peptide sequence (source
/// `SeqToValue`).
pub type SeqToValue = BTreeMap<String, f64>;

/// The score filter a feature or consensus map selects (source member
/// `better_`, a function pointer).
///
/// It is state: the source sets it only when retention times are collected
/// from a feature or consensus map, from `score_cutoff` and the score
/// direction of that map's first identification, and keeps it for later
/// calls. Retention times collected from a plain identification list use
/// whatever filter is current, which is no filter unless a map set one; so
/// `score_cutoff` has no effect on identification lists aligned before any
/// map. This is reproduced.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ScoreFilter {
    /// Every score passes (the source's initial value).
    #[default]
    All,
    /// A score passes when it is at least the minimum.
    AtLeast,
    /// A score passes when it is at most the minimum.
    AtMost,
}

impl ScoreFilter {
    fn passes(self, score: f64, minimum: f64) -> bool {
        match self {
            Self::All => true,
            Self::AtLeast => score >= minimum,
            Self::AtMost => score <= minimum,
        }
    }
}

/// What an alignment returns.
#[derive(Clone, Debug)]
pub struct IdentificationAlignment {
    /// One transformation per input, in input order: the internal reference
    /// gets `identity`, every other input a description of its
    /// correspondences (source median RT, reference RT, annotated with the
    /// sequence) with no model fitted yet, or `identity` when no reference
    /// RT information is left.
    pub transformations: Vec<TransformationDescription>,
    /// The lines the source writes to `OPENMS_LOG_INFO`, without line ends:
    /// any `Using score type` line of an automatically picked score type, then
    /// "Alignment based on:" and one line per input, framed by the source's
    /// leading and trailing blank lines as empty strings. `OPENMS_LOG_DEBUG`
    /// lines are not reproduced.
    pub info: Vec<String>,
    /// The lines the source writes to `OPENMS_LOG_WARN`.
    pub warnings: Vec<String>,
}

/// A map alignment algorithm based on peptide identifications from MS2
/// spectra (source `MapAlignmentAlgorithmIdentification`).
///
/// Peptide identifications are grouped by the sequence of their best-scoring
/// hit, and their retention times collected; only the best hit per spectrum
/// is taken as the identification. Groups with the same sequence in
/// different maps are points of correspondence between the maps and the
/// basis of the alignment.
///
/// Each map is aligned to a reference retention-time scale, which comes
/// either from a reference ([`set_reference`](Self::set_reference), or an
/// input chosen by index) or is computed as a consensus of the inputs: the
/// median over all maps of each group's median retention time. The median
/// retention time of each group in a map is then mapped onto the group's
/// reference retention time. The source goes on to smooth this mapping with
/// a model the caller fits; the returned descriptions carry the data points.
///
/// The object is stateful, as the source's is: a reference set with
/// [`set_reference`](Self::set_reference) is filtered by, and persists after,
/// each alignment; the score filter persists (see [`ScoreFilter`]); and the
/// effective `min_run_occur` and `min_score` are those of the last
/// alignment.
#[derive(Clone, Debug, PartialEq)]
pub struct MapAlignmentAlgorithmIdentification {
    /// Name of the score type to use for ranking and filtering, for
    /// identification data only (source `score_type`, default empty: picked
    /// automatically).
    pub score_type: String,
    /// Use only identifications above a score cut-off, `min_score` (source
    /// `score_cutoff`, default `false`).
    pub score_cutoff: bool,
    /// With `score_cutoff`, the minimum score for an identification to be
    /// considered (source `min_score`, default 0.05). Unless there are very
    /// few runs or identifications, the source advises increasing it to focus
    /// on more informative peptides.
    pub min_score: f64,
    /// Minimum number of runs, including the reference if any, in which a
    /// peptide must occur to be used (source `min_run_occur`, default 2,
    /// minimum 2). Values above the number of runs are lowered to it, with a
    /// warning.
    pub min_run_occur: i64,
    /// Maximum realistic RT difference of a peptide (median per run versus
    /// reference); peptides with higher shifts are outliers and not used
    /// (source `max_rt_shift`, default 0.5, minimum 0). 0 disables the filter;
    /// above 1 the value is in seconds; at most 1 it is a fraction of the
    /// range of the reference RT scale.
    pub max_rt_shift: f64,
    /// Also use unassigned peptide identifications of feature and consensus
    /// maps (source `use_unassigned_peptides`, default `true`).
    pub use_unassigned_peptides: bool,
    /// For feature and consensus maps, use the retention time of the feature
    /// centroid rather than of the peptide identification; of several
    /// identifications on one feature, only the one closest in RT to the
    /// centroid is used (source `use_feature_rt`, default `false`). Precludes
    /// `use_unassigned_peptides`.
    pub use_feature_rt: bool,
    /// For identification data with adducts, treat differently adducted
    /// variants of the same molecule as different (source `use_adducts`,
    /// default `true`).
    pub use_adducts: bool,
    state: State,
}

/// The source's protected members, set from the parameters at the points
/// where the source copies them.
#[derive(Clone, Debug, Default, PartialEq)]
struct State {
    reference_index: Option<usize>,
    reference: SeqToValue,
    min_run_occur: usize,
    use_feature_rt: bool,
    use_adducts: bool,
    min_score: f64,
    score_cutoff: bool,
    score_type: String,
    better: ScoreFilter,
    /// Info lines logged while collecting retention times.
    log: Vec<String>,
}

impl Default for MapAlignmentAlgorithmIdentification {
    fn default() -> Self {
        Self {
            score_type: String::new(),
            score_cutoff: false,
            min_score: 0.05,
            min_run_occur: 2,
            max_rt_shift: 0.5,
            use_unassigned_peptides: true,
            use_feature_rt: false,
            use_adducts: true,
            state: State::default(),
        }
    }
}

fn flag(value: bool) -> ParamValue {
    ParamValue::String(if value { "true" } else { "false" }.into())
}

impl MapAlignmentAlgorithmIdentification {
    /// The source parameter defaults, as the constructor sets them:
    /// `score_type`, `score_cutoff`, `min_score`, `min_run_occur`,
    /// `max_rt_shift`, `use_unassigned_peptides`, `use_feature_rt` and
    /// `use_adducts`, with their descriptions and restrictions, in
    /// declaration order.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters of this algorithm in the layout of
    /// [`MapAlignmentAlgorithmIdentification::defaults`] (source
    /// `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let boolean = |name: &str, value: bool, description: &str| {
            let mut e = entry(name, flag(value), description, false);
            e.valid_strings = vec!["true".into(), "false".into()];
            e
        };
        let mut occur = entry(
            "min_run_occur",
            ParamValue::Integer(self.min_run_occur),
            "Minimum number of runs (incl. reference, if any) in which a peptide must occur to be used for the alignment.\nUnless you have very few runs or identifications, increase this value to focus on more informative peptides.",
            false,
        );
        occur.min_int = 2;
        let mut shift = entry(
            "max_rt_shift",
            ParamValue::Float(self.max_rt_shift),
            "Maximum realistic RT difference for a peptide (median per run vs. reference). Peptides with higher shifts (outliers) are not used to compute the alignment.\nIf 0, no limit (disable filter); if > 1, the final value in seconds; if <= 1, taken as a fraction of the range of the reference RT scale.",
            false,
        );
        shift.min_float = 0.0;
        tree(vec![
            entry(
                "score_type",
                ParamValue::String(self.score_type.clone()),
                "Name of the score type to use for ranking and filtering (.oms input only). If left empty, a score type is picked automatically.",
                false,
            ),
            boolean(
                "score_cutoff",
                self.score_cutoff,
                "Use only IDs above a score cut-off (parameter 'min_score') for alignment?",
            ),
            entry(
                "min_score",
                ParamValue::Float(self.min_score),
                "If 'score_cutoff' is 'true': Minimum score for an ID to be considered.\nUnless you have very few runs or identifications, increase this value to focus on more informative peptides.",
                false,
            ),
            occur,
            shift,
            boolean(
                "use_unassigned_peptides",
                self.use_unassigned_peptides,
                "Should unassigned peptide identifications be used when computing an alignment of feature or consensus maps? If 'false', only peptide IDs assigned to features will be used.",
            ),
            boolean(
                "use_feature_rt",
                self.use_feature_rt,
                "When aligning feature or consensus maps, don't use the retention time of a peptide identification directly; instead, use the retention time of the centroid of the feature (apex of the elution profile) that the peptide was matched to. If different identifications are matched to one feature, only the peptide closest to the centroid in RT is used.\nPrecludes 'use_unassigned_peptides'.",
            ),
            boolean(
                "use_adducts",
                self.use_adducts,
                "If IDs contain adducts, treat differently adducted variants of the same molecule as different.",
            ),
        ])
    }

    /// Apply parameters, as source `setParameters`. Missing parameters take
    /// their defaults; unknown ones are ignored. The alignment state (the
    /// reference and score filter) is kept, as the source keeps its members.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction (the source's `Exception::InvalidParameter`);
    /// nothing changes then.
    pub fn set_parameters(&mut self, parameters: &Param) -> Result<Vec<String>> {
        let mut handler = DefaultParamHandler::new(MAP_ALIGNMENT_ALGORITHM_IDENTIFICATION_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        let state = self.state.clone();
        let (updated, warnings) = handler.set_parameters_with(parameters, |p| {
            let text = |key: &str| -> Result<String> { Ok(p.value(key)?.as_str()?.to_owned()) };
            let boolean = |key: &str| -> Result<bool> { p.value(key)?.to_bool() };
            Ok(Self {
                score_type: text("score_type")?,
                score_cutoff: boolean("score_cutoff")?,
                min_score: p.value("min_score")?.to_f64()?,
                min_run_occur: p.value("min_run_occur")?.to_i64()?,
                max_rt_shift: p.value("max_rt_shift")?.to_f64()?,
                use_unassigned_peptides: boolean("use_unassigned_peptides")?,
                use_feature_rt: boolean("use_feature_rt")?,
                use_adducts: boolean("use_adducts")?,
                state,
            })
        })?;
        *self = updated;
        Ok(warnings)
    }

    /// A new algorithm with parameters (source construction followed by
    /// `setParameters`).
    ///
    /// # Errors
    ///
    /// As [`MapAlignmentAlgorithmIdentification::set_parameters`].
    pub fn from_param(parameters: &Param) -> Result<Self> {
        let mut algorithm = Self::default();
        algorithm.set_parameters(parameters)?;
        Ok(algorithm)
    }

    /// The reference retention time per sequence, as currently stored.
    pub fn reference(&self) -> &SeqToValue {
        &self.state.reference
    }

    /// The current score filter (see [`ScoreFilter`]).
    pub fn score_filter(&self) -> ScoreFilter {
        self.state.better
    }

    /// Set a reference for the alignment (source `setReference`): the median
    /// retention time per sequence of `data` replaces the stored reference.
    /// Empty data clears the reference.
    ///
    /// `use_feature_rt`, `score_cutoff` and `score_type` are taken from the
    /// parameters here, but `min_score` is not: until the first alignment the
    /// source compares scores with 0 when collecting the reference, and so
    /// does this.
    ///
    /// Returns the lines the source writes to `OPENMS_LOG_INFO` meanwhile: for
    /// identification data whose score type is picked automatically,
    /// `Using score type: <name>`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when no retention time can be
    /// extracted from `data` (the source's `Exception::MissingInformation`),
    /// and the collection errors of [`AlignmentData`]. On error the
    /// reference is left empty, as in the source, which clears it first.
    pub fn set_reference<D: AlignmentData + ?Sized>(&mut self, data: &D) -> Result<Vec<String>> {
        self.state.log.clear();
        let loaded = self.load_reference(data);
        let log = std::mem::take(&mut self.state.log);
        loaded?;
        Ok(log)
    }

    /// `setReference` without draining the info lines.
    fn load_reference<D: AlignmentData + ?Sized>(&mut self, data: &D) -> Result<()> {
        self.state.reference.clear();
        if data.is_empty() {
            return Ok(());
        }
        self.state.use_feature_rt = self.use_feature_rt;
        self.state.score_cutoff = self.score_cutoff;
        self.state.score_type = self.score_type.clone();
        let mut rt_data = SeqToList::new();
        let sorted = data.retention_times(self, &mut rt_data)?;
        let reference = medians(&mut rt_data, sorted)?;
        if reference.is_empty() {
            return Err(Error::MissingInformation(
                "Could not extract retention time information from the reference file".into(),
            ));
        }
        self.state.reference = reference;
        Ok(())
    }

    /// Align inputs (source `align`), without progress reporting.
    ///
    /// # Errors
    ///
    /// As [`MapAlignmentAlgorithmIdentification::align_with_progress`].
    pub fn align<D: AlignmentData>(
        &mut self,
        data: &[D],
        reference_index: Option<usize>,
    ) -> Result<IdentificationAlignment> {
        self.align_with_progress(data, reference_index, &mut ProgressReporter::silent())
    }

    /// Align inputs to a reference retention-time scale (source `align`),
    /// reporting progress ("aligning maps", 0 to 3) as the source does.
    ///
    /// With `reference_index`, that input is the reference and gets the
    /// `identity` transformation; without it, a reference set by
    /// [`set_reference`](Self::set_reference) is used, or a consensus scale is
    /// computed. The steps, as in the source:
    ///
    /// 1. `min_run_occur` is lowered to the number of runs (inputs, plus an
    ///    external reference) with a warning when it exceeds it.
    /// 2. Retention times are collected per input and sequence, and their
    ///    median taken per input.
    /// 3. With a reference, sequences of it that occur in fewer than
    ///    `min_run_occur` runs (the reference counting as one) are removed
    ///    from it; without, the reference is the median of the per-input
    ///    medians of the sequences occurring in at least `min_run_occur` runs.
    /// 4. `max_rt_shift` of at most 1 is scaled by the reference RT range; 0
    ///    becomes unlimited.
    /// 5. Each input's data points are its sequences present in the reference
    ///    whose median is within `max_rt_shift` of the reference; the others
    ///    are counted as outliers.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a `reference_index` past the end
    /// (the source's `Exception::IndexOverflow`), and the errors of
    /// [`set_reference`](Self::set_reference) and of [`AlignmentData`];
    /// progress errors as well. A failed alignment can leave the reference
    /// cleared or filtered, as the source's state is after an exception.
    pub fn align_with_progress<D: AlignmentData>(
        &mut self,
        data: &[D],
        reference_index: Option<usize>,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<IdentificationAlignment> {
        let internal = reference_index.is_some();
        if internal {
            // A reference left over from a previous call is not an extra run.
            self.state.reference.clear();
        }
        let mut warnings = Vec::new();
        self.state.log.clear();
        self.check_parameters(data.len(), &mut warnings)?;
        progress.section(0, 3, "aligning maps", |progress| {
            self.state.reference_index = reference_index;
            if let Some(index) = reference_index {
                if index >= data.len() {
                    return Err(Error::InvalidValue(format!(
                        "reference index {index} is out of range for {} inputs",
                        data.len()
                    )));
                }
                self.load_reference(&data[index])?;
            }
            let mut rt_data = Vec::with_capacity(data.len());
            let mut all_sorted = true;
            for (i, input) in data.iter().enumerate() {
                if reference_index == Some(i) {
                    continue;
                }
                let mut run = SeqToList::new();
                all_sorted &= input.retention_times(self, &mut run)?;
                rt_data.push(run);
            }
            progress.set(1)?;
            let (transformations, summary) =
                self.compute_transformations(&mut rt_data, all_sorted, &mut warnings)?;
            let mut info = std::mem::take(&mut self.state.log);
            info.extend(summary);
            progress.set(2)?;
            progress.set(3)?;
            Ok(IdentificationAlignment {
                transformations,
                info,
                warnings,
            })
        })
    }

    /// Source `checkParameters_`.
    fn check_parameters(&mut self, runs: usize, warnings: &mut Vec<String>) -> Result<()> {
        let mut runs = runs;
        let min_run_occur = usize::try_from(self.min_run_occur)
            .map_err(|_| Error::InvalidValue("'min_run_occur' must not be negative".into()))?;
        self.state.min_run_occur = min_run_occur;
        if !self.state.reference.is_empty() {
            runs += 1;
        }
        self.state.use_feature_rt = self.use_feature_rt;
        if min_run_occur > runs {
            warnings.push(format!(
                "Warning: Value of parameter 'min_run_occur' (here: {min_run_occur}) is higher than the number of runs incl. reference (here: {runs}). Using {runs} instead."
            ));
            self.state.min_run_occur = runs;
        }
        self.state.score_cutoff = self.score_cutoff;
        if self.state.score_cutoff && self.state.score_type.is_empty() {
            self.state.score_type = self.score_type.clone();
        }
        self.state.min_score = self.min_score;
        self.state.use_adducts = self.use_adducts;
        Ok(())
    }

    /// Source `computeTransformations_`: returns the transformations and the
    /// info lines.
    // The run loop goes one past the last run, for a reference at the end, as
    // the source's `i < size + 1` does; an iterator would hide that slot.
    #[allow(clippy::needless_range_loop)]
    fn compute_transformations(
        &mut self,
        rt_data: &mut [SeqToList],
        sorted: bool,
        warnings: &mut Vec<String>,
    ) -> Result<(Vec<TransformationDescription>, Vec<String>)> {
        let size = rt_data.len();
        let medians_per_run = rt_data
            .iter_mut()
            .map(|run| medians(run, sorted))
            .collect::<Result<Vec<_>>>()?;
        let mut medians_per_seq = SeqToList::new();
        for run in &medians_per_run {
            for (sequence, &value) in run {
                medians_per_seq
                    .entry(sequence.clone())
                    .or_default()
                    .push(value);
            }
        }
        let min_run_occur = self.state.min_run_occur;
        let reference_given = !self.state.reference.is_empty();
        if reference_given {
            self.state.reference.retain(|sequence, _| {
                medians_per_seq
                    .get(sequence)
                    .is_some_and(|runs| runs.len() + 1 >= min_run_occur)
            });
        } else {
            medians_per_seq.retain(|_, runs| runs.len() >= min_run_occur);
            self.state.reference = medians(&mut medians_per_seq, false)?;
        }
        if self.state.reference.is_empty() {
            warnings.push("No reference RT information left after filtering!".into());
        }
        let mut max_rt_shift = self.max_rt_shift;
        if max_rt_shift <= 1.0 {
            let mut rt_min = f64::INFINITY;
            let mut rt_max = -rt_min;
            for &value in self.state.reference.values() {
                rt_min = rt_min.min(value);
                rt_max = rt_max.max(value);
            }
            max_rt_shift *= rt_max - rt_min;
        }
        if max_rt_shift == 0.0 {
            max_rt_shift = f64::MAX;
        }

        let mut info = vec![String::new(), "Alignment based on:".to_owned()];
        let mut transformations = Vec::with_capacity(size + 1);
        let mut offset = 0;
        for i in 0..=size {
            if self.state.reference_index == Some(i) {
                transformations.push(identity()?);
                info.push(format!(
                    "- {} data points for sample {} (reference)",
                    self.state.reference.len(),
                    i + 1
                ));
                offset = 1;
            }
            if i >= size {
                break;
            }
            if self.state.reference.is_empty() {
                transformations.push(identity()?);
                continue;
            }
            let mut points = Vec::new();
            let mut outliers = 0usize;
            for (sequence, &value) in &medians_per_run[i] {
                if let Some(&reference) = self.state.reference.get(sequence) {
                    if (value - reference).abs() <= max_rt_shift {
                        points.push(DataPoint::with_note(value, reference, sequence.clone()));
                    } else {
                        outliers += 1;
                    }
                }
            }
            let count = points.len();
            transformations.push(TransformationDescription::new(points)?);
            let mut line = format!("- {count} data points for sample {}", i + offset + 1);
            if outliers > 0 {
                line.push_str(&format!(" ({outliers} outliers removed)"));
            }
            info.push(line);
        }
        info.push(String::new());
        if !reference_given {
            self.state.reference.clear();
        }
        Ok((transformations, info))
    }

    /// Collect the retention times of the best hits of identifications
    /// (source `getRetentionTimes_(const PeptideIdentificationList&, ...)`),
    /// filtered with the current [`ScoreFilter`] against `min_score`.
    fn peptide_retention_times(
        &self,
        peptides: &[PeptideIdentification],
        rt_data: &mut SeqToList,
    ) -> Result<()> {
        for peptide in peptides {
            if let Some(best) = best_hit(peptide) {
                if self.state.better.passes(best.score, self.state.min_score) {
                    let rt = peptide.rt.ok_or_else(|| {
                        Error::MissingInformation(
                            "a peptide identification used for alignment has no retention time"
                                .into(),
                        )
                    })?;
                    rt_data
                        .entry(best.sequence.to_string())
                        .or_default()
                        .push(rt);
                }
            }
        }
        Ok(())
    }

    /// Source `getRetentionTimes_(const IsFCMap auto&, ...)` over the
    /// features' base data and the map's unassigned identifications.
    fn feature_retention_times<'a>(
        &mut self,
        features: impl ExactSizeIterator<Item = &'a BaseFeature> + Clone,
        unassigned: &[PeptideIdentification],
        rt_data: &mut SeqToList,
    ) -> Result<bool> {
        self.state.better = if !self.state.score_cutoff {
            ScoreFilter::All
        } else {
            // features[0].getPeptideIdentifications()[0], which the source
            // reads without checking that either exists.
            let first = features
                .clone()
                .next()
                .and_then(|f| f.peptide_identifications.first())
                .ok_or_else(|| {
                    Error::MissingInformation(
                        "score_cutoff needs an identification on the first feature to tell the score direction"
                            .into(),
                    )
                })?;
            if first.higher_score_better {
                ScoreFilter::AtLeast
            } else {
                ScoreFilter::AtMost
            }
        };
        for feature in features {
            if self.state.use_feature_rt {
                let mut sequence = String::new();
                let mut rt_distance = f64::MAX;
                let mut any_hit = false;
                for peptide in &feature.peptide_identifications {
                    if peptide.hits.is_empty() {
                        continue;
                    }
                    any_hit = true;
                    // The source's missing RT is NaN, which is never closer.
                    let Some(rt) = peptide.rt else { continue };
                    let distance = (rt - feature.rt).abs();
                    if distance < rt_distance {
                        if let Some(best) = best_hit(peptide) {
                            if self.state.better.passes(best.score, self.state.min_score) {
                                sequence = best.sequence.to_string();
                                rt_distance = distance;
                            }
                        }
                    }
                }
                // As in the source, a feature whose identifications all fail
                // the score filter is still recorded, under the empty
                // sequence.
                if any_hit {
                    rt_data.entry(sequence).or_default().push(feature.rt);
                }
            } else {
                self.peptide_retention_times(&feature.peptide_identifications, rt_data)?;
            }
        }
        if !self.state.use_feature_rt && self.use_unassigned_peptides {
            self.peptide_retention_times(unassigned, rt_data)?;
        }
        // Remove duplicates, which occur when an identification was assigned
        // to several features.
        for values in rt_data.values_mut() {
            values.sort_by(f64::total_cmp);
            values.dedup();
        }
        Ok(true)
    }

    /// Source `getRetentionTimes_(const IdentificationData&, ...)`.
    fn identification_retention_times(
        &mut self,
        data: &IdentificationData,
        rt_data: &mut SeqToList,
    ) -> Result<bool> {
        if data.observation_match_count() == 0 {
            return Ok(true);
        }
        let score = self.identification_score_type(data)?;
        let score_type = data.score_type(score)?.clone();
        for id in data.best_matches(score, false)? {
            let hit = data.observation_match(id)?;
            if self.state.score_cutoff {
                match hit.result.score(score) {
                    Some(value) if !score_type.is_better_score(self.state.min_score, value) => {}
                    _ => continue,
                }
            }
            let mut molecule = data.molecule_string(hit.identified_molecule)?;
            if self.state.use_adducts {
                if let Some(adduct) = hit.adduct {
                    molecule.push_str(&format!("+[{}]", data.adduct(adduct)?.name()));
                }
            }
            let rt = data.observation(hit.observation)?.rt.ok_or_else(|| {
                Error::MissingInformation(
                    "an observation used for alignment has no retention time".into(),
                )
            })?;
            rt_data.entry(molecule).or_default().push(rt);
        }
        Ok(false)
    }

    /// Source `handleIdDataScoreType_`: the named score type, or with no name
    /// the most recent score of the first observation match that has one
    /// (source `pickScoreType` over the observation matches), whose name then
    /// becomes the score type.
    fn identification_score_type(&mut self, data: &IdentificationData) -> Result<ScoreTypeId> {
        if self.state.score_type.is_empty() {
            for (_, hit) in data.observation_matches() {
                let picked = hit.result.most_recent_score(|step| {
                    let software = data.processing_step(step)?.software;
                    Ok(data
                        .processing_software(software)?
                        .assigned_scores
                        .as_slice())
                })?;
                if let Some((score, _)) = picked {
                    self.state.score_type = data.score_type(score)?.cv_term.name.clone();
                    self.state
                        .log
                        .push(format!("Using score type: {}", self.state.score_type));
                    return Ok(score);
                }
            }
            return Err(Error::MissingInformation("no scores found".into()));
        }
        data.find_score_type(&self.state.score_type)?
            .ok_or_else(|| {
                Error::MissingInformation(format!(
                    "score type '{}' not found",
                    self.state.score_type
                ))
            })
    }
}

/// The `identity` transformation of a reference or of a run without
/// reference information.
fn identity() -> Result<TransformationDescription> {
    let mut description = TransformationDescription::default();
    description.fit_model(ModelConfig::Identity)?;
    Ok(description)
}

/// Source `getBestScoringHit`: the first hit with the best score in the
/// identification's direction.
fn best_hit(peptide: &PeptideIdentification) -> Option<&crate::identification::PeptideHit> {
    let mut best: Option<&crate::identification::PeptideHit> = None;
    for hit in &peptide.hits {
        let better = match best {
            None => true,
            Some(current) => {
                if peptide.higher_score_better {
                    hit.score > current.score
                } else {
                    hit.score < current.score
                }
            }
        };
        if better {
            best = Some(hit);
        }
    }
    best
}

/// Source `computeMedians_`: the median retention time per sequence.
fn medians(rt_data: &mut SeqToList, _sorted: bool) -> Result<SeqToValue> {
    let mut result = SeqToValue::new();
    for (sequence, values) in rt_data.iter_mut() {
        // Sorting an already sorted list changes nothing, so the source's
        // `sorted` shortcut needs no separate path.
        result.insert(sequence.clone(), median(values)?);
    }
    Ok(result)
}

/// Data an identification-based alignment can collect retention times from
/// (the source's `align` and `setReference` templates, instantiated for
/// `FeatureMap`, `ConsensusMap`, `PeptideIdentificationList` and
/// `IdentificationData`).
pub trait AlignmentData {
    /// Whether there is nothing to collect (`empty()`).
    fn is_empty(&self) -> bool;

    /// Collect retention times per sequence into `rt_data`, returning whether
    /// each list is already sorted (source `getRetentionTimes_`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] for an identification used for
    /// alignment without a retention time (the source records NaN, which
    /// then decides medians and sorts arbitrarily), and as each
    /// implementation documents.
    fn retention_times(
        &self,
        algorithm: &mut MapAlignmentAlgorithmIdentification,
        rt_data: &mut SeqToList,
    ) -> Result<bool>;
}

impl AlignmentData for Vec<PeptideIdentification> {
    fn is_empty(&self) -> bool {
        <[PeptideIdentification]>::is_empty(self)
    }
    fn retention_times(
        &self,
        algorithm: &mut MapAlignmentAlgorithmIdentification,
        rt_data: &mut SeqToList,
    ) -> Result<bool> {
        algorithm.peptide_retention_times(self, rt_data)?;
        Ok(false)
    }
}

impl AlignmentData for FeatureMap {
    fn is_empty(&self) -> bool {
        self.features.is_empty()
    }
    /// # Errors
    ///
    /// With `score_cutoff`, returns [`Error::MissingInformation`] when the
    /// first feature has no identification, where the source reads past the
    /// end of an empty vector.
    fn retention_times(
        &self,
        algorithm: &mut MapAlignmentAlgorithmIdentification,
        rt_data: &mut SeqToList,
    ) -> Result<bool> {
        algorithm.feature_retention_times(
            self.features.iter().map(|f| &f.base),
            &self.unassigned_peptide_identifications,
            rt_data,
        )
    }
}

impl AlignmentData for ConsensusMap {
    fn is_empty(&self) -> bool {
        self.features.is_empty()
    }
    /// # Errors
    ///
    /// As for [`FeatureMap`].
    fn retention_times(
        &self,
        algorithm: &mut MapAlignmentAlgorithmIdentification,
        rt_data: &mut SeqToList,
    ) -> Result<bool> {
        algorithm.feature_retention_times(
            self.features.iter().map(|f| &f.base),
            &self.unassigned_peptide_identifications,
            rt_data,
        )
    }
}

impl AlignmentData for IdentificationData {
    fn is_empty(&self) -> bool {
        IdentificationData::is_empty(self)
    }
    /// The best match per observation under the score type (`score_type`, or
    /// picked automatically), optionally cut at `min_score`, keyed by the
    /// identified molecule (with `+[adduct]` under `use_adducts`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when there is no score to pick,
    /// or the named score type does not exist (the source's
    /// `Exception::MissingInformation`), and graph errors.
    fn retention_times(
        &self,
        algorithm: &mut MapAlignmentAlgorithmIdentification,
        rt_data: &mut SeqToList,
    ) -> Result<bool> {
        algorithm.identification_retention_times(self, rt_data)
    }
}

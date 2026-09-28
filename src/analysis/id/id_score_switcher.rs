// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `ANALYSIS/ID/IDScoreSwitcherAlgorithm.h`: switch the main score of
//! identifications to one stored as a hit meta value.
//!
//! The score types themselves (source `Scores`, which this class re-exports
//! as its deprecated `ScoreType` typedef) are [`ScoreType`] in
//! [`crate::analysis::scores`]. That module also holds an earlier native
//! switcher, [`crate::analysis::scores::ScoreSwitcher`], which refuses to
//! overwrite a differing backup score; [`IDScoreSwitcherAlgorithm`] follows
//! the source exactly and is the one the TOPP tools use.
//!
//! The source templates over "identification types" (`PeptideIdentification`,
//! `ProteinIdentification`); here they implement [`SwitchableIdentification`].

use crate::analysis::scores::{ScoreType, all_score_names};
use crate::identification::{PeptideIdentification, ProteinIdentification};
use crate::kernel::ConsensusMap;
use crate::metadata::{MetaInfo, MetaValue, MetaValueData};
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{Error, Result};

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// Relative tolerance under which a stored old score counts as equal to the
/// main score (source `tolerance_`).
pub const SCORE_TOLERANCE: f64 = 1e-6;

/// An identification whose hits carry a main score and meta values: the
/// source's `IdentificationType` template parameter.
pub trait SwitchableIdentification: Clone {
    /// The main score type (source `getScoreType`).
    fn score_type(&self) -> &str;
    /// Whether higher main scores are better (source `isHigherScoreBetter`).
    fn higher_score_better(&self) -> bool;
    /// Set the main score type and orientation.
    fn set_score_type(&mut self, score_type: &str, higher_better: bool);
    /// Meta values of the first hit, if any.
    fn first_hit_metadata(&self) -> Option<&MetaInfo>;
    /// Every hit's main score and meta values, with a short label for
    /// messages, in hit order.
    fn hit_scores_mut(&mut self) -> Vec<(&mut f64, &mut MetaInfo, String)>;
    /// Validate the whole identification.
    ///
    /// # Errors
    ///
    /// As the type's own `validate`.
    fn validate_identification(&self) -> Result<()>;
}

impl SwitchableIdentification for PeptideIdentification {
    fn score_type(&self) -> &str {
        &self.score_type
    }
    fn higher_score_better(&self) -> bool {
        self.higher_score_better
    }
    fn set_score_type(&mut self, score_type: &str, higher_better: bool) {
        self.score_type = score_type.into();
        self.higher_score_better = higher_better;
    }
    fn first_hit_metadata(&self) -> Option<&MetaInfo> {
        self.hits.first().map(|hit| &hit.metadata)
    }
    fn hit_scores_mut(&mut self) -> Vec<(&mut f64, &mut MetaInfo, String)> {
        self.hits
            .iter_mut()
            .map(|hit| {
                let label = format!("peptide hit '{}'", hit.sequence);
                (&mut hit.score, &mut hit.metadata, label)
            })
            .collect()
    }
    fn validate_identification(&self) -> Result<()> {
        self.validate()
    }
}

impl SwitchableIdentification for ProteinIdentification {
    fn score_type(&self) -> &str {
        &self.score_type
    }
    fn higher_score_better(&self) -> bool {
        self.higher_score_better
    }
    fn set_score_type(&mut self, score_type: &str, higher_better: bool) {
        self.score_type = score_type.into();
        self.higher_score_better = higher_better;
    }
    fn first_hit_metadata(&self) -> Option<&MetaInfo> {
        self.hits.first().map(|hit| &hit.metadata)
    }
    fn hit_scores_mut(&mut self) -> Vec<(&mut f64, &mut MetaInfo, String)> {
        self.hits
            .iter_mut()
            .map(|hit| {
                let label = format!("protein hit '{}'", hit.accession);
                (&mut hit.score, &mut hit.metadata, label)
            })
            .collect()
    }
    fn validate_identification(&self) -> Result<()> {
        self.validate()
    }
}

/// Result of [`IDScoreSwitcherAlgorithm::find_score_type`] (source
/// `ScoreSearchResult`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ScoreSearchResult {
    /// True if the main score is already of the requested type.
    pub is_main_score_type: bool,
    /// The main score name when `is_main_score_type`, else the name of the
    /// meta value holding a score of the requested type, else empty.
    pub score_name: String,
}

/// Score name, orientation and type found by
/// [`IDScoreSwitcherAlgorithm::determine_score_name_orientation_and_type`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScoreNameOrientation {
    /// The main score type name, e.g. `XTandem` or `Percolator_qvalue`.
    pub name: String,
    /// Whether higher main scores are better.
    pub higher_better: bool,
    /// The score category of `name`, when [`ScoreType::from_name`] knows it.
    /// The source leaves its output parameter unchanged otherwise.
    pub score_type: Option<ScoreType>,
}

/// Original and requested score of a switch, used to switch back (source
/// `IDSwitchResult`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IDSwitchResult {
    /// The main score name before the switch.
    pub original_score_name: String,
    /// Whether higher original scores are better.
    pub original_score_higher_better: bool,
    /// The type of the original score ([`ScoreType::Raw`] when unknown).
    pub original_score_type: ScoreType,
    /// Whether higher requested scores are better.
    pub requested_score_higher_better: bool,
    /// The type of the requested score.
    pub requested_score_type: ScoreType,
    /// The main score name after the switch: a search-engine score name
    /// (e.g. `X!Tandem_score`) or a score category (e.g. `PEP`).
    pub requested_score_name: String,
    /// Whether the main score was switched.
    pub score_switched: bool,
}

impl Default for IDSwitchResult {
    fn default() -> Self {
        Self {
            original_score_name: String::new(),
            original_score_higher_better: true,
            original_score_type: ScoreType::Raw,
            requested_score_higher_better: true,
            requested_score_type: ScoreType::Raw,
            requested_score_name: String::new(),
            score_switched: false,
        }
    }
}

/// Hits switched by one operation, and the warning the source logs when it
/// corrects the requested score orientation.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SwitchCount {
    /// Number of hits whose main score was replaced (the source `counter`).
    pub hits: usize,
    /// Warnings the source logs.
    pub warnings: Vec<String>,
}

/// Switch the main score of identifications to a score stored as a hit meta
/// value (source `IDScoreSwitcherAlgorithm`).
///
/// Switching supports raw scores, E-values, posterior (error) probabilities,
/// FDRs and q-values; the orientation of the new score is set and the old
/// score is kept as a meta value so no data is lost.
///
/// The fields are the source parameters. The source switches
/// ([`switch_to_general_score_type`](Self::switch_to_general_score_type))
/// store the score they chose back into these fields, as the source does into
/// its members.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct IDScoreSwitcherAlgorithm {
    /// `new_score`: name of the meta value to use as the new score.
    pub new_score: String,
    /// `new_score_orientation`: `Some(true)` for `higher_better`,
    /// `Some(false)` for `lower_better`, `None` for the source default `""`,
    /// which the source treats as lower-better.
    pub new_score_orientation: Option<bool>,
    /// `new_score_type`: name to use as the type of the new score; empty
    /// means the same as `new_score`.
    pub new_score_type: String,
    /// `old_score`: name for the meta value storing the old score; empty means
    /// the old score type.
    pub old_score: String,
    /// `proteins`: apply to protein scores instead of PSM scores. Read by the
    /// TOPP tool, not by this class.
    pub proteins: bool,
}

impl IDScoreSwitcherAlgorithm {
    /// The default parameters, in source declaration order.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters in the layout of [`IDScoreSwitcherAlgorithm::defaults`],
    /// carrying the current values (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut param = Param::new();
        param.set_value(
            "new_score",
            ParamValue::String(self.new_score.clone()),
            "Name of the meta value to use as the new score",
            &[],
        )?;
        let orientation = match self.new_score_orientation {
            None => "",
            Some(true) => "higher_better",
            Some(false) => "lower_better",
        };
        param.set_value(
            "new_score_orientation",
            ParamValue::String(orientation.into()),
            "Orientation of the new score (are higher or lower values better?)",
            &[],
        )?;
        param.set_valid_strings(
            "new_score_orientation",
            &["lower_better".to_owned(), "higher_better".to_owned()],
        )?;
        param.set_value(
            "new_score_type",
            ParamValue::String(self.new_score_type.clone()),
            "Name to use as the type of the new score (default: same as 'new_score')",
            &[],
        )?;
        param.set_value(
            "old_score",
            ParamValue::String(self.old_score.clone()),
            "Name to use for the meta value storing the old score (default: old score type)",
            &[],
        )?;
        param.set_value(
            "proteins",
            ParamValue::String(if self.proteins { "true" } else { "false" }.into()),
            "Apply to protein scores instead of PSM scores",
            &[],
        )?;
        param.set_valid_strings("proteins", &["true".to_owned(), "false".to_owned()])?;
        Ok(param)
    }

    /// Build from parameters, as source `setParameters` followed by
    /// `updateMembers_`. Missing parameters take their defaults; unknown ones
    /// are ignored (see [`IDScoreSwitcherAlgorithm::from_param_with_warnings`]).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction, including a `new_score_orientation` left at
    /// its default `""`, which is not one of its valid strings (source
    /// `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`IDScoreSwitcherAlgorithm::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`IDScoreSwitcherAlgorithm::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new("IDScoreSwitcherAlgorithm")?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |param| {
            let text = |key: &str| -> Result<String> { Ok(param.value(key)?.as_str()?.to_owned()) };
            Ok(Self {
                new_score: text("new_score")?,
                new_score_orientation: match text("new_score_orientation")?.as_str() {
                    "" => None,
                    other => Some(other == "higher_better"),
                },
                new_score_type: text("new_score_type")?,
                old_score: text("old_score")?,
                proteins: text("proteins")? == "true",
            })
        })
    }

    /// The effective orientation of the new score: the source compares the
    /// parameter with `"higher_better"`.
    pub fn higher_better(&self) -> bool {
        self.new_score_orientation == Some(true)
    }

    /// The effective new score type: `new_score_type`, or `new_score` when it
    /// is empty.
    pub fn effective_new_score_type(&self) -> &str {
        if self.new_score_type.is_empty() {
            &self.new_score
        } else {
            &self.new_score_type
        }
    }

    /// Whether `score_name` is of `score_type` (source `isScoreType`, which
    /// delegates to `Scores::isScoreType`): see [`ScoreType::matches`].
    pub fn is_score_type(score_name: &str, score_type: ScoreType) -> bool {
        score_type.matches(score_name)
    }

    /// Source `toScoreTypeEnum`: see [`ScoreType::parse`].
    ///
    /// # Errors
    ///
    /// As [`ScoreType::parse`].
    pub fn to_score_type_enum(score_type: &str) -> Result<ScoreType> {
        ScoreType::parse(score_type)
    }

    /// Source `isScoreTypeHigherBetter`: see [`ScoreType::higher_is_better`].
    pub fn is_score_type_higher_better(score_type: ScoreType) -> bool {
        score_type.higher_is_better()
    }

    /// Every score name used in OpenMS (source `getScoreNames`), e.g.
    /// `q-value`, `ln(hyperscore)`.
    pub fn score_names() -> Vec<&'static str> {
        all_score_names().collect()
    }

    /// Search `id` for a score of `score_type` (source `findScoreType`).
    ///
    /// If the main score is of that type, its name is returned with
    /// `is_main_score_type`. Otherwise the meta values of the *first hit only*
    /// are searched for each known name of the type, in
    /// [`ScoreType::names`] order, trying `<name>` before `<name>_score`. The
    /// name is empty when nothing is found.
    pub fn find_score_type<I: SwitchableIdentification>(
        id: &I,
        score_type: ScoreType,
    ) -> ScoreSearchResult {
        let main = id.score_type();
        if score_type.matches(main) {
            return ScoreSearchResult {
                is_main_score_type: true,
                score_name: main.to_owned(),
            };
        }
        let mut result = ScoreSearchResult::default();
        if let Some(metadata) = id.first_hit_metadata() {
            for name in score_type.names() {
                if metadata.contains_key(*name) {
                    result.score_name = (*name).to_owned();
                    break;
                }
                let suffixed = format!("{name}_score");
                if metadata.contains_key(&suffixed) {
                    result.score_name = suffixed;
                    break;
                }
            }
        }
        result
    }

    /// Switch the main score of every hit of `id` to the meta value
    /// `new_score` (source `switchScores(IdentificationType&, Size&)`),
    /// returning the number of hits switched.
    ///
    /// The old main score is kept as a meta value named `old_score`, or the
    /// old score type when that is empty. If that meta value already exists
    /// and differs from the main score by more than [`SCORE_TOLERANCE`]
    /// (relative, `|a - b| * 2 / (a + b)`), the old score goes to that name
    /// with a `~` appended instead, overwriting any value there; if it exists
    /// and agrees, nothing is stored. This guards e.g. a `q-value` main score
    /// next to a different `q-value` meta value. Finally the identification
    /// takes the effective new score type and orientation.
    ///
    /// The source notes that the difference test can add the backup to some
    /// hits and not to others, when a hit's scores agree by chance.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a hit lacks the `new_score` meta
    /// value (source `Exception::MissingInformation`), when that value or an
    /// existing old-score value is not numeric (source conversion error), or
    /// when the result fails validation. `id` is unchanged on error; the
    /// source leaves the hits before the failing one switched.
    pub fn switch_scores<I: SwitchableIdentification>(&self, id: &mut I) -> Result<usize> {
        let mut next = id.clone();
        let old_score_meta = if self.old_score.is_empty() {
            next.score_type().to_owned()
        } else {
            self.old_score.clone()
        };
        let mut count = 0;
        for (score, metadata, label) in next.hit_scores_mut() {
            let new_value = metadata
                .get(&self.new_score)
                .ok_or_else(|| {
                    bad(format!(
                        "Meta value '{}' not found for {label}",
                        self.new_score
                    ))
                })?
                .as_f64()?;
            match metadata.get(&old_score_meta) {
                Some(value) if !matches!(value.data(), MetaValueData::Empty) => {
                    let stored = value.as_f64()?;
                    if ((stored - *score) * 2.0 / (stored + *score)).abs() > SCORE_TOLERANCE {
                        metadata.insert(format!("{old_score_meta}~"), MetaValue::try_from(*score)?);
                    }
                }
                _ => {
                    metadata.insert(old_score_meta.clone(), MetaValue::try_from(*score)?);
                }
            }
            *score = new_value;
            count += 1;
        }
        next.set_score_type(self.effective_new_score_type(), self.higher_better());
        next.validate_identification()?;
        *id = next;
        Ok(count)
    }

    /// Adopt `score_name` as the new score for `score_type` (the shared part
    /// of the two `switchToGeneralScoreType` overloads).
    fn adopt(&mut self, score_name: String, score_type: ScoreType) -> Vec<String> {
        self.new_score_type = score_name
            .strip_suffix("_score")
            .unwrap_or(&score_name)
            .to_owned();
        self.new_score = score_name;
        let mut warnings = Vec::new();
        if self.higher_better() != score_type.higher_is_better() {
            warnings.push(
                "Requested score type does not match the expected score direction. Correcting!"
                    .to_owned(),
            );
            self.new_score_orientation = Some(score_type.higher_is_better());
        }
        warnings
    }

    /// Switch every identification of `ids` to a score of the general
    /// `score_type` (source `switchToGeneralScoreType` for vectors).
    ///
    /// Only the first identification is inspected
    /// ([`find_score_type`](Self::find_score_type)). If its main score is
    /// already of the type, nothing changes, on the assumption that all other
    /// identifications agree. Otherwise the score found in its first hit's
    /// meta values becomes `new_score`, with the `_score` suffix removed for
    /// `new_score_type`; the orientation is corrected to the type's with a
    /// warning; and every identification is switched.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the first identification has no
    /// score of the type ("First encountered ID does not have the requested
    /// score type."), or as [`switch_scores`](Self::switch_scores). `ids` and
    /// `self` are unchanged on error.
    pub fn switch_to_general_score_type<I: SwitchableIdentification>(
        &mut self,
        ids: &mut [I],
        score_type: ScoreType,
    ) -> Result<SwitchCount> {
        let Some(first) = ids.first() else {
            return Ok(SwitchCount::default());
        };
        let found = Self::find_score_type(first, score_type);
        if found.is_main_score_type {
            return Ok(SwitchCount::default());
        }
        if found.score_name.is_empty() {
            return Err(bad(
                "First encountered ID does not have the requested score type.",
            ));
        }
        let mut next_self = self.clone();
        let warnings = next_self.adopt(found.score_name, score_type);
        let mut next = ids.to_vec();
        let mut hits = 0;
        for id in &mut next {
            hits += next_self.switch_scores(id)?;
        }
        ids.clone_from_slice(&next);
        *self = next_self;
        Ok(SwitchCount { hits, warnings })
    }

    /// Switch the peptide identifications of a consensus map to a score of the
    /// general `score_type` (source `switchToGeneralScoreType(ConsensusMap&)`).
    ///
    /// The features are searched in order: at the first feature with
    /// identifications whose first identification has the type as its main
    /// score, nothing changes; otherwise the first feature whose first
    /// identification has such a score in its first hit's meta values decides
    /// the new score. Unassigned identifications are not searched. Then all
    /// assigned identifications, and the unassigned ones when
    /// `unassigned_peptides_too`, are switched.
    ///
    /// # Errors
    ///
    /// As [`switch_to_general_score_type`](Self::switch_to_general_score_type).
    pub fn switch_consensus_map_to_general_score_type(
        &mut self,
        map: &mut ConsensusMap,
        score_type: ScoreType,
        unassigned_peptides_too: bool,
    ) -> Result<SwitchCount> {
        let mut new_type = String::new();
        for feature in &map.features {
            if let Some(first) = feature.peptide_identifications.first() {
                let found = Self::find_score_type(first, score_type);
                if found.is_main_score_type {
                    return Ok(SwitchCount::default());
                }
                if !found.score_name.is_empty() {
                    new_type = found.score_name;
                    break;
                }
            }
        }
        if new_type.is_empty() {
            return Err(bad(
                "First encountered ID does not have the requested score type.",
            ));
        }
        let mut next_self = self.clone();
        let warnings = next_self.adopt(new_type, score_type);
        let hits = next_self.switch_map(map, unassigned_peptides_too)?;
        *self = next_self;
        Ok(SwitchCount { hits, warnings })
    }

    /// Source `applyFunctionOnPeptideIDs(switchScores, unassigned)`, atomic.
    fn switch_map(&self, map: &mut ConsensusMap, unassigned: bool) -> Result<usize> {
        let mut next = map.clone();
        let mut hits = 0;
        for feature in &mut next.features {
            for id in &mut feature.peptide_identifications {
                hits += self.switch_scores(id)?;
            }
        }
        if unassigned {
            for id in &mut next.unassigned_peptide_identifications {
                hits += self.switch_scores(id)?;
            }
        }
        *map = next;
        Ok(hits)
    }

    /// Switch the peptide identifications of a consensus map to `new_score`
    /// (source `switchScores(ConsensusMap&, Size&, bool)`).
    ///
    /// Nothing changes when the first feature with identifications already has
    /// `new_score` as the main score type of its first identification.
    ///
    /// # Errors
    ///
    /// As [`switch_scores`](Self::switch_scores); `map` is unchanged on error.
    pub fn switch_consensus_map_scores(
        &self,
        map: &mut ConsensusMap,
        unassigned_peptides_too: bool,
    ) -> Result<usize> {
        if let Some(first) = map
            .features
            .iter()
            .find_map(|feature| feature.peptide_identifications.first())
        {
            if first.score_type == self.new_score {
                return Ok(0);
            }
        }
        self.switch_map(map, unassigned_peptides_too)
    }

    /// Switch a list of peptide identifications to `new_score` (source
    /// `switchScores(PeptideIdentificationList&, Size&)`).
    ///
    /// Nothing changes when the list is empty or its first identification
    /// already has `new_score` as its main score type.
    ///
    /// # Errors
    ///
    /// As [`switch_scores`](Self::switch_scores); `ids` is unchanged on error.
    pub fn switch_peptide_scores(&self, ids: &mut [PeptideIdentification]) -> Result<usize> {
        match ids.first() {
            None => return Ok(0),
            Some(first) if first.score_type == self.new_score => return Ok(0),
            _ => {}
        }
        let mut next = ids.to_vec();
        let mut hits = 0;
        for id in &mut next {
            hits += self.switch_scores(id)?;
        }
        ids.clone_from_slice(&next);
        Ok(hits)
    }

    /// Main score name, orientation and type of a peptide identification list
    /// (source `determineScoreNameOrientationAndType(PeptideIdentificationList&)`).
    ///
    /// Only the first identification is inspected, on the assumption that all
    /// agree. `None` for an empty list, where the source leaves its output
    /// parameters unchanged.
    pub fn determine_score_name_orientation_and_type(
        ids: &[PeptideIdentification],
    ) -> Option<ScoreNameOrientation> {
        ids.first().map(|first| ScoreNameOrientation {
            name: first.score_type.clone(),
            higher_better: first.higher_score_better,
            score_type: ScoreType::from_name(&first.score_type),
        })
    }

    /// Main score name, orientation and type of a consensus map (source
    /// `determineScoreNameOrientationAndType(const ConsensusMap&, ...)`).
    ///
    /// The first identification of each feature with identifications is
    /// inspected in turn, stopping at the first whose score name has a known
    /// type; if none has, the name and orientation of the last one inspected
    /// are returned without a type. When no feature has identifications and
    /// `include_unassigned`, the first unassigned identification is used. With
    /// nothing to inspect the name is empty and higher is better.
    pub fn determine_consensus_map_score_name_orientation_and_type(
        map: &ConsensusMap,
        include_unassigned: bool,
    ) -> ScoreNameOrientation {
        let mut result = ScoreNameOrientation {
            name: String::new(),
            higher_better: true,
            score_type: None,
        };
        for feature in &map.features {
            if let Some(first) = feature.peptide_identifications.first() {
                result.name = first.score_type.clone();
                result.higher_better = first.higher_score_better;
                if let Some(found) = ScoreType::from_name(&result.name) {
                    result.score_type = Some(found);
                    return result;
                }
            }
        }
        if result.name.is_empty() && include_unassigned {
            if let Some(first) = map.unassigned_peptide_identifications.first() {
                result.name = first.score_type.clone();
                result.higher_better = first.higher_score_better;
                result.score_type = ScoreType::from_name(&result.name);
            }
        }
        result
    }

    /// The switcher the static source helpers configure: `new_score`,
    /// orientation, `proteins` false and `old_score` empty.
    fn configured(new_score: &str, higher_better: bool) -> Self {
        Self {
            new_score: new_score.to_owned(),
            new_score_orientation: Some(higher_better),
            ..Self::default()
        }
    }

    /// Switch a peptide identification list to the requested score type
    /// (source static `switchToScoreType(PeptideIdentificationList&,
    /// std::string)`).
    ///
    /// The original score name, orientation and type are recorded first (type
    /// [`ScoreType::Raw`] when unknown). An empty request keeps the main score.
    /// Otherwise the request is parsed as a category (e.g. `RAW`, `PEP`,
    /// `q-value`); if it differs from the original type, the list is switched
    /// with [`switch_to_general_score_type`](Self::switch_to_general_score_type)
    /// and `score_switched` is set. The requested name, orientation and type
    /// are then read back from the list.
    ///
    /// # Errors
    ///
    /// As [`ScoreType::parse`] and
    /// [`switch_to_general_score_type`](Self::switch_to_general_score_type);
    /// `ids` is unchanged on error.
    pub fn switch_to_score_type(
        ids: &mut [PeptideIdentification],
        requested_score_type: &str,
    ) -> Result<IDSwitchResult> {
        let mut result = IDSwitchResult::default();
        apply_found(
            Self::determine_score_name_orientation_and_type(ids),
            &mut result.original_score_name,
            &mut result.original_score_higher_better,
            &mut result.original_score_type,
        );
        result.requested_score_name = result.original_score_name.clone();
        result.requested_score_type = result.original_score_type;
        result.requested_score_higher_better = result.original_score_higher_better;
        if requested_score_type.is_empty() {
            return Ok(result);
        }
        result.requested_score_type = ScoreType::parse(requested_score_type)?;
        if result.requested_score_type != result.original_score_type {
            result.requested_score_higher_better = result.requested_score_type.higher_is_better();
            let mut switcher = Self::configured(
                &result.requested_score_name,
                result.requested_score_higher_better,
            );
            switcher.switch_to_general_score_type(ids, result.requested_score_type)?;
            result.score_switched = true;
        }
        apply_found(
            Self::determine_score_name_orientation_and_type(ids),
            &mut result.requested_score_name,
            &mut result.requested_score_higher_better,
            &mut result.requested_score_type,
        );
        Ok(result)
    }

    /// Switch a consensus map to the requested score type (source static
    /// `switchToScoreType(ConsensusMap&, std::string, bool)`); as
    /// [`switch_to_score_type`](Self::switch_to_score_type), using
    /// [`switch_consensus_map_to_general_score_type`](Self::switch_consensus_map_to_general_score_type).
    ///
    /// # Errors
    ///
    /// As [`switch_to_score_type`](Self::switch_to_score_type); `map` is
    /// unchanged on error.
    pub fn switch_consensus_map_to_score_type(
        map: &mut ConsensusMap,
        requested_score_type: &str,
        include_unassigned: bool,
    ) -> Result<IDSwitchResult> {
        let mut result = IDSwitchResult::default();
        let found =
            Self::determine_consensus_map_score_name_orientation_and_type(map, include_unassigned);
        result.original_score_name = found.name;
        result.original_score_higher_better = found.higher_better;
        if let Some(found) = found.score_type {
            result.original_score_type = found;
        }
        result.requested_score_name = result.original_score_name.clone();
        result.requested_score_type = result.original_score_type;
        result.requested_score_higher_better = result.original_score_higher_better;
        if requested_score_type.is_empty() {
            return Ok(result);
        }
        result.requested_score_type = ScoreType::parse(requested_score_type)?;
        if result.requested_score_type != result.original_score_type {
            result.requested_score_higher_better = result.requested_score_type.higher_is_better();
            let mut switcher = Self::configured(
                &result.requested_score_name,
                result.requested_score_higher_better,
            );
            switcher.switch_consensus_map_to_general_score_type(
                map,
                result.requested_score_type,
                include_unassigned,
            )?;
            result.score_switched = true;
        }
        let found =
            Self::determine_consensus_map_score_name_orientation_and_type(map, include_unassigned);
        result.requested_score_name = found.name;
        result.requested_score_higher_better = found.higher_better;
        if let Some(found) = found.score_type {
            result.requested_score_type = found;
        }
        Ok(result)
    }

    /// Switch a peptide identification list back to its original score after
    /// [`switch_to_score_type`](Self::switch_to_score_type) (source static
    /// `switchBackScoreType(PeptideIdentificationList&, IDSwitchResult)`);
    /// nothing happens unless `result.score_switched`.
    ///
    /// # Errors
    ///
    /// As [`switch_peptide_scores`](Self::switch_peptide_scores).
    pub fn switch_back_score_type(
        ids: &mut [PeptideIdentification],
        result: &IDSwitchResult,
    ) -> Result<usize> {
        if !result.score_switched {
            return Ok(0);
        }
        Self::configured(
            &result.original_score_name,
            result.original_score_higher_better,
        )
        .switch_peptide_scores(ids)
    }

    /// Switch a consensus map back to its original score (source static
    /// `switchBackScoreType(ConsensusMap&, IDSwitchResult, bool)`); nothing
    /// happens unless `result.score_switched`.
    ///
    /// # Errors
    ///
    /// As [`switch_consensus_map_scores`](Self::switch_consensus_map_scores).
    pub fn switch_back_consensus_map_score_type(
        map: &mut ConsensusMap,
        result: &IDSwitchResult,
        include_unassigned: bool,
    ) -> Result<usize> {
        if !result.score_switched {
            return Ok(0);
        }
        Self::configured(
            &result.original_score_name,
            result.original_score_higher_better,
        )
        .switch_consensus_map_scores(map, include_unassigned)
    }
}

/// Apply a list's [`ScoreNameOrientation`] to output fields the way the
/// source's output parameters are written: nothing for an empty list, and the
/// type only when known.
fn apply_found(
    found: Option<ScoreNameOrientation>,
    name: &mut String,
    higher_better: &mut bool,
    score_type: &mut ScoreType,
) {
    if let Some(found) = found {
        *name = found.name;
        *higher_better = found.higher_better;
        if let Some(found) = found.score_type {
            *score_type = found;
        }
    }
}

// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Identification score categories and transactional main-score replacement.

use crate::identification::{PeptideIdentification, ProteinIdentification};
use crate::kernel::{ConsensusMap, FeatureMap};
use crate::metadata::{MetaInfo, MetaValue, MetaValueData};
use crate::{Error, Result};

/// The six score categories and their source ordering.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ScoreType {
    Raw,
    RawEValue,
    PosteriorProbability,
    PosteriorErrorProbability,
    Fdr,
    QValue,
}
impl ScoreType {
    pub const ALL: [Self; 6] = [
        Self::Raw,
        Self::RawEValue,
        Self::PosteriorProbability,
        Self::PosteriorErrorProbability,
        Self::Fdr,
        Self::QValue,
    ];
    /// Exact, case-sensitive source names, in source set order.
    pub fn names(self) -> &'static [&'static str] {
        match self {
            Self::Raw => &[
                "MS:1001492",
                "Mascot",
                "OMSSA",
                "SEQUEST:xcorr",
                "XTandem",
                "hyperscore",
                "ln(hyperscore)",
                "mvh",
                "svm",
            ],
            Self::RawEValue => &[
                "E-Value",
                "MS:1002053",
                "MS:1002257",
                "SpecEValue",
                "evalue",
                "expect",
            ],
            Self::PosteriorProbability => &["Posterior Probability"],
            Self::PosteriorErrorProbability => &[
                "MS:1001493",
                "PEP",
                "Posterior Error Probability",
                "pep",
                "posterior_error_probability",
            ],
            Self::Fdr => &["FDR", "false discovery rate", "fdr"],
            Self::QValue => &["MS:1001491", "q-Value", "q-value", "qval", "qvalue"],
        }
    }
    pub fn higher_is_better(self) -> bool {
        matches!(self, Self::Raw | Self::PosteriorProbability)
    }
    /// Parse a category, not a search-engine score name. Source normalization
    /// strips the literal `_score` suffix, then ignores case, spaces, '-' and '_'.
    pub fn parse(category: &str) -> Result<Self> {
        let key: String = normalize_score_name(category)
            .chars()
            .filter(|c| !matches!(c, '-' | '_' | ' '))
            .flat_map(char::to_lowercase)
            .collect();
        match key.as_str() {
            "raw" => Ok(Self::Raw),
            "rawevalue" => Ok(Self::RawEValue),
            "pp" | "posteriorprobability" => Ok(Self::PosteriorProbability),
            "pep" | "posteriorerrorprobability" => Ok(Self::PosteriorErrorProbability),
            "fdr" | "falsediscoveryrate" => Ok(Self::Fdr),
            "qvalue" => Ok(Self::QValue),
            _ => Err(bad(&format!("unknown score category {category:?}"))),
        }
    }
    pub fn matches(self, name: &str) -> bool {
        self.names().contains(&normalize_score_name(name))
    }
    /// Exact lookup, corresponding to source findIDTypeByName.
    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|t| t.names().contains(&name))
    }
    pub fn from_normalized_name(name: &str) -> Option<Self> {
        Self::from_name(normalize_score_name(name))
    }
}
pub fn normalize_score_name(name: &str) -> &str {
    name.strip_suffix("_score").unwrap_or(name)
}
pub fn all_score_names() -> impl Iterator<Item = &'static str> {
    ScoreType::ALL
        .into_iter()
        .flat_map(|t| t.names().iter().copied())
}
fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScoreSearchResult {
    pub is_main_score: bool,
    pub name: String,
}
fn find_score(
    main: &str,
    first_hit: Option<&MetaInfo>,
    category: ScoreType,
) -> Option<ScoreSearchResult> {
    if category.matches(main) {
        return Some(ScoreSearchResult {
            is_main_score: true,
            name: main.into(),
        });
    }
    let metadata = first_hit?;
    for name in category.names() {
        for key in [(*name).to_owned(), format!("{name}_score")] {
            if metadata.contains_key(&key) {
                return Some(ScoreSearchResult {
                    is_main_score: false,
                    name: key,
                });
            }
        }
    }
    None
}
pub fn find_peptide_score(
    id: &PeptideIdentification,
    category: ScoreType,
) -> Option<ScoreSearchResult> {
    find_score(
        &id.score_type,
        id.hits.first().map(|h| &h.metadata),
        category,
    )
}
pub fn find_protein_score(
    id: &ProteinIdentification,
    category: ScoreType,
) -> Option<ScoreSearchResult> {
    find_score(
        &id.score_type,
        id.hits.first().map(|h| &h.metadata),
        category,
    )
}

/// Replace each main score from a numeric hit metadata value and retain the old
/// score. No implicit sorting or rank assignment occurs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScoreSwitcher {
    pub new_score: String,
    pub higher_score_better: bool,
    /// None uses new_score verbatim.
    pub new_score_type: Option<String>,
    /// None uses the record's previous score_type.
    pub old_score: Option<String>,
}
impl ScoreSwitcher {
    pub fn new(new_score: impl Into<String>, higher_score_better: bool) -> Self {
        Self {
            new_score: new_score.into(),
            higher_score_better,
            new_score_type: None,
            old_score: None,
        }
    }
    fn validate(&self) -> Result<()> {
        if self.new_score.is_empty()
            || self.new_score_type.as_ref().is_some_and(String::is_empty)
            || self.old_score.as_ref().is_some_and(String::is_empty)
        {
            return Err(bad("explicit score names must not be empty"));
        }
        Ok(())
    }
    fn new_type(&self) -> &str {
        self.new_score_type.as_deref().unwrap_or(&self.new_score)
    }
    fn replace_hit(&self, score: &mut f64, metadata: &mut MetaInfo, old_type: &str) -> Result<()> {
        let next = metadata
            .get(&self.new_score)
            .ok_or_else(|| bad(&format!("missing score metadata {:?}", self.new_score)))?
            .as_f64()?;
        let old_key = self.old_score.as_deref().unwrap_or(old_type);
        let destination = match metadata.get(old_key) {
            None => Some(old_key.to_owned()),
            Some(value) if matches!(value.data(), MetaValueData::Empty) => Some(old_key.to_owned()),
            Some(value) if different(value.as_f64()?, *score) => Some(format!("{old_key}~")),
            _ => None,
        };
        if let Some(key) = destination {
            if key == "target_decoy" {
                return Err(bad("score backup cannot replace target_decoy metadata"));
            }
            let mut keep_existing = false;
            if let Some(existing) = metadata.get(&key) {
                if !matches!(existing.data(), MetaValueData::Empty) {
                    if different(existing.as_f64()?, *score) {
                        return Err(bad(&format!(
                            "score backup {key:?} would overwrite a different value"
                        )));
                    }
                    keep_existing = true;
                }
            }
            if !keep_existing {
                metadata.insert(key, MetaValue::try_from(*score)?);
            }
        }
        *score = next;
        Ok(())
    }
    pub fn switch_peptides(&self, ids: &mut [PeptideIdentification]) -> Result<usize> {
        self.validate()?;
        let mut next = ids.to_vec();
        let mut count = 0;
        for id in &mut next {
            id.validate()?;
            for hit in &mut id.hits {
                self.replace_hit(&mut hit.score, &mut hit.metadata, &id.score_type)?;
                count += 1;
            }
            id.score_type = self.new_type().into();
            id.higher_score_better = self.higher_score_better;
        }
        ids.clone_from_slice(&next);
        Ok(count)
    }
    pub fn switch_proteins(&self, ids: &mut [ProteinIdentification]) -> Result<usize> {
        self.validate()?;
        let mut next = ids.to_vec();
        let mut count = 0;
        for id in &mut next {
            id.validate()?;
            for hit in &mut id.hits {
                self.replace_hit(&mut hit.score, &mut hit.metadata, &id.score_type)?;
                count += 1;
            }
            id.score_type = self.new_type().into();
            id.higher_score_better = self.higher_score_better;
        }
        ids.clone_from_slice(&next);
        Ok(count)
    }
    /// Top-level assigned peptides plus optional unassigned peptides, matching
    /// the source consensus-map traversal; protein scores are unchanged.
    pub fn switch_consensus_map(
        &self,
        map: &mut ConsensusMap,
        include_unassigned: bool,
    ) -> Result<usize> {
        self.validate()?;
        map.validate()?;
        let mut next = map.clone();
        let mut count = 0;
        for feature in &mut next.features {
            count += self.switch_peptides(&mut feature.peptide_identifications)?;
        }
        if include_unassigned {
            count += self.switch_peptides(&mut next.unassigned_peptide_identifications)?;
        }
        *map = next;
        Ok(count)
    }
    /// Native feature-map equivalent; subordinate records remain unchanged.
    pub fn switch_feature_map(
        &self,
        map: &mut FeatureMap,
        include_unassigned: bool,
    ) -> Result<usize> {
        self.validate()?;
        map.validate()?;
        let mut next = map.clone();
        let mut count = 0;
        for feature in &mut next.features {
            count += self.switch_peptides(&mut feature.peptide_identifications)?;
        }
        if include_unassigned {
            count += self.switch_peptides(&mut next.unassigned_peptide_identifications)?;
        }
        *map = next;
        Ok(count)
    }
}

// Stable symmetric relative comparison, including opposite signs and zero.
// This preserves the source 1e-6 criterion without intermediate overflow.
fn different(a: f64, b: f64) -> bool {
    let scale = a.abs().max(b.abs());
    if scale == 0.0 {
        return false;
    }
    let a = a / scale;
    let b = b / scale;
    (2.0 * (a - b) / (a + b)).abs() > 1e-6
}

/// Switch each record to an available category. Source lookup uses the first
/// hit; the subsequent numeric replacement checks every hit. Unlike C++, each
/// record is examined so an already-correct first record cannot hide later IDs.
pub fn switch_peptides_to_category(
    ids: &mut [PeptideIdentification],
    category: ScoreType,
) -> Result<usize> {
    let mut next = ids.to_vec();
    let mut count = 0;
    for id in &mut next {
        id.validate()?;
        let found = find_peptide_score(id, category)
            .ok_or_else(|| bad("requested score category is unavailable"))?;
        if found.is_main_score {
            continue;
        }
        let mut switcher = ScoreSwitcher::new(&found.name, category.higher_is_better());
        switcher.new_score_type = Some(normalize_score_name(&found.name).into());
        count += switcher.switch_peptides(std::slice::from_mut(id))?;
    }
    ids.clone_from_slice(&next);
    Ok(count)
}
pub fn switch_proteins_to_category(
    ids: &mut [ProteinIdentification],
    category: ScoreType,
) -> Result<usize> {
    let mut next = ids.to_vec();
    let mut count = 0;
    for id in &mut next {
        id.validate()?;
        let found = find_protein_score(id, category)
            .ok_or_else(|| bad("requested score category is unavailable"))?;
        if found.is_main_score {
            continue;
        }
        let mut switcher = ScoreSwitcher::new(&found.name, category.higher_is_better());
        switcher.new_score_type = Some(normalize_score_name(&found.name).into());
        count += switcher.switch_proteins(std::slice::from_mut(id))?;
    }
    ids.clone_from_slice(&next);
    Ok(count)
}

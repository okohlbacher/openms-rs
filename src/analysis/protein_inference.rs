// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Basic protein inference by aggregation of representative peptide scores.
//!
//! Uses the pinned BasicProteinInferenceAlgorithm and IDBoostGraph algorithms.
//! Inference is transactional; run identifiers are exact, unique associations.
//! See `docs/PROTEIN_INFERENCE_SUPPORT.md` for source differences and limits.

use super::id_filter::{keep_n_best_peptide_hits, remove_dangling_protein_references};
use super::scores::{ScoreSwitcher, ScoreType, find_peptide_score, normalize_score_name};
use crate::identification::{
    PeptideIdentification, ProteinGroup, ProteinIdentification, TargetDecoyType,
};
use crate::kernel::ConsensusMap;
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

/// Source aggregation methods. `Mean` is named `sum` in C++, which divides by count.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AggregationMethod {
    #[default]
    Best,
    /// Multiply positive contributions; skip zero and negative contributions.
    Product,
    Mean,
}
impl AggregationMethod {
    pub fn from_source_name(name: &str) -> Result<Self> {
        match name {
            "best" | "maximum" => Ok(Self::Best),
            "product" => Ok(Self::Product),
            "sum" => Ok(Self::Mean),
            _ => Err(bad("unknown protein score aggregation method")),
        }
    }
    pub fn source_name(self) -> &'static str {
        match self {
            Self::Best => "best",
            Self::Product => "product",
            Self::Mean => "sum",
        }
    }
}

/// Native configuration with source defaults and bounded working collections.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BasicProteinInference {
    pub aggregation: AggregationMethod,
    pub min_peptides_per_protein: usize,
    pub treat_charge_variants_separately: bool,
    pub treat_modification_variants_separately: bool,
    pub use_shared_peptides: bool,
    pub skip_count_annotation: bool,
    pub annotate_indistinguishable_groups: bool,
    pub greedy_group_resolution: bool,
    /// None uses the current main score. Source categories are RAW, PEP and q-value.
    pub score_type: Option<ScoreType>,
    /// Maximum total input peptide and protein hits, before filtering or cloning.
    pub max_input_hits: usize,
    /// Maximum total input peptide evidence records, before filtering or cloning.
    pub max_input_evidences: usize,
}
impl Default for BasicProteinInference {
    fn default() -> Self {
        Self {
            aggregation: AggregationMethod::Best,
            min_peptides_per_protein: 1,
            treat_charge_variants_separately: true,
            treat_modification_variants_separately: true,
            use_shared_peptides: true,
            skip_count_annotation: false,
            annotate_indistinguishable_groups: true,
            greedy_group_resolution: false,
            score_type: None,
            max_input_hits: 1_000_000,
            max_input_evidences: 5_000_000,
        }
    }
}

/// Counts describe the committed result. Spectrum-level records are retained.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ProteinInferenceReport {
    pub protein_runs: usize,
    pub peptide_hits_removed: usize,
    pub proteins_removed: usize,
    pub indistinguishable_groups: usize,
    /// Number of top PSMs whose shared evidence was reduced by greedy resolution.
    pub resolved_peptide_hits: usize,
}

impl BasicProteinInference {
    /// Infer all supplied runs independently and commit both slices on success.
    ///
    /// The best original-score candidate is retained per spectrum **before**
    /// selecting an optional score category, matching the source vector overload.
    /// Retained PSM scores are restored after inference; score backups remain in
    /// metadata. A one-element protein slice uses this vector-overload selection
    /// for one run; `run_single` has a distinct candidate-selection convention.
    /// Duplicate run IDs/accessions, unknown runs, incompatible scores, invalid
    /// probabilities and nonfinite retained protein scores are errors.
    pub fn run(
        &self,
        peptide_ids: &mut [PeptideIdentification],
        protein_ids: &mut [ProteinIdentification],
    ) -> Result<ProteinInferenceReport> {
        self.run_impl(peptide_ids, protein_ids, true)
    }

    /// Source single-object convention: switch scores first, sort candidates and
    /// use only the first for aggregation, but retain other candidates unless
    /// evidence cleanup removes them. Unrelated runs remain untouched.
    pub fn run_single(
        &self,
        peptide_ids: &mut [PeptideIdentification],
        protein_id: &mut ProteinIdentification,
    ) -> Result<ProteinInferenceReport> {
        let indices: Vec<_> = peptide_ids
            .iter()
            .enumerate()
            .filter_map(|(i, id)| (id.identifier == protein_id.identifier).then_some(i))
            .collect();
        self.check_input_size(
            indices.iter().map(|&i| &peptide_ids[i]),
            std::slice::from_ref(protein_id),
        )?;
        let mut peptides: Vec<_> = indices.iter().map(|&i| peptide_ids[i].clone()).collect();
        let report = self.run_impl(&mut peptides, std::slice::from_mut(protein_id), false)?;
        for (index, id) in indices.into_iter().zip(peptides) {
            peptide_ids[index] = id;
        }
        Ok(report)
    }

    /// Infer a union protein run from every top-level consensus-feature peptide
    /// identification, ignoring their run identifiers. Optionally include unassigned
    /// IDs. The output protein hits are sorted by inferred score, as in C++.
    /// Map-stored protein runs and excluded unassigned IDs remain untouched.
    pub fn run_consensus_map(
        &self,
        map: &mut ConsensusMap,
        protein_id: &mut ProteinIdentification,
        include_unassigned: bool,
    ) -> Result<ProteinInferenceReport> {
        map.validate()?;
        let selected = || {
            map.features
                .iter()
                .flat_map(|f| &f.peptide_identifications)
                .chain(if include_unassigned {
                    &map.unassigned_peptide_identifications[..]
                } else {
                    &[]
                })
        };
        self.check_input_size(selected(), std::slice::from_ref(protein_id))?;
        let identifiers: Vec<_> = selected().map(|id| id.identifier.clone()).collect();
        let mut peptides: Vec<_> = selected().cloned().collect();
        for id in &mut peptides {
            id.identifier.clone_from(&protein_id.identifier);
        }
        let mut protein = protein_id.clone();
        let report = self.run(&mut peptides, std::slice::from_mut(&mut protein))?;
        protein.sort()?;
        for (id, identifier) in peptides.iter_mut().zip(identifiers) {
            id.identifier = identifier;
        }
        let mut output = peptides.into_iter();
        for feature in &mut map.features {
            for id in &mut feature.peptide_identifications {
                *id = output.next().expect("same number of peptide IDs");
            }
        }
        if include_unassigned {
            for id in &mut map.unassigned_peptide_identifications {
                *id = output.next().expect("same number of peptide IDs");
            }
        }
        *protein_id = protein;
        Ok(report)
    }

    fn run_impl(
        &self,
        peptide_ids: &mut [PeptideIdentification],
        protein_ids: &mut [ProteinIdentification],
        top_one: bool,
    ) -> Result<ProteinInferenceReport> {
        self.preflight(peptide_ids, protein_ids)?;
        let mut peptides = peptide_ids.to_vec();
        let mut proteins = protein_ids.to_vec();
        if top_one {
            keep_n_best_peptide_hits(&mut peptides, 1)?;
        }
        let mut saved = Vec::with_capacity(peptides.len());
        for id in &mut peptides {
            let mut previous = self.switch_score(id)?;
            if let Some(previous) = &mut previous {
                let mut joined: Vec<_> = std::mem::take(&mut id.hits)
                    .into_iter()
                    .zip(std::mem::take(&mut previous.hits))
                    .collect();
                joined.sort_by(|a, b| {
                    if id.higher_score_better {
                        b.0.score.partial_cmp(&a.0.score).unwrap()
                    } else {
                        a.0.score.partial_cmp(&b.0.score).unwrap()
                    }
                });
                (id.hits, previous.hits) = joined.into_iter().unzip();
            } else {
                id.sort()?;
            }
            saved.push(previous);
        }
        let mut by_run: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
        for (index, id) in peptides.iter().enumerate() {
            by_run.entry(&id.identifier).or_default().push(index);
        }
        // Own keys before peptide mutation; indices remain stable even for empty IDs.
        let by_run: BTreeMap<String, Vec<usize>> = by_run
            .into_iter()
            .map(|(key, value)| (key.to_owned(), value))
            .collect();
        let mut report = ProteinInferenceReport {
            protein_runs: proteins.len(),
            ..Default::default()
        };
        for protein in &mut proteins {
            let indices = by_run
                .get(&protein.identifier)
                .map_or(&[][..], Vec::as_slice);
            report.resolved_peptide_hits += self.process_run(&mut peptides, indices, protein)?;
        }
        for (id, previous) in peptides.iter_mut().zip(saved) {
            if let Some(previous) = previous {
                // ScoreSwitcher preserves both scores and enforces backup collisions.
                for (hit, (score, backup_key)) in id.hits.iter_mut().zip(previous.hits) {
                    let mut record = PeptideIdentification {
                        score_type: id.score_type.clone(),
                        hits: vec![std::mem::take(hit)],
                        ..Default::default()
                    };
                    let mut restore = ScoreSwitcher::new(backup_key, previous.higher);
                    restore.new_score_type = Some(previous.score_type.clone());
                    restore.switch_peptides(std::slice::from_mut(&mut record))?;
                    *hit = record.hits.pop().expect("one restored hit");
                    hit.score = score;
                }
                id.score_type = previous.score_type;
                id.higher_score_better = previous.higher;
            }
            id.validate()?;
        }
        // Restoring first keeps sidecar scores aligned when cleanup removes candidates.
        if self.min_peptides_per_protein > 0 {
            remove_dangling_protein_references(&mut peptides, &proteins, true)?;
        }
        report.peptide_hits_removed = peptide_ids.iter().map(|id| id.hits.len()).sum::<usize>()
            - peptides.iter().map(|id| id.hits.len()).sum::<usize>();
        report.proteins_removed = protein_ids.iter().map(|id| id.hits.len()).sum::<usize>()
            - proteins.iter().map(|id| id.hits.len()).sum::<usize>();
        report.indistinguishable_groups = proteins
            .iter()
            .map(|id| id.indistinguishable_groups.len())
            .sum();
        peptide_ids.clone_from_slice(&peptides);
        protein_ids.clone_from_slice(&proteins);
        Ok(report)
    }

    fn preflight(
        &self,
        peptides: &[PeptideIdentification],
        proteins: &[ProteinIdentification],
    ) -> Result<()> {
        self.check_input_size(peptides.iter(), proteins)?;
        if self.score_type.is_some_and(|score| {
            !matches!(
                score,
                ScoreType::Raw | ScoreType::PosteriorErrorProbability | ScoreType::QValue
            )
        }) {
            return Err(bad("inference score category must be RAW, PEP or q-value"));
        }
        let mut runs = BTreeSet::new();
        for protein in proteins {
            protein.validate()?;
            if !runs.insert(protein.identifier.as_str()) {
                return Err(bad("protein inference requires unique run identifiers"));
            }
            let mut accessions = BTreeSet::new();
            for hit in &protein.hits {
                if hit.accession.is_empty() || !accessions.insert(hit.accession.as_str()) {
                    return Err(bad(
                        "protein inference requires nonempty unique accessions per run",
                    ));
                }
            }
        }
        for peptide in peptides {
            peptide.validate()?;
            if !runs.contains(peptide.identifier.as_str()) {
                return Err(bad(
                    "peptide identification references an unknown inference run",
                ));
            }
            for hit in &peptide.hits {
                if hit.sequence.is_empty() {
                    return Err(bad("protein inference requires nonempty peptide sequences"));
                }
            }
        }
        Ok(())
    }

    fn check_input_size<'a>(
        &self,
        peptides: impl Iterator<Item = &'a PeptideIdentification>,
        proteins: &[ProteinIdentification],
    ) -> Result<()> {
        if self.max_input_hits == 0 || self.max_input_evidences == 0 {
            return Err(bad("protein inference resource limits must be positive"));
        }
        let mut hits = 0usize;
        let mut evidences = 0usize;
        for protein in proteins {
            add_bounded(&mut hits, protein.hits.len(), self.max_input_hits)?;
        }
        for peptide in peptides {
            add_bounded(&mut hits, peptide.hits.len(), self.max_input_hits)?;
            for hit in &peptide.hits {
                add_bounded(
                    &mut evidences,
                    hit.evidences.len(),
                    self.max_input_evidences,
                )?;
            }
        }
        Ok(())
    }

    fn switch_score(&self, id: &mut PeptideIdentification) -> Result<Option<SavedScore>> {
        let Some(category) = self.score_type else {
            return Ok(None);
        };
        if id.hits.is_empty() {
            return Ok(None);
        }
        let found = find_peptide_score(id, category)
            .ok_or_else(|| bad("requested protein-inference score category is unavailable"))?;
        if found.is_main_score {
            return Ok(None);
        }
        if id.score_type.is_empty() {
            return Err(bad("score switching requires a named original main score"));
        }
        let mut previous = SavedScore {
            score_type: id.score_type.clone(),
            higher: id.higher_score_better,
            hits: id
                .hits
                .iter()
                .map(|hit| (hit.score, id.score_type.clone()))
                .collect(),
        };
        let mut switcher = ScoreSwitcher::new(&found.name, category.higher_is_better());
        switcher.new_score_type = Some(normalize_score_name(&found.name).into());
        switcher.switch_peptides(std::slice::from_mut(id))?;
        for (hit, (score, key)) in id.hits.iter().zip(&mut previous.hits) {
            if let Some(backup) = hit.metadata.get(&format!("{}~", previous.score_type)) {
                if backup.as_f64().ok() == Some(*score) {
                    key.push('~');
                }
            }
        }
        Ok(Some(previous))
    }

    fn process_run(
        &self,
        peptides: &mut [PeptideIdentification],
        indices: &[usize],
        protein: &mut ProteinIdentification,
    ) -> Result<usize> {
        let first = indices
            .iter()
            .map(|&i| &peptides[i])
            .find(|id| !id.hits.is_empty());
        let (score_type, higher) = first.map_or_else(
            || (protein.score_type.clone(), protein.higher_score_better),
            |id| (id.score_type.clone(), id.higher_score_better),
        );
        let pep_scores = ScoreType::from_normalized_name(&score_type)
            == Some(ScoreType::PosteriorErrorProbability);
        let probability = pep_scores
            || ScoreType::from_normalized_name(&score_type)
                == Some(ScoreType::PosteriorProbability);
        if probability && higher == pep_scores {
            return Err(bad(
                "probability score orientation is inconsistent with its type",
            ));
        }
        if self.greedy_group_resolution && !higher && !pep_scores {
            return Err(bad(
                "greedy protein resolution requires higher-better scores or PEP",
            ));
        }
        // Grouping by sequence then charge preserves the source's lowest-charge accession choice.
        let mut representatives: BTreeMap<String, BTreeMap<i32, usize>> = BTreeMap::new();
        for &index in indices {
            let id = &peptides[index];
            let Some(hit) = id.hits.first() else { continue };
            if id.score_type != score_type || id.higher_score_better != higher {
                return Err(bad(
                    "mixed score names or orientations within an inference run",
                ));
            }
            if probability && !(0.0..=1.0).contains(&hit.score) {
                return Err(bad("inference probabilities must lie in 0..=1"));
            }
            if !self.use_shared_peptides {
                let Some(value) = hit.metadata.get("protein_references") else {
                    continue;
                };
                if value.unit().is_some() {
                    return Err(bad("protein_references must be a unitless string"));
                }
                match value.as_str()? {
                    "non-unique" => continue,
                    "unique" | "unmatched" => {}
                    _ => return Err(bad("invalid protein_references annotation")),
                }
            }
            let sequence = if self.treat_modification_variants_separately {
                hit.sequence.to_string()
            } else {
                hit.sequence.as_str().into()
            };
            let charge = if self.treat_charge_variants_separately {
                hit.charge
            } else {
                0
            };
            let slot = representatives
                .entry(sequence)
                .or_default()
                .entry(charge)
                .or_insert(index);
            if better(hit.score, peptides[*slot].hits[0].score, higher) {
                *slot = index;
            }
        }
        let accession_indices: BTreeMap<_, _> = protein
            .hits
            .iter()
            .enumerate()
            .map(|(i, hit)| (hit.accession.as_str(), i))
            .collect();
        let mut counts = vec![0usize; protein.hits.len()];
        let mut scores = vec![None; protein.hits.len()];
        for charges in representatives.values() {
            let first = *charges
                .first_key_value()
                .expect("representative group is nonempty")
                .1;
            let accessions = peptides[first].hits[0].protein_accessions();
            for &index in charges.values() {
                if peptides[index].hits[0].protein_accessions() != accessions {
                    return Err(bad(
                        "charge variants of a peptide have inconsistent protein accessions",
                    ));
                }
                let hit_score = peptides[index].hits[0].score;
                let score = if pep_scores {
                    1.0 - hit_score
                } else {
                    hit_score
                };
                for accession in &accessions {
                    let Some(&i) = accession_indices.get(accession) else {
                        continue;
                    };
                    counts[i] += 1;
                    let old = scores[i];
                    scores[i] = Some(match self.aggregation {
                        AggregationMethod::Best => old.map_or(score, |old| {
                            if better(score, old, higher || pep_scores) {
                                score
                            } else {
                                old
                            }
                        }),
                        AggregationMethod::Product => {
                            old.unwrap_or(1.0) * if score > 0.0 { score } else { 1.0 }
                        }
                        AggregationMethod::Mean => old.unwrap_or(0.0) + score,
                    });
                    if !scores[i].is_some_and(f64::is_finite) {
                        return Err(bad("protein score aggregation overflow"));
                    }
                }
            }
        }
        protein.protein_groups.clear();
        protein.indistinguishable_groups.clear();
        protein.score_type = if pep_scores {
            "Posterior Probability".into()
        } else {
            score_type
        };
        protein.higher_score_better = higher || pep_scores;
        let mut undefined = BTreeSet::new();
        for (i, hit) in protein.hits.iter_mut().enumerate() {
            if !self.skip_count_annotation {
                hit.metadata.insert(
                    "nr_found_peptides".into(),
                    MetaValue::from(i64::try_from(counts[i]).map_err(|_| {
                        bad("protein peptide count exceeds metadata integer range")
                    })?),
                );
            }
            hit.score = match (self.aggregation, scores[i]) {
                (AggregationMethod::Mean, Some(sum)) => sum / counts[i] as f64,
                (_, Some(score)) => score,
                (AggregationMethod::Product, None) => 1.0,
                _ => {
                    undefined.insert(hit.accession.clone());
                    0.0
                }
            };
        }
        let mut i = 0usize;
        protein.hits.retain(|_| {
            let keep = counts[i] >= self.min_peptides_per_protein;
            i += 1;
            keep
        });
        let resolved = if self.annotate_indistinguishable_groups || self.greedy_group_resolution {
            self.group_and_resolve(peptides, indices, protein, &undefined)?
        } else {
            0
        };
        if protein
            .hits
            .iter()
            .any(|hit| undefined.contains(&hit.accession))
        {
            return Err(bad(
                "unreferenced protein would have a nonfinite score; require at least one peptide or resolve groups",
            ));
        }
        let metadata = &mut protein.search_parameters.metadata;
        metadata.insert("InferenceEngine".into(), "TOPPProteinInference".into());
        metadata.insert(
            "InferenceEngineVersion".into(),
            concat!("openms-rust-", env!("CARGO_PKG_VERSION")).into(),
        );
        metadata.insert(
            "TOPPProteinInference:aggregation_method".into(),
            self.aggregation.source_name().into(),
        );
        for (key, value) in [
            ("use_shared_peptides", self.use_shared_peptides),
            (
                "treat_charge_variants_separately",
                self.treat_charge_variants_separately,
            ),
            (
                "treat_modification_variants_separately",
                self.treat_modification_variants_separately,
            ),
        ] {
            metadata.insert(
                format!("TOPPProteinInference:{key}"),
                i64::from(value).into(),
            );
        }
        protein.validate()?;
        Ok(resolved)
    }

    fn group_and_resolve(
        &self,
        peptides: &mut [PeptideIdentification],
        indices: &[usize],
        protein: &mut ProteinIdentification,
        undefined: &BTreeSet<String>,
    ) -> Result<usize> {
        let accession_indices: BTreeMap<_, _> = protein
            .hits
            .iter()
            .enumerate()
            .map(|(i, hit)| (hit.accession.as_str(), i))
            .collect();
        let mut neighbors = vec![BTreeSet::new(); protein.hits.len()];
        for &i in indices {
            if let Some(hit) = peptides[i].hits.first() {
                for accession in hit.protein_accessions() {
                    if let Some(&p) = accession_indices.get(accession) {
                        neighbors[p].insert(i);
                    }
                }
            }
        }
        let mut same_neighbors: BTreeMap<BTreeSet<usize>, Vec<usize>> = BTreeMap::new();
        for (i, psms) in neighbors.into_iter().enumerate() {
            if !psms.is_empty() {
                if self.greedy_group_resolution
                    && self.aggregation == AggregationMethod::Mean
                    && undefined.contains(&protein.hits[i].accession)
                {
                    return Err(bad(
                        "unscored mean cannot participate in greedy protein resolution",
                    ));
                }
                same_neighbors.entry(psms).or_default().push(i);
            }
        }
        let mut groups: Vec<_> = same_neighbors.into_iter().collect();
        for (_, members) in &mut groups {
            members.sort_by(|&a, &b| protein.hits[a].accession.cmp(&protein.hits[b].accession));
        }
        groups.sort_by(|a, b| {
            protein.hits[a.1[0]]
                .accession
                .cmp(&protein.hits[b.1[0]].accession)
        });
        let mut resolved = 0usize;
        if self.greedy_group_resolution {
            let mut targets = Vec::with_capacity(groups.len());
            for (_, members) in &groups {
                let mut target = false;
                for &p in members {
                    match protein.hits[p].target_decoy_type()? {
                        TargetDecoyType::Target => target = true,
                        TargetDecoyType::Decoy => {}
                        _ => {
                            return Err(bad(
                                "greedy resolution requires protein target_decoy annotations",
                            ));
                        }
                    }
                }
                targets.push(target);
            }
            let mut parents: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
            for (g, (psms, _)) in groups.iter().enumerate() {
                for &i in psms {
                    parents.entry(i).or_default().push(g);
                }
            }
            let mut clusters: BTreeMap<Vec<usize>, Vec<usize>> = BTreeMap::new();
            for (i, group_ids) in parents {
                if group_ids.len() > 1 {
                    clusters.entry(group_ids).or_default().push(i);
                }
            }
            let mut clusters: Vec<_> = clusters.into_iter().collect();
            // Source uses hash iteration; native ordering is first PSM, then lexical accession.
            clusters.sort_by_key(|(_, psms)| psms[0]);
            for (parents, psms) in clusters {
                let mut best = parents[0];
                for &candidate in &parents[1..] {
                    // None carries source -infinity ordering without storing it in a record.
                    let candidate_hit = &protein.hits[groups[candidate].1[0]];
                    let previous_hit = &protein.hits[groups[best].1[0]];
                    let score = (!undefined.contains(&candidate_hit.accession))
                        .then_some(candidate_hit.score);
                    let previous = (!undefined.contains(&previous_hit.accession))
                        .then_some(previous_hit.score);
                    if score > previous
                        || (score == previous
                            && (targets[candidate], groups[candidate].0.len())
                                > (targets[best], groups[best].0.len()))
                    {
                        best = candidate;
                    }
                }
                let mut remove = BTreeSet::new();
                for loser in parents.into_iter().filter(|&g| g != best) {
                    for &p in &groups[loser].1 {
                        remove.insert(protein.hits[p].accession.as_str());
                    }
                    for i in &psms {
                        groups[loser].0.remove(i);
                    }
                }
                for i in psms {
                    let hit = &mut peptides[i].hits[0];
                    let before = hit.evidences.len();
                    hit.evidences
                        .retain(|e| !remove.contains(e.protein_accession.as_str()));
                    if hit.evidences.len() != before {
                        resolved += 1;
                    }
                }
            }
        }
        // The single-object source overload retains lower-ranked candidates;
        // those still count when removing unreferenced proteins after resolution.
        let retained: BTreeSet<_> = indices
            .iter()
            .flat_map(|&i| &peptides[i].hits)
            .flat_map(|hit| hit.protein_accessions().into_iter().map(str::to_owned))
            .collect();
        let mut grouped = BTreeSet::new();
        for (_, members) in groups {
            // Unresolved source groups use max(-1, scores); resolved groups use their score.
            let mut group = ProteinGroup {
                probability: if self.greedy_group_resolution {
                    protein.hits[members[0]].score
                } else {
                    -1.0
                },
                ..Default::default()
            };
            for p in members {
                let hit = &protein.hits[p];
                if self.greedy_group_resolution && !retained.contains(&hit.accession) {
                    continue;
                }
                grouped.insert(hit.accession.clone());
                group.accessions.push(hit.accession.clone());
                if !self.greedy_group_resolution {
                    group.probability = group.probability.max(hit.score);
                }
            }
            if self.annotate_indistinguishable_groups && !group.accessions.is_empty() {
                protein.indistinguishable_groups.push(group);
            }
        }
        if self.greedy_group_resolution {
            protein.hits.retain(|h| retained.contains(&h.accession));
            if self.annotate_indistinguishable_groups {
                for hit in &protein.hits {
                    if !grouped.contains(&hit.accession) {
                        protein.indistinguishable_groups.push(ProteinGroup {
                            probability: hit.score,
                            accessions: vec![hit.accession.clone()],
                            ..Default::default()
                        });
                    }
                }
            }
        }
        ProteinGroup::sort(&mut protein.indistinguishable_groups)?;
        Ok(resolved)
    }
}

struct SavedScore {
    score_type: String,
    higher: bool,
    hits: Vec<(f64, String)>,
}
fn better(score: f64, previous: f64, higher: bool) -> bool {
    if higher {
        score > previous
    } else {
        score < previous
    }
}
fn add_bounded(total: &mut usize, extra: usize, limit: usize) -> Result<()> {
    *total = total
        .checked_add(extra)
        .filter(|&sum| sum <= limit)
        .ok_or_else(|| bad("protein inference input exceeds resource limit"))?;
    Ok(())
}
fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

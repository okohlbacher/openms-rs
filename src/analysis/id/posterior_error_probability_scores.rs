// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The identification-facing members of the posterior error probability
//! model: extracting and transforming search-engine scores, and writing the
//! fitted probabilities back.
//!
//! Ports `extractAndTransformScores`, `updateScores`, `transformScore_` and
//! `getScore_` of `src/openms/include/OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h`
//! (`PosteriorErrorProbabilityModel.cpp:874-1164`, core SDK `bc9cc12`). The
//! numerical model is [`crate::math::posterior_error_probability`]; these
//! members need peptide and protein identifications, which the math module
//! may not depend on, so they live here, next to their consumer
//! (`IDPosteriorErrorProbability`). See
//! `docs/POSTERIOR_ERROR_PROBABILITY_SUPPORT.md`.
//!
//! # Example
//!
//! ```
//! use openms::analysis::id::posterior_error_probability_scores::transform_score;
//! use openms::identification::PeptideHit;
//!
//! let mut hit = PeptideHit::default();
//! hit.score = 1e-3;
//! // X! Tandem's E-value is the hit score here; -log10 transforms it.
//! assert!((transform_score("XTANDEM", &hit, "E-Value")? - 3.0).abs() < 1e-12);
//! # Ok::<(), openms::Error>(())
//! ```

use std::collections::{BTreeMap, BTreeSet};

use crate::identification::{PeptideHit, PeptideIdentification, ProteinIdentification};
use crate::math::posterior_error_probability::PosteriorErrorProbabilityModel;
use crate::metadata::meta_info_description::data_value_text;
use crate::metadata::{MetaValue, MetaValueData};
use crate::{Error, Result};

/// The search engines `extractAndTransformScores` looks for, in the source's
/// order (`PosteriorErrorProbabilityModel.cpp:966-969`), before upper-casing.
pub const SEARCH_ENGINES: [&str; 17] = [
    "XTandem",
    "OMSSA",
    "MASCOT",
    "SpectraST",
    "MyriMatch",
    "SimTandem",
    "MSGFPlus",
    "MS-GF+",
    "Comet",
    "MSFragger",
    "tide-search",
    "Sage",
    "SimpleSearchEngine",
    "ProSE",
    "OpenMS/ConsensusID_best",
    "OpenMS/ConsensusID_worst",
    "OpenMS/ConsensusID_average",
];

/// Scores per search engine (and charge), as `extractAndTransformScores`
/// returns them: the key is the upper-cased engine name, followed by `,` and
/// the charge when scores are split by charge; the value holds all scores,
/// the target scores and the decoy scores, in that order.
pub type EngineScores = BTreeMap<String, Vec<Vec<f64>>>;

/// The source's `getScore_`: the first requested score type that is either
/// the identification's current score type (then the hit's score), or a meta
/// value of the hit under that name or that name plus `_score`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when no requested type is found (the source prints
/// `actual_score_type` to standard output and throws
/// `Exception::UnableToFit`), and when the meta value found is not numeric
/// (the source's `DataValue` conversion to `double` throws
/// `ConversionError`).
pub fn get_score(
    requested_score_types: &[&str],
    hit: &PeptideHit,
    actual_score_type: &str,
) -> Result<f64> {
    for requested in requested_score_types {
        if actual_score_type == *requested {
            return Ok(hit.score);
        }
        if let Some(value) = hit.metadata.get(*requested) {
            return numeric(value, requested);
        }
        let with_suffix = format!("{requested}_score");
        if let Some(value) = hit.metadata.get(&with_suffix) {
            return numeric(value, &with_suffix);
        }
    }
    Err(Error::InvalidValue(format!(
        "Expected score type for search engine not found: None of the expected score types {} for search engine found (actual score type '{actual_score_type}')",
        requested_score_types.join(",")
    )))
}

/// `static_cast<double>(const DataValue&)`: integers and doubles convert,
/// anything else is a conversion error.
fn numeric(value: &MetaValue, name: &str) -> Result<f64> {
    match value.data() {
        MetaValueData::Float(number) => Ok(*number),
        // The source converts its SignedSize to double the same way.
        MetaValueData::Integer(number) => Ok(*number as f64),
        _ => Err(Error::InvalidValue(format!(
            "Could not convert non-numerical meta value '{name}' ('{}') to double",
            data_value_text(value)
        ))),
    }
}

/// The source's `transformScore_`: the per-engine transformation of a hit's
/// score into the model's discriminant score.
///
/// `engine` is the upper-cased engine name. `-log10` for OMSSA, X! Tandem,
/// Mascot, SimTandem, MS-GF+, Comet and MSFragger; `100 *` for SpectraST's
/// f-value; unchanged for MyriMatch's MVH and the hyperscore of
/// SimpleSearchEngine, ProSE and Sage. A Mascot hit with score 0 gives NaN
/// (source issue #740), which the callers skip.
///
/// # Errors
///
/// [`Error::InvalidValue`] for an engine without a transformation - which
/// includes `TIDE-SEARCH` and the three `OPENMS/CONSENSUSID_*` names that
/// [`SEARCH_ENGINES`] nevertheless lists, exactly as in the source (its
/// `UnableToFit`, "The chosen search engine is currently not supported") -
/// and as [`get_score`].
pub fn transform_score(engine: &str, hit: &PeptideHit, current_score_type: &str) -> Result<f64> {
    let get = |types: &[&str]| get_score(types, hit, current_score_type);
    match engine {
        "OMSSA" => Ok(-get(&["OMSSA"])?.log10()),
        "MYRIMATCH" => get(&["mvh"]),
        "XTANDEM" | "SIMTANDEM" => Ok(-get(&["E-Value"])?.log10()),
        "MASCOT" => {
            if hit.score == 0.0 {
                return Ok(f64::NAN);
            }
            Ok(-get(&["EValue", "expect"])?.log10())
        }
        "SPECTRAST" => Ok(100.0 * get(&["f-val"])?),
        "MSGFPLUS" | "MS-GF+" => Ok(-get(&["MS:1002053", "expect"])?.log10()),
        "COMET" => Ok(-get(&["MS:1002257", "expect"])?.log10()),
        "SIMPLESEARCHENGINE" | "PROSE" | "SAGE" => get(&["hyperscore", "ln(hyperscore)"]),
        "MSFRAGGER" => Ok(-get(&["expect"])?.log10()),
        _ => Err(Error::InvalidValue(format!(
            "No parameters for chosen search engine: the chosen search engine '{engine}' is currently not supported"
        ))),
    }
}

/// The source's `extractAndTransformScores`: the transformed scores of every
/// supported search engine's identifications, keyed by engine (and charge).
///
/// For each engine of [`SEARCH_ENGINES`] (upper-cased) and each protein run
/// whose upper-cased search engine equals it - for an `OpenMS/ConsensusID*`
/// run, the part before the first `:` of its `ConsensusIDBaseSearch` meta
/// value - every peptide identification of that run (same identifier) is
/// sorted by score, and either its first hit (`top_hits_only`) or every hit
/// is transformed with [`transform_score`]; NaN results are skipped. With
/// `split_charge` this is repeated per distinct hit charge, keeping only hits
/// of that charge. An engine (and charge) enters the map only with **more
/// than two** scores. With `top_hits_only` and `target_decoy_available`, a
/// top hit whose raw score is below `fdr_for_targets_smaller` is a target,
/// otherwise a decoy; without `top_hits_only` the target and decoy lists stay
/// empty, as in the source.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `split_charge` is set but no hit carries a
/// charge (the source's `Exception::Precondition`); as
/// [`PeptideIdentification::sort`] for non-finite scores; and as
/// [`transform_score`] - in particular for a run of an engine the source
/// lists but cannot transform (`tide-search`, the ConsensusID names).
pub fn extract_and_transform_scores(
    protein_ids: &[ProteinIdentification],
    peptide_ids: &[PeptideIdentification],
    split_charge: bool,
    top_hits_only: bool,
    target_decoy_available: bool,
    fdr_for_targets_smaller: f64,
) -> Result<EngineScores> {
    let mut charges = BTreeSet::new();
    if split_charge {
        for peptide in peptide_ids {
            for hit in &peptide.hits {
                charges.insert(hit.charge);
            }
        }
        if charges.is_empty() {
            return Err(Error::InvalidValue(
                "'split_charge' is set, but the list of charge states is empty".into(),
            ));
        }
    }
    // Without split_charge the source's do-while runs once and never reads
    // the charge.
    let passes: Vec<Option<i32>> = if split_charge {
        charges.into_iter().map(Some).collect()
    } else {
        vec![None]
    };
    let mut all_scores = EngineScores::new();
    for charge in passes {
        for supported in SEARCH_ENGINES {
            let supported = supported.to_ascii_uppercase();
            let (mut scores, mut target, mut decoy) = (Vec::new(), Vec::new(), Vec::new());
            for protein in protein_ids {
                if run_engine(protein) != supported {
                    continue;
                }
                for peptide in peptide_ids {
                    if protein.identifier != peptide.identifier {
                        continue;
                    }
                    let mut peptide = peptide.clone();
                    peptide.sort()?;
                    let wanted = |hit: &PeptideHit| charge.is_none_or(|c| hit.charge == c);
                    if top_hits_only {
                        let Some(top) = peptide.hits.first() else {
                            continue;
                        };
                        if !wanted(top) {
                            continue;
                        }
                        let score = transform_score(&supported, top, &peptide.score_type)?;
                        if score.is_nan() {
                            continue;
                        }
                        scores.push(score);
                        if target_decoy_available {
                            if top.score < fdr_for_targets_smaller {
                                target.push(score);
                            } else {
                                decoy.push(score);
                            }
                        }
                    } else {
                        for hit in peptide.hits.iter().filter(|hit| wanted(hit)) {
                            let score = transform_score(&supported, hit, &peptide.score_type)?;
                            if !score.is_nan() {
                                scores.push(score);
                            }
                        }
                    }
                }
            }
            if scores.len() > 2 {
                let key = match charge {
                    Some(c) => format!("{supported},{c}"),
                    None => supported.clone(),
                };
                // `std::map::insert` keeps an existing entry.
                all_scores.entry(key).or_insert(vec![scores, target, decoy]);
            }
        }
    }
    Ok(all_scores)
}

/// The upper-cased search engine of a run, resolving an
/// `OpenMS/ConsensusID*` run to the engine before the first `:` of its
/// `ConsensusIDBaseSearch` meta value (empty when absent).
fn run_engine(protein: &ProteinIdentification) -> String {
    let engine = if protein.search_engine.starts_with("OpenMS/ConsensusID") {
        let base = protein
            .metadata
            .get("ConsensusIDBaseSearch")
            .map(data_value_text)
            .unwrap_or_default();
        match base.find(':') {
            Some(end) => base[..end].to_owned(),
            None => base,
        }
    } else {
        protein.search_engine.clone()
    };
    engine.to_ascii_uppercase()
}

/// The two flags `updateScores` reports through its reference parameters.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ScoreUpdate {
    /// A computed probability fell outside `[0, 1]`: the fit is invalid.
    pub unable_to_fit_data: bool,
    /// Always `false`: the source's quality check is a TODO.
    pub data_might_not_be_well_fit: bool,
}

/// The source's `updateScores`: replace the scores of the hits of every run
/// of `search_engine` by the model's posterior error probability (or, with
/// `prob_correct`, by one minus it).
///
/// Runs match on the upper-cased search engine exactly - no ConsensusID
/// resolution here, unlike [`extract_and_transform_scores`], as in the
/// source. Each updated hit first stores its old score as the meta value
/// `<score type>_score`; a transformed score of NaN (a Mascot score of 0)
/// becomes probability 1. With `split_charge`, only hits of `charge` are
/// updated. The identification's score type is not changed.
///
/// # Errors
///
/// As [`transform_score`] and
/// [`PosteriorErrorProbabilityModel::compute_probability`]. Divergence: the
/// latter refuses an unfitted model and a non-finite score; the source
/// computes with whatever it holds, which for a transformed score of
/// `+inf`/`-inf` (an E-value of 0 or of infinity) is the limiting probability.
/// Identifications updated before the error keep their new scores, as they
/// would when the source throws.
pub fn update_scores(
    model: &PosteriorErrorProbabilityModel,
    search_engine: &str,
    charge: i32,
    prob_correct: bool,
    split_charge: bool,
    protein_ids: &[ProteinIdentification],
    peptide_ids: &mut [PeptideIdentification],
) -> Result<ScoreUpdate> {
    let engine = search_engine.to_ascii_uppercase();
    let mut update = ScoreUpdate::default();
    for protein in protein_ids {
        if protein.search_engine.to_ascii_uppercase() != engine {
            continue;
        }
        for peptide in peptide_ids.iter_mut() {
            if protein.identifier != peptide.identifier {
                continue;
            }
            let score_type = format!("{}_score", peptide.score_type);
            let mut hits = peptide.hits.clone();
            for hit in hits
                .iter_mut()
                .filter(|hit| !split_charge || hit.charge == charge)
            {
                hit.metadata
                    .insert(score_type.clone(), MetaValue::try_from(hit.score)?);
                let transformed = transform_score(&engine, hit, &peptide.score_type)?;
                let probability = if transformed.is_nan() {
                    1.0
                } else {
                    let probability = model.compute_probability(transformed)?;
                    if !(0.0..=1.0).contains(&probability) {
                        update.unable_to_fit_data = true;
                    }
                    probability
                };
                hit.score = if prob_correct {
                    1.0 - probability
                } else {
                    probability
                };
            }
            peptide.hits = hits;
        }
    }
    Ok(update)
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Percolator feature set and integration helper.
//!
//! Port of `ANALYSIS/ID/PercolatorFeatureSetHelper.h` and `.cpp` at
//! `bc9cc12`. These functions handle (compute, aggregate, integrate)
//! Percolator features: they calculate or extract search-engine specific
//! features for later use with `PercolatorAdapter`, and merge protein
//! identification runs from several search engines.
//!
//! The source is a class of static members; they become free functions.
//! Every function that mutates identifications works on a copy and commits
//! only on success, so an error leaves the input unchanged, where the source
//! throws halfway through a loop. `docs/PERCOLATOR_SUPPORT.md` maps every
//! member and lists the divergences.

use crate::comparison::Tolerance;
use crate::concept::constants::user_param::ISOTOPE_ERROR;
use crate::identification::{
    EnzymeTermSpecificity, PeptideHit, PeptideIdentification, ProteinHit, ProteinIdentification,
};
use crate::metadata::{MetaInfo, MetaValue, MetaValueData};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::BTreeSet;

/// Registers the MS-GF+ Percolator features and makes sure every hit has
/// them.
///
/// MS-GF+ does not always produce all scores, so the source focuses on the
/// main ones: `MS:1002049` (MS-GF:RawScore), `MS:1002050`
/// (MS-GF:DeNovoScore), `MS:1002052` (MS-GF:SpecEValue), `MS:1002053`
/// (MS-GF:EValue), then [`ISOTOPE_ERROR`], are appended to `feature_set` in
/// that order. Each hit lacking one of the four scores gets it as `0.0`.
/// MS-GF+ should be run with its `addFeatures` flag enabled.
///
/// Infallible, as in the source.
pub fn add_msgf_features(peptide_ids: &mut [PeptideIdentification], feature_set: &mut Vec<String>) {
    const SCORES: [&str; 4] = ["MS:1002049", "MS:1002050", "MS:1002052", "MS:1002053"];
    feature_set.extend(SCORES.iter().map(|s| (*s).to_owned()));
    feature_set.push(ISOTOPE_ERROR.to_owned());
    for identification in peptide_ids.iter_mut() {
        for hit in &mut identification.hits {
            for key in SCORES {
                hit.metadata
                    .entry(key.to_owned())
                    .or_insert_with(|| source_float(0.0));
            }
        }
    }
}

/// Registers and computes the X!Tandem Percolator features.
///
/// The ion series `a, b, c, x, y, z` are probed on the **first hit of the
/// first identification**: a series is used when that hit has both an
/// `<ion>_score` and an `<ion>_ions` meta value whose text is non-empty.
/// For each used series `XTANDEM:frac_ion_<ion>` is registered, then
/// `XTANDEM:deltascore`. Then, for the first hit of every identification:
/// `XTANDEM:deltascore = score - nextscore`, and for every used series
/// `XTANDEM:frac_ion_<ion> = <ion>_ions / unmodified sequence length`. Other
/// hits are not touched, as in the source.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `peptide_ids` is empty or an identification
/// has no hits (undefined behaviour in the source, which dereferences
/// `front()`), when a first hit lacks a numeric `nextscore` or a used
/// `<ion>_ions` (the source throws `Exception::ConversionError`), or when a
/// first hit has an empty sequence (the source divides by zero). On error
/// `peptide_ids` and `feature_set` are unchanged.
pub fn add_xtandem_features(
    peptide_ids: &mut [PeptideIdentification],
    feature_set: &mut Vec<String>,
) -> Result<()> {
    const IONS: [&str; 6] = ["a", "b", "c", "x", "y", "z"];
    let first = peptide_ids
        .first()
        .and_then(|id| id.hits.first())
        .ok_or_else(|| invalid("X!Tandem features need at least one identification with a hit"))?;
    let found: Vec<&str> = IONS
        .iter()
        .copied()
        .filter(|ion| {
            text_non_empty(first.metadata.get(&format!("{ion}_score")))
                && text_non_empty(first.metadata.get(&format!("{ion}_ions")))
        })
        .collect();
    let mut updated = Vec::with_capacity(peptide_ids.len());
    for identification in peptide_ids.iter() {
        let hit = identification
            .hits
            .first()
            .ok_or_else(|| invalid("X!Tandem features need a hit in every identification"))?;
        let mut meta = hit.metadata.clone();
        let delta = hit.score - numeric(&hit.metadata, "nextscore")?;
        meta.insert("XTANDEM:deltascore".into(), source_float(delta));
        let length = hit.sequence.len();
        if length == 0 && !found.is_empty() {
            return Err(invalid("X!Tandem ion fractions need a non-empty sequence"));
        }
        for ion in &found {
            let ions = numeric(&hit.metadata, &format!("{ion}_ions"))?;
            meta.insert(
                format!("XTANDEM:frac_ion_{ion}"),
                source_float(ions / length as f64),
            );
        }
        updated.push(meta);
    }
    for (identification, meta) in peptide_ids.iter_mut().zip(updated) {
        identification.hits[0].metadata = meta;
    }
    feature_set.extend(found.iter().map(|ion| format!("XTANDEM:frac_ion_{ion}")));
    feature_set.push("XTANDEM:deltascore".into());
    Ok(())
}

/// Registers the MSFragger Percolator features: `MS:1001330` (expect
/// score), `hyperscore`, `nextscore` and [`ISOTOPE_ERROR`], appended in that
/// order. Nothing is computed.
pub fn add_msfragger_features(extra_features: &mut Vec<String>) {
    extra_features.extend(
        ["MS:1001330", "hyperscore", "nextscore", ISOTOPE_ERROR]
            .iter()
            .map(|s| (*s).to_owned()),
    );
}

/// Collects and registers the andes search engine Percolator features.
///
/// andes annotates every PSM with its own numeric features as meta values
/// named with the `andes:` prefix (e.g. `andes:RankScore`). Rather than
/// hard-coding the feature names, every `andes:`-prefixed meta value that is
/// an integer or a float on any hit is collected, keeping andes the single
/// source of truth for its feature set. The names are appended in byte
/// order without duplicates, as the source's `std::set<std::string>`.
///
/// A key that is numeric on one hit and a string on another is collected,
/// as in the source.
pub fn add_andes_features(peptide_ids: &[PeptideIdentification], feature_set: &mut Vec<String>) {
    let mut andes = BTreeSet::new();
    for identification in peptide_ids {
        for hit in &identification.hits {
            for (key, value) in &hit.metadata {
                if key.starts_with("andes:")
                    && matches!(
                        value.data(),
                        MetaValueData::Integer(_) | MetaValueData::Float(_)
                    )
                {
                    andes.insert(key.clone());
                }
            }
        }
    }
    feature_set.extend(andes);
}

/// Registers and computes the Comet Percolator features.
///
/// Registers, in order: [`ISOTOPE_ERROR`], `COMET:deltaCn`,
/// `COMET:deltaLCn`, `COMET:lnExpect`, `MS:1002252` (XCorr, unchanged),
/// `MS:1002255` (Sp, unchanged), `COMET:lnNumSP`, `COMET:lnRankSP`,
/// `COMET:IonFrac`.
///
/// Per identification, `second` is the XCorr (`MS:1002252`) of the second
/// hit (`0` with fewer than two hits) and `worst` that of the last hit. Then
/// per hit, each value is written only when the hit does not already carry
/// it, except `COMET:lnExpect`, which is always overwritten:
///
/// - `COMET:deltaCn = (xcorr - second) / max(1, xcorr)`
/// - `COMET:deltaLCn = (xcorr - worst) / max(1, xcorr)`
/// - `COMET:lnExpect = ln(MS:1002257)` (the E-value)
/// - `COMET:lnNumSP = ln(max(1, num_matched_peptides))` when
///   `num_matched_peptides` exists, else the **unlogged** `MS:1002255`
///   (the source's fallback, marked "TODO: remove?")
/// - `COMET:lnRankSP = ln(max(1, MS:1002256))`
/// - `COMET:IonFrac = MS:1002258 / MS:1002259` (matched / total ions)
///
/// A logarithm of zero or a division by zero stores a non-finite float,
/// exactly as the source; such a value fails [`MetaValue::validate`], so a
/// writer refuses it rather than emitting it.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a meta value the formulas need is missing
/// or not numeric (the source throws `Exception::ConversionError`). On
/// error `peptide_ids` and `feature_set` are unchanged.
pub fn add_comet_features(
    peptide_ids: &mut [PeptideIdentification],
    feature_set: &mut Vec<String>,
) -> Result<()> {
    let mut updated: Vec<Vec<MetaInfo>> = Vec::with_capacity(peptide_ids.len());
    for identification in peptide_ids.iter() {
        let mut worst = 0.0;
        let mut second = 0.0;
        for (index, hit) in identification.hits.iter().enumerate() {
            let xcorr = numeric(&hit.metadata, "MS:1002252")?;
            worst = xcorr;
            if index == 1 {
                second = xcorr;
            }
        }
        let mut metas = Vec::with_capacity(identification.hits.len());
        for hit in &identification.hits {
            let source = &hit.metadata;
            let mut meta = source.clone();
            let xcorr = numeric(source, "MS:1002252")?;
            if !source.contains_key("COMET:deltaCn") {
                meta.insert(
                    "COMET:deltaCn".into(),
                    source_float((xcorr - second) / xcorr.max(1.0)),
                );
            }
            if !source.contains_key("COMET:deltaLCn") {
                meta.insert(
                    "COMET:deltaLCn".into(),
                    source_float((xcorr - worst) / xcorr.max(1.0)),
                );
            }
            let ln_expect = numeric(source, "MS:1002257")?.ln();
            meta.insert("COMET:lnExpect".into(), source_float(ln_expect));
            if !source.contains_key("COMET:lnNumSP") {
                let ln_num_sp = if source.contains_key("num_matched_peptides") {
                    numeric(source, "num_matched_peptides")?.max(1.0).ln()
                } else {
                    numeric(source, "MS:1002255")?
                };
                meta.insert("COMET:lnNumSP".into(), source_float(ln_num_sp));
            }
            if !source.contains_key("COMET:lnRankSP") {
                let ln_rank_sp = numeric(source, "MS:1002256")?.max(1.0).ln();
                meta.insert("COMET:lnRankSP".into(), source_float(ln_rank_sp));
            }
            if !source.contains_key("COMET:IonFrac") {
                let matched = numeric(source, "MS:1002258")?;
                let total = numeric(source, "MS:1002259")?;
                meta.insert("COMET:IonFrac".into(), source_float(matched / total));
            }
            metas.push(meta);
        }
        updated.push(metas);
    }
    for (identification, metas) in peptide_ids.iter_mut().zip(updated) {
        for (hit, meta) in identification.hits.iter_mut().zip(metas) {
            hit.metadata = meta;
        }
    }
    feature_set.extend(
        [
            ISOTOPE_ERROR,
            "COMET:deltaCn",
            "COMET:deltaLCn",
            "COMET:lnExpect",
            "MS:1002252",
            "MS:1002255",
            "COMET:lnNumSP",
            "COMET:lnRankSP",
            "COMET:IonFrac",
        ]
        .iter()
        .map(|s| (*s).to_owned()),
    );
    Ok(())
}

/// Registers and computes the Mascot Percolator features.
///
/// Registers `MS:1001171` (Mascot score, unchanged), `MASCOT:delta_score`
/// and `MASCOT:hasMod`, in that order. Each identification's hits are sorted
/// by score (stable, in the identification's score direction), then
/// `MASCOT:delta_score` is assigned by [`assign_delta_score`] from
/// `MS:1001171`, and `MASCOT:hasMod` is `1` for a modified sequence and `0`
/// otherwise (an integer, as the source's `bool` argument promotes to the
/// `int` `DataValue` constructor).
///
/// The source documents features 1–9 of Mascot Percolator's basic set
/// (mass, charge, mScore, dScore, deltaM, absDeltaM, isoDeltaM, uniquePeps,
/// mc) and 10–18 of its extended set (totInt, intMatchedTot,
/// relIntMatchedTot, binom, fragMassError, absFragMassError,
/// fracIonsMatched, seqCov, intMatched). Only mScore and dScore are
/// computed here; the generic ones are added by the `.pin` writer and the
/// extended set is not available from the Mascot adapter.
///
/// **Divergence:** the source computes both values on a *copy* of each
/// identification's hits and discards it, so only the sort survives
/// (defect candidate, see `docs/PERCOLATOR_SUPPORT.md`). This function
/// writes the values, since the registered features would otherwise be
/// missing on every hit; [`add_mascot_features_source`] reproduces the
/// source exactly.
///
/// # Errors
///
/// [`Error::InvalidValue`] when an identification fails validation (the
/// sort needs finite scores) or a hit lacks a numeric `MS:1001171` (the
/// source throws `Exception::ConversionError`). On error nothing changes.
pub fn add_mascot_features(
    peptide_ids: &mut [PeptideIdentification],
    feature_set: &mut Vec<String>,
) -> Result<()> {
    let mut updated = peptide_ids.to_vec();
    for identification in &mut updated {
        identification.sort()?;
        assign_delta_score(&mut identification.hits, "MS:1001171", "MASCOT:delta_score")?;
        for hit in &mut identification.hits {
            let has_mod = i64::from(hit.sequence.is_modified());
            hit.metadata
                .insert("MASCOT:hasMod".into(), MetaValue::from(has_mod));
        }
    }
    peptide_ids.clone_from_slice(&updated);
    register_mascot(feature_set);
    Ok(())
}

/// The source's `addMASCOTFeatures` exactly: registers the three Mascot
/// features and sorts each identification's hits by score, but writes no
/// meta value, because the source computes them on a discarded copy.
///
/// The discarded computation still converts `MS:1001171` on every hit, so
/// the source throws where that is missing or not numeric; that check is
/// kept.
///
/// # Errors
///
/// As [`add_mascot_features`]. On error nothing changes.
pub fn add_mascot_features_source(
    peptide_ids: &mut [PeptideIdentification],
    feature_set: &mut Vec<String>,
) -> Result<()> {
    let mut updated = peptide_ids.to_vec();
    for identification in &mut updated {
        identification.sort()?;
        let mut discarded = identification.hits.clone();
        assign_delta_score(&mut discarded, "MS:1001171", "MASCOT:delta_score")?;
    }
    peptide_ids.clone_from_slice(&updated);
    register_mascot(feature_set);
    Ok(())
}

fn register_mascot(feature_set: &mut Vec<String>) {
    feature_set.extend(
        ["MS:1001171", "MASCOT:delta_score", "MASCOT:hasMod"]
            .iter()
            .map(|s| (*s).to_owned()),
    );
}

/// Merges one protein identification run into another.
///
/// Takes over the search parameters of the first run, sets the search
/// engine to `"multiple"` once runs from differing search engines are
/// combined, and unions the protein hits by accession. Peptide hits
/// reference proteins by accession, so they need no update.
///
/// Uses the current local time for a newly created run; see
/// [`merge_multise_protein_ids_at`] for the details and errors.
pub fn merge_multise_protein_ids(
    all_protein_ids: &mut Vec<ProteinIdentification>,
    new_protein_ids: &mut [ProteinIdentification],
) -> Result<()> {
    merge_multise_protein_ids_at(
        all_protein_ids,
        new_protein_ids,
        chrono::Local::now().naive_local(),
    )
}

/// [`merge_multise_protein_ids`] with an explicit creation time, for
/// reproducible output.
///
/// Only the first run of each list takes part, as in the source. Steps:
///
/// 1. When `all_protein_ids` is empty, a run is created with date/time
///    `now`, identifier `TopPerc_<yyyy-MM-dd>`, and the search engine and
///    search parameters of the new run. Otherwise, when the engines differ,
///    the target's engine becomes `"multiple"`.
/// 2. The new run's hits are sorted by accession (stable here; the source's
///    `std::sort` leaves equal accessions in unspecified order). If the
///    target has no hits, the two hit lists are **swapped**, so the new run
///    is left holding the target's former (empty) list, as in the source.
///    Otherwise the target's hits become `std::set_union` of both lists
///    under accession order: for equal accessions the target's hit is kept,
///    and the new run keeps its sorted hits. The union assumes the target's
///    hits are sorted by accession, which holds after a first merge; for an
///    unsorted target it reproduces the source algorithm step by step.
/// 3. When the target's search-parameter metadata lacks `SE:<engine>`, the
///    new run's settings are recorded there: `SE:<engine>` = engine version,
///    and `<engine>:db`, `:db_version`, `:taxonomy`, `:charges`,
///    `:fixed_modifications` and `:variable_modifications` (comma-joined),
///    `:missed_cleavages`, `:fragment_mass_tolerance`,
///    `:fragment_mass_tolerance_unit` (`ppm`/`Da`),
///    `:precursor_mass_tolerance`, `:precursor_mass_tolerance_unit`,
///    `:digestion_enzyme` and `:enzyme_term_specificity` (the source's
///    `EnzymaticDigestion::Specificity` integer: none 0, semi 1, full 2,
///    unknown 3).
/// 4. The new run's primary MS run paths are appended to the target's.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `new_protein_ids` is empty (undefined
/// behaviour in the source) or a tolerance is not finite. On error both
/// lists are unchanged.
pub fn merge_multise_protein_ids_at(
    all_protein_ids: &mut Vec<ProteinIdentification>,
    new_protein_ids: &mut [ProteinIdentification],
    now: chrono::NaiveDateTime,
) -> Result<()> {
    let new_run = new_protein_ids
        .first()
        .ok_or_else(|| invalid("no protein identification run to merge"))?;
    let engine = new_run.search_engine.clone();
    let mut target = match all_protein_ids.first() {
        Some(existing) => {
            let mut run = existing.clone();
            if run.search_engine != engine {
                run.search_engine = "multiple".into();
            }
            run
        }
        None => ProteinIdentification {
            date_time: Some(now.format("%Y-%m-%dT%H:%M:%S").to_string()),
            identifier: format!("TopPerc_{}", now.format("%Y-%m-%d")),
            search_engine: engine.clone(),
            search_parameters: new_run.search_parameters.clone(),
            ..ProteinIdentification::default()
        },
    };
    let mut new_hits = new_run.hits.clone();
    new_hits.sort_by(accession_order);
    if target.hits.is_empty() {
        std::mem::swap(&mut target.hits, &mut new_hits);
    } else {
        target.hits = set_union_by_accession(&target.hits, &new_hits);
    }
    let marker = format!("SE:{engine}");
    if !target.search_parameters.metadata.contains_key(&marker) {
        let sp = &new_run.search_parameters;
        let meta = &mut target.search_parameters.metadata;
        let text = |value: &str| MetaValue::from(value);
        meta.insert(marker, text(&new_run.search_engine_version));
        meta.insert(format!("{engine}:db"), text(&sp.database));
        meta.insert(format!("{engine}:db_version"), text(&sp.database_version));
        meta.insert(format!("{engine}:taxonomy"), text(&sp.taxonomy));
        meta.insert(format!("{engine}:charges"), text(&sp.charges));
        meta.insert(
            format!("{engine}:fixed_modifications"),
            text(&sp.fixed_modifications.join(",")),
        );
        meta.insert(
            format!("{engine}:variable_modifications"),
            text(&sp.variable_modifications.join(",")),
        );
        meta.insert(
            format!("{engine}:missed_cleavages"),
            MetaValue::from(sp.missed_cleavages),
        );
        let (fragment, fragment_unit) = tolerance_parts(sp.fragment_tolerance);
        let (precursor, precursor_unit) = tolerance_parts(sp.precursor_tolerance);
        meta.insert(
            format!("{engine}:fragment_mass_tolerance"),
            MetaValue::try_from(fragment)?,
        );
        meta.insert(
            format!("{engine}:fragment_mass_tolerance_unit"),
            text(fragment_unit),
        );
        meta.insert(
            format!("{engine}:precursor_mass_tolerance"),
            MetaValue::try_from(precursor)?,
        );
        meta.insert(
            format!("{engine}:precursor_mass_tolerance_unit"),
            text(precursor_unit),
        );
        meta.insert(
            format!("{engine}:digestion_enzyme"),
            text(&sp.digestion_enzyme),
        );
        meta.insert(
            format!("{engine}:enzyme_term_specificity"),
            MetaValue::from(specificity_code(sp.enzyme_specificity)),
        );
    }
    target
        .primary_ms_run_paths
        .extend(new_run.primary_ms_run_paths.iter().cloned());

    // Commit.
    new_protein_ids[0].hits = new_hits;
    if all_protein_ids.is_empty() {
        all_protein_ids.push(target);
    } else {
        all_protein_ids[0] = target;
    }
    Ok(())
}

/// Checks requested extra Percolator features against the PSMs and removes
/// those that are unavailable.
///
/// A requested feature survives only if every hit in `psms` carries it as a
/// meta value. Surviving features keep their order. The removed names are
/// returned in their original order; the source logs a warning for each
/// ("A extra_feature requested (…) was not available - removed.").
pub fn check_extra_features(psms: &[PeptideHit], extra_features: &mut Vec<String>) -> Vec<String> {
    let (kept, removed): (Vec<String>, Vec<String>) = extra_features
        .drain(..)
        .partition(|feature| psms.iter().all(|hit| hit.metadata.contains_key(feature)));
    *extra_features = kept;
    removed
}

/// Rescales a fragment feature to penalize features calculated from few
/// ions, adapted from Percolator's `MSGFPlusReader` converter.
///
/// `value * (1 + 7)^2 / (1 + min(matched_main_ions, 7))^2`. Protected and
/// unused in the source; public here so the formula stays reachable.
///
/// A `matched_main_ions` of `-1` divides by zero, as in the source, and
/// yields a non-finite result.
pub fn rescale_fragment_feature(feature_value: f64, matched_main_ions: i32) -> f64 {
    const LIMIT: i64 = 7;
    let numerator = (1 + LIMIT) * (1 + LIMIT);
    let capped = 1 + i64::from(matched_main_ions).min(LIMIT);
    let denominator = capped * capped;
    feature_value * (numerator as f64 / denominator as f64)
}

/// Assigns the frequently occurring delta score feature.
///
/// For each hit but the last, `output_ref = score_ref(hit) -
/// score_ref(next hit)`; the last (or only) hit gets `0.0`. An empty slice is
/// left alone. Protected in the source; public here, since callers need it
/// for engines this helper does not cover.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a hit lacks `score_ref` or it is not an
/// integer or a float (the source's `double(DataValue)` throws
/// `Exception::ConversionError`; unlike [`add_comet_features`]' text
/// conversion, a numeric string is refused here too). On error `hits` is
/// unchanged.
pub fn assign_delta_score(
    hits: &mut [PeptideHit],
    score_ref: &str,
    output_ref: &str,
) -> Result<()> {
    if hits.is_empty() {
        return Ok(());
    }
    let scores = hits
        .iter()
        .map(|hit| typed_numeric(&hit.metadata, score_ref))
        .collect::<Result<Vec<_>>>()?;
    let count = hits.len();
    for (index, hit) in hits.iter_mut().enumerate() {
        let value = if index + 1 < count {
            scores[index] - scores[index + 1]
        } else {
            0.0
        };
        hit.metadata
            .insert(output_ref.to_owned(), source_float(value));
    }
    Ok(())
}

/// Accession order of protein hits, the source's `lq_ProteinHit`
/// comparator (byte-wise `std::string` less-than).
pub fn accession_order(a: &ProteinHit, b: &ProteinHit) -> Ordering {
    a.accession.as_bytes().cmp(b.accession.as_bytes())
}

/// `std::set_union` exactly as libstdc++ implements it, under
/// [`accession_order`].
fn set_union_by_accession(first: &[ProteinHit], second: &[ProteinHit]) -> Vec<ProteinHit> {
    let mut out = Vec::with_capacity(first.len() + second.len());
    let (mut i, mut j) = (0, 0);
    while i < first.len() && j < second.len() {
        match accession_order(&first[i], &second[j]) {
            Ordering::Less => {
                out.push(first[i].clone());
                i += 1;
            }
            Ordering::Greater => {
                out.push(second[j].clone());
                j += 1;
            }
            Ordering::Equal => {
                out.push(first[i].clone());
                i += 1;
                j += 1;
            }
        }
    }
    out.extend(first[i..].iter().cloned());
    out.extend(second[j..].iter().cloned());
    out
}

fn tolerance_parts(tolerance: Tolerance) -> (f64, &'static str) {
    match tolerance {
        Tolerance::Absolute(value) => (value, "Da"),
        Tolerance::Ppm(value) => (value, "ppm"),
    }
}

/// The source `EnzymaticDigestion::Specificity` integer.
fn specificity_code(value: EnzymeTermSpecificity) -> i64 {
    match value {
        EnzymeTermSpecificity::None => 0,
        EnzymeTermSpecificity::Semi => 1,
        EnzymeTermSpecificity::Full => 2,
        EnzymeTermSpecificity::Unknown => 3,
    }
}

/// `DataValue::toString()` is non-empty: present, and not an empty string.
fn text_non_empty(value: Option<&MetaValue>) -> bool {
    match value.map(MetaValue::data) {
        None | Some(MetaValueData::Empty) => false,
        Some(MetaValueData::String(text)) => !text.is_empty(),
        Some(_) => true,
    }
}

/// The source's `StringUtils::toDouble(getMetaValue(key).toString())`:
/// numbers pass through, text is parsed (surrounding whitespace and a
/// leading `+` allowed), anything else is an error.
fn numeric(meta: &MetaInfo, key: &str) -> Result<f64> {
    let value = meta
        .get(key)
        .ok_or_else(|| invalid(format!("meta value '{key}' is missing")))?;
    match value.data() {
        MetaValueData::Integer(number) => Ok(*number as f64),
        MetaValueData::Float(number) => Ok(*number),
        MetaValueData::String(text) => crate::concept::string_conversion::to_double(text)
            .map_err(|_| invalid(format!("meta value '{key}' is not numeric: '{text}'"))),
        _ => Err(invalid(format!("meta value '{key}' is not numeric"))),
    }
}

/// The source's `double(DataValue)`: only integers and floats convert.
fn typed_numeric(meta: &MetaInfo, key: &str) -> Result<f64> {
    match meta.get(key).map(MetaValue::data) {
        Some(MetaValueData::Integer(number)) => Ok(*number as f64),
        Some(MetaValueData::Float(number)) => Ok(*number),
        Some(_) => Err(invalid(format!("meta value '{key}' is not numeric"))),
        None => Err(invalid(format!("meta value '{key}' is missing"))),
    }
}

/// A float stored as the source stores it, non-finite included.
fn source_float(value: f64) -> MetaValue {
    MetaValue::source_float(value)
}

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Median normalisation of the maps of a consensus map.
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h`.
//! The accession/description filter the source shares with
//! `ConsensusMapNormalizerAlgorithmThreshold` (its `friend`) is
//! [`IdentificationFilter`] here. See `docs/CONSENSUS_MAP_NORMALIZER_SUPPORT.md`.
//!
//! The source class is a stateless bag of static functions whose copy
//! constructor and assignment are suppressed; here it is a unit struct with
//! associated functions, and nothing is lost because it holds no state.

use crate::concept::boost_regex::BoostRegex;
use crate::error::{Error, Result};
use crate::kernel::features::{ConsensusFeature, ConsensusMap, FeatureHandle};
use crate::math::statistic_functions;

/// Largest number of feature handles, summed over all consensus features, that
/// one normalisation call accepts. Checked before anything is allocated. The
/// source has no ceiling.
pub const MAX_HANDLES: usize = 1 << 27;

/// Whether maps are brought to a common median by scaling or by shifting
/// (source `ConsensusMapNormalizerAlgorithmMedian::NormalizationMethod`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NormalizationMethod {
    /// Multiply by `median(reference) / median(map)` (source `NM_SCALE`).
    Scale,
    /// Add `median(largest) - median(map)` (source `NM_SHIFT`). The source logs
    /// that shifting is not recommended for regular log-normal MS data; use it
    /// only when the intensities are known to be on an additive scale.
    Shift,
}

/// The accession and description regular expressions that select which
/// consensus features contribute to normalisation coefficients (source
/// `passesFilters_`, which recompiles both expressions for every feature).
///
/// An empty expression, or one that matches the empty string, selects every
/// feature, including features without any identification.
#[derive(Clone, Debug)]
pub struct IdentificationFilter {
    /// Compiled accession expression; `None` when it is empty or matches the
    /// empty string, which the source treats as "no filter".
    accession: Option<BoostRegex>,
    /// Compiled description expression, `None` under the same rule.
    description: Option<BoostRegex>,
}

impl IdentificationFilter {
    /// Compile the two filters. Both use Boost perl syntax with default
    /// options, as the source's `boost::regex(acc_filter)` does.
    ///
    /// # Errors
    ///
    /// Returns the [`BoostRegex::new`] error for a malformed expression. The
    /// source constructs the same `boost::regex` and lets its `regex_error`
    /// escape uncaught.
    pub fn new(accession_filter: &str, description_filter: &str) -> Result<Self> {
        Ok(Self {
            accession: compile_nontrivial(accession_filter)?,
            description: compile_nontrivial(description_filter)?,
        })
    }

    /// A filter that passes every feature, as two empty expressions do.
    pub fn pass_all() -> Self {
        Self {
            accession: None,
            description: None,
        }
    }

    /// Whether `feature` of `map` passes both filters (source `passesFilters_`).
    ///
    /// When both filters are trivial every feature passes, even one without an
    /// identification. Otherwise the feature passes when some protein accession
    /// of some peptide hit (empty accessions are ignored, as the source's
    /// `extractProteinAccessionsSet` does) matches the accession filter and,
    /// unless the description filter is trivial, the description of the first
    /// protein hit with that accession in some protein identification of `map`
    /// matches the description filter. Matching is a search, not a full match.
    ///
    /// Source `regex_search` receives `c_str()`, so text after an embedded NUL
    /// byte is invisible to it; this port searches the same prefix.
    ///
    /// # Errors
    ///
    /// Returns the [`BoostRegex::is_search_match`] error when the regex engine
    /// exceeds its backtracking budget, where Boost would throw.
    pub fn passes(&self, feature: &ConsensusFeature, map: &ConsensusMap) -> Result<bool> {
        if self.accession.is_none() && self.description.is_none() {
            return Ok(true);
        }
        for identification in &feature.peptide_identifications {
            for hit in &identification.hits {
                for accession in hit.protein_accessions() {
                    if let Some(regex) = &self.accession {
                        if !regex.is_search_match(c_str(accession))? {
                            continue;
                        }
                    }
                    let Some(description) = &self.description else {
                        return Ok(true);
                    };
                    for protein in &map.protein_identifications {
                        if let Some(protein_hit) = protein.find_hit(accession) {
                            if description.is_search_match(c_str(&protein_hit.description()))? {
                                return Ok(true);
                            }
                        }
                    }
                }
            }
        }
        Ok(false)
    }
}

/// `None` for an expression the source treats as absent: empty, or matching
/// the empty string.
fn compile_nontrivial(pattern: &str) -> Result<Option<BoostRegex>> {
    if pattern.is_empty() {
        return Ok(None);
    }
    let regex = BoostRegex::new(pattern)?;
    if regex.is_search_match(b"")? {
        Ok(None)
    } else {
        Ok(Some(regex))
    }
}

/// The bytes a C string built from `text` would expose: up to the first NUL.
fn c_str(text: &str) -> &[u8] {
    let bytes = text.as_bytes();
    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
    &bytes[..end]
}

/// The medians of all maps and the reference map (the result of source
/// `computeMedians`, whose return value is the reference index and whose
/// out-parameter is the median vector).
#[derive(Clone, Debug, PartialEq)]
pub struct Medians {
    /// Index of the map whose column header records the largest `size`; the
    /// first such map on ties, and map 0 when every size is 0.
    pub reference_map: usize,
    /// One median per map, or all `1.0` when some map had no passing handle.
    pub medians: Vec<f64>,
    /// Number of consensus features that passed the filters (the source logs
    /// it as "Using N/M consensus features").
    pub passing_features: usize,
    /// `false` when some map had no handle in a passing feature. The source
    /// then logs a warning, sets every median to 1.0 and returns 0, so the
    /// normalisation is a multiply-by-one; `reference_map` is 0 in that case,
    /// as the source's return value is.
    pub enough_features: bool,
}

/// Median normalisation (source `ConsensusMapNormalizerAlgorithmMedian`).
#[derive(Clone, Copy, Debug, Default)]
pub struct ConsensusMapNormalizerAlgorithmMedian;

impl ConsensusMapNormalizerAlgorithmMedian {
    /// Compute the median handle intensity of every map and the reference map
    /// (source `computeMedians`).
    ///
    /// The maps are the column headers `0..n`, `n` being the number of column
    /// headers. Only handles of features that pass `filter` contribute.
    ///
    /// The reference is the argmax of the stored column-header `size`, **not**
    /// of the number of handles present; with every size at 0 (as some linking
    /// paths leave them) it silently collapses to map 0. That is the source
    /// behaviour, pinned by its class test, and kept.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a column header in `0..n` is
    /// missing (source `Exception::ElementNotFound`), when a handle names a map
    /// index outside `0..n` (undefined behaviour in the source, which indexes a
    /// vector with it), when the handle count exceeds [`MAX_HANDLES`], and the
    /// [`IdentificationFilter::passes`] error.
    pub fn compute_medians(map: &ConsensusMap, filter: &IdentificationFilter) -> Result<Medians> {
        let number_of_maps = check_maps(map)?;
        let reference_map = largest_map(map, number_of_maps);
        let mut intensities = vec![Vec::<f64>::new(); number_of_maps];
        let mut passing_features = 0;
        for feature in &map.features {
            if !filter.passes(feature, map)? {
                continue;
            }
            passing_features += 1;
            for handle in feature.handles() {
                intensities[map_slot(handle)].push(f64::from(handle.intensity));
            }
        }
        if intensities.iter().any(Vec::is_empty) {
            return Ok(Medians {
                reference_map: 0,
                medians: vec![1.0; number_of_maps],
                passing_features,
                enough_features: false,
            });
        }
        let medians = intensities
            .iter_mut()
            .map(|values| statistic_functions::median(values))
            .collect::<Result<Vec<_>>>()?;
        Ok(Medians {
            reference_map,
            medians,
            passing_features,
            enough_features: true,
        })
    }

    /// Normalise the handle intensities of every map to a common median
    /// (source `normalizeMaps`).
    ///
    /// [`NormalizationMethod::Scale`] multiplies each handle of map `j` by
    /// `medians[reference] / medians[j]`, evaluated as
    /// `(intensity * medians[reference]) / medians[j]` in `f64` and stored as
    /// `f32`, as the source's `float` intensity is.
    ///
    /// [`NormalizationMethod::Shift`] ignores the reference and shifts towards
    /// the largest median instead: it adds `medians[max] - medians[j]`. The
    /// source searches the maximum starting from
    /// `std::numeric_limits<double>::min()`, the smallest *positive* double, so
    /// when no median exceeds it the target is map 0; this port keeps that.
    ///
    /// Only the handles are rewritten; the consensus features' own intensities
    /// are left as they are, as in the source.
    ///
    /// Returns the [`Medians`] used. When they report
    /// `enough_features == false`, which the source only logs, every median is
    /// 1.0: scaling then leaves every intensity unchanged, and shifting adds
    /// `1.0 - 1.0` in `f64`, which is exact except that an intensity below the
    /// rounding unit of `1.0` is absorbed, as in the source.
    ///
    /// # Errors
    ///
    /// As [`Self::compute_medians`]. Also [`Error::InvalidValue`] when a new
    /// intensity is not finite, which a zero median produces; the source
    /// stores the infinity or NaN. On any error the map is unchanged.
    pub fn normalize_maps(
        map: &mut ConsensusMap,
        method: NormalizationMethod,
        filter: &IdentificationFilter,
    ) -> Result<Medians> {
        let medians = Self::compute_medians(map, filter)?;
        let values = &medians.medians;
        let shift_target = {
            let mut max_median = f64::MIN_POSITIVE;
            let mut index = 0;
            for (i, &median) in values.iter().enumerate() {
                if median > max_median {
                    max_median = median;
                    index = i;
                }
            }
            index
        };
        rewrite_intensities(map, |handle| {
            let j = map_slot(handle);
            let intensity = f64::from(handle.intensity);
            match method {
                NormalizationMethod::Scale => {
                    // The source's `medians[index_of_largest_map]` with an
                    // empty median vector is only reached without handles.
                    intensity * values[medians.reference_map] / values[j]
                }
                NormalizationMethod::Shift => intensity + values[shift_target] - values[j],
            }
        })?;
        Ok(medians)
    }
}

/// Check the column headers `0..n` exist and every handle names one of them,
/// and bound the handle count. Returns `n`.
pub(crate) fn check_maps(map: &ConsensusMap) -> Result<usize> {
    let number_of_maps = map.column_headers.len();
    for index in 0..number_of_maps {
        if !map.column_headers.contains_key(&(index as u64)) {
            return Err(Error::InvalidValue(format!(
                "consensus map has {number_of_maps} column headers but none for map index {index}"
            )));
        }
    }
    let mut handles = 0usize;
    for feature in &map.features {
        handles = handles.saturating_add(feature.len());
        if handles > MAX_HANDLES {
            return Err(Error::InvalidValue(format!(
                "consensus map has more than {MAX_HANDLES} feature handles"
            )));
        }
        for handle in feature.handles() {
            if handle.map_index >= number_of_maps as u64 {
                return Err(Error::InvalidValue(format!(
                    "feature handle names map index {} but the consensus map has {number_of_maps} column headers",
                    handle.map_index
                )));
            }
        }
    }
    Ok(number_of_maps)
}

/// Index of the map whose column header has the largest `size`, first on
/// ties, starting from map 0 (source: `find(0)` then a strict `>` scan).
/// Callers have run [`check_maps`].
pub(crate) fn largest_map(map: &ConsensusMap, number_of_maps: usize) -> usize {
    let size = |index: usize| {
        map.column_headers
            .get(&(index as u64))
            .map_or(0, |header| header.size)
    };
    let mut best = 0;
    for index in 0..number_of_maps {
        if size(index) > size(best) {
            best = index;
        }
    }
    best
}

/// The handle's map index as a vector slot. Callers have run [`check_maps`],
/// which bounds it by the number of column headers.
pub(crate) fn map_slot(handle: &FeatureHandle) -> usize {
    handle.map_index as usize
}

/// Replace every handle intensity by `new_intensity(handle)`, narrowed to
/// `f32`. Transactional: every value is computed and checked before the map is
/// touched.
pub(crate) fn rewrite_intensities(
    map: &mut ConsensusMap,
    mut new_intensity: impl FnMut(&FeatureHandle) -> f64,
) -> Result<()> {
    let mut replaced = Vec::with_capacity(map.features.len());
    for feature in &map.features {
        let mut handles = feature.handles().to_vec();
        for handle in &mut handles {
            let value = new_intensity(handle) as f32;
            if !value.is_finite() {
                return Err(Error::InvalidValue(format!(
                    "normalisation produced a non-finite intensity for handle ({}, {})",
                    handle.map_index, handle.unique_id
                )));
            }
            handle.intensity = value;
            handle.validate()?;
        }
        replaced.push(handles);
    }
    for (feature, handles) in map.features.iter_mut().zip(replaced) {
        feature.set_handles(handles)?;
    }
    Ok(())
}

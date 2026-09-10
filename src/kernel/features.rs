// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Native feature and consensus containers with attached identification records.

use super::NumericRange;
use super::geometry::{BoundingBox2D, ConvexHull2D, Point2D};
use crate::format::FileType;
use crate::identification::{PeptideIdentification, ProteinIdentification};
use crate::metadata::{DataProcessing, MetaInfo, MetaValue, validate_meta};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::ops::{Deref, DerefMut};

/// Common measured values. ID zero denotes an unassigned ID, as in OpenMS.
/// Quality and intensity may be signed; width must be nonnegative when checked.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BaseFeature {
    pub rt: f64,
    pub mz: f64,
    pub intensity: f32,
    pub quality: f32,
    pub charge: i32,
    pub width: f32,
    pub unique_id: u64,
    pub peptide_identifications: Vec<PeptideIdentification>,
    pub metadata: MetaInfo,
}

impl BaseFeature {
    pub fn new(rt: f64, mz: f64, intensity: f32) -> Self {
        Self {
            rt,
            mz,
            intensity,
            ..Self::default()
        }
    }

    pub fn validate(&self) -> Result<()> {
        validate_values(self.rt, self.mz, self.intensity, self.width)?;
        finite(f64::from(self.quality), "feature quality")?;
        validate_meta(&self.metadata)?;
        for identification in &self.peptide_identifications {
            identification.validate()?;
        }
        Ok(())
    }

    /// Store width and the source featureXML FWHM metadata bridge together.
    pub fn set_width(&mut self, width: f32) -> Result<()> {
        finite(f64::from(width), "feature width")?;
        if width < 0.0 {
            return Err(Error::InvalidValue(
                "feature width must be nonnegative".into(),
            ));
        }
        self.metadata
            .insert("FWHM".into(), MetaValue::try_from(f64::from(width))?);
        self.width = width;
        Ok(())
    }

    pub const fn position(&self) -> Point2D {
        Point2D::new(self.rt, self.mz)
    }
}

/// A feature with per-dimension quality, mass-trace hulls and subordinate features.
/// Common fields are also accessible through `Deref`, e.g. `feature.rt`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Feature {
    pub base: BaseFeature,
    pub quality_rt: f32,
    pub quality_mz: f32,
    pub convex_hulls: Vec<ConvexHull2D>,
    pub subordinates: Vec<Feature>,
}

impl Deref for Feature {
    type Target = BaseFeature;
    fn deref(&self) -> &Self::Target {
        &self.base
    }
}
impl DerefMut for Feature {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.base
    }
}
impl From<BaseFeature> for Feature {
    fn from(base: BaseFeature) -> Self {
        Self {
            base,
            ..Self::default()
        }
    }
}

impl Feature {
    /// Checked operations reject deeper subordinate trees before recursive cloning.
    pub const MAX_SUBORDINATE_DEPTH: usize = 128;

    pub fn new(rt: f64, mz: f64, intensity: f32) -> Self {
        BaseFeature::new(rt, mz, intensity).into()
    }

    /// Iterative traversal avoids recursion during validation. Hull geometry is
    /// already validated by its constructors and has no mutable public internals.
    pub fn validate(&self) -> Result<()> {
        let mut pending = vec![(self, 0)];
        while let Some((feature, depth)) = pending.pop() {
            if depth > Self::MAX_SUBORDINATE_DEPTH {
                return Err(Error::InvalidValue(
                    "subordinate feature depth exceeds 128".into(),
                ));
            }
            feature.base.validate()?;
            finite(f64::from(feature.quality_rt), "RT quality")?;
            finite(f64::from(feature.quality_mz), "m/z quality")?;
            pending.extend(feature.subordinates.iter().map(|sub| (sub, depth + 1)));
        }
        Ok(())
    }

    /// One mass trace retains its exact hull; multiple traces use the rectangular
    /// bounding union, matching `Feature::getConvexHull` rather than a polygon union.
    pub fn convex_hull(&self) -> ConvexHull2D {
        if self.convex_hulls.len() == 1 {
            return self.convex_hulls[0].clone();
        }
        self.hull_bounding_box()
            .map_or_else(ConvexHull2D::new, ConvexHull2D::from_bounding_box)
    }

    pub fn hull_bounding_box(&self) -> Option<BoundingBox2D> {
        self.convex_hulls
            .iter()
            .filter_map(ConvexHull2D::bounding_box)
            .reduce(BoundingBox2D::union)
    }

    /// Membership in any individual mass-trace hull, excluding gaps between traces.
    pub fn encloses(&self, rt: f64, mz: f64) -> Result<bool> {
        let point = Point2D::new(rt, mz);
        point.validate()?;
        for hull in &self.convex_hulls {
            if hull.encloses(point)? {
                return Ok(true);
            }
        }
        Ok(false)
    }
}

/// An owned snapshot of a feature in a source map. Identity is the pair
/// `(map_index, unique_id)`; a consensus never stores that pair twice.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FeatureHandle {
    pub map_index: u64,
    pub unique_id: u64,
    pub rt: f64,
    pub mz: f64,
    pub intensity: f32,
    pub charge: i32,
    pub width: f32,
}

impl FeatureHandle {
    pub fn new(map_index: u64, feature: &BaseFeature) -> Self {
        Self {
            map_index,
            unique_id: feature.unique_id,
            rt: feature.rt,
            mz: feature.mz,
            intensity: feature.intensity,
            charge: feature.charge,
            width: feature.width,
        }
    }

    pub const fn key(&self) -> (u64, u64) {
        (self.map_index, self.unique_id)
    }
    pub fn validate(&self) -> Result<()> {
        validate_values(self.rt, self.mz, self.intensity, self.width)
    }
}

/// Independent dimension bounds. Empty containers have no ranges.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FeatureRanges {
    pub rt: Option<NumericRange>,
    pub mz: Option<NumericRange>,
    pub intensity: Option<NumericRange>,
}

impl FeatureRanges {
    fn add(&mut self, rt: f64, mz: f64, intensity: f32) {
        extend(&mut self.rt, rt);
        extend(&mut self.mz, mz);
        extend(&mut self.intensity, f64::from(intensity));
    }
}

/// Consensus summary and handles ordered by map index then unique ID.
/// Insertion does not implicitly recompute the summary.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ConsensusFeature {
    pub base: BaseFeature,
    handles: Vec<FeatureHandle>,
}

impl Deref for ConsensusFeature {
    type Target = BaseFeature;
    fn deref(&self) -> &Self::Target {
        &self.base
    }
}
impl DerefMut for ConsensusFeature {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.base
    }
}
impl From<BaseFeature> for ConsensusFeature {
    fn from(base: BaseFeature) -> Self {
        Self {
            base,
            handles: Vec::new(),
        }
    }
}

impl ConsensusFeature {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_feature(map_index: u64, feature: &BaseFeature) -> Result<Self> {
        feature.validate()?;
        let mut result = Self::from(feature.clone());
        result.insert(FeatureHandle::new(map_index, feature))?;
        Ok(result)
    }

    pub fn handles(&self) -> &[FeatureHandle] {
        &self.handles
    }
    pub fn len(&self) -> usize {
        self.handles.len()
    }
    pub fn is_empty(&self) -> bool {
        self.handles.is_empty()
    }
    /// Clears the handles only; measured summary and metadata are retained.
    pub fn clear(&mut self) {
        self.handles.clear();
    }

    pub fn validate(&self) -> Result<()> {
        self.base.validate()?;
        for handle in &self.handles {
            handle.validate()?;
        }
        Ok(())
    }

    /// Duplicate identities fail without changing existing data.
    pub fn insert(&mut self, handle: FeatureHandle) -> Result<()> {
        handle.validate()?;
        match self
            .handles
            .binary_search_by_key(&handle.key(), FeatureHandle::key)
        {
            Ok(_) => Err(Error::InvalidValue(format!(
                "duplicate consensus handle {:?}",
                handle.key()
            ))),
            Err(index) => {
                self.handles.insert(index, handle);
                Ok(())
            }
        }
    }

    /// Replaces all handles transactionally, rejecting duplicate identities.
    pub fn set_handles(&mut self, handles: Vec<FeatureHandle>) -> Result<()> {
        let mut handles = handles;
        for handle in &handles {
            handle.validate()?;
        }
        handles.sort_by_key(FeatureHandle::key);
        if handles
            .windows(2)
            .any(|pair| pair[0].key() == pair[1].key())
        {
            return Err(Error::InvalidValue(
                "duplicate consensus handle identity".into(),
            ));
        }
        self.handles = handles;
        Ok(())
    }

    /// Set union with another consensus, retaining this object's value when an
    /// identity exists in both, matching C++ insertion of a consensus feature.
    pub fn merge(&mut self, other: &Self) -> Result<()> {
        for handle in &other.handles {
            handle.validate()?;
        }
        let mut combined = self.handles.clone();
        combined.extend_from_slice(&other.handles);
        combined.sort_by_key(FeatureHandle::key);
        combined.dedup_by_key(|handle| handle.key());
        self.handles = combined;
        Ok(())
    }

    /// Bounds of handles only, excluding the stored consensus summary.
    pub fn handle_ranges(&self) -> FeatureRanges {
        let mut ranges = FeatureRanges::default();
        for handle in &self.handles {
            ranges.add(handle.rt, handle.mz, handle.intensity);
        }
        ranges
    }

    /// Arithmetic mean RT, m/z and intensity. Charge follows the source's
    /// running vote: most occurrences, then smaller absolute charge; unresolved
    /// sign ties retain the running winner in handle-key order.
    pub fn compute_consensus(&mut self) -> Result<()> {
        self.compute_mean(false)
    }

    /// Mean RT and intensity, minimum m/z, and the same charge vote as above.
    pub fn compute_monoisotopic_consensus(&mut self) -> Result<()> {
        self.compute_mean(true)
    }

    fn compute_mean(&mut self, monoisotopic: bool) -> Result<()> {
        self.require_handles()?;
        let mut rt = 0.0;
        let mut mz = if monoisotopic { f64::MAX } else { 0.0 };
        let mut intensity = 0.0;
        let mut occurrences = BTreeMap::<i32, usize>::new();
        let (mut charge, mut most) = (0i32, 0);
        for handle in &self.handles {
            rt += handle.rt;
            if monoisotopic {
                mz = mz.min(handle.mz);
            } else {
                mz += handle.mz;
            }
            intensity += f64::from(handle.intensity);
            let count = occurrences.entry(handle.charge).or_default();
            *count += 1;
            if *count > most
                || (*count == most && handle.charge.unsigned_abs() < charge.unsigned_abs())
            {
                charge = handle.charge;
                most = *count;
            }
        }
        let count = self.handles.len() as f64;
        self.assign_summary(
            rt / count,
            if monoisotopic { mz } else { mz / count },
            intensity / count,
            charge,
        )
    }

    /// Computes neutral mass in `mz`, mean RT and summed intensity. Source-map
    /// lookup uses unique IDs, ignoring handle map indices as in C++. The optional
    /// `dc_charge_adduct_mass` metadata value overrides `charge * proton_mass`.
    ///
    /// Unknown charge is rejected. Weighted mode requires nonnegative intensities
    /// and positive total intensity. All failures leave the summary unchanged.
    pub fn compute_decharge_consensus(
        &mut self,
        map: &FeatureMap,
        intensity_weighted: bool,
    ) -> Result<()> {
        self.require_handles()?;
        map.validate()?;
        let index = unique_index(map.features.iter().map(|feature| feature.unique_id))?;
        let intensity: f64 = self.handles.iter().map(|h| f64::from(h.intensity)).sum();
        if intensity_weighted
            && (intensity <= 0.0 || self.handles.iter().any(|h| h.intensity < 0.0))
        {
            return Err(Error::InvalidValue(
                "weighted consensus requires nonnegative intensities and positive total".into(),
            ));
        }
        let mut rt = 0.0;
        let mut mass = 0.0;
        for handle in &self.handles {
            if handle.charge == 0 {
                return Err(Error::InvalidValue(
                    "decharging requires a nonzero charge".into(),
                ));
            }
            let feature = index
                .get(&handle.unique_id)
                .map(|&i| &map.features[i])
                .ok_or_else(|| {
                    Error::InvalidValue(format!(
                        "source feature ID {} is absent or unassigned",
                        handle.unique_id
                    ))
                })?;
            let adduct = match feature.metadata.get("dc_charge_adduct_mass") {
                Some(value) => match value.as_str() {
                    Ok(text) => text
                        .parse::<f64>()
                        .map_err(|_| Error::InvalidValue("invalid dc_charge_adduct_mass".into()))?,
                    Err(_) => value.as_f64()?,
                },
                None => f64::from(handle.charge) * crate::chemistry::PROTON_MASS_U,
            };
            finite(adduct, "adduct mass")?;
            let weight = if intensity_weighted {
                f64::from(handle.intensity) / intensity
            } else {
                1.0 / self.len() as f64
            };
            rt += handle.rt * weight;
            mass += (handle.mz * f64::from(handle.charge.unsigned_abs()) - adduct) * weight;
        }
        self.assign_summary(rt, mass, intensity, 0)
    }

    fn require_handles(&self) -> Result<()> {
        if self.is_empty() {
            return Err(Error::InvalidValue(
                "cannot compute an empty consensus".into(),
            ));
        }
        for handle in &self.handles {
            handle.validate()?;
        }
        Ok(())
    }

    fn assign_summary(&mut self, rt: f64, mz: f64, intensity: f64, charge: i32) -> Result<()> {
        finite(rt, "consensus RT")?;
        finite(mz, "consensus m/z or mass")?;
        let intensity = intensity as f32;
        finite(f64::from(intensity), "consensus intensity")?;
        self.rt = rt;
        self.mz = mz;
        self.intensity = intensity;
        self.charge = charge;
        Ok(())
    }
}

/// Owned features; subordinate features stay attached to their parents.
#[derive(Clone, Debug, PartialEq)]
pub struct FeatureMap {
    pub features: Vec<Feature>,
    pub protein_identifications: Vec<ProteinIdentification>,
    pub unassigned_peptide_identifications: Vec<PeptideIdentification>,
    pub unique_id: u64,
    pub identifier: String,
    pub metadata: MetaInfo,
    pub data_processing: Vec<DataProcessing>,
    pub loaded_file_path: String,
    pub loaded_file_type: FileType,
}

impl Default for FeatureMap {
    fn default() -> Self {
        Self {
            features: Vec::new(),
            protein_identifications: Vec::new(),
            unassigned_peptide_identifications: Vec::new(),
            unique_id: 0,
            identifier: String::new(),
            metadata: MetaInfo::new(),
            data_processing: Vec::new(),
            loaded_file_path: String::new(),
            loaded_file_type: FileType::Unknown,
        }
    }
}

impl FeatureMap {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn from_features(features: Vec<Feature>) -> Self {
        Self {
            features,
            ..Self::default()
        }
    }
    pub fn len(&self) -> usize {
        self.features.len()
    }
    pub fn is_empty(&self) -> bool {
        self.features.is_empty()
    }

    pub fn validate(&self) -> Result<()> {
        for feature in &self.features {
            feature.validate()?;
        }
        unique_index(self.features.iter().map(|feature| feature.unique_id))?;
        validate_meta(&self.metadata)?;
        for processing in &self.data_processing {
            processing.validate()?;
        }
        for identification in &self.protein_identifications {
            identification.validate()?;
        }
        for identification in &self.unassigned_peptide_identifications {
            identification.validate()?;
        }
        Ok(())
    }

    /// Searches top-level features only. Zero is not a usable lookup key;
    /// duplicate assigned IDs are errors, not an arbitrary first match.
    pub fn unique_id_to_index(&self, unique_id: u64) -> Result<Option<usize>> {
        lookup_id(
            self.features.iter().map(|feature| feature.unique_id),
            unique_id,
        )
    }

    /// Top-level centroids and hulls contribute RT/m/z; only top-level feature
    /// intensities contribute intensity, matching `FeatureMap::updateRanges`.
    pub fn ranges(&self) -> Result<FeatureRanges> {
        self.validate()?;
        let mut ranges = FeatureRanges::default();
        for feature in &self.features {
            ranges.add(feature.rt, feature.mz, feature.intensity);
            if let Some(bounds) = feature.hull_bounding_box() {
                extend(&mut ranges.rt, bounds.rt_range().min);
                extend(&mut ranges.rt, bounds.rt_range().max);
                extend(&mut ranges.mz, bounds.mz_range().min);
                extend(&mut ranges.mz, bounds.mz_range().max);
            }
        }
        Ok(ranges)
    }

    /// Selects unique zero-based indices in caller order, transactionally.
    pub fn select(&mut self, indices: &[usize]) -> Result<()> {
        self.validate()?;
        select(&mut self.features, indices)
    }

    pub fn sort_by_position(&mut self) -> Result<()> {
        self.sort_by(|a, b| cmp_position(&a.base, &b.base))
    }
    pub fn sort_by_rt(&mut self) -> Result<()> {
        self.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap())
    }
    pub fn sort_by_mz(&mut self) -> Result<()> {
        self.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap())
    }
    pub fn sort_by_intensity(&mut self, reverse: bool) -> Result<()> {
        self.sort_by(|a, b| ordered(a.intensity, b.intensity, reverse))
    }
    pub fn sort_by_quality(&mut self, reverse: bool) -> Result<()> {
        self.sort_by(|a, b| ordered(a.quality, b.quality, reverse))
    }

    fn sort_by(
        &mut self,
        compare: impl FnMut(&Feature, &Feature) -> std::cmp::Ordering,
    ) -> Result<()> {
        self.validate()?;
        self.features.sort_by(compare);
        Ok(())
    }

    pub fn clear(&mut self, clear_metadata: bool) {
        if clear_metadata {
            *self = Self::default();
        } else {
            self.features.clear();
        }
    }
}

/// Description of a source-map column in a consensus map.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ColumnHeader {
    pub filename: String,
    pub label: String,
    pub size: usize,
    pub unique_id: u64,
    pub metadata: MetaInfo,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ConsensusMap {
    pub features: Vec<ConsensusFeature>,
    pub protein_identifications: Vec<ProteinIdentification>,
    pub unassigned_peptide_identifications: Vec<PeptideIdentification>,
    pub column_headers: BTreeMap<u64, ColumnHeader>,
    pub unique_id: u64,
    pub identifier: String,
    pub experiment_type: String,
    pub metadata: MetaInfo,
    pub data_processing: Vec<DataProcessing>,
    pub loaded_file_path: String,
    pub loaded_file_type: FileType,
}

impl Default for ConsensusMap {
    fn default() -> Self {
        Self {
            features: Vec::new(),
            column_headers: BTreeMap::new(),
            protein_identifications: Vec::new(),
            unassigned_peptide_identifications: Vec::new(),
            unique_id: 0,
            identifier: String::new(),
            experiment_type: "label-free".into(),
            metadata: MetaInfo::new(),
            data_processing: Vec::new(),
            loaded_file_path: String::new(),
            loaded_file_type: FileType::Unknown,
        }
    }
}

impl ConsensusMap {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn from_features(features: Vec<ConsensusFeature>) -> Self {
        Self {
            features,
            ..Self::default()
        }
    }
    pub fn len(&self) -> usize {
        self.features.len()
    }
    pub fn is_empty(&self) -> bool {
        self.features.is_empty()
    }

    /// Checks numeric values and assigned summary IDs. Column referential
    /// integrity is available separately through [`Self::validate_consistency`].
    pub fn validate(&self) -> Result<()> {
        if !matches!(
            self.experiment_type.as_str(),
            "label-free" | "labeled_MS1" | "labeled_MS2"
        ) {
            return Err(Error::InvalidValue(
                "unknown consensus experiment type".into(),
            ));
        }
        for header in self.column_headers.values() {
            validate_meta(&header.metadata)?;
        }
        for feature in &self.features {
            feature.validate()?;
        }
        unique_index(self.features.iter().map(|feature| feature.unique_id))?;
        validate_meta(&self.metadata)?;
        for processing in &self.data_processing {
            processing.validate()?;
        }
        for identification in &self.protein_identifications {
            identification.validate()?;
        }
        for identification in &self.unassigned_peptide_identifications {
            identification.validate()?;
        }
        Ok(())
    }

    /// Each handle must refer to a described map, and (filename,label) pairs
    /// must be unique. Header `size` is a count, not a bound on unique IDs.
    pub fn validate_consistency(&self) -> Result<()> {
        self.validate()?;
        let mut descriptions = BTreeSet::new();
        for header in self.column_headers.values() {
            if !descriptions.insert((&header.filename, &header.label)) {
                return Err(Error::InvalidValue(
                    "duplicate consensus column filename and label".into(),
                ));
            }
        }
        for feature in &self.features {
            for handle in feature.handles() {
                if !self.column_headers.contains_key(&handle.map_index) {
                    return Err(Error::InvalidValue(format!(
                        "consensus handle references missing map {}",
                        handle.map_index
                    )));
                }
            }
        }
        Ok(())
    }

    pub fn unique_id_to_index(&self, unique_id: u64) -> Result<Option<usize>> {
        lookup_id(
            self.features.iter().map(|feature| feature.unique_id),
            unique_id,
        )
    }

    /// Includes both summaries and all handles, as in `ConsensusMap::updateRanges`.
    pub fn ranges(&self) -> Result<FeatureRanges> {
        self.validate()?;
        let mut ranges = FeatureRanges::default();
        for feature in &self.features {
            ranges.add(feature.rt, feature.mz, feature.intensity);
            for handle in feature.handles() {
                ranges.add(handle.rt, handle.mz, handle.intensity);
            }
        }
        Ok(ranges)
    }

    pub fn select(&mut self, indices: &[usize]) -> Result<()> {
        self.validate()?;
        select(&mut self.features, indices)
    }
    pub fn sort_by_position(&mut self) -> Result<()> {
        self.sort_by(|a, b| cmp_position(&a.base, &b.base))
    }
    pub fn sort_by_rt(&mut self) -> Result<()> {
        self.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap())
    }
    pub fn sort_by_mz(&mut self) -> Result<()> {
        self.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap())
    }
    pub fn sort_by_intensity(&mut self, reverse: bool) -> Result<()> {
        self.sort_by(|a, b| ordered(a.intensity, b.intensity, reverse))
    }
    pub fn sort_by_quality(&mut self, reverse: bool) -> Result<()> {
        self.sort_by(|a, b| ordered(a.quality, b.quality, reverse))
    }
    /// Decreasing number of handles; ties retain input order.
    pub fn sort_by_size(&mut self) -> Result<()> {
        self.sort_by(|a, b| b.len().cmp(&a.len()))
    }
    /// Lexicographic order of complete handle identities, including unique IDs.
    pub fn sort_by_maps(&mut self) -> Result<()> {
        self.sort_by(|a, b| {
            a.handles
                .iter()
                .map(FeatureHandle::key)
                .cmp(b.handles.iter().map(FeatureHandle::key))
        })
    }

    fn sort_by(
        &mut self,
        compare: impl FnMut(&ConsensusFeature, &ConsensusFeature) -> std::cmp::Ordering,
    ) -> Result<()> {
        self.validate()?;
        self.features.sort_by(compare);
        Ok(())
    }

    pub fn clear(&mut self, clear_metadata: bool) {
        if clear_metadata {
            *self = Self::default();
        } else {
            self.features.clear();
        }
    }
}

fn finite(value: f64, name: &str) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::InvalidValue(format!("{name} must be finite")))
    }
}

fn validate_values(rt: f64, mz: f64, intensity: f32, width: f32) -> Result<()> {
    Point2D::new(rt, mz).validate()?;
    finite(f64::from(intensity), "feature intensity")?;
    finite(f64::from(width), "feature width")?;
    if width < 0.0 {
        return Err(Error::InvalidValue(
            "feature width must be nonnegative".into(),
        ));
    }
    Ok(())
}

fn extend(range: &mut Option<NumericRange>, value: f64) {
    match range {
        Some(range) => {
            range.min = range.min.min(value);
            range.max = range.max.max(value);
        }
        None => {
            *range = Some(NumericRange {
                min: value,
                max: value,
            })
        }
    }
}

fn unique_index(ids: impl Iterator<Item = u64>) -> Result<BTreeMap<u64, usize>> {
    let mut index = BTreeMap::new();
    for (i, id) in ids.enumerate() {
        if id != 0 && index.insert(id, i).is_some() {
            return Err(Error::InvalidValue(format!(
                "duplicate assigned feature ID {id}"
            )));
        }
    }
    Ok(index)
}

fn lookup_id(ids: impl Iterator<Item = u64>, unique_id: u64) -> Result<Option<usize>> {
    if unique_id == 0 {
        return Err(Error::InvalidValue(
            "zero is an unassigned feature ID".into(),
        ));
    }
    Ok(unique_index(ids)?.get(&unique_id).copied())
}

fn select<T: Clone>(values: &mut Vec<T>, indices: &[usize]) -> Result<()> {
    let mut seen = BTreeSet::new();
    for &index in indices {
        if index >= values.len() || !seen.insert(index) {
            return Err(Error::InvalidValue(
                "selection indices must be unique and in bounds".into(),
            ));
        }
    }
    *values = indices.iter().map(|&index| values[index].clone()).collect();
    Ok(())
}

fn cmp_position(a: &BaseFeature, b: &BaseFeature) -> std::cmp::Ordering {
    a.rt.partial_cmp(&b.rt)
        .unwrap()
        .then_with(|| a.mz.partial_cmp(&b.mz).unwrap())
}

fn ordered(a: f32, b: f32, reverse: bool) -> std::cmp::Ordering {
    let order = a.partial_cmp(&b).unwrap();
    if reverse { order.reverse() } else { order }
}

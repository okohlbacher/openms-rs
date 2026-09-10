// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Atomic application of retention-time transformations to native data containers.
//! Ordering is retained even when a transformation is nonmonotonic.

use super::transformations::TransformationDescription;
use crate::identification::PeptideIdentification;
use crate::kernel::MSExperiment;
use crate::kernel::features::{BaseFeature, ConsensusFeature, ConsensusMap, Feature, FeatureMap};
use crate::kernel::geometry::ConvexHull2D;
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::collections::BTreeMap;

/// Source `MapAlignmentTransformer` behavior with explicit resource limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MapAlignmentTransformer {
    pub store_original_rt: bool,
    /// Native extension. The source experiment overload leaves attached peptide
    /// RTs unchanged; feature/consensus peptide RTs are always transformed.
    pub transform_spectrum_identifications: bool,
    /// Total transformed coordinates, including hull points, per operation.
    pub max_rt_values: usize,
    /// Spectra, chromatograms, features, hulls and peptide records per operation.
    pub max_records: usize,
    /// Top-level features have depth zero; cannot exceed the kernel depth limit.
    pub max_subordinate_depth: usize,
}
impl Default for MapAlignmentTransformer {
    fn default() -> Self {
        Self {
            store_original_rt: false,
            transform_spectrum_identifications: false,
            max_rt_values: 1_000_000,
            max_records: 1_000_000,
            max_subordinate_depth: Feature::MAX_SUBORDINATE_DEPTH,
        }
    }
}
struct Budget {
    options: MapAlignmentTransformer,
    values: usize,
    records: usize,
}
impl Budget {
    fn new(options: MapAlignmentTransformer) -> Result<Self> {
        if options.max_rt_values == 0
            || options.max_records == 0
            || options.max_subordinate_depth > Feature::MAX_SUBORDINATE_DEPTH
        {
            return Err(bad("invalid alignment-transformer resource limits"));
        }
        Ok(Self {
            options,
            values: 0,
            records: 0,
        })
    }
    fn records(&mut self, amount: usize) -> Result<()> {
        self.records = self
            .records
            .checked_add(amount)
            .ok_or_else(|| bad("alignment record count overflow"))?;
        if self.records > self.options.max_records {
            return Err(bad("alignment record limit exceeded"));
        }
        Ok(())
    }
    fn value(&mut self, rt: f64, transformation: &TransformationDescription) -> Result<f64> {
        self.values = self
            .values
            .checked_add(1)
            .ok_or_else(|| bad("alignment coordinate count overflow"))?;
        if self.values > self.options.max_rt_values {
            return Err(bad("alignment coordinate limit exceeded"));
        }
        transformation.apply(rt)
    }
}
fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
struct PositionPlan {
    rt: f64,
    original: Option<String>,
}
impl PositionPlan {
    fn new(
        rt: f64,
        metadata: &BTreeMap<String, String>,
        transformation: &TransformationDescription,
        budget: &mut Budget,
    ) -> Result<Self> {
        let transformed = budget.value(rt, transformation)?;
        Ok(Self {
            rt: transformed,
            original: (budget.options.store_original_rt && !metadata.contains_key("original_RT"))
                .then(|| rt.to_string()),
        })
    }
    fn apply(self, rt: &mut f64, metadata: &mut BTreeMap<String, String>) {
        *rt = self.rt;
        if let Some(value) = self.original {
            metadata.insert("original_RT".into(), value);
        }
    }
}
struct PeptidePlan {
    rt: f64,
    original: Option<MetaValue>,
}
fn plan_peptides(
    ids: &[PeptideIdentification],
    transformation: &TransformationDescription,
    budget: &mut Budget,
) -> Result<Vec<Option<PeptidePlan>>> {
    budget.records(ids.len())?;
    ids.iter()
        .map(|id| {
            id.validate()?;
            match id.rt {
                None => Ok(None),
                Some(rt) => Ok(Some(PeptidePlan {
                    rt: budget.value(rt, transformation)?,
                    original: if budget.options.store_original_rt
                        && !id.metadata.contains_key("original_RT")
                    {
                        Some(MetaValue::try_from(rt)?)
                    } else {
                        None
                    },
                })),
            }
        })
        .collect()
}
fn apply_peptides(ids: &mut [PeptideIdentification], plans: Vec<Option<PeptidePlan>>) {
    for (id, plan) in ids.iter_mut().zip(plans) {
        if let Some(plan) = plan {
            id.rt = Some(plan.rt);
            if let Some(value) = plan.original {
                id.metadata.insert("original_RT".into(), value);
            }
        }
    }
}
struct BasePlan {
    position: PositionPlan,
    peptides: Vec<Option<PeptidePlan>>,
}
impl BasePlan {
    fn new(
        feature: &BaseFeature,
        transformation: &TransformationDescription,
        budget: &mut Budget,
    ) -> Result<Self> {
        Ok(Self {
            position: PositionPlan::new(feature.rt, &feature.metadata, transformation, budget)?,
            peptides: plan_peptides(&feature.peptide_identifications, transformation, budget)?,
        })
    }
    fn apply(self, feature: &mut BaseFeature) {
        self.position.apply(&mut feature.rt, &mut feature.metadata);
        apply_peptides(&mut feature.peptide_identifications, self.peptides);
    }
}
struct FeaturePlan {
    base: BasePlan,
    hulls: Vec<ConvexHull2D>,
    subordinates: Vec<FeaturePlan>,
}
impl FeaturePlan {
    fn new(
        feature: &Feature,
        transformation: &TransformationDescription,
        budget: &mut Budget,
        depth: usize,
    ) -> Result<Self> {
        if depth > budget.options.max_subordinate_depth {
            return Err(bad("alignment subordinate depth limit exceeded"));
        }
        budget.records(1)?;
        let base = BasePlan::new(&feature.base, transformation, budget)?;
        budget.records(feature.convex_hulls.len())?;
        let mut hulls = Vec::new();
        for hull in &feature.convex_hulls {
            // The source replaces scan envelopes with the transformed ordered
            // outline; preserve that representation rather than rebuilding a hull.
            let mut points = hull.hull_points();
            for point in &mut points {
                point.rt = budget.value(point.rt, transformation)?;
            }
            let mut transformed = ConvexHull2D::new();
            transformed.set_hull_points(&points)?;
            hulls.push(transformed);
        }
        let subordinates = feature
            .subordinates
            .iter()
            .map(|sub| Self::new(sub, transformation, budget, depth + 1))
            .collect::<Result<_>>()?;
        Ok(Self {
            base,
            hulls,
            subordinates,
        })
    }
    fn apply(self, feature: &mut Feature) {
        self.base.apply(&mut feature.base);
        feature.convex_hulls = self.hulls;
        for (sub, plan) in feature.subordinates.iter_mut().zip(self.subordinates) {
            plan.apply(sub);
        }
    }
}
fn plan_consensus(
    feature: &ConsensusFeature,
    transformation: &TransformationDescription,
    budget: &mut Budget,
) -> Result<ConsensusFeature> {
    budget.records(1)?;
    let base = BasePlan::new(&feature.base, transformation, budget)?;
    let handles = feature
        .handles()
        .iter()
        .map(|handle| {
            let mut handle = *handle;
            handle.rt = budget.value(handle.rt, transformation)?;
            Ok(handle)
        })
        .collect::<Result<_>>()?;
    let mut result = feature.clone();
    base.apply(&mut result.base);
    result.set_handles(handles)?;
    Ok(result)
}
impl MapAlignmentTransformer {
    /// Transform spectrum RTs and every chromatogram sample. Mass/intensity arrays
    /// retain their order and values. Ranges in the native kernel are computed
    /// from current data, so no stale range cache remains after transformation.
    pub fn transform_experiment(
        &self,
        experiment: &mut MSExperiment,
        transformation: &TransformationDescription,
    ) -> Result<()> {
        let mut budget = Budget::new(*self)?;
        budget.records(experiment.spectra.len())?;
        budget.records(experiment.chromatograms.len())?;
        experiment.validate()?;
        let mut spectra = Vec::new();
        for spectrum in &experiment.spectra {
            let position =
                PositionPlan::new(spectrum.rt, &spectrum.metadata, transformation, &mut budget)?;
            let ids = if self.transform_spectrum_identifications {
                Some(plan_peptides(
                    &spectrum.peptide_identifications,
                    transformation,
                    &mut budget,
                )?)
            } else {
                None
            };
            spectra.push((position, ids));
        }
        let mut chromatograms = Vec::new();
        for chromatogram in &experiment.chromatograms {
            let mut times = Vec::new();
            for peak in &chromatogram.peaks {
                times.push(budget.value(peak.rt, transformation)?);
            }
            // Kernel metadata currently stores strings. Display of a typed float
            // list supplies deterministic, round-trippable source-style syntax.
            let original =
                if self.store_original_rt && !chromatogram.metadata.contains_key("original_rt") {
                    Some(
                        MetaValue::try_from(
                            chromatogram.peaks.iter().map(|p| p.rt).collect::<Vec<_>>(),
                        )?
                        .to_string(),
                    )
                } else {
                    None
                };
            chromatograms.push((times, original));
        }
        for (spectrum, (position, ids)) in experiment.spectra.iter_mut().zip(spectra) {
            position.apply(&mut spectrum.rt, &mut spectrum.metadata);
            if let Some(ids) = ids {
                apply_peptides(&mut spectrum.peptide_identifications, ids);
            }
        }
        for (chromatogram, (times, original)) in
            experiment.chromatograms.iter_mut().zip(chromatograms)
        {
            for (peak, rt) in chromatogram.peaks.iter_mut().zip(times) {
                peak.rt = rt;
            }
            if let Some(original) = original {
                chromatogram.metadata.insert("original_rt".into(), original);
            }
        }
        Ok(())
    }
    /// Includes assigned peptide RTs, ordered hull points, all subordinate
    /// features, and unassigned peptide RTs. Existing original_RT is retained.
    pub fn transform_feature_map(
        &self,
        map: &mut FeatureMap,
        transformation: &TransformationDescription,
    ) -> Result<()> {
        let mut budget = Budget::new(*self)?;
        map.validate()?;
        let plans = map
            .features
            .iter()
            .map(|f| FeaturePlan::new(f, transformation, &mut budget, 0))
            .collect::<Result<Vec<_>>>()?;
        let ids = plan_peptides(
            &map.unassigned_peptide_identifications,
            transformation,
            &mut budget,
        )?;
        for (feature, plan) in map.features.iter_mut().zip(plans) {
            plan.apply(feature);
        }
        apply_peptides(&mut map.unassigned_peptide_identifications, ids);
        Ok(())
    }
    /// Includes every handle RT and assigned/unassigned peptide RT. Handles keep
    /// their `(map_index, unique_id)` identity and are not re-consensused.
    pub fn transform_consensus_map(
        &self,
        map: &mut ConsensusMap,
        transformation: &TransformationDescription,
    ) -> Result<()> {
        let mut budget = Budget::new(*self)?;
        map.validate()?;
        let features = map
            .features
            .iter()
            .map(|f| plan_consensus(f, transformation, &mut budget))
            .collect::<Result<_>>()?;
        let ids = plan_peptides(
            &map.unassigned_peptide_identifications,
            transformation,
            &mut budget,
        )?;
        map.features = features;
        apply_peptides(&mut map.unassigned_peptide_identifications, ids);
        Ok(())
    }
    /// Peptides without an RT remain unchanged and do not gain original_RT.
    pub fn transform_peptide_identifications(
        &self,
        ids: &mut [PeptideIdentification],
        transformation: &TransformationDescription,
    ) -> Result<()> {
        let mut budget = Budget::new(*self)?;
        let plans = plan_peptides(ids, transformation, &mut budget)?;
        apply_peptides(ids, plans);
        Ok(())
    }
    pub fn transform_feature(
        &self,
        feature: &mut Feature,
        transformation: &TransformationDescription,
    ) -> Result<()> {
        let mut budget = Budget::new(*self)?;
        feature.validate()?;
        let plan = FeaturePlan::new(feature, transformation, &mut budget, 0)?;
        plan.apply(feature);
        Ok(())
    }
    pub fn transform_consensus_feature(
        &self,
        feature: &mut ConsensusFeature,
        transformation: &TransformationDescription,
    ) -> Result<()> {
        let mut budget = Budget::new(*self)?;
        feature.validate()?;
        *feature = plan_consensus(feature, transformation, &mut budget)?;
        Ok(())
    }
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The features table of a FeatureMap as an Arrow table.
//!
//! Ports the features-table members of `FORMAT/FeatureMapArrowIO.h`
//! (`FeatureMapArrowIO.cpp`): `exportFeaturesToArrow` and
//! `importFeaturesFromArrow`, the [`FeatureSchema`](crate::format::arrow_schema_registry::FeatureSchema) table that
//! `features.parquet` of a `.featureparquet` bundle holds. The PSM table and
//! the bundle functions (`exportPSMsToArrow`, `importPSMsFromArrow`,
//! `exportToParquet`, `importFromParquet`) build and read their PSM rows with
//! `QPXFile`, which is not ported; see `docs/ARROW_PARQUET_SUPPORT.md`.

use crate::format::arrow_io_helpers::{
    self as helpers, get_double_value, get_float_value, get_int32_value, get_int64_value, is_null,
    list_array, list_element, list_range, meta_entries, metavalues_array, read_meta_values,
    struct_array,
};
use crate::format::arrow_schema_registry::FeatureSchema;
use crate::format::arrow_schema_registry::validation::{Mode, validate};
use crate::kernel::features::{BaseFeature, Feature, FeatureMap};
use crate::kernel::geometry::{ConvexHull2D, Point2D};
use crate::metadata::MetaValue;
use crate::{Error, Result};
use arrow_array::cast::AsArray;
use arrow_array::types::Float64Type;
use arrow_array::{
    Array, ArrayRef, Float32Array, Float64Array, Int32Array, Int64Array, RecordBatch,
};
use std::collections::HashMap;
use std::sync::Arc;

/// Import and export of a FeatureMap's features table.
///
/// Ports the features-table half of `OpenMS::FeatureMapArrowIO`, a class of
/// static members; each is an associated function here.
#[derive(Clone, Copy, Debug)]
pub struct FeatureMapArrowIO;

/// One exported row before the columns are built.
struct Row<'a> {
    feature: &'a Feature,
    parent: Option<i64>,
    depth: i32,
}

/// The features depth first, each followed by its subordinates, as the
/// source's recursive `appendFeature` visits them.
fn flatten(map: &FeatureMap) -> Result<Vec<Row<'_>>> {
    let mut rows = Vec::new();
    let mut pending: Vec<(&Feature, Option<i64>, i32)> = map
        .features
        .iter()
        .rev()
        .map(|feature| (feature, None, 0))
        .collect();
    while let Some((feature, parent, depth)) = pending.pop() {
        if usize::try_from(depth).unwrap_or(usize::MAX) > Feature::MAX_SUBORDINATE_DEPTH {
            return Err(Error::InvalidValue(format!(
                "FeatureMapArrowIO: subordinate feature depth exceeds {}",
                Feature::MAX_SUBORDINATE_DEPTH
            )));
        }
        rows.push(Row {
            feature,
            parent,
            depth,
        });
        // The source's static_cast<int64_t>(UInt64).
        let id = feature.unique_id as i64;
        pending.extend(
            feature
                .subordinates
                .iter()
                .rev()
                .map(|sub| (sub, Some(id), depth + 1)),
        );
    }
    Ok(rows)
}

/// The bounding box of all hulls, as the source's `getConvexHull()
/// .getBoundingBox()` computes it for a feature with hulls: the union of
/// each hull's bounding box. `None` when there is no hull point at all.
fn bounding_box(hulls: &[ConvexHull2D]) -> Option<[f64; 4]> {
    let mut bounds: Option<[f64; 4]> = None;
    for hull in hulls {
        for point in hull.hull_points() {
            let b = bounds.get_or_insert([point.rt, point.rt, point.mz, point.mz]);
            b[0] = b[0].min(point.rt);
            b[1] = b[1].max(point.rt);
            b[2] = b[2].min(point.mz);
            b[3] = b[3].max(point.mz);
        }
    }
    bounds
}

impl FeatureMapArrowIO {
    /// The features of `map` as an Arrow table, one row per feature and
    /// subordinate.
    ///
    /// Ports `exportFeaturesToArrow(feature_map)`. The features are listed
    /// depth first, each followed by its subordinates. Columns, in
    /// [`FeatureSchema`] order: `unique_id` (the `UInt64` id as `int64`,
    /// wrapping as the source's cast does), `parent_feature_id` (null for a
    /// top-level feature), `depth` (0 for top level), `rt`, `mz`,
    /// `intensity`, `charge`, `quality` (the overall quality), `quality_rt`,
    /// `quality_mz`, `width` (null when 0, the unset default), the bounding
    /// box of the hulls (`rt_bb_min`, `rt_bb_max`, `mz_bb_min`, `mz_bb_max`;
    /// null when the feature has no hull), `convex_hulls` (each hull's index
    /// and points) and `metavalues` (every metadata entry except `FWHM`,
    /// which `width` carries).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when subordinates nest deeper than
    /// [`Feature::MAX_SUBORDINATE_DEPTH`] (the source recurses without a
    /// bound) or a list column exceeds Arrow's 2^31-element limit.
    pub fn export_features_to_arrow(map: &FeatureMap) -> Result<RecordBatch> {
        let rows = flatten(map)?;
        let n = rows.len();
        let mut unique_id = Vec::with_capacity(n);
        let mut parent = Vec::with_capacity(n);
        let mut depth = Vec::with_capacity(n);
        let mut rt = Vec::with_capacity(n);
        let mut mz = Vec::with_capacity(n);
        let mut intensity = Vec::with_capacity(n);
        let mut charge = Vec::with_capacity(n);
        let mut quality = Vec::with_capacity(n);
        let mut quality_rt = Vec::with_capacity(n);
        let mut quality_mz = Vec::with_capacity(n);
        let mut width = Vec::with_capacity(n);
        let mut bb = [
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
            Vec::with_capacity(n),
        ];
        let mut hull_counts = Vec::with_capacity(n);
        let mut hull_index = Vec::new();
        let mut point_counts = Vec::new();
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        let mut metavalues = Vec::with_capacity(n);
        for row in &rows {
            let f = row.feature;
            unique_id.push(f.unique_id as i64);
            parent.push(row.parent);
            depth.push(row.depth);
            rt.push(f.rt);
            mz.push(f.mz);
            intensity.push(f.intensity);
            charge.push(f.charge);
            quality.push(f.quality);
            quality_rt.push(f.quality_rt);
            quality_mz.push(f.quality_mz);
            width.push((f.width != 0.0).then_some(f.width));
            let bounds = if f.convex_hulls.is_empty() {
                None
            } else {
                bounding_box(&f.convex_hulls)
            };
            for (column, index) in bb.iter_mut().zip(0..4) {
                column.push(bounds.map(|b| b[index]));
            }
            hull_counts.push(Some(f.convex_hulls.len()));
            for (index, hull) in f.convex_hulls.iter().enumerate() {
                hull_index.push(i32::try_from(index).map_err(|_| {
                    Error::InvalidValue("more than 2^31 convex hulls in one feature".into())
                })?);
                let points = hull.hull_points();
                point_counts.push(Some(points.len()));
                for point in points {
                    xs.push(point.rt);
                    ys.push(point.mz);
                }
            }
            metavalues.push(Some(meta_entries(&f.metadata, &["FWHM"])));
        }
        let hull_type = FeatureSchema::convex_hull_type();
        let hull_entry = list_element(&hull_type)?;
        let arrow_schema::DataType::Struct(hull_fields) = hull_entry else {
            return Err(Error::InvalidValue(format!(
                "{hull_entry} is not a struct type"
            )));
        };
        let points_type = hull_fields[1].data_type();
        let points = struct_array(
            list_element(points_type)?,
            vec![
                Arc::new(Float64Array::from(xs)) as ArrayRef,
                Arc::new(Float64Array::from(ys)),
            ],
        )?;
        let hulls = struct_array(
            hull_entry,
            vec![
                Arc::new(Int32Array::from(hull_index)) as ArrayRef,
                list_array(points_type, &point_counts, points)?,
            ],
        )?;
        let [bb0, bb1, bb2, bb3] = bb;
        let columns: Vec<ArrayRef> = vec![
            Arc::new(Int64Array::from(unique_id)),
            Arc::new(Int64Array::from(parent)),
            Arc::new(Int32Array::from(depth)),
            Arc::new(Float64Array::from(rt)),
            Arc::new(Float64Array::from(mz)),
            Arc::new(Float32Array::from(intensity)),
            Arc::new(Int32Array::from(charge)),
            Arc::new(Float32Array::from(quality)),
            Arc::new(Float32Array::from(quality_rt)),
            Arc::new(Float32Array::from(quality_mz)),
            Arc::new(Float32Array::from(width)),
            Arc::new(Float64Array::from(bb0)),
            Arc::new(Float64Array::from(bb1)),
            Arc::new(Float64Array::from(bb2)),
            Arc::new(Float64Array::from(bb3)),
            list_array(&hull_type, &hull_counts, hulls)?,
            metavalues_array(&FeatureSchema::metavalues_type(), &metavalues)?,
        ];
        let schema = FeatureSchema::schema();
        let table = RecordBatch::try_new(schema.clone(), columns).map_err(|error| {
            Error::InvalidValue(format!(
                "FeatureMapArrowIO: Failed to build the features table: {error}"
            ))
        })?;
        let result = validate(&table, &schema, Mode::Strict);
        if !result.valid {
            return Err(Error::InvalidValue(format!(
                "Schema validation failed: {result}"
            )));
        }
        Ok(table)
    }

    /// Append the features of a features table to `map`.
    ///
    /// Ports `importFeaturesFromArrow(table, feature_map)`. An empty table
    /// adds nothing. The table is checked in subset mode. Each row becomes a
    /// feature with its id, position, intensity, charge, qualities, width (set
    /// with its `FWHM` metadata, as the source's `setWidth`, when not null and
    /// not 0), hulls (in list order; `hull_index` is not consulted) and
    /// metadata. Subordinates are then attached to the feature whose
    /// `unique_id` their `parent_feature_id` names, deepest rows first and in
    /// row order within a depth, so each keeps its subtree; a row whose parent
    /// is not in the table becomes a top-level feature with a warning line,
    /// as the source logs it. Top-level features are appended in row order.
    ///
    /// When two rows share a `unique_id`, children attach to the later one,
    /// as the source's id map keeps the last. A row attached to a parent that
    /// was itself already attached elsewhere is lost in the source (it is
    /// moved into a moved-from feature); the port drops it too and adds a
    /// warning line.
    ///
    /// Returns the warning lines.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an incompatible schema, a column of the
    /// wrong type or a hull point that is not finite;
    /// [`Error::MissingInformation`] for a missing required column. `map` is
    /// unchanged on error.
    pub fn import_features_from_arrow(
        table: &RecordBatch,
        map: &mut FeatureMap,
    ) -> Result<Vec<String>> {
        let mut warnings = Vec::new();
        if table.num_rows() == 0 {
            return Ok(warnings);
        }
        let result = validate(table, &FeatureSchema::schema(), Mode::Subset);
        if !result.valid {
            return Err(Error::InvalidValue(format!(
                "Incompatible schema: {result}"
            )));
        }
        let required = |name: &str| {
            helpers::get_column(table, name, true).ok_or_else(|| {
                Error::MissingInformation(
                    "FeatureMapArrowIO: Missing one or more required columns".into(),
                )
            })
        };
        let unique_id = required(FeatureSchema::UNIQUE_ID)?;
        let parent_id = required(FeatureSchema::PARENT_FEATURE_ID)?;
        let depth = required(FeatureSchema::DEPTH)?;
        let rt = required(FeatureSchema::RT)?;
        let mz = required(FeatureSchema::MZ)?;
        let intensity = required(FeatureSchema::INTENSITY)?;
        let charge = required(FeatureSchema::CHARGE)?;
        let quality = required(FeatureSchema::QUALITY)?;
        let quality_rt = required(FeatureSchema::QUALITY_RT)?;
        let quality_mz = required(FeatureSchema::QUALITY_MZ)?;
        let width = helpers::get_column(table, FeatureSchema::WIDTH, false);
        let hulls = helpers::get_column(table, FeatureSchema::CONVEX_HULLS, false);
        let metavalues = helpers::get_column(table, FeatureSchema::METAVALUES, false);
        struct Entry {
            depth: i32,
            row: usize,
            unique_id: i64,
            parent: Option<i64>,
            feature: Option<Feature>,
        }
        let mut entries = Vec::with_capacity(table.num_rows());
        for row in 0..table.num_rows() {
            let id = get_int64_value(Some(unique_id.as_ref()), row, 0)?;
            let parent = if is_null(Some(parent_id.as_ref()), row) {
                None
            } else {
                Some(get_int64_value(Some(parent_id.as_ref()), row, 0)?)
            };
            let base = BaseFeature {
                // The source's static_cast<UInt64>(int64_t).
                unique_id: id as u64,
                rt: get_double_value(Some(rt.as_ref()), row, 0.0)?,
                mz: get_double_value(Some(mz.as_ref()), row, 0.0)?,
                intensity: get_float_value(Some(intensity.as_ref()), row, 0.0)?,
                charge: get_int32_value(Some(charge.as_ref()), row, 0)?,
                quality: get_float_value(Some(quality.as_ref()), row, 0.0)?,
                ..BaseFeature::default()
            };
            let mut feature = Feature {
                base,
                quality_rt: get_float_value(Some(quality_rt.as_ref()), row, 0.0)?,
                quality_mz: get_float_value(Some(quality_mz.as_ref()), row, 0.0)?,
                ..Feature::default()
            };
            let w = get_float_value(width.as_deref(), row, 0.0)?;
            if w != 0.0 {
                // BaseFeature::setWidth: the width and its FWHM metadata bridge.
                feature.width = w;
                feature
                    .metadata
                    .insert("FWHM".into(), MetaValue::source_float(f64::from(w)));
            }
            if let Some(column) = hulls.as_deref() {
                feature.convex_hulls = read_hulls(column, row)?;
            }
            read_meta_values(metavalues.as_deref(), row, &mut feature.metadata, &[])?;
            entries.push(Entry {
                depth: get_int32_value(Some(depth.as_ref()), row, 0)?,
                row,
                unique_id: id,
                parent,
                feature: Some(feature),
            });
        }
        // Deepest first, row order within a depth (the source's stable_sort).
        entries.sort_by_key(|entry| std::cmp::Reverse(entry.depth));
        let index: HashMap<i64, usize> = entries
            .iter()
            .enumerate()
            .map(|(i, entry)| (entry.unique_id, i))
            .collect();
        for i in 0..entries.len() {
            let Some(parent) = entries[i].parent else {
                continue;
            };
            match index.get(&parent).copied() {
                Some(at) => {
                    let child = entries[i].feature.take();
                    match (child, entries[at].feature.as_mut()) {
                        (Some(child), Some(target)) => target.subordinates.push(child),
                        (Some(_), None) => warnings.push(format!(
                            "FeatureMapArrowIO: feature {} was attached to feature {parent}, \
                             which had already been attached elsewhere; it is lost, as in \
                             the source",
                            entries[i].unique_id
                        )),
                        (None, _) => {}
                    }
                }
                None => {
                    warnings.push(format!(
                        "FeatureMapArrowIO: Could not find parent feature with id {parent} for \
                         feature {}. Will be added as top-level.",
                        entries[i].unique_id
                    ));
                    entries[i].parent = None;
                }
            }
        }
        let mut top: Vec<&mut Entry> = entries
            .iter_mut()
            .filter(|entry| entry.parent.is_none())
            .collect();
        top.sort_by_key(|entry| entry.row);
        map.features
            .extend(top.into_iter().filter_map(|entry| entry.feature.take()));
        Ok(warnings)
    }
}

/// The convex hulls at `row` of a `convex_hulls` column.
fn read_hulls(column: &dyn Array, row: usize) -> Result<Vec<ConvexHull2D>> {
    let type_error = || {
        Error::InvalidValue(format!(
            "FeatureMapArrowIO: convex_hulls has unexpected type {}",
            column.data_type()
        ))
    };
    let Some((values, start, end)) = list_range(column, row)? else {
        return Ok(Vec::new());
    };
    let entries = values
        .as_struct_opt()
        .filter(|s| s.num_columns() >= 2)
        .ok_or_else(type_error)?;
    let lists = entries.column(1);
    let mut hulls = Vec::with_capacity(end - start);
    for h in start..end {
        let mut points = Vec::new();
        if let Some((values, from, to)) = list_range(lists.as_ref(), h)? {
            let coordinates = values
                .as_struct_opt()
                .filter(|s| s.num_columns() >= 2)
                .ok_or_else(type_error)?;
            let x = coordinates
                .column(0)
                .as_primitive_opt::<Float64Type>()
                .ok_or_else(type_error)?;
            let y = coordinates
                .column(1)
                .as_primitive_opt::<Float64Type>()
                .ok_or_else(type_error)?;
            points.reserve(to - from);
            for p in from..to {
                points.push(Point2D::new(x.value(p), y.value(p)));
            }
        }
        let mut hull = ConvexHull2D::new();
        hull.set_hull_points(&points)?;
        hulls.push(hull);
    }
    Ok(hulls)
}

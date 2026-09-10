// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Native featureXML 1.9 transport, including legacy hulls and nested features.
//! Read filters are half-open and never affect writing. See FEATUREXML_SUPPORT.md.

use super::identification_xml::{self as xml, Node};
use super::{FileType, map_xml, path_io};
use crate::chemistry::ModificationsDB;
use crate::kernel::{ConvexHull2D, Feature, FeatureMap, Point2D};
use crate::metadata::{MetaInfo, MetaValueData};
use crate::{Error, Result};
use std::collections::BTreeSet;
use std::io::{BufRead, Write};
use std::ops::Range;
use std::path::Path;

pub const VERSION: &str = "1.9";

/// Passive source options. Empty/inverted ranges select no values. `size_only`
/// is retained by source options but has no effect on `read`; use `read_size`.
#[derive(Clone, Debug, PartialEq)]
pub struct FeatureFileOptions {
    pub load_convex_hulls: bool,
    pub load_subordinates: bool,
    pub metadata_only: bool,
    pub size_only: bool,
    pub rt_range: Option<Range<f64>>,
    pub mz_range: Option<Range<f64>>,
    pub intensity_range: Option<Range<f64>>,
}
impl Default for FeatureFileOptions {
    fn default() -> Self {
        Self {
            load_convex_hulls: true,
            load_subordinates: true,
            metadata_only: false,
            size_only: false,
            rt_range: None,
            mz_range: None,
            intensity_range: None,
        }
    }
}

/// Cumulative conversion limits. Depth counts subordinate feature levels.
#[derive(Clone, Copy, Debug)]
pub struct Limits {
    pub max_xml_bytes: u64,
    pub max_records: usize,
    pub max_list_items: usize,
    pub max_depth: usize,
    pub max_work: usize,
    pub max_payload_bytes: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_xml_bytes: 64 * 1024 * 1024,
            max_records: 1_000_000,
            max_list_items: 1_000_000,
            max_depth: 128,
            max_work: 50_000_000,
            max_payload_bytes: 256 * 1024 * 1024,
        }
    }
}
#[derive(Clone, Debug, Default)]
pub struct ReadOptions {
    pub feature_options: FeatureFileOptions,
    pub limits: Limits,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct WriteOptions {
    pub limits: Limits,
}

fn bad(text: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: text.into(),
    }
}
fn unsupported(text: impl Into<String>) -> Error {
    Error::Unsupported(text.into())
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b)
        .ok_or_else(|| bad("featureXML size overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b)
        .ok_or_else(|| bad("featureXML size overflow"))
}
fn xml_options(limits: Limits) -> xml::ReadOptions {
    xml::ReadOptions {
        max_xml_bytes: limits.max_xml_bytes,
        max_records: limits.max_records,
        max_list_items: limits.max_list_items,
    }
}
fn uid(text: &str) -> Result<u64> {
    map_xml::unique_id(text)
}
fn scalar(text: &str) -> Result<f32> {
    let value = xml::finite(text)? as f32;
    if !value.is_finite() {
        return Err(bad("featureXML value exceeds f32 range"));
    }
    Ok(value)
}
fn check(node: &Node, attrs: &[&str], children: &[&str], text: bool) -> Result<()> {
    // Source permits metadata wherever the current object is active. Unlike the
    // identification schema helper, geometry does not impose a child sort order.
    if !text && !node.text.trim().is_empty() {
        return Err(bad(format!("unexpected text in {}", node.name)));
    }
    for key in node.attrs.keys() {
        if !attrs.contains(&key.as_str()) {
            return Err(unsupported(format!("{} attribute {key}", node.name)));
        }
    }
    for child in &node.children {
        if !children.contains(&child.name.as_str()) {
            return Err(unsupported(format!("{} inside {}", child.name, node.name)));
        }
    }
    Ok(())
}

struct Work {
    remaining: usize,
    bytes: usize,
    records: usize,
    limits: Limits,
}
impl Work {
    fn new(limits: Limits) -> Self {
        Self {
            remaining: limits.max_work,
            bytes: limits.max_payload_bytes,
            records: 0,
            limits,
        }
    }
    fn consume(&mut self, count: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(count)
            .ok_or_else(|| bad("featureXML cumulative work limit exceeded"))?;
        Ok(())
    }
    fn allocate(&mut self, bytes: usize) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| bad("featureXML cumulative payload limit exceeded"))?;
        Ok(())
    }
    fn record(&mut self) -> Result<()> {
        self.consume(1)?;
        self.records = add(self.records, 1)?;
        if self.records > self.limits.max_records {
            return Err(bad("featureXML record limit exceeded"));
        }
        Ok(())
    }
    fn slots<T>(&mut self, count: usize) -> Result<()> {
        self.consume(count)?;
        self.allocate(mul(count, std::mem::size_of::<T>())?)
    }
    fn string(&mut self, text: &str) -> Result<()> {
        self.consume(text.len())?;
        self.allocate(mul(text.len(), 6)?)
    }
    fn node(&mut self, node: &Node, depth: usize) -> Result<()> {
        if depth > 2 * self.limits.max_depth.min(Feature::MAX_SUBORDINATE_DEPTH) + 16 {
            return Err(bad("featureXML nesting limit exceeded"));
        }
        self.record()?;
        self.allocate(std::mem::size_of::<Node>())?;
        self.string(&node.name)?;
        self.string(&node.text)?;
        self.slots::<(String, String)>(node.attrs.len())?;
        for (key, value) in &node.attrs {
            self.string(key)?;
            self.string(value)?;
        }
        for child in &node.children {
            self.node(child, depth + 1)?;
        }
        Ok(())
    }
    fn meta(&mut self, metadata: &MetaInfo) -> Result<()> {
        xml::measure_metadata(metadata, &mut self.remaining, &mut self.bytes)?;
        for value in metadata.values() {
            let count = match value.data() {
                MetaValueData::StringList(v) => v.len(),
                MetaValueData::IntegerList(v) => v.len(),
                MetaValueData::FloatList(v) => v.len(),
                _ => 0,
            };
            if count > self.limits.max_list_items {
                return Err(bad("featureXML metadata list limit exceeded"));
            }
            if let Some(unit) = value.unit() {
                self.string(unit.accession())?;
                self.string(unit.name())?;
                self.string(unit.cv_ref())?;
            }
            value.validate()?;
        }
        Ok(())
    }
    fn feature(&mut self, feature: &Feature, depth: usize) -> Result<()> {
        if depth > self.limits.max_depth.min(Feature::MAX_SUBORDINATE_DEPTH) {
            return Err(bad("featureXML subordinate depth exceeded"));
        }
        self.record()?;
        // Scalar text can be hundreds of bytes for finite subnormal f64 values;
        // nested IDs also grow with depth. Charge this before formatting either.
        self.allocate(std::mem::size_of::<Feature>() + 8 * std::mem::size_of::<Node>() + 4096)?;
        self.allocate(mul(depth + 1, 128)?)?;
        for value in [
            feature.rt,
            feature.mz,
            f64::from(feature.intensity),
            f64::from(feature.quality),
            f64::from(feature.quality_rt),
            f64::from(feature.quality_mz),
            f64::from(feature.width),
        ] {
            if !value.is_finite() {
                return Err(bad("nonfinite feature value"));
            }
        }
        if feature.width < 0.0 {
            return Err(bad("negative feature width"));
        }
        // Width has no XML element. Source reads FWHM into top-level width only.
        if depth != 0 && feature.width != 0.0 {
            return Err(unsupported(
                "featureXML cannot restore nonzero subordinate width",
            ));
        }
        if depth == 0 {
            let stored = feature
                .metadata
                .get("FWHM")
                .map(|v| v.as_f64())
                .transpose()?
                .unwrap_or(0.0) as f32;
            if stored != feature.width {
                return Err(bad("feature width must match FWHM metadata; use set_width"));
            }
        }
        self.meta(&feature.metadata)?;
        for hull in &feature.convex_hulls {
            // Both outline-only and scan-backed hulls expose a bounded size.
            let count = hull.point_count_bound();
            self.slots::<Point2D>(mul(count, 3)?)?;
            self.slots::<Node>(count)?;
            self.allocate(mul(count, 2048)?)?;
        }
        for child in &feature.subordinates {
            self.feature(child, depth + 1)?;
        }
        Ok(())
    }
}

pub fn read(input: impl BufRead) -> Result<FeatureMap> {
    read_with_options(input, &ReadOptions::default())
}
pub fn read_with_options(input: impl BufRead, options: &ReadOptions) -> Result<FeatureMap> {
    read_with_registry(input, options, ModificationsDB::global())
}
pub fn read_with_registry(
    input: impl BufRead,
    options: &ReadOptions,
    registry: &ModificationsDB,
) -> Result<FeatureMap> {
    Ok(read_document(input, options, registry, false)?.0)
}
/// Atomic replacement: neither parse nor conversion failure changes the target.
pub fn read_into(
    input: impl BufRead,
    target: &mut FeatureMap,
    options: &ReadOptions,
) -> Result<()> {
    let draft = read_with_options(input, options)?;
    *target = draft;
    Ok(())
}
pub fn read_size(input: impl BufRead, options: &ReadOptions) -> Result<usize> {
    Ok(read_document(input, options, ModificationsDB::global(), true)?.1)
}

fn read_document(
    input: impl BufRead,
    options: &ReadOptions,
    registry: &ModificationsDB,
    size_only: bool,
) -> Result<(FeatureMap, usize)> {
    let opts = xml_options(options.limits);
    let stop = options.feature_options.metadata_only || size_only;
    let mut work = Work::new(options.limits);
    let root = xml::parse_xml_with_budget(
        input,
        &opts,
        2 * options.limits.max_depth.min(Feature::MAX_SUBORDINATE_DEPTH) + 16,
        stop.then_some("featureList"),
        &mut work.remaining,
        &mut work.bytes,
    )?;
    if root.name != "featureMap" {
        return Err(bad("expected featureMap root"));
    }
    check(
        &root,
        &[
            "version",
            "document_id",
            "id",
            "unique_id",
            "xmlns:xsi",
            "xsi:noNamespaceSchemaLocation",
        ],
        &[
            "UserParam",
            "userParam",
            "dataProcessing",
            "IdentificationRun",
            "UnassignedPeptideIdentification",
            "featureList",
            "description",
        ],
        false,
    )?;
    xml::measure_node(&root, &mut work.remaining, &mut work.bytes)?;
    let mut map = FeatureMap {
        identifier: root.optional("document_id").unwrap_or("").into(),
        ..Default::default()
    };
    if let Some(id) = root.optional("id") {
        map.unique_id = uid(id)?;
    }
    if let Some(id) = root.optional("unique_id") {
        map.unique_id = uid(id)?;
    }
    let mut registry = xml::clone_registry(registry, &mut work.remaining, &mut work.bytes)?;
    let mut context = map_xml::ReadContext::default();
    for node in &root.children {
        match node.name.as_str() {
            "IdentificationRun" => map.protein_identifications.push(map_xml::read_run(
                node,
                &opts,
                &mut context,
                &mut registry,
                &mut work.remaining,
                &mut work.bytes,
            )?),
            "dataProcessing" => map
                .data_processing
                .push(map_xml::read_processing(node, &opts)?),
            _ => {}
        }
    }
    map.metadata = metadata(&root, &opts)?;
    let mut count = 0;
    let mut seen_list = false;
    for node in &root.children {
        match node.name.as_str() {
            "UnassignedPeptideIdentification" => {
                map.unassigned_peptide_identifications
                    .push(map_xml::read_peptide(
                        node,
                        &opts,
                        &context,
                        &registry,
                        &mut work.remaining,
                        &mut work.bytes,
                    )?)
            }
            "featureList" => {
                if seen_list {
                    return Err(bad("multiple featureList elements"));
                }
                seen_list = true;
                if options.feature_options.metadata_only {
                    break;
                }
                count = xml::number(node.get("count")?)?;
                if size_only {
                    break;
                }
                check(node, &["count"], &["feature"], false)?;
                work.slots::<Feature>(node.children.len())?;
                for child in &node.children {
                    let mut feature =
                        read_feature(child, options, &opts, &context, &registry, &mut work, 0)?;
                    if accepts(&feature, &options.feature_options)? {
                        if let Some(width) = feature.metadata.get("FWHM") {
                            let width = width.as_f64()? as f32;
                            feature.set_width(width)?;
                        }
                        map.features.push(feature);
                    }
                }
            }
            _ => {}
        }
    }
    if !seen_list {
        return Err(bad("featureMap requires featureList"));
    }
    Ok((map, count))
}
fn metadata(node: &Node, opts: &xml::ReadOptions) -> Result<MetaInfo> {
    // Normalize only the old spelling; all values use the shared typed codec.
    if node.children.iter().all(|n| n.name != "userParam") {
        return xml::read_meta(node, opts);
    }
    let mut holder = Node::new("metadata");
    for child in &node.children {
        if matches!(child.name.as_str(), "UserParam" | "userParam") {
            let mut child = child.clone();
            child.name = "UserParam".into();
            holder.children.push(child);
        }
    }
    xml::read_meta(&holder, opts)
}
fn accepts(feature: &Feature, options: &FeatureFileOptions) -> Result<bool> {
    for (value, range) in [
        (feature.rt, &options.rt_range),
        (feature.mz, &options.mz_range),
        (f64::from(feature.intensity), &options.intensity_range),
    ] {
        if let Some(range) = range {
            if !range.start.is_finite() || !range.end.is_finite() {
                return Err(bad("nonfinite feature filter range"));
            }
            if !range.contains(&value) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}
fn read_feature(
    node: &Node,
    options: &ReadOptions,
    opts: &xml::ReadOptions,
    context: &map_xml::ReadContext,
    registry: &ModificationsDB,
    work: &mut Work,
    depth: usize,
) -> Result<Feature> {
    if depth > options.limits.max_depth.min(Feature::MAX_SUBORDINATE_DEPTH) {
        return Err(bad("featureXML subordinate depth exceeded"));
    }
    check(
        node,
        &["id"],
        &[
            "position",
            "intensity",
            "quality",
            "overallquality",
            "charge",
            "convexhull",
            "subordinate",
            "PeptideIdentification",
            "UserParam",
            "userParam",
            "model",
            "description",
        ],
        false,
    )?;
    let mut feature = Feature::from(crate::kernel::BaseFeature {
        unique_id: uid(node.get("id")?)?,
        ..Default::default()
    });
    for child in &node.children {
        match child.name.as_str() {
            "position" | "quality" => {
                check(child, &["dim"], &[], true)?;
                let dim: usize = xml::number(child.get("dim")?)?;
                match (child.name.as_str(), dim) {
                    ("position", 0) => feature.rt = xml::finite(&child.text)?,
                    ("position", 1) => feature.mz = xml::finite(&child.text)?,
                    ("quality", 0) => feature.quality_rt = scalar(&child.text)?,
                    ("quality", 1) => feature.quality_mz = scalar(&child.text)?,
                    _ => return Err(bad("feature dimension must be zero or one")),
                }
            }
            "intensity" | "overallquality" | "charge" => {
                check(child, &[], &[], true)?;
                match child.name.as_str() {
                    "intensity" => feature.intensity = scalar(&child.text)?,
                    "overallquality" => feature.quality = scalar(&child.text)?,
                    _ => feature.charge = xml::number(&child.text)?,
                }
            }
            "convexhull" if options.feature_options.load_convex_hulls => {
                feature.convex_hulls.push(read_hull(child, work)?);
            }
            "subordinate" if options.feature_options.load_subordinates => {
                check(child, &[], &["feature"], false)?;
                work.slots::<Feature>(child.children.len())?;
                for sub in &child.children {
                    let value =
                        read_feature(sub, options, opts, context, registry, work, depth + 1)?;
                    if accepts(&value, &options.feature_options)? {
                        feature.subordinates.push(value);
                    }
                }
            }
            "PeptideIdentification" => {
                feature.peptide_identifications.push(map_xml::read_peptide(
                    child,
                    opts,
                    context,
                    registry,
                    &mut work.remaining,
                    &mut work.bytes,
                )?);
            }
            "UserParam" | "userParam" => {
                let mut holder = Node::new("metadata");
                let mut entry = child.clone();
                entry.name = "UserParam".into();
                holder.children.push(entry);
                let values = xml::read_meta(&holder, opts)?;
                feature.metadata.extend(values);
            }
            _ => {}
        }
    }
    Ok(feature)
}
fn read_hull(node: &Node, work: &mut Work) -> Result<ConvexHull2D> {
    check(node, &["nr"], &["pt", "hullpoint"], false)?;
    work.slots::<Point2D>(mul(node.children.len(), 2)?)?;
    let mut points = Vec::new();
    points
        .try_reserve_exact(node.children.len())
        .map_err(|_| bad("hull allocation failed"))?;
    for child in &node.children {
        let point = if child.name == "pt" {
            check(child, &["x", "y"], &[], false)?;
            Point2D::new(xml::finite(child.get("x")?)?, xml::finite(child.get("y")?)?)
        } else {
            check(child, &[], &["hposition"], false)?;
            let mut values = [0.0; 2];
            for position in &child.children {
                check(position, &["dim"], &[], true)?;
                let dim: usize = xml::number(position.get("dim")?)?;
                if dim > 1 {
                    return Err(bad("hull dimension must be zero or one"));
                }
                values[dim] = xml::finite(&position.text)?;
            }
            Point2D::new(values[0], values[1])
        };
        points.push(point);
    }
    let mut hull = ConvexHull2D::new();
    hull.set_hull_points(&points)?;
    Ok(hull)
}

pub fn write(output: impl Write, map: &FeatureMap) -> Result<()> {
    write_with_options(output, map, &WriteOptions::default())
}
pub fn write_with_options(
    output: impl Write,
    map: &FeatureMap,
    options: &WriteOptions,
) -> Result<()> {
    write_with_registry(output, map, options, ModificationsDB::global())
}
pub fn write_with_registry(
    mut output: impl Write,
    map: &FeatureMap,
    options: &WriteOptions,
    registry: &ModificationsDB,
) -> Result<()> {
    let bytes = encode(map, options, registry)?;
    output.write_all(&bytes)?;
    Ok(())
}
fn encode(map: &FeatureMap, options: &WriteOptions, registry: &ModificationsDB) -> Result<Vec<u8>> {
    let mut work = Work::new(options.limits);
    work.slots::<Feature>(map.features.len())?;
    work.string(&map.identifier)?;
    work.meta(&map.metadata)?;
    for processing in &map.data_processing {
        work.record()?;
        if !processing.software.cv_terms.is_empty()
            || !processing.software.cv_terms.metadata.is_empty()
        {
            return Err(unsupported("map XML software CV terms and metadata"));
        }
        work.string(&processing.software.name)?;
        work.string(&processing.software.version)?;
        work.meta(&processing.metadata)?;
    }
    let mut unique = BTreeSet::new();
    for feature in &map.features {
        if feature.unique_id != 0 && !unique.insert(feature.unique_id) {
            return Err(bad("duplicate assigned feature ID"));
        }
        work.feature(feature, 0)?;
    }
    for protein in &map.protein_identifications {
        if !protein.protein_groups.is_empty() || !protein.indistinguishable_groups.is_empty() {
            return Err(unsupported(
                "featureXML does not encode structured protein groups; use consensusXML or idXML",
            ));
        }
    }
    xml::measure_identifications(
        &map.protein_identifications,
        &map.unassigned_peptide_identifications,
        &mut work.remaining,
        &mut work.bytes,
    )?;
    let mut peptides = Vec::new();
    work.slots::<&crate::identification::PeptideIdentification>(
        map.unassigned_peptide_identifications.len(),
    )?;
    peptides.extend(map.unassigned_peptide_identifications.iter());
    let mut pending: Vec<_> = map.features.iter().collect();
    while let Some(feature) = pending.pop() {
        xml::measure_identifications(
            &[],
            &feature.peptide_identifications,
            &mut work.remaining,
            &mut work.bytes,
        )?;
        work.slots::<&crate::identification::PeptideIdentification>(
            feature.peptide_identifications.len(),
        )?;
        peptides.extend(feature.peptide_identifications.iter());
        pending.extend(feature.subordinates.iter().rev());
    }
    let definitions = super::modification_definitions::collect_iter_with_budget(
        &map.protein_identifications,
        peptides,
        registry,
        &mut work.remaining,
        &mut work.bytes,
    )?;
    let mut registry = xml::clone_registry(registry, &mut work.remaining, &mut work.bytes)?;
    let context = map_xml::WriteContext::new(&map.protein_identifications)?;
    let mut root = Node::new("featureMap");
    root.attr("version", VERSION);
    root.attr("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    root.attr("xsi:noNamespaceSchemaLocation", "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/FeatureXML_1_9.xsd");
    root.attr("id", format!("fm_{}", map.unique_id));
    if !map.identifier.is_empty() {
        root.attr("document_id", &map.identifier);
    }
    xml::write_meta(&mut root, &map.metadata)?;
    for processing in &map.data_processing {
        root.children.push(map_xml::write_processing(processing)?);
    }
    for protein in &map.protein_identifications {
        xml::measure_search(
            &protein.search_parameters,
            &mut work.remaining,
            &mut work.bytes,
        )?;
        let mut search = protein.search_parameters.clone();
        if let Some(records) = definitions.get(&protein.identifier) {
            super::modification_definitions::attach_with_budget(
                &mut search,
                records,
                &mut work.remaining,
                &mut work.bytes,
            )?;
        }
        super::modification_definitions::register_search_parameters_with_budget(
            &search,
            &mut registry,
            &mut work.remaining,
            &mut work.bytes,
        )?;
        root.children.push(map_xml::write_run(
            protein,
            &context,
            &search,
            &protein.metadata,
        )?);
    }
    for id in &map.unassigned_peptide_identifications {
        root.children.push(map_xml::write_peptide(
            id,
            "UnassignedPeptideIdentification",
            &context,
            &registry,
            &mut work.remaining,
            &mut work.bytes,
        )?);
    }
    let mut list = Node::new("featureList");
    list.attr("count", map.len());
    for feature in &map.features {
        list.children.push(write_feature(
            feature, "f_", &context, &registry, &mut work,
        )?);
    }
    root.children.push(list);
    work.node(&root, 0)?;
    xml::render(
        &root,
        &xml::WriteOptions {
            max_xml_bytes: usize::try_from(options.limits.max_xml_bytes)
                .map_err(|_| bad("XML byte limit overflow"))?,
            max_records: options.limits.max_records,
        },
    )
}
fn text_node(name: &str, value: impl ToString) -> Node {
    let mut n = Node::new(name);
    n.text = value.to_string();
    n
}
fn write_feature(
    feature: &Feature,
    prefix: &str,
    context: &map_xml::WriteContext,
    registry: &ModificationsDB,
    work: &mut Work,
) -> Result<Node> {
    let mut node = Node::new("feature");
    let id = format!("{prefix}{}", feature.unique_id);
    node.attr("id", &id);
    for (dim, value) in [feature.rt, feature.mz].into_iter().enumerate() {
        let mut n = text_node("position", value);
        n.attr("dim", dim);
        node.children.push(n);
    }
    node.children
        .push(text_node("intensity", feature.intensity));
    for (dim, value) in [feature.quality_rt, feature.quality_mz]
        .into_iter()
        .enumerate()
    {
        let mut n = text_node("quality", value);
        n.attr("dim", dim);
        node.children.push(n);
    }
    node.children
        .push(text_node("overallquality", feature.quality));
    node.children.push(text_node("charge", feature.charge));
    for (i, hull) in feature.convex_hulls.iter().enumerate() {
        let mut h = Node::new("convexhull");
        h.attr("nr", i);
        let mut hull = hull.clone();
        hull.compress();
        for point in hull.hull_points() {
            let mut p = Node::new("pt");
            p.attr("x", point.rt);
            p.attr("y", point.mz);
            h.children.push(p);
        }
        node.children.push(h);
    }
    if !feature.subordinates.is_empty() {
        let mut sub = Node::new("subordinate");
        let prefix = format!("{id}_");
        for f in &feature.subordinates {
            sub.children
                .push(write_feature(f, &prefix, context, registry, work)?);
        }
        node.children.push(sub);
    }
    for id in &feature.peptide_identifications {
        node.children.push(map_xml::write_peptide(
            id,
            "PeptideIdentification",
            context,
            registry,
            &mut work.remaining,
            &mut work.bytes,
        )?);
    }
    xml::write_meta(&mut node, &feature.metadata)?;
    Ok(node)
}

pub fn load(path: impl AsRef<Path>) -> Result<FeatureMap> {
    load_with_options(path, &ReadOptions::default())
}
pub fn load_with_options(path: impl AsRef<Path>, options: &ReadOptions) -> Result<FeatureMap> {
    let path = path.as_ref();
    let mut map = read_with_options(path_io::open(path)?, options)?;
    map.loaded_file_path = path
        .to_str()
        .ok_or_else(|| bad("loaded filename must be UTF-8"))?
        .into();
    map.loaded_file_type = FileType::FeatureXml;
    Ok(map)
}
pub fn load_into(
    path: impl AsRef<Path>,
    target: &mut FeatureMap,
    options: &ReadOptions,
) -> Result<()> {
    let draft = load_with_options(path, options)?;
    *target = draft;
    Ok(())
}
pub fn load_size(path: impl AsRef<Path>, options: &ReadOptions) -> Result<usize> {
    read_size(path_io::open(path.as_ref())?, options)
}
pub fn store(path: impl AsRef<Path>, map: &FeatureMap) -> Result<()> {
    store_with_options(path, map, &WriteOptions::default())
}
pub fn store_with_options(
    path: impl AsRef<Path>,
    map: &FeatureMap,
    options: &WriteOptions,
) -> Result<()> {
    let path = path.as_ref();
    if !super::file_types::has_valid_extension(
        path.to_string_lossy().as_ref(),
        FileType::FeatureXml,
    ) {
        return Err(bad("expected featureXML file extension"));
    }
    let bytes = encode(map, options, ModificationsDB::global())?;
    path_io::store(path, &bytes)
}

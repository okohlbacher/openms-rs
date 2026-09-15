// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Native featureXML 1.9 transport, including legacy hulls and nested features.
//! Read filters are half-open and never affect writing. See FEATUREXML_SUPPORT.md.

#[path = "featurexml_scaling.rs"]
mod scaling;
pub use scaling::{Allowance, InputScaling, OutputScaling};

use super::identification_xml::{self as xml, Detach, Node};
use super::{FileType, map_xml, path_io};
use crate::chemistry::ModificationsDB;
use crate::kernel::{ConvexHull2D, Feature, FeatureMap, Point2D};
use crate::metadata::{MetaInfo, MetaValueData};
use crate::{Error, Result};
use std::collections::BTreeSet;
use std::io::{BufRead, Write};
use std::ops::Range;
use std::path::Path;

/// The featureXML schema version this adapter reads and writes.
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

/// Absolute cumulative conversion ceilings. Depth counts subordinate feature
/// levels.
///
/// Every ceiling but `max_xml_bytes` and `max_depth` defaults to unbounded, so
/// the size-derived allowances in [`InputScaling`] and [`OutputScaling`] decide
/// on their own; a caller that sets one keeps it exactly, and the effective
/// ceiling is then the smaller of the two. `max_xml_bytes` is the one quantity
/// nothing can be derived from — it bounds the document itself, on input the
/// decoded bytes accepted from the stream and on output the rendered bytes —
/// so it keeps a large but finite default.
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
            max_xml_bytes: 8 * 1024 * 1024 * 1024,
            max_records: usize::MAX,
            max_list_items: usize::MAX,
            max_depth: 128,
            max_work: usize::MAX,
            max_payload_bytes: usize::MAX,
        }
    }
}
impl Limits {
    /// The fixed ceilings this adapter used before they became size-derived:
    /// 64 MiB of XML, one million elements and list items, 128 subordinate
    /// levels, 50 million work units and 256 MiB of payload.
    ///
    /// Combined with [`InputScaling::fixed`] or [`OutputScaling::fixed`] this
    /// reproduces the former behaviour exactly.
    #[must_use]
    pub const fn former() -> Self {
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
/// Everything one featureXML read is parameterised by: the source
/// `FeatureFileOptions`, the absolute ceilings, and their growth with the size
/// of the document.
#[derive(Clone, Debug, Default)]
pub struct ReadOptions {
    pub feature_options: FeatureFileOptions,
    pub limits: Limits,
    /// Growth of the cumulative ceilings with the decoded document size.
    pub scaling: InputScaling,
}
/// Everything one featureXML write is parameterised by: the absolute ceilings
/// and their growth with the size of the map.
#[derive(Clone, Copy, Debug, Default)]
pub struct WriteOptions {
    pub limits: Limits,
    /// Growth of the cumulative ceilings with the counted size of the map.
    pub scaling: OutputScaling,
}

/// The cumulative ceilings one read or write runs under, after the size-derived
/// allowances have been reconciled with the absolute [`Limits`].
#[derive(Clone, Copy, Debug)]
struct Ceilings {
    work: usize,
    payload_bytes: usize,
    records: usize,
    list_items: usize,
    depth: usize,
}
impl Ceilings {
    /// Reader ceilings earned by a document of `consumed` decoded bytes.
    fn read(limits: Limits, scaling: InputScaling, consumed: usize) -> Self {
        Self {
            work: scaling.work.capped(consumed, limits.max_work),
            payload_bytes: scaling
                .payload_bytes
                .capped(consumed, limits.max_payload_bytes),
            records: scaling.records.capped(consumed, limits.max_records),
            list_items: scaling.list_items.capped(consumed, limits.max_list_items),
            depth: limits.max_depth,
        }
    }
    /// Writer ceilings earned by a map of `units` counted parts.
    fn write(limits: Limits, scaling: OutputScaling, units: usize) -> Self {
        Self {
            work: scaling.work.capped(units, limits.max_work),
            payload_bytes: scaling
                .payload_bytes
                .capped(units, limits.max_payload_bytes),
            records: scaling.records.capped(units, limits.max_records),
            list_items: scaling.list_items.capped(units, limits.max_list_items),
            depth: limits.max_depth,
        }
    }
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
fn xml_options(limits: Limits, ceilings: Ceilings) -> xml::ReadOptions {
    xml::ReadOptions {
        max_xml_bytes: limits.max_xml_bytes,
        max_records: ceilings.records,
        max_list_items: ceilings.list_items,
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
    ceilings: Ceilings,
}
impl Work {
    fn new(ceilings: Ceilings) -> Self {
        Self {
            remaining: ceilings.work,
            bytes: ceilings.payload_bytes,
            records: 0,
            ceilings,
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
        if self.records > self.ceilings.records {
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
        if depth > 2 * self.ceilings.depth.min(Feature::MAX_SUBORDINATE_DEPTH) + 16 {
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
            if count > self.ceilings.list_items {
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
        if depth > self.ceilings.depth.min(Feature::MAX_SUBORDINATE_DEPTH) {
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

/// Read a featureXML document with the default options.
///
/// # Errors
///
/// Malformed or unsupported XML, a field this dialect cannot represent, or an
/// exceeded ceiling; see [`Limits`] and [`InputScaling`].
pub fn read(input: impl BufRead) -> Result<FeatureMap> {
    read_with_options(input, &ReadOptions::default())
}
/// Read a featureXML document, filtering and bounding it as `options` says.
///
/// # Errors
///
/// As [`read`].
pub fn read_with_options(input: impl BufRead, options: &ReadOptions) -> Result<FeatureMap> {
    read_with_registry(input, options, ModificationsDB::global())
}
/// Read a featureXML document, resolving modifications against `registry`
/// rather than the global database, which is left untouched.
///
/// # Errors
///
/// As [`read`], plus chemistry that `registry` cannot resolve.
pub fn read_with_registry(
    input: impl BufRead,
    options: &ReadOptions,
    registry: &ModificationsDB,
) -> Result<FeatureMap> {
    Ok(read_document(input, options, registry, false)?.0)
}
/// Read into `target`, replacing it only on success.
///
/// Atomic replacement: neither parse nor conversion failure changes the target.
///
/// # Errors
///
/// As [`read`].
pub fn read_into(
    input: impl BufRead,
    target: &mut FeatureMap,
    options: &ReadOptions,
) -> Result<()> {
    let draft = read_with_options(input, options)?;
    *target = draft;
    Ok(())
}
/// Return the declared `featureList/@count` without interpreting the features.
///
/// The document is read only through the opening `featureList` tag, so feature
/// payload is never decoded.
///
/// # Errors
///
/// As [`read`], for the prefix that is read.
pub fn read_size(input: impl BufRead, options: &ReadOptions) -> Result<usize> {
    Ok(read_document(input, options, ModificationsDB::global(), true)?.1)
}

/// Attributes and children source accepts on the `featureMap` root.
const ROOT_ATTRS: &[&str] = &[
    "version",
    "document_id",
    "id",
    "unique_id",
    "xmlns:xsi",
    "xsi:noNamespaceSchemaLocation",
];
const ROOT_CHILDREN: &[&str] = &[
    "UserParam",
    "userParam",
    "dataProcessing",
    "IdentificationRun",
    "UnassignedPeptideIdentification",
    "featureList",
    "description",
];

/// What the document header contributes to every feature it precedes.
struct Header {
    context: map_xml::ReadContext,
    registry: ModificationsDB,
}

/// Convert the `featureMap` header into `map` and return the identification
/// context its features are read against.
///
/// `root` carries the children that closed before the header was needed. In a
/// schema-valid document that is all of them but `featureList`, which the
/// FeatureXML 1.9 sequence places last.
fn read_header(
    root: &Node,
    map: &mut FeatureMap,
    options: &xml::ReadOptions,
    source: &ModificationsDB,
    work: &mut Work,
) -> Result<Header> {
    if root.name != "featureMap" {
        return Err(bad("expected featureMap root"));
    }
    check(root, ROOT_ATTRS, ROOT_CHILDREN, false)?;
    xml::measure_node(root, &mut work.remaining, &mut work.bytes)?;
    map.identifier = root.optional("document_id").unwrap_or("").into();
    if let Some(id) = root.optional("id") {
        map.unique_id = uid(id)?;
    }
    if let Some(id) = root.optional("unique_id") {
        map.unique_id = uid(id)?;
    }
    let mut registry = xml::clone_registry(source, &mut work.remaining, &mut work.bytes)?;
    let mut context = map_xml::ReadContext::default();
    for node in &root.children {
        match node.name.as_str() {
            "IdentificationRun" => map.protein_identifications.push(map_xml::read_run(
                node,
                options,
                &mut context,
                &mut registry,
                &mut work.remaining,
                &mut work.bytes,
            )?),
            "dataProcessing" => map
                .data_processing
                .push(map_xml::read_processing(node, options)?),
            _ => {}
        }
    }
    map.metadata = metadata(root, options)?;
    for node in &root.children {
        if node.name == "UnassignedPeptideIdentification" {
            map.unassigned_peptide_identifications
                .push(map_xml::read_peptide(
                    node,
                    options,
                    &context,
                    &registry,
                    &mut work.remaining,
                    &mut work.bytes,
                )?);
        }
    }
    Ok(Header { context, registry })
}

fn push_feature(
    map: &mut FeatureMap,
    mut feature: Feature,
    options: &FeatureFileOptions,
) -> Result<()> {
    if accepts(&feature, options)? {
        if let Some(width) = feature.metadata.get("FWHM") {
            let width = width.as_f64()? as f32;
            feature.set_width(width)?;
        }
        map.features.push(feature);
    }
    Ok(())
}

/// Converts `feature` elements as the parser hands them over, so that the tree
/// in memory is one feature rather than the whole `featureList`.
struct Streamer<'a> {
    options: &'a ReadOptions,
    xml: &'a xml::ReadOptions,
    source: &'a ModificationsDB,
    ceilings: Ceilings,
    records: usize,
    map: FeatureMap,
    header: Option<Header>,
}
impl Streamer<'_> {
    /// Charge one detached element against the shared budgets the parser holds.
    fn take(
        &mut self,
        root: &Node,
        node: Node,
        remaining: &mut usize,
        bytes: &mut usize,
    ) -> Result<()> {
        let mut work = Work {
            remaining: *remaining,
            bytes: *bytes,
            records: self.records,
            ceilings: self.ceilings,
        };
        let outcome = self.convert(root, &node, &mut work);
        *remaining = work.remaining;
        *bytes = work.bytes;
        self.records = work.records;
        outcome
    }
    fn convert(&mut self, root: &Node, node: &Node, work: &mut Work) -> Result<()> {
        if self.header.is_none() {
            self.header = Some(read_header(
                root,
                &mut self.map,
                self.xml,
                self.source,
                work,
            )?);
        }
        let Self {
            options,
            xml,
            map,
            header,
            ..
        } = self;
        let header = header
            .as_ref()
            .ok_or_else(|| bad("featureXML header unavailable"))?;
        work.slots::<Feature>(1)?;
        let feature = read_feature(
            node,
            options,
            xml,
            &header.context,
            &header.registry,
            work,
            0,
        )?;
        push_feature(map, feature, &options.feature_options)
    }
}

fn read_document(
    input: impl BufRead,
    options: &ReadOptions,
    registry: &ModificationsDB,
    size_only: bool,
) -> Result<(FeatureMap, usize)> {
    let limits = options.limits;
    // Ceilings that can admit nothing are refused before the input is touched.
    // They are the shared parser's, so they keep its message.
    if limits.max_records == 0 || limits.max_list_items == 0 {
        return Err(xml::bad("invalid identification XML limits"));
    }
    let prefix_only = options.feature_options.metadata_only || size_only;
    // The one ceiling nothing can be derived from bounds the decode; every
    // other ceiling is then earned by the bytes the decode actually produced.
    let byte_cap = usize::try_from(limits.max_xml_bytes).unwrap_or(usize::MAX);
    let text = xml::decode_document(input, byte_cap, prefix_only.then_some("featureList"))?;
    // Only `featureList` children are streamed, and the prefix reader stops at
    // its opening tag, so a document without one would be held whole either
    // way. The opening tag's spelling is exact in XML, which makes this a cheap
    // necessary condition to refuse on before any tree is built. A document
    // that has the tag only inside a comment still fails the same way below,
    // after the parse.
    if !text.contains("<featureList") {
        return Err(bad("featureMap requires featureList"));
    }
    let ceilings = Ceilings::read(limits, options.scaling, text.len());
    let opts = xml_options(limits, ceilings);
    let depth = 2 * ceilings.depth.min(Feature::MAX_SUBORDINATE_DEPTH) + 16;
    let mut remaining = ceilings.work;
    let mut bytes = ceilings.payload_bytes;
    let mut sink = Streamer {
        options,
        xml: &opts,
        source: registry,
        ceilings,
        records: 0,
        map: FeatureMap::default(),
        header: None,
    };
    let root = if prefix_only {
        xml::parse_text_with_budget(
            &text,
            &opts,
            depth,
            Some("featureList"),
            &mut remaining,
            &mut bytes,
            None,
        )?
    } else {
        let mut take = |root: &Node, node: Node, work: &mut usize, payload: &mut usize| {
            sink.take(root, node, work, payload)
        };
        xml::parse_text_with_budget(
            &text,
            &opts,
            depth,
            None,
            &mut remaining,
            &mut bytes,
            Some(Detach {
                container: "featureList",
                element: "feature",
                take: &mut take,
            }),
        )?
    };
    drop(text);
    let Streamer {
        records,
        mut map,
        header,
        ..
    } = sink;
    let mut work = Work {
        remaining,
        bytes,
        records,
        ceilings,
    };
    if header.is_none() {
        read_header(&root, &mut map, &opts, registry, &mut work)?;
    } else {
        check(&root, ROOT_ATTRS, ROOT_CHILDREN, false)?;
    }
    let mut count = 0;
    let mut seen_list = false;
    for (index, node) in root.children.iter().enumerate() {
        if node.name != "featureList" {
            continue;
        }
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
        // Detaching left only what is not a feature, which must be nothing.
        check(node, &["count"], &["feature"], false)?;
        if index + 1 != root.children.len() {
            return Err(unsupported("featureMap content after featureList"));
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

/// Write `map` as featureXML with the default options.
///
/// The complete document is prepared and validated before any byte reaches
/// `output`.
///
/// # Errors
///
/// A field this dialect cannot represent, a duplicate assigned feature ID, an
/// exceeded ceiling, or the stream's own I/O failure.
pub fn write(output: impl Write, map: &FeatureMap) -> Result<()> {
    write_with_options(output, map, &WriteOptions::default())
}
/// Write `map` as featureXML, bounded as `options` says.
///
/// # Errors
///
/// As [`write()`].
pub fn write_with_options(
    output: impl Write,
    map: &FeatureMap,
    options: &WriteOptions,
) -> Result<()> {
    write_with_registry(output, map, options, ModificationsDB::global())
}
/// Write `map` as featureXML, taking modification definitions from `registry`
/// rather than the global database.
///
/// # Errors
///
/// As [`write()`], plus chemistry `registry` cannot describe portably.
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
/// The counted size of one feature and its subtree: the parts the writer
/// charges for.
///
/// Hull point counts come from [`ConvexHull2D::point_count_bound`], so counting
/// is linear in the number of hulls rather than of points. Recursion stops at
/// the same subordinate depth writing does, so a cyclic or absurdly deep map
/// cannot make counting unbounded; such a map is refused by the writer itself.
fn feature_units(feature: &Feature, depth: usize) -> usize {
    let mut units = 1usize
        .saturating_add(feature.metadata.len())
        .saturating_add(feature.peptide_identifications.len());
    for id in &feature.peptide_identifications {
        units = units.saturating_add(id.hits.len());
    }
    for hull in &feature.convex_hulls {
        units = units.saturating_add(hull.point_count_bound());
    }
    if depth < Feature::MAX_SUBORDINATE_DEPTH {
        for child in &feature.subordinates {
            units = units.saturating_add(feature_units(child, depth + 1));
        }
    }
    units
}

/// The counted size of a whole map, which earns the writer's ceilings.
fn map_units(map: &FeatureMap) -> usize {
    let mut units = map
        .data_processing
        .len()
        .saturating_add(map.metadata.len())
        .saturating_add(map.protein_identifications.len())
        .saturating_add(map.unassigned_peptide_identifications.len());
    for protein in &map.protein_identifications {
        units = units.saturating_add(protein.hits.len());
    }
    for id in &map.unassigned_peptide_identifications {
        units = units.saturating_add(id.hits.len());
    }
    for feature in &map.features {
        units = units.saturating_add(feature_units(feature, 0));
    }
    units
}

fn encode(map: &FeatureMap, options: &WriteOptions, registry: &ModificationsDB) -> Result<Vec<u8>> {
    let ceilings = Ceilings::write(options.limits, options.scaling, map_units(map));
    let mut work = Work::new(ceilings);
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
            // A ceiling wider than the address space is the address space.
            max_xml_bytes: usize::try_from(options.limits.max_xml_bytes).unwrap_or(usize::MAX),
            max_records: ceilings.records,
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

/// Load a featureXML file, recording its path and type on the returned map.
///
/// Plain, gzip and bzip2 input are detected by content.
///
/// # Errors
///
/// The file's own I/O failure, or any error of [`read`].
pub fn load(path: impl AsRef<Path>) -> Result<FeatureMap> {
    load_with_options(path, &ReadOptions::default())
}
/// Load a featureXML file, filtering and bounding it as `options` says.
///
/// # Errors
///
/// As [`load`].
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
/// Load into `target`, replacing it only on success.
///
/// # Errors
///
/// As [`load`].
pub fn load_into(
    path: impl AsRef<Path>,
    target: &mut FeatureMap,
    options: &ReadOptions,
) -> Result<()> {
    let draft = load_with_options(path, options)?;
    *target = draft;
    Ok(())
}
/// Return a file's declared `featureList/@count` without reading its features.
///
/// # Errors
///
/// As [`load`], for the prefix that is read.
pub fn load_size(path: impl AsRef<Path>, options: &ReadOptions) -> Result<usize> {
    read_size(path_io::open(path.as_ref())?, options)
}
/// Store `map` at `path`, replacing the destination atomically.
///
/// The extension must be a featureXML one; `.gz` and `.bz2` suffixes select
/// output compression.
///
/// # Errors
///
/// An unexpected extension, the file's own I/O failure, or any error of
/// [`write()`].
pub fn store(path: impl AsRef<Path>, map: &FeatureMap) -> Result<()> {
    store_with_options(path, map, &WriteOptions::default())
}
/// Store `map` at `path`, bounded as `options` says.
///
/// # Errors
///
/// As [`store`].
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
